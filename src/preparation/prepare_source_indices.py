"""Rebuild source sample indices/counts without fitting any statistical model.

The default output is an audit staging directory, never the active input folder.
Pass --output-dir to choose another empty destination. Existing frozen input
files cannot be overwritten by this script.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import numpy as np
import pandas as pd
import h5py

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/'05_输入数据'
FROZEN = RAW/'派生输入'
MISSING = {'','nan','none','null','not provided','not applicable',
           'missing: not provided','missing: not applicable','unspecified',
           'unknown','not collected','missing: unknown'}
TAXA = ['Lacticaseibacillus','Lactobacillus','Leuconostoc',
        'Limosilactobacillus','Ligilactobacillus','Pediococcus']
ORDINAL = {'Never':0,'Rarely (a few times/month)':1,
           'Rarely (less than once/week)':1,'Occasionally (1-2 times/week)':2,
           'Regularly (3-5 times/week)':3,'Daily':4}

def clean_subject(s):
    s=s.astype('string').str.strip()
    return s.mask(s.str.lower().isin(MISSING)|s.str.lower().str.startswith('missing:',na=False))

def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def frame_equal(a,b):
    if not a.index.equals(b.index) or not a.columns.equals(b.columns):return False
    for col in a:
        if pd.api.types.is_numeric_dtype(a[col]) and pd.api.types.is_numeric_dtype(b[col]):
            if not np.allclose(a[col],b[col],rtol=1e-12,atol=1e-12,equal_nan=True):return False
        else:
            if not a[col].fillna('<NA>').astype(str).eq(b[col].fillna('<NA>').astype(str)).all():return False
    return True

def prepare(output):
    output=output.resolve()
    if output==FROZEN.resolve() and any((output/n).exists() for n in
        ['first_sample_index.tsv','pilot_analysis_design.tsv','lactobacillales_genus_counts_all_samples.tsv']):
        raise FileExistsError('Frozen inputs exist; rebuild to a separate staging directory and compare instead.')
    output.mkdir(parents=True,exist_ok=True)
    paths = [RAW/'01_原始输入'/n for n in ['4168_annotated_feature_table.tsv','4168_ASV_table.biom',
             'sample_information_10317_matched_4168.tsv','vioscreen_micromacro.tsv']]
    paths += [RAW/'02_模型输入'/n for n in ['analysis_matrix_stool_base.tsv','food_exposure_matrix_stool.tsv']]
    for p in paths:
        if not p.is_file():raise FileNotFoundError(p)
    manifest=[{'file':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha256(p)} for p in paths]
    meta=pd.read_csv(paths[2],sep='\t',dtype=str,keep_default_na=False).set_index('sample_name')
    base=pd.read_csv(paths[4],sep='\t',dtype={'sample_id':str}).set_index('sample_id')
    food=pd.read_csv(paths[5],sep='\t',dtype={'sample_id':str},
                     usecols=['sample_id','yogurt_all_types_except_frozen__freq_year']).set_index('sample_id')
    nut=pd.read_csv(paths[3],sep='\t',dtype={'#SampleID':str},
        usecols=['#SampleID','Energy_in_kcal','Total_Dietary_Fiber_in_g']).set_index('#SampleID')
    assert meta.index.is_unique and base.index.is_unique and food.index.is_unique and nut.index.is_unique
    assert set(base.index)==set(food.index) and set(base.index)<=set(meta.index)
    st=meta.loc[base.index].copy()
    assert st.body_site.eq('UBERON:feces').all()
    st['audit_date']=pd.to_datetime(st.collection_timestamp,errors='coerce',format='mixed')
    st['audit_subject']=clean_subject(st.host_subject_id)
    st['audit_sample']=st.index
    assert st.audit_subject.notna().all() and st.audit_date.notna().all()
    # The full timestamp (including time of day) is compared before sample ID.
    # Equal timestamps use lexicographic sample_id; all three keys are unique.
    st=st.sort_values(['audit_subject','audit_date','audit_sample'],na_position='last',kind='stable')
    assert not st.duplicated(['audit_subject','audit_date','audit_sample']).any()
    first=st.drop_duplicates('audit_subject',keep='first')
    first_index=first[['audit_subject','audit_date']].copy();first_index.index.name='sample_id'
    assert len(st)==3850 and len(first_index)==3204
    ids=first_index.index

    # Full annotated Lactobacillales counts are reconstructed from raw taxonomy.
    tax=pd.read_csv(paths[0],sep='\t',index_col=0)
    lab=tax.loc[tax.index.str.contains(r'(?:^|;)o__Lactobacillales(?:;|$)')]
    family=lab.index.to_series().str.extract(r'(?:^|;)f__([^;]+)',expand=False).fillna('unclassified_family')
    genus=lab.index.to_series().str.extract(r'(?:^|;)g__([^;]+)',expand=False)
    labels=genus.fillna(family+'_unclassified').to_numpy()
    counts=lab.groupby(labels).sum().T
    counts.index.name='sample_id'
    assert counts.shape==(4168,43) and all(t in counts for t in TAXA)
    assert np.equal(tax.to_numpy(),np.floor(tax.to_numpy())).all() and tax.ge(0).all().all()
    # BIOM participates in independent total-count validation; its ASV IDs do
    # not supply a taxonomic map and are not substituted for annotation labels.
    with h5py.File(paths[1],'r') as f:
        bm_ids=[i.decode() if isinstance(i,bytes) else str(i) for i in f['sample/ids'][:]]
        data=f['sample/matrix/data'][:];ptr=f['sample/matrix/indptr'][:]
        bm_total=pd.Series([data[ptr[i]:ptr[i+1]].sum() for i in range(len(bm_ids))],index=bm_ids)
    annotated=tax.sum(axis=0)
    count_check=pd.DataFrame({'annotated_total':annotated,'biom_total':bm_total.loc[annotated.index]})
    count_check['difference']=count_check.annotated_total-count_check.biom_total
    assert count_check.difference.eq(0).all()
    assert np.allclose(base.total_reads,annotated.loc[base.index])

    # Pure preparation section from the earlier pilot; no GLM code is run.
    b=base.loc[ids];f=food.loc[ids];m=meta.loc[ids];n=nut.loc[ids]
    d=pd.DataFrame(index=ids)
    d['age']=b.age_years_clean;d['bmi']=b.bmi_clean;d['sex']=b.sex_clean
    d['probiotic']=m.probiotic_frequency.map(ORDINAL)
    d['plant_ferment']=m.fermented_plant_frequency.map(ORDINAL)
    d['yogurt_freq']=f.iloc[:,0];d['yogurt']=np.log1p(d.yogurt_freq)
    for c in ['antibiotic_recent','ibd_binary','ibs_binary','collection_year_clean']:d[c]=b[c]
    d['country']=np.where(m.country_residence=='United States','US',np.where(m.country_residence=='United Kingdom','UK',np.where(m.country_residence.isin(['not provided','not applicable','missing: not provided','missing: not applicable'])|m.country_residence.eq('')|m.country_residence.isna(),'Unknown','Other')))
    d['energy']=np.log(n.Energy_in_kcal.where(n.Energy_in_kcal>0))
    d['fiber_density']=n.Total_Dietary_Fiber_in_g/n.Energy_in_kcal*1000
    d['depth']=b.log10_total_reads
    d['bowel_quality']=np.select([m.bowel_movement_quality.str.contains('Type 1 and 2',na=False),m.bowel_movement_quality.str.contains('Type 5, 6 and 7',na=False),m.bowel_movement_quality.str.contains('type 3 and 4',case=False,na=False)],['constipated','loose','normal'],default='unknown')
    d['bowel_frequency']=m.bowel_movement_frequency.where(m.bowel_movement_frequency.isin(['Less than one','One','Two','Three','Four','Five or more']),'unknown')
    d=d.replace([np.inf,-np.inf],np.nan)
    age_ok=d.age.ge(18)&d.age.le(100)
    required=['yogurt','probiotic','plant_ferment','antibiotic_recent','ibd_binary','ibs_binary']
    complete=d[required].notna().all(axis=1)
    use=age_ok&complete
    eligible_index=first_index.loc[use].copy()
    eligible_index['total_reads']=b.loc[use,'total_reads']
    eligible_index['in_formal_depth5000']=eligible_index.total_reads.ge(5000)
    eligible_index['in_association_female_male']=eligible_index.in_formal_depth5000 & b.loc[use,'sex_clean'].isin(['female','male'])
    assert len(eligible_index)==2877 and eligible_index.in_formal_depth5000.sum()==2760
    selection=first_index.copy()
    selection['age_18_100']=age_ok
    selection['required_six_fields_complete']=complete
    for c in required:selection[c+'_missing']=d[c].isna()
    selection['in_eligible_2877']=use
    selection['total_reads']=b.total_reads
    selection['in_formal_2760']=use & b.total_reads.ge(5000)
    selection['exclusion_reason']=np.select([~age_ok,~complete],['age_outside_18_to_100_or_missing','required_exposure_or_disease_history_missing'],default='included_eligible')

    # Preserve the historical file schema because agp_core.prepare reads its
    # index. Its standardized values are compatibility artifacts, not formal
    # model covariates. Formal code rebuilds values and recovers year from raw dates.
    d=d.loc[use].copy()
    for c in ['yogurt','probiotic','plant_ferment']:
        d[c]=(d[c]-d[c].mean())/d[c].std(ddof=0)
    for c in ['age','bmi','collection_year_clean','energy','fiber_density','depth']:
        if d[c].isna().any():d[c+'_missing']=d[c].isna().astype(int)
        d[c]=d[c].fillna(d[c].median())
        d[c]=(d[c]-d[c].mean())/d[c].std(ddof=0)

    comparisons=[]
    # Compare frozen references before writing any staged output.
    fp=FROZEN/'first_sample_index.tsv'
    if fp.exists():
        old=pd.read_csv(fp,sep='\t',index_col=0,dtype=str)
        old['audit_date']=pd.to_datetime(old.audit_date,format='mixed')
        ok=frame_equal(first_index,old)
        comparisons.append({'file':fp.name,'comparison':'all IDs in order, participant IDs and full timestamps','agrees':ok})
        assert ok,'Rebuilt first-sample index differs from frozen reference'
    fp=FROZEN/'pilot_analysis_design.tsv'
    if fp.exists():
        old=pd.read_csv(fp,sep='\t',index_col=0,dtype={'sample_id':str})
        ok=frame_equal(d,old)
        comparisons.append({'file':fp.name,'comparison':'full table including index order and standardized legacy values','agrees':ok})
        assert ok,'Rebuilt historical compatibility design differs from frozen reference'
    fp=FROZEN/'lactobacillales_genus_counts_all_samples.tsv'
    if fp.exists():
        old=pd.read_csv(fp,sep='\t',index_col=0,dtype={'sample_id':str})
        ok=frame_equal(counts,old)
        comparisons.append({'file':fp.name,'comparison':'4168 samples × 43 groups, all integer counts and labels','agrees':ok})
        assert ok,'Rebuilt count matrix differs from frozen reference'

    same_day=st.assign(audit_day=st.audit_date.dt.normalize()).groupby(['audit_subject','audit_day']).size()
    same_stamp=st.groupby(['audit_subject','audit_date']).size()
    sd=st.assign(audit_day=st.audit_date.dt.normalize())
    tie_audit=sd.loc[sd.duplicated(['audit_subject','audit_day'],keep=False),
                    ['audit_subject','collection_timestamp','audit_date','audit_sample']].copy()
    tie_audit['chosen_as_first_sample']=tie_audit.index.isin(ids)
    tie_audit['exact_timestamp_tie']=st.duplicated(['audit_subject','audit_date'],keep=False).loc[tie_audit.index]
    norm=st.audit_subject.str.replace(r"^\['([^']+)'\]$",r'\1',regex=True)
    legacy_year=pd.to_numeric(base.loc[ids,'collection_year_clean'],errors='coerce')
    dates=first_index.audit_date.dt.year
    year_audit=pd.DataFrame({'collection_timestamp_original':meta.loc[ids,'collection_timestamp'],
               'legacy_collection_year_clean':legacy_year,'mixed_parsed_year':dates,
               'in_eligible_2877':ids.isin(eligible_index.index)})
    year_audit=year_audit.loc[legacy_year.isna()|legacy_year.ne(dates)]
    summary={'status':'PASS','model_fits':0,'input_files':manifest,
             'stool_samples':len(st),'unique_first_sample_people':len(ids),
             'adult_age_18_100':int(age_ok.sum()),'adult_required_data_complete':int(use.sum()),
             'main_depth_ge5000':int(eligible_index.in_formal_depth5000.sum()),
             'association_female_male':int(eligible_index.in_association_female_male.sum()),
             'complete_mixed_dates':int(st.audit_date.notna().sum()),
             'date_with_slashes':int(st.collection_timestamp.str.contains('/').sum()),
             'same_person_same_day_groups':int(same_day.gt(1).sum()),
             'same_person_exact_timestamp_groups':int(same_stamp.gt(1).sum()),
             'same_person_exact_timestamp_samples':int(same_stamp.loc[same_stamp.gt(1)].sum()),
             'sorting_rule':'subject string ascending, full mixed-parsed timestamp ascending, sample_id lexicographically ascending; first per subject; missing dates last (none observed)',
             'subject_representation_normalization_unique_people':int(norm.nunique()),
             'subject_representation_normalization_changes_number_of_people':bool(norm.nunique()!=st.audit_subject.nunique()),
             'year_recovery_rows':len(year_audit),'year_recovery_eligible_rows':int(year_audit.in_eligible_2877.sum()),
             'lab_taxonomy_rows':len(lab),'lab_genus_or_unclassified_groups':counts.shape[1],
             'annotation_vs_BIOM_all4168_equal':True,'base_depth_matches_annotation':True,
             'frozen_comparisons':comparisons,
             'required_fields_missing_among_first_samples':selection[[c+'_missing' for c in required]].sum().to_dict(),
             'legacy_design_policy':'Rebuilt compatibility table retains old transformed columns, used only for its eligible index by formal code. Formal covariates always derive from raw inputs, including mixed date year.'}
    first_index.to_csv(output/'first_sample_index.tsv',sep='\t')
    st[['audit_subject','audit_date']].to_csv(output/'all_stool_sample_index.tsv',sep='\t',index_label='sample_id')
    counts.to_csv(output/'lactobacillales_genus_counts_all_samples.tsv',sep='\t')
    counts[TAXA].to_csv(output/'six_target_genus_counts_all_samples.tsv',sep='\t')
    d.to_csv(output/'pilot_analysis_design.tsv',sep='\t')
    eligible_index.to_csv(output/'eligible_first_sample_index.tsv',sep='\t')
    selection.to_csv(output/'first_sample_eligibility_audit.tsv',sep='\t')
    tie_audit.to_csv(output/'same_day_and_timestamp_sort_audit.tsv',sep='\t',index_label='sample_id')
    year_audit.to_csv(output/'collection_year_recovery_audit.tsv',sep='\t',index_label='sample_id')
    count_check.to_csv(output/'BIOM_annotation_total_agreement.tsv',sep='\t',index_label='sample_id')
    pd.DataFrame(comparisons).to_csv(output/'frozen_input_comparison.tsv',sep='\t',index=False)
    (output/'source_rebuild_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='input_files'},ensure_ascii=False,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'00_规范与审计/原始输入重建核验')
    args=parser.parse_args()
    prepare(args.output_dir)
