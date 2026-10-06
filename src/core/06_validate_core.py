"""Compare newly computed core results with submitted workbook and archive.

This reads old results only as references, after the model scripts have run.
Numerical tolerance: absolute 1e-8 + relative 1e-7. All p/q decisions at 0.05
must match separately; identifiers and non-numerical fields match exactly.
"""
import json
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from core_paths import ROOT, DATA, OUT, PRIMARY, TAXONOMY_INPUT

T=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
read=lambda path:pd.read_csv(path,sep='\t')
# Export original PC omnibus and six PC block tests from freshly fitted arrays.
a=np.load(OUT/'pc_joint_sandwich.npz');b=a['beta'];v=a['covariance']
counts=read(PRIMARY/'association_target_counts.tsv')
blocks=[]
for i,taxon in enumerate(T):
    ix=slice(i*30+20,(i+1)*30);beta=b[ix];cov=v[ix,ix]
    rank=int(np.linalg.matrix_rank(cov));stat=float(beta@np.linalg.pinv(cov)@beta)
    blocks.append(dict(taxon=taxon,n=len(counts),events=int(counts[taxon].ge(1).sum()),chi2=stat,df=rank,p=float(chi2.sf(stat,rank)),valid=True))
blocks=pd.DataFrame(blocks);blocks['q_BH_six']=multipletests(blocks.p,method='fdr_bh')[1]
blocks.to_csv(OUT/'refitted_primary_pc_blocks.tsv',sep='\t',index=False)
c=np.zeros((50,180))
for i in range(1,6):
    c[(i-1)*10:i*10,20:30]=-np.eye(10)
    c[(i-1)*10:i*10,i*30+20:(i+1)*30]=np.eye(10)
d=c@b;cv=c@v@c.T;rank=int(np.linalg.matrix_rank(cv));stat=float(d@np.linalg.pinv(cv)@d)
omnibus=pd.DataFrame([dict(comparison='all six target genera have equal 10-PC log-odds coefficient vectors',chi2=stat,df=rank,nominal_constraints=50,p=float(chi2.sf(stat,rank)),n_subjects=len(counts),same_subject_correlation='joint sandwich using six GLMs individual score contributions',valid=True,scope='explanatory relative-composition association heterogeneity, not prediction improvement or causal heterogeneity')])
omnibus.to_csv(OUT/'refitted_primary_pc_omnibus.tsv',sep='\t',index=False)
tax_omnibus=read(OUT/'refitted_taxonomy_pc_omnibus.tsv')
for release in ['R07','R06']:
    p=TAXONOMY_INPUT/release
    z=pd.read_csv(p/'standardized_clr.tsv.gz',sep='\t',index_col=0)
    pca=PCA(n_components=10,svd_solver='full').fit(z)
    mask=tax_omnibus.release.eq(release)
    tax_omnibus.loc[mask,'n_selected_genera']=z.shape[1]
    tax_omnibus.loc[mask,'PC_variance_explained']=float(pca.explained_variance_ratio_.sum())
    # Reconstruct the full named background genus universe from the archived
    # all-ASV assignments, using aggregate_full_asv_sintax.py's original rule.
    assignments=read(DATA/'输入/03_分析与结果/taxonomy_validation'/(release+'_all_asvs_genus_assignments.tsv')).fillna({'genus':'','family':'','order':'','domain':''})
    bacterial=assignments.domain.eq('Bacteria') & assignments.domain_bootstrap.ge(.8)
    forbidden=(assignments['order'].eq('Lactobacillales') & assignments.order_bootstrap.ge(.8)) | (assignments.family.eq('Bifidobacteriaceae') & assignments.family_bootstrap.ge(.8)) | (assignments.genus.isin(T) & assignments.genus_bootstrap.ge(.8))
    named=assignments.genus.ne('') & ~assignments.genus.str.contains('unclassified|uncultured',case=False,regex=True) & assignments.genus_bootstrap.ge(.8)
    families=assignments.family.where(assignments.family_bootstrap.ge(.8),'unclassified_family')
    labels='microbe::'+families+'|'+assignments.genus
    tax_omnibus.loc[mask,'n_background_genera']=labels[bacterial & ~forbidden & named].nunique()
tax_omnibus.to_csv(OUT/'refitted_taxonomy_pc_omnibus.tsv',sep='\t',index=False)

MAPPING={
 'C01':('refitted_936_associations.tsv',['taxon','feature']),
 'C02':('refitted_primary_pc_blocks.tsv',['taxon']),
 'C03':('refitted_primary_pc_omnibus.tsv',[]),
 'C04':('表S22c_全部15组主成分向量比较.tsv',['taxon_A','taxon_B']),
 'C05':('表S22e_全部156背景属整体差异.tsv',['feature']),
 'C06':('表S22f_全部2340项背景属系数比较.tsv',['feature','taxon_A','taxon_B']),
 'C07':('表S23_整体差异稳健性.tsv',['condition']),
 'C08':('表S24a_重点背景属分类支持与覆盖.tsv',['release','genus','support']),
 'C09':('表S24b_重新注释候选ASV明细.tsv',['release','query_id']),
 'C10':('表S24c_参考库属标签.tsv',['release','reference_genus']),
 'C11':('表S24g_重新注释15组整体模式比较.tsv',['release','taxon_A','taxon_B']),
 'C12':('表S24h_三项突出差异的阈值与分类复核.tsv',['scenario','feature']),
 'C13':('表S24i_三背景两目标属效应量.tsv',['scenario','feature','taxon']),
 'C14':('refitted_taxonomy_pc_blocks.tsv',['release','taxon']),
 'C15':('refitted_taxonomy_pc_omnibus.tsv',['release']),
}
results=[]; details=[]
def compare(actual,expected,label,keys=None,exclude=()):
    actual=actual.copy();expected=expected.copy()
    cols=[col for col in expected if col not in exclude]
    missing=[col for col in cols if col not in actual]
    assert not missing,(label,'missing columns',missing)
    assert len(actual)==len(expected),(label,'row count differs',len(actual),len(expected))
    if keys:
        for frame in [actual,expected]:
            assert not frame.duplicated(keys).any(),(label,'duplicate keys',keys)
        actual=actual.sort_values(keys).reset_index(drop=True)
        expected=expected.sort_values(keys).reset_index(drop=True)
    overall=True; maximum=0.;maximum_relative=0.; cells=0;failures=[]
    for col in cols:
        aa=actual[col];ee=expected[col]
        is_num=pd.api.types.is_numeric_dtype(ee) and not pd.api.types.is_bool_dtype(ee)
        if is_num:
            av=pd.to_numeric(aa).to_numpy(dtype=float);ev=pd.to_numeric(ee).to_numpy(dtype=float)
            ok=np.isclose(av,ev,atol=1e-8,rtol=1e-7,equal_nan=True)
            finite=np.isfinite(av)&np.isfinite(ev)
            diff=np.abs(av[finite]-ev[finite]);reldiff=diff/np.maximum(np.abs(ev[finite]),1e-300)
            err=float(diff.max()) if diff.size else 0.
            relerr=float(reldiff.max()) if reldiff.size else 0.
            maximum=max(maximum,err);maximum_relative=max(maximum_relative,relerr)
            decisions=True
            if col=='p' or col.startswith('q_') or col.startswith('p_'):
                decisions=bool(np.array_equal(av<.05,ev<.05))
            passed=bool(ok.all() and decisions)
            details.append(dict(table=label,column=col,n=len(aa),numeric=True,maximum_absolute_error=err,maximum_relative_error=relerr,significance_decisions_match=decisions,pass_=passed))
        else:
            av=aa.fillna('<MISSING>').astype(str);ev=ee.fillna('<MISSING>').astype(str)
            # Pandas may infer Excel integer IDs and TSV floats differently.
            ok=av.eq(ev).to_numpy();passed=bool(ok.all())
            details.append(dict(table=label,column=col,n=len(aa),numeric=False,maximum_absolute_error=None,maximum_relative_error=None,significance_decisions_match=None,pass_=passed))
        if not passed:
            failures.append(col)
            bad=np.flatnonzero(~ok)[:5]
            print('MISMATCH',label,col,[(int(j),str(aa.iloc[j]),str(ee.iloc[j])) for j in bad],flush=True)
        overall=overall and passed;cells+=len(aa)
    results.append(dict(table=label,rows=len(expected),compared_columns=len(cols),compared_cells=cells,max_absolute_error=maximum,max_relative_error=maximum_relative,pass_=overall,failed_columns=' | '.join(failures),excluded_columns=' | '.join(exclude)))
    return overall

xls=pd.ExcelFile(ROOT/'reference/submitted_results.xlsx')
for prefix,(file,keys) in MAPPING.items():
    sheet=next(s for s in xls.sheet_names if s.startswith(prefix+'_'))
    compare(read(OUT/file),pd.read_excel(xls,sheet),sheet,keys)

for archived in sorted((DATA/'新增分析').glob('*.tsv')):
    current=OUT/archived.name
    if not current.exists():
        continue
    # Both original and migrated scripts preserve row order for these outputs.
    compare(read(current),read(archived),'archive/'+archived.name)

summary=dict(status='pass' if all(row['pass_'] for row in results) else 'fail',
 tolerance=dict(absolute=1e-8,relative=1e-7,significance_threshold=.05),
 scope='Fresh GLM fits from newly rebuilt primary matrices and BIOM-reaggregated taxonomy matrices; taxonomy annotation outputs are supplied inputs, not rerun SINTAX classifications.',
 submission_tables_checked=list(MAPPING),submission_table_count=len(MAPPING),
 excluded_metadata={},
 results=results,
 model_counts=dict(primary_single_background=936,two_background_plus_PC_sensitivity=48,stacked_independent_checks=4,taxonomy_Veillonella=12,taxonomy_PC=12,key_pair_sensitivity=24,total_GLM_fits=1036),
 upstream_scope='Primary inputs are read from runs/preparation/project/03_分析与结果/association_models after source reconstruction. Taxonomy matrices are rebuilt in step 00 from BIOM and supplied labels.')
(ROOT/'reports/core_verification.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
pd.DataFrame(results).to_csv(ROOT/'reports/core_table_comparison.tsv',sep='\t',index=False)
pd.DataFrame(details).to_csv(ROOT/'reports/core_column_comparison.tsv',sep='\t',index=False)
print(pd.DataFrame(results).to_string(index=False),flush=True)
assert summary['status']=='pass','See reports/core_verification.json for mismatches'
print('ALL CORE COMPARISONS PASSED',flush=True)
