"""Reconstruct the cohort and model inputs from the six supplied source tables."""
from pathlib import Path
import os, sys, shutil, importlib.util, json, time, hashlib
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT/'data/archive'
SOURCE = ARCHIVE/'08_本轮方法复核/provenance/recovered_preparation/new_article/04_复现代码'
STAGE = ROOT/'runs/preparation/project'

def load(name, relative):
    dest = ROOT/'src/preparation'/relative
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.is_file():
        raise FileNotFoundError(f'Missing versioned preparation script: {dest}')
    spec = importlib.util.spec_from_file_location(name, dest)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def compare_frames(new, old):
    import pandas as pd
    import numpy as np
    a = pd.read_csv(new,sep='\t',index_col=0,dtype={'sample_id':str})
    b = pd.read_csv(old,sep='\t',index_col=0,dtype={'sample_id':str})
    assert a.index.equals(b.index) and a.columns.equals(b.columns), str(new)
    error = float(np.nanmax(np.abs(a.to_numpy(float)-b.to_numpy(float))))
    equal = bool(np.allclose(a,b,rtol=1e-10,atol=1e-10,equal_nan=True))
    assert equal, (new,error)
    return {'file':new.name,'shape':list(a.shape),'same_ids_and_columns':True,
            'numeric_equal':equal,'max_abs_error':error,'rtol':1e-10,'atol':1e-10}

def main():
    import pandas as pd
    start = time.time()
    raw = STAGE/'05_输入数据'
    for folder,names in {
        '01_原始输入':['4168_annotated_feature_table.tsv','4168_ASV_table.biom','sample_information_10317_matched_4168.tsv','vioscreen_micromacro.tsv'],
        '02_模型输入':['analysis_matrix_stool_base.tsv','food_exposure_matrix_stool.tsv']
    }.items():
        (raw/folder).mkdir(parents=True, exist_ok=True)
        for name in names:
            dest = raw/folder/name
            if not dest.exists():
                # A local hard link avoids a redundant 200 MB copy; scripts only read inputs.
                try: os.link(ARCHIVE/'02_统计输入/六份源数据'/name,dest)
                except OSError: shutil.copy2(ARCHIVE/'02_统计输入/六份源数据'/name,dest)
            assert hashlib.sha256(dest.read_bytes()).digest()==hashlib.sha256((ARCHIVE/'02_统计输入/六份源数据'/name).read_bytes()).digest(), f'Stale staged input: {name}'
    indices = load('cohort_source','prepare_source_indices.py')
    indices.ROOT, indices.RAW = STAGE, raw
    indices.FROZEN = ARCHIVE/'07_历史与诊断归档/既有扩展_移出展示/原始电子结果/目标属筛选来源'
    indices.prepare(raw/'派生输入')
    prep = load('source_predictors','agp_core.py')
    prep.ROOT, prep.R, prep.PREV = STAGE, raw, raw/'派生输入'
    prep.P = STAGE/'03_分析与结果/agp_models'
    prep.P.mkdir(parents=True,exist_ok=True)
    prep.prepare()
    model = load('association_preparation','association_models/fit_associations.py')
    model.PROJECT, model.SRC, model.RAW = STAGE, prep.P, raw/'01_原始输入'
    model.HERE = STAGE/'03_分析与结果/association_models'
    model.HERE.mkdir(parents=True,exist_ok=True)
    _,_,_,_,summary = model.prepare()
    (model.HERE/'association_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    comparisons = [compare_frames(model.HERE/name,ARCHIVE/'02_统计输入/冻结统计输入'/name)
        for name in ['association_covariate_design.tsv','association_pc_scores.tsv','association_standardized_clr.tsv','association_target_counts.tsv']]
    source_summary = json.loads((raw/'派生输入/source_rebuild_summary.json').read_text(encoding='utf-8'))
    n = source_summary
    flow = {'stool_records':n['stool_samples'],'first_samples':n['unique_first_sample_people'],
        'adult_samples':n['adult_age_18_100'],'complete_samples':n['adult_required_data_complete'],
        'descriptive_samples':n['main_depth_ge5000'],'association_samples':n['association_female_male']}
    flow.update(other_stools=flow['stool_records']-flow['first_samples'],
        age_excluded=flow['first_samples']-flow['adult_samples'],
        missing_excluded=flow['adult_samples']-flow['complete_samples'],
        depth_excluded=flow['complete_samples']-flow['descriptive_samples'],
        sex_excluded=flow['descriptive_samples']-flow['association_samples'])
    original_flow = json.loads((ARCHIVE/'11_投稿格式修订_20261006/图件复绘/source_data/fig1_flow.json').read_text(encoding='utf-8'))
    assert all(original_flow[k]==v for k,v in flow.items())
    # Display field list is metadata; all numerical flow values above are regenerated.
    flow['complete_required_fields'] = original_flow['complete_required_fields']
    cohort = ROOT/'runs/cohort'; cohort.mkdir(parents=True,exist_ok=True)
    (cohort/'fig1_flow.json').write_text(json.dumps(flow,ensure_ascii=False,indent=2),encoding='utf-8')
    eligible = pd.read_csv(raw/'派生输入/eligible_first_sample_index.tsv',sep='\t',index_col=0)
    ids = eligible.index[eligible.in_formal_depth5000]
    counts = pd.read_csv(raw/'派生输入/six_target_genus_counts_all_samples.tsv',sep='\t',index_col=0).loc[ids]
    labels = ['乳酪杆菌属','乳杆菌属','明串珠菌属','黏液乳杆菌属','联合乳杆菌属','片球菌属']
    rows = [dict(taxon=taxon,label=label,threshold=t,detected=int(counts[taxon].ge(t).sum()),n_total=len(ids),percent=100*counts[taxon].ge(t).mean())
        for taxon,label in zip(indices.TAXA,labels) for t in (1,3)]
    detection = pd.DataFrame(rows)
    old = pd.read_csv(ARCHIVE/'11_投稿格式修订_20261006/图件复绘/source_data/fig1_detection.csv')
    pd.testing.assert_frame_equal(detection,old,rtol=1e-10,atol=1e-10)
    detection.to_csv(cohort/'fig1_detection.csv',index=False)
    result={'status':'PASS','elapsed_seconds':time.time()-start,'source_rebuild':source_summary,
        'association_preparation':summary,'matrix_comparisons':comparisons,
        'figure1_flow_and_detection_match':True,'scope':'Rebuilt from six supplied processed source tables; does not reprocess FASTQ or reproduce upstream questionnaire matching.'}
    (ROOT/'reports/preparation_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__': main()
