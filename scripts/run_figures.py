"""Rebuild submitted Figures 1-3 from this checkout's newly calculated results.

No statistical models are fitted by this entry point. It requires outputs of the
statistical runs and records exactly which files are consumed. Figure 1 may use
the frozen submitted flow/count tables only with the explicit opt-in flag;
that limitation is then prominent in the provenance report.
"""
from pathlib import Path
import argparse, subprocess, sys, json, hashlib, shutil, csv

ROOT=Path(__file__).resolve().parents[1]
def absolute(p):return p if p.is_absolute() else ROOT/p
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def relative(p):
    try:return str(p.relative_to(ROOT)).replace('\\','/')
    except ValueError:return str(p)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--core-dir',type=Path,default=Path('runs/core'))
    ap.add_argument('--deepening-dir',type=Path,default=Path('runs/deepening/results'))
    ap.add_argument('--diet-dir',type=Path,default=Path('runs/deepening/diet_preprocessing'))
    ap.add_argument('--cohort-source-dir',type=Path,default=Path('runs/cohort'))
    ap.add_argument('--allow-frozen-figure1',action='store_true')
    ap.add_argument('--font-latin',type=Path)
    ap.add_argument('--font-cjk',type=Path)
    args=ap.parse_args()
    core,deep,diet,cohort=map(absolute,[args.core_dir,args.deepening_dir,args.diet_dir,args.cohort_source_dir])
    out=ROOT/'runs/figures';source=out/'source_data';source.mkdir(parents=True,exist_ok=True)
    sources={
        'C04':core/'表S22c_全部15组主成分向量比较.tsv',
        'C13':core/'表S24i_三背景两目标属效应量.tsv',
        'P02':diet/'更广泛饮食调整_全部属间差异.tsv',
        'S02':deep/'S02_四状态标准化概率.tsv',
        'A01':deep/'A01_检出后相对丰度.tsv'}
    missing=[relative(p) for p in sources.values() if not p.is_file()]
    if missing:raise FileNotFoundError('Statistical rerun outputs missing; refusing silent frozen-data substitution: '+str(missing))
    import pandas as pd
    tables={};provenance=[]
    for name,p in sources.items():
        rows=pd.read_csv(p,sep='\t').to_dict('records')
        for i,row in enumerate(rows):row.update({'_source':relative(p),'_row':i+2})
        tables[name]=rows
        provenance.append({'figure_panels':{'C04':['2a'],'C13':['2b'],'P02':['2c'],'S02':['3a'],'A01':['3b']}[name],'source':relative(p),'sha256':sha(p),'rows':len(rows),'kind':'current_statistical_run_output'})
    (source/'plot_tables.json').write_text(json.dumps(tables,ensure_ascii=False,indent=2),encoding='utf-8')
    frozen=ROOT/'src/figures/source_data'
    used_frozen=False
    for name in ['fig1_flow.json','fig1_detection.csv']:
        p=cohort/name
        kind='current_cohort_reconstruction_output'
        if not p.is_file():
            if not args.allow_frozen_figure1:raise FileNotFoundError(f'{relative(p)} missing. Supply independently recalculated cohort source tables, or explicitly --allow-frozen-figure1 and report that limitation.')
            p=frozen/name;used_frozen=True;kind='frozen_submitted_plot_input_NOT_new_cohort_reproduction'
        shutil.copy2(p,source/name)
        provenance.append({'figure_panels':['1a' if name.endswith('json') else '1b'],'source':relative(p),'sha256':sha(p),'kind':kind})
    report={'entry_point':'scripts/run_figures.py','figures_2_and_3_inputs':'Current statistical-run output tables, not submitted Excel values','figure_1_frozen_inputs_used':used_frozen,'claim_limit':'Rendering and plotted-value comparison only. Model-run validation is separately documented. No upstream sequencing reanalysis is claimed.','sources':provenance,'status':'inputs_prepared'}
    rp=ROOT/'reports/figures_provenance.json';rp.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    command=[sys.executable,str(ROOT/'src/figures/rebuild_figures.py'),'--source-dir',str(source),'--output-dir',str(out)]
    for key in ['font_latin','font_cjk']:
        if getattr(args,key):command+=['--'+key.replace('_','-'),str(absolute(getattr(args,key)))]
    with (out/'render.log').open('w',encoding='utf-8') as log:
        result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode==0:
        with (out/'export_validation.log').open('w',encoding='utf-8') as log:
            result=subprocess.run([sys.executable,str(ROOT/'src/figures/validate_exports.py'),'--output-dir',str(out)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    report['status']='rendered_and_numeric_assertions_passed' if result.returncode==0 else 'render_failed'
    report['returncode']=result.returncode
    if result.returncode==0:
        report['numeric_audit']=relative(out/'numeric_plot_audit.json')
        report['rendered_visual_QA']='pending human/model image inspection; automated bounds checks are not visual QA'
    rp.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    raise SystemExit(result.returncode)

if __name__=='__main__':main()
