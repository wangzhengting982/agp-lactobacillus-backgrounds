"""Re-fit dietary and extension analyses from archived input tables.

Usage: python scripts/run_deepening.py [--from-stage NAME] [--through-stage NAME]
All numerical libraries run with one thread. Logs and results are new files in
runs/deepening. Archived PC scores/results are never used as fitted outputs.
"""
from pathlib import Path
import os,sys,subprocess,time,json,argparse,shutil
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'src/deepening'
RUN=Path(os.environ.get('AGP_DEEPENING_OUTPUT',ROOT/'runs/deepening'))
RUN.mkdir(parents=True,exist_ok=True)
stages=['00a_inventory','00b_diet_pca','00c_metadata','01_prepare','02_diet_interactions','03_joint_states_abundance','04_within_person','05_asv','06_core_robustness','07_verify_statistics','08_BMI_strata_resolution','09_sex_contrasts','10_abundance_influence']
p=argparse.ArgumentParser();p.add_argument('--from-stage',choices=stages);p.add_argument('--through-stage',choices=stages);args=p.parse_args()
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']});env.update(PYTHONUTF8='1',PYTHONIOENCODING='utf-8')
chosen=stages[stages.index(args.from_stage) if args.from_stage else 0:(stages.index(args.through_stage)+1) if args.through_stage else len(stages)]
records=[]
for stage in chosen:
    t=time.time();log=RUN/(stage+'.log')
    print('START',stage,flush=True)
    with log.open('w',encoding='utf-8') as out:
        proc=subprocess.run([sys.executable,'-u',str(SRC/(stage+'.py'))],cwd=ROOT,env=env,stdout=out,stderr=subprocess.STDOUT)
    records.append(dict(stage=stage,seconds=round(time.time()-t,3),exit_code=proc.returncode,log=str(log.relative_to(ROOT))))
    (ROOT/'reports/deepening_run_status.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
    print('DONE',stage,'seconds',records[-1]['seconds'],'exit',proc.returncode,flush=True)
    if proc.returncode:
        print(log.read_text(encoding='utf-8')[-6000:],flush=True);sys.exit(proc.returncode)
    if stage=='00b_diet_pca':
        out=RUN/'results';out.mkdir(exist_ok=True)
        for old,new in [('更广泛饮食调整_全部六配对模型.tsv','P01_广泛饮食调整全部关联.tsv'),('更广泛饮食调整_全部属间差异.tsv','P02_广泛饮食调整全部差异.tsv')]:shutil.copyfile(RUN/'diet_preprocessing'/old,out/new)
print('All selected analyses completed.',flush=True)
