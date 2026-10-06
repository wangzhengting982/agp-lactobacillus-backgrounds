"""Fresh fits for current N01/N02 sensitivity tables; source archive is read-only."""
from pathlib import Path
import os, sys, subprocess, shutil, json, time, hashlib
ROOT=Path(__file__).resolve().parents[1]
A=ROOT/'data/archive'; OLD=A/'08_本轮方法复核'
RUN=ROOT/'runs/sensitivity'; VIEW=RUN/'input_view'; SRC=ROOT/'src/sensitivity'
RUN.mkdir(parents=True,exist_ok=True);SRC.mkdir(parents=True,exist_ok=True)
env=os.environ.copy()
env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']})
env.update(PYTHONUTF8='1',PYTHONIOENCODING='utf-8',AGP_AUDIT_SOURCE_PACKAGE=str(VIEW))

def overlay(src,dest):
    dest.parent.mkdir(parents=True,exist_ok=True)
    # View files may be hard links to immutable source input: detach before replacing.
    assert dest.resolve().is_relative_to(RUN.resolve())
    if dest.exists():dest.unlink()
    shutil.copy2(src,dest)

def prepare():
    for folder in ['01_核心统计','02_统计输入','03_饮食与菌群深化']:
        for src in (A/folder).rglob('*'):
            if src.is_file():
                dest=VIEW/src.relative_to(A);dest.parent.mkdir(parents=True,exist_ok=True)
                if not dest.exists():
                    try:os.link(src,dest)
                    except OSError:shutil.copy2(src,dest)
    rebuilt=ROOT/'runs/preparation/project/03_分析与结果/association_models'
    for name in ['association_covariate_design.tsv','association_pc_scores.tsv','association_standardized_clr.tsv','association_target_counts.tsv']:
        overlay(rebuilt/name,VIEW/'02_统计输入/冻结统计输入'/name)
    for src in (ROOT/'runs/deepening/diet_preprocessing').glob('*.tsv'):
        overlay(src,VIEW/'02_统计输入/饮食预处理'/src.name)
    for src in (ROOT/'runs/deepening/results').glob('*.tsv'):
        overlay(src,VIEW/'03_饮食与菌群深化/完整分析结果'/src.name)
    overlay(ROOT/'runs/preparation/project/05_输入数据/派生输入/lactobacillales_genus_counts_all_samples.tsv',VIEW/'02_统计输入/冻结统计输入/lactobacillales_genus_counts_all_samples.tsv')
    # Fresh classification-specific matrices are connected when available; their
    # upstream ASV classifications remain supplied data, not a new SINTAX run.
    taxonomy=ROOT/'runs/core/taxonomy_inputs'
    used_taxonomy=[]
    for label in ['R06','R07']:
        assert (taxonomy/label/'covariate_design.tsv').exists(), 'Run core taxonomy preparation first'
        for src in (taxonomy/label).glob('*.tsv*'):
            overlay(src,VIEW/'01_核心统计/输入/03_分析与结果/taxonomy_validation/association_sensitivity'/label/src.name)
            used_taxonomy.append(str(src.relative_to(ROOT)))
    changes=[]
    for relative in ['diversity/run_diversity_audit.py','provenance/01_intake.py','provenance/02_round_sensitivity.py','bloom/run_bloom_sensitivity.py','bloom/run_qc5000_sensitivity.py','within_person/run_within_audit.py']:
        old=OLD/relative;code=old.read_text(encoding='utf-8-sig');new=code
        if 'runtime_paths=os.environ' in new:
            start=new.index('runtime_paths=os.environ');end=new.index('\nimport numpy',start)
            new=new[:start]+'# Use only this independent Python environment.\n'+new[end:]
        new='\n'.join(line for line in new.split('\n') if not line.startswith('sys.path[:0]='))
        if relative.startswith('provenance/'):
            new=new.replace('P=Path(__file__).resolve().parent;',"P=Path(os.environ['AGP_AUDIT_OUTPUT']);P.mkdir(parents=True,exist_ok=True);")
        dest=SRC/relative;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(new,encoding='utf-8')
        changes.append({'file':relative,'source_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),
            'changes':'Runtime path injection removed; provenance output root configurable. Numerical formulas unchanged.'})
    for relative in ['diversity/fixed_plan.md','bloom/fixed_plan.txt','bloom/qc5000_fixed_plan.txt']:
        overlay(OLD/relative,RUN/relative)
    for name in ['AGP_Table_S1_official.xlsx','AGP_Table_S1_sample_round_long.tsv','AGP_BLOOM.fasta']:
        overlay(OLD/'provenance/sources'/name,RUN/'provenance/sources'/name)
    (ROOT/'reports/sensitivity_portability_changes.json').write_text(json.dumps({'changes':changes,'taxonomy_rebuilt_inputs':used_taxonomy},ensure_ascii=False,indent=2),encoding='utf-8')

def main():
    prepare();records=[]
    stages=[
        ('diversity','diversity/run_diversity_audit.py',['--source-root',VIEW,'--output-dir',RUN/'diversity']),
        ('intake','provenance/01_intake.py',[]),
        ('round','provenance/02_round_sensitivity.py',['--pooled']),
        ('bloom','bloom/run_bloom_sensitivity.py',['--source-root',VIEW,'--provenance-root',RUN/'provenance','--output-dir',RUN/'bloom']),
        ('qc5000','bloom/run_qc5000_sensitivity.py',['--source-root',VIEW,'--bloom-root',RUN/'bloom']),
        ('within_person','within_person/run_within_audit.py',['--source',VIEW,'--output',RUN/'within_person'])]
    env['AGP_AUDIT_OUTPUT']=str(RUN/'provenance')
    for name,script,args in stages:
        start=time.time();log=RUN/(name+'.log');print('START',name,flush=True)
        with log.open('w',encoding='utf-8') as out:
            process=subprocess.Popen([sys.executable,'-u',str(SRC/script),*map(str,args)],cwd=ROOT,env=env,stdout=out,stderr=subprocess.STDOUT)
            code=process.wait()
        records.append({'stage':name,'pid':process.pid,'seconds':time.time()-start,'exit_code':code,'log':str(log.relative_to(ROOT))})
        (ROOT/'reports/sensitivity_run_status.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
        print('DONE',name,code,flush=True)
        if code:
            print(log.read_text(encoding='utf-8')[-7000:]);sys.exit(code)

if __name__=='__main__':main()
