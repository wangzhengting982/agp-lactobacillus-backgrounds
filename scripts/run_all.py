"""Execute the complete current-submission reproduction, stopping at first failure."""
from pathlib import Path
import os,sys,subprocess,json,time,uuid,platform
ROOT=Path(__file__).resolve().parents[1]
env=os.environ.copy();env.update(PYTHONUTF8='1',PYTHONIOENCODING='utf-8')
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
stages=[
    'scripts/verify_inputs.py','scripts/run_preparation.py','scripts/run_core.py',
    'scripts/run_deepening.py','src/deepening/validate_deepening.py',
    'scripts/run_sensitivity.py','scripts/validate_sensitivity.py',
    'scripts/run_figures.py','scripts/validate_figure_pixels.py',
    'src/deepening/check_manuscript_claims.py','scripts/summarize_validation.py']
def main():
    assert sys.version_info[:2]==(3,12),'Use the tested Python 3.12 environment.'
    assert sys.prefix!=sys.base_prefix,'Use the independent project .venv.'
    record={'run_id':time.strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8],
        'interpreter':sys.executable,'python':sys.version,'platform':platform.platform(),
        'started_local':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'steps':[],'success':False}
    logs=ROOT/'reports/full_runs'/record['run_id'];logs.mkdir(parents=True,exist_ok=True)
    start=time.time()
    for script in stages:
        begin=time.time();log=logs/(Path(script).stem+'.log');print('START',script,flush=True)
        with log.open('w',encoding='utf-8') as output:
            process=subprocess.Popen([sys.executable,'-u',str(ROOT/script)],cwd=ROOT,env=env,stdout=output,stderr=subprocess.STDOUT)
            code=process.wait()
        record['steps'].append({'script':script,'pid':process.pid,'seconds':round(time.time()-begin,3),'exit_code':code,'log':str(log.relative_to(ROOT))})
        record['elapsed_seconds']=round(time.time()-start,3)
        (ROOT/'reports/full_run_status.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
        print('DONE',script,'exit',code,flush=True)
        if code:
            print(log.read_text(encoding='utf-8')[-7000:],flush=True)
            raise SystemExit(code)
    record.update(success=True,completed_local=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    value=json.dumps(record,ensure_ascii=False,indent=2)
    (ROOT/'reports/full_run_status.json').write_text(value,encoding='utf-8')
    (logs/'run_status.json').write_text(value,encoding='utf-8')
    print('COMPLETE: all current-submission reproduction checks passed.',flush=True)
if __name__=='__main__':main()
