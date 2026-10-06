"""Recompute the core model outputs and compare them with submitted results.

Run from any working directory using the project's independent interpreter:
    .venv/Scripts/python.exe scripts/run_core.py
Run scripts/run_preparation.py first. Primary inputs then come from
runs/preparation, taxonomy inputs are rebuilt from data/archive, and all core
output is written to runs/core. Archive model results are comparisons only.
"""
from pathlib import Path
import os
import sys
import subprocess
import json
import time
import argparse
import uuid
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--steps', nargs='*', help='Optional individual script name prefixes (e.g. 04 06). Default runs all steps.')
args=parser.parse_args()
env=os.environ.copy()
env.update(PYTHONUTF8='1',PYTHONUNBUFFERED='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
out=ROOT/'runs/core'; out.mkdir(parents=True,exist_ok=True)
reports=ROOT/'reports';reports.mkdir(exist_ok=True)
scripts=['00_rebuild_taxonomy_inputs.py','01_genus_contrasts.py','02_taxonomy_trace.py','02b_refit_taxonomy_pc.py','03_taxonomy_PC_contrasts.py','04_complete_background_differences.py','05_key_difference_sensitivity.py','06_validate_core.py']
if args.steps:
    scripts=[s for s in scripts if any(s.startswith(prefix) for prefix in args.steps)]
    assert scripts,'No matching steps'
run_id=time.strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8]
record=dict(run_id=run_id,interpreter=sys.executable,started_at=time.strftime('%Y-%m-%dT%H:%M:%S'),threads=1,parent_pid=os.getpid(),steps=[])
record_path=reports/('core_run_status.json' if not args.steps else 'core_selected_run_status.json')
history=reports/'core_run_history';history.mkdir(exist_ok=True)
def save_record():
    value=json.dumps(record,ensure_ascii=False,indent=2)
    record_path.write_text(value,encoding='utf8')
    (history/(run_id+'.json')).write_text(value,encoding='utf8')
for name in scripts:
    start=time.perf_counter()
    print('START',name,flush=True)
    with (out/(name+'.log')).open('w',encoding='utf8') as log:
        process=subprocess.Popen([sys.executable,str(ROOT/'src/core'/name)],cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf8',errors='replace')
        step=dict(script=name,pid=process.pid,started_at=time.strftime('%Y-%m-%dT%H:%M:%S'),status='running')
        record['steps'].append(step)
        save_record()
        for line in process.stdout:
            log.write(line);log.flush();print(line,end='',flush=True)
        code=process.wait()
    step.update(exit_code=code,seconds=round(time.perf_counter()-start,3),status='completed' if code==0 else 'failed')
    record['success']=code==0
    save_record()
    if code:
        raise SystemExit(code)
record['completed_at']=time.strftime('%Y-%m-%dT%H:%M:%S')
save_record()
print('CORE RUN COMPLETE',flush=True)
