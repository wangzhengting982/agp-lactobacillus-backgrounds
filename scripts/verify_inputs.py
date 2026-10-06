"""Fail before fitting if an archived input or submission reference has changed."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def verify():
    entries=json.loads((ROOT/'reports/input_manifest.json').read_text(encoding='utf-8'))
    unique={e['path']:e for e in entries}
    failed=[]
    for name,e in unique.items():
        p=ROOT/'data/archive'/name
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=e['sha256']:failed.append(name)
    references=json.loads((ROOT/'reference_manifest.json').read_text(encoding='utf-8'))
    for e in references:
        p=ROOT/e['path']
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=e['sha256']:failed.append(e['path'])
    result={'status':'PASS' if not failed else 'FAIL','unique_archive_files_checked':len(unique),'reference_files_checked':len(references),'changed_or_missing':failed}
    (ROOT/'reports/input_integrity.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if failed:raise SystemExit(1)
    return result
if __name__=='__main__':verify()
