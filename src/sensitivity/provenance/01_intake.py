from pathlib import Path
import sys, os, json, hashlib, urllib.request, datetime
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
# Use only this independent Python environment.

import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix
P=Path(os.environ['AGP_AUDIT_OUTPUT']);P.mkdir(parents=True,exist_ok=True); S=Path(os.environ.get('AGP_AUDIT_SOURCE_PACKAGE',str(P.parent/'source_package'))); (P/'sources').mkdir(exist_ok=True)
def write(name,o): (P/name).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=str),encoding='utf8')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fetch(url,name):
    req=urllib.request.Request(url,headers={'User-Agent':'AGP-provenance-research-audit/1.0'})
    with urllib.request.urlopen(req,timeout=40) as r:b=r.read()
    p=P/'sources'/name;p.write_bytes(b);return b
read=lambda p:pd.read_csv(p,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
X=read(S/'02_统计输入/冻结统计输入/association_covariate_design.tsv'); C=read(S/'02_统计输入/冻结统计输入/association_target_counts.tsv').loc[X.index]
official=P/'sources'/'AGP_Table_S1_official.xlsx'
rows=[];workbook=pd.ExcelFile(official)
for sheet in workbook.sheet_names:
    if not sheet.startswith('AG') or sheet=='AGP sample summary': continue
    f=pd.read_excel(workbook,sheet_name=sheet,dtype='string')
    # Same published columns used by the previous explicit Table S1 mapping.
    if '#SampleID' not in f.columns: print('SHEET_COLUMNS',sheet,list(f.columns)); continue
    for sample in f['#SampleID'].dropna(): rows.append({'sample_id':str(sample),'sequencing_round':sheet})
new=pd.DataFrame(rows)
old=pd.read_csv(official.parent/'AGP_Table_S1_sample_round_long.tsv',sep='\t',dtype=str)
if len(new)==0:
    write('table_s1_columns.json',{sheet:list(pd.read_excel(official,sheet_name=sheet,nrows=0).columns) for sheet in workbook.sheet_names if sheet.startswith('AG')})
    raise RuntimeError('Need exact column audit; no guessed sample column.')
new=new.drop_duplicates().sort_values(['sample_id','sequencing_round']); oldpair=old[['sample_id','sequencing_round']].drop_duplicates().sort_values(['sample_id','sequencing_round'])
assert set(map(tuple,new.to_numpy()))==set(map(tuple,oldpair.to_numpy())), 'Mapping does not match existing audited Table S1 long data'
new.to_csv(P/'round_mapping_all.tsv',sep='\t',index=False)
g=new.groupby('sample_id').sequencing_round.agg(list); m=pd.DataFrame(index=X.index)
m['n_rounds']=m.index.map(g.map(len)).fillna(0).astype(int)
m['rounds']=m.index.map(g.map(lambda v:'|'.join(v))).fillna('')
m['mapping_status']=np.select([m.n_rounds.eq(1),m.n_rounds.gt(1)],['unique','multiple'],default='unmapped')
for t in ['Lacticaseibacillus','Lactobacillus']:m[t+'_detected']=C[t].ge(1).astype(int)
m.to_csv(P/'round_mapping_primary2748.tsv',sep='\t')
summary=m.groupby('mapping_status').agg(n=('n_rounds','size'),Lacticaseibacillus_events=('Lacticaseibacillus_detected','sum'),Lactobacillus_events=('Lactobacillus_detected','sum'))
summary.to_csv(P/'round_mapping_coverage.tsv',sep='\t')
stats=m.loc[m.n_rounds.eq(1)].groupby('rounds').agg(n=('n_rounds','size'),Lacticaseibacillus_events=('Lacticaseibacillus_detected','sum'),Lactobacillus_events=('Lactobacillus_detected','sum'))
for t in ['Lacticaseibacillus','Lactobacillus']:stats[t+'_non_events']=stats.n-stats[t+'_events']
stats.to_csv(P/'round_strata_diagnostics.tsv',sep='\t')
write('round_mapping_audit.json',{'official_file':str(official),'sha256':sha(official),'workbook_sheets':workbook.sheet_names,'derived_long_exact_pairs_match':True,'primary_n':len(m),'unique_n':int(m.n_rounds.eq(1).sum()),'multiple_n':int(m.n_rounds.gt(1).sum()),'unmapped_n':int(m.n_rounds.eq(0).sum()),'unique_rounds':len(stats),'scope':'Published sequencing round proxy only; not exact Qiita prep/artifact or current ASV generating run.'})
print(summary.to_string());print(stats.to_string())
# Official historical bloom candidates: exact sequence-presence audit only.
url='https://raw.githubusercontent.com/biocore/American-Gut/68fd6d4b2fa6aeb5b4f5272c6f1006defe5b160e/data/AG/BLOOM.fasta'
ref=P/'sources'/'AGP_BLOOM.fasta'
b=ref.read_bytes() if ref.is_file() else fetch(url,'AGP_BLOOM.fasta')
assert hashlib.sha256(b).hexdigest()=='18e0a6b012a1d6f2109b8f068ba5a78febf6c5111958933adb7341f6d57d43fd','Historical candidate reference changed'
lines=b.decode().splitlines(); refs=[]
for line in lines:
    if line.startswith('>'):refs.append([line[1:],''])
    else:refs[-1][1]+=line.strip().upper()
windows={}
for label,seq in refs:
    for ss in [seq,seq.translate(str.maketrans('ACGT','TGCA'))[::-1]]:
        for j in range(len(ss)-99):windows.setdefault(ss[j:j+100],set()).add(label)
biom=S/'02_统计输入/六份源数据/4168_ASV_table.biom'
with h5py.File(biom) as h:
    decode=lambda v:v.decode() if isinstance(v,bytes) else str(v)
    seq=[decode(v) for v in h['observation/ids'][:]];ids=[decode(v) for v in h['sample/ids'][:]]
    a=h['observation/matrix']; mat=csr_matrix((a['data'][:],a['indices'][:],a['indptr'][:]),shape=(len(seq),len(ids)))
    attrs={k:decode(v) if isinstance(v,(str,bytes)) else np.asarray(v).tolist() for k,v in h.attrs.items()}
index=[ids.index(i) for i in X.index]; mi=mat[:,index];tot=np.asarray(mi.sum(axis=0)).ravel()
r07=pd.read_csv(S/'02_统计输入/冻结统计输入/R07_all_asvs_genus_assignments.tsv',sep='\t').set_index('biom_row')
hits=[]
for i,ss in enumerate(seq):
    if ss not in windows:continue
    count=mi.getrow(i).toarray().ravel();record={'biom_row':i,'sequence':ss,'historical_bloom_labels':'|'.join(sorted(windows[ss])),'reads_all4168':int(mat.getrow(i).sum()),'reads_primary2748':int(count.sum()),'positive_primary2748':int((count>0).sum()),'largest_sample_reads':int(count.max()),'R07_genus':str(r07.loc[i,'genus']),'R07_genus_bootstrap':float(r07.loc[i,'genus_bootstrap'])}
    hits.append(record)
pd.DataFrame(hits).to_csv(P/'bloom_exact_sequence_hits.tsv',sep='\t',index=False)
bloom=np.asarray(mi[[d['biom_row'] for d in hits]].sum(axis=0)).ravel() if hits else np.zeros(len(index))
pd.DataFrame({'sample_id':X.index,'bloom_candidate_reads':bloom,'all_reads':tot,'bloom_candidate_fraction':bloom/tot}).to_csv(P/'bloom_candidate_sample_burden.tsv',sep='\t',index=False)
write('bloom_sequence_audit.json',{'reference_url':url,'reference_sha256':sha(P/'sources'/'AGP_BLOOM.fasta'),'reference_sequences':len(refs),'retrieved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'matching_rule':'100 nt ASV exactly equals any 100 nt contiguous window in forward or reverse-complement historical bloom reference; no mismatches','matched_ASVs':len(hits),'primary2748_total_matched_reads':int(bloom.sum()),'primary2748_positive_people':int((bloom>0).sum()),'median_bloom_candidate_fraction':float(np.median(bloom/tot)),'max_bloom_candidate_fraction':float(np.max(bloom/tot)),'input_biom_sha256':sha(biom),'biom_attributes':attrs,'interpretation':'Candidate retention is a sequence-level observation, not proof of contamination, shipping growth, or exact upstream filtering history.'})
print('BLOOM_HITS',json.dumps(hits,ensure_ascii=False))
