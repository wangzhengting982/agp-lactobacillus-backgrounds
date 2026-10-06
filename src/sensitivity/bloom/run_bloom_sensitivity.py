"""R07 ASV-level historical-bloom-candidate sensitivity, not a contamination diagnosis."""
from pathlib import Path
import sys, os, json, argparse, hashlib, warnings
from datetime import datetime, timezone
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
p=argparse.ArgumentParser()
p.add_argument('--source-root',type=Path,required=True)
p.add_argument('--provenance-root',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--runtime',type=Path,action='append',default=[])
args=p.parse_args()
for rt in reversed(args.runtime):sys.path.insert(0,str(rt))
import numpy as np,pandas as pd,h5py,scipy,statsmodels
from scipy.sparse import csr_matrix,coo_matrix
from scipy.linalg import block_diag
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
sys.stdout.reconfigure(encoding='utf-8')
R=args.source_root.resolve(); O=args.output_dir.resolve(); P=args.provenance_root.resolve()
O.mkdir(parents=True,exist_ok=True)
Q=R/'01_核心统计/输入/03_分析与结果/taxonomy_validation/association_sensitivity/R07'
paths={
 'biom':R/'02_统计输入/六份源数据/4168_ASV_table.biom',
 'taxonomy':R/'02_统计输入/冻结统计输入/R07_all_asvs_genus_assignments.tsv',
 'old_X':Q/'covariate_design.tsv','old_Z':Q/'standardized_clr.tsv.gz','old_Y':Q/'target_counts.tsv',
 'old_contrasts':R/'01_核心统计/新增分析/表S24h_三项突出差异的阈值与分类复核.tsv',
 'food_pc':R/'02_统计输入/饮食预处理/food_PC_得分_内部复现.tsv',
 'nutrient_pc':R/'02_统计输入/饮食预处理/nutrient_PC_得分_内部复现.tsv',
 'candidates':P/'bloom_exact_sequence_hits.tsv',
 'bloom_fasta':P/'sources/AGP_BLOOM.fasta',
 'candidate_audit':P/'bloom_sequence_audit.json'}
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
def read(x):
 d=pd.read_csv(x,sep='\t',dtype={'sample_id':str}).set_index('sample_id');assert d.index.is_unique;return d
def save(rows,name):pd.DataFrame(rows).to_csv(O/name,sep='\t',index=False,encoding='utf-8-sig')
def dec(x):return x.decode() if isinstance(x,bytes) else str(x)
hashes=[dict(key=k,path=str(v),bytes=v.stat().st_size,sha256=sha(v))for k,v in paths.items()]
save(hashes,'input_hashes.tsv')
X0=read(paths['old_X']);Z0=read(paths['old_Z']).loc[X0.index];Y0=read(paths['old_Y']).loc[X0.index]
assert X0.shape==(2748,20)
PC=pd.concat([read(paths['food_pc']).iloc[:,:10],read(paths['nutrient_pc']).iloc[:,:10]],axis=1).loc[X0.index]
BG=['Anaerococcus','Finegoldia','Peptoniphilus_A'];TT=['Lacticaseibacillus','Lactobacillus']
SIX=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
CAND=[4246,14343,21871,38369,45119,59390,69046]
hits=pd.read_csv(paths['candidates'],sep='\t');assert sorted(hits.biom_row.tolist())==CAND
assign=pd.read_csv(paths['taxonomy'],sep='\t').sort_values('biom_row').reset_index(drop=True)
for col in ['genus','domain','order','family']:assign[col]=assign[col].fillna('')
assert assign.biom_row.tolist()==list(range(len(assign)))
with h5py.File(paths['biom']) as f:
 samples=[dec(x)for x in f['sample/ids'][:]]; seq=[dec(x)for x in f['observation/ids'][:]]
 m=f['observation/matrix'];allmat=csr_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=(len(seq),len(samples)))
 assert np.equal(allmat.data,np.round(allmat.data)).all()
 mat=allmat[:,[samples.index(s)for s in X0.index]].astype(np.int64)
assert len(seq)==len(assign)==95567
for row in hits.itertuples():assert seq[row.biom_row]==row.sequence
# Independently reproduce exact-window/strand candidate membership from the official fasta.
refs=[];buf=''
for line in paths['bloom_fasta'].read_text().splitlines():
 if line.startswith('>'):
  if buf:refs.append(buf.upper())
  buf=''
 else:buf+=line.strip()
if buf:refs.append(buf.upper())
windows=set()
for s in refs:
 for strand in [s,s.translate(str.maketrans('ACGT','TGCA'))[::-1]]:
  windows.update(strand[i:i+100]for i in range(len(strand)-99))
assert [i for i,s in enumerate(seq)if s.upper() in windows]==CAND
assert all(len(s)==100 for s in seq)
primary=.8
bacterial=assign.domain.eq('Bacteria')&assign.domain_bootstrap.ge(primary)
forbidden=(assign['order'].eq('Lactobacillales')&assign.order_bootstrap.ge(primary))|(assign.family.eq('Bifidobacteriaceae')&assign.family_bootstrap.ge(primary))|(assign.genus.isin(SIX)&assign.genus_bootstrap.ge(primary))
named=assign.genus.ne('')&~assign.genus.str.contains('unclassified|uncultured',case=False)&assign.genus_bootstrap.ge(primary)
background_mask=(bacterial&~forbidden&named).to_numpy()
depth_mask=(bacterial&~forbidden).to_numpy()
family=assign.family.where(assign.family_bootstrap.ge(primary),'unclassified_family')
labels=('microbe::'+family+'|'+assign.genus).to_numpy()
unique=sorted(set(labels[background_mask]));lookup={s:i for i,s in enumerate(unique)}
ix=np.flatnonzero(background_mask)
incidence=coo_matrix((np.ones(len(ix),dtype=np.int64),([lookup[labels[i]]for i in ix],ix)),shape=(len(unique),len(seq))).tocsr()
assert set(Z0.columns)<=set(unique)
matched={g:[c for c in Z0.columns if c.split('|')[-1]==g]for g in BG}
assert all(len(v)==1 for v in matched.values());matched={g:v[0]for g,v in matched.items()}
metadata=dict(started_utc=datetime.now(timezone.utc).isoformat(),n=2748,
 matched_features=matched,n_fixed_CLR_genera=len(Z0.columns),n_all_named_background_genera=len(unique),
 candidates=CAND,software={'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__,'h5py':h5py.__version__},
 sample_order_sha256=hashlib.sha256('\n'.join(X0.index).encode()).hexdigest(),script_sha256=sha(Path(__file__)),plan_sha256=sha(O/'fixed_plan.txt'))
(O/'run_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
units=[];scenarios={};composition=[];target_summaries=[];candidate_contributions=[]
for row in hits.itertuples():
 i=row.biom_row;candidate_contributions.append(dict(biom_row=i,genus=assign.loc[i,'genus'],genus_bootstrap=assign.loc[i,'genus_bootstrap'],
  reads=int(mat[i].sum()),positive_people=int((mat[i].toarray().ravel()>0).sum()),
  included_in_background=bool(background_mask[i]),included_in_depth=bool(depth_mask[i]),
  included_target=assign.loc[i,'genus'] if assign.loc[i,'genus']in SIX and assign.loc[i,'genus_bootstrap']>=.8 else '',
  background_label=labels[i]if background_mask[i]else ''))
save(candidate_contributions,'candidate_contributions.tsv')
for scenario in ['R07_baseline','R07_filtered']:
 keep=np.ones(len(seq),dtype=np.int64)
 if scenario=='R07_filtered':keep[CAND]=0
 mm=mat.multiply(keep[:,None]).tocsr();mm.eliminate_zeros()
 bg=pd.DataFrame((incidence@mm).T.toarray(),index=X0.index,columns=unique)
 den=bg.sum(axis=1);assert den.gt(0).all()
 log=np.log(bg[Z0.columns].div(den,axis=0)+.00005);clr=log.sub(log.mean(axis=1),axis=0)
 # Fixed feature set stays in CLR geometric mean even if a nonfocal genus becomes absent.
 z=clr.sub(clr.mean()).div(clr.std(ddof=0))
 assert np.isfinite(z[list(matched.values())]).all().all()
 depth=np.asarray(mm[depth_mask].sum(axis=0)).ravel()
 prob=bg.div(den,axis=0)
 shannon=-(prob*np.log(prob.where(prob>0,1))).sum(axis=1)
 raw=pd.DataFrame({'depth_nonlab':np.log10(1+depth),'nonlab_shannon':shannon,'nonlab_logrichness':np.log1p(bg.gt(0).sum(axis=1))},index=X0.index)
 x=X0.copy()
 for c in raw:
  x[c]=(raw[c]-raw[c].mean())/raw[c].std(ddof=0)
  units.append(dict(scenario=scenario,kind='rebuilt_covariate',variable=c,mean=raw[c].mean(),SD_ddof0=raw[c].std(ddof=0)))
 for g,c in matched.items():units.append(dict(scenario=scenario,kind='focal_CLR',variable=g,mean=clr[c].mean(),SD_ddof0=clr[c].std(ddof=0)))
 y=pd.DataFrame({g:np.asarray(mm[(assign.genus.eq(g)&assign.genus_bootstrap.ge(.8)).to_numpy()].sum(axis=0)).ravel()for g in SIX},index=X0.index)
 for g in TT:
  for th in [1,3]:target_summaries.append(dict(scenario=scenario,target=g,threshold=th,n=len(y),events=int(y[g].ge(th).sum()),total_target_reads=int(y[g].sum())))
 composition.append(dict(scenario=scenario,total_ASV_reads=int(mm.sum()),named_background_reads=int(den.sum()),depth_background_reads=int(depth.sum()),
  minimum_background_depth=int(depth.min()),minimum_named_background_depth=int(den.min()),
  n_zero_variance_fixed_CLR=int(clr.std(ddof=0).eq(0).sum()),n_fixed_CLR_genera=len(Z0.columns)))
 scenarios[scenario]=(x,z,y)
 y.to_csv(O/(scenario+'_target_counts.tsv'),sep='\t')
 pd.concat([x[list(raw)],z[list(matched.values())]],axis=1).to_csv(O/(scenario+'_rebuilt_variables.tsv'),sep='\t')
save(units,'standardization_units.tsv');save(composition,'composition_summary.tsv');save(target_summaries,'target_events.tsv')
xb,zb,yb=scenarios['R07_baseline']
checks=dict(target_max_abs_count_error=int(np.abs(yb[Y0.columns]-Y0).to_numpy().max()),
 CLR_max_abs_error=float(np.abs(zb[Z0.columns]-Z0).to_numpy().max()),
 covariate_max_abs_error=float(np.abs(xb-X0).to_numpy().max()),
 independent_official_candidate_match=True)
(O/'reconstruction_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
assert checks['target_max_abs_count_error']==0 and checks['CLR_max_abs_error']<1e-8 and checks['covariate_max_abs_error']<1e-8,checks
print('Reconstruction gate passed '+json.dumps(checks),flush=True)
old=pd.read_csv(paths['old_contrasts'],sep='\t')
models=[];contrasts=[];diags=[];independent=[];reproduction=[]
def fitpair(scenario,threshold,g,with_pc=True):
 x,z,y=scenarios[scenario];d=pd.concat([x,PC],axis=1)if with_pc else x.copy()
 d=d.assign(background=z[matched[g]]);a=d.to_numpy(float);assert np.linalg.matrix_rank(a)==a.shape[1]
 fits=[]
 for t in TT:
  yy=y[t].ge(threshold).to_numpy(float)
  with warnings.catch_warnings(record=True)as ws:
   warnings.simplefilter('always');f=sm.GLM(yy,a,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0')
  mu=np.asarray(f.fittedvalues);info=a.T@(a*(mu*(1-mu))[:,None]);eig=np.linalg.eigvalsh(info)
  u=(a*(yy-mu)[:,None])@np.linalg.inv(info)
  se_error=float(np.max(np.abs(np.sqrt((u*u).sum(axis=0))-f.bse)));score=float(np.max(np.abs(a.T@(yy-mu))))
  good=bool(f.converged and not ws and np.isfinite(f.params).all() and np.isfinite(f.cov_params()).all() and eig.min()>0 and se_error<1e-7 and score<1e-5)
  meta=dict(scenario=scenario,threshold=threshold,background=g,target=t,with_diet_PC=with_pc,n=len(y),events=int(yy.sum()),parameters=a.shape[1])
  diags.append(dict(**meta,valid=good,converged=bool(f.converged),rank=int(np.linalg.matrix_rank(a)),hessian_min_eigenvalue=float(eig.min()),hessian_condition=float(np.linalg.cond(info)),score_max=score,HC0_se_max_error=se_error,warnings=' | '.join(str(w.message)for w in ws)))
  beta=float(f.params[-1]);se=float(f.bse[-1]);models.append(dict(**meta,valid=good,beta=beta,se=se,OR=np.exp(beta),low=np.exp(beta-1.96*se),high=np.exp(beta+1.96*se),p=2*norm.sf(abs(beta/se))))
  assert good,diags[-1];fits.append((f,u))
 delta=float(fits[0][0].params[-1]-fits[1][0].params[-1]);se=float(np.sqrt(np.sum((fits[0][1][:,-1]-fits[1][1][:,-1])**2)))
 row=dict(scenario=scenario,threshold=threshold,background=g,with_diet_PC=with_pc,n=len(y),target_A=TT[0],target_B=TT[1],valid=True,delta_log_OR=delta,se_joint=se,ratio_OR=np.exp(delta),low=np.exp(delta-1.96*se),high=np.exp(delta+1.96*se),p=2*norm.sf(abs(delta/se)))
 contrasts.append(row)
 stacked=block_diag(a,a);ys=np.r_[y[TT[0]].ge(threshold).to_numpy(float),y[TT[1]].ge(threshold).to_numpy(float)]
 with warnings.catch_warnings(record=True)as ws:
  warnings.simplefilter('always');sf=sm.GLM(ys,stacked,family=sm.families.Binomial()).fit(method='newton',start_params=np.zeros(stacked.shape[1]),maxiter=200,tol=1e-10,cov_type='cluster',cov_kwds={'groups':np.tile(np.arange(len(y)),2),'use_correction':False})
 k=a.shape[1];co=np.zeros(2*k);co[k-1]=1;co[2*k-1]=-1
 sd=float(co@sf.params);ss=float(np.sqrt(co@sf.cov_params()@co));pe=float(np.max(np.abs(sf.params-np.r_[fits[0][0].params,fits[1][0].params])))
 passed=bool(sf.mle_retvals['converged']and not ws and abs(sd-delta)<1e-7 and abs(ss-se)<1e-7 and pe<1e-7)
 independent.append(dict(scenario=scenario,threshold=threshold,background=g,with_diet_PC=with_pc,passed=passed,delta_abs_error=abs(sd-delta),SE_abs_error=abs(ss-se),all_beta_max_error=pe,stacked_converged=bool(sf.mle_retvals['converged']),warnings=' | '.join(str(w.message)for w in ws)))
 assert passed,independent[-1]
 print(f'{scenario} threshold={threshold} dietPC={with_pc} {g}: ratio={row["ratio_OR"]:.9f},p={row["p"]:.4g}',flush=True)
 return row
def checkpoint():
 save(models,'all_model_results.tsv');save(contrasts,'all_contrast_results.tsv');save(diags,'model_diagnostics.tsv');save(independent,'independent_checks.tsv');save(reproduction,'original_R07_reproduction.tsv')
try:
 for g in BG:
  row=fitpair('R07_baseline',1,g,False)
  prior=old[old.scenario.eq('R07_reads1')&old.feature.str.endswith('|'+g)].iloc[0]
  errs={'delta_abs_error':abs(row['delta_log_OR']-prior.delta_log_OR),'SE_abs_error':abs(row['se_joint']-prior.SE_joint),'ratio_abs_error':abs(row['ratio_OR']-prior.ratio_OR)}
  reproduction.append(dict(background=g,passed=max(errs.values())<1e-8,**errs));checkpoint();assert reproduction[-1]['passed'],reproduction[-1]
 for scenario in ['R07_baseline','R07_filtered']:
  for threshold in [1,3]:
   for g in BG:fitpair(scenario,threshold,g,True);checkpoint()
 df=pd.DataFrame(contrasts)
 for _,idx in df.groupby(['scenario','threshold','with_diet_PC']).groups.items():df.loc[idx,'q_BH_selected3']=multipletests(df.loc[idx,'p'],method='fdr_bh')[1]
 df.to_csv(O/'all_contrast_results.tsv',sep='\t',index=False,encoding='utf-8-sig')
 df[df.with_diet_PC].to_csv(O/'bloom_contrast_results.tsv',sep='\t',index=False,encoding='utf-8-sig')
 pd.DataFrame(models).query('with_diet_PC').to_csv(O/'bloom_model_results.tsv',sep='\t',index=False,encoding='utf-8-sig')
 status=dict(status='success',finished_utc=datetime.now(timezone.utc).isoformat(),planned_models=24,planned_contrasts=12,baseline_gate_models=6,baseline_gate_contrasts=3,
  all_models_valid=all(r['valid']for r in diags),all_independent_checks_passed=all(r['passed']for r in independent),all_input_hashes_unchanged=all(sha(paths[r['key']])==r['sha256']for r in hashes))
 (O/'completion_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(status,ensure_ascii=False),flush=True)
 print(df[df.with_diet_PC][['scenario','threshold','background','ratio_OR','low','high','p','q_BH_selected3']].to_string(index=False),flush=True)
except Exception as e:
 checkpoint();(O/'completion_status.json').write_text(json.dumps(dict(status='failed',error=repr(e)),ensure_ascii=False,indent=2),encoding='utf-8');raise
