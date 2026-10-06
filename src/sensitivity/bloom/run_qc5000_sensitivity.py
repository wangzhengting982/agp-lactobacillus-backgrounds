"""Post-result fixed subset check; retains each scenario's original 2748-person scale."""
from pathlib import Path
import sys,os,argparse,json,hashlib,warnings
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
p=argparse.ArgumentParser();p.add_argument('--source-root',type=Path,required=True);p.add_argument('--bloom-root',type=Path,required=True);p.add_argument('--runtime',type=Path,action='append',default=[]);a=p.parse_args()
for r in reversed(a.runtime):sys.path.insert(0,str(r))
import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix
from scipy.linalg import block_diag
from scipy.stats import norm
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
sys.stdout.reconfigure(encoding='utf-8');R=a.source_root.resolve();O=a.bloom_root.resolve()
Q=R/'01_核心统计/输入/03_分析与结果/taxonomy_validation/association_sensitivity/R07'
P=R/'02_统计输入/饮食预处理'
read=lambda p:pd.read_csv(p,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
files={'biom':R/'02_统计输入/六份源数据/4168_ASV_table.biom','X':Q/'covariate_design.tsv',
 'foodPC':P/'food_PC_得分_内部复现.tsv','nutrientPC':P/'nutrient_PC_得分_内部复现.tsv'}
for s in ['R07_baseline','R07_filtered']:
 files[s+'_variables']=O/(s+'_rebuilt_variables.tsv');files[s+'_counts']=O/(s+'_target_counts.tsv')
hashes=[dict(key=k,path=str(v),sha256=sha(v))for k,v in files.items()]
pd.DataFrame(hashes).to_csv(O/'qc5000_input_hashes.tsv',sep='\t',index=False,encoding='utf-8-sig')
x0=read(files['X']);assert x0.shape==(2748,20)
with h5py.File(files['biom'])as h:
 samples=[s.decode()if isinstance(s,bytes)else str(s)for s in h['sample/ids'][:]];rr=h['observation/matrix']
 mat=csr_matrix((rr['data'][:],rr['indices'][:],rr['indptr'][:]),shape=(len(h['observation/ids']),len(samples)))[:,[samples.index(s)for s in x0.index]]
tot=np.asarray(mat.sum(axis=0)).ravel();drop=np.asarray(mat[[4246,14343,21871,38369,45119,59390,69046]].sum(axis=0)).ravel();remain=tot-drop
ix=x0.index[remain>=5000];assert len(ix)==2343
PC=pd.concat([read(files['foodPC']).iloc[:,:10],read(files['nutrientPC']).iloc[:,:10]],axis=1).loc[x0.index]
BG=['Anaerococcus','Finegoldia','Peptoniphilus_A'];TT=['Lacticaseibacillus','Lactobacillus'];CV=['depth_nonlab','nonlab_shannon','nonlab_logrichness']
scenarios={};events=[]
for scenario in ['R07_baseline','R07_filtered']:
 v=read(files[scenario+'_variables']).loc[x0.index];y=read(files[scenario+'_counts']).loc[x0.index]
 x=x0.copy();x[CV]=v[CV];x=pd.concat([x,PC],axis=1).loc[ix];z=v.drop(columns=CV).loc[ix];y=y.loc[ix]
 assert x.shape==(2343,40)and x.index.equals(z.index)and x.index.equals(y.index)
 scenarios[scenario]=(x,z,y)
 for th in [1,3]:
  for t in TT:events.append(dict(scenario=scenario,threshold=th,target=t,n=len(y),events=int(y[t].ge(th).sum()),target_reads=int(y[t].sum())))
pd.DataFrame(events).to_csv(O/'qc5000_target_events.tsv',sep='\t',index=False,encoding='utf-8-sig')
meta=dict(started_utc=datetime.now(timezone.utc).isoformat(),n_original=2748,n_retained=2343,n_excluded=405,
 subset_sample_order_sha256=hashlib.sha256('\n'.join(ix).encode()).hexdigest(),selection='filtered total ASV reads >=5000',
 minimum_retained_filtered_total_reads=int(remain[remain>=5000].min()),standardization='unchanged scenario-specific 2748-person standardized values',
 script_sha256=sha(Path(__file__)),plan_sha256=sha(O/'qc5000_fixed_plan.txt'))
(O/'qc5000_run_metadata.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
models=[];contrasts=[];diags=[];independent=[]
def save():
 for rows,name in [(models,'qc5000_model_results.tsv'),(contrasts,'qc5000_contrast_results.tsv'),(diags,'qc5000_diagnostics.tsv'),(independent,'qc5000_independent_checks.tsv')]:pd.DataFrame(rows).to_csv(O/name,sep='\t',index=False,encoding='utf-8-sig')
for scenario in ['R07_baseline','R07_filtered']:
 x,z,y=scenarios[scenario]
 for th in [1,3]:
  for g in BG:
   matched=[c for c in z if c.split('|')[-1]==g];assert len(matched)==1
   d=x.assign(background=z[matched[0]]);design=d.to_numpy(float);fits=[]
   for t in TT:
    yy=y[t].ge(th).to_numpy(float);m=dict(scenario=scenario,threshold=th,background=g,target=t,n=len(y),events=int(yy.sum()),parameters=design.shape[1])
    try:
     rank=int(np.linalg.matrix_rank(design));assert rank==design.shape[1]
     with warnings.catch_warnings(record=True)as ws:
      warnings.simplefilter('always');f=sm.GLM(yy,design,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0')
     mu=np.asarray(f.fittedvalues);info=design.T@(design*(mu*(1-mu))[:,None]);eig=np.linalg.eigvalsh(info)
     u=(design*(yy-mu)[:,None])@np.linalg.inv(info);serr=float(np.max(np.abs(np.sqrt((u*u).sum(0))-f.bse)));score=float(np.max(np.abs(design.T@(yy-mu))))
     valid=bool(f.converged and not ws and np.isfinite(f.params).all()and np.isfinite(f.cov_params()).all()and eig.min()>0 and serr<1e-7 and score<1e-5)
     diags.append(dict(**m,valid=valid,converged=bool(f.converged),rank=rank,hessian_min_eigenvalue=float(eig.min()),hessian_condition=float(np.linalg.cond(info)),score_max=score,HC0_se_max_error=serr,warnings=' | '.join(str(w.message)for w in ws)))
     b=float(f.params[-1]);se=float(f.bse[-1]);models.append(dict(**m,valid=valid,beta=b,se=se,OR=np.exp(b),low=np.exp(b-1.96*se),high=np.exp(b+1.96*se),p=2*norm.sf(abs(b/se))))
     fits.append((f,u)if valid else None)
    except Exception as e:
     models.append(dict(**m,valid=False,p=np.nan,error=repr(e)));diags.append(dict(**m,valid=False,error=repr(e)));fits.append(None)
   cm=dict(scenario=scenario,threshold=th,background=g,n=len(y),target_A=TT[0],target_B=TT[1])
   if any(f is None for f in fits):
    contrasts.append(dict(**cm,valid=False,p=np.nan,error='At least one single-outcome fit invalid; no rescue.'));save();continue
   delta=float(fits[0][0].params[-1]-fits[1][0].params[-1]);se=float(np.sqrt(np.sum((fits[0][1][:,-1]-fits[1][1][:,-1])**2)))
   cr=dict(**cm,valid=True,delta_log_OR=delta,se_joint=se,ratio_OR=np.exp(delta),low=np.exp(delta-1.96*se),high=np.exp(delta+1.96*se),p=2*norm.sf(abs(delta/se)))
   try:
    da=block_diag(design,design);ys=np.r_[y[TT[0]].ge(th).to_numpy(float),y[TT[1]].ge(th).to_numpy(float)]
    with warnings.catch_warnings(record=True)as ws:
     warnings.simplefilter('always');sf=sm.GLM(ys,da,family=sm.families.Binomial()).fit(method='newton',start_params=np.zeros(da.shape[1]),maxiter=200,tol=1e-10,cov_type='cluster',cov_kwds={'groups':np.tile(np.arange(len(y)),2),'use_correction':False})
    k=design.shape[1];co=np.zeros(2*k);co[k-1]=1;co[2*k-1]=-1
    sd=float(co@sf.params);ss=float(np.sqrt(co@sf.cov_params()@co));be=float(np.max(np.abs(sf.params-np.r_[fits[0][0].params,fits[1][0].params])))
    passed=bool(sf.mle_retvals['converged']and not ws and abs(sd-delta)<1e-7 and abs(ss-se)<1e-7 and be<1e-7)
    independent.append(dict(**cm,passed=passed,delta_abs_error=abs(sd-delta),SE_abs_error=abs(ss-se),all_beta_max_error=be,stacked_converged=bool(sf.mle_retvals['converged']),warnings=' | '.join(str(w.message)for w in ws)))
    if not passed:cr.update(valid=False,p=np.nan,error='Independent numerical verification failed; no rescue.')
   except Exception as e:
    independent.append(dict(**cm,passed=False,error=repr(e)));cr.update(valid=False,p=np.nan,error='Independent fit failed: '+repr(e))
   contrasts.append(cr);save();print(json.dumps(cr),flush=True)
c=pd.DataFrame(contrasts)
for _,idx in c.groupby(['scenario','threshold']).groups.items():
 q=multipletests(c.loc[idx,'p'].fillna(1),method='fdr_bh')[1];c.loc[idx,'q_BH_selected3']=np.where(c.loc[idx,'p'].notna(),q,np.nan)
c.to_csv(O/'qc5000_contrast_results.tsv',sep='\t',index=False,encoding='utf-8-sig')
diag=pd.DataFrame(diags);ind=pd.DataFrame(independent)
summary=dict(completed_utc=datetime.now(timezone.utc).isoformat(),n_models=len(models),valid_models=int(pd.DataFrame(models).valid.sum()),n_contrasts=len(c),valid_contrasts=int(c.valid.sum()),
 n_retained=2343,subset_sample_order_sha256=meta['subset_sample_order_sha256'],
 all_models_converged=bool(diag.get('converged',pd.Series(False)).fillna(False).all()),all_rank_full=bool((diag['rank']==diag.parameters).all()),
 minimum_Hessian_eigenvalue=float(diag.hessian_min_eigenvalue.min()),maximum_Hessian_condition=float(diag.hessian_condition.max()),maximum_abs_score=float(diag.score_max.max()),
 all_independent_passed=bool(ind.passed.all()),independent_delta_max_error=float(ind.delta_abs_error.max()),independent_SE_max_error=float(ind.SE_abs_error.max()),independent_all_beta_max_error=float(ind.all_beta_max_error.max()),
 all_inputs_unchanged=all(sha(Path(r['path']))==r['sha256']for r in hashes),
 each_specification=[dict(scenario=s,threshold=int(t),valid=int(v.valid.sum()),q_below05=int(v.q_BH_selected3.lt(.05).sum()))for(s,t),v in c.groupby(['scenario','threshold'])])
(O/'qc5000_validation_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
