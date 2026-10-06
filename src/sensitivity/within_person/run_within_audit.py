from pathlib import Path
import os, sys, json, hashlib, warnings, datetime, argparse
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import pandas as pd
import scipy, statsmodels
from scipy.stats import t as tdist, norm
from statsmodels.discrete.conditional_models import ConditionalLogit
from statsmodels.stats.multitest import multipletests

parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[1]/'source_package')
parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent)
args=parser.parse_args(); root=args.source; out=args.output;out.mkdir(parents=True,exist_ok=True)
start=datetime.datetime.now(datetime.timezone.utc).isoformat()
read=lambda p:pd.read_csv(p,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
inp=root/'02_统计输入'; frozen=inp/'冻结统计输入'
files={'samples':frozen/'纵向分析样本_内部索引.tsv','counts':frozen/'lactobacillales_genus_counts_all_samples.tsv','X':frozen/'association_covariate_design.tsv','Z':frozen/'association_standardized_clr.tsv','taxa':inp/'六份源数据'/'4168_annotated_feature_table.tsv','previous':root/'03_饮食与菌群深化'/'完整分析结果'/'L01_同人检出变化全部规格.tsv'}
hashes={k:{'path':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in files.items()}
(out/'run_started.json').write_text(json.dumps({'started_utc':start,'inputs':hashes,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__},ensure_ascii=False,indent=2),encoding='utf-8')
v=read(files['samples']);yy=read(files['counts']);X=read(files['X']);Z=read(files['Z']);old=pd.read_csv(files['previous'],sep='\t')
assert len(v)==664 and v.host_subject_id.nunique()==153 and len(X)==2748
BG=['Anaerococcus','Finegoldia','Peptoniphilus_A']; TT=['Lacticaseibacillus','Lactobacillus'];BC=['microbe::Peptoniphilaceae|'+g for g in BG]
ids=list(dict.fromkeys(list(X.index)+list(v.index)))
first=pd.read_csv(files['taxa'],sep='\t',nrows=0).columns[0]
tax=pd.read_csv(files['taxa'],sep='\t',index_col=0,usecols=[first]+ids)
forbidden=tax.index.str.contains(r'(?:^|;)o__Lactobacillales(?:;|$)|(?:^|;)f__Bifidobacteriaceae(?:;|$)')
ok=tax.index.str.startswith('d__Bacteria;')&~forbidden
ge=tax.index.to_series().str.extract(r'(?:^|;)g__([^;]+)',expand=False)
fa=tax.index.to_series().str.extract(r'(?:^|;)f__([^;]+)',expand=False).fillna('unclassified_family')
known=ge.notna()&~ge.fillna('').str.contains('unclassified|uncultured',case=False)
allbg=tax.loc[ok&known].groupby((fa+'|'+ge)[ok&known]).sum().T;allbg.columns=['microbe::'+c for c in allbg]
den=allbg.sum(axis=1);lg=np.log(allbg.reindex(columns=Z.columns).div(den,axis=0)+.00005)
clr=lg.sub(lg.mean(axis=1),axis=0);allz=clr.sub(clr.loc[X.index].mean()).div(clr.loc[X.index].std(ddof=0))
clr_error=float(np.abs(allz.loc[X.index,Z.columns]-Z).to_numpy().max());assert clr_error<1e-10
depth=np.log10(tax.loc[ok].sum(axis=0)+1).loc[v.index];v['depth']=(depth-depth.mean())/depth.std(ddof=0)
most=v.host_subject_id.value_counts().idxmax();specs=[(1,'all'),(1,'drop_most_sampled'),(3,'all')]
rows=[];diags=[]
for th,spec in specs:
 for g,c in zip(BG,BC):
  for target in TT:
   a=v.copy()
   if spec=='drop_most_sampled':a=a.loc[a.host_subject_id.ne(most)]
   a['y']=yy.loc[a.index,target].ge(th).astype(int)
   vary=a.groupby('host_subject_id').y.nunique().eq(2);a=a.loc[a.host_subject_id.isin(vary[vary].index)]
   d=pd.DataFrame({'background':allz.loc[a.index,c],'depth':a.depth,'time_year':a.days_since_first/365.25},index=a.index)
   G=int(a.host_subject_id.nunique());meta=dict(threshold=th,spec=spec,background=g,target=target,people=G,samples=len(a),events=int(a.y.sum()))
   try:
    assert np.isfinite(d.to_numpy()).all() and np.linalg.matrix_rank(d)==3
    with warnings.catch_warnings(record=True) as ws:
     warnings.simplefilter('always')
     model=ConditionalLogit(a.y,d,groups=a.host_subject_id)
     fit=model.fit(method='bfgs',maxiter=500,disp=False)
    params=np.asarray(fit.params).copy()
    initial_params=params.copy()
    # ConditionalLogit.fit in this runtime does not forward optimizer kwargs.
    # Explicit Newton refinement ensures a genuinely small absolute score.
    refinement_steps=0
    for iteration in range(20):
     score=model.score(params)
     if np.abs(score).max()<1e-7: break
     hh=1e-5
     hi=-np.column_stack([(model.score(params+np.eye(3)[j]*hh)-model.score(params-np.eye(3)[j]*hh))/(2*hh) for j in range(3)])
     hi=(hi+hi.T)/2
     step=np.linalg.solve(hi,score);ll=model.loglike(params);alpha=1.0
     while alpha>1e-8 and model.loglike(params+alpha*step)<ll-1e-10: alpha/=2
     assert alpha>1e-8, 'Newton line search failed'
     params=params+alpha*step;refinement_steps+=1
    score=model.score(params)
    scores=np.asarray([model.score_grp(j,params) for j in range(model._n_groups)])
    sum_error=float(np.abs(scores.sum(0)-score).max());assert sum_error<1e-8
    # Central finite differences, separate from statsmodels forward-difference Hessian.
    h=1e-5;info=-np.column_stack([(model.score(params+np.eye(3)[j]*h)-model.score(params-np.eye(3)[j]*h))/(2*h) for j in range(3)])
    info=(info+info.T)/2; eig=np.linalg.eigvalsh(info);assert eig.min()>0
    bread=np.linalg.inv(info);cov0=bread@scores.T@scores@bread;cov=cov0*G/(G-1)
    all_fd=np.asarray([[(model.loglike_grp(j,params+np.eye(3)[k]*h)-model.loglike_grp(j,params-np.eye(3)[k]*h))/(2*h) for k in range(3)] for j in range(model._n_groups)])
    grad_error=float(np.abs(all_fd-scores).max());assert grad_error<2e-5
    sm_info=-model.hessian(params);sm_info=(sm_info+sm_info.T)/2
    curvature_error=float(np.abs(np.diag(bread)**.5-np.diag(np.linalg.inv(sm_info))**.5).max());assert curvature_error<1e-4
    prior=old.loc[old.threshold.eq(th)&old.spec.eq(spec)&old.background.eq(g)&old.target.eq(target)].iloc[0]
    beta=float(params[0]);se=float(np.sqrt(cov[0,0]));crit=float(tdist.ppf(.975,G-1))
    prior_error=abs(beta-float(prior.beta));assert prior_error<.002
    score_max=float(np.abs(score).max());assert score_max<1e-7
    row=dict(**meta,valid=True,beta=beta,OR=float(np.exp(beta)),model_se=float(np.sqrt(bread[0,0])),cluster_se_uncorrected=float(np.sqrt(cov0[0,0])),cluster_se=se,low=float(np.exp(beta-crit*se)),high=float(np.exp(beta+crit*se)),p=float(2*tdist.sf(abs(beta/se),G-1)),reference='t(G-1), person-clustered sandwich, G/(G-1)',previous_OR=float(prior.OR),previous_q=float(prior.q_BH_6),beta_error_vs_archive=prior_error)
    diag=dict(**meta,valid=True,group_score_error=grad_error,sum_score_error=sum_error,score_max=score_max,information_min_eigenvalue=float(eig.min()),information_condition=float(np.linalg.cond(info)),model_se_check_error=curvature_error,refinement_steps=refinement_steps,refinement_max_parameter_change=float(np.abs(params-initial_params).max()),warnings=[str(w.message) for w in ws])
   except Exception as e:
    row=dict(**meta,valid=False,p=np.nan,error=repr(e));diag=dict(**meta,valid=False,error=repr(e))
   rows.append(row);diags.append(diag)
   pd.DataFrame(rows).to_csv(out/'within_results_partial.csv',index=False,encoding='utf-8-sig')
   print(json.dumps({k:row[k] for k in ['threshold','spec','background','target','valid','p']},ensure_ascii=False),flush=True)
rr=pd.DataFrame(rows)
for th,spec in specs:
 k=rr.threshold.eq(th)&rr.spec.eq(spec)
 vals=rr.loc[k,'p'].fillna(1).to_numpy();q=multipletests(vals,method='fdr_bh')[1]
 rr.loc[k,'q_BH_6']=np.where(rr.loc[k,'p'].notna(),q,np.nan)
rr.to_csv(out/'within_person_cluster_results.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(diags).to_csv(out/'within_person_diagnostics.csv',index=False,encoding='utf-8-sig')
summary={'started_utc':start,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'CLR_rebuild_max_error':clr_error,'n_models':len(rr),'valid_models':int(rr.valid.sum()),'max_absolute_score':float(pd.DataFrame(diags).score_max.max()),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'inputs':hashes,'specification_counts':[{ 'threshold':th,'spec':spec,'significant_cluster_q':int((rr.loc[rr.threshold.eq(th)&rr.spec.eq(spec),'q_BH_6']<.05).sum())} for th,spec in specs]}
(out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
