"""Exploratory attribution of genus-pattern differences across the entire 156-background family.
The extension and complete correction families are declared before these fits, following the
omnibus/two-background analyses. No significant subset replaces the full comparison family.
"""
from pathlib import Path
from core_paths import DATA, OUT, PRIMARY, TAXONOMY_INPUT
import sys,os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd,statsmodels.api as sm,itertools,json,warnings
from scipy.stats import chi2,norm
from statsmodels.stats.multitest import multipletests
R=DATA;O=OUT;P=(DATA/'输入/03_分析与结果/association_models')
read=lambda n:pd.read_csv(PRIMARY/n,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
X=read('association_covariate_design.tsv');Z=read('association_standardized_clr.tsv').loc[X.index];Y=read('association_target_counts.tsv').loc[X.index].ge(1).astype(float)
T=Y.columns.tolist();assert len(T)==6
old=pd.read_csv(P/'single_genus_all_associations.tsv',sep='\t').set_index(['feature','taxon'])
rows=[];global_rows=[];quality=[];refitted=[]
for fi,feature in enumerate(Z):
    a=np.column_stack([X,Z[feature]]);k=a.shape[1];b=[];u=[]
    for t in T:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always')
            f=sm.GLM(Y[t].to_numpy(),a,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0')
        assert f.converged and not ws
        pr=np.asarray(f.fittedvalues);infl=(a*(Y[t].to_numpy()-pr)[:,None])@np.linalg.inv(a.T@(a*(pr*(1-pr))[:,None]))
        assert abs(f.params[-1]-old.loc[(feature,t),'beta_log_OR_per_SD_CLR'])<1e-8
        assert abs(f.bse[-1]-old.loc[(feature,t),'SE_HC0'])<1e-8
        b.append(f.params[-1]);u.append(infl[:,-1]);quality.append(dict(feature=feature,taxon=t,reproduced=True,beta_absolute_error=float(abs(f.params[-1]-old.loc[(feature,t),'beta_log_OR_per_SD_CLR'])),se_absolute_error=float(abs(f.bse[-1]-old.loc[(feature,t),'SE_HC0']))))
        refitted.append(dict(taxon=t,feature=feature,n=len(Y),events=int(Y[t].sum()),beta_log_OR_per_SD_CLR=float(f.params[-1]),SE_HC0=float(f.bse[-1]),OR=float(np.exp(f.params[-1])),CI_low=float(np.exp(f.params[-1]-1.95996398454*f.bse[-1])),CI_high=float(np.exp(f.params[-1]+1.95996398454*f.bse[-1])),p=float(f.pvalues[-1]),valid=bool(f.converged)))
    b=np.array(b);u=np.column_stack(u);v=u.T@u;c=np.column_stack([-np.ones(5),np.eye(5)]);d=c@b;cv=c@v@c.T;s=float(d@np.linalg.solve(cv,d))
    global_rows.append(dict(feature=feature,chi2=s,df=5,p=float(chi2.sf(s,5))))
    for i,j in itertools.combinations(range(6),2):
        d=b[i]-b[j];se=np.sqrt(v[i,i]+v[j,j]-2*v[i,j]);p=float(2*norm.sf(abs(d/se)))
        rows.append(dict(feature=feature,taxon_A=T[i],taxon_B=T[j],beta_A=b[i],beta_B=b[j],delta_log_OR=d,SE_joint=se,ratio_OR=np.exp(d),CI_low=np.exp(d-1.96*se),CI_high=np.exp(d+1.96*se),p=p))
    if fi%25==0:print('background',fi+1,'/156',flush=True)
res=pd.DataFrame(rows);glob=pd.DataFrame(global_rows)
res['q_BH_2340']=multipletests(res.p,method='fdr_bh')[1];res['p_Holm_2340']=multipletests(res.p,method='holm')[1]
glob['q_BH_156']=multipletests(glob.p,method='fdr_bh')[1]
res.to_csv(O/'表S22f_全部2340项背景属系数比较.tsv',sep='\t',index=False);glob.to_csv(O/'表S22e_全部156背景属整体差异.tsv',sep='\t',index=False)
pd.DataFrame(quality).to_csv(O/'936原模型逐项重现.tsv',sep='\t',index=False)
print('SIGNIFICANT GLOBAL',glob.q_BH_156.lt(.05).sum(),'PAIRS',res.q_BH_2340.lt(.05).sum(),flush=True)
print(res.sort_values('q_BH_2340').head(15).to_string(index=False),flush=True)

refitted=pd.DataFrame(refitted);refitted['q_BH_all']=multipletests(refitted.p,method='fdr_bh')[1]
refitted.to_csv(O/'refitted_936_associations.tsv',sep='\t',index=False)
