from pathlib import Path
import os,sys,json,warnings,itertools
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
# Dependencies are supplied by the project environment.
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.stats import norm,chi2
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[2]
ARCH=ROOT/'data/archive'
RUN=Path(os.environ.get('AGP_DEEPENING_OUTPUT',ROOT/'runs/deepening'))
PREP=Path(os.environ.get('AGP_PREPARATION_ROOT',ROOT/'runs/preparation/project'))
FROZEN=ARCH/'02_统计输入/冻结统计输入'
A=Path(os.environ.get('AGP_ASSOC_ROOT',PREP/'03_分析与结果/association_models'))
D=RUN/'diet_preprocessing';D.mkdir(parents=True,exist_ok=True)
RAW=ARCH/'02_统计输入/六份源数据'
R=RUN;O=RUN/'results';O.mkdir(parents=True,exist_ok=True)

read=lambda p:pd.read_csv(p,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
X=read(A/'association_covariate_design.tsv');Z=read(A/'association_standardized_clr.tsv').loc[X.index];Y=read(A/'association_target_counts.tsv').loc[X.index]
T=Y.columns.tolist();TT=['Lacticaseibacillus','Lactobacillus'];BG=['Anaerococcus','Finegoldia','Peptoniphilus_A'];BC=['microbe::Peptoniphilaceae|'+g for g in BG]
PC=pd.concat([read(D/'food_PC_得分_内部复现.tsv').iloc[:,:10],read(D/'nutrient_PC_得分_内部复现.tsv').iloc[:,:10]],axis=1).loc[X.index]
XD=pd.concat([X,PC],axis=1)
def scale(x):return (x-x.mean())/x.std(ddof=0)
def save(df,name):pd.DataFrame(df).to_csv(O/(name+'.tsv'),sep='\t',index=False)
def bh(df,column='p',out='q_BH'):
 df=pd.DataFrame(df).copy();v=pd.to_numeric(df[column],errors='coerce');good=v.notna()&np.isfinite(v)
 # Keep non-estimable attempted tests in the declared family with p=1, but leave their reported q missing.
 vals=multipletests(v.fillna(1).to_numpy(),method='fdr_bh')[1];df[out]=np.where(good,vals,np.nan);return df
def glm(d,y):
 a=d.to_numpy(float);yy=np.asarray(y,float)
 try:
  assert np.isfinite(a).all() and np.isfinite(yy).all() and np.linalg.matrix_rank(a)==a.shape[1]
  with warnings.catch_warnings(record=True) as ws:
   warnings.simplefilter('always');f=sm.GLM(yy,a,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0')
  good=f.converged and np.isfinite(f.params).all() and np.isfinite(f.bse).all() and not ws
  if not good:return None,{'valid':False,'error':';'.join(str(w.message) for w in ws) or 'nonconverged/nonfinite'}
  mu=np.asarray(f.fittedvalues);u=(a*(yy-mu)[:,None])@np.linalg.inv(a.T@(a*(mu*(1-mu))[:,None]))
  assert np.max(np.abs(np.sqrt((u*u).sum(axis=0))-f.bse))<1e-6
  return {'fit':f,'b':np.asarray(f.params),'v':np.asarray(f.cov_params()),'u':u,'cols':list(d.columns),'d':d},{'valid':True,'n':len(d),'events':int(yy.sum()),'parameters':a.shape[1]}
 except Exception as e:return None,{'valid':False,'error':repr(e),'n':len(d),'events':int(yy.sum())}
def effect(f,term,meta):
 if f is None:return dict(**meta,p=np.nan)
 j=f['cols'].index(term);b=f['b'][j];se=np.sqrt(f['v'][j,j]);return dict(**meta,beta=b,se=se,OR=np.exp(b),low=np.exp(b-1.96*se),high=np.exp(b+1.96*se),p=2*norm.sf(abs(b/se)))
def contrast(fs,term,meta):
 if any(f is None for f in fs):return dict(**meta,valid=False,p=np.nan)
 j=[f['cols'].index(term) for f in fs];d=fs[0]['b'][j[0]]-fs[1]['b'][j[1]];se=np.sqrt(np.sum((fs[0]['u'][:,j[0]]-fs[1]['u'][:,j[1]])**2))
 return dict(**meta,valid=True,delta=d,se=se,ratio_OR=np.exp(d),low=np.exp(d-1.96*se),high=np.exp(d+1.96*se),p=2*norm.sf(abs(d/se)))
