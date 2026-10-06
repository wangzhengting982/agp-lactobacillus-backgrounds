from pathlib import Path
import sys,os,json,warnings,re
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
# Dependencies are supplied by the project environment.
import pandas as pd,numpy as np,statsmodels.api as sm
from sklearn.decomposition import PCA
from scipy.stats import norm
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

P=A;O=D
read=lambda p:pd.read_csv(P/p,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
X=read('association_covariate_design.tsv');Z=read('association_standardized_clr.tsv').loc[X.index];Y=read('association_target_counts.tsv').loc[X.index]
old=pd.read_csv(FROZEN/'single_genus_all_associations.tsv',sep='\t').set_index(['feature','taxon'])
F=pd.read_pickle(O/'_food_primary.pkl').loc[X.index];N=pd.read_pickle(O/'_nutrient_primary.pkl').loc[X.index]
E=pd.to_numeric(N.Energy_in_kcal);summary={'scope':'Exploratory dietary adjustment of three previously selected contrasts; not an independent replication.','preprocessing':{},'fits':[],'energy_restriction':{}}
f=F[[c for c in F if c.endswith('__freq_year')]].copy()
for c in f:
 flag=c.replace('__freq_year','__inconsistent')
 if flag in F:f.loc[F[flag].eq(1),c]=np.nan
n=N.drop(columns=['survey_id'],errors='ignore').apply(pd.to_numeric,errors='coerce')
nc=[c for c in n if '_in_' in c and not c.startswith('Energy_') and '/' not in c]
excluded=[c for c in n if c not in nc]
n=n[nc].copy()
summary['nutrient_excluded_units_or_energy']=excluded
def pcs(d,label,density=False):
 d=d.replace([np.inf,-np.inf],np.nan)
 include=(d.notna().mean().ge(.9)&d.gt(0).mean().ge(.05)&d.nunique().gt(1))
 omitted=d.columns[~include].tolist();d=d.loc[:,include]
 dup=d.T.duplicated();omitted+=d.columns[dup].tolist();d=d.loc[:,~dup]
 d=d.fillna(d.median())
 if density:
  d=d.div(E,axis=0)*1000;positive_scale=d.where(d.gt(0)).median();d=np.log1p(d.div(positive_scale))
 else:d=np.log1p(d)
 lo=d.quantile(.005);hi=d.quantile(.995);d=d.clip(lo,hi,axis=1)
 sd=d.std(ddof=0);good=sd.gt(1e-12);omitted+=d.columns[~good].tolist();d=d.loc[:,good];sd=sd[good]
 a=(d-d.mean())/sd
 p=PCA(n_components=20,svd_solver='full');scores=p.fit_transform(a)
 sc=pd.DataFrame(scores,index=X.index,columns=[label+str(i+1) for i in range(20)])
 sc=sc/sc.std(ddof=0)
 summary['preprocessing'][label]={'used_columns':a.columns.tolist(),'omitted_columns':omitted,'n_features':len(a.columns),'explained_variance_cumulative':{k:float(p.explained_variance_ratio_[:k].sum()) for k in [5,10,20]},'finite':bool(np.isfinite(a).all().all())}
 pd.DataFrame(p.components_.T,index=a.columns,columns=sc.columns).rename_axis('feature').to_csv(O/(label+'_全部载荷.tsv'),sep='\t')
 sc.rename_axis('sample_id').to_csv(O/(label+'_得分_内部复现.tsv'),sep='\t')
 return sc
fp=pcs(f,'food_PC');np_=pcs(n,'nutrient_PC',True)
energykeep=E.between(500,5000)
summary['energy_restriction']={'retained':int(energykeep.sum()),'excluded':int((~energykeep).sum()),'below_500':int(E.lt(500).sum()),'above_5000':int(E.gt(5000).sum()),'events_retained':Y.loc[energykeep,['Lacticaseibacillus','Lactobacillus']].ge(1).sum().to_dict()}
empty=pd.DataFrame(index=X.index)
specs={'original':empty,'food10':fp.iloc[:,:10],'nutrient10':np_.iloc[:,:10],'food10_nutrient10':pd.concat([fp.iloc[:,:10],np_.iloc[:,:10]],axis=1),'food5_nutrient5':pd.concat([fp.iloc[:,:5],np_.iloc[:,:5]],axis=1),'food20_nutrient20':pd.concat([fp,np_],axis=1),'energy_restricted_original':empty,'energy_restricted_food10_nutrient10':pd.concat([fp.iloc[:,:10],np_.iloc[:,:10]],axis=1)}
genera=['Anaerococcus','Finegoldia','Peptoniphilus_A'];targets=['Lacticaseibacillus','Lactobacillus'];rows=[];contrasts=[]
for spec,extra in specs.items():
 mask=energykeep if spec.startswith('energy_restricted') else pd.Series(True,index=X.index)
 for g in genera:
  feat='microbe::Peptoniphilaceae|'+g;d=pd.concat([X,extra,Z[[feat]]],axis=1).loc[mask];a=d.to_numpy(dtype=float)
  assert np.linalg.matrix_rank(a)==a.shape[1]
  u=[];bs=[]
  for t in targets:
   yy=Y.loc[d.index,t].ge(1).astype(float).to_numpy()
   with warnings.catch_warnings(record=True) as ws:
    warnings.simplefilter('always');fit=sm.GLM(yy,a,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0')
   valid=bool(fit.converged and np.isfinite(fit.params).all() and np.isfinite(fit.bse).all() and not ws)
   assert valid, (spec,g,t,[str(w.message) for w in ws])
   mu=np.asarray(fit.fittedvalues);info=a.T@(a*(mu*(1-mu))[:,None]);infl=(a*(yy-mu)[:,None])@np.linalg.inv(info)
   b=float(fit.params[-1]);se=float(fit.bse[-1]);assert abs(np.sqrt(np.sum(infl[:,-1]**2))-se)<1e-7
   u.append(infl[:,-1]);bs.append(b)
   row=dict(spec=spec,background=g,target=t,n=len(d),events=int(yy.sum()),parameters=a.shape[1],beta=b,se=se,OR=float(np.exp(b)),low=float(np.exp(b-1.96*se)),high=float(np.exp(b+1.96*se)),p=float(fit.pvalues[-1]),valid=valid)
   if spec=='original':
    row['original_beta_error']=abs(b-old.loc[(feat,t),'beta_log_OR_per_SD_CLR']);row['original_se_error']=abs(se-old.loc[(feat,t),'SE_HC0'])
    assert row['original_beta_error']<1e-8 and row['original_se_error']<1e-8
   rows.append(row)
  dif=bs[0]-bs[1];sed=float(np.sqrt(np.sum((u[0]-u[1])**2)))
  contrasts.append(dict(spec=spec,background=g,n=len(d),delta_log_OR=dif,se_joint=sed,ratio_OR=float(np.exp(dif)),low=float(np.exp(dif-1.96*sed)),high=float(np.exp(dif+1.96*sed)),p=float(2*norm.sf(abs(dif/sed)))))
 print('completed',spec,flush=True)
rr=pd.DataFrame(rows);cc=pd.DataFrame(contrasts)
for spec in specs:
 k=cc.spec.eq(spec);cc.loc[k,'q_BH_within_selected3']=multipletests(cc.loc[k,'p'],method='fdr_bh')[1]
rr.to_csv(O/'更广泛饮食调整_全部六配对模型.tsv',sep='\t',index=False);cc.to_csv(O/'更广泛饮食调整_全部属间差异.tsv',sep='\t',index=False)
summary['all_models_valid']=bool(rr.valid.all());summary['model_count']=len(rr)
summary['max_original_reproduction_error']=float(rr[['original_beta_error','original_se_error']].max().max())
summary['primary_contrasts']=cc.loc[cc.spec.eq('food10_nutrient10')].to_dict('records')
(O/'饮食预检验汇总.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k!='preprocessing'},ensure_ascii=False,indent=2),flush=True)
print('PCA',json.dumps({k:{kk:vv for kk,vv in v.items() if kk in ['n_features','explained_variance_cumulative']} for k,v in summary['preprocessing'].items()},ensure_ascii=False),flush=True)
print(cc.to_string(index=False),flush=True)
