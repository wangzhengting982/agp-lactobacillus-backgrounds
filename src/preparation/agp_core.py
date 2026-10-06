from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import json, hashlib, warnings, sys, time
import numpy as np
import pandas as pd
from scipy.special import logit
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import log_loss, roc_auc_score, average_precision_score, brier_score_loss
from sklearn.exceptions import ConvergenceWarning
from joblib import Parallel, delayed, parallel_config
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'03_分析与结果/agp_models'
PREV=ROOT/'05_输入数据/派生输入'
R=ROOT/'05_输入数据'
TAXA=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
H=['age','bmi','antibiotic','ibd','ibs','year','depth_nonlab','energy','fiber_density']
D=['yogurt','plant_ferment','probiotic']
S=['nonlab_shannon','nonlab_logrichness']
CATS=['sex','country']

class MicrobeCLR(TransformerMixin, BaseEstimator):
 def __init__(self,min_prevalence=.10,pseudocount=.00005):
  self.min_prevalence=min_prevalence;self.pseudocount=pseudocount
 def fit(self,X,y=None):
  self.feature_names_in_=np.array(X.columns,dtype=object)
  a=np.asarray(X,dtype=float);self.keep_=np.mean(a>0,axis=0)>=self.min_prevalence
  if self.keep_.sum()<2:raise ValueError('Too few eligible background genera')
  return self
 def transform(self,X):
  a=np.asarray(X,dtype=float)
  rel=a[:,self.keep_]/np.maximum(a.sum(axis=1,keepdims=True),1)
  z=np.log(rel+self.pseudocount)
  return z-z.mean(axis=1,keepdims=True)
 def get_feature_names_out(self,input_features=None):return self.feature_names_in_[self.keep_]

def prepare():
 ix=pd.read_csv(PREV/'pilot_analysis_design.tsv',sep='\t',index_col=0,dtype={'sample_id':str}).index
 assert len(ix)==2877 and ix.is_unique
 b=pd.read_csv(R/'02_模型输入/analysis_matrix_stool_base.tsv',sep='\t',dtype={'sample_id':str}).set_index('sample_id').loc[ix]
 m=pd.read_csv(R/'01_原始输入/sample_information_10317_matched_4168.tsv',sep='\t',dtype=str).set_index('sample_name').loc[ix]
 f=pd.read_csv(R/'02_模型输入/food_exposure_matrix_stool.tsv',sep='\t',usecols=['sample_id','yogurt_all_types_except_frozen__freq_year'],dtype={'sample_id':str}).set_index('sample_id').loc[ix]
 n=pd.read_csv(R/'01_原始输入/vioscreen_micromacro.tsv',sep='\t',usecols=['#SampleID','Energy_in_kcal','Total_Dietary_Fiber_in_g'],dtype={'#SampleID':str}).set_index('#SampleID').loc[ix]
 tax=pd.read_csv(R/'01_原始输入/4168_annotated_feature_table.tsv',sep='\t',index_col=0).loc[:,ix]
 forbidden=tax.index.str.contains(r'(?:^|;)o__Lactobacillales(?:;|$)|(?:^|;)f__Bifidobacteriaceae(?:;|$)')
 bacterial=tax.index.str.startswith('d__Bacteria;')
 ok=bacterial&~forbidden
 genera=tax.index.to_series().str.extract(r'(?:^|;)g__([^;]+)',expand=False)
 families=tax.index.to_series().str.extract(r'(?:^|;)f__([^;]+)',expand=False).fillna('unclassified_family')
 labels=families+'|'+genera
 known=genera.notna()&~genera.fillna('').str.contains('unclassified|uncultured',case=False)
 bg=tax.loc[ok&known].groupby(labels[ok&known]).sum().T
 bg.columns=['microbe::'+str(x) for x in bg.columns]
 nonlab_depth=tax.loc[ok].sum(axis=0)
 assert not any(any(t in c for t in TAXA+['Bifidobacterium','Streptococcus','Lactococcus']) for c in bg.columns)
 excl=pd.DataFrame({'taxonomy':tax.index,'excluded_order_or_family':forbidden,'bacterial':bacterial,'named_genus':known.to_numpy(),'used_for_background':(ok&known).to_numpy()})
 excl.to_csv(P/'taxonomy_feature_audit.tsv',sep='\t',index=False)
 x=pd.DataFrame(index=ix)
 x['age']=b.age_years_clean;x['bmi']=b.bmi_clean;x['sex']=b.sex_clean.fillna('unknown').astype(str)
 for old,new in [('antibiotic_recent','antibiotic'),('ibd_binary','ibd'),('ibs_binary','ibs'),('collection_year_clean','year')]:x[new]=b[old]
 x['year']=pd.to_datetime(m.collection_timestamp,format='mixed',errors='coerce').dt.year
 x['country']=np.where(m.country_residence=='United States','US',np.where(m.country_residence=='United Kingdom','UK',np.where(m.country_residence.isna()|m.country_residence.isin(['not provided','not applicable','missing: not provided','missing: not applicable']),'Unknown','Other')))
 x['depth_nonlab']=np.log10(nonlab_depth+1);x['energy']=np.log(n.Energy_in_kcal.where(n.Energy_in_kcal>0));x['fiber_density']=n.Total_Dietary_Fiber_in_g/n.Energy_in_kcal*1000
 om={'Never':0,'Rarely (a few times/month)':1,'Rarely (less than once/week)':1,'Occasionally (1-2 times/week)':2,'Regularly (3-5 times/week)':3,'Daily':4}
 x['yogurt']=np.log1p(f.iloc[:,0]);x['plant_ferment']=m.fermented_plant_frequency.map(om);x['probiotic']=m.probiotic_frequency.map(om)
 x['bowel_quality']=np.select([m.bowel_movement_quality.str.contains('Type 1 and 2',na=False),m.bowel_movement_quality.str.contains('Type 5, 6 and 7',na=False),m.bowel_movement_quality.str.contains('type 3 and 4',case=False,na=False)],['constipated','loose','normal'],default='unknown')
 x['bowel_frequency']=m.bowel_movement_frequency.where(m.bowel_movement_frequency.isin(['Less than one','One','Two','Three','Four','Five or more']),'unknown')
 pp=bg.div(bg.sum(axis=1).replace(0,np.nan),axis=0).fillna(0)
 x[S[0]]=-(pp*np.log(pp.where(pp>0,1))).sum(axis=1);x[S[1]]=np.log1p((bg>0).sum(axis=1))
 x['total_reads_qc']=b.total_reads
 x=x.replace([np.inf,-np.inf],np.nan)
 assert x[D].notna().all().all()
 x=pd.concat([x,bg],axis=1)
 y=pd.read_csv(PREV/'lactobacillales_genus_counts_all_samples.tsv',sep='\t',index_col=0,dtype={'sample_id':str}).loc[ix,TAXA]
 x.to_pickle(P/'raw_predictors.pkl');y.to_csv(P/'target_genus_counts.tsv',sep='\t')
 x[H+D+S+CATS+['bowel_quality','bowel_frequency','total_reads_qc']].to_csv(P/'raw_covariates.tsv',sep='\t')
 bg.to_csv(P/'background_genus_counts.tsv',sep='\t')
 summary={'n':len(x),'background_named_genera':bg.shape[1],'background_genera_global_prevalence_ge10pct':int(((bg>0).mean()>=.1).sum()),'excluded_taxonomy_rows':int(forbidden.sum()),'countries':x.country.value_counts().to_dict(),'events_ge1':y.ge(1).sum().to_dict(),'events_ge3':y.ge(3).sum().to_dict(),'depth_ge5000_n':int(x.total_reads_qc.ge(5000).sum()),'nonlab_depth_zero':int(nonlab_depth.eq(0).sum()),'covariate_missing':x[H+D+CATS].isna().sum().to_dict(),'sample_index_unique':bool(x.index.is_unique)}
 (P/'prepared_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

def model_pipeline(x,model,with_bowel=False):
 num=H.copy();cats=CATS.copy()
 if 'D' in model:num+=D
 if 'S' in model:num+=S
 if with_bowel:cats+=['bowel_quality','bowel_frequency']
 numpipe=Pipeline([('impute',SimpleImputer(strategy='median',add_indicator=True,keep_empty_features=True)),('scale',StandardScaler())])
 blocks=[('covariates',numpipe,num),('categories',OneHotEncoder(handle_unknown='ignore',sparse_output=False),cats)]
 if 'E' in model:
  cols=[c for c in x if c.startswith('microbe::')]
  blocks.append(('ecology',Pipeline([('clr',MicrobeCLR()),('scale',StandardScaler())]),cols))
 return Pipeline([('features',ColumnTransformer(blocks,remainder='drop',sparse_threshold=0)),('logistic',LogisticRegression(C=1,solver='lbfgs',max_iter=1200,tol=1e-6))])

def metrics(y,p):
 p=np.clip(p,1e-10,1-1e-10)
 return {'n':len(y),'events':int(np.sum(y)),'prevalence':float(np.mean(y)),'log_loss':log_loss(y,p,labels=[0,1]),'AUROC':roc_auc_score(y,p),'AP':average_precision_score(y,p),'Brier':brier_score_loss(y,p)}


if __name__=="__main__":
 if len(sys.argv)==1 or sys.argv[1]=="prepare":prepare()
 else:raise SystemExit("Use formal_agp_models.py for all formal model fits.")
