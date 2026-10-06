from pathlib import Path
import sys,os,json,hashlib,re
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
# Dependencies are supplied by the project environment.
import numpy as np,pandas as pd,h5py
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
X=pd.read_csv(P/'association_covariate_design.tsv',sep='\t',dtype={'sample_id':str}).set_index('sample_id');ix=X.index
S={'scope':'Inventory all six sources and variable availability; no inferential diet association scan.','primary_n':len(ix),'files':{}}
files=list(RAW.iterdir())
for p in files:
 if not p.is_file():continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
 S['files'][p.name]={'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()}
F=pd.read_csv(RAW/'food_exposure_matrix_stool.tsv',sep='\t',dtype={'sample_id':str},low_memory=False).set_index('sample_id')
B=pd.read_csv(RAW/'analysis_matrix_stool_base.tsv',sep='\t',dtype={'sample_id':str},low_memory=False).set_index('sample_id')
M=pd.read_csv(RAW/'sample_information_10317_matched_4168.tsv',sep='\t',dtype=str,low_memory=False).set_index('sample_name')
N=pd.read_csv(RAW/'vioscreen_micromacro.tsv',sep='\t',dtype={'#SampleID':str,'survey_id':str},low_memory=False).set_index('#SampleID')
for name,d in [('food_exposure_matrix_stool.tsv',F),('analysis_matrix_stool_base.tsv',B),('sample_information_10317_matched_4168.tsv',M),('vioscreen_micromacro.tsv',N)]:
 S['files'][name].update(rows=len(d),columns_including_id=len(d.columns)+1,unique_index=bool(d.index.is_unique),primary_ids_matched=int(ix.isin(d.index).sum()))
print('nutrient mapping',S['files']['vioscreen_micromacro.tsv'],flush=True)
assert F.index.is_unique and B.index.is_unique and M.index.is_unique
assert ix.isin(F.index).all() and ix.isin(B.index).all() and ix.isin(M.index).all()
m=M.loc[ix];f=F.loc[ix]
S['identity']={'primary_people':int(m.host_subject_id.nunique()),'all_people':int(M.host_subject_id.nunique())}
if ix.isin(N.index).all() and N.index.is_unique:
 n=N.loc[ix];S['nutrient_join']='exact_sample_id'
elif 'survey_id' in m and N.survey_id.is_unique and m.survey_id.isin(N.survey_id).all():
 n=N.set_index('survey_id').loc[m.survey_id];n.index=ix;S['nutrient_join']='exact_survey_id'
else:raise ValueError('Nutrient identity unresolved; do not guess or deduplicate arbitrarily.')
if 'survey_id' in n:
 S['identity']['primary_survey_ids']=int(n.survey_id.nunique())
 joined=M[['host_subject_id']].join(N[['survey_id']],how='left',validate='one_to_one')
 S['identity']['all_metadata_sample_ids_matched_to_nutrient_file']=int(joined.survey_id.notna().sum())
 S['identity']['people_with_multiple_survey_ids']=int(joined.groupby('host_subject_id').survey_id.nunique().gt(1).sum())
missing={'','not provided','not applicable','missing','unspecified','nan','unknown','not collected','restricted access','na','n/a'}
def audit(d,kind):
 rows=[]
 for c in d:
  s=d[c];v=pd.to_numeric(s,errors='coerce');valid=v.notna()&np.isfinite(v)
  rows.append({'kind':kind,'feature':c,'n':len(s),'valid_numeric':int(valid.sum()),'nonzero':int(v.loc[valid].ne(0).sum()),'negative':int(v.loc[valid].lt(0).sum()),'nunique_numeric':int(v.loc[valid].nunique()),'min':v.min(),'median':v.median(),'max':v.max(),'numeric_missing':int((~valid).sum())})
 return pd.DataFrame(rows)
foods=[c for c in f if c.endswith('__freq_year')];doses=[c for c in f if c.endswith('__dose_g_day')]
fa=audit(f[foods+doses],'food');fa.to_csv(O/'食物频次与摄入量_字段审计.tsv',sep='\t',index=False)
na=audit(n.drop(columns=['survey_id'],errors='ignore'),'nutrient');na.to_csv(O/'营养指标_字段审计.tsv',sep='\t',index=False)
S['food']={'frequency_fields':len(foods),'dose_fields':len(doses),'all_2748_complete_frequency':int(fa.loc[fa.feature.isin(foods),'numeric_missing'].eq(0).sum()),'frequency_nonconstant':int(fa.loc[fa.feature.isin(foods),'nunique_numeric'].gt(1).sum()),'frequency_nonzero_ge5pct':int(((fa.feature.isin(foods))&(fa.nonzero.ge(.05*len(ix)))&fa.valid_numeric.ge(.9*len(ix))).sum()),'negative_entries':int(fa.negative.sum())}
flags=[c for c in f if c.endswith('__inconsistent')]
flagmat=f[flags].apply(pd.to_numeric,errors='coerce')
S['food'].update(inconsistent_cells=int(flagmat.eq(1).sum().sum()),persons_any_inconsistent=int(flagmat.eq(1).any(axis=1).sum()),inconsistent_foods=int(flagmat.eq(1).any().sum()))
S['nutrients']={'numeric_fields':len(na),'complete_fields':int(na.numeric_missing.eq(0).sum()),'nonconstant_fields':int(na.nunique_numeric.gt(1).sum()),'negative_entries':int(na.negative.sum()),'energy_quantiles':pd.to_numeric(n['Energy_in_kcal']).quantile([0,.01,.25,.5,.75,.99,1]).to_dict()}
S['nutrients']['exact_duplicate_columns']=int(n.drop(columns='survey_id',errors='ignore').T.duplicated().sum())
nr=[]
for c in m:
 s=m[c].astype('string').str.strip();ok=s.notna()&~s.str.lower().isin(missing)
 vc=s.loc[ok].value_counts()
 nr.append({'feature':c,'valid_n':int(ok.sum()),'distinct_nonmissing':len(vc),'largest_category_n':int(vc.iloc[0]) if len(vc) else 0,'diet_or_health':bool(re.search('diet|ferment|yogurt|probiotic|milk|cheese|fiber|veget|fruit|antibiotic|ibd|ibs|bowel|medic|acid|pills|smok|alcohol|exercise',c,re.I))})
ma=pd.DataFrame(nr);ma.to_csv(O/'全部问卷字段_可用性.tsv',sep='\t',index=False)
S['metadata']={'fields':len(ma),'entirely_missing':int(ma.valid_n.eq(0).sum()),'nonconstant_with_90pct_available':int((ma.valid_n.ge(.9*len(ix))&ma.distinct_nonmissing.gt(1)).sum())}
ma.loc[ma.diet_or_health&ma.distinct_nonmissing.gt(1)].to_csv(O/'饮食与健康问卷_可用性.tsv',sep='\t',index=False)
tax=pd.read_csv(RAW/'4168_annotated_feature_table.tsv',sep='\t',index_col=0)
S['files']['4168_annotated_feature_table.tsv'].update(taxa_rows=len(tax),sample_columns=len(tax.columns),primary_ids_matched=int(ix.isin(tax.columns).sum()))
S['taxonomy']={'total_reads_primary_quantiles':tax.loc[:,ix].sum().quantile([0,.25,.5,.75,1]).to_dict(),'base_total_exactly_matches':bool(np.array_equal(tax.loc[:,B.index].sum().to_numpy(),B.total_reads.to_numpy()))}
with h5py.File(RAW/'4168_ASV_table.biom','r') as h:
 def dec(a):return [x.decode() if isinstance(x,bytes) else str(x) for x in a]
 sid=dec(h['sample/ids'][:]);aid=dec(h['observation/ids'][:]);data=h['sample/matrix/data'][:];indptr=h['sample/matrix/indptr'][:]
 sums=np.array([data[indptr[i]:indptr[i+1]].sum() for i in range(len(sid))]);ss=pd.Series(sums,index=sid)
 S['biom']={'asvs':len(aid),'samples':len(sid),'primary_ids_matched':int(ix.isin(sid).sum()),'asv_id_lengths':pd.Series([len(a) for a in aid]).value_counts().to_dict(),'observation_groups':list(h['observation'].keys()),'sample_groups':list(h['sample'].keys()),'total_count_matches_annotated':bool(np.array_equal(ss.loc[tax.columns].to_numpy(),tax.sum().to_numpy())),'nonzero_entries':len(data)}
f[foods+doses+flags].to_pickle(O/'_food_primary.pkl');n.to_pickle(O/'_nutrient_primary.pkl')
(O/'全部源数据核查.json').write_text(json.dumps(S,ensure_ascii=False,indent=2,default=lambda v:int(v) if isinstance(v,np.integer) else float(v)),encoding='utf-8')
print(json.dumps(S,ensure_ascii=False,indent=2,default=lambda v:int(v) if isinstance(v,np.integer) else float(v)),flush=True)
