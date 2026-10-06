from pathlib import Path
import sys,json,warnings,itertools
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common import RAW,X,XD,Z,Y,TT,BG,BC,glm,effect,contrast,bh,scale,np,pd,sm,norm,ARCH,RUN
import patsy,h5py
from scipy.sparse import csr_matrix
R=RUN;O=RUN/'sex_results';O.mkdir(parents=True,exist_ok=True)
M=pd.read_csv(O/'_additional_metadata.tsv',sep='\t',dtype={'sample_id':str}).set_index('sample_id').loc[X.index]
drug=pd.DataFrame({'hormonal':M.contraceptive.str.startswith('Yes').astype(float),'contraceptive_missing':M.contraceptive.eq('not provided').astype(float)},index=X.index)
spl=patsy.dmatrix('cr(age, df=4, constraints="center") - 1',{'age':X.age},return_type='dataframe');spl.index=X.index;spl.columns=['age_spline_'+str(i) for i in range(4)]
schemes={'primary':(XD,pd.Series(True,index=X.index)), 'not_pregnant_hormonal_adjusted':(pd.concat([XD,drug],axis=1),M.pregnant.eq('No')), 'no_pregnancy_no_hormonal':(XD,M.pregnant.eq('No')&M.contraceptive.eq('No')), 'age_spline4':(pd.concat([XD.drop(columns='age'),spl],axis=1),pd.Series(True,index=X.index))}
rows=[];difs=[];inter=[];diagn=[];fits_primary={}
def lin(f,v,meta,kind='OR'):
    b=float(v@f['b']);se=float(np.sqrt(v@f['v']@v));return dict(**meta,beta=b,se=se,estimate=np.exp(b),low=np.exp(b-1.96*se),high=np.exp(b+1.96*se),p=2*norm.sf(abs(b/se)))
def diffvec(fs,vecs,meta):
    if any(f is None for f in fs):return dict(**meta,valid=False,p=np.nan)
    b=vecs[0]@fs[0]['b']-vecs[1]@fs[1]['b'];u=fs[0]['u']@vecs[0]-fs[1]['u']@vecs[1];se=np.sqrt(u@u)
    return dict(**meta,valid=True,beta=b,se=se,ratio_OR=np.exp(b),low=np.exp(b-1.96*se),high=np.exp(b+1.96*se),p=2*norm.sf(abs(b/se)))
for th in [1,3]:
 for scheme,(base,mask) in schemes.items():
  for g,c in zip(BG,BC):
   d=base.loc[mask].copy();d=d.loc[:,d.nunique().gt(1)|d.columns.to_series(index=d.columns).eq('const')];d['background']=Z.loc[d.index,c];d['sex_interaction']=X.loc[d.index,'sex_male']*d.background
   fs=[]
   for t in TT:
    f,m=glm(d,Y.loc[d.index,t].ge(th));fs.append(f);diagn.append(dict(threshold=th,scheme=scheme,background=g,target=t,**m))
    if f is not None:
     if th==1 and scheme=='primary':fits_primary[(g,t)]=f
     inter.append(effect(f,'sex_interaction',dict(threshold=th,scheme=scheme,background=g,target=t,**m)))
     for sex in [0,1]:
      v=np.zeros(len(d.columns));v[d.columns.get_loc('background')]=1;v[d.columns.get_loc('sex_interaction')]=sex
      rows.append(lin(f,v,dict(threshold=th,scheme=scheme,background=g,target=t,sex='female' if sex==0 else 'male',n=int(X.loc[d.index,'sex_male'].eq(sex).sum()),events=int(Y.loc[d.index[X.loc[d.index,'sex_male'].eq(sex)],t].ge(th).sum()))))
   for typ in ['female','male','sex_difference']:
    v=np.zeros(len(d.columns))
    if typ!='sex_difference':v[d.columns.get_loc('background')]=1
    if typ!='female':v[d.columns.get_loc('sex_interaction')]=1
    difs.append(diffvec(fs,[v,v],dict(threshold=th,scheme=scheme,background=g,comparison=typ,n=len(d))))
  print('completed',th,scheme,flush=True)
dd=pd.DataFrame(difs);ii=pd.DataFrame(inter)
for th in [1,3]:
 for scheme in schemes:
  for typ in ['female','male','sex_difference']:
   k=dd.threshold.eq(th)&dd.scheme.eq(scheme)&dd.comparison.eq(typ);dd.loc[k,'q_BH_3']=bh(dd.loc[k]).q_BH.to_numpy()
  k=ii.threshold.eq(th)&ii.scheme.eq(scheme);ii.loc[k,'q_BH_6']=bh(ii.loc[k]).q_BH.to_numpy()
pd.DataFrame(rows).to_csv(O/'H01_性别内目标属关联.tsv',sep='\t',index=False);dd.to_csv(O/'H02_属间差异及性别直接比较.tsv',sep='\t',index=False);ii.to_csv(O/'H03_原性别交互复核.tsv',sep='\t',index=False);pd.DataFrame(diagn).to_csv(O/'H04_模型诊断.tsv',sep='\t',index=False)
orig=pd.read_csv(RUN/'results/R03_性别背景交互.tsv',sep='\t');new=ii.query('threshold==1 and scheme=="primary"');comp=orig.merge(new,on=['background','target'],suffixes=('_old','_new'));assert np.allclose(comp.beta_old,comp.beta_new,atol=1e-9) and np.allclose(comp.se_old,comp.se_new,atol=1e-9)
# Independent stacked participant-clustered reconstruction for every primary background.
vcheck=[]
for g in BG:
 f1,f2=[fits_primary[(g,t)] for t in TT];a=f1['d'].to_numpy();n,k=a.shape;stack=np.zeros((2*n,2*k));stack[:n,:k]=a;stack[n:,k:]=a;ys=np.r_[Y[TT[0]].ge(1),Y[TT[1]].ge(1)].astype(float)
 fit=sm.GLM(ys,stack,family=sm.families.Binomial()).fit(maxiter=150,cov_type='cluster',cov_kwds={'groups':np.r_[np.arange(n),np.arange(n)],'use_correction':False})
 vec=np.zeros(2*k);j=f1['cols'].index('sex_interaction');vec[j]=1;vec[k+j]=-1;b=float(vec@fit.params);se=float(np.sqrt(vec@fit.cov_params()@vec));ref=dd.query('threshold==1 and scheme=="primary" and comparison=="sex_difference" and background==@g').iloc[0]
 assert abs(b-ref.beta)<1e-7 and abs(se-ref.se)<1e-7;vcheck.append({'background':g,'beta_error':abs(b-ref.beta),'se_error':abs(se-ref.se)})
(O/'H05_独立统计核验.json').write_text(json.dumps(vcheck,ensure_ascii=False,indent=2),encoding='utf8')
# Fixed target ASVs, assessed without source or strain inference.
with h5py.File(RAW/'4168_ASV_table.biom','r') as h:
 sid=[x.decode() for x in h['sample/ids'][:]];grp=h['observation/matrix'];mat=csr_matrix((grp['data'][:],grp['indices'][:],grp['indptr'][:]),shape=(95567,len(sid)))[:,pd.Index(sid).get_indexer(X.index)]
CY={r:pd.Series(np.asarray(mat[r].toarray()).ravel(),index=X.index) for r in [40974,90943]};ar=[];ac=[];ad=[]
for th in [1,3]:
 for g,c in zip(BG,BC):
  d=XD.assign(background=Z[c]);d['sex_interaction']=X.sex_male*Z[c];fs=[]
  for r in CY:
   yy=CY[r].ge(th);event=[int(yy[X.sex_male.eq(s)].sum()) for s in [0,1]];eligible=min(event)>=10
   f,m=glm(d,yy) if eligible else (None,{'valid':False,'error':'fewer than 10 ASV detections in one sex'})
   fs.append(f);ad.append(dict(threshold=th,background=g,ASV=r,female_events=event[0],male_events=event[1],**m));ar.append(effect(f,'sex_interaction',dict(threshold=th,background=g,ASV=r,female_events=event[0],male_events=event[1],**m)))
  for typ in ['female','male','sex_difference']:
   v=np.zeros(len(d.columns))
   if typ!='sex_difference':v[d.columns.get_loc('background')]=1
   if typ!='female':v[d.columns.get_loc('sex_interaction')]=1
   ac.append(diffvec(fs,[v,v],dict(threshold=th,background=g,comparison=typ,ASV_A=40974,ASV_B=90943)))
ar=pd.DataFrame(ar);ac=pd.DataFrame(ac)
for th in [1,3]:
 k=ar.threshold.eq(th);ar.loc[k,'q_BH_6']=bh(ar.loc[k]).q_BH.to_numpy()
 for typ in ['female','male','sex_difference']:
  k=ac.threshold.eq(th)&ac.comparison.eq(typ);ac.loc[k,'q_BH_3']=bh(ac.loc[k]).q_BH.to_numpy()
ar.to_csv(O/'H06_固定ASV性别交互.tsv',sep='\t',index=False);ac.to_csv(O/'H07_固定ASV直接比较.tsv',sep='\t',index=False);pd.DataFrame(ad).to_csv(O/'H08_ASV模型诊断.tsv',sep='\t',index=False)
print('PRIMARY',dd.query('scheme=="primary" and threshold==1').to_string(index=False),flush=True)
print('SENSITIVITY',dd.query('comparison=="sex_difference"').to_string(index=False),flush=True)
print('ASV',ar.to_string(index=False),flush=True)
