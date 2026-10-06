from common import *
from statsmodels.discrete.conditional_models import ConditionalLogit
v=read(FROZEN/'纵向分析样本_内部索引.tsv')
yy=read(PREP/'05_输入数据/派生输入/lactobacillales_genus_counts_all_samples.tsv')
z=pd.read_pickle(O/'_all_core_clr.pkl');depth=pd.read_pickle(O/'_all_depth_family.pkl').depth_nonlab
assert len(v)==664 and v.host_subject_id.nunique()==153
v['depth']=scale(depth.loc[v.index]);rows=[]
specs=[(1,s) for s in ['all','drop_most_sampled','first_last','first365days']]+[(3,'all')]
for th,spec in specs:
 for g,c in zip(BG,BC):
  for t in TT:
   a=v.copy()
   if spec=='drop_most_sampled':a=a.loc[a.host_subject_id.ne(a.host_subject_id.value_counts().idxmax())]
   if spec=='first_last':a=a.loc[a.groupby('host_subject_id').cumcount().eq(0)|a.groupby('host_subject_id').cumcount(ascending=False).eq(0)]
   if spec=='first365days':a=a.loc[a.days_since_first.le(365)]
   a['y']=yy.loc[a.index,t].ge(th).astype(int);vary=a.groupby('host_subject_id').y.nunique().eq(2);a=a.loc[a.host_subject_id.isin(vary[vary].index)]
   d=pd.DataFrame({'background':z.loc[a.index,c],'depth':a.depth,'time_year':a.days_since_first/365.25},index=a.index)
   meta=dict(threshold=th,spec=spec,background=g,target=t,people=a.host_subject_id.nunique(),samples=len(a),events=int(a.y.sum()))
   try:
    with warnings.catch_warnings(record=True) as ws:
     warnings.simplefilter('always');f=ConditionalLogit(a.y,d,groups=a.host_subject_id).fit(method='bfgs',maxiter=500,disp=False)
    b=float(f.params['background']);se=float(f.bse['background']);valid=np.isfinite([b,se]).all() and not ws;assert valid,[str(w.message) for w in ws]
    rows.append(dict(**meta,valid=True,beta=b,se=se,OR=np.exp(b),low=np.exp(b-1.96*se),high=np.exp(b+1.96*se),p=2*norm.sf(abs(b/se))))
   except Exception as e:rows.append(dict(**meta,valid=False,error=repr(e),p=np.nan))
 print('within completed',th,spec,flush=True)
rr=pd.DataFrame(rows)
for th,spec in specs:
 k=rr.threshold.eq(th)&rr.spec.eq(spec);rr.loc[k,'q_BH_6']=bh(rr.loc[k]).q_BH.to_numpy()
save(rr,'L01_同人检出变化全部规格');print(rr.to_string(index=False),flush=True)
