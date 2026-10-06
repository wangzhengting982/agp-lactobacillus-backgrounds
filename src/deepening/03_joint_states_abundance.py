from common import *
from scipy.linalg import qr
dep=pd.read_pickle(O/'_all_depth_family.pkl').loc[X.index];family=dep.family_balance
rows=[];con=[];joint=[]
for th in [1,3]:
 for mode in ['three_joint','whole_family']:
  extra=Z[BC].rename(columns=dict(zip(BC,BG))) if mode=='three_joint' else pd.DataFrame({'family_balance':family})
  d=pd.concat([XD,extra],axis=1);fs=[]
  for t in TT:
   f,m=glm(d,Y[t].ge(th));fs.append(f)
   for term in extra:rows.append(effect(f,term,dict(threshold=th,mode=mode,background=term,target=t,**m)))
  for term in extra:con.append(contrast(fs,term,dict(threshold=th,mode=mode,background=term)))
  if all(f is not None for f in fs):
   jj=[list(d.columns).index(c) for c in extra];delta=fs[0]['b'][jj]-fs[1]['b'][jj];u=fs[0]['u'][:,jj]-fs[1]['u'][:,jj];v=u.T@u;s=float(delta@np.linalg.solve(v,delta));joint.append(dict(threshold=th,mode=mode,chi2=s,df=len(jj),p=chi2.sf(s,len(jj))))
rr=pd.DataFrame(rows);cc=pd.DataFrame(con)
for th in [1,3]:
 for mode in ['three_joint','whole_family']:
  k=rr.threshold.eq(th)&rr['mode'].eq(mode);rr.loc[k,'q_BH']=bh(rr.loc[k]).q_BH.to_numpy()
  k=cc.threshold.eq(th)&cc['mode'].eq(mode);cc.loc[k,'q_BH']=bh(cc.loc[k]).q_BH.to_numpy()
save(rr,'J01_联合与菌科关联');save(cc,'J02_联合与菌科直接差异');save(joint,'J03_联合整体差异')
print('JOINT',cc.to_string(index=False),flush=True)
# Four outcomes retain the complete cohort. statsmodels orders outcome-major coefficient blocks.
state_rows=[];prob=[];state_diag=[]
for th in [1,3]:
 y=Y[TT[0]].ge(th).astype(int)+2*Y[TT[1]].ge(th).astype(int)
 for g,z in [(g,Z[c]) for g,c in zip(BG,BC)]+[('family_balance',family)]:
  d=XD.assign(background=z);a=d.to_numpy(float)
  try:
   with warnings.catch_warnings(record=True) as ws:
    warnings.simplefilter('always');f=sm.MNLogit(y.to_numpy(),a).fit(method='newton',maxiter=200,disp=False,cov_type='HC0')
   valid=bool(f.mle_retvals['converged'] and np.isfinite(f.params).all() and np.isfinite(f.cov_params()).all() and not ws)
   assert valid,[str(w.message) for w in ws]
   k=a.shape[1];params=np.asarray(f.params);v=np.asarray(f.cov_params());b=params[-1,:]
   for s1,s2 in [(2,1),(3,1),(3,2)]:
    q1=(s1-1)*k+k-1;q2=(s2-1)*k+k-1;delta=b[s1-1]-b[s2-1];se=np.sqrt(v[q1,q1]+v[q2,q2]-2*v[q1,q2]);state_rows.append(dict(threshold=th,background=g,state_A=s1,state_B=s2,n=len(d),valid=True,beta=delta,se=se,RRR=np.exp(delta),low=np.exp(delta-1.96*se),high=np.exp(delta+1.96*se),p=2*norm.sf(abs(delta/se))))
   for label,level in [('p10',float(z.quantile(.1))),('p90',float(z.quantile(.9)))]:
    ap=a.copy();ap[:,-1]=level;pr=np.asarray(f.predict(ap));assert np.allclose(pr.sum(axis=1),1)
    for s in range(4):
     # Gradient of sample-averaged multinomial probability, outcome-major order.
     grad=np.concatenate([(ap*(pr[:,s]*((1 if s==j else 0)-pr[:,j]))[:,None]).mean(axis=0) for j in [1,2,3]])
     se=np.sqrt(grad@v@grad);mean=pr[:,s].mean();prob.append(dict(threshold=th,background=g,level=label,standardized_value=level,state=s,probability=mean,low=max(0,mean-1.96*se),high=min(1,mean+1.96*se)))
   state_diag.append(dict(threshold=th,background=g,valid=True,n=len(y),counts=json.dumps(y.value_counts().to_dict()),parameters=3*k))
  except Exception as e:
   state_diag.append(dict(threshold=th,background=g,valid=False,error=repr(e)))
   for s1,s2 in [(2,1),(3,1),(3,2)]:state_rows.append(dict(threshold=th,background=g,state_A=s1,state_B=s2,valid=False,p=np.nan))
rr=pd.DataFrame(state_rows)
for th in [1,3]:
 for fam in [False,True]:
  k=rr.threshold.eq(th)&rr.background.eq('family_balance').eq(fam);rr.loc[k,'q_BH']=bh(rr.loc[k]).q_BH.to_numpy()
save(rr,'S01_四种检出状态直接比较');save(prob,'S02_四状态标准化概率');save(state_diag,'S03_多项模型诊断');print('STATES',rr.to_string(index=False),flush=True)
# Positive abundance: selected-positive conditional estimand, no absolute-abundance claim.
ar=[];di=[]
for th in [1,3]:
 for g,c in zip(BG,BC):
  for t in TT:
   ix=Y.index[Y[t].ge(th)];d=XD.assign(background=Z[c]).loc[ix];yy=np.log(Y.loc[ix,t]/dep.loc[ix,'named_depth'])
   # Drop constant-zero columns and use QR if required, preserving background as last column.
   d=d.loc[:,d.nunique().gt(1)|d.columns.to_series(index=d.columns).eq('const')]
   try:
    a=d.to_numpy(float);f0=sm.OLS(yy,a).fit();h=f0.get_influence().hat_matrix_diag;keep=h<1-1e-10
    removed=int((~keep).sum());d=d.loc[keep];yy=yy.loc[d.index];d=d.loc[:,d.nunique().gt(1)|d.columns.to_series(index=d.columns).eq('const')];a=d.to_numpy(float)
    if np.linalg.matrix_rank(a)<a.shape[1]:raise ValueError('rank deficiency remains')
    f=sm.OLS(yy,a).fit(cov_type='HC3');j=list(d.columns).index('background');b=float(f.params.iloc[j]);se=float(f.bse.iloc[j]);valid=np.isfinite([b,se]).all();assert valid
    # Singleton removal leaves estimable background coefficient unchanged.
    assert abs(b-f0.params.iloc[-1])<1e-7
    ar.append(dict(threshold=th,background=g,target=t,eligible=len(ix),n=len(d),singletons_removed=removed,valid=True,beta=b,se=se,ratio_geomean=np.exp(b),low=np.exp(b-1.96*se),high=np.exp(b+1.96*se),p=2*norm.sf(abs(b/se))))
   except Exception as e:ar.append(dict(threshold=th,background=g,target=t,eligible=len(ix),valid=False,error=repr(e),p=np.nan))
rr=pd.DataFrame(ar)
for th in [1,3]:
 k=rr.threshold.eq(th);rr.loc[k,'q_BH_6']=bh(rr.loc[k]).q_BH.to_numpy()
save(rr,'A01_检出后相对丰度');print('ABUNDANCE',rr.to_string(index=False),flush=True)
