from common import *
F=pd.read_pickle(O/'_diet12.pkl');xx=pd.concat([X.drop(columns=['yogurt','plant_ferment','probiotic','fiber_density']),F],axis=1);rows=[];fits={};diagnostics=[]
for t in T:
 f,m=glm(xx,Y[t].ge(1));fits[t]=f;diagnostics.append(dict(target=t,**m))
 for c in F:rows.append(effect(f,c,dict(target=t,exposure=c,**m)))
rr=bh(rows);save(rr,'D02_全部72项摄入关联');cc=[];gg=[]
for c in F:
 for ta,tb in itertools.combinations(T,2):cc.append(contrast([fits[ta],fits[tb]],c,dict(exposure=c,target_A=ta,target_B=tb)))
 if all(fits[t] is not None for t in T):
  j=list(xx.columns).index(c);b=np.array([fits[t]['b'][j] for t in T]);u=np.column_stack([fits[t]['u'][:,j] for t in T]);C=np.column_stack([-np.ones(5),np.eye(5)]);delta=C@b;vv=C@(u.T@u)@C.T;s=float(delta@np.linalg.solve(vv,delta));gg.append(dict(exposure=c,chi2=s,df=5,p=chi2.sf(s,5)))
 else:gg.append(dict(exposure=c,p=np.nan))
save(bh(cc),'D03_全部180项摄入属间比较');save(bh(gg),'D04_12项摄入六属整体比较');save(diagnostics,'D02_模型诊断')
rows=[];pred=[]
for th in [1,3]:
 for intake in ['yogurt','plant_ferment','probiotic']:
  for g,c in zip(BG,BC):
   d=XD.assign(background=Z[c]);d['interaction']=X[intake]*Z[c]
   for t in TT:
    f,m=glm(d,Y[t].ge(th));rows.append(effect(f,'interaction',dict(threshold=th,intake=intake,background=g,target=t,**m)))
rr=pd.DataFrame(rows)
for th in [1,3]:
 ix=rr.threshold.eq(th);rr.loc[ix,'q_BH_18']=bh(rr.loc[ix],out='q').q.to_numpy()
save(rr,'D05_18项摄入背景交互及阈值复核')
print('diet',pd.DataFrame(rows).valid.value_counts().to_dict(),flush=True);print('diet associations',pd.read_csv(O/'D02_全部72项摄入关联.tsv',sep='\t').query('q_BH<0.05').to_string(index=False),flush=True);print('interactions',rr.to_string(index=False),flush=True)
