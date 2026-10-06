from common import *
checks=[]
def ck(name,value,detail=None):checks.append(dict(check=name,passed=bool(value),detail=detail));assert value,(name,detail)
# An independent outcome-reference change validates multinomial coefficient and covariance indexing.
y=Y[TT[0]].ge(1).astype(int)+2*Y[TT[1]].ge(1).astype(int);d=XD.assign(background=Z[BC[0]])
f=sm.MNLogit(y.map({0:1,1:0,2:2,3:3}).to_numpy(),d.to_numpy()).fit(method='newton',maxiter=200,disp=False,cov_type='HC0')
s=pd.read_csv(O/'S01_四种检出状态直接比较.tsv',sep='\t');ref=s.query('threshold==1 and background=="Anaerococcus" and state_A==2 and state_B==1').iloc[0]
ck('多项模型重新指定参照后系数一致',abs(f.params[-1,1]-ref.beta)<1e-8)
ck('多项模型重新指定参照后稳健标准误一致',abs(f.bse[-1,1]-ref.se)<1e-8)
# Independent stacked participant-clustered GLM checks the three-background difference covariance.
extra=Z[BC].rename(columns=dict(zip(BC,BG)));d=pd.concat([XD,extra],axis=1);a=d.to_numpy();k=a.shape[1];aa=np.block([[a,np.zeros_like(a)],[np.zeros_like(a),a]]);yy=np.r_[Y[TT[0]].ge(1),Y[TT[1]].ge(1)].astype(float)
f=sm.GLM(yy,aa,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='cluster',cov_kwds={'groups':np.tile(np.arange(len(X)),2),'use_correction':False})
v=np.asarray(f.cov_params());params=np.asarray(f.params);C=np.zeros((3,2*k))
for i in range(3):C[i,k-3+i]=1;C[i,2*k-3+i]=-1
de=C@params;vv=C@v@C.T;stat=float(de@np.linalg.solve(vv,de));j=pd.read_csv(O/'J03_联合整体差异.tsv',sep='\t').query('threshold==1 and mode=="three_joint"').iloc[0]
ck('堆叠聚类模型独立复核联合统计量',abs(stat-j.chi2)<1e-6,{'stacked':stat,'separate':j.chi2})
for file,count in [('D02_全部72项摄入关联',72),('D03_全部180项摄入属间比较',180),('D04_12项摄入六属整体比较',12),('D05_18项摄入背景交互及阈值复核',36),('S01_四种检出状态直接比较',24),('A01_检出后相对丰度',12),('L01_同人检出变化全部规格',30),('V02_ASV完整关联',26),('V03_背景ASV属间差异',7)]:
 df=pd.read_csv(O/(file+'.tsv'),sep='\t');ck(file+'完整条目',len(df)==count);ck(file+'全部模型有效',df.valid.all() if 'valid' in df else True)
ck('三背景联合估计不冒充三个独立差异',pd.read_csv(O/'J02_联合与菌科直接差异.tsv',sep='\t').query('threshold==1 and mode=="three_joint"').q_BH.ge(.05).all())
(O/'统计独立核验.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8');print('statistical checks',len(checks),'passed',flush=True)
