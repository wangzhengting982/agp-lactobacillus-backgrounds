from common import *
# Declared after diagnostic identification, before these complete-BMI fits. Retain original failures.
base=read(RAW/'analysis_matrix_stool_base.tsv').loc[X.index];F=pd.read_pickle(D/'_food_primary.pkl').loc[X.index]
notice='原男性及酸奶零摄入分层的部分模型中，BMI缺失指示系数趋向负无穷，其稳健标准误未通过数值核对。原设定保留为未通过核验；追加两分层的BMI非缺失分析，两目标采用相同人群，原变量尺度不变。'
(O/'分层数值问题与追加规则.txt').write_text(notice,encoding='utf-8');rows=[];con=[]
for spec,mask in [('male_BMI_observed',base.sex_clean.eq('male')),('no_yogurt_BMI_observed',F.yogurt_all_types_except_frozen__freq_year.eq(0))]:
 mask=mask&X.bmi_missing.eq(0);d0=XD.loc[mask];d0=d0.loc[:,d0.nunique().gt(1)|d0.columns.to_series(index=d0.columns).eq('const')]
 for g,c in zip(BG,BC):
  d=d0.assign(background=Z.loc[mask,c]);fs=[]
  for t in TT:
   f,m=glm(d,Y.loc[d.index,t].ge(1));fs.append(f);rows.append(effect(f,'background',dict(spec=spec,background=g,target=t,**m)))
  con.append(contrast(fs,'background',dict(spec=spec,background=g,n=len(d))))
rr=pd.DataFrame(con)
for spec in rr.spec.unique():
 k=rr.spec.eq(spec);rr.loc[k,'q_BH_3']=bh(rr.loc[k]).q_BH.to_numpy()
save(rows,'R04_BMI非缺失分层关联');save(rr,'R05_BMI非缺失分层属间差异');print(rr.to_string(index=False),flush=True)
