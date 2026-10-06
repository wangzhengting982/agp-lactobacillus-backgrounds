from common import *
base=read(RAW/'analysis_matrix_stool_base.tsv').loc[X.index]
f=pd.read_pickle(D/'_food_primary.pkl').loc[X.index];yy0=f['yogurt_all_types_except_frozen__freq_year'].eq(0)
tax=pd.read_csv(RAW/'4168_annotated_feature_table.tsv',sep='\t',index_col=0).loc[:,X.index]
ge=tax.index.to_series().str.extract(r'(?:^|;)g__([^;]+)',expand=False);fa=tax.index.to_series().str.extract(r'(?:^|;)f__([^;]+)',expand=False).fillna('unclassified_family')
ok=tax.index.str.startswith('d__Bacteria;')&~tax.index.str.contains(r'(?:^|;)o__Lactobacillales(?:;|$)|(?:^|;)f__Bifidobacteriaceae(?:;|$)')
known=ge.notna()&~ge.fillna('').str.contains('unclassified|uncultured',case=False)
bg=tax.loc[ok&known].groupby((fa+'|'+ge)[ok&known]).sum().T;bg.columns=['microbe::'+c for c in bg]
den=bg.sum(axis=1);versions={}
for pc in [5e-6,5e-4]:
 lg=np.log(bg[Z.columns].div(den,axis=0)+pc);z=lg.sub(lg.mean(axis=1),axis=0);versions['pseudo_'+str(pc)]=z[BC].apply(scale)
versions['background_presence']=bg[BC].ge(1).astype(float)
den2=den-bg[[c for c in bg if c.startswith('microbe::Peptoniphilaceae|')]].sum(axis=1)
versions['reference_excludes_family']=np.log(bg[BC].div(den2,axis=0)+.00005).apply(scale)
bow=base[['bowel_movement_quality','bowel_movement_frequency']].fillna('unknown').apply(lambda s:s.astype(str).str.strip().str.lower());bd=pd.get_dummies(bow,drop_first=True,dtype=float);bd.columns=['bowel_'+str(i) for i in range(len(bd.columns))]
conditions={'bowel_adjustment':pd.Series(True,index=X.index),'female':base.sex_clean.eq('female'),'male':base.sex_clean.eq('male'),'United_States':base.country_raw.eq('United States'),'United_Kingdom':base.country_raw.eq('United Kingdom'),'no_antibiotic_no_IBD_IBS':base.antibiotic_recent.eq(0)&base.ibd_binary.eq(0)&base.ibs_binary.eq(0),'no_yogurt':yy0,'no_probiotic':base.probiotic_frequency_raw.eq('Never')}
rows=[];contr=[]
for spec in list(versions)+list(conditions):
 mask=conditions.get(spec,pd.Series(True,index=X.index));z=versions.get(spec,Z[BC]);xx=XD.copy()
 if spec=='bowel_adjustment':xx=pd.concat([xx,bd],axis=1)
 xx=xx.loc[mask];xx=xx.loc[:,xx.nunique().gt(1)|xx.columns.to_series(index=xx.columns).eq('const')]
 for g,c in zip(BG,BC):
  d=xx.assign(background=z.loc[mask,c]);fs=[]
  for t in TT:
   fit,m=glm(d,Y.loc[d.index,t].ge(1));fs.append(fit);rows.append(effect(fit,'background',dict(spec=spec,background=g,target=t,**m)))
  contr.append(contrast(fs,'background',dict(spec=spec,background=g,n=len(d))))
rr=pd.DataFrame(rows);cc=pd.DataFrame(contr)
for spec in cc.spec.unique():
 k=cc.spec.eq(spec);cc.loc[k,'q_BH_3']=bh(cc.loc[k]).q_BH.to_numpy()
save(rr,'R01_核心关联全部配套复核');save(cc,'R02_核心差异全部配套复核')
ir=[]
for g,c in zip(BG,BC):
 d=XD.assign(background=Z[c]);d['sex_interaction']=X.sex_male*Z[c]
 for t in TT:
  fit,m=glm(d,Y[t].ge(1));ir.append(effect(fit,'sex_interaction',dict(background=g,target=t,**m)))
save(bh(ir),'R03_性别背景交互');print('VALID',rr.valid.value_counts().to_dict(),flush=True);print(cc.to_string(index=False),flush=True);print(bh(ir).to_string(index=False),flush=True)
