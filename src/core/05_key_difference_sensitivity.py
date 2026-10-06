from pathlib import Path
from core_paths import DATA, OUT, PRIMARY, TAXONOMY_INPUT
import os,sys
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd,statsmodels.api as sm,json
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
R=DATA;O=OUT;P=(DATA/'输入/03_分析与结果/association_models');Q=TAXONOMY_INPUT
scan=pd.read_csv(O/'表S22f_全部2340项背景属系数比较.tsv',sep='\t')
sel=scan.sort_values(['q_BH_2340','p','feature']).head(3)
assert sel.taxon_A.eq('Lacticaseibacillus').all() and sel.taxon_B.eq('Lactobacillus').all()
sel.to_csv(O/'展示配对选择规则_完整家族q最小3项.tsv',sep='\t',index=False)
selected=sel.feature.tolist();taxa=['Lacticaseibacillus','Lactobacillus'];rows=[];co=[]
for scenario in ['original_reads1','original_reads3','R07_reads1','R06_reads1']:
    if scenario.startswith('original'):
        xp=PRIMARY/'association_covariate_design.tsv';zp=PRIMARY/'association_standardized_clr.tsv';yp=PRIMARY/'association_target_counts.tsv'
    else:
        p=Q/scenario.split('_')[0];xp=p/'covariate_design.tsv';zp=p/'standardized_clr.tsv.gz';yp=p/'target_counts.tsv'
    rd=lambda p:pd.read_csv(p,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
    x=rd(xp);z=rd(zp).loc[x.index];y=rd(yp).loc[x.index,taxa].ge(3 if scenario.endswith('reads3') else 1).astype(float)
    for feature in selected:
        genus=feature.split('|')[-1];matched=[f for f in z.columns if f.split('|')[-1]==genus]
        if len(matched)!=1:
            rows.append(dict(scenario=scenario,feature=feature,status='not uniquely matched',p=np.nan));continue
        a=np.column_stack([x,z[matched[0]]]);bs=[];us=[]
        for t in taxa:
            f=sm.GLM(y[t].to_numpy(),a,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0');assert f.converged
            p=np.asarray(f.fittedvalues);u=(a*(y[t].to_numpy()-p)[:,None])@np.linalg.inv(a.T@(a*(p*(1-p))[:,None]))
            bs.append(f.params[-1]);us.append(u[:,-1]);se=f.bse[-1]
            co.append(dict(scenario=scenario,feature=feature,matched_feature=matched[0],taxon=t,OR=np.exp(f.params[-1]),CI_low=np.exp(f.params[-1]-1.96*se),CI_high=np.exp(f.params[-1]+1.96*se)))
        delta=bs[0]-bs[1];se=np.sqrt(np.sum((us[0]-us[1])**2))
        rows.append(dict(scenario=scenario,feature=feature,matched_feature=matched[0],status='estimated',beta_A=bs[0],beta_B=bs[1],delta_log_OR=delta,SE_joint=se,ratio_OR=np.exp(delta),CI_low=np.exp(delta-1.96*se),CI_high=np.exp(delta+1.96*se),p=float(2*norm.sf(abs(delta/se)))))
out=pd.DataFrame(rows)
for scenario in out.scenario.unique():
    ii=out.scenario.eq(scenario)
    if scenario=='original_reads1':
        out.loc[ii,'q_original_2340']=out.loc[ii,'feature'].map(sel.set_index('feature').q_BH_2340)
    else:
        out.loc[ii,'q_BH_3']=multipletests(out.loc[ii,'p'].fillna(1),method='fdr_bh')[1]
out.to_csv(O/'表S24h_三项突出差异的阈值与分类复核.tsv',sep='\t',index=False)
pd.DataFrame(co).to_csv(O/'表S24i_三背景两目标属效应量.tsv',sep='\t',index=False)
print(out.to_string(index=False))
