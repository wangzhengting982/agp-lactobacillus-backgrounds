from pathlib import Path
from core_paths import DATA, OUT, PRIMARY, TAXONOMY_INPUT
import sys
import numpy as np,pandas as pd,itertools
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
R=DATA;O=OUT;P=(DATA/'输入/03_分析与结果/taxonomy_validation/association_sensitivity')
T=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
rows=[]
for release in ['R07','R06']:
    a=np.load(O/'taxonomy_refits'/release/'pc_joint_sandwich.npz');b=a['beta'];v=a['covariance'];assert b.shape==(180,)
    for i,j in itertools.combinations(range(6),2):
        c=np.zeros((10,180));c[:,i*30+20:(i+1)*30]=np.eye(10);c[:,j*30+20:(j+1)*30]=-np.eye(10)
        d=c@b;cv=c@v@c.T;assert np.linalg.matrix_rank(cv)==10;s=float(d@np.linalg.solve(cv,d))
        rows.append(dict(release=release,taxon_A=T[i],taxon_B=T[j],chi2=s,df=10,p=float(chi2.sf(s,10))))
out=pd.DataFrame(rows)
for release in ['R07','R06']:
    idx=out.release.eq(release);out.loc[idx,'q_BH_15']=multipletests(out.loc[idx,'p'],method='fdr_bh')[1];out.loc[idx,'p_Holm_15']=multipletests(out.loc[idx,'p'],method='holm')[1]
out.to_csv(O/'表S24g_重新注释15组整体模式比较.tsv',sep='\t',index=False)
print(out.to_string(index=False))
