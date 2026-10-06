from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent))
from common import X,XD,Z,Y,TT,BG,BC,np,pd,sm,bh,RUN
R=RUN;O=RUN/'sex_results';O.mkdir(parents=True,exist_ok=True);(R/'修订记录').mkdir(exist_ok=True)
# Declared before fitting: remove exactly one maximum-count positive participant per target;
# retain every background, the frozen design and all six tests. No alternate deletion rules.
(R/'修订记录/高读数影响复核范围.txt').write_text('在v39丰度结果已知后，针对原始注释表各目标属计数最高者，分别剔除1人；并列时按sample_id升序取首人。其余保留原40列饮食扩展协变量、背景CLR及HC3估计。主阈值>=1的2目标×3背景共6项统一BH，不追加其他删样方案。',encoding='utf8')
depth=pd.read_pickle(RUN/'results/_all_depth_family.pkl')
depth['named_depth']=depth['named_depth'].astype(float)
rows=[]
for t in TT:
 eligible=Y.index[Y[t].ge(1)];drop=Y.loc[eligible].sort_index()[t].idxmax();ids=eligible[eligible!=drop]
 for g,c in zip(BG,BC):
  d=XD.loc[ids].copy();d=d.loc[:,d.nunique().gt(1)|d.columns.to_series(index=d.columns).eq('const')];d['background']=Z.loc[ids,c]
  yy=np.log(Y.loc[ids,t]/depth.loc[ids,'named_depth']);f=sm.OLS(yy,d).fit(cov_type='HC3');v=float(f.params['background']);se=float(f.bse['background']);ci=f.conf_int().loc['background']
  valid=np.isfinite(f.params).all() and np.isfinite(f.bse).all() and np.linalg.matrix_rank(d)==d.shape[1]
  rows.append(dict(target=t,background=g,n=len(ids),dropped_id=drop,dropped_count=float(Y.loc[drop,t]),fraction_of_target_reads=float(Y.loc[drop,t]/Y[t].sum()),valid=valid,beta=v,se=se,ratio_geomean=np.exp(v),low=np.exp(ci[0]),high=np.exp(ci[1]),p=f.pvalues['background']))
out=bh(rows);out.to_csv(O/'A02_最高读数样本剔除复核.tsv',sep='\t',index=False);print(out.drop(columns='dropped_id').to_string(index=False),flush=True)
