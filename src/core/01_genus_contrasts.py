from pathlib import Path
from core_paths import DATA, OUT, PRIMARY, TAXONOMY_INPUT
import os,sys
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd,statsmodels.api as sm,json,hashlib,warnings,itertools
from scipy.stats import chi2,norm
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
R=DATA; O=OUT;O.mkdir(exist_ok=True)
P=(DATA/'输入/03_分析与结果/association_models')
T=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
files={x:(PRIMARY if x.startswith('association_') else P)/x for x in ['association_covariate_design.tsv','association_standardized_clr.tsv','association_target_counts.tsv','association_pc_scores.tsv','single_genus_all_associations.tsv','pc_joint_sandwich.npz']}
hashes={str(p):sha(p) for p in files.values()}
read=lambda n:pd.read_csv(files[n],sep='\t',dtype={'sample_id':str}).set_index('sample_id')
X=read('association_covariate_design.tsv');Z=read('association_standardized_clr.tsv').loc[X.index];C=read('association_target_counts.tsv').loc[X.index,T]
Y=C.ge(1).astype(float)
assert X.shape==(2748,20) and Z.shape==(2748,156)
diag=[];check=[];coeff=[];omnibus=[];contrasts=[]

def fit_joint(a,y,label):
    bet=[];inf=[]
    for t in y.columns:
        with warnings.catch_warnings(record=True) as ws,threadpool_limits(limits=1):
            warnings.simplefilter('always')
            f=sm.GLM(y[t].to_numpy(),a,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0')
        v=np.asarray(f.cov_params());pr=np.asarray(f.fittedvalues)
        influence=(a*(y[t].to_numpy()-pr)[:,None])@np.linalg.inv(a.T@(a*(pr*(1-pr))[:,None]))
        err=float(np.max(np.abs(influence.T@influence-v)))
        valid=bool(f.converged and np.isfinite(v).all() and np.isfinite(f.params).all() and not ws)
        diag.append(dict(analysis=label,taxon=t,n=len(y),events=int(y[t].sum()),parameters=a.shape[1],rank=int(np.linalg.matrix_rank(a)),valid=valid,warnings=' | '.join(str(w.message) for w in ws),cov_reconstruction_error=err))
        assert valid and err<1e-7,(label,t,diag[-1])
        bet.append(f.params);inf.append(influence)
    b=np.concatenate(bet);u=np.column_stack(inf);v=u.T@u
    return b,v,u

def equal_vectors(b,v,groups,k,last):
    mat=np.zeros(((len(groups)-1)*last,len(b)))
    for h,g in enumerate(groups[1:]):
        for j in range(last):
            mat[h*last+j,groups[0]*k+k-last+j]=-1
            mat[h*last+j,g*k+k-last+j]=1
    d=mat@b;cv=mat@v@mat.T;rank=np.linalg.matrix_rank(cv)
    assert rank==len(d)
    stat=float(d@np.linalg.solve(cv,d))
    return dict(chi2=stat,df=int(rank),p=float(chi2.sf(stat,rank)))

base=pd.read_csv(P/'single_genus_all_associations.tsv',sep='\t')
for cutoff in [1,3]:
    y=C.ge(cutoff).astype(float)
    for genus in ['Oscillibacter','Veillonella']:
        feature=next(c for c in Z if c.endswith('|'+genus));a=np.column_stack([X,Z[feature]])
        b,v,u=fit_joint(a,y,f'{genus}_reads{cutoff}');k=a.shape[1]
        omnibus.append(dict(background=genus,cutoff=cutoff,**equal_vectors(b,v,list(range(6)),k,1)))
        for i,t in enumerate(T):
            j=(i+1)*k-1;se=np.sqrt(v[j,j]);beta=b[j]
            coeff.append(dict(background=genus,cutoff=cutoff,taxon=t,beta=beta,se=se,OR=np.exp(beta),low=np.exp(beta-1.96*se),high=np.exp(beta+1.96*se)))
            if cutoff==1:
                old=base[(base.taxon==t)&(base.feature==feature)].iloc[0]
                assert abs(beta-old.beta_log_OR_per_SD_CLR)<1e-8 and abs(se-old.SE_HC0)<1e-8
        for i,j in itertools.combinations(range(6),2):
            ii=(i+1)*k-1;jj=(j+1)*k-1;d=b[ii]-b[jj];se=np.sqrt(v[ii,ii]+v[jj,jj]-2*v[ii,jj])
            contrasts.append(dict(background=genus,cutoff=cutoff,taxon_A=T[i],taxon_B=T[j],beta_A=b[ii],beta_B=b[jj],delta_log_OR=d,se=se,ratio_OR=np.exp(d),CI_low=np.exp(d-1.95996398454*se),CI_high=np.exp(d+1.95996398454*se),p=2*norm.sf(abs(d/se))))
        # Independent stacked GLM verifies all six coefficients and same-person covariance.
        biga=np.zeros((len(y)*6,k*6))
        for i in range(6):biga[i*len(y):(i+1)*len(y),i*k:(i+1)*k]=a
        with threadpool_limits(limits=1):
            stacked=sm.GLM(y.to_numpy().T.ravel(),biga,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='cluster',cov_kwds={'groups':np.tile(np.arange(len(y)),6),'use_correction':False})
        be=float(np.max(np.abs(stacked.params-b)));ve=float(np.max(np.abs(stacked.cov_params()-v)))
        assert be<1e-7 and ve<1e-7
        check.append(dict(background=genus,cutoff=cutoff,beta_max_difference=be,cov_max_difference=ve))
        np.savez_compressed(O/f'{genus}_reads{cutoff}_joint.npz',beta=b,covariance=v,influence=u)
        print('pair contrasts',genus,cutoff,omnibus[-1],flush=True)

co=pd.DataFrame(contrasts);om=pd.DataFrame(omnibus)
for cutoff in [1,3]:
    ix=co.cutoff.eq(cutoff);co.loc[ix,'q_BH_30']=multipletests(co.loc[ix,'p'],method='fdr_bh')[1];co.loc[ix,'p_Holm_30']=multipletests(co.loc[ix,'p'],method='holm')[1]
    ix=om.cutoff.eq(cutoff);om.loc[ix,'q_BH_2']=multipletests(om.loc[ix,'p'],method='fdr_bh')[1]
co.to_csv(O/'表S22b_两背景全部属间比较.tsv',sep='\t',index=False);om.to_csv(O/'表S22a_两背景整体差异.tsv',sep='\t',index=False)
pd.DataFrame(coeff).to_csv(O/'表S22d_两背景各属关联系数.tsv',sep='\t',index=False)

sens=[];pcpairs=[];pcfits={}
for pc_n in [5,10,15]:
    with threadpool_limits(limits=1):pc=PCA(n_components=pc_n,svd_solver='full').fit(Z);sc=pc.transform(Z)
    sc=(sc-sc.mean(axis=0))/sc.std(axis=0,ddof=0);a=np.column_stack([X,sc]);k=a.shape[1]
    b,v,u=fit_joint(a,Y,f'PC{pc_n}_reads1');pcfits[pc_n]=(b,v,u,k)
    result=equal_vectors(b,v,list(range(6)),k,pc_n)
    sens.append(dict(condition=f'PC{pc_n}_all6',cutoff=1,PCs=pc_n,target_count=6,excluded='none',variance_explained=float(pc.explained_variance_ratio_.sum()),**result))
    if pc_n==10:
        np.savez_compressed(O/'pc_joint_sandwich.npz',beta=b,covariance=v,influence=u)
        old=np.load(P/'pc_joint_sandwich.npz');expected=json.loads((P/'pc_direct_cross_genus_heterogeneity.json').read_text(encoding='utf8'))
        assert abs(result['chi2']-expected['chi2'])<1e-6
        for i,j in itertools.combinations(range(6),2):
            pcpairs.append(dict(taxon_A=T[i],taxon_B=T[j],**equal_vectors(b,v,[i,j],k,10)))
        for drop in range(6):
            sens.append(dict(condition='leave_out_'+T[drop],cutoff=1,PCs=10,target_count=5,excluded=T[drop],variance_explained=float(pc.explained_variance_ratio_.sum()),**equal_vectors(b,v,[i for i in range(6) if i!=drop],k,10)))
        sens.append(dict(condition='two_most_prevalent',cutoff=1,PCs=10,target_count=2,excluded='four_less_prevalent',variance_explained=float(pc.explained_variance_ratio_.sum()),**equal_vectors(b,v,[0,1],k,10)))
        b3,v3,u3=fit_joint(a,C.ge(3).astype(float),'PC10_reads3')
        sens.append(dict(condition='reads3_all6',cutoff=3,PCs=10,target_count=6,excluded='none',variance_explained=float(pc.explained_variance_ratio_.sum()),**equal_vectors(b3,v3,list(range(6)),k,10)))
    print('PC sensitivity',pc_n,result,flush=True)
pcpairs=pd.DataFrame(pcpairs);pcpairs['q_BH_15']=multipletests(pcpairs.p,method='fdr_bh')[1];pcpairs['p_Holm_15']=multipletests(pcpairs.p,method='holm')[1]
pcpairs.to_csv(O/'表S22c_全部15组主成分向量比较.tsv',sep='\t',index=False)
pd.DataFrame(sens).to_csv(O/'表S23_整体差异稳健性.tsv',sep='\t',index=False)
pd.DataFrame(diag).to_csv(O/'模型数值诊断.tsv',sep='\t',index=False)
pd.DataFrame(check).to_csv(O/'堆叠模型交叉核验.tsv',sep='\t',index=False)
assert all(sha(Path(p))==h for p,h in hashes.items())
(O/'分析与核验摘要.json').write_text(json.dumps(dict(input_sha256=hashes,original_12_estimates_reproduced=True,original_PC_test_reproduced=True,stacked_checks=check,model_fits=len(diag),all_valid=True,significant_pairs=co.loc[(co.cutoff==1)&co.q_BH_30.lt(.05)].to_dict('records'),omnibus=om.to_dict('records'),PC_comparisons=pcpairs.to_dict('records'),sensitivity=sens),ensure_ascii=False,indent=2),encoding='utf8')
print('COMPLETE',len(diag),co.loc[(co.cutoff==1)&co.q_BH_30.lt(.05)].to_string(index=False),flush=True)
