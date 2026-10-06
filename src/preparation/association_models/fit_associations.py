from pathlib import Path
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import json, warnings, hashlib, sys
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2
from scipy.linalg import block_diag
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed, parallel_config
from threadpoolctl import threadpool_limits

PROJECT=Path(__file__).resolve().parents[2]
HERE=PROJECT/'03_分析与结果/association_models'
SRC=PROJECT/'03_分析与结果/agp_models'
RAW=PROJECT/'05_输入数据/01_原始输入'
TAXA=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
NUM=['age','bmi','antibiotic','ibd','ibs','year','depth_nonlab','energy','fiber_density','yogurt','plant_ferment','probiotic','nonlab_shannon','nonlab_logrichness']

def prepare():
    cov=pd.read_csv(SRC/'raw_covariates.tsv',sep='\t',dtype={'sample_id':str}).set_index('sample_id')
    bg=pd.read_csv(SRC/'background_genus_counts.tsv',sep='\t',dtype={'sample_id':str}).set_index('sample_id').loc[cov.index]
    targets=pd.read_csv(SRC/'target_genus_counts.tsv',sep='\t',dtype={'sample_id':str}).set_index('sample_id').loc[cov.index,TAXA]
    all_ids=cov.index.copy()
    depth_eligible=cov.total_reads_qc.ge(5000)
    ix=cov.index[cov.sex.isin(['female','male'])&depth_eligible]
    cov=cov.loc[ix].copy(); bg=bg.loc[ix]; targets=targets.loc[ix]
    keep=bg.gt(0).mean().ge(.1)
    meta=pd.read_csv(RAW/'sample_information_10317_matched_4168.tsv',sep='\t',dtype=str,usecols=['sample_name','collection_timestamp']).set_index('sample_name').loc[ix]
    cov['year']=pd.to_datetime(meta.collection_timestamp,format='mixed',errors='coerce').dt.year
    assert cov.year.notna().all()
    medians=cov[NUM].median()
    num=cov[NUM].fillna(medians)
    X=(num-num.mean())/num.std(ddof=0)
    X['bmi_missing']=cov.bmi.isna().astype(float)
    X['sex_male']=cov.sex.eq('male').astype(float)
    for c in ['US','UK','Unknown']:
        X['country_'+c]=cov.country.eq(c).astype(float)
    X=sm.add_constant(X)
    assert np.linalg.matrix_rank(X.to_numpy())==X.shape[1]
    rel=bg.loc[:,keep].div(bg.sum(axis=1),axis=0)
    lg=np.log(rel+.00005); clr=lg.sub(lg.mean(axis=1),axis=0)
    clr_std=clr.sub(clr.mean()).div(clr.std(ddof=0))
    pca=PCA(n_components=10,svd_solver='full')
    scores=pca.fit_transform(clr_std)
    scores=(scores-scores.mean(axis=0))/scores.std(axis=0,ddof=0)
    pcs=pd.DataFrame(scores,index=ix,columns=[f'PC{i+1:02d}' for i in range(10)])
    var=pd.DataFrame({'PC':pcs.columns,'explained_variance_ratio':pca.explained_variance_ratio_,'cumulative_explained_variance_ratio':np.cumsum(pca.explained_variance_ratio_)})
    var.to_csv(HERE/'pc_explained_variance.tsv',sep='\t',index=False)
    pd.DataFrame(pca.components_.T,index=clr_std.columns,columns=pcs.columns).to_csv(HERE/'pc_loadings.tsv',sep='\t')
    X.to_csv(HERE/'association_covariate_design.tsv',sep='\t')
    pcs.to_csv(HERE/'association_pc_scores.tsv',sep='\t')
    clr_std.to_csv(HERE/'association_standardized_clr.tsv',sep='\t')
    targets.to_csv(HERE/'association_target_counts.tsv',sep='\t')
    summary={'n':len(ix),'original_eligible_n':len(all_ids),'depth_5000_eligible_n':int(depth_eligible.sum()),'excluded_unknown_or_other_sex_after_depth_QC':int(depth_eligible.sum())-len(ix),'events':targets.ge(1).sum().to_dict(),'background_genera':int(keep.sum()),'PCs':10,'PC_variance_explained':float(pca.explained_variance_ratio_.sum()),'baseline_parameters':X.shape[1],'pc_model_parameters':X.shape[1]+10,'BMI_missing':int(cov.bmi.isna().sum()),'BMI_imputation_median':float(medians.bmi),'year_source':'raw collection_timestamp parsed with format=mixed','country_reference':'Other','sex_reference':'female','CLR_pseudocount_relative_abundance':.00005,'primary_status':'exploratory explanatory association following known pilot results'}
    return X,pcs,clr_std,targets.ge(1).astype(int),summary

def fit_glm(X,y,meta):
    a=np.asarray(X,float); yy=np.asarray(y,float)
    with threadpool_limits(limits=1),warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        fit=sm.GLM(yy,a,family=sm.families.Binomial()).fit(maxiter=150,tol=1e-10,cov_type='HC0')
    warning_text=' | '.join(str(w.message) for w in ws)
    audit={**meta,'n':len(yy),'events':int(yy.sum()),'parameters':a.shape[1],'events_per_parameter':float(yy.sum()/a.shape[1]),'design_rank':int(np.linalg.matrix_rank(a)),'design_condition_number':float(np.linalg.cond(a)),'converged':bool(fit.converged),'warnings':warning_text,'max_abs_beta':float(np.max(np.abs(fit.params))),'min_probability':float(fit.fittedvalues.min()),'max_probability':float(fit.fittedvalues.max()),'covariance_rank':int(np.linalg.matrix_rank(fit.cov_params()))}
    valid=fit.converged and audit['design_rank']==a.shape[1] and np.isfinite(fit.params).all() and np.isfinite(fit.cov_params()).all() and not any('separation' in str(w.message).lower() for w in ws)
    audit['valid']=bool(valid)
    return fit,audit

def single_fit(X,z,y,taxon,feature):
    design=np.column_stack([X,z])
    fit,audit=fit_glm(design,y,{'analysis':'single_genus','taxon':taxon,'feature':feature})
    beta=fit.params[-1];se=np.sqrt(fit.cov_params()[-1,-1]);p=float(fit.pvalues[-1]) if audit['valid'] else 1.
    result={'taxon':taxon,'feature':feature,'n':len(y),'events':int(np.sum(y)),'beta_log_OR_per_SD_CLR':float(beta),'SE_HC0':float(se),'OR':float(np.exp(beta)),'CI_low':float(np.exp(beta-1.95996398454*se)),'CI_high':float(np.exp(beta+1.95996398454*se)),'p':p,'valid':audit['valid']}
    return result,audit

def main():
    if '--from-prepared' in sys.argv:
        read=lambda name:pd.read_csv(HERE/name,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
        X=read('association_covariate_design.tsv')
        pcs=read('association_pc_scores.tsv').loc[X.index]
        clr=read('association_standardized_clr.tsv').loc[X.index]
        y=read('association_target_counts.tsv').loc[X.index,TAXA].ge(1).astype(int)
        summary=json.loads((HERE/'association_summary.json').read_text(encoding='utf-8'))
    else:
        X,pcs,clr,y,summary=prepare()
    Xpc=pd.concat([X,pcs],axis=1); a=Xpc.to_numpy(); k=a.shape[1]
    pc_results=[];coeff=[];audit=[];influences=[];betas=[]
    for taxon in TAXA:
        fit,qa=fit_glm(Xpc,y[taxon],{'analysis':'PC_block','taxon':taxon,'feature':'PC01-PC10'})
        audit.append(qa); beta=fit.params[-10:];cov=np.asarray(fit.cov_params())[-10:,-10:]
        rank=np.linalg.matrix_rank(cov); stat=float(beta@np.linalg.pinv(cov)@beta)
        pc_results.append({'taxon':taxon,'n':len(y),'events':int(y[taxon].sum()),'chi2':stat,'df':int(rank),'p':float(chi2.sf(stat,rank)) if qa['valid'] else 1.,'valid':qa['valid']})
        for j,c in enumerate(Xpc.columns):
            b=fit.params[j];se=np.sqrt(fit.cov_params()[j,j])
            coeff.append({'taxon':taxon,'term':c,'beta':b,'SE_HC0':se,'OR':float(np.exp(b)),'CI_low':float(np.exp(b-1.95996398454*se)),'CI_high':float(np.exp(b+1.95996398454*se)),'p':float(fit.pvalues[j])})
        prob=np.asarray(fit.fittedvalues);bread=a.T@(a*(prob*(1-prob))[:,None]);bread_inv=np.linalg.inv(bread)
        influence=(a*(y[taxon].to_numpy()-prob)[:,None])@bread_inv
        assert np.allclose(influence.T@influence,np.asarray(fit.cov_params()),rtol=2e-5,atol=1e-7)
        influences.append(influence);betas.append(np.asarray(fit.params))
    pc_table=pd.DataFrame(pc_results);pc_table['q_BH_six']=multipletests(pc_table.p,method='fdr_bh')[1]
    pc_table.to_csv(HERE/'pc_block_associations.tsv',sep='\t',index=False)
    pd.DataFrame(coeff).to_csv(HERE/'pc_model_all_coefficients.tsv',sep='\t',index=False)
    infl=np.column_stack(influences);joint_cov=infl.T@infl; joint_beta=np.concatenate(betas)
    R=np.zeros((50,6*k))
    for t in range(1,6):
        for j in range(10):
            R[(t-1)*10+j,k-10+j]=-1
            R[(t-1)*10+j,t*k+k-10+j]=1
    diff=R@joint_beta;V=R@joint_cov@R.T;rank=np.linalg.matrix_rank(V)
    stat=float(diff@np.linalg.pinv(V)@diff)
    heterogeneity={'comparison':'all six target genera have equal 10-PC log-odds coefficient vectors','chi2':stat,'df':int(rank),'nominal_constraints':50,'p':float(chi2.sf(stat,rank)),'n_subjects':len(y),'same_subject_correlation':'joint sandwich using six GLMs individual score contributions','valid':all(q['valid'] for q in audit),'scope':'explanatory relative-composition association heterogeneity, not prediction improvement or causal heterogeneity'}
    (HERE/'pc_direct_cross_genus_heterogeneity.json').write_text(json.dumps(heterogeneity,ensure_ascii=False,indent=2),encoding='utf-8')
    np.savez_compressed(HERE/'pc_joint_sandwich.npz',covariance=joint_cov,beta=joint_beta,contrast_matrix=R,contrast_covariance=V)
    print('PC BLOCK',pc_table.to_string(index=False),json.dumps(heterogeneity),flush=True)
    tasks=[(X.to_numpy(),clr[feature].to_numpy(),y[taxon].to_numpy(),taxon,feature) for taxon in TAXA for feature in clr.columns]
    results=[]
    with parallel_config(backend='loky',n_jobs=2,inner_max_num_threads=1):
        for i,(res,qa) in enumerate(Parallel(return_as='generator_unordered')(delayed(single_fit)(*t) for t in tasks),1):
            results.append(res);audit.append(qa)
            if i%100==0: print('SINGLE GENUS',i,'/',len(tasks),flush=True)
    allres=pd.DataFrame(results).sort_values(['taxon','feature']);allres['q_BH_all']=multipletests(allres.p,method='fdr_bh')[1]
    allres.to_csv(HERE/'single_genus_all_associations.tsv',sep='\t',index=False)
    allres.pivot(index='feature',columns='taxon',values='beta_log_OR_per_SD_CLR').to_csv(HERE/'heatmap_complete_beta.tsv',sep='\t')
    allres.pivot(index='feature',columns='taxon',values='q_BH_all').to_csv(HERE/'heatmap_complete_q.tsv',sep='\t')
    pd.DataFrame(audit).to_csv(HERE/'model_quality_audit.tsv',sep='\t',index=False)
    summary['GLM_fits']=len(audit);summary['invalid_models']=sum(not q['valid'] for q in audit);summary['warning_models']=sum(bool(q['warnings']) for q in audit)
    summary['single_genus_total_tests']=len(allres)
    summary['single_genus_q_lt_05_counts']=allres.loc[allres.q_BH_all.lt(.05)].groupby('taxon').size().to_dict()
    (HERE/'association_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
    print(allres.loc[allres.q_BH_all.lt(.05)].sort_values('q_BH_all').head(30).to_string(index=False),flush=True)

if __name__=='__main__':main()
