"""Refit taxonomy PC models from rebuilt input matrices, without using old fits.

Statistical settings match archived taxonomy_association_sensitivity.py, lines
115-170. Old coefficients are read only after the current fits for comparison.
Upstream SINTAX reclassification is outside this entry point's scope.
"""
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '1'
import json
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.decomposition import PCA
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
from core_paths import DATA, OUT, PRIMARY, TAXONOMY_INPUT

TAXA = ['Lacticaseibacillus', 'Lactobacillus', 'Leuconostoc',
        'Limosilactobacillus', 'Ligilactobacillus', 'Pediococcus']
pc_results, omnibus_results, audits = [], [], []
for release in ['R07', 'R06']:
    source = TAXONOMY_INPUT / release
    reference = DATA / '输入/03_分析与结果/taxonomy_validation/association_sensitivity' / release
    dest = OUT / 'taxonomy_refits' / release
    dest.mkdir(parents=True, exist_ok=True)
    def read(name):
        return pd.read_csv(source / name, sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    x = read('covariate_design.tsv')
    z = read('standardized_clr.tsv.gz').loc[x.index]
    y = read('target_counts.tsv').loc[x.index, TAXA].ge(1).astype(int)
    assert x.shape == (2748, 20) and z.shape[0] == 2748
    pca = PCA(n_components=10, svd_solver='full')
    with threadpool_limits(limits=1):
        scores = pca.fit_transform(z)
    scores = (scores - scores.mean(axis=0)) / scores.std(axis=0, ddof=0)
    a = np.column_stack([x, scores]); k = a.shape[1]
    assert np.linalg.matrix_rank(a) == k == 30
    betas, influences, local = [], [], []
    for taxon in TAXA:
        yy = y[taxon].to_numpy(dtype=float)
        with threadpool_limits(limits=1), warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always')
            fit = sm.GLM(yy, a, family=sm.families.Binomial()).fit(maxiter=150, tol=1e-10, cov_type='HC0')
        covariance = np.asarray(fit.cov_params())
        valid = bool(fit.converged and np.isfinite(fit.params).all() and np.isfinite(covariance).all()
                     and not any('separation' in str(w.message).lower() for w in ws))
        assert valid
        beta = np.asarray(fit.params)[-10:]; pc_cov = covariance[-10:, -10:]
        rank = int(np.linalg.matrix_rank(pc_cov))
        stat = float(beta @ np.linalg.pinv(pc_cov) @ beta)
        local.append(dict(release=release,taxon=taxon,n=len(y),events=int(yy.sum()),chi2=stat,df=rank,p=float(chi2.sf(stat,rank)),valid=valid))
        fitted = np.asarray(fit.fittedvalues)
        influence = (a * (yy-fitted)[:,None]) @ np.linalg.inv(a.T @ (a * (fitted*(1-fitted))[:,None]))
        reconstruction_error = float(np.max(np.abs(influence.T @ influence - covariance)))
        assert np.allclose(influence.T @ influence, covariance, rtol=2e-5, atol=1e-7)
        audits.append(dict(release=release,taxon=taxon,converged=bool(fit.converged),warnings=' | '.join(str(w.message) for w in ws),cov_reconstruction_max_absolute_error=reconstruction_error))
        influences.append(influence); betas.append(np.asarray(fit.params))
    local = pd.DataFrame(local); local['q_BH_six'] = multipletests(local.p,method='fdr_bh')[1]
    pc_results.extend(local.to_dict('records'))
    joint_beta=np.concatenate(betas); u=np.column_stack(influences); joint_cov=u.T@u
    contrast=np.zeros((50,180))
    for t in range(1,6):
        for j in range(10):
            contrast[(t-1)*10+j,20+j]=-1
            contrast[(t-1)*10+j,t*k+20+j]=1
    difference=contrast@joint_beta; V=contrast@joint_cov@contrast.T
    rank=int(np.linalg.matrix_rank(V)); stat=float(difference@np.linalg.pinv(V)@difference)
    assert rank == 50
    omnibus_results.append(dict(release=release,comparison='all six target genera have equal 10-PC log-odds coefficient vectors',chi2=stat,df=rank,nominal_constraints=50,p=float(chi2.sf(stat,rank)),n_subjects=len(y),same_subject_correlation='joint sandwich using six GLMs individual score contributions',valid=True,scope='explanatory relative-composition association heterogeneity, not prediction improvement or causal heterogeneity'))
    np.savez_compressed(dest/'pc_joint_sandwich.npz',beta=joint_beta,covariance=joint_cov)
    old=np.load(reference/'pc_joint_sandwich.npz')
    errors=dict(beta_max_absolute_error=float(np.max(np.abs(joint_beta-old['beta']))),covariance_max_absolute_error=float(np.max(np.abs(joint_cov-old['covariance']))))
    (dest/'comparison_to_archive.json').write_text(json.dumps(errors,indent=2),encoding='utf8')
    assert np.allclose(joint_beta,old['beta'],rtol=1e-7,atol=1e-8)
    assert np.allclose(joint_cov,old['covariance'],rtol=1e-7,atol=1e-8)
    print(release,'PC refit complete',errors,flush=True)
pd.DataFrame(pc_results).to_csv(OUT/'refitted_taxonomy_pc_blocks.tsv',sep='\t',index=False)
pd.DataFrame(omnibus_results).to_csv(OUT/'refitted_taxonomy_pc_omnibus.tsv',sep='\t',index=False)
pd.DataFrame(audits).to_csv(OUT/'taxonomy_pc_model_diagnostics.tsv',sep='\t',index=False)
