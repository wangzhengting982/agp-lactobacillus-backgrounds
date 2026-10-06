"""Read-only targeted sensitivity audit from the archived local source package.

Usage: python run_diversity_audit.py --source-root PATH --output-dir PATH
Dependencies: numpy, pandas, scipy, statsmodels. Existing local runtimes may be
passed as repeated --runtime PATH flags. No original files are modified.
"""
from pathlib import Path
import argparse, os, sys, json, hashlib, warnings, platform
from datetime import datetime, timezone
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[key] = '1'
parser = argparse.ArgumentParser()
parser.add_argument('--source-root', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
parser.add_argument('--runtime', type=Path, action='append', default=[])
args = parser.parse_args()
for path in reversed(args.runtime):
    sys.path.insert(0, str(path.resolve()))
import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm
from scipy.linalg import block_diag
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
sys.stdout.reconfigure(encoding='utf-8')
ROOT = args.source_root.resolve()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
BG = ['Anaerococcus', 'Finegoldia', 'Peptoniphilus_A']
TARGETS = ['Lacticaseibacillus', 'Lactobacillus']
DROP = ['nonlab_shannon', 'nonlab_logrichness']
paths = {
    'covariates': ROOT/'02_统计输入/冻结统计输入/association_covariate_design.tsv',
    'background': ROOT/'02_统计输入/冻结统计输入/association_standardized_clr.tsv',
    'targets': ROOT/'02_统计输入/冻结统计输入/association_target_counts.tsv',
    'food_pc': ROOT/'02_统计输入/饮食预处理/food_PC_得分_内部复现.tsv',
    'nutrient_pc': ROOT/'02_统计输入/饮食预处理/nutrient_PC_得分_内部复现.tsv',
    'previous_contrasts': ROOT/'03_饮食与菌群深化/完整分析结果/P02_广泛饮食调整全部差异.tsv',
    'previous_models': ROOT/'03_饮食与菌群深化/完整分析结果/P01_广泛饮食调整全部关联.tsv',
}
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):
    d = pd.read_csv(path, sep='\t', dtype={'sample_id': str}).set_index('sample_id')
    assert d.index.is_unique, str(path)
    return d
def write(rows, name):
    pd.DataFrame(rows).to_csv(OUT/name, sep='\t', index=False, encoding='utf-8-sig')
manifest = [{'key': key, 'relative_path': str(path.relative_to(ROOT)),
             'bytes': path.stat().st_size, 'sha256': sha(path)} for key, path in paths.items()]
write(manifest, 'input_hashes.tsv')
X = read(paths['covariates'])
Z = read(paths['background']).loc[X.index]
Y = read(paths['targets']).loc[X.index]
PC = pd.concat([read(paths['food_pc']).iloc[:, :10],
                read(paths['nutrient_pc']).iloc[:, :10]], axis=1).loc[X.index]
XD = pd.concat([X, PC], axis=1)
assert X.shape == (2748, 20) and XD.shape == (2748, 40)
assert XD.columns.is_unique and all(c in XD for c in DROP)
assert all(d.index.equals(X.index) for d in [Z, Y, PC])
assert all(np.isfinite(d.to_numpy(float)).all() for d in [XD, Z, Y])
assert list(PC.columns) == ([f'food_PC{i}' for i in range(1,11)] +
                            [f'nutrient_PC{i}' for i in range(1,11)])
previous_c = pd.read_csv(paths['previous_contrasts'], sep='\t')
previous_m = pd.read_csv(paths['previous_models'], sep='\t')
models, contrasts, diagnostics, independent, reproduction = [], [], [], [], []
started = datetime.now(timezone.utc).isoformat()
environment = {
    'started_utc': started, 'python': sys.version, 'platform': platform.platform(),
    'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__,
    'statsmodels': statsmodels.__version__, 'runtime_paths': [str(p) for p in args.runtime],
    'n': len(X), 'base_columns': list(XD.columns), 'dropped_columns': DROP,
    'sample_order_sha256': hashlib.sha256('\n'.join(X.index).encode()).hexdigest(),
    'script_sha256': sha(Path(__file__)),
    'plan_sha256': sha(OUT/'fixed_plan.md') if (OUT/'fixed_plan.md').exists() else None,
}
(OUT/'run_metadata.json').write_text(json.dumps(environment, ensure_ascii=False, indent=2), encoding='utf-8')

def fit_one(design, response, meta):
    a, yy = design.to_numpy(float), np.asarray(response, float)
    assert np.linalg.matrix_rank(a) == a.shape[1]
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        f = sm.GLM(yy, a, family=sm.families.Binomial()).fit(
            maxiter=200, tol=1e-10, cov_type='HC0')
    mu = np.asarray(f.fittedvalues)
    h = a.T @ (a * (mu*(1-mu))[:,None])
    eig = np.linalg.eigvalsh(h)
    u = (a * (yy-mu)[:,None]) @ np.linalg.inv(h)
    se_error = float(np.max(np.abs(np.sqrt(np.sum(u*u,axis=0))-f.bse)))
    score_max = float(np.max(np.abs(a.T@(yy-mu))))
    valid = bool(f.converged and not ws and np.isfinite(f.params).all()
                 and np.isfinite(f.cov_params()).all() and eig.min()>0
                 and se_error < 1e-7 and score_max < 1e-5)
    diagnostics.append(dict(**meta, n=len(yy), events=int(yy.sum()),
        parameters=a.shape[1], rank=int(np.linalg.matrix_rank(a)), valid=valid,
        converged=bool(f.converged), hessian_min_eigenvalue=float(eig.min()),
        hessian_condition=float(np.linalg.cond(h)), max_abs_score=score_max,
        max_HC0_SE_error=se_error, warning=' | '.join(str(w.message) for w in ws),
        min_fitted_probability=float(mu.min()), max_fitted_probability=float(mu.max())))
    beta, se = float(f.params[-1]), float(f.bse[-1])
    models.append(dict(**meta, n=len(yy), events=int(yy.sum()), parameters=a.shape[1],
        valid=valid, beta=beta, se_HC0=se, OR=np.exp(beta), low=np.exp(beta-1.96*se),
        high=np.exp(beta+1.96*se), p=2*norm.sf(abs(beta/se))))
    assert valid, diagnostics[-1]
    return f, u

def run_pair(threshold, spec, background):
    xx = XD if spec == 'with_diversity' else XD.drop(columns=DROP)
    d = xx.assign(background=Z['microbe::Peptoniphilaceae|'+background])
    a = d.to_numpy(float)
    fs=[]
    for target in TARGETS:
        fs.append(fit_one(d, Y[target].ge(threshold),
                         dict(threshold=threshold, spec=spec, background=background, target=target)))
    delta=float(fs[0][0].params[-1]-fs[1][0].params[-1])
    uu=fs[0][1][:,-1]-fs[1][1][:,-1]
    se=float(np.sqrt(np.sum(uu*uu)))
    row=dict(threshold=threshold, spec=spec, background=background, n=len(d),
        target_A=TARGETS[0], target_B=TARGETS[1], valid=True, delta_log_OR=delta,
        se_joint=se, ratio_OR=np.exp(delta), low=np.exp(delta-1.96*se),
        high=np.exp(delta+1.96*se), p=2*norm.sf(abs(delta/se)))
    contrasts.append(row)
    # A distinct optimizer and statsmodels' grouped-score covariance implementation.
    stacked=block_diag(a, a)
    ys=np.r_[Y[TARGETS[0]].ge(threshold).to_numpy(float),
             Y[TARGETS[1]].ge(threshold).to_numpy(float)]
    groups=np.tile(np.arange(len(d)),2)
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        sf=sm.GLM(ys, stacked, family=sm.families.Binomial()).fit(
            method='newton', start_params=np.zeros(stacked.shape[1]), maxiter=200,
            tol=1e-10, cov_type='cluster',
            cov_kwds={'groups': groups, 'use_correction': False})
    k=a.shape[1]; c=np.zeros(2*k);c[k-1]=1;c[2*k-1]=-1
    sd=float(c@sf.params); ss=float(np.sqrt(c@sf.cov_params()@c))
    pe=float(np.max(np.abs(sf.params-np.r_[fs[0][0].params,fs[1][0].params])))
    converged=bool(sf.mle_retvals.get('converged',False))
    ok=bool(converged and not ws and abs(sd-delta)<1e-7 and abs(ss-se)<1e-7 and pe<1e-7)
    independent.append(dict(threshold=threshold,spec=spec,background=background,
        stacked_converged=converged,stacked_rank=int(np.linalg.matrix_rank(stacked)),
        stacked_parameters=2*k,stacked_delta=sd,stacked_se=ss,
        delta_absolute_error=abs(sd-delta),se_absolute_error=abs(ss-se),
        max_all_coefficients_error=pe,passed=ok,
        warning=' | '.join(str(w.message) for w in ws)))
    assert ok, independent[-1]
    print(f'Completed threshold={threshold}, {spec}, {background}: '
          f'OR ratio={row["ratio_OR"]:.9f}, p={row["p"]:.3g}', flush=True)

def checkpoint():
    write(models,'model_results.tsv');write(contrasts,'contrast_results.tsv')
    write(diagnostics,'model_diagnostics.tsv');write(independent,'independent_checks.tsv')
    write(reproduction,'original_reproduction.tsv')

try:
    # Gate: all three original comparisons must match before running any new specification.
    for bg in BG:
        run_pair(1,'with_diversity',bg)
        old=previous_c[(previous_c.spec=='food10_nutrient10')&(previous_c.background==bg)].iloc[0]
        new=contrasts[-1]
        errors={field:abs(float(new[field])-float(old[field]))
                for field in ['delta_log_OR','se_joint','ratio_OR','low','high']}
        model_errors=[]
        for m in models[-2:]:
            oldm=previous_m[(previous_m.spec=='food10_nutrient10')&
                           (previous_m.background==bg)&(previous_m.target==m['target'])].iloc[0]
            model_errors.append(abs(m['OR']-float(oldm.OR)))
        ok=max([*errors.values(),*model_errors])<1e-8
        reproduction.append(dict(background=bg,passed=ok,
            **{k+'_abs_error':v for k,v in errors.items()},max_model_OR_error=max(model_errors)))
        checkpoint()
        assert ok, reproduction[-1]
    for threshold,spec in [(1,'without_diversity'),(3,'with_diversity'),(3,'without_diversity')]:
        for bg in BG:
            run_pair(threshold,spec,bg);checkpoint()
    cdf=pd.DataFrame(contrasts)
    for _, idx in cdf.groupby(['threshold','spec']).groups.items():
        cdf.loc[idx,'q_BH_selected3']=multipletests(cdf.loc[idx,'p'],method='fdr_bh')[1]
    cdf.to_csv(OUT/'contrast_results.tsv',sep='\t',index=False,encoding='utf-8-sig')
    side=cdf.pivot(index=['threshold','background'],columns='spec',values=['ratio_OR','low','high','p','q_BH_selected3','delta_log_OR'])
    side.columns=['_'.join(c) for c in side.columns]
    side['ratio_OR_relative_change_percent']=(side['ratio_OR_without_diversity']/side['ratio_OR_with_diversity']-1)*100
    side.reset_index().to_csv(OUT/'specification_comparison.tsv',sep='\t',index=False,encoding='utf-8-sig')
    status={'status':'success','finished_utc':datetime.now(timezone.utc).isoformat(),
            'models':len(models),'contrasts':len(contrasts),'all_model_valid':all(d['valid'] for d in diagnostics),
            'all_reproduction_passed':all(r['passed'] for r in reproduction),
            'all_independent_checks_passed':all(r['passed'] for r in independent)}
    (OUT/'completion_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(status,ensure_ascii=False),flush=True)
    print(cdf[['threshold','spec','background','ratio_OR','low','high','p','q_BH_selected3']].to_string(index=False),flush=True)
except Exception as e:
    checkpoint()
    (OUT/'completion_status.json').write_text(json.dumps({'status':'failed','error':repr(e),'finished_utc':datetime.now(timezone.utc).isoformat()},ensure_ascii=False,indent=2),encoding='utf-8')
    raise
