"""Only the three fixed manuscript contrasts; never alter original files."""
from pathlib import Path
import sys,os,json,warnings,hashlib,argparse
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
# Use only this independent Python environment.

import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
P=Path(os.environ['AGP_AUDIT_OUTPUT']);P.mkdir(parents=True,exist_ok=True);S=Path(os.environ.get('AGP_AUDIT_SOURCE_PACKAGE',str(P.parent/'source_package')))
read=lambda p:pd.read_csv(p,sep='\t',dtype={'sample_id':str}).set_index('sample_id')
X=read(S/'02_统计输入/冻结统计输入/association_covariate_design.tsv')
Z=read(S/'02_统计输入/冻结统计输入/association_standardized_clr.tsv').loc[X.index]
C=read(S/'02_统计输入/冻结统计输入/association_target_counts.tsv').loc[X.index]
PC=pd.concat([read(S/'02_统计输入/饮食预处理'/f'{n}_PC_得分_内部复现.tsv').iloc[:,:10] for n in ['food','nutrient']],axis=1).loc[X.index]
XD=pd.concat([X,PC],axis=1);assert XD.shape==(2748,40)
M=read(P/'round_mapping_primary2748.tsv');ix=M.index[M.n_rounds.eq(1)];D=XD.loc[ix];Y=C.loc[ix,['Lacticaseibacillus','Lactobacillus']].ge(1).astype(float)
parser=argparse.ArgumentParser();parser.add_argument('--pooled',action='store_true');args=parser.parse_args()
schemes={'same_subset_baseline':D,'all_round_fixed_effects':pd.concat([D,pd.get_dummies(M.loc[ix,'rounds'],prefix='round',drop_first=True,dtype=float)],axis=1)}
if args.pooled:
    counts=M.loc[ix,'rounds'].value_counts();pooled=M.loc[ix,'rounds'].where(M.loc[ix,'rounds'].isin(counts[counts.ge(20)].index),'other_rounds_n_lt20')
    schemes['sample_count_pooled_round_fixed_effects']=pd.concat([D,pd.get_dummies(pooled,prefix='round',drop_first=True,dtype=float)],axis=1)
    table=pd.crosstab(pooled,Y.Lacticaseibacillus).add_prefix('Lacticaseibacillus_').join(pd.crosstab(pooled,Y.Lactobacillus).add_prefix('Lactobacillus_'))
    table.to_csv(P/'pooled_round_strata_diagnostics.tsv',sep='\t')
def independent_columns(d):
    cols=[];drops=[];rank=0
    for c in d:
        a=d[cols+[c]].to_numpy(float);r=np.linalg.matrix_rank(a)
        if r>rank:cols.append(c);rank=r
        else:drops.append({'column':c,'reason':'all_zero' if np.all(d[c].to_numpy()==0) else 'exact_linear_dependence_with_previous_columns'})
    return d[cols],drops
effects=[];contrasts=[];diagnostics=[];level_diag=[];columns={};coefficients=[]
for scheme,design0 in schemes.items():
    design,drops=independent_columns(design0);columns[scheme]={'input_columns':list(design0.columns),'retained_columns':list(design.columns),'removed_columns':drops}
    for c in design:
        vals=np.unique(design[c]);
        if c=='const' or len(vals)!=2:continue
        for t in Y:
            for v in vals:
                yy=Y.loc[design[c].eq(v),t];level_diag.append({'scheme':scheme,'term':c,'level':float(v),'target':t,'n':len(yy),'events':int(yy.sum()),'non_events':int(len(yy)-yy.sum()),'zero_event_or_nonevent_level':bool(yy.sum() in [0,len(yy)])})
    for genus in ['Anaerococcus','Finegoldia','Peptoniphilus_A']:
        bg='microbe::Peptoniphilaceae|'+genus;a=np.column_stack([design,Z.loc[ix,bg]]);names=list(design.columns)+['background'];fits={}
        for target in Y:
            yy=Y[target].to_numpy();audit={'scheme':scheme,'background':genus,'target':target,'n':len(yy),'events':int(yy.sum()),'parameters':a.shape[1],'rank':int(np.linalg.matrix_rank(a)),'condition_number':float(np.linalg.cond(a))}
            try:
                with warnings.catch_warnings(record=True) as ws,threadpool_limits(limits=1):
                    warnings.simplefilter('always');f=sm.GLM(yy,a,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0')
                mu=np.asarray(f.fittedvalues);h=a.T@(a*(mu*(1-mu))[:,None]);u=(a*(yy-mu)[:,None])@np.linalg.inv(h);v=np.asarray(f.cov_params());err=float(np.max(np.abs(u.T@u-v)))
                zero_rounds=[]
                if scheme!='same_subset_baseline':
                    lev=M.loc[ix,'rounds'] if scheme=='all_round_fixed_effects' else pooled
                    for r in sorted(lev.unique()):
                        yb=yy[lev.eq(r).to_numpy()]
                        if yb.sum() in [0,len(yb)]:zero_rounds.append(str(r))
                warning_text=' | '.join(str(w.message) for w in ws)
                valid=bool(f.converged and np.isfinite(f.params).all() and np.isfinite(v).all() and np.diag(v).min()>=0 and np.linalg.matrix_rank(a)==a.shape[1] and not zero_rounds and not any('separation' in str(w.message).lower() for w in ws) and err<1e-6)
                audit.update(converged=bool(f.converged),iterations=f.fit_history['iteration'],warnings=warning_text,max_abs_beta=float(np.max(np.abs(f.params))),min_fitted_probability=float(mu.min()),max_fitted_probability=float(mu.max()),zero_event_rounds=zero_rounds,negative_HC0_variances=int((np.diag(v)<0).sum()),minimum_HC0_variance=float(np.diag(v).min()),HC0_reconstruction_error=err,valid=valid)
                fits[target]=(np.asarray(f.params),u,v,valid)
                for k,name in enumerate(names):coefficients.append({'scheme':scheme,'background':genus,'target':target,'term':name,'beta':float(f.params[k]),'SE_HC0':float(np.sqrt(v[k,k])) if v[k,k]>=0 else np.nan,'model_valid':valid})
                b=float(f.params[-1]);se=float(np.sqrt(v[-1,-1]));effects.append({'scheme':scheme,'background':genus,'target':target,'beta':b,'SE_HC0':se,'OR':float(np.exp(b)),'low':float(np.exp(b-1.95996398454*se)),'high':float(np.exp(b+1.95996398454*se)),'valid':valid})
            except Exception as e:audit.update(valid=False,error=repr(e))
            diagnostics.append(audit)
        if len(fits)==2:
            a1,u1,v1,ok1=fits['Lacticaseibacillus'];a2,u2,v2,ok2=fits['Lactobacillus'];delta=a1[-1]-a2[-1];se=np.sqrt(np.sum((u1[:,-1]-u2[:,-1])**2));valid=ok1 and ok2
            contrasts.append({'scheme':scheme,'background':genus,'n':len(ix),'delta_logOR':float(delta),'SE':float(se),'ratio_OR':float(np.exp(delta)),'low':float(np.exp(delta-1.95996398454*se)),'high':float(np.exp(delta+1.95996398454*se)),'p':float(2*norm.sf(abs(delta/se))) if valid else np.nan,'valid':valid})
        else:contrasts.append({'scheme':scheme,'background':genus,'n':len(ix),'p':np.nan,'valid':False})
        print(scheme,genus,'valid',contrasts[-1]['valid'],flush=True)
co=pd.DataFrame(contrasts)
for scheme in schemes:
    sel=co.scheme.eq(scheme);ps=co.loc[sel,'p'];q=multipletests(ps.fillna(1),method='fdr_bh')[1];co.loc[sel,'q_BH_3']=np.where(ps.notna(),q,np.nan)
base=co.loc[co.scheme.eq('same_subset_baseline')].set_index('background')
if 'delta_logOR' in co:co['change_delta_logOR_from_same_subset']=co.apply(lambda r:r.get('delta_logOR',np.nan)-base.loc[r.background].get('delta_logOR',np.nan),axis=1)
co.loc[~co.valid].assign(diagnostic_only=True).to_csv(P/'raw_failed_estimates.tsv',sep='\t',index=False)
for c in ['delta_logOR','SE','ratio_OR','low','high','p','q_BH_3','change_delta_logOR_from_same_subset']:
    if c in co:co.loc[~co.valid,c]=np.nan
co.to_csv(P/'round_contrasts.tsv',sep='\t',index=False)
ef=pd.DataFrame(effects);ef.loc[~ef.valid].assign(diagnostic_only=True).to_csv(P/'raw_failed_target_effects.tsv',sep='\t',index=False)
for c in ['beta','SE_HC0','OR','low','high']:ef.loc[~ef.valid,c]=np.nan
ef.to_csv(P/'round_target_effects.tsv',sep='\t',index=False)
pd.DataFrame(coefficients).to_csv(P/'round_all_coefficients.tsv',sep='\t',index=False)
pd.DataFrame(level_diag).to_csv(P/'round_binary_level_diagnostics.tsv',sep='\t',index=False)
(P/'round_model_diagnostics.json').write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2),encoding='utf8')
(P/'round_design_columns.json').write_text(json.dumps(columns,ensure_ascii=False,indent=2),encoding='utf8')
(P/'round_runtime.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'statsmodels':sm.__version__,'included_n':len(ix),'schemes':list(schemes),'note':'Separate per-scheme 3-test families. Invalid unpenalized batch models are not inferential evidence.'},ensure_ascii=False,indent=2),encoding='utf8')
print(co.to_string(index=False),flush=True)
