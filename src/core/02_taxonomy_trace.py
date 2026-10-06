from pathlib import Path
from core_paths import DATA, OUT, PRIMARY, TAXONOMY_INPUT
import sys,os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pandas as pd,numpy as np,h5py,json,hashlib,re
from scipy.sparse import csr_matrix
R=DATA;O=OUT;O.mkdir(exist_ok=True)
P=(DATA/'输入/03_分析与结果/taxonomy_validation');RAW=(DATA/'输入/05_输入数据/01_原始输入');A=(DATA/'输入/03_分析与结果/association_models')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
ix=pd.read_csv(PRIMARY/'association_covariate_design.tsv',sep='\t',dtype={'sample_id':str}).sample_id.tolist()
with h5py.File(RAW/'4168_ASV_table.biom') as h:
    samples=[x.decode() if isinstance(x,bytes) else str(x) for x in h['sample/ids'][:]]
    seq=[x.decode() if isinstance(x,bytes) else str(x) for x in h['observation/ids'][:]]
    m=h['observation/matrix'];mat=csr_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=(len(seq),len(samples))).astype(np.int64)
    keys=list(h['observation/metadata'].keys()) if 'observation/metadata' in h else []
pos=[samples.index(s) for s in ix];mi=mat[:,pos]
origin=pd.read_csv(RAW/'4168_annotated_feature_table.tsv',sep='\t',index_col=0).loc[:,ix]
seqidx=pd.read_csv(P/'asv_sequence_index.tsv',sep='\t');assert seqidx.asv_sequence.tolist()==seq
rows=[];samplecols={};asvout=[];refrows=[]
for genus in ['Oscillibacter','Veillonella']:
    inds=origin.index.str.contains('g__'+genus+';',regex=False);orig=origin.loc[inds].sum(axis=0).astype(int)
    samplecols['original_'+genus]=orig
    for release in ['R07','R06']:
        a=pd.read_csv(P/(release+'_all_asvs_genus_assignments.tsv'),sep='\t')
        inv=pd.read_csv(P/'references'/(release+'_genus_inventory.tsv'),sep='\t')
        for _,r in inv.iterrows():
            if genus in str(r.to_dict()):refrows.append(dict(release=release,**r.to_dict()))
        for support in [.5,.7,.8,.9]:
            mask=a.genus.eq(genus)&a.genus_bootstrap.ge(support)&a.domain.eq('Bacteria')&a.domain_bootstrap.ge(.8)
            c=np.asarray(mi[mask.to_numpy()].sum(axis=0)).ravel();present=c>0
            rows.append(dict(release=release,genus=genus,support=support,assigned_ASVs=int(mask.sum()),ASVs_observed_2748=int(np.asarray(mi[mask.to_numpy()].sum(axis=1)).ravel().astype(bool).sum()),reads=int(c.sum()),positive=int(present.sum()),prevalence=float(present.mean()),passes_10pct_filter=bool(present.mean()>=.1),original_positive=int(orig.gt(0).sum()),original_reads=int(orig.sum()),same_positive=int((present&orig.gt(0).to_numpy()).sum()),read_ratio_to_original=float(c.sum()/orig.sum()),scope='aggregate comparison; original per-ASV mapping unavailable'))
            samplecols[f'{release}_{genus}_{support}']=pd.Series(c,index=ix)
        subset=a.loc[a.genus.eq(genus)].copy();subset['release']=release
        subset['sequence']=seqidx.loc[subset.biom_row,'asv_sequence'].to_numpy()
        subset['reads_2748']=np.asarray(mi[subset.biom_row.to_numpy()].sum(axis=1)).ravel()
        asvout.append(subset)
df=pd.DataFrame(rows);df.to_csv(O/'表S24a_重点背景属分类支持与覆盖.tsv',sep='\t',index=False)
pd.concat(asvout).to_csv(O/'表S24b_重新注释候选ASV明细.tsv',sep='\t',index=False)
pd.DataFrame(samplecols).rename_axis('sample_id').to_csv(O/'分类计数逐样本.tsv',sep='\t')
pd.DataFrame(refrows).to_csv(O/'表S24c_参考库属标签.tsv',sep='\t',index=False)
# Direct high-confidence classification sensitivity for Veillonella uses existing independently rebuilt CLR/design.
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import chi2,norm
import itertools
taxa=['Lacticaseibacillus','Lactobacillus','Leuconostoc','Limosilactobacillus','Ligilactobacillus','Pediococcus']
vr=[];vc=[];vo=[]
for release in ['R07','R06']:
    sub=TAXONOMY_INPUT/release
    def read(n):return pd.read_csv(sub/n,sep='\t',dtype={'sample_id':str}).set_index('sample_id').loc[ix]
    x=read('covariate_design.tsv');z=read('standardized_clr.tsv.gz');y=read('target_counts.tsv')[taxa].ge(1).astype(float)
    feat=[c for c in z if c.endswith('|Veillonella')];assert len(feat)==1
    design=np.column_stack([x,z[feat[0]]]);k=design.shape[1];bs=[];us=[]
    for t in taxa:
        f=sm.GLM(y[t],design,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-10,cov_type='HC0')
        assert f.converged and np.isfinite(f.cov_params()).all().all()
        mu=np.asarray(f.fittedvalues);u=(design*(y[t].to_numpy()-mu)[:,None])@np.linalg.inv(design.T@(design*(mu*(1-mu))[:,None]))
        bs.append(np.asarray(f.params));us.append(u);b=float(f.params.iloc[-1]);se=float(f.bse.iloc[-1])
        vr.append(dict(release=release,taxon=t,beta=b,OR=float(np.exp(b)),low=float(np.exp(b-1.96*se)),high=float(np.exp(b+1.96*se)),p=float(f.pvalues.iloc[-1])))
    b=np.concatenate(bs);u=np.column_stack(us);v=u.T@u
    con=np.zeros((5,len(b)))
    for j in range(1,6):con[j-1,j*k+k-1]=1;con[j-1,k-1]=-1
    d=con@b;cv=con@v@con.T;stat=float(d@np.linalg.solve(cv,d));vo.append(dict(release=release,chi2=stat,df=5,p=float(chi2.sf(stat,5))))
    for i,j in itertools.combinations(range(6),2):
        a1=(i+1)*k-1;a2=(j+1)*k-1;d=b[a1]-b[a2];se=np.sqrt(v[a1,a1]+v[a2,a2]-2*v[a1,a2])
        vc.append(dict(release=release,taxon_A=taxa[i],taxon_B=taxa[j],delta_log_OR=d,ratio_OR=np.exp(d),low=np.exp(d-1.96*se),high=np.exp(d+1.96*se),p=float(2*norm.sf(abs(d/se)))))
vr=pd.DataFrame(vr);vc=pd.DataFrame(vc)
for release in ['R07','R06']:
    ii=vr.release.eq(release);vr.loc[ii,'q_BH_6']=multipletests(vr.loc[ii,'p'],method='fdr_bh')[1]
    ii=vc.release.eq(release);vc.loc[ii,'q_BH_15']=multipletests(vc.loc[ii,'p'],method='fdr_bh')[1]
vr.to_csv(O/'表S24d_韦荣氏球菌属重新注释关联.tsv',sep='\t',index=False);vc.to_csv(O/'表S24e_韦荣氏球菌属重新注释属间比较.tsv',sep='\t',index=False);pd.DataFrame(vo).to_csv(O/'表S24f_韦荣氏球菌属重新注释整体差异.tsv',sep='\t',index=False)
meta=dict(original_observation_metadata_keys=keys,original_ASV_taxonomy_map_available=False,primary_threshold=.8,classification_method='same SINTAX full R07/R06 runs; no new classification or selected label merging',reason='Oscillibacter exists in both reference and assignment outputs but high-confidence sample prevalence fails original 10% filter; no original per-ASV assignment allows exact old-to-new lineage tracing',original_Oscillibacter_taxonomy_rows=origin.index[origin.index.str.contains('g__Oscillibacter;',regex=False)].tolist(),input_sha256={str(p):sha(p) for p in [RAW/'4168_ASV_table.biom',RAW/'4168_annotated_feature_table.tsv',P/'R07_all_asvs_genus_assignments.tsv',P/'R06_all_asvs_genus_assignments.tsv']},notes='Veillonella direct comparisons were added after observing annotation coverage but before fitting these supplementary models; release-specific BH families 6 and 15. They test annotation sensitivity in same AGP participants, not external replication.')
(O/'分类追踪核查.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf8')
print(df.to_string(index=False),flush=True);print('Veillonella contrasts',vc.loc[vc.q_BH_15.lt(.05)].to_string(index=False),flush=True)
