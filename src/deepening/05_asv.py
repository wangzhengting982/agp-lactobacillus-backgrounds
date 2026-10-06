from common import *
import h5py
from scipy.sparse import csr_matrix
P=ARCH/'02_统计输入/冻结统计输入';a=pd.read_csv(P/'R07_all_asvs_genus_assignments.tsv',sep='\t');b=pd.read_csv(P/'R06_all_asvs_genus_assignments.tsv',sep='\t')
with h5py.File(RAW/'4168_ASV_table.biom','r') as h:
 dec=lambda ar:[x.decode() if isinstance(x,bytes) else str(x) for x in ar]
 sid=dec(h['sample/ids'][:]);seq=dec(h['observation/ids'][:]);m=h['observation/matrix'];mat=csr_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=(len(seq),len(sid)))[:,pd.Index(sid).get_indexer(X.index)]
assert np.array_equal(a.biom_row,np.arange(len(seq))) and np.array_equal(b.biom_row,np.arange(len(seq)))
cons=a.genus.eq(b.genus)&a.genus_bootstrap.ge(.8)&b.genus_bootstrap.ge(.8)&a.domain.eq('Bacteria')&b.domain.eq('Bacteria')&a.domain_bootstrap.ge(.8)&b.domain_bootstrap.ge(.8)
prev=np.asarray((mat>0).sum(axis=1)).ravel()/len(X);take=cons&prev.__ge__(.05)&a.genus.isin(BG+TT)
candidates=a.loc[take,['biom_row','genus','genus_bootstrap']].copy();candidates['R06_support']=b.loc[take,'genus_bootstrap'].to_numpy();candidates['positive_people']=(prev[take]*len(X)).round().astype(int);candidates['sequence']=[seq[i] for i in candidates.biom_row];save(candidates,'V01_ASV候选完整列表')
den=pd.read_pickle(O/'_all_depth_family.pkl').loc[X.index,'named_depth'];rows=[];con=[]
for _,r in candidates.iterrows():
 counts=np.asarray(mat[int(r.biom_row)].toarray()).ravel();afs=[]
 if r.genus in BG:
  zz=scale(np.log1p(pd.Series(counts,index=X.index)/den*1e5));d=XD.assign(background_ASV=zz)
  for t in TT:
   f,m=glm(d,Y[t].ge(1));afs.append(f);rows.append(effect(f,'background_ASV',dict(direction='background_ASV',ASV_row=int(r.biom_row),ASV_genus=r.genus,target=t,**m)))
  con.append(contrast(afs,'background_ASV',dict(ASV_row=int(r.biom_row),ASV_genus=r.genus)))
 else:
  for g,c in zip(BG,BC):
   d=XD.assign(background=Z[c]);f,m=glm(d,counts>=1);rows.append(effect(f,'background',dict(direction='target_ASV',ASV_row=int(r.biom_row),ASV_genus=r.genus,target=r.genus,background=g,**m)))
rr=pd.DataFrame(rows)
for direction in ['background_ASV','target_ASV']:
 k=rr.direction.eq(direction);rr.loc[k,'q_BH']=bh(rr.loc[k]).q_BH.to_numpy()
save(rr,'V02_ASV完整关联');save(bh(con),'V03_背景ASV属间差异');print('CANDIDATES',candidates.groupby('genus').size().to_dict(),flush=True);print(rr.to_string(index=False),flush=True)
