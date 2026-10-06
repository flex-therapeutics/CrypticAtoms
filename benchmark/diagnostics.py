"""Coordinate diagnostics using the manuscript's published definitions."""
import numpy as np
import biotite.structure as st
from scipy.spatial import cKDTree
from .metrics import match_atoms, _kabsch_align_residue, _hungarian_pair_within_residue
from scipy.optimize import linear_sum_assignment
HEAVY=dict(ALA=5,ARG=11,ASN=8,ASP=8,CYS=6,GLN=9,GLU=9,GLY=4,HIS=10,ILE=8,LEU=8,LYS=9,MET=8,PHE=11,PRO=7,SER=6,THR=7,TRP=14,TYR=12,VAL=7)
def residue(a,r):
    return a[(a.res_id==r)&(a.atom_name!='OXT')]
def complete(a):
    return len(a)>0 and len(set(a.res_name))==1 and len(a)==HEAVY.get(str(a.res_name[0]),-1) and len(set(a.atom_name))==len(a) and {'N','CA','C','O'}<=set(a.atom_name)
def maximum_displacement(a,b):
    if not complete(a) or not complete(b) or set(a.res_name)!=set(b.res_name) or set(a.atom_name)!=set(b.atom_name):
        return np.nan
    ai,bi=_hungarian_pair_within_residue(a,b)
    if len(ai)!=len(a) or len(bi)!=len(b):return np.nan
    return float(np.linalg.norm(a.coord[ai]-b.coord[bi],axis=1).max())
def stability(sample_in_apo, sample_in_holo, apo, holo, holo_in_apo, pocket):
    stable=[r for r in pocket if maximum_displacement(residue(holo_in_apo,r),residue(apo,r))<2]
    da=np.array([maximum_displacement(residue(sample_in_apo,r),residue(apo,r)) for r in stable])
    dh=np.array([maximum_displacement(residue(sample_in_holo,r),residue(holo,r)) for r in stable])
    ok=np.isfinite(da)&np.isfinite(dh)
    return dict(n_stable=len(stable),n_evaluable_stable=int(ok.sum()),stable_coverage_percent=100*ok.sum()/len(stable) if stable else np.nan,preserved_apo_percent=100*np.mean(da[ok]<2) if ok.any() else np.nan)
def pocket_clashes(sample, apo, holo, pocket):
    def keys(a):return list(zip(map(int,a.res_id),map(str,a.atom_name)))
    adj={}
    for ref in [apo,holo]:
        kk=keys(ref)
        def connect(x,y):
            adj.setdefault(x,set()).add(y);adj.setdefault(y,set()).add(x)
        for i,j,_ in st.connect_via_residue_names(ref,inter_residue=False).as_array():connect(kk[i],kk[j])
        idx={k:i for i,k in enumerate(kk)}
        for rid in sorted(set(ref.res_id)):
            x,y=(int(rid),'C'),(int(rid)+1,'N')
            if x in idx and y in idx and np.linalg.norm(ref.coord[idx[x]]-ref.coord[idx[y]])<1.8:connect(x,y)
        sg=np.where((ref.res_name=='CYS')&(ref.atom_name=='SG'))[0]
        for ii,i in enumerate(sg):
            for j in sg[ii+1:]:
                if np.linalg.norm(ref.coord[i]-ref.coord[j])<2.3:connect(kk[i],kk[j])
    excluded=set()
    for x,ns in adj.items():
        for y in ns|set().union(*(adj.get(z,set()) for z in ns)):excluded.add(frozenset((x,y)))
    p=sample[np.isin(sample.res_id,pocket)];kk=keys(p)
    assert len(set(kk))==len(kk),'duplicate atom keys'
    hits=[(i,j) for i,j in cKDTree(p.coord).query_pairs(2.) if kk[i][0]!=kk[j][0] and frozenset((kk[i],kk[j])) not in excluded] if len(p) else []
    residues=set(p.res_id);clashing={kk[i][0] for pair in hits for i in pair}
    return dict(n_pocket_residues=len(residues),n_clashing_residues=len(clashing),residue_clash_percent=100*len(clashing)/len(residues) if residues else np.nan,n_clash_pairs=len(hits))

def backbone_sidechain_sums(a,b,moving):
 bb={"N","CA","C","O"}
 vals=[]
 for rid in moving:
  aa=a[(a.res_id==rid)&(a.atom_name!='OXT')];hh=b[(b.res_id==rid)&(b.atom_name!='OXT')];assert set(aa.atom_name)==set(hh.atom_name),(rid,'incomplete atoms');assert set(aa.res_name)==set(hh.res_name),(rid,'residue mismatch')
  ab,hb=match_atoms(aa[np.isin(aa.atom_name,list(bb))],hh[np.isin(hh.atom_name,list(bb))]);bsq=((ab.coord-hb.coord)**2).sum(1);al=_kabsch_align_residue(aa,hh);sc=[]
  for el in set(aa.element):
   ai=np.where((aa.element==el)&~np.isin(aa.atom_name,list(bb)))[0];hi=np.where((hh.element==el)&~np.isin(hh.atom_name,list(bb)))[0];assert len(ai)==len(hi)
   if len(ai):
    ii,jj=linear_sum_assignment(np.linalg.norm(al.coord[ai,None,:]-hh.coord[hi][None,:,:],axis=-1));sc.extend(((aa.coord[ai[ii]]-hh.coord[hi[jj]])**2).sum(1).tolist())
  vals.append(dict(res_id=int(rid),res_name=str(hh.res_name[0]),bb_sum=float(bsq.sum()),bb_n=len(bsq),sc_sum=float(sum(sc)),sc_n=len(sc)))
 return vals
