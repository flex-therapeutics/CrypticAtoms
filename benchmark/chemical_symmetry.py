"""Candidate aa_mov implementation: complete heavy atoms, chemical symmetry only.
Independent of benchmark/metrics.py; original benchmark results remain unchanged.
Coordinates must already be aligned on the declared pocket C-alpha set.
"""
import numpy as np
SWAPS={'ASP':[('OD1','OD2')],'GLU':[('OE1','OE2')],'ARG':[('NH1','NH2')],
       'PHE':[('CD1','CD2'),('CE1','CE2')],'TYR':[('CD1','CD2'),('CE1','CE2')],
       'VAL':[('CG1','CG2')],'LEU':[('CD1','CD2')]}
def matched_residue(a,b,rid):
 a=a[a.res_id==rid];b=b[b.res_id==rid]
 if not len(a) or not len(b):raise ValueError(f'missing residue {rid}')
 if len(set(a.res_name))!=1 or set(a.res_name)!=set(b.res_name):raise ValueError(f'residue identity mismatch {rid}')
 an={str(n):i for i,n in enumerate(a.atom_name) if str(n)!='OXT'};bn={str(n):i for i,n in enumerate(b.atom_name) if str(n)!='OXT'}
 if len(an)!=len(a) or len(bn)!=len(b):raise ValueError(f'duplicate or terminal atom at moving residue {rid}')
 if set(an)!=set(bn):raise ValueError(f'incomplete heavy atoms at {rid}: {set(an)^set(bn)}')
 names=sorted(an);ai=np.array([an[n] for n in names]);base=[bn[n] for n in names];options=[base]
 swaps=SWAPS.get(str(a.res_name[0]),[])
 if swaps:
  mapping={n:n for n in names}
  for x,y in swaps:mapping[x]=y;mapping[y]=x
  options.append([bn[mapping[n]] for n in names])
 best=min(options,key=lambda bi:np.square(a.coord[ai]-b.coord[bi]).sum())
 return a.coord[ai],b.coord[best]
def displacement(a,b,residues):
 return {r:float(np.linalg.norm(np.subtract(*matched_residue(a,b,r)),axis=1).max()) for r in residues}
def rmsd(a,b,residues):
 deltas=[np.subtract(*matched_residue(a,b,r)) for r in residues]
 return float(np.sqrt(np.mean(np.sum(np.concatenate(deltas)**2,axis=1))))
