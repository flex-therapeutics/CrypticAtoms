"""Generated-ligand diagnostics; requires requirements-geometry.txt.

Apo2Mol connectivity must be reconstructed before use. Unparseable graphs
remain missing; coordinates alone do not establish chemical validity.
"""
import numpy as np
from scipy.spatial.distance import cdist
from rdkit import Chem, DataStructs
from rdkit.Chem import rdDistGeom, rdFingerprintGenerator

def radius_of_gyration(xyz):
    xyz=np.asarray(xyz)
    return float(np.sqrt(np.mean(np.sum((xyz-xyz.mean(0))**2,axis=1))))

def tanimoto(molecule, reference):
    fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
    return float(DataStructs.TanimotoSimilarity(fp.GetFingerprint(molecule),fp.GetFingerprint(reference)))

def ligand_geometry(molecule, protein_xyz, protein_elements):
    """Heavy-atom molecule with a conformer and the generated protein coordinates.

    Keeps the original 1.25 upper-bound, 0.70 lower-bound and 0.75 vdW criteria.
    A graph failure must be recorded separately, not passed as an empty molecule.
    """
    m=Chem.RemoveHs(molecule);Chem.SanitizeMol(m)
    xyz=m.GetConformer().GetPositions();n=len(xyz)
    if not n:raise ValueError('No ligand heavy atoms')
    pt=Chem.GetPeriodicTable();vr=np.array([pt.GetRvdw(a.GetAtomicNum()) for a in m.GetAtoms()])
    pr=np.array([pt.GetRvdw(pt.GetAtomicNumber(str(e).title())) for e in protein_elements])
    dist=cdist(xyz,np.asarray(protein_xyz));hits=dist<.75*(vr[:,None]+pr[None,:])
    bounds=rdDistGeom.GetMoleculeBoundsMatrix(m,set15bounds=True,scaleVDW=True,doTriangleSmoothing=True,useMacrocycle14config=False)
    bonds=[tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in m.GetBonds()]
    top=Chem.GetDistanceMatrix(m);dd=cdist(xyz,xyz)
    stretched=[dd[i,j]>1.25*bounds[i,j] for i,j in bonds]
    clashes=[(i,j) for i in range(n) for j in range(i+1,n) if top[i,j]>2 and dd[i,j]<.70*bounds[j,i]]
    atoms={a for pair in clashes for a in pair};fragments=len(Chem.GetMolFrags(m))
    return dict(n_atoms=n,n_bonds=len(bonds),n_stretched_bonds=int(sum(stretched)),stretched_bond_percent=100*np.mean(stretched) if bonds else np.nan,intraligand_clash_percent=100*len(atoms)/n,n_internal_clash_pairs=len(clashes),fragments=fragments,bond_issue_percent=100*int(fragments>1 or any(stretched)),protein_ligand_clash_percent=100*hits.any(axis=1).mean(),protein_ligand_clash_pairs=int(hits.sum()),ligand_rg=radius_of_gyration(xyz))
