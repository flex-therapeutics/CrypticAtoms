from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.spatial.distance import cdist
from rdkit import Chem,RDLogger
from rdkit.Chem import rdDistGeom
RDLogger.DisableLog('rdApp.*')
r=Path(__file__).parent;w=r.parent;d=pd.read_csv(r/'per_sample.csv');scores=pd.read_csv(r/'selected_scores.csv');man=pd.read_csv(w/'validity/input_manifest.csv');pt=Chem.GetPeriodicTable();rows=[]
def protein(p):
 xyz=[];radii=[]
 for l in Path(p).read_text().splitlines():
  if l.startswith('ATOM  ') and l[16] in [' ','A']:
   e=l[76:78].strip() or l[12:16].strip()[0]
   if e in ['H','D']:continue
   xyz.append([float(l[30:38]),float(l[38:46]),float(l[46:54])]);radii.append(pt.GetRvdw(pt.GetAtomicNumber(e.title())))
 assert xyz,p
 return np.array(xyz),np.array(radii)
for q in d[d.has_coordinates].itertuples():
 folder='apo2mol' if q.method=='apo2mol' else 'sage-flex';stem=w/'ligand_export_20260917/bundle'/q.target/folder/f'sample_{q.sample_idx:03d}';rec=dict(method=q.method,target=q.target,sample_idx=q.sample_idx)
 if q.method=='apo2mol':
  ls=stem.with_suffix('.xyz').read_text().splitlines()[2:];x=np.array([[float(v) for v in l.split()[1:4]] for l in ls if l.strip()]);z=[pt.GetAtomicNumber(l.split()[0]) for l in ls if l.strip()]
  sq=scores[(scores.plot_method==q.method)&(scores.target==q.target)&(scores.sample_idx==q.sample_idx)].iloc[0];p=next(Path(str(sq[k])) for k in ['structure_path','struct_path'] if Path(str(sq[k])).is_file())
 else:
  raw=next(Chem.MaeMolSupplier(str(stem.with_suffix('.mae')),sanitize=False,removeHs=False));ids=[a.GetIdx() for a in raw.GetAtoms() if a.GetAtomicNum()>1];x=raw.GetConformer().GetPositions()[ids];z=[raw.GetAtomWithIdx(i).GetAtomicNum() for i in ids];p=Path(man[(man.method==q.method)&(man.target==q.target)&(man.sample_idx==q.sample_idx)].iloc[0].local_path)
 px,pr=protein(p);vr=np.array([pt.GetRvdw(int(a)) for a in z]);dist=cdist(x,px);hits=dist < .75*(vr[:,None]+pr[None,:]);rec.update(n_atoms=len(x),protein_ligand_clash_percent=100*hits.any(axis=1).mean(),protein_ligand_clash_pairs=int(hits.sum()),min_protein_distance=float(dist.min()),protein_path=str(p))
 if q.chemical_valid:
  try:
   m=Chem.SDMolSupplier(str(stem.with_suffix('.sdf')),removeHs=True)[0];assert m is not None;xyz=m.GetConformer().GetPositions();assert len(xyz)==len(x)
   bounds=rdDistGeom.GetMoleculeBoundsMatrix(m,set15bounds=True,scaleVDW=True,doTriangleSmoothing=True,useMacrocycle14config=False);bonds=[tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in m.GetBonds()];top=Chem.GetDistanceMatrix(m);dd=cdist(xyz,xyz);nb=np.array([dd[i,j]>1.25*bounds[i,j] for i,j in bonds]);pairs=np.array([(i,j) for i in range(len(x)) for j in range(i+1,len(x)) if top[i,j]>2],dtype=int).reshape(-1,2);cl=[(i,j) for i,j in pairs if dd[i,j]<.70*bounds[j,i]];atoms={a for ij in cl for a in ij};rec.update(n_bonds=len(bonds),n_stretched_bonds=int(nb.sum()),stretched_bond_percent=100*nb.mean() if len(nb) else np.nan,intraligand_clash_percent=100*len(atoms)/len(x),n_internal_clash_pairs=len(cl),fragments=len(Chem.GetMolFrags(m)),graph_status='ok')
  except Exception as e:rec.update(graph_status='error',error=str(e))
 else:rec['graph_status']='unparseable'
 rows.append(rec)
out=pd.DataFrame(rows);out['bond_issue_percent']=np.where(out.graph_status.eq('ok'),100*((out.fragments>1)|(out.n_stretched_bonds>0)),np.nan);out.to_csv(r/'geometry_per_sample.csv',index=False);cols=['bond_issue_percent','stretched_bond_percent','intraligand_clash_percent','protein_ligand_clash_percent'];summary=out.groupby(['method','target'])[cols].mean().reset_index();summary.to_csv(r/'geometry_per_target.csv',index=False);print(summary.groupby('method')[cols].agg(['median','max']).to_string());print(out.groupby(['method','graph_status']).size());print('fragments',out.groupby('method').fragments.max())
(r/'geometry_definitions.json').write_text(json.dumps(dict(bond_issue_percent='percentage of parseable ligands with more than one connected component or at least one overstretched bond',stretched_bonds='distance > 1.25 times RDKit distance-geometry upper bound; percentage of ligand heavy-atom bonds',within_ligand_clash='nonbonded heavy-atom pairs separated by >2 graph edges with distance <0.70 times RDKit DG lower bound; percentage of ligand heavy atoms in at least one such pair',protein_ligand_clash='distance <0.75 times sum of RDKit vdW radii; percentage of ligand heavy atoms contacting protein heavy atoms',aggregation='mean within each target, then show all 21 target means',missing_graphs='unavailable, never interpreted as passing',apo2mol_limitation='reconstructed connectivity is conditional on coordinates and cannot identify all absent bonds'),indent=2))
