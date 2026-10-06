from pathlib import Path
import json,numpy as np,pandas as pd
from rdkit import Chem,DataStructs,RDLogger,rdBase
from rdkit.Chem import Descriptors,rdFingerprintGenerator
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
r=Path(__file__).resolve().parent;w=r.parent;bundle=w/'ligand_export_20260917/bundle';src=w.parent/'cryptic-pockets-benchmark'
fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
targets=json.loads((w/'fig2/targets.json').read_text());scores=pd.read_csv(r/'selected_scores.csv');assert not scores.duplicated(['plot_method','target','sample_idx']).any()
def xyzpdb(p):
 return np.array([[float(x[30:38]),float(x[38:46]),float(x[46:54])] for x in p.read_text().splitlines() if x.startswith(('ATOM  ','HETATM')) and x[76:78].strip() not in ['H','D']])
def rg(x):return float(np.sqrt(np.mean(np.sum((x-x.mean(0))**2,axis=1))))
refs={};refrows=[]
for t in targets:
 name=t['name'];ccd=t['ligands']['drug_ccd'];p=src/'docking/ligands'/f'{ccd}_ideal.sdf';m=Chem.SDMolSupplier(str(p),removeHs=True)[0];assert m is not None,(name,p)
 x=xyzpdb(r/'references'/f'{name}.pdb');assert len(x)>0
 refs[name]=(m,fp.GetFingerprint(m),rg(x));refrows.append(dict(target=name,ccd=ccd,heavy_atoms=m.GetNumHeavyAtoms(),crystal_heavy_atoms=len(x),molecular_weight=Descriptors.MolWt(m),ligand_rg=rg(x)))
refdf=pd.DataFrame(refrows).set_index('target');rows=[]
for t in targets:
 name=t['name'];refm,reffp,refrg=refs[name];rr=refdf.loc[name]
 for method,folder in [('apo2mol','apo2mol'),('molgen_ns1.0','sage-flex')]:
  for i in range(30):
   stem=bundle/name/folder/f'sample_{i:03d}';row=dict(target=name,method=method,sample_idx=i,has_coordinates=False,chemical_valid=False,graph_source='Apo2Mol official reconstruction' if method=='apo2mol' else 'native Maestro bonds',error='')
   x=None;m=None
   try:
    if method=='apo2mol':
     lines=stem.with_suffix('.xyz').read_text().splitlines();x=np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:] if l.strip()]);assert len(x)==int(lines[0])
    else:
     raw=next(Chem.MaeMolSupplier(str(stem.with_suffix('.mae')),sanitize=False,removeHs=False));assert raw is not None
     ids=[a.GetIdx() for a in raw.GetAtoms() if a.GetAtomicNum()>1];x=raw.GetConformer().GetPositions()[ids]
    assert len(x)>0 and np.isfinite(x).all();row.update(has_coordinates=True,heavy_atoms=len(x),heavy_atom_ratio=len(x)/rr.heavy_atoms,ligand_rg=rg(x),rg_ratio=rg(x)/refrg)
    if stem.with_suffix('.sdf').exists():
     m=Chem.SDMolSupplier(str(stem.with_suffix('.sdf')),sanitize=True,removeHs=True)[0]
    if m is not None:
     Chem.SanitizeMol(m);nfrag=len(Chem.GetMolFrags(m));row.update(chemical_valid=True,fragments=nfrag,molecular_weight=Descriptors.MolWt(m),mw_ratio=Descriptors.MolWt(m)/rr.molecular_weight,tanimoto=DataStructs.TanimotoSimilarity(fp.GetFingerprint(m),reffp),smiles=Chem.MolToSmiles(m))
     assert m.GetNumHeavyAtoms()==len(x),(name,method,i,'atom mismatch')
    else:row['error']='no sanitizable molecular graph'
   except Exception as e:row['error']=str(e)
   rows.append(row)
d=pd.DataFrame(rows);metric=scores[['plot_method','target','sample_idx','aa_rmsd_moving_vs_apo','aa_rmsd_moving_vs_holo','tier1']].rename(columns={'plot_method':'method','aa_rmsd_moving_vs_apo':'Dsa','aa_rmsd_moving_vs_holo':'Dsh'})
d=d.merge(metric,on=['method','target','sample_idx'],how='left',validate='one_to_one');d.to_csv(r/'per_sample.csv',index=False);refdf.to_csv(r/'references.csv')
cols=['heavy_atom_ratio','mw_ratio','tanimoto','rg_ratio'];summary=d.groupby(['method','target'])[cols].median().reset_index();summary.to_csv(r/'per_target.csv',index=False)
corr=[]
for (method,target),g in d.groupby(['method','target']):
 q=g[['ligand_rg','Dsa']].dropna();rho=float(spearmanr(q.ligand_rg,q.Dsa).statistic) if len(q)>=10 and q.ligand_rg.nunique()>1 and q.Dsa.nunique()>1 else np.nan
 corr.append(dict(method=method,target=target,n=len(q),spearman_rg_Dsa=rho))
pd.DataFrame(corr).to_csv(r/'size_motion_correlations.csv',index=False)
counts=d.groupby('method').agg(n_expected=('sample_idx','size'),n_coordinates=('has_coordinates','sum'),n_sanitizable=('chemical_valid','sum'));counts.to_csv(r/'validity_counts.csv');print(counts.to_string());print(summary.groupby('method')[cols].median().to_string());print(pd.DataFrame(corr).groupby('method').spearman_rg_Dsa.agg(['median','min','max']).to_string());print('ref atom count differences',refdf[refdf.heavy_atoms.ne(refdf.crystal_heavy_atoms)].to_string())
(r/'analysis_config.json').write_text(json.dumps(dict(rdkit_version=rdBase.rdkitVersion,morgan_radius=2,fingerprint_bits=2048,use_chirality=False,reference_chemistry='CCD ideal SDF',reference_extent='crystal ligand in approved Figure 2 session',coordinate_radius='unweighted heavy-atom RMS distance to centroid',fingerprint_filter='sanitizable molecular graph; no repairs',summary_unit='median per reference pair',correlation_unit='within target and method, across samples',apo2mol_atom_count='reference atom count supplied to generator'),indent=2))
