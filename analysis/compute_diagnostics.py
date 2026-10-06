"""Aggregate supplied stable-residue, generated-ligand and Figure 4 diagnostics."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parent
s=pd.read_csv(R/'stable_per_sample_2p0.csv');cols=['preserved_'+x+'_percent' for x in ['apo','holo','both','either']]
s.groupby(['method','target']).agg(n_samples=('sample_idx','size'),n_scored=('preserved_apo_percent','count'),n_stable=('n_stable','first'),mean_coverage=('coverage_percent','mean'),min_coverage=('coverage_percent','min'),**{x:(x,'mean') for x in cols}).reset_index().to_csv(R/'stable_per_target_2p0.csv',index=False)
g=pd.read_csv(R/'generated_ligands_geometry_per_sample.csv');cols=['bond_issue_percent','stretched_bond_percent','intraligand_clash_percent','protein_ligand_clash_percent']
g.groupby(['method','target'])[cols].mean().reset_index().to_csv(R/'generated_ligands_geometry_per_target.csv',index=False)
d=pd.read_csv(R/'generated_ligands_per_sample.csv');cols=['heavy_atom_ratio','mw_ratio','tanimoto','rg_ratio']
d.groupby(['method','target'])[cols].median().reset_index().to_csv(R/'generated_ligands_per_target.csv',index=False)
corr=[]
for (method,target),g in d.groupby(['method','target']):
 q=g[['ligand_rg','Dsa']].dropna();rho=float(spearmanr(q.ligand_rg,q.Dsa).statistic) if len(q)>=10 and q.ligand_rg.nunique()>1 and q.Dsa.nunique()>1 else np.nan
 corr.append(dict(method=method,target=target,n=len(q),spearman_rg_Dsa=rho))
pd.DataFrame(corr).to_csv(R/'generated_ligands_size_motion_correlations.csv',index=False)
d=pd.read_csv(R/'physics_ml_cases_per_residue.csv');q=d.groupby(['target','method','sample_idx'])[['bb_sum','bb_n','sc_sum','sc_n']].sum();q['bb_rmsd']=np.sqrt(q.bb_sum/q.bb_n);q['sc_rmsd']=np.sqrt(q.sc_sum/q.sc_n)
q.reset_index().to_csv(R/'recomputed_physics_ml_cases.csv',index=False)
