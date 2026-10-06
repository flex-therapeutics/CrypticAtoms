"""Reproduce manuscript tables in an output copy; fail on changed results."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reproduced');args=p.parse_args();out=args.output.resolve()
if out==ROOT or ROOT.is_relative_to(out):raise SystemExit('Output must not overwrite the source repository')
for d in ['analysis','controls']:
 shutil.copytree(ROOT/d,out/d,dirs_exist_ok=True)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
commands=[['compute_recovery.py'],['compute_reference_chronology.py'],['compute_docking_recovery.py'],['compute_diagnostics.py'],['compute_sensitivity.py'],['review_statistics.py','--include-controls','--no-figure'],['review_docking_statistics.py']]
for name,*options in commands:
 print('Running',name,flush=True);subprocess.run([sys.executable,str(out/'analysis'/name),*options],check=True,env=env)
checked=[]
for pattern in ['review_statistics*.csv','review_docking_statistics*.csv','joint_apo_holo*.csv','cutoff_recovery*.csv','cumulative_tier1.csv','docking_recovery*.csv','docking_nonimproving_references.csv','stable_per_target_2p0.csv','sensitivity_summary.csv','generated_ligands_geometry_per_target.csv','generated_ligands_per_target.csv','generated_ligands_size_motion_correlations.csv']:
 for ref in (ROOT/'analysis').glob(pattern):
  result=out/'analysis'/ref.name
  pd.testing.assert_frame_equal(pd.read_csv(ref),pd.read_csv(result),check_dtype=False,rtol=1e-7,atol=1e-9);checked.append(ref.name)
x=pd.concat([pd.read_csv(ROOT/'analysis/structural_per_entry.csv'),pd.read_csv(ROOT/'controls/structural_per_entry.csv')]).set_index(['method','target']).sort_index()
y=pd.read_csv(out/'analysis/recomputed_structural_per_entry.csv').set_index(['method','target']).sort_index()
assert x.index.equals(y.index)
for c in ['n_valid','tier1','tier2','tier3','n_tier1','map_rg','sigma_x','sigma_y']:np.testing.assert_allclose(x[c],y[c],rtol=1e-7,atol=1e-9)
x=pd.read_csv(ROOT/'analysis/physics_ml_cases_per_sample.csv').set_index(['target','method','sample_idx']).sort_index();y=pd.read_csv(out/'analysis/recomputed_physics_ml_cases.csv').set_index(['target','method','sample_idx']).sort_index();assert x.index.equals(y.index)
for c in ['bb_rmsd','sc_rmsd']:np.testing.assert_allclose(x[c],y[c],rtol=1e-7,atol=1e-9)
report=dict(status='passed',tables_compared=sorted(set(checked)),structural_method_target_cells=len(pd.read_csv(out/'analysis/recomputed_structural_per_entry.csv')),figure4_diagnostic_samples=len(y),rtol=1e-7,atol=1e-9)
(out/'validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
