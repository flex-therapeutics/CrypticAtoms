"""Recompute tiers, joint apo/holo recovery and pocket-map spread from scores."""
from pathlib import Path
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent
s=pd.concat([pd.read_csv(R/'selected_structural_scores.csv'),pd.read_csv(R.parent/'controls/selected_structural_scores.csv')],ignore_index=True)
assert not s.duplicated(['plot_method','target','sample_idx']).any()
meta=pd.read_csv(R/'reference_metadata.csv').set_index('target');rows=[];joint=[];fractions=[]
for (method,target),g in s.groupby(['plot_method','target']):
 a=g.aa_rmsd_moving_vs_apo.to_numpy();h=g.aa_rmsd_moving_vs_holo.to_numpy();D=float(meta.loc[target,'D'])
 assert np.isfinite(a).all() and np.isfinite(h).all() and D>0
 x=(a*a+D*D-h*h)/(2*D);y=np.sqrt(np.maximum(0,a*a-x*x))
 row=dict(method=method,target=target,n_valid=len(g),map_rg=float(np.sqrt(np.mean((x-x.mean())**2+(y-y.mean())**2))),sigma_x=float(np.std(x,ddof=1)),sigma_y=float(np.std(y,ddof=1)))
 for tier,cut in [(1,np.inf),(2,1.5),(3,1.)]:
  apo=(a<D)&(a<h)&(a<cut);holo=(h<D)&(h<a)&(h<cut)
  row[f'tier{tier}']=int(holo.any());row[f'n_tier{tier}']=int(holo.sum())
  fractions.append(dict(method=method,target=target,tier=tier,n_success=int(holo.sum()),budget=30,success_percent=100*holo.sum()/30))
  joint.append(dict(method=method,target=target,protein=target.split('__')[0],tier=tier,n_samples=len(g),apo_recovered=int(apo.any()),holo_recovered=int(holo.any()),both_recovered=int(apo.any() and holo.any()),apo_samples=','.join(map(str,g.loc[apo,'sample_idx'].astype(int))),holo_samples=','.join(map(str,g.loc[holo,'sample_idx'].astype(int)))))
 rows.append(row)
pd.DataFrame(rows).to_csv(R/'recomputed_structural_per_entry.csv',index=False)
pd.DataFrame(fractions).to_csv(R/'recomputed_success_fractions.csv',index=False)
r=pd.DataFrame(joint);r.to_csv(R/'joint_apo_holo_per_pair.csv',index=False)
r.groupby(['method','tier'])[['apo_recovered','holo_recovered','both_recovered']].sum().to_csv(R/'joint_apo_holo_counts.csv')
