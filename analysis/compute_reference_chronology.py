"""Reproduce release-date curves and source-date comparisons without plotting."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent
s=pd.read_csv(R/'structural_per_entry.csv');meta=pd.read_csv(R/'reference_metadata.csv')
meta=meta.sort_values(['holo_release','holo_pdb','target']).reset_index(drop=True)
curves=[]
for slug,*_ in json.loads((R/'methods.json').read_text()):
 vals=s[s.method==slug].set_index('target').loc[meta.target,'tier1'].to_numpy()
 curves.append(pd.DataFrame(dict(method=slug,target=meta.target,holo_release=meta.holo_release,tier1=vals,cumulative=np.cumsum(vals))))
spec=[('boltz2','2023-06-01','holo_release','lt','PDB release cutoff'),('esmfold2','2021-09-30','holo_release','le','2021 checkpoint; September end boundary'),('bioemu_refine','2023-11-23','holo_release','le','Documented PDB snapshot; not full training cutoff'),('rfd3_fix7','2024-12-16','holo_deposition','lt','Public training configuration deposition cutoff'),('bbflow_refine','2019-12-04','holo_release','le','Latest PDB release in public ATLAS training list'),('molgen_ns1.0','2023-03-29','holo_release','le','Latest PDB release in supplied SAGE v1 prepared dataset')]
dates=pd.read_csv(R/'holo_dates.csv');rows=[];targets=[]
for method,cutoff,field,op,basis in spec:
 q=s[s.method.eq(method)].merge(dates,on='target',validate='one_to_one');assert len(q)==21
 q['cutoff']=cutoff;q['date_field']=field;q['basis']=basis;q['earlier_operator']=op
 dt=pd.to_datetime(q[field]);early=dt.lt(pd.Timestamp(cutoff)) if op=='lt' else dt.le(pd.Timestamp(cutoff))
 q['group']=np.where(early,'Earlier','Later');targets.append(q)
 for group in ['Earlier','Later']:
  part=q[q.group.eq(group)];n=len(part);k=int(part.tier1.sum())
  rows.append(dict(method=method,cutoff=cutoff,date_field=field,basis=basis,earlier_operator=op,group=group,recovered=k,n_pairs=n,recovery_percent=100*k/n))
pd.concat(curves).to_csv(R/'cumulative_tier1.csv',index=False)
pd.DataFrame(rows).to_csv(R/'cutoff_recovery_summary.csv',index=False)
pd.concat(targets).to_csv(R/'cutoff_recovery_per_target.csv',index=False)
