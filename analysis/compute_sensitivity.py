"""Summarize the finalized paired chemical-correspondence sensitivity audit."""
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parent
v=pd.read_csv(R/'sensitivity_paired.csv');summary=[]
for method,g in v.groupby('method'):
 rec=dict(method=method,n_compared=len(g),targets_compared=g.target.nunique(),median_abs_dh_change=abs(g.chemical_dh-g.stored_dh).median(),max_abs_dh_change=abs(g.chemical_dh-g.stored_dh).max())
 for t in [1,2,3]:
  rec[f'sample_flips_t{t}']=int((g[f'old_t{t}']!=g[f'new_t{t}']).sum());p=g.groupby('target')[[f'old_t{t}',f'new_t{t}']].any();rec[f'entry_flips_t{t}']=int((p.iloc[:,0]!=p.iloc[:,1]).sum());rec[f'old_entries_t{t}']=int(p.iloc[:,0].sum());rec[f'new_entries_t{t}']=int(p.iloc[:,1].sum())
 summary.append(rec)
pd.DataFrame(summary).to_csv(R/'sensitivity_summary.csv',index=False)
