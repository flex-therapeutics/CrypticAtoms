"""Reproduce Figure 3f and its threshold sensitivity from filtered docking data.

Run from any directory with Python, NumPy, and pandas installed. Inputs and
outputs reside beside this script. No docking or new sample selection occurs.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
ORDER = ['boltz2', 'esmfold2', 'bioemu_refine', 'confornet_nomsa',
         'bbflow_refine', 'rfd3_fix7', 'apo2mol', 'molgen_ns1.0',
         'prody_rfd3', 'rosetta_repack', 'rosetta_fastrelax']
dock = pd.read_csv(ROOT / 'clash_filtered_docking_per_target.csv')
refs = pd.read_csv(ROOT / 'docking_references.csv')
assert dock.n_pending.eq(0).all(), 'Docking must be complete'
assert not dock.duplicated(['method', 'target']).any()
assert set(dock.method) == set(ORDER)
for method, group in dock.groupby('method'):
    assert len(group) == len(refs) and set(group.target) == set(refs.target)
    assert np.allclose(group.holo_minus_apo.to_numpy(),
                       refs.set_index('target').loc[group.target, 'delta'].to_numpy())

# Lower Vina scores are preferable. No recovery fraction is defined when
# the holo reference does not improve on apo. Missing eligible samples are
# unsuccessful, rather than omitted, on the common reference-improving set.
dock['reference_improves'] = dock.holo_minus_apo.lt(0)
dock['eligible_sample'] = dock.n_docked.gt(0) & dock.delta_apo.notna()
dock['improvement_fraction'] = np.where(
    dock.reference_improves & dock.eligible_sample,
    dock.delta_apo / dock.holo_minus_apo, np.nan)
rows = []
for threshold in [0.5, 0.75, 1.0]:
    key = f'recovered_{int(threshold * 100)}'
    dock[key] = (dock.reference_improves & dock.eligible_sample &
                 dock.delta_apo.le(threshold * dock.holo_minus_apo))
    for method in ORDER:
        group = dock[dock.method.eq(method) & dock.reference_improves]
        rows.append(dict(method=method, threshold=threshold,
                         n_reference_pairs=len(group),
                         n_eligible=int(group.eligible_sample.sum()),
                         n_no_eligible=int((~group.eligible_sample).sum()),
                         n_recovered=int(group[key].sum())))
summary = pd.DataFrame(rows)
summary.to_csv(ROOT / 'docking_recovery_counts.csv', index=False)
dock.to_csv(ROOT / 'docking_recovery_per_target.csv', index=False)
refs[refs.delta.ge(0)].to_csv(ROOT / 'docking_nonimproving_references.csv', index=False)
print(summary[summary.threshold.eq(0.75)].to_string(index=False))
