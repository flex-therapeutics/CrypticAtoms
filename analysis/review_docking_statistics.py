"""Paired docking inference with shared-protein clustering; no new docking runs."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import itertools
import json
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

A = Path(__file__).resolve().parent
SEED = 153829
B = 20000


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum(1., np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    out = np.empty_like(p)
    out[order] = adjusted
    return out


def inference(values, seed):
    # All targets sharing the prefix refer to the same protein; paired holo
    # references remain in their protein's cluster throughout resampling.
    v = values.dropna().astype(float)
    v = v[np.isfinite(v)]
    clusters = v.groupby(v.index.str.split('__').str[0])
    sums = clusters.sum().to_numpy()
    sizes = clusters.size().to_numpy()
    k, n = len(sums), len(v)
    assert 1 <= k <= 17
    # Exact protein-cluster sign flips: every pair in one protein changes sign
    # together. Statistic is the pair-weighted mean, as in the plotted dataset.
    integers = np.arange(2**k, dtype=np.uint32)
    null = np.zeros(2**k)
    for j, val in enumerate(sums):
        null += np.where((integers >> j) & 1, 1., -1.) * val
    observed = abs(sums.sum())
    p = float(np.mean(np.abs(null) >= observed - 1e-12))
    rng = np.random.default_rng(seed)
    sample = rng.integers(0, k, (B, k))
    boot = sums[sample].sum(axis=1) / sizes[sample].sum(axis=1)
    lo, hi = np.quantile(boot, [.025, .975])
    # Pair-level Wilcoxon is sensitivity only because shared proteins correlate.
    w = wilcoxon(v.to_numpy(), alternative='two-sided', zero_method='wilcox', method='auto') if np.any(v != 0) else None
    return dict(n_targets=n, n_proteins=k, mean_difference=float(v.mean()),
                median_difference=float(v.median()), mean_ci_low=float(lo),
                mean_ci_high=float(hi), cluster_signflip_p=p,
                wilcoxon_p=float(w.pvalue) if w else 1.,
                negative_pairs=int((v < 0).sum()), zero_pairs=int((v == 0).sum()),
                positive_pairs=int((v > 0).sum()), n_sign_patterns=2**k)


def main():
    import runpy
    runpy.run_path(str(A / "filtered_docking_statistics.py"), run_name="__main__")

if __name__ == "__main__": main()
