"""Cluster bootstrap for Phase C (Phase C spec 4). Needs numpy.

The unit is the MMseqs2 cluster: a resample draws clusters with replacement, and every member
of a drawn cluster enters the resample as often as its cluster was drawn. The result is a
weight matrix W (B x n): W[b, i] = number of times protein i is in resample b. Metrics read W
(metrics.py), so one W serves every candidate, variant and stratum of a test set: the
differences are paired. The models are not refitted inside the bootstrap; the intervals
describe sampling of the test set only.
"""

import hashlib

import numpy as np

N_RESAMPLES = 2000
CI = (2.5, 97.5)


def seed_for(base_seed: int, name: str) -> int:
    """A fixed seed per test set: the same name and base seed give the same resamples."""
    digest = hashlib.sha256(f"{base_seed}:{name}".encode()).hexdigest()
    return int(digest[:16], 16)


def cluster_weights(cluster_ids, n_resamples: int, seed: int) -> np.ndarray:
    """Weight matrix (n_resamples x len(cluster_ids)), int32.

    Clusters are sorted by id before drawing, so the row order of the proteins does not change
    which clusters a resample holds."""
    ids = np.asarray(cluster_ids, dtype=object)
    if len(ids) == 0:
        return np.zeros((n_resamples, 0), dtype=np.int32)
    uniq, inverse = np.unique(ids.astype(str), return_inverse=True)
    k = len(uniq)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, k, size=(n_resamples, k))
    offset = (np.arange(n_resamples) * k)[:, None]
    counts = np.bincount((draws + offset).ravel(), minlength=n_resamples * k)
    counts = counts.reshape(n_resamples, k).astype(np.int32)
    return counts[:, inverse]


def interval(values) -> dict:
    """Percentile interval over the resamples where the metric is defined (not NaN)."""
    v = np.asarray(values, dtype=np.float64)
    ok = v[~np.isnan(v)]
    if len(ok) == 0:
        return {"lo": None, "hi": None, "n_defined": 0}
    lo, hi = np.percentile(ok, CI)
    return {"lo": float(lo), "hi": float(hi), "n_defined": int(len(ok))}


def half_width(ci: dict):
    """(hi - lo) / 2, or None when the interval is not defined."""
    if ci["lo"] is None or ci["hi"] is None:
        return None
    return (ci["hi"] - ci["lo"]) / 2.0
