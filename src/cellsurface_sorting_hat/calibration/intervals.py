"""Confidence intervals for sensitivity and specificity.

``wilson`` is for simple counts. ``cluster_bootstrap`` resamples homology clusters (not single
proteins), so related proteins do not make the interval too narrow. Models are not refitted; the
interval describes sampling of the calibration set only.
"""

import math
from numbers import Integral

import numpy as np


def _count(value, name):
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer (got {value!r})")
    return int(value)


def wilson(k, n, z=1.959964):
    """Wilson score interval for k successes in n trials: ``(value, lo, hi)``; None when n is 0."""
    n = _count(n, "n")
    k = _count(k, "k")
    if not isinstance(z, int | float) or isinstance(z, bool) or not math.isfinite(z) or z <= 0:
        raise ValueError(f"z must be a positive finite number (got {z!r})")
    if n < 0:
        raise ValueError(f"n must not be negative (n={n})")
    if n == 0:
        return None
    if not 0 <= k <= n:
        raise ValueError(f"k must be between 0 and n (k={k}, n={n})")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return p, max(0.0, centre - half), min(1.0, centre + half)


def _rates(y, call):
    pos, neg = y == 1, y == 0
    sens = (call & pos).sum() / pos.sum() if pos.sum() else np.nan
    spec = (~call & neg).sum() / neg.sum() if neg.sum() else np.nan
    return sens, spec


def _clean_inputs(y, call, clusters):
    y_raw, call_raw = np.asarray(y), np.asarray(call)
    if y_raw.ndim != 1 or call_raw.ndim != 1:
        raise ValueError("y and call must be one-dimensional")
    clusters = list(clusters)
    if not (len(y_raw) == len(call_raw) == len(clusters)):
        raise ValueError("y, call and clusters must have the same length")
    if len(y_raw) == 0:
        raise ValueError("input is empty: need at least one labelled protein")
    if y_raw.dtype.kind not in "biuf":
        raise ValueError("y must be numeric: 1 for a positive, 0 for a negative")
    if y_raw.dtype.kind == "f" and not np.isfinite(y_raw).all():
        raise ValueError("y holds a non-finite value")
    if not set(np.unique(y_raw).tolist()) <= {0, 1}:
        raise ValueError("y must hold 0 and 1 only")
    if call_raw.dtype.kind not in "biu" or not set(np.unique(call_raw).tolist()) <= {0, 1}:
        raise ValueError("call must be boolean (True/False or 0/1)")
    for c in clusters:
        if isinstance(c, bool) or not isinstance(c, str | int | np.integer) or c == "":
            raise ValueError(f"every protein needs a cluster label (got {c!r})")
    return y_raw.astype(int), call_raw.astype(bool), np.asarray([str(c) for c in clusters])


def cluster_bootstrap(y, call, clusters, n_boot=2000, seed=1):
    """Sensitivity and specificity with a 95% interval that respects clusters.

    ``y`` is 1 for a known positive and 0 for a known negative; ``call`` is True when the module
    called the protein. The interval is the widest of two: the percentile interval of a cluster
    bootstrap, and the Wilson interval on the number of independent clusters of that class (the
    effective sample size). The second guards against two failures of the bootstrap: zero width when
    every resample gives the same value (all positives called), and too narrow an interval when there
    are few clusters. Returns ``{"sensitivity", "specificity", "n_clusters_pos", "n_clusters_neg"}``;
    a rate is None when its class is empty. Each rate is ``{"value", "lo", "hi"}``, ``lo <= value <= hi``.
    """
    if isinstance(n_boot, bool) or not isinstance(n_boot, Integral) or n_boot < 1:
        raise ValueError(f"n_boot must be a positive integer (got {n_boot!r})")
    y, call, clusters = _clean_inputs(y, call, clusters)
    uniq, inverse = np.unique(clusters, return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(len(uniq))]
    point = _rates(y, call)
    rng = np.random.default_rng(seed)
    boot = np.full((n_boot, 2), np.nan)
    for b in range(n_boot):
        draw = rng.integers(0, len(uniq), size=len(uniq))
        idx = np.concatenate([members[i] for i in draw])
        boot[b] = _rates(y[idx], call[idx])
    n_clusters = {1: len(set(clusters[y == 1])), 0: len(set(clusters[y == 0]))}
    out = {"n_clusters_pos": n_clusters[1], "n_clusters_neg": n_clusters[0]}
    for j, (name, label) in enumerate((("sensitivity", 1), ("specificity", 0))):
        value = point[j]
        if np.isnan(value):
            out[name] = None
            continue
        col = boot[:, j][~np.isnan(boot[:, j])]
        lo, hi = np.percentile(col, [2.5, 97.5]) if len(col) else (value, value)
        nc = n_clusters[label]
        _, wlo, whi = wilson(int(round(value * nc)), nc)
        out[name] = {
            "value": float(value),
            "lo": float(min(lo, wlo, value)),
            "hi": float(max(hi, whi, value)),
        }
    return out
