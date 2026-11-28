"""Metrics and paired tests for the ablation power study.

Metric choice (decision D-C, revised 2026-11-27 after review). The primary metric is the
WITHIN-FOLD STRATIFIED CONCORDANCE, not pooled average precision:

  * Pooling out-of-fold scores across folds is incoherent here, not merely uncalibrated. Measured on
    adhesion_v1, length-only AUC is 0.912 in Saccharomycetes and 0.327 in Eurotiomycetes -- two
    folds learn "long => adhesin" and the third learns the opposite. 47% of the pairs in a pooled
    AUC compare scores from models that disagree about the sign of their main feature.
  * Stratified concordance simply drops those cross-fold pairs. It is prevalence-invariant, which
    matters because fold prevalence runs 0.36 / 0.65 / 0.31.
  * Average precision is still reported, descriptively, always beside its fold's prevalence. It is
    prevalence-dependent and has no clean within-fold-stratified form.

Test choice. The keep/remove decision uses a PAIRED PERMUTATION TEST by cluster, not a bootstrap
confidence interval: it is exact under the exchangeability null and does not depend on the shape of
the bootstrap distribution. The bootstrap is kept for the effect size and its interval, and it
REFITS inside each resample -- a bootstrap that conditions on the fitted models carries test-set
variance only and was measured to be anti-conservative (9.5% rejection at a nominal 5%).
"""

import numpy as np

MIN_PERMUTATIONS = 50


class StatsError(Exception):
    pass


def _auc(y, score):
    """Mann-Whitney concordance for one group. Ties count a half. None if a class is absent."""
    pos = score[y == 1]
    neg = score[y == 0]
    if pos.size == 0 or neg.size == 0:
        return None, 0
    order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    ranks = np.empty(order.size, dtype=float)
    ranks[order] = np.arange(1, order.size + 1, dtype=float)
    # average ranks over ties so that a tie contributes exactly a half pair
    combined = np.concatenate([pos, neg])
    _, inverse, counts = np.unique(combined, return_inverse=True, return_counts=True)
    sums = np.zeros(counts.size)
    np.add.at(sums, inverse, ranks)
    ranks = (sums / counts)[inverse]
    npos, nneg = pos.size, neg.size
    concordant = ranks[:npos].sum() - npos * (npos + 1) / 2
    return concordant / (npos * nneg), npos * nneg


def stratified_concordance(y, score, fold):
    """Pair-weighted within-fold concordance: concordant within-fold pairs / all within-fold pairs.

    Cross-fold pairs are never compared, so three separately fitted models cannot contradict
    each other on one axis.
    """
    y = np.asarray(y)
    score = np.asarray(score, dtype=float)
    fold = np.asarray(fold)
    total_pairs = 0.0
    total_concordant = 0.0
    for f in np.unique(fold):
        sel = fold == f
        auc, pairs = _auc(y[sel], score[sel])
        if auc is None:
            continue
        total_concordant += auc * pairs
        total_pairs += pairs
    if total_pairs == 0:
        return float("nan")
    return total_concordant / total_pairs


def average_precision(y, score):
    """Average precision, the step-wise sum used by sklearn (no trapezoid interpolation)."""
    y = np.asarray(y)
    score = np.asarray(score, dtype=float)
    if y.sum() == 0:
        return float("nan")
    order = np.argsort(-score, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y)
    precision = tp / np.arange(1, y.size + 1)
    return float((precision * y).sum() / y.sum())


def per_fold_report(y, score, fold):
    """AP and AUC per fold, each beside its prevalence, so a prevalence artefact is visible."""
    y = np.asarray(y)
    score = np.asarray(score, dtype=float)
    fold = np.asarray(fold)
    out = {}
    for f in np.unique(fold):
        sel = fold == f
        auc, pairs = _auc(y[sel], score[sel])
        prev = float(y[sel].mean())
        ap = average_precision(y[sel], score[sel])
        out[str(f)] = {
            "n": int(sel.sum()),
            "prevalence": prev,
            "auc": None if auc is None else float(auc),
            "ap": ap,
            # AP normalised against its own prevalence, which reorders the folds
            "ap_normalised": None if prev >= 1 else float((ap - prev) / (1 - prev)),
            "pairs": int(pairs),
        }
    return out


def paired_permutation_test(
    y, score_a, score_b, cluster, fold, n=10000, seed=0, alternative="greater"
):
    """Paired permutation test on the difference in stratified concordance (b minus a).

    For each CLUSTER independently, with probability one half the two models' scores are swapped
    for every member of that cluster. Clusters are the exchangeable unit because a paralogue family
    is one observation, not nine.

    Returns (observed_difference, p_value).
    """
    if n < MIN_PERMUTATIONS:
        raise StatsError(f"n must be at least {MIN_PERMUTATIONS}, got {n}")
    y = np.asarray(y)
    a = np.asarray(score_a, dtype=float)
    b = np.asarray(score_b, dtype=float)
    fold = np.asarray(fold)
    cluster = np.asarray(cluster)

    observed = stratified_concordance(y, b, fold) - stratified_concordance(y, a, fold)

    _, cluster_index = np.unique(cluster, return_inverse=True)
    n_clusters = cluster_index.max() + 1
    rng = np.random.default_rng(seed)

    at_least_as_extreme = 0
    for _ in range(n):
        swap = rng.random(n_clusters) < 0.5
        mask = swap[cluster_index]
        pa = np.where(mask, b, a)
        pb = np.where(mask, a, b)
        diff = stratified_concordance(y, pb, fold) - stratified_concordance(y, pa, fold)
        if alternative == "greater":
            if diff >= observed:
                at_least_as_extreme += 1
        elif alternative == "two-sided":
            if abs(diff) >= abs(observed):
                at_least_as_extreme += 1
        else:
            raise StatsError(f"unknown alternative {alternative!r}")
    # add-one correction: a permutation p is never exactly zero
    p = (at_least_as_extreme + 1) / (n + 1)
    return float(observed), float(p)


def cluster_bootstrap_indices(cluster, rng):
    """Resample clusters with replacement and return the member indices they carry."""
    cluster = np.asarray(cluster)
    uniq, inverse = np.unique(cluster, return_inverse=True)
    members = [np.flatnonzero(inverse == k) for k in range(uniq.size)]
    picked = rng.integers(0, uniq.size, size=uniq.size)
    return np.concatenate([members[k] for k in picked])


def percentile_interval(samples, alpha=0.05):
    samples = np.asarray([s for s in samples if np.isfinite(s)], dtype=float)
    if samples.size == 0:
        return float("nan"), float("nan")
    lo = float(np.percentile(samples, 100 * alpha / 2))
    hi = float(np.percentile(samples, 100 * (1 - alpha / 2)))
    return lo, hi


def bca_interval(samples, observed, jackknife, alpha=0.05):
    """Bias-corrected and accelerated interval.

    `jackknife` holds the leave-one-cluster-out statistics, which give the acceleration. Falls back
    to the percentile interval when the spread is degenerate.
    """
    from scipy.stats import norm

    samples = np.asarray([s for s in samples if np.isfinite(s)], dtype=float)
    jackknife = np.asarray([s for s in jackknife if np.isfinite(s)], dtype=float)
    if samples.size == 0:
        return float("nan"), float("nan")

    proportion = float(np.mean(samples < observed))
    if proportion <= 0 or proportion >= 1:
        return percentile_interval(samples, alpha)
    z0 = norm.ppf(proportion)

    centred = jackknife.mean() - jackknife
    denom = 6 * (np.sum(centred**2) ** 1.5)
    acc = 0.0 if denom == 0 else float(np.sum(centred**3) / denom)

    out = []
    for q in (alpha / 2, 1 - alpha / 2):
        z = norm.ppf(q)
        adjusted = z0 + (z0 + z) / (1 - acc * (z0 + z))
        out.append(float(np.percentile(samples, 100 * norm.cdf(adjusted))))
    return out[0], out[1]
