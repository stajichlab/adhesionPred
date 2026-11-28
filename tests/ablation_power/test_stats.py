"""Tests for analysis/ablation_power/stats.py.

The statistics are the part that can be wrong while looking right, so every one of these checks a
value that was worked out independently of the implementation -- by hand, or against scipy/sklearn.
"""

import numpy as np
import pytest
import stats


# --------------------------------------------------------------- concordance
def test_concordance_of_a_perfect_ranking_is_one():
    y = np.array([0, 0, 1, 1])
    assert stats.stratified_concordance(y, np.array([1.0, 2, 3, 4]), np.array(["a"] * 4)) == 1.0


def test_concordance_of_a_reversed_ranking_is_zero():
    y = np.array([0, 0, 1, 1])
    assert stats.stratified_concordance(y, np.array([4.0, 3, 2, 1]), np.array(["a"] * 4)) == 0.0


def test_a_tie_counts_as_half_a_pair():
    # one positive, one negative, identical scores -> exactly 0.5 by hand
    y = np.array([0, 1])
    got = stats.stratified_concordance(y, np.array([1.0, 1.0]), np.array(["a", "a"]))
    assert got == pytest.approx(0.5)


def test_concordance_matches_sklearn_roc_auc_within_a_single_fold():
    roc = pytest.importorskip("sklearn.metrics").roc_auc_score
    rng = np.random.default_rng(0)
    for _ in range(5):
        y = rng.integers(0, 2, 40)
        if len(np.unique(y)) < 2:
            continue
        s = rng.standard_normal(40)
        got = stats.stratified_concordance(y, s, np.array(["only"] * 40))
        assert got == pytest.approx(roc(y, s))


def test_cross_fold_pairs_are_never_compared():
    """Two folds whose scores are on opposite scales: each is perfect, so the total must be 1.0.

    A pooled AUC would be far below 1.0 here. This is the whole point of stratifying.
    """
    y = np.array([0, 1, 0, 1])
    score = np.array([0.0, 1.0, 100.0, 101.0])
    fold = np.array(["a", "a", "b", "b"])
    assert stats.stratified_concordance(y, score, fold) == 1.0


def test_concordance_is_pair_weighted_across_folds():
    # fold a: 2x2 = 4 pairs, perfect. fold b: 1x1 = 1 pair, reversed.
    y = np.array([0, 0, 1, 1, 0, 1])
    score = np.array([0.0, 1, 2, 3, 1.0, 0.0])
    fold = np.array(["a"] * 4 + ["b"] * 2)
    assert stats.stratified_concordance(y, score, fold) == pytest.approx(4 / 5)


def test_a_fold_with_one_class_contributes_nothing():
    y = np.array([0, 1, 1, 1])
    score = np.array([0.0, 1, 2, 3])
    fold = np.array(["a", "a", "b", "b"])
    # fold b has no negative, so only fold a's single perfect pair counts
    assert stats.stratified_concordance(y, score, fold) == 1.0


# --------------------------------------------------------------- average precision
def test_average_precision_matches_sklearn():
    sk = pytest.importorskip("sklearn.metrics").average_precision_score
    rng = np.random.default_rng(1)
    for _ in range(5):
        y = rng.integers(0, 2, 30)
        if y.sum() == 0:
            continue
        s = rng.standard_normal(30)
        assert stats.average_precision(y, s) == pytest.approx(sk(y, s))


def test_average_precision_by_hand():
    # ranked: pos, neg, pos -> precision at the two hits is 1/1 and 2/3
    y = np.array([1, 0, 1])
    s = np.array([3.0, 2.0, 1.0])
    assert stats.average_precision(y, s) == pytest.approx((1.0 + 2 / 3) / 2)


# --------------------------------------------------------------- per-fold report
def test_per_fold_report_exposes_prevalence_and_normalised_ap():
    y = np.array([1, 1, 0, 1, 0, 0, 0, 0])
    s = np.array([4.0, 3, 2, 1, 0.5, 0.4, 0.3, 0.2])
    fold = np.array(["a"] * 4 + ["b"] * 4)
    got = stats.per_fold_report(y, s, fold)
    assert got["a"]["prevalence"] == pytest.approx(0.75)
    assert got["b"]["prevalence"] == pytest.approx(0.0)
    assert got["a"]["pairs"] == 3  # 3 positives x 1 negative


# --------------------------------------------------------------- permutation test
def test_identical_models_give_a_large_p_value():
    rng = np.random.default_rng(2)
    y = rng.integers(0, 2, 40)
    s = rng.standard_normal(40)
    observed, p = stats.paired_permutation_test(
        y, s, s.copy(), np.arange(40), np.array(["a"] * 40), n=200, seed=0
    )
    assert observed == pytest.approx(0.0)
    assert p > 0.2


def test_a_perfectly_predictive_component_gives_a_small_p_value():
    y = np.array([0] * 15 + [1] * 15)
    rng = np.random.default_rng(3)
    worse = rng.standard_normal(30)
    better = y.astype(float) + 0.01 * rng.standard_normal(30)
    _, p = stats.paired_permutation_test(
        y, worse, better, np.arange(30), np.array(["a"] * 30), n=500, seed=0
    )
    assert p < 0.05


def test_the_p_value_is_never_zero():
    """The add-one correction matters: a permutation p of exactly 0 is not a p-value."""
    y = np.array([0] * 10 + [1] * 10)
    better = y.astype(float)
    worse = np.zeros(20)
    _, p = stats.paired_permutation_test(
        y, worse, better, np.arange(20), np.array(["a"] * 20), n=100, seed=0
    )
    assert p > 0


def test_permutation_swaps_whole_clusters_together():
    """With one cluster the permutation distribution has exactly two outcomes, -1 and +1.

    If scores were swapped per PROTEIN instead of per cluster, intermediate differences would
    appear and both p-values below would move off these exact values. That is what this pins.
    """
    y = np.array([0, 0, 1, 1])
    good = np.array([0.0, 0, 1, 1])  # ranks positives above negatives
    bad = np.array([1.0, 1, 0, 0])  # the reverse
    cluster = np.array(["one"] * 4)
    fold = np.array(["f"] * 4)

    # component is WORSE: observed is the minimum of the two outcomes, so every permutation
    # is "at least as extreme" and the one-sided p is 1.0
    observed, p = stats.paired_permutation_test(y, good, bad, cluster, fold, n=200, seed=0)
    assert observed == pytest.approx(-1.0)
    assert p == pytest.approx(1.0)

    # component is BETTER: only the no-swap half of the draws reaches the observed value,
    # so the p-value sits near one half -- the most a single cluster can ever provide
    observed, p = stats.paired_permutation_test(y, bad, good, cluster, fold, n=2000, seed=0)
    assert observed == pytest.approx(1.0)
    assert 0.4 < p < 0.6


def test_permutation_is_deterministic_for_a_fixed_seed():
    rng = np.random.default_rng(4)
    y = rng.integers(0, 2, 30)
    a, b = rng.standard_normal(30), rng.standard_normal(30)
    cluster, fold = np.arange(30), np.array(["f"] * 30)
    first = stats.paired_permutation_test(y, a, b, cluster, fold, n=200, seed=7)
    second = stats.paired_permutation_test(y, a, b, cluster, fold, n=200, seed=7)
    assert first == second


def test_permutation_rejects_too_few_draws():
    with pytest.raises(stats.StatsError):
        stats.paired_permutation_test(
            np.array([0, 1]),
            np.array([0.0, 1]),
            np.array([0.0, 1]),
            np.array([0, 1]),
            np.array(["f", "f"]),
            n=5,
        )


def test_unknown_alternative_is_refused():
    with pytest.raises(stats.StatsError):
        stats.paired_permutation_test(
            np.array([0, 1]),
            np.array([0.0, 1]),
            np.array([0.0, 1]),
            np.array([0, 1]),
            np.array(["f", "f"]),
            n=100,
            alternative="sideways",
        )


# --------------------------------------------------------------- bootstrap helpers
def test_cluster_bootstrap_keeps_cluster_members_together():
    cluster = np.array(["a", "a", "a", "b", "c"])
    rng = np.random.default_rng(0)
    for _ in range(25):
        idx = stats.cluster_bootstrap_indices(cluster, rng)
        picked = cluster[idx]
        # every time cluster 'a' appears it brings all three of its members
        assert np.sum(picked == "a") % 3 == 0


def test_cluster_bootstrap_returns_one_row_per_resampled_member():
    cluster = np.array(["a", "a", "b"])
    rng = np.random.default_rng(1)
    sizes = {len(stats.cluster_bootstrap_indices(cluster, rng)) for _ in range(30)}
    assert sizes <= {2, 3, 4, 6}  # two clusters drawn, each of size 2 or 1


def test_percentile_interval_brackets_the_bulk():
    samples = np.linspace(0, 1, 1001)
    lo, hi = stats.percentile_interval(samples, alpha=0.05)
    assert lo == pytest.approx(0.025, abs=1e-3)
    assert hi == pytest.approx(0.975, abs=1e-3)


def test_percentile_interval_ignores_non_finite_samples():
    lo, hi = stats.percentile_interval([0.0, 1.0, np.nan, np.inf], alpha=0.5)
    assert np.isfinite(lo) and np.isfinite(hi)


def test_bca_matches_percentile_for_a_symmetric_unbiased_sample():
    pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(5)
    samples = rng.standard_normal(4000)
    jack = rng.standard_normal(80)
    observed = float(np.median(samples))
    bca = stats.bca_interval(samples, observed, jack, alpha=0.05)
    pct = stats.percentile_interval(samples, alpha=0.05)
    assert bca[0] == pytest.approx(pct[0], abs=0.25)
    assert bca[1] == pytest.approx(pct[1], abs=0.25)
