"""metrics.py: hand-computed values, ties, weights, undefined cases (spec 4, section 6)."""

import math

import pytest

np = pytest.importorskip("numpy")
sk = pytest.importorskip("sklearn.metrics")

import metrics as m  # noqa: E402

# 10 proteins; scores with two tie groups (0.8: pos + neg; 0.5: neg + pos).
Y = np.array([1, 1, 0, 1, 0, 0, 1, 0, 0, 0], dtype=bool)
S = np.array([0.9, 0.8, 0.8, 0.7, 0.6, 0.5, 0.5, 0.4, 0.3, 0.1])
CALL = S >= 0.6  # rows 0-4
ONE = np.ones(10)


def test_metrics_hand_computed():
    # P = 4, N = 6. Called: rows 0-4 -> TP 3 (rows 0, 1, 3), FP 2 (rows 2, 4).
    assert m.recall(ONE, Y, CALL)[0] == 0.75
    assert m.precision(ONE, Y, CALL)[0] == 0.6
    assert m.fpr(ONE, Y, CALL)[0] == pytest.approx(2 / 6, abs=1e-15)
    # ROC-AUC by pairs: 6 + 5.5 + 5 + 3.5 = 20 of 24 pairs (ties count 0.5)
    assert m.roc_auc(ONE, Y, S)[0] == pytest.approx(20 / 24, abs=1e-15)
    # AP: 0.25 * (1 + 2/3 + 3/4 + 4/7) = 251/336
    assert m.pr_auc(ONE, Y, S)[0] == pytest.approx(251 / 336, abs=1e-15)


def test_precision_at_recall_reads_the_step_curve_with_ties():
    # thresholds (TP, FP): 0.9 (1,0) 0.8 (2,1) 0.7 (3,1) 0.6 (3,2) 0.5 (4,3) ...
    # recall >= 0.5 first at 0.8 (precision 2/3), but 0.7 gives 3/4: the highest wins
    assert m.precision_at_recall(ONE, Y, S, 0.5)[0] == pytest.approx(0.75, abs=1e-15)
    assert m.precision_at_recall(ONE, Y, S, 0.8)[0] == pytest.approx(4 / 7, abs=1e-15)
    assert m.precision_at_recall(ONE, Y, S, 0.9)[0] == pytest.approx(4 / 7, abs=1e-15)


def test_recall_at_fpr():
    assert m.recall_at_fpr(ONE, Y, S, 0.01)[0] == 0.25  # only the 0.9 threshold has FPR 0
    assert m.recall_at_fpr(ONE, Y, S, 0.17)[0] == 0.75  # FPR 1/6 at 0.8 and at 0.7
    assert m.recall_at_fpr(ONE, Y, S, 1.0)[0] == 1.0
    # review I-3: the top-scored protein is a negative, so every non-empty call set has FPR > 0;
    # the empty call set (recall 0, FPR 0) is the only one within FPR 0.01: the result is 0.0
    y = np.array([0, 1, 1, 0, 0], dtype=bool)
    s = np.array([0.9, 0.8, 0.7, 0.2, 0.1])
    got = m.recall_at_fpr(np.ones(5), y, s, 0.01)[0]
    assert got == 0.0 and not math.isnan(got)
    assert m.scored_summary(np.ones(5), y, s)["recall_at_fpr_0.01"][0] == 0.0


def test_fpr_at_recall_inside_a_stratum():
    nsec = np.zeros(10, dtype=bool)
    nsec[[4, 5]] = True  # two negatives with scores 0.6 and 0.5
    all_neg = ~Y
    assert m.fpr_at_recall(ONE, Y, S, 0.75, all_neg)[0] == pytest.approx(1 / 6, abs=1e-15)
    assert m.fpr_at_recall(ONE, Y, S, 0.75, nsec)[0] == 0.0  # threshold 0.7 calls neither
    assert m.fpr_at_recall(ONE, Y, S, 1.0, nsec)[0] == 1.0  # threshold 0.5 calls both


def test_brier_and_prevalence():
    assert m.brier(np.ones(2), [1, 0], [0.9, 0.1])[0] == pytest.approx(0.01, abs=1e-15)
    # recall 0.8, FPR 0.1 at prevalence 0.05: 0.04 / (0.04 + 0.095)
    got = m.precision_at_prevalence(0.8, 0.1, 0.05)
    assert got == pytest.approx(0.04 / 0.135, abs=1e-15)


def test_weights_equal_duplicated_rows():
    w = np.array([2, 1, 1, 3, 1, 0, 1, 2, 1, 1], dtype=float)
    idx = np.repeat(np.arange(10), w.astype(int))
    for f in (m.roc_auc, m.pr_auc):
        assert f(w, Y, S)[0] == pytest.approx(f(np.ones(len(idx)), Y[idx], S[idx])[0], abs=1e-12)
    assert m.recall(w, Y, CALL)[0] == pytest.approx(
        m.recall(np.ones(len(idx)), Y[idx], CALL[idx])[0]
    )


def test_weighted_values_match_sklearn():
    rng = np.random.default_rng(3)
    for _ in range(20):
        n = 60
        y = rng.random(n) < 0.3
        s = np.round(rng.random(n), 1)  # many ties
        w = rng.integers(0, 4, size=n).astype(float)
        if w[y].sum() == 0 or w[~y].sum() == 0:
            continue
        assert m.roc_auc(w, y, s)[0] == pytest.approx(sk.roc_auc_score(y, s, sample_weight=w))
        assert m.pr_auc(w, y, s)[0] == pytest.approx(
            sk.average_precision_score(y, s, sample_weight=w)
        )


def test_rows_of_w_are_independent():
    W = np.stack([ONE, np.r_[np.ones(5), np.zeros(5)]])
    got = m.recall(W, Y, CALL)
    assert got[0] == 0.75 and got[1] == 1.0  # second row: positives 0, 1, 3, all called


def test_undefined_metrics_are_nan():
    # Review Focus 1: a stratum with no positives or no negatives gives NaN, not an error
    neg_only = np.zeros(4, dtype=bool)
    s = np.array([0.1, 0.2, 0.3, 0.4])
    assert math.isnan(m.recall(np.ones(4), neg_only, s > 0.2)[0])
    assert m.fpr(np.ones(4), neg_only, s > 0.2)[0] == 0.5
    for f in (m.roc_auc, m.pr_auc):
        assert math.isnan(f(np.ones(4), neg_only, s)[0])
    assert math.isnan(m.precision_at_recall(np.ones(4), neg_only, s, 0.8)[0])
    pos_only = np.ones(4, dtype=bool)
    assert math.isnan(m.fpr(np.ones(4), pos_only, s > 0.2)[0])
    assert math.isnan(m.roc_auc(np.ones(4), pos_only, s)[0])
    assert m.pr_auc(np.ones(4), pos_only, s)[0] == 1.0
    assert math.isnan(m.precision(np.ones(4), pos_only, np.zeros(4, dtype=bool))[0])
    assert math.isnan(m.recall(np.zeros((1, 0)), np.zeros(0, bool), np.zeros(0, bool))[0])


def test_non_finite_scores_are_refused():
    # Review Focus 3: a NaN score must not be sorted silently
    with pytest.raises(ValueError, match="finite"):
        m.roc_auc(ONE, Y, np.r_[S[:9], np.nan])


def test_level_can_differ_per_resample():
    W = np.stack([ONE, ONE])
    got = m.precision_at_recall(W, Y, S, np.array([0.5, 0.8]))
    assert got == pytest.approx([0.75, 4 / 7])


def test_scored_summary_equals_the_single_functions():
    W = np.stack([ONE, np.r_[np.ones(5), 2 * np.ones(5)]])
    got = m.scored_summary(W, Y, S, (0.8, 0.9), 0.01)
    assert set(got) == {"roc_auc", "pr_auc", "precision_at_recall_0.8",
                        "precision_at_recall_0.9", "recall_at_fpr_0.01"}  # fmt: skip
    assert got["roc_auc"] == pytest.approx(m.roc_auc(W, Y, S))
    assert got["pr_auc"] == pytest.approx(m.pr_auc(W, Y, S))
    assert got["precision_at_recall_0.9"] == pytest.approx(m.precision_at_recall(W, Y, S, 0.9))
    assert got["recall_at_fpr_0.01"] == pytest.approx(m.recall_at_fpr(W, Y, S, 0.01))
