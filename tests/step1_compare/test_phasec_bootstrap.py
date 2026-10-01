"""bootstrap.py: whole clusters move together, seeds, pairing, empty resamples (spec 4, 6)."""

import pytest

np = pytest.importorskip("numpy")

import bootstrap as bs  # noqa: E402
import metrics as m  # noqa: E402

CLUSTERS = ["c1", "c1", "c1", "c2", "c3", "c3", "c4", "c5", "c5", "c5"]


def test_bootstrap_resamples_whole_clusters():
    W = bs.cluster_weights(CLUSTERS, 500, seed=1)
    assert W.shape == (500, 10) and W.dtype == np.int32
    for members in ([0, 1, 2], [4, 5], [7, 8, 9]):
        assert (W[:, members] == W[:, members[:1]]).all()
    # each resample draws as many clusters as there are (5)
    first = [0, 3, 4, 6, 7]  # one member per cluster
    assert (W[:, first].sum(axis=1) == 5).all()


def test_bootstrap_seeded():
    a = bs.cluster_weights(CLUSTERS, 50, seed=11)
    b = bs.cluster_weights(CLUSTERS, 50, seed=11)
    c = bs.cluster_weights(CLUSTERS, 50, seed=12)
    assert (a == b).all() and not (a == c).all()
    assert bs.seed_for(20261001, "S1:all") == bs.seed_for(20261001, "S1:all")
    assert bs.seed_for(20261001, "S1:all") != bs.seed_for(20261001, "S1:Scer_SGD")


def test_row_order_does_not_change_the_cluster_draws():
    perm = [9, 0, 5, 3, 1, 8, 2, 7, 4, 6]
    a = bs.cluster_weights(CLUSTERS, 30, seed=5)
    b = bs.cluster_weights([CLUSTERS[i] for i in perm], 30, seed=5)
    assert (a[:, perm] == b).all()


def test_paired_bootstrap_uses_same_indices():
    y = np.array([1, 1, 0, 1, 0, 0, 1, 0, 0, 0], dtype=bool)
    s = np.array([0.9, 0.8, 0.8, 0.7, 0.6, 0.5, 0.5, 0.4, 0.3, 0.1])
    W = bs.cluster_weights(CLUSTERS, 400, seed=2)
    a = m.roc_auc(W, y, s)
    b = m.roc_auc(W, y, s)  # a second candidate with the same scores
    diff = bs.interval(a - b)
    assert diff["lo"] == 0.0 and diff["hi"] == 0.0  # paired: identical candidates differ by 0
    unpaired = m.roc_auc(bs.cluster_weights(CLUSTERS, 400, seed=3), y, s)
    assert bs.interval(a - unpaired)["hi"] > 0  # different indices do not cancel


def test_interval_skips_undefined_resamples():
    # Review Focus 2: a resample without positives has NaN recall; the interval uses the rest
    vals = np.array([np.nan, 0.5, 0.6, 0.7, np.nan])
    ci = bs.interval(vals)
    assert (
        ci["n_defined"] == 3
        and ci["lo"] == pytest.approx(0.505)
        and ci["hi"] == pytest.approx(0.695)
    )
    assert bs.interval([np.nan, np.nan]) == {"lo": None, "hi": None, "n_defined": 0}
    assert bs.half_width({"lo": None, "hi": None, "n_defined": 0}) is None
    assert bs.half_width({"lo": 0.6, "hi": 0.8, "n_defined": 9}) == pytest.approx(0.1)


def test_all_positives_in_one_cluster_gives_some_empty_resamples():
    y = np.array([1, 1, 0, 0, 0, 0], dtype=bool)
    clusters = ["p", "p", "a", "b", "c", "d"]
    W = bs.cluster_weights(clusters, 300, seed=4)
    rec = m.recall(W, y, np.array([1, 0, 0, 0, 0, 0], dtype=bool))
    assert np.isnan(rec).any() and not np.isnan(rec).all()
    ci = bs.interval(rec)
    assert 0 < ci["n_defined"] < 300 and ci["lo"] == ci["hi"] == 0.5


def test_empty_test_set_gives_an_empty_matrix():
    assert bs.cluster_weights([], 7, seed=1).shape == (7, 0)
