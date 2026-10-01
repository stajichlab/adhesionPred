"""models.py and universe.py: rows, nested protocol, learnability, chance level (spec 3.3, 6)."""

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import evalio  # noqa: E402
import models  # noqa: E402
import universe  # noqa: E402
from metrics import roc_auc  # noqa: E402
from sklearn.model_selection import StratifiedGroupKFold  # noqa: E402

M8, M35 = universe.MODELS


def toy_universe(n, dim, seed, y=None, signal=0.0, lengths=None, cterm=None):
    """Universe with random features; embedding column 0 carries `signal` for y == 1."""
    rng = np.random.default_rng(seed)
    y = np.zeros(n, dtype=bool) if y is None else np.asarray(y, dtype=bool)
    lengths = np.full(n, 300) if lengths is None else np.asarray(lengths)
    crow = np.full(n, -1)
    long = np.nonzero(lengths > 1022)[0]
    crow[long] = np.arange(len(long))
    nterm = rng.normal(size=(n, dim))
    nterm[:, 0] += signal * y
    u = universe.Universe(
        hashes=[f"h{i:05d}" for i in range(n)],
        length=lengths,
        sp=rng.random(n) < 0.5,
        sp_prob=rng.random(n),
        rank=rng.integers(0, 4, size=n).astype(np.int8),
        gpi_prob=rng.random(n),
        st=rng.random(n) * 0.4,
        comp=rng.dirichlet(np.ones(20), size=n),
        cterm_row=crow,
    )
    for model in universe.MODELS:
        u.nterm[model] = nterm
        u.cterm[model] = rng.normal(size=(len(long), dim)) if cterm is None else cterm
    return u


def test_cterm_variant_row_selection():
    # proteins of 1,022 aa (no C-terminal row) and 1,023 aa (C-terminal row 0)
    cterm = np.full((1, 4), 7.0)
    u = toy_universe(2, 4, seed=1, lengths=[1022, 1023], cterm=cterm)
    assert list(u.cterm_row) == [-1, 0]
    for name in ("M8-C", "M35-C"):
        X = universe.features(u, name, [0, 1])
        assert (X[0] == u.nterm[M8][0]).all()
        assert (X[1] == 7.0).all()
    for name in ("M8", "M35"):
        assert (universe.features(u, name, [1]) == u.nterm[M8][1]).all()


def test_cterm_rows_must_match_lengths():
    with pytest.raises(evalio.StopError, match="disagree"):
        universe.check_cterm_rows(np.array([1023]), np.array([-1]))
    with pytest.raises(evalio.StopError, match="disagree"):
        universe.check_cterm_rows(np.array([1022]), np.array([0]))


def test_feature_matrices_of_the_baselines_and_h():
    u = toy_universe(3, 4, seed=2, lengths=[100, 200, 400])
    assert universe.features(u, "B0", [0, 2])[:, 0] == pytest.approx(np.log([100, 400]))
    assert universe.features(u, "B1", [1]).shape == (1, 21)
    h = universe.features(u, "H", [0], "M35")
    assert h.shape == (1, 7) and h[0, 4] == u.sp_prob[0] and h[0, 6] == u.st[0]


def test_youden_threshold_and_platt():
    y = np.array([1, 1, 0, 1, 0, 0], dtype=bool)
    s = np.array([3.0, 2.0, 2.0, 1.0, 0.0, -1.0])
    # thresholds 3: J=1/3; 2: 2/3-1/3=1/3; 1: 1-1/3=2/3; 0: 1-2/3 -> threshold 1.0
    assert models.youden_threshold(y, s) == 1.0
    a, b = models.fit_platt(np.r_[y, y], np.r_[s, s])
    p = models.platt(s, a, b)
    assert a > 0 and (np.diff(p[[5, 4, 3, 0]]) > 0).all()


def test_inner_folds_stop_without_enough_positive_clusters():
    # Review Focus 5: 2 positive clusters cannot fill 3 inner folds
    y = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0], dtype=bool)
    groups = ["a", "a", "b", "c", "d", "e", "f", "g", "h"]
    with pytest.raises(models.ModelError, match="2 positive and 6 negative clusters"):
        models.inner_folds(y, groups, seed=1)


def test_inner_folds_keep_clusters_whole():
    rng = np.random.default_rng(0)
    y = rng.random(90) < 0.3
    groups = np.array([f"c{i // 3}" for i in range(90)])
    for train, test in models.inner_folds(y, groups, seed=5):
        assert not set(groups[train]) & set(groups[test])
        assert y[test].any() and not y[test].all()


def test_planted_signal_is_recovered():
    n = 400
    rng = np.random.default_rng(11)
    y = rng.random(n) < 0.25
    u = toy_universe(n, 10, seed=12, y=y, signal=2.5)
    train, test = np.arange(0, 300), np.arange(300, 400)
    params, scores = models.fit_unit(
        u, train, y[train], [f"c{i}" for i in train], test, 3, ("M8", "H")
    )
    assert roc_auc(np.ones(100), y[test], scores["M8"]["score"])[0] > 0.9
    assert params["H"]["h_variant"] in models.H_VARIANTS
    assert params["M8"]["C"] in models.C_GRID and params["M8"]["n_train"] == 300


def _outer_oof(u, y, groups, candidates, seed):
    oof = np.empty(len(y))
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    for train, test in cv.split(np.zeros(len(y)), y, groups):
        _, s = models.fit_unit(u, train, y[train], groups[train], test, seed, candidates)
        oof[test] = s[candidates[0]]["score"]
    return oof


def test_shuffled_labels_give_chance_auc():
    # Labels drawn at random, independent of the features. 600 proteins x 300 dimensions: the
    # model can memorise its training rows, so a test fold that leaked into training would give
    # an AUC far above 0.5 (0.97 in the prototype when every fold trained on all rows).
    n = 600
    rng = np.random.default_rng(21)
    y = rng.random(n) < 0.3
    u = toy_universe(n, 300, seed=121)
    groups = np.array([f"c{i}" for i in range(n)])
    oof = _outer_oof(u, y, groups, ("M8",), seed=evalio.SEED)
    auc = roc_auc(np.ones(n), y, oof)[0]
    assert 0.45 <= auc <= 0.55, auc


def test_rules_in_fit_unit_use_training_rows_only():
    n = 60
    rng = np.random.default_rng(4)
    y = rng.random(n) < 0.4
    u = toy_universe(n, 3, seed=5, y=y)
    train = np.arange(40)
    params, scores = models.fit_unit(
        u, train, y[train], [f"c{i}" for i in train], np.arange(40, 60), 1, ("R0", "R2")
    )
    assert params["R0"]["g"] is None and params["R2"]["g"] in (
        "highly_probable",
        "probable",
        "weakly",
    )
    assert scores["R2"]["score"] is None and scores["R2"]["call"].dtype == bool
