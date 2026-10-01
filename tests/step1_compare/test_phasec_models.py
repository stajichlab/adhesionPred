"""models.py and universe.py: rows, nested protocol, learnability, chance level (spec 3.3, 6)."""

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import evalio  # noqa: E402
import models  # noqa: E402
import rules  # noqa: E402
import universe  # noqa: E402
from metrics import pr_auc, roc_auc  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
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


def test_rule_settings_come_from_training_rows_only():
    n = 60
    rng = np.random.default_rng(4)
    y = rng.random(n) < 0.4
    u = toy_universe(n, 3, seed=5, y=y)
    train, test = np.arange(40), np.arange(40, 60)
    groups = [f"c{i}" for i in train]
    params, scores = models.fit_unit(u, train, y[train], groups, test, 1, ("R0", "R2"))
    assert params["R0"]["g"] is None and params["R2"]["g"] in (
        "highly_probable",
        "probable",
        "weakly",
    )
    assert scores["R2"]["score"] is None and scores["R2"]["call"].dtype == bool
    # g and t equal fit_rule on the training rows
    want = rules.fit_rule("R2", y[train], u.sp[train], u.rank[train], u.st[train])
    assert (params["R2"]["g"], params["R2"]["t"]) == (want["g"], want["t"])
    # test rows are not an input: changing their features or labels changes nothing
    u2 = toy_universe(n, 3, seed=5, y=y)
    r2 = np.random.default_rng(99)
    u2.sp[test] = r2.random(20) < 0.5
    u2.st[test] = r2.random(20) * 0.4
    u2.rank[test] = r2.integers(0, 4, size=20).astype(np.int8)
    y_flipped = y.copy()
    y_flipped[test] = ~y_flipped[test]
    p2, _ = models.fit_unit(u2, train, y_flipped[train], groups, test, 1, ("R0", "R2"))
    assert p2 == params


# --- oracle, ties, in-sample versus inner out-of-fold, scaling, boundaries -------------------


def _setup(signal, seed=31, dim=150, n=300, n_train=240):
    rng = np.random.default_rng(seed)
    y = rng.random(n) < 0.3
    u = toy_universe(n, dim, seed=seed + 1, y=y, signal=signal)
    train, test = np.arange(n_train), np.arange(n_train, n)
    groups = np.array([f"c{i // 2}" for i in train])
    return u, y, train, test, groups


def distinct_universe(n, dim, seed, y, signal):
    """Universe whose four embedding windows differ; half of the proteins are over 1,022 aa."""
    y = np.asarray(y, dtype=bool)
    lengths = np.where(np.arange(n) % 2 == 0, 1100, 300)
    u = toy_universe(n, dim, seed=seed, y=y, lengths=lengths)
    rng = np.random.default_rng(seed + 1000)
    long = np.nonzero(lengths > 1022)[0]
    for model in universe.MODELS:
        u.nterm[model] = rng.normal(size=(n, dim))
        u.nterm[model][:, 0] += signal * y
        u.cterm[model] = rng.normal(size=(len(long), dim))
        u.cterm[model][:, 0] += signal * y[long]
    return u


def _ap(y, score):
    return float(pr_auc(np.ones(len(y)), y, score)[0])


def _oracle_choice(u, cand, train, y, folds, variants=(None,)):
    """Strict > over the ascending grid (variants in order): the first best setting wins."""
    best = None
    for v in variants:
        X = universe.features(u, cand, train, v)
        for C in models.C_GRID:
            ap = _ap(y, models.oof_scores(X, y, folds, C, {}))
            if best is None or ap > best[0]:
                best = (ap, v, C)
    return best[1], best[2]


def test_fit_unit_equals_an_independent_recomputation():
    u, y, train, test, groups = _setup(1.5)
    # a mean shift in the test rows: a scaler that saw them would change every score
    for model in universe.MODELS:
        u.nterm[model][test, 1] += 40.0
    yt = y[train]
    params, scores = models.fit_unit(u, train, yt, groups, test, 3, ("M8", "H"))
    folds = models.inner_folds(yt, groups, 3)
    _, C = _oracle_choice(u, "M8", train, yt, folds)
    p = params["M8"]
    assert p["C"] == C
    X = universe.features(u, "M8", train)
    oof = models.oof_scores(X, yt, folds, C, {})
    assert p["threshold"] == models.youden_threshold(yt, oof)
    # unweighted Platt fit with an intercept: the mean probability equals the prevalence
    prob = models.platt(oof, p["platt_a"], p["platt_b"])
    assert prob.mean() == pytest.approx(yt.mean(), abs=1e-3)
    ref = models.make_lr(C).fit(X, yt).decision_function(universe.features(u, "M8", test))
    assert np.allclose(scores["M8"]["score"], ref)
    assert (scores["M8"]["call"] == (ref >= p["threshold"])).all()
    assert np.allclose(scores["M8"]["prob"], models.platt(ref, p["platt_a"], p["platt_b"]))
    variant, C = _oracle_choice(u, "H", train, yt, folds, models.H_VARIANTS)
    assert (params["H"]["h_variant"], params["H"]["C"]) == (variant, C)
    Xh = universe.features(u, "H", train, variant)
    refh = models.make_lr(C).fit(Xh, yt).decision_function(universe.features(u, "H", test, variant))
    assert np.allclose(scores["H"]["score"], refh)


def test_inner_scaling_is_per_fold():
    # Oracle with an explicit per-fold scaler: it must equal oof_scores on shifted features.
    u, y, train, _, groups = _setup(1.0, seed=7, dim=20)
    X = universe.features(u, "M8", train).copy()
    X[:, 2] += np.where(np.arange(len(train)) < 120, 0.0, 500.0)
    yt = y[train]
    folds = models.inner_folds(yt, groups, 3)
    got = models.oof_scores(X, yt, folds, 0.1, {})
    want = np.empty(len(yt))
    for tr, te in folds:
        mu, sd = X[tr].mean(axis=0), X[tr].std(axis=0)
        lr = LogisticRegression(C=0.1, class_weight="balanced", max_iter=5000)
        lr.fit((X[tr] - mu) / sd, yt[tr])
        want[te] = lr.decision_function((X[te] - mu) / sd)
    assert np.allclose(got, want, atol=1e-6)


def test_ties_go_to_the_smaller_c_and_the_first_variant():
    # every window separates the classes perfectly: all inner PR-AUC values are 1.0
    n = 300
    y = np.random.default_rng(31).random(n) < 0.3
    u = distinct_universe(n, 5, seed=32, y=y, signal=30.0)
    train, test = np.arange(240), np.arange(240, 300)
    groups = np.array([f"c{i // 2}" for i in train])
    params, _ = models.fit_unit(u, train, y[train], groups, test, 3, ("M8", "H"))
    assert set(params["M8"]["inner_pr_auc"].values()) == {1.0}
    assert params["M8"]["C"] == 0.001
    assert set(params["H"]["inner_pr_auc"].values()) == {1.0}
    assert (params["H"]["h_variant"], params["H"]["C"]) == ("M8", 0.001)


def test_youden_tie_goes_to_the_higher_threshold():
    y = np.array([1, 0, 1, 0], dtype=bool)
    s = np.array([4.0, 3.0, 2.0, 1.0])  # J = 0.5 at threshold 4 and at threshold 2
    assert models.youden_threshold(y, s) == 4.0


def test_choice_follows_inner_out_of_fold_not_in_sample_scores():
    # 80 rows, 100 noise dimensions: the training rows are fitted almost perfectly in sample,
    # so the in-sample best setting differs from the inner out-of-fold best (checked below).
    n = 80
    y = np.random.default_rng(8).random(n) < 0.3
    u = distinct_universe(n, 100, seed=58, y=y, signal=0.0)
    train, groups = np.arange(n), np.array([f"c{i}" for i in range(n)])
    params, _ = models.fit_unit(u, train, y, groups, np.arange(10), 5, ("M8", "H"))
    folds = models.inner_folds(y, groups, 5)
    _, C = _oracle_choice(u, "M8", train, y, folds)
    variant, Ch = _oracle_choice(u, "H", train, y, folds, models.H_VARIANTS)
    assert params["M8"]["C"] == C
    assert (params["H"]["h_variant"], params["H"]["C"]) == (variant, Ch)

    def in_sample(cand, v, c):
        X = universe.features(u, cand, train, v)
        return _ap(y, models.make_lr(c).fit(X, y).decision_function(X))

    m8 = {c: in_sample("M8", None, c) for c in models.C_GRID}
    h = {(v, c): in_sample("H", v, c) for v in models.H_VARIANTS for c in models.C_GRID}
    best_m8 = next(c for c in models.C_GRID if m8[c] == max(m8.values()))
    best_h = next(k for k in h if h[k] == max(h.values()))
    assert best_m8 != C and best_h[0] != variant  # the fixture is informative


def test_cterm_rows_with_several_long_proteins():
    cterm = np.array([[10.0] * 4, [20.0] * 4, [30.0] * 4])
    lengths = [1022, 1500, 1023, 2000, 50]
    u = toy_universe(5, 4, seed=3, lengths=lengths, cterm=cterm)
    assert list(u.cterm_row) == [-1, 0, 1, 2, -1]
    X = universe.features(u, "M8-C", [3, 0, 2, 1, 4])
    assert (X[0] == 30.0).all() and (X[1] == u.nterm[M8][0]).all()
    assert (X[2] == 20.0).all() and (X[3] == 10.0).all() and (X[4] == u.nterm[M8][4]).all()


def test_call_is_greater_or_equal_at_the_threshold(monkeypatch):
    u, y, train, test, groups = _setup(1.5)
    _, first = models.fit_unit(u, train, y[train], groups, test, 3, ("M8",))
    s = first["M8"]["score"]
    target = float(np.sort(s)[len(s) // 2])
    monkeypatch.setattr(models, "youden_threshold", lambda yy, ss: target)
    params, second = models.fit_unit(u, train, y[train], groups, test, 3, ("M8",))
    assert params["M8"]["threshold"] == target
    assert (second["M8"]["call"] == (second["M8"]["score"] >= target)).all()
    assert second["M8"]["call"][np.argmax(second["M8"]["score"] == target)]


def test_composition_divides_by_the_full_length():
    c = universe.composition("AAXX")
    assert c[0] == 0.5 and sum(c) == pytest.approx(0.5)
    assert "full length" in universe.composition.__doc__
    with pytest.raises(models.ModelError, match="empty"):
        universe.composition("")


def test_fit_unit_stops_on_unequal_lengths():
    u = toy_universe(20, 3, seed=1)
    y = np.arange(10) % 2 == 0
    with pytest.raises(models.ModelError, match="must agree"):
        models.fit_unit(u, np.arange(10), y, ["a"] * 9, np.arange(10, 20), 1, ("B0",))
    with pytest.raises(models.ModelError, match="must agree"):
        models.fit_unit(u, np.arange(10), y[:8], ["a"] * 10, np.arange(10, 20), 1, ("B0",))


def test_universe_rejects_duplicate_hashes():
    u = toy_universe(3, 2, seed=1)
    with pytest.raises(models.ModelError, match="more than once"):
        universe.Universe(
            hashes=["a", "a", "b"],
            length=u.length,
            sp=u.sp,
            sp_prob=u.sp_prob,
            rank=u.rank,
            gpi_prob=u.gpi_prob,
            st=u.st,
            comp=u.comp,
            cterm_row=u.cterm_row,
        )


def test_h_needs_a_variant():
    u = toy_universe(3, 2, seed=1)
    with pytest.raises(models.ModelError, match="h_variant"):
        universe.features(u, "H", [0])
    with pytest.raises(models.ModelError, match="h_variant"):
        universe.features(u, "H", [0], "B0")
