"""Tests for analysis/ablation_power/model.py and simulate.py."""

import benchmark_data as bd
import model
import numpy as np
import pytest
import simulate
from test_benchmark_data import REAL, make_benchmark

pytest.importorskip("sklearn")


@pytest.fixture
def toy(tmp_path):
    """Two folds with a WEAK composition signal.

    The classes must overlap: if the baseline already reaches a concordance of 1.0 there is no
    headroom, every delta is exactly zero, and a test that a real signal is detected would fail
    for a reason that has nothing to do with the code.
    """
    rng = np.random.default_rng(0)
    proteins, clusters, seqs = [], [], {}
    for clade in ("Alpha", "Beta"):
        for i in range(16):
            tier = "positive" if i % 2 == 0 else "hard_negative"
            acc = f"{clade}{i}"
            # alanine fraction overlaps heavily between the classes
            frac = rng.normal(0.30 if tier == "positive" else 0.22, 0.10)
            frac = float(np.clip(frac, 0.02, 0.6))
            length = int(rng.normal(90, 25))
            n_ala = int(length * frac)
            seq = "A" * n_ala + "".join(rng.choice(list("CDEFGHIKL"), max(length - n_ala, 5)))
            proteins.append(
                {"accession": acc, "tier": tier, "clade_class": clade, "length": len(seq)}
            )
            clusters.append({"accession": acc, "cluster_representative": acc})
            seqs[acc] = seq
    folds = [
        {
            "fold": c,
            "level": "class",
            "members": c,
            "n_positive": 8,
            "n_hard_negative": 8,
            "scored": "yes",
        }
        for c in ("Alpha", "Beta")
    ]
    return bd.load(make_benchmark(tmp_path, proteins, clusters, folds, seqs))


# --------------------------------------------------------------- model
def test_loco_scores_cover_every_scorable_protein(toy):
    scores, _ = model.loco_scores(toy, "composition_length")
    assert np.isfinite(scores).all()
    assert len(scores) == len(toy)


def test_loco_leaves_unscorable_proteins_as_nan(tmp_path):
    proteins, clusters, seqs = [], [], {}
    for clade in ("Alpha", "Beta", "Tiny"):
        for i in range(6):
            tier = "positive" if i % 2 == 0 else "hard_negative"
            acc = f"{clade}{i}"
            proteins.append({"accession": acc, "tier": tier, "clade_class": clade, "length": 50})
            clusters.append({"accession": acc, "cluster_representative": acc})
            seqs[acc] = "A" * 40 if tier == "positive" else "C" * 20
    folds = [
        {
            "fold": c,
            "level": "class",
            "members": c,
            "n_positive": 3,
            "n_hard_negative": 3,
            "scored": "yes" if c != "Tiny" else "no",
        }
        for c in ("Alpha", "Beta", "Tiny")
    ]
    bench = bd.load(make_benchmark(tmp_path, proteins, clusters, folds, seqs))
    scores, _ = model.loco_scores(bench, "composition_length")
    tiny = np.array([f == "Tiny" for f in bench.fold])
    assert np.isnan(scores[tiny]).all()
    assert np.isfinite(scores[~tiny]).all()


def test_a_component_with_the_wrong_number_of_rows_is_refused(toy):
    with pytest.raises(model.ModelError, match="rows"):
        model.loco_scores(toy, "length", component=np.zeros(3))


def test_a_one_dimensional_component_is_accepted_as_a_column(toy):
    scores, _ = model.loco_scores(toy, "length", component=np.zeros(len(toy)))
    assert np.isfinite(scores).all()


def test_choose_bar_returns_a_member_of_the_declared_family(toy):
    name, scores, table, _c = model.choose_bar(toy, tune=False)
    assert name in bd.FEATURE_SETS
    assert set(table) == set(bd.FEATURE_SETS)
    assert np.isfinite(scores).all()


def test_choose_bar_picks_the_best_of_the_family(toy):
    name, _scores, table, _c = model.choose_bar(toy, tune=False)
    best = max(table, key=lambda k: table[k]["stratified_concordance"])
    assert name == best


# --------------------------------------------------------------- planted signals
def test_homogeneous_signal_separates_the_classes():
    y = np.array([0] * 400 + [1] * 400)
    x = simulate.plant_homogeneous(y, 2.0, np.random.default_rng(0))
    assert x[y == 1, 0].mean() - x[y == 0, 0].mean() == pytest.approx(2.0, abs=0.25)


def test_a_zero_effect_signal_is_pure_noise():
    y = np.array([0] * 400 + [1] * 400)
    x = simulate.plant_homogeneous(y, 0.0, np.random.default_rng(1))
    assert x[y == 1, 0].mean() - x[y == 0, 0].mean() == pytest.approx(0.0, abs=0.25)


def test_heterogeneous_signal_is_weaker_outside_the_strong_fold():
    y = np.array([1] * 600)
    fold = np.array(["S"] * 300 + ["W"] * 300)
    x = simulate.plant_heterogeneous(y, fold, 3.0, np.random.default_rng(2), "S", attenuation=1 / 3)
    strong = x[fold == "S", 0].mean()
    weak = x[fold == "W", 0].mean()
    assert strong == pytest.approx(3.0, abs=0.3)
    assert weak == pytest.approx(1.0, abs=0.3)


def test_heterogeneous_needs_a_strong_fold_named():
    with pytest.raises(simulate.SimulationError, match="strong_fold"):
        simulate.power_curve(None, "length", None, "heterogeneous", [1.0], replicates=1)


def test_binary_signal_respects_its_sensitivity_and_fpr():
    y = np.array([1] * 2000 + [0] * 2000)
    x = simulate.plant_binary(y, 0.4, 0.05, np.random.default_rng(3))[:, 0]
    assert x[y == 1].mean() == pytest.approx(0.4, abs=0.04)
    assert x[y == 0].mean() == pytest.approx(0.05, abs=0.02)


def test_an_unknown_signal_shape_is_refused(toy):
    name, scores, _t, _c = model.choose_bar(toy, tune=False)
    with pytest.raises(simulate.SimulationError, match="unknown shape"):
        simulate.power_curve(toy, name, scores, "spiral", [1.0], replicates=1, n_perm=50)


# --------------------------------------------------------------- the decision rule end to end
def test_a_strong_planted_signal_is_detected(toy):
    name, scores, _t, _c = model.choose_bar(toy, tune=False)
    component = simulate.plant_homogeneous(toy.y, 6.0, np.random.default_rng(4))
    delta, p, rejected = simulate.one_replicate(toy, name, scores, component, n_perm=500, seed=0)
    assert delta > 0
    assert rejected, f"a d=6 signal should be detected, got p={p}"


def test_pure_noise_is_usually_not_detected(toy):
    name, scores, _t, _c = model.choose_bar(toy, tune=False)
    rng = np.random.default_rng(5)
    rejects = sum(
        simulate.one_replicate(
            toy, name, scores, rng.standard_normal((len(toy), 1)), n_perm=200, seed=i
        )[2]
        for i in range(20)
    )
    assert rejects <= 5, f"{rejects}/20 noise components rejected; the rule is badly calibrated"


def test_type_i_rate_reports_the_fields_the_report_needs(toy):
    name, scores, _t, _c = model.choose_bar(toy, tune=False)
    got = simulate.type_i_rate(toy, name, scores, dims=1, replicates=5, n_perm=100)
    assert set(got) >= {
        "dims",
        "rejection_rate",
        "mean_delta_under_null",
        "fraction_delta_positive",
    }
    assert 0.0 <= got["rejection_rate"] <= 1.0


# --------------------------------------------------------------- MDE interpolation
def test_minimum_detectable_effect_interpolates_between_grid_points():
    rows = [
        {"effect": 1.0, "power": 0.6, "mean_delta": 0.10},
        {"effect": 2.0, "power": 1.0, "mean_delta": 0.20},
    ]
    effect, delta = simulate.minimum_detectable_effect(rows, 0.8)
    assert effect == pytest.approx(1.5)
    assert delta == pytest.approx(0.15)


def test_minimum_detectable_effect_is_none_when_the_grid_never_reaches_target():
    rows = [{"effect": e, "power": 0.1, "mean_delta": 0.01} for e in (1.0, 2.0, 3.0)]
    assert simulate.minimum_detectable_effect(rows, 0.8) is None


def test_minimum_detectable_effect_takes_the_first_point_if_already_powered():
    rows = [{"effect": 0.5, "power": 0.95, "mean_delta": 0.3}]
    assert simulate.minimum_detectable_effect(rows, 0.8)[0] == pytest.approx(0.5)


# --------------------------------------------------------------- the real benchmark
@pytest.mark.skipif(not REAL.exists(), reason="adhesion_v1 is not present")
def test_on_the_real_benchmark_log_length_alone_beats_the_21_feature_model():
    """Finding 5: the bar must be the best of the family, not composition+length by default."""
    bench = bd.load(REAL)
    _name, _scores, table, _c = model.choose_bar(bench, tune=False)
    assert (
        table["length"]["stratified_concordance"]
        > table["composition_length"]["stratified_concordance"]
    ), table


def test_bar_and_component_are_fitted_at_the_same_regularisation(toy):
    """Regression: the bar was tuned while components ran at a fixed C=1.0.

    Their difference then measured the change in C, not the component, and a pure-noise component
    appeared to cost 0.15 concordance on the real benchmark. choose_bar returns the per-fold C it
    settled on so every later fit can be pinned to it.
    """
    name, scores, _table, c_per_fold = model.choose_bar(toy, tune=True)
    assert c_per_fold, "choose_bar must report the C it chose, per fold"

    # refitting the bar at its own reported C must reproduce the bar exactly
    repeat, used = model.loco_scores(toy, name, c=c_per_fold)
    assert used == c_per_fold
    np.testing.assert_allclose(repeat, scores, atol=1e-12)


def test_a_zero_component_barely_moves_the_score_when_c_is_pinned(toy):
    """A constant column carries no information, so the delta must be ~0, not a C artefact."""
    name, scores, _table, c_per_fold = model.choose_bar(toy, tune=True)
    delta, _p, _rejected = simulate.one_replicate(
        toy, name, scores, np.zeros((len(toy), 1)), c=c_per_fold, n_perm=100, seed=0
    )
    assert abs(delta) < 0.02, f"a constant component moved the score by {delta}"


# --------------------------------------------------------------- parallelism
def test_components_are_reproducible_from_the_task_seed_alone(toy):
    """Per-task seeding is what makes a run independent of execution order and worker count."""
    task = {"shape": "homogeneous", "effect": 1.5, "replicate": 3, "seed": 42}
    first = simulate.build_component(toy, task)
    second = simulate.build_component(toy, task)
    np.testing.assert_array_equal(first, second)


def test_different_replicates_get_different_components(toy):
    a = simulate.build_component(toy, {"shape": "homogeneous", "effect": 1.0, "seed": 1})
    b = simulate.build_component(toy, {"shape": "homogeneous", "effect": 1.0, "seed": 2})
    assert not np.array_equal(a, b)


def test_task_builders_cover_the_whole_grid():
    noise = simulate.noise_tasks([1, 5], replicates=3, seed=0)
    assert len(noise) == 6
    assert {t["dims"] for t in noise} == {1, 5}
    assert len({t["seed"] for t in noise}) == 6, "every task needs its own seed"

    signal = simulate.signal_tasks("binary", [0.1, 0.2], replicates=4, seed=0, fpr=0.05)
    assert len(signal) == 8
    assert len({t["seed"] for t in signal}) == 8


def test_build_component_refuses_an_unknown_shape(toy):
    with pytest.raises(simulate.SimulationError, match="unknown shape"):
        simulate.build_component(toy, {"shape": "spiral", "effect": 1.0, "seed": 0})


def test_build_component_refuses_heterogeneous_without_a_strong_fold(toy):
    with pytest.raises(simulate.SimulationError, match="strong_fold"):
        simulate.build_component(toy, {"shape": "heterogeneous", "effect": 1.0, "seed": 0})


def test_serial_and_parallel_give_identical_results(toy):
    """The whole point: parallelism must change the runtime and nothing else."""
    name, scores, _t, c = model.choose_bar(toy, tune=False)
    tasks = simulate.signal_tasks("homogeneous", [0.5, 2.0], replicates=3, seed=0)
    serial = simulate.run_tasks(toy, name, scores, tasks, bar_c=c, n_perm=100, workers=1)
    parallel = simulate.run_tasks(toy, name, scores, tasks, bar_c=c, n_perm=100, workers=2)
    assert [r["delta"] for r in serial] == [r["delta"] for r in parallel]
    assert [r["rejected"] for r in serial] == [r["rejected"] for r in parallel]


def test_power_curve_is_identical_serial_and_parallel(toy):
    name, scores, _t, c = model.choose_bar(toy, tune=False)
    kw = {"replicates": 3, "c": c, "n_perm": 100, "seed": 0}
    serial = simulate.power_curve(toy, name, scores, "homogeneous", [1.0, 3.0], workers=1, **kw)
    parallel = simulate.power_curve(toy, name, scores, "homogeneous", [1.0, 3.0], workers=2, **kw)
    assert serial == parallel


def test_power_curve_returns_one_row_per_grid_point_in_order(toy):
    name, scores, _t, c = model.choose_bar(toy, tune=False)
    rows = simulate.power_curve(
        toy, name, scores, "homogeneous", [2.0, 0.5], replicates=2, c=c, n_perm=100
    )
    assert [r["effect"] for r in rows] == [0.5, 2.0]
    assert all(r["replicates"] == 2 for r in rows)
