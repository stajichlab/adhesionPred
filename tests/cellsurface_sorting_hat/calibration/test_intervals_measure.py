import json

import numpy as np
import pytest

from cellsurface_sorting_hat.calibration.intervals import cluster_bootstrap, wilson
from cellsurface_sorting_hat.calibration.measure import (
    build_measure,
    make_entry,
    status_from_measure,
    write_status_source,
)
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source


def test_wilson_matches_known_values():
    v, lo, hi = wilson(8, 10)
    assert (
        v == 0.8 and lo == pytest.approx(0.4902, abs=1e-3) and hi == pytest.approx(0.9433, abs=1e-3)
    )
    assert wilson(0, 10)[1] == 0.0 and wilson(10, 10)[2] == pytest.approx(1.0)
    assert wilson(0, 0) is None


def test_bootstrap_point_estimates_and_ordering():
    y = np.array([1] * 30 + [0] * 70)
    call = np.array([True] * 24 + [False] * 6 + [True] * 7 + [False] * 63)
    clusters = [f"c{i}" for i in range(100)]
    r = cluster_bootstrap(y, call, clusters, n_boot=500, seed=3)
    assert r["sensitivity"]["value"] == pytest.approx(0.8) and r["specificity"][
        "value"
    ] == pytest.approx(0.9)
    for k in (r["sensitivity"], r["specificity"]):
        assert 0 <= k["lo"] <= k["value"] <= k["hi"] <= 1
    assert r == cluster_bootstrap(y, call, clusters, n_boot=500, seed=3)  # same seed, same result


def test_clusters_widen_the_interval_compared_with_single_proteins():
    y = np.array([1] * 40 + [0] * 40)
    call = np.array(([True] * 20 + [False] * 20) * 2)
    single = cluster_bootstrap(y, call, [str(i) for i in range(80)], n_boot=800, seed=1)[
        "sensitivity"
    ]
    # 8 clusters of 5 positives that are all called or all missed: far less independent information
    pos_clusters = [f"p{i // 5}" for i in range(40)]
    call2 = np.array([(i // 5) % 2 == 0 for i in range(40)] + [False] * 40)
    grouped = cluster_bootstrap(
        y, call2, pos_clusters + [f"n{i}" for i in range(40)], n_boot=800, seed=1
    )["sensitivity"]
    assert (grouped["hi"] - grouped["lo"]) > (single["hi"] - single["lo"])


def test_bootstrap_without_negatives_has_no_specificity():
    r = cluster_bootstrap([1, 1, 1], [True, False, True], ["a", "b", "c"], n_boot=100)
    assert r["specificity"] is None and r["sensitivity"] is not None


@pytest.mark.parametrize(
    "args", [([1, 0], [True], ["a", "b"]), ([2, 0], [True, False], ["a", "b"])]
)
def test_bootstrap_refuses_bad_input(args):
    with pytest.raises(ValueError):
        cluster_bootstrap(*args)


def test_status_rule_follows_phase_c_and_needs_negatives_for_an_estimate():
    def m(n_pos, sens, n_neg=None, spec=None):
        out = {"n_pos": n_pos, "sensitivity": sens}
        if spec:
            out.update(n_neg=n_neg, specificity=spec)
        return out

    tight = {"value": 0.6, "lo": 0.55, "hi": 0.65}
    wide = {"value": 0.8, "lo": 0.5, "hi": 0.95}
    sp_ok = {"value": 0.96, "lo": 0.95, "hi": 0.97}
    sp_wide = {"value": 0.9, "lo": 0.6, "hi": 0.99}
    assert status_from_measure(m(232, tight, 4244, sp_ok)) == "estimated"
    assert status_from_measure(m(25, wide, 100, sp_ok)) == "smoke"  # sensitivity half-width 0.225
    assert status_from_measure(m(16, tight, 100, sp_ok)) == "smoke"  # fewer than 20 positives
    assert status_from_measure(m(232, tight, 15, sp_ok)) == "smoke"  # fewer than 20 negatives
    assert status_from_measure(m(232, tight, 4244, sp_wide)) == "smoke"  # specificity too wide
    assert status_from_measure(m(232, tight)) == "smoke"  # never tested on negatives
    # negatives were counted but no specificity was measured: still not an estimate
    assert status_from_measure({"n_pos": 232, "n_neg": 500, "sensitivity": tight}) == "smoke"
    assert status_from_measure({"n_pos": 0}) == "unvalidated"


def test_a_leakage_cap_limits_an_estimate_to_smoke():
    measure = {
        "n_pos": 232,
        "n_neg": 4244,
        "sensitivity": {"value": 0.6, "lo": 0.55, "hi": 0.65},
        "specificity": {"value": 0.96, "lo": 0.95, "hi": 0.97},
    }
    assert make_entry([1], measure)["status"] == "estimated"
    assert make_entry([1], measure, cap="smoke")["status"] == "smoke"


def test_build_measure_counts_and_notes():
    m = build_measure(
        "set",
        "truth.tsv",
        [1, 1, 0, 0],
        [True, False, False, False],
        ["a", "b", "c", "d"],
        notes="n",
        n_boot=50,
    )
    assert (m["n_pos"], m["n_neg"], m["calibration_set"], m["notes"]) == (2, 2, "set", "n")


def test_write_status_source_uses_the_module_identity_and_validates(tmp_path):
    write_module(
        tmp_path, ModuleSpec("pfam_adhesion", "1", {"a": 1}), [], [{"id": "A", "state": "ok"}]
    )
    m = build_measure("set", "t", [1, 0], [True, False], ["a", "b"], n_boot=20)
    path = write_status_source(tmp_path, "pfam_adhesion", [make_entry([40], m)])
    rec = load_status_source(path)
    run = json.loads((tmp_path / "modules" / "pfam_adhesion.json").read_text())
    assert rec.identity.params_hash == run["params_hash"] and rec.entries[0].taxa == (40,)
    with pytest.raises(ValueError, match="appears in two entries"):
        write_status_source(
            tmp_path, "pfam_adhesion", [make_entry([40], m), make_entry([40, 41], m)]
        )


def test_all_positives_called_gives_a_width_not_a_zero_width_interval():
    y = [1] * 20 + [0] * 20
    r = cluster_bootstrap(y, [True] * 20 + [False] * 20, [f"c{i}" for i in range(40)], n_boot=300)
    assert (
        r["sensitivity"]["value"] == 1.0 and r["sensitivity"]["lo"] < 0.9
    )  # Wilson on 20 clusters
    assert r["specificity"]["value"] == 1.0 and r["specificity"]["lo"] < 0.9
    assert (r["n_clusters_pos"], r["n_clusters_neg"]) == (20, 20)


def test_one_cluster_per_class_cannot_give_an_estimate():
    y = [1] * 30 + [0] * 30
    call = [True] * 24 + [False] * 6 + [False] * 28 + [True] * 2
    clusters = ["pos"] * 30 + ["neg"] * 30
    m = build_measure("s", "t", y, call, clusters, n_boot=200)
    assert (m["n_clusters_pos"], m["n_clusters_neg"]) == (1, 1)
    assert status_from_measure(m) == "smoke"  # many proteins, almost no independent clusters
    assert (m["sensitivity"]["hi"] - m["sensitivity"]["lo"]) / 2 > 0.10


def test_wilson_refuses_k_above_n():
    with pytest.raises(ValueError):
        wilson(5, 3)


def test_a_tight_interval_with_few_recorded_clusters_is_not_an_estimate():
    tight = {"value": 0.6, "lo": 0.55, "hi": 0.65}
    base = {
        "n_pos": 300,
        "n_neg": 3000,
        "sensitivity": tight,
        "specificity": {"value": 0.96, "lo": 0.95, "hi": 0.97},
    }
    assert status_from_measure({**base, "n_clusters_pos": 40, "n_clusters_neg": 300}) == "estimated"
    assert status_from_measure({**base, "n_clusters_pos": 5, "n_clusters_neg": 300}) == "smoke"
    assert status_from_measure({**base, "n_clusters_pos": 40, "n_clusters_neg": 3}) == "smoke"


# ---- hardening tests added for review edge cases -------------------------------------------


def _measure(
    n_pos=100,
    n_neg=100,
    sens=(0.5, 0.6, 0.7),
    spec=(0.9, 0.95, 0.99),
    cp=None,
    cn=None,
):
    out = {
        "n_pos": n_pos,
        "n_neg": n_neg,
        "sensitivity": dict(zip(("lo", "value", "hi"), sens, strict=True)),
        "specificity": dict(zip(("lo", "value", "hi"), spec, strict=True)),
    }
    if cp is not None:
        out["n_clusters_pos"] = cp
    if cn is not None:
        out["n_clusters_neg"] = cn
    return out


@pytest.mark.parametrize(
    "n_pos,n_neg,expected",
    [(20, 20, "estimated"), (19, 20, "smoke"), (20, 19, "smoke"), (21, 21, "estimated")],
)
def test_count_boundaries_19_vs_20(n_pos, n_neg, expected):
    assert status_from_measure(_measure(n_pos, n_neg)) == expected


@pytest.mark.parametrize(
    "cp,cn,expected",
    [(20, 20, "estimated"), (19, 20, "smoke"), (20, 19, "smoke"), (0, 20, "smoke")],
)
def test_cluster_boundaries_19_vs_20(cp, cn, expected):
    assert status_from_measure(_measure(cp=cp, cn=cn)) == expected


@pytest.mark.parametrize(
    "lo,hi,expected",
    [
        (0.5, 0.7, "estimated"),  # half-width 0.10 exactly: allowed (<=)
        (0.6, 0.8, "estimated"),  # 0.8 - 0.6 is 0.20000000000000007 in floating point
        (0.3, 0.5, "estimated"),
        (0.5, 0.7001, "smoke"),  # 0.10005
        (0.4999, 0.7, "smoke"),
    ],
)
def test_half_width_boundary_for_both_rates(lo, hi, expected):
    value = (lo + hi) / 2
    assert status_from_measure(_measure(sens=(lo, value, hi))) == expected
    assert status_from_measure(_measure(spec=(lo, value, hi))) == expected


def test_status_thresholds_can_be_changed_by_keyword_only_for_callers_that_ask():
    m = _measure(n_pos=19, n_neg=19)
    assert status_from_measure(m) == "smoke"
    assert status_from_measure(m, min_positives=19, min_negatives=19) == "estimated"


def test_unvalidated_cases():
    assert status_from_measure({}) == "unvalidated"
    assert status_from_measure({"n_pos": 5}) == "unvalidated"  # no sensitivity
    m = _measure()
    m["n_pos"] = 0
    assert status_from_measure(m) == "unvalidated"  # sensitivity without any positive


def test_endpoints_at_zero_and_one_are_valid_rates():
    m = _measure(sens=(0.9, 1.0, 1.0), spec=(0.95, 1.0, 1.0))
    assert status_from_measure(m) == "estimated"
    m = _measure(sens=(0.0, 0.0, 0.05), spec=(0.0, 0.0, 0.05))
    assert status_from_measure(m) == "estimated"  # tight but poor: width rule only


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -0.1, 1.1])
def test_status_refuses_bad_rate_values(bad):
    m = _measure()
    m["sensitivity"]["hi"] = bad
    with pytest.raises(ValueError, match="sensitivity"):
        status_from_measure(m)


@pytest.mark.parametrize("bad", [-1, 2.5, float("nan"), True, "3"])
def test_status_refuses_bad_counts(bad):
    m = _measure()
    m["n_pos"] = bad
    with pytest.raises(ValueError, match="n_pos"):
        status_from_measure(m)


def test_status_refuses_unordered_interval():
    m = _measure(sens=(0.7, 0.6, 0.5))
    with pytest.raises(ValueError, match="sensitivity"):
        status_from_measure(m)


def test_cap_takes_the_weaker_label_and_never_strengthens():
    est = _measure()
    smoke = _measure(n_pos=5)
    assert make_entry([1], est, cap="estimated")["status"] == "estimated"
    assert make_entry([1], est, cap="smoke")["status"] == "smoke"
    assert make_entry([1], est, cap="unvalidated")["status"] == "unvalidated"
    assert make_entry([1], smoke, cap="estimated")["status"] == "smoke"  # cap cannot raise
    assert make_entry([1], smoke, cap="smoke")["status"] == "smoke"
    assert make_entry([1], smoke, cap="unvalidated")["status"] == "unvalidated"
    with pytest.raises(ValueError):
        make_entry([1], est, cap="strong")


def test_never_tested_on_negatives_is_at_most_smoke_even_with_huge_counts():
    m = {
        "n_pos": 10_000,
        "n_neg": 10_000,
        "sensitivity": {"value": 0.6, "lo": 0.59, "hi": 0.61},
    }
    assert status_from_measure(m) == "smoke"


def test_make_entry_refuses_bad_taxa():
    with pytest.raises(ValueError):
        make_entry([True], _measure())
    with pytest.raises(ValueError):
        make_entry(["x"], _measure())
    with pytest.raises(ValueError):
        make_entry([-5], _measure())


@pytest.mark.parametrize("k,n", [(-1, 5), (6, 5), (1.5, 5), (float("nan"), 5)])
def test_wilson_refuses_bad_k(k, n):
    with pytest.raises(ValueError):
        wilson(k, n)


def test_wilson_refuses_non_finite_or_negative_n_and_z():
    with pytest.raises(ValueError):
        wilson(0, float("nan"))
    with pytest.raises(ValueError):
        wilson(0, -1)
    with pytest.raises(ValueError):
        wilson(1, 2, z=float("nan"))


def test_wilson_endpoints_and_symmetry():
    v, lo, hi = wilson(0, 1)
    assert (v, lo) == (0.0, 0.0) and 0 < hi < 1
    v, lo, hi = wilson(1, 1)
    assert v == 1.0 and hi == pytest.approx(1.0) and 0 < lo < 1
    _, lo1, hi1 = wilson(3, 10)
    _, lo2, hi2 = wilson(7, 10)
    assert lo1 == pytest.approx(1 - hi2) and hi1 == pytest.approx(1 - lo2)


def test_bootstrap_zero_positives_and_zero_negatives():
    r = cluster_bootstrap([0, 0, 0], [True, False, False], ["a", "b", "c"], n_boot=50)
    assert r["sensitivity"] is None and r["specificity"] is not None
    assert (r["n_clusters_pos"], r["n_clusters_neg"]) == (0, 3)
    m = build_measure("s", "t", [0, 0], [False, False], ["a", "b"], n_boot=20)
    assert "sensitivity" not in m and m["n_pos"] == 0
    assert status_from_measure(m) == "unvalidated"


def test_bootstrap_refuses_empty_input():
    with pytest.raises(ValueError, match="empty"):
        cluster_bootstrap([], [], [])


def test_bootstrap_single_item_each_class_is_wide():
    r = cluster_bootstrap([1, 0], [True, False], ["a", "b"], n_boot=50)
    for k in ("sensitivity", "specificity"):
        assert r[k]["value"] == 1.0 and r[k]["lo"] < 0.5
        assert 0 <= r[k]["lo"] <= r[k]["value"] <= r[k]["hi"] <= 1


def test_all_same_labels_and_all_called_or_none_called():
    n = 40
    clusters = [f"c{i}" for i in range(n)]
    r = cluster_bootstrap([1] * n, [False] * n, clusters, n_boot=100)
    assert r["sensitivity"]["value"] == 0.0 and r["sensitivity"]["lo"] == 0.0
    assert 0 < r["sensitivity"]["hi"] < 0.2  # Wilson on 40 clusters, not a zero-width interval
    r = cluster_bootstrap([0] * n, [True] * n, clusters, n_boot=100)
    assert r["specificity"]["value"] == 0.0 and r["specificity"]["hi"] > 0


def test_interval_endpoints_stay_within_zero_and_one():
    y = [1] * 25 + [0] * 25
    call = [True] * 25 + [False] * 25
    r = cluster_bootstrap(y, call, [f"c{i}" for i in range(50)], n_boot=100)
    for k in (r["sensitivity"], r["specificity"]):
        assert 0.0 <= k["lo"] <= k["value"] <= k["hi"] <= 1.0
    # the result must pass the repository validator
    m = build_measure("s", "t", y, call, [f"c{i}" for i in range(50)], n_boot=100)
    assert status_from_measure(m) == "estimated"


def test_seed_determinism_and_dependence():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 2, 200)
    call = rng.random(200) < 0.5
    clusters = [f"c{i % 37}" for i in range(200)]
    a = cluster_bootstrap(y, call, clusters, n_boot=300, seed=5)
    assert a == cluster_bootstrap(y, call, clusters, n_boot=300, seed=5)
    b = cluster_bootstrap(y, call, clusters, n_boot=300, seed=6)
    # the point estimates do not depend on the seed; the resampled intervals may
    assert a["sensitivity"]["value"] == b["sensitivity"]["value"]
    # default seed is fixed, so two default calls agree
    assert cluster_bootstrap(y, call, clusters, n_boot=50) == cluster_bootstrap(
        y, call, clusters, n_boot=50
    )


def test_default_resample_count_is_2000_and_seed_is_1(monkeypatch):
    import inspect

    sig = inspect.signature(cluster_bootstrap)
    assert sig.parameters["n_boot"].default == 2000 and sig.parameters["seed"].default == 1
    sig = inspect.signature(build_measure)
    assert sig.parameters["n_boot"].default == 2000 and sig.parameters["seed"].default == 1


def test_interval_is_the_widest_of_bootstrap_and_wilson():
    # 5 positive clusters of 8 proteins each: Wilson on 5 clusters is wide
    y = [1] * 40 + [0] * 40
    call = [True] * 40 + [False] * 40
    clusters = [f"p{i // 8}" for i in range(40)] + [f"n{i}" for i in range(40)]
    r = cluster_bootstrap(y, call, clusters, n_boot=200)
    _, wlo, whi = wilson(5, 5)
    assert r["sensitivity"]["lo"] <= wlo + 1e-12  # at least as wide as Wilson(5/5)
    assert r["sensitivity"]["hi"] >= whi - 1e-12
    assert r["sensitivity"]["lo"] < 0.6


@pytest.mark.parametrize(
    "y,call",
    [
        ([1, 0.5], [True, False]),
        ([1, float("nan")], [True, False]),
        ([1, -1], [True, False]),
        ([1, 0], [1, 2]),
        ([1, 0], [True, float("nan")]),
        ([1, 0], ["yes", "no"]),
    ],
)
def test_bootstrap_refuses_non_binary_or_non_finite_labels_and_calls(y, call):
    with pytest.raises(ValueError):
        cluster_bootstrap(y, call, ["a", "b"], n_boot=10)


@pytest.mark.parametrize("n_boot", [0, -3, 1.5, True])
def test_bootstrap_refuses_bad_n_boot(n_boot):
    with pytest.raises(ValueError, match="n_boot"):
        cluster_bootstrap([1, 0], [True, False], ["a", "b"], n_boot=n_boot)


def test_bootstrap_refuses_missing_cluster_labels():
    with pytest.raises(ValueError, match="cluster"):
        cluster_bootstrap([1, 0], [True, False], ["a", None])
    with pytest.raises(ValueError, match="cluster"):
        cluster_bootstrap([1, 0], [True, False], ["a", ""])


def test_bootstrap_accepts_call_as_0_1_ints_and_numpy_bool():
    a = cluster_bootstrap([1, 0], [1, 0], ["a", "b"], n_boot=10)
    b = cluster_bootstrap([1, 0], np.array([True, False]), ["a", "b"], n_boot=10)
    assert a == b


def test_write_status_source_validates_before_writing(tmp_path):
    write_module(
        tmp_path, ModuleSpec("pfam_adhesion", "1", {"a": 1}), [], [{"id": "A", "state": "ok"}]
    )
    m = build_measure("set", "t", [1, 0], [True, False], ["a", "b"], n_boot=20)
    good = make_entry([40], m)
    target = tmp_path / "status" / "pfam_adhesion.json"
    bad_status = {**good, "taxa": [41], "status": "excellent"}
    overclaim = {**good, "taxa": [42], "status": "estimated"}  # 1 positive cannot be estimated
    bad_measure = {**good, "taxa": [43], "measure": {**m, "n_pos": -1}}
    for entries in ([bad_status], [overclaim], [bad_measure], [good, {**good, "taxa": [40]}]):
        with pytest.raises(ValueError):
            write_status_source(tmp_path, "pfam_adhesion", entries)
        assert not target.exists()
    with pytest.raises(FileNotFoundError):
        write_status_source(tmp_path, "no_such_module", [good])
    assert not (tmp_path / "status" / "no_such_module.json").exists()


def test_write_status_source_refuses_nan_in_a_measure(tmp_path):
    write_module(
        tmp_path, ModuleSpec("pfam_adhesion", "1", {"a": 1}), [], [{"id": "A", "state": "ok"}]
    )
    m = build_measure("set", "t", [1, 0], [True, False], ["a", "b"], n_boot=20)
    m["notes"] = "x"
    entry = make_entry([40], m)
    entry["measure"]["sensitivity"]["value"] = float("nan")
    with pytest.raises(ValueError):
        write_status_source(tmp_path, "pfam_adhesion", [entry])
    assert not (tmp_path / "status" / "pfam_adhesion.json").exists()


def test_write_status_source_keeps_a_weaker_status_than_the_measure_allows(tmp_path):
    write_module(
        tmp_path, ModuleSpec("pfam_adhesion", "1", {"a": 1}), [], [{"id": "A", "state": "ok"}]
    )
    y = [1] * 25 + [0] * 25
    call = [True] * 25 + [False] * 25
    m = build_measure("set", "t", y, call, [f"c{i}" for i in range(50)], n_boot=50)
    entry = make_entry([40], m, cap="smoke")  # leakage cap
    path = write_status_source(tmp_path, "pfam_adhesion", [entry])
    assert load_status_source(path).entries[0].status == "smoke"
