"""rules.py: R0-R2 truth table, cut-off values, Youden fit (spec 3.3, section 6)."""

import itertools
import math

import pytest

np = pytest.importorskip("numpy")

import rules  # noqa: E402

GPI = ("highly_probable", "probable", "weakly", "none", "too_short")


def test_rule_truth_table():
    # every SP x GPI class x Ser+Thr combination, with Ser+Thr exactly at the cut-off t = 0.25
    st_values = (0.2499, 0.25, 0.30)
    for sp, gpi, st in itertools.product((True, False), GPI, st_values):
        rank = rules.gpi_rank([gpi])
        for g in rules.G_VALUES:
            gpi_ok = rules.GPI_RANK[gpi] >= rules.GPI_RANK[g]
            assert rules.rule_call("R0", [sp], rank, [st])[0] == sp
            assert rules.rule_call("R1", [sp], rank, [st], g)[0] == (sp and gpi_ok)
            want = sp and (gpi_ok or st >= 0.25)
            assert rules.rule_call("R2", [sp], rank, [st], g, 0.25)[0] == want, (sp, gpi, st, g)


def test_cut_off_values_parse_equal_to_the_grid():
    # ser_thr_frac is stored with 6 decimals; "0.300000" must reach t = 0.30
    for t in rules.T_VALUES:
        st = float(f"{t:.6f}")
        assert rules.rule_call("R2", [True], rules.gpi_rank(["none"]), [st], "weakly", t)[0]
        below = float(f"{t - 0.000001:.6f}")
        assert not rules.rule_call("R2", [True], rules.gpi_rank(["none"]), [below], "weakly", t)[0]
    assert rules.T_VALUES[2] == 0.30 and 0.2 + 0.05 * 2 != 0.30  # why the grid is literal


def test_t_grid_is_0_20_to_0_40():
    # owner decision 2026-10-01: the grid starts at 0.20 (at 0.10, R2 equals R0 on the real data)
    assert rules.T_VALUES == (0.20, 0.25, 0.30, 0.35, 0.40)
    # positives at Ser+Thr 0.15, negatives at 0.05: a grid with 0.10 or 0.15 would give J = 1
    y = np.array([1] * 4 + [0] * 4, dtype=bool)
    rank = rules.gpi_rank(["none"] * 8)
    st = np.array([0.15] * 4 + [0.05] * 4)
    fit = rules.fit_rule("R2", y, np.ones(8, dtype=bool), rank, st)
    assert fit == {"g": "highly_probable", "t": 0.40, "j": 0.0}  # no t of the grid calls them


def test_too_short_is_never_gpi():
    rank = rules.gpi_rank(["too_short"])
    assert not rules.rule_call("R1", [True], rank, [0.0], "weakly")[0]


def test_unknown_gpi_call_is_refused():
    with pytest.raises(ValueError, match="unknown gpi_call"):
        rules.gpi_rank(["maybe"])


def test_youden():
    y = [1, 1, 0, 0]
    assert rules.youden(y, [1, 0, 1, 0]) == 0.0
    assert rules.youden(y, [1, 1, 0, 0]) == 1.0
    assert math.isnan(rules.youden([0, 0], [1, 0]))


def test_fit_rule_finds_the_planted_cut_and_breaks_ties_to_the_stricter_rule():
    # positives: SP, no GPI, Ser+Thr 0.32; negatives: SP, no GPI, Ser+Thr 0.20
    y = np.array([1] * 5 + [0] * 5, dtype=bool)
    sp = np.ones(10, dtype=bool)
    rank = rules.gpi_rank(["none"] * 10)
    st = np.array([0.32] * 5 + [0.20] * 5)
    fit = rules.fit_rule("R2", y, sp, rank, st)
    # t = 0.25 and t = 0.30 both give J = 1; the grid runs from 0.40 down, so 0.30 wins
    assert fit == {"g": "highly_probable", "t": 0.30, "j": 1.0}
    assert rules.fit_rule("R0", y, sp, rank, st)["j"] == 0.0


def test_fit_rule_g_for_r1():
    y = np.array([1, 1, 0, 0], dtype=bool)
    rank = rules.gpi_rank(["probable", "weakly", "weakly", "none"])
    fit = rules.fit_rule("R1", y, np.ones(4, dtype=bool), rank, np.zeros(4))
    # g = probable: J = 0.5; g = weakly: J = 1 - 0.5 = 0.5; tie -> probable (stricter)
    assert fit["g"] == "probable" and fit["j"] == 0.5


def test_fit_rule_refuses_rows_without_a_class():
    with pytest.raises(ValueError, match="lack a class"):
        rules.fit_rule(
            "R2", np.ones(3, bool), np.ones(3, bool), rules.gpi_rank(["none"] * 3), np.zeros(3)
        )
