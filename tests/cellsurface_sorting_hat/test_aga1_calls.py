"""The six AGA1-class calls of categories.yaml on a synthetic work directory (every combination of the three-valued inputs)."""

import itertools

import pytest
from engine_helpers import R0, base_modules, ok, run_eval, table

from cellsurface_sorting_hat.engine import call_eligible, load_config, reads_of_call

CALLED, NOT, UNK = "called", "not_called", "not_assessable"


def kleene_and(*values):
    if NOT in values:
        return NOT
    return UNK if UNK in values else CALLED


def modules(gpi, run, sp, repeat14, model="not_called"):
    mods = base_modules(step1=sp, repeat14=table("repeat14", {"P": ok(call=repeat14)}))
    mods["gpi"] = table("gpi", {"P": {"state": UNK} if gpi == UNK else ok(call=gpi)})
    mods["disorder"] = table(
        "disorder", {"P": {"state": "error"} if run is None else ok(longest_run=str(run))}
    )
    mods["aga1_model"] = table("aga1_model", {"P": ok(call=model)})
    return mods


def values(mods):
    res = run_eval(mods)
    got = {}
    for (protein, call, variant), r in res.items():
        if protein == "P" and variant in ("", "R0"):
            got[call] = r.value
    return got


RUNS = {
    250: CALLED,
    200: CALLED,
    199: NOT,
    0: NOT,
    None: UNK,
}  # 200 is the cut and counts; None = no usable disorder row


@pytest.mark.parametrize(
    "gpi, run, sp, rep",
    list(itertools.product([CALLED, NOT, UNK], [250, 199, None], [CALLED, NOT], [CALLED, NOT])),
)
def test_the_calls_follow_kleene_logic(gpi, run, sp, rep):
    got = values(modules(gpi, run, sp, rep))
    long_run = RUNS[run]
    assert got["gpi_anchored_protein"] == gpi
    assert got["long_disordered_protein"] == long_run
    both = kleene_and(gpi, long_run)
    assert got["gpi_long_disorder_protein"] == both
    broad = kleene_and(sp, both)
    assert got["ser_thr_gpi_wall_candidate"] == broad
    assert got["flocculin_like_candidate"] == kleene_and(broad, rep)


@pytest.mark.parametrize("run, expected", [(200, CALLED), (199, NOT)])
def test_the_disorder_cut_is_inclusive_and_read_from_thresholds(run, expected):
    assert load_config().thresholds["disorder_run_min"] == 200
    assert values(modules(CALLED, run, CALLED, CALLED))["long_disordered_protein"] == expected


@pytest.mark.parametrize(
    "model, sp, expected",
    [(CALLED, CALLED, CALLED), (CALLED, NOT, NOT), (NOT, CALLED, NOT), (NOT, NOT, NOT)],
)
def test_the_model_call_needs_the_model_and_a_signal_peptide(model, sp, expected):
    assert values(modules(CALLED, 250, sp, NOT, model=model))["model_surface_candidate"] == expected


def test_an_absent_aga1_module_makes_its_calls_not_assessable_not_zero():
    mods = modules(CALLED, 250, CALLED, CALLED)
    for name in ("gpi", "disorder", "aga1_model"):
        del mods[name]
    got = values(mods)
    for call in (
        "gpi_anchored_protein",
        "long_disordered_protein",
        "gpi_long_disorder_protein",
        "model_surface_candidate",
    ):
        assert got[call] == UNK
    assert got["ser_thr_gpi_wall_candidate"] == UNK and got["flocculin_like_candidate"] == UNK


def test_call_structure_matches_what_the_engine_can_measure():
    cfg = load_config()
    by = {c["name"]: c for c in cfg.calls}
    assert not by["gpi_long_disorder_protein"].get("per_variant")  # ungated, no ref
    assert (
        call_eligible(cfg, "gpi_long_disorder_protein")[0] is True
    )  # reads two modules: a call status file can be measured
    assert sorted(reads_of_call(cfg, "gpi_long_disorder_protein")) == ["disorder", "gpi"]
    for name in ("gpi_anchored_protein", "long_disordered_protein"):
        assert (
            call_eligible(cfg, name)[0] is False
        )  # they read one module: measured with --module, not as a call
        assert len(reads_of_call(cfg, name)) == 1
    for name in (
        "ser_thr_gpi_wall_candidate",
        "flocculin_like_candidate",
        "model_surface_candidate",
    ):
        assert by[name]["per_variant"] is True
        assert call_eligible(cfg, name) == (
            False,
            "contains ref",
        )  # derived: the status is the weakest of the leaf calls
    narrow = reads_of_call(cfg, "flocculin_like_candidate", "R0")
    assert "repeat14" in narrow and "repeat02" not in narrow  # detector 14 only, as measured


def test_the_aga1_calls_leave_the_existing_labels_and_calls_unchanged():
    got = values(modules(CALLED, 250, CALLED, NOT))  # GPI and a long run, but no repeat
    assert got["ser_thr_gpi_wall_candidate"] == CALLED
    assert got["tandem_repeat_protein"] == NOT and got["cell_wall_adhesion_candidate"] == NOT
    assert (
        got["other_surface_no_mechanism"] == CALLED
    )  # no existing mechanism evidence: the existing label rule is not changed here
    assert R0 == "step1_rule@R0"
