"""Engine behaviour of PR 3 at the level of ``evaluate``.

Tests that fail on ``main`` carry ``xfail(strict=True)`` with the task id. The task that lands the
behaviour removes the mark in the same commit.
"""

from engine_helpers import base_modules, config_with, ok, run_eval, table

NA = "not_applicable"
ESTIMATED = lambda module, taxon: ("estimated", "x")  # noqa: E731


def antigen_na():
    return table("antigen_lookup", {"P": {"state": "not_applicable"}})


# ---- Task 3.2: unknown mechanism input, not_applicable -----------------------------------------


def test_probe_p2_an_unknown_mechanism_makes_other_surface_no_mechanism_not_assessable():
    mods = base_modules(step1="called", pfam_hsba=table("pfam_hsba", {"P": {"state": "error"}}))
    got = run_eval(mods)
    r = got[("P", "other_surface_no_mechanism", "R0")]
    assert r.value == "not_assessable"
    assert r.other_basis == "hsba_domain"
    assert got[("P", "surface_attachment_candidate", "R0")].value == "not_assessable"


def test_a_called_mechanism_wins_over_another_unknown_one():
    mods = base_modules(
        step1="called",
        repeat02=table("repeat02", {"P": ok(call="called")}),
        pfam_hsba=table("pfam_hsba", {"P": {"state": "error"}}),
    )
    assert run_eval(mods)[("P", "other_surface_no_mechanism", "R0")].value == "not_called"


def test_probe_p3_not_applicable_antigen_rank():
    mods = base_modules(step1="called", antigen_lookup=antigen_na())
    got = run_eval(mods, status_of=ESTIMATED)
    assert got[("P", "cocci_specificity_rank_top15", "")].value == NA
    sero = got[("P", "serodiagnostic_marker_candidate", "R0")]
    assert sero.value == "not_called"
    assert sero.status_basis == "not_applicable: antigen_lookup"
    assert sero.status == "unvalidated"  # an NA leaf contributes no status
    assert got[("P", "other_surface_no_mechanism", "R0")].value == "called"


def test_a_not_applicable_record_has_no_status_and_names_the_module():
    got = run_eval(base_modules(step1="called", antigen_lookup=antigen_na()), status_of=ESTIMATED)
    r = got[("P", "cocci_specificity_rank_top15", "")]
    assert r.status == "unvalidated"
    assert r.status_basis == "not_applicable: antigen_lookup"


def test_a_secreted_protein_outside_rs_is_never_a_serodiagnostic_candidate():
    """FD-A1: dropping NA from the operands would make this call true for every secreted protein."""
    for step1 in ("called", "not_called"):
        mods = base_modules(step1=step1, antigen_lookup=antigen_na())
        assert run_eval(mods)[("P", "serodiagnostic_marker_candidate", "R0")].value != "called"


def test_not_of_a_not_applicable_call_is_not_applicable(tmp_path):
    cfg = config_with(
        tmp_path,
        extra_calls=[
            {"name": "not_rank", "expr": {"not": {"ref": "cocci_specificity_rank_top15"}}},
        ],
    )
    got = run_eval(base_modules(step1="called", antigen_lookup=antigen_na()), cfg=cfg)
    assert got[("P", "not_rank", "")].value == NA


def test_call_hash_changes_with_the_engine_semantics():
    from cellsurface_sorting_hat import engine

    # value of call_hash(cfg, "tandem_repeat_protein") on main at 14a3011, ENGINE_SEMANTICS "1"
    old = "94228eccdcd64a1065c0750f1ad680d9a5895619057c2fad6734ddd0e365afdd"
    got = engine.call_hash(engine.load_config(), "tandem_repeat_protein")
    assert (engine.ENGINE_SEMANTICS == "1") == (got == old)  # the hash follows the semantics


def test_engine_semantics_were_bumped():
    from cellsurface_sorting_hat import engine

    assert engine.ENGINE_SEMANTICS != "1"


# ---- Task 3.4: the narrow call reads its own flag ------------------------------------------------


def _adhesion(hit, narrow):
    return table("pfam_adhesion", {"P": ok(hit=hit, hit_narrow=narrow)})


def test_a_hit_in_a_family_that_is_not_narrow_feeds_only_the_broad_call():
    got = run_eval(base_modules(step1="called", pfam_adhesion=_adhesion("1", "0")))
    assert got[("P", "wall_family_domain", "")].value == "called"
    assert got[("P", "wall_family_narrow_domain", "")].value == "not_called"
    assert got[("P", "cell_wall_adhesion_candidate", "R0")].value == "not_called"
    assert got[("P", "surface_attachment_candidate", "R0")].value == "called"
    assert got[("P", "other_surface_no_mechanism", "R0")].value == "not_called"


def test_a_hit_in_a_narrow_family_feeds_both_calls():
    got = run_eval(base_modules(step1="called", pfam_adhesion=_adhesion("1", "1")))
    assert got[("P", "wall_family_narrow_domain", "")].value == "called"
    assert got[("P", "cell_wall_adhesion_candidate", "R0")].value == "called"
    assert got[("P", "surface_attachment_candidate", "R0")].value == "called"


def test_the_narrow_call_needs_a_signal_peptide():
    got = run_eval(base_modules(step1="not_called", pfam_adhesion=_adhesion("1", "1")))
    assert got[("P", "wall_family_narrow_domain", "")].value == "called"
    assert got[("P", "cell_wall_adhesion_candidate", "R0")].value == "not_called"


# ---- Task 3.5: derived marker -------------------------------------------------------------------


def test_composite_and_other_records_name_their_leaf_calls():
    mods = base_modules(step1="called", repeat02=table("repeat02", {"P": ok(call="called")}))
    got = run_eval(mods)
    cwac = got[("P", "cell_wall_adhesion_candidate", "R0")]
    assert cwac.value == "called"
    assert cwac.derived_from == "signal_peptide_protein,tandem_repeat_protein"
    sac = got[("P", "surface_attachment_candidate", "R0")]
    assert "tandem_repeat_protein" in sac.derived_from.split(",")
    other = run_eval(base_modules(step1="called"))[("P", "other_surface_no_mechanism", "R0")]
    assert other.value == "called"
    assert "signal_peptide_protein" in other.derived_from.split(",")


def test_leaf_records_have_no_derived_marker():
    got = run_eval(
        base_modules(step1="called", repeat02=table("repeat02", {"P": ok(call="called")}))
    )
    for call in ("signal_peptide_protein", "tandem_repeat_protein", "wall_family_domain"):
        variant = "R0" if call == "signal_peptide_protein" else ""
        assert got[("P", call, variant)].derived_from == ""
