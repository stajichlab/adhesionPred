"""Rules, gating, status and the `other_*` formulas."""

from pathlib import Path

import pytest
import yaml

from cellsurface_sorting_hat import engine
from cellsurface_sorting_hat.engine import (
    ConfigError,
    ModuleTable,
    collect_evidence,
    evaluate,
    load_config,
    referenced_modules,
)

R0 = "step1_rule@R0"
R2 = "step1_rule@R2"


def table(name, rows):
    return ModuleTable(name, {k: dict(v) for k, v in rows.items()})


def run(modules, status_of=None, ids=("P",), taxon=40, measured_call_of=None, call_status_of=None):
    cfg = load_config()
    status_of = status_of or (lambda module, t: ("unvalidated", "x"))
    records = evaluate(
        cfg,
        list(ids),
        dict.fromkeys(ids, taxon),
        modules,
        status_of,
        measured_call_of,
        call_status_of,
    )
    return {(r.protein, r.call, r.variant): r for r in records}


def ok(**fields):
    return {"state": "ok", **fields}


def base_modules(step1="not_called", **over):
    mods = {
        R0: table(R0, {"P": ok(call=step1)}),
        "repeat02": table("repeat02", {"P": ok(call="not_called")}),
        "repeat14": table("repeat14", {"P": ok(call="not_called")}),
        "pfam_adhesion": table("pfam_adhesion", {"P": ok(hit="0")}),
        "pfam_hydrophobin": table("pfam_hydrophobin", {"P": ok(hit="0")}),
        "pfam_hsba": table("pfam_hsba", {"P": ok(hit="0")}),
        "pfam_allergen": table("pfam_allergen", {"P": ok(hit="0")}),
        "antigen_lookup": table("antigen_lookup", {"P": ok(percentile="50")}),
        "allergen_homology": table("allergen_homology", {"P": ok(identity="0", coverage="0")}),
    }
    mods.update(over)
    return mods


def test_packaged_config_loads_and_has_the_spec_calls():
    names = [c["name"] for c in load_config().calls]
    assert names == [
        "signal_peptide_protein",
        "tandem_repeat_protein",
        "wall_family_domain",
        "hydrophobin_domain",
        "hsba_domain",
        "cell_wall_adhesion_candidate",
        "cocci_specificity_rank_top15",
        "serodiagnostic_marker_candidate",
        "iuis_allergen_similarity",
        "iuis_allergen_homolog",
        "other_not_surface",
        "other_surface_no_mechanism",
    ]


def test_only_available_step1_variants_get_per_variant_records():
    res = run(base_modules())
    variants = {v for (_, call, v) in res if call == "signal_peptide_protein"}
    assert variants == {"R0"}


def test_a_true_or_input_wins_over_an_unknown_input():
    mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})
    r = run(mods)[("P", "tandem_repeat_protein", "")]
    assert r.value == "called"


def test_unknown_or_input_with_a_false_other_is_unknown():
    mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
    assert run(mods)[("P", "tandem_repeat_protein", "")].value == "not_assessable"


def test_false_surface_makes_gated_antigen_false_even_when_antigen_is_unknown():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
    res = run(mods)
    assert res[("P", "cocci_specificity_rank_top15", "")].value == "not_assessable"
    assert res[("P", "serodiagnostic_marker_candidate", "R0")].value == "not_called"


def test_antigen_threshold_is_inclusive():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="15")}))
    assert run(mods)[("P", "cocci_specificity_rank_top15", "")].value == "called"
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="15.01")}))
    assert run(mods)[("P", "cocci_specificity_rank_top15", "")].value == "not_called"


def test_empty_or_nan_number_is_unknown():
    for bad in ("", "nan", "abc"):
        mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile=bad)}))
        assert run(mods)[("P", "cocci_specificity_rank_top15", "")].value == "not_assessable", bad


@pytest.mark.parametrize(
    "identity,coverage,pfam,expected",
    [
        ("82", "90", "0", "called"),
        ("69.9", "90", "0", "not_called"),
        ("82", "79", "0", "not_called"),
        ("40", "85", "1", "called"),
        ("0", "0", "0", "not_called"),
    ],
)
def test_allergen_two_tier_call(identity, coverage, pfam, expected):
    mods = base_modules(
        allergen_homology=table(
            "allergen_homology", {"P": ok(identity=identity, coverage=coverage)}
        ),
        pfam_allergen=table("pfam_allergen", {"P": ok(hit=pfam)}),
    )
    assert run(mods)[("P", "iuis_allergen_homolog", "")].value == expected


def test_allergen_needs_no_surface_call():
    mods = base_modules(step1="not_called")
    mods["allergen_homology"] = table("allergen_homology", {"P": ok(identity="90", coverage="95")})
    res = run(mods)
    assert res[("P", "signal_peptide_protein", "R0")].value == "not_called"
    assert res[("P", "iuis_allergen_homolog", "")].value == "called"


def test_other_not_surface_leaves_out_unknown_calls_and_names_them():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
    r = run(mods)[("P", "other_not_surface", "R0")]
    assert (r.value, r.other_basis) == ("called", "cocci_specificity_rank_top15")


def test_other_not_surface_is_false_when_a_mechanism_call_is_true():
    mods = base_modules(repeat02=table("repeat02", {"P": ok(call="called")}))
    assert run(mods)[("P", "other_not_surface", "R0")].value == "not_called"


def test_an_allergen_flag_does_not_change_the_other_calls():
    mods = base_modules(
        allergen_homology=table("allergen_homology", {"P": ok(identity="90", coverage="95")})
    )
    res = run(mods)
    assert res[("P", "iuis_allergen_homolog", "")].value == "called"
    assert res[("P", "other_not_surface", "R0")].value == "called"


def test_other_surface_no_mechanism():
    res = run(base_modules(step1="called"))
    assert res[("P", "other_surface_no_mechanism", "R0")].value == "called"
    assert res[("P", "other_not_surface", "R0")].value == "not_called"


def test_other_is_unknown_when_the_surface_call_is_unknown():
    mods = base_modules()
    mods[R0] = table(R0, {"P": {"state": "error"}})
    res = run(mods)
    assert res[("P", "other_not_surface", "R0")].value == "not_assessable"
    assert res[("P", "other_surface_no_mechanism", "R0")].value == "not_assessable"


def test_other_can_be_true_in_one_variant_and_unknown_in_another():
    mods = base_modules()
    mods[R2] = table(R2, {"P": {"state": "error"}})
    res = run(mods)
    assert res[("P", "other_not_surface", "R0")].value == "called"
    assert res[("P", "other_not_surface", "R2")].value == "not_assessable"


def test_status_is_the_weakest_of_the_modules_that_decided_the_result():
    statuses = {
        R0: ("estimated", "taxon:40"),
        "repeat02": ("smoke", "taxon:30"),
        "repeat14": ("unvalidated", "no status_source"),
    }
    mods = base_modules(step1="called")
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})  # repeat14 stays not_called
    res = run(mods, status_of=lambda m, t: statuses.get(m, ("unvalidated", "")))
    # tandem_repeat_protein is true because of repeat02 only: status smoke, not unvalidated
    assert res[("P", "tandem_repeat_protein", "")].status == "smoke"
    # the gated call is true because of repeat02 and step 1: weakest is smoke
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert (r.value, r.status) == ("called", "smoke")
    assert "step1_rule@R0:taxon:40" in r.status_basis


def test_invalid_protein_rows_make_every_call_unknown():
    mods = base_modules()
    for name in mods:
        mods[name] = table(name, {"P": {"state": "na_invalid"}})
    res = run(mods)
    assert {r.value for r in res.values()} == {"not_assessable"}


def _config_with(tmp_path, mutate):
    packaged = Path(engine.__file__).with_name("categories.yaml")
    data = yaml.safe_load(packaged.read_text())
    mutate(data)
    path = tmp_path / "c.yaml"
    path.write_text(yaml.safe_dump(data))
    return path


@pytest.mark.parametrize(
    "mutate,message",
    [
        (lambda d: d.update(default_gate="nope"), "default_gate"),
        (lambda d: d["calls"].append(dict(d["calls"][0])), "duplicate call"),
        (lambda d: d["calls"][1].update(expr={"ref": "later_call"}), "unknown or later"),
        (lambda d: d["calls"][1].update(expr={"bogus": 1}), "bad node"),
        (lambda d: d["calls"][6]["expr"]["test"].update(op="~"), "bad operator"),
        (lambda d: d["calls"][6]["expr"]["test"].update(value="$nope"), "unknown threshold"),
        (
            lambda d: d["calls"][1].update(expr={"ref": "signal_peptide_protein"}),
            "ungated call cannot use",
        ),
    ],
)
def test_bad_config_is_refused(tmp_path, mutate, message):
    with pytest.raises(ConfigError, match=message):
        load_config(_config_with(tmp_path, mutate))


def test_an_unknown_result_has_status_unvalidated_and_no_basis():
    # R0 is in state error; repeat02 is called with status smoke. The gated call is unknown.
    mods = base_modules(step1="called")
    mods[R0] = table(R0, {"P": {"state": "error"}})
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})
    res = run(mods, status_of=lambda m, t: ("smoke", "t"))
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert (r.value, r.status, r.status_basis) == ("not_assessable", "unvalidated", "")
    # an OR with a false input and an unknown input is unknown and also has no deciding module
    mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
    r = run(mods, status_of=lambda m, t: ("estimated", "t"))[("P", "tandem_repeat_protein", "")]
    assert (r.value, r.status) == ("not_assessable", "unvalidated")


@pytest.mark.parametrize(
    "identity,length,expected",
    [
        ("35", "80", "called"),
        ("34.9", "120", "not_called"),
        ("60", "79", "not_called"),
        ("", "90", "not_assessable"),
    ],
)
def test_iuis_allergen_similarity_is_the_35_percent_80_aa_evidence_tier(identity, length, expected):
    mods = base_modules(
        allergen_homology=table(
            "allergen_homology", {"P": ok(identity=identity, aligned_length=length, coverage="0")}
        )
    )
    assert run(mods)[("P", "iuis_allergen_similarity", "")].value == expected


def test_evidence_rows_are_collected_only_for_ok_rows_with_a_value():
    cfg = load_config()
    mods = {
        "antigen_lookup": table(
            "antigen_lookup",
            {
                "A": ok(percentile="7.1", antigenicity=""),
                "B": {"state": "not_applicable", "percentile": "9"},
            },
        )
    }
    rows = collect_evidence(cfg, ["A", "B"], mods)
    assert rows == [("A", "antigen_lookup", "percentile", "7.1")]


def test_referenced_modules_expand_the_step1_variants():
    names = referenced_modules(load_config())
    assert "step1_rule@R0" in names and "step1_ml@card" in names
    assert {
        "repeat02",
        "repeat14",
        "pfam_adhesion",
        "pfam_hydrophobin",
        "pfam_hsba",
        "pfam_allergen",
        "antigen_lookup",
        "allergen_homology",
    } <= set(names)
    assert not any("{step1}" in n for n in names)


@pytest.mark.parametrize(
    "mutate,message",
    [
        (
            lambda d: d["calls"][-2].update(surface="tandem_repeat_protein"),
            "needs an earlier per_variant surface",
        ),
        (
            lambda d: d["calls"][1].update(expr={"call": "{step1}"}),
            "only allowed in a per_variant call",
        ),
        (lambda d: d.update(evidence=["no_dot"]), "evidence entry"),
    ],
)
def test_more_bad_config_is_refused(tmp_path, mutate, message):
    with pytest.raises(ConfigError, match=message):
        load_config(_config_with(tmp_path, mutate))


@pytest.mark.parametrize("bad", ["inf", "-inf", "Infinity"])
def test_a_non_finite_number_is_unknown(bad):
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile=bad)}))
    assert run(mods)[("P", "cocci_specificity_rank_top15", "")].value == "not_assessable"


def test_a_false_and_takes_its_status_from_the_false_inputs_only():
    # surface (R0) is false with status smoke; antigen is true with status estimated
    statuses = {R0: ("smoke", "t"), "antigen_lookup": ("estimated", "t")}
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="5")}))
    res = run(mods, status_of=lambda m, t: statuses.get(m, ("unvalidated", "")))
    r = res[("P", "serodiagnostic_marker_candidate", "R0")]
    assert (r.value, r.status) == ("not_called", "smoke")
    assert r.status_basis == "step1_rule@R0:t"


def _status_fn(mapping):
    return lambda m, t: (mapping[m], "t") if m in mapping else ("unvalidated", "")


def _basis_modules(record):
    return {item.split(":", 1)[0] for item in record.status_basis.split(";") if item}


def test_other_true_status_uses_surface_and_not_called_mechanisms_only():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
    statuses = _status_fn(
        {
            R0: "estimated",
            "repeat02": "smoke",
            "repeat14": "estimated",
            "pfam_adhesion": "estimated",
            "pfam_hydrophobin": "estimated",
            "pfam_hsba": "estimated",
            "allergen_homology": "estimated",
            "pfam_allergen": "estimated",
        }
    )
    r = run(mods, status_of=statuses)[("P", "other_not_surface", "R0")]
    assert r.value == "called"
    assert r.status == "smoke"
    assert _basis_modules(r) == {
        R0,
        "repeat02",
        "repeat14",
        "pfam_adhesion",
        "pfam_hydrophobin",
        "pfam_hsba",
    }  # antigen_lookup is unknown and does not contribute; allergen flags are not mechanisms


def test_other_false_from_surface_mismatch_uses_surface_only():
    statuses = _status_fn({R0: "estimated", "repeat02": "smoke"})
    r = run(base_modules(step1="called"), status_of=statuses)[("P", "other_not_surface", "R0")]
    assert r.value == "not_called"
    assert (r.status, _basis_modules(r)) == ("estimated", {R0})


def test_other_false_from_one_called_mechanism_uses_that_mechanism_only():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="5")}))
    statuses = _status_fn({R0: "smoke", "antigen_lookup": "estimated"})
    r = run(mods, status_of=statuses)[("P", "other_not_surface", "R0")]
    assert r.value == "not_called"
    assert (r.status, _basis_modules(r)) == ("estimated", {"antigen_lookup"})


def test_other_false_from_two_called_mechanisms_takes_the_weakest_and_lists_both():
    mods = base_modules(
        antigen_lookup=table("antigen_lookup", {"P": ok(percentile="5")}),
        repeat02=table("repeat02", {"P": ok(call="called")}),
    )
    statuses = _status_fn({R0: "estimated", "antigen_lookup": "estimated", "repeat02": "smoke"})
    r = run(mods, status_of=statuses)[("P", "other_not_surface", "R0")]
    assert r.value == "not_called"
    assert r.status == "smoke"
    assert _basis_modules(r) == {"antigen_lookup", "repeat02"}


def test_or_true_with_two_true_inputs_takes_the_weakest_and_lists_both():
    mods = base_modules()
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})
    mods["repeat14"] = table("repeat14", {"P": ok(call="called")})
    statuses = _status_fn({"repeat02": "estimated", "repeat14": "smoke"})
    r = run(mods, status_of=statuses)[("P", "tandem_repeat_protein", "")]
    assert (r.value, r.status) == ("called", "smoke")
    assert _basis_modules(r) == {"repeat02", "repeat14"}


def test_or_false_takes_status_from_all_inputs():
    statuses = _status_fn({"repeat02": "estimated", "repeat14": "smoke"})
    r = run(base_modules(), status_of=statuses)[("P", "tandem_repeat_protein", "")]
    assert (r.value, r.status) == ("not_called", "smoke")
    assert _basis_modules(r) == {"repeat02", "repeat14"}


# ---- a status counts for the call on which the module was measured (decision of 2026-10-06) ----


def _measured(mapping):
    """``measured_call_of`` from {module: call}; a module not listed has an entry without ``call=``."""
    return lambda module, taxon: mapping.get(module)


def _all_estimated(module, taxon):
    return ("estimated", f"taxon:{taxon}")


def _allergen_modules(identity="100"):
    return base_modules(
        allergen_homology=table(
            "allergen_homology",
            {"P": ok(identity=identity, coverage="100", aligned_length="1000")},
        )
    )


def test_a_module_measured_on_call_a_is_unvalidated_for_call_b_that_reads_it():
    mods = _allergen_modules()
    res = run(
        mods,
        _all_estimated,
        measured_call_of=_measured({"allergen_homology": "iuis_allergen_similarity"}),
    )
    similarity = res[("P", "iuis_allergen_similarity", "")]
    homolog = res[("P", "iuis_allergen_homolog", "")]
    assert (similarity.value, similarity.status) == ("called", "estimated")
    assert (homolog.value, homolog.status) == ("called", "unvalidated")
    assert "allergen_homology:module measured on call iuis_allergen_similarity" in (
        homolog.status_basis
    )


def test_the_measured_call_is_the_one_that_keeps_its_status():
    mods = _allergen_modules()
    res = run(
        mods,
        _all_estimated,
        measured_call_of=_measured({"allergen_homology": "iuis_allergen_homolog"}),
    )
    assert res[("P", "iuis_allergen_homolog", "")].status == "estimated"
    assert res[("P", "iuis_allergen_similarity", "")].status == "unvalidated"


def test_a_legacy_entry_without_call_keeps_the_status_in_every_call():
    mods = _allergen_modules()
    res = run(mods, _all_estimated, measured_call_of=_measured({}))
    assert res[("P", "iuis_allergen_similarity", "")].status == "estimated"
    assert res[("P", "iuis_allergen_homolog", "")].status == "estimated"
    res = run(mods, _all_estimated)  # no measured_call_of at all
    assert res[("P", "iuis_allergen_homolog", "")].status == "estimated"


def _composite_modules():
    # R0 called and pfam_adhesion hit: cell_wall_adhesion_candidate is called through wall_family_domain
    return base_modules(step1="called", pfam_adhesion=table("pfam_adhesion", {"P": ok(hit="1")}))


def test_a_composite_call_follows_the_leaf_call_through_ref():
    measured = {R0: "signal_peptide_protein", "pfam_adhesion": "wall_family_domain"}
    res = run(_composite_modules(), _all_estimated, measured_call_of=_measured(measured))
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert (r.value, r.status) == ("called", "estimated")
    # pfam_adhesion measured on another call: the composite call loses it, the leaf call keeps R0
    measured["pfam_adhesion"] = "tandem_repeat_protein"
    res = run(_composite_modules(), _all_estimated, measured_call_of=_measured(measured))
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert (r.value, r.status) == ("called", "unvalidated")
    assert "pfam_adhesion:module measured on call tandem_repeat_protein" in r.status_basis
    assert res[("P", "signal_peptide_protein", "R0")].status == "estimated"
    assert res[("P", "wall_family_domain", "")].status == "unvalidated"


def test_a_step_1_module_measured_on_a_composite_call_does_not_count_for_its_own_leaf_call():
    measured = {R0: "cell_wall_adhesion_candidate", "pfam_adhesion": "wall_family_domain"}
    res = run(_composite_modules(), _all_estimated, measured_call_of=_measured(measured))
    assert res[("P", "signal_peptide_protein", "R0")].status == "unvalidated"
    assert res[("P", "cell_wall_adhesion_candidate", "R0")].status == "unvalidated"


def test_the_weakest_deciding_module_rule_still_holds_with_measured_calls():
    mods = base_modules()
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})
    mods["repeat14"] = table("repeat14", {"P": ok(call="called")})
    statuses = _status_fn({"repeat02": "estimated", "repeat14": "smoke"})
    measured = {"repeat02": "tandem_repeat_protein", "repeat14": "tandem_repeat_protein"}
    r = run(mods, status_of=statuses, measured_call_of=_measured(measured))[
        ("P", "tandem_repeat_protein", "")
    ]
    assert (r.value, r.status) == ("called", "smoke")


# ---- call status hook (per-call status files) ----------------------------------------------------

REPEAT = "tandem_repeat_protein"


def _hook(table):
    """A ``call_status_of`` hook from {(call, variant): result}; it records what it was asked."""
    asked = []

    def hook(call, variant, taxon):
        asked.append((call, variant, taxon))
        return table.get((call, variant))

    hook.asked = asked
    return hook


def _fields(res):
    return {k: (r.value, r.status, r.status_basis, r.other_basis) for k, r in res.items()}


def test_a_hook_that_never_applies_changes_nothing():
    mods = base_modules(step1="called", repeat02=table("repeat02", {"P": ok(call="called")}))
    plain = run(mods, _all_estimated, measured_call_of=_measured({R0: "signal_peptide_protein"}))
    hooked = run(
        mods,
        _all_estimated,
        measured_call_of=_measured({R0: "signal_peptide_protein"}),
        call_status_of=_hook({}),
    )
    assert _fields(hooked) == _fields(plain)


def test_a_call_status_replaces_the_module_statuses_of_that_call():
    mods = base_modules(repeat02=table("repeat02", {"P": ok(call="called")}))
    hook = _hook({(REPEAT, ""): ("smoke", "call:tandem_repeat_protein:taxon:40")})
    r = run(mods, _all_estimated, call_status_of=hook)[("P", REPEAT, "")]
    assert (r.value, r.status, r.status_basis) == (
        "called",
        "smoke",
        "call:tandem_repeat_protein:taxon:40",
    )
    stronger = _hook({(REPEAT, ""): ("estimated", "call:tandem_repeat_protein:taxon:40")})
    r = run(mods, lambda m, t: ("unvalidated", "x"), call_status_of=stronger)[("P", REPEAT, "")]
    assert r.status == "estimated"  # the call was measured; the modules alone were unvalidated


def test_an_unknown_value_never_takes_a_call_status():
    bad = {"P": {"state": "error"}}
    mods = base_modules(repeat02=table("repeat02", bad), repeat14=table("repeat14", bad))
    hook = _hook({(REPEAT, ""): ("estimated", "call:tandem_repeat_protein:taxon:40")})
    r = run(mods, _all_estimated, call_status_of=hook)[("P", REPEAT, "")]
    assert (r.value, r.status, r.status_basis) == ("not_assessable", "unvalidated", "")


def test_a_composite_call_takes_the_weakest_status_over_its_measured_leaves():
    mods = base_modules(step1="called", repeat02=table("repeat02", {"P": ok(call="called")}))
    hook = _hook({(REPEAT, ""): ("smoke", "call:tandem_repeat_protein:taxon:40")})
    res = run(
        mods,
        _all_estimated,
        measured_call_of=_measured({R0: "signal_peptide_protein"}),
        call_status_of=hook,
    )
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert (r.value, r.status) == ("called", "smoke")
    assert r.status_basis == "call:tandem_repeat_protein:taxon:40;step1_rule@R0:taxon:40"
    # the R0 module measured on another call does not count for the R0 leaf
    res = run(
        mods,
        _all_estimated,
        measured_call_of=_measured({R0: "something_else"}),
        call_status_of=hook,
    )
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert r.status == "unvalidated"
    assert r.status_basis == (
        "call:tandem_repeat_protein:taxon:40;step1_rule@R0:module measured on call something_else"
    )


def test_a_false_or_takes_its_status_from_all_its_false_inputs():
    mods = base_modules()  # both repeat detectors not_called
    hook = _hook({(REPEAT, ""): ("smoke", "call:tandem_repeat_protein:taxon:40")})
    r = run(mods, _all_estimated, call_status_of=hook)[("P", REPEAT, "")]
    assert (r.value, r.status, r.status_basis) == (
        "not_called",
        "smoke",
        "call:tandem_repeat_protein:taxon:40",
    )


def test_a_false_and_takes_the_call_status_of_the_false_leaf_and_the_module_status_of_the_other():
    mods = base_modules(step1="called")  # candidate is false because repeat and domain are false
    hook = _hook({(REPEAT, ""): ("smoke", "call:tandem_repeat_protein:taxon:40")})
    res = run(
        mods,
        lambda m, t: ("estimated", "x") if m != "pfam_adhesion" else ("unvalidated", "x"),
        measured_call_of=_measured({}),
        call_status_of=hook,
    )
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert r.value == "not_called"
    assert r.status == "unvalidated"  # the weakest of smoke (repeat call) and unvalidated (pfam)
    # items follow the sorted (module, leaf) pairs; the call item stands where its first module would
    assert r.status_basis == "pfam_adhesion:x;call:tandem_repeat_protein:taxon:40"


def test_other_surface_no_mechanism_takes_the_call_status_of_a_called_mechanism():
    mods = base_modules(step1="called", repeat02=table("repeat02", {"P": ok(call="called")}))
    hook = _hook({(REPEAT, ""): ("smoke", "call:tandem_repeat_protein:taxon:40")})
    res = run(
        mods,
        _all_estimated,
        measured_call_of=_measured({R0: "signal_peptide_protein"}),
        call_status_of=hook,
    )
    r = res[("P", "other_surface_no_mechanism", "R0")]
    assert (r.value, r.status) == ("not_called", "smoke")
    assert r.status_basis == "call:tandem_repeat_protein:taxon:40"


def test_the_variant_of_a_leaf_is_the_label_only_for_a_per_variant_leaf():
    mods = base_modules(
        step1="called",
        **{R2: table(R2, {"P": ok(call="called")})},
        repeat02=table("repeat02", {"P": ok(call="called")}),
    )
    hook = _hook(
        {("signal_peptide_protein", "R0"): ("smoke", "call:signal_peptide_protein:taxon:40")}
    )
    res = run(mods, _all_estimated, measured_call_of=_measured({}), call_status_of=hook)
    assert res[("P", "signal_peptide_protein", "R0")].status == "smoke"
    assert res[("P", "signal_peptide_protein", "R2")].status == "estimated"  # no file for R2
    candidate_r2 = res[("P", "cell_wall_adhesion_candidate", "R2")]
    assert "call:signal_peptide_protein" not in candidate_r2.status_basis
    asked = set(hook.asked)
    assert ("signal_peptide_protein", "R2", 40) in asked
    assert (REPEAT, "", 40) in asked and (REPEAT, "R0", 40) not in asked  # plain leaf: no variant


def test_a_stale_call_status_falls_back_to_the_module_status_and_says_why():
    mods = base_modules(repeat02=table("repeat02", {"P": ok(call="called")}))
    hook = _hook({(REPEAT, ""): (None, "reads differ")})
    r = run(mods, _all_estimated, call_status_of=hook)[("P", REPEAT, "")]
    assert r.status == "estimated"
    assert r.status_basis == "repeat02:call status stale (reads differ); taxon:40"


def test_hydrophobin_hit_is_not_a_wall_family_hit_and_counts_as_mechanism():
    mods = base_modules(
        step1="called", pfam_hydrophobin=table("pfam_hydrophobin", {"P": ok(hit="1")})
    )
    res = run(mods)
    assert res[("P", "hydrophobin_domain", "")].value == "called"
    assert res[("P", "wall_family_domain", "")].value == "not_called"
    assert res[("P", "cell_wall_adhesion_candidate", "R0")].value == "not_called"
    assert res[("P", "other_surface_no_mechanism", "R0")].value == "not_called"


def test_hsba_hit_is_a_separate_call_and_counts_as_mechanism():
    mods = base_modules(step1="called", pfam_hsba=table("pfam_hsba", {"P": ok(hit="1")}))
    res = run(mods)
    assert res[("P", "hsba_domain", "")].value == "called"
    assert res[("P", "hydrophobin_domain", "")].value == "not_called"
    assert res[("P", "other_surface_no_mechanism", "R0")].value == "not_called"
