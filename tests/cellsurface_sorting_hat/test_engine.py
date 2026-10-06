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


def run(modules, status_of=None, ids=("P",), taxon=40):
    cfg = load_config()
    status_of = status_of or (lambda module, t: ("unvalidated", "x"))
    records = evaluate(cfg, list(ids), dict.fromkeys(ids, taxon), modules, status_of)
    return {(r.protein, r.call, r.variant): r for r in records}


def ok(**fields):
    return {"state": "ok", **fields}


def base_modules(step1="not_called", **over):
    mods = {
        R0: table(R0, {"P": ok(call=step1)}),
        "repeat02": table("repeat02", {"P": ok(call="not_called")}),
        "repeat14": table("repeat14", {"P": ok(call="not_called")}),
        "pfam_adhesion": table("pfam_adhesion", {"P": ok(hit="0")}),
        "pfam_allergen": table("pfam_allergen", {"P": ok(hit="0")}),
        "antigen_lookup": table("antigen_lookup", {"P": ok(percentile="50")}),
        "allergen_homology": table("allergen_homology", {"P": ok(identity="0", coverage="0")}),
    }
    mods.update(over)
    return mods


def test_packaged_config_loads_and_has_the_spec_calls():
    names = [c["name"] for c in load_config().calls]
    assert names == [
        "surface_glycoprotein",
        "adhesion_repeat",
        "adhesion_domain",
        "cell_wall_adhesion_ungated",
        "cell_wall_adhesion_candidate",
        "antigen_candidate",
        "antigen_candidate_surface",
        "allergen_homolog_hit",
        "allergen_candidate",
        "other_not_surface",
        "other_surface_no_mechanism",
    ]


def test_only_available_step1_variants_get_per_variant_records():
    res = run(base_modules())
    variants = {v for (_, call, v) in res if call == "surface_glycoprotein"}
    assert variants == {"R0"}


def test_a_true_or_input_wins_over_an_unknown_input():
    mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})
    r = run(mods)[("P", "adhesion_repeat", "")]
    assert r.value == "called"


def test_unknown_or_input_with_a_false_other_is_unknown():
    mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
    assert run(mods)[("P", "adhesion_repeat", "")].value == "not_assessable"


def test_false_surface_makes_gated_antigen_false_even_when_antigen_is_unknown():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
    res = run(mods)
    assert res[("P", "antigen_candidate", "")].value == "not_assessable"
    assert res[("P", "antigen_candidate_surface", "R0")].value == "not_called"


def test_antigen_threshold_is_inclusive():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="15")}))
    assert run(mods)[("P", "antigen_candidate", "")].value == "called"
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="15.01")}))
    assert run(mods)[("P", "antigen_candidate", "")].value == "not_called"


def test_empty_or_nan_number_is_unknown():
    for bad in ("", "nan", "abc"):
        mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile=bad)}))
        assert run(mods)[("P", "antigen_candidate", "")].value == "not_assessable", bad


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
    assert run(mods)[("P", "allergen_candidate", "")].value == expected


def test_allergen_needs_no_surface_call():
    mods = base_modules(step1="not_called")
    mods["allergen_homology"] = table("allergen_homology", {"P": ok(identity="90", coverage="95")})
    res = run(mods)
    assert res[("P", "surface_glycoprotein", "R0")].value == "not_called"
    assert res[("P", "allergen_candidate", "")].value == "called"


def test_other_not_surface_leaves_out_unknown_calls_and_names_them():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
    r = run(mods)[("P", "other_not_surface", "R0")]
    assert (r.value, r.other_basis) == ("called", "antigen_candidate")


def test_other_not_surface_is_false_when_a_mechanism_call_is_true():
    mods = base_modules(
        allergen_homology=table("allergen_homology", {"P": ok(identity="90", coverage="95")})
    )
    assert run(mods)[("P", "other_not_surface", "R0")].value == "not_called"


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
    # adhesion_repeat is true because of repeat02 only: status smoke, not unvalidated
    assert res[("P", "adhesion_repeat", "")].status == "smoke"
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
        (lambda d: d["calls"][5]["expr"]["test"].update(op="~"), "bad operator"),
        (lambda d: d["calls"][5]["expr"]["test"].update(value="$nope"), "unknown threshold"),
        (
            lambda d: d["calls"][1].update(expr={"ref": "surface_glycoprotein"}),
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
    r = run(mods, status_of=lambda m, t: ("estimated", "t"))[("P", "adhesion_repeat", "")]
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
def test_allergen_homolog_hit_is_the_35_percent_80_aa_evidence_tier(identity, length, expected):
    mods = base_modules(
        allergen_homology=table(
            "allergen_homology", {"P": ok(identity=identity, aligned_length=length, coverage="0")}
        )
    )
    assert run(mods)[("P", "allergen_homolog_hit", "")].value == expected


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
        "pfam_allergen",
        "antigen_lookup",
        "allergen_homology",
    } <= set(names)
    assert not any("{step1}" in n for n in names)


@pytest.mark.parametrize(
    "mutate,message",
    [
        (
            lambda d: d["calls"][-2].update(surface="adhesion_repeat"),
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
    assert run(mods)[("P", "antigen_candidate", "")].value == "not_assessable"


def test_a_false_and_takes_its_status_from_the_false_inputs_only():
    # surface (R0) is false with status smoke; antigen is true with status estimated
    statuses = {R0: ("smoke", "t"), "antigen_lookup": ("estimated", "t")}
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="5")}))
    res = run(mods, status_of=lambda m, t: statuses.get(m, ("unvalidated", "")))
    r = res[("P", "antigen_candidate_surface", "R0")]
    assert (r.value, r.status) == ("not_called", "smoke")
    assert r.status_basis == "step1_rule@R0:t"
