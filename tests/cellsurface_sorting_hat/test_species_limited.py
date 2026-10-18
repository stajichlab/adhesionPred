"""Issue #93: a module that is limited to some taxa (``species_limited`` in categories.yaml).

An absent module table gives ``not_applicable`` for a protein outside the taxa of the module, and
stays ``not_assessable`` inside them. A present table decides by its rows, as before.
Tests that fail before the change carry ``xfail(strict=True)``; the commit that lands the
behaviour removes the mark.
"""

import copy
import gzip
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from engine_helpers import base_modules, ok, table
from engine_helpers import run_eval as _run_eval
from test_cli import build, read_long

import cellsurface_sorting_hat.engine as engine_module
from cellsurface_sorting_hat.cli import main
from cellsurface_sorting_hat.engine import ConfigError, evaluate, load_config
from cellsurface_sorting_hat.modules import lookups

RS = 246410
OUTSIDE = 40
NA = "not_applicable"
NOT_ASSESSABLE = "not_assessable"
RANK = ("P", "cocci_specificity_rank_top15", "")


def run_eval(modules, **kw):
    """``run_eval`` with a config that limits antigen_lookup to RS, whatever the packaged file says."""
    if "cfg" not in kw:
        with tempfile.TemporaryDirectory() as folder:
            kw["cfg"] = _config(Path(folder), {"antigen_lookup": [RS]})
    return _run_eval(modules, **kw)


def no_antigen(step1="called", **over):
    mods = base_modules(step1=step1, **over)
    del mods["antigen_lookup"]
    return mods


def test_absent_module_outside_the_taxa_is_not_applicable():
    got = run_eval(no_antigen(), taxon=OUTSIDE)
    assert got[RANK].value == NA
    assert got[RANK].status_basis == "not_applicable: antigen_lookup"
    assert got[("P", "other_surface_no_mechanism", "R0")].value == "called"
    assert got[("P", "other_surface_no_mechanism", "R0")].other_basis == ""
    assert got[("P", "other_not_surface", "R0")].value == "not_called"


def test_absent_module_outside_the_taxa_gives_a_label_for_a_not_surface_protein():
    got = run_eval(no_antigen(step1="not_called"), taxon=OUTSIDE)
    assert got[("P", "other_not_surface", "R0")].value == "called"
    assert got[("P", "other_surface_no_mechanism", "R0")].value == "not_called"


def test_absent_module_inside_the_taxa_stays_not_assessable():
    got = run_eval(no_antigen(), taxon=RS)
    assert got[RANK].value == NOT_ASSESSABLE
    assert got[("P", "other_surface_no_mechanism", "R0")].value == NOT_ASSESSABLE
    assert got[("P", "other_surface_no_mechanism", "R0")].other_basis == (
        "cocci_specificity_rank_top15"
    )


@pytest.mark.parametrize("taxon", [RS, OUTSIDE])
def test_a_present_table_decides_by_its_rows(taxon):
    top = base_modules(
        step1="called", antigen_lookup=table("antigen_lookup", {"P": ok(percentile="5")})
    )
    assert run_eval(top, taxon=taxon)[RANK].value == "called"
    low = base_modules(step1="called")
    assert run_eval(low, taxon=taxon)[RANK].value == "not_called"
    na_row = base_modules(
        step1="called", antigen_lookup=table("antigen_lookup", {"P": {"state": NA}})
    )
    assert run_eval(na_row, taxon=taxon)[RANK].value == NA
    error = base_modules(
        step1="called", antigen_lookup=table("antigen_lookup", {"P": {"state": "error"}})
    )
    assert run_eval(error, taxon=taxon)[RANK].value == NOT_ASSESSABLE


@pytest.mark.parametrize(
    "step1",
    ["called", "not_called"],
)
def test_no_serodiagnostic_call_for_a_secreted_protein_outside_the_taxa(step1):
    """FD-A1: the absent module must not turn into a mass call."""
    got = run_eval(no_antigen(step1=step1), taxon=OUTSIDE)
    assert got[("P", "serodiagnostic_marker_candidate", "R0")].value == "not_called"


def test_a_not_called_record_with_a_not_applicable_leaf_has_a_basis():
    got = run_eval(no_antigen(), taxon=OUTSIDE)
    sero = got[("P", "serodiagnostic_marker_candidate", "R0")]
    assert sero.value == "not_called"
    assert sero.status_basis == "not_applicable: antigen_lookup"
    assert sero.status == "unvalidated"


def test_matching_is_exact_not_by_lineage():
    # 41 is a toy taxon below the Onygenales; only the exact taxon of the list counts
    assert run_eval(no_antigen(), taxon=41)[RANK].value == NA


# ---- through the CLI, with a taxon map -----------------------------------------------------------


def _nodes(tmp_path):
    from conftest import TOY_NODES

    nodes = dict(TOY_NODES)
    nodes[RS] = 31
    path = tmp_path / "nodes_rs.dmp"
    path.write_text("".join(f"{t}\t|\t{p}\t|\tno rank\t|\n" for t, p in nodes.items()))
    return path


def test_mixed_taxa_get_a_decision_per_protein(tmp_path, write_module):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    for suffix in (".tsv.gz", ".json"):
        (wd / "modules" / f"antigen_lookup{suffix}").unlink()
    taxon_map.write_text("SOW1\t246410\nDUP1\t40\nENZ1\t40\nBAD1\t40\nSTAR1\t246410\n")
    out = tmp_path / "out"
    argv = [
        "--fasta",
        str(fasta),
        "--taxon-map",
        str(taxon_map),
        "--taxdump",
        str(_nodes(tmp_path)),
    ]
    argv += ["--workdir", str(wd), "--out", str(out)]
    assert main(argv) == 0
    with gzip.open(out / "calls.long.tsv.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"), strict=False)) for line in fh]
    value = {(r["protein"], r["call"], r["variant"]): r["value"] for r in rows}
    assert value[("SOW1", "cocci_specificity_rank_top15", "")] == NOT_ASSESSABLE
    assert value[("STAR1", "cocci_specificity_rank_top15", "")] == NOT_ASSESSABLE
    assert value[("DUP1", "cocci_specificity_rank_top15", "")] == NA
    assert value[("ENZ1", "cocci_specificity_rank_top15", "")] == NA
    assert value[("DUP1", "serodiagnostic_marker_candidate", "R0")] == "not_called"


# ---- config validation ---------------------------------------------------------------------------


def _config(tmp_path, species_limited="absent"):
    data = yaml.safe_load(
        Path(engine_module.__file__).with_name("categories.yaml").read_text(encoding="utf-8")
    )
    data = copy.deepcopy(data)
    if species_limited != "absent":
        data["species_limited"] = species_limited
    else:
        data.pop("species_limited", None)
    path = tmp_path / "categories.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return load_config(path)


def test_the_packaged_config_limits_the_antigen_module_to_rs():
    assert load_config().species_limited == {"antigen_lookup": [RS]}


def test_the_key_is_optional(tmp_path):
    assert _config(tmp_path).species_limited == {}
    assert _config(tmp_path, {"antigen_lookup": [RS, 5501]}).species_limited == {
        "antigen_lookup": [RS, 5501]
    }


@pytest.mark.parametrize(
    ("value", "message"),
    [
        ({"no_such_module": [RS]}, "no_such_module"),
        ({"antigen_lookup": []}, "non-empty list"),
        ({"antigen_lookup": ["246410"]}, "integer"),
        ({"antigen_lookup": [1.5]}, "integer"),
        ({"antigen_lookup": [True]}, "integer"),
        ({"antigen_lookup": RS}, "non-empty list"),
        ([RS], "mapping"),
        ("antigen_lookup", "mapping"),
        ({"antigen_lookup": [0]}, "at least 2"),
        ({"antigen_lookup": [-5]}, "at least 2"),
        ({"antigen_lookup": [1]}, "at least 2"),
        ({"antigen_lookup": [RS, RS]}, "duplicate"),
        ({"step1_rule@R0": [RS]}, "step 1"),
    ],
)
def test_a_bad_species_limited_block_is_a_config_error(tmp_path, value, message):
    with pytest.raises(ConfigError, match=message):
        _config(tmp_path, value)


# ---- invalid proteins: absent == run, also for a protein that is not valid ---------------------


def _values_of(out, protein):
    return {k[1:]: r["value"] for k, r in read_long(out).items() if k[0] == protein}


def _cli(tmp_path, write_module, nodes, run_antigen):
    tmp_path.mkdir(exist_ok=True)
    fasta, taxon_map, wd = build(tmp_path, write_module)
    if not run_antigen:
        for suffix in (".tsv.gz", ".json"):
            (wd / "modules" / f"antigen_lookup{suffix}").unlink()
    out = tmp_path / "out"
    argv = ["--fasta", str(fasta), "--taxon-map", str(taxon_map), "--taxdump", str(nodes)]
    argv += ["--workdir", str(wd), "--out", str(out)]
    assert main(argv) == 0
    return out


@pytest.mark.xfail(strict=True, reason="issue #93 review: invalid proteins")
def test_an_invalid_protein_gives_the_same_leaf_results_absent_and_run(tmp_path, write_module):
    nodes = _nodes(tmp_path)
    absent = _cli(tmp_path / "a", write_module, nodes, run_antigen=False)
    run = _cli(tmp_path / "b", write_module, nodes, run_antigen=True)
    assert _values_of(absent, "BAD1") == _values_of(run, "BAD1")
    assert _values_of(absent, "BAD1")[("serodiagnostic_marker_candidate", "R0")] == NOT_ASSESSABLE
    # a valid protein outside the list is still not_applicable when the table is absent
    assert _values_of(absent, "ENZ1")[("cocci_specificity_rank_top15", "")] == NA


@pytest.mark.xfail(strict=True, reason="issue #93 review: invalid proteins")
def test_evaluate_skips_the_rule_for_the_invalid_ids():
    mods = no_antigen()
    cfg = load_config()
    kw = {"invalid_ids": frozenset({"P"})}
    records = evaluate(cfg, ["P"], {"P": OUTSIDE}, mods, lambda m, t: ("unvalidated", "x"), **kw)
    got = {(r.call, r.variant): r for r in records}
    assert got[("cocci_specificity_rank_top15", "")].value == NOT_ASSESSABLE
    assert got[("serodiagnostic_marker_candidate", "R0")].value == NOT_ASSESSABLE
    valid = evaluate(cfg, ["P"], {"P": OUTSIDE}, mods, lambda m, t: ("unvalidated", "x"))
    assert {(r.call, r.variant): r for r in valid}[("cocci_specificity_rank_top15", "")].value == NA


# ---- two sources of truth for the applicable taxa ------------------------------------------------


def _mismatch_run(tmp_path, write_module, module_taxa):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    meta_path = wd / "modules" / "antigen_lookup.json"
    meta = json.loads(meta_path.read_text())
    if module_taxa is not None:
        meta["params"] = {"applicable_taxa": module_taxa}
    meta_path.write_text(json.dumps(meta))
    out = tmp_path / "out"
    argv = ["--fasta", str(fasta), "--taxon-map", str(taxon_map), "--taxdump"]
    argv += [str(_nodes(tmp_path)), "--workdir", str(wd), "--out", str(out)]
    assert main(argv) == 0
    return out


@pytest.mark.xfail(strict=True, reason="issue #93 review: two lists")
def test_a_module_list_that_differs_from_the_yaml_list_is_a_warning(tmp_path, write_module):
    out = _mismatch_run(tmp_path, write_module, [RS, 5501])
    report = (out / "report.md").read_text()
    assert "antigen_lookup" in report
    assert "[246410, 5501]" in report
    assert "[246410]" in report
    info = json.loads((out / "run.json").read_text())
    assert info["applicable_taxa_mismatch"] == [
        {"module": "antigen_lookup", "module_taxa": [RS, 5501], "config_taxa": [RS]}
    ]


@pytest.mark.parametrize("module_taxa", [None, [RS]])
def test_no_warning_when_the_lists_agree_or_the_record_has_none(
    tmp_path, write_module, module_taxa
):
    out = _mismatch_run(tmp_path, write_module, module_taxa)
    assert "applicable taxa" not in (out / "report.md").read_text()
    assert not json.loads((out / "run.json").read_text()).get("applicable_taxa_mismatch")


# ---- a child taxon of the listed taxon is not in the list ----------------------------------------

CHILD = 999999


def _child_nodes(tmp_path):
    from conftest import TOY_NODES

    nodes = dict(TOY_NODES)
    nodes[RS] = 31
    nodes[CHILD] = RS
    path = tmp_path / "nodes_child.dmp"
    path.write_text("".join(f"{t}\t|\t{p}\t|\tno rank\t|\n" for t, p in nodes.items()))
    return path


def _child_run(tmp_path, write_module, run_antigen):
    tmp_path.mkdir(exist_ok=True)
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text("".join(f"{pid}\t{CHILD}\n" for pid in TAXA_IDS))
    for suffix in (".tsv.gz", ".json"):
        (wd / "modules" / f"antigen_lookup{suffix}").unlink()
    if run_antigen:
        proteins = [SimpleNamespace(id=pid, state="ok") for pid in VALID_IDS]
        rows = lookups.antigen_rows(proteins, dict.fromkeys(VALID_IDS, CHILD), {}, {}, {RS})
        rows = [{"percentile": "", **r} for r in rows]  # the module writes every column
        write_module(wd, "antigen_lookup", rows, meta={"params": {"applicable_taxa": [RS]}})
    out = tmp_path / "out"
    argv = ["--fasta", str(fasta), "--taxon-map", str(taxon_map), "--taxdump"]
    argv += [str(_child_nodes(tmp_path)), "--workdir", str(wd), "--out", str(out)]
    assert main(argv) == 0
    return out


TAXA_IDS = ["SOW1", "ENZ1", "DUP1", "BAD1", "STAR1"]
VALID_IDS = ["SOW1", "ENZ1", "DUP1", "STAR1"]


def test_a_child_taxon_of_the_listed_taxon_is_not_in_the_list(tmp_path, write_module):
    absent = _child_run(tmp_path / "a", write_module, run_antigen=False)
    run = _child_run(tmp_path / "b", write_module, run_antigen=True)
    for pid in VALID_IDS:
        a, r = _values_of(absent, pid), _values_of(run, pid)
        assert a == r, pid
        assert a[("cocci_specificity_rank_top15", "")] == NA
