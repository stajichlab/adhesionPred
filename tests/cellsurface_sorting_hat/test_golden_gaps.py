"""Golden gap cases (D10): seven proteins that the old rules left without a label or a clear call.

The expected rows in ``golden/gaps/`` were written by hand from the spec (D4, D5, D14, D18) before
the engine was changed. They are not regenerated: a difference is a defect to explain.

  G1 no signal peptide, repeat02 called          G5 signal peptide, antigen_lookup not applicable
  G2 signal peptide, hydrophobin hit only        G6 signal peptide, hit in a family that is not narrow
  G3 signal peptide, HsbA hit only               G7 signal peptide, antigen rank in the top 15 only
  G4 signal peptide, pfam_hsba row in error
"""

import csv
import gzip
from pathlib import Path

import pytest

from cellsurface_sorting_hat.cli import main

GAPS = Path(__file__).parent / "golden" / "gaps"
IDS = [f"G{i}" for i in range(1, 8)]
SEQS = ["MKTAYIAKQRQISFVKSHFSRQ", "MNLLPQWERTYIPASDFG", "MSTNPKPQRKTKRNTNRRPQDV", "MAAGHHKLLWQQ"]


def _fasta():
    return "".join(f">{pid}\n{SEQS[i % 4]}{'ACDEFGHIK'[i]}\n" for i, pid in enumerate(IDS))


def _rows(fn):
    return [fn(pid) for pid in IDS]


def build(tmp_path, write_module):
    fasta = tmp_path / "gaps.faa"
    fasta.write_text(_fasta())
    wd = tmp_path / "wd"

    def ok(pid, **fields):
        return {"id": pid, "state": "ok", **fields}

    def mod(name, fn):
        write_module(wd, name, _rows(fn))

    mod("step1_rule@R0", lambda p: ok(p, call="not_called" if p == "G1" else "called"))
    mod("repeat02", lambda p: ok(p, call="called" if p == "G1" else "not_called"))
    mod("repeat14", lambda p: ok(p, call="not_called"))
    mod("pfam_adhesion", lambda p: ok(p, hit="1" if p == "G6" else "0", hit_narrow="0"))
    mod("pfam_hydrophobin", lambda p: ok(p, hit="1" if p == "G2" else "0"))
    mod(
        "pfam_hsba",
        lambda p: {"id": p, "state": "error"}
        if p == "G4"
        else ok(p, hit="1" if p == "G3" else "0"),
    )
    mod("pfam_allergen", lambda p: ok(p, hit="0"))
    mod(
        "antigen_lookup",
        lambda p: ok(p, percentile="7.1") if p == "G7" else {"id": p, "state": "not_applicable"},
    )
    mod("allergen_homology", lambda p: ok(p, identity="0", coverage="0", aligned_length="100"))
    return fasta, wd


@pytest.fixture
def gap_run(tmp_path, write_module, nodes_dmp):
    fasta, wd = build(tmp_path, write_module)
    out = tmp_path / "out"
    argv = ["--fasta", str(fasta), "--taxon", "40", "--taxdump", str(nodes_dmp)]
    argv += ["--workdir", str(wd), "--out", str(out)]
    assert main(argv) == 0
    return out


def _tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read(out, name):
    with gzip.open(out / name, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _long(out):
    return {(r["protein"], r["call"], r["variant"]): r for r in _read(out, "calls.long.tsv.gz")}


def test_gap_call_values_and_other_basis(gap_run):
    got = _long(gap_run)
    expected = _tsv(GAPS / "expected_calls.tsv")
    assert len(expected) == 77
    for e in expected:
        r = got[(e["protein"], e["call"], e["variant"])]
        assert (r["value"], r["other_basis"]) == (e["value"], e["other_basis"]), e


def test_gap_status_basis_names_the_not_applicable_module(gap_run):
    got = _long(gap_run)
    for e in _tsv(GAPS / "expected_basis.tsv"):
        r = got[(e["protein"], e["call"], e["variant"])]
        assert r["status_basis"] == e["status_basis"], e


def test_gap_primary_labels(gap_run):
    wide = {r["protein"]: r for r in _read(gap_run, "calls.wide.tsv.gz")}
    expected = _tsv(GAPS / "expected_wide.tsv")
    assert [e["protein"] for e in expected] == IDS
    for e in expected:
        r = wide[e["protein"]]
        assert (r["primary_label"], r["primary_basis"]) == (e["primary_label"], e["primary_basis"])


def test_a_secreted_protein_outside_rs_is_not_a_serodiagnostic_candidate(gap_run):
    """G5: antigen_lookup is not applicable. Dropping NA from the operands would call it (FD-A1)."""
    got = _long(gap_run)
    for pid in ("G1", "G2", "G3", "G4", "G5", "G6"):
        assert got[(pid, "serodiagnostic_marker_candidate", "R0")]["value"] != "called"
    assert got[("G7", "serodiagnostic_marker_candidate", "R0")]["value"] == "called"


# Task 4.2: golden extensions. The expected rows in ``golden/extensions/`` were written by hand from
# categories.yaml and the spec before the cases were run. Do not regenerate them from the engine.
#
#   X1 signal peptide, X2 no signal peptide: antigen_lookup, cys_rich and expression all write
#      not_applicable (the three modules that can)
#   X3 signal peptide, two unknown mechanisms          X4 no signal peptide, one called, two unknown
#   X5 no signal peptide, two unknown, none called
#   V1, V2: two step 1 variants (R0, R1) that disagree, for the per_variant calls
EXT = Path(__file__).parent / "golden" / "extensions"
XIDS = ["X1", "X2", "X3", "X4", "X5"]
VIDS = ["V1", "V2"]
NA_MODULES = ("antigen_lookup", "cys_rich", "expression")


def _write_fasta(path, ids):
    path.write_text(
        "".join(f">{pid}\n{SEQS[i % 4]}{'ACDEFGHIK'[i]}\n" for i, pid in enumerate(ids))
    )


def _na_row(pid, module):
    """A not_applicable row; the table still has the columns of the module (empty values)."""
    cols = {"antigen_lookup": ["percentile"], "cys_rich": ["tier"], "expression": ["log2fc"]}
    return {"id": pid, "state": "not_applicable", **dict.fromkeys(cols[module], "")}


def _run_cli(tmp_path, fasta, wd, nodes_dmp):
    out = tmp_path / "out"
    argv = ["--fasta", str(fasta), "--taxon", "40", "--taxdump", str(nodes_dmp)]
    argv += ["--workdir", str(wd), "--out", str(out)]
    assert main(argv) == 0
    return out


@pytest.fixture
def na_run(tmp_path, write_module, nodes_dmp):
    fasta = tmp_path / "x.faa"
    _write_fasta(fasta, XIDS)
    wd = tmp_path / "wd"
    unknown = {"X3", "X4", "X5"}

    def rows(fn):
        return [fn(p) for p in XIDS]

    def ok(pid, **fields):
        return {"id": pid, "state": "ok", **fields}

    signal = {"X1": "called", "X3": "called"}
    write_module(wd, "step1_rule@R0", rows(lambda p: ok(p, call=signal.get(p, "not_called"))))
    write_module(
        wd, "repeat02", rows(lambda p: ok(p, call="called" if p == "X4" else "not_called"))
    )
    write_module(wd, "repeat14", rows(lambda p: ok(p, call="not_called")))
    write_module(wd, "pfam_adhesion", rows(lambda p: ok(p, hit="0", hit_narrow="0")))
    for mod in ("pfam_hydrophobin", "pfam_hsba"):
        write_module(
            wd,
            mod,
            rows(lambda p: {"id": p, "state": "error"} if p in unknown else ok(p, hit="0")),
        )
    write_module(wd, "pfam_allergen", rows(lambda p: ok(p, hit="0")))
    write_module(
        wd,
        "allergen_homology",
        rows(lambda p: ok(p, identity="0", coverage="0", aligned_length="1")),
    )
    for mod in NA_MODULES:
        write_module(wd, mod, rows(lambda p, mod=mod: _na_row(p, mod)))
    return _run_cli(tmp_path, fasta, wd, nodes_dmp)


@pytest.fixture
def variant_run(tmp_path, write_module, nodes_dmp):
    fasta = tmp_path / "v.faa"
    _write_fasta(fasta, VIDS)
    wd = tmp_path / "wd"

    def rows(fn):
        return [fn(p) for p in VIDS]

    def ok(pid, **fields):
        return {"id": pid, "state": "ok", **fields}

    write_module(
        wd, "step1_rule@R0", rows(lambda p: ok(p, call="called" if p == "V1" else "not_called"))
    )
    write_module(
        wd, "step1_rule@R1", rows(lambda p: ok(p, call="called" if p == "V2" else "not_called"))
    )
    write_module(wd, "repeat02", rows(lambda p: ok(p, call="not_called")))
    write_module(wd, "repeat14", rows(lambda p: ok(p, call="not_called")))
    write_module(wd, "pfam_adhesion", rows(lambda p: ok(p, hit="1", hit_narrow="1")))
    write_module(wd, "pfam_hydrophobin", rows(lambda p: ok(p, hit="0")))
    write_module(wd, "pfam_hsba", rows(lambda p: ok(p, hit="0")))
    write_module(wd, "pfam_allergen", rows(lambda p: ok(p, hit="0")))
    write_module(
        wd,
        "allergen_homology",
        rows(lambda p: ok(p, identity="0", coverage="0", aligned_length="1")),
    )
    write_module(wd, "antigen_lookup", rows(lambda p: _na_row(p, "antigen_lookup")))
    return _run_cli(tmp_path, fasta, wd, nodes_dmp)


def _check_calls(out, expected_file, n_rows):
    got = _long(out)
    expected = _tsv(EXT / expected_file)
    assert len(expected) == n_rows
    for e in expected:
        r = got[(e["protein"], e["call"], e["variant"])]
        assert (r["value"], r["other_basis"]) == (e["value"], e["other_basis"]), e


def test_na_at_every_module_gives_the_hand_written_calls(na_run):
    _check_calls(na_run, "expected_calls.tsv", 55)


def test_na_at_every_module_gives_the_hand_written_labels(na_run):
    wide = {r["protein"]: r for r in _read(na_run, "calls.wide.tsv.gz")}
    for e in [r for r in _tsv(EXT / "expected_wide.tsv") if r["protein"] in XIDS]:
        got = wide[e["protein"]]
        assert (got["primary_label"], got["primary_basis"]) == (
            e["primary_label"],
            e["primary_basis"],
        ), e


def test_na_modules_that_no_call_reads_add_no_evidence(na_run):
    """cys_rich and expression write not_applicable; they are evidence only, and an NA row has none."""
    evidence = _read(na_run, "evidence.tsv.gz")
    assert [r for r in evidence if r["field"].split(".")[0] in NA_MODULES] == []


def test_several_unknown_mechanisms_are_all_named(na_run):
    got = _long(na_run)
    for pid in ("X3", "X5"):
        call = "other_surface_no_mechanism" if pid == "X3" else "other_not_surface"
        r = got[(pid, call, "R0")]
        assert r["value"] == "not_assessable"
        assert r["other_basis"] == "hydrophobin_domain,hsba_domain"
    # one called mechanism outranks the unknown ones (X4)
    assert got[("X4", "other_not_surface", "R0")]["value"] == "not_called"


def test_per_variant_calls_follow_their_own_step1_variant(variant_run):
    _check_calls(variant_run, "expected_variants.tsv", 36)


def test_the_primary_label_reads_the_default_gate_only(variant_run):
    wide = {r["protein"]: r for r in _read(variant_run, "calls.wide.tsv.gz")}
    for e in [r for r in _tsv(EXT / "expected_wide.tsv") if r["protein"] in VIDS]:
        got = wide[e["protein"]]
        assert (got["primary_label"], got["primary_basis"]) == (
            e["primary_label"],
            e["primary_basis"],
        ), e
