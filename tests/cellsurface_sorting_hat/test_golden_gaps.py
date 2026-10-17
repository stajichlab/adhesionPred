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


@pytest.mark.xfail(strict=True, reason="Tasks 3.2 and 3.4: not_applicable, narrow call")
def test_gap_call_values_and_other_basis(gap_run):
    got = _long(gap_run)
    expected = _tsv(GAPS / "expected_calls.tsv")
    assert len(expected) == 77
    for e in expected:
        r = got[(e["protein"], e["call"], e["variant"])]
        assert (r["value"], r["other_basis"]) == (e["value"], e["other_basis"]), e


@pytest.mark.xfail(strict=True, reason="Task 3.2: not_applicable basis")
def test_gap_status_basis_names_the_not_applicable_module(gap_run):
    got = _long(gap_run)
    for e in _tsv(GAPS / "expected_basis.tsv"):
        r = got[(e["protein"], e["call"], e["variant"])]
        assert r["status_basis"] == e["status_basis"], e


@pytest.mark.xfail(strict=True, reason="Task 3.7: primary_label")
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
