"""The committed run record of the PR 6 family-specificity runs (jobs 30423745 and 30426504)."""

import csv
import hashlib
import importlib.util
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_family_specificity_results_2026-10-24"
REPORT = ROOT / "docs/reports/2026-10-24-pr6-family-specificity-results.md"
SCRIPT = ROOT / "analysis/wall_fusion_truth/pr6_family_specificity_results_summary.py"


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _fam(run):
    return {r["pfam_acc"]: r for r in _tsv(REC / run / "family_specificity.tsv")}


def test_folder_and_report_exist_and_avoid_ignored_names():
    assert REC.is_dir() and REPORT.is_file()
    assert (REC / "README.md").is_file()
    assert not any(f.suffix == ".csv" or f.suffix == ".log" for f in REC.rglob("*"))
    assert not list(REC.rglob("domtbl_*"))
    assert REC.name not in {"build", "results", "input", "output", "dist"}


def test_no_committed_file_holds_a_sequence():
    seq = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{40,}")
    for f in sorted(REC.rglob("*")):
        if f.is_file():
            for line in f.read_text().splitlines():
                assert not seq.search(line), f
    assert not list(REC.rglob("*.faa")) and not list(REC.rglob("*.hmm"))


def test_row_counts_and_circular_recall_and_leakage():
    assert len(_fam("run1")) == 30
    assert len(_fam("run2")) == 3
    for run in ("run1", "run2"):
        for r in _fam(run).values():
            assert r["own_recall_circular"] == "1.000" and r["own_hits"] == r["n_own_members"]
            assert r["leakage"] == "tuned_on_truth"
            assert (r["n_neg_n1"], r["n_neg_n2"], r["kept_n_proteins"]) == ("137", "24", "32491")


def test_sha256_lists_match_committed_files():
    lines = (REC / "SHA256SUMS.txt").read_text().splitlines()
    assert len(lines) == 15
    for line in lines:
        digest, name = line.split("  ", 1)
        assert _sha(REC / name) == digest, name
    hpcc = dict(
        reversed(x.split("  ", 1)) for x in (REC / "hpcc/hpcc_sha256.txt").read_text().splitlines()
    )
    for run, d in (("run1", "family_specificity"), ("run2", "family_specificity_extras_out")):
        for f in (
            "family_specificity",
            "negatives_hits",
            "positives_hits",
            "proteome_hits",
            "provenance",
        ):
            assert _sha(REC / run / f"{f}.tsv") == hpcc[f"{d}/{f}.tsv"], (run, f)
    assert (
        _sha(REC / "hpcc/pr6_family_specificity.30423745.log.txt")
        == hpcc["logs/pr6_family_specificity.30423745.log"]
    )
    assert (
        _sha(REC / "hpcc/pr6_family_specificity_extras.30426504.log.txt")
        == hpcc["logs/pr6_family_specificity_extras.30426504.log"]
    )
    readme = (REC / "README.md").read_text()
    assert len(hpcc) == 12
    for digest in hpcc.values():
        assert digest in readme, digest


def test_logs_pins_and_inputs():
    one = (REC / "hpcc/pr6_family_specificity.30423745.log.txt").read_text()
    two = (REC / "hpcc/pr6_family_specificity_extras.30426504.log.txt").read_text()
    assert "pin ok: HEAD e36c6e4a02031f30e1c64a74f581ee0f7389e8d7" in one
    assert "pin ok: HEAD f56277354f174f0a852f86cf7fc61234c84995e0" in two
    assert "30 families, 132 positives, 161 negatives" in one
    assert "3 families, 132 positives, 161 negatives" in two
    assert "30 fetched, 30 with a GA line, 0 without" in one
    assert "3 fetched, 3 with a GA line, 0 without" in two
    for run in ("run1", "run2"):
        prov = (REC / run / "provenance.tsv").read_text()
        assert "kept_total\t32491" in prov
        assert "tier_changes\tnone" in prov and "panel_changes\tnone" in prov
    assert "tier_cap\tsmoke" in (REC / "run2/provenance.tsv").read_text()
    assert "pfam_release: 38.2" in (REC / "run1/provenance.txt").read_text()


def test_run1_negatives_and_headline_rates():
    fam = _fam("run1")
    neg = _tsv(REC / "run1/negatives_hits.tsv")
    assert sorted((r["pfam_acc"], r["accession"]) for r in neg) == [
        ("PF00082", "J9VDC2"),
        ("PF02055", "Q4WBR2"),
        ("PF02055", "Q8J0I9"),
    ]
    assert sum(int(r["neg_hits_n1"]) + int(r["neg_hits_n2"]) for r in fam.values()) == 3
    above = {
        a: (int(r["kept_hits_total"]), r["kept_per10k"])
        for a, r in fam.items()
        if float(r["kept_per10k"]) > 5
    }
    assert above == {
        "PF00082": (37, "11.39"),
        "PF00534": (23, "7.08"),
        "PF13641": (20, "6.16"),
        "PF02366": (18, "5.54"),
        "PF00982": (17, "5.23"),
        "PF16192": (17, "5.23"),
    }
    assert fam["PF13641"]["other_positive_hits_wall"] == "6"
    assert fam["PF02366"]["other_positive_hits_wall"] == "1"
    hy = fam["PF11765"]
    assert hy["in_family_table"] == "yes" and hy["kept_hits_total"] == "14"
    assert hy["Calb_SC5314_hits"] == "12"
    rep1 = fam["PF29797"]
    assert rep1["kept_hits_total"] == "0" and rep1["own_hits"] == "1"
    assert rep1["other_positive_hits_wall"] == rep1["other_positive_hits_adhesion"] == "0"


def test_run2_families_are_under_five_per_10000():
    fam = _fam("run2")
    assert {a: (r["name"], r["kept_hits_total"], r["kept_per10k"]) for a, r in fam.items()} == {
        "PF27678": ("BIG1_C", "1", "0.31"),
        "PF09346": ("SMI1_KNR4", "5", "1.54"),
        "PF02544": ("Steroid_dh", "10", "3.08"),
    }
    for r in fam.values():
        assert float(r["kept_per10k"]) < 5
        assert r["neg_hits_n1"] == r["neg_hits_n2"] == "0"
        assert r["other_positive_hits_adhesion"] == r["other_positive_hits_wall"] == "0"
        assert r["other_positive_hits_fusion"] == "0"
    assert _tsv(REC / "run2/negatives_hits.tsv") == []


def test_summary_table_is_deterministic_and_excludes_own_members():
    spec = importlib.util.spec_from_file_location("results_summary", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    rows = mod.summarize(REC)
    assert len(rows) == 33
    committed = _tsv(REC / "summary.tsv")
    assert [{k: str(v) for k, v in r.items()} for r in rows] == committed
    by = {r["pfam_acc"]: r for r in committed}
    assert (by["PF02366"]["own_in_proteomes"], by["PF02366"]["kept_per10k_excl_own"]) == (
        "2",
        "4.92",
    )
    assert by["PF00082"]["kept_per10k_excl_own"] == "11.08"
    assert all(r["own_in_proteomes"] == "0" for r in committed if r["run"] == "run2")


def test_family_table_and_family_measure_unchanged_from_origin_main():
    for rel in ("data/sorting_hat/family_table.tsv", "data/sorting_hat/family_measure.tsv"):
        head = subprocess.run(
            ["git", "-C", str(ROOT), "diff", "--quiet", "origin/main", "--", rel],
            capture_output=True,
        )
        if head.returncode > 1:  # no origin/main in this checkout (for example a detached clone)
            continue
        assert head.returncode == 0, rel


def test_report_states_the_main_numbers_and_limits():
    text = REPORT.read_text()
    for needle in (
        "11.39",
        "7.08",
        "6.16",
        "5.54",
        "5.23",
        "J9VDC2",
        "Q4WBR2",
        "Q8J0I9",
        "0.31",
        "1.54",
        "3.08",
        "circular",
        "NOT established",
        "owner-signed",
    ):
        assert needle in text, needle
