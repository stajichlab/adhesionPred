"""The committed run record of the PR 6 family-specificity runs (jobs 30423745 and 30426504)."""

import csv
import hashlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_family_specificity_results_2026-10-24"
REPORT = ROOT / "docs/reports/2026-10-24-pr6-family-specificity-results.md"
MEMBERS = ROOT / "docs/reports/data/kiss_cluster_split_2026-10-22/members_split.tsv"
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
    assert len(lines) == 18
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


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _helpers():
    sys.path.insert(0, str(SCRIPT.parent))
    try:
        return _load(SCRIPT, "results_summary"), importlib.import_module("pr6_kiss_proteome_map")
    finally:
        sys.path.remove(str(SCRIPT.parent))


def _summary():
    return {(r["run"], r["pfam_acc"]): r for r in _tsv(REC / "summary.tsv")}


def test_every_kiss_positive_is_mapped_or_listed_unmapped():
    members = [r for r in _tsv(MEMBERS) if r["class"] == "positive"]
    rows = _tsv(REC / "kiss_positive_proteome_map.tsv")
    assert len(members) == 132
    assert [r["accession"] for r in rows] == [m["accession"] for m in members]
    assert {r["status"] for r in rows} <= {"mapped", "not_in_kept_proteomes", "unmapped"}
    for r in rows:
        if r["status"] == "mapped":
            assert r["proteome"] and r["proteome_id"], r
        if r["status"] == "unmapped":
            assert r["basis"], r  # listed explicitly, never silently dropped
    by_status = {
        s: sum(r["status"] == s for r in rows)
        for s in ("mapped", "not_in_kept_proteomes", "unmapped")
    }
    assert by_status == {"mapped": 85, "not_in_kept_proteomes": 47, "unmapped": 0}
    by_proteome = {p: sum(r["proteome"] == p for r in rows) for p in {r["proteome"] for r in rows}}
    assert by_proteome == {"": 47, "Afum_Af293_UniProt": 26, "Calb_SC5314": 35, "Scer_S288C": 24}
    assert not [r for r in rows if r["proteome"] == "Cimm_RS"]  # no Coccidioides positive
    assert len({(r["proteome"], r["proteome_id"]) for r in rows if r["status"] == "mapped"}) == 85


def test_map_is_reproducible_from_the_committed_uniprot_export():
    _, pm = _helpers()
    spec = importlib.util.spec_from_file_location(
        "build_map", SCRIPT.parent / "pr6_build_kiss_positive_proteome_map.py"
    )
    sys.path.insert(0, str(SCRIPT.parent))
    try:
        build = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(build)
    finally:
        sys.path.remove(str(SCRIPT.parent))
    rows = build.build(REC, MEMBERS)
    assert [{k: str(v) for k, v in r.items()} for r in rows] == _tsv(
        REC / "kiss_positive_proteome_map.tsv"
    )
    assert len(_tsv(REC / "kiss_positives_uniprot_2026_03.tsv")) == 112


def test_id_rules_and_named_cases():
    _, pm = _helpers()
    assert pm.calb_a22_id("CAALFM_C206100WA") == "C2_06100W_A"
    assert pm.calb_a22_id("CAALFM_CR00610WA orf19.1") == "CR_00610W_A"
    assert pm.calb_a22_id("CaO19.7472") is None
    assert pm.scer_systematic_id("YNL238W") == "YNL238W"
    assert pm.scer_systematic_id("N1122") is None
    assert pm.normalize_proteome_id("Cimm_RS", "XP_001239171.1") == pm.normalize_proteome_id(
        "Cimm_RS", "XP_001239171.2"
    )
    assert pm.normalize_proteome_id("Afum_Af293_UniProt", "tr|Q4WGG2|Q4WGG2_ASPFU") == "Q4WGG2"
    by = {r["accession"]: r for r in _tsv(REC / "kiss_positive_proteome_map.tsv")}
    assert by["Q5AAL9"]["proteome_id"] == "CR_00610W_A"  # IFF4
    assert by["P13134"]["proteome_id"] == "YNL238W"  # KEX2
    assert by["A0A1D8PEG3"]["proteome_id"] == "C1_08990C_A"
    assert by["Q59X23"]["proteome_id"] == "C2_06100W_A"
    assert by["P29465"]["proteome_id"] == "YBR023C"  # CHS3
    assert by["jgi:Neucr2:9773"]["status"] == "not_in_kept_proteomes"
    assert all(
        by[a]["status"] == "not_in_kept_proteomes"
        for a in by
        if by[a]["organism_id"] in {"2059318", "5270"}
    )


def test_summary_table_is_deterministic_and_counts_own_members_by_proteome_id():
    mod, _ = _helpers()
    rows = mod.summarize(REC)
    assert len(rows) == 33
    committed = _tsv(REC / "summary.tsv")
    assert [{k: str(v) for k, v in r.items()} for r in rows] == committed
    assert list(committed[0]) == mod.COLUMNS
    s = _summary()
    pmt = s[("run1", "PF02366")]
    assert (pmt["own_in_proteomes"], pmt["kept_hits_excl_own"], pmt["kept_per10k_excl_own"]) == (
        "5",
        "13",
        "4.00",
    )
    assert s[("run1", "PF16192")]["own_in_proteomes"] == "5"
    s8 = s[("run1", "PF00082")]
    assert (s8["own_in_proteomes"], s8["kept_hits_excl_own"], s8["kept_per10k_excl_own"]) == (
        "3",
        "34",
        "10.46",
    )
    assert s[("run1", "PF11765")]["own_in_proteomes"] == "1"  # IFF4 is one of its 14 hits
    big1 = s[("run2", "PF27678")]
    assert (big1["own_in_proteomes"], big1["kept_hits_excl_own"], big1["non_kiss_hits"]) == (
        "1",
        "0",
        "0",
    )
    g1 = s[("run1", "PF00534")]
    assert (g1["kept_hits_total"], g1["own_in_proteomes"], g1["kept_hits_excl_own"]) == (
        "23",
        "0",
        "23",
    )
    assert g1["kept_per10k_excl_own"] == "7.08"
    for r in committed:
        assert int(r["own_in_proteomes"]) + int(r["other_kiss_positive_in_proteomes"]) + int(
            r["non_kiss_hits"]
        ) == int(r["kept_hits_total"])
        assert int(r["own_in_proteomes"]) <= int(r["n_own_members"])
        assert r["own_mapped_not_hit"] == "0"
        assert (r["over_limit_excl_own"] == "yes") == (float(r["kept_per10k_excl_own"]) > 5)


def test_other_kiss_positive_hits_and_the_corrected_families_above_five():
    s = _summary()
    others = {k[1]: int(r["other_kiss_positive_in_proteomes"]) for k, r in s.items()}
    assert {a: n for a, n in others.items() if n} == {"PF13641": 4}
    gt = s[("run1", "PF13641")]
    assert (gt["kept_hits_excl_own"], gt["kept_per10k_excl_own"]) == ("20", "6.16")
    assert (gt["non_kiss_hits"], gt["non_kiss_per10k"]) == ("16", "4.92")
    above = sorted(a for (run, a), r in s.items() if r["over_limit_excl_own"] == "yes")
    assert above == ["PF00082", "PF00534", "PF13641"]


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
        "Correction",
        "CR_00610W_A",
        "capped at `smoke`",
        "NOT established",
        "owner-signed",
    ):
        assert needle in text, needle
