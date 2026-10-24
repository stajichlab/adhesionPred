"""Record test for docs/reports/data/aga1_prediction_step_2026-10-24/ (real SignalP 6 and PredGPI outputs, 850 proteins).

The table, the raw files and the README hashes are real data; the fixtures of build_record_table below are SYNTHETIC.
"""

import csv
import hashlib
import importlib.util
import re
import subprocess
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/aga1_prediction_step_2026-10-24"
TABLE = REC / "aga1_signalp_predgpi.tsv"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "analysis/aga1_class" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


brt = _load("build_record_table")


@pytest.fixture(scope="module")
def rows():
    with TABLE.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_table_matches_manifest_and_proteins(rows):
    manifest = (ROOT / "analysis/aga1_class/manifest_accessions.txt").read_text().split()
    assert [r["accession"] for r in rows] == manifest and len(rows) == 850
    with (ROOT / "docs/reports/data/aga1_class_scoping/proteins.tsv").open(newline="") as fh:
        prot = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for r in rows:
        assert r["role"] == prot[r["accession"]]["role"]
    assert not any("disorder" in c for c in rows[0])  # the disorder run is not part of this record


def _by_group(rows, col, value):
    n, hit = Counter(), Counter()
    for r in rows:
        n[r["role_group"]] += 1
        hit[r["role_group"]] += r[col] == value
    return n, hit


def test_signalp_counts_by_role(rows):
    n, sp = _by_group(rows, "signal_peptide_pred", "yes")
    assert sum(sp.values()) == 619 and Counter(r["signal_peptide_pred"] for r in rows)["no"] == 231
    assert (n["class_positive_A"], sp["class_positive_A"]) == (35, 35)
    assert (n["hard_negative"], sp["hard_negative"]) == (172, 160)
    assert (n["background_assumed_negative"], sp["background_assumed_negative"]) == (389, 303)
    assert (n["class_positive_B"], sp["class_positive_B"]) == (32, 29)
    assert (n["kiss_row_positive"], sp["kiss_row_positive"]) == (70, 20)
    _, uni = _by_group(rows, "uniprot_signal_feature", "yes")
    assert (uni["class_positive_A"], uni["hard_negative"], uni["background_assumed_negative"]) == (
        35,
        163,
        389,
    )
    assert uni["kiss_row_positive"] == 0 and all(
        r["uniprot_signal_feature"] == "" for r in rows if r["role_group"] == "kiss_row_positive"
    )


def test_predgpi_counts_by_role(rows):
    n, gp = _by_group(rows, "gpi_pred", "yes")
    _, uni = _by_group(rows, "uniprot_gpi_feature", "yes")
    expect = {  # group: (predicted GPI, UniProt GPI)
        "class_positive_A": (27, 24), "class_positive_B": (28, 11), "hard_negative": (53, 35),
        "background_assumed_negative": (99, 95), "adhesin_other_mechanism": (20, 9),
        "excluded_surface_other_adhesion_phenotype": (16, 7),
    }  # fmt: skip
    for g, (p, u) in expect.items():
        assert (gp[g], uni[g]) == (p, u), g
    only = sum(
        r["gpi_pred"] == "yes" and r["uniprot_gpi_feature"] == "no"
        for r in rows
        if r["role_group"] == "class_positive_A"
    )
    assert only == 4
    assert Counter(r["gpi_pred"] for r in rows) == {"no": 597, "yes": 250, "NA": 3}
    assert Counter(r["gpi_confidence_class"] for r in rows) == {
        "NA": 600,
        "high": 232,
        "probable": 8,
        "weak": 10,
    }
    assert n["class_positive_A"] == 35


def test_reference_proteins_and_short_proteins(rows):
    by = {r["accession"]: r for r in rows}
    for acc, omega in (("P32323", "700"), ("Q01589", "307"), ("P47179", "1141")):
        r = by[acc]
        assert (
            r["gpi_pred"],
            r["gpi_omega_site"],
            r["gpi_confidence_class"],
            r["gpi_fpr_upper_bound"],
        ) == ("yes", omega, "high", "0.0015")
    assert by["P32334"]["gpi_pred"] == "no" and by["P32334"]["gpi_omega_site"] == "NA"
    for acc in ("Q3E837", "Q8TGS4", "Q8TGS5"):
        assert by[acc]["gpi_pred"] == "NA" and "length<=40" in by[acc]["gpi_note"]
    assert sum(r["gpi_pred"] == "NA" for r in rows) == 3
    for r in rows:  # no invented false positive rate; bound only with a GPI call
        assert (r["gpi_fpr_upper_bound"] != "NA") == (r["gpi_pred"] == "yes")


def _readme_hashes():
    text = (REC / "README.md").read_text()
    found = re.findall(r"^([0-9a-f]{64})  (\S+)$", text, flags=re.M)
    return {path: sha for sha, path in found}


def test_readme_hashes_match_committed_files():
    h = _readme_hashes()
    assert set(h) == {
        "raw/predgpi/predgpi_raw.gff3", "raw/predgpi/provenance.txt", "raw/signalp/prediction_results.txt",
        "raw/signalp/provenance.txt", "aga1_signalp_predgpi.tsv",
    }  # fmt: skip
    for rel, sha in h.items():
        assert hashlib.sha256((REC / rel).read_bytes()).hexdigest() == sha, rel
    assert (
        "708384fbdbe4a2b24714c8547a12ef23671faff7e1bd75fd5c9cf8475ac63a0e"
        in (REC / "README.md").read_text()
    )
    assert (
        "708384fbdbe4a2b24714c8547a12ef23671faff7e1bd75fd5c9cf8475ac63a0e"
        in (REC / "raw/signalp/provenance.txt").read_text()
    )
    assert "424f5e8861dfd4c0cb56fd7533c5e69610120086" in (REC / "README.md").read_text()


def test_raw_files_are_the_documented_shape():
    g = (REC / "raw/predgpi/predgpi_raw.gff3").read_text().splitlines()
    assert len(g) == 850 and all(len(ln.split("\t")) == 9 for ln in g)
    s = [
        ln
        for ln in (REC / "raw/signalp/prediction_results.txt").read_text().splitlines()
        if not ln.startswith("#")
    ]
    assert len(s) == 850
    assert not any(
        line.startswith(">")
        for p in REC.rglob("*")
        if p.is_file()
        for line in p.read_text().splitlines()
    )  # no FASTA


def test_record_paths_are_not_gitignored():
    files = [p for p in REC.rglob("*") if p.is_file()]
    assert len(files) == 6
    r = subprocess.run(["git", "check-ignore", *map(str, files)], cwd=ROOT, capture_output=True)
    assert r.returncode == 1, (
        "a record file matches .gitignore"
    )  # 1 = none ignored (0 = at least one ignored)
    assert not any(
        part in ("build", "results", "input", "output", "dist")
        for p in files
        for part in p.relative_to(ROOT).parts
    )


# ---- build_record_table, SYNTHETIC fixtures ----
def test_build_record_table_roles_and_failures():
    prot = [
        {"accession": a, "gene": "g", "organism": "o", "clade_coarse": "c", "role": role, "pos_tier": tier, "curated_class": "x",
         "uniprot_gpi_feature": "no", "uniprot_signal_feature": "yes"}
        for a, role, tier in (("S1", "class_positive", "A"), ("S2", "class_positive", "B"), ("S3", "hard_negative", ""))
    ]  # fmt: skip
    pred = [{"accession": a, "signal_peptide_pred": "yes"} for a in ("S1", "S2", "S3")]
    out = brt.build(pred, prot)
    assert [r["role_group"] for r in out] == [
        "class_positive_A",
        "class_positive_B",
        "hard_negative",
    ]
    assert list(out[0])[:5] == ["accession", "gene", "organism", "clade_coarse", "role_group"]
    with pytest.raises(SystemExit, match="same accessions"):
        brt.build(pred[:2], prot)
    with pytest.raises(SystemExit, match="same accessions"):
        brt.build([*pred, {"accession": "S4"}], prot)
    bad = [{**prot[0], "pos_tier": ""}, *prot[1:]]
    with pytest.raises(SystemExit, match="without tier"):
        brt.build(pred, bad)
    with pytest.raises(SystemExit, match="duplicate"):
        brt.build(pred, [*prot, prot[0]])
