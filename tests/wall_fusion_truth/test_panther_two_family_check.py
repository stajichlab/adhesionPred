"""The committed record of the PANTHER two-family hand check: tables exist, columns, no sequence lines, both families present."""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/panther_two_family_check"
COLUMNS = {
    "family_totals.tsv": [
        "panther_family",
        "family_name",
        "n_uniprot",
        "n_interpro",
        "domain_of_life",
        "kingdom",
        "n",
    ],
    "subfamily_sizes.tsv": [
        "panther_family",
        "subfamily",
        "subfamily_name",
        "n_all",
        "n_fungal",
        "fungal_share",
    ],
    "fungal_clades_by_subfamily.tsv": ["panther_family", "subfamily", "fine_clade", "n_fungal"],
    "domain_counts.tsv": [
        "panther_family",
        "scope",
        "database",
        "entry",
        "entry_name",
        "n_with_entry",
        "fraction",
    ],
    "reviewed_members.tsv": [
        "panther_family",
        "subfamily",
        "accession",
        "gene",
        "protein_name",
        "organism",
    ],
    "kiss_in_family.tsv": ["panther_family", "subfamily", "accession", "kiss_cluster", "gene"],
    "kiss_mapping.tsv": [
        "kiss_cluster",
        "gene",
        "accession",
        "panther_family",
        "subfamily",
        "pfam_scan_class",
    ],
    "pfam_vs_family_counts.tsv": [
        "panther_family",
        "pfam",
        "n_fungal_in_family",
        "n_fungal_with_pfam",
        "n_fungal_both",
    ],
}


def _rows(name):
    with open(REC / name, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_tables_exist_with_columns():
    for name, cols in COLUMNS.items():
        rows = _rows(name)
        assert rows, name
        assert set(cols) <= set(rows[0]), name


def test_no_sequence_lines():
    seq = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{40,}$")
    files = [f for f in sorted(REC.rglob("*")) if f.suffix in {".tsv", ".md", ".py"}]
    assert len(files) >= 8
    for f in files:
        assert not any(seq.match(ln.strip()) for ln in f.read_text().splitlines()), f.name


def test_both_families_present():
    for name in (
        "family_totals.tsv",
        "subfamily_sizes.tsv",
        "domain_counts.tsv",
        "reviewed_members.tsv",
    ):
        assert {"PTHR46430", "PTHR35192"} <= {r["panther_family"] for r in _rows(name)}, name


def test_known_counts():
    sub = {r["subfamily"]: r for r in _rows("subfamily_sizes.tsv")}
    assert set(sub) >= {"PTHR46430:SF1", "PTHR46430:SF2", "PTHR46430:SF3", "PTHR35192:SF2"}
    assert sub["PTHR35192:SF2"]["n_all"] == sub["PTHR35192:SF2"]["n_fungal"]
    fam = {(r["panther_family"], r["kingdom"]): int(r["n"]) for r in _rows("family_totals.tsv")}
    assert fam[("PTHR46430", "Fungi")] == 4935
    kiss = {r["kiss_cluster"]: r for r in _rows("kiss_mapping.tsv")}
    assert set(kiss) == {"KS0065", "KS0084", "KS0056"}
    assert kiss["KS0065"]["subfamily"] == "PTHR46430:SF2"
    assert kiss["KS0084"]["subfamily"] == "PTHR46430:SF1"
    assert kiss["KS0056"]["panther_family"] == "PTHR35192"
