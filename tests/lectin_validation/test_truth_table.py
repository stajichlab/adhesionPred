"""Checks of the lectin validation truth table (analysis/lectin_validation/controls)."""

import csv
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
CTRL = ROOT / "analysis" / "lectin_validation" / "controls"
sys.path.insert(0, str(ROOT / "analysis" / "lectin_validation"))

import lectin_metrics as lm  # noqa: E402

REQUIRED = [
    "accession", "gene", "species", "taxon_id", "label", "stratum", "evidence_level", "source",
    "quote", "tuned", "leakage", "curator", "date", "notes",
]  # fmt: skip


@pytest.fixture(scope="module")
def rows():
    with open(CTRL / "controls.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_required_columns(rows):
    for col in REQUIRED:
        assert col in rows[0], col


def test_no_duplicate_accessions(rows):
    accs = [r["accession"] for r in rows]
    assert len(accs) == len(set(accs))


def test_every_row_has_source_and_quote(rows):
    for r in rows:
        assert r["source"].strip(), r["accession"]
        assert r["quote"].strip(), r["accession"]


def test_leakage_values_valid(rows):
    for r in rows:
        assert r["leakage"] in lm.VALID_LEAKAGE, (r["accession"], r["leakage"])


def test_label_and_class_values_valid(rows):
    for r in rows:
        assert r["label"] in ("positive", "negative")
        if r["label"] == "positive":
            assert r["label_class"] in lm.POSITIVE_CLASSES, r["accession"]
            assert r["stratum"] == r["label_class"]
            assert r["tier"] in ("T1", "T2"), r["accession"]
            assert r["h_flag"] in ("H0", "H1a", "H1b", "H2"), r["accession"]
        else:
            assert r["label_class"] in lm.NEGATIVE_CLASSES, r["accession"]


def test_curator_tuned_date(rows):
    for r in rows:
        assert r["curator"].strip()
        assert r["tuned"] in ("yes", "no")
        assert len(r["date"]) == 10


def test_reasoned_negatives_have_a_reason(rows):
    for r in rows:
        if r["label"] == "negative" and r["label_class"] != "assumed":
            assert "Reason:" in r["notes"], r["accession"]
            assert r["negative_evidence"] != "assumed"


def test_assumed_rows_are_marked(rows):
    for r in rows:
        if r["label_class"] == "assumed":
            assert r["evidence_level"] == "assumed_no_reason"
            assert "ASSUMED" in r["notes"]
        else:
            assert r["evidence_level"] != "assumed_no_reason"


def test_sequences_and_clusters_cover_every_row(rows):
    ids = []
    with open(CTRL / "sequences.faa") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip())
    assert sorted(ids) == sorted(r["accession"] for r in rows)
    clu = dict(line.rstrip("\n").split("\t") for line in open(CTRL / "clusters.tsv"))
    assert set(clu) == {r["accession"] for r in rows}
    for r in rows:
        assert clu[r["accession"]] == r["cluster30"]


# Rows whose UniProt organism is above species rank (disclosed in unverified_notes.md, item 9).
RANK_ALLOW_LIST = {"A0A1U7Q1Z0": "genus"}


def test_taxon_ids_species_names_and_ranks_match_taxonomy(rows):
    names = "/srv/projects/db/taxonomy/names.dmp"
    nodes = "/srv/projects/db/taxonomy/nodes.dmp"
    if not (os.path.exists(names) and os.path.exists(nodes)):
        pytest.skip("taxonomy files not available")
    want = {r["taxon_id"]: r["species"] for r in rows}
    found, rank = {}, {}
    with open(names) as fh:
        for line in fh:
            f = line.split("\t|\t")
            if f[0] in want and "scientific name" in line:
                found[f[0]] = f[1]
    with open(nodes) as fh:
        for line in fh:
            f = line.split("\t|\t")
            if f[0] in want:
                rank[f[0]] = f[2]
    for t, sp in want.items():
        assert found.get(t) == sp, (t, sp, found.get(t))
    for r in rows:
        expected = RANK_ALLOW_LIST.get(r["accession"], "species")
        assert rank[r["taxon_id"]] == expected, (r["accession"], rank[r["taxon_id"]])


def test_every_positive_has_a_valid_selection_route(rows):
    valid = {"unilectin3d", "literature_or_keyword", "model_gated_pfam_xref"}
    for r in rows:
        if r["label"] == "positive":
            assert r["selection_route"] in valid, r["accession"]
        else:
            assert r["selection_route"].strip(), r["accession"]


def test_no_identical_sequences_under_two_accessions():
    seqs, cur = {}, None
    with open(CTRL / "sequences.faa") as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].strip()
                seqs[cur] = ""
            else:
                seqs[cur] += line.strip()
    assert len(set(seqs.values())) == len(seqs)


def test_in_reference_rows_name_a_seed_hit(rows):
    for r in rows:
        if r["leakage"] == "in_reference":
            assert r["seed_hits"], r["accession"]
        if r["leakage"] == "none":
            assert not r["seed_hits"], r["accession"]


def test_positive_and_negative_sets_are_disjoint(rows):
    pos = {r["accession"] for r in rows if r["label"] == "positive"}
    neg = {r["accession"] for r in rows if r["label"] == "negative"}
    assert not pos & neg
