"""Structure of the PR 6 candidate-group table (a proposal made from annotation text; ids and text only)."""

import csv
import re
from collections import Counter
from pathlib import Path

TABLE = (
    Path(__file__).resolve().parents[2]
    / "docs/superpowers/specs/2026-10-22-pr6-custom-hmms-data/candidate_groups.tsv"
)
SLICES = {"adhesion": 12, "wall": 28, "fusion": 4}
STRATA = {"shared", "clade_only"}
SEED = {"seed-like", "seed-independent"}


def rows():
    with open(TABLE, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_every_missed_cluster_appears_once_with_its_slice_count():
    data = rows()
    ids = [r["cluster_id"] for r in data]
    assert len(ids) == len(set(ids)) == 44
    assert dict(Counter(r["slice"] for r in data)) == SLICES


def test_every_cluster_has_exactly_one_group_of_its_own_slice():
    for r in rows():
        assert r["group"][0] == {"adhesion": "A", "wall": "W", "fusion": "F"}[r["slice"]], r[
            "cluster_id"
        ]
    assert len({r["group"] for r in rows()}) == 12


def test_values_are_from_the_closed_vocabularies():
    for r in rows():
        assert r["stratum_fine"] in STRATA
        assert (
            r["seed_class_slice"] == "seed-independent"
        )  # the table lists clusters the slice's own models miss
        assert r["seed_class_all39"] in SEED
        assert r["jgi_accept"], r["cluster_id"]  # no blank: a locus-tag row says so


def test_the_table_holds_no_sequence():
    for line in TABLE.read_text().splitlines():
        for field in line.split("\t"):
            assert not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{40,}", field)
