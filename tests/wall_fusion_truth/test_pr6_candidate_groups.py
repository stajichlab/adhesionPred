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


SPEC = TABLE.parents[1] / "2026-10-22-pr6-custom-hmms-design.md"
INPUTS = Path(__file__).resolve().parents[2] / "scripts/sorting_hat/refresh_3_10/inputs.tsv"
KEPT_PROTEOMES = ("Afum_Af293_UniProt", "Calb_SC5314", "Cimm_RS", "Scer_S288C")


def test_the_spec_records_the_owner_decisions_and_no_open_proposals():
    text = SPEC.read_text()
    assert "## 7. Owner decisions (final, 2026-10-22)" in text
    assert "Open questions for the owner" not in text
    assert "PROPOSAL" not in text
    for number in ("| 1 |", "| 9b |", "| 12 |"):
        assert number in text


def test_the_extra_call_limit_in_the_spec_is_computed_from_the_kept_proteomes():
    with open(INPUTS, newline="") as handle:
        sizes = {
            r["prot"]: int(r["n_proteins"])
            for r in csv.DictReader(handle, delimiter="\t")
            if r["n_proteins"].isdigit()
        }
    total = sum(sizes[p] for p in KEPT_PROTEOMES)
    text = SPEC.read_text()
    assert f"{total:,} proteins" in text
    assert f"{total / 10000 * 5:.1f}" in text
    for removed in ("Afum_Af293_Fungi5k", "Bder_ER3"):
        assert removed in text and removed not in KEPT_PROTEOMES
