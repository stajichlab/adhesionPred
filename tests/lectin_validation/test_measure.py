"""Tests of the grouped numbers of measure.py on a toy truth table.

The toy table has one cluster (C1) that holds a seed-similar (`in_reference`) row and a `none` row.
Review 1 found that the report counted such a cluster in both leakage groups. These tests keep the
cluster-level split correct.
"""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2] / "analysis" / "lectin_validation"
sys.path.insert(0, str(ROOT))

import lectin_metrics as lm  # noqa: E402

COLS = [
    "accession", "label", "label_class", "cluster30", "species", "taxon_id", "leakage", "seed_hits",
    "unilectin_class", "tier", "phylum", "negative_evidence", "protein_name", "length", "selection_route",
]  # fmt: skip


def _row(acc, label, cls, clu, leak, seed="", route="unilectin3d", ucls="AAL-like (PropLec6A)"):
    return {
        "accession": acc, "label": label, "label_class": cls, "cluster30": clu, "species": "Sp " + clu,
        "taxon_id": "1", "leakage": leak, "seed_hits": seed, "unilectin_class": ucls, "tier": "T2",
        "phylum": "Ascomycota", "negative_evidence": "", "protein_name": "x", "length": "100",
        "selection_route": route,
    }  # fmt: skip


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    d = tmp_path_factory.mktemp("toy")
    rows = [
        _row("a", "positive", "structure_backed", "C1", "in_reference", "PF07938:80.0"),
        _row("b", "positive", "structure_backed", "C1", "none"),
        _row("c", "positive", "structure_backed", "C2", "none", route="model_gated_pfam_xref"),
        _row("d", "positive", "structure_backed", "C3", "none"),
        _row("n1", "negative", "n1_shared", "C4", "none", route="shared_set", ucls=""),
        _row("as1", "negative", "assumed", "C5", "none", route="random_sample_filtered_by_pfam_xref", ucls=""),
    ]  # fmt: skip
    with open(d / "truth.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(d / "hits.tsv", "w") as fh:
        fh.write("accession\tmodel\tmodel_name\tscore\tevalue\tbest_domain_score\n")
        fh.write("a\tPF07938\tFungal_lectin\t100.0\t1e-30\t100\n")
        fh.write("c\tPF07938\tFungal_lectin\t40.0\t1e-10\t40\n")
    out = d / "res"
    subprocess.run(
        [sys.executable, str(ROOT / "measure.py"), "--truth", str(d / "truth.tsv"), "--hits", str(d / "hits.tsv"), "--out-dir", str(out)],
        check=True,
        cwd=str(d),
    )  # fmt: skip

    def read(name):
        with open(out / name) as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    return read


def test_cluster_level_leakage_split_does_not_double_count(results):
    rows = results("leakage_split.tsv")
    cl = {
        r["leakage"]: r
        for r in rows
        if r["group"] == "ALL_PRIMARY_POSITIVES" and r["basis"].startswith("cluster_level")
    }
    sim, ind = cl["cluster_has_in_reference_member"], cl["cluster_all_none"]
    assert (sim["n_clu"], sim["k_clu"]) == ("1", "1")
    assert (ind["n_clu"], ind["k_clu"]) == ("2", "1")
    assert int(sim["n_clu"]) + int(ind["n_clu"]) == 3  # clusters C1, C2, C3 once each


def test_row_level_split_is_labelled_and_differs(results):
    rows = results("leakage_split.tsv")
    r = next(
        x
        for x in rows
        if x["group"] == "ALL_PRIMARY_POSITIVES"
        and x["leakage"] == "none"
        and x["basis"].startswith("row_level")
    )
    assert r["n_clu"] == "3"  # C1 is counted here and in the in_reference group


def test_family_recall_seed_independent_is_cluster_level(results):
    rows = results("family_recall.tsv")
    ind = next(
        r
        for r in rows
        if r["model"] == "PF07938" and r["split"].startswith("cluster has no such member")
    )
    sim = next(
        r for r in rows if r["model"] == "PF07938" and r["split"].startswith("cluster has a member")
    )
    assert (sim["n_clu"], sim["k_clu"]) == ("1", "1")
    assert (ind["n_clu"], ind["k_clu"]) == ("2", "1")


def test_independent_route_group(results):
    rows = results("recall_by_group.tsv")
    ind = next(
        r
        for r in rows
        if r["group"] == "ALL_PRIMARY_independent_route" and r["model_set"] == "union15"
    )
    gated = next(
        r
        for r in rows
        if r["group"] == "ALL_PRIMARY_model_gated_route" and r["model_set"] == "union15"
    )
    assert (ind["n_seq"], ind["k_seq"]) == ("3", "1")
    assert (gated["n_seq"], gated["k_seq"]) == ("1", "1")


def test_assumed_is_not_in_false_hit_table(results):
    rows = results("false_hits.tsv")
    assert all(r["group"] != "assumed" for r in rows)
    cc = results("construction_check.tsv")
    assert {r["group"] for r in cc} == {"assumed_construction_check"}


def test_floors_apply_to_pf07938_score(results):
    rows = results("floors.tsv")
    r100 = next(
        r
        for r in rows
        if r["floor"] == "score>=100"
        and r["group"] == "ALL_PRIMARY_POSITIVES"
        and r["model_set"] == "PF07938_only"
    )
    r50 = next(
        r
        for r in rows
        if r["floor"] == "score>=50"
        and r["group"] == "ALL_PRIMARY_POSITIVES"
        and r["model_set"] == "PF07938_only"
    )
    assert r100["k_seq"] == "1"  # only a (100.0); c (40.0) is below the floor
    assert r50["k_seq"] == "1"


def test_strata_shared_clusters_file(results):
    rows = results("strata_shared_clusters.tsv")
    assert all(
        r["cluster"].startswith("#") for r in rows
    )  # no real cluster spans two strata in the toy


def test_cluster_seed_split_unit():
    members = [("a", "C1"), ("b", "C1"), ("c", "C2")]
    sim, ind = lm.cluster_seed_split(members, {"a"})
    assert sim == [("a", "C1"), ("b", "C1")]
    assert ind == [("c", "C2")]
    sim, ind = lm.cluster_seed_split(members, set())
    assert sim == [] and len(ind) == 3


def test_identity_bin():
    assert lm.identity_bin(None) == "no seed hit (<30%)"
    assert lm.identity_bin(100.0) == "98% or more"
    assert lm.identity_bin(98.0) == "98% or more"
    assert lm.identity_bin(95.0) == "90 to 97.9%"
    assert lm.identity_bin(60.0) == "50 to 89.9%"
    assert lm.identity_bin(30.0) == "30 to 49.9%"
