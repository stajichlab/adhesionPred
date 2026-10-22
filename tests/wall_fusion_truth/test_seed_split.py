"""Tests of seed_split.py (plan Task 5.3): a cluster with any member at 30% or more is seed-like."""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "seed_split", ROOT / "analysis/wall_fusion_truth/seed_split.py"
)
ss = importlib.util.module_from_spec(spec)
sys.modules["seed_split"] = ss
spec.loader.exec_module(ss)


def write(tmp_path, members, clusters, hits):
    (tmp_path / "members.tsv").write_text(
        "accession\tclass\tslice\n" + "".join(f"{a}\t{c}\t{s}\n" for a, c, s in members)
    )
    (tmp_path / "clu.tsv").write_text("".join(f"{r}\t{m}\n" for r, m in clusters))
    (tmp_path / "hits.m8").write_text(
        "".join(f"{q}\t{t}\t{p}\t100\t0.9\t{tc}\t{e}\t200\n" for q, t, p, tc, e in hits)
    )
    return [
        "--members",
        str(tmp_path / "members.tsv"),
        "--clusters",
        str(tmp_path / "clu.tsv"),
        "--seed-hits",
        str(tmp_path / "hits.m8"),
    ]


MEMBERS = [
    ("A1", "positive", "adhesion"),
    ("A2", "positive", "adhesion"),  # cluster of A1: one member at 35%
    ("B1", "positive", "wall"),  # cluster of B1: 20% only
    ("C1", "positive", "wall"),
    ("N1", "negative", ""),  # a mixed cluster, no hit
    ("D1", "positive", "adhesion"),  # exactly 30%
    ("E1", "positive", "fusion"),  # 40% but target coverage 0.4
]
CLUSTERS = [
    ("A1", "A1"),
    ("A1", "A2"),
    ("B1", "B1"),
    ("C1", "C1"),
    ("C1", "N1"),
    ("D1", "D1"),
    ("E1", "E1"),
]
HITS = [
    ("A2", "PF00001|seqx", 35.0, 0.8, 1e-10),
    ("B1", "PF00002|seqy", 20.0, 0.9, 1e-5),
    ("D1", "PF00003|seqz", 30.0, 0.5, 1e-3),
    ("E1", "PF00004|seqw", 40.0, 0.4, 1e-9),
]


def read(path):
    rows = [ln.split("\t") for ln in Path(path).read_text().splitlines()]
    return [dict(zip(rows[0], r, strict=True)) for r in rows[1:]]


def test_35_percent_is_seed_like_and_20_percent_is_seed_independent(tmp_path):
    out = tmp_path / "out"
    assert ss.main([*write(tmp_path, MEMBERS, CLUSTERS, HITS), "--out", str(out)]) == 0
    cl = {r["cluster_id"]: r for r in read(out / "clusters_split.tsv")}
    mem = {r["accession"]: r for r in read(out / "members_split.tsv")}
    # a cluster with any member at 30% or more is seed-like, whichever member it is
    assert mem["A1"]["seed_class"] == mem["A2"]["seed_class"] == "seed-like"
    assert mem["A2"]["member_seed_hit"] == "yes" and mem["A1"]["member_seed_hit"] == "no"
    assert mem["A2"]["best_identity"] == "35.0" and mem["A2"]["nearest_model"] == "PF00001"
    assert (
        mem["A2"]["leakage_suggested"] == "in_reference"
        and mem["A1"]["leakage_suggested"] == "none"
    )
    assert mem["B1"]["seed_class"] == "seed-independent"  # 20% is below the threshold
    assert mem["D1"]["seed_class"] == "seed-like"  # 30% counts (30% or more)
    assert mem["E1"]["seed_class"] == "seed-independent"  # target coverage 0.4 is below 0.5
    assert mem["C1"]["seed_class"] == "seed-independent"
    mixed = [r for r in cl.values() if r["mixed_pos_neg"] == "yes"]
    assert len(mixed) == 1 and mixed[0]["classes"] == "negative,positive"
    assert len(cl) == 5


def test_summary_is_per_class_and_slice_never_pooled(tmp_path):
    out = tmp_path / "out"
    ss.main([*write(tmp_path, MEMBERS, CLUSTERS, HITS), "--out", str(out)])
    s = {(r["class"], r["slice"]): r for r in read(out / "summary.tsv")}
    assert s[("positive", "adhesion")]["clusters_seed_like"] == "2"
    assert s[("positive", "wall")]["clusters_seed_independent"] == "2"
    assert s[("positive", "fusion")]["clusters"] == "1"
    assert s[("negative", "")]["clusters"] == "1"
    assert ("positive", "") not in s  # no pooled line


def test_identity_as_a_fraction_is_read_as_a_percentage(tmp_path):
    hits = [(q, t, p / 100, tc, e) for q, t, p, tc, e in HITS]
    out = tmp_path / "out"
    ss.main([*write(tmp_path, MEMBERS, CLUSTERS, hits), "--out", str(out)])
    mem = {r["accession"]: r for r in read(out / "members_split.tsv")}
    assert mem["A2"]["seed_class"] == "seed-like" and mem["B1"]["seed_class"] == "seed-independent"


def test_clusters_tsv_is_the_common_rules_file(tmp_path):
    out = tmp_path / "out"
    ss.main([*write(tmp_path, MEMBERS, CLUSTERS, HITS), "--out", str(out)])
    lines = (out / "clusters.tsv").read_text().splitlines()
    assert lines[0] == "accession\tcluster_id" and len(lines) == 1 + len(MEMBERS)


def test_a_member_that_was_not_clustered_stops_the_script(tmp_path):
    args = write(tmp_path, MEMBERS, CLUSTERS[:-1], HITS[:-1])
    with pytest.raises(ValueError, match="clusters and members differ"):
        ss.main([*args, "--out", str(tmp_path / "out")])


def test_a_hit_for_an_unknown_sequence_stops_the_script(tmp_path):
    args = write(tmp_path, MEMBERS, CLUSTERS, [*HITS, ("ZZ", "PF00001|s", 50.0, 0.9, 1e-9)])
    with pytest.raises(ValueError, match="not members"):
        ss.main([*args, "--out", str(tmp_path / "out")])


def test_a_bad_class_stops_the_script(tmp_path):
    args = write(tmp_path, [("A1", "weird", "wall")], [("A1", "A1")], [])
    with pytest.raises(ValueError, match="is not one of"):
        ss.main([*args, "--out", str(tmp_path / "out")])


def test_clade_and_stratum_join_from_a_control_set(tmp_path):
    cdir = tmp_path / "set"
    cdir.mkdir()
    (cdir / "controls.tsv").write_text(
        "accession\tclade_fine\tclade_coarse\tcluster_stratum_fine\tcluster_stratum_coarse\t"
        "kiss_cluster_number\tjgi_accept\n"
        "A1\tSaccharomycotina\tSaccharomycotina\tshared\tclade_only\t5\trow_verified\n"
        "A2\tSordariomycetes\tPezizomycotina\tshared\tclade_only\t5\tportal_validated\n"
    )
    out = tmp_path / "out"
    ss.main([*write(tmp_path, MEMBERS, CLUSTERS, HITS), "--controls", str(cdir), "--out", str(out)])
    first = next(r for r in read(out / "clusters_split.tsv") if "A1" in r["members"])
    assert first["clades_fine"] == "Saccharomycotina,Sordariomycetes"
    assert first["strata_fine"] == "shared" and first["strata_coarse"] == "clade_only"
    assert first["jgi_accept"] == "portal_validated,row_verified"
