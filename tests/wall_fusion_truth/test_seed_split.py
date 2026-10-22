"""Tests of seed_split.py (plan Task 5.3): a cluster with any member at 30% or more is seed-like.

Two classes are written: seed_class_slice (the models of the row's own slice) and seed_class_all39.
"""

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

# members: accession, class, slice, overlap_with
MEMBERS = [
    ("A1", "positive", "adhesion", ""),
    ("A2", "positive", "adhesion", ""),  # cluster of A1: A2 hits a family model at 35%
    (
        "W1",
        "positive",
        "wall",
        "",
    ),  # hits a FAMILY model at 35%: not a wall-slice hit, an all39 hit
    ("B1", "positive", "wall", ""),  # 20% only
    ("C1", "positive", "wall", ""),
    ("N1", "negative", "", ""),  # a mixed cluster with C1, no hit
    ("D1", "positive", "adhesion", ""),  # exactly 30% to a family model
    ("E1", "positive", "fusion", ""),  # 40% to a wall model but target coverage 0.4
    ("F1", "positive", "fusion", ""),  # 31% to a wall model
    ("N2", "negative", "", ""),  # a cluster of only negatives; hits a wall model at 40%
    ("X1", "positive", "wall", ""),
    (
        "X1#neg",
        "negative",
        "",
        "X1",
    ),  # a protein that is a control and a negative: one mixed cluster
]
CLUSTERS = [
    ("A1", "A1"), ("A1", "A2"), ("W1", "W1"), ("B1", "B1"), ("C1", "C1"), ("C1", "N1"),
    ("D1", "D1"), ("E1", "E1"), ("F1", "F1"), ("N2", "N2"), ("X1", "X1"), ("X1", "X1#neg"),
]  # fmt: skip
# query, target, identity, target coverage, E-value
HITS = [
    ("A2", "PFF1|seqx", 35.0, 0.8, 1e-10),
    ("W1", "PFF1|seqy", 35.0, 0.8, 1e-10),
    ("B1", "PFW1|seqz", 20.0, 0.9, 1e-5),
    ("D1", "PFF2|seqz", 30.0, 0.5, 1e-3),
    ("E1", "PFW2|seqw", 40.0, 0.4, 1e-9),
    ("F1", "PFW2|seqw", 31.0, 0.9, 1e-9),
    ("N2", "PFW1|seqv", 40.0, 0.9, 1e-9),
]


def write(tmp_path, members=MEMBERS, clusters=CLUSTERS, hits=HITS):
    (tmp_path / "members.tsv").write_text(
        "accession\tclass\tslice\toverlap_with\n" + "".join("\t".join(m) + "\n" for m in members)
    )
    (tmp_path / "clu.tsv").write_text("".join(f"{r}\t{m}\n" for r, m in clusters))
    (tmp_path / "hits.m8").write_text(
        "".join(f"{q}\t{t}\t{p}\t100\t0.9\t{tc}\t{e}\t200\n" for q, t, p, tc, e in hits)
    )
    (tmp_path / "wall.tsv").write_text("id\tname\nPFW1\ta\nPFW2\tb\n")
    (tmp_path / "family.tsv").write_text(
        "pfam_acc\tname\tactive\nPFF1\tf1\tyes\nPFF2\tf2\tyes\nPFF3\tf3\tno\n"
    )
    return [
        "--members", str(tmp_path / "members.tsv"), "--clusters", str(tmp_path / "clu.tsv"),
        "--seed-hits", str(tmp_path / "hits.m8"), "--wall-models", str(tmp_path / "wall.tsv"),
        "--family-table", str(tmp_path / "family.tsv"),
    ]  # fmt: skip


def read(path):
    rows = [ln.split("\t") for ln in Path(path).read_text().splitlines()]
    return [dict(zip(rows[0], r, strict=True)) for r in rows[1:]]


def run(tmp_path, **kw):
    out = tmp_path / "out"
    assert ss.main([*write(tmp_path, **kw), "--out", str(out)]) == 0
    mem = {r["accession"]: r for r in read(out / "members_split.tsv")}
    return out, mem


def test_35_percent_is_seed_like_and_20_percent_is_seed_independent_in_both_columns(tmp_path):
    _, mem = run(tmp_path)
    # a cluster with any member at 30% or more is seed-like, whichever member it is
    for a in ("A1", "A2"):
        assert mem[a]["seed_class_slice"] == mem[a]["seed_class_all39"] == "seed-like"
    assert mem["A2"]["member_hit_slice"] == "yes" and mem["A1"]["member_hit_slice"] == "no"
    assert mem["A2"]["best_identity_slice"] == "35.0" and mem["A2"]["nearest_model_slice"] == "PFF1"
    assert (
        mem["A2"]["leakage_suggested"] == "in_reference"
        and mem["A1"]["leakage_suggested"] == "none"
    )
    assert (
        mem["B1"]["seed_class_slice"] == mem["B1"]["seed_class_all39"] == "seed-independent"
    )  # 20%
    assert (
        mem["D1"]["seed_class_slice"] == mem["D1"]["seed_class_all39"] == "seed-like"
    )  # exactly 30%
    assert mem["F1"]["seed_class_slice"] == mem["F1"]["seed_class_all39"] == "seed-like"  # 31%
    # 40% but target coverage 0.4: no hit at all
    assert mem["E1"]["seed_class_slice"] == mem["E1"]["seed_class_all39"] == "seed-independent"


def test_a_wall_row_is_judged_against_wall_models_only_in_the_slice_column(tmp_path):
    _, mem = run(tmp_path)
    # W1 hits a family model: not a hit for a wall row (slice column), a hit for all 39
    assert mem["W1"]["seed_class_slice"] == "seed-independent"
    assert mem["W1"]["seed_class_all39"] == "seed-like"
    assert mem["W1"]["member_hit_slice"] == "no" and mem["W1"]["member_hit_all39"] == "yes"
    assert mem["W1"]["nearest_model_all39"] == "PFF1" and mem["W1"]["nearest_model_slice"] == ""
    assert mem["W1"]["leakage_suggested"] == "none"  # a control: judged by its slice models


def test_a_cluster_of_only_negatives_has_no_slice_class(tmp_path):
    _, mem = run(tmp_path)
    assert mem["N2"]["seed_class_slice"] == "not_applicable"
    assert mem["N2"]["seed_class_all39"] == "seed-like"
    assert mem["N2"]["leakage_suggested"] == "in_reference"  # a negative is judged by all 39


def test_mixed_clusters_are_flagged_including_the_overlap_copy(tmp_path):
    out, mem = run(tmp_path)
    cl = {r["cluster_id"]: r for r in read(out / "clusters_split.tsv")}
    mixed = {r["cluster_id"] for r in cl.values() if r["mixed_pos_neg"] == "yes"}
    assert mixed == {mem["C1"]["cluster_id"], mem["X1"]["cluster_id"]}
    assert mem["C1"]["in_mixed_cluster"] == mem["N1"]["in_mixed_cluster"] == "yes"
    assert mem["A1"]["in_mixed_cluster"] == "no"
    x = cl[mem["X1"]["cluster_id"]]
    assert (
        x["has_overlap_copy"] == "yes" and cl[mem["C1"]["cluster_id"]]["has_overlap_copy"] == "no"
    )
    assert mem["X1#neg"]["overlap_with"] == "X1"
    # a mixed cluster is judged by its control members (C1 has no hit)
    assert x["seed_class_slice"] == "seed-independent"


def test_summary_is_per_class_and_slice_never_pooled(tmp_path):
    out, _ = run(tmp_path)
    s = {(r["class"], r["slice"]): r for r in read(out / "summary.tsv")}
    assert s[("positive", "adhesion")]["slice_seed_like"] == "2"
    assert s[("positive", "adhesion")]["all39_seed_like"] == "2"
    assert s[("positive", "wall")]["slice_seed_like"] == "0"
    assert s[("positive", "wall")]["all39_seed_like"] == "1"  # W1
    assert s[("positive", "fusion")]["clusters"] == "2"
    assert s[("negative", "")]["clusters"] == "3"
    assert ("positive", "") not in s  # no pooled line


def test_identity_as_a_fraction_is_read_as_a_percentage(tmp_path):
    hits = [(q, t, p / 100, tc, e) for q, t, p, tc, e in HITS]
    _, mem = run(tmp_path, hits=hits)
    assert (
        mem["A2"]["seed_class_slice"] == "seed-like"
        and mem["B1"]["seed_class_slice"] == "seed-independent"
    )


def test_clusters_tsv_is_the_common_rules_file(tmp_path):
    out, _ = run(tmp_path)
    lines = (out / "clusters.tsv").read_text().splitlines()
    assert lines[0] == "accession\tcluster_id" and len(lines) == 1 + len(MEMBERS)


def test_a_member_that_was_not_clustered_stops_the_script(tmp_path):
    args = write(tmp_path, clusters=CLUSTERS[:-1])
    with pytest.raises(ValueError, match="clusters and members differ"):
        ss.main([*args, "--out", str(tmp_path / "out")])


def test_a_hit_for_an_unknown_sequence_stops_the_script(tmp_path):
    args = write(tmp_path, hits=[*HITS, ("ZZ", "PFW1|s", 50.0, 0.9, 1e-9)])
    with pytest.raises(ValueError, match="not members"):
        ss.main([*args, "--out", str(tmp_path / "out")])


def test_a_hit_to_a_model_in_neither_list_stops_the_script(tmp_path):
    args = write(tmp_path, hits=[*HITS, ("A1", "PF99999|s", 50.0, 0.9, 1e-9)])
    with pytest.raises(ValueError, match="neither model list"):
        ss.main([*args, "--out", str(tmp_path / "out")])


def test_a_bad_class_or_slice_stops_the_script(tmp_path):
    args = write(tmp_path, [("A1", "weird", "wall", "")], [("A1", "A1")], [])
    with pytest.raises(ValueError, match="is not one of"):
        ss.main([*args, "--out", str(tmp_path / "out")])
    args = write(tmp_path, [("A1", "positive", "odd", "")], [("A1", "A1")], [])
    with pytest.raises(ValueError, match="slice 'odd'"):
        ss.main([*args, "--out", str(tmp_path / "out")])


def test_clade_stratum_join_reads_both_files_and_lists_the_family_reps(tmp_path):
    cdir = tmp_path / "set"
    cdir.mkdir()
    head = (
        "accession\tclade_fine\tclade_coarse\tcluster_stratum_fine\tcluster_stratum_coarse\t"
        "kiss_cluster_number\tjgi_accept\n"
    )
    (cdir / "controls.tsv").write_text(
        head + "A1\tSaccharomycotina\tSaccharomycotina\tshared\tclade_only\t5\trow_verified\n"
    )
    (cdir / "controls_local_sequence.tsv").write_text(
        head + "A2\tSordariomycetes\tPezizomycotina\tshared\tclade_only\t5\tportal_validated\n"
    )
    (cdir / "family_reps.tsv").write_text(
        "accession\tslice\tkiss_cluster_number\nP1\twall\t7\nP2\twall\t7\nP3\twall\t8\nP4\tfusion\t9\n"
    )
    out = tmp_path / "out"
    ss.main([*write(tmp_path), "--controls", str(cdir), "--out", str(out)])
    first = next(r for r in read(out / "clusters_split.tsv") if "A1" in r["members"].split(","))
    assert first["clades_fine"] == "Saccharomycotina,Sordariomycetes"  # the local file was read too
    assert first["strata_fine"] == "shared" and first["strata_coarse"] == "clade_only"
    assert first["jgi_accept"] == "portal_validated,row_verified"
    s = {(r["class"], r["slice"]): r for r in read(out / "summary.tsv")}
    assert s[("family_rep (not clustered)", "wall")]["clusters"] == "2"  # distinct Kiss clusters
    assert s[("family_rep (not clustered)", "fusion")]["slice_seed_like"] == ""


def test_the_model_sets_are_the_21_and_the_18_of_the_repository():
    sets = ss.read_model_sets()
    assert len(sets["wall_fusion"]) == 21 and len(sets["family"]) == 18
    assert not sets["wall_fusion"] & sets["family"]
