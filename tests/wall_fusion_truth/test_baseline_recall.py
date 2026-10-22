"""baseline_recall.py (Task 5.4 draft): hand-counted fixture, no data, no HPCC."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "baseline_recall", ROOT / "analysis/wall_fusion_truth/baseline_recall.py"
)
br = importlib.util.module_from_spec(spec)
spec.loader.exec_module(br)

CTRL_COLS = [
    "accession",
    "gene",
    "species",
    "stratum",
    "cluster_id",
    "clade_fine",
    "clade_coarse",
    "cluster_stratum_fine",
    "cluster_stratum_coarse",
    "jgi_accept",
]
MEM_COLS = ["accession", "cluster_id", "class", "slice", "in_mixed_cluster"]
CL_COLS = [
    "cluster_id",
    "seed_class_slice",
    "seed_class_all39",
    "strata_fine",
    "best_identity_all39",
]


def tsv(path, cols, rows):
    path.write_text(
        "\t".join(cols) + "\n" + "".join("\t".join(r.get(c, "") for c in cols) + "\n" for r in rows)
    )


def dom(rows, ga=True, ok=True):
    head = "# target name accession tlen query name accession qlen E-value\n"
    opt = f"# Option settings: hmmsearch {'--cut_ga ' if ga else ''}--domtblout x\n"
    body = "".join(
        f"{t} - 100 {m} {acc} 90 1e-20 50.0 0.1 1 1 1e-20 1e-20 50.0 0.1 1 90 1 90 1 90 0.9 d\n"
        for t, m, acc in rows
    )
    return head + opt + body + ("# [ok]\n" if ok else "")


def ctrl(acc, gene, sl, cl, clade, coarse, stratum, accept="row_verified"):
    return {
        "accession": acc, "gene": gene, "species": "Sp " + clade, "stratum": sl,
        "cluster_id": cl, "clade_fine": clade, "clade_coarse": coarse,
        "cluster_stratum_fine": stratum, "cluster_stratum_coarse": stratum,
        "jgi_accept": accept,
    }  # fmt: skip


@pytest.fixture
def fx(tmp_path):
    (tmp_path / "ctl").mkdir()
    tsv(
        tmp_path / "wall.tsv",
        ["id", "set"],
        [
            {"id": "PW0001", "set": "wall_remodeling"},
            {"id": "PW0002", "set": "wall_synthesis"},
            {"id": "PX0001", "set": "cell_fusion_candidate_models"},
        ],
    )
    tsv(
        tmp_path / "fam.tsv",
        ["pfam_acc", "active"],
        [{"pfam_acc": "PF0001", "active": "yes"}, {"pfam_acc": "PF0002", "active": "yes"}]
        + [{"pfam_acc": "PF0003", "active": "no"}],
    )
    controls = [
        ctrl("adh1", "g1", "adhesion", "C1", "Saccharomycotina", "Sacc", "clade_only"),
        ctrl("adh2", "g2", "adhesion", "C1", "Saccharomycotina", "Sacc", "clade_only"),
        ctrl("adh3", "g3", "adhesion", "C2", "Sordariomycetes", "Pez", "clade_only"),
        ctrl("wall1", "g4", "wall", "C3", "Eurotiomycetes", "Pez", "shared"),
        ctrl("wall2", "g5", "wall", "C4", "Eurotiomycetes", "Pez", "clade_only"),
        ctrl("fus1", "g6", "fusion", "C5", "Sordariomycetes", "Pez", "clade_only"),
    ]
    tsv(tmp_path / "ctl/controls.tsv", CTRL_COLS, controls[:-1])
    tsv(tmp_path / "ctl/controls_local_sequence.tsv", CTRL_COLS, controls[-1:])
    tsv(
        tmp_path / "ctl/family_reps.tsv",
        ["accession", "slice"],
        [{"accession": "rep1", "slice": "wall"}, {"accession": "rep2", "slice": "adhesion"}],
    )
    members = [
        {"accession": c["accession"], "cluster_id": c["cluster_id"], "class": "positive",
         "slice": c["stratum"], "in_mixed_cluster": "no"} for c in controls
    ]  # fmt: skip
    members += [
        {"accession": "neg1", "cluster_id": "N1", "class": "negative", "slice": "",
         "in_mixed_cluster": "no"},
        {"accession": "neg2", "cluster_id": "C3", "class": "negative", "slice": "",
         "in_mixed_cluster": "yes"},
    ]  # fmt: skip
    tsv(tmp_path / "members.tsv", MEM_COLS, members)
    sc = {
        "C1": ("seed-like", "seed-like"),
        "C2": ("seed-independent", "seed-independent"),
        "C3": ("seed-independent", "seed-like"),
        "C4": ("seed-independent", "seed-independent"),
        "C5": ("seed-independent", "seed-independent"),
        "N1": ("not_applicable", "seed-independent"),
    }
    tsv(
        tmp_path / "clusters.tsv",
        CL_COLS,
        [
            {"cluster_id": k, "seed_class_slice": v[0], "seed_class_all39": v[1],
             "strata_fine": "clade_only", "best_identity_all39": "12.0"} for k, v in sc.items()
        ],
    )  # fmt: skip
    (tmp_path / "dom.txt").write_text(
        dom(
            [
                ("adh1", "famA", "PF0001.24"),
                ("wall1", "wallA", "PW0001.10"),
                (
                    "wall1",
                    "famA",
                    "PF0001.24",
                ),  # a wall row hit by an adhesion family: false-hit panel
                ("fus1", "fusA", "PX0001.3"),
                ("neg1", "famA", "PF0001.24"),
                ("neg2", "wallA", "PW0001.10"),
            ]
        )
    )
    (tmp_path / "dom_rep.txt").write_text(dom([("rep1", "wallA", "PW0002.1")]))
    return tmp_path


def run(fx, extra=()):
    out = fx / "out"
    argv = ["--domtbl", str(fx / "dom.txt"), "--members", str(fx / "members.tsv"),
            "--clusters", str(fx / "clusters.tsv"), "--controls", str(fx / "ctl"),
            "--wall-models", str(fx / "wall.tsv"), "--family-table", str(fx / "fam.tsv"),
            "--out", str(out), *extra]  # fmt: skip
    assert br.main(argv) == 0
    return out


def table(path):
    rows = [ln.split("\t") for ln in path.read_text().splitlines()]
    return [dict(zip(rows[0], r, strict=True)) for r in rows[1:]]


def test_recall_equals_the_hand_count(fx):
    out = run(fx)
    t = {
        (r["slice"], r["seed_class_kind"], r["seed_class"]): r
        for r in table(out / "recall_by_slice.tsv")
    }
    a = t[("adhesion", "slice", "seed-like")]  # C1: adh1 hit, adh2 not
    assert (a["n_rows"], a["k_rows"], a["n_clusters"], a["k_clusters_any"], a["k_clusters_all"]) == (
        "2", "1", "1", "1", "0",
    )  # fmt: skip
    assert a["status"] == "descriptive"  # 1 cluster, below the floor of 5
    b = t[("adhesion", "slice", "seed-independent")]  # C2: adh3 not hit
    assert (b["n_rows"], b["k_rows"], b["n_clusters"], b["k_clusters_any"]) == ("1", "0", "1", "0")
    # wall1 is hit by a wall model: recalled; wall2 is not. C3 differs between the two seed classes.
    assert t[("wall", "slice", "seed-independent")]["k_rows"] == "1"
    assert t[("wall", "all39", "seed-like")]["n_rows"] == "1"
    assert t[("wall", "all39", "seed-independent")]["n_rows"] == "1"


def test_no_pooled_line_and_seed_classes_stay_apart(fx):
    out = run(fx)
    for name in ("recall_by_slice.tsv", "recall_by_clade.tsv"):
        for r in table(out / name):
            assert r["slice"] in br.SLICES
            assert r["seed_class"] in br.SEED_CLASSES
            assert r["seed_class_kind"] in br.KINDS
    kinds = {r["seed_class_kind"] for r in table(out / "recall_by_slice.tsv")}
    assert kinds == {"slice", "all39"}


def test_fusion_is_descriptive_and_every_number_has_its_n(fx):
    out = run(fx)
    fus = [r for r in table(out / "recall_by_slice.tsv") if r["slice"] == "fusion"]
    assert {r["status"] for r in fus} <= {"descriptive", "none"}
    for r in table(out / "recall_by_clade.tsv"):
        assert r["n_rows"].isdigit() and r["n_clusters"].isdigit()
        if r["slice"] == "fusion":
            assert r["status"] == "descriptive"
    members = {r["accession"]: r for r in table(out / "members_hits.tsv")}
    assert members["fus1"]["hit_fusion_models"] == "yes"


def test_clade_table_has_strata_and_floor_status(fx):
    out = run(fx)
    rows = [
        r
        for r in table(out / "recall_by_clade.tsv")
        if r["level"] == "fine" and r["slice"] == "wall" and r["clade"] == "Eurotiomycetes"
    ]
    cells = {(r["stratum"], r["seed_class_kind"], r["seed_class"]): r for r in rows}
    shared = cells[("shared", "slice", "seed-independent")]  # wall1, cluster C3
    assert (shared["n_rows"], shared["k_rows"], shared["n_clusters"]) == ("1", "1", "1")
    only = cells[("clade_only", "slice", "seed-independent")]  # wall2, cluster C4
    assert (only["n_rows"], only["k_rows"]) == ("1", "0")
    assert all(r["status"] != "estimated" for r in rows)  # no cell has 5 clusters


def test_wall_rows_are_a_false_hit_panel_and_never_a_negative(fx):
    out = run(fx)
    panel = {r["family_model"]: r for r in table(out / "wall_rows_false_hits.tsv")}
    assert panel["PF0001"]["n_wall_rows"] == "2" and panel["PF0001"]["n_hit_rows"] == "1"
    assert panel["PF0002"]["n_hit_rows"] == "0"
    assert "never negatives" in panel["PF0001"]["note"]
    neg = table(out / "negatives_false_hits.tsv")
    # the negative counts come from the class-negative members only: wall1 is not among them
    assert all("wall1" not in r["hit_accessions"] for r in neg)
    cell_ = {(r["model"], r["cluster_composition"]): r for r in neg}
    assert cell_[("PF0001", "negative_only")]["n_negatives"] == "1"
    assert cell_[("PF0001", "negative_only")]["n_hit"] == "1"
    assert cell_[("PW0001", "mixed")]["n_hit"] == "1"
    assert cell_[("PW0001", "negative_only")]["n_hit"] == "0"


def test_the_clusters_pfam_misses_are_listed(fx):
    out = run(fx)
    missed = {r["cluster_id"]: r for r in table(out / "clusters_missed.tsv")}
    assert set(missed) == {"C2", "C4"}  # C1 and C3 have a hit, C5 (fusion) has a fusion-model hit
    assert missed["C2"]["slice"] == "adhesion" and missed["C2"]["genes"] == "g3"
    assert missed["C4"]["seed_class_slice"] == "seed-independent"


def test_family_reps_are_a_flagged_descriptive_group_only(fx):
    out = run(fx, ["--domtbl-famrep", str(fx / "dom_rep.txt")])
    rows = {r["slice"]: r for r in table(out / "family_rep_descriptive.tsv")}
    assert rows["wall"]["k_hit_own_models"] == "1" and rows["adhesion"]["k_hit_own_models"] == "0"
    for r in rows.values():
        assert r["group"] == "family_rep" and r["counts_toward_floor"] == "no"
        assert r["may_be_test_control"] == "no"
    # a representative is not a control: it appears in no recall table
    members = {r["accession"] for r in table(out / "members_hits.tsv")}
    assert not ({"rep1", "rep2"} & members)


def test_model_without_ga_is_reported_not_folded_in(fx):
    (fx / "noga.txt").write_text("PF0002\n")
    (fx / "dom_noga.txt").write_text(dom([("adh2", "famB", "PF0002.5")], ga=False))
    out = run(
        fx, ["--models-no-ga", str(fx / "noga.txt"), "--domtbl-noga", str(fx / "dom_noga.txt")]
    )
    rows = table(out / "no_ga_models.tsv")
    assert rows[0]["model"] == "PF0002" and rows[0]["n_sequences_hit"] == "1"
    assert "not in any recall number" in rows[0]["note"]
    members = {r["accession"]: r for r in table(out / "members_hits.tsv")}
    assert members["adh2"]["hit_own_slice_models"] == "no"  # the E-value hit did not count


def test_writes_only_its_tables_and_no_call_or_status(fx):
    out = run(fx)
    assert {p.name for p in out.iterdir()} == set(br.OUTPUTS)
    assert not any(
        p.name.startswith(("call", "status")) or "status" in p.name for p in out.iterdir()
    )
    assert not any(ln.startswith("import cellsurface_sorting_hat") for ln in
                   (ROOT / "analysis/wall_fusion_truth/baseline_recall.py").read_text().splitlines())  # fmt: skip


def test_deterministic_and_refuses_a_non_empty_out(fx):
    out = run(fx)
    first = {p.name: p.read_bytes() for p in out.iterdir()}
    (fx / "out2").mkdir()
    argv = ["--domtbl", str(fx / "dom.txt"), "--members", str(fx / "members.tsv"),
            "--clusters", str(fx / "clusters.tsv"), "--controls", str(fx / "ctl"),
            "--wall-models", str(fx / "wall.tsv"), "--family-table", str(fx / "fam.tsv"),
            "--out", str(fx / "out2")]  # fmt: skip
    assert br.main(argv) == 0
    assert first == {p.name: p.read_bytes() for p in (fx / "out2").iterdir()}
    with pytest.raises(SystemExit):
        br.main(argv)  # out2 is no longer empty


def test_domtbl_must_be_finished_and_made_with_cut_ga(tmp_path):
    p = tmp_path / "d.txt"
    p.write_text(dom([("a", "m", "PF0001.1")], ok=False))
    with pytest.raises(ValueError, match="trailer"):
        br.parse_domtbl(p)
    p.write_text(dom([("a", "m", "PF0001.1")], ga=False))
    with pytest.raises(ValueError, match="cut_ga"):
        br.parse_domtbl(p)
    assert br.parse_domtbl(p, require_cut_ga=False) == {"a": {"PF0001"}}


def test_a_hit_for_an_unknown_sequence_stops_the_run(fx):
    (fx / "dom.txt").write_text(dom([("ghost", "famA", "PF0001.24")]))
    with pytest.raises(SystemExit, match="not in members"):
        run(fx)


def test_wilson_interval_is_sane():
    lo, hi = br.wilson(5, 10)
    assert 0.2 < lo < 0.5 < hi < 0.8
    assert br.wilson(0, 0)[0] != br.wilson(0, 0)[0]  # nan
