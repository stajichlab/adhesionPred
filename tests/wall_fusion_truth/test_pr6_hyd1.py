"""PR 6 hyd1_like drafts: members, units, training-rule boundaries, manifests, the freeze guard, comparators, the pre-check and static checks of the
sbatch scripts. Fixtures only; the committed drift record is read for the member counts."""

import csv
import importlib.util
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
DRIFT = ROOT / "docs/reports/data/pr6_drift_check_2026-10-22"
PROPOSED = WF / "pr6_hyd1_freeze_PROPOSED.json"


def load(name):
    spec = importlib.util.spec_from_file_location(name, WF / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


mem = load("pr6_hyd1_members")
bld = load("pr6_hyd1_build")
rep = load("pr6_hyd1_report")
pre = load("pr6_hyd1_negatives_precheck")
e33b = sys.modules["pr6_ecm33_build"]


def drift(hit, clade, qcov="0.9", cls="supported_homolog", query="A0ACE4KFX7"):
    return {"query": query, "hit": hit, "clade_fine": clade, "query_cover": qcov, "class": cls}


def member(i, clade, query=False):
    return {
        "hit_id": "A0ACE4KFX7" if query else f"UniRef50_{i}",
        "query_accession": "A0ACE4KFX7" if query else "",
        "clade_fine": clade,
        "is_query": "yes" if query else "no",
    }


def fixture_members(spec):
    """spec: {clade: n}; the query is added to Eurotiomycetes when it is in spec."""
    rows = []
    for clade, n in spec.items():
        for i in range(n):
            rows.append(
                member(f"{clade[:3]}{i}", clade, query=(clade == "Eurotiomycetes" and i == 0))
            )
    return rows


# ------------------------------------------------------------------------------------------------ members


def test_member_selection_query_cover_boundary_and_no_clade():
    d = [
        drift("UniRef50_a", "Sordariomycetes", "0.500"),  # kept: >= 0.5
        drift("UniRef50_b", "Sordariomycetes", "0.499"),  # dropped: query coverage below 0.5
        drift("UniRef50_c", "(no clade)", "0.99"),  # supported, not a clade: excluded list
        drift("UniRef50_d", "", "0.99"),  # same
        drift("UniRef50_e", "Dothideomycetes", "0.9", cls="unsupported"),
        drift("UniRef50_f", "Dothideomycetes", "0.9", cls="bias_driven"),
        drift("UniRef50_g", "Dothideomycetes", "0.9", query="P32323"),  # another query
    ]
    members, excluded = mem.select_members(d, [])
    assert {m["hit_id"] for m in members} == {"A0ACE4KFX7", "UniRef50_a"}
    assert [m["is_query"] for m in members if m["hit_id"] == "A0ACE4KFX7"] == ["yes"]
    assert sorted(e["hit_id"] for e in excluded) == ["UniRef50_c", "UniRef50_d"]
    assert "(no clade)" not in {m["clade_fine"] for m in members}


def test_real_drift_record_gives_the_documented_members():
    members, excluded = mem.select_members(
        mem.read_tsv(DRIFT / "drift_hits.tsv"), mem.read_tsv(DRIFT / "hits_manifest.tsv")
    )
    hits = [m for m in members if m["is_query"] != "yes"]
    assert len(hits) == 29 and len(excluded) == 3 and len(members) == 30
    by = {}
    for m in members:
        by[m["clade_fine"]] = by.get(m["clade_fine"], 0) + 1
    assert by == {
        "Sordariomycetes": 19,
        "Dothideomycetes": 7,
        "Leotiomycetes": 3,
        "Eurotiomycetes": 1,
    }
    assert all(float(m["query_cover"]) >= 0.5 for m in members)
    # 27 of the 32 supported hits are reached only by iteration; 5 are in the one-pass set (3 named + 2 without a clade)
    assert sum(1 for m in hits if m["in_phmmer"] == "yes") == 3
    assert sum(1 for e in excluded if e["in_phmmer"] == "yes") == 2
    assert all(m["taxid"] for m in members)


def test_query_clade_constant_matches_the_kiss_control_set():
    with open(ROOT / "data/controls/wall-fusion-kiss2019/controls.tsv", encoding="utf-8-sig") as f:
        row = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}[mem.QUERY]
    assert row["clade_fine"] == mem.QUERY_CLADE and row["gene"] == "hyd1_like"
    assert row["stratum"] == "adhesion" and row["label"] == "positive"


# ------------------------------------------------------------------------------------------------ units and training rule


def test_real_units_under_both_readings():
    members, _ = mem.select_members(
        mem.read_tsv(DRIFT / "drift_hits.tsv"), mem.read_tsv(DRIFT / "hits_manifest.tsv")
    )
    tot = {u["unit"]: u["status"] for u in mem.unit_table(members, "total")}
    per = {u["unit"]: u["status"] for u in mem.unit_table(members, "per_clade")}
    assert tot == {
        "clade:Sordariomycetes": "unit",
        "clade:Dothideomycetes": "unit",
        "clade:Leotiomycetes": "too_few_members",
        "clade:Eurotiomycetes": "kiss_single",
    }
    # per_clade: holding out Sordariomycetes leaves only Dothideomycetes with >= 5 sequences, and the reverse
    assert (
        per["clade:Sordariomycetes"] == "not_buildable"
        and per["clade:Dothideomycetes"] == "not_buildable"
    )
    assert (
        per["clade:Eurotiomycetes"] == "kiss_single"
        and per["clade:Leotiomycetes"] == "too_few_members"
    )
    with pytest.raises(ValueError):
        mem.unit_table(members, "bogus")


@pytest.mark.parametrize(
    "outside,total,per",
    [
        ({"A": 5, "B": 1}, True, False),  # 6 sequences in 2 clades; only one clade has 5
        ({"A": 4, "B": 0}, False, False),
        ({"A": 3, "B": 2}, True, False),  # exactly 5 in 2 clades
        ({"A": 2, "B": 2}, False, False),  # 4 sequences: below 5
        ({"A": 9}, False, False),  # one clade only
        ({"A": 5, "B": 5}, True, True),
        ({"A": 5, "B": 4, "C": 5}, True, True),
        ({}, False, False),
    ],
)
def test_training_rule_boundaries(outside, total, per):
    outside = {k: v for k, v in outside.items() if v}
    assert mem.training_rules(outside) == (total, per)


def test_unit_needs_five_members_and_kiss_single_is_never_a_rule_unit():
    members = fixture_members({"A": 5, "B": 4, "Eurotiomycetes": 1, "C": 6})
    st = {u["clade_fine"]: u["status"] for u in mem.unit_table(members)}
    assert st["A"] == "unit" and st["C"] == "unit" and st["B"] == "too_few_members"
    assert st["Eurotiomycetes"] == "kiss_single"
    # nothing outside: not buildable, whatever the unit size
    only = fixture_members({"A": 6})
    assert mem.unit_table(only)[0]["status"] == "not_buildable"


def test_members_main_writes_ids_only(tmp_path):
    out = tmp_path / "o"
    rc = mem.main(
        [
            "--drift-hits",
            str(DRIFT / "drift_hits.tsv"),
            "--hits-manifest",
            str(DRIFT / "hits_manifest.tsv"),
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    ids = (out / "member_ids.txt").read_text().split()
    assert (
        len(ids) == 29 and "A0ACE4KFX7" not in ids and all(i.startswith("UniRef50_") for i in ids)
    )
    text = " ".join(
        (out / n).read_text() for n in ("members.tsv", "excluded.tsv", "units.tsv", "summary.txt")
    )
    assert not re.search(r"[ACDEFGHIKLMNPQRSTVWY]{40}", text)  # no sequence


# ------------------------------------------------------------------------------------------------ manifests


def clusters_for(members, groups):
    """groups: list of lists of fasta ids that share a cluster; every other id is its own cluster."""
    cl = {bld.fasta_id(m): bld.fasta_id(m) for m in members}
    for g in groups:
        for i in g:
            cl[i] = g[0]
    return cl


def test_manifest_held_out_never_in_training_clade_and_kiss_units():
    members = fixture_members({"A": 6, "B": 7, "Eurotiomycetes": 1})
    cl = clusters_for(members, [])
    for unit, n in (("clade:A", 6), ("clade:Eurotiomycetes", 1)):
        man = bld.make_manifest(members, unit, cl, "none")
        train = {t["id"] for t in man["training"]}
        assert not train & set(man["held_out_ids"]) and man["n_held_out"] == n
        assert man["status"] == "buildable"
    kiss = bld.make_manifest(members, "clade:Eurotiomycetes", cl, "none")
    assert kiss["held_out_ids"] == ["A0ACE4KFX7"] and all(
        t["id"] != "A0ACE4KFX7" for t in kiss["training"]
    )
    # the manifest passes the alignment check only if the held-out id and the neighbours are absent
    ids = [t["id"] for t in kiss["training"]]
    assert bld.check_alignment_ids(kiss, ids) == []
    assert bld.check_alignment_ids(kiss, ids + ["A0ACE4KFX7"])


def test_cluster30_removes_neighbours_and_can_empty_the_training_set():
    members = fixture_members({"A": 5, "B": 6, "Eurotiomycetes": 1})
    a0, b0 = (
        bld.fasta_id(members[0]),
        bld.fasta_id([m for m in members if m["clade_fine"] == "B"][0]),
    )
    cl = clusters_for(members, [[b0, a0]])
    man = bld.make_manifest(members, "clade:A", cl, "cluster30")
    assert b0 in man["excluded_30pct_neighbours"] and b0 not in {t["id"] for t in man["training"]}
    none = bld.make_manifest(members, "clade:A", cl, "none")
    assert b0 in {t["id"] for t in none["training"]} and none["n_excluded_neighbours"] == 0
    # everything outside A clusters with A: nothing left to train on
    allc = clusters_for(members, [[bld.fasta_id(m) for m in members]])
    dead = bld.make_manifest(members, "clade:A", allc, "cluster30")
    assert dead["status"] == "not_buildable" and dead["n_training"] == 0
    with pytest.raises(ValueError):
        bld.make_manifest(
            members, "clade:A", None, "cluster30"
        )  # cluster30 needs the cluster table
    with pytest.raises(ValueError):
        bld.make_manifest(members, "clade:Nope", cl, "none")
    with pytest.raises(ValueError):
        bld.make_manifest(members, "bogus:A", cl, "none")


def test_cluster_unit_holds_out_the_whole_cluster_and_has_one_column():
    members = fixture_members({"A": 5, "B": 6, "C": 5, "Eurotiomycetes": 1})
    ids = [bld.fasta_id(m) for m in members]
    grp = ids[:3]  # three members of A form one cluster
    cl = clusters_for(members, [grp])
    man = bld.make_manifest(members, f"cluster:{grp[0]}", cl, "cluster30")
    assert set(man["held_out_ids"]) == set(grp) and man["n_excluded_neighbours"] == 0
    same = bld.make_manifest(members, f"cluster:{grp[0]}", cl, "none")
    assert {t["id"] for t in same["training"]} == {t["id"] for t in man["training"]}
    rows = bld.cluster_units(members, cl)
    big = {r["cluster"]: r for r in rows}
    assert big[grp[0]]["n_members"] == 3 and big[grp[0]]["status"] == "too_few_members"
    q = [r for r in rows if r["has_query"] == "yes"][0]
    assert q["status"] == "kiss_single"


def test_task_list_dedupes_the_kiss_unit_and_gives_cluster_units_one_column():
    units = [
        {"unit": "clade:S", "status": "unit"},
        {"unit": "clade:E", "status": "kiss_single"},
        {"unit": "clade:L", "status": "too_few_members"},
    ]
    clus = [
        {"unit": "cluster:q", "status": "kiss_single"},  # the query alone: same unit as clade:E
        {"unit": "cluster:x", "status": "unit"},
        {"unit": "cluster:y", "status": "too_few_members"},
    ]
    t = bld.task_list(units, clus, "cluster30", "none")
    assert t == [
        ("clade:S", "cluster30"),
        ("clade:S", "none"),
        ("clade:E", "cluster30"),
        ("clade:E", "none"),
        ("cluster:x", "cluster30"),
        ("all", "none"),
    ]
    # without a kiss clade unit the query's own cluster is kept
    t2 = bld.task_list([{"unit": "clade:S", "status": "unit"}], clus, "cluster30", "none")
    assert ("cluster:q", "cluster30") in t2


# ------------------------------------------------------------------------------------------------ freeze


def frozen_fixture():
    fz = json.loads(PROPOSED.read_text())
    fz["frozen"] = True
    return fz


def test_the_proposed_freeze_is_a_proposal_and_matches_the_decisions():
    fz = json.loads(PROPOSED.read_text())
    assert fz["frozen"] is False and fz["owner_signature"] is None
    assert fz["extra_calls_per_10000"] == 5 and fz["cutoff_floor"] == 0.0
    rows = list(
        csv.DictReader(open(ROOT / "scripts/sorting_hat/refresh_3_10/inputs.tsv"), delimiter="\t")
    )
    sizes = {r["prot"]: int(r["n_proteins"]) for r in rows if r["prot"] in fz["kept_proteomes"]}
    assert sizes == fz["kept_proteome_sizes"] and sum(sizes.values()) == 32491
    assert math.isclose(e33b.extra_call_limit_total(sizes, 5), 16.2455)
    for k in e33b.FREEZE_KEYS + bld.EXTRA_FREEZE_KEYS:
        assert k in fz
    assert fz["comparators"]["jackhmmer_seed_ceiling"].startswith("jackhmmer")
    assert fz["heldout_truth"]["external_set"] is None  # an owner decision, not made here
    assert fz["evalue_db_size"] == 70198728  # the round 2 provenance
    prov = (ROOT / "docs/reports/data/pr6_gap_uniref50_2026-10-22/provenance.txt").read_text()
    assert "uniprot_fasta_sequences: 70198728" in prov


def test_freeze_guard_refuses_missing_proposed_and_incomplete(tmp_path):
    with pytest.raises(SystemExit):
        bld.load_freeze(tmp_path / "none.json")
    with pytest.raises(SystemExit):
        bld.load_freeze(PROPOSED)  # frozen false
    bad = frozen_fixture()
    del bad["comparators"]
    f = tmp_path / "f.json"
    f.write_text(json.dumps(bad))
    with pytest.raises(SystemExit):
        bld.load_freeze(f)
    bad = frozen_fixture()
    bad["training_rule"]["reading"] = "other"
    f.write_text(json.dumps(bad))
    with pytest.raises(SystemExit):
        bld.load_freeze(f)
    good = frozen_fixture()
    f.write_text(json.dumps(good))
    assert bld.load_freeze(f)["frozen"] is True


def test_no_manifest_task_list_or_cutoff_before_the_freeze(tmp_path):
    members = tmp_path / "m.tsv"
    members.write_text(
        "hit_id\tquery_accession\tclade_fine\tis_query\nA0ACE4KFX7\tA0ACE4KFX7\tEurotiomycetes\tyes\n"
    )
    for argv in (
        [
            "manifest",
            "--members",
            str(members),
            "--unit",
            "clade:Eurotiomycetes",
            "--mode",
            "none",
            "--freeze",
            str(PROPOSED),
            "--out",
            str(tmp_path / "o"),
        ],
        [
            "tasks",
            "--units",
            str(members),
            "--units-cluster",
            str(members),
            "--freeze",
            str(PROPOSED),
        ],
        [
            "cutoff",
            "--freeze",
            str(PROPOSED),
            "--proteome-hits",
            str(members),
            "--labelled",
            str(members),
            "--neg-groups",
            str(members),
            "--neg-hits",
            str(members),
            "--sizes",
            str(members),
            "--out",
            str(tmp_path / "c.json"),
        ],
    ):
        with pytest.raises(SystemExit):
            bld.main(argv)
    assert not (tmp_path / "o").exists() and not (tmp_path / "c.json").exists()


def test_reused_code_is_the_ecm33_code_and_the_cutoff_takes_no_positive_score():
    assert bld.lowest_cutoff is e33b.lowest_cutoff and bld.run_cutoff is e33b.run_cutoff
    import inspect

    assert "positive" not in " ".join(inspect.signature(bld.lowest_cutoff).parameters)


# ------------------------------------------------------------------------------------------------ comparators and report


def test_unit_row_counts_u_without_the_seed_ceiling():
    held = [f"h{i}" for i in range(10)]
    row = rep.unit_row(
        "clade:S",
        held,
        custom=set(held[:8]),
        pfam=set(held[:1]),
        relaxed={"h1"},
        seed={"h2"},
        jack=set(held),
        min_u=3,
        kiss=False,
    )
    assert (
        row["u_indep"] == 7 and row["u_with_ceiling"] == 0
    )  # the ceiling comparator would make u = 0
    assert row["custom_recovers_of_u"] == 5 and row["ceil_u_half"] == 4 and row["rule_ok"] == "yes"
    assert (
        row["testable_u"] == "yes"
        and row["label"] == "descriptive"
        and row["recall_custom"] == "8/10"
    )
    low = rep.unit_row("clade:S", held, set(held[:3]), set(), set(), set(), set(), 3, False)
    assert low["u_indep"] == 10 and low["rule_ok"] == "no"  # needs ceil(10/2) = 5 of the 10, has 3
    boundary = rep.unit_row("x", held, set(), set(held[:8]), set(), set(), set(), 3, False)
    assert boundary["u_indep"] == 2 and boundary["rule_ok"] == "not_testable"  # u < min_u


def test_kiss_unit_row_seed_in_unit_and_partial_label():
    row = rep.unit_row(
        "clade:Eurotiomycetes", ["A0ACE4KFX7"], {"A0ACE4KFX7"}, set(), set(), None, set(), 3, True
    )
    assert row["phmmer_seed_hit"] == "seed_in_unit" and row["label"] == "partial"
    assert row["u_indep"] == 1 and row["rule_ok"] == "not_testable"  # n = 1: never a rule unit
    assert row["recall_custom"] == "1/1"


def write_unit_dir(root, unit, mode, held, fz_sha, hits, seed=True):
    d = root / unit.replace(":", "_") / mode
    d.mkdir(parents=True)
    man = {
        "unit": unit,
        "mode": mode,
        "status": "buildable",
        "n_held_out": len(held),
        "held_out_ids": held,
    }
    (d / "manifest.json").write_text(json.dumps(man))
    (d / "cutoff.json").write_text(
        json.dumps(
            {
                "freeze_sha256": fz_sha,
                "positives_used": False,
                "cutoff": 20.0,
                "extra_calls_total": 4,
                "extra_call_limit_total": 16.2455,
            }
        )
    )
    (d / "hits_custom.tsv").write_text("id\tscore\n" + "".join(f"{i}\t{s}\n" for i, s in hits))
    (d / "hits_pfam_ga.tsv").write_text("id\n")
    (d / "hits_relaxed_panel.tsv").write_text("id\n")
    if seed:
        (d / "hits_phmmer_seed.tsv").write_text("id\tevalue\n")
        (d / "hits_jackhmmer_seed.tsv").write_text(
            "id\tevalue\n" + "".join(f"{i}\t1e-9\n" for i in held)
        )
    return d


def test_report_refuses_unfrozen_or_other_freeze_and_runs_end_to_end(tmp_path):
    fz = tmp_path / "fz.json"
    fz.write_text(json.dumps(frozen_fixture()))
    other = tmp_path / "other.json"
    other.write_text(json.dumps({**frozen_fixture(), "min_u": 4}))
    held = ["a", "b", "c", "d", "e"]
    res = tmp_path / "res"
    write_unit_dir(
        res, "clade:S", "cluster30", held, rep.sha256_file(fz), [("a", 90), ("b", 80), ("c", 10)]
    )
    write_unit_dir(
        res,
        "clade:Eurotiomycetes",
        "cluster30",
        ["A0ACE4KFX7"],
        rep.sha256_file(fz),
        [("A0ACE4KFX7", 70)],
        seed=False,
    )
    with pytest.raises(SystemExit):
        rep.main(["--results", str(res), "--freeze", str(PROPOSED), "--out", str(tmp_path / "x")])
    assert (
        rep.main(["--results", str(res), "--freeze", str(other), "--out", str(tmp_path / "y")]) == 1
    )  # made with another freeze
    assert (
        rep.main(["--results", str(res), "--freeze", str(fz), "--out", str(tmp_path / "out")]) == 0
    )
    rows = list(csv.DictReader(open(tmp_path / "out/recall_by_unit.tsv"), delimiter="\t"))
    s = [r for r in rows if r["unit"] == "clade:S"][0]
    assert (
        s["custom_hit"] == "2"
        and s["u_indep"] == "5"
        and s["u_with_ceiling"] == "0"
        and s["label"] == "descriptive"
    )
    k = [r for r in rows if r["unit"] == "clade:Eurotiomycetes"][0]
    assert k["label"] == "partial" and k["phmmer_seed_hit"] == "seed_in_unit"
    assert not any(r["label"] == "estimated" for r in rows)


# ------------------------------------------------------------------------------------------------ pre-check


def tbl_line(target, evalue, score, query):
    cols = [
        target,
        "-",
        query,
        "-",
        f"{evalue:g}",
        f"{score}",
        "0.1",
        f"{evalue:g}",
        f"{score}",
        "0.1",
    ]
    cols += ["1.0", "1", "0", "0", "1", "1", "1", "0", "desc"]
    return " ".join(cols)


def write_tbl(path, hits):
    path.write_text("#h\n" + "\n".join(tbl_line(*h) for h in hits) + "\n# [ok]\n")
    return path


def test_precheck_lists_family_members_and_panel_hits_separately(tmp_path):
    ctl = tmp_path / "controls.tsv"
    ctl.write_text(
        "accession\tgene\tspecies\tstratum\nN1\tg1\tsp1\tcell_wall\nN2\tg2\tsp2\tcell_wall\nN3\tg3\tsp3\tother\nN4\tg4\tsp4\tother\n"
    )
    fa = tmp_path / "n.faa"
    fa.write_text(">N1\nAAAA\n>N2\nCCCC\n>N3\nDDDD\n>N4\nEEEE\n")
    strict = write_tbl(tmp_path / "s.tbl", [("N1", 1e-8, 40.0, "A0ACE4KFX7")])
    weak = write_tbl(
        tmp_path / "w.tbl", [("N1", 1e-8, 40.0, "A0ACE4KFX7"), ("N2", 1e-4, 22.0, "A0ACE4KFX7")]
    )
    panel = write_tbl(
        tmp_path / "p.tbl", [("N1", 1e-9, 30.0, "Hydrophobin"), ("N3", 1e-6, 25.0, "HsbA")]
    )
    rel = write_tbl(
        tmp_path / "r.tbl", [("T1", 1e-6, 30.0, "A0ACE4KFX7"), ("T2", 1e-2, 11.0, "A0ACE4KFX7")]
    )
    out = tmp_path / "out"
    rc = pre.main(
        [
            "summarize",
            "--strict",
            str(strict),
            "--weak",
            str(weak),
            "--panel",
            str(panel),
            "--controls",
            str(ctl),
            "--negatives-fasta",
            str(fa),
            "--related",
            str(rel),
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    fam = list(csv.DictReader(open(out / "negatives_family_table.tsv"), delimiter="\t"))
    assert [r["accession"] for r in fam] == ["N1"]
    pan = list(csv.DictReader(open(out / "negatives_panel_table.tsv"), delimiter="\t"))
    assert {(r["accession"], r["model"], r["also_phmmer_family"]) for r in pan} == {
        ("N1", "Hydrophobin", "yes"),
        ("N3", "HsbA", "no"),
    }
    assert json.loads(
        "{" + (out / "proposed_negatives_ignored.json").read_text().strip().rstrip(",") + "}"
    ) == {"negatives_ignored": ["N1"]}
    assert [r["target"] for r in csv.DictReader(open(out / "related.tsv"), delimiter="\t")] == [
        "T1"
    ]  # T2 is above E 1e-3
    assert (
        "N3" in (out / "summary.txt").read_text()
        and "LOUD" not in (out / "summary.txt").read_text()
    )
    again = pre.main(
        [
            "summarize",
            "--strict",
            str(strict),
            "--weak",
            str(weak),
            "--panel",
            str(panel),
            "--controls",
            str(ctl),
            "--negatives-fasta",
            str(fa),
            "--out",
            str(out),
        ]
    )
    assert again == 1  # nothing is overwritten


def test_shared_negatives_have_no_hydrophobin_named_entry_but_mp1_has_hsba():
    """What reading controls.tsv can say (the job decides the rest): no name or protein_name contains hydrophobin/rodlet; B6QUR1 carries Pfam HsbA."""
    with open(ROOT / "data/controls/hard-negatives-shared/controls.tsv", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 161
    named = [
        r["accession"]
        for r in rows
        if re.search(r"hydrophobin|rodlet", r["gene"] + " " + r["protein_name"], re.I)
    ]
    assert named == []
    mp1 = [r for r in rows if r["accession"] == "B6QUR1"][0]
    assert "PF12296" in mp1["notes"] or "PF12296" in mp1["family_domain_pfam"]


# ------------------------------------------------------------------------------------------------ static checks


SCRIPTS = ["pr6_hyd1_extract.sbatch", "pr6_hyd1_build.sbatch", "pr6_hyd1_negatives_precheck.sbatch"]


@pytest.mark.parametrize("name", SCRIPTS)
def test_sbatch_is_valid_bash_marked_draft_and_submits_nothing(name):
    p = WF / name
    assert subprocess.run(["bash", "-n", str(p)]).returncode == 0
    text = p.read_text()
    assert "DRAFT, NOT SUBMITTED" in text and "set -euo pipefail" in text
    body = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))
    assert not re.search(r"(^|[;&|]\s*)(sbatch|srun|ssh|git push)\b|rm -rf /", body)
    assert "SCRATCH" in text and "/tmp" not in body


def test_build_job_orders_freeze_manifest_cutoff_before_the_held_out_searches():
    t = (WF / "pr6_hyd1_build.sbatch").read_text()
    i_freeze = t.index("grep -q '\"frozen\": true'")
    i_manifest = t.index("$BUILD manifest")
    i_cut = t.index("$BUILD cutoff")
    i_guard = t.index("the held-out step does not start")
    i_held = t.index('seqkit grep -f "$W/unit/heldout_ids.txt"')
    assert i_freeze < i_manifest < i_cut < i_guard < i_held
    assert "-Z" in t and "heldout.faa" in t
    assert "hydrophobin_relaxed.hmm" in t and "--nobias --cut_ga" in t
    assert 'IS_KISS" = no' in t  # the seed comparators are skipped in the Kiss unit
    assert "positive" not in t.split("# 4. The cutoff")[1].split("# 5.")[0].replace(
        "no positive score", ""
    )


def test_scripts_import_no_call_status_or_tier_code_and_write_no_committed_table():
    for n in ("members", "build", "report", "negatives_precheck"):
        src = (WF / f"pr6_hyd1_{n}.py").read_text()
        assert not re.search(
            r"^\s*(from|import)\s+(cellsurface_sorting_hat|adhesion_predict)", src, re.M
        )
        assert not re.search(
            r"data/sorting_hat|\bfamily_table\.tsv|family_measure",
            re.sub(r'""".*?"""', "", src, flags=re.S),
        )
