"""PR 6 ecm33 drafts: members and units, the training manifest, the independence check, the cutoff rule, the freeze guard, the report, and
static checks of the two sbatch scripts. Fixtures only; the real committed round 2 table is read for the member counts."""

import csv
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
REC = ROOT / "docs/reports/data/pr6_gap_uniref50_2026-10-22"


def load(name, path=None):
    spec = importlib.util.spec_from_file_location(name, path or WF / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


mem = load("pr6_ecm33_members")
bld = load("pr6_ecm33_build")
rep = load("pr6_ecm33_report")
relaxed = load("hydrophobin_relaxed", ROOT / "analysis/hydrophobin_truth/relaxed.py")

FREEZE_PROPOSED = WF / "pr6_ecm33_freeze_PROPOSED.json"


def hit(
    acc,
    clade,
    where="outside_clade",
    score="100",
    ev="1e-30",
    method="phmmer",
    query="Q5AGC4",
    taxid="1",
):
    return {
        "query": query,
        "method": method,
        "hit_accession": acc,
        "evalue": ev,
        "score": score,
        "taxid": taxid,
        "clade_fine": clade,
        "where": where,
    }


# --------------------------------------------------------------------------------------------------- members and units


def test_members_select_conservative_fungal_hits_with_a_clade_only():
    hits = [
        hit("UniRef50_Q5AGC4", "Saccharomycotina", "self"),
        hit("UniRef50_A", "Eurotiomycetes"),
        hit("UniRef50_B", "", "outside_clade"),  # fungal, no fine clade
        hit("UniRef50_C", "", "non_fungal"),
        hit("UniRef50_D", "Sordariomycetes", "inside_clade"),
    ]
    members, excluded = mem.select_members(hits)
    assert [m["hit_id"] for m in members] == ["UniRef50_D", "UniRef50_A", "UniRef50_Q5AGC4"] or len(
        members
    ) == 3
    assert {m["hit_id"] for m in members} == {"UniRef50_Q5AGC4", "UniRef50_A", "UniRef50_D"}
    assert [e["hit_id"] for e in excluded] == ["UniRef50_B"]
    assert "(no clade)" not in {m["clade_fine"] for m in members}
    # the non-fungal hit is never a member and never listed as an excluded fungal hit
    assert all(m["hit_id"] != "UniRef50_C" for m in members + excluded)


def test_the_no_clade_bucket_is_not_a_clade():
    hits = [hit(f"UniRef50_N{i}", "(no clade)") for i in range(5)] + [
        hit("UniRef50_E1", "Eurotiomycetes")
    ]
    members, excluded = mem.select_members(hits)
    assert len(members) == 1 and len(excluded) == 5
    units = mem.unit_table(members)
    assert [u["clade_fine"] for u in units] == ["Eurotiomycetes"]


def members_of(counts):
    return [{"clade_fine": c, "hit_id": f"{c}_{i}"} for c, n in counts.items() for i in range(n)]


def test_unit_boundaries_min_members_and_training_rule():
    # 5 members is testable, 4 is not; the rule needs 5 outside sequences in 2 clades
    u = {r["clade_fine"]: r for r in mem.unit_table(members_of({"A": 5, "B": 4, "C": 1}))}
    assert u["A"]["testable"] == "yes" and u["B"]["testable"] == "no"
    assert (
        u["A"]["n_outside_seq"] == 5
        and u["A"]["n_outside_clades"] == 2
        and u["A"]["training_rule_met"] == "yes"
    )
    assert u["A"]["status"] == "unit"
    assert (
        u["B"]["status"] == "too_few_members"
    )  # rule met (6 sequences in 2 clades), clade too small
    # exactly 4 outside sequences: not buildable
    v = {r["clade_fine"]: r for r in mem.unit_table(members_of({"A": 6, "B": 3, "C": 1}))}
    assert (
        v["A"]["n_outside_seq"] == 4
        and v["A"]["training_rule_met"] == "no"
        and v["A"]["status"] == "not_buildable"
    )
    # 5 outside sequences in ONE clade: not enough clades
    w = {r["clade_fine"]: r for r in mem.unit_table(members_of({"A": 6, "B": 5}))}
    assert w["A"]["training_rule_met"] == "no"


def test_real_round2_table_gives_the_documented_members():
    members, excluded = mem.select_members(mem.read_hits(REC / "gap_homologs.tsv.gz"))
    assert len(members) == 481 and len({m["clade_fine"] for m in members}) == 19
    assert len(excluded) == 5
    units = mem.unit_table(members)
    assert sum(1 for u in units if u["status"] == "unit") == 11
    assert sum(1 for u in units if u["status"] == "too_few_members") == 8
    sacch = next(u for u in units if u["clade_fine"] == "Saccharomycotina")
    assert (
        sacch["n_members"] == 107 and sacch["status"] == "unit"
    )  # 106 inside-clade hits plus the query itself


def test_members_main_writes_no_sequence(tmp_path):
    assert mem.main(["--homologs", str(REC / "gap_homologs.tsv.gz"), "--out", str(tmp_path)]) == 0
    for f in tmp_path.iterdir():
        assert not re.search(r"[ACDEFGHIKLMNPQRSTVWY]{40,}", f.read_text()), f.name
    ids = (tmp_path / "member_ids.txt").read_text().split()
    assert (
        "UniRef50_Q5AGC4" not in ids and len(ids) == 480
    )  # the query is read from the Kiss control set, not from UniRef50


# --------------------------------------------------------------------------------------------------- manifest


def clustered_members():
    ms = []
    for clade, ids in {
        "X": ["x1", "x2"],
        "Y": ["y1", "y2", "y3"],
        "Z": ["z1", "z2", "z3", "z4"],
    }.items():
        ms += [
            {"clade_fine": clade, "hit_id": i, "query_accession": "", "is_query": "no"} for i in ids
        ]
    return ms


def test_manifest_held_out_never_in_training():
    clusters = {i["hit_id"]: i["hit_id"] for i in clustered_members()}  # every member alone
    man = bld.make_manifest(clustered_members(), "X", clusters, "cluster30")
    assert set(man["held_out_ids"]) == {"x1", "x2"}
    assert not ({t["id"] for t in man["training"]} & set(man["held_out_ids"]))
    assert man["n_training"] == 7 and man["status"] == "buildable"


def test_manifest_excludes_30pct_neighbours_of_the_held_out_members():
    clusters = {i["hit_id"]: i["hit_id"] for i in clustered_members()}
    clusters["y1"] = "x1"  # y1 clusters with the held-out x1
    clusters["z4"] = "x2"
    man = bld.make_manifest(clustered_members(), "X", clusters, "cluster30")
    assert man["excluded_30pct_neighbours"] == ["y1", "z4"]
    train = {t["id"] for t in man["training"]}
    assert not (train & {"y1", "z4", "x1", "x2"})
    none = bld.make_manifest(clustered_members(), "X", clusters, "none")
    assert none["excluded_30pct_neighbours"] == [] and none["n_training"] == 7


def test_manifest_not_buildable_when_the_independence_rule_empties_the_training_set():
    clusters = {
        i["hit_id"]: "x1" for i in clustered_members()
    }  # one conserved cluster spanning every clade
    man = bld.make_manifest(clustered_members(), "X", clusters, "cluster30")
    assert (
        man["status"] == "not_buildable"
        and man["n_training"] == 0
        and "minimum 5 in 2" in man["reason"]
    )


def test_manifest_all_has_no_held_out_unit_and_unknown_unit_fails():
    man = bld.make_manifest(clustered_members(), "all", None, "none")
    assert man["held_out_ids"] == [] and man["n_training"] == 9
    with pytest.raises(ValueError):
        bld.make_manifest(clustered_members(), "Nope", {}, "none")


def test_query_member_uses_its_accession_as_the_fasta_id():
    m = {
        "hit_id": "UniRef50_Q5AGC4",
        "query_accession": "Q5AGC4",
        "is_query": "yes",
        "clade_fine": "S",
    }
    assert bld.fasta_id(m) == "Q5AGC4"
    assert (
        bld.fasta_id({"hit_id": "UniRef50_X", "query_accession": "", "is_query": "no"})
        == "UniRef50_X"
    )


def test_alignment_check_flags_held_out_neighbours_and_gaps():
    clusters = {i["hit_id"]: i["hit_id"] for i in clustered_members()}
    clusters["y1"] = "x1"
    man = bld.make_manifest(clustered_members(), "X", clusters, "cluster30")
    ok = [t["id"] for t in man["training"]]
    assert bld.check_alignment_ids(man, ok) == []
    assert any("held-out id x1" in p for p in bld.check_alignment_ids(man, ok + ["x1"]))
    assert any("neighbour y1" in p for p in bld.check_alignment_ids(man, ok + ["y1"]))
    assert any("missing" in p for p in bld.check_alignment_ids(man, ok[1:]))
    assert any("not in the manifest" in p for p in bld.check_alignment_ids(man, ok + ["q9"]))


def test_occupancy_check_boundaries():
    aln = ["ACDEFGHIKL", "ACDEFGHIK-", "ACDEFGH---", "AC--------"]
    r = bld.occupancy_check(aln, occupancy=0.5, min_columns=9, min_median_coverage=0.6)
    assert (
        r["n_columns"] == 9
    )  # columns 0 to 8 hold a residue in at least 2 of 4 sequences (the limit is inclusive)
    assert r["ok"] and r["median_coverage"] == pytest.approx((7 / 9 + 1.0) / 2)
    assert not bld.occupancy_check(aln, 0.5, 10, 0.6)["ok"]
    assert not bld.occupancy_check(aln, 0.5, 9, 0.95)["ok"]
    assert not bld.occupancy_check([], 0.5, 1, 0.1)["ok"]


# --------------------------------------------------------------------------------------------------- cutoff


def test_cutoff_rule_equals_the_hydrophobin_rule_and_takes_no_positive_score():
    cands = [
        {"id": f"p{i}", "proteome": "P1" if i % 2 else "P2", "score": float(s)}
        for i, s in enumerate([3, 8, 12, 15, 21, 40])
    ]
    sizes = {"P1": 20000, "P2": 10000}
    hn = {"G1": [0.0, 5.0, 30.0, 0.0], "G2": [0.0] * 9 + [14.0]}
    for per10k, rate in ((5, 0.05), (1, 0.3), (0, 1.0)):
        a = bld.lowest_cutoff(cands, sizes, hn, per10k, rate, 0.0, ["G2"])
        b = relaxed.lowest_cutoff(cands, sizes, hn, per10k, rate, 0.0, ["G2"])
        assert a == b
    params = bld.lowest_cutoff.__code__.co_varnames[: bld.lowest_cutoff.__code__.co_argcount]
    assert "positives" not in " ".join(params) and "pos" not in params


def test_extra_call_limit_arithmetic_over_the_kept_proteomes():
    sizes = {"Afum_Af293_UniProt": 9647, "Calb_SC5314": 6212, "Cimm_RS": 9910, "Scer_S288C": 6722}
    assert sum(sizes.values()) == 32491
    assert bld.extra_call_limit_total(sizes, 5) == pytest.approx(16.2455)
    assert math.floor(bld.extra_call_limit_total(sizes, 5) * 10) / 10 == 16.2


def frozen_fixture():
    fz = json.loads(FREEZE_PROPOSED.read_text())
    fz["frozen"] = True
    return fz


def test_run_cutoff_excludes_labelled_proteins_and_ignored_negatives():
    fz = frozen_fixture()
    sizes = {"Afum_Af293_UniProt": 9647, "Calb_SC5314": 6212, "Cimm_RS": 9910, "Scer_S288C": 6722}
    hits = [
        {"proteome": "Calb_SC5314", "id": "ecm33_true", "score": "400"},
        {"proteome": "Calb_SC5314", "id": "junk", "score": "12"},
    ]
    neg_groups = {"P38248": "gpi_wall_enzyme", "n1": "gpi_wall_enzyme", "n2": "gpi_wall_enzyme"}
    neg_hits = {"P38248": 350.0}  # the true Ecm33 among the negatives
    res = bld.run_cutoff(fz, hits, {"ecm33_true"}, neg_groups, neg_hits, sizes)
    assert res["positives_used"] is False
    assert "P38248" in res["negatives_ignored"]
    # the true Ecm33 (labelled) and P38248 (ignored) do not push the cutoff up to 400 or 350: only the unlabelled 12 matters
    assert res["cutoff"] <= 12.0 + 1e-9
    assert res["extra_calls_by_proteome"]["Calb_SC5314"] <= 1


def test_freeze_guard_refuses_missing_unfrozen_and_incomplete_files(tmp_path):
    with pytest.raises(SystemExit):
        bld.load_freeze(tmp_path / "nope.json")
    with pytest.raises(SystemExit):
        bld.load_freeze(FREEZE_PROPOSED)  # frozen: false
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"frozen": True}))
    with pytest.raises(SystemExit):
        bld.load_freeze(bad)
    good = tmp_path / "good.json"
    good.write_text(json.dumps(frozen_fixture()))
    assert bld.load_freeze(good)["extra_calls_per_10000"] == 5


def test_the_proposed_freeze_matches_the_inputs_and_keeps_the_decided_items():
    fz = json.loads(FREEZE_PROPOSED.read_text())
    assert (
        fz["frozen"] is False and fz["owner_signature"] is None
    )  # the OWNER freezes it, not this work
    assert (
        fz["extra_calls_per_10000"] == 5
        and fz["hard_negative_rate"] == 0.05
        and fz["cutoff_floor"] == 0.0
    )
    rows = list(
        csv.DictReader(open(ROOT / "scripts/sorting_hat/refresh_3_10/inputs.tsv"), delimiter="\t")
    )
    sizes = {r["prot"]: int(r["n_proteins"]) for r in rows if r["prot"] in fz["kept_proteomes"]}
    assert sizes == fz["kept_proteome_sizes"] and sum(sizes.values()) == 32491
    # decided 2026-10-22: P38248 left out of the rate, both independence modes, min_u 3
    assert fz["negatives_ignored"] == ["P38248"]
    assert fz["independence_mode"] == "cluster30" and fz["independence_second_column"] == "none"
    assert fz["min_u"] == 3
    assert set(fz["decided_by_owner_2026-10-22"]) == {
        "negatives_ignored",
        "independence_mode",
        "independence_second_column",
        "min_u",
    }
    assert "OWNER DECISION NEEDED" not in FREEZE_PROPOSED.read_text()
    # the unchanged proposals stay proposals and keep their values
    assert set(fz["proposals_unchanged"]) == {
        "cutoff_floor",
        "extra_calls_per_10000",
        "hard_negative_rate",
        "min_unit_members",
        "alignment_check",
        "training_rule",
    }
    assert fz["alignment_check"] == {
        "occupancy": 0.5,
        "min_columns": 100,
        "min_median_coverage": 0.6,
        "note": fz["alignment_check"]["note"],
    }
    assert fz["min_unit_members"] == 5
    for k in bld.FREEZE_KEYS:
        assert k in fz


def test_freeze_requires_two_different_independence_modes(tmp_path):
    bad = frozen_fixture()
    bad["independence_second_column"] = "cluster30"
    f = tmp_path / "f.json"
    f.write_text(json.dumps(bad))
    with pytest.raises(SystemExit):
        bld.load_freeze(f)
    bad["independence_second_column"] = "weird"
    f.write_text(json.dumps(bad))
    with pytest.raises(SystemExit):
        bld.load_freeze(f)


def test_both_modes_keep_the_held_out_unit_out_of_the_alignment():
    """The held-out-absent-from-alignment manifest check holds for BOTH modes (cluster30 headline, none second column)."""
    clusters = {i["hit_id"]: i["hit_id"] for i in clustered_members()}
    clusters["y1"] = "x1"
    for mode in ("cluster30", "none"):
        man = bld.make_manifest(clustered_members(), "X", clusters, mode)
        ok = [t["id"] for t in man["training"]]
        assert bld.check_alignment_ids(man, ok) == []
        assert any("held-out id x1" in p for p in bld.check_alignment_ids(man, ok + ["x1"]))
        assert not ({"x1", "x2"} & set(ok))
    c30 = bld.make_manifest(clustered_members(), "X", clusters, "cluster30")
    none = bld.make_manifest(clustered_members(), "X", clusters, "none")
    assert c30["n_training"] == none["n_training"] - 1 and c30["excluded_30pct_neighbours"] == [
        "y1"
    ]
    assert any(
        "neighbour y1" in p
        for p in bld.check_alignment_ids(c30, [t["id"] for t in c30["training"]] + ["y1"])
    )
    assert (
        bld.check_alignment_ids(none, [t["id"] for t in none["training"]]) == []
    )  # y1 is legitimate training data in mode none


def test_array_size_is_two_tasks_per_unit_plus_the_all_model():
    members, _ = mem.select_members(mem.read_hits(REC / "gap_homologs.tsv.gz"))
    n_units = sum(1 for u in mem.unit_table(members) if u["status"] == "unit")
    assert n_units == 11 and 2 * n_units + 1 == 23  # tasks 0..22


def test_p38248_is_a_shared_negative_and_a_true_ecm33():
    rows = list(
        csv.DictReader(
            open(ROOT / "data/controls/hard-negatives-shared/controls.tsv"), delimiter="\t"
        )
    )
    p = next(r for r in rows if r["accession"] == "P38248")
    assert p["label"] == "negative" and p["gene"] == "ECM33"


# --------------------------------------------------------------------------------------------------- report


def test_unit_row_arithmetic():
    held = [f"h{i}" for i in range(10)]
    r = rep.unit_row(
        "U",
        held,
        custom={"h0", "h1", "h2", "h3", "h4", "h5"},
        pfam={"h0"},
        relaxed={"h1"},
        it={"h2"},
        min_u=3,
    )
    # comparators hit h0,h1,h2; missed by all = h3..h9 (7); custom recovers h3,h4,h5 of them
    assert r["u_missed_by_every_comparator"] == 7 and r["custom_recovers_of_u"] == 3
    assert r["ceil_u_half"] == 4 and r["rule_ok"] == "no" and r["testable_u"] == "yes"
    assert r["status"] == "descriptive"
    z = rep.unit_row("U", held, custom=set(held), pfam=set(held), relaxed=set(), it=set(), min_u=3)
    assert (
        z["u_missed_by_every_comparator"] == 0 and z["rule_ok"] == "not_testable"
    )  # the hydrophobin lesson: u = 0


def write_unit(
    d,
    fz_sha,
    held,
    custom_scores,
    cutoff=50.0,
    status="buildable",
    n_excl=0,
    n_train=300,
    p38248=None,
    neg_extra=None,
):
    d.mkdir(parents=True)
    (d / "manifest.json").write_text(
        json.dumps(
            {
                "status": status,
                "held_out_ids": held,
                "n_held_out": len(held),
                "n_excluded_neighbours": n_excl,
                "n_training": n_train,
            }
        )
    )
    base = {
        "freeze_sha256": fz_sha,
        "positives_used": False,
        "extra_calls_total": 2,
        "extra_call_limit_total": 16.2455,
        "extra_calls_by_proteome": {"P": 2},
        "hard_negative_call_rates": {"G": 0.0},
        "hard_negative_rate_limit": 0.05,
    }
    (d / "cutoff.json").write_text(json.dumps({**base, "cutoff": cutoff}))
    (d / "relaxed_cutoff.json").write_text(json.dumps({**base, "cutoff": 10.0}))
    (d / "hits_custom.tsv").write_text(
        "id\tscore\n" + "".join(f"{i}\t{s}\n" for i, s in custom_scores.items())
    )
    (d / "hits_pfam_ga.tsv").write_text("id\n")
    (d / "hits_pfam_relaxed.tsv").write_text("id\tscore\n")
    (d / "hits_iter.tsv").write_text("id\tevalue\n")
    (d / "neg_hits.tsv").write_text(
        "accession\tscore\n"
        + (f"P38248\t{p38248}\n" if p38248 is not None else "")
        + "".join(f"{a}\t{v}\n" for a, v in (neg_extra or {}).items())
    )
    (d / "neg_hits_relaxed.tsv").write_text("accession\tscore\n")


def units_file(tmp_path, units):
    p = tmp_path / "units.tsv"
    p.write_text(
        "clade_fine\tn_members\tn_outside_seq\tn_outside_clades\ttraining_rule_met\ttestable\tstatus\n"
        + "".join(f"{u}\t9\t30\t5\tyes\tyes\tunit\n" for u in units)
    )
    return p


def freeze_file(tmp_path):
    fz = tmp_path / "fz.json"
    fz.write_text(json.dumps(frozen_fixture()))
    return fz


def test_report_refuses_without_a_frozen_freeze_or_with_a_different_one(tmp_path):
    held = ["h0", "h1", "h2", "h3"]
    fz = freeze_file(tmp_path)
    fz_sha = rep.sha256_file(fz)
    res = tmp_path / "res"
    for mode in ("cluster30", "none"):
        write_unit(res / "U" / mode, fz_sha, held, {"h0": 60.0, "h1": 70.0})
    args = [
        "--results",
        str(res),
        "--units",
        str(units_file(tmp_path, ["U"])),
        "--out",
        str(tmp_path / "out"),
    ]
    with pytest.raises(SystemExit):
        rep.main([*args, "--freeze", str(tmp_path / "missing.json")])
    unfrozen = tmp_path / "unfrozen.json"
    unfrozen.write_text(FREEZE_PROPOSED.read_text())
    with pytest.raises(SystemExit):
        rep.main([*args, "--freeze", str(unfrozen)])
    assert rep.main([*args, "--freeze", str(fz)]) == 0
    other = tmp_path / "other.json"  # a cutoff made with another freeze file is refused
    other.write_text(json.dumps({**frozen_fixture(), "min_u": 7}))
    assert rep.main([*args, "--freeze", str(other)]) == 1


def test_report_end_to_end_both_modes_side_by_side_descriptive(tmp_path):
    held = ["h0", "h1", "h2", "h3"]
    fz = freeze_file(tmp_path)
    sha = rep.sha256_file(fz)
    res = tmp_path / "res"
    write_unit(
        res / "U" / "cluster30",
        sha,
        held,
        {"h0": 60.0, "h1": 70.0, "h2": 20.0},
        n_excl=40,
        n_train=260,
        p38248=350.0,
    )
    write_unit(
        res / "U" / "none",
        sha,
        held,
        {"h0": 60.0, "h1": 70.0, "h2": 90.0, "h3": 99.0},
        n_excl=0,
        n_train=300,
        p38248=20.0,
    )
    out = tmp_path / "out"
    assert (
        rep.main(
            [
                "--results",
                str(res),
                "--units",
                str(units_file(tmp_path, ["U"])),
                "--freeze",
                str(fz),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    rows = list(csv.DictReader(open(out / "recall_by_unit.tsv"), delimiter="\t"))
    assert [(r["independence_mode"], r["role"]) for r in rows] == [
        ("cluster30", "headline"),
        ("none", "descriptive"),
    ]
    assert [r["custom_hit"] for r in rows] == ["2", "4"] and all(
        r["status"] == "descriptive" for r in rows
    )
    wide = next(csv.DictReader(open(out / "side_by_side.tsv"), delimiter="\t"))
    assert (
        wide["n_excluded_neighbours_headline"] == "40"
        and wide["n_training_headline"] == "260"
        and wide["n_training_second"] == "300"
    )
    assert (
        wide["custom_hit_headline"] == "2"
        and wide["custom_hit_second"] == "4"
        and wide["second_is_same_training_set"] == "no"
    )
    assert sorted(p.name for p in out.iterdir()) == [
        "extra_calls.tsv",
        "hard_negative_calls.tsv",
        "ignored_negatives_line.tsv",
        "recall_by_unit.tsv",
        "side_by_side.tsv",
        "summary.txt",
    ]
    assert not any(c in (out / "recall_by_unit.tsv").read_text() for c in ("estimated", "smoke"))
    # the freeze lists only P38248: one line per mode; hit above the cutoff of 50 in the headline, below it in the second column
    pl = list(csv.DictReader(open(out / "ignored_negatives_line.tsv"), delimiter="\t"))
    assert [(r["independence_mode"], r["custom_hit"]) for r in pl] == [
        ("cluster30", "yes"),
        ("none", "no"),
    ]
    # and it is not in any hard-negative rate row
    assert all(
        r["group"] != "P38248"
        for r in csv.DictReader(open(out / "hard_negative_calls.tsv"), delimiter="\t")
    )


def test_report_writes_one_line_per_ignored_negative(tmp_path):
    """negatives_ignored with two accessions (P38248 and Q4WNS8): each has its own line per mode, none is in a rate row."""
    held = ["h0", "h1", "h2", "h3"]
    fzd = frozen_fixture()
    fzd["negatives_ignored"] = ["P38248", "Q4WNS8"]
    fz = tmp_path / "fz2.json"
    fz.write_text(json.dumps(fzd))
    sha = rep.sha256_file(fz)
    res = tmp_path / "res"
    write_unit(
        res / "U" / "cluster30",
        sha,
        held,
        {"h0": 60.0},
        p38248=350.0,
        neg_extra={"Q4WNS8": 20.0},
    )
    write_unit(
        res / "U" / "none",
        sha,
        held,
        {"h0": 60.0},
        p38248=20.0,
        neg_extra={"Q4WNS8": 400.0},
    )
    out = tmp_path / "out"
    assert (
        rep.main(
            [
                "--results",
                str(res),
                "--units",
                str(units_file(tmp_path, ["U"])),
                "--freeze",
                str(fz),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    lines = list(csv.DictReader(open(out / "ignored_negatives_line.tsv"), delimiter="\t"))
    assert [(r["independence_mode"], r["accession"], r["custom_hit"]) for r in lines] == [
        ("cluster30", "P38248", "yes"),
        ("cluster30", "Q4WNS8", "no"),
        ("none", "P38248", "no"),
        ("none", "Q4WNS8", "yes"),
    ]
    summary = (out / "summary.txt").read_text()
    assert "P38248, Q4WNS8" in summary and "ignored_negatives_line.tsv" in summary
    # the rate rows are per stratum: no accession of the ignored list appears in them
    assert not any(
        r["group"] in ("P38248", "Q4WNS8")
        for r in csv.DictReader(open(out / "hard_negative_calls.tsv"), delimiter="\t")
    )


def test_run_cutoff_leaves_every_ignored_negative_out_of_the_rate():
    fz = frozen_fixture()
    fz["negatives_ignored"] = ["P38248", "Q4WNS8"]
    sizes = {"Afum_Af293_UniProt": 9647, "Calb_SC5314": 6212, "Cimm_RS": 9910, "Scer_S288C": 6722}
    hits = [{"proteome": "Calb_SC5314", "id": "junk", "score": "12"}]
    groups = {a: "gpi_wall_enzyme" for a in ("P38248", "Q4WNS8", "n1", "n2")}
    res = bld.run_cutoff(fz, hits, set(), groups, {"P38248": 350.0, "Q4WNS8": 300.0}, sizes)
    assert res["negatives_ignored"] == ["P38248", "Q4WNS8"]
    assert res["cutoff"] <= 12.0 + 1e-9  # neither true Ecm33 pushes the cutoff up


def test_report_second_mode_not_rebuilt_repeats_the_headline(tmp_path):
    held = ["h0", "h1", "h2", "h3"]
    fz = freeze_file(tmp_path)
    sha = rep.sha256_file(fz)
    res = tmp_path / "res"
    write_unit(res / "U" / "cluster30", sha, held, {"h0": 60.0}, n_excl=0, n_train=300)
    (res / "U" / "none").mkdir(parents=True)
    (res / "U" / "none" / "SAME_AS_HEADLINE").write_text("same training set\n")
    out = tmp_path / "out"
    assert (
        rep.main(
            [
                "--results",
                str(res),
                "--units",
                str(units_file(tmp_path, ["U"])),
                "--freeze",
                str(fz),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    wide = next(csv.DictReader(open(out / "side_by_side.tsv"), delimiter="\t"))
    assert (
        wide["second_is_same_training_set"] == "yes"
        and wide["custom_hit_headline"] == wide["custom_hit_second"] == "1"
    )


def test_report_not_buildable_headline_is_listed_not_scored(tmp_path):
    held = ["h0", "h1", "h2", "h3"]
    fz = freeze_file(tmp_path)
    sha = rep.sha256_file(fz)
    res = tmp_path / "res"
    d = res / "U" / "cluster30"
    d.mkdir(parents=True)
    (d / "manifest.json").write_text(
        json.dumps(
            {
                "status": "not_buildable",
                "held_out_ids": ["h"],
                "n_held_out": 1,
                "n_training": 0,
                "n_excluded_neighbours": 9,
            }
        )
    )
    write_unit(
        res / "U" / "none", sha, held, {"h0": 60.0, "h1": 70.0}, n_excl=0, n_train=300
    )  # the lighter rule still builds
    out = tmp_path / "out"
    assert (
        rep.main(
            [
                "--results",
                str(res),
                "--units",
                str(units_file(tmp_path, ["U"])),
                "--freeze",
                str(fz),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    rows = list(csv.DictReader(open(out / "recall_by_unit.tsv"), delimiter="\t"))
    assert rows[0]["rule_ok"] == "not_buildable" and rows[0]["role"] == "headline"
    assert rows[1]["custom_hit"] == "2" and rows[1]["role"] == "descriptive"


def test_scripts_import_no_call_status_or_tier_code_and_do_not_write_the_committed_tables():
    for name in ("pr6_ecm33_members.py", "pr6_ecm33_build.py", "pr6_ecm33_report.py"):
        text = (WF / name).read_text()
        assert "cellsurface_sorting_hat" not in text.replace("data/sorting_hat", "")
        assert (
            "family_measure.tsv" not in re.sub(r'"""[\s\S]*?"""', "", text)
            or "does not touch" in text
        )


# --------------------------------------------------------------------------------------------------- static checks of the jobs


def code_lines(text):
    return "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))


def test_extract_job_static_checks():
    s = (WF / "pr6_ecm33_extract.sbatch").read_text()
    assert (
        "DRAFT, NOT SUBMITTED" in s
        and "${SCRATCH:?" in s
        and "/tmp" not in code_lines(s).replace("$TMP", "")
    )
    assert "#SBATCH -p batch" in s and "seqkit grep" in s and "28781193344" in s
    assert 'fail "$OUT exists' in s and "check_pin" in s
    assert code_lines(s).count("seqkit grep") == 1  # UniRef50 is streamed once


def test_build_job_orders_cutoff_before_the_held_out_searches():
    s = (WF / "pr6_ecm33_build.sbatch").read_text()
    assert "DRAFT, NOT SUBMITTED" in s and "${SCRATCH:?" in s and "--array" in s
    assert "mafft --localpair --maxiterate 1000" in s
    assert "grep -q '\"frozen\": true'" in s
    i_cut = s.index("$BUILD cutoff")
    i_held = s.index(
        'hmmsearch -T "$FLOOR" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/h_custom.tbl"'
    )
    assert i_cut < i_held
    assert "no cutoff: the held-out step does not start" in s
    assert s.index("check-alignment") < s.index("hmmbuild")
    assert s.index("mmseqs easy-cluster") < s.index("$BUILD manifest")
    assert "--min-seq-id 0.3 -c 0.5" in s
    # both independence modes: a (unit, mode) task list, the headline manifest first, the second column skipped when identical
    assert "tasks.list" in s and "HEAD_MODE" in s and "SECOND_MODE" in s and "SAME_AS_HEADLINE" in s
    assert 'UDIR="$OUT/$UNIT/$MODE"' in s and "printf 'all\\tnone\\n'" in s
    assert s.count("$BUILD manifest") == 2 and "--mode" in s
    assert "/tmp" not in code_lines(s).replace("$TMP", "")
    assert "family_measure" not in s and "family_table" not in s.replace("data/sorting_hat", "")


def test_sbatch_files_are_valid_bash():
    import subprocess

    for n in ("pr6_ecm33_extract.sbatch", "pr6_ecm33_build.sbatch"):
        assert subprocess.run(["bash", "-n", str(WF / n)], capture_output=True).returncode == 0
