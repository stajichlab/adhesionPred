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


def test_the_proposed_freeze_matches_the_inputs_and_is_a_proposal():
    fz = json.loads(FREEZE_PROPOSED.read_text())
    assert fz["frozen"] is False and fz["owner_signature"] is None
    assert fz["extra_calls_per_10000"] == 5 and fz["hard_negative_rate"] == 0.05
    rows = list(
        csv.DictReader(open(ROOT / "scripts/sorting_hat/refresh_3_10/inputs.tsv"), delimiter="\t")
    )
    sizes = {r["prot"]: int(r["n_proteins"]) for r in rows if r["prot"] in fz["kept_proteomes"]}
    assert sizes == fz["kept_proteome_sizes"] and sum(sizes.values()) == 32491
    assert (
        "P38248" in fz["negatives_ignored"]
        and "OWNER DECISION NEEDED" in fz["negatives_ignored_note"]
    )
    for k in bld.FREEZE_KEYS:
        assert k in fz


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


def write_unit(d, fz_sha, held, custom_scores, cutoff=50.0, status="buildable"):
    d.mkdir(parents=True)
    (d / "manifest.json").write_text(
        json.dumps({"status": status, "held_out_ids": held, "n_held_out": len(held)})
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


def units_file(tmp_path, units):
    p = tmp_path / "units.tsv"
    p.write_text(
        "clade_fine\tn_members\tn_outside_seq\tn_outside_clades\ttraining_rule_met\ttestable\tstatus\n"
        + "".join(f"{u}\t9\t30\t5\tyes\tyes\tunit\n" for u in units)
    )
    return p


def test_report_refuses_without_a_frozen_freeze_or_with_a_different_one(tmp_path):
    held = ["h0", "h1", "h2", "h3"]
    fz = tmp_path / "fz.json"
    fz.write_text(json.dumps(frozen_fixture()))
    fz_sha = rep.sha256_file(fz)
    res = tmp_path / "res"
    write_unit(res / "U", fz_sha, held, {"h0": 60.0, "h1": 70.0})
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
    # a cutoff made with another freeze file: refused
    other = tmp_path / "other.json"
    other.write_text(json.dumps({**frozen_fixture(), "min_u": 7}))
    assert rep.main([*args, "--freeze", str(other)]) == 1


def test_report_end_to_end_is_descriptive_and_writes_no_call_or_tier(tmp_path):
    held = ["h0", "h1", "h2", "h3"]
    fz = tmp_path / "fz.json"
    fz.write_text(json.dumps(frozen_fixture()))
    res = tmp_path / "res"
    write_unit(
        res / "U", rep.sha256_file(fz), held, {"h0": 60.0, "h1": 70.0, "h2": 20.0}
    )  # h2 is under the cutoff of 50
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
    row = next(csv.DictReader(open(out / "recall_by_unit.tsv"), delimiter="\t"))
    assert row["custom_hit"] == "2" and row["n_held_out"] == "4" and row["status"] == "descriptive"
    assert row["u_missed_by_every_comparator"] == "4"
    assert sorted(p.name for p in out.iterdir()) == [
        "extra_calls.tsv",
        "hard_negative_calls.tsv",
        "recall_by_unit.tsv",
        "summary.txt",
    ]
    assert not any(c in (out / "recall_by_unit.tsv").read_text() for c in ("estimated", "smoke"))
    assert (
        next(csv.DictReader(open(out / "extra_calls.tsv"), delimiter="\t"))["within_limit"] == "yes"
    )


def test_report_not_buildable_unit_is_listed_not_scored(tmp_path):
    fz = tmp_path / "fz.json"
    fz.write_text(json.dumps(frozen_fixture()))
    d = tmp_path / "res" / "U"
    d.mkdir(parents=True)
    (d / "manifest.json").write_text(
        json.dumps({"status": "not_buildable", "held_out_ids": ["h"], "n_held_out": 1})
    )
    out = tmp_path / "out"
    assert (
        rep.main(
            [
                "--results",
                str(tmp_path / "res"),
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
    assert (
        next(csv.DictReader(open(out / "recall_by_unit.tsv"), delimiter="\t"))["rule_ok"]
        == "not_buildable"
    )


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
    assert "/tmp" not in code_lines(s).replace("$TMP", "")
    assert "family_measure" not in s and "family_table" not in s.replace("data/sorting_hat", "")


def test_sbatch_files_are_valid_bash():
    import subprocess

    for n in ("pr6_ecm33_extract.sbatch", "pr6_ecm33_build.sbatch"):
        assert subprocess.run(["bash", "-n", str(WF / n)], capture_output=True).returncode == 0
