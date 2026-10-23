"""Independent-set evaluation of the frozen custom Ecm33 model (drafts): test-set selection (real candidates.tsv), the cluster30 separation, u and
the rule at their boundaries, the freeze guard, the sha256 discipline of the report, static checks of the two sbatch scripts and of the proposed
freeze. Fixtures only except for the committed tables; nothing is run on HPCC."""

import csv
import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
CAND = ROOT / "docs/reports/data/ecm33_independent_set_candidates/candidates.tsv"
BLD = ROOT / "docs/reports/data/pr6_ecm33_build_2026-10-23"
MODEL = ROOT / "data/custom_models/ecm33_all_2026-10-23"


def load(name):
    spec = importlib.util.spec_from_file_location(name, WF / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


mem = load("pr6_ecm33_eval_members")
rep = load("pr6_ecm33_eval_report")
PROPOSED = WF / "pr6_ecm33_eval_freeze_PROPOSED.json"


def tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


# ------------------------------------------------------------------------------------------------------ selection (real table)


def test_selection_counts_are_recounted_from_candidates():
    rows = mem.read_tsv(CAND)
    prot = mem.select_proteins(rows)
    clu = mem.cluster_table(prot)
    cn = mem.counts(prot, clu)
    assert (cn["n_proteins_all"], cn["n_clusters_all"]) == (725, 187)
    assert (cn["n_proteins_pe"], cn["n_clusters_pe"]) == (563, 109)
    assert cn["tiers_proteins"] == {"T2_named": 321, "T3_panther_only": 404}
    assert cn["subfamily_proteins"] == {"SF3": 718, "SF12": 7}
    assert cn["subfamily_clusters"]["SF12"] == 4
    assert sum(cn["clades_clusters"].values()) == 187
    # every selected protein is outside every training cluster
    assert all(p["overlap_with_training"] == "new_cluster" for p in prot)


def test_purity_check_counts_164_of_481():
    rows, summ = mem.purity(mem.read_tsv(CAND), mem.read_tsv(BLD / "members/members.tsv"))
    assert summ["n_members"] == 481
    assert summ["n_not_in_panther_cluster"] == 164


def test_cluster_attribute_rules_on_a_fixture():
    def p(acc, cl, clade, tier, sf, e):
        return {
            "accession": acc, "uniref50": cl, "clade_fine": clade, "tier": tier,
            "panther_subfamily": sf, "panther_evalue": e, "overlap_with_training": "new_cluster",
        }  # fmt: skip

    prot = [
        p("a", "UniRef50_X", "Zed", "T3_panther_only", "SF3", "1e-5"),
        p("b", "UniRef50_X", "Alpha", "T2_named", "SF12", "1e-12"),
        p(
            "c", "UniRef50_Y", "Alpha", "T3_panther_only", "SF3", "1e-10"
        ),  # exactly at the cut: strong
        p("d", "UniRef50_Z", "Alpha", "T3_panther_only", "SF3", "1.1e-10"),
    ]
    cl = {c["cluster_id"]: c for c in mem.cluster_table(prot)}
    x = cl["UniRef50_X"]
    assert (x["clade_fine"], x["clade_multi"]) == ("Alpha", "yes")  # tie: alphabetical
    assert (x["tier"], x["tier_multi"]) == ("T2_named", "yes")  # best tier
    assert x["subfamily"] == "SF12" and x["headline_pe"] == "yes" and x["n_proteins_pe"] == 1
    assert cl["UniRef50_Y"]["headline_pe"] == "yes"
    assert cl["UniRef50_Z"]["headline_pe"] == "no"


# ------------------------------------------------------------------------------------------------------------- separation


def test_separation_excludes_a_test_sequence_that_clusters_with_a_training_member():
    cluster_of = {"T1": "T1", "T2": "T1", "S1": "T1", "S2": "S2", "S3": "S2"}
    rows = {
        r["cluster_id"]: r for r in mem.separation(cluster_of, ["T1", "T2"], ["S1", "S2", "S3"])
    }
    assert rows["S1"]["separated"] == "no" and rows["S1"]["training_ids"] == "T1,T2"
    assert rows["S2"]["separated"] == "yes" and rows["S3"]["separated"] == "yes"


def test_separation_refuses_a_missing_or_shared_id():
    with pytest.raises(SystemExit):
        mem.separation(
            {"T1": "T1"}, ["T1"], ["S1"]
        )  # S1 was never clustered: it would pass silently
    with pytest.raises(SystemExit):
        mem.separation({"T1": "T1"}, ["T1"], ["T1"])


# ------------------------------------------------------------------------------------------------------------ u and the rule


def hits_with(custom, others=()):
    h = {c: set() for c in rep.COMPARATORS}
    h["custom"] = set(custom)
    for name, ids in others:
        h[name] = set(ids)
    return h


@pytest.mark.parametrize(
    "u,rec,testable,ok",
    [
        (2, 2, "no", "not_testable"),  # below min_u: never testable, even with every u recovered
        (3, 1, "yes", "no"),  # ceil(3/2) = 2 needed
        (3, 2, "yes", "yes"),
        (4, 2, "yes", "yes"),  # ceil(4/2) = 2
        (4, 1, "yes", "no"),
        (5, 3, "yes", "yes"),  # ceil(5/2) = 3
        (5, 2, "yes", "no"),
    ],
)
def test_u_and_the_rule_at_the_boundaries(u, rec, testable, ok):
    ids = {f"c{i}" for i in range(u + 4)}
    in_u = sorted(ids)[:u]
    comparator_hit = sorted(ids)[u:]  # these four are found by some comparator
    h = hits_with(
        in_u[:rec] + comparator_hit,
        [("pfam_ga", comparator_hit[:2]), ("jackhmmer_seed_only", comparator_hit[2:])],
    )
    g = rep.group_row("headline", "all", "all", ids, h, 3, 10)
    assert g["u"] == u and g["custom_recovers_of_u"] == rec
    assert g["ceil_u_half"] == (math.ceil(u / 2) if u else 0)
    assert (g["testable_u"], g["rule_ok"]) == (testable, ok)


@pytest.mark.parametrize(
    "k,n,lo,hi",
    [
        (0, 3, 0.0, 0.708),
        (3, 3, 0.292, 1.0),
        (2, 3, 0.094, 0.992),
        (1, 2, 0.013, 0.987),
        (5, 5, 0.478, 1.0),
    ],
)
def test_clopper_pearson_known_values(k, n, lo, hi):
    got = rep.clopper_pearson(k, n)
    assert round(got[0], 3) == lo and round(got[1], 3) == hi
    if k == 0:
        assert got[0] == 0.0
    if k == n:
        assert got[1] == 1.0
    assert rep.fmt_ci(k, n) == f"{lo:.3f}-{hi:.3f}"
    assert rep.fmt_ci(0, 0) == "NA"


def test_verdict_applies_only_to_the_whole_headline_row():
    ids = {f"c{i}" for i in range(6)}
    h = hits_with({"c0", "c1"})
    assert rep.group_row("headline", "all", "all", ids, h, 3, 1)["rule_ok"] == "no"  # u 6 needs 3
    for col, kind in (("headline", "clade"), ("descriptive_1", "all"), ("descriptive_2", "tier")):
        assert rep.group_row(col, kind, "x", ids, h, 3, 1)["rule_ok"] == "descriptive"


def test_u_counts_only_clusters_missed_by_every_one_of_the_five_comparators():
    ids = {"a", "b", "c", "d", "e", "f"}
    h = hits_with(
        set(),
        [
            ("pfam_ga", {"a"}),
            ("pfam_relaxed", {"b"}),
            ("iteration_alignment_seeded", {"c"}),
            ("phmmer_seed_only", {"d"}),
            ("jackhmmer_seed_only", {"e"}),
        ],
    )
    g = rep.group_row("headline", "all", "all", ids, h, 3, 6)
    assert g["u"] == 1 and g["any_comparator"] == 5  # only f; the custom model is not a comparator


def test_threshold_application_is_inclusive_and_fixed():
    rows = [{"id": "a", "score": "29.6"}, {"id": "b", "score": "29.59"}]
    assert set(rep.by_score(rows, 29.6)) == {"a"}
    ev = [{"id": "a", "evalue": "1e-5"}, {"id": "b", "evalue": "1.1e-5"}]
    assert set(rep.by_evalue(ev)) == {"a"}


def test_groups_split_headline_from_descriptive_and_keep_sf12_apart():
    def c(i, clade, tier, sf, pe, n=2):
        return {
            "cluster_id": i, "clade_fine": clade, "tier": tier, "subfamily": sf, "n_proteins": n,
            "best_panther_evalue": "1e-20" if pe else "1e-3", "headline_pe": "yes" if pe else "no",
        }  # fmt: skip

    cl = [
        c("a", "A", "T2_named", "SF3", True),
        c("b", "A", "T3_panther_only", "SF12", True),
        c("c", "B", "T3_panther_only", "SF3", True),
        c("d", "B", "T3_panther_only", "SF3", False),
    ]
    sep = {"a": True, "b": True, "c": False, "d": True}
    rows = rep.build_groups(cl, sep, hits_with({"a", "b"}), 3)
    pick = {(r["column"], r["group_type"], r["group"]): r for r in rows}
    assert (
        pick[("headline", "all", "all")]["n_clusters"] == 2
    )  # c is not separated, d is not strong
    assert pick[("descriptive_1", "all", "all")]["n_clusters"] == 3
    assert pick[("descriptive_2", "all", "all")]["n_clusters"] == 4
    assert pick[("headline", "subfamily", "SF12")]["n_clusters"] == 1
    assert ("headline", "clade", "B") not in pick
    assert all(r["status"] in ("headline", "descriptive") for r in rows)


# ------------------------------------------------------------------------------------------------------ freeze and sha256


def test_load_frozen_refuses_missing_unfrozen_and_unsigned(tmp_path):
    with pytest.raises(SystemExit):
        rep.load_frozen(tmp_path / "nope.json")
    with pytest.raises(SystemExit):
        rep.load_frozen(PROPOSED)  # frozen: false
    fz = json.load(open(PROPOSED))
    fz["frozen"] = True
    p = tmp_path / "f.json"
    p.write_text(json.dumps(fz))
    with pytest.raises(SystemExit):
        rep.load_frozen(p)  # no owner_signature
    fz["owner_signature"] = "Test Owner, 2026-10-24"
    p.write_text(json.dumps(fz))
    assert rep.load_frozen(p)["u"]["min_u"] == 3


def make_results(tmp_path, fz_path, hits=None, tamper=False, freeze_sha=None):
    """A complete small result directory + inputs for rep.main, with a 4-cluster fixture freeze."""
    fz = json.load(open(PROPOSED))
    fz.update(frozen=True, owner_signature="Test Owner")
    fz["test_set"]["expected_clusters"] = 4
    fz["test_set"]["expected_proteins"] = 8
    cj = MODEL / "cutoff.json"
    fz_path.write_text(json.dumps(fz))
    res = tmp_path / "res"
    res.mkdir()
    hits = hits or {}
    names = rep.HIT_FILES
    for k, fn in names.items():
        rows = hits.get(k, [])
        if k == "pfam_ga":
            tsv(res / fn, ["id"], [[r] for r in rows])
        elif k in ("custom", "pfam_relaxed"):
            tsv(res / fn, ["id", "score"], rows)
        else:
            tsv(res / fn, ["id", "evalue"], rows)
    files = {fn: rep.sha256_file(res / fn) for fn in names.values()}
    if tamper:
        (res / "hits_custom.tsv").write_text("id\tscore\nUniRef50_A\t99\n")
    (res / "run.json").write_text(
        json.dumps(
            {"freeze_sha256": freeze_sha or rep.sha256_file(fz_path), "commit": "x", "files": files}
        )
    )
    clusters = tmp_path / "clusters.tsv"
    tsv(
        clusters,
        mem.CLUSTER_HEADER,
        [
            [f"UniRef50_{x}", "Alpha", "no", "T2_named", "no", "SF3", 2, "1e-20", "yes", 2]
            for x in "ABCD"
        ],
    )
    sep = tmp_path / "sep.tsv"
    tsv(sep, mem.SEP_HEADER, [[f"UniRef50_{x}", "r", 0, "", "yes"] for x in "ABCD"])
    return res, clusters, sep, cj


def run_main(tmp_path, res, clusters, sep, cj, fz_path):
    out = tmp_path / "out"
    rc = rep.main(
        [
            "--results",
            str(res),
            "--clusters",
            str(clusters),
            "--separation",
            str(sep),
            "--freeze",
            str(fz_path),
            "--cutoff-json",
            str(cj),
            "--out",
            str(out),
        ]
    )
    return rc, out


def test_report_end_to_end_on_a_fixture(tmp_path):
    fz = tmp_path / "freeze.json"
    res, clusters, sep, cj = make_results(
        tmp_path,
        fz,
        hits={"custom": [["UniRef50_A", "50"], ["UniRef50_B", "29.5"]], "pfam_ga": ["UniRef50_C"]},
    )
    rc, out = run_main(tmp_path, res, clusters, sep, cj, fz)
    assert rc == 0
    g = {
        (r["column"], r["group_type"], r["group"]): r
        for r in csv.DictReader(open(out / "groups.tsv"), delimiter="\t")
    }
    top = g[("headline", "all", "all")]
    assert (
        top["custom"] == "1" and top["pfam_ga"] == "1" and top["u"] == "3"
    )  # B is below 29.6; A, C, D... u = A, B, D
    assert top["custom_recovers_of_u"] == "1" and top["rule_ok"] == "no"  # u 3 needs 2
    assert top["ci_recovered_of_u"] == rep.fmt_ci(1, 3) and top["ci_recall_custom"] == rep.fmt_ci(
        1, 4
    )
    summ = (out / "summary.txt").read_text()
    assert "VERDICT (whole headline set, min_u 3): FAIL" in summ
    assert "wide row is not a result" in summ
    assert all("CI" in ln for ln in summ.splitlines() if ln.startswith("["))
    assert sum(ln.startswith("[") for ln in summ.splitlines()) == len(g)
    # no call, status, tier or family_measure file is written
    assert sorted(p.name for p in out.iterdir()) == [
        "clusters.tsv",
        "groups.tsv",
        "not_separated.tsv",
        "summary.txt",
        "u_clusters.tsv",
    ]


def test_report_refuses_a_result_with_a_different_freeze_sha(tmp_path):
    fz = tmp_path / "freeze.json"
    res, clusters, sep, cj = make_results(tmp_path, fz, freeze_sha="0" * 64)
    rc, out = run_main(tmp_path, res, clusters, sep, cj, fz)
    assert rc == 1 and not out.exists()


def test_report_refuses_a_result_file_that_changed_after_the_run(tmp_path):
    fz = tmp_path / "freeze.json"
    res, clusters, sep, cj = make_results(tmp_path, fz, tamper=True)
    rc, out = run_main(tmp_path, res, clusters, sep, cj, fz)
    assert rc == 1 and not out.exists()


def test_report_refuses_without_a_signed_freeze(tmp_path):
    res, clusters, sep, cj = make_results(tmp_path, tmp_path / "freeze.json")
    with pytest.raises(SystemExit):
        run_main(tmp_path, res, clusters, sep, cj, PROPOSED)  # the proposal is never accepted
    with pytest.raises(SystemExit):
        run_main(tmp_path, res, clusters, sep, cj, tmp_path / "absent.json")


def test_report_refuses_a_hit_that_is_not_a_test_cluster(tmp_path):
    fz = tmp_path / "freeze.json"
    res, clusters, sep, cj = make_results(tmp_path, fz, hits={"custom": [["UniRef50_ZZZ", "99"]]})
    rc, _ = run_main(tmp_path, res, clusters, sep, cj, fz)
    assert rc == 1


# --------------------------------------------------------------------------------------------------- the proposed freeze


def test_proposed_freeze_matches_the_committed_model_and_tables():
    fz = json.load(open(PROPOSED))
    assert fz["frozen"] is False and fz["owner_signature"] == ""
    cj = json.load(open(MODEL / "cutoff.json"))
    assert fz["custom_model"]["cutoff_bits"] == cj["cutoff"] == 29.6
    sums = dict(reversed(line.split()) for line in open(MODEL / "SHA256SUMS"))
    assert fz["custom_model"]["sha256"] == sums["model.hmm"] == rep.sha256_file(MODEL / "model.hmm")
    assert fz["custom_model"]["cutoff_json_sha256"] == rep.sha256_file(MODEL / "cutoff.json")
    rel = json.load(open(BLD / "ecm33_build/all/none/relaxed_cutoff.json"))
    assert fz["searches"]["relaxed_pfam_cutoff_bits"] == rel["cutoff"] == 15.0
    assert fz["searches"]["evalue_db_size"] == 70198728
    prov = (BLD / "ecm33_build/all/none/provenance.txt").read_text()
    assert fz["separation"]["training_members_faa_sha256"] in prov
    assert fz["u"]["min_u"] == 3 and fz["report_only"] is True
    assert fz["u"]["comparators"] == [
        "pfam_ga", "pfam_relaxed", "iteration_alignment_seeded", "phmmer_seed_only", "jackhmmer_seed_only"
    ]  # fmt: skip
    ts = fz["test_set"]
    prot = mem.select_proteins(mem.read_tsv(CAND))
    cn = mem.counts(prot, mem.cluster_table(prot), ts["headline_panther_evalue"])
    assert (
        (ts["expected_proteins"], ts["expected_clusters"])
        == (cn["n_proteins_all"], cn["n_clusters_all"])
        == (725, 187)
    )
    assert (
        (ts["expected_headline_pe_proteins"], ts["expected_headline_pe_clusters"])
        == (cn["n_proteins_pe"], cn["n_clusters_pe"])
        == (563, 109)
    )
    assert fz["proposals"] and fz["decided_by_owner_2026-10-22_and_23"]


def test_proposed_freeze_records_the_owner_decisions_of_2026_10_23():
    fz = json.load(open(PROPOSED))
    assert fz["frozen"] is False
    assert len(fz["decided_by_owner_2026-10-23"]) == 5
    assert "-Z 70198728" in fz["searches"]["iteration_alignment_seeded"]
    assert "NO -Z" not in json.dumps(fz)
    assert not any("without -Z" in p for p in fz["proposals"])
    assert "Clopper-Pearson" in fz["u"]["interval"] and fz["u"]["min_u"] == 3


def test_the_signed_freeze_does_not_exist_yet():
    assert not (WF / "pr6_ecm33_eval_freeze.json").exists()  # only the owner creates it


# ------------------------------------------------------------------------------------------------- scripts and no sequences


@pytest.mark.parametrize("name", ["pr6_ecm33_eval_extract.sbatch", "pr6_ecm33_eval_search.sbatch"])
def test_sbatch_static_checks(name):
    text = (WF / name).read_text()
    assert "DRAFT, NOT SUBMITTED" in text
    assert "-p batch" in text and "--time=" in text and "--mem=" in text and "NOT MEASURED" in text
    assert "pr6_ecm33_eval_freeze.json" in text and '"frozen": true' in text
    assert (
        "git push" not in text and "\nsbatch " not in text
    )  # the submit command only appears in comments
    if shutil.which("bash"):
        assert subprocess.run(["bash", "-n", str(WF / name)]).returncode == 0


def test_search_script_has_the_frozen_options():
    t = (WF / "pr6_ecm33_eval_search.sbatch").read_text()
    assert (
        '-Z "$EVALDB"' in t and t.count('-Z "$EVALDB"') == 3
    )  # phmmer, jackhmmer seed-only and the alignment-seeded jackhmmer (owner decision 2026-10-23)
    assert "70198728" in t
    assert (
        'jackhmmer -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5 -Z "$EVALDB" --cpu "$THREADS" --noali -o /dev/null --tblout "$W/iter.tbl" "$W/train.sto"'
        in t
    )
    for line in t.splitlines():  # every iterative or phmmer command line carries -Z
        if line.startswith(("jackhmmer ", "phmmer ")):
            assert '-Z "$EVALDB"' in line, line
    assert (
        "hmmsearch --cut_ga" in t and "-T 29.6" not in t
    )  # the report applies the frozen cutoff from the freeze
    assert "training members are NOT in any target" in t


def test_extract_script_uses_the_freeze_separation_command():
    t = (WF / "pr6_ecm33_eval_extract.sbatch").read_text()
    assert "--min-seq-id 0.3 -c 0.5" in t and "MMseqs2/17-b804f" in t and "seqkit grep" in t


def test_committed_outputs_carry_no_sequences(tmp_path):
    out = tmp_path / "sel"
    assert mem.main(["select", "--candidates", str(CAND), "--out", str(out)]) == 0
    for f in out.iterdir():
        for line in f.read_text().splitlines():
            assert not line.startswith(">")
            assert not re.search(r"[A-Z]{40,}", line), f"{f.name} holds a sequence-like string"
    # no FASTA is added to the repository by this work
    assert not list(WF.glob("*eval*.fa*")) and not list(WF.glob("*eval*.faa"))
    assert (out / "test_cluster_ids.txt").read_text().count("\n") == 187
    assert (out / "headline_pe_cluster_ids.txt").read_text().count("\n") == 109
