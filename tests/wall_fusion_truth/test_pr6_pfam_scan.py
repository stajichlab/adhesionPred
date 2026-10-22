"""PR 6 step 6.0 drafts: pr6_pfam_scan_summary.py (fixtures, hand counts) and pr6_pfam_scan.sbatch (static checks)."""

import csv
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
SBATCH = WF / "pr6_pfam_scan.sbatch"
BASE = ROOT / "docs/reports/data/kiss_baseline_2026-10-22"
REC = ROOT / "docs/reports/data/kiss_cluster_split_2026-10-22"
GROUPS = ROOT / "docs/superpowers/specs/2026-10-22-pr6-custom-hmms-data/candidate_groups.tsv"


def load():
    spec = importlib.util.spec_from_file_location(
        "pr6_pfam_scan_summary", WF / "pr6_pfam_scan_summary.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pr6_pfam_scan_summary"] = mod
    spec.loader.exec_module(mod)
    return mod


MISSED_HEADER = (
    "slice\tcluster_id\tn_members\taccessions\tgenes\tspecies\tclades_fine\tstrata_fine\t"
    "seed_class_slice\tseed_class_all39\tbest_identity_all39\tjgi_accept\n"
)


def missed_row(slice_, cid, accs, gene="g"):
    return (
        f"{slice_}\t{cid}\t{len(accs)}\t{','.join(accs)}\t{gene}\tSpecies x\tClade\tclade_only\t"
        "seed-independent\tseed-independent\t\trow_verified\n"
    )


def dom(model, acc, tlen, seq, qlen, ievalue, score, a_from, a_to):
    """One hmmscan --domtblout line (target = model, query = sequence)."""
    return (
        f"{model} {acc} {tlen} {seq} - {qlen} 1e-20 {score} 0.0 1 1 1e-20 {ievalue} {score} 0.0 "
        f"1 {tlen} {a_from} {a_to} {a_from} {a_to} 0.99 description of the model\n"
    )


def domtbl(rows, option):
    return f"# hmmscan :: search sequence(s) against a profile database\n# Option settings: {option}\n#\n{''.join(rows)}# [ok]\n"


@pytest.fixture()
def fx(tmp_path):
    (tmp_path / "missed.tsv").write_text(
        MISSED_HEADER
        + missed_row("wall", "KS1", ["M1"])
        + missed_row("wall", "KS2", ["M2"])
        + missed_row("adhesion", "KS3", ["M3"])
        + missed_row("fusion", "KS4", ["M4a", "M4b"])
        + missed_row("adhesion", "KS5", ["M5a", "M5b"])
    )
    (tmp_path / "groups.tsv").write_text(
        "cluster_id\tslice\tgroup\n"
        "KS1\twall\tW1\nKS2\twall\tW1\nKS3\tadhesion\tA1\nKS4\tfusion\tF1\nKS5\tadhesion\tA1\n"
    )
    (tmp_path / "wall.tsv").write_text("id\tname\tset\nPF00001\tWallOne\twall_remodeling\n")
    (tmp_path / "family.tsv").write_text(
        "pfam_acc\tname\tactive\nPF05730\tCFEM\tyes\nPF99999\tInactive\tno\n"
    )
    ga = domtbl(
        [
            dom("CFEM", "PF05730.20", 80, "M1", 100, 1e-30, 90.0, 5, 70),  # a panel model
            dom("Other1", "PF90001.5", 200, "M2", 100, 1e-10, 40.0, 1, 30),
            dom("Other1", "PF90001.5", 200, "M2", 100, 1e-8, 35.0, 41, 70),  # union 60 of 100
            dom("Other2", "PF90002.7", 100, "M3", 100, 1e-6, 30.0, 1, 30),  # cover 0.30
            dom("Other4", "PF90004.1", 100, "M5a", 100, 1e-6, 30.0, 1, 20),  # cover 0.20
            dom("Other4", "PF90004.1", 100, "M5b", 100, 1e-9, 50.0, 1, 60),  # cover 0.60
        ],
        "hmmscan --cut_ga --domtblout x db q",
    )
    weak = domtbl(
        [
            dom("Other1", "PF90001.5", 200, "M2", 100, 1e-10, 40.0, 1, 30),
            dom("Weak3", "PF90003.2", 100, "M4b", 100, 5e-4, 12.0, 10, 50),  # no GA hit: a weak hit
            dom("CFEM", "PF05730.20", 80, "M1", 100, 1e-30, 90.0, 5, 70),
        ],
        "hmmscan -E 1e-3 --domE 1e-3 --domtblout y db q",
    )
    (tmp_path / "ga.txt").write_text(ga)
    (tmp_path / "weak.txt").write_text(weak)
    return tmp_path


def run(mod, fx, out="out", weak=True, extra=None):
    args = [
        "--domtbl-ga",
        str(fx / "ga.txt"),
        "--missed",
        str(fx / "missed.tsv"),
        "--groups",
        str(fx / "groups.tsv"),
        "--wall-models",
        str(fx / "wall.tsv"),
        "--family-table",
        str(fx / "family.tsv"),
        "--out",
        str(fx / out),
    ]
    if weak:
        args += ["--domtbl-weak", str(fx / "weak.txt")]
    return mod.main(args + (extra or []))


def tsv(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_classes_proposals_and_counts_match_the_hand_count(fx):
    mod = load()
    assert run(mod, fx) == 0
    rows = {r["cluster_id"]: r for r in tsv(fx / "out/scan_per_cluster.tsv")}
    assert {k: (v["class"], v["proposal"]) for k, v in rows.items()} == {
        "KS1": ("hit_panel_model", "already_in_panel"),
        "KS2": ("hit_other_pfam", "add_to_panel_proposal"),
        "KS3": ("hit_other_pfam", "partial_domain_for_owner"),
        "KS4": ("no_pfam_hit", "custom_target"),
        "KS5": ("hit_other_pfam", "add_to_panel_proposal"),  # the best member counts
    }
    assert "frac=0.600" in rows["KS2"]["other_pfam_families"]
    assert "frac=0.300" in rows["KS3"]["other_pfam_families"]
    assert (
        "frac=0.600" in rows["KS5"]["other_pfam_families"]
        and "seq=M5b" in rows["KS5"]["other_pfam_families"]
    )
    counts = {(r["slice"], r["group"]): r for r in tsv(fx / "out/scan_counts.tsv")}
    assert counts[("wall", "ALL")]["n_clusters"] == "2"
    assert (
        counts[("wall", "ALL")]["hit_panel_model"],
        counts[("wall", "ALL")]["hit_other_pfam"],
    ) == ("1", "1")
    assert (
        counts[("adhesion", "ALL")]["hit_other_pfam"],
        counts[("adhesion", "ALL")]["no_pfam_hit"],
    ) == ("2", "0")
    assert counts[("fusion", "F1")]["no_pfam_hit"] == "1"
    # only the cluster without a GA hit is a custom target; a weak hit does not change its class
    targets = tsv(fx / "out/custom_targets.tsv")
    assert [t["cluster_id"] for t in targets] == ["KS4"]
    assert targets[0]["weak_families"] == "PF90003:Weak3"
    # the weak hit of M2 on a family that has a GA hit for the same member is not listed as weak
    assert rows["KS2"]["weak_families"] == ""
    prop = tsv(fx / "out/panel_proposals.tsv")
    assert [(p["pfam_acc"], p["n_clusters"]) for p in prop] == [("PF90001", "1"), ("PF90004", "1")]
    assert "proposal_for_owner" in prop[0]["status"]


def test_a_panel_hit_is_reported_loudly(fx, capsys):
    mod = load()
    run(mod, fx)
    err = capsys.readouterr().err
    assert "REPORT LOUD: 1 missed cluster(s) are hit by a panel model" in err
    assert "REPORT LOUD" in (fx / "out/summary.txt").read_text()
    assert [r["cluster_id"] for r in tsv(fx / "out/panel_hits.tsv")] == ["KS1"]


def test_no_loud_report_without_a_panel_hit(fx, capsys):
    mod = load()
    lines = [ln for ln in (fx / "ga.txt").read_text().splitlines(keepends=True) if "CFEM" not in ln]
    (fx / "ga.txt").write_text("".join(lines))
    run(mod, fx)
    assert "REPORT LOUD" not in capsys.readouterr().err
    assert "REPORT LOUD" not in (fx / "out/summary.txt").read_text()
    assert tsv(fx / "out/panel_hits.tsv") == []


def test_covered_fraction_is_the_union_of_the_intervals():
    mod = load()
    assert mod.covered_fraction([(1, 30), (41, 70)], 100) == pytest.approx(0.6)
    assert mod.covered_fraction([(1, 50), (30, 60)], 100) == pytest.approx(
        0.6
    )  # overlap counted once
    assert mod.covered_fraction([(1, 10), (11, 20)], 100) == pytest.approx(0.2)  # adjacent
    assert mod.covered_fraction([(1, 200)], 100) == 1.0 and mod.covered_fraction([(1, 5)], 0) == 0.0


@pytest.mark.parametrize(
    "bad",
    ["no_trailer", "ga_without_cut_ga", "hmmsearch_orientation", "unknown_query"],
)
def test_ga_table_refusals(fx, bad):
    mod = load()
    text = (fx / "ga.txt").read_text()
    if bad == "no_trailer":
        text = text.replace("# [ok]\n", "")
    elif bad == "ga_without_cut_ga":
        text = text.replace("--cut_ga", "-E 1e-5")
    elif bad == "hmmsearch_orientation":
        text = text.replace(" PF90002.7 ", " - ")
    else:
        text = text.replace(" M3 ", " NOTAMEMBER ")
    (fx / "ga.txt").write_text(text)
    with pytest.raises((ValueError, SystemExit)):
        run(mod, fx)


def test_weak_table_made_with_cut_ga_is_refused(fx):
    mod = load()
    (fx / "weak.txt").write_text(
        (fx / "weak.txt").read_text().replace("hmmscan -E", "hmmscan --cut_ga -E")
    )
    with pytest.raises(ValueError, match="must not be made with --cut_ga"):
        run(mod, fx)


def test_missing_group_and_existing_output_are_refused(fx):
    mod = load()
    (fx / "groups.tsv").write_text("cluster_id\tslice\tgroup\nKS1\twall\tW1\n")
    with pytest.raises(SystemExit, match="no function group"):
        run(mod, fx)
    (fx / "groups.tsv").write_text(
        "cluster_id\tslice\tgroup\n" + "".join(f"KS{i}\twall\tW1\n" for i in range(1, 6))
    )
    (fx / "out").mkdir()
    (fx / "out/keep.txt").write_text("x")
    with pytest.raises(SystemExit, match="not empty"):
        run(mod, fx)


def test_outputs_are_deterministic_and_hold_no_sequences(fx):
    mod = load()
    run(mod, fx, out="a")
    run(mod, fx, out="b")
    names = sorted(p.name for p in (fx / "a").iterdir())
    assert names == sorted(mod.OUTPUTS)
    for n in names:
        assert (fx / "a" / n).read_bytes() == (fx / "b" / n).read_bytes()
        assert not re.search(r"^[ACDEFGHIKLMNPQRSTVWY]{30,}$", (fx / "a" / n).read_text(), re.M)


def test_sbatch_follows_the_job_conventions():
    assert subprocess.run(["bash", "-n", str(SBATCH)], capture_output=True).returncode == 0
    text = SBATCH.read_text()
    code = "\n".join(
        ln
        for ln in text.splitlines()
        if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")
    )
    assert text.startswith("#!/bin/bash -l")
    assert "DRAFT, NOT SUBMITTED" in text.splitlines()[1]
    assert "set -euo pipefail" in code and "#SBATCH -p short" in code
    assert '"${SCRATCH:?' in code
    assert "/tmp" not in "\n".join(ln for ln in code.splitlines() if not ln.startswith("#SBATCH"))
    assert "BASH_SOURCE" not in code
    assert "module load hmmer/3.4" in code
    assert "hmmscan --cut_ga" in code and "-E 1e-3 --domE 1e-3" in code
    searches = [ln.strip() for ln in code.splitlines() if ln.strip().startswith("hmmscan ")]
    assert (
        len(searches) == 2 and "--cut_ga" in searches[0] and "--cut_ga" not in searches[1]
    )  # the weak table has no GA cut-off
    assert '[ ! -e "$OUT" ]' in code and "check_pin" in code and "sidecar_ok" in code
    assert (
        'Pfam${PFAM_RELEASE}/"' in code
        and "4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e" in code
    )
    for ext in ("h3f", "h3i", "h3m", "h3p"):
        assert ext in code
    assert (
        "hard_negative" not in text
    )  # no negative is scanned: this step is about the 44 missed clusters only


def test_the_pinned_member_count_matches_the_committed_record():
    rows = tsv(BASE / "clusters_missed.tsv")
    members = {a for r in rows for a in r["accessions"].split(",")}
    n = re.search(r"^N_MISSED_MEMBERS=(\d+)", SBATCH.read_text(), re.M)
    assert len(rows) == 44 and len(members) == int(n.group(1)) == 66
    in_split = {r["accession"]: r for r in tsv(REC / "members_split.tsv")}
    assert members <= set(in_split) and not any("#" in a for a in members)
    # every missed member is a positive control, never a negative copy
    assert all(in_split[a]["class"] == "positive" for a in members)
    groups = {r["cluster_id"] for r in tsv(GROUPS)}
    assert {r["cluster_id"] for r in rows} <= groups
