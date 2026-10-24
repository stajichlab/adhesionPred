"""The committed run record of the PR 6 ecm33 independent-set evaluation."""

import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_ecm33_independent_eval_2026-10-24"
REPORT = ROOT / "docs/reports/2026-10-24-pr6-ecm33-independent-eval.md"
FREEZE = ROOT / "analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json"
FREEZE_SHA = "22e6727879f9847b166b3ccdc75ceb702a87589b504983c5208589ee1dffc270"


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_no_committed_file_holds_a_sequence():
    seq = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{40,}")
    files = [f for f in sorted(REC.rglob("*")) if f.is_file()]
    assert len(files) == 25
    for f in files:
        assert f.suffix in {".tsv", ".txt", ".md", ".json", ".sha256"}, f
        for line in f.read_text().splitlines():
            assert not seq.search(line), f
    assert not list(REC.rglob("*.faa")) and not list(REC.rglob("*.hmm"))


def test_run_json_carries_the_signed_freeze_hash_and_hit_file_hashes():
    assert _sha(FREEZE) == FREEZE_SHA
    run = json.loads((REC / "search/run.json").read_text())
    assert run["freeze_sha256"] == FREEZE_SHA
    assert run["commit"] == "a8d0876afc2bf9a56f52c87d616f7cc8feb848f9"
    assert len(run["files"]) == 6
    for name, digest in run["files"].items():
        assert _sha(REC / "search" / name) == digest, name
    ext = json.loads((REC / "extract/run.json").read_text())
    assert ext["freeze_sha256"] == FREEZE_SHA
    assert ext["reps_fasta_sha256"].startswith("493d386e31a87e45")
    assert FREEZE_SHA in (REC / "report/summary.txt").read_text()


def test_counts_187_109_and_66_separated():
    counts = json.loads((REC / "extract/counts.json").read_text())
    assert counts["n_clusters_all"] == 187
    assert counts["n_clusters_pe"] == 109
    rows = [line.split("\t") for line in (REC / "report/clusters.tsv").read_text().splitlines()[1:]]
    header = (REC / "report/clusters.tsv").read_text().splitlines()[0].split("\t")
    sep = header.index("separated")
    pe = header.index("headline_pe")
    assert len(rows) == 187
    assert sum(r[sep] == "yes" for r in rows) == 66
    assert sum(r[pe] == "yes" for r in rows) == 109
    assert sum(r[pe] == "yes" and r[sep] == "yes" for r in rows) == 15
    assert (
        "test clusters: 187 (at the PANTHER E cut: 109); separated: 66"
        in (REC / "extract/provenance.txt").read_text()
    )


def test_summary_verdict_headline_and_recall_rows():
    text = (REC / "report/summary.txt").read_text()
    assert "not testable (u < 3): u = 1" in text
    head = next(x for x in text.splitlines() if x.startswith("[headline] all=all"))
    assert "n 15 clusters" in head
    assert "recall custom 14/15 (CI 0.681-0.998)" in head
    assert "jackhmmer_seed_only 14/15" in head
    assert "phmmer_seed_only 5/15 (CI 0.118-0.616)" in head
    assert "pfam_ga 1/15" in head and "pfam_relaxed 1/15" in head
    d2 = next(x for x in text.splitlines() if x.startswith("[descriptive_2] all=all"))
    assert "custom 185/187" in d2 and "jackhmmer_seed_only 185/187" in d2
    assert "phmmer_seed_only 124/187" in d2 and "pfam_ga 28/187" in d2
    assert "estimated" not in text.replace("not estimated", "")


def test_report_states_the_verdict_and_the_reading():
    text = REPORT.read_text()
    assert "not testable (u < 3): u = 1" in text
    assert "UniRef50_A0A9P5VCA1" in text
    assert "14/15" in text and "5/15" in text and "108/109" in text
    assert "51 of 187" in text and "no recall advantage" in text
    assert FREEZE_SHA in text
    u = (REC / "report/u_clusters.tsv").read_text()
    assert "UniRef50_A0A9P5VCA1" in u and u.count("\n") == 3
