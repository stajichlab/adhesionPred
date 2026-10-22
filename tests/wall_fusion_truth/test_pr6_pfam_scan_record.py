"""The committed run record of PR 6 step 6.0 (job 30235946): no sequence, the pin, the class counts, the two panel hits."""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/pr6_pfam_scan_2026-10-22"
PIN = "ebb29a9a4ceaee9caf0453f74ed411b2508cce82"
NO_HIT = {"KS0017", "KS0080", "KS0107", "KS0086"}


def _rows(name):
    with open(REC / name, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_no_committed_table_holds_a_sequence():
    seq = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]{40,}$")
    files = [f for f in sorted(REC.glob("*")) if f.suffix in {".tsv", ".txt", ".md"}]
    assert files
    for f in files:
        assert not any(seq.match(line.strip()) for line in f.read_text().splitlines()), f.name


def test_provenance_carries_the_pin_and_the_pfam_release():
    text = (REC / "provenance.txt").read_text()
    assert f"commit: {PIN} (pinned: {PIN})" in text
    assert "pfam_release: 38.2" in text
    assert "pfam_sha256_same_as_earlier_jobs: yes" in text
    assert "queries: 66 sequences" in text


def test_44_clusters_split_into_2_38_and_4():
    rows = _rows("scan_per_cluster.tsv")
    assert len(rows) == 44 and len({r["cluster_id"] for r in rows}) == 44
    counts = {c: sum(r["class"] == c for r in rows) for c in {r["class"] for r in rows}}
    assert counts == {"hit_panel_model": 2, "hit_other_pfam": 38, "no_pfam_hit": 4}
    by_slice = {s: sum(r["slice"] == s for r in rows) for s in ("adhesion", "wall", "fusion")}
    assert by_slice == {"adhesion": 12, "wall": 28, "fusion": 4}
    total = _rows("scan_counts.tsv")
    for r in total:
        if r["group"] == "ALL":
            n = int(r["n_clusters"])
            assert n == int(r["hit_panel_model"]) + int(r["hit_other_pfam"]) + int(r["no_pfam_hit"])


def test_no_pfam_hit_set_and_panel_hits_are_the_named_clusters():
    rows = _rows("scan_per_cluster.tsv")
    assert {r["cluster_id"] for r in rows if r["class"] == "no_pfam_hit"} == NO_HIT
    assert {r["cluster_id"] for r in _rows("custom_targets.tsv")} == NO_HIT
    panel = {r["cluster_id"]: r["panel_models"] for r in _rows("panel_hits.tsv")}
    assert panel == {"KS0082": "PF22799:PIR1-like_C", "KS0116": "PF05730:CFEM"}


def test_proposals_are_14_clusters_in_12_families_and_the_rest_partial():
    rows = _rows("scan_per_cluster.tsv")
    other = [r for r in rows if r["class"] == "hit_other_pfam"]
    prop = [r for r in other if r["proposal"] == "add_to_panel_proposal"]
    part = [r for r in other if r["proposal"] == "partial_domain_for_owner"]
    assert (len(prop), len(part)) == (14, 24)
    fam = _rows("panel_proposals.tsv")
    assert len(fam) == 12 and sum(int(r["n_clusters"]) for r in fam) == 14
    assert all(r["status"].startswith("proposal_for_owner") for r in fam)
