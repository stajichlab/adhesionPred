"""Model card, curated-ortholog count and margin report of the custom Ecm33 `all` model (2026-10-24)."""

import csv
import hashlib
import json
import re
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "data/custom_models/ecm33_all_2026-10-23"
CARD = MODEL / "MODEL_CARD.md"
COUNT = ROOT / "docs/reports/2026-10-24-ecm33-curated-ortholog-count.md"
MARGIN = ROOT / "docs/reports/2026-10-24-ecm33-margin-report.md"
ORTH = ROOT / "docs/reports/data/ecm33_curated_orthologs"
BUILD = ROOT / "docs/reports/data/pr6_ecm33_build_2026-10-23"
EVAL = ROOT / "docs/reports/data/pr6_ecm33_independent_eval_2026-10-24"
MODEL_SHA = "9a6e2c812540b0b1892f081611bf4044ec674bc524d68f5c6fc3d16001c50937"


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_documents_exist_with_their_headings():
    card = CARD.read_text()
    for h in (
        "## 1. Name and version",
        "## 2. What it is for, and what it is not",
        "## 3. How it was built",
        "## 4. Cutoff rule and file hashes",
        "## 5. How to run it",
        "## 6. What the evaluations showed",
        "## 7. Known limits",
        '## 8. What "additional development" could mean',
        "## 9. Owner decision",
    ):
        assert h in card, h
    count = COUNT.read_text()
    for h in (
        "## 1. Answer in short",
        "## 3. Swiss-Prot reviewed entries",
        "## 4. Literature-named proteins",
        "## 5. Conclusion in numbers",
    ):
        assert h in count, h
    margin = MARGIN.read_text()
    for h in (
        "## 2. Summary",
        "## 3. Held-out members per unit",
        "## 5. False hits",
        "## 6. Shared hard negatives",
        "## 7. What this does and does not say",
    ):
        assert h in margin, h
    assert "does NOT change the cutoff" in margin and "new signed freeze" in margin
    assert "MODEL_CARD.md" in (MODEL / "README.md").read_text()


def test_card_hash_lines_equal_sha256sums_and_model_is_unchanged():
    sums = (MODEL / "SHA256SUMS").read_text().strip().splitlines()
    assert len(sums) == 3
    card = CARD.read_text()
    for line in sums:
        assert line in card, line
        digest, name = line.split()
        assert _sha(MODEL / name) == digest
    assert _sha(MODEL / "model.hmm") == MODEL_SHA
    freeze = ROOT / "analysis/wall_fusion_truth/pr6_ecm33_freeze.json"
    cut = json.loads((MODEL / "cutoff.json").read_text())
    assert cut["cutoff"] == 29.6
    assert _sha(freeze) == cut["freeze_sha256"] and cut["freeze_sha256"] in card


def test_card_command_matches_the_report_script():
    src = (ROOT / "analysis/wall_fusion_report/wall_fusion_report.py").read_text()
    assert '["-T", "%g" % cutoff_bits]' in src and '"--noali"' in src and '"--tblout"' in src
    assert "hmmsearch -T 29.6 --noali --cpu N --tblout hits.tbl" in CARD.read_text()


def test_card_numbers_recompute():
    members = _tsv(BUILD / "members/members.tsv")
    assert len(members) == 481 and len({m["clade_fine"] for m in members}) == 19
    side = _tsv(BUILD / "report/side_by_side.tsv")
    assert sum(int(r["n_held_out"]) for r in side) == 469
    assert sum(int(r["custom_hit_headline"]) for r in side) == 469
    assert all(r["u_headline"] == "0" and r["u_second"] == "0" for r in side)
    summary = (EVAL / "report/summary.txt").read_text()
    assert "not testable (u < 3): u = 1" in summary
    assert "recall custom 14/15" in summary and "jackhmmer_seed_only 14/15" in summary
    assert "n 187 clusters" in summary and "recall custom 185/187" in summary
    card = CARD.read_text()
    for s in ("481 members in 19 clades", "14/15", "28 of 187", "469 of 469", "u = 1"):
        assert s in card, s


def test_ortholog_table_and_count():
    cand = [
        r
        for r in _tsv(ROOT / "docs/reports/data/ecm33_independent_set_candidates/candidates.tsv")
        if r["tier"] == "T1_swissprot"
    ]
    assert len(cand) == 13
    rows = _tsv(ORTH / "curated_orthologs.tsv")
    swiss = [r for r in rows if r["set"] == "swissprot_PTHR31018"]
    assert sorted(r["accession"] for r in swiss) == sorted(r["accession"] for r in cand)
    assert len({r["uniref50"] for r in swiss}) == 9
    assert (
        len({r["organism"].split(" (")[0] for r in swiss}) == 5
        or len({r["organism"] for r in swiss}) >= 5
    )
    assert len({r["clade_fine"] for r in swiss}) == 3
    members = {m["hit_id"] for m in _tsv(BUILD / "members/members.tsv")}
    assert all(r["uniref50"] in members for r in swiss)
    assert all(r["uniref50_is_training_member"] == "yes" for r in swiss)
    lit = [r for r in rows if r["set"] == "literature_named"]
    assert len(lit) == 8 and len({r["accession"] for r in lit}) == 8
    assert all(r["independent"] == "no" for r in rows)
    sep = {r["cluster_id"]: r["separated"] for r in _tsv(EVAL / "extract/separation.tsv")}
    n1 = [r for r in rows if r["accession"] == "N1RSC8"][0]
    assert sep[n1["uniref50"]] == "no"
    assert sum(1 for v in sep.values() if v == "yes") == 66
    tc = {r["cluster_id"]: r for r in _tsv(EVAL / "extract/test_clusters.tsv")}
    sepd = [tc[c] for c in tc if sep[c] == "yes"]
    assert len(sepd) == 66 and sum(r["tier"] == "T2_named" for r in sepd) == 12
    assert not any(r["tier"] == "T1_swissprot" for r in tc.values())
    text = COUNT.read_text()
    assert "**Zero.**" in text and "13** reviewed entries" in text
    # ids only: no sequence-like strings
    seq = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{40,}")
    for f in ORTH.rglob("*"):
        if f.is_file():
            assert not seq.search(f.read_text()), f


def test_margin_tables_recompute_from_committed_tables(tmp_path):
    sys.path.insert(0, str(ROOT / "analysis/wall_fusion_truth"))
    try:
        import ecm33_margin_tables as m
    finally:
        sys.path.pop(0)
    rows = m.heldout_margins()
    assert len(rows) == 22
    c30 = [r for r in rows if r["mode"] == "cluster30"]
    assert min(r["min_score"] for r in c30) == 86.9
    assert min(r["min_margin"] for r in rows) == 41.4
    # independent recompute of one row from the raw file
    scores = [
        float(r["score"])
        for r in _tsv(BUILD / "ecm33_build/Eurotiomycetes/cluster30/hits_custom.tsv")
    ]
    assert min(scores) == 86.9 and statistics.median(scores) == 230.2
    cut = json.loads((MODEL / "cutoff.json").read_text())
    cutoff = cut["cutoff"]
    ph = _tsv(BUILD / "ecm33_build/all/none/proteome_hits.tsv")
    above = [r for r in ph if float(r["score"]) >= cutoff]
    assert len(above) == 14
    afum = sorted(float(r["score"]) for r in above if r["proteome"] == "Afum_Af293_UniProt")
    assert (
        afum[:4] == [29.6, 33.8, 34.1, 42.7]
        and len(afum) - 2 == cut["extra_calls_by_proteome"]["Afum_Af293_UniProt"]
    )
    assert sum(cut["extra_calls_by_proteome"].values()) == 4
    below = [
        float(r["score"])
        for r in ph
        if r["proteome"] == "Afum_Af293_UniProt" and float(r["score"]) < cutoff
    ]
    assert max(below) == 29.5
    assert sum(cutoff - 10 <= float(r["score"]) < cutoff for r in ph) == 19
    assert all(v == 0.0 for v in cut["hard_negative_call_rates"].values())
    neg = {
        r["accession"]: float(r["score"]) for r in _tsv(BUILD / "ecm33_build/all/none/neg_hits.tsv")
    }
    assert neg["P38248"] == 387.1 and neg["Q4WNS8"] == 444.2 and neg["G2X962"] == 21.1
    cutoff2, sc, n = m.eval_scores()
    assert (len(sc), n, sc[0]) == (185, 187, 41.7)
    assert sum(s < cutoff + 10 for s in sc) == 0 and sum(s < cutoff + 20 for s in sc) == 1
    assert sum(s < cutoff + 45 for s in sc) == 9
    assert round(statistics.median(sc), 1) == 262.3 and sc[-1] == 477.1
    missing = {r["cluster_id"] for r in _tsv(EVAL / "extract/test_clusters.tsv")} - {
        r["id"] for r in _tsv(EVAL / "search/hits_custom.tsv")
    }
    assert missing == {"UniRef50_A0A9P5VCA1", "UniRef50_A0A9W6T7D3"}
    # the committed margin tables equal a fresh run
    out = tmp_path / "o"
    m.main(out)
    for name in ("heldout_margins.tsv", "proteome_margins.tsv", "eval_histogram.txt"):
        assert (out / name).read_text() == (
            ROOT / "docs/reports/data/ecm33_margin_report_2026-10-24" / name
        ).read_text()
    text = MARGIN.read_text()
    assert (
        ROOT / "docs/reports/data/ecm33_margin_report_2026-10-24/eval_histogram.txt"
    ).read_text().strip() in text
