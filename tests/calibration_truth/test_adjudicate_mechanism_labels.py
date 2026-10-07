import csv
import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "analysis/calibration_truth/adjudicate_mechanism_labels.py"
)
spec = importlib.util.spec_from_file_location("adjudicate_mechanism_labels", SCRIPT)
adj = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adj)

COLS = [
    "accession",
    "gene",
    "organism",
    "cls",
    "evidence_level",
    "cluster_id",
    "mechanism_label",
    "label_confidence",
    "mechanism_source_pmid",
    "mechanism_quote",
    "evidence_basis",
]


def decide(votes):
    return adj.decide(votes)


def test_highest_tier_wins_over_lower_tiers():
    votes = [("A", "2b-i", 1), ("B", "2a", 3), ("UniProt", "2a", 2)]
    assert decide(votes) == ("2a", "stated_in_paper", "higher_tier_wins", "no")


def test_agreement_is_recorded_as_agreement():
    assert decide([("A", "2a", 1), ("B", "2a", 3)])[2] == "agreement"


def test_two_labels_at_the_same_top_tier_go_to_an_expert():
    label, tier, rule, expert = decide([("A", "2a", 1), ("B", "other", 1)])
    assert (label, rule, expert) == ("unresolved", "conflict_at_same_tier", "yes")


def test_a_lower_tier_conflict_does_not_need_an_expert():
    assert decide([("A", "2a", 1), ("B", "other", 3)])[3] == "no"


def test_no_votes_is_no_claim():
    assert decide([]) == ("unknown", "none", "no_claim", "no")


def row(acc, label, basis, cluster="c1"):
    return {
        "accession": acc,
        "gene": acc,
        "organism": "x",
        "cls": "adhesin",
        "evidence_level": "E1",
        "cluster_id": cluster,
        "mechanism_label": label,
        "label_confidence": "high",
        "mechanism_source_pmid": "1",
        "mechanism_quote": "q",
        "evidence_basis": basis,
    }


def write(path, rows, cols=COLS):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def run(tmp_path, a, b, uni, monkeypatch):
    pa, pb, pu = tmp_path / "a.tsv", tmp_path / "b.tsv", tmp_path / "u.tsv"
    write(pa, a)
    write(pb, b)
    write(pu, uni, ["accession", "n_repeat_features"])
    out = tmp_path / "out"
    monkeypatch.setattr(
        sys,
        "argv",
        ["x", "--agent-a", str(pa), "--agent-b", str(pb), "--uniprot", str(pu), "--out", str(out)],
    )
    adj.main()
    return out


def test_an_other_vote_is_not_evidence_against_repeats(tmp_path, monkeypatch):
    # P1: agent A says 2a by family; agent B says `other` from a paper; UniProt records repeats.
    a = [row("P1", "2a", "family_inference")]
    b = [row("P1", "other", "stated_in_paper")]
    out = run(tmp_path, a, b, [{"accession": "P1", "n_repeat_features": "5"}], monkeypatch)
    rec = list(csv.DictReader(open(out / "curation_table.consensus.tsv"), delimiter="\t"))[0]
    assert rec["final_label"] == "other"  # the mechanism question follows the paper
    assert rec["repeat_evidence_tier"] == "database_annotation"  # the repeat question does not
    assert "2a clusters" not in (out / "adjudication_summary.md").read_text()


def test_cluster_counts_are_cumulative_by_tier(tmp_path, monkeypatch):
    a = [
        row("P1", "2a", "stated_in_paper", "c1"),
        row("P2", "2a", "family_inference", "c2"),
        row("P3", "unknown", "background_knowledge", "c3"),
    ]
    b = [
        row("P1", "2a", "stated_in_paper", "c1"),
        row("P2", "unknown", "family_inference", "c2"),
        row("P3", "unknown", "background_knowledge", "c3"),
    ]
    uni = [{"accession": "P3", "n_repeat_features": "2"}]
    out = run(tmp_path, a, b, uni, monkeypatch)
    text = (out / "adjudication_summary.md").read_text()
    assert "| paper_stated | 1 | 20 |" in text
    assert "| + database_annotation (UniProt repeat features) | 2 | 20 |" in text
    assert "| + family_inference | 3 | 20 |" in text


def test_one_repeat_feature_is_not_enough(tmp_path, monkeypatch):
    a = [row("P1", "unknown", "background_knowledge")]
    b = [row("P1", "unknown", "background_knowledge")]
    out = run(tmp_path, a, b, [{"accession": "P1", "n_repeat_features": "1"}], monkeypatch)
    rec = list(csv.DictReader(open(out / "curation_table.consensus.tsv"), delimiter="\t"))[0]
    assert rec["repeat_evidence_tier"] == "none"


def test_a_bad_basis_is_refused(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="evidence_basis"):
        run(
            tmp_path,
            [row("P1", "2a", "guess")],
            [row("P1", "2a", "stated_in_paper")],
            [],
            monkeypatch,
        )


def test_different_row_sets_are_refused(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="different rows"):
        run(
            tmp_path,
            [row("P1", "2a", "stated_in_paper")],
            [row("P2", "2a", "stated_in_paper")],
            [],
            monkeypatch,
        )
