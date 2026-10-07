import csv
import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2] / "analysis/calibration_truth/compare_mechanism_labels.py"
)
spec = importlib.util.spec_from_file_location("compare_mechanism_labels", SCRIPT)
cml = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cml)

COLS = [
    "accession",
    "gene",
    "cls",
    "evidence_level",
    "cluster_id",
    "mechanism_label",
    "label_confidence",
    "mechanism_source_pmid",
    "mechanism_quote",
    "evidence_basis",
]


def write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLS})


def row(acc, cluster, label, conf="high", quote="tandem repeats", basis="stated_in_paper"):
    return {
        "accession": acc,
        "gene": acc,
        "cls": "adhesin",
        "evidence_level": "E1",
        "cluster_id": cluster,
        "mechanism_label": label,
        "label_confidence": conf,
        "mechanism_source_pmid": "1",
        "mechanism_quote": quote,
        "evidence_basis": basis,
    }


def run(tmp_path, a, b, monkeypatch):
    pa, pb = tmp_path / "a.tsv", tmp_path / "b.tsv"
    write(pa, a)
    write(pb, b)
    out = tmp_path / "out"
    monkeypatch.setattr(
        sys, "argv", ["x", "--agent-a", str(pa), "--agent-b", str(pb), "--out", str(out)]
    )
    cml.main()
    return out


def test_kappa_known_values():
    assert cml.cohen_kappa([("a", "a"), ("b", "b")]) == 1.0
    # observed 0.5, expected 0.5 -> kappa 0
    assert cml.cohen_kappa([("a", "a"), ("a", "b"), ("b", "a"), ("b", "b")]) == pytest.approx(0.0)
    assert cml.cohen_kappa([("a", "a"), ("a", "a")]) is None  # undefined when expected is 1


def test_quote_check_needs_words_that_fit_the_label():
    assert cml.quote_check("2a", "20 tandem repeats") == "pass"
    assert cml.quote_check("2a", "GPI-anchored cell wall protein") == "fail"
    assert cml.quote_check("2c", "a class I hydrophobin") == "pass"
    assert cml.quote_check("other", "anything") == "na"


def test_summary_counts_agreement_and_strict_clusters(tmp_path, monkeypatch):
    a = [
        row("P1", "c1", "2a"),
        row("P2", "c2", "2a", quote="cell wall protein"),  # 2a, quote fails
        row("P3", "c3", "2a"),
        row("P4", "c4", "other", conf="low", quote="adhesion"),
    ]
    b = [
        row("P1", "c1", "2a"),
        row("P2", "c2", "2a", quote="cell wall protein"),
        row("P3", "c3", "other", quote="adhesion"),  # disagrees on 2a
        row("P4", "c4", "other", conf="low", quote="adhesion"),
    ]
    out = run(tmp_path, a, b, monkeypatch)
    text = (out / "summary.md").read_text()
    assert "3 of 4 agree" in text  # P1, P2, P4
    # strict rule: only c1 (both 2a and both quotes pass)
    assert (
        "| Both agree on 2a and both quotes pass the keyword check (strictest) | 1 | 20 |" in text
    )
    assert "| Both agents label at least one row of the cluster 2a | 2 | 20 |" in text
    assert "| Either agent labels at least one row of the cluster 2a (loosest) | 3 | 20 |" in text
    queue = list(csv.DictReader(open(out / "review_queue.tsv"), delimiter="\t"))
    reasons = {r["accession"]: r["queue_reasons"] for r in queue}
    assert "2a_vs_not_2a" in reasons["P3"]
    assert "quote_check_fail_A" in reasons["P2"]
    assert "low_confidence" in reasons["P4"]
    assert "P1" not in reasons  # agreed, high, quote passes: no review needed


def test_extra_columns_are_carried_and_counted(tmp_path, monkeypatch):
    a = [row("P1", "c1", "2a", basis="stated_in_paper")]
    b = [row("P1", "c1", "2a", basis="family_inference")]
    out = run(tmp_path, a, b, monkeypatch)
    comp = list(csv.DictReader(open(out / "comparison.tsv"), delimiter="\t"))[0]
    assert comp["evidence_basis_A"] == "stated_in_paper"
    assert comp["evidence_basis_B"] == "family_inference"
    assert "`family_inference`: 0 / 1" in (out / "summary.md").read_text()


def test_different_row_sets_are_refused(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="different rows"):
        run(tmp_path, [row("P1", "c1", "2a")], [row("P2", "c1", "2a")], monkeypatch)


def test_a_label_outside_the_class_list_is_refused(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="not one of"):
        run(tmp_path, [row("P1", "c1", "repeat")], [row("P1", "c1", "2a")], monkeypatch)


def test_an_empty_label_is_refused(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="not one of"):
        run(tmp_path, [row("P1", "c1", "")], [row("P1", "c1", "2a")], monkeypatch)


def test_a_missing_required_column_is_refused(tmp_path):
    p = tmp_path / "a.tsv"
    p.write_text("accession\tcls\nP1\tadhesin\n")
    with pytest.raises(ValueError, match="missing columns"):
        cml.read_table(p)
