"""Record test for docs/reports/data/aga1_disorder_2026-10-24/ (real AIUPred and IUPred3 results, 850 proteins)."""

import csv
import hashlib
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
REC = ROOT / "docs/reports/data/aga1_disorder_2026-10-24"
TOOLS = ("aiupred", "iupred3")


def _rows(tool):
    with (REC / f"aga1_signalp_predgpi_{tool}.tsv").open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.mark.parametrize("tool", TOOLS)
def test_table_covers_manifest_and_names_the_tool(tool):
    rows = _rows(tool)
    manifest = (ROOT / "analysis/aga1_class/manifest_accessions.txt").read_text().split()
    assert [r["accession"] for r in rows] == manifest and len(rows) == 850
    assert {r["disorder_tool"] for r in rows} == {tool}
    assert sum(int(r["disorder_n_res"]) for r in rows) == 445655
    assert all(0.0 <= float(r["disorder_frac"]) <= 1.0 for r in rows)


def test_signalp_and_gpi_columns_match_the_earlier_record():
    earlier = ROOT / "docs/reports/data/aga1_prediction_step_2026-10-24/aga1_signalp_predgpi.tsv"
    with earlier.open(newline="") as fh:
        base = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for tool in TOOLS:
        for r in _rows(tool):
            assert all(r[k] == v for k, v in base[r["accession"]].items()), (tool, r["accession"])


def test_group_medians_of_the_report():
    def median_run(tool, group):
        runs = sorted(
            int(r["disorder_longest_run"]) for r in _rows(tool) if r["role_group"] == group
        )
        mid = len(runs) // 2
        return runs[mid] if len(runs) % 2 else (runs[mid - 1] + runs[mid]) / 2

    assert median_run("aiupred", "class_positive_A") == 302
    assert median_run("aiupred", "hard_negative") == 59
    assert median_run("aiupred", "background_assumed_negative") == 16
    assert median_run("iupred3", "class_positive_A") == 191
    assert median_run("iupred3", "hard_negative") == 26


def test_candidate_rule_counts_of_the_report():
    ai = _rows("aiupred")
    iu = {r["accession"]: r for r in _rows("iupred3")}

    def count(group, cut, tool_rows):
        return sum(
            r["signal_peptide_pred"] == "yes"
            and r["gpi_pred"] == "yes"
            and int(tool_rows[r["accession"]]["disorder_longest_run"]) >= cut
            for r in ai
            if r["role_group"] == group
        )

    a = {r["accession"]: r for r in ai}
    assert count("class_positive_A", 200, a) == 23 and count("hard_negative", 200, a) == 8
    assert count("background_assumed_negative", 200, a) == 7
    assert count("class_positive_A", 100, iu) == 21 and count("hard_negative", 100, iu) == 14


def test_readme_hashes_match_committed_files():
    text = (REC / "README.md").read_text()
    found = {p: s for s, p in re.findall(r"^([0-9a-f]{64})  (\S+)$", text, flags=re.M)}
    committed = {p: s for p, s in found.items() if (REC / p).exists()}
    assert set(committed) == {
        "aga1_signalp_predgpi_aiupred.tsv", "aga1_signalp_predgpi_iupred3.tsv",
        "raw_provenance/aiupred_provenance.txt", "raw_provenance/iupred3_provenance.txt",
    }  # fmt: skip
    for rel, sha in committed.items():
        assert hashlib.sha256((REC / rel).read_bytes()).hexdigest() == sha, rel
    assert set(found) - set(committed) == {"aiupred/disorder_raw.txt", "iupred3/disorder_raw.txt"}
    assert "708384fbdbe4a2b24714c8547a12ef23671faff7e1bd75fd5c9cf8475ac63a0e" in text
    assert "424f5e8861dfd4c0cb56fd7533c5e69610120086" in text
