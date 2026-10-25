"""Tests of build_hand_check_sheet.py (SYNTHETIC rows) and of the committed hand-check record."""

import csv
import hashlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "analysis/aga1_class"
spec = importlib.util.spec_from_file_location(
    "build_hand_check_sheet", D / "build_hand_check_sheet.py"
)
bh = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bh)
REC = ROOT / "docs/reports/data/aga1_hand_check_2026-10-26"
SAMPLE_REC = ROOT / "docs/reports/data/aga1_candidate_tiers_2026-10-26/hand_check_sample.tsv"


def sample():
    return pd.DataFrame({
        "accession": ["A", "B"], "tier": ["1", "2"], "sample": ["all tier 1", "sampled tier 2"], "species": ["S1", "S2"],
        "gene": ["g1", ""], "entry_name": ["a", "b"], "protein_name": ["Flo11 domain", "Unchar"], "length": ["500", "300"],
        "tier_reason": ["r1", "r2"], "sp_yes": ["1", "1"], "gpi_yes": ["1", "0"], "disorder_longest_run": ["300", "10"],
        "rep_coverage": ["0.4", "0"], "rep_n_copies": ["5", "0"], "broad_flag": ["1", "0"], "narrow_flag": ["1", "0"],
        "n_models_call": ["3", "2"], "score_logreg": ["0.9", "0.6"], "score_rf": ["0.8", "0.5"], "score_gboost": ["0.7", "0.4"],
        "mean_rank_species": ["0.99", "0.95"], "pfam": ["PF1", ""], "interpro": ["", ""], "domains": ["Flo11", ""],
        "keywords": ["Signal", ""], "subcellular_location": ["", ""], "uniprot_release": ["2026_03", "2026_03"],
    })  # fmt: skip


def feats(accs):
    return pd.DataFrame({"accession": accs, **{c: ["0.1"] * len(accs) for c in bh.FEATURE_COLS}})


def lead(acc, pmids="none", conf="low"):
    return {
        "accession": acc, "gene": "g", "family_hint": "", "uniprot_function_or_similarity": "none", "go_cellular_component": "none",
        "pmids": pmids, "literature_lead": "no literature found", "wall_or_enzyme_flags": "none", "lead_confidence": conf, "problems": "none",
    }  # fmt: skip


def test_build_joins_features_and_makes_an_empty_verdict_template():
    sheet = bh.build(sample(), [feats(["A"]), feats(["B"])])
    assert (
        list(sheet.accession) == ["A", "B"]
        and "frac_ST" in sheet.columns
        and "tier" in sheet.columns
    )
    t = bh.verdict_template(sheet)
    assert list(t.columns) == bh.VERDICT_COLS and (t["verdict"] == "").all()


def test_build_failures():
    with pytest.raises(bh.SheetError, match="no feature row"):
        bh.build(sample(), [feats(["A"])])
    with pytest.raises(bh.SheetError, match="more than one"):
        bh.build(sample(), [feats(["A", "B"]), feats(["A"])])


def test_check_verdicts():
    s = sample()
    t = bh.verdict_template(bh.build(s, [feats(["A", "B"])]))
    assert bh.check(t, s) == 2
    t.loc[0, "verdict"] = "supported"
    assert bh.check(t, s) == 1
    bad = t.copy()
    bad.loc[1, "verdict"] = "adhesin"
    with pytest.raises(bh.SheetError, match="outside the vocabulary"):
        bh.check(bad, s)
    with pytest.raises(bh.SheetError, match="missing"):
        bh.check(t.iloc[:1], s)
    extra = pd.concat(
        [t, pd.DataFrame([{**t.iloc[0].to_dict(), "accession": "Z"}])], ignore_index=True
    )
    with pytest.raises(bh.SheetError, match="not in the sample"):
        bh.check(extra, s)
    with pytest.raises(bh.SheetError, match="duplicate"):
        bh.check(pd.concat([t, t.iloc[:1]], ignore_index=True), s)


def test_merge_leads_and_render():
    ev = bh.build(sample(), [feats(["A", "B"])])
    merged = bh.merge_leads(
        ev, [pd.DataFrame([lead("A", "123;456", "high")]), pd.DataFrame([lead("B")])]
    )
    md = bh.render_sheets(merged)
    assert (
        md.count("\n## ") == 2
        and "PMIDs: 123;456." in md
        and "NOT verdicts" in md
        and "Verdict:" in md
    )


@pytest.mark.parametrize(
    "rows, msg",
    [
        ([lead("A")], "no lead row"),
        ([lead("A"), lead("B"), lead("C")], "not in the sample"),
        ([lead("A"), lead("A"), lead("B")], "duplicate"),
        ([lead("A", pmids="PMID123"), lead("B")], "pmids must be numbers"),
        ([lead("A", conf="certain"), lead("B")], "lead_confidence"),
    ],
)
def test_merge_leads_failures(rows, msg):
    ev = bh.build(sample(), [feats(["A", "B"])])
    with pytest.raises(bh.SheetError, match=msg):
        bh.merge_leads(ev, [pd.DataFrame(rows)])


def test_cli_refuses_existing_out_dir(tmp_path):
    s, f = tmp_path / "s.tsv", tmp_path / "f.tsv"
    sample().to_csv(s, sep="\t", index=False)
    feats(["A", "B"]).to_csv(f, sep="\t", index=False)
    r = subprocess.run(
        [
            sys.executable,
            "-I",
            str(D / "build_hand_check_sheet.py"),
            "build",
            "--sample",
            str(s),
            "--features",
            str(f),
            "--out-dir",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 1 and r.stderr.startswith("FATAL:")


def test_committed_record():
    text = (REC / "README.md").read_text()
    found = {p: s for s, p in re.findall(r"^([0-9a-f]{64})  (\S+)$", text, flags=re.M)}
    committed = {
        p.name for p in REC.iterdir() if p.is_file() and p.name not in ("README.md", "verdicts.tsv")
    }
    assert committed <= set(found)
    for p in committed:
        assert hashlib.sha256((REC / p).read_bytes()).hexdigest() == found[p], p
    sample_df = bh.rd(SAMPLE_REC)
    assert len(sample_df) == 97
    v = bh.rd(REC / "verdicts.tsv")
    bh.check(
        v, sample_df
    )  # every sample protein is present, no unknown accession, verdicts inside the vocabulary
    with (REC / "evidence_auto.tsv").open(newline="") as fh:
        assert [r["accession"] for r in csv.DictReader(fh, delimiter="\t")] == list(
            sample_df["accession"]
        )
