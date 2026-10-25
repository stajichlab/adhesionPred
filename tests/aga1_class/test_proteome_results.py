"""Tests of merge_proteome_chunks.py and proteome_summary.py (SYNTHETIC data) and of the committed proteome-test record."""

import csv
import hashlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "analysis/aga1_class"


def load(name):
    spec = importlib.util.spec_from_file_location(name, D / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mc = load("merge_proteome_chunks")
ps = load("proteome_summary")
REC = ROOT / "docs/reports/data/aga1_proteome_test_2026-10-25"

SIG = "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"


def make_results(tmp_path, commit="abc123", skew=None):
    chunks = tmp_path / "chunks.tsv"
    chunks.write_text("chunk\tspecies\tn_proteins\tsha256\nS_001\tS\t1\th1\nS_002\tS\t1\th2\n")
    res = tmp_path / "res"
    for tool, name, body in (
        ("signalp", "prediction_results.txt", lambda c: SIG + f"{c}\tOTHER\t1\t0\t\n"),
        ("predgpi", "predgpi_raw.gff3", lambda c: f"{c}\tPredGPI\tChain\t1\t9\t1.0\t.\t.\te\n"),
        ("aiupred", "disorder_raw.txt", lambda c: f"# banner\n#>{c}\n1\tA\t0.1\n"),
    ):
        for i, c in enumerate(("S_001", "S_002"), 1):
            d = res / tool / c
            d.mkdir(parents=True)
            (d / name).write_text(body(f"P{i}"))
            sha = f"h{i}" if skew != (tool, c) else "bad"
            (d / "provenance.txt").write_text(
                f"commit: {commit} (pinned: {commit})\ntool: {tool}\nchunk: {c} (1 proteins)\nchunk sha256 {sha}\n"
            )
    return chunks, res


def test_merge_concatenates_in_chunk_order_and_keeps_one_header(tmp_path):
    chunks, res = make_results(tmp_path)
    out = tmp_path / "merged"
    lines = mc.merge(chunks, res, out)
    sig = (out / "signalp/prediction_results.txt").read_text().splitlines()
    assert (
        sig[0].startswith("# SignalP")
        and sig[1].startswith("# ID")
        and sig[2:] == ["P1\tOTHER\t1\t0\t", "P2\tOTHER\t1\t0\t"]
    )
    assert (out / "predgpi/predgpi_raw.gff3").read_text().count("PredGPI") == 2
    assert (out / "aiupred/disorder_raw.txt").read_text().count("#>P") == 2
    assert lines[0] == "commit: abc123" and "merge_provenance.txt" in " ".join(
        p.name for p in out.iterdir()
    )
    sha = hashlib.sha256((out / "aiupred/disorder_raw.txt").read_bytes()).hexdigest()
    assert sha in (out / "merge_provenance.txt").read_text()


@pytest.mark.parametrize("what", ["exists", "missing", "hash", "commit"])
def test_merge_failures(tmp_path, what):
    chunks, res = make_results(tmp_path, skew=("predgpi", "S_002") if what == "hash" else None)
    out = tmp_path / "merged"
    if what == "exists":
        out.mkdir()
    if what == "missing":
        (res / "signalp/S_002/prediction_results.txt").unlink()
    if what == "commit":
        (res / "aiupred/S_001/provenance.txt").write_text(
            "commit: zzz (pinned: zzz)\ntool: aiupred\nchunk sha256 h1\n"
        )
    with pytest.raises(mc.MergeError):
        mc.merge(chunks, res, out)


def test_cli_fatal_exit(tmp_path):
    chunks, res = make_results(tmp_path)
    r = subprocess.run(
        [
            sys.executable,
            "-I",
            str(D / "merge_proteome_chunks.py"),
            "--chunks",
            str(chunks),
            "--results",
            str(res),
            "--out",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 1 and r.stderr.startswith("FATAL:")


def synth():
    feat = pd.DataFrame({
        "accession": ["A", "B", "C", "D"], "gene": ["", "", "", ""], "length": [900, 800, 300, 500],
        "sp_yes": [1, 1, 0, 1], "gpi_yes": [1, 1, 0, 0], "gpi_missing": [0, 0, 0, 0],
        "disorder_longest_run": [300, 250, 400, 10], "stalk_run_len": [300, 0, 0, 0],
        "rep_coverage": [0.5, 0.0, 0.0, 0.0], "rep_n_copies": [10, 0, 0, 0], "ctail_hydro": [2.0, 1.0, 0.0, 0.0],
    })  # fmt: skip
    prot = pd.DataFrame(
        {"accession": list("ABCD"), "species": ["S1", "S1", "S2", "S2"], "entry_name": list("abcd")}
    )
    sc = {
        m: pd.DataFrame({"accession": list("ABCD"), f"score_{m}": [0.9, 0.8, 0.95, 0.1]})
        for m in ("logreg", "rf")
    }
    oof = pd.DataFrame(
        {
            "role": ["hard_negative"] * 10 + ["class_positive"],
            "score_logreg": np.linspace(0, 0.9, 11),
            "score_rf": np.linspace(0, 0.9, 11),
        }
    )
    train = pd.DataFrame({"accession": ["A"], "role": ["class_positive"], "pos_tier": ["A"]})
    names = {
        "A": "Flocculation protein FLO1",
        "B": "Uncharacterized protein",
        "C": "Verprolin",
        "D": "Adhesin AWP1",
    }
    return feat, sc, oof, train, prot, names


def test_summary_flags_calls_and_known_roles():
    feat, sc, oof, train, prot, names = synth()
    df, thr, models = ps.build(feat, sc, oof, train, prot, names, 0.2)
    byacc = df.set_index("accession")
    assert byacc.at["A", "broad_flag"] == 1 and byacc.at["A", "narrow_flag"] == 1
    assert byacc.at["B", "broad_flag"] == 1 and byacc.at["B", "narrow_flag"] == 0  # no repeat call
    assert (
        byacc.at["C", "broad_flag"] == 0 and byacc.at["D", "broad_flag"] == 0
    )  # C has no SP/GPI, D no GPI
    assert thr["logreg"] == pytest.approx(
        np.quantile(np.linspace(0, 0.9, 11)[:10], 0.8)
    )  # negatives only
    assert byacc.at["A", "in_850_role"] == "class_positive" and byacc.at["B", "in_850_role"] == ""
    s = ps.summarise(df, models).set_index("species")
    assert (
        s.at["S1", "broad_flag"] == 2
        and s.at["S1", "novel_broad"] == 1
        and s.at["S1", "tierA_broad"] == 1
    )
    rec = ps.name_recovery(df, models).set_index(["species", "scope"])
    assert (
        rec.at[("S2", "adhesin_named"), "n_proteins"] == 1
        and rec.at[("S1", "adhesin_named"), "broad"] == 1
    )
    assert (
        rec.at[("S2", "other"), "2of3_models_and_SP"] == 0
    )  # C has high scores but no signal peptide


def test_summary_rejects_mismatches():
    feat, sc, oof, train, prot, names = synth()
    with pytest.raises(ps.SummaryError, match="different accessions"):
        ps.build(feat, sc, oof, train, prot.iloc[:3], names, 0.05)
    with pytest.raises(ps.SummaryError, match="lacks columns"):
        ps.build(feat.drop(columns=["sp_yes"]), sc, oof, train, prot, names, 0.05)
    with pytest.raises(ps.SummaryError, match="lacks score_rf"):
        ps.build(feat, sc, oof.drop(columns=["score_rf"]), train, prot, names, 0.05)


def test_committed_record_hashes_and_headline_counts():
    text = (REC / "README.md").read_text()
    found = {p: s for s, p in re.findall(r"^([0-9a-f]{64})  (\S+)$", text, flags=re.M)}
    committed = {p for p in found if (REC / p).is_file()}
    assert committed == {p.name for p in REC.iterdir() if p.is_file() and p.name != "README.md"}
    for p in committed:
        assert hashlib.sha256((REC / p).read_bytes()).hexdigest() == found[p], p
    with (REC / "summary_by_species.tsv").open(newline="") as fh:
        rows = {r["species"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert sum(int(r["n_proteins"]) for r in rows.values()) == 27811
    assert (
        int(rows["Cgla_CBS138"]["broad_flag"]),
        int(rows["Calb_SC5314"]["broad_flag"]),
        int(rows["Scer_S288C"]["broad_flag"]),
    ) == (25, 25, 11)
    assert "87ee2da4a08b3218013d6631beb0ef7997118afc" in text
