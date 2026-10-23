"""Tests for the custom-model (non-Pfam, -T) extension of wall_fusion_report.py.

A stub hmmsearch reads hits from FASTA headers: PFnnnnn=SCORE for a Pfam hit (needs
--cut_ga), ECM=SCORE for a hit of the custom model (needs -T, forbids --cut_ga).
"""

# ruff: noqa: UP031
import hashlib
import importlib.util
import json
import re
import stat
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis/wall_fusion_report/wall_fusion_report.py"
MODELS_TSV = ROOT / "analysis/wall_fusion_report/wall_fusion_models.tsv"
CUSTOM_TSV = ROOT / "analysis/wall_fusion_report/wall_fusion_custom_models.tsv"
MODEL_DIR = ROOT / "data/custom_models/ecm33_all_2026-10-23"
FAKE_PFAM_FIXTURE = Path(__file__).resolve().parent / "data/fake_pfam_models.tsv"
spec = importlib.util.spec_from_file_location("wall_fusion_report_c", SCRIPT)
wf = importlib.util.module_from_spec(spec)
sys.modules["wall_fusion_report_c"] = wf
spec.loader.exec_module(wf)

STUB = """#!{python}
import re, sys
args = sys.argv[1:]
if args and args[0] == "-h":
    print("# HMMER 3.4 (stub)")
    sys.exit(0)
tbl = args[args.index("--tblout") + 1]
target = args[-1]
custom = "-T" in args
if custom:
    assert "--cut_ga" not in args
    cutoff = float(args[args.index("-T") + 1])
else:
    assert "--cut_ga" in args and "-T" not in args
rows = []
for line in open(target):
    if not line.startswith(">"):
        continue
    words = line[1:].split()
    for w in words[1:]:
        m = re.match(r"(PF\\d+|ECM)=([0-9.]+)$", w)
        if not m:
            continue
        sc = float(m.group(2))
        if custom and m.group(1) == "ECM" and sc >= cutoff:
            rows.append((words[0], "ecm33_all", "-", sc))
        if not custom and m.group(1) != "ECM":
            rows.append((words[0], "NAME", m.group(1) + ".1", sc))
with open(tbl, "w") as out:
    out.write("# stub\\n")
    for prot, name, acc, sc in rows:
        out.write("%s - %s %s 1e-10 %.1f 0.0 1e-10 %.1f 1.0 1 1 1 1 1 1 desc\\n" % (prot, name, acc, sc, sc))
"""


def fake_pfam(path):
    lines = [ln for ln in FAKE_PFAM_FIXTURE.read_text().splitlines() if not ln.startswith("#")]
    blocks = []
    for ln in lines[1:]:
        acc, name = ln.split("\t")
        blocks.append(
            "HMMER3/f\nNAME  %s\nACC   %s.1\nDESC  fake\nGA    25.0 20.0;\nHMM x\n//\n"
            % (name, acc)
        )
    Path(path).write_text("".join(blocks))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_custom_table(path, model, digest=None, cutoff="29.6", **over):
    cols = list(wf.CUSTOM_COLUMNS)
    row = {
        "id": "ecm33_test",
        "name": "ecm33_all",
        "model_path": str(model),
        "cutoff_bits": cutoff,
        "sha256": digest or sha(model),
        "kind": "custom",
        "use": "report_only",
        "validation": "unvalidated_beyond_circular_holdout",
        "supporting_record": "rec",
    }
    row.update(over)
    Path(path).write_text("\t".join(cols) + "\n" + "\t".join(row[c] for c in cols) + "\n")


@pytest.fixture
def env(tmp_path):
    stub = tmp_path / "hmmsearch"
    stub.write_text(STUB.format(python=sys.executable))
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    pfam = tmp_path / "Pfam38.2" / "Pfam-A.hmm"
    pfam.parent.mkdir()
    fake_pfam(pfam)
    model = tmp_path / "model.hmm"
    model.write_text("HMMER3/f\nNAME  ecm33_all\nLENG  1\n//\n")
    return tmp_path, str(stub), str(pfam), model


def run(env, fasta_records, table=None, extra=()):
    tmp, stub, pfam, model = env
    fa = tmp / "prot1.faa"
    fa.write_text("".join(">%s\nMKV\n" % r for r in fasta_records))
    if table is None:
        table = tmp / "custom.tsv"
        write_custom_table(table, model)
    out = tmp / "out"
    argv = [
        "--out-dir", str(out), "--models-tsv", str(MODELS_TSV), "--pfam-hmm", pfam,
        "--hmmsearch", stub, "--tmpdir", str(tmp / "scratch"),
        "--custom-models-tsv", str(table), *extra, str(fa),
    ]  # fmt: skip
    return wf.main(argv), out


def rows(path):
    lines = [ln for ln in Path(path).read_text().splitlines() if not ln.startswith("#")]
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:]]


# ------------------------------------------------------------------ table and command


def test_committed_table_parses_and_matches_committed_model():
    (row,) = wf.load_custom_models(CUSTOM_TSV)
    assert Path(row["model_path"]) == MODEL_DIR / "model.hmm"
    assert row["sha256"] == sha(MODEL_DIR / "model.hmm")
    assert row["sha256"].startswith("9a6e2c812540b0b1")
    cut = json.loads((MODEL_DIR / "cutoff.json").read_text())
    assert row["cutoff"] == cut["cutoff"] == 29.6
    assert (row["kind"], row["use"], row["validation"]) == (
        "custom",
        "report_only",
        "unvalidated_beyond_circular_holdout",
    )
    wf.verify_custom_model(row)
    sums = (MODEL_DIR / "SHA256SUMS").read_text()
    for name in ("model.hmm", "cutoff.json", "provenance.txt"):
        assert "%s  %s" % (sha(MODEL_DIR / name), name) in sums
    assert cut["freeze_sha256"] in (MODEL_DIR / "provenance.txt").read_text()


def test_sbatch_pins_the_same_sha():
    text = (ROOT / "analysis/wall_fusion_report/submit_wall_fusion_report.sbatch").read_text()
    pinned = re.search(r"ECM33_PINNED_SHA256=([0-9a-f]{64})", text).group(1)
    assert pinned == sha(MODEL_DIR / "model.hmm")
    assert "Refusing" in text


@pytest.mark.parametrize(
    "bad",
    [
        {"kind": "pfam"},
        {"use": "call"},
        {"validation": "validated"},
        {"cutoff_bits": "high"},
        {"sha256": "abc"},
        {"id": "PF12454"},
    ],
)
def test_bad_custom_table_rejected(tmp_path, bad):
    model = tmp_path / "m.hmm"
    model.write_text("NAME  ecm33_all\n//\n")
    t = tmp_path / "t.tsv"
    write_custom_table(t, model, **bad)
    with pytest.raises(wf.ReportError):
        wf.load_custom_models(t)


def test_command_construction():
    pf = wf.hmmsearch_command("hs", "m.hmm", "t.faa", "o.tbl", 2)
    assert "--cut_ga" in pf and "-T" not in pf
    cu = wf.hmmsearch_command("hs", "m.hmm", "t.faa", "o.tbl", 2, 29.6)
    assert cu[1:3] == ["-T", "29.6"] and "--cut_ga" not in cu


# ------------------------------------------------------------------ runs


def test_sha_mismatch_refused(env):
    tmp, _stub, _pfam, model = env
    t = tmp / "bad.tsv"
    write_custom_table(t, model, digest="0" * 64)
    with pytest.raises(wf.ReportError, match="sha256"):
        run(env, ["P1 ECM=100.0"], table=t)
    assert not (tmp / "out").exists() or not list((tmp / "out").glob("*.md"))


def test_name_mismatch_refused(env):
    tmp, _stub, _pfam, model = env
    t = tmp / "bad.tsv"
    write_custom_table(t, model, name="other")
    with pytest.raises(wf.ReportError, match="NAME"):
        run(env, ["P1 ECM=100.0"], table=t)


def test_custom_hits_kept_separate_from_pfam(env):
    rc, out = run(env, ["P1 ECM=100.0", "P2 ECM=10.0 PF03198=250.0", "P3 PF12454=90.0"])
    assert rc == 0
    custom = rows(out / "wall_fusion_custom_hits.tsv")
    assert [r["protein"] for r in custom] == ["P1"]  # P2 is below the 29.6 cutoff
    assert custom[0]["model_sha256"] == sha(env[3]) and custom[0]["cutoff_bits"] == "29.6"
    assert custom[0]["use"] == "report_only"
    counts = rows(out / "wall_fusion_custom_counts.tsv")
    assert counts[0]["n_hit_proteins"] == "1" and counts[0]["custom_model"] == "ecm33_test"
    pfam_hits = rows(out / "wall_fusion_hits.tsv")
    assert {r["protein"] for r in pfam_hits} == {"P2", "P3"}
    assert all(r["model_acc"].startswith("PF") for r in pfam_hits)
    assert not any(
        "ecm33_test" in ln for ln in (out / "wall_fusion_hits.tsv").read_text().splitlines()
    )
    for r in rows(out / "wall_fusion_counts.tsv"):
        assert not r["model_acc"].startswith("ecm33")
    md = (out / "wall_fusion_report.md").read_text()
    assert "## Custom models" in md and sha(env[3]) in md and "29.6" in md
    assert "unvalidated_beyond_circular_holdout" not in md.split("## Custom models")[0]
    for t in ("custom_hits", "custom_counts"):
        assert (out / ("wall_fusion_%s.tsv.NOTICE" % t)).exists()


def test_without_option_outputs_unchanged(env):
    tmp, stub, pfam, _model = env
    fa = tmp / "prot1.faa"
    fa.write_text(">P1 PF03198=250.0\nMKV\n")
    out = tmp / "out"
    assert (
        wf.main(
            [
                "--out-dir",
                str(out),
                "--models-tsv",
                str(MODELS_TSV),
                "--pfam-hmm",
                pfam,
                "--hmmsearch",
                stub,
                "--tmpdir",
                str(tmp / "s"),
                str(fa),
            ]
        )  # fmt: skip
        == 0
    )
    names = sorted(p.name for p in out.iterdir())
    assert names == [
        "wall_fusion_counts.tsv",
        "wall_fusion_counts.tsv.NOTICE",
        "wall_fusion_hits.tsv",
        "wall_fusion_hits.tsv.NOTICE",
        "wall_fusion_report.md",
    ]
    assert "Custom models" not in (out / "wall_fusion_report.md").read_text()
