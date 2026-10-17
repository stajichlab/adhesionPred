"""Tests for analysis/wall_fusion_report/wall_fusion_report.py.

The tests use a stub hmmsearch and a tiny fake Pfam file. They do not need HMMER or Pfam.
The stub reads hits from the FASTA headers: a header word ACC=SCORE (for example
PF03198=250.0) gives a hit of that model with that sequence score. A domain score equal
to the sequence score is written, unless the word DOM=VALUE is also in the header.
"""

# ruff: noqa: UP031
import importlib.util
import stat
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis/wall_fusion_report/wall_fusion_report.py"
MODELS_TSV = ROOT / "analysis/wall_fusion_report/wall_fusion_models.tsv"
spec = importlib.util.spec_from_file_location("wall_fusion_report", SCRIPT)
wf = importlib.util.module_from_spec(spec)
sys.modules["wall_fusion_report"] = wf
spec.loader.exec_module(wf)

# IDs that docs/reports/cell-wall-fusion/05-pfam-id-verification.md Table A marks MATCH
# and that the reports tie to wall synthesis or remodeling. Committed copy: no network.
VERIFIED_05_WALL = {
    "PF03198": "Glyco_hydro_72",
    "PF00722": "Glyco_hydro_16",
    "PF00332": "Glyco_hydro_17",
    "PF03663": "Glyco_hydro_76",
    "PF00704": "Glyco_hydro_18",
    "PF00150": "Cellulase",
    "PF01644": "Chitin_synth_1",
    "PF03142": "Chitin_synth_2",
    "PF08407": "Chitin_synth_1N",
    "PF02364": "Glucan_synthase",
    "PF10342": "Kre9_KNH",
    "PF03856": "SUN",
    "PF12454": "Ecm33",
    "PF10287": "YJL171C_Tos1_C",
    "PF10290": "YJL171C_Tos1_N",
    "PF07983": "X8",
}
# Fusion IDs are not in 05. Their Pfam names come from the 07 scan of Pfam-A.hmm 38.2.
VERIFIED_07_FUSION = {"PF27758": "PRM1", "PF12351": "Fig1"}
# IDs that must never enter: 05 UNVERIFIABLE or not tied to the sets, and SH3_1 (Fus1).
EXCLUDED = {
    "PF01522",  # CE4 mapping UNVERIFIABLE
    "PF10333",  # Pig-X claim UNVERIFIABLE
    "PF00187",  # CBM18 link UNVERIFIABLE, lectin report
    "PF01607",  # Avr4 link UNVERIFIABLE, lectin report
    "PF00018",  # SH3_1: not a fusion family
}

STUB = """#!{python}
import re, sys
args = sys.argv[1:]
if args and args[0] == "-h":
    print("# hmmsearch :: search profile(s) against a sequence database")
    print("# HMMER 3.4 (stub)")
    sys.exit(0)
tbl = args[args.index("--tblout") + 1]
target = args[-1]
assert "--cut_ga" in args
rows = []
for line in open(target):
    if not line.startswith(">"):
        continue
    words = line[1:].split()
    dom = None
    for w in words[1:]:
        if w.startswith("DOM="):
            dom = float(w[4:])
    for w in words[1:]:
        m = re.match(r"(PF\\d+)=([0-9.]+)$", w)
        if m:
            sc = float(m.group(2))
            rows.append((words[0], m.group(1), sc, sc if dom is None else dom))
with open(tbl, "w") as out:
    out.write("# stub\\n")
    for prot, acc, sc, dsc in rows:
        out.write("%s - NAME %s.1 1e-10 %.1f 0.0 1e-10 %.1f 1.0 1 1 1 1 1 1 desc\\n" % (prot, acc, sc, dsc))
"""


def fake_pfam(path):
    models = wf.load_models(str(MODELS_TSV))
    blocks = []
    for m in models:
        blocks.append(
            "HMMER3/f\nNAME  %s\nACC   %s.1\nDESC  fake\nGA    25.0 20.0;\nHMM x\n//\n"
            % (m["name"], m["id"])
        )
    blocks.append("HMMER3/f\nNAME  Other\nACC   PF99999.1\nGA    1.0 1.0;\n//\n")
    Path(path).write_text("".join(blocks))


@pytest.fixture
def env(tmp_path):
    stub = tmp_path / "hmmsearch"
    stub.write_text(STUB.format(python=sys.executable))
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    pfam = tmp_path / "Pfam38.2" / "Pfam-A.hmm"
    pfam.parent.mkdir()
    fake_pfam(pfam)
    return tmp_path, str(stub), str(pfam)


def write_fasta(path, records):
    Path(path).write_text("".join(">%s\nMKV\nLLA\n" % r for r in records))


def run(env, fastas, extra=()):
    tmp, stub, pfam = env
    out = tmp / "out"
    argv = [
        "--out-dir",
        str(out),
        "--models-tsv",
        str(MODELS_TSV),
        "--pfam-hmm",
        pfam,
        "--hmmsearch",
        stub,
        "--tmpdir",
        str(tmp / "scratch"),
        *extra,
        *[str(f) for f in fastas],
    ]
    assert wf.main(argv) == 0
    return out


def read_tsv(path):
    lines = [ln for ln in Path(path).read_text().splitlines() if not ln.startswith("#")]
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:]]


# ------------------------------------------------------------------ model table


def test_model_tsv_only_verified_ids():
    models = wf.load_models(str(MODELS_TSV))
    allowed = {**VERIFIED_05_WALL, **VERIFIED_07_FUSION}
    ids = {m["id"] for m in models}
    assert ids == set(allowed)
    for m in models:
        assert m["name"] == allowed[m["id"]]
    assert not ids & EXCLUDED
    for m in models:
        assert m["supporting_report"].strip()
        assert m["verification_status"].strip()


def test_model_sets_and_status_text():
    models = wf.load_models(str(MODELS_TSV))
    by_set = {}
    for m in models:
        by_set.setdefault(m["set"], set()).add(m["id"])
    assert by_set["cell_fusion"] == set(VERIFIED_07_FUSION)
    assert by_set["wall_remodeling"] == set(VERIFIED_05_WALL)
    for m in models:
        if m["id"] in VERIFIED_05_WALL:
            assert m["verification_status"].startswith("05 MATCH")
        else:
            assert m["verification_status"].startswith("NOT in 05")


def test_bad_model_tsv_stops(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("id\tname\n")
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))
    hdr = "\t".join(wf.MODEL_COLUMNS) + "\n"
    p.write_text(hdr + "PF1\tx\tx\twall_remodeling\ta\tb\n")
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))
    p.write_text(hdr + "PF03198\tx\tx\tadhesion\ta\tb\n")
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))
    p.write_text(hdr + "PF03198\tx\tx\twall_remodeling\ta\tb\n" * 2)
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))


def test_set_selection():
    wall = wf.load_models(str(MODELS_TSV), "wall_remodeling")
    fus = wf.load_models(str(MODELS_TSV), "cell_fusion")
    assert {m["set"] for m in wall} == {"wall_remodeling"}
    assert {m["id"] for m in fus} == {"PF27758", "PF12351"}
    with pytest.raises(wf.ReportError):
        wf.load_models(str(MODELS_TSV), "adhesion")


# ------------------------------------------------------------------ runs


def test_counts_hits_and_label(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    b = tmp / "specB.faa.gz"
    write_fasta(a, ["p1 PF03198=250.0", "p2 PF03198=40.0 PF27758=30.0", "p3"])
    import gzip

    with gzip.open(b, "wt") as fh:
        fh.write(">q1 PF12351=50.0\nMKV\n>q2\nMKV\n")
    out = run(env, [a, b])
    counts = read_tsv(out / "wall_fusion_counts.tsv")
    c = {(r["proteome"], r["model_acc"]): r for r in counts}
    n_models = len(wf.load_models(str(MODELS_TSV)))
    assert len(counts) == 2 * n_models
    assert c[("specA", "PF03198")]["n_hit_proteins"] == "2"
    assert c[("specA", "PF03198")]["best_score"] == "250.0"
    assert c[("specA", "PF03198")]["n_proteins"] == "3"
    assert c[("specA", "PF03198")]["family"] == "GH72"
    assert c[("specB", "PF03198")]["best_score"] == "NA"
    assert c[("specB", "PF12351")]["n_hit_proteins"] == "1"
    assert c[("specB", "PF12351")]["set"] == "cell_fusion"
    hits = read_tsv(out / "wall_fusion_hits.tsv")
    assert len(hits) == 4
    labels = {(h["protein"], h["model_acc"]): h["label"] for h in hits}
    assert labels[("p1", "PF03198")] == "wall_remodeling_family_hit: GH72"
    assert labels[("q1", "PF12351")] == "cell_fusion_family_hit: Fig1"
    for h in hits:
        assert h["supporting_report"]


def test_no_call_guarantee(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF03198=250.0 PF27758=90.0"])
    out = run(env, [a])
    assert sorted(p.name for p in out.iterdir()) == [
        "wall_fusion_counts.tsv",
        "wall_fusion_hits.tsv",
        "wall_fusion_report.md",
    ]
    for name in ("wall_fusion_counts.tsv", "wall_fusion_hits.tsv"):
        text = (out / name).read_text()
        header = [ln for ln in text.splitlines() if not ln.startswith("#")][0].split("\t")
        for col in header:
            assert "call" not in col.lower()
            assert "status" not in col.lower()
        assert "NOT evidence of adhesion" in text.splitlines()[0]
        assert "false-positive" in text.splitlines()[1]
    md = (out / "wall_fusion_report.md").read_text()
    assert "makes no call" in md
    assert "NOT evidence of adhesion" in md
    assert "6 of 13" in md
    assert "SH3_1" in md


def test_set_option_limits_models(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF03198=250.0 PF27758=90.0"])
    out = run(env, [a], ["--set", "cell_fusion"])
    hits = read_tsv(out / "wall_fusion_hits.tsv")
    assert [h["model_acc"] for h in hits] == ["PF27758"]
    counts = read_tsv(out / "wall_fusion_counts.tsv")
    assert {r["set"] for r in counts} == {"cell_fusion"}


def test_best_score_kept_for_repeated_id(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF03198=100.0", "p1 PF03198=300.0"])
    out = run(env, [a])
    hits = read_tsv(out / "wall_fusion_hits.tsv")
    assert len(hits) == 1
    assert hits[0]["seq_score"] == "300.0"


def test_domain_below_ga_is_flagged(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF03198=250.0 DOM=5.0", "p2 PF03198=250.0"])
    out = run(env, [a])
    h = {r["protein"]: r for r in read_tsv(out / "wall_fusion_hits.tsv")}
    assert h["p1"]["domain_ge_ga"] == "no"
    assert h["p2"]["domain_ge_ga"] == "yes"


def test_missing_model_stops(env):
    tmp, stub, _pfam = env
    empty = tmp / "empty.hmm"
    empty.write_text("HMMER3/f\nNAME  Other\nACC   PF99999.1\nGA    1.0 1.0;\n//\n")
    a = tmp / "specA.faa"
    write_fasta(a, ["p1"])
    with pytest.raises(wf.ReportError):
        run((tmp, stub, str(empty)), [a])


def test_needs_tmpdir_or_scratch(env, monkeypatch):
    tmp, stub, pfam = env
    monkeypatch.delenv("SCRATCH", raising=False)
    a = tmp / "specA.faa"
    write_fasta(a, ["p1"])
    argv = [
        "--out-dir", str(tmp / "o"), "--models-tsv", str(MODELS_TSV),
        "--pfam-hmm", pfam, "--hmmsearch", stub, str(a),
    ]  # fmt: skip
    with pytest.raises(wf.ReportError):
        wf.main(argv)


def test_duplicate_names_stop(env):
    tmp = env[0]
    (tmp / "d1").mkdir()
    (tmp / "d2").mkdir()
    a, b = tmp / "d1" / "sp.faa", tmp / "d2" / "sp.fasta"
    write_fasta(a, ["p1"])
    write_fasta(b, ["p2"])
    with pytest.raises(wf.ReportError):
        run(env, [a, b])


def test_empty_fasta_does_not_abort(env, capsys):
    tmp = env[0]
    a, b = tmp / "specA.faa", tmp / "empty.faa"
    write_fasta(a, ["p1 PF03198=90.0"])
    b.write_text("")
    out = run(env, [a, b])
    assert "no sequences" in capsys.readouterr().err
    counts = read_tsv(out / "wall_fusion_counts.tsv")
    assert {r["n_proteins"] for r in counts if r["proteome"] == "empty"} == {"0"}


def test_proteome_name():
    assert wf.proteome_name("/x/Sc.S288C.faa.gz") == "Sc.S288C"
    assert wf.proteome_name("a.fasta.zst") == "a"


def test_script_names_no_config_files():
    """The script lives in analysis/ only. It must not name status or config files to write."""
    text = SCRIPT.read_text()
    assert "family_table.tsv" not in text.split('"""', 2)[2]
    assert "categories.yaml" not in text.split('"""', 2)[2]
