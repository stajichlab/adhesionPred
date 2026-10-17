"""Tests for analysis/wall_fusion_report/wall_fusion_report.py.

The tests use a stub hmmsearch and a tiny fake Pfam file. They do not need HMMER or Pfam.
The stub reads hits from the FASTA headers: a header word ACC=SCORE (for example
PF03198=250.0) gives a hit of that model with that sequence score. A domain score equal
to the sequence score is written, unless the word DOM=VALUE is also in the header.
"""

# ruff: noqa: UP031
import importlib.util
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FAKE_PFAM_FIXTURE = Path(__file__).resolve().parent / "data/fake_pfam_models.tsv"
SCRIPT = ROOT / "analysis/wall_fusion_report/wall_fusion_report.py"
MODELS_TSV = ROOT / "analysis/wall_fusion_report/wall_fusion_models.tsv"
spec = importlib.util.spec_from_file_location("wall_fusion_report", SCRIPT)
wf = importlib.util.module_from_spec(spec)
sys.modules["wall_fusion_report"] = wf
spec.loader.exec_module(wf)

VERIFIED_TSV = Path(__file__).resolve().parent / "data/verified_pfam_ids.tsv"
DOC05 = ROOT / "docs/reports/cell-wall-fusion/05-pfam-id-verification.md"
SCAN_HITS = ROOT / "analysis/cell_wall_check/full_pfam_ga_hits.tsv"


def read_verified():
    """Committed list extracted from reports 05 and 07: {id: row}."""
    lines = [ln for ln in VERIFIED_TSV.read_text().splitlines() if not ln.startswith("#")]
    head = lines[0].split("\t")
    rows = [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:]]
    assert len({r["id"] for r in rows}) == len(rows)
    return {r["id"]: r for r in rows}


VERIFIED = read_verified()
# IDs 05 verified or 07 measured that stay out on purpose (scope, not verification).
EXCLUDED = {
    "PF00026",  # generic aspartyl protease
    "PF10333",  # Pga1, an ER GPI enzyme; Pig-X claim UNVERIFIABLE
    "PF00187",  # CBM18, lectin report
    "PF01607",  # CBM14, lectin report
    "PF01476",  # LysM, lectin report
    "PF00018",  # SH3_1: not a fusion family (07 row 26)
    "PF01822",  # WSC: sensor
    "PF06423",
    "PF05024",
    "PF04113",
    "PF04114",
    "PF06432",  # GPI anchor synthesis
}
ENZYME_SETS = {"wall_remodeling", "wall_synthesis"}

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


def fake_pfam(path, rename=None, drop_ga=None):
    """Write the fake Pfam file from its own fixture, not from the model TSV."""
    lines = [ln for ln in FAKE_PFAM_FIXTURE.read_text().splitlines() if not ln.startswith("#")]
    blocks = []
    for ln in lines[1:]:
        acc, name = ln.split("\t")
        name = (rename or {}).get(acc, name)
        ga = "" if acc == drop_ga else "GA    25.0 20.0;\n"
        blocks.append("HMMER3/f\nNAME  %s\nACC   %s.1\nDESC  fake\n%sHMM x\n//\n" % (name, acc, ga))
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


def test_verified_list_matches_report_05_and_scan():
    """The committed list is re-derived from the 05 table and the 07 scan table."""
    parsed = {}
    for line in DOC05.read_text().splitlines():
        if line.startswith("| PF"):
            c = [x.strip() for x in line.strip().strip("|").split("|")]
            parsed[c[0]] = c[2].split(",")[0].strip()
    from_05 = {k: v for k, v in VERIFIED.items() if v["verification_level"] == "05_interpro"}
    assert {k: v["name"] for k, v in from_05.items()} == parsed
    assert len(parsed) == 53
    scan = {}
    for line in SCAN_HITS.read_text().splitlines()[1:]:
        f = line.split("\t")
        scan[f[1]] = f[2]
    only_07 = {k: v for k, v in VERIFIED.items() if v["verification_level"] == "07_hmmscan_only"}
    assert set(only_07) == {"PF27758", "PF12351", "PF29664", "PF05390"}
    for acc, row in only_07.items():
        assert row["supporting_report"] == "07_model_hit_check"
        assert scan[acc] == row["name"]
        assert acc not in parsed
    for row in from_05.values():
        assert row["supporting_report"] == "05_pfam_id_verification"


def test_model_tsv_follows_inclusion_rule():
    """Every model is verified in 05 or measured in 07, with the level and name of the list."""
    models = wf.load_models(str(MODELS_TSV))
    ids = {m["id"] for m in models}
    assert not ids & EXCLUDED
    assert ids <= set(VERIFIED)
    for m in models:
        v = VERIFIED[m["id"]]
        assert m["name"] == v["name"]
        assert m["verification_level"] == v["verification_level"]
        assert m["supporting_report"].strip()
        assert m["verification_status"].strip()
        if v["verification_level"] == "05_interpro":
            assert m["verification_status"].startswith("05 MATCH")
        else:
            assert m["verification_status"].startswith("NOT in 05")
    # symmetric rule: both 07-only wall models and both fusion models are in
    assert {"PF29664", "PF05390", "PF27758", "PF12351"} <= ids


def test_sets_split_enzymes_from_non_enzymes():
    models = wf.load_models(str(MODELS_TSV))
    by_set = {}
    for m in models:
        by_set.setdefault(m["set"], set()).add(m["id"])
    assert set(by_set) == set(wf.SETS)
    assert by_set["wall_remodeling"] == {
        "PF03198", "PF00722", "PF00332", "PF03663", "PF00704", "PF00150", "PF01522", "PF29664",
    }  # fmt: skip
    assert by_set["wall_synthesis"] == {"PF01644", "PF03142", "PF08407", "PF02364"}
    assert by_set["wall_structural_other"] == {
        "PF10342", "PF05390", "PF03856", "PF12454", "PF10287", "PF10290", "PF07983",
    }  # fmt: skip
    assert by_set["cell_fusion_candidate_models"] == {"PF27758", "PF12351"}
    notes = {m["id"]: m["verification_status"] for m in models}
    assert "NO hit" in notes["PF12454"]
    assert "unclear" in notes["PF03856"].lower()
    assert "CE4 mapping is UNVERIFIABLE" in notes["PF01522"]
    assert "UNVERIFIABLE" not in notes["PF01522"].replace("CE4 mapping is UNVERIFIABLE", "")
    assert "475.9" in notes["PF27758"] and "475.5" in notes["PF27758"]


def test_bad_model_tsv_stops(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("id\tname\n")
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))
    hdr = "\t".join(wf.MODEL_COLUMNS) + "\n"
    p.write_text(hdr + "PF1\tx\tx\twall_remodeling\ta\t05_interpro\tb\n")
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))
    p.write_text(hdr + "PF03198\tx\tx\tadhesion\ta\t05_interpro\tb\n")
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))
    p.write_text(hdr + "PF03198\tx\tx\twall_remodeling\ta\tguess\tb\n")
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))
    p.write_text(hdr + "PF03198\tx\tx\twall_remodeling\ta\t05_interpro\tb\n" * 2)
    with pytest.raises(wf.ReportError):
        wf.load_models(str(p))


def test_set_selection():
    wall = wf.load_models(str(MODELS_TSV), "wall_remodeling")
    fus = wf.load_models(str(MODELS_TSV), "cell_fusion_candidate_models")
    both = wf.load_models(str(MODELS_TSV), "wall_synthesis,wall_structural_other")
    assert {m["set"] for m in wall} == {"wall_remodeling"}
    assert {m["id"] for m in fus} == {"PF27758", "PF12351"}
    assert {m["set"] for m in both} == {"wall_synthesis", "wall_structural_other"}
    for bad in ("adhesion", "cell_fusion", "wall_remodeling,x"):
        with pytest.raises(wf.ReportError):
            wf.load_models(str(MODELS_TSV), bad)


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
    assert c[("specB", "PF12351")]["set"] == "cell_fusion_candidate_models"
    hits = read_tsv(out / "wall_fusion_hits.tsv")
    assert len(hits) == 4
    labels = {(h["protein"], h["model_acc"]): h["label"] for h in hits}
    assert labels[("p1", "PF03198")] == "GH72 (PF03198) domain hit"
    assert labels[("q1", "PF12351")] == "Fig1 (PF12351) domain hit"
    sets = {(h["protein"], h["model_acc"]): h["set"] for h in hits}
    assert sets[("p1", "PF03198")] == "wall_remodeling"
    assert sets[("q1", "PF12351")] == "cell_fusion_candidate_models"
    for h in hits:
        assert h["supporting_report"]
        assert h["verification_level"] in ("05_interpro", "07_hmmscan_only")
        assert "fusion" not in h["label"] and "wall" not in h["label"]


def test_set_totals_do_not_mix_enzymes_and_non_enzymes(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF07983=90.0", "p2 PF03198=250.0", "p3 PF01644=80.0", "p4 PF27758=90.0"])
    out = run(env, [a])
    md = (out / "wall_fusion_report.md").read_text()
    table = md.split("## Per proteome and per set")[1].split("## Notes")[0]
    head = [ln for ln in table.splitlines() if ln.startswith("| proteome")][0]
    cols = [x.strip() for x in head.strip("|").split("|")]
    row = [ln for ln in table.splitlines() if ln.startswith("| specA")][0]
    vals = dict(zip(cols, [x.strip() for x in row.strip("|").split("|")], strict=True))
    assert vals["wall_remodeling"] == "1"
    assert vals["wall_synthesis"] == "1"
    assert vals["wall_structural_other"] == "1"
    assert vals["cell_fusion_candidate_models"] == "1"
    assert "total" not in head.lower()


def test_no_call_guarantee(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF03198=250.0 PF27758=90.0"])
    out = run(env, [a])
    assert sorted(p.name for p in out.iterdir()) == [
        "wall_fusion_counts.tsv",
        "wall_fusion_counts.tsv.NOTICE",
        "wall_fusion_hits.tsv",
        "wall_fusion_hits.tsv.NOTICE",
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


def test_notice_survives_comment_read(env):
    """A pandas comment='#' read drops the header lines; the sidecar file keeps the notice."""
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF03198=250.0"])
    out = run(env, [a])
    for name in ("wall_fusion_counts.tsv", "wall_fusion_hits.tsv"):
        rows = read_tsv(out / name)  # same filter as comment="#"
        assert all("NOT evidence" not in v for r in rows for v in r.values())
        notice = (out / (name + ".NOTICE")).read_text()
        assert "makes no call" in notice and "NOT evidence of adhesion" in notice
        assert "comment=" in notice


def test_set_option_limits_models(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF03198=250.0 PF27758=90.0"])
    out = run(env, [a], ["--set", "cell_fusion_candidate_models"])
    hits = read_tsv(out / "wall_fusion_hits.tsv")
    assert [h["model_acc"] for h in hits] == ["PF27758"]
    counts = read_tsv(out / "wall_fusion_counts.tsv")
    assert {r["set"] for r in counts} == {"cell_fusion_candidate_models"}


def test_duplicate_id_inside_one_proteome_warns(env, capsys):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["dup PF03198=100.0", "dup PF03198=300.0", "u1"])
    out = run(env, [a])
    err = capsys.readouterr().err
    assert "specA: 1 records repeat an ID inside the file" in err
    hits = read_tsv(out / "wall_fusion_hits.tsv")
    assert len(hits) == 1
    assert hits[0]["seq_score"] == "300.0"
    counts = {r["model_acc"]: r for r in read_tsv(out / "wall_fusion_counts.tsv")}
    assert counts["PF03198"]["n_proteins"] == "2"  # unique IDs, not 3 records
    assert counts["PF03198"]["n_hit_proteins"] == "1"
    md = (out / "wall_fusion_report.md").read_text()
    assert "repeat an ID inside the file" in md
    assert "counts unique protein IDs" in md


def test_duplicate_ids_across_proteomes_warn(env, capsys):
    tmp = env[0]
    a, b = tmp / "specA.faa", tmp / "specB.faa"
    write_fasta(a, ["same PF03198=250.0", "u1"])
    write_fasta(b, ["same PF03198=150.0", "u2"])
    out = run(env, [a, b])
    assert "WARNING: 1 protein IDs occur in more than one input proteome" in capsys.readouterr().err
    assert "**Warning:**" in (out / "wall_fusion_report.md").read_text()
    assert len(read_tsv(out / "wall_fusion_hits.tsv")) == 2


def test_no_duplicate_warning_when_ids_differ(env, capsys):
    tmp = env[0]
    a, b = tmp / "specA.faa", tmp / "specB.faa"
    write_fasta(a, ["a1"])
    write_fasta(b, ["b1"])
    run(env, [a, b])
    assert "WARNING" not in capsys.readouterr().err


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not installed")
def test_zst_input(env):
    tmp = env[0]
    plain = tmp / "specZ.faa"
    write_fasta(plain, ["z1 PF03198=99.0", "z1 PF03198=120.0", "z2"])
    subprocess.run(["zstd", "-q", str(plain), "-o", str(tmp / "specZ.faa.zst")], check=True)
    plain.unlink()
    out = run(env, [tmp / "specZ.faa.zst"])
    counts = {r["model_acc"]: r for r in read_tsv(out / "wall_fusion_counts.tsv")}
    assert counts["PF03198"]["n_hit_proteins"] == "1"
    assert counts["PF03198"]["n_proteins"] == "2"
    assert counts["PF03198"]["best_score"] == "120.0"


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
    with pytest.raises(wf.ReportError, match="not found"):
        run((tmp, stub, str(empty)), [a])


def test_model_name_mismatch_stops(env):
    tmp, stub, pfam = env
    fake_pfam(pfam, rename={"PF03198": "Glyco_hydro_72_renamed", "PF27758": "Other_name"})
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF03198=250.0"])
    with pytest.raises(wf.ReportError) as exc:
        run(env, [a])
    msg = str(exc.value)
    assert "PF03198" in msg and "Glyco_hydro_72_renamed" in msg and "PF27758" in msg
    assert not (tmp / "out" / "wall_fusion_hits.tsv").exists()


def test_model_without_ga_line_stops_with_exit_code_2(env):
    tmp, stub, pfam = env
    fake_pfam(pfam, drop_ga="PF03198")
    a = tmp / "specA.faa"
    write_fasta(a, ["p1"])
    with pytest.raises(wf.ReportError, match="no GA line"):
        run(env, [a])
    res = subprocess.run(
        [
            sys.executable, str(SCRIPT), "--out-dir", str(tmp / "o2"),
            "--models-tsv", str(MODELS_TSV), "--pfam-hmm", pfam, "--hmmsearch", stub,
            "--tmpdir", str(tmp / "scratch"), str(a),
        ],
        capture_output=True,
        text=True,
    )  # fmt: skip
    assert res.returncode == 2
    assert "no GA line" in res.stderr


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
