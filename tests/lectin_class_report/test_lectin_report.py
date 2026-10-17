"""Tests for analysis/lectin_class_report/lectin_report.py.

The tests use a stub hmmsearch and a tiny fake Pfam file. They do not need HMMER or Pfam.
The stub reads hits from the FASTA headers: a header word ACC=SCORE (for example
PF07938=250.0) gives a hit of that model with that sequence score. A domain score equal
to the sequence score is written, unless the word DOM=VALUE is also in the header.
"""

# ruff: noqa: UP031
import gzip
import importlib.util
import shutil
import stat
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/lectin_class_report/lectin_report.py"
spec = importlib.util.spec_from_file_location("lectin_report", SCRIPT)
lr = importlib.util.module_from_spec(spec)
sys.modules["lectin_report"] = lr
spec.loader.exec_module(lr)

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
    blocks = []
    for acc, *_ in lr.MODELS:
        blocks.append(
            "HMMER3/f\nNAME  Model_%s\nACC   %s.1\nDESC  fake\nGA    25.0 20.0;\nHMM x\n//\n"
            % (acc, acc)
        )
    # one extra model that must not be extracted
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
    text = "".join(">%s\nMKV\nLLA\n" % r for r in records)
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)


def run(env, fastas, extra=()):
    tmp, stub, pfam = env
    out = tmp / "out"
    argv = [
        "--out-dir",
        str(out),
        "--pfam-hmm",
        pfam,
        "--hmmsearch",
        stub,
        "--tmpdir",
        str(tmp / "scratch"),
        *extra,
        *[str(f) for f in fastas],
    ]
    assert lr.main(argv) == 0
    return out


def read_tsv(path):
    lines = Path(path).read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"), strict=False)) for ln in lines[1:]]


def test_counts_hits_and_label(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    b = tmp / "specB.faa.gz"
    write_fasta(a, ["p1 PF07938=250.0", "p2 PF07938=40.0 PF01476=30.0", "p3"])
    write_fasta(b, ["q1 PF01476=50.0", "q2"])
    out = run(env, [a, b])
    counts = read_tsv(out / "lectin_report_counts.tsv")
    c = {(r["proteome"], r["model_acc"]): r for r in counts}
    assert len(counts) == 2 * len(lr.MODELS)
    assert c[("specA", "PF07938")]["n_hit_proteins"] == "2"
    assert c[("specA", "PF07938")]["best_score"] == "250.0"
    assert c[("specA", "PF07938")]["n_proteins"] == "3"
    assert c[("specA", "PF01476")]["n_hit_proteins"] == "1"
    assert c[("specB", "PF07938")]["n_hit_proteins"] == "0"
    assert c[("specB", "PF07938")]["best_score"] == "NA"
    assert c[("specB", "PF01476")]["best_score"] == "50.0"
    assert c[("specB", "PF01476")]["n_proteins"] == "2"
    hits = read_tsv(out / "lectin_report_hits.tsv")
    assert len(hits) == 4
    for h in hits:
        if h["model_acc"] == "PF07938":
            assert h["label"] == "PF07938 (AAL-like / FleA-family) domain hit"
    assert "signal_peptide" not in hits[0]


def test_report_text(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF07938=250.0"])
    out = run(env, [a])
    md = (out / "lectin_report.md").read_text()
    assert "PF07938 (AAL-like / FleA-family) domain hit" in md
    assert lr.FIXED_TEXT in md
    assert "does not show that the protein is on the cell surface or secreted" in md
    assert "Lectin-like families (not called)" in md
    assert "paracoccin" in md.lower()
    assert "Lectin_leg-like" in md  # M3 note
    assert "makes no call" in md
    assert "(a2) group" not in md  # no signal peptide data
    assert not list(out.glob("*status*"))


def test_a2_group_only_with_signalp(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(
        a,
        [
            "sp|P11111|NAME1 PF07938=250.0",  # SP yes (matched on accession)
            "p2 PF07938=120.0",  # SP no
            "p3 PF07938=300.0",  # no SP data
            "p4 PF07938=60.0",  # below the floor, SP yes
        ],
    )
    sp = tmp / "sp.txt"
    sp.write_text("# SignalP-6.0\n# ID\tPrediction\nP11111\tSP\t0.9\np2\tOTHER\t1.0\np4\tSP\t0.9\n")
    out = run(env, [a], ["--signalp", str(sp)])
    md = (out / "lectin_report.md").read_text()
    assert "## (a2) group: PF07938 score >= 100 and a signal peptide" in md
    assert "1 with a signal peptide, 1 without, 1 with no signal peptide data" in md
    assert "| specA | sp|P11111|NAME1 | 250.0 |" in md
    assert "p4" not in md.split("## (a2) group")[1].split("## Notes")[0]
    hits = read_tsv(out / "lectin_report_hits.tsv")
    col = {h["protein"]: h["signal_peptide"] for h in hits}
    assert col == {"sp|P11111|NAME1": "yes", "p2": "no", "p3": "NA", "p4": "yes"}


def test_domain_below_ga_is_flagged(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF07938=250.0 DOM=10.0", "p2 PF07938=250.0"])
    out = run(env, [a])
    counts = {r["model_acc"]: r for r in read_tsv(out / "lectin_report_counts.tsv")}
    assert counts["PF07938"]["n_hit_proteins"] == "2"
    assert counts["PF07938"]["n_hit_domain_ge_ga"] == "1"


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd not installed")
def test_zst_input(env):
    import subprocess

    tmp = env[0]
    plain = tmp / "specZ.faa"
    write_fasta(plain, ["z1 PF07938=99.0"])
    subprocess.run(["zstd", "-q", str(plain), "-o", str(tmp / "specZ.faa.zst")], check=True)
    plain.unlink()
    out = run(env, [tmp / "specZ.faa.zst"])
    counts = {
        (r["proteome"], r["model_acc"]): r for r in read_tsv(out / "lectin_report_counts.tsv")
    }
    assert counts[("specZ", "PF07938")]["n_hit_proteins"] == "1"


def test_needs_tmpdir_or_scratch(env, monkeypatch):
    tmp, stub, pfam = env
    a = tmp / "specA.faa"
    write_fasta(a, ["p1"])
    monkeypatch.delenv("SCRATCH", raising=False)
    with pytest.raises(lr.ReportError):
        lr.main(["--out-dir", str(tmp / "o"), "--pfam-hmm", pfam, "--hmmsearch", stub, str(a)])


def test_missing_model_stops(env):
    tmp, stub, pfam = env
    Path(pfam).write_text("HMMER3/f\nNAME  Other\nACC   PF99999.1\nGA    1.0 1.0;\n//\n")
    a = tmp / "specA.faa"
    write_fasta(a, ["p1"])
    with pytest.raises(lr.ReportError):
        run(env, [a])


def test_duplicate_names_stop(env):
    tmp = env[0]
    (tmp / "d1").mkdir()
    (tmp / "d2").mkdir()
    write_fasta(tmp / "d1/x.faa", ["p1"])
    write_fasta(tmp / "d2/x.fasta", ["p1"])
    with pytest.raises(lr.ReportError):
        run(env, [tmp / "d1/x.faa", tmp / "d2/x.fasta"])


def test_signalp6_file_with_header_in_id_column(tmp_path):
    sp = tmp_path / "sp.txt"
    sp.write_text(
        "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 1\n"
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        "p1 Genus_species_strain\tSP\t0.0002\t0.9998\tCS pos: 21-22. Pr: 0.98\n"
        "p2 Genus_other\tOTHER\t1.0\t0.0\t\n"
    )
    assert lr.read_signalp(str(sp)) == {"p1": True, "p2": False}


def test_name_list_with_signalp(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF07938=250.0", "p2 PF07938=300.0", "p3 PF07938=70.0"])
    sp = tmp / "sp.txt"
    sp.write_text("p1\tSP\np2\tSP\np3\tOTHER\n")
    names = tmp / "names.txt"
    names.write_text("# six\np1\tGenus one\np2\tGenus two\npX\tnot there\n")
    out = run(env, [a], ["--signalp", str(sp), "--name-list", str(names)])
    md = (out / "lectin_report.md").read_text()
    a2 = md.split("## (a2) group")[1].split("## Notes")[0]
    assert "### Named proteins" in a2
    assert (
        "| p1 | Genus one | specA | 250.0 | yes | PF07938 (AAL-like / FleA-family) domain hit |"
        in a2
    )
    assert "| pX | not there | not in the input proteomes, or no PF07938 hit |" in a2
    assert "No paper on the function of these proteins was found" in a2
    assert "Kerr 2016" in a2
    assert "Orthology of these proteins to FleA is not checked" in a2


def test_name_list_without_signalp(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF07938=250.0"])
    names = tmp / "names.txt"
    names.write_text("p1\tGenus one\n")
    out = run(env, [a], ["--name-list", str(names)])
    md = (out / "lectin_report.md").read_text()
    assert "## Named proteins" in md
    assert "(a2) group: PF07938" not in md  # the (a2) section needs signal peptide data
    assert "| p1 | Genus one | specA | 250.0 | NA |" in md


def test_duplicate_ids_warn_not_fail(env, capsys):
    tmp = env[0]
    a = tmp / "specA.faa"
    b = tmp / "specB.faa"
    write_fasta(a, ["same PF07938=250.0", "u1"])
    write_fasta(b, ["same PF07938=150.0", "u2"])
    sp = tmp / "sp.txt"
    sp.write_text("same\n")
    out = run(env, [a, b], ["--signalp", str(sp)])
    err = capsys.readouterr().err
    assert "WARNING: 1 protein IDs occur in more than one input proteome" in err
    assert "Warning:" in (out / "lectin_report.md").read_text()


def test_no_duplicate_warning_when_ids_differ(env, capsys):
    tmp = env[0]
    a = tmp / "specA.faa"
    b = tmp / "specB.faa"
    write_fasta(a, ["a1"])
    write_fasta(b, ["b1"])
    run(env, [a, b])
    assert "WARNING" not in capsys.readouterr().err


@pytest.mark.parametrize(
    "text",
    [
        # SignalP 4.1 short format: a number in column 2
        "# name\tCmax\tpos\tSmax\tpos\tSmean\tD\t?\tDmaxCut\tNetworks-used\np1\t0.110\t19\t0.2\t5\t0.1\t0.15\tN\t0.45\tSignalP-noTM\n",
        # yes/no table
        "p1\tyes\np2\tno\n",
        # SignalP 4 or 5 style words that are not in the accepted set
        "p1\tSP(Sec/SPI)\t0.9\n",
        "p1\tNO_SP\t0.9\n",
        # lower case word
        "p1\tsp\n",
        # space separated two columns
        "p1 SP\n",
        # tab line without a prediction word
        "p1\t\n",
    ],
)
def test_signalp_bad_formats_stop(tmp_path, text):
    f = tmp_path / "sp.txt"
    f.write_text(text)
    with pytest.raises(lr.ReportError) as exc:
        lr.read_signalp(str(f))
    assert "line" in str(exc.value)


def test_signalp_crlf_and_all_words(tmp_path):
    f = tmp_path / "sp.txt"
    f.write_bytes(
        b"# SignalP-6.0\r\n"
        b"a\r\nb\r\n"
        b"c hdr\tSP\t0.9\t0.1\tCS pos: 21-22. Pr: 0.98\r\n"
        b"d\tOTHER\t1.0\r\n"
        b"e\tLIPO\r\nf\tTAT\r\ng\tTATLIPO\r\nh\tPILIN\r\n"
    )
    assert lr.read_signalp(str(f)) == {
        "a": True, "b": True, "c": True, "d": False,
        "e": True, "f": True, "g": True, "h": True,
    }  # fmt: skip
    assert lr.lookup_sp(lr.read_signalp(str(f)), "b") is True


def test_bad_signalp_file_stops_the_run(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1 PF07938=250.0"])
    sp = tmp / "sp.txt"
    sp.write_text("p1\tyes\n")
    with pytest.raises(lr.ReportError):
        run(env, [a], ["--signalp", str(sp)])


def test_endemic_note_text(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["p1"])
    md = (run(env, [a]) / "lectin_report.md").read_text()
    assert "found no fungal lectin literature" not in md
    assert "no fungal lectin gene with host-binding data" in md
    assert "PMID 41805922" in md
    assert "CFP28" in md and "PMID 25023485" in md
    assert "Paracoccin (Paracoccidioides) is a known gap" in md


def test_duplicate_id_inside_one_proteome(env, capsys):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["dup PF07938=100.0", "dup PF07938=250.0", "u1"])
    out = run(env, [a])
    err = capsys.readouterr().err
    assert "specA: 1 records repeat an ID inside the file" in err
    counts = {r["model_acc"]: r for r in read_tsv(out / "lectin_report_counts.tsv")}
    assert counts["PF07938"]["n_hit_proteins"] == "1"
    assert counts["PF07938"]["best_score"] == "250.0"
    hits = read_tsv(out / "lectin_report_hits.tsv")
    assert len(hits) == 1
    assert "repeat an ID inside the file" in (out / "lectin_report.md").read_text()


def test_empty_fasta_does_not_abort(env, capsys):
    tmp = env[0]
    a = tmp / "specA.faa"
    e = tmp / "specEmpty.faa"
    write_fasta(a, ["p1 PF07938=250.0"])
    e.write_text("")
    out = run(env, [e, a])
    assert "specEmpty: the FASTA file has no sequences" in capsys.readouterr().err
    counts = {
        (r["proteome"], r["model_acc"]): r for r in read_tsv(out / "lectin_report_counts.tsv")
    }
    assert counts[("specEmpty", "PF07938")]["n_proteins"] == "0"
    assert counts[("specEmpty", "PF07938")]["n_hit_proteins"] == "0"
    assert counts[("specA", "PF07938")]["n_hit_proteins"] == "1"
    assert "no sequences" in (out / "lectin_report.md").read_text()


def test_name_list_bare_accession_and_crlf(env):
    tmp = env[0]
    a = tmp / "specA.faa"
    write_fasta(a, ["sp|P11111|NAME1 PF07938=250.0", "p2 PF07938=150.0"])
    names = tmp / "names.txt"
    names.write_bytes(b"P11111\tbare accession\r\nsp|P11111|NAME1\tfull ID\r\np2\tplain\r\n")
    out = run(env, [a], ["--name-list", str(names)])
    md = (out / "lectin_report.md").read_text()
    assert "| P11111 | bare accession | specA | 250.0 |" in md
    assert "| sp|P11111|NAME1 | full ID | specA | 250.0 |" in md
    assert "| p2 | plain | specA | 150.0 |" in md


def test_proteome_name():
    assert lr.proteome_name("/a/b/Sp_one.proteins.fa") == "Sp_one"
    assert lr.proteome_name("UP000000560.fasta.gz") == "UP000000560"
    assert lr.proteome_name("x.faa.zst") == "x"
