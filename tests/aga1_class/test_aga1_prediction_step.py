"""Tests of the AGA1-class prediction step drafts (analysis/aga1_class/).

ALL tool outputs below are SYNTHETIC, hand-written, tiny, and only follow the documented or assumed layouts (UNVERIFIED until the
first real output is seen). No real sequence or prediction is used.
"""

import csv
import importlib.util
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "analysis/aga1_class"


def load(name):
    spec = importlib.util.spec_from_file_location(name, D / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


conv = load("convert_predictions")
stage = load("stage_inputs")

SIGNALP = (
    "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: SYNTHETIC\n"
    "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
    "P1 synthetic desc\tSP\t0.0100\t0.9900\tCS pos: 22-23. Pr: 0.9500\n"
    "P2\tOTHER\t0.9990\t0.0010\t\n"
    "P3\tOTHER\t0.8000\t0.2000\t\n"
)
PREDGPI = (
    "id\tprediction\tomega_site\tfp_rate\tscore\n"
    "P1\tGPI-anchored\t300\t0.0100\t12.5\n"
    "P2\tNot GPI-anchored\t\t55.0\t\n"
    "P3\tGPI-anchored\tNA\t0.5\t3.0\n"
)
DISORDER = (
    "# SYNTHETIC per-residue scores\n"
    ">P1 x\n1 M 0.1\n2 S 0.6\n3 T 0.7\n4 S 0.8\n5 A 0.2\n"
    ">P2\n1 M 0.0\n2 A 0.0\n"
    ">P3\n1 M 0.5\n2 A 0.5\n3 G 0.4\n"
)


@pytest.fixture
def files(tmp_path):
    def w(name, text):
        p = tmp_path / name
        p.write_text(text)
        return str(p)

    return {
        "manifest": w("m.txt", "P1\nP2\nP3\n"),
        "signalp": w("sp.txt", SIGNALP),
        "predgpi": w("gpi.tsv", PREDGPI),
        "disorder": w("dis.txt", DISORDER),
        "lengths": w("len.tsv", "accession\tlength\nP1\t5\nP2\t2\nP3\t3\n"),
        "tmp": tmp_path,
    }


def run(f, **over):
    a = {k: f[k] for k in ("manifest", "signalp", "predgpi", "disorder")}
    a.update(over)
    return conv.convert(
        a["manifest"],
        a["signalp"],
        a["predgpi"],
        a["disorder"],
        "aiupred",
        0.5,
        over.get("lengths"),
    )


def test_convert_values(files):
    rows = {r["accession"]: r for r in run(files)}
    assert list(rows) == ["P1", "P2", "P3"]
    p1 = rows["P1"]
    assert (p1["signal_peptide_pred"], p1["sp_cleavage_pos"], p1["sp_cleavage_prob"]) == (
        "yes",
        22,
        0.95,
    )
    assert (p1["gpi_pred"], p1["gpi_omega_site"], p1["gpi_fp_rate"], p1["gpi_score"]) == (
        "yes",
        300,
        0.01,
        12.5,
    )
    assert p1["disorder_n_res"] == 5 and p1["disorder_longest_run"] == 3
    assert p1["disorder_mean"] == "0.4800" and p1["disorder_frac"] == "0.6000"
    assert rows["P2"]["signal_peptide_pred"] == "no" and rows["P2"]["sp_cleavage_pos"] == "NA"
    assert rows["P2"]["gpi_omega_site"] == "NA" and rows["P2"]["disorder_longest_run"] == 0
    assert (
        rows["P3"]["gpi_omega_site"] == "NA" and rows["P3"]["disorder_longest_run"] == 2
    )  # score == threshold counts


def test_cli_writes_table_and_checks_lengths(files):
    out = files["tmp"] / "t.tsv"
    args = ["--manifest", files["manifest"], "--signalp", files["signalp"], "--predgpi", files["predgpi"],
            "--disorder", files["disorder"], "--disorder-tool", "aiupred", "--lengths", files["lengths"], "--out", str(out)]  # fmt: skip
    assert conv.main(args) == 0
    rows = list(csv.reader(out.open(), delimiter="\t"))
    assert rows[0] == conv.OUT_COLS and len(rows) == 4
    assert conv.main(args) == 1  # refuses an existing output


def edit(f, key, old, new):
    p = Path(f[key])
    assert old in p.read_text()
    p.write_text(p.read_text().replace(old, new))


def test_fail_missing_accession(files, capsys):
    edit(files, "signalp", "P3\tOTHER\t0.8000\t0.2000\t\n", "")
    with pytest.raises(conv.ConvertError, match="signalp: 1 manifest accessions missing"):
        run(files)
    out = files["tmp"] / "t.tsv"
    rc = conv.main(["--manifest", files["manifest"], "--signalp", files["signalp"], "--predgpi", files["predgpi"],
                    "--disorder", files["disorder"], "--disorder-tool", "aiupred", "--out", str(out)])  # fmt: skip
    assert rc == 1 and not out.exists() and "FATAL" in capsys.readouterr().err


def test_fail_extra_accession(files):
    edit(files, "predgpi", "P3\t", "P9\t")
    with pytest.raises(conv.ConvertError, match="predgpi"):
        run(files)


@pytest.mark.parametrize(
    "key,old,new,msg",
    [
        ("signalp", "P2\tOTHER", "P1\tOTHER", "duplicate ID P1"),
        ("predgpi", "P2\tNot", "P1\tNot", "duplicate id P1"),
        ("disorder", ">P2\n", ">P1\n", "duplicate record P1"),
        ("signalp", "0.9900", "abc", "SP probability is not a number"),
        ("signalp", "CS pos: 22-23. Pr: 0.9500", "", "without a parsable"),
        ("signalp", "Eukarya", "Gram-positive", "not a Eukarya run"),
        ("signalp", "P2\tOTHER", "P2\tLIPO", "unexpected prediction"),
        ("signalp", "P2\tOTHER\t0.9990\t0.0010\t", "P2\tOTHER\t0.9990", "expected 5 TAB"),
        ("predgpi", "GPI-anchored\t300", "maybe\t300", "unexpected prediction"),
        ("predgpi", "0.0100", "nan", "not finite"),
        ("predgpi", "id\tprediction", "id\tcall", "header must start with"),
        ("disorder", "2 S 0.6", "2 S x", "not a number"),
        ("disorder", "3 T 0.7", "4 T 0.7", "out of order"),
        ("disorder", "2 S 0.6", "2 S 1.6", "outside [0, 1]"),
        ("disorder", "2 S 0.6", "2 S", "expected 'pos aa score'"),
    ],
)
def test_fail_unparseable(files, key, old, new, msg):
    edit(files, key, old, new)
    with pytest.raises(conv.ConvertError, match=re.escape(msg)):
        run(files)


@pytest.mark.parametrize("key", ["signalp", "predgpi", "disorder", "manifest"])
def test_fail_empty_file(files, key):
    Path(files[key]).write_text("")
    with pytest.raises(conv.ConvertError, match="is empty"):
        run(files)


def test_fail_header_only_files(files):
    Path(files["predgpi"]).write_text("id\tprediction\tomega_site\tfp_rate\n")
    with pytest.raises(conv.ConvertError, match="manifest accessions missing"):
        run(files)
    Path(files["signalp"]).write_text("# SignalP-6.0\tOrganism: Eukarya\n")
    with pytest.raises(conv.ConvertError, match="no prediction rows"):
        run(files)


def test_fail_length_mismatch(files):
    edit(files, "lengths", "P1\t5", "P1\t6")
    with pytest.raises(conv.ConvertError, match="differ in length"):
        run(files, lengths=files["lengths"])


def test_missing_file(files):
    with pytest.raises(conv.ConvertError, match="does not exist"):
        run(files, signalp=str(files["tmp"] / "nope.txt"))


# ---- staging manifest ----
def test_manifest_matches_proteins_tsv():
    acc = stage.accessions_from_proteins()
    text = (D / "manifest_accessions.txt").read_text()
    assert text == "\n".join(acc) + "\n"
    assert len(acc) == 850 == len(set(acc)) and acc == sorted(acc)
    assert all(re.fullmatch(r"[A-Za-z0-9_.]+", a) for a in acc)
    assert stage.main([]) == 0


def test_manifest_holds_no_sequence():
    for ln in (D / "manifest_accessions.txt").read_text().splitlines():
        assert len(ln) < 20


def test_build_fasta_synthetic_and_failures():
    fasta, rows = stage.build(["A1", "B2"], {"A1": "ma" * 40 + "*", "B2": "MXK"}, {"A1": "s"})
    lines = fasta.splitlines()
    assert lines[0] == ">A1" and len(lines[1]) == 60 and lines[2] == "A" * 0 + "MA" * 10
    assert rows == [("A1", 80, 0, "s"), ("B2", 3, 1, "")]
    with pytest.raises(SystemExit, match="no sequence"):
        stage.build(["A1", "Z"], {"A1": "MK"}, {})
    with pytest.raises(SystemExit, match="empty sequence"):
        stage.build(["A1"], {"A1": "*"}, {})


def test_stage_refuses_output_inside_repo_or_existing(tmp_path):
    with pytest.raises(SystemExit, match="outside the repo"):
        stage.main(["--out", str(ROOT / "analysis/aga1_class/staged_here")])
    with pytest.raises(SystemExit, match="exists"):
        stage.main(["--out", str(tmp_path)])


# ---- sbatch drafts ----
SBATCH = ["signalp.sbatch", "predgpi.sbatch", "disorder.sbatch"]


@pytest.mark.parametrize("name", SBATCH)
def test_sbatch_syntax_and_guards(name):
    p = D / name
    assert subprocess.run(["bash", "-n", str(p)], capture_output=True).returncode == 0
    t = p.read_text()
    assert "DRAFT, NOT SUBMITTED" in t
    assert "check_pin()" in t and 'git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD' in t
    assert "${SCRATCH:?" in t
    assert '[ ! -e "$OUT" ] || fail' in t
    assert "#SBATCH -p short" in t
    hours = re.search(r"#SBATCH --time=(\d+):", t)
    assert hours and int(hours.group(1)) <= 2
    assert "-o /path/to/logs/%x.%j.log" in t and "log directory must exist" in t
    assert "provenance.txt" in t


@pytest.mark.parametrize("name", SBATCH)
def test_no_placeholders_left(name):
    t = (D / name).read_text()
    assert "PLACEHOLDER" not in t and "CMD_CONFIRMED" not in t


def test_predgpi_command_and_home():
    t = (D / "predgpi.sbatch").read_text()
    assert 'export PREDGPI_HOME="$PREDGPI_DIR"' in t
    assert (
        '[ -f "$PREDGPI_HOME/predgpi.py" ] || fail' in t and '[ -d "$PREDGPI_HOME" ] || fail' in t
    )
    assert 'predgpi.py -f "$FA" -o "$TMP/predgpi_raw.gff3" -m gff3' in t


def test_disorder_commands_and_guards():
    t = (D / "disorder.sbatch").read_text()
    assert '"$AIUPRED_PY" -i "$FA" -o' in t
    assert "--force-cpu" in t and "ALLOW_CPU" in t and "import torch, numpy, scipy" in t
    assert t.index("ALLOW_CPU:-no") < t.index('"$TMP/disorder_raw.txt" --force-cpu')
    assert "iupred3.py" in t and " long" in t and "import scipy" in t
    assert "#SBATCH --gres" not in t  # the GPU partition is named by the owner on the command line


def test_signalp_records_versions():
    t = (D / "signalp.sbatch").read_text()
    assert "signalp6 --version" in t and "--organism eukarya" in t and "module -t list" in t
