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
    "P4\tSP\t0.1\t0.9\tCS pos: 19-20. Pr: 0.5000\n"
    "P5\tOTHER\t0.9\t0.1\t\n"
)
# SYNTHETIC gff3 in the layout of the real PredGPI output (P1 5 aa would be too short in reality; lengths below are synthetic too)
PREDGPI = (
    "P1\tPredGPI\tGPI-anchor\t60\t60\t1.0\t.\t.\tOntology_term=GO:0046658;evidence=ECO:0000256\n"
    "P2\tPredGPI\tChain\t1\t50\t1.0\t.\t.\tevidence=ECO:0000256\n"
    "P3\tPredGPI\tChain\t1\t39\t1.0\t.\t.\tevidence=ECO:0000256\n"
    "P4\tPredGPI\tGPI-anchor\t70\t70\t0.55\t.\t.\tOntology_term=GO:0046658;evidence=ECO:0000256\n"
    "P5\tPredGPI\tGPI-anchor\t80\t80\t0.7\t.\t.\tOntology_term=GO:0046658;evidence=ECO:0000256\n"
)
LENGTHS = "accession\tlength\nP1\t100\nP2\t50\nP3\t39\nP4\t90\nP5\t95\n"


def _rec(name, head, n):
    rows = [f"{i} A {sc}" for i, sc in enumerate(head, 1)] + [
        f"{i} A 0.0" for i in range(len(head) + 1, n + 1)
    ]
    return f">{name} x\n" + "\n".join(rows) + "\n"


DISORDER = (
    "# SYNTHETIC per-residue scores\n"
    + _rec("P1", [0.1, 0.6, 0.7, 0.8, 0.2], 100)
    .replace("2 A 0.6", "2 S 0.6")
    .replace("3 A 0.7", "3 T 0.7")
    + _rec("P2", [], 50)
    + _rec("P3", [0.5, 0.5, 0.4], 39)
    + _rec("P4", [], 90)
    + _rec("P5", [], 95)
)


@pytest.fixture
def files(tmp_path):
    def w(name, text):
        p = tmp_path / name
        p.write_text(text)
        return str(p)

    return {
        "manifest": w("m.txt", "P1\nP2\nP3\nP4\nP5\n"),
        "signalp": w("sp.txt", SIGNALP),
        "predgpi": w("gpi.tsv", PREDGPI),
        "disorder": w("dis.txt", DISORDER),
        "lengths": w("len.tsv", LENGTHS),
        "tmp": tmp_path,
    }


def run(f, disorder=True, **over):
    a = {k: f[k] for k in ("manifest", "signalp", "predgpi", "disorder", "lengths")}
    a.update(over)
    return conv.convert(
        a["manifest"],
        a["signalp"],
        a["predgpi"],
        a["disorder"] if disorder else None,
        "aiupred",
        0.5,
        a["lengths"],
    )


def test_convert_values(files):
    rows = {r["accession"]: r for r in run(files)}
    assert list(rows) == ["P1", "P2", "P3", "P4", "P5"]
    p1 = rows["P1"]
    assert (p1["signal_peptide_pred"], p1["sp_cleavage_pos"], p1["sp_cleavage_prob"]) == (
        "yes",
        22,
        0.95,
    )
    assert (p1["gpi_pred"], p1["gpi_omega_site"], p1["gpi_confidence_class"], p1["gpi_fpr_upper_bound"]) == (
        "yes", 60, "high", 0.0015,
    )  # fmt: skip
    assert p1["disorder_n_res"] == 100 and p1["disorder_longest_run"] == 3
    assert p1["disorder_mean"] == "0.0240" and p1["disorder_frac"] == "0.0300"
    assert rows["P2"]["signal_peptide_pred"] == "no" and rows["P2"]["sp_cleavage_pos"] == "NA"
    assert (
        rows["P2"]["gpi_pred"],
        rows["P2"]["gpi_omega_site"],
        rows["P2"]["gpi_confidence_class"],
    ) == ("no", "NA", "NA")
    assert rows["P2"]["gpi_fpr_upper_bound"] == "NA" and rows["P2"]["disorder_longest_run"] == 0
    assert rows["P3"]["disorder_longest_run"] == 2  # score == threshold counts
    assert (rows["P4"]["gpi_confidence_class"], rows["P4"]["gpi_fpr_upper_bound"]) == ("weak", 0.01)
    assert (rows["P5"]["gpi_confidence_class"], rows["P5"]["gpi_fpr_upper_bound"]) == (
        "probable",
        0.005,
    )
    assert not any(
        "fp_rate" in c or c == "gpi_score" for c in rows["P1"]
    )  # no invented false positive rate


def test_short_protein_is_na_not_no(files):
    rows = {r["accession"]: r for r in run(files)}
    # P3 is 39 aa: PredGPI never scores it, its Chain row is a silent "no" in the tool, so it must be NA here
    assert rows["P3"]["gpi_pred"] == "NA" and "length<=40" in rows["P3"]["gpi_note"]
    assert rows["P2"]["gpi_note"] == ""  # 50 aa: scored, a real "no"
    # exactly 40 aa is still too short, 41 is scored
    assert conv.PREDGPI_MIN_LEN == 41
    one = {"P": 40}, {"Q": 41}
    for lens, want in zip(one, ("NA", "no"), strict=True):
        acc = next(iter(lens))
        p = files["tmp"] / f"g_{acc}.gff3"
        p.write_text(f"{acc}\tPredGPI\tChain\t1\t{lens[acc]}\t1.0\t.\t.\te=1\n")
        assert conv.parse_predgpi(str(p), lens)[acc]["pred"] == want


def test_no_disorder_flag_leaves_columns_out(files):
    out = files["tmp"] / "t.tsv"
    args = ["--manifest", files["manifest"], "--signalp", files["signalp"], "--predgpi", files["predgpi"],
            "--lengths", files["lengths"], "--no-disorder", "--out", str(out)]  # fmt: skip
    assert conv.main(args) == 0
    head = out.read_text().splitlines()[0].split("\t")
    assert head == [c for c in conv.OUT_COLS if not c.startswith("disorder_")]
    assert "gpi_confidence_class" in head and "gpi_fp_rate" not in head


@pytest.mark.parametrize(
    "extra",
    [
        [],
        ["--no-disorder", "--disorder", "x", "--disorder-tool", "aiupred"],
        ["--disorder", "x"],
        ["--no-disorder", "--disorder-tool", "aiupred"],
    ],
)
def test_disorder_block_must_be_explicit(files, extra):
    args = ["--manifest", files["manifest"], "--signalp", files["signalp"], "--predgpi", files["predgpi"],
            "--lengths", files["lengths"], "--out", str(files["tmp"] / "t.tsv"), *extra]  # fmt: skip
    with pytest.raises(SystemExit) as e:
        conv.main(args)
    assert e.value.code == 2  # argparse error: never a silent default


def test_cli_writes_table_and_checks_lengths(files):
    out = files["tmp"] / "t.tsv"
    args = ["--manifest", files["manifest"], "--signalp", files["signalp"], "--predgpi", files["predgpi"],
            "--disorder", files["disorder"], "--disorder-tool", "aiupred", "--lengths", files["lengths"], "--out", str(out)]  # fmt: skip
    assert conv.main(args) == 0
    rows = list(csv.reader(out.open(), delimiter="\t"))
    assert rows[0] == conv.OUT_COLS and len(rows) == 6
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
                    "--disorder", files["disorder"], "--disorder-tool", "aiupred", "--lengths", files["lengths"], "--out", str(out)])  # fmt: skip
    assert rc == 1 and not out.exists() and "FATAL" in capsys.readouterr().err


def test_fail_extra_accession(files):
    edit(files, "predgpi", "P3\t", "P9\t")
    with pytest.raises(conv.ConvertError, match="predgpi"):
        run(files)


@pytest.mark.parametrize(
    "key,old,new,msg",
    [
        ("signalp", "P2\tOTHER", "P1\tOTHER", "duplicate ID P1"),
        ("predgpi", "P2\tPredGPI", "P1\tPredGPI", "duplicate ID P1"),
        ("disorder", ">P2 x\n", ">P1\n", "duplicate record P1"),
        ("signalp", "0.9900", "abc", "SP probability is not a number"),
        ("signalp", "CS pos: 22-23. Pr: 0.9500", "", "without a parsable"),
        ("signalp", "Eukarya", "Gram-positive", "not a Eukarya run"),
        ("signalp", "P2\tOTHER", "P2\tLIPO", "unexpected prediction"),
        ("signalp", "P2\tOTHER\t0.9990\t0.0010\t", "P2\tOTHER\t0.9990", "expected 5 TAB"),
        ("predgpi", "GPI-anchor\t60\t60", "Mystery\t60\t60", "unknown feature type"),
        ("predgpi", "GPI-anchor\t60\t60", "GPI-anchor\t60\t61", "start 60 != end 61"),
        ("predgpi", "GPI-anchor\t60\t60", "GPI-anchor\t160\t160", "outside the protein"),
        ("predgpi", "60\t60\t1.0", "60\t60\t0.9", "unknown GPI-anchor confidence score"),
        (
            "predgpi",
            "P2\tPredGPI\tChain\t1\t50",
            "P2\tPredGPI\tChain\t1\t49",
            "does not span the protein",
        ),
        (
            "predgpi",
            "P2\tPredGPI\tChain\t1\t50",
            "P2\tPredGPI\tChain\t2\t50",
            "does not span the protein",
        ),
        ("predgpi", "P2\tPredGPI\tChain", "P2\tOther\tChain", "expected PredGPI"),
        (
            "predgpi",
            "P2\tPredGPI\tChain\t1\t50\t1.0\t.\t.\t",
            "P2\tPredGPI\tChain\t1\t50\t1.0\t.\t",
            "expected 9 TAB",
        ),
        ("predgpi", "P2\tPredGPI\tChain\t1\t50", "P2\tPredGPI\tChain\tone\t50", "not integers"),
        (
            "predgpi",
            "P3\tPredGPI\tChain\t1\t39",
            "P3\tPredGPI\tGPI-anchor\t30\t30",
            "PredGPI should not score it",
        ),
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


@pytest.mark.parametrize("key", ["signalp", "predgpi", "disorder", "manifest", "lengths"])
def test_fail_empty_file(files, key):
    Path(files[key]).write_text("")
    with pytest.raises(conv.ConvertError, match="is empty"):
        run(files)


def test_fail_missing_gpi_row_and_no_rows(files):
    edit(files, "predgpi", "P2\tPredGPI\tChain\t1\t50\t1.0\t.\t.\tevidence=ECO:0000256\n", "")
    with pytest.raises(conv.ConvertError, match="predgpi: 1 manifest accessions missing"):
        run(files)
    Path(files["signalp"]).write_text("# SignalP-6.0\tOrganism: Eukarya\n")
    with pytest.raises(conv.ConvertError, match="no prediction rows"):
        run(files)


def test_fail_lengths_missing_accession(files):
    edit(files, "lengths", "P5\t95\n", "")
    with pytest.raises(conv.ConvertError, match="lengths: 1 manifest accessions missing"):
        run(files)


def test_fail_length_mismatch(files):
    edit(files, "lengths", "P1\t100", "P1\t101")
    with pytest.raises(conv.ConvertError, match="differ in length"):
        run(files)


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


def test_aiupred_hash_header_layout_is_read(files):
    """AIUPred 2.x writes a banner, then per record '#>ID', score rows and a blank line (real layout of 2026-10-24)."""
    plain = run(files)
    banner = "#  AIUPred (synthetic banner)\n# Position\tResidue\tDisorder\n"
    lines = Path(files["disorder"]).read_text().splitlines()
    body = "".join(
        ("#" + ln + "\n") if ln.startswith(">") else (ln + "\n")
        for ln in lines
        if not ln.startswith("# SYNTH")
    )
    Path(files["disorder"]).write_text(banner + body.replace("\n#>", "\n\n#>"))
    assert "#>P1" in Path(files["disorder"]).read_text()
    assert run(files) == plain


def test_small_overshoot_of_the_score_range_is_accepted(files):
    p = Path(files["disorder"])
    p.write_text(p.read_text().replace("2 S 0.6", "2 S 1.0096").replace("3 T 0.7", "3 T -0.0123"))
    assert run(files)[0]["disorder_n_res"] == 100
