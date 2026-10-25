"""Tests of stage_proteomes.py and proteome_array.sbatch (SYNTHETIC FASTA files only; no real sequence)."""

import hashlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "analysis/aga1_class"
spec = importlib.util.spec_from_file_location("stage_proteomes", D / "stage_proteomes.py")
sp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sp)

FA1 = ">sp|P00001|AAA_YEAST Alpha OS=Yeast GN=AAA PE=1\nMKSTA\nXUKS*\n>tr|A0A000|A0A000_YEAST Beta OS=Yeast\nMSSS\n"
FA2 = ">tr|Q99999|Q99999_CAND Gamma OS=Candida GN=GGG\nMTTTT\n"


@pytest.fixture
def fastas(tmp_path):
    a, b = tmp_path / "a.faa", tmp_path / "b.faa"
    a.write_text(FA1)
    b.write_text(FA2)
    return a, b


def test_stage_chunks_and_tables(fastas, tmp_path):
    out = tmp_path / "out"
    chunks, prots = sp.stage([("Sp1", fastas[0]), ("Sp2", fastas[1])], out, chunk_size=1)
    assert [c[0] for c in chunks] == ["Sp1_001", "Sp1_002", "Sp2_001"]
    assert (
        out / "chunks/Sp1_001.faa"
    ).read_text() == ">P00001\nMKSTAXUKS\n"  # bare accession, trailing * removed
    assert prots[0] == (
        "P00001",
        "Sp1",
        "Sp1_001",
        9,
        2,
        "AAA_YEAST",
        "AAA",
    )  # X and U counted, not changed
    assert prots[1][6] == "" and prots[2][0] == "Q99999"
    assert (out / "manifest.txt").read_text().split() == ["A0A000", "P00001", "Q99999"]
    for line in (out / "SHA256SUMS").read_text().splitlines():
        sha, name = line.split("  ")
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == sha
    for name, _species, _n, sha in (
        r.split("\t") for r in (out / "chunks.tsv").read_text().splitlines()[1:]
    ):
        assert hashlib.sha256((out / "chunks" / f"{name}.faa").read_bytes()).hexdigest() == sha


def test_stage_failures(fastas, tmp_path):
    with pytest.raises(sp.StageError, match="exists"):
        sp.stage([("S", fastas[0])], tmp_path)
    with pytest.raises(sp.StageError, match="given twice"):
        sp.stage([("S", fastas[0]), ("S", fastas[1])], tmp_path / "o1")
    dup = tmp_path / "dup.faa"
    dup.write_text(FA2.replace("Q99999", "P00001"))
    with pytest.raises(sp.StageError, match="occurs in"):
        sp.stage([("S1", fastas[0]), ("S2", dup)], tmp_path / "o2")
    for text, msg in (
        ("", "no records"),
        (">sp|P1|X_Y\n", "without a sequence"),
        (">bad header\nMK\n", "UniProt accession"),
        ("MK\n", "sequence before"),
    ):
        bad = tmp_path / "bad.faa"
        bad.write_text(text)
        with pytest.raises(sp.StageError, match=msg):
            sp.stage([("S", bad)], tmp_path / "ob")


def test_cli_refuses_existing_out_and_reports_fatal(fastas, tmp_path):
    out = tmp_path / "cli"
    ok = subprocess.run(
        [sys.executable, "-I", str(D / "stage_proteomes.py"), "--out", str(out), f"S={fastas[0]}"],
        capture_output=True,
        text=True,
    )
    assert ok.returncode == 0 and "2 proteins" in ok.stdout
    again = subprocess.run(
        [sys.executable, "-I", str(D / "stage_proteomes.py"), "--out", str(out), f"S={fastas[0]}"],
        capture_output=True,
        text=True,
    )
    assert again.returncode == 1 and again.stderr.startswith("FATAL:")


def test_array_sbatch_guards():
    p = D / "proteome_array.sbatch"
    assert subprocess.run(["bash", "-n", str(p)], capture_output=True).returncode == 0
    t = p.read_text()
    assert "check_pin()" in t and 'git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD' in t
    assert "${SCRATCH:?" in t and "${SLURM_ARRAY_TASK_ID:?" in t
    assert (
        '[ ! -e "$D" ] || fail' in t and t.count('[ ! -e "$D" ] || fail') == 2
    )  # before the run and before the copy
    assert 'sha256sum < "$FA"' in t and "chunks.tsv" in t
    assert "#SBATCH -p short" in t
    hours = re.search(r"#SBATCH --time=(\d+):", t)
    assert hours and int(hours.group(1)) <= 2
    assert (
        "no CUDA device" in t and "--force-cpu" not in t
    )  # AIUPred never falls back to the CPU here
    assert 'predgpi.py -f "$FA" -o "$TMP/predgpi_raw.gff3" -m gff3' in t
    assert 'signalp6 --fastafile "$FA" --organism eukarya' in t and "--mode fast" in t
    assert "aiupred.py" in t and "-g 0" in t
    for needle in ("PLACEHOLDER", "CMD_CONFIRMED", "DRAFT, NOT SUBMITTED"):
        assert needle not in t
