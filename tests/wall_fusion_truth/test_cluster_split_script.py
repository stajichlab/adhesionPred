"""Static checks of cluster_split.sh (a draft: it is never run in a test) and of fetch_seeds_models.py."""

import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SH = ROOT / "analysis/wall_fusion_truth/cluster_split.sh"


def test_script_parses_and_follows_the_job_conventions():
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    text = SH.read_text()
    code = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))
    assert text.startswith("#!/bin/bash -l")
    assert "set -euo pipefail" in code
    assert "#SBATCH -p short" in text
    assert '"${SCRATCH:?' in code
    assert "/tmp" not in code and "BASH_SOURCE" not in code
    assert "DRAFT, NOT SUBMITTED" in text.splitlines()[1]
    assert "easy-cluster" in code and "--min-seq-id 0.3 -c 0.5" in code
    assert "module load MMseqs2/17-b804f" in code
    assert '[ ! -e "$OUT" ]' in code  # an existing output is never overwritten
    assert "check_pin" in code and "sidecar_ok" in code


def test_script_names_no_negative_source():
    text = SH.read_text()
    assert "hard_negative_seeds" not in text and "hard-negatives-shared" not in text


def load_fetch():
    spec = importlib.util.spec_from_file_location(
        "fetch_seeds_models", ROOT / "analysis/wall_fusion_truth/fetch_seeds_models.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fetch_seeds_models"] = mod
    spec.loader.exec_module(mod)
    return mod


def test_parse_stockholm_ungaps_and_joins_blocks():
    fm = load_fetch()
    txt = "#=GF AC   PF99999.1\n#=GF ID   Test_fam\nseq1/1-9  AC-DE.F\nseq2/1-9  MKL--V\nseq1/1-9  GHI\n//\n"
    ver, name, seqs = fm.parse_stockholm(txt)
    assert (ver, name) == ("PF99999.1", "Test_fam")
    assert seqs == {"seq1/1-9": "ACDEFGHI", "seq2/1-9": "MKLV"}


def test_model_list_is_the_wall_models_and_the_active_family_rows():
    import csv

    fm = load_fetch()
    with fm.WALL_MODELS.open(newline="") as fh:
        wall = [r["id"] for r in csv.DictReader(fh, delimiter="\t")]
    with fm.FAMILY_TABLE.open(newline="") as fh:
        fam = [r["pfam_acc"] for r in csv.DictReader(fh, delimiter="\t") if r["active"] == "yes"]
    got = fm.model_list()
    assert got == sorted(set(wall) | set(fam))
    assert len(wall) == 21 and len(fam) == 18 and len(got) == len(set(wall) | set(fam))
