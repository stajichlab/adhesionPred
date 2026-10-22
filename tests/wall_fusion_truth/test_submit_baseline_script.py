"""Static checks of submit_baseline.sbatch (a draft: it is never run in a test) and its GA-split step."""

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SH = ROOT / "analysis/wall_fusion_truth/submit_baseline.sbatch"


def code():
    return "\n".join(ln for ln in SH.read_text().splitlines() if not ln.lstrip().startswith("#"))


def test_script_parses_and_follows_the_job_conventions():
    assert subprocess.run(["bash", "-n", str(SH)], capture_output=True).returncode == 0
    text = SH.read_text()
    c = code()
    assert text.startswith("#!/bin/bash -l")
    assert "DRAFT, NOT SUBMITTED" in text.splitlines()[1]
    assert "set -euo pipefail" in c
    assert "#SBATCH -p short" in text
    assert '"${SCRATCH:?' in c
    assert "/tmp" not in c and "BASH_SOURCE" not in c
    assert "check_pin" in c and "sidecar_ok" in c
    assert '[ ! -e "$OUT" ]' in c  # an existing output is never overwritten
    assert "module load hmmer/3.4" in c


def test_pfam_path_release_check_and_hash_are_recorded_not_fatal():
    c = code()
    assert "/srv/projects/db/pfam/current/Pfam-A.hmm" in c
    assert 'readlink -f "$PFAM_HMM"' in c
    assert '*"Pfam${PFAM_RELEASE}/"*' in c and "named *Pfam${PFAM_RELEASE}" in c
    # the release name stops the job; a different sha256 is only reported
    assert "REPORT: sha256" in c
    assert "4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e" in c
    assert "pfam_sha256_same_as_earlier_jobs" in c
    sha_block = c[c.index("HMM_SHA=") : c.index("/usr/bin/python3.12 -I -")]
    assert "exit" not in sha_block and "fail" not in sha_block


def test_models_are_fetched_by_name_and_checked_against_the_lists():
    c = code()
    assert "hmmfetch -f" in c and "N_EXP" in c and "-eq 39" in c
    assert "hmmsearch --cut_ga" in c
    assert "wall_fusion_models.tsv" in c and "family_table.tsv" in c
    assert "hard_negative" not in SH.read_text()
    assert "# [ok]" in SH.read_text() or "'# \\[ok\\]'" in c


def test_family_reps_are_searched_apart_and_nothing_is_called():
    c = code()
    assert "family_rep_sequences.faa" in c and "domtbl_famrep.txt" in c
    assert "calls" not in c and "status/" not in c and "write_status" not in c


def awk_program() -> str:
    text = SH.read_text()
    start = text.index("awk -v ga=")
    body = text[start : text.index('\' "$TMP/models.hmm"', start)]
    return body.split("'", 1)[1]


def test_a_model_without_ga_is_split_off_and_listed(tmp_path):
    blocks = (
        "HMMER3/f\nNAME  A\nACC   PF00001.5\nGA    27.0 27.0;\nHMM   x\n//\n"
        "HMMER3/f\nNAME  B\nACC   PF00002.7\nHMM   y\n//\n"
        "HMMER3/f\nNAME  C\nACC   PF00003.1\nGA    30.0 30.0;\nHMM   z\n//\n"
    )
    (tmp_path / "models.hmm").write_text(blocks)
    prog = awk_program()
    subprocess.run(
        ["awk", "-v", f"ga={tmp_path}/ga.hmm", "-v", f"noga={tmp_path}/noga.hmm",
         "-v", f"list={tmp_path}/noga.txt", prog, str(tmp_path / "models.hmm")],
        check=True,
    )  # fmt: skip
    assert (tmp_path / "ga.hmm").read_text().count("NAME") == 2
    assert (tmp_path / "noga.hmm").read_text().count("NAME") == 1
    assert (tmp_path / "noga.txt").read_text().split() == ["PF00002"]
    assert re.search(r"NAME\s+B", (tmp_path / "noga.hmm").read_text())
