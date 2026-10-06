"""The shipped family table, the Phase C species table and the job scripts."""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest

from cellsurface_sorting_hat.modules.pfam import load_family_table

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = sorted((ROOT / "scripts" / "sorting_hat").glob("*"))
SH = ROOT / "scripts" / "sorting_hat"


def test_shipped_family_table_loads_and_starts_with_every_family_inactive():
    families = load_family_table(ROOT / "data" / "sorting_hat" / "family_table.tsv")
    assert len(families) == 15
    assert not any(
        f.active for f in families
    )  # a family is made active only after its specificity test
    assert {f.module for f in families} == {"pfam_adhesion", "pfam_allergen"}


def test_phasec_species_table_has_one_row_per_species():
    with open(ROOT / "data" / "sorting_hat" / "phasec_set_species.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len({r["scientific_name"] for r in rows}) == len(rows) == 6


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_scripts_are_valid_bash_and_follow_the_site_rules(script):
    text = script.read_text()
    assert subprocess.run(["bash", "-n", str(script)], capture_output=True).returncode == 0
    assert "BASH_SOURCE" not in text  # breaks under sbatch
    if script.suffix == ".sbatch":
        assert "SCRATCH:?" in text and "#SBATCH" in text and "set -euo pipefail" in text
        assert "/tmp" not in text


def test_there_is_a_script_for_every_job_that_submit_modules_starts():
    names = {p.stem for p in SCRIPTS if p.suffix == ".sbatch"}
    text = (ROOT / "scripts" / "sorting_hat" / "submit_modules.sh").read_text()
    for job in text.split("for job in ", 1)[1].split(";", 1)[0].split():
        assert job in names


@pytest.mark.skipif(shutil.which("shellcheck") is None, reason="shellcheck is not installed")
@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_shellcheck(script):
    assert (
        subprocess.run(["shellcheck", "-S", "warning", str(script)], capture_output=True).returncode
        == 0
    )


def test_the_pfam_job_checks_model_names_against_the_family_table():
    text = (ROOT / "scripts" / "sorting_hat" / "pfam_hmmsearch.sbatch").read_text()
    assert "differ from the family table" in text and "hmmfetch -f" in text and "--cut_ga" in text


def test_the_pfam_job_fetches_by_model_name_and_the_table_names_are_unique():
    """hmmfetch -f on the pressed Pfam database needs the NAME; an accession without a version is refused."""
    text = (ROOT / "scripts" / "sorting_hat" / "pfam_hmmsearch.sbatch").read_text()
    assert "cut -f2 | sort -u" in text and "cut -f1 | sort -u" not in text
    families = load_family_table(ROOT / "data" / "sorting_hat" / "family_table.tsv")
    names = [f.name for f in families]
    assert len(names) == len(set(names))
    assert "Pfam${PFAM_RELEASE}/" in text  # the release must be a whole directory name


def test_every_sbatch_script_cleans_its_scratch_directory_on_exit():
    for script in SCRIPTS:
        if script.suffix == ".sbatch":
            assert "trap 'rm -rf \"$TMP\"' EXIT" in script.read_text(), script.name


# --- checks added for the items deferred from earlier reviews ---


def test_blast_job_and_its_documented_module_command_use_the_same_blast_settings():
    """The module command requires --evalue, --seg and --max-target-seqs; they must match blastp."""
    text = (SH / "blast_allergen.sbatch").read_text()
    assert "-max_target_seqs 200" in text and "-evalue 1 " in text and "-seg no" in text
    assert (
        "cellsurface_sorting_hat_module allergen" in text
        and "--max-target-seqs 200 --evalue 1 --seg no" in text
    )
    assert "must exceed the database size" in text  # leave-species-out recall depends on it
    assert "must exceed the database size" in (SH / "submit_modules.sh").read_text()


def test_blast_job_writes_a_completion_marker_and_fails_when_there_is_no_output():
    text = (SH / "blast_allergen.sbatch").read_text()
    assert "blast.tsv.done" in text and "queries=" in text
    assert "no output file" in text and "FATAL" in text


def test_repeat_job_documents_module_commands_that_pass_the_detector_script():
    text = (SH / "repeats.sbatch").read_text()
    for name, script in (
        ("repeat02", "02_repeat_profile.py"),
        ("repeat14", "14_repeat_detect_general.py"),
    ):
        lines = [ln for ln in text.splitlines() if f"cellsurface_sorting_hat_module {name}" in ln]
        assert lines and all(
            f"--script $PROJ_ROOT/analysis/cocci_repeats/{script}" in ln for ln in lines
        )


def test_pfam_job_runs_hmmsearch_with_cut_ga_and_domtblout_and_documents_the_pfam_command():
    text = (SH / "pfam_hmmsearch.sbatch").read_text()
    lines = [ln for ln in text.splitlines() if ln.startswith("hmmsearch ")]
    assert len(lines) == 1 and "--cut_ga" in lines[0] and "--domtblout" in lines[0]
    assert "# [ok]" in text  # the domain table trailer the module command requires


def test_tmhmm_job_runs_in_the_scratch_directory():
    text = (SH / "tmhmm.sbatch").read_text()
    assert text.index('cd "$TMP"') < text.index("tmhmm -short")
