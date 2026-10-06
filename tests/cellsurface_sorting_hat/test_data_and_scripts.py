"""The shipped family table, the Phase C species table and the job scripts."""

import csv
import os
import re
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
# Every check reads the CODE lines only: a comment that names a flag must not satisfy a test.


def code(name_or_text):
    text = (SH / name_or_text).read_text() if (SH / name_or_text).exists() else name_or_text
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))


def line_with(text, start):
    lines = [ln for ln in text.splitlines() if ln.lstrip().startswith(start)]
    assert len(lines) == 1, (start, lines)
    return lines[0]


def check_blast(text):
    c = code(text)
    blastp = c[c.index("blastp -query") :].split("-out ", 1)[0]
    for flag in ("-max_target_seqs 200", "-evalue 1 ", "-seg no"):
        assert flag in blastp, flag
    assert 'ge "$CAP"' in c and "CAP=200" in c  # the database must be smaller than the cap
    assert 'mv "$OUT/.tmp.blast.tsv" "$OUT/blast.tsv"' in c  # tmp copy, then atomic mv
    assert 'rm -f "$OUT/blast.tsv.done"' in c
    assert "blast.tsv.done" in c and "queries=" in c and "no output file" in c


def check_pfam(text):
    c = code(text)
    hmmsearch = line_with(c, "hmmsearch ")
    assert "--cut_ga" in hmmsearch and "--domtblout" in hmmsearch
    tail = line_with(c, "tail -n 3")
    assert "'# \\[ok\\]'" in tail and "FATAL" in tail
    assert "hmmfetch -f" in c and "grep -c '^NAME'" in c and "hmmfetch found" in c


def check_python_calls(text):
    for ln in code(text).splitlines():
        if re.search(r"(^|[\s;&|(])python3?(\s|$)", ln) or ".py" in ln:
            assert "/usr/bin/python3.12 " in ln, ln


def test_blast_job_and_its_documented_module_command_use_the_same_blast_settings():
    """The module command requires --evalue, --seg and --max-target-seqs; they must match blastp."""
    text = (SH / "blast_allergen.sbatch").read_text()
    check_blast(text)
    # the documented module command (a comment) carries the same values
    assert "--max-target-seqs 200 --evalue 1 --seg no" in text
    assert "must exceed the database size" in text  # leave-species-out recall depends on it
    assert "111 sequences" in text
    assert "must exceed the database size" in (SH / "submit_modules.sh").read_text()


def test_pfam_job_runs_hmmsearch_with_cut_ga_and_checks_the_ok_trailer():
    check_pfam((SH / "pfam_hmmsearch.sbatch").read_text())


def test_python_is_always_called_as_usr_bin_python3_12():
    for script in SCRIPTS:
        check_python_calls(script.read_text())
    assert "/usr/bin/python3.12" in code("repeats.sbatch")


def test_gpu_is_requested_only_by_the_signalp_job():
    for script in SCRIPTS:
        has = "#SBATCH --gres=gpu" in script.read_text()
        assert has == (script.name == "signalp_gpu.sbatch"), script.name


def test_submit_modules_starts_every_job_and_each_has_a_script():
    text = code("submit_modules.sh")
    jobs = text.split("for job in ", 1)[1].split(";", 1)[0].split()
    assert jobs == ["signalp_gpu", "pfam_hmmsearch", "repeats", "blast_allergen", "tmhmm"]
    assert {p.stem for p in SCRIPTS if p.suffix == ".sbatch"} == set(jobs)


def test_blast_job_writes_a_completion_marker_and_fails_when_there_is_no_output():
    c = code("blast_allergen.sbatch")
    assert "blast.tsv.done" in c and "queries=" in c
    assert "no output file" in c and "FATAL" in c


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


def test_repeat_job_passes_the_detector_min_len_of_the_module_code():
    from cellsurface_sorting_hat.modules.repeats import DETECTOR_MIN_LEN

    for ln in code("repeats.sbatch").splitlines():
        if ".py" in ln:
            assert f"--min-len {DETECTOR_MIN_LEN} " in ln, ln


def test_signalp_job_records_module_and_version_without_swallowing_errors():
    c = code("signalp_gpu.sbatch")
    assert "module=%s" in c and "signalp6 --version" in c and "|| true" not in c
    assert 'mv "$WORKDIR/raw/signalp/.tmp.version.txt"' in c


def test_fetch_script_does_not_claim_provenance_and_cleans_tmp_files_on_failure():
    text = (SH / "fetch_proteomes.sh").read_text()
    assert "and write provenance" not in text.split("set -euo", 1)[0]
    assert "trap 'rm -f \"$DEST\"/.tmp.*' EXIT" in code("fetch_proteomes.sh")


def test_tmhmm_job_makes_paths_absolute_before_it_changes_directory():
    c = code("tmhmm.sbatch")
    assert c.index("realpath") < c.index('OUT="') < c.index('cd "$TMP"') < c.index("tmhmm -short")
    assert 'FASTA="$(realpath "$FASTA")"' in c and 'WORKDIR="$(realpath -m "$WORKDIR")"' in c


@pytest.mark.parametrize(
    "script,check,needle",
    [
        ("blast_allergen.sbatch", check_blast, "-max_target_seqs 200"),
        ("blast_allergen.sbatch", check_blast, "-evalue 1 "),
        ("blast_allergen.sbatch", check_blast, "-seg no"),
        ("blast_allergen.sbatch", check_blast, 'if [ "$NDB" -ge "$CAP" ]'),
        ("blast_allergen.sbatch", check_blast, 'rm -f "$OUT/blast.tsv.done"'),
        ("blast_allergen.sbatch", check_blast, 'mv "$OUT/.tmp.blast.tsv"'),
        ("pfam_hmmsearch.sbatch", check_pfam, "hmmsearch --cut_ga"),
        ("pfam_hmmsearch.sbatch", check_pfam, "tail -n 3"),
        ("pfam_hmmsearch.sbatch", check_pfam, "hmmfetch found"),
        ("repeats.sbatch", check_python_calls, "/usr/bin/python3.12 "),
    ],
)
def test_each_new_check_fails_when_its_line_is_removed(script, check, needle):
    """Mutate an in-memory copy of the script: delete or weaken the line that has the needle."""
    text = (SH / script).read_text()
    check(text)  # passes on the real script
    mutated = []
    for ln in text.splitlines():
        if needle in ln and not ln.lstrip().startswith("#"):
            if check is check_python_calls:
                ln = ln.replace("/usr/bin/python3.12 ", "python3 ")
            else:
                continue
        mutated.append(ln)
    assert mutated != text.splitlines()
    with pytest.raises((AssertionError, ValueError)):
        check("\n".join(mutated))


def test_comments_alone_do_not_satisfy_the_blast_and_pfam_checks():
    for script, check in (
        ("blast_allergen.sbatch", check_blast),
        ("pfam_hmmsearch.sbatch", check_pfam),
    ):
        text = (SH / script).read_text()
        commented = "\n".join(
            "# " + ln if not ln.startswith("#") else ln for ln in text.splitlines()
        )
        with pytest.raises((AssertionError, ValueError)):
            check(commented)


# --- dry runs with stub tools (nothing is submitted, no real tool runs) ---


def _stub(bindir, name, body):
    path = bindir / name
    path.write_text("#!/bin/bash\n" + body)
    path.chmod(0o755)


def _env(tmp_path, **extra):
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    scratch = tmp_path / "scratch"
    scratch.mkdir(exist_ok=True)
    env = {k: v for k, v in os.environ.items() if k not in ("FASTA", "WORKDIR")}
    env.update(PATH=f"{bindir}:{env['PATH']}", SCRATCH=str(scratch), **extra)
    return bindir, env


def _run_script(script, env, cwd):
    # "module" is a shell function on the cluster; a stub function stands in for it
    cmd = f"module() {{ :; }}; export -f module; bash {SH / script}"
    return subprocess.run(["bash", "-c", cmd], env=env, cwd=cwd, capture_output=True, text=True)


def _blast_fixture(tmp_path, n_db):
    (tmp_path / "db.faa").write_text("".join(f">a{i}\nMKV\n" for i in range(n_db)))
    (tmp_path / "q.faa").write_text(">q1\nMKV\n>q2\nMKV\n")
    bindir, env = _env(
        tmp_path,
        FASTA=str(tmp_path / "q.faa"),
        WORKDIR=str(tmp_path / "w"),
        ALLERGEN_FASTA=str(tmp_path / "db.faa"),
    )
    calls = tmp_path / "calls.txt"
    _stub(bindir, "makeblastdb", f"echo makeblastdb >> {calls}\n")
    _stub(
        bindir,
        "blastp",
        f'if [ "$1" = -version ]; then echo "blastp: 2.14.0+"; exit 0; fi\n'
        f"echo blastp >> {calls}\n"
        'while [ $# -gt 0 ]; do [ "$1" = -out ] && out="$2"; shift; done\n'
        'printf "q1\\ta0\\t50\\t3\\t3\\t3\\t30\\t0.1\\n" > "$out"\n',
    )
    return env, calls


def test_blast_job_stops_before_blastp_when_the_database_reaches_the_cap(tmp_path):
    env, calls = _blast_fixture(tmp_path, 200)
    r = _run_script("blast_allergen.sbatch", env, tmp_path)
    assert r.returncode != 0 and "must exceed the database size" in r.stderr
    assert not calls.exists()  # blastp and makeblastdb were not called
    assert not (tmp_path / "w" / "raw" / "allergen" / "blast.tsv").exists()


def test_blast_job_with_a_small_database_runs_blastp_and_writes_marker_and_removes_old_one(
    tmp_path,
):
    env, calls = _blast_fixture(tmp_path, 111)
    out = tmp_path / "w" / "raw" / "allergen"
    out.mkdir(parents=True)
    (out / "blast.tsv.done").write_text("old")
    r = _run_script("blast_allergen.sbatch", env, tmp_path)
    assert r.returncode == 0, r.stderr
    assert "blastp" in calls.read_text()
    done = (out / "blast.tsv.done").read_text()
    assert "queries=2" in done and "blast=blastp: 2.14.0+" in done
    assert (out / "blast.tsv").stat().st_size > 0
    assert not list(out.glob(".tmp.*"))


def test_blast_job_fails_when_blastp_fails_and_leaves_no_marker(tmp_path):
    env, _ = _blast_fixture(tmp_path, 111)
    _stub(tmp_path / "bin", "blastp", "exit 3\n")
    r = _run_script("blast_allergen.sbatch", env, tmp_path)
    assert r.returncode != 0
    assert not (tmp_path / "w" / "raw" / "allergen" / "blast.tsv.done").exists()


def test_tmhmm_job_works_with_relative_paths(tmp_path):
    (tmp_path / "p.faa").write_text(">q1\nMKV\n")
    bindir, env = _env(tmp_path, FASTA="p.faa", WORKDIR="work")  # relative on purpose
    _stub(
        bindir,
        "tmhmm",
        'test -s "$2" || exit 9\necho "q1 len=3 ExpAA=0.1 First60=0.0 PredHel=0 Topology=o"\n',
    )
    r = _run_script("tmhmm.sbatch", env, tmp_path)
    assert r.returncode == 0, r.stderr
    lines = (tmp_path / "work" / "raw" / "tmhmm" / "tmhmm.tsv").read_text().splitlines()
    assert lines[0].startswith("protein_id") and lines[1].split("\t")[:2] == ["q1", "3"]
    assert not list((tmp_path / "scratch").glob("csh_tmhmm.*"))  # the trap removed it


def test_submit_modules_passes_absolute_paths_to_sbatch_and_supports_dry_run(tmp_path):
    for n in ("p.faa", "fam.tsv", "db.faa"):
        (tmp_path / n).write_text("x\n")
    bindir, env = _env(
        tmp_path,
        PROJ_ROOT=str(ROOT),
        FASTA="p.faa",
        WORKDIR="work",
        FAMILY_TABLE="fam.tsv",
        ALLERGEN_FASTA="db.faa",
        PFAM_RELEASE="38.2",
    )
    rec = tmp_path / "rec.txt"
    _stub(
        bindir,
        "sbatch",
        f'echo "$FASTA $WORKDIR $FAMILY_TABLE $ALLERGEN_FASTA $PROJ_ROOT" >> {rec}\necho 123\n',
    )
    r = subprocess.run(
        ["bash", str(SH / "submit_modules.sh")],
        env=env,
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    rows = rec.read_text().splitlines()
    assert len(rows) == 5
    for row in rows:
        assert all(f.startswith("/") for f in row.split())
    dry = subprocess.run(
        ["bash", str(SH / "submit_modules.sh")],
        env={**env, "DRY_RUN": "1"},
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert dry.returncode == 0 and dry.stdout.count("sbatch --parsable") == 5
    assert len(rec.read_text().splitlines()) == 5  # the dry run did not call sbatch
    assert str(tmp_path / "work") in dry.stdout
