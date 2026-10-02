"""Static checks of the Phase B SLURM scripts and a J1 smoke run with a stub SignalP."""

import gzip
import json
import os
import subprocess

import paths
import pytest
import truth_table
from conftest import load_script

JOBS = paths.STEP1_DIR / "jobs"
SCRIPTS = sorted(JOBS.glob("*.sh"))
STUB = paths.STEP1_DIR.parents[1] / "tests/step1_compare/fixtures/phaseb/stub_signalp"


def test_three_job_scripts_exist():
    assert [s.name for s in SCRIPTS] == ["j0_pilot.sh", "j1_features.sh", "j2_embed.sh"]


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_job_script_conventions(script):
    text = script.read_text()
    assert "BASH_SOURCE" not in text
    assert '"${SCRATCH:?' in text
    assert ': "${PROJ_ROOT:?' in text and ': "${STEP1_WORKDIR:?' in text
    assert "#SBATCH -p exfab" in text and "#SBATCH --gres=gpu:1" in text
    assert "#SBATCH --time=" in text
    assert text.startswith("#!/bin/bash -l\n")
    assert "set -euo pipefail" in text
    subprocess.run(["bash", "-n", str(script)], check=True)


def test_j1_stub_signalp_is_executable():
    assert os.access(STUB / "bin" / "signalp6", os.X_OK)


needs_modules = pytest.mark.skipif(
    "BASH_FUNC_module%%" not in os.environ
    or not os.path.exists("/opt/linux/rocky/8.x/x86_64/modules/predgpi/202001"),
    reason="needs the HPCC module system and predgpi/202001",
)


def _tiny_work(tmp_path):
    work = tmp_path / "work"
    (work / "downloads").mkdir(parents=True)
    seqs = {
        "G1": "MKLLSVLALLLAAGSAQAS" + "ST" * 40 + "GAAAGLLSLLAALLF",
        "G2": "MSEEKKQ" * 12,
        "G3": "MKV" * 15,
    }
    (work / "downloads" / "p.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    sets = tmp_path / "sets.tsv"
    sets.write_text("set_id\tkind\tlocation\tnote\nP\tdownload\tp.fasta\t\n")
    prepare = load_script("05_prepare_sequences")
    argv = ["--sets", str(sets), "--work-dir", str(work), "--input-dir", str(work / "downloads")]
    assert prepare.main(argv) == 0
    return work


def _run_j1(work, tmp_path, **extra):
    env = {
        **os.environ,
        "SCRATCH": str(tmp_path / "scratch"),
        "PROJ_ROOT": str(paths.STEP1_DIR.parents[1]),
        "STEP1_WORKDIR": str(work),
        "J1_PARTS": "2",
        "J1_SIGNALP_MODULE": str(STUB / "modules" / "signalp" / "6-gpu"),
        **extra,
    }
    (tmp_path / "scratch").mkdir(exist_ok=True)
    script = JOBS / "j1_features.sh"
    return subprocess.run(["bash", str(script)], env=env, capture_output=True, text=True)


@needs_modules
def test_j1_runs_resumes_and_feeds_07(tmp_path):
    work = _tiny_work(tmp_path)
    first = _run_j1(work, tmp_path)
    assert first.returncode == 0, first.stdout[-2000:] + first.stderr[-2000:]
    out = work / "phaseb"
    # Part membership: record k of N (row order) goes to part floor((k - 1) * parts / N).
    # With 3 records and 2 parts: rows 0 and 1 in part_000, row 2 in part_001.
    hashes = [r["seq_sha256"] for r in truth_table.read_tsv(out / "unique_sequences.tsv.gz")]
    expected = {"part_000": hashes[:2], "part_001": hashes[2:]}
    for part, want in expected.items():
        sp = gzip.decompress((out / "signalp" / part / "prediction_results.txt.gz").read_bytes())
        sp_ids = [x.split("\t")[0] for x in sp.decode().splitlines() if not x.startswith("#")]
        gpi = gzip.decompress((out / "predgpi" / f"{part}.tsv.gz").read_bytes())
        gpi_ids = [x.split("\t")[0] for x in gpi.decode().splitlines()[1:]]
        assert sp_ids == want and gpi_ids == want, part
    again = _run_j1(work, tmp_path)
    assert again.returncode == 0 and again.stdout.count("skip") == 4
    header = "source_id\tgene_id\tlabel\tsubset\tstratum\td8_class\thomology_only\trole\n"
    (work / "truth_set_triaged.tsv.gz").write_bytes(gzip.compress(header.encode()))
    run = json.dumps({"all_sources": True, "truth_set_sha256": "0" * 64})
    (work / "d8_run.json").write_text(run)
    (work / "sequence_run.json").write_text(run)  # 07 checks that both name the same truth set
    prepare = json.loads((out / "prepare_run.json").read_text())
    prepare["truth_set_sha256"] = "0" * 64  # 05 copies this from sequence_run.json
    (out / "prepare_run.json").write_text(json.dumps(prepare))
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work), "--sets", str(work.parent / "sets.tsv")]) == 0
    cov = truth_table.read_tsv(out / "feature_coverage.tsv")
    assert cov == [
        {
            "set_id": "P",
            "members": "3",
            "unique_sequences": "3",
            "no_signalp": "0",
            "no_predgpi": "0",
            "over_1022": "0",
            "gpi_too_short": "0",
        }
    ]


@needs_modules
def test_j1_signalp_failure_leaves_no_done_marker(tmp_path):
    work = _tiny_work(tmp_path)
    run = _run_j1(work, tmp_path, STUB_SIGNALP_FAIL="1")
    assert run.returncode != 0
    assert not list((work / "phaseb" / "signalp").glob("part_*/prediction_results.txt.gz"))
    assert "stub signalp6 failure" in run.stderr
    assert "SignalP failed for part_000" in run.stderr


def test_j2_script_passes_the_planned_job_count():
    text = (JOBS / "j2_embed.sh").read_text()
    assert 'JOB_COUNT="${J2_JOB_COUNT:-1}"' in text
    assert '--job-count "$JOB_COUNT"' in text and '--job-index "$JOB_INDEX"' in text


def test_j0_keeps_throughput_when_the_diff_step_fails(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    (tmp_path / "scratch").mkdir()
    fake_py = tmp_path / "fake_python"
    fake_py.write_text(
        "#!/bin/bash\n"
        'out=""; while [ $# -gt 0 ]; do [ "$1" = --out ] && out=$2; shift; done\n'
        'case "$out" in\n'
        '  *throughput.json) echo "{}" > "$out" ;;\n'
        "  *) echo 'stub diff failure' >&2; exit 1 ;;\n"
        "esac\n"
    )
    fake_py.chmod(0o755)
    env = {
        **os.environ,
        "SCRATCH": str(tmp_path / "scratch"),
        "PROJ_ROOT": str(paths.STEP1_DIR.parents[1]),
        "STEP1_WORKDIR": str(work),
        "STEP1_ENV_PY": str(fake_py),
        "PATH": f"{STUB / 'bin'}:{os.environ['PATH']}",
    }
    run = subprocess.run(
        ["bash", str(JOBS / "j0_pilot.sh")], env=env, capture_output=True, text=True
    )
    assert run.returncode != 0
    assert (work / "phaseb" / "j0" / "throughput.json").exists()
    assert (work / "phaseb" / "j0" / "nvidia_smi.csv").exists()
    assert not (work / "phaseb" / "j0" / "gpu_cpu_diff.json").exists()


def test_j1_stops_on_empty_input(tmp_path):
    work = tmp_path / "work"
    (work / "phaseb").mkdir(parents=True)
    (work / "phaseb" / "unique_sequences.fasta.gz").write_bytes(gzip.compress(b""))
    run = _run_j1(work, tmp_path)
    assert run.returncode == 2
    assert "STOP: no sequences" in run.stderr


@needs_modules
def test_j1_stops_without_a_gpu(tmp_path):
    work = _tiny_work(tmp_path)
    run = _run_j1(work, tmp_path, STUB_NVIDIA_FAIL="1")
    assert run.returncode != 0
    assert "No devices were found" in run.stderr
    assert not list((work / "phaseb" / "signalp").glob("part_*/prediction_results.txt.gz"))
