"""Static checks of the Phase C shell scripts and a C1 smoke run with the stub MMseqs2."""

import json
import os
import re
import subprocess

import paths
import pytest

PHASEC = paths.STEP1_DIR / "phasec"
SCRIPTS = sorted(PHASEC.glob("*.sh"))
# node features of AVX2 CPUs that sbatch accepts on UCR HPCC (epyc nodes: ryzen, amd, milan)
AVX2_FEATURES = ("ryzen", "milan", "genoa", "rome")
NON_AVX2_MMSEQS = "/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs"


def needs_avx2(text: str) -> bool:
    """The script loads the MMseqs2 module, names the AVX2 build, or is a SLURM script that
    runs mmseqs."""
    loads = re.search(r"module\s+load\b[^\n]*MMseqs2", text) or "17-b804f-avx2" in text
    return bool(loads) or ("#SBATCH" in text and "mmseqs" in text.lower())


def has_avx2_constraint(text: str) -> bool:
    found = re.findall(r"^#SBATCH\s+--constraint=(\S+)", text, flags=re.M)
    return any(f in AVX2_FEATURES for c in found for f in re.split(r"[&|,]", c))


def uses_non_avx2_binary(text: str) -> bool:
    return re.search(rf'STEP1_MMSEQS="?{re.escape(NON_AVX2_MMSEQS)}"?', text) is not None


def test_two_shell_scripts_exist():
    assert [s.name for s in SCRIPTS] == ["09_cluster_and_split.sh", "c1_evaluate.sh"]


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_shell_script_conventions(script):
    text = script.read_text()
    assert text.startswith("#!/bin/bash -l\n")
    assert "set -euo pipefail" in text
    assert "BASH_SOURCE" not in text
    assert '"${SCRATCH:?' in text
    assert ': "${PROJ_ROOT:?' in text and ': "${STEP1_WORKDIR:?' in text
    subprocess.run(["bash", "-n", str(script)], check=True)


def test_c1_requests_cpu_on_an_avx2_partition_with_a_time_limit():
    text = (PHASEC / "c1_evaluate.sh").read_text()
    assert "#SBATCH -p epyc" in text and "#SBATCH --time=" in text
    assert "#SBATCH --constraint=ryzen" in text  # ruling C-15
    assert """trap 'cp "$TMP/wall.txt" "$WALL_OUT"' EXIT""" in text  # review M-1
    assert "--gres" not in text  # no GPU (spec 5)
    assert "OMP_NUM_THREADS=1" in text
    assert 'STEPS="${C1_STEPS:-09 10 11}"' in text
    assert "wall_seconds=" in text


def test_avx2_tools_have_a_cpu_constraint():
    # owner rule: a job that runs an AVX2 tool must request a node feature that has AVX2
    scripts = sorted(paths.STEP1_DIR.rglob("*.sh"))
    flagged = {s.relative_to(paths.STEP1_DIR).as_posix(): s.read_text() for s in scripts
               if needs_avx2(s.read_text())}  # fmt: skip
    assert {"phasec/c1_evaluate.sh", "phasec/09_cluster_and_split.sh"} <= set(flagged)
    # the Phase B jobs run on exfab GPU nodes and use no MMseqs2
    assert not set(flagged) & {"jobs/j0_pilot.sh", "jobs/j1_features.sh", "jobs/j2_embed.sh"}
    bad = [
        n for n, t in flagged.items() if not has_avx2_constraint(t) and not uses_non_avx2_binary(t)
    ]
    assert bad == [], f"no #SBATCH --constraint with an AVX2 feature: {bad}"


def test_avx2_check_flags_a_planted_script():
    head = "#!/bin/bash -l\n#SBATCH -p epyc\n"
    plain = head + 'module load "${M:-MMseqs2/17-b804f}"\nmmseqs easy-cluster in out tmp\n'
    assert needs_avx2(plain) and not has_avx2_constraint(plain)
    assert has_avx2_constraint(plain.replace(head, head + "#SBATCH --constraint=ryzen\n"))
    assert not has_avx2_constraint(plain.replace(head, head + "#SBATCH --constraint=intel\n"))
    assert uses_non_avx2_binary(f"STEP1_MMSEQS={NON_AVX2_MMSEQS}\n")
    assert not needs_avx2("#!/bin/bash -l\n#SBATCH -p exfab\nsignalp6 --help\n")


def _env(tmp_path, fx, **extra):
    import phasec_fixture as pf

    (tmp_path / "scratch").mkdir(exist_ok=True)
    return {
        **os.environ,
        "SCRATCH": str(tmp_path / "scratch"),
        "PROJ_ROOT": str(paths.STEP1_DIR.parents[1]),
        "STEP1_WORKDIR": str(fx["work"]),
        "STEP1_MMSEQS": str(pf.STUB_MMSEQS),
        "STEP1_ENV_PY": os.environ.get("STEP1_ENV_PY", __import__("sys").executable),
        "SLURM_CPUS_PER_TASK": "2",
        "PHASEC_SPECIES": str(fx["species"]),
        "PHASEC_SETS": str(fx["sets"]),
        "PHASEC_N_RESAMPLES": "20",
        "PHASEC_CANDIDATES": "B1,R2,M8",
        **extra,
    }


def test_c1_smoke_run_with_the_stub(tmp_path):
    pytest.importorskip("sklearn")
    import phasec_fixture as pf
    from conftest import load_phasec

    fx = pf.make_work(tmp_path)
    assert load_phasec("08_build_eval_tables").main(pf.build_argv(fx)) == 0
    run = subprocess.run(["bash", str(PHASEC / "c1_evaluate.sh")], env=_env(tmp_path, fx),
                         capture_output=True, text=True)  # fmt: skip
    assert run.returncode == 0, run.stdout[-2000:] + run.stderr[-2000:]
    for step in ("09", "10", "11"):
        assert f"C1 step {step} wall_seconds=" in run.stdout
    out = fx["work"] / "phasec"
    m = json.loads((out / "metrics.json").read_text())
    assert m["settings"]["n_resamples"] == 20 and m["settings"]["candidates"] == ["B1", "R2", "M8"]
    assert list((out / "logs").glob("wall.*.txt"))


def test_c1_stops_when_mmseqs_fails(tmp_path):
    pytest.importorskip("sklearn")
    import phasec_fixture as pf
    from conftest import load_phasec

    fx = pf.make_work(tmp_path)
    assert load_phasec("08_build_eval_tables").main(pf.build_argv(fx)) == 0
    env = _env(tmp_path, fx, STUB_MMSEQS_FAIL="easy-cluster")
    run = subprocess.run(["bash", str(PHASEC / "c1_evaluate.sh")], env=env,
                         capture_output=True, text=True)  # fmt: skip
    assert run.returncode == 2 and "STOP: mmseqs easy-cluster exited with 1" in run.stderr
    assert "C1 step 09 wall_seconds" not in run.stdout
    assert not (fx["work"] / "phasec" / "split_members.tsv.gz").exists()


def test_c1_keeps_the_wall_times_after_a_failed_step(tmp_path):
    # review M-1: step 09 passes, step 10 stops; wall.<job>.txt must still hold the 09 line
    pytest.importorskip("sklearn")
    import phasec_fixture as pf
    from conftest import load_phasec

    fx = pf.make_work(tmp_path)
    assert load_phasec("08_build_eval_tables").main(pf.build_argv(fx)) == 0
    env = _env(tmp_path, fx, SLURM_JOB_ID="m1test", PHASEC_CANDIDATES="M8,XGB")
    run = subprocess.run(["bash", str(PHASEC / "c1_evaluate.sh")], env=env,
                         capture_output=True, text=True)  # fmt: skip
    assert run.returncode == 2 and "unknown candidates ['XGB']" in run.stderr
    assert "C1 step 09 wall_seconds=" in run.stdout and "C1 done" not in run.stdout
    wall = (fx["work"] / "phasec" / "logs" / "wall.m1test.txt").read_text()
    assert wall.startswith("C1 step 09 wall_seconds=") and "step 10" not in wall
