"""COLUMNS.md and README.md document the Phase B outputs exactly as the code writes them.

Table sections: the SET of names in the first cell of each table row must equal the code
tuple (rows that say "As in X" still list their names, so the set is complete; a stale name
fails too). Key sections: the bullet list after a marker such as `**Keys:**` must equal the set
of keys in a real output. The text after the first back-ticked name of a bullet is free (it
holds value tokens such as a schema string) and is never read as a key.
"""

import json
import re

import chunk_plan
import paths
import pytest
import seqsets
from conftest import load_script
from test_chunk_plan import _j0, _plan_inputs, _write_j0
from test_phaseb_features import _work
from test_phaseb_prepare import _phaseb_inputs

COLUMNS = (paths.STEP1_DIR / "COLUMNS.md").read_text()
README = (paths.STEP1_DIR / "README.md").read_text()
JOBS = paths.STEP1_DIR / "jobs"
KEYS = "**Keys:**"


def section(heading: str, text: str = COLUMNS) -> str:
    """Text of the COLUMNS.md section whose heading starts with `heading`."""
    parts = re.split(r"^## ", text, flags=re.M)
    hit = [p for p in parts[1:] if p.split()[0] == heading]
    assert len(hit) == 1, f"COLUMNS.md needs exactly one heading {heading}"
    return hit[0]


def table_names(heading: str, text: str = COLUMNS) -> set[str]:
    """Column names in the first cell of every table row of the section."""
    names: list[str] = []
    rows = [ln for ln in section(heading, text).splitlines() if ln.startswith("|")]
    assert len(rows) > 2, f"section {heading} has no table"
    for line in rows[2:]:  # skip the header row and the |---| row
        cell = line.strip("|").split("|")[0]
        names += [n.strip().strip("`").strip() for n in cell.split(",") if n.strip()]
    assert len(names) == len(set(names)), f"section {heading} names a column twice"
    return set(names)


def key_names(heading: str, marker: str = KEYS, text: str = COLUMNS) -> set[str]:
    """Names in the bullet list that follows `marker` in the section."""
    lines = section(heading, text).splitlines()
    assert marker in lines, f"section {heading} has no line {marker!r}"
    names: list[str] = []
    for line in lines[lines.index(marker) + 1 :]:
        if not line.strip():
            if names:
                break
            continue
        match = re.match(r"- `([^`]+)`", line)
        assert match, f"section {heading}: {line!r} is not a bullet that starts with a key"
        names.append(match.group(1))
    assert names and len(names) == len(set(names)), f"section {heading}: empty or duplicate list"
    return set(names)


def test_phaseb_outputs_have_columns_md_headings():
    names = []
    for script in ("05_prepare_sequences", "06_plan_embedding", "07_build_features"):
        names += load_script(script).OUTPUT_NAMES
    names += [
        "j0/throughput.json",
        "j0/gpu_cpu_diff.json",
        "j0/nvidia_smi.csv",
        "signalp/part_NNN/",
        "predgpi/part_NNN.tsv.gz",
        "emb/<model>/<chunk_id>.npy",
        "emb/<model>.nterm.npy",
        "emb/chunk_manifest.tsv",
        "emb/embedding_run.json",
    ]
    headings = [line[3:].split()[0] for line in COLUMNS.splitlines() if line.startswith("## ")]
    missing = [n for n in names if n not in headings]
    assert not missing, f"COLUMNS.md has no heading for {missing}"


def test_table_columns_equal_the_code_constants():
    import predgpi_scores

    f07 = load_script("07_build_features")
    expected = {
        "sequence_members.tsv.gz": seqsets.MEMBER_COLUMNS,
        "unique_sequences.tsv.gz": seqsets.UNIQUE_COLUMNS,
        "chunk_plan.tsv": chunk_plan.PLAN_COLUMNS,
        "chunk_members.tsv.gz": chunk_plan.MEMBER_COLUMNS,
        "predgpi/part_NNN.tsv.gz": predgpi_scores.COLUMNS,
        "features_unique.tsv.gz": f07.UNIQUE_FEATURE_COLUMNS,
        "features.tsv.gz": f07.MEMBER_FEATURE_COLUMNS,
        "feature_coverage.tsv": f07.COVERAGE_COLUMNS,
    }
    for heading, columns in expected.items():
        assert table_names(heading) == set(columns), heading


def test_chunk_manifest_columns_equal_the_code_constant():
    pytest.importorskip("numpy")
    import assemble_embeddings

    got = table_names("emb/chunk_manifest.tsv")
    assert got == set(assemble_embeddings.MANIFEST_COLUMNS)


def test_features_run_keys_equal_the_code(tmp_path, fixtures_dir):
    f07 = load_script("07_build_features")
    work, out = _work(tmp_path, fixtures_dir)
    assert f07.main(["--work-dir", str(work)]) == 0
    log = json.loads((out / "features_run.json").read_text())
    assert key_names("features_run.json") == set(log)
    assert key_names("features_run.json", "**Keys of `input_sha256`:**") == set(log["input_sha256"])


def test_prepare_run_keys_equal_the_code(tmp_path, monkeypatch):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 0
    log = json.loads((work / "phaseb" / "prepare_run.json").read_text())
    assert key_names("prepare_run.json") == set(log)
    one = next(iter(log["inputs"].values()))
    assert key_names("prepare_run.json", "**Keys of each `inputs` object:**") == set(one)


def test_job_plan_keys_equal_the_code_for_both_rate_sources(tmp_path):
    plan06 = load_script("06_plan_embedding")
    out, _ = _plan_inputs(tmp_path / "assumed", [50, 400, 1100, 1300, 900])
    assert plan06.main(["--work-dir", str(out.parent), "--rate", "2000"]) == 0
    assumed = json.loads((out / "job_plan.json").read_text())
    assert assumed["rate_source"] == "assumed"
    base = key_names("job_plan.json")
    assert base == set(assumed)

    out, _ = _plan_inputs(tmp_path / "j0", [50, 400, 1100, 1300, 900])
    _write_j0(out, _j0())
    assert plan06.main(["--work-dir", str(out.parent), "--chunk-residues", "1000"]) == 0
    from_j0 = json.loads((out / "job_plan.json").read_text())
    assert from_j0["rate_source"] == "J0"
    added = key_names("job_plan.json", "**Keys added when `rate_source` is `J0`:**")
    assert base | added == set(from_j0)
    assert not base & added


def test_plan_jobs_keys_are_in_the_documented_job_plan():
    plan = chunk_plan.plan_jobs({"m": 1000.0}, {"m": 1.0}, 5000, [5000])
    assert set(plan) <= key_names("job_plan.json")


# --- parser self-tests: a stale or missing name must make the comparison fail ---


def test_table_parser_sees_a_missing_row_and_a_stale_name():
    text = COLUMNS
    row = next(ln for ln in text.splitlines() if ln.startswith("| residues |"))
    assert "residues" not in table_names("chunk_plan.tsv", text.replace(row + "\n", ""))
    stale = text.replace(row, row + "\n| old_column | Gone. |")
    assert table_names("chunk_plan.tsv", stale) != set(chunk_plan.PLAN_COLUMNS)
    assert table_names("chunk_plan.tsv") == set(chunk_plan.PLAN_COLUMNS)


def test_key_parser_sees_a_missing_key_and_an_extra_key():
    line = next(ln for ln in COLUMNS.splitlines() if ln.startswith("- `longest_job_seconds`"))
    got = key_names("job_plan.json")
    assert "longest_job_seconds" in got
    assert "longest_job_seconds" not in key_names(
        "job_plan.json", text=COLUMNS.replace(line + "\n", "")
    )
    extra = COLUMNS.replace(line, line + "\n- `not_written_by_the_code`")
    assert "not_written_by_the_code" in key_names("job_plan.json", text=extra)


def test_readme_job_rules_match_the_code():
    plan06 = load_script("06_plan_embedding")
    assert f"{plan06.MIN_RATE:,.0f} residues/s" in README
    assert "J1_PARTS" in README and "delete\n  `phaseb/signalp` and `phaseb/predgpi`" in README
    assert "J1_PARTS" in (JOBS / "j1_features.sh").read_text()
    assert "J2_JOB_COUNT" in README and "J2_JOB_COUNT" in (JOBS / "j2_embed.sh").read_text()
    assert "--array=0-<n_jobs-1>" in README and "--time=<time_minutes>" in README
    assert "--device cuda" in (JOBS / "j0_pilot.sh").read_text()
    assert "sequence_run.json" in README and "d8_run.json" in README
    assert "07 does not need J2 or assemble." in README


def test_readme_has_sbatch_examples_and_defines_env_py():
    assert "ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python" in README
    assert "PY=/usr/bin/python3.12" in README
    for job in ("j0_pilot.sh", "j1_features.sh"):
        assert f'"$S1/jobs/{job}"' in README
    assert README.count('--export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR"') >= 2
    assert "export PROJ_ROOT=" in README and "export STEP1_WORKDIR=" in README


def test_nvidia_smi_csv_is_documented_as_the_script_writes_it():
    text = (JOBS / "j0_pilot.sh").read_text()
    query = "nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv"
    assert query in text and "nvidia_smi.csv" in text
    assert f"`{query}`" in section("j0/nvidia_smi.csv")
    assert "phaseb/j0/nvidia_smi.csv" in README


@pytest.mark.parametrize("script", sorted(p.name for p in JOBS.glob("*.sh")))
def test_readme_names_every_job_script(script):
    assert f"jobs/{script}" in README
