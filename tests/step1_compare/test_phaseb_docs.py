"""COLUMNS.md and README.md document the Phase B outputs as the code writes them."""

import json
import re

import chunk_plan
import paths
import pytest
from conftest import load_script
from test_phaseb_features import _work

COLUMNS = (paths.STEP1_DIR / "COLUMNS.md").read_text()
README = (paths.STEP1_DIR / "README.md").read_text()
JOBS = paths.STEP1_DIR / "jobs"


def _section(heading: str) -> str:
    """Text of the COLUMNS.md section whose heading starts with `heading`."""
    parts = re.split(r"^## ", COLUMNS, flags=re.M)
    hit = [p for p in parts[1:] if p.split()[0] == heading]
    assert len(hit) == 1, f"COLUMNS.md needs exactly one heading {heading}"
    return hit[0]


def _check_names(heading: str, names) -> None:
    text = _section(heading)
    missing = [n for n in names if f"`{n}`" not in text and not re.search(rf"\b{n}\b", text)]
    assert not missing, f"section {heading} does not name {missing}"


def test_phaseb_outputs_have_columns_md_headings():
    names = []
    for script in ("05_prepare_sequences", "06_plan_embedding", "07_build_features"):
        names += load_script(script).OUTPUT_NAMES
    names += [
        "j0/throughput.json",
        "j0/gpu_cpu_diff.json",
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


def test_header_constants_are_documented():
    pytest.importorskip("numpy")
    import assemble_embeddings
    import predgpi_scores
    import seqsets

    f07 = load_script("07_build_features")
    _check_names("sequence_members.tsv.gz", seqsets.MEMBER_COLUMNS)
    _check_names("unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS)
    _check_names("chunk_plan.tsv", chunk_plan.PLAN_COLUMNS)
    _check_names("chunk_members.tsv.gz", chunk_plan.MEMBER_COLUMNS)
    _check_names("predgpi/part_NNN.tsv.gz", predgpi_scores.COLUMNS)
    _check_names("features_unique.tsv.gz", f07.UNIQUE_FEATURE_COLUMNS)
    _check_names("features.tsv.gz", f07.MEMBER_FEATURE_COLUMNS)
    _check_names("feature_coverage.tsv", f07.COVERAGE_COLUMNS)
    _check_names("emb/chunk_manifest.tsv", assemble_embeddings.MANIFEST_COLUMNS)


def test_json_keys_are_documented(tmp_path, fixtures_dir):
    f07 = load_script("07_build_features")
    work, out = _work(tmp_path, fixtures_dir)
    assert f07.main(["--work-dir", str(work)]) == 0
    log = json.loads((out / "features_run.json").read_text())
    _check_names("features_run.json", log)
    _check_names("features_run.json", log["input_sha256"])
    plan = chunk_plan.plan_jobs({"m": 1000.0}, {"m": 1.0}, 5000, [5000])
    _check_names("job_plan.json", plan)


def test_j0_json_keys_are_documented():
    pytest.importorskip("torch")
    import gpu_cpu_diff
    import throughput_pilot

    for module, heading in (
        (throughput_pilot, "j0/throughput.json"),
        (gpu_cpu_diff, "j0/gpu_cpu_diff.json"),
    ):
        assert f"`{module.SCHEMA}`" in _section(heading)
    best = throughput_pilot.best_runs(
        [
            {
                "model": "m",
                "status": "ok",
                "batch_size": 8,
                "residues_per_s": 5.0,
                "proteins_per_s": 1.0,
                "seconds": 1.0,
            }
        ]
    )
    _check_names("j0/throughput.json", best["m"])


def test_sidecar_keys_are_documented():
    text = (JOBS / "embed_chunks.py").read_text()
    block = text[text.index("meta = {") :].split("}", 1)[0]
    keys = re.findall(r'"([a-z_0-9]+)":', block)
    assert "repr_layer" in keys
    _check_names("emb/<model>/<chunk_id>.npy", keys + ["chunk_id", "model", "shape", "dtype"])


def test_readme_job_rules_match_the_code():
    plan06 = load_script("06_plan_embedding")
    assert f"{plan06.MIN_RATE:,.0f} residues/s" in README
    assert "J1_PARTS" in README and "delete\n  `phaseb/signalp` and `phaseb/predgpi`" in README
    assert "J1_PARTS" in (JOBS / "j1_features.sh").read_text()
    assert "J2_JOB_COUNT" in README and "J2_JOB_COUNT" in (JOBS / "j2_embed.sh").read_text()
    assert "--array=0-<n_jobs-1>" in README and "--time=<time_minutes>" in README
    assert "--device cuda" in (JOBS / "j0_pilot.sh").read_text()
    assert "sequence_run.json" in README and "d8_run.json" in README


@pytest.mark.parametrize("script", sorted(p.name for p in JOBS.glob("*.sh")))
def test_readme_names_every_job_script(script):
    assert f"jobs/{script}" in README
