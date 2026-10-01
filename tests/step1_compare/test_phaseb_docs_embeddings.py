"""COLUMNS.md keys of the outputs that need torch and fair-esm (CPU runs on tiny inputs)."""

import json

import pytest

pytest.importorskip("torch")
pytest.importorskip("esm")

import assemble_embeddings  # noqa: E402
import gpu_cpu_diff  # noqa: E402
import throughput_pilot  # noqa: E402
from phaseb_fixture import MODEL, make_plan, run_cpu  # noqa: E402
from test_phaseb_docs import key_names  # noqa: E402
from test_phaseb_pilot import _truth  # noqa: E402


def test_embedding_run_and_sidecar_keys_equal_the_code(tmp_path):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    assert assemble_embeddings.main(["--work-dir", str(work), "--models", MODEL]) == 0
    log = json.loads((out / "emb" / "embedding_run.json").read_text())
    assert key_names("emb/embedding_run.json") == set(log)
    window = log["models"][MODEL]["nterm"]
    assert key_names("emb/embedding_run.json", "**Keys of each window object in `models`:**") == (
        set(window)
    )
    sidecar = json.loads(next((out / "emb" / MODEL).glob("*.json")).read_text())
    assert key_names("emb/<model>/<chunk_id>.npy", "**Keys of the sidecar:**") == set(sidecar)


def test_j0_json_keys_equal_the_code(tmp_path):
    _truth(tmp_path)
    argv = ["--work-dir", str(tmp_path), "--device", "cpu", "--models", MODEL]
    assert throughput_pilot.main(argv + ["--n", "20", "--batch-sizes", "2"]) == 0
    rec = json.loads((tmp_path / "phaseb" / "j0" / "throughput.json").read_text())
    assert key_names("j0/throughput.json") == set(rec)
    assert rec["runs"][0]["status"] == "ok"
    assert key_names("j0/throughput.json", "**Keys of each `runs` object with status `ok`:**") == (
        set(rec["runs"][0])
    )
    assert key_names("j0/throughput.json", "**Keys of each `best` object:**") == set(
        rec["best"][MODEL]
    )
    assert f"`{throughput_pilot.SCHEMA}`" in _line("j0/throughput.json", "schema")

    out = tmp_path / "diff.json"
    argv = ["--work-dir", str(tmp_path), "--out", str(out), "--models", MODEL]
    argv += ["--device-a", "cpu", "--device-b", "cpu", "--n", "10", "--n-long", "1"]
    assert gpu_cpu_diff.main(argv) == 0
    diff = json.loads(out.read_text())
    assert key_names("j0/gpu_cpu_diff.json") == set(diff)
    assert key_names("j0/gpu_cpu_diff.json", "**Keys of each `models` object:**") == set(
        diff["models"][MODEL]
    )
    assert f"`{gpu_cpu_diff.SCHEMA}`" in _line("j0/gpu_cpu_diff.json", "schema")


def _line(heading, key):
    from test_phaseb_docs import section

    return next(ln for ln in section(heading).splitlines() if ln.startswith(f"- `{key}`"))
