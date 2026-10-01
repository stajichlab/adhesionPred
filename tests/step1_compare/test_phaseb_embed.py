"""J2 runner on CPU with ESM-2 8M (needs torch, fair-esm and the 8M weights)."""

import os
import subprocess
import sys

import pytest

np = pytest.importorskip("numpy")
torch = pytest.importorskip("torch")
pytest.importorskip("esm")

import chunk_plan  # noqa: E402
import embed_chunks  # noqa: E402
import embed_store  # noqa: E402
from conftest import STEP1_DIR  # noqa: E402
from phaseb_fixture import LONG, MODEL, change_row0_sequence, make_plan, run_cpu  # noqa: E402

CPU = torch.device("cpu")


def test_job_embeds_every_chunk_member_once(tmp_path):
    work = tmp_path / "w"
    out, unique = make_plan(work)
    result = run_cpu(work, tmp_path / "scratch")
    chunks = chunk_plan.read_plan(out)
    assert result == {"done": len(chunks), "skipped": 0, "stopped": False}
    for c in chunks:
        arr, meta = embed_store.load_chunk(out / "emb", MODEL, c.chunk_id)
        assert arr.shape == (len(c.rows), 320) and arr.dtype == np.float32
        assert np.isfinite(arr).all()
        assert meta["members_sha256"] == c.members_sha256 and meta["repr_layer"] == 6


def test_rerun_skips_finished_chunks_after_hash_check(tmp_path):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    first = run_cpu(work, tmp_path / "scratch")
    again = run_cpu(work, tmp_path / "scratch")
    assert again == {"done": 0, "skipped": first["done"], "stopped": False}


def test_corrupted_chunk_is_recomputed_with_the_same_hash(tmp_path):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    chunk_id = chunk_plan.read_plan(out)[0].chunk_id
    npy, js = embed_store.chunk_paths(out / "emb", MODEL, chunk_id)
    good = embed_store.load_chunk(out / "emb", MODEL, chunk_id)[1]["array_sha256"]
    np.save(npy, np.zeros_like(np.load(npy)))
    assert not embed_store.chunk_is_done(out / "emb", MODEL, chunk_id, "x", 1)
    again = run_cpu(work, tmp_path / "scratch")
    assert again["done"] == 1
    # determinism on CPU: the recomputed chunk is bit-identical to the first run
    assert embed_store.load_chunk(out / "emb", MODEL, chunk_id)[1]["array_sha256"] == good


def test_kill_and_resume_gives_the_same_arrays(tmp_path):
    work_a, work_b = tmp_path / "a", tmp_path / "b"
    out_a, _ = make_plan(work_a)
    out_b, _ = make_plan(work_b)
    run_cpu(work_b, tmp_path / "scratch")  # uninterrupted reference
    script = STEP1_DIR / "jobs" / "embed_chunks.py"
    env = {
        **os.environ,
        "PYTHONPATH": os.pathsep.join(
            [str(STEP1_DIR.parents[1] / "src"), str(STEP1_DIR), str(STEP1_DIR / "jobs")]
        ),
    }
    base = [sys.executable, str(script), "--work-dir", str(work_a), "--models", MODEL]
    base += ["--device", "cpu", "--batch-size", "2", "--scratch-dir", str(tmp_path / "s")]
    killed = subprocess.run(base + ["--stop-after-chunks", "1"], env=env, capture_output=True)
    assert killed.returncode == 3, killed.stderr.decode()[-2000:]
    resumed = subprocess.run(base, env=env, capture_output=True, text=True)
    assert resumed.returncode == 0, resumed.stderr[-2000:]
    assert '"skipped": 1' in resumed.stdout
    for c in chunk_plan.read_plan(out_a):
        a = embed_store.load_chunk(out_a / "emb", MODEL, c.chunk_id)[1]["array_sha256"]
        b = embed_store.load_chunk(out_b / "emb", MODEL, c.chunk_id)[1]["array_sha256"]
        assert a == b, c.chunk_id


def test_chunk_rows_do_not_depend_on_batch_size(tmp_path):
    seqs = ["MKTLLVAGLLSSAAFA", "MSTTSSTTSTPSSTSA" * 4, LONG[-1022:]]
    one = embed_chunks.embed_window_sequences(seqs, MODEL, 1, CPU, 6)
    three = embed_chunks.embed_window_sequences(seqs, MODEL, 3, CPU, 6)
    np.testing.assert_allclose(one, three, atol=1e-5)


def test_cterm_chunk_row_is_the_embedding_of_the_last_1022_residues(tmp_path):
    from surface_glyco.embeddings import get_esm_embeddings

    work = tmp_path / "w"
    out, unique = make_plan(work, chunk_residues=5000)
    run_cpu(work, tmp_path / "scratch")
    cterm = [c for c in chunk_plan.read_plan(out) if c.window == "cterm"]
    assert len(cterm) == 1 and len(cterm[0].rows) == 1
    arr, _ = embed_store.load_chunk(out / "emb", MODEL, cterm[0].chunk_id)
    direct, _ = get_esm_embeddings([{"id": "x", "sequence": LONG[-1022:]}], MODEL, 1, CPU)
    np.testing.assert_allclose(arr[0], direct[0], atol=1e-5)
    nterm_direct, _ = get_esm_embeddings([{"id": "x", "sequence": LONG}], MODEL, 1, CPU)
    assert not np.allclose(arr[0], nterm_direct[0], atol=1e-3)


def test_a_skipped_sequence_stops_the_chunk(monkeypatch):
    import surface_glyco.embeddings as emb

    def drop_last(records, **kw):
        n = len(records) - 1
        return np.zeros((n, 320), dtype=np.float32), [r["id"] for r in records[:n]], list(range(n))

    monkeypatch.setattr(emb, "get_esm_embeddings", drop_last)
    with pytest.raises(embed_chunks.EmbedError, match="not embedded"):
        embed_chunks.embed_window_sequences(["MKV", "MKL"], MODEL, 2, CPU, 6)


def test_nan_result_is_refused(tmp_path):
    arr = np.full((2, 4), np.nan, dtype=np.float32)
    with pytest.raises(ValueError, match="NaN or inf"):
        embed_store.save_chunk(tmp_path / "p", tmp_path / "s", MODEL, "nterm_0000", arr, {})
    assert not (tmp_path / "p").exists()


@pytest.mark.skipif(torch.cuda.is_available(), reason="needs a host without a GPU")
def test_cuda_request_without_gpu_stops(tmp_path, capsys):
    work = tmp_path / "w"
    make_plan(work)
    argv = ["--work-dir", str(work), "--device", "cuda", "--scratch-dir", str(tmp_path / "s")]
    assert embed_chunks.main(argv) == 2
    assert "cuda" in capsys.readouterr().err


def _cli(work, tmp_path, *extra):
    argv = ["--work-dir", str(work), "--models", MODEL, "--device", "cpu", "--batch-size", "2"]
    return embed_chunks.main(argv + ["--scratch-dir", str(tmp_path / "s"), *extra])


def test_stale_plan_stops_the_job(tmp_path, capsys):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    change_row0_sequence(out)
    assert _cli(work, tmp_path) == 2
    assert "re-run 06_plan_embedding.py" in capsys.readouterr().err
    assert not (out / "emb").exists()


def test_chunk_json_with_other_members_is_recomputed(tmp_path):
    import json

    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    chunk = chunk_plan.read_plan(out)[0]
    _, js = embed_store.chunk_paths(out / "emb", MODEL, chunk.chunk_id)
    meta = json.loads(js.read_text())
    meta["members_sha256"] = "0" * 64  # array and hash stay valid; only the members differ
    js.write_text(json.dumps(meta))
    assert embed_store.load_chunk(out / "emb", MODEL, chunk.chunk_id)  # passes its hash check
    done = embed_store.chunk_is_done
    assert not done(out / "emb", MODEL, chunk.chunk_id, chunk.members_sha256, len(chunk.rows))
    again = run_cpu(work, tmp_path / "scratch")
    assert again["done"] == 1
    meta = json.loads(js.read_text())
    assert meta["members_sha256"] == chunk.members_sha256


def test_job_count_must_match_the_plan(tmp_path, capsys):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    assert _cli(work, tmp_path, "--job-count", "2", "--job-index", "0") == 2
    assert "J2_JOB_COUNT" in capsys.readouterr().err
    assert not (out / "emb").exists()
