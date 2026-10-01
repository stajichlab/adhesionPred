"""Assembly of the J2 chunks into matrices (CPU, ESM-2 8M)."""

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("torch")
pytest.importorskip("esm")

import assemble_embeddings  # noqa: E402
import chunk_plan  # noqa: E402
import embed_store  # noqa: E402
import truth_table  # noqa: E402
from phaseb_fixture import MODEL, change_row0_sequence, make_plan, run_cpu  # noqa: E402


def test_assemble_builds_row_ordered_matrices(tmp_path):
    work = tmp_path / "w"
    out, unique = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    assert assemble_embeddings.main(["--work-dir", str(work), "--models", MODEL]) == 0
    nterm = np.load(out / "emb" / f"{MODEL}.nterm.npy")
    cterm = np.load(out / "emb" / f"{MODEL}.cterm.npy")
    assert nterm.shape == (len(unique), 320) and nterm.dtype == np.float32
    assert cterm.shape == (1, 320)
    for c in chunk_plan.read_plan(out):
        arr, _ = embed_store.load_chunk(out / "emb", MODEL, c.chunk_id)
        for pos, row in enumerate(c.rows):
            target = nterm[row] if c.window == "nterm" else cterm[int(unique[row]["cterm_row"])]
            np.testing.assert_array_equal(target, arr[pos])
    manifest = truth_table.read_tsv(out / "emb" / "chunk_manifest.tsv")
    assert [m["chunk_id"] for m in manifest] == [c.chunk_id for c in chunk_plan.read_plan(out)]
    for m in manifest:
        meta = embed_store.load_chunk(out / "emb", MODEL, m["chunk_id"])[1]
        assert m["array_sha256"] == meta["array_sha256"]
    full_c = assemble_embeddings.window_matrix(nterm, cterm, [u["cterm_row"] for u in unique])
    long_row = next(int(u["row"]) for u in unique if u["cterm_row"] != "")
    np.testing.assert_array_equal(full_c[long_row], cterm[0])
    short_rows = [int(u["row"]) for u in unique if u["cterm_row"] == ""]
    np.testing.assert_array_equal(full_c[short_rows], nterm[short_rows])


def test_assemble_stops_on_a_missing_chunk(tmp_path, capsys):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    chunk_id = chunk_plan.read_plan(out)[-1].chunk_id
    embed_store.chunk_paths(out / "emb", MODEL, chunk_id)[1].unlink()
    assert assemble_embeddings.main(["--work-dir", str(work), "--models", MODEL]) == 2
    assert chunk_id in capsys.readouterr().err
    assert not (out / "emb" / f"{MODEL}.nterm.npy").exists()


def test_assemble_stops_on_a_stale_plan(tmp_path, capsys):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    change_row0_sequence(out)
    assert assemble_embeddings.main(["--work-dir", str(work), "--models", MODEL]) == 2
    assert "re-run 06_plan_embedding.py" in capsys.readouterr().err
    assert not (out / "emb" / f"{MODEL}.nterm.npy").exists()
    assert not (out / "emb" / "embedding_run.json").exists()
