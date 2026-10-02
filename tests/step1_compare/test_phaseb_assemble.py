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


# --- fix round 1 (coordinator rulings) ---------------------------------------------------


def _custom_plan(work, seqs, chunk_residues):
    """Like phaseb_fixture.make_plan, but for a given list of sequences."""
    import seqhash
    import seqsets
    from conftest import load_script

    out = work / "phaseb"
    out.mkdir(parents=True)
    unique = seqsets.unique_rows({seqhash.seq_sha256(s): s for s in seqs})
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    plan = load_script("06_plan_embedding")
    argv = ["--work-dir", str(work), "--rate", "100", "--models", MODEL]
    assert plan.main(argv + ["--chunk-residues", str(chunk_residues)]) == 0
    return out, unique


def _assemble(work, *models):
    return assemble_embeddings.main(["--work-dir", str(work), "--models", *(models or [MODEL])])


def test_stray_tmp_files_are_never_read(tmp_path):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    assert _assemble(work) == 0
    emb = out / "emb"
    base = {w: (emb / f"{MODEL}.{w}.npy").read_bytes() for w in ("nterm", "cterm")}
    stray = emb / MODEL / ".tmp.nterm_0000.npy"
    for payload in (b"not an npy file", None):
        if payload is None:
            with open(stray, "wb") as handle:
                np.save(handle, np.full((3, 320), np.nan, dtype=np.float32))
        else:
            stray.write_bytes(payload)
        assert _assemble(work) == 0
        for w in ("nterm", "cterm"):
            mat = np.load(emb / f"{MODEL}.{w}.npy")
            assert mat.dtype == np.float32 and np.isfinite(mat).all()
            assert (emb / f"{MODEL}.{w}.npy").read_bytes() == base[w]


def test_each_cterm_row_holds_the_embedding_of_its_own_last_1022_residues(tmp_path):
    import seqwindow
    import torch
    from embed_chunks import embed_window_sequences

    longs = [
        "MKLSTA" + "STPSSTSAGN" * 130,
        "MQQ" + "AVLGSTNPEK" * 125,
        "MRR" + "TTSSEDPKQV" * 140,
        "MSN" + "GGSTPAEEKL" * 120,
    ]
    seqs = ["MKTLLVAGLLSSAAFA", longs[0], "MSTTSSTTSTPSSTSA" * 4, *longs[1:]]
    work = tmp_path / "w"
    out, unique = _custom_plan(work, seqs, chunk_residues=2100)
    run_cpu(work, tmp_path / "scratch")
    assert _assemble(work) == 0
    cterm = np.load(out / "emb" / f"{MODEL}.cterm.npy")
    long_rows = [u for u in unique if u["cterm_row"] != ""]
    assert len(long_rows) == 4 and cterm.shape[0] == 4
    # Two cterm chunks: the position inside a chunk differs from cterm_row for the later one.
    assert len([c for c in chunk_plan.read_plan(out) if c.window == "cterm"]) == 2
    for u in long_rows:
        window = seqwindow.window(u["sequence"], "cterm")
        assert len(window) == 1022
        direct = embed_window_sequences([window], MODEL, 1, torch.device("cpu"), 6)[0]
        np.testing.assert_allclose(cterm[int(u["cterm_row"])], direct, atol=1e-5, rtol=0)


def test_missing_chunk_of_the_second_model_stops_with_no_output(tmp_path, capsys, monkeypatch):
    import json
    import shutil

    import embed_constants

    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    second = "esm2_copy"  # assembly only reads files; no weights are needed for a copy
    shutil.copytree(out / "emb" / MODEL, out / "emb" / second)
    monkeypatch.setitem(embed_constants.MODEL_DIM, second, 320)  # the copy is a 320-wide model
    for js in (out / "emb" / second).glob("*.json"):
        js.write_text(json.dumps({**json.loads(js.read_text()), "model": second}))
    chunk_id = chunk_plan.read_plan(out)[-1].chunk_id
    embed_store.chunk_paths(out / "emb", second, chunk_id)[0].unlink()
    assert _assemble(work, MODEL, second) == 2
    err = capsys.readouterr().err
    assert "STOP:" in err and second in err and chunk_id in err
    assert sorted(p.name for p in (out / "emb").iterdir()) == sorted([MODEL, second])


def test_chunk_file_of_another_window_or_chunk_id_stops(tmp_path, capsys):
    import shutil

    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    plan = chunk_plan.read_plan(out)
    cterm = next(c for c in plan if c.window == "cterm")
    twin = next(c for c in plan if c.window == "nterm" and c.hashes == cterm.hashes)
    assert twin.members_sha256 == cterm.members_sha256
    for src, dst in zip(
        embed_store.chunk_paths(out / "emb", MODEL, twin.chunk_id),
        embed_store.chunk_paths(out / "emb", MODEL, cterm.chunk_id),
        strict=True,
    ):
        shutil.copyfile(src, dst)
    assert _assemble(work) == 2
    err = capsys.readouterr().err
    assert cterm.chunk_id in err and "window" in err
    assert not (out / "emb" / f"{MODEL}.nterm.npy").exists()


def test_manifest_has_members_sha256_in_plan_order(tmp_path):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    assert _assemble(work) == 0
    plan = chunk_plan.read_plan(out)
    manifest = truth_table.read_tsv(out / "emb" / "chunk_manifest.tsv")
    assert [(m["chunk_id"], m["window"], m["n_seqs"], m["members_sha256"]) for m in manifest] == [
        (c.chunk_id, c.window, str(len(c.rows)), c.members_sha256) for c in plan
    ]


def test_empty_plan_stops_with_a_clear_message(tmp_path, capsys):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    header = (out / "chunk_plan.tsv").read_text().splitlines()[0]
    (out / "chunk_plan.tsv").write_text(header + "\n")
    truth_table.write_tsv(out / "chunk_members.tsv.gz", chunk_plan.MEMBER_COLUMNS, [])
    assert _assemble(work) == 2
    err = capsys.readouterr().err
    assert err.startswith("STOP:") and "no chunks" in err and "06_plan_embedding.py" in err
    assert not (out / "emb").exists()


@pytest.mark.parametrize(
    ("change", "text"),
    [
        ({"model": "esm2_t12_35M_UR50D"}, "records model"),
        ({"repr_layer": 12}, "repr_layer"),
    ],
)
def test_assembly_stops_on_a_chunk_of_another_model_or_layer(tmp_path, capsys, change, text):
    import json

    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    chunk = chunk_plan.read_plan(out)[0]
    _, js = embed_store.chunk_paths(out / "emb", MODEL, chunk.chunk_id)
    meta = json.loads(js.read_text())
    meta.update(change)
    js.write_text(json.dumps(meta))
    assert _assemble(work) == 2
    assert text in capsys.readouterr().err
    assert not (out / "emb" / f"{MODEL}.nterm.npy").exists()


def test_assembly_stops_on_a_chunk_with_the_wrong_width(tmp_path, capsys):
    work = tmp_path / "w"
    out, _ = make_plan(work)
    run_cpu(work, tmp_path / "scratch")
    chunk = chunk_plan.read_plan(out)[0]
    npy, js = embed_store.chunk_paths(out / "emb", MODEL, chunk.chunk_id)
    arr = np.load(npy)[:, :100]
    meta = embed_store.load_chunk(out / "emb", MODEL, chunk.chunk_id)[1]
    embed_store.save_chunk(out / "emb", tmp_path / "s2", MODEL, chunk.chunk_id, arr, meta)
    assert _assemble(work) == 2
    assert "width" in capsys.readouterr().err
    assert not (out / "emb" / f"{MODEL}.nterm.npy").exists()
