"""J0 pilot and the GPU-CPU difference harness, run on CPU with 20 sequences (ESM-2 8M)."""

import json

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("esm")

import gpu_cpu_diff  # noqa: E402
import seqhash  # noqa: E402
import throughput_pilot  # noqa: E402
import truth_table  # noqa: E402

MODEL = "esm2_t6_8M_UR50D"
COLUMNS = ("source_id", "gene_id", "label", "fasta_id", "length", "seq_sha256", "sequence")


def _truth(tmp_path, n=30):
    rows = []
    for i in range(n):
        seq = ("MKTLLVAGLLSSAAFA" + "ACDEFGHIKLMNPQRSTVWY"[i % 20] * (5 + 7 * i))[: 40 + 9 * i]
        if i == 0:
            seq = "MKLSTA" + "STPSSTSA" * 140  # one protein longer than 1,022 aa
        rows.append(
            {
                "source_id": "Scer_SGD",
                "gene_id": f"G{i}",
                "label": "N-int",
                "fasta_id": f"G{i}",
                "length": str(len(seq)),
                "seq_sha256": seqhash.seq_sha256(seq),
                "sequence": seq,
            }
        )
    rows.append({**rows[1], "source_id": "Calb_CGD", "gene_id": "C1"})  # duplicate sequence
    path = tmp_path / "truth_sequences.tsv.gz"
    truth_table.write_tsv(path, COLUMNS, rows)
    return path, rows


def test_sample_is_unique_fixed_by_seed_and_capped(tmp_path):
    _, rows = _truth(tmp_path)
    a = throughput_pilot.sample_sequences(rows, 20, 7)
    b = throughput_pilot.sample_sequences(list(reversed(rows)), 20, 7)
    assert a == b
    assert len({k for k, _ in a}) == 20
    assert max(len(s) for _, s in a) <= 1022
    assert throughput_pilot.sample_sequences(rows, 20, 8) != a
    with pytest.raises(throughput_pilot.PilotError):
        throughput_pilot.sample_sequences(rows, 31, 7)  # 30 unique sequences only


def test_pilot_writes_the_throughput_schema_on_cpu(tmp_path):
    truth_path, _ = _truth(tmp_path)
    out = tmp_path / "phaseb" / "j0" / "throughput.json"
    argv = ["--work-dir", str(tmp_path), "--device", "cpu", "--models", MODEL]
    assert throughput_pilot.main(argv + ["--n", "20", "--batch-sizes", "2,4"]) == 0
    rec = json.loads(out.read_text())
    assert rec["schema"] == throughput_pilot.SCHEMA
    assert rec["device"] == "cpu" and rec["gpu_name"] == "cpu"
    assert rec["n_proteins"] == 20 and rec["seed"] == 20261001
    assert len(rec["sample_sha256"]) == 64 and len(rec["truth_sequences_sha256"]) == 64
    assert [r["batch_size"] for r in rec["runs"]] == [2, 4]
    for r in rec["runs"]:
        assert r["status"] == "ok" and r["dim"] == 320 and r["batch_failures"] == 0
        assert r["proteins_per_s"] > 0 and r["residues_per_s"] > 0
        assert r["peak_mem_bytes"] is None
    assert rec["best"][MODEL]["batch_size"] in (2, 4)
    assert rec["model_load_s"][MODEL] >= 0
    for key in ("torch", "torch_cuda", "esm", "git_commit", "python", "residues"):
        assert key in rec


def test_pilot_output_feeds_the_plan_script(tmp_path):
    from conftest import load_script

    truth_path, rows = _truth(tmp_path)
    argv = ["--work-dir", str(tmp_path), "--device", "cpu", "--models", MODEL]
    assert throughput_pilot.main(argv + ["--n", "20", "--batch-sizes", "4"]) == 0
    import seqsets

    unique = seqsets.unique_rows({r["seq_sha256"]: r["sequence"] for r in rows})
    out = tmp_path / "phaseb"
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    plan = load_script("06_plan_embedding")
    # 06 refuses a CPU pilot (device and rate floor), so the raw CPU record must be rejected.
    assert plan.main(["--work-dir", str(tmp_path), "--models", MODEL]) == 2
    # Stand-in for a GPU run: change only the device and the rate fields; keep the rest as written.
    j0_path = out / "j0" / "throughput.json"
    j0 = json.loads(j0_path.read_text())
    j0["device"] = "cuda"
    j0["best"][MODEL]["residues_per_s"] = 50000.0
    for r in j0["runs"]:
        r["residues_per_s"] = 50000.0
    j0_path.write_text(json.dumps(j0))
    assert plan.main(["--work-dir", str(tmp_path), "--models", MODEL]) == 0
    job = json.loads((out / "job_plan.json").read_text())
    assert job["rate_source"] == "J0" and job["batch_size"][MODEL] == 4


def test_gpu_cpu_diff_schema_with_cpu_on_both_sides(tmp_path):
    _truth(tmp_path)
    out = tmp_path / "diff.json"
    argv = ["--work-dir", str(tmp_path), "--out", str(out), "--models", MODEL]
    argv += ["--device-a", "cpu", "--device-b", "cpu", "--n", "10", "--n-long", "1"]
    assert gpu_cpu_diff.main(argv) == 0
    rec = json.loads(out.read_text())
    assert rec["schema"] == gpu_cpu_diff.SCHEMA and rec["n_windows"] == 11
    m = rec["models"][MODEL]
    assert m["max_abs_diff"] == 0.0 and m["repeat_identical_on_a"] is True
    assert m["min_cosine"] == pytest.approx(1.0, abs=1e-6)


@pytest.mark.skipif(torch.cuda.is_available(), reason="needs a host without a GPU")
def test_gpu_cpu_diff_stops_without_a_gpu(tmp_path, capsys):
    _truth(tmp_path)
    assert gpu_cpu_diff.main(["--work-dir", str(tmp_path), "--models", MODEL]) == 2
    assert "cuda" in capsys.readouterr().err


def test_batch_failures_are_counted_and_never_best(tmp_path, monkeypatch):
    import surface_glyco.embeddings as emb

    real = emb._embed_batch

    def fail_big_batches(model, alphabet, bc, batch, layer, device):
        if len(batch) > 2:
            raise RuntimeError("CUDA out of memory (simulated)")
        return real(model, alphabet, bc, batch, layer, device)

    monkeypatch.setattr(emb, "_embed_batch", fail_big_batches)
    _truth(tmp_path)
    out = tmp_path / "phaseb" / "j0" / "throughput.json"
    argv = ["--work-dir", str(tmp_path), "--device", "cpu", "--models", MODEL]
    assert throughput_pilot.main(argv + ["--n", "8", "--batch-sizes", "2,4"]) == 0
    rec = json.loads(out.read_text())
    by_size = {r["batch_size"]: r for r in rec["runs"]}
    assert by_size[2]["status"] == "ok" and by_size[2]["batch_failures"] == 0
    assert by_size[4]["status"] == "batch_failures" and by_size[4]["batch_failures"] == 2
    assert rec["best"][MODEL]["batch_size"] == 2
