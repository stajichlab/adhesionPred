#!/usr/bin/env python
"""J0: measure ESM-2 embedding throughput on a fixed sample of truth sequences.

Samples --n unique sequences (by seq_sha256) from truth_sequences.tsv.gz with a fixed seed,
cuts the N-terminal window (first 1,022 residues), and times
surface_glyco.embeddings.get_esm_embeddings for every model and batch size. Writes one JSON
file (schema SCHEMA) with proteins/s, residues/s, peak GPU memory, the GPU model and, per
model, the fastest batch size (`best`). 06_plan_embedding.py reads `best` and `model_load_s`.

get_esm_embeddings catches every batch error (out of memory included), prints
"Warning: batch failed" and retries the batch one sequence at a time, so no exception reaches
this script. The pilot captures that output and records `batch_failures` per run. A run with
batch failures gets status `batch_failures` (its time includes the slow retries) and is never
chosen as `best`. A run that lost sequences gets status `skipped_<n>`. STOP (exit 2): the
input is missing, --device cuda without a GPU, or no model has a successful run.
"""

import argparse
import contextlib
import io
import json
import random
import sys
import time
from pathlib import Path

import chunk_plan
import manifest
import paths
import runinfo
import seqwindow
import truth_table

SCHEMA = "step1-phaseb-j0-throughput/2"
BATCH_FAILED = "Warning: batch failed"  # printed by surface_glyco.embeddings.get_esm_embeddings
MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")


class PilotError(RuntimeError):
    """The pilot cannot produce a usable throughput record."""


def sample_sequences(truth_rows, n: int, seed: int) -> list[tuple[str, str]]:
    """(seq_sha256, N-terminal window) for n unique sequences; same input and seed, same set."""
    unique = {r["seq_sha256"]: r["sequence"] for r in truth_rows}
    keys = sorted(unique)
    if n > len(keys):
        raise PilotError(f"--n {n} is larger than the {len(keys)} unique truth sequences")
    chosen = sorted(random.Random(seed).sample(keys, n))
    return [(k, seqwindow.nterm_window(unique[k])) for k in chosen]


def _sync(torch, device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize()


def time_model(model: str, sample, batch_sizes, device) -> tuple[list[dict], float]:
    import torch

    from surface_glyco import embeddings

    start = time.time()
    embeddings.get_cached_model(model, device)
    _sync(torch, device)
    load_s = time.time() - start
    records = [{"id": k, "sequence": s} for k, s in sample]
    residues = sum(len(s) for _, s in sample)
    embeddings.get_esm_embeddings(records[:4], model, 4, device)  # warm-up, not timed
    runs = []
    for size in batch_sizes:
        if device.type == "cuda":
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats(device)
        run = {"model": model, "batch_size": size}
        captured = io.StringIO()
        _sync(torch, device)
        t0 = time.time()
        with contextlib.redirect_stdout(captured):
            arr, _, kept = embeddings.get_esm_embeddings(
                records, model, size, device, return_indices=True
            )
        _sync(torch, device)
        secs = time.time() - t0
        sys.stdout.write(captured.getvalue())
        run["batch_failures"] = captured.getvalue().count(BATCH_FAILED)
        if len(kept) != len(records):
            runs.append({**run, "status": f"skipped_{len(records) - len(kept)}"})
            continue
        run.update(
            {
                "status": "ok" if run["batch_failures"] == 0 else "batch_failures",
                "seconds": round(secs, 3),
                "proteins_per_s": round(len(records) / secs, 2),
                "residues_per_s": round(residues / secs, 1),
                "peak_mem_bytes": int(torch.cuda.max_memory_allocated(device))
                if device.type == "cuda"
                else None,
                "dim": int(arr.shape[1]),
            }
        )
        runs.append(run)
    return runs, round(load_s, 3)


def best_runs(runs) -> dict:
    best = {}
    for r in runs:
        if r["status"] != "ok":
            continue
        cur = best.get(r["model"])
        if cur is None or r["residues_per_s"] > cur["residues_per_s"]:
            best[r["model"]] = {
                k: r[k] for k in ("batch_size", "residues_per_s", "proteins_per_s", "seconds")
            }
    return best


def run(truth_path: Path, out_path: Path, models, batch_sizes, n, seed, device) -> dict:
    import esm
    import torch

    sample = sample_sequences(truth_table.read_tsv(truth_path), n, seed)
    runs, load_s = [], {}
    for model in models:
        got, load_s[model] = time_model(model, sample, batch_sizes, device)
        runs += got
    best = best_runs(runs)
    if not best:
        raise PilotError("no model has a successful run")
    record = {
        "schema": SCHEMA,
        "device": device.type,
        "gpu_name": torch.cuda.get_device_name(device) if device.type == "cuda" else "cpu",
        "torch": torch.__version__,
        "torch_cuda": torch.version.cuda,
        "esm": getattr(esm, "__version__", "unknown"),
        "n_proteins": len(sample),
        "residues": sum(len(s) for _, s in sample),
        "seed": seed,
        "sample_sha256": chunk_plan.members_digest([k for k, _ in sample]),
        "truth_sequences_sha256": manifest.sha256_file(truth_path),
        "batch_sizes": list(batch_sizes),
        "model_load_s": load_s,
        "runs": runs,
        "best": best,
        "git_commit": runinfo.git_commit(),
        "python": runinfo.python_version(),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writers = {out_path.name: lambda p: runinfo.write_json(p, record)}
    runinfo.atomic_write_all(out_path.parent, writers)
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--out", default=None, help="default: phaseb/j0/throughput.json")
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    parser.add_argument("--batch-sizes", default="8,16,32,64")
    parser.add_argument("--n", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    out = Path(args.out) if args.out else work / "phaseb" / "j0" / "throughput.json"
    try:
        import torch

        if args.device == "cuda" and not torch.cuda.is_available():
            raise PilotError("--device cuda but torch.cuda.is_available() is False")
        sizes = [int(x) for x in args.batch_sizes.split(",")]
        record = run(
            work / "truth_sequences.tsv.gz",
            out,
            args.models,
            sizes,
            args.n,
            args.seed,
            torch.device(args.device),
        )
    except (PilotError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(record["best"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
