#!/usr/bin/env python
"""J2: embed the chunks of one job with ESM-2 (layer 6, residue-mean pooling). Resumable.

Run with the adhesionPred conda env and
PYTHONPATH=$PROJ_ROOT/src:$PROJ_ROOT/analysis/step1_compare:$PROJ_ROOT/analysis/step1_compare/jobs.
Pooling is done by
surface_glyco.embeddings.get_esm_embeddings; this script only cuts windows, groups chunks,
checks the result and stores it. Reads $STEP1_WORKDIR/phaseb/unique_sequences.tsv.gz,
chunk_plan.tsv, chunk_members.tsv.gz and job_plan.json (batch size per model).

For each model and each chunk of this job (chunk k goes to job k mod --job-count):
- skip the chunk if phaseb/emb/<model>/<chunk_id>.npy and .json exist, the JSON names the
  same members_sha256, and the array SHA-256 matches the JSON;
- else embed it, write it under --scratch-dir, and copy it to phaseb/emb/<model>/.

STOP (exit 2): the plan files disagree; the plan was made from another
unique_sequences.tsv.gz or a row holds another sequence than planned (re-run 06); --job-count
differs from n_jobs in job_plan.json; --device cuda without a GPU; a sequence that
get_esm_embeddings skipped; NaN or inf in a result. Exit 3: --stop-after-chunks was reached
(test hook that simulates a killed job).
"""

import argparse
import json
import sys
import time
from pathlib import Path

import chunk_plan
import embed_store
import numpy as np
import paths
import seqwindow
import truth_table

MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")
REPR_LAYER = 6


class EmbedError(RuntimeError):
    """A chunk cannot be embedded completely and correctly."""


def embed_window_sequences(seqs: list[str], model: str, batch_size: int, device, layer: int):
    """float32 array (len(seqs), dim) from get_esm_embeddings, in input order."""
    from surface_glyco.embeddings import get_esm_embeddings

    records = [{"id": str(i), "sequence": s} for i, s in enumerate(seqs)]
    arr, _, kept = get_esm_embeddings(
        records,
        model_name=model,
        batch_size=batch_size,
        device=device,
        repr_layer=layer,
        return_indices=True,
    )
    if list(kept) != list(range(len(seqs))):
        missing = sorted(set(range(len(seqs))) - set(kept))
        raise EmbedError(f"{len(missing)} sequences were not embedded (positions {missing[:5]})")
    return np.asarray(arr, dtype=np.float32)


def run_job(
    work: Path,
    models,
    job_index: int,
    job_count: int,
    device,
    scratch: Path,
    batch_size: int | None = None,
    stop_after: int | None = None,
) -> dict:
    out = Path(work) / "phaseb"
    publish = out / "emb"
    chunks = {c.chunk_id: c for c in chunk_plan.read_plan(out)}
    job_plan = json.loads((out / "job_plan.json").read_text())
    if job_count != job_plan["n_jobs"]:
        raise EmbedError(
            f"--job-count {job_count} differs from n_jobs {job_plan['n_jobs']} in job_plan.json; "
            "submit J2 with J2_JOB_COUNT=n_jobs"
        )
    mine = chunk_plan.chunks_for_job(list(chunks), job_index, job_count)
    seq_by_row = {
        int(r["row"]): r["sequence"] for r in truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    }
    chunk_plan.check_plan_is_current(out, list(chunks.values()), seq_by_row)
    done = skipped = 0
    for model in models:
        size = batch_size or int(job_plan["batch_size"][model])
        for chunk_id in mine:
            chunk = chunks[chunk_id]
            if embed_store.chunk_is_done(
                publish, model, chunk_id, chunk.members_sha256, len(chunk.rows)
            ):
                skipped += 1
                print(f"skip {model}/{chunk_id} (done, hash verified)", flush=True)
                continue
            seqs = [seqwindow.window(seq_by_row[r], chunk.window) for r in chunk.rows]
            start = time.time()
            arr = embed_window_sequences(seqs, model, size, device, REPR_LAYER)
            meta = {
                "window": chunk.window,
                "members_sha256": chunk.members_sha256,
                "repr_layer": REPR_LAYER,
                "batch_size": size,
                "device": str(device),
                "seconds": round(time.time() - start, 2),
            }
            embed_store.save_chunk(publish, scratch, model, chunk_id, arr, meta)
            done += 1
            print(f"done {model}/{chunk_id} n={len(seqs)} s={meta['seconds']}", flush=True)
            if stop_after is not None and done >= stop_after:
                return {"done": done, "skipped": skipped, "stopped": True}
    return {"done": done, "skipped": skipped, "stopped": False}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    parser.add_argument("--job-index", type=int, default=0)
    parser.add_argument("--job-count", type=int, default=1)
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    parser.add_argument("--scratch-dir", required=True, help="node-local, e.g. $SCRATCH/emb")
    parser.add_argument("--batch-size", type=int, default=None, help="default: job_plan.json")
    parser.add_argument("--stop-after-chunks", type=int, default=None, help="test hook")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        import torch

        if args.device == "cuda" and not torch.cuda.is_available():
            raise EmbedError("--device cuda but torch.cuda.is_available() is False")
        result = run_job(
            work,
            args.models,
            args.job_index,
            args.job_count,
            torch.device(args.device),
            Path(args.scratch_dir),
            args.batch_size,
            args.stop_after_chunks,
        )
    except (EmbedError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result), flush=True)
    return 3 if result["stopped"] else 0


if __name__ == "__main__":
    sys.exit(main())
