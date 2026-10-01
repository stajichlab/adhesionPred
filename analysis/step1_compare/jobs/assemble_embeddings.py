#!/usr/bin/env python
"""After J2: check every chunk and build one matrix per model and window. Needs numpy.

Reads phaseb/unique_sequences.tsv.gz, the chunk plan and phaseb/emb/<model>/<chunk_id>.*.
Writes to phaseb/emb/:

  <model>.nterm.npy        float32, shape (unique sequences, dim); row = `row` column
  <model>.cterm.npy        float32, shape (sequences > 1,022 aa, dim); row = `cterm_row`
  chunk_manifest.tsv       one row per model and chunk: members and array SHA-256
  embedding_run.json       shapes, matrix SHA-256 values, chunk count, plan inputs

A sequence of 1,022 aa or less has no C-terminal row: its C-terminal window is the whole
sequence, so the M8-C and M35-C candidates use its `nterm` row (see `window_matrix`).
STOP (exit 2, no output): the plan was made from another unique_sequences.tsv.gz or a row
holds another sequence than planned (re-run 06); a chunk is missing, fails its hash, or holds other members than the
plan; a matrix row is filled twice or not at all; NaN or inf; a row count that is not the
unique sequence count.
"""

import argparse
import sys
from pathlib import Path

import chunk_plan
import embed_store
import manifest
import numpy as np
import paths
import runinfo
import truth_table

MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")
MANIFEST_COLUMNS = ("model", "chunk_id", "window", "n_seqs", "members_sha256", "array_sha256")


class AssembleError(ValueError):
    """The chunks do not form complete matrices."""


def window_matrix(nterm: np.ndarray, cterm: np.ndarray, cterm_rows: list[str]) -> np.ndarray:
    """Full C-terminal matrix (one row per unique sequence) for the M8-C / M35-C candidates."""
    out = nterm.copy()
    for row, crow in enumerate(cterm_rows):
        if crow != "":
            out[row] = cterm[int(crow)]
    return out


def assemble_model(publish: Path, model: str, chunks, unique, manifest_rows: list) -> dict:
    n_rows = {"nterm": len(unique), "cterm": sum(r["cterm_row"] != "" for r in unique)}
    target = {int(r["row"]): int(r["cterm_row"]) if r["cterm_row"] != "" else None for r in unique}
    mats: dict[str, np.ndarray] = {}
    filled = {w: np.zeros(n, dtype=bool) for w, n in n_rows.items()}
    for chunk in chunks:
        try:
            arr, meta = embed_store.load_chunk(publish, model, chunk.chunk_id)
        except (OSError, ValueError, KeyError) as exc:
            raise AssembleError(f"{model}/{chunk.chunk_id}: {exc}") from exc
        if meta.get("members_sha256") != chunk.members_sha256 or arr.shape[0] != len(chunk.rows):
            raise AssembleError(f"{model}/{chunk.chunk_id}: members differ from the plan")
        manifest_rows.append(
            {
                "model": model,
                "chunk_id": chunk.chunk_id,
                "window": chunk.window,
                "n_seqs": str(len(chunk.rows)),
                "members_sha256": chunk.members_sha256,
                "array_sha256": meta["array_sha256"],
            }
        )
        w = chunk.window
        if w not in mats:
            mats[w] = np.zeros((n_rows[w], arr.shape[1]), dtype=embed_store.DTYPE)
        idx = [r if w == "nterm" else target[r] for r in chunk.rows]
        if any(i is None for i in idx) or filled[w][idx].any():
            raise AssembleError(f"{model}/{chunk.chunk_id}: a row is filled twice or has no slot")
        mats[w][idx] = arr
        filled[w][idx] = True
    for w, n in n_rows.items():
        if n and (w not in mats or not filled[w].all()):
            raise AssembleError(f"{model}: {int((~filled[w]).sum())} {w} rows have no embedding")
        if w not in mats:
            mats[w] = np.zeros((0, mats["nterm"].shape[1]), dtype=embed_store.DTYPE)
        if not np.isfinite(mats[w]).all():
            raise AssembleError(f"{model}: {w} matrix has NaN or inf")
    return mats


def run(work: Path, models) -> dict:
    out = Path(work) / "phaseb"
    publish = out / "emb"
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    chunks = chunk_plan.read_plan(out)
    seq_by_row = {int(r["row"]): r["sequence"] for r in unique}
    chunk_plan.check_plan_is_current(out, chunks, seq_by_row)
    log = {"models": {}, "unique_sequences": len(unique), "chunks": len(chunks)}
    writers, manifest_rows = {}, []
    for model in models:
        mats = assemble_model(publish, model, chunks, unique, manifest_rows)
        log["models"][model] = {
            w: {
                "shape": list(m.shape),
                "dtype": str(m.dtype),
                "array_sha256": embed_store.array_sha256(m),
            }
            for w, m in mats.items()
        }
        for w, m in mats.items():
            writers[f"{model}.{w}.npy"] = lambda p, m=m: _save(p, m)
    log.update(
        {
            "unique_sequences_sha256": _sha(out / "unique_sequences.tsv.gz"),
            "chunk_plan_sha256": _sha(out / "chunk_plan.tsv"),
            "git_commit": runinfo.git_commit(),
            "python": runinfo.python_version(),
        }
    )
    writers["chunk_manifest.tsv"] = lambda p: truth_table.write_tsv(
        p, MANIFEST_COLUMNS, manifest_rows
    )
    writers["embedding_run.json"] = lambda p: runinfo.write_json(p, log)
    runinfo.atomic_write_all(publish, writers)
    return log


def _save(path: Path, arr: np.ndarray) -> None:
    with open(path, "wb") as handle:  # np.save on a handle keeps the temp name unchanged
        np.save(handle, arr)


def _sha(path: Path) -> str:
    return manifest.sha256_file(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        log = run(work, args.models)
    except (AssembleError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for model, mats in log["models"].items():
        print(model, {w: m["shape"] for w, m in mats.items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
