#!/usr/bin/env python
"""Spec 7 harness: difference between ESM-2 embeddings computed on two devices.

Samples --n unique truth sequences (fixed seed, the J0 sampler), adds the C-terminal window of
the first --n-long sequences longer than 1,022 aa, and embeds all windows with
surface_glyco.embeddings.get_esm_embeddings on --device-a and --device-b (default cuda and
cpu). It also embeds the windows twice on --device-a to check run-to-run repeatability.
Writes one JSON file per call: per model the maximum and mean absolute difference, the largest
difference relative to the largest absolute value, the lowest cosine similarity, and whether
the repeat on device A is bit-identical. No threshold is applied: the numbers are reported.

STOP (exit 2): a cuda device is requested and none is available; a sequence was not embedded.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import chunk_plan
import numpy as np
import paths
import runinfo
import seqwindow
import throughput_pilot
import truth_table

SCHEMA = "step1-phaseb-gpu-cpu-diff/1"
MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")


def windows_for(truth_rows, n: int, n_long: int, seed: int) -> list[str]:
    sample = throughput_pilot.sample_sequences(truth_rows, n, seed)
    unique = {r["seq_sha256"]: r["sequence"] for r in truth_rows}
    long_keys = sorted(k for k, s in unique.items() if seqwindow.needs_cterm(s))[:n_long]
    return [s for _, s in sample] + [seqwindow.cterm_window(unique[k]) for k in long_keys]


def embed(seqs, model, device, batch_size) -> np.ndarray:
    from surface_glyco.embeddings import get_esm_embeddings

    records = [{"id": str(i), "sequence": s} for i, s in enumerate(seqs)]
    arr, _, kept = get_esm_embeddings(records, model, batch_size, device, return_indices=True)
    if len(kept) != len(seqs):
        raise RuntimeError(f"{model} on {device}: {len(seqs) - len(kept)} sequences skipped")
    return np.asarray(arr, dtype=np.float32)


def compare(a: np.ndarray, b: np.ndarray) -> dict:
    diff = np.abs(a.astype(np.float64) - b.astype(np.float64))
    na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    cos = (a * b).sum(axis=1) / np.maximum(na * nb, 1e-12)
    return {
        "max_abs_diff": float(diff.max()),
        "mean_abs_diff": float(diff.mean()),
        "max_rel_diff": float(diff.max() / max(float(np.abs(b).max()), 1e-12)),
        "min_cosine": float(cos.min()),
    }


def run(truth_path, out_path, models, dev_a, dev_b, n, n_long, seed, batch_size) -> dict:
    import torch

    seqs = windows_for(truth_table.read_tsv(truth_path), n, n_long, seed)
    result = {
        "schema": SCHEMA,
        "device_a": str(dev_a),
        "device_b": str(dev_b),
        "gpu_name": torch.cuda.get_device_name(dev_a) if dev_a.type == "cuda" else "cpu",
        "torch": torch.__version__,
        "n_windows": len(seqs),
        "windows_sha256": chunk_plan.members_digest(
            [hashlib.sha256(s.encode()).hexdigest() for s in seqs]
        ),
        "batch_size": batch_size,
        "models": {},
        "git_commit": runinfo.git_commit(),
    }
    for model in models:
        a1 = embed(seqs, model, dev_a, batch_size)
        a2 = embed(seqs, model, dev_a, batch_size)
        b = embed(seqs, model, dev_b, batch_size)
        result["models"][model] = {
            **compare(a1, b),
            "repeat_identical_on_a": bool(np.array_equal(a1, a2)),
            "repeat_max_abs_diff_on_a": float(np.abs(a1 - a2).max()),
        }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writers = {out_path.name: lambda p: runinfo.write_json(p, result)}
    runinfo.atomic_write_all(out_path.parent, writers)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--out", default=None, help="default: phaseb/j0/gpu_cpu_diff.json")
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    parser.add_argument("--device-a", default="cuda")
    parser.add_argument("--device-b", default="cpu")
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--n-long", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--batch-size", type=int, default=8)
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    out = Path(args.out) if args.out else work / "phaseb" / "j0" / "gpu_cpu_diff.json"
    try:
        import torch

        devices = [torch.device(args.device_a), torch.device(args.device_b)]
        if any(d.type == "cuda" for d in devices) and not torch.cuda.is_available():
            raise RuntimeError("a cuda device was requested but torch.cuda.is_available() is False")
        result = run(
            work / "truth_sequences.tsv.gz",
            out,
            args.models,
            devices[0],
            devices[1],
            args.n,
            args.n_long,
            args.seed,
            args.batch_size,
        )
    except (RuntimeError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result["models"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
