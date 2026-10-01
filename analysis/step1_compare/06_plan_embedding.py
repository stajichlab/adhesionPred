#!/usr/bin/env python3
"""Phase B step 2: cut the unique sequences into embedding chunks and size the J2 jobs.

Reads $STEP1_WORKDIR/phaseb/unique_sequences.tsv.gz and the J0 throughput JSON
($STEP1_WORKDIR/phaseb/j0/throughput.json). Writes to $STEP1_WORKDIR/phaseb/:

  chunk_plan.tsv          one row per chunk (chunk_id, window, n_seqs, residues, members_sha256)
  chunk_members.tsv.gz    one row per chunk member (chunk_id, position, row, seq_sha256)
  job_plan.json           the formula inputs (rates, load times, residues), chunk_residues,
                          n_jobs, time_minutes, batch size per model

The formula is in chunk_plan.py. --rate gives one assumed residues/s rate for every model
instead of J0; the plan then records `rate_source: assumed`. Use it only for a dry plan.
STOP (exit 2, no output): no throughput JSON and no --rate; a model without a successful J0 run.
"""

import argparse
import json
import sys
from pathlib import Path

import chunk_plan
import manifest
import paths
import runinfo
import truth_table

OUTPUT_NAMES = ("chunk_plan.tsv", "chunk_members.tsv.gz", "job_plan.json")
MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")


class PlanError(ValueError):
    """The throughput input is missing or incomplete."""


def rates_from_j0(j0: dict, models) -> tuple[dict, dict, dict]:
    rates, load_s, batch = {}, {}, {}
    for m in models:
        best = j0.get("best", {}).get(m)
        if not best or best.get("residues_per_s", 0) <= 0:
            raise PlanError(f"J0 throughput has no successful run for {m}")
        rates[m] = float(best["residues_per_s"])
        batch[m] = int(best["batch_size"])
        load_s[m] = float(j0.get("model_load_s", {}).get(m, 0.0))
    return rates, load_s, batch


def run(unique, rates, load_s, batch, out_dir: Path, chunk_residues=None, provenance=None):
    size = chunk_residues or chunk_plan.chunk_residues_from_rate(min(rates.values()))
    chunks = []
    for window in ("nterm", "cterm"):
        chunks += chunk_plan.build_chunks(unique, window, size)
    residues = sum(c.residues for c in chunks)
    plan = {
        "models": list(rates),
        "rates_residues_per_s": rates,
        "model_load_s": load_s,
        "batch_size": batch,
        "chunk_residues": size,
        "chunks": len(chunks),
        "residues_per_model": residues,
        "unique_sequences": len(unique),
        **chunk_plan.plan_jobs(rates, load_s, residues),
    }
    plan.update(provenance or {})
    plan_rows = [
        {
            "chunk_id": c.chunk_id,
            "window": c.window,
            "n_seqs": str(len(c.rows)),
            "residues": str(c.residues),
            "members_sha256": c.members_sha256,
        }
        for c in chunks
    ]
    member_rows = [
        {"chunk_id": c.chunk_id, "position": str(i), "row": str(r), "seq_sha256": h}
        for c in chunks
        for i, (r, h) in enumerate(zip(c.rows, c.hashes, strict=True))
    ]
    runinfo.atomic_write_all(
        out_dir,
        {
            "chunk_plan.tsv": lambda p: truth_table.write_tsv(
                p, chunk_plan.PLAN_COLUMNS, plan_rows
            ),
            "chunk_members.tsv.gz": lambda p: truth_table.write_tsv(
                p, chunk_plan.MEMBER_COLUMNS, member_rows
            ),
            "job_plan.json": lambda p: runinfo.write_json(p, plan),
        },
    )
    return chunks, plan


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--throughput", default=None, help="default: phaseb/j0/throughput.json")
    parser.add_argument("--rate", type=float, default=None, help="assumed residues/s (dry plan)")
    parser.add_argument("--chunk-residues", type=int, default=None)
    parser.add_argument("--models", nargs="*", default=list(MODELS))
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    out = work / "phaseb"
    j0_path = Path(args.throughput) if args.throughput else out / "j0" / "throughput.json"
    try:
        unique_path = out / "unique_sequences.tsv.gz"
        unique = truth_table.read_tsv(unique_path)
        if args.rate is not None:
            rates = {m: args.rate for m in args.models}
            load_s, batch = {m: 0.0 for m in args.models}, {m: 16 for m in args.models}
            source = {"rate_source": "assumed", "throughput_sha256": ""}
        else:
            if not j0_path.exists():
                raise PlanError(f"{j0_path} not found; run J0 first or give --rate for a dry plan")
            rates, load_s, batch = rates_from_j0(json.loads(j0_path.read_text()), args.models)
            source = {"rate_source": "J0", "throughput_sha256": manifest.sha256_file(j0_path)}
        provenance = {
            **source,
            "unique_sequences_sha256": manifest.sha256_file(unique_path),
            "git_commit": runinfo.git_commit(),
            "python": runinfo.python_version(),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        chunks, plan = run(unique, rates, load_s, batch, out, args.chunk_residues, provenance)
    except (PlanError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(
        f"chunks={plan['chunks']} chunk_residues={plan['chunk_residues']} "
        f"n_jobs={plan['n_jobs']} time_minutes={plan['time_minutes']} "
        f"rate_source={plan['rate_source']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
