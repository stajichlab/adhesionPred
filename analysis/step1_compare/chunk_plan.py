"""Embedding chunks and job sizing for Phase B (J2). Standard library only.

A chunk is the unit of resumability: a fixed, ordered list of unique sequences for one window.
Entries are sorted by (window length, seq_sha256), so a chunk holds sequences of similar length
(little padding) and its membership depends only on the unique set and `chunk_residues`.

A job is a group of chunks. Job count and chunk size come from the measured J0 throughput:

  chunk_residues = max(MIN_CHUNK, floor(min_m(r_m) * CHUNK_SECONDS / 10,000) * 10,000)
  T_total        = sum_m (R / r_m + load_m)          R = residues of all windows to embed
  n_jobs         = max(1, ceil(T_total * SAFETY / TARGET_SECONDS))
  time_minutes   = ceil((T_total / n_jobs * 1.5 + 600) / 60)

r_m is residues per second of model m at its best batch size, and load_m its load time, both
from the J0 throughput JSON. Chunk k of the plan goes to job k mod n_jobs.
"""

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import manifest
import seqhash
import seqwindow
import truth_table

CHUNK_SECONDS = 300  # at most about 5 min of work is lost when a job is killed
MIN_CHUNK = 10_000
SAFETY = 1.25
TARGET_SECONDS = 4500  # 1.25 h: the middle of the 1 to 1.5 h job size rule
PLAN_COLUMNS = ("chunk_id", "window", "n_seqs", "residues", "members_sha256")
MEMBER_COLUMNS = ("chunk_id", "position", "row", "seq_sha256")


@dataclass(slots=True)
class Chunk:
    chunk_id: str
    window: str
    rows: list[int] = field(default_factory=list)
    hashes: list[str] = field(default_factory=list)
    residues: int = 0

    @property
    def members_sha256(self) -> str:
        return members_digest(self.hashes)


def members_digest(hashes) -> str:
    return hashlib.sha256("\n".join(hashes).encode("ascii")).hexdigest()


def window_entries(unique_rows, window: str) -> list[tuple[int, str, int]]:
    """(window length, seq_sha256, row) for every sequence that has this window."""
    out = []
    for r in unique_rows:
        seq = r["sequence"]
        if window == "cterm" and not seqwindow.needs_cterm(seq):
            continue
        out.append((len(seqwindow.window(seq, window)), r["seq_sha256"], int(r["row"])))
    return sorted(out)


def build_chunks(unique_rows, window: str, chunk_residues: int) -> list[Chunk]:
    if chunk_residues < 1:
        raise ValueError("chunk_residues must be positive")
    chunks: list[Chunk] = []
    current = None
    for length, digest, row in window_entries(unique_rows, window):
        if current is None or (current.rows and current.residues + length > chunk_residues):
            current = Chunk(f"{window}_{len(chunks):04d}", window)
            chunks.append(current)
        current.rows.append(row)
        current.hashes.append(digest)
        current.residues += length
    return chunks


def chunk_residues_from_rate(min_rate: float) -> int:
    if min_rate <= 0:
        raise ValueError("throughput must be positive")
    return max(MIN_CHUNK, int(min_rate * CHUNK_SECONDS // 10_000) * 10_000)


def plan_jobs(rates: dict[str, float], load_s: dict[str, float], residues: int) -> dict:
    """Job count and wall time for embedding `residues` residues with every model in `rates`."""
    if not rates:
        raise ValueError("no model rates")
    total = sum(residues / rates[m] + load_s.get(m, 0.0) for m in rates)
    n_jobs = max(1, math.ceil(total * SAFETY / TARGET_SECONDS))
    per_job = total / n_jobs
    return {
        "total_seconds": round(total, 1),
        "n_jobs": n_jobs,
        "seconds_per_job": round(per_job, 1),
        "time_minutes": math.ceil((per_job * 1.5 + 600) / 60),
    }


def chunks_for_job(chunk_ids: list[str], job_index: int, job_count: int) -> list[str]:
    if not 0 <= job_index < job_count:
        raise ValueError(f"job index {job_index} outside 0..{job_count - 1}")
    return [c for k, c in enumerate(chunk_ids) if k % job_count == job_index]


def read_plan(out_dir) -> list[Chunk]:
    """Chunks from chunk_plan.tsv and chunk_members.tsv.gz, checked against each other."""
    out_dir = Path(out_dir)
    plan = truth_table.read_tsv(out_dir / "chunk_plan.tsv")
    chunks = {r["chunk_id"]: Chunk(r["chunk_id"], r["window"]) for r in plan}
    for m in truth_table.read_tsv(out_dir / "chunk_members.tsv.gz"):
        chunk = chunks.get(m["chunk_id"])
        if chunk is None or int(m["position"]) != len(chunk.rows):
            raise ValueError(f"chunk_members.tsv.gz: unexpected row {m}")
        chunk.rows.append(int(m["row"]))
        chunk.hashes.append(m["seq_sha256"])
    for r in plan:
        chunk = chunks[r["chunk_id"]]
        if len(chunk.rows) != int(r["n_seqs"]) or chunk.members_sha256 != r["members_sha256"]:
            raise ValueError(f"chunk {r['chunk_id']}: members differ from chunk_plan.tsv")
        chunk.residues = int(r["residues"])
    return [chunks[r["chunk_id"]] for r in plan]


class StalePlanError(ValueError):
    """The chunk plan was made from another unique_sequences.tsv.gz."""


def check_plan_is_current(out_dir, chunks: list[Chunk], seq_by_row: dict[int, str]) -> None:
    """Raise StalePlanError unless the plan matches the current unique sequences.

    Two checks: the SHA-256 of unique_sequences.tsv.gz equals `unique_sequences_sha256` in
    job_plan.json, and the hash of the sequence in each member's row equals the member's
    planned seq_sha256."""
    out_dir = Path(out_dir)
    advice = "re-run 06_plan_embedding.py (and delete the old phaseb/emb/ chunks)"
    plan = json.loads((out_dir / "job_plan.json").read_text())
    current = manifest.sha256_file(out_dir / "unique_sequences.tsv.gz")
    if plan.get("unique_sequences_sha256") != current:
        raise StalePlanError(
            f"job_plan.json was made from another unique_sequences.tsv.gz "
            f"(recorded {plan.get('unique_sequences_sha256')!r}, current {current}); {advice}"
        )
    for chunk in chunks:
        for row, digest in zip(chunk.rows, chunk.hashes, strict=True):
            seq = seq_by_row.get(row)
            if seq is None or seqhash.seq_sha256(seq) != digest:
                raise StalePlanError(
                    f"chunk {chunk.chunk_id}: row {row} does not hold the planned sequence "
                    f"{digest[:12]}; {advice}"
                )
