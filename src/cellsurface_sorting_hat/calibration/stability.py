"""Stability of calls between two runs (for example two strains of one species).

Proteins are matched by sequence sha256 (``proteins.tsv.gz``), not by ID. Only proteins whose
sequence is identical in both runs are compared; the others are counted. A difference in calls for
an identical sequence comes from the taxon, the status or a module, not from the sequence.
"""

import csv
import gzip
from collections import Counter, defaultdict


def _read(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _calls(out_dir):
    calls = defaultdict(dict)
    for r in _read(f"{out_dir}/calls.long.tsv.gz"):
        calls[r["protein"]][(r["call"], r["variant"])] = r["value"]
    return calls


def compare_runs(out_a, out_b):
    """Return ``{"n_a", "n_b", "n_identical_sequences", "per_call": {(call, variant): Counter}}``."""
    pa = {r["id"]: r["sha256"] for r in _read(f"{out_a}/proteins.tsv.gz")}
    pb = {r["id"]: r["sha256"] for r in _read(f"{out_b}/proteins.tsv.gz")}
    by_sha_b = defaultdict(list)
    for pid, sha in pb.items():
        by_sha_b[sha].append(pid)
    ca, cb = _calls(out_a), _calls(out_b)
    per_call, n_identical, unique = defaultdict(Counter), 0, set()
    for pid, sha in pa.items():
        if sha not in by_sha_b:
            continue
        n_identical += 1
        unique.add(sha)
        other = sorted(by_sha_b[sha])[0]  # identical sequences have identical calls inside one run
        for key, value in ca[pid].items():
            per_call[key]["agree" if cb[other].get(key) == value else "differ"] += 1
    return {
        "n_a": len(pa),
        "n_b": len(pb),
        "n_identical_sequences": n_identical,  # proteins of run A, counted by ID
        "n_unique_identical_sequences": len(unique),
        "per_call": dict(per_call),
    }
