#!/usr/bin/env python3.12
"""Make the T1.6 fixture: a small set of sequences with planted relationships.

Planted relationships (identity is approximate: each mutated position gets a different residue):
  P_B1  positive, train         base sequence B1
  P_B2  positive, train         base sequence B2 (unrelated to B1)
  P_B3  positive, val           base sequence B3
  P_B4  positive, test          base sequence B4
  P_B1m positive, test          B1 mutated to ~85% identity  -> planted cross-split violation with P_B1
  N_c1  candidate               B1 mutated to ~65% identity  -> should pair with P_B1 at T_C 0.6, not 0.7
  N_c2  candidate               B2 mutated to ~35% identity  -> should not pair at T_C 0.6
  N_c3  candidate               a 30-residue stretch of B3 inside random sequence -> fails the coverage rule
  N_c4  candidate               unrelated random sequence -> no hits

Writes positives.faa, candidates.faa, split.faa (prefixed train_/val_/test_) and expected.json.
Run: /usr/bin/python3.12 make_fixture.py <outdir>
"""

import json
import os
import random
import sys

AA = "ACDEFGHIKLMNPQRSTVWY"
SEED = 20261010
LEN = 150


def rand_seq(rng, n):
    return "".join(rng.choice(AA) for _ in range(n))


def mutate(rng, seq, identity):
    """Replace each position with a different residue with probability 1 - identity."""
    out = []
    for c in seq:
        if rng.random() < identity:
            out.append(c)
        else:
            out.append(rng.choice([a for a in AA if a != c]))
    return "".join(out)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(SEED)
    b1, b2, b3, b4 = (rand_seq(rng, LEN) for _ in range(4))

    positives = [
        ("P_B1", b1),
        ("P_B2", b2),
        ("P_B3", b3),
        ("P_B4", b4),
        ("P_B1m", mutate(rng, b1, 0.85)),
    ]
    motif = b3[60:90]
    candidates = [
        ("N_c1", mutate(rng, b1, 0.65)),
        ("N_c2", mutate(rng, b2, 0.35)),
        ("N_c3", rand_seq(rng, 60) + motif + rand_seq(rng, 60)),
        ("N_c4", rand_seq(rng, LEN)),
    ]
    split_of = {"P_B1": "train", "P_B2": "train", "P_B3": "val", "P_B4": "test", "P_B1m": "test"}
    split_records = [(f"{split_of[pid]}_{pid}", seq) for pid, seq in positives]
    split_records += [
        ("train_N_c1", candidates[0][1])
    ]  # a negative in train, for the cross-split check only

    write_fasta(os.path.join(outdir, "positives.faa"), positives)
    write_fasta(os.path.join(outdir, "candidates.faa"), candidates)
    write_fasta(os.path.join(outdir, "split.faa"), split_records)
    expected = {
        "pairs_at_T_C_0.6": {"N_c1": "P_B1"},
        "pairs_at_T_C_0.7": {
            "N_c1": "P_B1"
        },  # N_c1 measured at 0.74 identity; the planted 0.65 is approximate
        "fails_coverage": ["N_c3"],
        "no_hits": ["N_c4"],
        "cross_split_violation": ["P_B1", "P_B1m"],
        "cross_split_max_fident_min": 0.75,
        "cluster_same": [["P_B1", "P_B1m"], ["P_B1", "N_c1"]],
        "cluster_different": [["P_B1", "P_B2"]],
    }
    with open(os.path.join(outdir, "expected.json"), "w") as fh:
        json.dump(expected, fh, indent=2)
    print(f"wrote fixture to {outdir}")


if __name__ == "__main__":
    main(sys.argv[1])
