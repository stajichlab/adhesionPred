#!/usr/bin/env python3.12
"""B6: nearest-neighbour shortcut baseline (spec section 5, owner stop margin O1).

For every validation and test protein, the score is its best fident (MMseqs2 fraction identity) to any TRAINING positive:
  - positives: from pos_vs_pos (positive vs positive); self-hits excluded
  - weak negatives: from the candidate-vs-positive search (full_pairs); candidates with no hit score 0
Identity is not coverage-filtered, consistent with a plain nearest-neighbour score.

AUROC is computed on validation and on test (positive = 1, weak negative = 0), by the Mann-Whitney rank formula.
Stop rule (owner decision O1): AUROC above 0.60 on test halts the arm for owner review before any embedding model.

Inputs: --splits splits/splits.tsv, --pos-pairs pos_vs_pos/pos_vs_pos.tsv.zst, --cand-pairs full_pairs/pairs_candidates_vs_positives.tsv.zst
Outputs (in --outdir): b6_scores.tsv, b6.md
"""

import argparse
import csv
import datetime
import os
import subprocess


def zst_rows(path):
    out = subprocess.run(
        ["zstd", "-dc", path], capture_output=True, text=True, check=True
    ).stdout.split("\n")
    h = out[0].split("\t")
    return [dict(zip(h, ln.split("\t"), strict=False)) for ln in out[1:] if ln]


def auroc(pos, neg):
    """Mann-Whitney AUROC: probability that a random positive scores above a random negative (ties count half)."""
    if not pos or not neg:
        return None
    wins = 0.0
    for p in pos:
        for n in neg:
            wins += 1.0 if p > n else (0.5 if p == n else 0.0)
    return wins / (len(pos) * len(neg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", required=True)
    ap.add_argument("--pos-pairs", required=True)
    ap.add_argument("--cand-pairs", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = list(csv.DictReader(open(args.splits), delimiter="\t"))
    train = {r["id"] for r in rows if r["split"] == "train" and r["class"] == "positive"}

    best = {}
    for r in zst_rows(args.pos_pairs):
        q, t = r["query"], r["target"]
        if q == t or t not in train:
            continue
        f = float(r["fident"])
        if f > best.get(q, 0.0):
            best[q] = f
    for r in zst_rows(args.cand_pairs):
        q, t = r["query"], r["target"]
        if t not in train:
            continue
        f = float(r["fident"])
        if f > best.get(q, 0.0):
            best[q] = f

    scores = []
    for r in rows:
        if r["split"] not in ("validation", "test"):
            continue
        s = best.get(r["id"], 0.0)
        scores.append((r["split"], r["class"], r["id"], s))

    with open(os.path.join(args.outdir, "b6_scores.tsv"), "w") as fh:
        fh.write("split\tclass\tid\tb6_best_identity_to_train_positive\n")
        for sp, cl, i, s in scores:
            fh.write(f"{sp}\t{cl}\t{i}\t{s:.4f}\n")

    res = {}
    for sp in ("validation", "test"):
        pos = [s for x, c, i, s in scores if x == sp and c == "positive"]
        neg = [s for x, c, i, s in scores if x == sp and c == "weak_negative"]
        res[sp] = (len(pos), len(neg), auroc(pos, neg))
    today = datetime.date.today().isoformat()
    with open(os.path.join(args.outdir, "b6.md"), "w") as fh:
        fh.write("# B6 nearest-neighbour baseline\n\n")
        fh.write(
            f"- Run date: {today}\n- Score: best fident to any training positive; no coverage filter.\n\n"
        )
        fh.write(
            "| Split | Positives | Weak negatives | AUROC | Stop (AUROC > 0.60) |\n|---|---|---|---|---|\n"
        )
        for sp in ("validation", "test"):
            npos, nneg, a = res[sp]
            stop = "n/a" if a is None else ("HALT" if a > 0.60 else "pass")
            fh.write(
                f"| {sp} | {npos} | {nneg} | {'n/a' if a is None else f'{a:.3f}'} | {stop} |\n"
            )
        fh.write(
            "\nNegatives are length-matched, not similarity-matched. A high B6 AUROC means similarity to "
            "training positives separates the classes, so the model result must be read against it.\n"
        )
    for sp in ("validation", "test"):
        print(sp, res[sp])


if __name__ == "__main__":
    main()
