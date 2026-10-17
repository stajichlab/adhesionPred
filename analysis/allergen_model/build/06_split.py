#!/usr/bin/env python3.12
"""Cluster-level train/validation/test split with length-matched weak negatives (spec 3.2 and 4, pairing dropped 2026-10-10).

Inputs:
  --clusters  full_cluster/clusters.tsv   accession<TAB>cluster_id (IDs carry the P_ or N_ prefix from the FASTA)
  --positives positives/positives.tsv     accession, length, ...
  --candidates candidates/candidates.tsv  accession, genus, length, ...
  --outdir    output folder
  --seed      fixed seed (default 20261010)

Method:
  1. Positive clusters (clusters containing at least one positive) are shuffled with the seed and assigned to
     train 70% / validation 15% / test 15% by positive count. A cluster is never split across splits.
  2. Length-matched sampling: positives and candidates are binned by length (50-residue bins). For each bin, as many
     candidates as positives are sampled (1:1), at random with the seed. If a bin has fewer candidates, all are used.
  3. A sampled negative takes the split of its cluster if the cluster holds a positive. If not, its cluster is assigned
     a split by the same shuffled 70/15/15 rule, so all sampled members of one cluster share a split.

Outputs (in --outdir): splits.tsv (id, class, accession, split, cluster, length, length_bin),
  counts.md (per split and class, bins, and the 500-negative planning floor check).
"""

import argparse
import collections
import csv
import datetime
import os
import random
import sys

BIN = 50


def strip(i):
    return i.split("_", 1)[1] if "_" in i and i[:2] in ("P_", "N_") else i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--positives", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=20261010)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    clu = {}
    for ln in open(args.clusters):
        a, c = ln.rstrip("\n").split("\t")
        clu[strip(a)] = c

    pos = list(csv.DictReader(open(args.positives, encoding="utf-8"), delimiter="\t"))
    cand = list(csv.DictReader(open(args.candidates, encoding="utf-8"), delimiter="\t"))
    pos_len = {r["accession"]: int(r["length"]) for r in pos}
    cand_len = {r["accession"]: int(r["length"]) for r in cand}
    for a in pos_len:
        if a not in clu:
            sys.exit(f"FATAL: positive {a} has no cluster")
    for a in cand_len:
        if a not in clu:
            sys.exit(f"FATAL: candidate {a} has no cluster")

    # 1. Split positive clusters, by positive count.
    pos_clusters = collections.Counter(clu[a] for a in pos_len)
    order = sorted(pos_clusters)
    rng.shuffle(order)
    total = sum(pos_clusters.values())
    split_of_cluster = {}
    acc_n = 0
    for c in order:
        frac = acc_n / total
        if frac < 0.70:
            s = "train"
        elif frac < 0.85:
            s = "validation"
        else:
            s = "test"
        split_of_cluster[c] = s
        acc_n += pos_clusters[c]

    # 2. Length-matched sampling, 1:1 per length bin.
    pos_bins = collections.Counter(pos_len[a] // BIN for a in pos_len)
    cand_by_bin = collections.defaultdict(list)
    for a, L in cand_len.items():
        cand_by_bin[L // BIN].append(a)
    sampled = []
    short = {}
    for b in sorted(pos_bins):
        pool = sorted(cand_by_bin.get(b, []))
        k = pos_bins[b]
        take = min(k, len(pool))
        if take < k:
            short[b] = (k, len(pool))
        sampled += rng.sample(pool, take)

    # 3. Negative split by cluster.
    neg_cluster_split = {}
    for a in sampled:
        c = clu[a]
        if c in split_of_cluster:
            continue
        if c not in neg_cluster_split:
            r = rng.random()
            neg_cluster_split[c] = "train" if r < 0.70 else ("validation" if r < 0.85 else "test")

    rows = []
    for a in pos_len:
        c = clu[a]
        rows.append(
            [f"P_{a}", "positive", a, split_of_cluster[c], c, pos_len[a], pos_len[a] // BIN]
        )
    for a in sampled:
        c = clu[a]
        s = split_of_cluster.get(c) or neg_cluster_split[c]
        rows.append([f"N_{a}", "weak_negative", a, s, c, cand_len[a], cand_len[a] // BIN])

    with open(os.path.join(args.outdir, "splits.tsv"), "w") as fh:
        fh.write("id\tclass\taccession\tsplit\tcluster\tlength\tlength_bin\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    cnt = collections.Counter((r[3], r[1]) for r in rows)
    floor = 500
    today = datetime.date.today().isoformat()
    with open(os.path.join(args.outdir, "counts.md"), "w") as fh:
        fh.write("# Cluster split with length-matched weak negatives\n\n")
        fh.write(f"- Run date: {today}; seed {args.seed}; length bin {BIN} aa; 1:1 per bin\n")
        fh.write(f"- Positive clusters: {len(pos_clusters)}; sampled negatives: {len(sampled)}\n")
        fh.write(
            f"- Bins where candidates were fewer than positives: {short if short else 'none'}\n\n"
        )
        fh.write(
            "| Split | Positives | Weak negatives | Negatives meet the 500 planning floor |\n|---|---|---|---|\n"
        )
        for s in ("train", "validation", "test"):
            p = cnt[(s, "positive")]
            n = cnt[(s, "weak_negative")]
            fh.write(f"| {s} | {p} | {n} | {'yes' if n >= floor else 'NO'} |\n")
        fh.write(
            "\nThe 500 floor is a planning value (spec section 3.2). Below it, the split is not feasible for a "
            "specificity claim; AUROC is still reported and labelled.\n"
        )
    print("splits:", {f"{s}/{c}": v for (s, c), v in cnt.items()})


if __name__ == "__main__":
    main()
