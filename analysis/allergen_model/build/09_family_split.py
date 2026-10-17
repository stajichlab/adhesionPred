#!/usr/bin/env python3.12
"""Assign the family-matched negatives to the existing train/validation/test split, by cluster.

Clusters come from the family-matched run (family_matched/cluster/clu_cluster.tsv, positives and family negatives together).
Rule: a cluster that holds a positive takes that positive's split (positives are split by the earlier cluster split,
splits/splits.tsv). A cluster with no positive is assigned by seeded random draw, 70/15/15. A family negative is never
placed in a split other than its cluster's split.

Output: family_splits.tsv (id, class=family_negative, accession, split, cluster, family_cluster_rep).
"""

import csv
import os
import random
import sys

SEED = 20261010
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))


def main():
    clu_path = os.path.join(ROOT, "family_matched", "cluster", "clu_cluster.tsv")
    fam_path = os.path.join(ROOT, "family_matched", "family_negatives.tsv")
    split_path = os.path.join(ROOT, "splits", "splits.tsv")
    out_path = os.path.join(ROOT, "splits", "family_splits.tsv")

    member_rep = {}
    for ln in open(clu_path):
        rep, mem = ln.rstrip("\n").split("\t")
        member_rep[mem] = rep

    pos_split = {}
    for r in csv.DictReader(open(split_path), delimiter="\t"):
        if r["class"] == "positive":
            pos_split[r["accession"]] = r["split"]

    fam = [r["accession"] for r in csv.DictReader(open(fam_path), delimiter="\t")]
    rep_has_pos = {}
    rep_members = {}
    for a in fam:
        rep = member_rep.get(f"F_{a}") or member_rep.get(f"P_{a}") or member_rep.get(a)
        if rep is None:
            sys.exit(f"FATAL: family negative {a} is not in the family clustering")
        rep_members.setdefault(rep, []).append(a)
    for mem, rep in member_rep.items():
        if mem.startswith("P_") and mem[2:] in pos_split:
            rep_has_pos.setdefault(rep, set()).add(pos_split[mem[2:]])

    rng = random.Random(SEED)
    rand_split = {}
    for rep in sorted(rep_members):
        if rep in rep_has_pos:
            if len(rep_has_pos[rep]) > 1:
                sys.exit(f"FATAL: cluster {rep} holds positives from more than one split")
            continue
        x = rng.random()
        rand_split[rep] = "train" if x < 0.70 else ("validation" if x < 0.85 else "test")

    with open(out_path, "w") as fh:
        fh.write("id\tclass\taccession\tsplit\tcluster\tfamily_cluster_rep\n")
        for rep, accs in sorted(rep_members.items()):
            split = next(iter(rep_has_pos[rep])) if rep in rep_has_pos else rand_split[rep]
            for a in accs:
                fh.write(f"F_{a}\tfamily_negative\t{a}\t{split}\t{rep}\t{rep}\n")
    print(
        f"family negatives={len(fam)} clusters={len(rep_members)} with_positive={sum(1 for r in rep_members if r in rep_has_pos)}"
    )


if __name__ == "__main__":
    main()
