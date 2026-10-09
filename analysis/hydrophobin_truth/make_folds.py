#!/usr/bin/python3.12
"""Task L2: leave-one-cluster-out folds over the T2 clusters, and the seeded tuning/test split of proteome species groups.

Folds come from clusters_positives.tsv restricted to T2 (28 clusters). The Pfam-missed clusters (the 9 T2 entries without a
Pfam hydrophobin-class cross-reference) are flagged. The 12 proteomes form 10 species groups (the three A. fumigatus strains are
one group); groups are split 5 and 5 with at least two groups with truth positives in each part. The split is made once.
Usage: make_folds.py --dir analysis/hydrophobin_truth
"""

import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path

SEED = 20261008
PFAM_MISSED = "A0A0A2KI06 A0A6V8R0V1 G9MKB2 I1RDP9 I1RDQ0 O94217 Q00367 Q0KKA0 Q7Z9L4".split()
GROUP_OF = {
    "Afum_Af293": "Afum",
    "Afum_A1163": "Afum",
    "Afum_W72310": "Afum",
    "Scer_S288C": "Scer",
    "Calb_SC5314": "Calb",
    "Cimm_RS": "Cimm",
    "Bder_ER3": "Bder",
    "Bbas_ARSEF2860": "Bbas",
    "Fful_Race5": "Fful",
    "Fgra_PH-1": "Fgra",
    "Pexp_MD-8": "Pexp",
    "Post_PC9": "Post",
}


def loco_folds(clusters, pfam_missed):
    by = {}
    for p, c in clusters.items():
        by.setdefault(c, []).append(p)
    folds = []
    for c in sorted(by):
        test = sorted(by[c])
        train = sorted(p for p in clusters if clusters[p] != c)
        folds.append(
            {
                "held_out_cluster": c,
                "test": test,
                "train": train,
                "pfam_missed": int(bool(set(test) & set(pfam_missed))),
            }
        )
    return folds


def split_species_groups(groups, seed, truth_groups, n_tuning=5):
    names = sorted(groups)
    rng = random.Random(seed)
    while True:
        rng.shuffle(names)
        tuning = set(names[:n_tuning])
        test = set(names[n_tuning:])
        ok = all(len([g for g in part if g in truth_groups]) >= 2 for part in (tuning, test))
        if ok:
            return {g: ("tuning" if g in tuning else "test") for g in groups}


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    out_split = d / "proteome_split.tsv"
    if out_split.exists():
        raise SystemExit(f"{out_split} exists: the split is made once")
    cl = {
        r["id"]: r["cluster"]
        for r in csv.DictReader(open(d / "clusters_positives.tsv"), delimiter="\t")
        if r["tier"] == "T2"
    }
    folds = loco_folds(cl, set(PFAM_MISSED))
    with open(d / "folds.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["fold", "held_out_cluster", "n_test", "n_train", "pfam_missed", "test_members"])
        for i, f in enumerate(folds, 1):
            w.writerow(
                [
                    i,
                    f["held_out_cluster"],
                    len(f["test"]),
                    len(f["train"]),
                    f["pfam_missed"],
                    ",".join(f["test"]),
                ]
            )
    runs = list(csv.DictReader(open(d / "run_list.tsv"), delimiter="\t"))
    names = [r["name"] for r in runs]
    groups = {}
    for n in names:
        groups.setdefault(GROUP_OF[n], []).append(n)
    taxa_truth = {
        g
        for g in groups
        if any((d / "calibration" / n / "truth_all.tsv").exists() for n in groups[g])
    }
    part = split_species_groups(groups, SEED, taxa_truth)
    with open(out_split, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["proteome", "group", "has_truth", "part", "seed"])
        for n in names:
            g = GROUP_OF[n]
            w.writerow([n, g, int(g in taxa_truth), part[g], SEED])
    meta = {
        "seed": SEED,
        "clusters_positives_sha256": sha(d / "clusters_positives.tsv"),
        "n_folds": len(folds),
    }
    json.dump(meta, open(d / "folds_meta.json", "w"), indent=1)
    print(f"{len(folds)} folds; {sum(f['pfam_missed'] for f in folds)} flagged Pfam-missed")
    print(
        "proteome split:",
        {p: sorted(g for g in groups if part[g] == p) for p in ("tuning", "test")},
    )
    print("groups with truth:", sorted(taxa_truth))
    return 0


if __name__ == "__main__":
    sys.exit(main())
