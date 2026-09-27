#!/usr/bin/env python
"""Stage-2 proof of concept: can adhesins be told apart from OTHER surface glycoproteins?

Stage 1 ("is this a cell-surface glycoprotein?") is close to solved; the shipped model
already does it. The open question this script answers is stage 2: within the surface
population, can a classifier separate curated adhesins from non-adhesive surface proteins?

Positives  : data/curated/surface/surface.tsv adhesion_status == adhesin
Negatives  : non_adhesin (curated) + non_adhesin_putative (enzyme-annotated)
Unknown    : excluded (they are unlabeled, not negatives)

Both CV schemes use MMseqs2 30%-identity clusters as groups, so a protein family never
spans a fold. leave-genome-out additionally tests transfer to an unseen species.

Feature sets compared:
  esmc_full    ESM C 300M, mean over all residues
  esmc_nterm   ESM C 300M, mean over the first 300 residues (adhesion domains are N-terminal)
  esmc_both    concatenation of the two
  aa_comp      20 amino-acid frequencies + log length (the baseline that saturated stage 1)
  arch         architecture only: GPI, signal peptide, log length

Usage:
    python analysis/model_review/stage2_proof_of_concept.py \
        --embeddings out_surface/esmc_300m/Surface.npz --clusters clu_cluster.tsv
"""

import argparse
import csv
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]
AA = "ACDEFGHIKLMNPQRSTVWY"


def aa_features(seq):
    n = max(len(seq), 1)
    return [seq.count(a) / n for a in AA] + [np.log10(n)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--embeddings", required=True, help="npz from embed_pilot.py (ids, full, nterm)")
    ap.add_argument("--clusters", required=True, help="MMseqs2 easy-cluster *_cluster.tsv")
    ap.add_argument("--fasta", help="FASTA of the same proteins (for aa_comp features)")
    ap.add_argument("--curated-negatives-only", action="store_true",
                    help="use only the 63 curated non-adhesins, excluding enzyme-annotated "
                         "proposals (the harder and more honest test)")
    args = ap.parse_args()

    labels = {}
    meta = {}

    with open(REPO / "data" / "curated" / "surface" / "surface.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            st = r["adhesion_status"]
            if st == "adhesin":
                labels[r["accession"]] = 1
            elif st == "non_adhesin":
                labels[r["accession"]] = 0
            elif st == "non_adhesin_putative" and not args.curated_negatives_only:
                labels[r["accession"]] = 0
            meta[r["accession"]] = r

    d = np.load(args.embeddings, allow_pickle=True)
    ids = [str(i).split("|")[1] if "|" in str(i) else str(i) for i in d["ids"]]
    idx = {a: i for i, a in enumerate(ids)}

    clusters = {}
    for line in open(args.clusters):
        rep, mem = line.rstrip("\n").split("\t")
        clusters[mem.split("|")[1] if "|" in mem else mem] = rep

    seqs = {}
    if args.fasta:
        acc = None
        for line in open(args.fasta):
            if line.startswith(">"):
                acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
                seqs[acc] = ""
            else:
                seqs[acc] += line.strip()

    keep = [a for a in labels if a in idx]
    y = np.array([labels[a] for a in keep])
    rows = np.array([idx[a] for a in keep])
    groups = np.array([hash(clusters.get(a, a)) for a in keep])
    genomes = np.array([meta[a]["genome"] for a in keep])

    X = {
        "esmc_full": d["full"][rows].astype(np.float32),
        "esmc_nterm": d["nterm"][rows].astype(np.float32),
        "esmc_both": np.hstack([d["full"][rows], d["nterm"][rows]]).astype(np.float32),
        "arch": np.array([[meta[a]["gpi_anchor"] == "yes", meta[a]["signal_peptide"] == "yes",
                           np.log10(max(int(meta[a]["length"]), 1))] for a in keep], dtype=np.float32),
    }
    if seqs:
        X["aa_comp"] = np.array([aa_features(seqs.get(a, "")) for a in keep], dtype=np.float32)

    print(f"{y.sum()} adhesins vs {(y == 0).sum()} non-adhesins, "
          f"{len(set(groups))} homology clusters, {len(set(genomes))} genomes\n")

    def clf():
        return make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced"))

    def evaluate(Xf, splits):
        p = np.full(len(y), np.nan)
        for tr, te in splits:
            if y[tr].sum() == 0 or y[te].sum() == 0:
                continue
            p[te] = clf().fit(Xf[tr], y[tr]).predict_proba(Xf[te])[:, 1]
        ok = ~np.isnan(p)
        order = np.argsort(-p[ok])
        top = y[ok][order][: int(y[ok].sum())]
        return (roc_auc_score(y[ok], p[ok]), average_precision_score(y[ok], p[ok]),
                top.sum() / max(y[ok].sum(), 1), ok.sum())

    cv = list(StratifiedGroupKFold(5, shuffle=True, random_state=0).split(y, y, groups))
    lgo = [(np.where(genomes != g)[0], np.where(genomes == g)[0]) for g in sorted(set(genomes))]

    print(f"{'features':<12}{'CV scheme':<22}{'ROC-AUC':>9}{'PR-AUC':>9}{'P@k':>7}{'n':>7}")
    for name, Xf in X.items():
        for scheme, splits in [("homology-grouped", cv), ("leave-genome-out", lgo)]:
            a, pr, pk, n = evaluate(Xf, splits)
            print(f"{name:<12}{scheme:<22}{a:>9.3f}{pr:>9.3f}{pk:>7.2f}{n:>7}")

    print("\nP@k = precision among the top-k ranked proteins, k = number of true adhesins.")
    print("PR-AUC is the metric to watch: the classes are heavily imbalanced.")


if __name__ == "__main__":
    main()
