#!/usr/bin/env python
"""Does a stage-2 adhesin model transfer between fungal subphyla?

Three tests, in increasing strictness:
  A) leave-one-genome-out, scored per held-out genome
  B) leave-one-CLADE-out: train on other subphyla entirely
  C) size-matched control: the cross-clade training set is tiny, so this asks whether a
     training set of the SAME size drawn from WITHIN Saccharomycotina also collapses.
     If it does not, the collapse in (B) is about clade, not about sample size.

Usage:
    python analysis/model_review/stage2_clade_transfer.py --embeddings <Surface.npz>
"""

import argparse
import csv
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]
CLADE = {
    "Saccharomyces cerevisiae S288C": "Saccharomycotina",
    "Candida albicans SC5314": "Saccharomycotina",
    "Nakaseomyces glabratus CBS138": "Saccharomycotina",
    "Candidozyma auris B8441": "Saccharomycotina",
    "Schizosaccharomyces pombe 972h-": "Taphrinomycotina",
    "Aspergillus fumigatus Af293": "Pezizomycotina",
}


def clf():
    return make_pipeline(
        StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced")
    )


def score(Xtr, ytr, Xte, yte):
    p = clf().fit(Xtr, ytr).predict_proba(Xte)[:, 1]
    return roc_auc_score(yte, p), average_precision_score(yte, p)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--resamples", type=int, default=200)
    args = ap.parse_args()

    lab, meta = {}, {}
    with open(REPO / "data" / "curated" / "surface" / "surface.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            meta[r["accession"]] = r
            if r["adhesion_status"] == "adhesin":
                lab[r["accession"]] = 1
            elif r["adhesion_status"] == "non_adhesin":
                lab[r["accession"]] = 0

    d = np.load(args.embeddings, allow_pickle=True)
    ids = [str(i).split("|")[1] if "|" in str(i) else str(i) for i in d["ids"]]
    idx = {a: i for i, a in enumerate(ids)}
    keep = [a for a in lab if a in idx]
    y = np.array([lab[a] for a in keep])
    X = d["full"][[idx[a] for a in keep]].astype(np.float32)
    g = np.array([meta[a]["genome"] for a in keep])
    cl = np.array([CLADE.get(x, "?") for x in g])

    print("A) leave-one-genome-out")
    print(f"{'held-out genome':<34}{'clade':<19}{'pos':>4}{'neg':>5}{'ROC':>8}{'PR':>8}")
    for gg in sorted(set(g)):
        te = g == gg
        if y[te].sum() == 0 or (y[te] == 0).sum() == 0:
            print(
                f"{gg[:33]:<34}{CLADE.get(gg,'?'):<19}{int(y[te].sum()):>4}"
                f"{int((y[te]==0).sum()):>5}{'n/a':>8}{'(one class)':>12}"
            )
            continue
        a, p = score(X[~te], y[~te], X[te], y[te])
        print(
            f"{gg[:33]:<34}{CLADE.get(gg,'?'):<19}{int(y[te].sum()):>4}"
            f"{int((y[te]==0).sum()):>5}{a:>8.3f}{p:>8.3f}"
        )

    print("\nB) leave-one-clade-out")
    cross = {}
    for c in sorted(set(cl)):
        te = cl == c
        if y[te].sum() == 0 or (y[te] == 0).sum() == 0:
            print(f"{c:<22} n/a (only one class present)")
            continue
        a, p = score(X[~te], y[~te], X[te], y[te])
        cross[c] = (a, p)
        print(
            f"{c:<22}{int(y[te].sum()):>4} pos {int((y[te]==0).sum()):>4} neg"
            f"   ROC {a:.3f}  PR {p:.3f}"
        )

    if "Saccharomycotina" not in cross:
        return
    n_pos = int((y[cl != "Saccharomycotina"] == 1).sum())
    n_neg = int((y[cl != "Saccharomycotina"] == 0).sum())
    print(
        f"\nC) size-matched control: train on {n_pos} pos + {n_neg} neg drawn from WITHIN "
        f"Saccharomycotina, {args.resamples} resamples"
    )
    m = cl == "Saccharomycotina"
    Xs, ys = X[m], y[m]
    rng = np.random.default_rng(0)
    pos, neg = np.where(ys == 1)[0], np.where(ys == 0)[0]
    rocs, prs = [], []
    for _ in range(args.resamples):
        tr = np.concatenate(
            [rng.choice(pos, n_pos, replace=False), rng.choice(neg, n_neg, replace=False)]
        )
        te = np.setdiff1d(np.arange(len(ys)), tr)
        a, p = score(Xs[tr], ys[tr], Xs[te], ys[te])
        rocs.append(a)
        prs.append(p)
    ca = cross["Saccharomycotina"][0]
    print(
        f"   within-clade ROC {np.mean(rocs):.3f} (5-95%: {np.percentile(rocs,5):.3f}-"
        f"{np.percentile(rocs,95):.3f})   PR {np.mean(prs):.3f}"
    )
    print(f"   cross-clade  ROC {ca:.3f}")
    print(f"   resamples scoring <= cross-clade: {np.mean(np.array(rocs) <= ca):.1%}")


if __name__ == "__main__":
    main()
