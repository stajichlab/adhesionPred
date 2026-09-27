#!/usr/bin/env python
"""Does a protein language model find adhesins that HMM/domain annotation cannot?

This is the question that decides whether the PLM earns its place at all: if a Pfam rule
matches it, the cheaper and more interpretable method wins.

Compares, under homology-grouped CV on the curated labels:
  HMM rule        does the protein carry a known adhesin-family Pfam domain (binary)
  all Pfam        one-hot over every Pfam domain seen (the generous domain-based baseline)
  ESM C 300M      mean-pooled embedding
  ESM C + Pfam    both, to test whether domains add anything the embedding lacks

and then re-scores every model on the SUBSET WHERE DOMAIN ANNOTATION IS BLIND -- proteins
carrying no adhesin-family domain at all. That subset is where a PLM has to earn its keep.

Note on family choice: PF10528 is GLEYA, NOT Flo11 (which is PF10182) -- an easy and costly
mix-up. CPL1-like (PF21671) is listed because it is the Basidiomycota adhesin-family
candidate, but its adhesion role is unconfirmed; it does not affect these numbers, as the
curated label set contains no Basidiomycota proteins.

Usage:
    python analysis/model_review/esm_vs_hmm.py --embeddings <Surface.npz> --clusters <clu_cluster.tsv>
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
ADHESIN_PFAM = {
    "PF00624": "Flocculin",
    "PF13928": "Flocculin_t3",
    "PF05792": "Candida_ALS",
    "PF11766": "Candida_ALS_N",
    "PF11765": "Hyr1",
    "PF10182": "Flo11",
    "PF10528": "GLEYA",
    "PF07691": "PA14",
    "PF21671": "CPL1-like (unconfirmed)",
    "PF11763": "DIPSY",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--clusters", required=True)
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
    clus = {}
    for line in open(args.clusters):
        rep, mem = line.rstrip("\n").split("\t")
        clus[mem.split("|")[1] if "|" in mem else mem] = rep

    keep = [a for a in lab if a in idx]
    y = np.array([lab[a] for a in keep])
    pfsets = [{x for x in (meta[a]["pfam"] or "").split(";") if x} for a in keep]
    vocab = sorted({p for s in pfsets for p in s})
    cid = {c: i for i, c in enumerate(sorted({clus.get(a, a) for a in keep}))}
    groups = np.array([cid[clus.get(a, a)] for a in keep])

    emb = d["full"][[idx[a] for a in keep]].astype(np.float32)
    onehot = np.array([[p in s for p in vocab] for s in pfsets], dtype=float)
    feats = {
        "HMM rule (adhesin-family Pfam)": np.array(
            [[len(s & set(ADHESIN_PFAM)) > 0] for s in pfsets], dtype=float
        ),
        "all Pfam domains (one-hot)": onehot,
        "ESM C 300M": emb,
        "ESM C + Pfam": np.hstack([emb, onehot]),
    }

    def clf():
        return make_pipeline(
            StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced")
        )

    cv = list(StratifiedGroupKFold(5, shuffle=True, random_state=0).split(y, y, groups))
    blind = np.array([len(s & set(ADHESIN_PFAM)) == 0 for s in pfsets])

    print(
        f"{y.sum()} adhesins vs {(y == 0).sum()} curated non-adhesins, {len(cid)} homology clusters"
    )
    print(
        f"domain-blind subset: n={blind.sum()} "
        f"({int(y[blind].sum())} adhesins, {int((y[blind] == 0).sum())} non-adhesins)\n"
    )
    print(f"{'features':<34}{'ROC':>7}{'PR':>7}{'ROC|blind':>11}{'PR|blind':>10}")
    for name, Xf in feats.items():
        p = np.full(len(y), np.nan)
        for tr, te in cv:
            p[te] = clf().fit(Xf[tr], y[tr]).predict_proba(Xf[te])[:, 1]
        print(
            f"{name:<34}{roc_auc_score(y, p):>7.3f}{average_precision_score(y, p):>7.3f}"
            f"{roc_auc_score(y[blind], p[blind]):>11.3f}"
            f"{average_precision_score(y[blind], p[blind]):>10.3f}"
        )


if __name__ == "__main__":
    main()
