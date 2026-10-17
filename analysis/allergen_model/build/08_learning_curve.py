#!/usr/bin/env python3.12
"""Learning curve for M1 (ESM-2 mean embedding + logistic regression), spec sections 5 and 6.

Training subsets are drawn by CLUSTER (not by protein), so a cluster is in or out of a subset as a whole.
Fractions: 25%, 50%, 75%, 100% of the training clusters. Five random draws per fraction (100% has one, since it is the full set).
Evaluation: validation and test AUROC. Confidence intervals: 95% cluster bootstrap (1000 resamples of the evaluation clusters).
Reporting: every number carries the cluster count and the protein count behind it.

Inputs:
  --embeddings  npz with ids (P_<acc> or N_<acc>) and matrix
  --splits      splits.tsv (id, class, accession, split, cluster, length, length_bin)
  --outdir      output folder
"""

import argparse
import datetime
import os
import random

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

FRACTIONS = [0.25, 0.50, 0.75, 1.00]
DRAWS = 5
BOOT = 1000


def load_splits(path):
    import csv

    rows = list(csv.DictReader(open(path), delimiter="\t"))
    return rows


def ci95(vals):
    v = np.array(vals)
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def cluster_bootstrap_auc(y, s, clusters, rng):
    """95% CI of AUROC by resampling evaluation clusters with replacement."""
    uniq = sorted(set(clusters))
    idx_by = {c: [] for c in uniq}
    for i, c in enumerate(clusters):
        idx_by[c].append(i)
    out = []
    for _ in range(BOOT):
        pick = [rng.choice(uniq) for _ in uniq]
        ii = [i for c in pick for i in idx_by[c]]
        yy = y[ii]
        if len(set(yy)) < 2:
            continue
        out.append(roc_auc_score(yy, s[ii]))
    return ci95(out) if out else (None, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--extra-splits", default="", help="family_splits.tsv: added to the rows")
    ap.add_argument(
        "--train-negatives",
        default="weak_negative",
        help="comma list of classes used for training negatives",
    )
    ap.add_argument(
        "--eval-negatives",
        default="weak_negative",
        help="comma list of classes used for evaluation negatives",
    )
    args = ap.parse_args()
    train_neg = set(args.train_negatives.split(","))
    eval_neg = set(args.eval_negatives.split(","))
    os.makedirs(args.outdir, exist_ok=True)

    z = np.load(args.embeddings)
    ids = [str(i) for i in z["ids"]]
    X = np.asarray(z["matrix"], dtype=float)
    pos_in = {i: k for k, i in enumerate(ids)}
    rows = load_splits(args.splits)
    if args.extra_splits:
        rows += load_splits(args.extra_splits)
    missing = [r["id"] for r in rows if r["id"] not in pos_in]
    if missing:
        raise SystemExit(f"FATAL: {len(missing)} split rows have no embedding, e.g. {missing[:3]}")
    rows = [
        r for r in rows if r["id"] in pos_in
    ]  # embeddings may cover more proteins than the splits use

    def subset(split, negs):
        rs = [
            r
            for r in rows
            if r["split"] == split and (r["class"] == "positive" or r["class"] in negs)
        ]
        idx = np.array([pos_in[r["id"]] for r in rs])
        y = np.array([1 if r["class"] == "positive" else 0 for r in rs])
        cl = [r["cluster"] for r in rs]
        return idx, y, cl

    tr_idx, tr_y, tr_cl = subset("train", train_neg)
    va_idx, va_y, va_cl = subset("validation", eval_neg)
    te_idx, te_y, te_cl = subset("test", eval_neg)
    train_clusters = sorted(set(tr_cl))
    rng = random.Random(20261010)
    nrng = random.Random(7)

    results = []
    for f in FRACTIONS:
        draws = 1 if f == 1.0 else DRAWS
        va_auc, te_auc, nclu, npos, nneg = [], [], [], [], []
        for _ in range(draws):
            if f == 1.0:
                chosen = set(train_clusters)
            else:
                k = max(2, round(f * len(train_clusters)))
                chosen = set(rng.sample(train_clusters, k))
            m = np.array([c in chosen for c in tr_cl])
            Xt, yt = X[tr_idx[m]], tr_y[m]
            nclu.append(len(chosen))
            npos.append(int(yt.sum()))
            nneg.append(int(len(yt) - yt.sum()))
            clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
            clf.fit(Xt, yt)
            va_auc.append(roc_auc_score(va_y, clf.decision_function(X[va_idx])))
            te_auc.append(roc_auc_score(te_y, clf.decision_function(X[te_idx])))
        results.append(
            {
                "fraction": f,
                "draws": draws,
                "train_clusters": nclu,
                "train_pos": npos,
                "train_neg": nneg,
                "val_auc": va_auc,
                "test_auc": te_auc,
            }
        )

    # Confidence interval for the full-training test AUROC (single model), by cluster bootstrap.
    full = results[-1]
    clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
    clf.fit(X[tr_idx], tr_y)
    te_score = clf.decision_function(X[te_idx])
    te_ci = cluster_bootstrap_auc(te_y, te_score, te_cl, nrng)
    va_score = clf.decision_function(X[va_idx])
    va_ci = cluster_bootstrap_auc(va_y, va_score, va_cl, nrng)

    today = datetime.date.today().isoformat()
    with open(os.path.join(args.outdir, "learning_curve.md"), "w") as fh:
        fh.write("# Learning curve, M1 (ESM-2 650M mean embedding + logistic regression)\n\n")
        fh.write(
            f"- Training negatives: {sorted(train_neg)}; evaluation negatives: {sorted(eval_neg)}\n"
        )
        fh.write(f"- Run date: {today}\n")
        fh.write(
            f"- Training clusters available: {len(train_clusters)} "
            f"({len(tr_idx)} proteins: {int(tr_y.sum())} positives, {int(len(tr_y) - tr_y.sum())} negatives)\n"
        )
        fh.write(
            f"- Validation: {len(va_idx)} proteins in {len(set(va_cl))} clusters; "
            f"test: {len(te_idx)} proteins in {len(set(te_cl))} clusters\n"
        )
        fh.write(f"- Draws per fraction: {DRAWS} (100% has one, the full training set)\n\n")
        fh.write(
            "| Train fraction | Draws | Train clusters (mean) | Train positives (mean) | Validation AUROC (mean, sd) | Test AUROC (mean, sd) |\n"
        )
        fh.write("|---|---|---|---|---|---|\n")
        for r in results:
            vm, vs = (
                float(np.mean(r["val_auc"])),
                float(np.std(r["val_auc"], ddof=1)) if len(r["val_auc"]) > 1 else 0.0,
            )
            tm, ts = (
                float(np.mean(r["test_auc"])),
                float(np.std(r["test_auc"], ddof=1)) if len(r["test_auc"]) > 1 else 0.0,
            )
            fh.write(
                f"| {int(r['fraction'] * 100)}% | {r['draws']} | {np.mean(r['train_clusters']):.1f} | "
                f"{np.mean(r['train_pos']):.1f} | {vm:.3f} ({vs:.3f}) | {tm:.3f} ({ts:.3f}) |\n"
            )
        fh.write("\n## Full training set, with 95% cluster-bootstrap intervals\n\n")
        fh.write(
            f"- Test AUROC {full['test_auc'][0]:.3f}, 95% CI {te_ci[0]:.3f} to {te_ci[1]:.3f}, "
            f"from {len(set(te_cl))} test clusters ({len(te_idx)} proteins)\n"
        )
        fh.write(
            f"- Validation AUROC {full['val_auc'][0]:.3f}, 95% CI {va_ci[0]:.3f} to {va_ci[1]:.3f}, "
            f"from {len(set(va_cl))} validation clusters ({len(va_idx)} proteins)\n"
        )
        fh.write(
            "\nNo result here is a specificity claim. The split is below the 500-negative planning floor "
            "(spec section 3.2).\n"
        )
    print(open(os.path.join(args.outdir, "learning_curve.md")).read())


if __name__ == "__main__":
    main()
