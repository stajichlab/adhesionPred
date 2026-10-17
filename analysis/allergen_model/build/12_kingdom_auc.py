#!/usr/bin/env python3.12
"""Within-kingdom AUROC for the all-kingdom model (plan T4, reporting). Trains on the full all-kingdom training split,
scores validation and test, and reports AUROC for plant and animal groups separately, with cluster counts and 95%
cluster-bootstrap intervals."""

import csv
import random

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ROOT = "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/allkingdom/"
z = np.load(ROOT + "embed/embeddings.npz")
ids = [str(i) for i in z["ids"]]
X = np.asarray(z["matrix"], float)
pos = {i: k for k, i in enumerate(ids)}
group = {}
for r in csv.DictReader(
    open(ROOT + "../allergenonline_all/allergenonline_fungi.tsv"), delimiter="\t"
):
    if r["new"] == "yes" and r["uniprot"]:
        group["P_" + r["uniprot"]] = (
            "plant" if r["type"] in ("Aero Plant", "Food Plant", "Gliadin") else "animal"
        )
for r in csv.DictReader(open(ROOT + "kingdom_negatives.tsv"), delimiter="\t"):
    group["K_" + r["accession"]] = r["group"]
for r in csv.DictReader(open(ROOT + "family/family_negatives.tsv"), delimiter="\t"):
    group["K_" + r["accession"]] = r["group"]
rows = list(csv.DictReader(open(ROOT + "unified/allkingdom_splits.tsv"), delimiter="\t"))


def sub(split):
    return [r for r in rows if r["split"] == split]


tr = sub("train")
clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
clf.fit(X[[pos[r["id"]] for r in tr]], [1 if r["class"] == "positive" else 0 for r in tr])
rng = random.Random(7)
out = []
for split in ("validation", "test"):
    rs = sub(split)
    for g in ("plant", "animal"):
        gr = [r for r in rs if group.get(r["id"]) == g]
        y = np.array([1 if r["class"] == "positive" else 0 for r in gr])
        s = clf.decision_function(X[[pos[r["id"]] for r in gr]])
        cl = [r["cluster"] for r in gr]
        auc = roc_auc_score(y, s)
        uniq = sorted(set(cl))
        idx = {c: [] for c in uniq}
        for k, c in enumerate(cl):
            idx[c].append(k)
        boots = []
        for _ in range(1000):
            pick = [rng.choice(uniq) for _ in uniq]
            ii = [k for c in pick for k in idx[c]]
            if len(set(y[ii])) < 2:
                continue
            boots.append(roc_auc_score(y[ii], s[ii]))
        lo, hi = np.percentile(boots, [2.5, 97.5])
        out.append((split, g, int(y.sum()), int((1 - y).sum()), len(uniq), auc, lo, hi))
with open(ROOT + "kingdom_auc.tsv", "w") as f:
    f.write("split\tgroup\tpositives\tnegatives\tclusters\tauroc\tci_low\tci_high\n")
    for o in out:
        f.write("\t".join(str(round(x, 3)) if isinstance(x, float) else str(x) for x in o) + "\n")
for o in out:
    print(o)
