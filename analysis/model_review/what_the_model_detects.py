#!/usr/bin/env python
"""What sequence feature actually drives the model's adhesin calls?

Motivation: the curated Eurotiomycetes labels showed the model finds BAD1, SOWgp and CspA
(all repeat-rich) but completely misses rodA, CalA and gp43. If that split is driven by
tandem-repeat content rather than by clade or by training-set composition, then the limit is
the label definition -- "adhesin" spans several physical mechanisms -- and not model capacity.

Correlates out-of-fold model score against three sequence properties for the curated
adhesins, under the same homology-grouped CV used elsewhere.

Usage:
    python analysis/model_review/what_the_model_detects.py --embeddings <Surface.npz> [...]
        --clusters <clu_cluster.tsv>
"""

import argparse
import csv
import urllib.parse
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]
import sys  # noqa: E402

sys.path.insert(0, str(REPO / "analysis" / "curation"))
from curation_lib import get  # noqa: E402


def low_complexity(seq):
    """Fraction of residues in the three most common amino acids."""
    if not seq:
        return 0.0
    counts = sorted(Counter(seq).values(), reverse=True)
    return sum(counts[:3]) / len(seq)


def repeat_coverage(seq, k=8, min_occurrences=3):
    """Fraction of the protein covered by any k-mer occurring >= min_occurrences times."""
    if len(seq) < k * min_occurrences:
        return 0.0
    counts = Counter(seq[i : i + k] for i in range(len(seq) - k + 1))
    repeated = [m for m, n in counts.items() if n >= min_occurrences]
    return min(1.0, sum(counts[m] * k for m in repeated) / len(seq))


def fetch_sequences(accessions):
    seqs, acc = {}, None
    accessions = sorted(accessions)
    for i in range(0, len(accessions), 90):
        q = urllib.parse.urlencode(
            {"accessions": ",".join(accessions[i : i + 90]), "format": "fasta"}
        )
        for line in get(
            f"https://rest.uniprot.org/uniprotkb/accessions?{q}", "text/plain"
        ).splitlines():
            if line.startswith(">"):
                acc = line.split("|")[1]
                seqs[acc] = ""
            elif acc:
                seqs[acc] += line.strip()
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--embeddings", nargs="+", required=True)
    ap.add_argument("--clusters", required=True)
    args = ap.parse_args()

    with open(REPO / "data" / "curated" / "surface" / "surface.tsv") as f:
        meta = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}

    emb = {}
    for path in args.embeddings:
        d = np.load(path, allow_pickle=True)
        for i, a in enumerate(d["ids"]):
            k = str(a).split("|")[1] if "|" in str(a) else str(a)
            emb[k] = d["full"][i].astype(np.float32)

    lab = {
        a: (1 if r["adhesion_status"] == "adhesin" else 0)
        for a, r in meta.items()
        if r["adhesion_status"] in ("adhesin", "non_adhesin") and a in emb
    }
    accs = sorted(lab)
    seqs = fetch_sequences(accs)

    clus = {}
    for line in open(args.clusters):
        rep, mem = line.rstrip("\n").split("\t")
        clus[mem.split("|")[1] if "|" in mem else mem] = rep
    cid = {c: i for i, c in enumerate(sorted({clus.get(a, a) for a in accs}))}

    y = np.array([lab[a] for a in accs])
    groups = np.array([cid[clus.get(a, a)] for a in accs])
    X = np.vstack([emb[a] for a in accs])
    p = np.full(len(y), np.nan)
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=0).split(X, y, groups):
        m = make_pipeline(
            StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced")
        ).fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[:, 1]

    pos = [i for i, a in enumerate(accs) if y[i] == 1 and a in seqs]
    props = {
        "low-complexity frac (top-3 aa)": np.array([low_complexity(seqs[accs[i]]) for i in pos]),
        "8-mer repeat coverage": np.array([repeat_coverage(seqs[accs[i]]) for i in pos]),
        "length": np.array([len(seqs[accs[i]]) for i in pos], dtype=float),
    }
    score = p[pos]
    print(f"Among {len(pos)} curated adhesins, does model score track sequence properties?")
    for name, v in props.items():
        rho, pv = spearmanr(v, score)
        print(f"   {name:<32} Spearman rho={rho:+.3f}  p={pv:.2g}")

    found = score > 0.5
    print(
        f"\n   FOUND  (p>0.5, n={found.sum()}):  "
        + "  ".join(f"{k.split()[0]} {np.median(v[found]):.3f}" for k, v in props.items())
    )
    print(
        f"   MISSED (p<0.5, n={(~found).sum()}):  "
        + "  ".join(f"{k.split()[0]} {np.median(v[~found]):.3f}" for k, v in props.items())
    )
    for name, v in props.items():
        if found.sum() and (~found).sum():
            _, pv = mannwhitneyu(v[found], v[~found], alternative="greater")
            print(f"   {name:<32} found > missed, Mann-Whitney p={pv:.3g}")

    print("\n   adhesins the model misses:")
    order = sorted(range(len(pos)), key=lambda i: score[i])
    for i in order[:14]:
        if score[i] >= 0.5:
            break
        a = accs[pos[i]]
        print(
            f"      p={score[i]:.3f} repeat={props['8-mer repeat coverage'][i]:.2f} "
            f"len={int(props['length'][i]):>5}  {(meta[a]['gene'] or a)[:10]:<11}"
            f"{meta[a]['genome'][:36]}"
        )


if __name__ == "__main__":
    main()
