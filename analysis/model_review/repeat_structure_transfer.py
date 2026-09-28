#!/usr/bin/env python
"""Class 2a: can composition-agnostic REPEAT-STRUCTURE features transfer across clades?

Motivation (docs/TOOL-ARCHITECTURE.md §2): SOWgp and BAD1 sit in the FLO11 repeat class
(repeat coverage 1.00) but with inverted composition -- Pro/Cys-rich rather than Ser/Thr-rich.
If class 2a is to become clade-general rather than Saccharomycotina-only, it has to key on
repeat *structure* (period, copy number, periodicity) rather than on which residues repeat.

Trains on Saccharomycotina adhesins only, then asks whether the Onygenales/Eurotiales repeat
proteins are recovered.
"""

import argparse
import csv
import sys
import urllib.parse
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "analysis" / "curation"))
from curation_lib import get  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
PERIODS = (5, 8, 12, 20, 30)
LAGS = (5, 10, 20, 40)
SACCH = ("Saccharomyces", "Candida albicans", "Nakaseomyces", "Candidozyma")
# cross-clade repeat proteins: same structural class, different composition
TEST = {
    "Q8NK60": "SOWgp58",
    "Q8NK61": "SOWgp66",
    "Q96V71": "SOWgp82",
    "A4D962": "BAD1",
    "Q4WXC4": "CspA",
}


def repeat_desc(s):
    """Repeat descriptors that never reference WHICH residues repeat."""
    out = {}
    best_cov, best_k, best_copies = 0.0, 0, 0
    for k in PERIODS:
        cov, maxcopy = 0.0, 0
        if len(s) >= k * 3:
            c = Counter(s[i : i + k] for i in range(len(s) - k + 1))
            rep = {m: n for m, n in c.items() if n >= 3}
            if rep:
                cov = min(1.0, sum(c[m] * k for m in rep) / len(s))
                maxcopy = max(rep.values())
        out[f"cov_{k}"] = cov
        out[f"maxcopy_{k}"] = min(maxcopy, 50) / 50
        if cov > best_cov:
            best_cov, best_k, best_copies = cov, k, maxcopy
    out["best_cov"] = best_cov
    out["best_period"] = best_k / max(PERIODS)
    out["best_copies"] = min(best_copies, 50) / 50
    for d in LAGS:
        out[f"autocorr_{d}"] = (
            sum(1 for i in range(len(s) - d) if s[i] == s[i + d]) / (len(s) - d)
            if len(s) > d
            else 0.0
        )
    cnt = np.array([s.count(a) for a in AA], dtype=float)
    p = cnt / max(cnt.sum(), 1)
    out["entropy"] = float(-(p[p > 0] * np.log2(p[p > 0])).sum()) / 4.32
    out["loglen"] = np.log10(max(len(s), 1)) / 4
    return out


def comp_desc(s):
    n = max(len(s), 1)
    d = {f"f_{a}": s.count(a) / n for a in AA}
    d["loglen"] = np.log10(n) / 4
    return d


def fetch(accs):
    seqs, accs = {}, sorted(accs)
    for i in range(0, len(accs), 90):
        q = urllib.parse.urlencode({"accessions": ",".join(accs[i : i + 90]), "format": "fasta"})
        for blk in get(f"https://rest.uniprot.org/uniprotkb/accessions?{q}", "text/plain").split(
            ">"
        ):
            if blk.strip():
                h, *rest = blk.split("\n")
                seqs[h.split("|")[1]] = "".join(rest)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.parse_args()
    with open(REPO / "data" / "curated" / "surface" / "surface.tsv") as f:
        meta = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}
    lab = {
        a: (1 if r["adhesion_status"] == "adhesin" else 0)
        for a, r in meta.items()
        if r["adhesion_status"] in ("adhesin", "non_adhesin")
    }
    train = [a for a in lab if any(k in meta[a]["genome"] for k in SACCH)]
    seqs = fetch(set(train) | set(TEST))
    train = [a for a in train if a in seqs]
    y = np.array([lab[a] for a in train])
    print(f"train: {int(y.sum())} Saccharomycotina adhesins, {int((y == 0).sum())} non-adhesins")
    print("test: Onygenales/Eurotiales repeat proteins with INVERTED composition\n")

    for name, fn in [
        ("composition (20 aa freq + length)", comp_desc),
        ("repeat STRUCTURE (composition-agnostic)", repeat_desc),
    ]:
        keys = sorted(fn("ACDEFGHIKLMNPQRSTVWY" * 10))
        X = np.array([[fn(seqs[a]).get(k, 0.0) for k in keys] for a in train])
        m = make_pipeline(
            StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced")
        ).fit(X, y)
        found = 0
        print(f"{name}:")
        for acc, nm in TEST.items():
            if acc not in seqs:
                print(f"   {nm:<10} (sequence not retrieved)")
                continue
            p = m.predict_proba(np.array([[fn(seqs[acc]).get(k, 0.0) for k in keys]]))[0, 1]
            found += p > 0.5
            print(f"   {nm:<10} p={p:.3f}  {'FOUND' if p > 0.5 else 'missed'}")
        print(f"   -> {found}/{len(TEST)} recovered\n")


if __name__ == "__main__":
    main()
