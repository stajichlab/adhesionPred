#!/usr/bin/env python
"""Status report: what can the stage-2 model do in a target clade, and what is validatable?

Reproduces the tables in docs/model-review/STATUS.md. Trains stage 2 on Saccharomycotina
only and applies it cold to the Eurotiomycetes genomes (Aspergillus + Coccidioides), then
scores the characterized surface proteins of those clades individually.

The point of the report is the label-coverage column: a call rate that looks plausible means
nothing when the target clade has no labels to check it against.

Usage:
    python analysis/model_review/clade_status_report.py --embeddings <Surface.npz>
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]
TRAIN_CLADE = {
    "Saccharomyces cerevisiae S288C",
    "Candida albicans SC5314",
    "Nakaseomyces glabratus CBS138",
    "Candidozyma auris B8441",
}
TARGET = {
    "Aspergillus fumigatus Af293",
    "Coccidioides immitis RS",
    "Coccidioides posadasii C735 delta SOWgp",
}
CHARACTERIZED = {
    "Q4WXJ1": "CalA, A. fumigatus invasin, binds integrin a5b1 (PMID 27841851)",
    "Q4WXC4": "CspA, A. fumigatus repeat-rich GPI cell-wall protein (PMID 20656913)",
    "Q96V71": "SOWgp82, Coccidioides spherule outer wall glycoprotein",
    "A0A0E1RVD3": "Ag2/PRA proline-rich antigen, Coccidioides",
    "A0ACG8DAS1": "Gel1, Coccidioides 1,3-beta-glucanosyltransferase (enzyme)",
}
STATUSES = ["adhesin", "non_adhesin", "non_adhesin_putative", "unknown"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--embeddings", required=True)
    args = ap.parse_args()

    with open(REPO / "data" / "curated" / "surface" / "surface.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    meta = {r["accession"]: r for r in rows}

    print("1. LABEL COVERAGE")
    counts = Counter((r["genome"], r["adhesion_status"]) for r in rows)
    print(f"{'genome':<42}" + "".join(f"{s[:19]:>21}" for s in STATUSES) + "   labeled%")
    for g in sorted({r["genome"] for r in rows}):
        n = sum(counts[(g, s)] for s in STATUSES)
        lab = counts[(g, "adhesin")] + counts[(g, "non_adhesin")]
        mark = "  <-- target" if g in TARGET else ""
        print(
            f"{g[:41]:<42}"
            + "".join(f"{counts[(g, s)]:>21}" for s in STATUSES)
            + f"{lab / n:>10.1%}{mark}"
        )

    d = np.load(args.embeddings, allow_pickle=True)
    ids = [str(i).split("|")[1] if "|" in str(i) else str(i) for i in d["ids"]]
    idx = {a: i for i, a in enumerate(ids)}

    tr = [
        a
        for a, r in meta.items()
        if r["genome"] in TRAIN_CLADE
        and r["adhesion_status"] in ("adhesin", "non_adhesin")
        and a in idx
    ]
    ytr = np.array([1 if meta[a]["adhesion_status"] == "adhesin" else 0 for a in tr])
    Xtr = d["full"][[idx[a] for a in tr]].astype(np.float32)
    model = make_pipeline(
        StandardScaler(), LogisticRegression(max_iter=5000, class_weight="balanced")
    ).fit(Xtr, ytr)

    print(
        f"\n2. SACCHAROMYCOTINA-TRAINED MODEL ({int(ytr.sum())} pos, {int((ytr == 0).sum())} neg) "
        f"APPLIED COLD"
    )
    for g in sorted(TARGET):
        sub = [a for a, r in meta.items() if r["genome"] == g and a in idx]
        p = model.predict_proba(d["full"][[idx[a] for a in sub]].astype(np.float32))[:, 1]
        lab = [meta[a]["adhesion_status"] for a in sub]
        gold = sum(1 for x in lab if x in ("adhesin", "non_adhesin"))
        print(
            f"\n  {g}: {len(sub)} surface proteins, {(p > 0.5).sum()} called "
            f"({(p > 0.5).mean():.1%}), {gold} with ground truth"
        )
        for s in STATUSES:
            k = [i for i, x in enumerate(lab) if x == s]
            if k:
                print(
                    f"     {s:<22} n={len(k):<4} called {int((p[k] > 0.5).sum()):<4} "
                    f"median p={np.median(p[k]):.3f}"
                )

    print("\n3. CHARACTERIZED PROTEINS OF THESE CLADES")
    for acc, desc in CHARACTERIZED.items():
        if acc not in idx:
            print(f"   {'n/a':>7}  (not in embedding set)  {acc}  {desc[:58]}")
            continue
        p = model.predict_proba(d["full"][[idx[acc]]].astype(np.float32))[:, 1][0]
        r = meta.get(acc, {})
        print(
            f"   p={p:.3f}  {'CALLED' if p > 0.5 else 'missed'}  {acc}  "
            f"len={r.get('length', '?'):>5} GPI={r.get('gpi_anchor', '?'):<4} {desc[:58]}"
        )


if __name__ == "__main__":
    main()
