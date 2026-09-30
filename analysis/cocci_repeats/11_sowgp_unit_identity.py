#!/usr/bin/python3.12
"""Unit-against-unit difference matrices: which SOWgp unit was duplicated?

Input : sowgp_repeat_viz.copies.tsv (09)
Output:
  sowgp_unit_identity.png   top: within-allele unit x unit aa differences for the modal
                            (most strains) regular allele of each species x unit count;
                            bottom: each unit of the expanded modal alleles vs the units
                            of the same species' modal 3-unit allele
  sowgp_unit_identity.tsv   every regular distinct allele: pairwise unit differences

Differences = mismatches + gap columns in a global alignment of the two 47 aa units
(sowgp_units.unit_diff). Two identical units (0) in an allele are the simplest
signature of a recent duplication, but identity alone does not say which copy is the
source, and gene conversion also homogenizes units.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import sowgp_units as su  # noqa: E402

CAP = 6  # colour saturates at >= CAP differences; numbers are always printed
# one-hue blue ramp (validated default), dark = identical
CMAP = ListedColormap(["#104281", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb", "#eef4fc"])

H = su.distinct_alleles(HERE / "sowgp_repeat_viz.copies.tsv")
H = H[H.regular].reset_index(drop=True)


def dmat(us_a, us_b):
    return np.array([[su.unit_diff(a, b) for b in us_b] for a in us_a])


def short(c):
    return {"terminal": "term"}.get(c, c)


def heat(ax, M, rlab, clab, title):
    ax.imshow(np.minimum(M, CAP), cmap=CMAP, vmin=0, vmax=CAP)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(
                j,
                i,
                str(M[i, j]),
                ha="center",
                va="center",
                fontsize=8,
                color="white" if M[i, j] <= 2 else su.INK,
            )
    ax.set_xticks(range(len(clab)))
    ax.set_xticklabels(clab, rotation=45, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(rlab)))
    ax.set_yticklabels(rlab, fontsize=6.5)
    ax.set_title(title, fontsize=8, loc="left")
    for s in ax.spines.values():
        s.set_visible(False)


# ---- TSV: all regular alleles ----
rows = []
for _, r in H.iterrows():
    M = dmat(r.unit_seqs, r.unit_seqs)
    for i in range(len(M)):
        for j in range(i + 1, len(M)):
            rows.append(
                {
                    "allele": r.allele,
                    "species": r.species,
                    "n_strains": r.n,
                    "n_units": r.n_units,
                    "unit_i": i + 1,
                    "unit_j": j + 1,
                    "class_i": r.classes[i],
                    "class_j": r.classes[j],
                    "aa_diff": M[i, j],
                }
            )
pd.DataFrame(rows).to_csv(HERE / "sowgp_unit_identity.tsv", sep="\t", index=False)

# ---- figure ----
modal = {k: g.loc[g.n.idxmax()] for k, g in H[H.n_units >= 3].groupby(["species", "n_units"])}
ks = sorted({k for _, k in modal})
species = list(su.SPECIES_INK)
fig, axes = plt.subplots(4, len(ks), figsize=(3.3 * len(ks), 13.5))
for si, sp in enumerate(species):
    base = modal.get((sp, 3))
    for ki, k in enumerate(ks):
        ax_w, ax_b = axes[si, ki], axes[2 + si, ki]
        r = modal.get((sp, k))
        if r is None:
            ax_w.axis("off")
            ax_b.axis("off")
            continue
        lab = [f"{i + 1} {short(c)}" for i, c in enumerate(r.classes)]
        nall = int((H.species == sp).mul(H.n_units == k).sum())
        heat(
            ax_w,
            dmat(r.unit_seqs, r.unit_seqs),
            lab,
            lab,
            f"{sp} {k} units: {r.allele}\n{r.n} strain(s); {nall} distinct allele(s)",
        )
        if k == 3 or base is None:
            ax_b.axis("off")
            continue
        blab = [f"{i + 1} {short(c)}" for i, c in enumerate(base.classes)]
        heat(
            ax_b,
            dmat(r.unit_seqs, base.unit_seqs),
            lab,
            blab,
            f"{sp} {k}-unit {r.allele} (rows)\nvs 3-unit {base.allele} (cols)",
        )
fig.text(
    0.01,
    0.985,
    "Within-allele unit × unit aa differences " "(modal allele per species and unit count)",
    fontsize=11,
    va="top",
)
fig.text(
    0.01,
    0.5,
    "Expanded allele units vs the modal 3-unit allele of the same species",
    fontsize=11,
    va="top",
)
fig.text(
    0.01,
    0.005,
    f"Cells: aa differences (mismatch + gap) between 47 aa units; colour "
    f"saturates at {CAP}. 0 = identical units. Irregular-spacing alleles excluded.",
    fontsize=8,
    color=su.MUTED,
)
fig.tight_layout(rect=(0, 0.015, 1, 0.975), h_pad=3)
fig.savefig(HERE / "sowgp_unit_identity.png", dpi=170)
print(f"{len(H)} regular alleles; wrote sowgp_unit_identity.png/.tsv", file=sys.stderr)
