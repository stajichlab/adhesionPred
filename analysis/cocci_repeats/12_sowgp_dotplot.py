#!/usr/bin/python3.12
"""Self dot-plots of one SOWgp protein per species x unit count.

Input : sowgp_repeat_viz.copies.tsv (09)
Output: sowgp_dotplot.png

Each panel is the modal (most strains) regular allele for that species and unit count.
A dot at (i, j) means the windows of --window aa starting at i and j share at least
--min-match identical residues (ungapped). Off-diagonal lines at 47 aa spacing are the
tandem repeat; their number and length show the unit count. Grey marks anchored
units (PTDCYGDC): a square per unit on the diagonal, lines at unit starts.
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import sowgp_units as su  # noqa: E402


def dot_matrix(s, w, k):
    a = np.frombuffer(s.encode(), dtype=np.uint8)
    eq = (a[:, None] == a[None, :]).astype(np.int16)
    n = len(s) - w + 1
    # windowed diagonal sums via cumulative sum along diagonals
    M = np.zeros((n, n), dtype=np.int16)
    for d in range(w):
        M += eq[d : d + n, d : d + n]
    return M >= k


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--window", type=int, default=10)
    ap.add_argument("--min-match", type=int, default=7)
    args = ap.parse_args()

    H = su.distinct_alleles(HERE / "sowgp_repeat_viz.copies.tsv")
    H = H[H.regular & (H.n_units >= 3)]
    reps = [g.loc[g.n.idxmax()] for _, g in H.groupby(["species", "n_units"])]
    ncol = max(sum(r.species == sp for r in reps) for sp in su.SPECIES_INK)
    fig, axes = plt.subplots(2, ncol, figsize=(3.6 * ncol, 7.8))
    Lmax = max(len(r.seq) for r in reps)
    for si, sp in enumerate(su.SPECIES_INK):
        rs = [r for r in reps if r.species == sp]
        for ci in range(ncol):
            ax = axes[si, ci]
            if ci >= len(rs):
                ax.axis("off")
                continue
            r = rs[ci]
            M = dot_matrix(r.seq, args.window, args.min_match)
            for p, _ in r.units:  # one light square per anchored unit, on the diagonal
                ax.add_patch(
                    plt.Rectangle((p, p), su.PERIOD, su.PERIOD, color=su.GRID, lw=0, zorder=0)
                )
                ax.axvline(p, color=su.GRID, lw=0.5, zorder=0)
                ax.axhline(p, color=su.GRID, lw=0.5, zorder=0)
            yy, xx = np.nonzero(M)
            ax.scatter(xx, yy, s=0.4, color=su.SPECIES_INK[sp], lw=0, zorder=2)
            ax.set_xlim(0, Lmax)
            ax.set_ylim(Lmax, 0)
            ax.set_aspect("equal")
            ax.set_title(
                f"{sp} {r.n_units} units, {r.length} aa\n{r.allele}, " f"{r.n} strain(s)",
                fontsize=8,
                loc="left",
            )
            ax.tick_params(labelsize=6.5)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
    fig.suptitle(
        f"SOWgp self dot-plots (window {args.window} aa, >= {args.min_match} "
        "identical; grey squares/lines = anchored 47 aa units)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(HERE / "sowgp_dotplot.png", dpi=170)
    print(f"wrote sowgp_dotplot.png ({len(reps)} panels)", file=sys.stderr)


if __name__ == "__main__":
    main()
