#!/usr/bin/python3.12
"""SOWgp repeat unit map (MVR-style) and length staircase, one row per distinct sequence.

Input : sowgp_repeat_viz.copies.tsv (09_sowgp_repeat_viz.py; full-length copies >= 250 aa)
Output: sowgp_unit_map.png, sowgp_unit_map.tsv (one row per distinct sequence)

(a) each 47 aa unit anchored on PTDCYGDC drawn at its protein position, coloured by
    unit class (sowgp_units.unit_class). A circle marks a unit that differs from the
    modal sequence of its class *within the same species* (so fixed species
    differences, e.g. immitis PPPP vs posadasii PPPPP, are not flagged). Hatching marks
    a unit whose anchor motif is mutated (inferred from 2-period spacing).
(b) protein length vs integer unit count; dotted line = base length + 47 aa per unit.
"""

import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import sowgp_units as su  # noqa: E402

H = su.distinct_alleles(HERE / "sowgp_repeat_viz.copies.tsv")

# modal sequence per (species, class), weighted by strain count
mode = {}
for _, r in H.iterrows():
    for u, c in zip(r.unit_seqs, r.classes, strict=False):
        mode.setdefault((r.species, c), Counter())[u] += r.n
mode = {k: v.most_common(1)[0][0] for k, v in mode.items()}

fig, (ax, bx) = plt.subplots(
    1, 2, figsize=(15, 0.16 * len(H) + 2.5), gridspec_kw={"width_ratios": [3, 1.1]}
)
for i, r in H.iterrows():
    y = len(H) - 1 - i
    L = len(r.seq)
    first = r.units[0][0]
    last_end = r.units[-1][0] + len(r.unit_seqs[-1])
    ax.plot([0, first], [y, y], color="#c9c8c1", lw=3, solid_capstyle="butt")
    ax.plot([last_end, L], [y, y], color="#c9c8c1", lw=3, solid_capstyle="butt")
    for (p, inferred), u, c in zip(r.units, r.unit_seqs, r.classes, strict=False):
        known = c in su.CLASS_COLOR
        ax.barh(
            y,
            len(u) - 2,
            left=p + 1,
            height=0.78,
            color=su.CLASS_COLOR.get(c, su.OTHER_COLOR),
            edgecolor=su.INK if (not known or inferred) else "none",
            lw=0.6,
            hatch="////" if inferred else None,
        )
        if u != mode.get((r.species, c)):
            ax.plot(p + len(u) / 2, y, marker="o", ms=3, color="white", mec=su.INK, mew=0.6)
    lab = f"{r.n:>3}×" + ("" if r.regular else "  irregular spacing")
    ax.text(
        L + 6,
        y,
        lab,
        va="center",
        fontsize=6.5,
        color=su.INK if r.regular else su.SPECIES_INK["posadasii"],
    )
ax.set_yticks([])
ax.set_ylim(-1, len(H))
ax.set_xlabel("position in protein (aa)")
ax.set_title(
    f"(a) Repeat units per distinct SOWgp protein sequence ({len(H)} rows, "
    f"{H.n.sum()} strains), anchored on {su.ANCHOR}\n"
    "label = strains with that exact sequence; ○ = differs from same-species "
    "modal unit of its class; hatched = anchor motif mutated; "
    "outlined empty = class not resolved",
    fontsize=9,
    loc="left",
)
for sp in su.SPECIES_INK:
    ys = [len(H) - 1 - i for i in H.index[H.species == sp]]
    ax.annotate(
        sp,
        xy=(-8, (min(ys) + max(ys)) / 2),
        ha="right",
        va="center",
        rotation=90,
        fontsize=10,
        color=su.SPECIES_INK[sp],
        annotation_clip=False,
    )
    ax.axhline(min(ys) - 0.5, color=su.GRID, lw=0.8)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in su.CLASS_COLOR.values()]
ax.legend(
    handles,
    [f"unit {k}" if k != "terminal" else "terminal unit" for k in su.CLASS_COLOR],
    loc="upper left",
    bbox_to_anchor=(0, -0.06 * 80 / len(H)),
    ncol=5,
    fontsize=8,
    frameon=False,
)

# (b) length staircase
agg = H.groupby(["species", "n_units", "length"])["n"].sum().reset_index()
off = {"immitis": -0.12, "posadasii": 0.12}
for sp, g in agg.groupby("species"):
    bx.scatter(
        g.n_units + off[sp],
        g.length,
        s=12 + g.n * 5,
        color=su.SPECIES_INK[sp],
        alpha=0.85,
        edgecolor="white",
        lw=0.8,
        label=sp,
        zorder=3,
    )
    for _, r in g[g.n >= 3].iterrows():
        right = sp == "posadasii"
        bx.text(
            r.n_units + off[sp] + (0.12 if right else -0.12),
            r.length,
            f"{r.length} aa (n={r.n})",
            fontsize=6.5,
            va="center",
            ha="left" if right else "right",
            color=su.MUTED,
        )
# base length = most common length among 3-unit alleles of each species
for sp in su.SPECIES_INK:
    g = agg[(agg.species == sp) & (agg.n_units == 3)]
    base = int(g.loc[g.n.idxmax(), "length"])
    xs = range(3, int(agg.n_units.max()) + 1)
    bx.plot(
        [x + off[sp] for x in xs],
        [base + su.PERIOD * (x - 3) for x in xs],
        color=su.SPECIES_INK[sp],
        lw=1,
        ls=":",
        zorder=2,
    )
bx.set_xticks(range(int(agg.n_units.min()), int(agg.n_units.max()) + 1))
bx.set_xlabel("repeat units (anchor count)")
bx.set_ylabel("protein length (aa)")
bx.set_title(
    f"(b) Length = base + {su.PERIOD} aa × units\n(dotted: expected step from the "
    "modal 3-unit length;\npoints off the line = indels or gene-model differences)",
    fontsize=9,
    loc="left",
)
bx.grid(axis="y", color=su.GRID, lw=0.6)
for s in ("top", "right"):
    bx.spines[s].set_visible(False)
bx.legend(frameon=False, fontsize=8, loc="upper left")

fig.tight_layout()
fig.savefig(HERE / "sowgp_unit_map.png", dpi=170, bbox_inches="tight")
fig.savefig(HERE / "sowgp_unit_map.pdf", bbox_inches="tight")
(
    H.assign(
        classes=H.classes.str.join(","),
        unit_starts=H.units.map(lambda us: ",".join(str(p) for p, _ in us)),
        anchor_inferred=H.units.map(lambda us: sum(i for _, i in us)),
    )
    .drop(columns=["seq", "units", "unit_seqs"])
    .to_csv(HERE / "sowgp_unit_map.tsv", sep="\t", index=False)
)
print(
    f"{len(H)} distinct sequences, {H.n.sum()} strains; wrote sowgp_unit_map.png/.tsv",
    file=sys.stderr,
)
