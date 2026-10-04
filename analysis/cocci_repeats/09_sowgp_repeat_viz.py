#!/usr/bin/python3.12
"""Visualize the SOWgp tandem repeat array and its variability in full-length copies.

Question: how variable is the SOWgp repeat? Copies below ~250 aa are fragments of the
array (split gene models / collapsed short-read assemblies) and are excluded
(--min-len 250).

Repeat units are phase-anchored on the PTDCYGDC motif (sowgp_units.py): one unit per
anchor, 47 aa. An earlier version divided the profiler's repeat region into equal-length
pieces; those pieces drifted out of phase, and the fractional n_copies it reported
depended on where the region edges were called. n_units (anchor count) is an integer.

Outputs (prefix default: sowgp_repeat_viz):
  <prefix>.png             figure: architecture, unit-count, per-species unit logos, raster
  <prefix>.units.tsv       one row per repeat unit (anchored unit sequence, aligned form)
  <prefix>.copies.tsv      one row per full-length copy (profiler region + n_units)

Deps: matplotlib (Agg), Bio, pandas. Each unit is projected onto the modal internal unit
(47 columns) by pairwise global alignment. MAFFT --localpair was tried and rejected: the
terminal units, which run into the C-terminus, spread the alignment over ~72 columns.
"""

import argparse
import importlib.util
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import sowgp_units as su  # noqa: E402

_spec = importlib.util.spec_from_file_location("rp", str(HERE / "02_repeat_profile.py"))
rp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rp)

# residue colours by chemistry (standard logo scheme)
AA_COLOR = {}
for aas, c in (
    ("DE", "#c0392b"),
    ("KRH", "#2a78d6"),
    ("STNQ", "#1baf7a"),
    ("GP", "#eda100"),
    ("C", "#8e44ad"),
    ("AVLIMFWY", "#222220"),
):
    for a in aas:
        AA_COLOR[a] = c


def read_fasta(path):
    with open(path) as fh:
        name, buf = None, []
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def project_to_reference(seqs, ref):
    """Globally align each unit to `ref`; return strings in ref coordinates.

    Output length = len(ref). Residues inserted relative to ref are dropped; ref
    positions with no aligned residue are '-'. Units are already phase-anchored, so
    this only absorbs the 1 aa indels (e.g. PPPP vs PPPPP).
    """
    from Bio.Align import PairwiseAligner

    pa = PairwiseAligner()
    pa.mode = "global"
    pa.match_score, pa.mismatch_score = 2, -1
    pa.open_gap_score, pa.extend_gap_score = -4, -1
    pa.end_gap_score = 0  # terminal units end early
    out = []
    for s in seqs:
        a = pa.align(ref, s)[0]
        col = ["-"] * len(ref)
        for (r0, r1), (q0, _q1) in zip(*a.aligned, strict=False):
            for k in range(r1 - r0):
                col[r0 + k] = s[q0 + k]
        out.append("".join(col))
    return out


def draw_letter(ax, ch, x, y, h, color):
    """Draw one glyph scaled to width 0.9 and height h (data units) at (x, y)."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import PathPatch
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D

    tp = TextPath((0, 0), ch, size=1, prop=FontProperties(family="DejaVu Sans", weight="bold"))
    bb = tp.get_extents()
    tr = (
        Affine2D()
        .translate(-bb.x0, -bb.y0)
        .scale(0.9 / max(bb.width, 1e-6), h / max(bb.height, 1e-6))
        .translate(x - 0.45, y)
    )
    ax.add_patch(PathPatch(tr.transform_path(tp), color=color, lw=0))


# more saturated than su.SPECIES_INK so filled blocks stay visible at 79 rows
VIVID = {"immitis": "#1f6feb", "posadasii": "#e8501a"}


def panel_title(ax, title, size=13):
    """Left-aligned bold title; the leading "(x)" panel label is upper-cased."""
    label, _, rest = title.partition(" ")
    ax.set_title(f"{label.upper()} {rest}", fontsize=size, fontweight="bold", loc="left")


def logo(ax, aln, title):
    """Information-content logo (bits) from aligned strings; gap columns >50% dropped."""
    L = len(aln[0])
    cols = [j for j in range(L) if sum(a[j] == "-" for a in aln) <= 0.5 * len(aln)]
    for x, j in enumerate(cols):
        c = Counter(a[j] for a in aln if a[j] != "-")
        tot = sum(c.values())
        freqs = {k: v / tot for k, v in c.items()}
        H = -sum(f * math.log2(f) for f in freqs.values())
        ic = max(0.0, math.log2(20) - H)
        y = 0.0
        for ch, f in sorted(freqs.items(), key=lambda kv: kv[1]):
            h = f * ic
            if h > 0.02:
                draw_letter(ax, ch, x, y, h, AA_COLOR.get(ch, "#6b6a64"))
            y += h
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(0, math.log2(20))
    ax.set_ylabel("bits")
    ax.set_xticks(range(0, len(cols), 5))
    ax.set_xticklabels([str(i + 1) for i in range(0, len(cols), 5)], fontsize=7)
    panel_title(ax, title)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return cols


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pangenome-tsv", default=str(HERE / "sowgp_pangenome.tsv"))
    ap.add_argument("--pangenome-fa", default=str(HERE / "sowgp_pangenome.fa"))
    ap.add_argument("--seed-fa", default=str(HERE / "sowgp_seed.fa"))
    ap.add_argument("--prefix", default=str(HERE / "sowgp_repeat_viz"))
    ap.add_argument("--min-len", type=int, default=250, help="drop fragments below this")
    args = ap.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(args.pangenome_tsv, sep="\t")
    df = df[df["family"] == "SOWgp"].copy()
    seqs = {r.id.split("|")[0]: str(r.seq) for r in SeqIO.parse(args.pangenome_fa, "fasta")}
    df["seq"] = df["protein"].map(lambda p: seqs.get(p, ""))
    df = df[df["seq"].str.len() > 0].copy()
    df["species"] = df["species"].fillna("posadasii").replace({"Coccidioides": "posadasii"})

    # reference/long-read seeds: add as labeled rows, deduped by sequence
    seed_rows = []
    seen = set(df["seq"].str.upper())
    for pid, s in read_fasta(args.seed_fa):
        if len(s) >= args.min_len and s.upper() not in seen:
            species = "posadasii" if "82" in pid or "Silveira" in pid else "immitis"
            seed_rows.append(
                {
                    "protein": pid,
                    "strain": pid,
                    "length": len(s),
                    "species": species,
                    "source": "reference",
                    "seq": s.upper(),
                }
            )
            seen.add(s.upper())
    sdf = pd.DataFrame(seed_rows).assign(family="SOWgp", best_ident_sowgp=100.0)
    df = pd.concat([df.assign(source="pangenome"), sdf], ignore_index=True)

    df = df[(df["length"] >= args.min_len) & (df["seq"].str.len() >= args.min_len)].copy()
    df = df[~df["seq"].str.contains("X", na=False)]
    print(f"full-length copies retained: {len(df)}", file=sys.stderr)

    prof = pd.DataFrame([rp.best_repeat(s) for s in df["seq"]])
    df = pd.concat([df.reset_index(drop=True), prof.reset_index(drop=True)], axis=1)
    df["n_units"] = [len(su.anchor_units(s)) for s in df["seq"]]
    df["units_regular"] = [su.is_regular(su.anchor_units(s)) for s in df["seq"]]

    # ---- phase-anchored unit extraction ----
    unit_rows = []
    for _, row in df.iterrows():
        us = su.anchor_units(row["seq"])
        for i, ((p, inferred), u) in enumerate(zip(us, su.unit_seqs(row["seq"], us), strict=False)):
            unit_rows.append(
                {
                    "protein": row["protein"],
                    "strain": row["strain"],
                    "species": row["species"],
                    "source": row["source"],
                    "length": row["length"],
                    "unit_start": p,
                    "anchor_inferred": inferred,
                    "n_units": len(us),
                    "unit_idx": i,
                    "unit_class": su.unit_class(u),
                    "unit_seq": u,
                }
            )
    udf = pd.DataFrame(unit_rows)
    df.to_csv(args.prefix + ".copies.tsv", sep="\t", index=False)
    print(f"units: {len(udf)} from {udf['protein'].nunique()} copies", file=sys.stderr)

    internal = udf[(udf.unit_class != "terminal") & (udf.unit_seq.str.len() == su.PERIOD)]
    ref = internal.unit_seq.value_counts().index[0]
    print(f"reference unit (modal internal): {ref}", file=sys.stderr)
    aln = project_to_reference(list(udf["unit_seq"]), ref)
    aln_method = "projected onto modal unit"
    udf["aln"] = aln
    L = len(aln[0])
    cons_col = []
    for j in range(L):
        c = Counter(a[j] for a in aln if a[j] != "-")
        cons_col.append(c.most_common(1)[0][0] if c else "-")
    consensus = "".join(x for x in cons_col if x != "-")
    ids = []
    for a in aln:
        n = sum(1 for j in range(L) if a[j] != "-")
        s = sum(1 for j in range(L) if a[j] != "-" and a[j] == ref[j])
        ids.append(s / n if n else 0.0)
    udf["ref_identity"] = ids  # identity to the modal internal unit (ref)
    udf.to_csv(args.prefix + ".units.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.3, 0.55, 0.55], hspace=0.7, wspace=0.18)
    fig.suptitle(
        f"SOWgp repeat array variability (full-length copies >= {args.min_len} aa; "
        f"units anchored on {su.ANCHOR}; aligned: {aln_method})",
        fontsize=16,
        fontweight="bold",
    )

    # (a) architecture, collapsed to distinct sequences; units drawn at their positions
    ax = fig.add_subplot(gs[0, 0])
    hap = df.groupby("seq").agg(n=("strain", "size"), species=("species", "first")).reset_index()
    hap["n_units"] = [len(su.anchor_units(s)) for s in hap.seq]
    hap["L"] = hap.seq.str.len()
    hap = hap.sort_values(["species", "n_units", "L"]).reset_index(drop=True)
    for i, r in hap.iterrows():
        y = len(hap) - 1 - i
        us = su.anchor_units(r.seq)
        ax.plot([0, r.L], [y, y], color="#8a897f", lw=1, solid_capstyle="butt")
        for p, _ in us:
            ax.barh(
                y,
                su.PERIOD - 2,
                left=p + 1,
                height=1.0,
                color=VIVID[r.species],
                linewidth=0,
            )
    ax.set_yticks([])
    ax.set_xlabel("position in protein (aa)")
    panel_title(
        ax, f"(a) {len(hap)} distinct sequences ({len(df)} copies);\nblocks = anchored 47 aa units"
    )
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    # (b) integer unit count by species
    bx = fig.add_subplot(gs[0, 1])
    ks = range(int(df.n_units.min()), int(df.n_units.max()) + 1)
    w = 0.38
    for k, sp in enumerate(("immitis", "posadasii")):
        cnt = df[df.species == sp].n_units.value_counts()
        vals = [cnt.get(x, 0) for x in ks]
        xs = [x + (k - 0.5) * (w + 0.02) for x in ks]
        bx.bar(xs, vals, width=w, color=VIVID[sp], label=f"{sp} (n={sum(vals)})")
        for x, v in zip(xs, vals, strict=False):
            if v:
                bx.text(x, v + 0.8, str(v), ha="center", fontsize=7, color=su.MUTED)
    bx.set_xticks(list(ks))
    bx.set_xlabel("repeat units (PTDCYGDC anchor count)")
    bx.set_ylabel("copies (strains)")
    bx.legend(frameon=False, fontsize=13)
    bx.grid(axis="y", color=su.GRID, lw=0.6)
    bx.set_axisbelow(True)
    for s in ("top", "right"):
        bx.spines[s].set_visible(False)
    panel_title(bx, "(b) unit-count distribution")

    # (c) per-species information-content logos of internal (non-terminal) units
    for k, sp in enumerate(("immitis", "posadasii")):
        lx = fig.add_subplot(gs[1 + k, 0])
        m = (udf.species == sp) & (udf.unit_class != "terminal")
        sub = list(udf.loc[m, "aln"])
        logo(
            lx,
            sub,
            f"(c{k + 1}) {sp}: internal units (n={len(sub)}),\nposition within aligned unit",
        )

    # (d) raster of aligned units vs the modal internal unit
    dx = fig.add_subplot(gs[1:, 1])
    order = udf.sort_values(["species", "unit_class", "n_units", "protein", "unit_idx"]).index
    mat = np.zeros((len(order), L))
    for r, i in enumerate(order):
        a = aln[i]
        for j in range(L):
            mat[r, j] = 0.5 if a[j] == "-" else (1.0 if a[j] == ref[j] else 0.0)
    from matplotlib.colors import ListedColormap

    dx.imshow(
        mat,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(["#eb6834", "#e6e5df", "#256abf"]),
        vmin=0,
        vmax=1,
    )
    # species / class boundaries
    lab = udf.loc[order, ["species", "unit_class"]].agg(" ".join, axis=1).tolist()
    edges = [0] + [r for r in range(1, len(lab)) if lab[r] != lab[r - 1]] + [len(lab)]
    for a, b in zip(edges, edges[1:], strict=False):
        dx.axhline(a - 0.5, color="white", lw=1)
        if b - a >= 8:
            dx.text(
                L - 0.2,
                (a + b) / 2,
                f" {lab[a]} ({b - a})",
                va="center",
                fontsize=6.5,
                color=su.INK,
            )
    dx.set_yticks([])
    dx.set_xlabel("aligned unit position")
    panel_title(
        dx,
        "(d) aligned units vs modal internal unit:\nblue = same, orange = differs, grey = gap\n"
        "rows grouped by species and unit class; groups < 8 unlabelled",
    )

    plots = Path(args.prefix).parent / "plots"
    plots.mkdir(exist_ok=True)
    stem = plots / Path(args.prefix).name
    out = f"{stem}.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    print(f"wrote {out}", file=sys.stderr)
    print(f"column consensus (all units incl. terminal): {consensus}", file=sys.stderr)


if __name__ == "__main__":
    main()
