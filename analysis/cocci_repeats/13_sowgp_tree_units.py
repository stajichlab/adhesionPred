#!/usr/bin/python3.12
"""SOWgp repeat unit arrays mapped onto the Coccidioides strain phylogeny.

Question: did each unit count (3/4/5/6) arise once, or many times independently?

Inputs:
  --tree   IQ-TREE partitioned CDS ML tree, 488 taxa, with bootstrap support
           (PopGenomics/2025_All_Cocci/Phylogeny/.../Cocci_cds.488taxa_ascomycota.fa.part.aicc.contree)
  sowgp_repeat_viz.copies.tsv (09)  full-length copies (>= 250 aa), one per strain
  sowgp_pangenome.tsv (05)          all SOWgp-family gene models, incl. fragments

Outputs:
  sowgp_tree_units.png / .pdf   cladogram (nodes with support < --min-support collapsed),
                                species-by-name strip, per-strain mini unit map
  sowgp_tree_units.tsv          per-tip status, unit count, unit classes
  sowgp_tree_units.stats.tsv    Fitch-Hartigan parsimony changes per species vs a
                                within-species permutation null, for unit count and
                                (control) the flank haplotype outside the repeat array

The tree is midpoint-rooted. Within-species branches are very short and many
within-species nodes have low support, so the per-strain order within a species
is poorly resolved. The parsimony test uses the full (uncollapsed) ML tree; strains
without a full-length copy are wildcards (they add no cost).
"""

import argparse
import random
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import sowgp_units as su  # noqa: E402

TREE = Path(
    "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/"
    "Phylogeny/results/msa_filter_cds_ascomycota-buildtree/"
    "Cocci_cds.488taxa_ascomycota.fa.part.aicc.contree"
)
# gene-model strain names that differ from tree tip names
NAME_FIX = {
    "CimmitisRS_FungiDB": "CimmitisRS",
    "CimmitisWA211_FungiDB": "CimmitisWA211",
    "CposadasiiSilveira2022_FungiDB": "CposadasiiSilveira2022",
}


def name_species(n):
    return "posadasii" if ("posadasii" in n or n.startswith("Cpos")) else "immitis"


def fitch_changes(tree, state):
    """Hartigan/Fitch parsimony length for an unordered character. Missing tips = wildcard."""
    cost = 0

    def post(cl):
        nonlocal cost
        if cl.is_terminal():
            s = state.get(cl.name)
            return None if s is None else {s}
        sets = [x for x in (post(c) for c in cl.clades) if x is not None]
        if not sets:
            return None
        cnt = {}
        for s in sets:
            for v in s:
                cnt[v] = cnt.get(v, 0) + 1
        m = max(cnt.values())
        cost += len(sets) - m
        return {v for v, c in cnt.items() if c == m}

    post(tree.root)
    return cost


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tree", default=str(TREE))
    ap.add_argument("--min-support", type=float, default=70)
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    # ---- per-strain status ----
    cp = pd.read_csv(HERE / "sowgp_repeat_viz.copies.tsv", sep="\t")
    cp["tip"] = cp.strain.str.replace(".proteins", "", regex=False).replace(NAME_FIX)
    pg = pd.read_csv(HERE / "sowgp_pangenome.tsv", sep="\t")
    pg = pg[pg.family == "SOWgp"]
    pg_tips = set(pg.strain.str.replace(".proteins", "", regex=False).replace(NAME_FIX))
    full = {r.tip: r.seq for r in cp.itertuples()}

    tree = Phylo.read(args.tree, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    tips = [c.name for c in tree.get_terminals()]
    missing = sorted(set(full) - set(tips))
    if missing:
        print(f"full-length copies not in tree (not plotted): {missing}", file=sys.stderr)

    # clade species = majority species-by-name of the root child containing the tip
    clade_sp = {}
    for k in tree.root.clades:
        ts = [c.name for c in k.get_terminals()]
        maj = pd.Series([name_species(x) for x in ts]).mode()[0]
        clade_sp.update({x: maj for x in ts})

    rows = []
    for t in tips:
        seq = full.get(t)
        if seq is not None:
            us = su.anchor_units(seq)
            classes = [su.unit_class(u) for u in su.unit_seqs(seq, us)]
            status, n_units = "full-length", len(us)
        else:
            classes, us = [], []
            status = "fragment only" if t in pg_tips else "no SOWgp gene model"
            n_units = np.nan
        rows.append(
            {
                "tip": t,
                "name_species": name_species(t),
                "clade_species": clade_sp[t],
                "status": status,
                "n_units": n_units,
                "classes": ",".join(classes),
                "unit_starts": ",".join(str(p) for p, _ in us),
                "length": len(seq) if seq else np.nan,
            }
        )
    D = pd.DataFrame(rows)
    D.to_csv(HERE / "sowgp_tree_units.tsv", sep="\t", index=False)
    conflict = D[D.name_species != D.clade_species]
    print("name/clade species conflicts:", list(conflict.tip), file=sys.stderr)
    print(D.status.value_counts().to_string(), file=sys.stderr)

    # ---- parsimony test (full ML tree) ----
    # trait 'n_units' is the question; trait 'flank_haplotype' (protein sequence outside
    # the repeat array) is a control: if flanks cluster on the tree but unit count does
    # not, the tree resolves this locus and unit count changes recurrently.
    def flank(seq):
        us = su.anchor_units(seq)
        end = us[-1][0] + len(su.unit_seqs(seq, us)[-1])
        return seq[: us[0][0]] + "|" + seq[end:]

    flank_of = {t: flank(s) for t, s in full.items()}
    D["flank_haplotype"] = D.tip.map(
        lambda t: pd.factorize(pd.Series(list(flank_of.values())))[1].get_loc(flank_of[t])
        if t in flank_of
        else np.nan
    )
    D.to_csv(HERE / "sowgp_tree_units.tsv", sep="\t", index=False)
    rng = random.Random(args.seed)
    stats = []
    for sp in su.SPECIES_INK:
        sub = D[(D.clade_species == sp) & (D.status == "full-length")]
        for trait in ("n_units", "flank_haplotype"):
            state = dict(zip(sub.tip, sub[trait].astype(int), strict=False))
            obs = fitch_changes(tree, state)
            vals = list(state.values())
            null = []
            for _ in range(args.perm):
                rng.shuffle(vals)
                null.append(fitch_changes(tree, dict(zip(state.keys(), vals, strict=False))))
            null = np.array(null)
            p = (np.sum(null <= obs) + 1) / (len(null) + 1)
            stats.append(
                {
                    "species": sp,
                    "trait": trait,
                    "n_tips": len(state),
                    "n_states": len(set(vals)),
                    "state_counts": dict(pd.Series(vals).value_counts().sort_index()),
                    "min_possible": len(set(vals)) - 1,
                    "observed_changes": obs,
                    "null_mean": round(null.mean(), 1),
                    "null_5pct": np.percentile(null, 5),
                    "p_le_obs": round(p, 4),
                    "perms": args.perm,
                }
            )
    S = pd.DataFrame(stats)
    S.to_csv(HERE / "sowgp_tree_units.stats.tsv", sep="\t", index=False)
    print(S.drop(columns=["state_counts"]).to_string(index=False), file=sys.stderr)

    # ---- collapsed cladogram for display ----
    ct = Phylo.read(args.tree, "newick")
    ct.root_at_midpoint()
    for c in list(ct.get_nonterminals()):
        if c is not ct.root and c.confidence is not None and c.confidence < args.min_support:
            ct.collapse(c)
    ct.ladderize()
    ctips = [c.name for c in ct.get_terminals()]
    y = {n: i for i, n in enumerate(ctips)}
    depth = ct.depths(unit_branch_lengths=True)
    maxd = max(depth.values())
    Dd = D.set_index("tip")

    def ypos(cl):
        if cl.is_terminal():
            return y[cl.name]
        return np.mean([ypos(c) for c in cl.clades])

    N = len(ctips)
    fig_h = max(12, N * 0.052)
    fig, (tx, sx, mx) = plt.subplots(
        1,
        3,
        figsize=(13, fig_h),
        sharey=True,
        gridspec_kw={"width_ratios": [1.3, 0.08, 2.2], "wspace": 0.02},
    )
    ycache = {}

    def draw(cl):
        yc = ycache.setdefault(id(cl), ypos(cl))
        x = depth[cl]
        if not cl.is_terminal():
            ys = [draw(c) for c in cl.clades]
            tx.plot([x, x], [min(ys), max(ys)], color=su.MUTED, lw=0.5)
            for c, yy in zip(cl.clades, ys, strict=False):
                tx.plot(
                    [x, depth[c] if not c.is_terminal() else maxd],
                    [yy, yy],
                    color=su.MUTED if not c.is_terminal() else su.GRID,
                    lw=0.5,
                )
        return yc

    draw(ct.root)
    # support labels on the root children only (the species split)
    for c in ct.root.clades:
        if c.confidence is not None:
            tx.text(
                depth[c],
                ycache[id(c)],
                f"{c.confidence:.0f}",
                fontsize=6,
                ha="right",
                va="bottom",
                color=su.INK,
            )
    tx.set_xlim(-0.5, maxd + 0.3)
    tx.axis("off")
    tx.set_title(
        f"CDS ML tree, midpoint-rooted;\nnodes with support < " f"{args.min_support:.0f} collapsed",
        fontsize=8,
        loc="left",
    )

    # species-by-name strip; conflicts marked
    for t in ctips:
        r = Dd.loc[t]
        sx.barh(y[t], 1, color=su.SPECIES_INK[r.name_species], height=1.0, lw=0)
        if r.name_species != r.clade_species:
            mx.text(
                6.3,
                y[t],
                f"◀ name says {r.name_species}, tree places it in " f"{r.clade_species}: {t}",
                fontsize=6,
                va="center",
                color=su.INK,
            )
    sx.set_xlim(0, 1)
    sx.axis("off")
    sx.set_title("name", fontsize=7)

    # mini unit maps
    for t in ctips:
        r = Dd.loc[t]
        yy = y[t]
        if r.status == "full-length":
            classes = r.classes.split(",")
            for i, c in enumerate(classes):
                mx.barh(
                    yy,
                    0.9,
                    left=i,
                    height=0.8,
                    color=su.CLASS_COLOR.get(c, su.OTHER_COLOR),
                    edgecolor=su.INK if c not in su.CLASS_COLOR else "none",
                    lw=0.3,
                )
        elif r.status == "fragment only":
            mx.plot([0, 2.9], [yy, yy], color="#c9c8c1", lw=0.8)
        # 'no SOWgp gene model': left blank
    mx.set_xlim(-0.2, 8.5)
    mx.set_ylim(N - 0.5, -0.5)
    mx.set_xticks(range(0, 6))
    mx.set_xticklabels([f"u{i + 1}" for i in range(6)], fontsize=7)
    mx.xaxis.tick_top()
    for s in ("right", "bottom", "left"):
        mx.spines[s].set_visible(False)
    mx.tick_params(left=False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in su.CLASS_COLOR.values()]
    handles.append(plt.Line2D([0], [0], color="#c9c8c1", lw=2))
    labels = [f"unit {k}" if k != "terminal" else "terminal unit" for k in su.CLASS_COLOR]
    labels.append("fragment only (< 250 aa); blank = no SOWgp gene model")
    mx.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0, -0.002),
        ncol=3,
        fontsize=7,
        frameon=False,
    )
    stat_txt = "\n".join(
        f"{r.species} {r.trait}: {r.observed_changes} changes on ML tree, "
        f"permutation null mean {r.null_mean} (P = {r.p_le_obs}; "
        f"{r.n_tips} strains, {r.n_states} states)"
        for r in S.itertuples()
    )
    mx.set_title(
        "Repeat units per strain (full-length copy, in array order)\n" + stat_txt,
        fontsize=7.5,
        loc="left",
    )

    fig.savefig(HERE / "sowgp_tree_units.png", dpi=150, bbox_inches="tight")
    fig.savefig(HERE / "sowgp_tree_units.pdf", bbox_inches="tight")
    print("wrote sowgp_tree_units.png/.pdf/.tsv/.stats.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
