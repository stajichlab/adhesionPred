#!/usr/bin/env python3
"""Tables and figure for the PF28404 family tree.

Inputs (work directory): proteomes.tsv, hits.tsv, tips.tsv, tree/pf28404.treefile (IQ-TREE,
UFBoot support on inner nodes). The tree is rooted at the midpoint. Outputs go to --out-dir:
  copy_number_onygenales.tsv   copies per Onygenales proteome (all hits, any length)
  coccidioides_clades.tsv      per Coccidioides RS paralog: closest clade holding a non-Coccidioides tip
  nearest_anchor.tsv           each non-Coccidioides Onygenales tip: nearest RS paralog by patristic distance
  pf28404_tree.png/.pdf        the tree
The nearest-anchor assignment is a distance heuristic, not an orthology inference.
"""

import argparse
import collections
import csv
from pathlib import Path

from Bio import Phylo

ANCHORS = {"A_PRA3": "CIMG_02492", "B_07303": "CIMG_07303", "C_05560": "CIMG_05560",
           "D_07843": "CIMG_07843"}  # fmt: skip
COCCI_PREFIX = ("Cimmitis", "Coccidioides", "Cposadasii", "CposadasiiSilveira")


def tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--unit-genomes", help="cocci_repeats unit_genomes.tsv (genus column)")
    a = ap.parse_args()
    work, out = Path(a.work), Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    prot = [p for p in tsv(work / "proteomes.tsv") if p["order"] == "Onygenales"]
    hits = tsv(work / "hits.tsv")
    tips = {r["tip"]: r for r in tsv(work / "tips.tsv")}
    genus_of = {}
    if a.unit_genomes:
        genus_of = {r["label"]: r["genus"] for r in tsv(a.unit_genomes) if r["genus"]}

    def genus(label):
        if label in genus_of:
            return genus_of[label]
        g = label.split("_")[0]
        return "Coccidioides" if g.startswith(("Cimmitis", "Cposadasii", "CposadasiiS")) else g

    n = collections.Counter(h["label"] for h in hits if h["order"] == "Onygenales")
    with open(out / "copy_number_onygenales.tsv", "w") as fh:
        fh.write("proteome\tgenus\tpf28404_proteins\n")
        for p in sorted(prot, key=lambda p: (genus(p["label"]), p["label"])):
            fh.write(f"{p['label']}\t{genus(p['label'])}\t{n.get(p['label'], 0)}\n")

    tree = Phylo.read(work / "tree" / "pf28404.treefile", "newick")
    tree.root_at_midpoint()
    iscoc = lambda t: tips[t]["label"].startswith(COCCI_PREFIX)  # noqa: E731
    anchor_tip = {
        k: next(t for t in tips if t.endswith("|" + v + "-t26_1-p1")) for k, v in ANCHORS.items()
    }
    with open(out / "coccidioides_clades.tsv", "w") as fh:
        fh.write("anchor\tclade_tips\tcoccidioides_tips\tufboot\tnon_coccidioides_tips\n")
        for k, t in anchor_tip.items():
            for c in reversed(tree.get_path(t)[:-1]):
                ts = [x.name for x in c.get_terminals()]
                non = [x for x in ts if not iscoc(x)]
                if non:
                    names = ";".join(tips[x]["label"] for x in non)
                    fh.write(f"{k}\t{len(ts)}\t{len(ts) - len(non)}\t{c.confidence}\t{names}\n")
                    break
    with open(out / "nearest_anchor.tsv", "w") as fh:
        fh.write("proteome\tprotein_id\tnearest_anchor\tdistance\tgap_to_second\n")
        for t, r in tips.items():
            if r["order"] != "Onygenales" or iscoc(t):
                continue
            d = {k: tree.distance(t, a_) for k, a_ in anchor_tip.items()}
            ds = sorted(d.values())
            k = min(d, key=d.get)
            fh.write(f"{r['label']}\t{r['protein_id']}\t{k}\t{d[k]:.3f}\t{ds[1] - ds[0]:.3f}\n")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def colour(t):
        r = tips[t]
        if r["order"] != "Onygenales":
            return "#8a8a8a"
        return "#c0392b" if iscoc(t) else "#1f5fa8"

    fig, ax = plt.subplots(figsize=(9, 24))
    for cl in tree.get_terminals():
        cl.name = cl.name  # keep names for drawing
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: "", show_confidence=False)
    depths = tree.depths()
    ys = {cl: i + 1 for i, cl in enumerate(tree.get_terminals())}
    for cl in tree.get_terminals():
        lab = tips[cl.name]["label"].replace("_", " ")[:34]
        ax.text(
            depths[cl] * 1.005, ys[cl], f" {lab}", fontsize=3.6, va="center", color=colour(cl.name)
        )
    ax.set_title(
        "PF28404 family (217 proteins), IQ-TREE Q.pfam+I+R5, midpoint root\n"
        "red: Coccidioides; blue: other Onygenales; grey: outgroup orders",
        fontsize=8,
    )
    ax.set_xlabel("substitutions per site")
    ax.set_ylabel("")
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(out / "pf28404_tree.png", dpi=170)
    fig.savefig(out / "pf28404_tree.pdf")
    print("wrote", ", ".join(sorted(p.name for p in out.iterdir())))


if __name__ == "__main__":
    main()
