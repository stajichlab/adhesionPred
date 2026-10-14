#!/usr/bin/python3.12
"""Second proteome per target genus and the within-Aspergillus supplement (spec 5.3 and section 6.7).

second_per_genus: the next proteome by the select_genera.py rule among the proteomes not already chosen.
Supplement: up to N Aspergillus proteomes (one per species, passing the rule), chosen by farthest-point sampling on the patristic
distance of the Fungi_5k species tree, starting from the proteome already chosen for the genus. Clades: average-linkage
clustering of the patristic distances cut at CLADES groups (written in the output).
Usage: select_supplement.py --selection selection_discovery.tsv --out-second S.tsv --out-supp T.tsv
"""

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("select_genera", HERE / "select_genera.py")
sg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sg)


def second_per_genus(rows, first, known_taxa, merged=None, hybrids=frozenset()):
    chosen = {r["_name"] for r in first}
    rest = [r for r in rows if r["_name"] not in chosen]
    return sg.pick_per_genus(rest, known_taxa, merged, hybrids)


def farthest_points(names, dist, start, n):
    out = [start]
    mind = {x: dist(start, x) for x in names}
    while len(out) < n and len(out) < len(names):
        nxt = max((x for x in names if x not in out), key=lambda x: (mind[x], x))
        out.append(nxt)
        for x in names:
            mind[x] = min(mind[x], dist(nxt, x))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--out-second", required=True)
    ap.add_argument("--out-supp", required=True)
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--clades", type=int, default=8)
    a = ap.parse_args()
    F = sg.ROOT
    sel = list(csv.DictReader(open(a.selection), delimiter="\t"))
    rows = list(csv.DictReader(open(F / "samples.csv")))
    prefix = {}
    for line in open(F / "proteinfile2prefix.tsv"):
        name, tag = line.rstrip("\n").split("\t")[:2]
        prefix[tag] = name
    counts = {}
    for line in open(F / "protein_counts.tsv"):
        f = line.rstrip("\n").split(":")
        counts[f[0].replace(".proteins.fa", "")] = int(f[-1])
    for r in rows:
        name = prefix.get(r["LOCUSTAG"])
        r["_name"] = name
        r["_n"] = counts.get(name, 0) if name else 0
        if name and r["_n"] < 100 and (F / "input" / f"{name}.proteins.fa").exists():
            r["_n"] = sg.fasta_count(F / "input" / f"{name}.proteins.fa")
    rows = [r for r in rows if r["_name"] and (F / "input" / f"{r['_name']}.proteins.fa").exists()]
    known = {line.split("\t|\t")[0] for line in open("_workdir/sorting_hat/taxdump/nodes.dmp")}
    merged = {}
    for line in open("_workdir/sorting_hat/taxdump/merged.dmp"):
        f = line.split("\t|\t")
        merged[f[0]] = f[1].split("\t")[0]
    hybrids = {
        r["name"] for r in csv.DictReader(open(HERE / "hybrids.tsv"), delimiter="\t") if r["name"]
    }
    chosen = {r["name"] for r in sel}
    tgt_genera = {r["genus"] for r in sel if r["group"] in ("onygenales", "aspergillaceae")}
    first = [r for r in rows if r["_name"] in chosen and r["GENUS"] in tgt_genera]
    pool = [r for r in rows if r["GENUS"] in tgt_genera]
    second = second_per_genus(pool, first, known, merged, hybrids)
    group = {r["genus"]: r["group"] for r in sel}
    with open(a.out_second, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["group", "name", "fasta", "taxon", "order", "family", "genus", "asmid", "proteins"]
        )
        for r in second:
            w.writerow(
                [
                    group[r["GENUS"]],
                    r["_name"],
                    str(F / "input" / f"{r['_name']}.proteins.fa"),
                    r["NCBI_TAXONID"],
                    r["ORDER"],
                    r["FAMILY"],
                    r["GENUS"],
                    r["ASMID"],
                    r["_n"],
                ]
            )
    # supplement
    from Bio import Phylo
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    asp_all = [
        r
        for r in rows
        if r["GENUS"] == "Aspergillus"
        and r["_name"] not in hybrids
        and r["_n"] >= sg.MIN_PROTEINS
        and merged.get(r["NCBI_TAXONID"], r["NCBI_TAXONID"]) in known
    ]
    by_species = {}
    for r in sorted(asp_all, key=lambda r: (r["SPECIES"], r["_name"])):
        by_species.setdefault(r["SPECIES"], r)
    start = next(r["name"] for r in sel if r["genus"] == "Aspergillus")
    start_species = next(r["SPECIES"] for r in rows if r["_name"] == start)
    by_species[start_species] = next(r for r in rows if r["_name"] == start)
    names = [r["_name"] for r in by_species.values()]
    sys.setrecursionlimit(1_000_000)
    tree = Phylo.read(F / "Phylogeny/fungi_tree/final_tree.nw", "newick")
    tips = {t.name: t for t in tree.get_terminals() if t.name in set(names)}
    names = [n for n in names if n in tips]
    sub = Phylo.BaseTree.Tree(root=tree.common_ancestor(list(tips.values())))
    cache = {}

    def dist(x, y):
        if x == y:
            return 0.0
        k = (x, y) if x < y else (y, x)
        if k not in cache:
            cache[k] = sub.distance(tips[x], tips[y])
        return cache[k]

    pick = farthest_points(sorted(names), dist, start=start, n=a.n)
    import numpy as np

    mat = np.array([[dist(x, y) for y in pick] for x in pick])
    clade = fcluster(
        linkage(squareform(mat, checks=False), "average"), a.clades, criterion="maxclust"
    )
    cname = {r["_name"]: r for r in by_species.values()}
    with open(a.out_supp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "group",
                "name",
                "fasta",
                "taxon",
                "order",
                "family",
                "genus",
                "asmid",
                "proteins",
                "species",
                "clade",
            ]
        )
        for n, c in zip(pick, clade, strict=True):
            r = cname[n]
            w.writerow(
                [
                    "aspergillus_supplement",
                    n,
                    str(F / "input" / f"{n}.proteins.fa"),
                    merged.get(r["NCBI_TAXONID"], r["NCBI_TAXONID"]),
                    r["ORDER"],
                    r["FAMILY"],
                    r["GENUS"],
                    r["ASMID"],
                    r["_n"],
                    r["SPECIES"],
                    int(c),
                ]
            )
    print(
        len(second),
        "second proteomes;",
        len(pick),
        "supplement species;",
        len(set(clade)),
        "clades",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
