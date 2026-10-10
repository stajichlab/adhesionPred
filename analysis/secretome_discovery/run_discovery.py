#!/usr/bin/python3.12
"""Driver of the secretome discovery pass (B4): genus tables, family tests per pool and target, gains and losses on the tree.

Reads selection_discovery.tsv and <pools>/<name>.pools.tsv.gz (pools.py). Writes, per pool (S, U) and target (onygenales,
aspergillaceae): <out>/<pool>_<target>_families.tsv (all families tested, with q) and <out>/<pool>_<target>_depleted.tsv.
Usage: run_discovery.py --selection S.tsv --pools DIR --tree final_tree.nw --out DIR
"""

import argparse
import csv
import gzip
import io
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stats  # noqa: E402

TARGETS = ("onygenales", "aspergillaceae")


def genus_counts(path, pool):
    """(pool size, {family: proteins with the family}) for one proteome; pool 'S' = all rows, 'U' = rows with pool U."""
    size, fam = 0, Counter()
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if pool == "U" and r["pool"] != "U":
                continue
            size += 1
            for f in filter(None, r["pfam"].split(";")):
                fam[f] += 1
    return size, dict(fam)


def family_rows(genus):
    """genus: {name: (group, pool size, {family: hits})} -> {family: [row per genus]}."""
    fams = sorted({f for _, _, h in genus.values() for f in h})
    return {
        f: [
            {"genus": g, "group": grp, "pool": size, "hits": hits.get(f, 0)}
            for g, (grp, size, hits) in sorted(genus.items())
        ]
        for f in fams
    }


def prune_newick(newick, keep):
    from Bio import Phylo

    tree = Phylo.read(io.StringIO(newick), "newick")
    sys.setrecursionlimit(1_000_000)

    def walk(node):
        if node.is_terminal():
            return node if node.name in keep else None
        kids = [k for k in (walk(c) for c in node.clades) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kid = kids[0]
            kid.branch_length = (kid.branch_length or 0) + (node.branch_length or 0)
            return kid
        node.clades = kids
        return node

    root = walk(tree.root)
    tree.root = root
    out = io.StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--pools", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-genera", type=int, default=3)
    a = ap.parse_args()
    sel = list(csv.DictReader(open(a.selection), delimiter="\t"))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for pool in ("S", "U"):
        data = {}
        for r in sel:
            f = Path(a.pools) / f"{r['name']}.pools.tsv.gz"
            if not f.exists():
                continue
            size, fam = genus_counts(f, pool)
            data[r["genus"]] = (r["group"], size, fam)
        for target in TARGETS:
            sub = {
                g: (("target" if grp == target else "background"), s, h)
                for g, (grp, s, h) in data.items()
                if grp in (target, "background")
            }
            res = stats.report(family_rows(sub), min_genera=a.min_genera)
            for kind in ("enriched", "depleted"):
                rows = res[kind]
                p = out / f"{pool}_{target}_{kind}.tsv"
                with open(p, "w", newline="") as fh:
                    cols = [
                        "family",
                        "t_present",
                        "n_target",
                        "b_present",
                        "n_background",
                        "odds_ratio",
                        "fisher_p",
                        "fisher_q",
                        "logit_p",
                        "t_copies_median",
                        "b_copies_median",
                        "mwu_p",
                    ]
                    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
                    w.writeheader()
                    w.writerows(rows)
            print(
                pool,
                target,
                "tested",
                res["n_tested"],
                "enriched",
                len(res["enriched"]),
                "depleted",
                len(res["depleted"]),
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
