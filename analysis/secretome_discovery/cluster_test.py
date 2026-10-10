#!/usr/bin/python3.12
"""Pfam-free cluster test of the secretome discovery pass (spec 5.2).

Units are MMseqs2 clusters of results/Fungi5K_cluster.tsv.gz (representative and member protein IDs, no size filter).
A cluster is "Pfam-free" if no member in the chosen proteomes has a Pfam-A hit (Fungi_5k results/function/pfam). Presence of a genus
= at least one pool protein of the genus's chosen proteome in the cluster. The file is read once, grouped by representative.
Usage: cluster_test.py --selection S --pools DIR --out DIR [--cluster-file ...]
"""

import argparse
import csv
import gzip
import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stats  # noqa: E402

RES = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/results")
TARGETS = ("onygenales", "aspergillaceae")


def scan_clusters(path, tag_genus, pool, pfam):
    """{pool name: {cluster rep: {genus: n pool proteins}}} for Pfam-free clusters with at least one pool protein."""
    out = {p: {} for p in pool}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        rows = (line.rstrip("\n").split("\t") for line in fh)
        for rep, grp in itertools.groupby(rows, key=lambda r: r[0]):
            members = [r[1] for r in grp]
            chosen = [x for x in members if x.split("_", 1)[0] in tag_genus]
            if not chosen or any(x in pfam for x in chosen):
                continue
            for p, prot in pool.items():
                hit = [x for x in chosen if x in prot]
                if hit:
                    d = {}
                    for x in hit:
                        g = tag_genus[x.split("_", 1)[0]]
                        d[g] = d.get(g, 0) + 1
                    out[p][rep] = d
    return out


def test_clusters(clusters, genera, pool_size, min_genera=3):
    """clusters: {id: {genus: count}}; genera: {genus: 'target'|'background'}; returns stats.report-like dict."""
    from scipy import stats as sst

    nt = sum(1 for g in genera.values() if g == "target")
    nb = len(genera) - nt
    ps = {}
    for cid, d in clusters.items():
        tp = sum(1 for g in d if genera.get(g) == "target")
        bp = sum(1 for g in d if genera.get(g) == "background")
        _, p = sst.fisher_exact([[tp, nt - tp], [bp, nb - bp]], alternative="two-sided")
        ps[cid] = (tp, bp, p)
    ids = sorted(ps)
    q = dict(zip(ids, stats.bh([ps[i][2] for i in ids]), strict=True)) if ids else {}
    enriched, depleted = [], []
    for i in ids:
        tp, bp, p = ps[i]
        if tp / nt > bp / nb and tp >= min_genera:
            rows = [
                {"genus": g, "group": grp, "pool": pool_size[g], "hits": clusters[i].get(g, 0)}
                for g, grp in genera.items()
            ]
            s = stats.family_stats(rows)
            s.update(family=i, fisher_q=q[i])
            enriched.append(s)
        elif tp / nt < bp / nb:
            depleted.append(
                {"family": i, "t_present": tp, "b_present": bp, "fisher_p": p, "fisher_q": q[i]}
            )
    enriched = [s for s in enriched if s["fisher_q"] < 0.05]
    depleted = [s for s in depleted if s["fisher_q"] < 0.05]
    return {
        "enriched": sorted(enriched, key=lambda s: s["fisher_p"]),
        "depleted": depleted,
        "n_tested": len(ids),
    }


def load_pool_members(selection, pools_dir):
    pool = {"S": set(), "U": set()}
    for r in selection:
        with gzip.open(Path(pools_dir) / f"{r['name']}.pools.tsv.gz", "rt") as fh:
            for p in csv.DictReader(fh, delimiter="\t"):
                pool["S"].add(p["protein"])
                if p["pool"] == "U":
                    pool["U"].add(p["protein"])
    return pool


def load_pfam_hit_proteins(selection):
    hit = set()
    for r in selection:
        with gzip.open(RES / "function" / "pfam" / f"{r['name']}.domtblout.gz", "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    hit.add(line.split(None, 1)[0])
    return hit


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--pools", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--prefix", default="/bigdata/stajichlab/shared/projects/Fungi_5k/proteinfile2prefix.tsv"
    )
    ap.add_argument("--samples", default="/bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv")
    ap.add_argument("--cluster-file", default=str(RES / "Fungi5K_cluster.tsv.gz"))
    a = ap.parse_args()
    sel = list(csv.DictReader(open(a.selection), delimiter="\t"))
    name2tag = {}
    for line in open(a.prefix):
        name, tag = line.rstrip("\n").split("\t")[:2]
        name2tag[name] = tag
    tag_genus = {name2tag[r["name"]]: r["genus"] for r in sel}
    pool = load_pool_members(sel, a.pools)
    pfam = load_pfam_hit_proteins(sel)
    print(
        "pool S",
        len(pool["S"]),
        "pool U",
        len(pool["U"]),
        "proteins with a Pfam hit",
        len(pfam),
        flush=True,
    )
    found = scan_clusters(a.cluster_file, tag_genus, pool, pfam)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    psize = {}
    for r in sel:
        for pname in pool:
            with gzip.open(Path(a.pools) / f"{r['name']}.pools.tsv.gz", "rt") as fh:
                n = sum(
                    1
                    for p in csv.DictReader(fh, delimiter="\t")
                    if pname == "S" or p["pool"] == "U"
                )
            psize[(pname, r["genus"])] = n
    group = {r["genus"]: r["group"] for r in sel}
    for pname in pool:
        for target in TARGETS:
            genera = {
                g: ("target" if grp == target else "background")
                for g, grp in group.items()
                if grp in (target, "background")
            }
            ps = {g: psize[(pname, g)] for g in genera}
            clusters = {
                c: {g: n for g, n in d.items() if g in genera} for c, d in found[pname].items()
            }
            clusters = {c: d for c, d in clusters.items() if d}
            res = test_clusters(clusters, genera, ps)
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
            for kind in ("enriched", "depleted"):
                with open(out / f"clusters_{pname}_{target}_{kind}.tsv", "w", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
                    w.writeheader()
                    w.writerows(res[kind])
            print(
                pname,
                target,
                "Pfam-free clusters tested",
                res["n_tested"],
                "enriched q<0.05",
                len(res["enriched"]),
                flush=True,
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
