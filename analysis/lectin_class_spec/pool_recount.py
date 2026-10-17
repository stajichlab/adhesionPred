#!/usr/bin/env python3.12
"""Recount lectin-like Pfam families in the discovery pools (fact F5 of the lectin class spec).

Usage: pool_recount.py POOLS_DIR SELECTION_TSV OUT_TSV
POOLS_DIR      the local directory _workdir/discovery/pools (653 files <name>.pools.tsv.gz; not in git)
SELECTION_TSV  analysis/fungi_scan/selection_discovery.tsv (columns group, name, ...; one proteome per genus)
Each pool file has one row per pool S protein. Column `pool` is S or U (U is a subset of S).
Column `pfam` is a ';' list of Pfam accessions (Pfam-A release 37.x, unversioned).
Output: one row per family with proteins in pool S and pool U, and the number of genera
(= proteomes) that have at least one pool S protein, split by group.
"""

import collections
import csv
import gzip
import sys

pools_dir, sel_tsv, out_tsv = sys.argv[1:4]
FAMILIES = {
    "PF01419": "Jacalin",
    "PF00652": "Ricin_B_lectin",
    "PF18647": "Fungal_lectin_2",
    "PF14200": "RicinB_lectin_2",
    "PF13385": "Laminin_G_3",
    "PF22803": "GBD_Y3",
    "PF07938": "Fungal_lectin",
    "PF21691": "LDL",
    "PF08881": "CVNH",
    "PF01476": "LysM",
    "PF03388": "Lectin_leg-like",
    "PF00262": "Calreticulin",
    "PF00337": "Gal-bind_lectin",
    "PF07367": "FB_lectin",
    "PF21595": "CCL2-like",
    "PF11721": "Malectin",
    "PF07691": "PA14",
    "PF10528": "GLEYA",
}
sel = list(csv.DictReader(open(sel_tsv), delimiter="\t"))
n_s = collections.Counter()
n_u = collections.Counter()
genera = collections.defaultdict(lambda: collections.Counter())
tot_s = tot_u = 0
for r in sel:
    path = f"{pools_dir}/{r['name']}.pools.tsv.gz"
    seen = set()
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tot_s += 1
            is_u = row["pool"] == "U"
            tot_u += is_u
            fams = {x for x in row["pfam"].split(";") if x}
            for f in fams & FAMILIES.keys():
                n_s[f] += 1
                n_u[f] += is_u
                seen.add(f)
    for f in seen:
        genera[f][r["group"]] += 1
groups = ["onygenales", "aspergillaceae", "background"]
with open(out_tsv, "w") as o:
    o.write(
        "pfam\tname\tpool_S\tpool_U\t"
        + "\t".join(f"genera_{g}" for g in groups)
        + "\tgenera_total\n"
    )
    for f, name in FAMILIES.items():
        g = [genera[f][x] for x in groups]
        o.write(f"{f}\t{name}\t{n_s[f]}\t{n_u[f]}\t" + "\t".join(map(str, g)) + f"\t{sum(g)}\n")
print(f"proteomes {len(sel)}, pool S proteins {tot_s}, pool U proteins {tot_u}")
