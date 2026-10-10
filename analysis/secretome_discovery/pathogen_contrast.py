#!/usr/bin/python3.12
"""Descriptive contrast inside Onygenales (spec 5.3): genera with mammal-pathogenic species against the other Onygenales genera.

No test is run (27 genera). For each family given (default: the Onygenales enriched families of pool U), the number of genera with the
family in the strict pathogen group (mammal_pathogen), the broad group (adds mammal_pathogen_mixed) and in all other Onygenales genera.
Usage: pathogen_contrast.py --selection S --pools DIR --groups onygenales_host_groups.tsv --families families_U_onygenales.tsv --out F
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path


def genus_families(path, pool):
    fams = set()
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if pool == "U" and r["pool"] != "U":
                continue
            fams.update(filter(None, r["pfam"].split(";")))
    return fams


def contrast(presence, groups, families):
    """presence: {genus: set of families}; groups: {genus: group}; returns rows of counts."""
    strict = {g for g, x in groups.items() if x == "mammal_pathogen"}
    broad = strict | {g for g, x in groups.items() if x == "mammal_pathogen_mixed"}
    rows = []
    for f in families:
        s = sum(1 for g in strict if f in presence.get(g, ()))
        b = sum(1 for g in broad if f in presence.get(g, ()))
        others = [g for g in presence if g not in strict]
        o = sum(1 for g in others if f in presence[g])
        ob = [g for g in presence if g not in broad]
        o2 = sum(1 for g in ob if f in presence[g])
        rows.append(
            {
                "family": f,
                "strict_present": s,
                "strict_n": len(strict),
                "others_present": o,
                "others_n": len(others),
                "broad_present": b,
                "broad_n": len(broad),
                "others_broad_present": o2,
                "others_broad_n": len(ob),
            }
        )
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--pools", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pool", default="U")
    a = ap.parse_args()
    sel = [
        r for r in csv.DictReader(open(a.selection), delimiter="\t") if r["group"] == "onygenales"
    ]
    groups = {r["genus"]: r["group"] for r in csv.DictReader(open(a.groups), delimiter="\t")}
    missing = [r["genus"] for r in sel if r["genus"] not in groups]
    if missing:
        sys.exit(f"genera without a group: {missing}")
    presence = {
        r["genus"]: genus_families(Path(a.pools) / f"{r['name']}.pools.tsv.gz", a.pool) for r in sel
    }
    fams = [r["family"] for r in csv.DictReader(open(a.families), delimiter="\t")]
    rows = contrast(presence, {r["genus"]: groups[r["genus"]] for r in sel}, fams)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(
        len(rows),
        "families; strict group",
        rows[0]["strict_n"],
        "genera; broad",
        rows[0]["broad_n"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
