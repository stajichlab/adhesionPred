#!/usr/bin/python3.12
"""Class breakdown per proteome from the sorting hat calls (calls.long.tsv.gz of each run).

One row per proteome: proteins, and the number called for each call. The three `other_*` and
`cell_wall_adhesion_candidate` rows partition the proteins that R0 can assess (see the report of a
run). Evidence calls overlap. `wall_family_domain = not_called` means "no hit in an ACTIVE family".
Only PA14 is active (2026-10-07), so it does not mean "no wall domain".

Usage: class_breakdown.py --work WORKDIR --out FILE PROTEOME [PROTEOME ...]
"""

import argparse
import collections
import csv
import gzip
import sys

COLUMNS = [
    ("signal_peptide_protein", "R0", "signal peptide (R0)"),
    ("cell_wall_adhesion_candidate", "R0", "cell wall adhesion candidate"),
    ("other_surface_no_mechanism", "R0", "surface, no mechanism evidence"),
    ("other_not_surface", "R0", "not surface (no signal peptide)"),
    ("tandem_repeat_protein", "", "tandem repeat evidence"),
    ("wall_family_domain", "", "wall family domain (active families)"),
    ("iuis_allergen_similarity", "", "allergen similarity (35%, 80 aa)"),
    ("iuis_allergen_homolog", "", "allergen homolog (70%, 80% cov)"),
    ("serodiagnostic_marker_candidate", "R0", "serodiagnostic marker candidate"),
    ("cocci_specificity_rank_top15", "", "antigen rank top 15% (C. immitis RS only)"),
]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("proteomes", nargs="+")
    a = ap.parse_args()
    out = []
    for p in a.proteomes:
        c = collections.defaultdict(collections.Counter)
        proteins = set()
        with gzip.open(f"{a.work}/{p}/out/calls.long.tsv.gz", "rt") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                proteins.add(r["protein"])
                c[(r["call"], r["variant"])][r["value"]] += 1
        row = {"proteome": p, "proteins": len(proteins)}
        for call, var, label in COLUMNS:
            k = c[(call, var)]
            row[label] = f"{k['called']}" + (
                f" ({k['not_assessable']} n/a)" if k["not_assessable"] else ""
            )
        out.append(row)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    for r in out:
        print(r["proteome"], r)
    return 0


if __name__ == "__main__":
    sys.exit(main())
