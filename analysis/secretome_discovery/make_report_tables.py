#!/usr/bin/python3.12
"""Report tables of the discovery pass (B5): families with q < 0.05 with Pfam names, stability and tree counts; candidate protein lists.

Usage: make_report_tables.py --results DIR --pools DIR --selection S --names pfam_names.tsv --out DIR
The role column is "no stated role" for every family (the spec requires roles to be written before the statistics are seen;
that did not happen, so no roles are assigned here). The Pfam description is shown instead.
"""

import argparse
import csv
import gzip
import shutil
import sys
from pathlib import Path


def load_names(path):
    names = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        names[f[0]] = (f[1], f[2] if len(f) > 2 else "")
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--pools", required=True)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--names", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res, out = Path(a.results), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    names = load_names(a.names)
    sel = list(csv.DictReader(open(a.selection), delimiter="\t"))
    for pool in ("S", "U"):
        for target in ("onygenales", "aspergillaceae"):
            stab = {}
            sp = res / f"stability_{pool}_{target}.tsv"
            for r in csv.DictReader(open(sp), delimiter="\t"):
                stab[r["family"]] = r["n_seeds_q_lt_0.05"]
            rows = [
                r
                for r in csv.DictReader(open(res / f"{pool}_{target}_enriched.tsv"), delimiter="\t")
                if float(r["fisher_q"]) < 0.05
            ]
            cols = [
                "family",
                "name",
                "description",
                "role",
                "t_present",
                "n_target",
                "b_present",
                "n_background",
                "odds_ratio",
                "fisher_q",
                "logit_p",
                "t_copies_median",
                "b_copies_median",
                "mwu_p",
                "gains_in_target",
                "gains_elsewhere",
                "losses",
                "seeds_q_lt_0.05_of_5",
            ]
            with open(out / f"families_{pool}_{target}.tsv", "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(cols)
                for r in rows:
                    n = names.get(r["family"], ("?", ""))
                    w.writerow(
                        [
                            r["family"],
                            n[0],
                            n[1],
                            "no stated role",
                            r["t_present"],
                            r["n_target"],
                            r["b_present"],
                            r["n_background"],
                            r["odds_ratio"],
                            r["fisher_q"],
                            r["logit_p"],
                            r["t_copies_median"],
                            r["b_copies_median"],
                            r["mwu_p"],
                            r["gains_in_target"],
                            r["gains_elsewhere"],
                            r["losses"],
                            stab.get(r["family"], "0"),
                        ]
                    )
            if pool == "U" and target == "onygenales":
                fams = {r["family"] for r in rows}
            if pool == "U":
                shutil.copy(sp, out / sp.name)
    tg = [r for r in sel if r["group"] == "onygenales"]
    with gzip.open(out / "candidates_U_onygenales.tsv.gz", "wt") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "family",
                "name",
                "proteome",
                "genus",
                "protein",
                "length",
                "pct_cys",
                "pct_pro",
                "pct_st",
                "sequons_per_100aa",
                "gpi",
                "all_pfam",
            ]
        )
        for r in tg:
            with gzip.open(Path(a.pools) / f"{r['name']}.pools.tsv.gz", "rt") as ph:
                for p in csv.DictReader(ph, delimiter="\t"):
                    if p["pool"] != "U":
                        continue
                    pf = set(filter(None, p["pfam"].split(";")))
                    for f in sorted(pf & fams):
                        w.writerow(
                            [
                                f,
                                names.get(f, ("?",))[0],
                                r["name"],
                                r["genus"],
                                p["protein"],
                                p["length"],
                                p["pct_cys"],
                                p["pct_pro"],
                                p["pct_st"],
                                p["sequons_per_100aa"],
                                p["gpi"],
                                ";".join(sorted(pf)),
                            ]
                        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
