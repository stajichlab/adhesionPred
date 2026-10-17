#!/usr/bin/python3.12
"""Stability of the family results under 5 random subsamples of 150 background genera (spec section 3 and 7, B4).

Subsamples are stratified by CLASS (samples.csv) with allocation proportional to class size, at least one genus per class;
genera with an empty CLASS form their own stratum. Writes <out>/stability_<pool>_<target>.tsv: families with q < 0.05 in all
5 subsamples and the full-background q. Usage: stability.py --selection S --pools DIR --samples samples.csv --out DIR
"""

import argparse
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_discovery as rd  # noqa: E402
import stats  # noqa: E402


def subsample(by_class, n, rng):
    total = sum(len(v) for v in by_class.values())
    out = []
    for _, genera in sorted(by_class.items()):
        k = max(1, round(n * len(genera) / total))
        out += rng.sample(sorted(genera), min(k, len(genera)))
    return set(out)


def run_target(data, target, genera):
    sub = {
        g: ("target" if grp == target else "background", s, h)
        for g, (grp, s, h) in data.items()
        if grp == target or (grp == "background" and g in genera)
    }
    res = stats.report(rd.family_rows(sub))
    return {s["family"]: s["fisher_q"] for s in res["enriched"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--pools", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--seeds", type=int, default=5)
    a = ap.parse_args()
    sel = list(csv.DictReader(open(a.selection), delimiter="\t"))
    cls = {}
    for r in csv.DictReader(open(a.samples)):
        cls.setdefault(r["GENUS"], r["CLASS"])
    by_class = defaultdict(list)
    for r in sel:
        if r["group"] == "background":
            by_class[cls.get(r["genus"], "")].append(r["genus"])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for pool in ("S", "U"):
        data = {}
        for r in sel:
            f = Path(a.pools) / f"{r['name']}.pools.tsv.gz"
            size, fam = rd.genus_counts(f, pool)
            data[r["genus"]] = (r["group"], size, fam)
        for target in rd.TARGETS:
            full = run_target(
                data, target, {g for g, (grp, _, _) in data.items() if grp == "background"}
            )
            runs = [
                run_target(data, target, subsample(by_class, a.n, random.Random(seed)))
                for seed in range(a.seeds)
            ]
            fams = sorted(
                {f for q in runs for f in q if q[f] < 0.05},
                key=lambda f: full.get(f, 1.0),
            )
            with open(out / f"stability_{pool}_{target}.tsv", "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(
                    ["family", "full_q"]
                    + [f"q_seed{i}" for i in range(a.seeds)]
                    + ["n_seeds_q_lt_0.05"]
                )
                for f in fams:
                    qs = [q.get(f, 1.0) for q in runs]
                    w.writerow([f, full.get(f, 1.0)] + qs + [sum(1 for x in qs if x < 0.05)])
            all5 = sum(1 for f in fams if all(q.get(f, 1.0) < 0.05 for q in runs))
            print(
                pool,
                target,
                "q<0.05 in >=1 seed:",
                len(fams),
                "in all",
                a.seeds,
                ":",
                all5,
                "full-background q<0.05:",
                sum(1 for v in full.values() if v < 0.05),
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
