#!/usr/bin/python3.12
"""Aspergillus supplement and second-proteome check of the secretome discovery pass (spec 5.3 and section 6.7).

clade: per family and per clade of the Aspergillus supplement species: species with the family and median copies per 1,000 pool proteins.
second: presence concordance between the first and the second proteome of a genus, and the family test re-run with the second proteome
substituted for the first in every target genus that has one.
Usage: extra_checks.py clade|second ...  (see main)
"""

import argparse
import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_discovery as rd  # noqa: E402
import stats  # noqa: E402

median = statistics.median


def clade_table(species, fams, pools, families):
    """species: {name: clade}; fams: {name: {family: hits}}; pools: {name: pool size}."""
    rows = []
    for f in families:
        for c in sorted(set(species.values())):
            names = [n for n, cl in species.items() if cl == c]
            hits = [fams.get(n, {}).get(f, 0) for n in names]
            copies = [1000 * h / pools[n] for h, n in zip(hits, names, strict=True)]
            rows.append(
                {
                    "family": f,
                    "clade": c,
                    "n_species": len(names),
                    "n_present": sum(1 for h in hits if h > 0),
                    "copies_per_1000_median": median(copies),
                }
            )
    return rows


def concordance(first, second, families):
    """first, second: {genus: set of families}; returns {family: (genera with the same presence call, genera compared)}."""
    genera = sorted(set(first) & set(second))
    return {
        f: (sum(1 for g in genera if (f in first[g]) == (f in second[g])), len(genera))
        for f in families
    }


def substitute(data, second):
    out = dict(data)
    for g, (size, fam) in second.items():
        if g in out:
            out[g] = (out[g][0], size, fam)
    return out


def _counts(pools_dir, name, pool):
    return rd.genus_counts(Path(pools_dir) / f"{name}.pools.tsv.gz", pool)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("clade")
    c.add_argument("--supplement", required=True)
    c.add_argument("--pools", required=True)
    c.add_argument("--families", nargs="+", required=True)
    c.add_argument("--pool", default="U")
    c.add_argument("--out", required=True)
    s = sub.add_parser("second")
    s.add_argument("--selection", required=True)
    s.add_argument("--second", required=True)
    s.add_argument("--pools", required=True)
    s.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "clade":
        sup = list(csv.DictReader(open(a.supplement), delimiter="\t"))
        species = {r["name"]: r["clade"] for r in sup}
        fams, pools = {}, {}
        for r in sup:
            p, f = _counts(a.pools, r["name"], a.pool)
            pools[r["name"]], fams[r["name"]] = p, f
        rows = clade_table(species, fams, pools, a.families)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(len(sup), "species;", len(set(species.values())), "clades;", len(rows), "rows")
        return 0
    sel = list(csv.DictReader(open(a.selection), delimiter="\t"))
    sec = list(csv.DictReader(open(a.second), delimiter="\t"))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for pool in ("S", "U"):
        data = {r["genus"]: (r["group"], *_counts(a.pools, r["name"], pool)) for r in sel}
        second = {r["genus"]: _counts(a.pools, r["name"], pool) for r in sec}
        first_sets = {g: set(data[g][2]) for g in second}
        second_sets = {g: set(second[g][1]) for g in second}
        fam_all = sorted({f for s_ in (first_sets | second_sets).values() for f in s_})
        conc = concordance(first_sets, second_sets, fam_all)
        agree = sum(v[0] for v in conc.values())
        total = sum(v[1] for v in conc.values())
        print(
            pool,
            "presence concordance over",
            len(second),
            "genera and",
            len(fam_all),
            "families:",
            f"{agree}/{total} = {agree / total:.3f}",
        )
        data2 = substitute(data, second)
        for target in rd.TARGETS:
            for label, d in (("first", data), ("mixed_second", data2)):
                sub_d = {
                    g: (("target" if grp == target else "background"), sz, h)
                    for g, (grp, sz, h) in d.items()
                    if grp in (target, "background")
                }
                res = stats.report(rd.family_rows(sub_d))
                names = sorted(x["family"] for x in res["enriched"] if x["fisher_q"] < 0.05)
                with open(out / f"second_{pool}_{target}_{label}.txt", "w") as fh:
                    fh.write("\n".join(names) + "\n")
                print(pool, target, label, "enriched q<0.05:", len(names))
    return 0


if __name__ == "__main__":
    sys.exit(main())
