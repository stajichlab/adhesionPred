#!/usr/bin/env python3
"""Pick the proteomes searched for the PF28404 (ARB_05178) family.

Set: every Fungi_5k Eurotiomycetes proteome (Onygenales, Eurotiales, Chaetothyriales, ...), a
fixed-seed sample of other Pezizomycotina classes as an outgroup check, and the seven
Coccidioides reference proteomes of analysis/cys_candidates. Coccidioides pangenome strains are
not searched: strain copies add no new family members for this question.

Output TSV: label, path, class, order, source. Reads plain or gzip input.
"""

import argparse
import csv
import gzip
import random
import sys

OUTGROUP_CLASSES = ["Sordariomycetes", "Dothideomycetes", "Leotiomycetes", "Pezizomycetes",
                    "Lecanoromycetes"]  # fmt: skip


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--unit-genomes", required=True, help="analysis/cocci_repeats/unit_genomes.tsv")
    ap.add_argument(
        "--cocci-proteomes", required=True, help="cys_candidates proteomes.tsv (3 columns)"
    )
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-outgroup-class", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20261002)
    a = ap.parse_args()
    with opener(a.unit_genomes) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["source"] == "fungi5k"]
    chosen = [r for r in rows if r["class"] == "Eurotiomycetes"]
    rng = random.Random(a.seed)
    for cls in OUTGROUP_CLASSES:
        pool = sorted((r for r in rows if r["class"] == cls), key=lambda r: r["label"])
        chosen += rng.sample(pool, min(a.per_outgroup_class, len(pool)))
    out, seen = [], set()
    with opener(a.cocci_proteomes) as fh:
        for name, path, _sp in csv.reader(fh, delimiter="\t"):
            out.append((name, path, "Eurotiomycetes", "Onygenales", "cocci_reference"))
            seen.add(path)
    for r in chosen:
        if r["path"] not in seen:
            out.append((r["label"], r["path"], r["class"], r["order"], "fungi5k"))
            seen.add(r["path"])
    if len({o[0] for o in out}) != len(out):
        sys.exit("STOP: duplicate proteome labels")
    with open(a.out, "w") as fh:
        fh.write("label\tpath\tclass\torder\tsource\n")
        for o in out:
            fh.write("\t".join(o) + "\n")
    print(f"{len(out)} proteomes written to {a.out}")


if __name__ == "__main__":
    main()
