#!/usr/bin/python3.12
"""Descriptive summary of the sorting hat run on the scan proteomes (no truth: counts and candidate lists, not measurements).

Reads <workdir>/out/calls.long.tsv.gz and the Pfam module tables of each proteome in selection.tsv and run_list.tsv.
Writes proteome_summary.tsv, family_summary.tsv, attachment_basis.tsv and attachment_candidates.tsv in --out.
Usage: summarize_scan.py --dir analysis/fungi_scan --out docs/reports/data/sorting_hat/fungi_scan
"""

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

CALLS = [
    ("signal_peptide_protein", "R0"),
    ("tandem_repeat_protein", ""),
    ("wall_family_domain", ""),
    ("hydrophobin_domain", ""),
    ("hsba_domain", ""),
    ("cell_wall_adhesion_candidate", "R0"),
    ("surface_attachment_candidate", "R0"),
    ("iuis_allergen_similarity", ""),
    ("iuis_allergen_homolog", ""),
    ("other_not_surface", "R0"),
    ("other_surface_no_mechanism", "R0"),
]


def per_1000(n, size):
    return 0.0 if size == 0 else 1000.0 * n / size


def count_calls(long_rows):
    c = defaultdict(Counter)
    for r in long_rows:
        c[(r["call"], r["variant"])][r["value"]] += 1
    return c


def basis_counts(long_rows, call):
    return Counter(
        r["other_basis"]
        for r in long_rows
        if r["call"] == call and r["value"] == "called" and r["other_basis"]
    )


def family_counts(module_rows):
    c = Counter()
    for r in module_rows:
        if r.get("hit") == "1":
            for f in filter(None, r.get("families", "").split(",")):
                c[f] += 1
    return dict(c)


def read_tsv_gz(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d, out = Path(a.dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    sel = {r["name"]: r for r in csv.DictReader(open(d / "selection.tsv"), delimiter="\t")}
    runs = list(csv.DictReader(open(d / "run_list.tsv"), delimiter="\t"))
    prow, frows, brows, cand = [], [], [], []
    for r in runs:
        name, w = r["name"], Path(r["workdir"])
        longp = w / "out/calls.long.tsv.gz"
        if not longp.exists():
            print("missing", name)
            continue
        rows = read_tsv_gz(longp)
        c = count_calls(rows)
        n = len({x["protein"] for x in rows})
        row = {
            "proteome": name,
            "order": sel[name]["order"],
            "genus": sel[name]["genus"],
            "species": sel[name]["species"],
            "proteins": n,
        }
        for call, variant in CALLS:
            k = c[(call, variant)]
            row[call + (f"[{variant}]" if variant else "")] = k["called"]
            if k["not_assessable"]:
                row[call + "_n/a"] = k["not_assessable"]
        prow.append(row)
        for mod in ("pfam_adhesion", "pfam_hydrophobin", "pfam_hsba"):
            mp = w / f"modules/{mod}.tsv.gz"
            if mp.exists():
                for fam, k in family_counts(read_tsv_gz(mp)).items():
                    frows.append({"proteome": name, "module": mod, "family": fam, "proteins": k})
        for basis, k in basis_counts(rows, "surface_attachment_candidate").items():
            brows.append({"proteome": name, "basis": basis, "proteins": k})
        seqlen = {}
        for line in open(r["fasta"]):
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqlen[cur] = 0
            else:
                seqlen[cur] += len(line.strip())
        for x in rows:
            if x["call"] == "surface_attachment_candidate" and x["value"] == "called":
                cand.append(
                    {
                        "proteome": name,
                        "protein": x["protein"],
                        "length": seqlen.get(x["protein"], ""),
                        "basis": x["other_basis"],
                    }
                )
    for fname, rs in (
        ("proteome_summary.tsv", prow),
        ("family_summary.tsv", frows),
        ("attachment_basis.tsv", brows),
        ("attachment_candidates.tsv", cand),
    ):
        if not rs:
            continue
        keys = list(dict.fromkeys(k for r in rs for k in r))
        with open(out / fname, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
            w.writeheader()
            w.writerows(rs)
    print(len(prow), "proteomes summarised;", len(cand), "attachment candidates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
