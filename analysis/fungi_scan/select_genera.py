#!/usr/bin/python3.12
"""Choose one proteome per named genus for the secretome discovery pass (B1 of docs/superpowers/specs/2026-10-09-secretome-discovery-pass-design.md).

Groups: onygenales and aspergillaceae (targets), background (all other Pezizomycotina genera). Other subphyla are not used.
Rule per genus: drop fragmentary proteomes (under MIN_PROTEINS) and rows with an empty genus; drop proteomes whose taxon ID is not in the
taxonomy dump (a merged ID is mapped through merged.dmp first); drop known hybrids named in hybrids.tsv; prefer a RefSeq (GCF_) assembly;
take the proteome nearest the genus median protein count (ties by name). BUSCO completeness does not exist in Fungi_5k.
Usage: select_genera.py --out selection_discovery.tsv [--nodes ... --merged ... --hybrids hybrids.tsv]
"""

import argparse
import csv
import statistics
import sys
from pathlib import Path

ROOT = Path("/bigdata/stajichlab/shared/projects/Fungi_5k")
MIN_PROTEINS = 4000


def fasta_count(path):
    """Number of records of a FASTA file (protein_counts.tsv lists 1 for 16 files that hold thousands)."""
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def group_of(r):
    if r["SUBPHYLUM"] != "Pezizomycotina":
        return "out"
    if r["ORDER"] == "Onygenales":
        return "onygenales"
    if r["FAMILY"] == "Aspergillaceae":
        return "aspergillaceae"
    return "background"


def pick_per_genus(rows, known_taxa, merged=None, hybrids=frozenset(), min_proteins=MIN_PROTEINS):
    merged = merged or {}
    by = {}
    for r in rows:
        if not r["GENUS"] or r["_n"] < min_proteins or r["_name"] in hybrids:
            continue
        taxon = merged.get(r["NCBI_TAXONID"], r["NCBI_TAXONID"])
        if taxon not in known_taxa:
            continue
        r = dict(r, NCBI_TAXONID=taxon)
        by.setdefault(r["GENUS"], []).append(r)
    out = []
    for genus in sorted(by):
        cands = by[genus]
        refseq = [r for r in cands if r["ASMID"].startswith("GCF_")]
        pool = refseq or cands
        med = statistics.median(r["_n"] for r in cands)
        out.append(min(pool, key=lambda r: (abs(r["_n"] - med), r["_name"])))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nodes", default="_workdir/sorting_hat/taxdump/nodes.dmp")
    ap.add_argument("--merged", default="_workdir/sorting_hat/taxdump/merged.dmp")
    ap.add_argument("--hybrids", default=str(Path(__file__).with_name("hybrids.tsv")))
    a = ap.parse_args()
    rows = list(csv.DictReader(open(ROOT / "samples.csv")))
    prefix = {}
    for line in open(ROOT / "proteinfile2prefix.tsv"):
        name, tag = line.rstrip("\n").split("\t")[:2]
        prefix[tag] = name
    counts = {}
    for line in open(ROOT / "protein_counts.tsv"):
        f = line.rstrip("\n").split(":")
        counts[f[0].replace(".proteins.fa", "")] = int(f[-1])
    for r in rows:
        name = prefix.get(r["LOCUSTAG"])
        r["_name"] = name
        r["_n"] = counts.get(name, 0) if name else 0
        if name and r["_n"] < 100 and (ROOT / "input" / f"{name}.proteins.fa").exists():
            r["_n"] = fasta_count(ROOT / "input" / f"{name}.proteins.fa")
    rows = [
        r for r in rows if r["_name"] and (ROOT / "input" / f"{r['_name']}.proteins.fa").exists()
    ]
    known = {line.split("\t|\t")[0] for line in open(a.nodes)}
    merged = {}
    for line in open(a.merged):
        f = line.split("\t|\t")
        merged[f[0]] = f[1].split("\t")[0]
    hybrids = {r["name"] for r in csv.DictReader(open(a.hybrids), delimiter="\t") if r["name"]}
    for g in ("onygenales", "aspergillaceae", "background"):
        grp = [r for r in rows if group_of(r) == g]
        picked = pick_per_genus(grp, known, merged, hybrids)
        named = {r["GENUS"] for r in grp if r["GENUS"]}
        print(g, "genera with a passing proteome:", len(picked), "of", len(named), "named")
        print("  failing:", ", ".join(sorted(named - {r["GENUS"] for r in picked})) or "none")
        mode = "w" if g == "onygenales" else "a"
        with open(a.out, mode, newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            if mode == "w":
                w.writerow(
                    [
                        "group",
                        "name",
                        "fasta",
                        "taxon",
                        "order",
                        "family",
                        "genus",
                        "asmid",
                        "proteins",
                    ]
                )
            for r in picked:
                w.writerow(
                    [
                        g,
                        r["_name"],
                        str(ROOT / "input" / f"{r['_name']}.proteins.fa"),
                        r["NCBI_TAXONID"],
                        r["ORDER"],
                        r["FAMILY"],
                        r["GENUS"],
                        r["ASMID"],
                        r["_n"],
                    ]
                )
    return 0


if __name__ == "__main__":
    sys.exit(main())
