#!/usr/bin/python3.12
"""Evidence table for the lectin class spec (docs/superpowers/specs/2026-10-10-lectin-class-design.md, section 2).

Steps run on 2026-10-10 (login node, no SLURM):
1. Query the public UniLectin search API (POST https://unilectin.unige.ch/api/search-module, module unilectin3d,
   rows 100) for the phrase "Fungal and yeast lectins" and for every class, species and family name that came back,
   plus a few clade words. Keep documents with origin "Fungal and yeast lectins". The API returns at most 100
   documents per call and ignores `start`, so the union is a lower bound (197 structures found; the phrase query
   reports numFound 202). Saved as fungal_docs.json (--docs).
2. Fetch UniProt sequences for the distinct accessions (rest.uniprot.org, accessions endpoint). One accession is a
   UniParc ID (UPI000150AC21) and has no UniProt FASTA here.
3. `mmseqs easy-cluster seq.faa clu tmp --min-seq-id 0.3 -c 0.5` (MMseqs2 17-b804f) gives clu_cluster.tsv (--clusters).
4. `hmmscan --cut_ga` of the sequences against Pfam-A 38.2 gives scan.domtbl (--domtbl).
This script joins the three outputs into one TSV.

Usage: unilectin_fungal_table.py --docs fungal_docs.json --clusters clu_cluster.tsv --domtbl scan.domtbl --out OUT.tsv
"""

import argparse
import csv
import json
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    docs = json.load(open(a.docs))
    cl = {}
    for line in open(a.clusters):
        rep, mem = line.rstrip("\n").split("\t")
        cl[mem] = rep
    hits = defaultdict(set)
    for line in open(a.domtbl):
        if line.startswith("#"):
            continue
        f = line.split()
        acc = f[3].split("|")[1] if "|" in f[3] else f[3]
        hits[acc].add(f[1].split(".")[0] + ":" + f[0])
    by = {}
    n_struct = defaultdict(int)
    for x in docs:
        n_struct[x["uniprot"]] += 1
        by.setdefault(x["uniprot"], x)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "uniprot",
                "species",
                "unilectin_class",
                "unilectin_family",
                "protein_name",
                "n_pdb_structures",
                "cluster30",
                "pfam_ga_hits",
            ]
        )
        for acc, x in sorted(by.items(), key=lambda t: (t[1]["class"], t[1]["species"], t[0])):
            w.writerow(
                [
                    acc,
                    x["species"].strip(),
                    x["class"],
                    x["family"].strip(),
                    x["protein_name"],
                    n_struct[acc],
                    cl.get(acc, "NA"),
                    ",".join(sorted(hits.get(acc, []))) or "none",
                ]
            )


if __name__ == "__main__":
    main()
