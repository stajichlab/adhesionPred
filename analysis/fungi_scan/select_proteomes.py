#!/usr/bin/python3.12
"""Choose proteomes from Fungi_5k for the descriptive scan of cell wall protein calls (2026-10-09).

Source: /bigdata/stajichlab/shared/projects/Fungi_5k (samples.csv, protein_counts.tsv, proteinfile2prefix.tsv, input/*.proteins.fa).
One proteome per species within a genus, skipping fragmentary proteomes (fewer than MIN_PROTEINS), deterministic by (species, strain).
The taxon ID is NCBI_TAXONID from samples.csv and must exist in the taxonomy dump used by the tool.
Usage: select_proteomes.py --out selection.tsv [--nodes _workdir/sorting_hat/taxdump/nodes.dmp]
"""

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path("/bigdata/stajichlab/shared/projects/Fungi_5k")
MIN_PROTEINS = 4000
# (genus, species list or None, number of species)
TARGETS = [
    # Onygenales: the group that holds Coccidioides, Histoplasma, Blastomyces and the dermatophytes
    ("Coccidioides", None, 2),
    ("Histoplasma", None, 2),
    ("Blastomyces", None, 2),
    ("Paracoccidioides", None, 2),
    ("Emergomyces", None, 2),
    ("Emmonsia", None, 1),
    ("Uncinocarpus", None, 1),
    ("Trichophyton", ["rubrum", "benhamiae", "tonsurans"], 3),
    ("Arthroderma", None, 2),
    ("Microsporum", None, 2),
    ("Ascosphaera", None, 1),
    ("Malbranchea", None, 1),
    ("Amauroascus", None, 1),
    ("Chrysosporium", None, 1),
    ("Nannizziopsis", None, 1),
    ("Ophidiomyces", None, 1),
    ("Onygena", None, 1),
    # Aspergillaceae
    (
        "Aspergillus",
        [
            "fumigatus",
            "nidulans",
            "niger",
            "flavus",
            "terreus",
            "oryzae",
            "clavatus",
            "lentulus",
            "udagawae",
            "fischeri",
            "versicolor",
            "sydowii",
        ],
        12,
    ),
    ("Penicillium", None, 2),
    ("Talaromyces", None, 2),
    # yeasts
    ("Candida", ["albicans", "tropicalis", "parapsilosis"], 3),
    ("Nakaseomyces", None, 1),
    ("Candidozyma", None, 1),
    ("Saccharomyces", ["cerevisiae"], 1),
    ("Kluyveromyces", None, 1),
    # Basidiomycota
    ("Cryptococcus", None, 2),
    ("Ustilago", None, 1),
    ("Malassezia", None, 1),
    ("Pleurotus", None, 1),
    ("Schizophyllum", None, 1),
    # other Ascomycota and Mucoromycota
    ("Neurospora", ["crassa"], 1),
    ("Fusarium", ["graminearum"], 1),
    ("Pyricularia", None, 1),
    ("Rhizopus", None, 1),
]


def choose(rows, targets, known_taxa=None, min_proteins=MIN_PROTEINS):
    out = []
    for genus, species, n in targets:
        cands = [r for r in rows if r["GENUS"] == genus and r.get("_n", 10**9) >= min_proteins]
        if known_taxa is not None:
            cands = [r for r in cands if r["NCBI_TAXONID"] in known_taxa]
        if species:
            cands = [r for r in cands if r["SPECIES"] in {f"{genus} {x}" for x in species}]
        cands.sort(key=lambda r: (r["SPECIES"], r["STRAIN"], r["LOCUSTAG"]))
        seen, picked = set(), []
        for r in cands:
            if r["SPECIES"] in seen:
                continue
            seen.add(r["SPECIES"])
            picked.append(r)
            if len(picked) == n:
                break
        out += picked
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nodes", default="_workdir/sorting_hat/taxdump/nodes.dmp")
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
    rows = [
        r for r in rows if r["_name"] and (ROOT / "input" / f"{r['_name']}.proteins.fa").exists()
    ]
    known = {line.split("\t|\t")[0] for line in open(a.nodes)}
    picked = choose(rows, TARGETS, known_taxa=known)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["name", "fasta", "taxon", "order", "genus", "species", "strain", "proteins"])
        for r in picked:
            w.writerow(
                [
                    r["_name"],
                    str(ROOT / "input" / f"{r['_name']}.proteins.fa"),
                    r["NCBI_TAXONID"],
                    r["ORDER"],
                    r["GENUS"],
                    r["SPECIES"],
                    r["STRAIN"],
                    r["_n"],
                ]
            )
    print(
        len(picked),
        "proteomes;",
        sum(1 for r in picked if r["ORDER"] == "Onygenales"),
        "Onygenales",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
