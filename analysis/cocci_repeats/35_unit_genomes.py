#!/usr/bin/python3.12
"""Genome list for the SOWgp unit-distribution search (36-39).

Writes unit_genomes.tsv, one row per proteome FASTA:
  source    fungi5k | cocci_pangenome | cocci_longread
  path      absolute path of the FASTA
  label     file name without .proteins.fa
  locustag  protein-id prefix (Fungi_5k only)
  species, genus, family, order, class   taxonomy (Fungi_5k samples.csv; Cocci rows are
            set by hand: Onygenales / Coccidioides)

Fungi_5k file names are Genus_species_STRAIN and do not follow SPECIESIN, so each file is
matched to samples.csv by the locus-tag prefix of its first protein id. A file with no match
is reported and kept with empty taxonomy. Fungi_5k is ~5,800 files, so run on a compute node:
    srun -p short -c 2 --mem 4G -t 15 /usr/bin/python3.12 35_unit_genomes.py
"""

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "_common"))
import paths  # noqa: E402

sys.path.insert(0, str(HERE))
import importlib  # noqa: E402

fam = importlib.import_module("30_anchor_family_search")

FIELDS = ["source", "path", "label", "locustag", "species", "genus", "family", "order", "class"]


def first_prefix(path):
    with open(path) as fh:
        line = fh.readline()
    if not line.startswith(">"):
        return ""
    return line[1:].split()[0].split("_")[0]


def main():
    samples = {r["LOCUSTAG"]: r for r in csv.DictReader(open(paths.FUNGI5K_SAMPLES))}
    rows, unmatched = [], 0
    for f in sorted(Path(paths.FUNGI5K_INPUT).glob("*.proteins.fa")):
        tag = first_prefix(f)
        s = samples.get(tag)
        if s is None:
            unmatched += 1
            s = {}
        rows.append(
            {
                "source": "fungi5k",
                "path": str(f),
                "label": fam.label_of(f),
                "locustag": tag,
                "species": s.get("SPECIES", ""),
                "genus": s.get("GENUS", ""),
                "family": s.get("FAMILY", ""),
                "order": s.get("ORDER", ""),
                "class": s.get("CLASS", ""),
            }
        )
    print(f"Fungi_5k: {len(rows)} files, {unmatched} without a samples.csv match")
    for src, lst in (
        ("cocci_longread", fam.resolve_proteomes(["longread"])),
        ("cocci_pangenome", fam.resolve_proteomes(["pangenome"])),
    ):
        for f in lst:
            lab = fam.label_of(f)
            sp = (
                "Coccidioides immitis"
                if "immitis" in lab.lower()
                or lab.lower().startswith("cimm")
                or lab.startswith("Ci")
                else "Coccidioides posadasii"
                if "posadasii" in lab.lower() or lab.startswith("Cpos")
                else "Coccidioides sp."
            )
            rows.append(
                {
                    "source": src,
                    "path": str(f),
                    "label": lab,
                    "locustag": "",
                    "species": sp,
                    "genus": "Coccidioides",
                    "family": "Onygenaceae",
                    "order": "Onygenales",
                    "class": "Eurotiomycetes",
                }
            )
    with open(HERE / "unit_genomes.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    ony = [r for r in rows if r["source"] == "fungi5k" and r["order"] == "Onygenales"]
    print(f"total rows {len(rows)}; Fungi_5k Onygenales {len(ony)}")


if __name__ == "__main__":
    main()
