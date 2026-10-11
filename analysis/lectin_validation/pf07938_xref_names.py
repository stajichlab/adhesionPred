#!/usr/bin/env python3.12
"""Count fungal UniProt proteins with a PF07938 cross-reference by name class.

Input: the TSV of `data/raw/uniprot_fungi_pf07938_xref.tsv.gz` (UniProt stream, taxonomy 4751,
xref:pfam-PF07938; columns Entry, Reviewed, Protein names, Organism, Length, Pfam, Protein existence).
Name rule (case-insensitive): a name is "lectin or uncharacterized" when it matches
  lectin | fucose | FleA | uncharacterized | genomic scaffold
Every other name is a "non-lectin name". No paper is used. Names are UniProt automatic annotation.
Usage: pf07938_xref_names.py IN_TSV_GZ OUT_TSV
"""

import csv
import gzip
import re
import sys
from collections import Counter

LECTIN_OR_UNCHAR = re.compile(r"lectin|fucose|FleA|uncharacterized|genomic scaffold", re.I)
rows = list(csv.DictReader(gzip.open(sys.argv[1], "rt"), delimiter="\t"))
other = [r for r in rows if not LECTIN_OR_UNCHAR.search(r["Protein names"])]
multi = [r for r in rows if len([p for p in r["Pfam"].split(";") if p]) > 1]
reviewed = [r for r in rows if r["Reviewed"] == "reviewed"]
pe = Counter(r["Protein existence"] for r in other)
names = Counter(re.sub(r" \(.*", "", r["Protein names"]) for r in other)
with open(sys.argv[2], "w") as fh:
    fh.write("item\tvalue\n")
    fh.write(f"entries\t{len(rows)}\nreviewed\t{len(reviewed)}\n")
    fh.write(
        f"non_lectin_names\t{len(other)}\nentries_with_more_than_one_pfam_xref\t{len(multi)}\n"
    )
    for k, v in sorted(pe.items()):
        fh.write(f"non_lectin_names_protein_existence:{k}\t{v}\n")
    for k, v in names.most_common():
        fh.write(f"name:{k}\t{v}\n")
print(len(rows), len(other), dict(pe), len(multi), len(reviewed))
