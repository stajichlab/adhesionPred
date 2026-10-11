#!/usr/bin/env python3.12
"""Write the list of accessions to fetch from UniProt (candidates for the truth set).

Union of: UniLectin3D fungal accessions (56), Swiss-Prot fungal entries with the keyword
Lectin or a UniProt Pfam cross-reference to a lectin-like (M1/M2) or house-keeping (M3) model,
paracoccin C1G842, and the hard-negative entries named in the spec (4.3).
Usage: make_acc_list.py SWISSPROT_TSV_GZ UNILECTIN_TSV OUT_TXT
"""

import csv
import gzip
import sys

PF = {
    "PF07938",
    "PF01476",
    "PF08881",
    "PF07367",
    "PF21691",
    "PF00337",
    "PF21595",
    "PF18647",
    "PF22803",
    "PF01419",
    "PF00652",
    "PF14200",
    "PF13385",
    "PF01607",
    "PF00187",
    "PF03388",
    "PF00262",
    "PF07915",
    "PF02157",
    "PF11721",
}
EXTRA = ["C1G842", "Q4WVZ3", "Q4WER9", "Q4WRR0", "J9VLH0", "G4N906", "Q00363"]

accs = []
for r in csv.DictReader(gzip.open(sys.argv[1], "rt"), delimiter="\t"):
    pf = {x for x in r["Pfam"].split(";") if x}
    if pf & PF or "Lectin" in r["Keywords"]:
        accs.append(r["Entry"])
for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
    if not r["uniprot"].startswith("UPI"):
        accs.append(r["uniprot"])
accs += EXTRA
seen = []
for a in accs:
    if a not in seen:
        seen.append(a)
open(sys.argv[3], "w").write("\n".join(seen) + "\n")
print(len(seen))
