#!/usr/bin/env python3.12
"""Write the N1 rows of the shared hard-negative set as a FASTA.

Usage: n1_fasta.py CONTROLS_TSV SEQUENCES_FAA OUT_FAA
N1 = evidence_level N1 (owner decision 2026-10-09, Q5: only N1 counts toward the 20-cluster floor).
"""

import csv
import sys

controls, faa, out = sys.argv[1:4]
n1 = {
    r["accession"]
    for r in csv.DictReader(open(controls), delimiter="\t")
    if r["evidence_level"] == "N1"
}
keep = False
written = 0
with open(faa) as fin, open(out, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            keep = line[1:].split()[0] in n1
            written += keep
        if keep:
            fout.write(line)
print(f"N1 accessions {len(n1)}, sequences written {written}")
