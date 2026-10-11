#!/usr/bin/env python3.12
"""Copy the final truth files into controls/ (COMMON-RULES file names) and merge the rejected lists.

Usage: assemble_controls.py WORK_OUT_DIR CONTROLS_DIR REJECTED_MANUAL_TSV
Writes controls.tsv, sequences.faa, clusters.tsv, rejected.tsv in CONTROLS_DIR.
"""

import csv
import os
import shutil
import sys

work, dest, manual = sys.argv[1:4]
os.makedirs(dest, exist_ok=True)
shutil.copy(work + "/truth_final.tsv", dest + "/controls.tsv")
shutil.copy(work + "/sequences.faa", dest + "/sequences.faa")
shutil.copy(work + "/clusters.tsv", dest + "/clusters.tsv")
rows = list(csv.DictReader(open(work + "/rejected_build.tsv"), delimiter="\t"))
rows += list(csv.DictReader(open(manual), delimiter="\t"))
with open(dest + "/rejected.tsv", "w", newline="") as fh:
    w = csv.DictWriter(
        fh,
        fieldnames=["accession", "candidate_class", "reason"],
        delimiter="\t",
        lineterminator="\n",
    )
    w.writeheader()
    w.writerows(rows)
print(len(rows), "rejected lines")
