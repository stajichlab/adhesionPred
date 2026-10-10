#!/usr/bin/env python3.12
"""Join SignalP 6 prediction_results.txt with the UniLectin3D fungal table.

Usage: signalp6_join.py PREDICTION_RESULTS_TXT UNILECTIN_TABLE_TSV OUT_TSV
Output columns: uniprot, species, unilectin_class, pfam_ga_hits, signalp6_prediction, p_sp, cs_position.
Rows of PREDICTION_RESULTS_TXT without a table row (for example C1G842, paracoccin) are kept with empty table fields.
"""

import csv
import sys

pred_txt, table_tsv, out = sys.argv[1:4]
tab = {r["uniprot"]: r for r in csv.DictReader(open(table_tsv), delimiter="\t")}
rows = []
for line in open(pred_txt):
    if line.startswith("#"):
        continue
    p = line.rstrip("\n").split("\t")
    acc, pred, p_sp = p[0], p[1], p[3]
    cs = p[4] if len(p) > 4 else ""
    t = tab.get(acc, {})
    rows.append(
        [
            acc,
            t.get("species", ""),
            t.get("unilectin_class", ""),
            t.get("pfam_ga_hits", ""),
            pred,
            p_sp,
            cs,
        ]
    )
with open(out, "w") as fh:
    fh.write(
        "uniprot\tspecies\tunilectin_class\tpfam_ga_hits\tsignalp6_prediction\tp_sp\tcs_position\n"
    )
    for r in rows:
        fh.write("\t".join(r) + "\n")
n_sp = sum(1 for r in rows if r[4] != "OTHER")
print(f"{len(rows)} sequences, {n_sp} with a SignalP 6 signal peptide")
