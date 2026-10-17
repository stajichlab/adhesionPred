#!/usr/bin/env python3.12
"""Combine hmmsearch --tblout files (one per proteome) into one TSV.

Usage: combine_tbl.py SCAN_DIR OUT_TSV
Columns: proteome, protein, model, model_acc, seq_evalue, seq_score.
"""

import glob
import os
import sys

scan_dir, out = sys.argv[1:3]
with open(out, "w") as fh:
    fh.write("proteome\tprotein\tmodel\tmodel_acc\tseq_evalue\tseq_score\n")
    for f in sorted(glob.glob(os.path.join(scan_dir, "*.tbl"))):
        b = os.path.basename(f)[:-4]
        for line in open(f):
            if line.startswith("#"):
                continue
            p = line.split()
            fh.write("\t".join([b, p[0], p[2], p[3], p[4], p[5]]) + "\n")
