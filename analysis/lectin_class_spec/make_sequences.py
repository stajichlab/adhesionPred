#!/usr/bin/env python3.12
"""Make sequences.faa (header = UniProt accession) from a UniProt FASTA download.

Usage: make_sequences.py RAW_UNIPROT_FAA OUT_FAA
The raw file comes from:
  https://rest.uniprot.org/uniprotkb/accessions?accessions=<comma list of the accessions in
  unilectin_fungal_sequences.tsv, without the UniParc row>&format=fasta
UniProt data is licensed CC BY 4.0. Fetched 2026-10-10 (revision 1) and re-fetched for revision 2.
"""

import sys

raw, out = sys.argv[1:3]
n = 0
with open(raw) as fin, open(out, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            acc = line[1:].split("|")[1]
            fout.write(f">{acc}\n")
            n += 1
        else:
            fout.write(line)
print(f"{n} sequences")
