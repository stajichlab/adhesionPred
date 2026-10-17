#!/usr/bin/env python3.12
"""Fetch the Pfam seed alignments of the models from the InterPro API and write one FASTA.

Usage: fetch_seeds.py OUT_DIR
Writes OUT_DIR/seeds.faa.gz (header: ACC|seqname) and OUT_DIR/seeds_summary.tsv (accession,
version from the file, number of seed sequences). The seed set is the leakage reference:
a truth sequence with 30% identity or more to a seed sequence of a model is `in_reference`
for that model. InterPro serves the current Pfam release (compare the version column with
the Pfam 38.2 file used for the scan).
"""

import gzip
import io
import sys
import time
import urllib.request

MODELS = [
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
]
out = sys.argv[1]
fa = gzip.open(out + "/seeds.faa.gz", "wt")
summ = open(out + "/seeds_summary.tsv", "w")
summ.write("accession\tversion\tname\tn_seed_sequences\n")
for acc in MODELS:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/?annotation=alignment:seed&download"
    for _attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                txt = gzip.open(io.BytesIO(r.read()), "rt").read()
            break
        except Exception as exc:
            print("retry", acc, exc, file=sys.stderr)
            time.sleep(3)
    else:
        raise SystemExit("failed " + acc)
    ver = name = ""
    seqs = {}
    for line in txt.splitlines():
        if line.startswith("#=GF AC"):
            ver = line.split()[2]
        elif line.startswith("#=GF ID"):
            name = line.split()[2]
        elif line and not line.startswith("#") and not line.startswith("//"):
            f = line.split()
            if len(f) == 2:
                seqs[f[0]] = seqs.get(f[0], "") + f[1]
    for k, s in seqs.items():
        s = "".join(c for c in s if c.isalpha()).upper()
        fa.write(f">{acc}|{k}\n{s}\n")
    summ.write(f"{acc}\t{ver}\t{name}\t{len(seqs)}\n")
    print(acc, ver, name, len(seqs))
fa.close()
summ.close()
