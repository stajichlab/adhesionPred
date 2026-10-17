# ruff: noqa
"""Fetch UniProt FASTA for accessions (first column of a file). Usage: fetch_fasta.py LIST OUT.faa"""

import subprocess
import sys
import time
import urllib.parse

accs = [l.split("\t")[0].strip() for l in open(sys.argv[1]) if l.strip()]
seqs = {}
for i in range(0, len(accs), 40):
    q = " OR ".join(f"accession:{a}" for a in accs[i : i + 40])
    url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=" + urllib.parse.quote(q)
    text = subprocess.run(
        ["curl", "-s", "-m", "120", url], capture_output=True, text=True, check=True
    ).stdout
    cur = None
    for line in text.splitlines():
        if line.startswith(">"):
            cur = line.split("|")[1]
            seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
    time.sleep(0.3)
with open(sys.argv[2], "w") as f:
    for a in accs:
        if a in seqs:
            f.write(f">{a}\n{''.join(seqs[a])}\n")
missing = [a for a in accs if a not in seqs]
print(len(seqs), "fetched; missing:", " ".join(missing))
