# ruff: noqa
"""Write the 61 evidence-sheet proteins to a FASTA (header: proteome__id) and per-proteome FASTAs.
Usage: extract61.py <hydrophobin_truth dir> <outdir>"""

import csv
import sys
from pathlib import Path

d, out = Path(sys.argv[1]), Path(sys.argv[2])
runs = {r["name"]: r["fasta"] for r in csv.DictReader(open(d / "run_list.tsv"), delimiter="\t")}
rows = list(csv.DictReader(open(d / "evidence_sheets.tsv"), delimiter="\t"))


def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


cache = {}
n = 0
with open(out / "q61.faa", "w") as fh:
    for r in rows:
        p = r["proteome"]
        if p not in cache:
            cache[p] = read_fasta(runs[p])
        s = cache[p][r["id"]].rstrip("*")
        assert len(s) == int(r["length"]), (r["id"], len(s), r["length"])
        fh.write(f">{p}__{r['id']}\n{s}\n")
        n += 1
print(n, "sequences")
