#!/usr/bin/env python3.12
"""For positives missed by the 15 models: (a) hits to the 20 models below GA (E <= 10),
(b) all Pfam-A families that hit at the gathering cutoff (hmmscan --cut_ga).

Usage: scan_misses.py MISSES_TSV SEQUENCES_FAA MODELS_HMM OUT_DIR [--pfam-hmm PATH] [--tmpdir DIR]
Needs hmmsearch and hmmscan on PATH (module load hmmer/3.4).
"""

import argparse
import csv
import os
import subprocess

DEFAULT_PFAM = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"

ap = argparse.ArgumentParser()
ap.add_argument("misses")
ap.add_argument("seqs")
ap.add_argument("models")
ap.add_argument("out_dir")
ap.add_argument("--pfam-hmm", default=DEFAULT_PFAM)
ap.add_argument("--tmpdir", required=True)
a = ap.parse_args()
os.makedirs(a.tmpdir, exist_ok=True)
os.makedirs(a.out_dir, exist_ok=True)
want = [r["accession"] for r in csv.DictReader(open(a.misses), delimiter="\t")]
seqs, cur = {}, None
for line in open(a.seqs):
    if line.startswith(">"):
        cur = line[1:].split()[0]
        seqs[cur] = []
    elif cur:
        seqs[cur].append(line.strip())
fa = a.tmpdir + "/misses.faa"
with open(fa, "w") as fh:
    for k in want:
        fh.write(f">{k}\n{''.join(seqs[k])}\n")
# GA of each of the 20 models
ga, name = {}, {}
acc = None
for line in open(a.models):
    if line.startswith("ACC"):
        acc = line.split()[1].split(".")[0]
    if line.startswith("NAME"):
        nm = line.split()[1]
    if line.startswith("GA"):
        ga[acc] = float(line.split()[1])
        name[acc] = nm
subprocess.run(
    [
        "hmmsearch",
        "-E",
        "10",
        "--domE",
        "10",
        "--noali",
        "--tblout",
        a.tmpdir + "/sub.tbl",
        "-o",
        "/dev/null",
        a.models,
        fa,
    ],
    check=True,
)
with open(a.out_dir + "/miss_subga.tsv", "w") as out:
    out.write("accession\tmodel\tscore\tga\tevalue\n")
    for line in open(a.tmpdir + "/sub.tbl"):
        if line.startswith("#"):
            continue
        f = line.split()
        m = f[3].split(".")[0]
        out.write(f"{f[0]}\t{m}\t{f[5]}\t{ga[m]}\t{f[4]}\n")
subprocess.run(
    [
        "hmmscan",
        "--cut_ga",
        "--noali",
        "--tblout",
        a.tmpdir + "/pf.tbl",
        "-o",
        "/dev/null",
        a.pfam_hmm,
        fa,
    ],
    check=True,
)
with open(a.out_dir + "/miss_pfamscan.tsv", "w") as out:
    out.write("accession\tfamily\tscore\n")
    for line in open(a.tmpdir + "/pf.tbl"):
        if line.startswith("#"):
            continue
        f = line.split()
        out.write(f"{f[2]}\t{f[0]}\t{f[5]}\n")
