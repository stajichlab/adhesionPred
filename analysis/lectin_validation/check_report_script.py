#!/usr/bin/env python3.12
"""Run analysis/lectin_class_report/lectin_report.py on the truth set (one FASTA per class) and compare
its per-protein hits with the hmmsearch of scan_models.py (same Pfam file, --cut_ga).

Usage (hmmsearch on PATH): check_report_script.py REPORT_SCRIPT TRUTH_TSV SEQUENCES_FAA SCAN_HITS_TSV WORK_DIR OUT_TSV
Compares, for each class (treated as a proteome) and each of the 15 M1/M2 models:
  proteins with a hit (script) vs (own scan), and the best score per protein.
Writes OUT_TSV with one row per class and model and a final line with the number of mismatches.
"""

import csv
import os
import subprocess
import sys
from collections import defaultdict

M12 = [
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
]
script, truth_path, seq_path, scan_path, work, out_path = sys.argv[1:7]
os.makedirs(work, exist_ok=True)
truth = list(csv.DictReader(open(truth_path), delimiter="\t"))
cls = {r["accession"]: r["label_class"] for r in truth}
seqs, cur = {}, None
for line in open(seq_path):
    if line.startswith(">"):
        cur = line[1:].split()[0]
        seqs[cur] = []
    elif cur:
        seqs[cur].append(line.strip())
files = {}
for c in sorted(set(cls.values())):
    p = os.path.join(work, f"class_{c}.faa")
    with open(p, "w") as fh:
        for a, k in cls.items():
            if k == c:
                fh.write(f">{a}\n{''.join(seqs[a])}\n")
    files[c] = p
subprocess.run(
    [
        sys.executable,
        script,
        "--out-dir",
        os.path.join(work, "report"),
        "--tmpdir",
        work,
        "--cpu",
        "4",
        *files.values(),
    ],
    check=True,
)
# script output: hits file
hits_path = os.path.join(work, "report", "lectin_report_hits.tsv")
rows = list(csv.DictReader(open(hits_path), delimiter="\t"))
cols = list(rows[0].keys())
# find columns
pcol, mcol, scol, ycol = "protein", "model_acc", "seq_score", "proteome"
for need in (pcol, mcol, scol, ycol):
    if need not in cols:
        raise SystemExit("column missing in report hits file: " + need)
script_hits = defaultdict(dict)
for r in rows:
    m = r[mcol].split(".")[0]
    script_hits[(r[ycol].removeprefix("class_"), m)][r[pcol]] = float(r[scol])
own = defaultdict(dict)
for r in csv.DictReader(open(scan_path), delimiter="\t"):
    if r["model"] in M12:
        own[(cls[r["accession"]], r["model"])][r["accession"]] = float(r["score"])
n_mis = 0
with open(out_path, "w") as out:
    out.write(
        "class\tmodel\tscript_proteins_with_hit\town_proteins_with_hit\tsame_set\tmax_abs_score_diff\n"
    )
    for c in sorted(files):
        for m in M12:
            s = script_hits.get((c, m), {})
            o = own.get((c, m), {})
            same = set(s) == set(o)
            diff = max([abs(s[k] - o[k]) for k in s if k in o] or [0.0])
            if not same or diff > 0.05:
                n_mis += 1
            out.write(f"{c}\t{m}\t{len(s)}\t{len(o)}\t{'yes' if same else 'NO'}\t{diff:.2f}\n")
    out.write(f"#mismatches\t{n_mis}\n")
print("mismatches", n_mis)
