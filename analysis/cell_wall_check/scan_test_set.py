#!/usr/bin/env python3.12
"""hmmsearch --cut_ga of the sorting-hat family-table models (plus context models) on test_set.faa.
Usage: scan_test_set.py FAMILY_TABLE FASTA OUT_PREFIX [--pfam-hmm P] [--hmmsearch P]
Measure only. Writes OUT_PREFIX.models.hmm, .domtbl, .tbl, .hits.tsv (one row per domain), .versions.tsv.
Context models are NOT in the family table; they are only there to name what else matches."""
import argparse, csv, subprocess
DEFAULT_PFAM = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"
CONTEXT = ["PF03198","PF00722","PF00332","PF03663","PF00150","PF00704","PF03856","PF12454","PF10342","PF07983","PF01822","PF00026"]
ap = argparse.ArgumentParser()
ap.add_argument("table"); ap.add_argument("fasta"); ap.add_argument("prefix")
ap.add_argument("--pfam-hmm", default=DEFAULT_PFAM); ap.add_argument("--hmmsearch", default="hmmsearch")
a = ap.parse_args()
fam = {r["pfam_acc"]: r for r in csv.DictReader(open(a.table), delimiter="\t")}
want = list(fam) + [c for c in CONTEXT if c not in fam]
ws = set(want); got = {}; buf = []; keep = None; ver = None
for line in open(a.pfam_hmm):
    buf.append(line)
    if line.startswith("ACC"):
        ver = line.split()[1]; k = ver.split(".")[0]; keep = k if k in ws else None
    if line.startswith("//"):
        if keep: got[keep] = (ver, "".join(buf))
        buf, keep = [], None
miss = ws - set(got)
if miss: raise SystemExit("models not found: " + ",".join(sorted(miss)))
names = {}
with open(a.prefix + ".models.hmm", "w") as fh, open(a.prefix + ".versions.tsv", "w") as vh:
    for k in want:
        fh.write(got[k][1])
        nm = [l.split()[1] for l in got[k][1].splitlines() if l.startswith("NAME")][0]
        names[k] = nm
        vh.write(f"{k}\t{got[k][0]}\t{nm}\t{'family_table' if k in fam else 'context'}\n")
subprocess.run([a.hmmsearch, "--cut_ga", "--noali", "--tblout", a.prefix + ".tbl", "--domtblout",
                a.prefix + ".domtbl", "-o", "/dev/null", a.prefix + ".models.hmm", a.fasta], check=True)
with open(a.prefix + ".hits.tsv", "w") as out:
    out.write("protein\tmodel\tmodel_name\tsource\tseq_score\tseq_evalue\tdom_i\tdom_score\tdom_ievalue\tenv_from\tenv_to\n")
    for l in open(a.prefix + ".domtbl"):
        if l.startswith("#"): continue
        f = l.split()
        k = f[4].split(".")[0]
        out.write("\t".join([f[0], k, f[3], "family_table" if k in fam else "context", f[7], f[6], f[9], f[13], f[12], f[19], f[20]]) + "\n")
