#!/usr/bin/env python3
"""Cluster all proteins of proteins.tsv locally with MMseqs2 (30 percent identity, 50 percent coverage; the setting
of analysis/calibration_truth/c1_truth_count.py) to get a leave-cluster-out grouping that spans every source table
(the source tables' own cluster ids use different namespaces). Local CPU, seconds. Writes clusters.tsv (ids only).
Needs `mmseqs` on PATH. The FASTA and MMseqs files stay in the cache outside the repo.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import CACHE, OUT, all_sequences, read_tsv, write_tsv  # noqa: E402


def main():
    seqs, _ = all_sequences()
    ids = [r["accession"] for r in read_tsv(OUT / "proteins.tsv")]
    work = CACHE / "mmseqs"
    work.mkdir(parents=True, exist_ok=True)
    fa = work / "all.faa"
    with open(fa, "w") as fh:
        for a in ids:
            fh.write(f">{a}\n{seqs[a]}\n")
    subprocess.run(
        ["mmseqs", "easy-cluster", str(fa), str(work / "clu"), str(work / "tmp"), "--min-seq-id", "0.3", "-c", "0.5", "--threads", "2", "-v", "1"],
        check=True,
    )  # fmt: skip
    rep = {}
    members = {}
    for ln in open(work / "clu_cluster.tsv"):
        r, m = ln.split()
        members.setdefault(r, []).append(m)
    # stable cluster names: C0001.. ordered by smallest member accession
    order = sorted(members, key=lambda r: min(members[r]))
    for i, r in enumerate(order, 1):
        for m in members[r]:
            rep[m] = (f"MC{i:04d}", len(members[r]))
    rows = [{"accession": a, "mmseqs_cluster": rep[a][0], "cluster_size": rep[a][1]} for a in ids]
    write_tsv(OUT / "clusters.tsv", rows, ["accession", "mmseqs_cluster", "cluster_size"])
    print(len(rows), "proteins in", len(members), "clusters")


if __name__ == "__main__":
    main()
