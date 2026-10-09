#!/usr/bin/python3.12
"""Task L1: cluster the Jensen 2010 proteins (LP) with the T2 proteins, choose the v2 reserve, and report the
cysteine gaps of the 8 proteins whose stated pattern was not found in the resolved sequence.

The joint clustering here is used ONLY for the v2 reserve. Folds use clusters_positives.tsv.
Usage: lp_clusters.py --dir analysis/hydrophobin_truth   (needs mmseqs on PATH)
"""

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

UNVERIFIED = [
    "An07g03340",
    "An08g09880",
    "AN6401.2",
    "ACLA_001890",
    "ACLA_048810",
    "ACLA_072820",
    "ACLA_018290",
    "ACLA_007980",
]


def reserve(clusters, lp, t12, unverified):
    """LP clusters with no T1/T2 member and no unverified protein."""
    by = {}
    for p, c in clusters.items():
        by.setdefault(c, set()).add(p)
    keep = set()
    for c, members in by.items():
        if not members & lp:
            continue
        if members & t12 or members & unverified:
            continue
        keep.add(c)
    return keep


def actual_gaps(seq):
    pos = [i for i, a in enumerate(seq.upper()) if a == "C"]
    return [b - a - 1 for a, b in zip(pos, pos[1:], strict=False)]


def stated_gaps(pattern):
    m = re.fullmatch(r"CN\{(\d+)\}CCN\{(\d+)\}CN\{(\d+)\}CN\{(\d+)\}CCN\{(\d+)\}C", pattern or "")
    return [int(x) for x in m.groups()] if m else None


def read_fasta(path):
    out, k, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if k:
                out[k] = "".join(buf)
            k, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if k:
        out[k] = "".join(buf)
    return out


def cluster(seqs, tmp):
    names = list(seqs)
    faa = Path(tmp) / "in.faa"
    with open(faa, "w") as f:
        for i, k in enumerate(names):
            f.write(f">s{i}\n{seqs[k]}\n")
    pre = Path(tmp) / "clu"
    subprocess.run(
        [
            "mmseqs",
            "easy-cluster",
            str(faa),
            str(pre),
            str(Path(tmp) / "w"),
            "--min-seq-id",
            "0.3",
            "-c",
            "0.5",
        ],
        check=True,
        capture_output=True,
    )
    cl = {}
    for line in open(f"{pre}_cluster.tsv"):
        rep, mem = line.split()
        cl[names[int(mem[1:])]] = names[int(rep[1:])]
    return cl


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    lit = d / "literature"
    lp = read_fasta(lit / "jensen_sequences.faa")
    lp = {k.split("|")[1]: v for k, v in lp.items()}
    truth = [
        r for r in csv.DictReader(open(d / "truth_all.tsv"), delimiter="\t") if r["tier"] == "T2"
    ]
    t2 = {"T2|" + r["accession"]: r["sequence"] for r in truth}
    seqs = {**{"LP|" + k: v for k, v in lp.items()}, **t2}
    with tempfile.TemporaryDirectory() as tmp:
        cl = cluster(seqs, tmp)
    lp_ids = {"LP|" + k for k in lp}
    t_ids = set(t2)
    unv = {"LP|" + u for u in UNVERIFIED}
    keep = reserve(cl, lp_ids, t_ids, unv)
    rows = {
        r["protein_id"]: r
        for r in csv.DictReader(open(lit / "protein_lists.tsv"), delimiter="\t")
        if r["pmid"] == "21182770"
    }
    by = {}
    for p, c in cl.items():
        by.setdefault(c, set()).add(p)
    with open(lit / "lp_table.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_id", "species", "cluster", "cluster_has_t2", "unverified", "reserved"])
        for k in lp:
            c = cl["LP|" + k]
            w.writerow(
                [
                    k,
                    rows[k]["species"],
                    c,
                    int(bool(by[c] & t_ids)),
                    int("LP|" + k in unv),
                    int(c in keep),
                ]
            )
    n_lp_clusters = len({cl[p] for p in lp_ids})
    n_with_t2 = len({cl[p] for p in lp_ids if by[cl[p]] & t_ids})
    reserved = sorted(p for p in lp_ids if cl[p] in keep)
    with open(lit / "v2_reserve.tsv", "w") as fh:
        fh.write("protein_id\tcluster\n")
        for p in reserved:
            fh.write(f"{p[3:]}\t{cl[p]}\n")
    print(
        f"LP proteins {len(lp)}; LP clusters {n_lp_clusters}; with a T2 member {n_with_t2}; reserved clusters {len(keep)} ({len(reserved)} proteins)"
    )
    print(
        "gap check of the unverified proteins (stated gaps C1-C2, C3-C4, C4-C5, C5-C6, C7-C8 vs all consecutive gaps in the resolved sequence)"
    )
    with open(lit / "unverified_gap_report.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "protein_id",
                "stated_pattern",
                "stated_gaps",
                "n_cys",
                "length",
                "actual_consecutive_gaps",
            ]
        )
        for u in UNVERIFIED:
            s = lp[u]
            w.writerow(
                [
                    u,
                    rows[u]["cys_pattern_stated"],
                    stated_gaps(rows[u]["cys_pattern_stated"]),
                    s.upper().count("C"),
                    len(s),
                    actual_gaps(s),
                ]
            )
            print(
                u,
                rows[u]["cys_pattern_stated"],
                "| Cys",
                s.upper().count("C"),
                "len",
                len(s),
                "| actual",
                actual_gaps(s),
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
