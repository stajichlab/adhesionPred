#!/usr/bin/env python3.12
"""Run hmmsearch --cut_ga of the 20 models (15 M1/M2 + 5 M3) on a FASTA file.

Usage: scan_models.py FASTA OUT_PREFIX [--pfam-hmm PATH] [--hmmsearch PATH] [--cpu N]
Writes OUT_PREFIX.models.hmm (the 20 models), OUT_PREFIX.tbl (hmmsearch --tblout, GA applied)
and OUT_PREFIX.hits.tsv (accession, model, model_name, sequence score, E-value, best-domain score).
A protein with no hit has no row. The tblout sequence score is the score used for score floors.
"""

import argparse
import subprocess
import sys

DEFAULT_PFAM = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"
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
M3 = ["PF03388", "PF00262", "PF07915", "PF02157", "PF11721"]


def extract(pfam, accs, out):
    want = set(accs)
    got = {}
    buf = []
    keep = None
    with open(pfam) as fh:
        for line in fh:
            buf.append(line)
            if line.startswith("ACC"):
                a = line.split()[1].split(".")[0]
                keep = a if a in want else None
                ver = line.split()[1]
            if line.startswith("//"):
                if keep:
                    got[keep] = (ver, "".join(buf))
                buf, keep = [], None
    missing = want - set(got)
    if missing:
        raise SystemExit("models not found: " + ",".join(sorted(missing)))
    with open(out, "w") as fh:
        for a in accs:
            fh.write(got[a][1])
    return {a: got[a][0] for a in accs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("prefix")
    ap.add_argument("--pfam-hmm", default=DEFAULT_PFAM)
    ap.add_argument("--hmmsearch", default="hmmsearch")
    ap.add_argument("--cpu", type=int, default=4)
    a = ap.parse_args()
    vers = extract(a.pfam_hmm, M12 + M3, a.prefix + ".models.hmm")
    with open(a.prefix + ".model_versions.tsv", "w") as fh:
        for k, v in vers.items():
            fh.write(f"{k}\t{v}\n")
    subprocess.run(
        [
            a.hmmsearch,
            "--cut_ga",
            "--noali",
            "--cpu",
            str(a.cpu),
            "--tblout",
            a.prefix + ".tbl",
            "-o",
            "/dev/null",
            a.prefix + ".models.hmm",
            a.fasta,
        ],
        check=True,
    )
    names = {}
    with open(a.prefix + ".models.hmm") as fh:
        n = None
        for line in fh:
            if line.startswith("NAME"):
                n = line.split()[1]
            if line.startswith("ACC"):
                names[line.split()[1].split(".")[0]] = n
    with open(a.prefix + ".hits.tsv", "w") as out:
        out.write("accession\tmodel\tmodel_name\tscore\tevalue\tbest_domain_score\n")
        for line in open(a.prefix + ".tbl"):
            if line.startswith("#"):
                continue
            f = line.split()
            prot = f[0].split("|")[1] if "|" in f[0] else f[0]
            acc = f[3].split(".")[0]
            out.write(f"{prot}\t{acc}\t{names[acc]}\t{f[5]}\t{f[4]}\t{f[8]}\n")
    print("done", file=sys.stderr)


main()
