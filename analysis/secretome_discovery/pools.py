#!/usr/bin/python3.12
"""Pools of the secretome discovery pass (B3 of docs/superpowers/specs/2026-10-09-secretome-discovery-pass-design.md, section 5.1).

Pool S: signal_peptide_protein[R0] called, no mature TM helix, no C-terminal HDEL or KDEL.
Pool U: pool S and other_surface_no_mechanism[R0] called.
Per proteome, writes <out>/<name>.pools.tsv.gz: one row per pool S protein with the pool, length, composition, the GPI call
(Fungi_5k PredGPI GPI-anchor feature) and the Pfam-A families (Fungi_5k full Pfam-A, release 37.x; unversioned accessions).
Usage: pools.py --name NAME --fasta FASTA --workdir DIR --out OUTDIR
"""

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

RES = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/results/function")
SEQUON = re.compile(r"(?=N[^P][ST])")


def read_fasta(path):
    seqs, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def composition(seq):
    n = len(seq)
    pct = lambda letters: round(100 * sum(seq.count(x) for x in letters) / n, 2) if n else 0.0  # noqa: E731
    return {
        "length": n,
        "pct_cys": pct("C"),
        "pct_pro": pct("P"),
        "pct_st": pct("ST"),
        "sequons_per_100aa": round(100 * len(SEQUON.findall(seq)) / n, 2) if n else 0.0,
    }


def er_retained(seq):
    return seq.endswith("HDEL") or seq.endswith("KDEL")


def _lines(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        yield from fh


def read_pfam(path):
    out = {}
    for line in _lines(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        out.setdefault(f[0], set()).add(f[4].split(".")[0])
    return out


def read_gpi(path):
    out = set()
    for line in _lines(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2 and f[2] == "GPI-anchor":
            out.add(f[0].split()[0])
    return out


def read_calls(path):
    """(proteins with signal_peptide_protein[R0] called, proteins with other_surface_no_mechanism[R0] called)."""
    sp, nm = set(), set()
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["value"] != "called":
                continue
            if r["call"] == "signal_peptide_protein":
                sp.add(r["protein"])
            elif r["call"] == "other_surface_no_mechanism":
                nm.add(r["protein"])
    return sp, nm


def read_tm(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["state"] == "ok":
                out[r["id"]] = int(r["n_tm_mature"])
    return out


def pools(seqs, sp, no_mech, n_tm_mature):
    s = {p for p in sp if p in seqs and n_tm_mature.get(p, 0) == 0 and not er_retained(seqs[p])}
    return s, s & no_mech


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    wd = Path(a.workdir)
    seqs = read_fasta(a.fasta)
    sp, nm = read_calls(wd / "out" / "calls.long.tsv.gz")
    tm = read_tm(wd / "modules" / "tm.tsv.gz")
    s, u = pools(seqs, sp, nm, tm)
    pfam = read_pfam(RES / "pfam" / f"{a.name}.domtblout.gz")
    gpi = read_gpi(RES / "predgpi" / f"{a.name}.predgpi.gff3.gz")
    Path(a.out).mkdir(parents=True, exist_ok=True)
    cols = [
        "protein",
        "pool",
        "length",
        "pct_cys",
        "pct_pro",
        "pct_st",
        "sequons_per_100aa",
        "gpi",
        "pfam",
    ]
    with gzip.open(Path(a.out) / f"{a.name}.pools.tsv.gz", "wt") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for p in sorted(s):
            c = composition(seqs[p])
            w.writerow(
                [
                    p,
                    "U" if p in u else "S",
                    c["length"],
                    c["pct_cys"],
                    c["pct_pro"],
                    c["pct_st"],
                    c["sequons_per_100aa"],
                    int(p in gpi),
                    ";".join(sorted(pfam.get(p, ()))),
                ]
            )
    print(f"{a.name}: pool S {len(s)}, pool U {len(u)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
