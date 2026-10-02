#!/usr/bin/env python3
"""Collect PF28404 hmmsearch hits (domtbl per proteome) into hits.tsv and hits.faa.

hits.tsv: one row per protein (best domain), with proteome class and order. hits.faa: the
protein sequences, header `<label>|<protein_id>`. Reads plain or gzip FASTA and domtbl files.
A hit is any domtbl line (the search used --cut_ga, so every line passed the model's gathering
threshold). STOP (exit 2) when a domtbl is missing for a proteome in proteomes.tsv.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with opener(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].split()[0], []
            else:
                buf.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True)
    a = ap.parse_args()
    work = Path(a.work)
    with open(work / "proteomes.tsv") as fh:
        prot = list(csv.DictReader(fh, delimiter="\t"))
    rows, fasta = [], []
    for p in prot:
        dom = work / "domtbl" / f"{p['label']}.domtbl"
        if not dom.exists():
            sys.exit(f"STOP: missing domtbl for {p['label']}")
        best = {}
        with open(dom) as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                f = ln.split()
                tid, tlen = f[0], int(f[2])
                ie, sc = float(f[12]), float(f[13])
                hf, ht = int(f[15]), int(f[16])
                if tid not in best or ie < best[tid][0]:
                    best[tid] = (ie, sc, hf, ht, tlen)
        if not best:
            continue
        seqs = read_fasta(p["path"])
        for tid, (ie, sc, hf, ht, tlen) in sorted(best.items()):
            s = seqs.get(tid)
            if s is None:
                sys.exit(f"STOP: {tid} not in {p['path']}")
            rows.append(
                [p["label"], tid, p["class"], p["order"], tlen, ie, sc, hf, ht, p["source"]]
            )
            fasta.append((f"{p['label']}|{tid}", s.rstrip("*")))
    with open(work / "hits.tsv", "w") as fh:
        fh.write(
            "label\tprotein_id\tclass\torder\tlength\ti_evalue\tscore\thmm_from\thmm_to\tsource\n"
        )
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(work / "hits.faa", "w") as fh:
        for n, s in fasta:
            fh.write(f">{n}\n{s}\n")
    print(f"{len(rows)} proteins in {len({r[0] for r in rows})} of {len(prot)} proteomes")


if __name__ == "__main__":
    main()
