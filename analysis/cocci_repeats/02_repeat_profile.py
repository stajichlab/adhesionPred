#!/usr/bin/env python3
"""Profile tandem-repeat structure, composition and size for every protein in a proteome.

Built to curate more class-2a examples (docs/TOOL-ARCHITECTURE.md): repeat/avidity surface
proteins. The transfer test in analysis/model_review/repeat_structure_transfer.py showed that
the current 2a model misses SOWgp because the training positives are long and Ser/Thr-rich
while SOWgp is short and Pro/Cys-rich. Finding more short, Pro/Cys-rich repeat proteins in
Coccidioides is the curation step that should close that gap.

LONG-READ ASSEMBLIES MATTER. Tandem arrays collapse or fragment in short-read assemblies, so
repeat copy number measured on those is unreliable. This is run on the UArizona long-read
strains alongside the pangenome references, and the two are reported separately.

Repeat detection is periodicity-based and composition-agnostic: for each candidate period p,
score the fraction of positions where s[i] == s[i+p], find the best p, then delimit the
maximal region sustaining that periodicity. Reports period, copy number, region extent and
the repeat unit, so candidates can be inspected by eye rather than trusted blindly.
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
MIN_PERIOD, MAX_PERIOD = 4, 80
WINDOW = 3  # smoothing half-window for the periodicity track


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def periodicity(s, p):
    """Fraction of positions identical to the residue p ahead."""
    n = len(s) - p
    if n <= 0:
        return 0.0
    return sum(1 for i in range(n) if s[i] == s[i + p]) / n


def best_repeat(s):
    """Best tandem period and the extent of the region sustaining it."""
    if len(s) < MIN_PERIOD * 3:
        return {
            "period": 0,
            "score": 0.0,
            "start": 0,
            "end": 0,
            "n_copies": 0,
            "coverage": 0.0,
            "unit": "",
        }
    best_p, best_s = 0, 0.0
    for p in range(MIN_PERIOD, min(MAX_PERIOD, len(s) // 3) + 1):
        sc = periodicity(s, p)
        if sc > best_s:
            best_p, best_s = p, sc
    if best_p == 0 or best_s < 0.3:
        return {
            "period": 0,
            "score": round(best_s, 3),
            "start": 0,
            "end": 0,
            "n_copies": 0,
            "coverage": 0.0,
            "unit": "",
        }
    # delimit the region: smoothed per-position match track at the best period
    p = best_p
    match = [1 if s[i] == s[i + p] else 0 for i in range(len(s) - p)]
    sm = [
        sum(match[max(0, i - WINDOW) : i + WINDOW + 1])
        / len(match[max(0, i - WINDOW) : i + WINDOW + 1])
        for i in range(len(match))
    ]
    inside, start, best_span = False, 0, (0, 0)
    for i, v in enumerate(sm + [0.0]):
        if v >= 0.5 and not inside:
            inside, start = True, i
        elif v < 0.5 and inside:
            inside = False
            if i - start > best_span[1] - best_span[0]:
                best_span = (start, i)
    a, b = best_span
    ext = (b - a) + p
    return {
        "period": p,
        "score": round(best_s, 3),
        "start": a,
        "end": a + ext,
        "n_copies": round(ext / p, 1) if p else 0,
        "coverage": round(ext / len(s), 3),
        "unit": s[a : a + p],
    }


def composition(s):
    n = max(len(s), 1)
    c = Counter(s)
    return {
        "pct_ser_thr": round(100 * (c["S"] + c["T"]) / n, 1),
        "pct_pro": round(100 * c["P"] / n, 1),
        "pct_cys": round(100 * c["C"] / n, 1),
        "pct_gly_ala": round(100 * (c["G"] + c["A"]) / n, 1),
        "pct_charged": round(100 * sum(c[x] for x in "DEKR") / n, 1),
        "top3_frac": round(sum(sorted(c.values(), reverse=True)[:3]) / n, 3),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("fasta", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=80)
    args = ap.parse_args()

    rows = []
    for fa in args.fasta:
        strain = Path(fa).name
        for suf in (".proteins.fa", ".fasta", ".fa"):
            strain = strain.replace(suf, "")
        n = 0
        for pid, seq in read_fasta(fa):
            seq = seq.rstrip("*")
            if len(seq) < args.min_len:
                continue
            r = best_repeat(seq)
            rows.append(
                {
                    "strain": strain,
                    "protein": pid,
                    "length": len(seq),
                    **{f"rep_{k}": v for k, v in r.items()},
                    **composition(seq),
                }
            )
            n += 1
        print(f"  {strain}: {n} proteins profiled", file=sys.stderr)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    rep = [r for r in rows if r["rep_coverage"] >= 0.3 and r["rep_n_copies"] >= 3]
    print(f"\nwrote {len(rows)} rows to {args.out}", file=sys.stderr)
    print(
        f"{len(rep)} proteins with a substantial tandem repeat (coverage >=0.3, >=3 copies)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
