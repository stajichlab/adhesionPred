#!/usr/bin/env python3
"""Compare sorting hat calls between two A. fumigatus proteomes by best reciprocal BLAST hit.

Plan 2, Task 14 step 2. Pairs proteins of run A (reference) and run B by best reciprocal hit
(highest bitscore in each direction, mutual). Pairs with byte-identical sequences are reported
separately: they test determinism of the software only. For the other pairs the script bins
identity and coverage, builds the 2x2 table and Cohen's kappa for signal_peptide_protein[R0],
and counts agreement of every other call. Named cases are listed with their calls.

Input files may be plain or .gz. Usage:
  rbh_stability.py --a-fasta A.faa --b-fasta B.faa --a-calls A/out/calls.long.tsv.gz \
    --b-calls B/out/calls.long.tsv.gz --b-vs-a B_vs_A.tsv --a-vs-b A_vs_B.tsv --out-prefix OUT
BLAST tables use -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send bitscore evalue".
"""

import argparse
import collections
import csv
import gzip
import hashlib
import sys

NAMED = {  # accession in the Af293 UniProt proteome -> label
    "Q4WXC4": "CspA",
    "P41746": "RodA",
    "Q4WXJ1": "CalA",
    "P79017": "Asp f 2",
    "O60024": "Asp f 4",
    "O42799": "Asp f 7",
    "O60022": "Asp f 15",
}


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with opener(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).rstrip("*")
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).rstrip("*")
    return seqs


def best_hits(path):
    best = {}
    with opener(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            q, s, pid, ln, qlen, slen = (
                row[0],
                row[1],
                float(row[2]),
                int(row[3]),
                int(row[4]),
                int(row[5]),
            )
            bits = float(row[10])
            if q not in best or bits > best[q][0]:
                best[q] = (bits, s, pid, ln, qlen, slen)
    return best


def read_calls(path):
    calls = collections.defaultdict(dict)
    with opener(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            key = r["call"] + ("[" + r["variant"] + "]" if r["variant"] else "")
            calls[r["protein"]][key] = r["value"]
    return calls


def kappa(table):
    n = sum(table.values())
    if n == 0:
        return float("nan")
    labels = sorted({k[0] for k in table} | {k[1] for k in table})
    po = sum(table.get((x, x), 0) for x in labels) / n
    pe = sum(
        (sum(table.get((x, y), 0) for y in labels) / n)
        * (sum(table.get((y, x), 0) for y in labels) / n)
        for x in labels
    )
    return float("nan") if pe == 1 else (po - pe) / (1 - pe)


def bin_of(pid, cov):
    if pid >= 99 and cov >= 0.95:
        return "ge99_cov95"
    if pid >= 95 and cov >= 0.9:
        return "95_99_cov90"
    if pid >= 90 and cov >= 0.8:
        return "90_95_cov80"
    return "lt90_or_lowcov"


def short(pid):
    return pid.split("|")[1] if "|" in pid else pid


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--a-fasta", required=True)
    ap.add_argument("--b-fasta", required=True)
    ap.add_argument("--a-calls", required=True)
    ap.add_argument("--b-calls", required=True)
    ap.add_argument("--b-vs-a", required=True, help="BLAST of B queries against A")
    ap.add_argument("--a-vs-b", required=True, help="BLAST of A queries against B")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    fa, fb = read_fasta(a.a_fasta), read_fasta(a.b_fasta)
    ca, cb = read_calls(a.a_calls), read_calls(a.b_calls)
    sha_a = collections.defaultdict(list)
    for k, s in fa.items():
        sha_a[hashlib.sha256(s.encode()).hexdigest()].append(k)
    ba, ab = best_hits(a.b_vs_a), best_hits(a.a_vs_b)

    pairs, used_a = [], set()
    for qb, (_, sa, pid, ln, qlen, slen) in ba.items():
        back = ab.get(sa)
        if back is None or back[1] != qb:
            continue
        cov = ln / max(qlen, slen)
        ident = (
            hashlib.sha256(fb[qb].encode()).hexdigest() in sha_a
            and sa in sha_a[hashlib.sha256(fb[qb].encode()).hexdigest()]
        )
        pairs.append((sa, qb, pid, cov, ident))
        used_a.add(sa)
    out = []
    out.append(
        f"A proteins {len(fa)}, B proteins {len(fb)}, reciprocal best hit pairs {len(pairs)}"
    )
    out.append(
        f"B proteins with no reciprocal best hit {len(fb) - len(pairs)}; A proteins with none {len(fa) - len(used_a)}"
    )
    n_ident = sum(p[4] for p in pairs)
    out.append(
        f"pairs with identical sequences {n_ident}; non-identical pairs {len(pairs) - n_ident}"
    )
    out.append("")
    calls = sorted({k for d in ca.values() for k in d})
    for label, sel in (
        ("identical sequences", [p for p in pairs if p[4]]),
        ("non-identical pairs", [p for p in pairs if not p[4]]),
    ):
        out.append(f"## {label} (n={len(sel)})")
        bins = collections.Counter(bin_of(p[2], p[3]) for p in sel)
        out.append(
            "identity/coverage bins: " + ", ".join(f"{k}={v}" for k, v in sorted(bins.items()))
        )
        for c in calls:
            t = collections.Counter()
            for sa, qb, *_ in sel:
                va, vb = ca.get(sa, {}).get(c), cb.get(qb, {}).get(c)
                if va is not None and vb is not None:
                    t[(va, vb)] += 1
            n = sum(t.values())
            agree = sum(v for (x, y), v in t.items() if x == y)
            out.append(
                f"  {c}: n={n} agree={agree} kappa={kappa(t):.3f} table={dict(sorted(t.items()))}"
            )
        out.append("")
    out.append(
        "## named cases (A accession, B protein, identity %, coverage, identical, calls A -> B)"
    )
    pair_by_a = {p[0]: p for p in pairs}
    for acc, label in NAMED.items():
        sa = next((k for k in fa if short(k) == acc), None)
        p = pair_by_a.get(sa) if sa else None
        if p is None:
            out.append(
                f"  {label} ({acc}): no reciprocal best hit"
                if sa
                else f"  {label} ({acc}): not in A"
            )
            continue
        diffs = [
            f"{c}: {ca[sa].get(c)} -> {cb.get(p[1], {}).get(c)}"
            for c in calls
            if ca[sa].get(c) != cb.get(p[1], {}).get(c)
        ]
        out.append(
            f"  {label} ({acc}): {p[1]} {p[2]:.1f}% cov {p[3]:.2f} identical={p[4]}; "
            + ("; ".join(diffs) or "all calls equal")
        )
    text = "\n".join(out)
    print(text)
    with open(a.out_prefix + ".summary.txt", "w") as f:
        f.write(text + "\n")
    with open(a.out_prefix + ".pairs.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["a", "b", "pident", "coverage", "identical", "bin"] + calls)
        for sa, qb, pid, cov, ident in pairs:
            w.writerow(
                [sa, qb, f"{pid:.2f}", f"{cov:.3f}", int(ident), bin_of(pid, cov)]
                + [
                    (
                        "="
                        if ca.get(sa, {}).get(c) == cb.get(qb, {}).get(c)
                        else f"{ca.get(sa, {}).get(c)}>{cb.get(qb, {}).get(c)}"
                    )
                    for c in calls
                ]
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
