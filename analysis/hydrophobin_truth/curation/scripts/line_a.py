# ruff: noqa
"""Line A (doublet architecture) per the curation rubric.

Rule implemented: P = positions of all cysteines. For each window W of k consecutive cysteines (k = 8, 9, 10),
every choice of 8 cysteines that keeps W[0] as C1 and W[-1] as C8 (6 chosen from the inner k-2) is tested.
A choice qualifies if gap(C2,C3) = 0 and gap(C6,C7) = 0 (gap = residues between). For k = 8 this is just the
8 consecutive cysteines. Line A = yes if any choice qualifies AND r0 == called AND length <= 300.
Usage: line_a.py <q.faa> <evidence_sheets.tsv | '-'>  ('-' = test mode: r0 assumed called)
"""

import csv
import itertools
import sys


def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def doublet_windows(seq):
    """Return list of (k, positions_1based_of_8, gaps) for every qualifying 8-cysteine choice."""
    P = [i for i, a in enumerate(seq.upper()) if a == "C"]
    hits = []
    for k in (8, 9, 10):
        for i in range(0, len(P) - k + 1):
            W = P[i : i + k]
            inner = W[1:-1]
            for mid in itertools.combinations(inner, 6):
                c = [W[0], *mid, W[-1]]
                gaps = [b - a - 1 for a, b in zip(c, c[1:])]
                if gaps[1] == 0 and gaps[5] == 0:
                    hits.append((k, [x + 1 for x in c], gaps))
    return P, hits


def describe(seq):
    P, hits = doublet_windows(seq)
    allgaps = [b - a - 1 for a, b in zip(P, P[1:])]
    ncys = len(P)
    if not hits:
        return (
            False,
            f"n_cys={ncys}; all C-C gaps={','.join(map(str, allgaps))}; windows k=8..10 tested; none with C2C3=0 and C6C7=0",
        )
    # report the first (smallest k, earliest) hit, and the number of qualifying choices
    k, pos, gaps = hits[0]
    ks = sorted({h[0] for h in hits})
    return True, (
        f"n_cys={ncys}; first qualifying window k={k}: C1-C8 at {','.join(map(str, pos))}; gaps={','.join(map(str, gaps))}; "
        f"qualifying 8-cys choices={len(hits)} (window sizes {','.join(map(str, ks))})"
    )


if __name__ == "__main__":
    seqs = read_fasta(sys.argv[1])
    sheet = {}
    if sys.argv[2] != "-":
        for r in csv.DictReader(open(sys.argv[2]), delimiter="\t"):
            sheet[f"{r['proteome']}__{r['id']}"] = r
    print("key\tline_A\tA_detail")
    for name, s in seqs.items():
        motif, det = describe(s)
        r0 = sheet[name]["r0"] if name in sheet else "called(assumed, test mode)"
        L = len(s)
        ok = motif and r0.startswith("called") and L <= 300
        det = f"{det}; r0={r0}; length={L}"
        if motif and not ok:
            reasons = []
            if not r0.startswith("called"):
                reasons.append("r0 not called")
            if L > 300:
                reasons.append("length>300")
            det += "; motif found but fails: " + ",".join(reasons)
        print(f"{name}\t{'yes' if ok else 'no'}\t{det}")
