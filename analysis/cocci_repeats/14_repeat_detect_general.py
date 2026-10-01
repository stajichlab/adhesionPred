#!/usr/bin/python3.12
"""Family-agnostic tandem-repeat detector. Successor to 02_repeat_profile.py.

02_repeat_profile.py scores a period p by the fraction of positions with s[i] == s[i+p].
Exact identity fails once repeat units diverge, and the region edges it calls are
conservative, so its fractional n_copies runs low (README section 4).

This script changes three things.

1. SIMILARITY, NOT IDENTITY. A position counts as a match when BLOSUM62(s[i], s[i+p]) >= 1,
   i.e. identical or a conservative substitution.

2. A COMPOSITION-FREE SIGNIFICANCE TEST. Similarity matching raises the false-positive rate
   on low-complexity tracts: measured on the benchmark, a plain similarity detector called a
   repeat in 19.0% of random 2-4 letter tracts, against 0.7% for exact matching.
   Correcting the match rate for the sequence's own composition was tried first and was not
   enough (10.7%), and correcting it for the composition of the called region cut the
   controls but also threw out the real Ser/Thr-rich FLO/ALS-type candidates - a Ser/Thr
   adhesin repeat and a random Ser/Thr tract have the same composition, so no
   composition-based cut can separate them.
   What separates them is that a real array matches at ONE period while a random tract
   matches equally at every period. The statistic used is therefore a z-score: the match rate
   at the chosen period against the match rate at all periods that are not its multiples or
   divisors. It is applied twice, over the whole protein and again inside the called region,
   and both must reach Z_MIN. A composition-corrected score is still computed and reported,
   but it does not gate the call.

3. INTEGER COPY COUNTS WITHOUT A MOTIF ANCHOR. sowgp_units.py counts copies by anchoring on
   the literal motif PTDCYGDC, which only works for SOWgp. Here a consensus unit is derived
   from the detected period, turned into a log-odds PSSM, and scanned back over the whole
   protein. Unit starts are chained at ~p spacing and accepted on the PSSM score. The PSSM is
   rebuilt from the accepted units and the scan repeated. The number of accepted starts is an
   integer copy count, and the array can grow past the edges the periodicity track called.

Also reported: `z_seq`, `z_region`, `region_score`, and `unit_period`, a unit-of-units period (SOWgp has identical units at positions
1/3/4 of the posadasii 5-unit allele), and `region_entropy`, the Wootton-Federhen complexity
of the repeat region in bits, for filtering low-complexity calls by hand.

Output columns are a superset of 02_repeat_profile.py's, so 03_repeat_surface_candidates.py
can read this file unchanged.

Usage:
  ./14_repeat_detect_general.py <fasta...> --out <tsv>            # similarity mode (default)
  ./14_repeat_detect_general.py <fasta...> --out <tsv> --mode exact
"""

import argparse
import csv
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from Bio.Align import substitution_matrices

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
NAA = len(AA)
OTHER = NAA  # slot 20: X, B, Z, U, * and anything else. Never matches.

MIN_PERIOD, MAX_PERIOD = 4, 80
WINDOW = 3  # smoothing half-window, as in 02
MIN_SCORE = 0.15  # cheap pre-filter on the background-corrected periodicity
Z_MIN = 4.0  # the real test: the best period must beat the other periods by 4 SD
REGION_THRESHOLD = 0.5  # background-corrected per-position track cut for the region
PERIOD_DIVISOR_FRAC = 0.9  # a divisor of the best period replaces it at >= this fraction
PSSM_ACCEPT_FRAC = 0.25  # accept a unit at >= this fraction of the median seed-unit score
PSSM_ROUNDS = 2
PSEUDOCOUNT = 1.0  # Dirichlet weight, spread over the background composition

# --- substitution matrix -> 21x21 boolean "similar" table -------------------------------
_B62 = substitution_matrices.load("BLOSUM62")
SIM = np.zeros((NAA + 1, NAA + 1), dtype=bool)
for _a in AA:
    for _b in AA:
        SIM[AA_IDX[_a], AA_IDX[_b]] = _B62[_a, _b] >= 1
EXACT = np.zeros((NAA + 1, NAA + 1), dtype=bool)
EXACT[np.arange(NAA), np.arange(NAA)] = True

EMPTY = {
    "period": 0,
    "score": 0.0,
    "start": 0,
    "end": 0,
    "n_copies": 0,
    "coverage": 0.0,
    "unit": "",
    "n_units": 0,
    "unit_period": 0,
    "unit_identity": 0.0,
    "region_entropy": 0.0,
    "raw_score": 0.0,
    "expected": 0.0,
    "region_score": 0.0,
    "z_seq": 0.0,
    "z_region": 0.0,
}


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


def encode(s):
    return np.array([AA_IDX.get(c, OTHER) for c in s], dtype=np.int8)


def background(idx):
    """Amino-acid frequencies of this sequence, over the 20 standard residues."""
    f = np.bincount(idx[idx < NAA], minlength=NAA).astype(float)
    t = f.sum()
    return f / t if t else np.full(NAA, 1.0 / NAA)


def expected_rate(freq, table):
    """Chance rate of a `table` match between two residues drawn from `freq`."""
    return float(freq @ table[:NAA, :NAA].astype(float) @ freq)


def wf_entropy(s):
    """Wootton-Federhen complexity of a string, in bits (max log2(20) = 4.32)."""
    n = len(s)
    if n == 0:
        return 0.0
    c = Counter(s)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


# --- stage 1: periodicity ---------------------------------------------------------------


def period_scan(idx, table, exp):
    """Background-corrected periodicity score for every candidate period."""
    n = len(idx)
    hi = min(MAX_PERIOD, n // 3)
    out = {}
    denom = 1.0 - exp
    for p in range(MIN_PERIOD, hi + 1):
        raw = float(table[idx[:-p], idx[p:]].mean())
        out[p] = ((raw - exp) / denom if denom > 0 else 0.0, raw)
    return out


def raw_scan(idx, table):
    """Raw match rate at every candidate period."""
    hi = min(MAX_PERIOD, len(idx) // 3)
    return {p: float(table[idx[:-p], idx[p:]].mean()) for p in range(MIN_PERIOD, hi + 1)}


def period_z(raw, p):
    """How far the match rate at p stands above the match rate at unrelated periods.

    This is the statistic that decides whether a call is a repeat. The composition-corrected
    score alone is not enough: a Ser/Thr-rich adhesin repeat and a random Ser/Thr tract have
    the same composition, so any composition-based cut either keeps both or drops both. They
    differ in that the adhesin repeat matches at ONE period and the random tract matches
    equally at every period. Periods that are multiples or divisors of p are excluded from
    the comparison set, because a real period p also raises them.

    Returns None when too few unrelated periods are left to estimate a spread.
    """
    rel = [v for q, v in raw.items() if not (q % p in (0, 1) or (q + 1) % p == 0 or p % q == 0)]
    if len(rel) < 5:
        return None
    mu, sd = float(np.mean(rel)), float(np.std(rel))
    if sd < 0.01:
        return 0.0  # every period matches equally: a low-complexity tract, not a repeat
    return (raw[p] - mu) / sd


def z_scan(raw):
    return {p: (period_z(raw, p) or 0.0) for p in raw}


def pick_period(scores):
    """Best period, replaced by its smallest well-scoring divisor.

    A true period p also scores highly at 2p, 3p... The argmax can land on a multiple.
    Only proper divisors of the argmax are allowed to replace it, so an unrelated small
    period with an accidentally high score cannot win.
    """
    if not scores:
        return 0, 0.0, 0.0
    best_p = max(scores, key=lambda p: scores[p][0])
    best = scores[best_p][0]
    for d in range(MIN_PERIOD, best_p):
        if best_p % d == 0 and d in scores and scores[d][0] >= PERIOD_DIVISOR_FRAC * best:
            best_p = d
            break
    return best_p, scores[best_p][0], scores[best_p][1]


def periodic_region(idx, p, table, exp):
    """Longest run where the smoothed, background-corrected match track stays high."""
    m = table[idx[:-p], idx[p:]].astype(float)
    k = 2 * WINDOW + 1
    csum = np.concatenate([[0.0], np.cumsum(m)])
    lo = np.maximum(0, np.arange(len(m)) - WINDOW)
    hi = np.minimum(len(m), np.arange(len(m)) + WINDOW + 1)
    sm = (csum[hi] - csum[lo]) / (hi - lo)
    del k
    denom = 1.0 - exp
    norm = (sm - exp) / denom if denom > 0 else sm * 0.0
    inside, start, span = False, 0, (0, 0)
    for i, v in enumerate(list(norm) + [-1.0]):
        if v >= REGION_THRESHOLD and not inside:
            inside, start = True, i
        elif v < REGION_THRESHOLD and inside:
            inside = False
            if i - start > span[1] - span[0]:
                span = (start, i)
    return span


# --- stage 2: consensus PSSM and integer copy count -------------------------------------


def build_pssm(idx, starts, p, freq):
    """Log-odds PSSM over p columns from the residues at the given unit starts."""
    counts = np.zeros((p, NAA))
    for st in starts:
        w = idx[st : st + p]
        for j, a in enumerate(w):
            if a < NAA:
                counts[j, a] += 1
    bg = np.where(freq > 0, freq, 1e-6)
    prob = (counts + PSEUDOCOUNT * bg) / (counts.sum(axis=1, keepdims=True) + PSEUDOCOUNT)
    pssm = np.log(prob / bg)
    return np.hstack([pssm, np.zeros((p, 1))])  # column for OTHER scores 0


def pssm_track(idx, pssm, p):
    """PSSM score of every window of length p."""
    n = len(idx) - p + 1
    if n <= 0:
        return np.zeros(0)
    win = np.lib.stride_tricks.sliding_window_view(idx, p)
    return pssm[np.arange(p)[None, :], win].sum(axis=1)


def chain_units(track, p, seed_starts, thr):
    """Walk forwards and backwards at ~p spacing, taking the best window in each step."""
    slack = max(2, int(round(0.25 * p)))
    anchor = int(seed_starts[int(np.argmax(track[seed_starts]))])
    picks = [anchor]
    for direction in (1, -1):
        cur = anchor
        while True:
            lo = cur + direction * p - slack
            hi = cur + direction * p + slack + 1
            lo, hi = max(0, lo), min(len(track), hi)
            if hi - lo <= 0:
                break
            j = lo + int(np.argmax(track[lo:hi]))
            if track[j] < thr or abs(j - cur) < p // 2:
                break
            picks.append(j)
            cur = j
    return sorted(set(picks))


def _tile(a, p, ext, n):
    t = [a + i * p for i in range(max(1, ext // p)) if 0 <= a + i * p <= n - p]
    if not t:
        t = [a] if 0 <= a <= n - p else [max(0, n - p)]
    return t


def _round(idx, p, starts, freq):
    """One PSSM build + scan + chain. Returns (starts, total score, ok)."""
    pssm = build_pssm(idx, starts, p, freq)
    track = pssm_track(idx, pssm, p)
    if track.size == 0:
        return starts, 0.0, False
    seed = [s for s in starts if s < track.size] or [int(np.argmax(track))]
    med = float(np.median(track[seed]))
    if med <= 0:
        return starts, 0.0, False
    new = chain_units(track, p, seed, PSSM_ACCEPT_FRAC * med)
    return new, float(track[new].sum()), True


def refine_units(idx, p, span, freq):
    """Integer unit starts from a consensus PSSM, seeded on the periodicity region.

    The tiling phase is searched. The periodicity track puts the region edge wherever the
    smoothed match rate crosses the cut, which need not be a unit boundary. A phase that is
    out of register splits the first and last units across two windows, so a real terminal
    copy scores badly and is dropped. All p phases are tried and the one giving the most
    accepted units is kept (ties broken on total score).
    """
    a, b = span
    ext = (b - a) + p
    n = len(idx)
    best = None
    for delta in range(p):
        st = (
            _tile(a + delta - p, p, ext + p, n)
            if a + delta - p >= 0
            else _tile(a + delta, p, ext, n)
        )
        cand, sc, ok = _round(idx, p, st, freq)
        if not ok:
            continue
        key = (len(cand), sc)
        if best is None or key > best[0]:
            best = (key, cand)
    starts = best[1] if best else _tile(a, p, ext, n)
    for _ in range(PSSM_ROUNDS):
        pssm = build_pssm(idx, starts, p, freq)
        track = pssm_track(idx, pssm, p)
        if track.size == 0:
            break
        seed = [s for s in starts if s < track.size] or [int(np.argmax(track))]
        seed_scores = np.sort(track[seed])
        med = float(np.median(seed_scores))
        thr = PSSM_ACCEPT_FRAC * med
        if med <= 0:
            return starts, 0.0
        new = chain_units(track, p, seed, thr)
        if new == starts:
            break
        starts = new
    pssm = build_pssm(idx, starts, p, freq)
    track = pssm_track(idx, pssm, p)
    mean_score = float(np.mean(track[[s for s in starts if s < track.size]])) if track.size else 0.0
    return starts, mean_score


def unit_level_period(seq, starts, p):
    """Smallest lag k where units i and i+k are more alike than units in general.

    Detects a unit-of-units: SOWgp's posadasii 5-unit allele has identical units at 1/3/4.
    """
    units = [seq[s : s + p] for s in starts if s + p <= len(seq)]
    n = len(units)
    if n < 3:
        return 0, 0.0

    def ident(a, b):
        return sum(x == y for x, y in zip(a, b, strict=False)) / max(len(a), len(b), 1)

    allpairs = [ident(units[i], units[j]) for i in range(n) for j in range(i + 1, n)]
    base = float(np.mean(allpairs)) if allpairs else 0.0
    best_k, best_v = 0, 0.0
    for k in range(1, n // 2 + 1):
        v = float(np.mean([ident(units[i], units[i + k]) for i in range(n - k)]))
        if v > best_v:
            best_k, best_v = k, v
    if best_v <= base + 0.05 or best_v < 0.8:
        return 0, round(base, 3)
    return best_k, round(best_v, 3)


# --- the detector -----------------------------------------------------------------------


def detect(seq, mode="similarity"):
    s = seq.rstrip("*")
    if len(s) < MIN_PERIOD * 3:
        return dict(EMPTY)
    table = SIM if mode == "similarity" else EXACT
    idx = encode(s)
    freq = background(idx)
    exp = expected_rate(freq, table)
    scores = period_scan(idx, table, exp)
    if not scores:
        return dict(EMPTY)
    raw_all = {q: v[1] for q, v in scores.items()}
    # rank periods by the z-score, not by the composition-corrected rate: the z-score is the
    # statistic the call is made on, so the period chosen should be the one that maximises it
    zs = z_scan(raw_all)
    p, _, _ = pick_period({q: (zs[q], raw_all[q]) for q in scores})
    sc, raw = scores[p][0], scores[p][1]
    zseq = zs[p]
    if p == 0 or zseq < Z_MIN or sc < MIN_SCORE:
        out = dict(EMPTY)
        out["score"] = round(sc, 3)
        out["expected"] = round(exp, 3)
        out["z_seq"] = round(zseq, 2)
        return out
    span = periodic_region(idx, p, table, exp)
    if span == (0, 0):
        out = dict(EMPTY)
        out["score"] = round(sc, 3)
        out["expected"] = round(exp, 3)
        out["z_seq"] = round(zseq, 2)
        return out
    starts, _ = refine_units(idx, p, span, freq)
    if not starts:
        return dict(EMPTY)
    a = starts[0]
    end = min(len(s), starts[-1] + p)
    ext = end - a
    # Re-test inside the called region. The whole-protein test can be carried by a periodic
    # stretch that the region call then fails to reproduce, and a low-complexity tract sits
    # inside a protein of ordinary composition, so its own statistics are the honest ones.
    sub = idx[a:end]
    rraw = raw_scan(sub, table)
    zreg = period_z(rraw, p) if p in rraw else None
    if zreg is None:
        zreg = zseq  # region too short to hold enough periods; fall back to the whole protein
    fr = background(sub)
    er = expected_rate(fr, table)
    rscore = 0.0
    if ext > p and er < 1.0:
        rscore = (float(table[sub[:-p], sub[p:]].mean()) - er) / (1.0 - er)
    if zreg < Z_MIN:
        out = dict(EMPTY)
        out["score"] = round(sc, 3)
        out["expected"] = round(exp, 3)
        out["region_score"] = round(rscore, 3)
        out["z_seq"] = round(zseq, 2)
        out["z_region"] = round(zreg, 2)
        return out
    uk, uid = unit_level_period(s, starts, p)
    return {
        "period": p,
        "score": round(sc, 3),
        "start": a,
        "end": end,
        "n_copies": len(starts),
        "coverage": round(ext / len(s), 3),
        "unit": s[a : a + p],
        "n_units": len(starts),
        "unit_period": uk,
        "unit_identity": uid,
        "region_entropy": round(wf_entropy(s[a:end]), 3),
        "raw_score": round(raw, 3),
        "expected": round(exp, 3),
        "region_score": round(rscore, 3),
        "z_seq": round(zseq, 2),
        "z_region": round(zreg, 2),
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
    ap.add_argument("--mode", choices=["similarity", "exact"], default="similarity")
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
            r = detect(seq, mode=args.mode)
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
    print(f"{len(rep)} proteins with a substantial tandem repeat", file=sys.stderr)


if __name__ == "__main__":
    main()
