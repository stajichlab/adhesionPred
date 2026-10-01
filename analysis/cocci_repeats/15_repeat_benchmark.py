#!/usr/bin/python3.12
"""Build a tandem-repeat benchmark with known period and copy number.

Why: `02_repeat_profile.py` matches residues exactly (s[i] == s[i+p]). That must fail once
repeat units diverge, but nobody had measured where. This builds the set that measures it.

Four parts.

1. `sowgp` - real, low divergence. Distinct full-length SOWgp alleles from
   `sowgp_repeat_viz.copies.tsv`. Truth period 47 aa; truth copy number is the anchored unit
   count from `sowgp_units.py` (an integer, validated against the published SOWgp58/66/82
   copy numbers). Only alleles with regular 47 +-1 aa anchor spacing are used.

2. `synthetic` - the divergence axis. A donor unit is tiled N times, then every copy is
   mutated independently to a target identity, and the array is padded with non-repetitive
   flanks. Substitutions follow a BLOSUM62-informed model: residue a is replaced by b != a
   with probability proportional to exp(LAMBDA * BLOSUM62[a][b]). The target identity sets
   the per-site mutation rate; the identity actually realised is measured after generation
   (`true_identity`) and is what the evaluation plots.

   THIS MEASURES THE DETECTOR'S RESPONSE TO DIVERGENCE UNDER AN ASSUMED MUTATION MODEL.
   It is not a claim about how natural repeat arrays evolve. Real arrays diverge by
   concerted evolution, slippage, indels and gene conversion, none of which are modelled
   here. Indels in particular are absent, so the numbers below are an upper bound on
   performance for a given identity.

   Donors: `native` = real repeat units taken from `class2a_candidates.tsv` (periods 14-47).
   `background` = units drawn from the *C. immitis* RS proteome amino-acid composition
   (periods 12-75, so the period axis reaches past what the real set covers).

3. `negative` - proteins with no tandem repeat. Three kinds:
   `random` (background composition only), `lowcomplexity` (a 2-4 letter tract embedded in
   background flanks - the obvious false-positive mode), and `shuffled` (a synthetic
   positive shuffled, so composition is kept and periodicity destroyed).

4. `real_candidate` - the 41 class-2a candidates from `03_repeat_surface_candidates.py`.
   These have NO independent ground truth. They were selected by the detector under test, so
   recall on them is circular. They are carried only to see which calls the new detector
   keeps and which it drops, and are excluded from every rate reported as recall.

Outputs: `repeat_benchmark.fa` and `repeat_benchmark_truth.tsv`.

Usage: ./15_repeat_benchmark.py [--prefix repeat_benchmark] [--seed 20260929]
"""

import argparse
import csv
import math
import random
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import sowgp_units as su  # noqa: E402

sys.path.insert(0, str(HERE.parent / "_common"))
from paths import COCCI_PANGENOME  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"
LAMBDA = 0.5  # BLOSUM62 weight in the substitution model
COPY_NUMBERS = (3, 5, 8, 15)
IDENTITY_TARGETS = tuple(round(1.00 - 0.05 * i, 2) for i in range(15))  # 1.00 .. 0.30
REPLICATES = 4
BACKGROUND_PERIODS = (12, 24, 36, 48, 60, 75)
NATIVE_PERIODS = (14, 17, 20, 34, 39, 47)
N_NEG_RANDOM = 300
N_NEG_LOWCOMPLEXITY = 300
N_NEG_SHUFFLED = 200
LC_ALPHABETS = ("ST", "QN", "PAS", "GS", "SAT", "EK")

_B62 = substitution_matrices.load("BLOSUM62")
REF_FASTA = COCCI_PANGENOME / "input_run2" / "CimmitisRS_FungiDB.fasta"


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


def proteome_background(path):
    c = Counter()
    for _, s in read_fasta(path):
        c.update(x for x in s if x in AA)
    t = sum(c.values())
    return np.array([c[a] / t for a in AA])


def sub_table():
    """P(b | a) for a substitution away from a, proportional to exp(LAMBDA*BLOSUM62)."""
    m = np.zeros((20, 20))
    for i, a in enumerate(AA):
        w = np.array([math.exp(LAMBDA * _B62[a, b]) if b != a else 0.0 for b in AA])
        m[i] = w / w.sum()
    return m


SUB = sub_table()
AAI = {a: i for i, a in enumerate(AA)}


def mutate(unit, rate, rng):
    out = []
    for ch in unit:
        if ch in AAI and rng.random() < rate:
            out.append(rng.choices(AA, weights=SUB[AAI[ch]])[0])
        else:
            out.append(ch)
    return "".join(out)


def mean_pairwise_identity(units):
    n = len(units)
    vals = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = units[i], units[j]
            k = min(len(a), len(b))
            vals.append(sum(x == y for x, y in zip(a[:k], b[:k])) / k)
    return float(np.mean(vals)) if vals else 1.0


def draw(freq, n, rng):
    return "".join(rng.choices(AA, weights=freq, k=n))


def build_synthetic(donors, freq, rng, rows, seqs):
    for donor_class, period, unit in donors:
        for ncop in COPY_NUMBERS:
            for target in IDENTITY_TARGETS:
                # pairwise identity between two independently mutated copies is about
                # (1-rate)^2, so rate = 1 - sqrt(target). The realised value is measured.
                rate = 1.0 - math.sqrt(max(target, 0.0))
                for rep in range(REPLICATES):
                    units = [mutate(unit, rate, rng) for _ in range(ncop)]
                    ident = mean_pairwise_identity(units)
                    # flanks scale with the array, so the true repeat coverage always sits
                    # in about 0.5-0.8 and the 0.25 coverage threshold is reachable for
                    # every cell. With fixed 60-150 aa flanks a 3-copy, 12 aa-period array
                    # could not pass that threshold however well it was detected, which
                    # would have confounded the divergence measurement with array length.
                    alen = ncop * len(unit)
                    lo, hi = max(30, int(0.12 * alen)), min(300, max(60, int(0.45 * alen)))
                    lf = rng.randint(lo, max(lo, hi))
                    rf = rng.randint(lo, max(lo, hi))
                    left, right = draw(freq, lf, rng), draw(freq, rf, rng)
                    seq = left + "".join(units) + right
                    sid = f"syn_{donor_class}_p{period}_n{ncop}_i{int(target * 100):03d}_r{rep}"
                    seqs.append((sid, seq))
                    rows.append(
                        {
                            "seq_id": sid,
                            "part": "synthetic",
                            "donor_class": donor_class,
                            "true_period": period,
                            "true_copies": ncop,
                            "target_identity": target,
                            "true_identity": round(ident, 4),
                            "array_start": lf,
                            "array_end": lf + sum(len(u) for u in units),
                            "length": len(seq),
                            "has_repeat": 1,
                            "note": "",
                        }
                    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--prefix", default=str(HERE / "repeat_benchmark"))
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--copies-tsv", default=str(HERE / "sowgp_repeat_viz.copies.tsv"))
    ap.add_argument("--candidates", default=str(HERE / "class2a_candidates.tsv"))
    args = ap.parse_args()

    rng = random.Random(args.seed)
    rows, seqs = [], []

    print("reading proteome background composition ...", file=sys.stderr)
    freq = proteome_background(REF_FASTA)

    # --- part 1: real SOWgp alleles, truth from the PTDCYGDC anchors -------------------
    H = su.distinct_alleles(args.copies_tsv)
    n_sowgp = 0
    for _, r in H.iterrows():
        if not r.regular:
            continue
        sid = f"sowgp_{r.allele}"
        seqs.append((sid, r.seq))
        rows.append(
            {
                "seq_id": sid,
                "part": "sowgp",
                "donor_class": r.species,
                "true_period": su.PERIOD,
                "true_copies": int(r.n_units),
                "target_identity": "",
                "true_identity": round(mean_pairwise_identity(list(r.unit_seqs)), 4),
                "array_start": int(r.units[0][0]),
                "array_end": int(r.units[-1][0]) + su.PERIOD,
                "length": len(r.seq),
                "has_repeat": 1,
                "note": f"{r.n} strains",
            }
        )
        n_sowgp += 1
    print(f"  sowgp alleles (regular spacing): {n_sowgp}", file=sys.stderr)

    # --- part 2: synthetic -------------------------------------------------------------
    cand = pd.read_csv(args.candidates, sep="\t")
    donors = []
    for p in NATIVE_PERIODS:
        sub = cand[cand.rep_period == p]
        if sub.empty:
            continue
        u = str(sub.iloc[0].rep_unit)
        u = "".join(c if c in AA else rng.choice(AA) for c in u)
        donors.append(("native", len(u), u))
    for p in BACKGROUND_PERIODS:
        donors.append(("background", p, draw(freq, p, rng)))
    print(f"  donors: {[(c, p) for c, p, _ in donors]}", file=sys.stderr)
    build_synthetic(donors, freq, rng, rows, seqs)
    n_syn = sum(1 for r in rows if r["part"] == "synthetic")
    print(f"  synthetic sequences: {n_syn}", file=sys.stderr)

    # --- part 3: negatives -------------------------------------------------------------
    for i in range(N_NEG_RANDOM):
        s = draw(freq, rng.randint(300, 600), rng)
        seqs.append((f"neg_random_{i:03d}", s))
        rows.append(
            {
                "seq_id": f"neg_random_{i:03d}",
                "part": "negative",
                "donor_class": "random",
                "true_period": 0,
                "true_copies": 0,
                "target_identity": "",
                "true_identity": "",
                "array_start": 0,
                "array_end": 0,
                "length": len(s),
                "has_repeat": 0,
                "note": "",
            }
        )
    for i in range(N_NEG_LOWCOMPLEXITY):
        alpha = rng.choice(LC_ALPHABETS)
        tract = "".join(rng.choices(alpha, k=rng.randint(100, 250)))
        s = draw(freq, rng.randint(60, 150), rng) + tract + draw(freq, rng.randint(60, 150), rng)
        seqs.append((f"neg_lowcomplexity_{i:03d}", s))
        rows.append(
            {
                "seq_id": f"neg_lowcomplexity_{i:03d}",
                "part": "negative",
                "donor_class": "lowcomplexity",
                "true_period": 0,
                "true_copies": 0,
                "target_identity": "",
                "true_identity": "",
                "array_start": 0,
                "array_end": 0,
                "length": len(s),
                "has_repeat": 0,
                "note": alpha,
            }
        )
    syn_ids = [s for s in seqs if s[0].startswith("syn_")]
    for i, (_, s) in enumerate(rng.sample(syn_ids, N_NEG_SHUFFLED)):
        lst = list(s)
        rng.shuffle(lst)
        sh = "".join(lst)
        seqs.append((f"neg_shuffled_{i:03d}", sh))
        rows.append(
            {
                "seq_id": f"neg_shuffled_{i:03d}",
                "part": "negative",
                "donor_class": "shuffled",
                "true_period": 0,
                "true_copies": 0,
                "target_identity": "",
                "true_identity": "",
                "array_start": 0,
                "array_end": 0,
                "length": len(sh),
                "has_repeat": 0,
                "note": "",
            }
        )

    # --- part 4: real class-2a candidates, no independent truth ------------------------
    want = {}
    for _, r in cand.iterrows():
        want.setdefault(r.strain, set()).add(r.protein)
    src = {
        "CimmitisRS_FungiDB": COCCI_PANGENOME / "input_run2" / "CimmitisRS_FungiDB.fasta",
        "CposadasiiSilveira2022_FungiDB": COCCI_PANGENOME
        / "input_run2"
        / "CposadasiiSilveira2022_FungiDB.fasta",
    }
    from paths import COCCI_LONGREAD  # noqa: E402

    for d in sorted(Path(COCCI_LONGREAD).iterdir()):
        for fa in d.glob("*.proteins.fa"):
            src[fa.name.replace(".proteins.fa", "")] = fa
    n_cand = 0
    for strain, ids in want.items():
        fa = src.get(strain)
        if not fa or not Path(fa).exists():
            print(f"  WARNING: no fasta for {strain}", file=sys.stderr)
            continue
        for pid, s in read_fasta(fa):
            if pid in ids:
                s = s.rstrip("*")
                sid = f"cand_{strain}_{pid}"
                seqs.append((sid, s))
                rows.append(
                    {
                        "seq_id": sid,
                        "part": "real_candidate",
                        "donor_class": strain,
                        "true_period": "",
                        "true_copies": "",
                        "target_identity": "",
                        "true_identity": "",
                        "array_start": "",
                        "array_end": "",
                        "length": len(s),
                        "has_repeat": "",
                        "note": pid,
                    }
                )
                n_cand += 1
    print(f"  real class-2a candidates: {n_cand}", file=sys.stderr)

    fa_out = f"{args.prefix}.fa"
    tsv_out = f"{args.prefix}_truth.tsv"
    with open(fa_out, "w") as fh:
        for sid, s in seqs:
            fh.write(f">{sid}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    with open(tsv_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(seqs)} sequences to {fa_out}", file=sys.stderr)
    print(f"wrote {len(rows)} truth rows to {tsv_out}", file=sys.stderr)
    print(Counter(r["part"] for r in rows), file=sys.stderr)


if __name__ == "__main__":
    main()
