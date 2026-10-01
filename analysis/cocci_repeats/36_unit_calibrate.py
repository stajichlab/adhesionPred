#!/usr/bin/python3.12
"""Calibrate the SOWgp unit HMM (34) before it is run on 5,800 proteomes.

Two questions, both answered with hmmsearch on the same HMM the real search uses.

NULL. What is the highest domain score the HMM gives to proteins that cannot carry a SOWgp
unit? Four Onygenales proteomes (Coccidioides RS, Blastomyces ER-3, Histoplasma G217B,
Uncinocarpus 1704) are shuffled per protein (keeps length and composition, destroys order and
homology). Zero shuffled hits at E <= 1000 (chance expectation 0.5) supports calling a hit at E <= 1 at Z = 6e7. It is not a full null: a
per-protein shuffle under-states the chance rate for real convergent sequence (the LAAKIS
lesson in 30_anchor_family_search.py), so hits above the floor still need a look at the
alignment.

SENSITIVITY. How diverged can a unit be and still score above the floor? Arrays of
--copies units are built from the modal unit at fixed identity to it (changed positions drawn
from BLOSUM62-conditional substitution probabilities), embedded in random flanks. Two flank
types: natural composition, and Pro-rich (P = 25%). No indels are simulated, so the floor on
detectable divergence found here is OPTIMISTIC for real diverged units, which also carry
indels.

Needs hmmer/3.4 on PATH and the output of 34, 35. Run:  sbatch 36_unit_calibrate.sh
Writes unit_calibration_null.tsv and unit_calibration_sensitivity.tsv.
"""

import argparse
import math
import os
import random
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import importlib  # noqa: E402

fam = importlib.import_module("30_anchor_family_search")

AA = "ACDEFGHIKLMNPQRSTVWY"
NAT = dict(
    zip(
        AA,
        [
            8.25,
            1.37,
            5.45,
            6.75,
            3.86,
            7.07,
            2.27,
            5.96,
            5.84,
            9.66,
            2.42,
            4.06,
            4.70,
            3.93,
            5.53,
            6.56,
            5.34,
            6.87,
            1.08,
            2.92,
        ],
        strict=False,
    )
)
IDS = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]


def hmmsearch(hmm, fasta, out, cpu=8, z=60_000_000):
    subprocess.run(
        [
            "hmmsearch",
            "--cpu",
            str(cpu),
            "-Z",
            str(z),
            "--domZ",
            str(z),
            "-E",
            "1000",
            "--domE",
            "1000",
            "--noali",
            "--domtblout",
            str(out),
            "-o",
            os.devnull,
            str(hmm),
            str(fasta),
        ],
        check=True,
    )


def read_domtbl(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        rows.append(
            {
                "target": f[0],
                "tlen": int(f[2]),
                "fullscore": float(f[5]),
                "ievalue": float(f[12]),
                "score": float(f[13]),
                "ali_from": int(f[17]),
                "ali_to": int(f[18]),
            }
        )
    return rows


def blosum_sub_table():
    from Bio.Align import substitution_matrices

    B = substitution_matrices.load("BLOSUM62")
    bg = np.array([NAT[a] for a in AA])
    bg /= bg.sum()
    tab = {}
    for a in AA:
        p = np.array(
            [bg[j] * math.exp(0.3466 * B[a][b]) if b != a else 0.0 for j, b in enumerate(AA)]
        )
        tab[a] = p / p.sum()
    return tab, bg


def diverge(unit, ident, rng, tab):
    n = len(unit)
    keep = set(rng.sample(range(n), round(ident * n)))
    out = []
    for i, a in enumerate(unit):
        out.append(a if i in keep else rng.choices(AA, weights=tab[a])[0])
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--hmm", default=str(HERE / "sowgp_unit.hmm"))
    ap.add_argument("--n", type=int, default=300, help="simulated arrays per identity level")
    ap.add_argument("--copies", type=int, default=4)
    ap.add_argument(
        "--floor-e",
        type=float,
        default=1.0,
        help="domain i-evalue (Z = 6e7 proteins) at or below which a unit counts as found",
    )
    args = ap.parse_args()
    scratch = Path(os.environ["SCRATCH"]) / "unitcal"
    scratch.mkdir(parents=True, exist_ok=True)
    rng = random.Random(20260930)

    # ---- null
    gen = {}
    import csv

    for r in csv.DictReader(open(HERE / "unit_genomes.tsv"), delimiter="\t"):
        if r["source"] == "fungi5k" and r["label"] in (
            "Coccidioides_immitis_RS",
            "Blastomyces_dermatitidis_ER-3",
            "Histoplasma_ohiense_G217B",
            "Uncinocarpus_reesii_1704",
        ):
            gen[r["label"]] = r["path"]
    seqs = {}
    for lab, p in gen.items():
        for k, s in fam.read_fasta(p).items():
            seqs[f"{lab}|{k}"] = s
    ids = list(seqs)
    sh = fam.shuffled_copy([seqs[i] for i in ids])
    with open(scratch / "null.fa", "w") as fh:
        for i, s in zip(ids, sh, strict=False):
            fh.write(f">{i}\n{s}\n")
    hmmsearch(args.hmm, scratch / "null.fa", scratch / "null.domtbl")
    nd = read_domtbl(scratch / "null.domtbl")
    sc = sorted((d["score"] for d in nd), reverse=True)
    print(
        f"NULL: {len(ids)} shuffled proteins from {len(gen)} proteomes, {len(nd)} domain hits "
        f"at E<=1000 (Z=6e7 proteins)"
    )
    print("top 10 null domain scores:", [round(x, 1) for x in sc[:10]])
    print(
        "expected by chance at E<=1000: %.2f hits in this many proteins" % (1000 * len(ids) / 6e7)
    )
    with open(HERE / "unit_calibration_null.tsv", "w") as fh:
        fh.write("rank\tscore\n")
        for i, x in enumerate(sc[:200]):
            fh.write(f"{i+1}\t{x}\n")
    floor_e = args.floor_e

    # ---- sensitivity
    modal = "".join(
        ln.strip() for ln in open(HERE / "sowgp_unit_modal.fa") if not ln.startswith(">")
    )
    tab, bg = blosum_sub_table()
    pro = {a: (25.0 if a == "P" else NAT[a] * 75 / (100 - NAT["P"])) for a in AA}
    flank_w = {"natural": [NAT[a] for a in AA], "pro_rich": [pro[a] for a in AA]}
    with open(scratch / "sim.fa", "w") as fh:
        for flank, w in flank_w.items():
            for ident in IDS:
                for j in range(args.n):
                    units = [diverge(modal, ident, rng, tab) for _ in range(args.copies)]
                    f1 = "".join(rng.choices(AA, weights=w, k=rng.randint(40, 120)))
                    f2 = "".join(rng.choices(AA, weights=w, k=rng.randint(40, 120)))
                    fh.write(f">{flank}|{ident}|{j}\n{f1}{''.join(units)}{f2}\n")
        # true-negative controls with no unit at all
        for flank, w in flank_w.items():
            for j in range(args.n):
                s = "".join(rng.choices(AA, weights=w, k=args.copies * 47 + 160))
                fh.write(f">{flank}|none|{j}\n{s}\n")
    hmmsearch(args.hmm, scratch / "sim.fa", scratch / "sim.domtbl")
    by = {}
    for d in read_domtbl(scratch / "sim.domtbl"):
        by.setdefault(d["target"], []).append((d["score"], d["ievalue"]))
    print(
        f"\nSENSITIVITY: {args.copies} copies per array, {args.n} arrays per level; floor = domain "
        f"E <= {floor_e}"
    )
    hdr = f"{'flank':<9}{'identity':>9}{'arrays':>8}{'>=1 unit':>10}{'all copies':>12}{'median n':>10}{'median best':>12}"
    print(hdr)
    out = [
        "flank\tidentity\tarrays\tfrac_ge1_unit_above_floor\tfrac_all_copies_found\tmedian_units_found\tmedian_best_score"
    ]
    for flank in flank_w:
        for ident in IDS + ["none"]:
            n1 = nall = 0
            ns, best = [], []
            for j in range(args.n):
                dd = by.get(f"{flank}|{ident}|{j}", [])
                sc_ = [x for x, _ in dd]
                n_above = sum(1 for _, e in dd if e <= floor_e)
                n1 += n_above >= 1
                nall += n_above >= args.copies
                ns.append(n_above)
                best.append(max(sc_) if sc_ else 0.0)
            line = (
                f"{flank:<9}{str(ident):>9}{args.n:>8}{n1/args.n:>10.3f}{nall/args.n:>12.3f}"
                f"{np.median(ns):>10.1f}{np.median(best):>12.1f}"
            )
            print(line)
            out.append(
                f"{flank}\t{ident}\t{args.n}\t{n1/args.n:.3f}\t{nall/args.n:.3f}\t"
                f"{np.median(ns):.1f}\t{np.median(best):.1f}"
            )
    (HERE / "unit_calibration_sensitivity.tsv").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
