#!/usr/bin/python3.12
"""Build the protein sets for the Pfam confirmation test of the 14 detector calls.

Four sets, all drawn from the same seven proteomes 02 and 14 were run on:

  gained      called by 14 but not by 02 (123)
  shared      called by both (66)             <- the control for "are gained calls worse?"
  lost        called by 02 but not by 14 (17)
  background  called by neither, length-matched to the gained set   <- the null

"Called" means the 03_repeat_surface_candidates.py thresholds: rep_period > 0,
coverage >= 0.25, copies >= 2.5. This is the same rule 18_repeat_real_compare.py uses, so
the gained/shared/lost counts here reproduce that script exactly.

The background is drawn per gained protein: --bg-per-call proteins from the SAME strain
whose length is within --bg-tol of that gained protein's length, sampled without
replacement from the proteins neither detector called. Length matching matters because
long proteins carry more Pfam domains of every kind, and the gained set is long-biased.

Outputs `pfam_sets.fa` (one record per protein, id = strain|protein) and `pfam_sets.tsv`
(set label, length, and both detectors' period/copies/coverage).

Usage: /usr/bin/python3.12 21_pfam_sets.py [--bg-per-call 10] [--seed 1]
"""

import argparse
import random
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
COV, COP = 0.25, 2.5

LR = Path(
    "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/" "UArizona_strains/For_Marc"
)
PAN = Path(
    "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/"
    "2025_All_Cocci/Pangenome/input_run2"
)

# strain label (as written in the detector TSVs) -> proteome fasta
PROTEOMES = {
    "Coccidioides_immitis_CiB10637": LR / "CiB10637/Coccidioides_immitis_CiB10637.proteins.fa",
    "Coccidioides_immitis_CiB10992": LR / "CiB10992/Coccidioides_immitis_CiB10992.proteins.fa",
    "Coccidioides_immitis_VFC140": LR / "VFC140/Coccidioides_immitis_VFC140.proteins.fa",
    "Coccidioides_posadasii_Cpos1038": LR / "Cpos1038/Coccidioides_posadasii_Cpos1038.proteins.fa",
    "Coccidioides_posadasii_Cpos3700": LR / "Cpos3700/Coccidioides_posadasii_Cpos3700.proteins.fa",
    "CimmitisRS_FungiDB": PAN / "CimmitisRS_FungiDB.fasta",
    "CposadasiiSilveira2022_FungiDB": PAN / "CposadasiiSilveira2022_FungiDB.fasta",
}


def read_fasta(path):
    """id -> sequence. The id is the first whitespace-delimited token, as 14 records it."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def load(paths):
    d = pd.concat([pd.read_csv(p, sep="\t") for p in paths], ignore_index=True)
    d["key"] = d.strain.astype(str) + "|" + d.protein.astype(str)
    return d


def called(d):
    return set(d[(d.rep_coverage >= COV) & (d.rep_n_copies >= COP) & (d.rep_period > 0)].key)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--bg-per-call", type=int, default=10, help="background proteins sampled per gained protein"
    )
    ap.add_argument(
        "--bg-tol",
        type=float,
        default=0.15,
        help="background length must be within this fraction of the gained length",
    )
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-fa", default=str(HERE / "pfam_sets.fa"))
    ap.add_argument("--out-tsv", default=str(HERE / "pfam_sets.tsv"))
    args = ap.parse_args()
    rng = random.Random(args.seed)

    old = load([HERE / "repeat_profile_longread.tsv", HERE / "repeat_profile_reference.tsv"])
    new = load([HERE / "repeat_general_longread.tsv", HERE / "repeat_general_reference.tsv"])
    both = set(old.key) & set(new.key)
    o, n = called(old) & both, called(new) & both
    label = {}
    for k in o & n:
        label[k] = "shared"
    for k in n - o:
        label[k] = "gained"
    for k in o - n:
        label[k] = "lost"
    print(
        f"shared {sum(v == 'shared' for v in label.values())}  "
        f"gained {sum(v == 'gained' for v in label.values())}  "
        f"lost {sum(v == 'lost' for v in label.values())}",
        file=sys.stderr,
    )

    # per-protein length, from the new detector table (it profiled every protein)
    meta = new.set_index("key")[["strain", "protein", "length"]]
    newi = new.set_index("key")
    oldi = old.set_index("key")

    # background pool: everything neither detector called
    uncalled = sorted(both - o - n)
    print(f"uncalled pool {len(uncalled)}", file=sys.stderr)
    by_strain = {}
    for k in uncalled:
        by_strain.setdefault(meta.at[k, "strain"], []).append(k)

    bg, used = [], set()
    for k in sorted(k for k, v in label.items() if v == "gained"):
        L = meta.at[k, "length"]
        lo, hi = L * (1 - args.bg_tol), L * (1 + args.bg_tol)
        pool = [
            c
            for c in by_strain.get(meta.at[k, "strain"], [])
            if c not in used and lo <= meta.at[c, "length"] <= hi
        ]
        rng.shuffle(pool)
        take = pool[: args.bg_per_call]
        used.update(take)
        bg.extend(take)
    print(
        f"background drawn {len(bg)} "
        f"(target {args.bg_per_call} x {sum(v == 'gained' for v in label.values())})",
        file=sys.stderr,
    )
    for k in bg:
        label[k] = "background"

    fastas = {s: read_fasta(p) for s, p in PROTEOMES.items()}
    rows, missing = [], 0
    with open(args.out_fa, "w") as fh:
        for k in sorted(label):
            strain, prot = k.split("|", 1)
            seq = fastas[strain].get(prot)
            if seq is None:
                missing += 1
                continue
            seq = seq.rstrip("*")
            fh.write(f">{k}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
            rows.append(
                {
                    "key": k,
                    "strain": strain,
                    "protein": prot,
                    "set": label[k],
                    "length": len(seq),
                    "period_old": oldi.at[k, "rep_period"] if k in oldi.index else 0,
                    "copies_old": oldi.at[k, "rep_n_copies"] if k in oldi.index else 0,
                    "coverage_old": oldi.at[k, "rep_coverage"] if k in oldi.index else 0,
                    "period_new": newi.at[k, "rep_period"],
                    "copies_new": newi.at[k, "rep_n_copies"],
                    "coverage_new": newi.at[k, "rep_coverage"],
                    "region_entropy_new": newi.at[k, "rep_region_entropy"],
                    "unit_new": newi.at[k, "rep_unit"],
                }
            )
    if missing:
        print(f"WARNING: {missing} proteins not found in the proteome fastas", file=sys.stderr)
    df = pd.DataFrame(rows)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(
        df.groupby("set").agg(n=("key", "size"), median_len=("length", "median")), file=sys.stderr
    )
    print(f"wrote {args.out_fa} and {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
