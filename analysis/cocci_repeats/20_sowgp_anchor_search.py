#!/usr/bin/env python3.12
"""Search proteomes directly for the SOWgp anchor motif, and map each hit onto SOWgp58.

Written for the 2026-09-30 correction to REPORT_2026-09-29_sowgp_repeat_structure.md.

Why this exists: section 5 of that report stated that VFC140, Cpos1038 and Cpos3700 have a
SOWgp locus but NO protein call. That was wrong. The search there went through the 496-genome
pangenome orthogroups, and the five UArizona long-read proteomes were only ever profiled by
02_repeat_profile.py -- never searched for the anchor directly. All three strains do have
SOWgp protein models.

This script does the direct search that was missing:
  1. find every protein containing the PTDCYGDC anchor (the same motif sowgp_units.py uses)
  2. locally align each hit to every sequence in sowgp_seed.fa
  3. report which span of SOWgp58 each hit covers, and at what identity

A hit covering the whole reference is a (possibly short) allele. Two hits covering
OVERLAPPING halves of the reference, with consecutive locus tags, indicate one gene called as
two -- a split gene model, which is the failure mode the original report predicted but did not
observe directly.

Run (from this directory):
    /usr/bin/python3.12 20_sowgp_anchor_search.py
    /usr/bin/python3.12 20_sowgp_anchor_search.py --out sowgp_anchor_hits.tsv

Needs biopython. Site paths come from ../_common/paths.py. Takes about 20 s.
"""

import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "_common"))

ANCHOR = "PTDCYGDC"
REFERENCE = "SOWgp58_Cocci_immitis_published"


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def proteome_paths():
    """The same seven proteomes 01_signalp.sh and 02_repeat_profile.py use."""
    import paths  # noqa: E402  -- from ../_common

    lr = Path(paths.COCCI_LONGREAD)
    pan = Path(paths.COCCI_PANGENOME) / "input_run2"
    out = sorted(lr.glob("*/*.proteins.fa"))
    for f in ("CimmitisRS_FungiDB.fasta", "CposadasiiSilveira2022_FungiDB.fasta"):
        if (pan / f).exists():
            out.append(pan / f)
    return out


def aligner():
    al = Align.PairwiseAligner()
    al.mode = "local"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    return al


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--seed",
        default=str(HERE / "sowgp_seed.fa"),
        help="SOWgp seed FASTA (gunzip sowgp_seed.fa.gz first if absent)",
    )
    ap.add_argument("--out", default=None, help="optional TSV output")
    args = ap.parse_args()

    seed_path = Path(args.seed)
    if not seed_path.exists():
        sys.exit(f"{seed_path} not found -- run: gunzip -k {seed_path}.gz")
    seed = read_fasta(seed_path)
    if REFERENCE not in seed:
        sys.exit(f"{REFERENCE} not in {seed_path}")
    ref = seed[REFERENCE]
    al = aligner()

    print(f"anchor: {ANCHOR}")
    print(f"reference: {REFERENCE}, {len(ref)} aa\n")
    header = f"{'proteome':<32} {'protein':<24} {'len':>5} {'anch':>4} {'maps to SOWgp58':>18} {'id%':>6}"
    print(header)
    print("-" * len(header))

    rows = []
    for path in proteome_paths():
        label = path.name.replace(".proteins.fa", "").replace(".fasta", "")
        hits = [(n, s) for n, s in read_fasta(path).items() if ANCHOR in s]
        if not hits:
            print(f"{label:<32} {'-- no anchored protein --'}")
            continue
        for name, seq in sorted(hits):
            a = al.align(seq, ref)[0]
            rs, re_ = a.aligned[1][0][0], a.aligned[1][-1][1]
            ident = sum(1 for x, y in zip(a[0], a[1], strict=False) if x == y and x != "-")
            alen = sum(1 for x, y in zip(a[0], a[1], strict=False) if x != "-" and y != "-")
            pct = 100 * ident / max(alen, 1)
            span = f"{rs + 1}-{re_}"
            print(
                f"{label:<32} {name:<24} {len(seq):>5} {seq.count(ANCHOR):>4} "
                f"{span:>18} {pct:>6.1f}"
            )
            rows.append((label, name, len(seq), seq.count(ANCHOR), rs + 1, re_, round(pct, 1)))

    print("\nReading the result:")
    print("  one hit spanning most of the reference   -> an allele (short if few anchors)")
    print("  two hits covering OVERLAPPING halves,")
    print("  with consecutive locus tags              -> one gene called as two (split model)")

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("proteome\tprotein\tlength\tn_anchor\tref_start\tref_end\tidentity_pct\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
