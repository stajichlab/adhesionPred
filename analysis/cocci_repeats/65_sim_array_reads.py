#!/usr/bin/python3.12
"""Simulated reads from SOWgp alleles with a KNOWN number of repeat units.

Tests the u/4 depth model of REPORT_2026-10-01_sowgp_depth_vs_repeats.md, section 4.3:
if reads from an allele with u units are mapped to RS (4 units), is array/flank depth u/4?

Subcommands
  build  Write the allele windows (FASTA) and a truth table.
         (a) RS edits: the RS window GG704914:950,001-990,000 with one or two 197 nt genomic
             periods deleted or duplicated. A period is one 141 nt exon (exon 4, 5 or 6 of
             CIMG_04613-t26_1) plus the 56 nt intron after it. Removing or copying a whole
             period keeps the reading frame and the splice structure. Each edit is checked:
             the spliced CDS is translated and the PTDCYGDC anchors are counted (must equal u).
         (b) Real alleles from long-read assemblies, +/- 20 kb around the locus:
             Silveira 2022 (posadasii, 4 units), CiB10637 (immitis, 5), Cpos1038 (posadasii, 5),
             VFC140 (2). Unit counts are from longread_locus_summary.tsv (DNA units) and the
             2026-09-29 report (Silveira QVM09276, 4 anchored units).
  sim    Simulate paired reads from one allele: uniform fragment start, normal fragment
         length, read length L, substitution errors at a fixed rate, quality 30.

Run by 65_sim_array_depth.sh (SLURM). Use /usr/bin/python3.12.
"""

import argparse
import gzip
import random
import sys

RS_FA = (
    "/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/"
    "Genotyping/all_C_immitis_ref_RS/genome/FungiDB-68_CimmitisRS_Genome.fasta"
)
SILV_FA = (
    "/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/"
    "Genotyping/all_C_immitis_ref_RS/genome/FungiDB-68_CposadasiiSilveira2022_Genome.fasta"
)
LR = "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc"

WIN = (950001, 990000)  # RS window, 1-based inclusive
# CIMG_04613-t26_1 CDS segments, 1-based inclusive (GFF)
CDS = [
    (970129, 970235),
    (970316, 970423),
    (970477, 970599),
    (970656, 970796),
    (970853, 970993),
    (971050, 971190),
    (971247, 971460),
]
# 197 nt genomic periods: exon k start .. next exon start - 1  (k = exon 4, 5, 6)
PERIODS = {"A": (970656, 970852), "B": (970853, 971049), "C": (971050, 971246)}

# real alleles: (id, fasta, contig, centre, species, units)
REAL = [
    ("real_Silveira2022_u4", SILV_FA, "CP075069", 4603200, "posadasii", 4),
    (
        "real_CiB10637_u5",
        f"{LR}/CiB10637/Coccidioides_immitis_CiB10637.scaffolds.fa",
        "scaffold_2",
        4398448,
        "immitis",
        5,
    ),
    (
        "real_Cpos1038_u5",
        f"{LR}/Cpos1038/Coccidioides_posadasii_Cpos1038.scaffolds.fa",
        "scaffold_2",
        2243000,
        "posadasii",
        5,
    ),
    (
        "real_VFC140_u2",
        f"{LR}/VFC140/Coccidioides_immitis_VFC140.scaffolds.fa",
        "scaffold_2",
        4654525,
        "VFC140",
        2,
    ),
]
CODON = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_contig(path, name):
    seq, keep = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep:
                    break
                keep = line[1:].split()[0] == name
            elif keep:
                seq.append(line.strip())
    if not seq:
        sys.exit(f"contig {name} not found in {path}")
    return "".join(seq).upper()


def translate(s):
    return "".join(CODON.get(s[i : i + 3], "X") for i in range(0, len(s) - 2, 3))


def edit_rs(contig, ops):
    """ops: list of ('del'|'dup', period). Returns (window seq, spliced CDS)."""
    # work on a list of genomic pieces in RS coordinates so CDS can be re-spliced
    lo, hi = WIN
    order = []  # list of (start, end) RS intervals in output order
    cur = lo
    events = sorted(ops, key=lambda o: PERIODS[o[1]][0])
    for kind, p in events:
        s, e = PERIODS[p]
        if kind == "del":
            order.append((cur, s - 1))
            cur = e + 1
        else:  # dup: keep the period, then a copy of it
            order.append((cur, e))
            order.append((s, e))
            cur = e + 1
    order.append((cur, hi))
    win = "".join(contig[s - 1 : e] for s, e in order)
    cds = ""
    for s, e in order:
        for cs, ce in CDS:
            a, b = max(s, cs), min(e, ce)
            if a <= b:
                cds += contig[a - 1 : b]
    return win, cds


def build(a):
    rs = read_contig(RS_FA, "GG704914")
    designs = {
        2: [[("del", "A"), ("del", "B")], [("del", "B"), ("del", "C")]],
        3: [[("del", "A")], [("del", "B")], [("del", "C")]],
        4: [[]],
        5: [[("dup", "A")], [("dup", "B")], [("dup", "C")]],
        6: [
            [("dup", "A"), ("dup", "B")],
            [("dup", "B"), ("dup", "C")],
            [("dup", "A"), ("dup", "C")],
        ],
    }
    rows = []
    with open(a.out + ".fa", "w") as fa:
        for u, dl in designs.items():
            for ops in dl:
                win, cds = edit_rs(rs, ops)
                prot = translate(cds)
                n_anchor = prot.count("PTDCYGDC")
                tag = "_".join(f"{k}{p}" for k, p in ops) or "RS"
                aid = f"rsedit_u{u}_{tag}"
                ok = n_anchor == u and prot.endswith("*") and prot.count("*") == 1
                print(
                    f"{aid}: window {len(win)} nt, CDS {len(cds)} nt, protein {len(prot) - 1} aa, "
                    f"anchors {n_anchor}, valid={ok}",
                    file=sys.stderr,
                )
                if not ok:
                    sys.exit(f"edit {aid} failed validation")
                fa.write(f">{aid}\n{win}\n")
                rows.append((aid, "RS edit", "immitis", u, len(win)))
        for aid, path, ctg, centre, sp, u in REAL:
            c = read_contig(path, ctg)
            w = c[max(0, centre - 20000) : centre + 20000]
            fa.write(f">{aid}\n{w}\n")
            rows.append((aid, f"long-read {ctg}:{centre}", sp, u, len(w)))
            print(f"{aid}: window {len(w)} nt from {path}", file=sys.stderr)
    with open(a.out + ".tsv", "w") as f:
        f.write("allele\tsource\tspecies\ttrue_units\twindow_len\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


def sim(a):
    seq = None
    with open(a.fasta) as fh:
        name = None
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            elif name == a.allele:
                seq = line.strip()
                break
    if seq is None:
        sys.exit(f"allele {a.allele} not in {a.fasta}")
    rng = random.Random(a.seed)
    n_pairs = int(round(a.cov * len(seq) / (2 * a.len)))
    q = chr(33 + 30) * a.len
    bases = "ACGT"

    def mutate(r):
        r = list(r)
        for i in range(len(r)):
            if rng.random() < a.err:
                r[i] = rng.choice([b for b in bases if b != r[i]])
        return "".join(r)

    with gzip.open(a.out + "_R1.fq.gz", "wt") as f1, gzip.open(a.out + "_R2.fq.gz", "wt") as f2:
        for i in range(n_pairs):
            while True:
                fl = int(round(rng.gauss(a.ins, a.ins_sd)))
                if a.len <= fl < len(seq):
                    break
            s = rng.randrange(0, len(seq) - fl + 1)
            frag = seq[s : s + fl]
            if rng.random() < 0.5:
                frag = frag.translate(COMP)[::-1]
            r1 = mutate(frag[: a.len])
            r2 = mutate(frag[-a.len :].translate(COMP)[::-1])
            f1.write(f"@p{i}/1\n{r1}\n+\n{q}\n")
            f2.write(f"@p{i}/2\n{r2}\n+\n{q}\n")
    print(f"{a.allele} L={a.len} cov={a.cov} pairs={n_pairs}", file=sys.stderr)


ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
sub = ap.add_subparsers(dest="cmd", required=True)
b = sub.add_parser("build")
b.add_argument("--out", default="sim_alleles")
s = sub.add_parser("sim")
s.add_argument("--fasta", required=True)
s.add_argument("--allele", required=True)
s.add_argument("--len", type=int, required=True)
s.add_argument("--cov", type=float, required=True)
s.add_argument("--ins", type=float, required=True)
s.add_argument("--ins-sd", type=float, default=120)
s.add_argument("--err", type=float, default=0.002)
s.add_argument("--seed", type=int, required=True)
s.add_argument("--out", required=True)
a = ap.parse_args()
build(a) if a.cmd == "build" else sim(a)
