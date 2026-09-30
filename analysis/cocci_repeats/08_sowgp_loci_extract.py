#!/usr/bin/python3.12
"""Extract SOWgp-like genomic loci in long-read assemblies where the gene call was missing.

Reads the unsegmented tblastn output (sowgp_tblastn_noseg_<strain>.out), clusters HSP
subject coordinate spans into loci per contig, extracts a window around each locus, and
translates all 6 reading frames. Writes ORF candidates to sowgp_loci_orfs.fa.

The point is to distinguish:
  - locus absent from assembly       -> no locus
  - locus present but gene model missed it / locus fragmented -> ORFs of partial length
Independent of whether the repeat is "called" by the predictor.
"""

import os
import re
import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

HERE = os.path.dirname(os.path.abspath(__file__))
LONGREAD = "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc"
STRAINS = ["VFC140", "Cpos1038", "Cpos3700"]
WINDOW = 2000
MIN_ORF = 50
COLS = [
    "qseqid",
    "qlen",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def loci_from_hits(df):
    """Cluster subject HSPs into loci (contig + overlapping interval).

    Returns list of (sseqid, a, b) genomic intervals.
    """
    out = []
    for sseqid, g in df.groupby("sseqid"):
        ivs = sorted((min(a, b), max(a, b)) for a, b in g[["sstart", "send"]].values)
        cur = list(ivs[0])
        for a, b in ivs[1:]:
            if a <= cur[1] + 500:
                cur[1] = max(cur[1], b)
            else:
                out.append((sseqid, cur[0], cur[1]))
                cur = [a, b]
        out.append((sseqid, cur[0], cur[1]))
    return out


def six_frame_orfs(seq):
    """Return ORFs >= MIN_ORF aa in all 6 frames, with frame label."""
    out = []
    for strand, s in ((1, seq), (-1, seq.reverse_complement())):
        for frame in range(3):
            aa = str(Seq(s[frame:]).translate(to_stop=False))
            # split at stops
            for seg in re.split(r"\*", aa):
                for mm in re.finditer(r"M[^*]{0,}", seg):
                    if mm.group() and len(mm.group()) >= MIN_ORF:
                        out.append((strand, frame, mm.start(), mm.group()))
    return out


def main():
    orf_recs = []
    for strain in STRAINS:
        out = os.path.join(HERE, f"sowgp_tblastn_noseg_{strain}.out")
        df = pd.read_csv(out, sep="\t", header=None, names=COLS)
        fasta = os.path.join(LONGREAD, strain, f"Coccidioides_*_{strain}.scaffolds.fa")
        import glob

        path = glob.glob(fasta)[0]
        contigs = {r.id: r.seq for r in SeqIO.parse(path, "fasta")}
        for sseqid, a, b in loci_from_hits(df):
            lo, hi = max(0, a - WINDOW), min(len(contigs[sseqid]), b + WINDOW)
            window = contigs[sseqid][lo:hi]
            for strand, frame, start, orf in six_frame_orfs(window):
                info = f"{strain}|{sseqid}|{lo}|{hi}|f{'+' if strand > 0 else '-'}{frame}|{start}"
                orf_recs.append(
                    SeqIO.SeqRecord(Seq(orf), id=info, description=f"locus len={len(orf)}")
                )
        print(f"{strain}: {len(orf_recs)} ORFs so far", file=sys.stderr)
    with open(os.path.join(HERE, "sowgp_loci_orfs.fa"), "w") as fh:
        SeqIO.write(orf_recs, fh, "fasta")
    print(f"wrote sowgp_loci_orfs.fa with {len(orf_recs)} ORF candidates")


if __name__ == "__main__":
    main()
