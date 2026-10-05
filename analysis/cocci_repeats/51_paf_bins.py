#!/usr/bin/python3.12
"""Reduce one minimap2 PAF (assembly -> reference) to per-bin reference coverage.

Reads PAF on stdin. For every reference contig, builds a depth track from all alignment blocks
of at least --min-aln bp, then reports, per fixed-size bin of the reference:

    covered   fraction of the bin's bases with depth >= 1
    depth     mean number of assembly alignments over the bin

`covered` answers "is this stretch of the reference present in the assembly at all".
`depth` above 1 flags a stretch several reference copies map onto, which is how a collapsed
repeat looks: one assembled contig standing in for many reference copies.

Caveat that matters when reading the output: a bin that is uncovered is not proof that the
strain lacks the sequence. It can be a real deletion in that strain, or divergence too high to
align, or sequence the assembler dropped. Telling those apart needs the reads. This script only
measures what the assembly contains relative to one reference.

Usage:  minimap2 ... | 51_paf_bins.py STRAIN SPECIES REFNAME BINSIZE REFLENS > strain.bins.tsv
REFLENS is a two-column TSV (contig, length) for EVERY reference contig. It is required so that a
reference contig with no alignments at all is reported as uncovered, not silently dropped.
Output columns: strain species ref contig bin_start bin_end covered depth
"""

import sys

import numpy as np

strain, species, ref, binsize, reflens = (
    sys.argv[1],
    sys.argv[2],
    sys.argv[3],
    int(sys.argv[4]),
    sys.argv[5],
)
MIN_ALN = 500

tlen = {}
diff = {}
for line in open(reflens):
    c, n = line.split()
    tlen[c] = int(n)
    diff[c] = np.zeros(int(n) + 1, dtype=np.int32)
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    t, tl, ts, te, alnlen = f[5], int(f[6]), int(f[7]), int(f[8]), int(f[10])
    if alnlen < MIN_ALN:
        continue
    diff[t][ts] += 1
    diff[t][te] -= 1

for t in sorted(tlen):
    depth = np.cumsum(diff[t])[: tlen[t]]
    for b in range(0, tlen[t], binsize):
        seg = depth[b : b + binsize]
        sys.stdout.write(
            f"{strain}\t{species}\t{ref}\t{t}\t{b}\t{b + len(seg)}\t"
            f"{(seg >= 1).mean():.4f}\t{seg.mean():.3f}\n"
        )
