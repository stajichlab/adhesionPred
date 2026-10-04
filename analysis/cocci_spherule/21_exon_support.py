#!/usr/bin/env python3
"""Compare RefSeq gene models with RNA-seq alignments (STAR, no annotation given) for chosen genes.

For each gene (longest protein's mRNA from the RefSeq GFF3) it reports, from the pooled spherule
samples (48 h and 8 d, two replicates each):
  - mean depth and fraction of exon bases with depth >= MIN_DEPTH, in exons, introns and flanks
  - annotated introns with >= MIN_JUNC unique junction reads (STAR SJ.out.tab, summed)
  - novel junctions inside the gene span with >= NOVEL_JUNC unique reads
  - a call (display thresholds, not validated)
Flank bases that overlap an exon of any other gene are left out. Uses `samtools depth`.
Usage: 21_exon_support.py --gff <genomic.gff.gz> --star <dir> --table <table.tsv.gz>
                          --genes <genes.tsv> --out <tsv> [--controls up]
Only unique alignments (MAPQ 255) are counted. Standard library only, plus the samtools module.
"""

import argparse
import bisect
import csv
import gzip
import statistics as st
import subprocess
from collections import defaultdict
from pathlib import Path

MIN_DEPTH = 5
MIN_JUNC = 3
NOVEL_JUNC = 5
FLANK = 300
SAMPLES = ("Spherule48H.r1", "Spherule48H.r2", "Spherule8D.r1", "Spherule8D.r2")


def read_gff(path):
    """Return {xp: (seqid, strand, [(start, end) exons], gene)} and exon intervals per seqid."""
    mrna_gene, xp_mrna, exons = {}, {}, defaultdict(list)
    strand, seqid = {}, {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "mRNA":
                mrna_gene[a["ID"]] = a.get("Parent", "").replace("gene-", "")
                strand[a["ID"]], seqid[a["ID"]] = f[6], f[0]
            elif f[2] == "exon":
                exons[a["Parent"]].append((int(f[3]), int(f[4])))
            elif f[2] == "CDS" and a.get("Name", "").startswith("XP_"):
                xp_mrna[a["Name"]] = a["Parent"]
    models = {}
    for xp, m in xp_mrna.items():
        if m in exons:
            models[xp] = (seqid[m], strand[m], sorted(exons[m]), mrna_gene[m])
    all_exons = defaultdict(list)
    for m, ex in exons.items():
        for s, e in ex:
            all_exons[seqid[m]].append((s, e, mrna_gene[m]))
    for k in all_exons:
        all_exons[k].sort()
    return models, all_exons


def read_junctions(star_dir):
    """Sum unique-read counts per (seqid, start, end) over the spherule samples."""
    junc = defaultdict(int)
    for s in SAMPLES:
        with open(Path(star_dir) / s / "SJ.out.tab") as fh:
            for line in fh:
                f = line.split("\t")
                junc[(f[0], int(f[1]), int(f[2]))] += int(f[6])
    by_chr = defaultdict(list)
    for (c, s, e), n in junc.items():
        by_chr[c].append((s, e, n))
    for c in by_chr:
        by_chr[c].sort()
    return junc, by_chr


def depth(star_dir, seqid, start, end):
    bams = [str(Path(star_dir) / s / "Aligned.sortedByCoord.out.bam") for s in SAMPLES]
    out = subprocess.run(
        ["samtools", "depth", "-a", "-Q", "255", "-r", f"{seqid}:{max(start, 1)}-{end}", *bams],
        check=True, capture_output=True, text=True,
    ).stdout  # fmt: skip
    d = {}
    for line in out.splitlines():
        f = line.split("\t")
        d[int(f[1])] = sum(int(x) for x in f[2:])
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--star", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--genes", required=True, help="TSV with a gene_id column")
    ap.add_argument("--out", required=True)
    ap.add_argument("--controls", choices=["none", "up"], default="up")
    a = ap.parse_args()

    models, all_exons = read_gff(a.gff)
    table = {r["gene_id"]: r for r in csv.DictReader(gzip.open(a.table, "rt"), delimiter="\t")}
    targets = [r["gene_id"] for r in csv.DictReader(open(a.genes), delimiter="\t")]
    groups = dict.fromkeys(targets, "followup")
    if a.controls == "up":
        for g, r in table.items():
            if r["flag_spherule_up"] == "yes" and g not in groups:
                groups[g] = "other_spherule_up"
    for g, r in table.items():
        if r["anchor"] and g not in groups:
            groups[g] = "anchor"
    junc, by_chr = read_junctions(a.star)

    cols = ["gene_id", "group", "protein_id", "seqid", "strand", "gene_start", "gene_end", "n_exons",
            "exon_bases", "mean_depth_exon", "frac_exon_ge_min", "mean_depth_intron", "mean_depth_flank",
            "n_introns", "introns_supported", "novel_junctions", "min_intron_junction_reads", "call"]  # fmt: skip
    rows = []
    for g, grp in groups.items():
        xp = table[g]["protein_id"]
        if xp not in models:
            continue
        seqid, strand, exons, _ = models[xp]
        gs, ge = exons[0][0], exons[-1][1]
        dep = depth(a.star, seqid, gs - FLANK, ge + FLANK)
        exon_pos = {p for s, e in exons for p in range(s, e + 1)}
        other = [
            (s, e)
            for s, e, gn in all_exons[seqid]
            if gn != g and e >= gs - FLANK and s <= ge + FLANK
        ]
        other_pos = {p for s, e in other for p in range(s, e + 1)}
        intron_pos = [p for p in range(gs, ge + 1) if p not in exon_pos and p not in other_pos]
        flank_pos = [
            p
            for p in list(range(gs - FLANK, gs)) + list(range(ge + 1, ge + FLANK + 1))
            if p > 0 and p not in other_pos
        ]
        ed = [dep.get(p, 0) for p in sorted(exon_pos)]
        introns = [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]
        jr = [junc.get((seqid, s, e), 0) for s, e in introns]
        ann = set(introns)
        starts = [x[0] for x in by_chr[seqid]]
        lo = bisect.bisect_left(starts, gs)
        novel = 0
        for s, e, n in by_chr[seqid][lo:]:
            if s > ge:
                break
            if e <= ge and (s, e) not in ann and n >= NOVEL_JUNC:
                novel += 1
        frac = sum(1 for x in ed if x >= MIN_DEPTH) / len(ed)
        mean_ed = st.mean(ed)
        mean_fl = st.mean([dep.get(p, 0) for p in flank_pos]) if flank_pos else 0.0
        if mean_ed < 1:
            call = "no_coverage"
        elif frac < 0.9:
            call = "partial_exon_coverage"
        elif introns and min(jr) < MIN_JUNC:
            call = "intron_unsupported"
        elif novel:
            call = "novel_junction_in_gene"
        elif mean_fl > 0.5 * mean_ed:
            call = "flank_expression"
        else:
            call = "supported"
        rows.append({
            "gene_id": g, "group": grp, "protein_id": xp, "seqid": seqid, "strand": strand,
            "gene_start": gs, "gene_end": ge, "n_exons": len(exons), "exon_bases": len(ed),
            "mean_depth_exon": f"{mean_ed:.1f}", "frac_exon_ge_min": f"{frac:.3f}",
            "mean_depth_intron": f"{st.mean([dep.get(p, 0) for p in intron_pos]):.1f}" if intron_pos else "",
            "mean_depth_flank": f"{mean_fl:.1f}", "n_introns": len(introns),
            "introns_supported": sum(1 for n in jr if n >= MIN_JUNC), "novel_junctions": novel,
            "min_intron_junction_reads": min(jr) if jr else "", "call": call,
        })  # fmt: skip
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    by = defaultdict(lambda: defaultdict(int))
    for r in rows:
        by[r["group"]][r["call"]] += 1
    for grp, c in by.items():
        print(grp, sum(c.values()), dict(c))


if __name__ == "__main__":
    main()
