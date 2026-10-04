#!/usr/bin/python3.12
"""Parse 62_asm_locus_check.sh output: is the RS SOWgp locus in each strain's ASSEMBLY?

Query (target in the PAF) = RS GG704914:968,094-973,690. Query coordinates, 0-based half-open:
  gene  CIMG_04613   RS 970,094-971,690  -> [2000, 3597)
  array (aa 84-270)  RS 970,511-971,295  -> [2417, 3202)   785 nt = 4 genomic periods

Per strain and assembly type (asm = AAFTF contigs; ann = the scaffolds the gene models were
predicted on) the script reports:
  gene_cov, array_cov   fraction of RS positions covered by any minimap2 alignment
  array_spanned         one alignment covers array +/- 50 nt
  array_asm_len         assembly bases aligned inside the array (from the cs tag), if spanned
  array_n               'n' bases of the assembly inside the array (scaffold gaps), if spanned
  units_asm             4 + (array_asm_len - 785) / 197, rounded; only if spanned and array_n == 0
  blast_gene_cov        same as gene_cov but from blastn HSPs (pident >= 85, length >= 100)
  class                 contiguous   gene and array on one alignment, no N in the array
                        gap_in_array array covered by no single alignment, or N bases in it
                        partial      gene covered < 90% but > 10%
                        absent       gene covered <= 10% by both minimap2 and blastn
Joins pangenome status (sowgp_tree_units.tsv), the anchored search
(sowgp_anchored_pangenome.tsv) and depth (CIMG_04613.coverage.tsv).

Inputs : asm_locus_raw.tar.gz (from 62), sowgp_tree_units.tsv, sowgp_anchored_pangenome.tsv,
         CIMG_04613.coverage.tsv
Outputs: asm_locus_summary.tsv (one row per strain), stdout summary
Run    : /usr/bin/python3.12 63_asm_locus_parse.py
"""

import re
import tarfile
from collections import defaultdict

import numpy as np
import pandas as pd

GENE = (2000, 3597)
ARRAY = (2417, 3202)
PERIOD = 197
PAD = 50


def parse_cs(cs, ts):
    """Yield (op, tpos_start, tlen, qlen, qbases) along the target."""
    t = ts
    for m in re.finditer(r"(:\d+)|(\*[a-z][a-z])|(\+[a-z]+)|(-[a-z]+)", cs):
        tok = m.group(0)
        if tok[0] == ":":
            n = int(tok[1:])
            yield ("m", t, n, n, "")
            t += n
        elif tok[0] == "*":
            yield ("x", t, 1, 1, tok[2])
            t += 1
        elif tok[0] == "+":
            yield ("i", t, 0, len(tok) - 1, tok[1:])
        else:
            yield ("d", t, len(tok) - 1, 0, "")
            t += len(tok) - 1


def region_stats(cs, ts, a, b):
    """Assembly bases aligned inside target [a, b), and 'n' bases among them."""
    qn = nn = 0
    for op, t, tl, ql, qb in parse_cs(cs, ts):
        if op == "m":
            qn += max(0, min(t + tl, b) - max(t, a))
        elif op == "x":
            if a <= t < b:
                qn += 1
                nn += qb == "n"
        elif op == "i":
            if a < t < b:
                qn += ql
                nn += qb.count("n")
    return qn, nn


def cover(intervals, a, b):
    pos = set()
    for s, e in intervals:
        pos.update(range(max(s, a), min(e, b)))
    return len(pos) / (b - a)


raw = defaultdict(dict)
mito_len = {}
with tarfile.open("asm_locus_raw.tar.gz") as tf:
    for m in tf.getmembers():
        if not m.isfile():
            continue
        name = m.name.split("/")[-1]
        mm = re.match(r"(.+)\.(asm|ann|mito)\.(paf|blastn\.tsv)$", name)
        if mm:
            raw[(mm.group(1), mm.group(2))][mm.group(3)] = tf.extractfile(m).read().decode()
        ml = re.match(r"(.+)\.mito\.len$", name)
        if ml:
            mito_len[ml.group(1)] = int(tf.extractfile(m).read().decode().strip())

rows = []
for (strain, kind), d in sorted(raw.items()):
    alns = []
    for line in d.get("paf", "").splitlines():
        f = line.split("\t")
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), "")
        alns.append(
            {"contig": f[0], "ts": int(f[7]), "te": int(f[8]), "cs": cs, "mapq": int(f[11])}
        )
    hsps = []
    for line in d.get("blastn.tsv", "").splitlines():
        f = line.split("\t")
        if float(f[2]) >= 85 and int(f[3]) >= 100:
            hsps.append((int(f[6]) - 1, int(f[7])))
    gene_cov = cover([(x["ts"], x["te"]) for x in alns], *GENE)
    array_cov = cover([(x["ts"], x["te"]) for x in alns], *ARRAY)
    blast_cov = cover(hsps, *GENE)
    span = [x for x in alns if x["ts"] <= ARRAY[0] - PAD and x["te"] >= ARRAY[1] + PAD]
    gene_one = any(x["ts"] <= GENE[0] and x["te"] >= GENE[1] for x in alns)
    alen = an = units = None
    if span:
        best = max(span, key=lambda x: x["te"] - x["ts"])
        alen, an = region_stats(best["cs"], best["ts"], *ARRAY)
        if an == 0:
            units = 4 + round((alen - (ARRAY[1] - ARRAY[0])) / PERIOD)
    gene_contigs = {x["contig"] for x in alns if x["te"] > GENE[0] and x["ts"] < GENE[1]}
    if gene_cov <= 0.10 and blast_cov <= 0.10:
        cls = "absent"
    elif span and an == 0 and gene_one:
        cls = "contiguous"
    elif not span or (an or 0) > 0:
        cls = "gap_in_array"
    else:
        cls = "partial"
    rows.append(
        {
            "strain": strain,
            "kind": kind,
            "n_aln": len(alns),
            "gene_contigs": len(gene_contigs),
            "gene_cov": round(gene_cov, 3),
            "array_cov": round(array_cov, 3),
            "blast_gene_cov": round(blast_cov, 3),
            "gene_one_aln": gene_one,
            "array_spanned": bool(span),
            "array_asm_len": alen,
            "array_n": an,
            "units_asm": units,
            "class": cls,
        }
    )

R = pd.DataFrame(rows)
W = R.pivot(index="strain", columns="kind")
W.columns = [f"{k}_{c}" for c, k in W.columns]
W = W.reset_index()

tu = pd.read_csv("sowgp_tree_units.tsv", sep="\t")
tu["strain"] = tu.tip.str.replace(r"^Coccidioides_(immitis|posadasii)_", "", regex=True)
ap = pd.read_csv("sowgp_anchored_pangenome.tsv", sep="\t")
ap = ap[ap.member.astype(str) == "True"]
ap["strain"] = ap.proteome.str.replace(r"^Coccidioides_(immitis|posadasii)_", "", regex=True)
anch = set(ap.strain)
cov = pd.read_csv("CIMG_04613.coverage.tsv", sep="\t")
W = W.merge(tu[["strain", "clade_species", "status", "n_units"]], on="strain", how="left")
W["status"] = W.status.fillna("not in tree")
W["anchored_hit"] = W.strain.isin(anch)
W = W.merge(cov[["strain", "ratio", "genome_mean"]], on="strain", how="left")
W["mito_bin_len"] = W.strain.map(mito_len)
W["locus_in_mito_bin"] = W.get("mito_gene_cov", pd.Series(0, index=W.index)).fillna(0) > 0.5
W.to_csv("asm_locus_summary.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
print(
    f"strains: {len(W)}  (asm rows {W.asm_class.notna().sum()}, ann rows {W.ann_class.notna().sum()}, "
    f"mito rows {W.mito_class.notna().sum()})"
)
for k in ("asm", "ann"):
    print(f"\n== {k}: class by pangenome status")
    print(pd.crosstab(W.status, W[f"{k}_class"], margins=True).to_string())
print("\n== locus (gene >50% covered) in the AAFTF mitochondrial bin, by status")
print(pd.crosstab(W.status, W.locus_in_mito_bin, margins=True).to_string())
print(
    "\nmitochondrial bin size (bp) over all strains: "
    f"median {W.mito_bin_len.median():.0f}, IQR {W.mito_bin_len.quantile(0.25):.0f}-{W.mito_bin_len.quantile(0.75):.0f}, "
    f"max {W.mito_bin_len.max():.0f}; strains with bin > 200 kb: {(W.mito_bin_len > 200000).sum()}"
)
nm = W[W.status == "no SOWgp gene model"].copy()
print(f"\n== the {len(nm)} no-model strains: anchored-search hit in {nm.anchored_hit.sum()}")
print(
    nm[
        [
            "strain",
            "clade_species",
            "anchored_hit",
            "ratio",
            "asm_class",
            "asm_gene_contigs",
            "asm_array_n",
            "ann_class",
            "ann_array_n",
            "locus_in_mito_bin",
            "mito_bin_len",
        ]
    ].to_string(index=False)
)
nm["why"] = np.select(
    [
        nm.locus_in_mito_bin & (nm.asm_class == "absent"),
        nm.asm_array_spanned.astype(bool) & (nm.asm_array_n.fillna(0) > 0),
        ~nm.asm_array_spanned.astype(bool) & (nm.asm_gene_cov > 0),
    ],
    [
        "locus contig moved to the mitochondrial bin",
        "N gap inside the array (one AAFTF scaffold)",
        "AAFTF contigs break in the array (N gap after scaffolding)",
    ],
    "other",
)
print("\nwhy the 39 have no gene model (assembly view):")
print(nm.why.value_counts().to_string())
print(
    "\nN bases in the array, annotated scaffolds, no-model strains: "
    f"{nm.ann_array_n.describe()[['count', 'min', '50%', 'max']].to_dict()}"
)
for s in ("full-length", "fragment only", "no SOWgp gene model"):
    x = W[W.status == s]
    print(
        f"{s:22s} n={len(x):3d} contiguous gap-free array in ann scaffolds: {(x.ann_class == 'contiguous').sum()} "
        f"({(x.ann_class == 'contiguous').mean():.0%})"
    )
fl = W[(W.status == "full-length") & W.asm_units_asm.notna()]
print("\n== check of units_asm against gene-model unit count, full-length strains (asm)")
print(pd.crosstab(fl.n_units, fl.asm_units_asm).to_string())
fl = W[(W.status == "full-length") & W.ann_units_asm.notna()]
print("\n== same, ann")
print(pd.crosstab(fl.n_units, fl.ann_units_asm).to_string())

# read length (64) against array assembly class: do short-read libraries break the array more often?
rl = pd.read_csv("cram_readlen.tsv", sep="\t")[["strain", "modal_len"]]
X = W.merge(rl, on="strain", how="inner")
X["rl_class"] = pd.cut(
    X.modal_len, [0, 101, 126, 151, 251, 301], labels=["<=101", "126", "150", "251", "301"]
)
print("\n== annotated-scaffold array class by read-length class (tree strains)")
print(
    pd.crosstab(
        X[X.status != "not in tree"].rl_class, X[X.status != "not in tree"].ann_class, margins=True
    ).to_string()
)
print("\n== pangenome status by read-length class (tree strains)")
print(
    pd.crosstab(
        X[X.status != "not in tree"].rl_class, X[X.status != "not in tree"].status, margins=True
    ).to_string()
)
