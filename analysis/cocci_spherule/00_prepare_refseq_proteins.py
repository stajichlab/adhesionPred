#!/usr/bin/env python3
"""Write the C. immitis RS RefSeq proteins (GCF_000149335.2) with a protein-to-gene map.

RefSeq FASTA headers carry XP_ accessions. The CIMG_ locus tag comes from the GFF3 (CDS
Name -> Parent mRNA -> locus_tag). The FASTA ID is the XP_ accession, so the step 1
feature job (analysis/step1_compare/jobs/j1_features.sh) can run on it unchanged.
Usage: 00_prepare_refseq_proteins.py <ref_dir> <step1_workdir>
Outputs: <step1_workdir>/phaseb/unique_sequences.fasta.gz, <ref_dir>/protein_map.tsv
"""

import csv
import gzip
import re
import sys
from pathlib import Path

ref = Path(sys.argv[1])
work = Path(sys.argv[2])
pre = "GCF_000149335.2_ASM14933v2_"

xp_gene = {}
with gzip.open(ref / f"{pre}genomic.gff.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        lt = attrs.get("locus_tag")
        xp = attrs.get("Name") or attrs.get("protein_id")
        if lt and xp and xp.startswith("XP_"):
            xp_gene[xp] = lt

records = []
name, product, seq = None, "", []
with gzip.open(ref / f"{pre}protein.faa.gz", "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            if name:
                records.append((name, product, "".join(seq)))
            name = line[1:].split()[0]
            m = re.search(r"\s(.*?)\s\[", line)
            product = m.group(1) if m else ""
            seq = []
        else:
            seq.append(line.strip())
    if name:
        records.append((name, product, "".join(seq)))

missing = [n for n, _, _ in records if n not in xp_gene]
if missing:
    sys.exit(f"STOP: {len(missing)} proteins have no locus tag in the GFF3, first: {missing[:3]}")

out = work / "phaseb"
out.mkdir(parents=True, exist_ok=True)
with gzip.open(out / "unique_sequences.fasta.gz", "wt") as fh:
    for n, _, s in records:
        fh.write(f">{n}\n{s}\n")
with open(ref / "protein_map.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["protein_id", "gene_id", "product", "length"])
    for n, p, s in records:
        w.writerow([n, xp_gene[n], p, len(s)])
genes = {xp_gene[n] for n, _, _ in records}
print(f"{len(records)} proteins, {len(genes)} genes")
