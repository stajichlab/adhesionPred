#!/usr/bin/env python3
"""Select the follow-up gene set and the proteome list for the ortholog search.

Gene set S = (A) specific extreme genes without a signal peptide
           + (B) Cys-rich spherule-up genes without a signal peptide that pass the specificity flag.
Definitions are the ones in docs/reports/2026-10-03-cocci-spherule-surface-table.md.
Search set = all Onygenales proteomes in Fungi5k, plus a fixed list of outgroup genera (one
proteome per genus, the first by name) so that "outside Onygenales" has a fixed meaning.
Outputs also: up_proteins.fa (all spherule-up genes).
Usage: 10_select_followup.py --table <table.tsv.gz> --seqs <unique_sequences.fasta.gz> --outdir <dir>
Outputs: followup_genes.tsv, followup_proteins.fa, search_proteomes.tsv
"""

import argparse
import csv
import gzip
import os
import re
from pathlib import Path

import duckdb

DB = "/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb"
INPUT = "/bigdata/stajichlab/shared/projects/Fungi_5k/input"
# One proteome per genus outside Onygenales (first file by name). Fixed list, chosen to span
# Eurotiomycetes and other Pezizomycotina, and two yeasts, so a hit can be placed in the tree.
OUTGROUP_GENERA = [
    "Aspergillus", "Penicillium", "Talaromyces", "Monascus", "Exophiala", "Cladophialophora",
    "Fonsecaea", "Neurospora", "Fusarium", "Trichoderma", "Magnaporthe", "Botrytis",
    "Sclerotinia", "Saccharomyces", "Candida", "Cryptococcus", "Ustilago",
]  # fmt: skip


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--seqs", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(gzip.open(a.table, "rt"), delimiter="\t"))
    nz = lambda x: float(x) if x else None  # noqa: E731
    sel = []
    for r in rows:
        up = r["flag_spherule_up"] == "yes"
        ext = r["flag_spherule_extreme"] == "yes"
        nosp = r["signalp_call"] != "SP"
        spec = r["flag_cocci_specific"] == "yes"
        cys_rich = nz(r["cys_frac"]) >= 0.0392  # genome-wide 95th percentile, see report section 6
        tags = []
        if ext and spec and nosp:
            tags.append("A_extreme_specific_noSP")
        if up and cys_rich and nosp and spec:
            tags.append("B_cysrich_up_specific_noSP")
        if tags:
            sel.append((r, ";".join(tags)))
    cols = ["gene_id", "protein_id", "set", "product", "length", "log2fc_48h", "padj_48h",
            "tpm_spherule48h", "signalp_sp_prob", "cys_frac", "pro_frac", "prevalence"]  # fmt: skip
    with open(out / "followup_genes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r, t in sel:
            w.writerow([r.get(c, "") if c != "set" else t for c in cols])
    seqs, name, ch = {}, None, []
    with gzip.open(a.seqs, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(ch)
                name, ch = line[1:].split()[0], []
            else:
                ch.append(line.strip())
        if name:
            seqs[name] = "".join(ch)
    with open(out / "up_proteins.fa", "w") as fh:  # all spherule-up genes, for Pfam and Phobius
        for r in rows:
            if r["flag_spherule_up"] == "yes":
                fh.write(f">{r['gene_id']}|{r['protein_id']}\n{seqs[r['protein_id']]}\n")
    with open(out / "followup_proteins.fa", "w") as fh:
        for r, _ in sel:
            fh.write(f">{r['gene_id']}|{r['protein_id']}\n{seqs[r['protein_id']]}\n")

    con = duckdb.connect(DB, read_only=True)
    sp = con.execute('select SPECIESIN,STRAIN,"ORDER",GENUS,LOCUSTAG from species').fetchall()
    files = sorted(f for f in os.listdir(INPUT) if f.endswith(".proteins.fa"))
    by_key = {}
    for f in files:
        by_key.setdefault(norm(f[: -len(".proteins.fa")]), f)
    chosen, unmapped, seen_genus = [], [], set()
    for spin, strain, order, genus, lt in sp:
        is_ony = order == "Onygenales"
        if not is_ony and (genus not in OUTGROUP_GENERA or genus in seen_genus):
            continue
        key = norm(f"{spin}_{strain}")
        f = by_key.get(key)
        if f is None:  # fall back to any file of the same species name
            pre = norm(str(spin))
            cands = [v for k, v in by_key.items() if k.startswith(pre)]
            f = sorted(cands)[0] if cands else None
        if f is None:
            unmapped.append((spin, strain, order))
            continue
        if not is_ony:
            seen_genus.add(genus)
        chosen.append((f"{INPUT}/{f}", spin, strain, order or "", genus or "", lt))
    # C. posadasii references from the pangenome inputs (not in the Fungi5k name map)
    pg = "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/input"
    for f, sp_name, strain in (
        (
            "Coccidioides_posadasii_UCSF_Cp_Silveira.proteins.fa",
            "Coccidioides posadasii",
            "Silveira",
        ),
        ("Coccidioides_posadasii_C735-TKO.proteins.fa", "Coccidioides posadasii", "C735-TKO"),
    ):
        chosen.append((f"{pg}/{f}", sp_name, strain, "Onygenales", "Coccidioides", "pangenome"))
    with open(out / "search_proteomes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "species", "strain", "order", "genus", "locustag"])
        w.writerows(chosen)
    n_ony = sum(1 for c in chosen if c[3] == "Onygenales")
    print(
        f"genes: {len(sel)} ({sum(1 for _, t in sel if 'A_' in t)} in A, "
        f"{sum(1 for _, t in sel if 'B_' in t)} in B)"
    )
    print(
        f"proteomes: {len(chosen)} ({n_ony} Onygenales); unmapped: {len(unmapped)} {unmapped[:5]}"
    )
    _ = os.path.getsize  # keep os import used


if __name__ == "__main__":
    main()
