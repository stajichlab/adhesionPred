#!/usr/bin/env python3
"""One row per C. immitis RS gene: spherule expression, surface features, specificity.

Joins, all keyed on the CIMG_ gene ID of the RefSeq / FungiDB-46 annotation (see README):
  - DESeq2 (01_deseq2_rs1.R) and RS1 kallisto TPM       -> expression
  - SignalP 6, PredGPI, TMHMM, Ser/Thr/Cys/Pro fraction  -> surface features
  - cocci_antigens/cocci_antigen_ranking.tsv             -> specificity (pangenome, 488 proteomes)
One protein per gene: the longest one. `n_proteins` says how many isoforms the gene has.

Usage: 02_build_table.py --work <_workdir/cocci_spherule> --out <table.tsv>
Needs only the standard library.
"""

import argparse
import csv
import gzip
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "analysis" / "step1_compare"))
import feature_parsers as fp  # noqa: E402

TPM = "/bigdata/stajichlab/shared/projects/Coccidioides/csSeq/UCSD_SpheruleMycelium/reports/RS1_kallisto.TPM.csv"
RANKING = REPO / "analysis" / "cocci_antigens" / "cocci_antigen_ranking.tsv"

# Thresholds for the flags. They are display filters, not validated cut-offs.
UP_LOG2FC = 2.0
UP_PADJ = 0.05
EXTREME_LOG2FC = 6.0


def read_fasta_gz(path):
    seqs, name, chunks = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fmt(x, nd=3):
    return "" if x is None else (f"{x:.{nd}g}" if isinstance(x, float) else str(x))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    work = Path(args.work)
    s1 = work / "s1" / "phaseb"

    pmap = list(csv.DictReader(open(work / "ref" / "protein_map.tsv"), delimiter="\t"))
    seqs = read_fasta_gz(s1 / "unique_sequences.fasta.gz")
    by_gene = defaultdict(list)
    for r in pmap:
        by_gene[r["gene_id"]].append(r)

    sp = {}
    for part in sorted((s1 / "signalp").glob("part_*/prediction_results.txt.gz")):
        sp.update(fp.parse_signalp(part))
    gpi = {}
    for part in sorted((s1 / "predgpi").glob("part_*.tsv.gz")):
        gpi.update(fp.parse_predgpi_scores(part))
    for name, d in (("SignalP", sp), ("PredGPI", gpi)):
        missing = [p["protein_id"] for p in pmap if p["protein_id"] not in d]
        if missing:
            sys.exit(f"STOP: {len(missing)} proteins have no {name} call, first {missing[:3]}")
    tm = {}
    with open(work / "tmhmm.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tm[r["protein_id"]] = r

    de = {
        r["gene"]: r
        for r in csv.DictReader(open(work / "deseq2" / "deseq2_rs1_gene.tsv"), delimiter="\t")
    }

    tpm = defaultdict(lambda: [0.0] * 6)
    with open(TPM) as fh:
        rd = csv.reader(fh)
        next(rd)
        for row in rd:
            g = row[0].strip('"').split("-t26")[0]
            for i, v in enumerate(row[1:7]):
                tpm[g][i] += float(v)

    rank = {}
    with open(RANKING) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["protein"].split("-t26")[0]
            if g not in rank or int(r["rank"]) < int(rank[g]["rank"]):
                rank[g] = r

    cols = [
        "gene_id", "product", "n_proteins", "protein_id", "length",
        "tpm_mycelia", "tpm_spherule48h", "tpm_spherule8d",
        "log2fc_48h", "padj_48h", "log2fc_8d", "padj_8d",
        "signalp_call", "signalp_sp_prob", "predgpi_call", "predgpi_fpr", "tm_helices",
        "ser_thr_frac", "cys_frac", "pro_frac",
        "orthogroup", "prevalence", "copy_number_mean", "max_fungal_crossreact_pid",
        "n_confounder_genera", "human_homolog_pid", "anchor",
        "flag_spherule_up", "flag_spherule_extreme", "flag_signal_peptide", "flag_gpi_call",
        "flag_tm", "flag_cocci_specific",
    ]  # fmt: skip
    rows = []
    for g in sorted(by_gene):
        prots = sorted(by_gene[g], key=lambda r: -int(r["length"]))
        p = prots[0]
        pid = p["protein_id"]
        s = seqs[pid]
        L = len(s)
        d = de.get(g, {})
        t = tpm.get(g)
        l48, p48 = num(d.get("log2FoldChange.Spherule48H")), num(d.get("padj.Spherule48H"))
        l8, p8 = num(d.get("log2FoldChange.Spherule8D")), num(d.get("padj.Spherule8D"))
        rk = rank.get(g, {})
        sc, gc = sp[pid], gpi[pid]
        tmh = int(tm[pid]["pred_hel"]) if pid in tm else None
        xr, nconf = num(rk.get("max_fungal_crossreact_pid")), num(rk.get("n_confounder_genera"))
        hum = num(rk.get("human_homolog_pid"))
        up = l48 is not None and p48 is not None and l48 >= UP_LOG2FC and p48 < UP_PADJ
        row = {
            "gene_id": g, "product": p["product"], "n_proteins": len(prots), "protein_id": pid, "length": L,
            "tpm_mycelia": fmt(st.mean(t[0:2]), 4) if t else "",
            "tpm_spherule48h": fmt(st.mean(t[2:4]), 4) if t else "",
            "tpm_spherule8d": fmt(st.mean(t[4:6]), 4) if t else "",
            "log2fc_48h": fmt(l48), "padj_48h": fmt(p48), "log2fc_8d": fmt(l8), "padj_8d": fmt(p8),
            "signalp_call": sc.prediction, "signalp_sp_prob": fmt(sc.sp_prob),
            "predgpi_call": gc.call, "predgpi_fpr": fmt(gc.fpr), "tm_helices": fmt(tmh),
            "ser_thr_frac": fmt(fp.ser_thr_fraction(s)), "cys_frac": fmt(s.count("C") / L),
            "pro_frac": fmt(s.count("P") / L),
            "orthogroup": rk.get("orthogroup", ""), "prevalence": rk.get("prevalence", ""),
            "copy_number_mean": rk.get("copy_number_mean", ""),
            "max_fungal_crossreact_pid": rk.get("max_fungal_crossreact_pid", ""),
            "n_confounder_genera": rk.get("n_confounder_genera", ""),
            "human_homolog_pid": rk.get("human_homolog_pid", ""), "anchor": rk.get("anchor", ""),
            "flag_spherule_up": "yes" if up else "no",
            "flag_spherule_extreme": "yes" if up and l48 >= EXTREME_LOG2FC else "no",
            "flag_signal_peptide": "yes" if sc.prediction == "SP" else "no",
            "flag_gpi_call": "yes" if gc.call not in ("none", "too_short") else "no",
            "flag_tm": "" if tmh is None else ("yes" if tmh >= 1 else "no"),
            "flag_cocci_specific": (
                "" if not rk else "yes" if (xr == 0 and nconf == 0 and hum == 0) else "no"
            ),
        }  # fmt: skip
        rows.append(row)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n = len(rows)
    count = lambda k, v="yes": sum(1 for r in rows if r[k] == v)  # noqa: E731
    print(
        f"{n} genes; in DESeq2: {sum(1 for r in rows if r['log2fc_48h'])}; in ranking: {sum(1 for r in rows if r['orthogroup'])}"
    )
    print(
        f"spherule_up {count('flag_spherule_up')}, extreme {count('flag_spherule_extreme')}, SP {count('flag_signal_peptide')}, GPI {count('flag_gpi_call')}, TM>=1 {count('flag_tm')}, cocci_specific {count('flag_cocci_specific')}"
    )


if __name__ == "__main__":
    main()
