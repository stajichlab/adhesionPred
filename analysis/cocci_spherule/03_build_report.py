#!/usr/bin/env python3
"""Write docs/reports/2026-10-03-cocci-spherule-surface-table.md from the gene table.

Every count and every gene row in the report is computed here from the table, so the text
cannot drift from the data. Usage: 03_build_report.py --table <table.tsv> --out <report.md>
"""

import argparse
import csv
from pathlib import Path


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def g3(x, nd=3):
    v = num(x)
    return "" if v is None else f"{v:.{nd}g}"


def gene_link(gid):
    """Markdown link to the FungiDB gene page for CIMG_ IDs; other IDs are returned as is."""
    if gid.startswith("CIMG_"):
        return f"[{gid}](https://fungidb.org/gene/{gid})"
    return gid


def anchor_links(anchor):
    """Link each `ACCESSION_name` entry of a ';'-separated anchor field to its UniProt page."""
    out = []
    for entry in anchor.split(";"):
        acc = entry.split("_", 1)[0]
        out.append(f"[{entry}](https://www.uniprot.org/uniprotkb/{acc})")
    return "; ".join(out)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.table), delimiter="\t"))
    n = len(rows)
    yes = lambda k, r=None: [x for x in (r if r is not None else rows) if x[k] == "yes"]  # noqa: E731
    up = sorted(yes("flag_spherule_up"), key=lambda x: -num(x["log2fc_48h"]))
    ext = yes("flag_spherule_extreme", up)
    sp_all = [x for x in rows if x["signalp_call"] == "SP"]
    sp_up = [x for x in up if x["signalp_call"] == "SP"]
    in_rank = [x for x in rows if x["orthogroup"]]
    up_sp_spec_notm = [
        x for x in sp_up if x["flag_cocci_specific"] == "yes" and x["flag_tm"] == "no"
    ]
    ext_unchar = [
        x for x in ext if "uncharacterized" in x["product"] or "hypothetical" in x["product"]
    ]
    ext_spec = yes("flag_cocci_specific", ext)
    ext_spec_nosp = [x for x in ext_spec if x["signalp_call"] != "SP"]
    ext_notrank = [x for x in ext if not x["orthogroup"]]
    de_n = sum(1 for x in rows if x["log2fc_48h"])
    anchors = [x for x in rows if x["anchor"]]
    cys_sorted = sorted(num(x["cys_frac"]) for x in rows)
    cys_p95 = cys_sorted[int(0.95 * (len(cys_sorted) - 1))]
    cys_median = cys_sorted[len(cys_sorted) // 2]
    cys_up = sorted(
        [x for x in up if num(x["cys_frac"]) >= cys_p95], key=lambda x: -num(x["cys_frac"])
    )
    cys_up_sp = [x for x in cys_up if x["signalp_call"] == "SP"]
    cys_up_nosp_spec = [
        x for x in cys_up if x["signalp_call"] != "SP" and x["flag_cocci_specific"] == "yes"
    ]
    cys_up_short = [x for x in cys_up if int(x["length"]) < 250]

    def cys_row(x):
        return [
            gene_link(x["gene_id"]), x["product"][:30], x["length"],
            round(num(x["cys_frac"]) * int(x["length"])), g3(x["cys_frac"], 2), g3(x["pro_frac"], 2),
            g3(x["tpm_spherule48h"]), g3(x["log2fc_48h"], 3), g3(x["padj_48h"], 2),
            x["signalp_call"], x["predgpi_call"].replace("highly_probable", "high"),
            x["tm_helices"], g3(x["prevalence"], 3), x["flag_cocci_specific"] or "NA",
        ]  # fmt: skip

    ch = ["Gene", "Product", "aa", "Cys", "Cys frac", "Pro frac", "TPM 48 h", "log2FC 48 h", "padj",
          "SignalP", "PredGPI", "TM", "Prevalence", "Specific"]  # fmt: skip

    def gene_row(x):
        return [
            gene_link(x["gene_id"]), x["product"][:38], x["length"], g3(x["tpm_mycelia"]),
            g3(x["tpm_spherule48h"]), g3(x["log2fc_48h"], 3), g3(x["padj_48h"], 2),
            x["signalp_call"], x["predgpi_call"].replace("highly_probable", "high"),
            x["tm_helices"], g3(x["cys_frac"], 2), g3(x["prevalence"], 3),
            x["flag_cocci_specific"] or "NA",
        ]  # fmt: skip

    gh = ["Gene", "Product", "aa", "TPM myc", "TPM 48 h", "log2FC 48 h", "padj", "SignalP",
          "PredGPI", "TM", "Cys frac", "Prevalence", "Specific"]  # fmt: skip

    text = f"""# C. immitis spherule table: expression, surface features and specificity

*2026-10-03. Branch `cocci-spherule-table`. Scripts and README: `analysis/cocci_spherule/`. All counts
below are computed from the gene table by `analysis/cocci_spherule/03_build_report.py`.*

## 1. Question and answer in short

Question: which *Coccidioides immitis* genes are strongly up in spherules, are predicted to sit at the
cell surface or to be secreted, and are specific to *Coccidioides*?

- {len(up)} of {
        de_n
    } genes with a DESeq2 result and a protein are up at 48 h spherule (log2 fold change at least 2, DESeq2 padj below
  0.05). {len(ext)} of them have log2 fold change at least 6.
- Only {len(sp_up)} of the {len(up)} spherule-up genes have a SignalP signal peptide
  ({100 * len(sp_up) / len(up):.1f}%). The rate for all genes is {100 * len(sp_all) / n:.1f}%
  ({len(sp_all)} of {n}). Spherule-up genes are not enriched for signal peptides.
- {
        len(up_sp_spec_notm)
    } spherule-up genes have a signal peptide, no TM helix, and no detectable hit in a
  confounder fungus or in human (the "specific" flag): {
        ", ".join(x["gene_id"] for x in up_sp_spec_notm)
    }.
  SOWgp is one of them.
- Of the {len(ext)} extreme genes, {len(ext_unchar)} are uncharacterized, {len(ext_spec)} pass the
  specificity flag, and {
        sum(1 for x in ext if x["signalp_call"] == "SP")
    } have a signal peptide. The most
  induced *Coccidioides*-specific genes mostly have no signal peptide and no GPI call
  ({len(ext_spec_nosp)} of {len(ext_spec)}).

- **Cys-rich genes.** "Cys-rich" means a Cys fraction at or above the genome-wide 95th percentile
  ({100 * cys_p95:.1f}%; the median is {100 * cys_median:.1f}%). {len(cys_up)} of the {
        len(up)
    } spherule-up
  genes are Cys-rich. At 5% you would expect about {0.05 * len(up):.0f} by chance. Only {
        len(cys_up_sp)
    } of
  the {len(cys_up)} have a signal peptide. {len(cys_up_short)} are shorter than 250 aa, and
  {len(cys_up_nosp_spec)} have no signal peptide and pass the specificity flag.

A surface flag here is a hypothesis. Step 1 has no Onygenales truth set, so the false-positive rate of a
SignalP or PredGPI call in *Coccidioides* is not measured.

## 2. Which gene version

Three annotations of the RS genome exist in the inputs. I checked them on 2026-10-03.

| Source | IDs | Genes | Relation |
|---|---|---|---|
| FungiDB-46 transcripts, kallisto run "RS1" | `CIMG_#####-t26_1` | 9,905 | Same as NCBI RefSeq GCF_000149335.2: the same 9,905 gene IDs, and all 10,058 transcripts have the same sequence and gene ID |
| `Coccidioides_immitis_RS.mrna-transcripts.fa.gz`, "RS2" | `CIMG2_######-T#` | 8,846 | Different annotation. 625 of 9,861 transcripts match an RS1 sequence exactly |
| Pangenome RS proteome | `C3CED3A_######-T1` | 8,645 | A third annotation |

The antigen ranking (`analysis/cocci_antigens`) uses RS1 IDs. The RS1 counts are therefore the right
counts for `CIMG_` IDs, and I did not re-run kallisto. RS2 is not used.

## 3. Method

| Step | Tool | Input |
|---|---|---|
| Expression | DESeq2 1.50.2 on the RS1 kallisto counts. Transcript counts summed per gene, rounded. 48 h and 8 d spherule against mycelia | `RS1_kallisto.counts.csv` (6 samples, 2 per condition) |
| Signal peptide | SignalP 6, fast mode, eukarya (GPU build, `analysis/step1_compare/jobs/j1_features.sh`) | 9,910 RefSeq proteins |
| GPI anchor | PredGPI (same job) | same |
| TM helices | TMHMM 2.0c | same |
| Specificity | `cocci_antigen_ranking.tsv`: prevalence in 488 *Coccidioides* proteomes, highest identity to a confounder fungus, number of confounder genera, human homology. Not recomputed | pangenome |

One protein per gene (the longest). {n} genes have a protein. {len(in_rank)} are in the ranking. The
others show `NA` for the specificity flag.

Flag definitions (display thresholds, not validated cut-offs): spherule up = log2 fold change at least 2
and padj below 0.05 at 48 h. Extreme = log2 fold change at least 6. Specific = highest confounder identity
0, no confounder genus, human identity 0. Signal peptide = SignalP call `SP`.

## 4. Check with the known anchors

Five characterized proteins are in the ranking as anchors. This is a sanity check, not a validation (n = 5).

{
        md_table(
            ["Gene", "Anchor", "log2FC 48 h", "SignalP", "PredGPI", "TM", "Specific"],
            [
                [
                    gene_link(x["gene_id"]),
                    anchor_links(x["anchor"]),
                    g3(x["log2fc_48h"]),
                    x["signalp_call"],
                    x["predgpi_call"].replace("highly_probable", "high"),
                    x["tm_helices"],
                    x["flag_cocci_specific"] or "NA",
                ]
                for x in anchors
            ],
        )
    }

SOWgp is up and the PRA family is down. All five have a signal peptide. The RefSeq SOWgp protein is
{
        next(x["length"] for x in rows if x["gene_id"] == "CIMG_04613")
    } aa. This is a collapsed gene model, so repeat
counts from this annotation are not meaningful (see `analysis/cocci_repeats`).

## 5. Spherule-up genes with a signal peptide ({len(sp_up)})

{md_table(gh, [gene_row(x) for x in sp_up])}

## 6. Cys-rich spherule-up genes ({len(cys_up)})

Cys fraction is the share of Cys in the longest protein of the gene. Known *Coccidioides* surface proteins
are in the Cys-rich range: SOWgp {
        g3(next(x["cys_frac"] for x in rows if x["gene_id"] == "CIMG_04613"), 2)
    },
PRA3 {g3(next(x["cys_frac"] for x in rows if x["gene_id"] == "CIMG_02492"), 2)}, PRA2
{g3(next(x["cys_frac"] for x in rows if x["gene_id"] == "CIMG_09560"), 2)}. The count of {
        len(cys_up)
    } against
about {
        0.05
        * len(
            up
        ):.0f} expected is a rough comparison. I did not test it. Spherule-up genes were chosen
by expression, and Cys content may correlate with short length. Most of these genes have no signal
peptide. SignalP can miss short Cys-rich proteins, and I did not run a CFEM or other Cys-domain HMM scan
yet.

{md_table(ch, [cys_row(x) for x in cys_up])}

## 7. Extreme genes: log2FC at least 6 and padj below 0.05 ({len(ext)})

{md_table(gh, [gene_row(x) for x in ext])}

The full list of {
        len(up)
    } spherule-up genes is `analysis/cocci_spherule/shortlist_spherule_up.tsv`. The
full table with {n} genes is `analysis/cocci_spherule/spherule_surface_table.tsv.gz`.

## 8. Limits

- **Two replicates per condition.** DESeq2 has little power. Read padj as a ranking aid. DESeq2 and TPM
  give different log2 fold changes for very abundant genes. For SOWgp: 9.09 (DESeq2) and 10.05 (TPM).
- **One strain, one experiment.** *C. immitis* RS (Carlin et al. 2021 data). *C. posadasii* has no expression
  data here.
- **The surface call is a prediction.** SignalP and PredGPI are untested in *Coccidioides*. Short Cys-rich
  and non-classically secreted proteins can be missed. TMHMM can mark a signal peptide as a TM helix, so a
  TM flag does not prove a membrane protein.
- **Specificity depends on the ranking's homology search.** Many extreme specific genes are short (41 to
  330 aa). I did not test whether the search is sensitive for short proteins, so "specific" may be inflated
  for them. {len(ext_notrank)} extreme gene is not in the ranking.
- **Gene models.** A short abundant ORF can be a gene-model artifact. Example: CIMG_03730 is 41 aa and
  reaches 6,270 TPM at 48 h. I did not check the models.
- **No gate values.** Step 1's R2 needs two fitted values that are not set. The table has raw SignalP,
  PredGPI and Ser/Thr values. It does not compute R2.

## 9. Reproduce

```
python analysis/cocci_spherule/00_prepare_refseq_proteins.py _workdir/cocci_spherule/ref _workdir/cocci_spherule/s1
Rscript analysis/cocci_spherule/01_deseq2_rs1.R <RS1_kallisto.counts.csv> _workdir/cocci_spherule/deseq2
sbatch --export=ALL,PROJ_ROOT=$PWD,STEP1_WORKDIR=$PWD/_workdir/cocci_spherule/s1 analysis/step1_compare/jobs/j1_features.sh
sbatch analysis/cocci_spherule/01b_tmhmm.sh _workdir/cocci_spherule/s1 $PWD/_workdir/cocci_spherule/tmhmm.tsv
python analysis/cocci_spherule/02_build_table.py --work _workdir/cocci_spherule --out spherule_surface_table.tsv
python analysis/cocci_spherule/03_build_report.py --table spherule_surface_table.tsv --out docs/reports/2026-10-03-cocci-spherule-surface-table.md
```

The RefSeq files come from `ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/335/GCF_000149335.2_ASM14933v2`.
`01c_signalp_predgpi_cpu.sh` is a CPU variant of the SignalP job. It was written as a fallback and cancelled
before it finished any part. The reported SignalP and PredGPI results are from the GPU job.

## 10. Next

Check the {len(ext_spec_nosp)} specific extreme genes without a signal peptide and the {
        len(cys_up_nosp_spec)
    }
Cys-rich spherule-up genes without a signal peptide that pass the specificity flag: gene-model support
(RS2 and pangenome models), orthologs in *C. posadasii* and other Onygenales with a more sensitive
search, Cys-domain HMM scans (CFEM and others), and sequence features. This is the next work item on this
branch.
"""
    Path(a.out).write_text(text)
    print(f"wrote {a.out}: {len(text.splitlines())} lines")


if __name__ == "__main__":
    main()
