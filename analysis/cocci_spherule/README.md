# *Coccidioides immitis* spherule table

Question: which *C. immitis* genes are (a) strongly up in spherules, (b) predicted to sit at the cell
surface or be secreted, and (c) specific to *Coccidioides*? This table puts the three kinds of evidence
side by side, one row per gene. It does not rank genes by a single score, and it validates nothing.

## Which gene version

Three annotations of the RS genome exist in this repository's inputs. The checks below were run on
2026-10-03.

| Source | IDs | Genes | Relation to the others |
|---|---|---|---|
| FungiDB-46 transcripts (`UCSD_SpheruleMycelium/data/CimmitisRS_Broad.mRNA.fasta`, kallisto run "RS1") | `CIMG_#####-t26_1` | 9,905 | Identical to NCBI RefSeq GCF_000149335.2: same 9,905 gene IDs, and all 10,058 transcripts have the same sequence and gene ID |
| `Coccidioides_immitis_RS.mrna-transcripts.fa.gz` ("RS2") | `CIMG2_######-T#` | 8,846 | A different annotation. 625 of its 9,861 transcripts match an RS1 sequence exactly |
| Pangenome RS proteome (`Pangenome/input/Coccidioides_immitis_RS.proteins.fa`) | `C3CED3A_######-T1` | 8,645 | A third annotation |

The antigen ranking (`analysis/cocci_antigens`) uses `CIMG_#####-t26_1-p1` proteins, which are RS1. So
RS1, the RefSeq annotation and the ranking agree, and the existing RS1 kallisto counts are the right ones
for the `CIMG_` IDs. No re-quantification was needed for the ID question. RS2 is not used here.

## Method

| Step | Script | Output (in `_workdir/cocci_spherule/`) |
|---|---|---|
| RefSeq proteins and protein-to-gene map | `00_prepare_refseq_proteins.py` | `ref/protein_map.tsv`, `s1/phaseb/unique_sequences.fasta.gz` |
| Differential expression, DESeq2 on the RS1 kallisto counts, 48 h spherule and 8 d spherule against mycelia | `01_deseq2_rs1.R` | `deseq2/deseq2_rs1_gene.tsv` |
| SignalP 6 (fast, eukarya) and PredGPI | `analysis/step1_compare/jobs/j1_features.sh` (reused) | `s1/phaseb/signalp/`, `s1/phaseb/predgpi/` |
| TMHMM 2.0c | `01b_tmhmm.sh` | `tmhmm.tsv` |
| Join | `02_build_table.py` | `spherule_surface_table.tsv` |

Specificity columns come from `analysis/cocci_antigens/cocci_antigen_ranking.tsv` (prevalence in 488
*Coccidioides* proteomes, copy number, highest identity to a confounder fungus, human homology). I did not
recompute them.

## Limits

- **Two replicates per condition** (`samples.csv`). DESeq2 has little power. Read `padj` as a ranking aid.
  Library composition matters: SOWgp is a large share of the spherule reads, and DESeq2 normalised counts
  and TPM give different log2 fold changes for it (9.09 against 10.05 at 48 h).
- **One strain** (*C. immitis* RS), one experiment (Carlin et al. 2021 data). *C. posadasii* has no
  expression data here.
- **The surface call is a prediction.** Step 1 has no Onygenales truth, so the false-positive rate of a
  signal-peptide or GPI call in *Coccidioides* is unmeasured. A flag is a hypothesis.
- **No gate values.** Step 1's R2 needs two fitted values that the owner has not set, so the table has raw
  SignalP, PredGPI and Ser/Thr values and a yes/no for a SignalP call (R0) and a PredGPI call. It does not
  compute R2.
- One protein per gene (the longest). `n_proteins` shows genes with isoforms.
- The Fungi5k database holds SignalP calls for a different annotation of RS (7,630 proteins, 371 with a
  signal peptide in its `signalp` table). They are not used here, because the gene models differ.
- Flags use display thresholds, not validated cut-offs: spherule up means log2 fold change at least 2
  and `padj` below 0.05 at 48 h; extreme means at least 6.

## Results

Report with the gene tables: `docs/reports/2026-10-03-cocci-spherule-surface-table.md` (built by
`03_build_report.py`). Files in this directory:

- `spherule_surface_table.tsv.gz`: 9,757 genes, 33 columns.
- `shortlist_spherule_up.tsv`: the 564 spherule-up genes (log2 fold change at least 2, padj below 0.05).
- Cys-rich means a Cys fraction at or above the genome-wide 95th percentile (3.9%). The report lists the
  Cys-rich spherule-up genes in section 6.
- `01c_signalp_predgpi_cpu.sh` is a CPU fallback for the SignalP and PredGPI job. It was cancelled before it
  finished any part. The results come from the GPU job (`step1_compare/jobs/j1_features.sh`).

## Follow-up (branch `cocci-extreme-genes`)

Report: `docs/reports/2026-10-03-cocci-spherule-followup.md`. Scripts `10_select_followup.py`,
`11_followup_job.sh` (phmmer, diamond, hmmscan Pfam-A, Phobius), `12_summarize_followup.py`,
`13_build_followup_report.py`. Files: `followup_genes.tsv` (45 genes), `search_proteomes.tsv`,
`followup_genes_summary.tsv`, `phobius_spherule_up.tsv`. The raw search outputs are in
`_workdir/cocci_spherule/followup/results/` (git-ignored).
