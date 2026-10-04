# C. immitis spherule table: expression, surface features and specificity

*2026-10-03. Branch `cocci-spherule-table`. Scripts and README: `analysis/cocci_spherule/`. All counts
below are computed from the gene table by `analysis/cocci_spherule/03_build_report.py`.*

## 1. Question and answer in short

Question: which *Coccidioides immitis* genes are strongly up in spherules, are predicted to sit at the
cell surface or to be secreted, and are specific to *Coccidioides*?

- 564 of 9630 genes with a DESeq2 result and a protein are up at 48 h spherule (log2 fold change at least 2, DESeq2 padj below
  0.05). 36 of them have log2 fold change at least 6.
- Only 23 of the 564 spherule-up genes have a SignalP signal peptide
  (4.1%). The rate for all genes is 4.7%
  (458 of 9757). Spherule-up genes are not enriched for signal peptides.
- 4 spherule-up genes have a signal peptide, no TM helix, and no detectable hit in a
  confounder fungus or in human (the "specific" flag): CIMG_04613, CIMG_12717, CIMG_01980, CIMG_00509.
  SOWgp is one of them.
- Of the 36 extreme genes, 28 are uncharacterized, 27 pass the
  specificity flag, and 2 have a signal peptide. The most
  induced *Coccidioides*-specific genes mostly have no signal peptide and no GPI call
  (26 of 27).

- **Cys-rich genes.** "Cys-rich" means a Cys fraction at or above the genome-wide 95th percentile
  (3.9%; the median is 1.3%). 41 of the 564 spherule-up
  genes are Cys-rich. At 5% you would expect about 28 by chance. Only 4 of
  the 41 have a signal peptide. 35 are shorter than 250 aa, and
  25 have no signal peptide and pass the specificity flag.

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

One protein per gene (the longest). 9757 genes have a protein. 8986 are in the ranking. The
others show `NA` for the specificity flag.

Flag definitions (display thresholds, not validated cut-offs): spherule up = log2 fold change at least 2
and padj below 0.05 at 48 h. Extreme = log2 fold change at least 6. Specific = highest confounder identity
0, no confounder genus, human identity 0. Signal peptide = SignalP call `SP`.

## 4. Check with the known anchors

Five characterized proteins are in the ranking as anchors. This is a sanity check, not a validation (n = 5).

| Gene | Anchor | log2FC 48 h | SignalP | PredGPI | TM | Specific |
|---|---|---|---|---|---|---|
| [CIMG_02492](https://fungidb.org/gene/CIMG_02492) | [Q2TVJ9_PRA3](https://www.uniprot.org/uniprotkb/Q2TVJ9) | -3.11 | SP | high | 0 | yes |
| [CIMG_02795](https://fungidb.org/gene/CIMG_02795) | [P0CB51_CiX1_CFantigen](https://www.uniprot.org/uniprotkb/P0CB51); [Q1E3R8_CTS1_CFantigen](https://www.uniprot.org/uniprotkb/Q1E3R8) | 1.13 | SP | none | 0 | no |
| [CIMG_04613](https://fungidb.org/gene/CIMG_04613) | [Q96V71_SOWgp82](https://www.uniprot.org/uniprotkb/Q96V71); [Q8NK60_SOWgp58](https://www.uniprot.org/uniprotkb/Q8NK60); [Q8NK61_SOWgp66](https://www.uniprot.org/uniprotkb/Q8NK61) | 9.09 | SP | high | 0 | yes |
| [CIMG_09560](https://fungidb.org/gene/CIMG_09560) | [Q6K1L8_PRA2](https://www.uniprot.org/uniprotkb/Q6K1L8) | -4.86 | SP | none | 0 | no |
| [CIMG_09696](https://fungidb.org/gene/CIMG_09696) | [A0A0E1RVD3_PRA_Cimm_RS](https://www.uniprot.org/uniprotkb/A0A0E1RVD3); [Q12295_Ag2_PRA](https://www.uniprot.org/uniprotkb/Q12295) | -2.79 | SP | high | 0 | no |

SOWgp is up and the PRA family is down. All five have a signal peptide. The RefSeq SOWgp protein is
324 aa. This is a collapsed gene model, so repeat
counts from this annotation are not meaningful (see `analysis/cocci_repeats`).

## 5. Spherule-up genes with a signal peptide (23)

| Gene | Product | aa | TPM myc | TPM 48 h | log2FC 48 h | padj | SignalP | PredGPI | TM | Cys frac | Prevalence | Specific |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [CIMG_04613](https://fungidb.org/gene/CIMG_04613) | uncharacterized protein CIMG_04613 | 324 | 13.1 | 1.5e+04 | 9.09 | 0 | SP | high | 0 | 0.065 | 0.92 | yes |
| [CIMG_10032](https://fungidb.org/gene/CIMG_10032) | expression library immunization antige | 224 | 18 | 5.14e+03 | 7.16 | 0 | SP | high | 0 | 0.018 | 0.996 | no |
| [CIMG_12717](https://fungidb.org/gene/CIMG_12717) | uncharacterized protein CIMG_12717 | 208 | 0 | 5.26 | 5.24 | 0.0028 | SP | none | 0 | 0.038 | 0.002 | yes |
| [CIMG_06994](https://fungidb.org/gene/CIMG_06994) | Cu-Zn superoxide dismutase | 248 | 17.3 | 1.14e+03 | 5.09 | 5.4e-236 | SP | high | 0 | 0.024 | 0.994 | no |
| [CIMG_07738](https://fungidb.org/gene/CIMG_07738) | uncharacterized protein CIMG_07738 | 209 | 8.98 | 559 | 4.98 | 5.2e-157 | SP | high | 0 | 0.019 | 0.998 | no |
| [CIMG_01980](https://fungidb.org/gene/CIMG_01980) | uncharacterized protein CIMG_01980 | 164 | 4.45 | 151 | 4.01 | 1.6e-40 | SP | none | 0 | 0.037 | 0.975 | yes |
| [CIMG_04143](https://fungidb.org/gene/CIMG_04143) | aspartyl protease 4 | 494 | 4.71 | 98.8 | 3.45 | 2.8e-63 | SP | high | 0 | 0.0081 | 0.998 | no |
| [CIMG_04351](https://fungidb.org/gene/CIMG_04351) | cellobiose dehydrogenase | 421 | 3.01 | 54 | 3.19 | 1.3e-33 | SP | none | 5 | 0.0071 | 0.992 | no |
| [CIMG_05585](https://fungidb.org/gene/CIMG_05585) | uncharacterized protein CIMG_05585 | 452 | 4.88 | 87.2 | 3.19 | 9.9e-57 | SP | none | 0 | 0.0089 | 0.978 | no |
| [CIMG_00509](https://fungidb.org/gene/CIMG_00509) | uncharacterized protein CIMG_00509 | 99 | 5.66 | 108 | 3.15 | 7.6e-29 | SP | none | 0 | 0.081 | 0.764 | yes |
| [CIMG_05057](https://fungidb.org/gene/CIMG_05057) | uncharacterized protein CIMG_05057 | 305 | 4.68 | 82.8 | 3.13 | 1.3e-42 | SP | none | 0 | 0.0033 | 0.988 | no |
| [CIMG_09471](https://fungidb.org/gene/CIMG_09471) | tyrosinase central domain-containing p | 414 | 1.08 | 18.7 | 3.05 | 5.6e-17 | SP | none | 0 | 0.015 | 0.992 | no |
| [CIMG_07502](https://fungidb.org/gene/CIMG_07502) | uncharacterized protein CIMG_07502 | 191 | 30.3 | 459 | 2.98 | 1.4e-63 | SP | high | 0 | 0.058 | 0.996 | no |
| [CIMG_03314](https://fungidb.org/gene/CIMG_03314) | alpha-mannosidase 1 | 519 | 11.2 | 159 | 2.83 | 1.2e-64 | SP | none | 1 | 0.0096 | 0.986 | no |
| [CIMG_07583](https://fungidb.org/gene/CIMG_07583) | tyrosinase | 658 | 41.8 | 487 | 2.54 | 3.7e-82 | SP | none | 0 | 0.0061 | 0.988 | no |
| [CIMG_03617](https://fungidb.org/gene/CIMG_03617) | uncharacterized protein CIMG_03617 | 76 | 1.27 | 14.6 | 2.52 | 3.2e-05 | SP | none | 1 | 0.053 |  | NA |
| [CIMG_00659](https://fungidb.org/gene/CIMG_00659) | uncharacterized protein CIMG_00659 | 237 | 4.83 | 50.3 | 2.43 | 9.6e-28 | SP | high | 0 | 0.013 | 0.992 | no |
| [CIMG_08973](https://fungidb.org/gene/CIMG_08973) | FRE family ferric-chelate reductase | 778 | 2.95 | 29.3 | 2.35 | 8e-21 | SP | none | 6 | 0.015 | 0.996 | no |
| [CIMG_12068](https://fungidb.org/gene/CIMG_12068) | multicopper oxidase | 636 | 4.78 | 45.8 | 2.32 | 4.7e-31 | SP | none | 0 | 0.013 | 0.99 | no |
| [CIMG_12472](https://fungidb.org/gene/CIMG_12472) | recyclin-1 | 214 | 5.09 | 51.5 | 2.32 | 6.1e-15 | SP | none | 0 | 0.019 | 0.996 | no |
| [CIMG_00364](https://fungidb.org/gene/CIMG_00364) | uncharacterized protein CIMG_00364 | 524 | 16.5 | 156 | 2.29 | 3.6e-42 | SP | high | 5 | 0.0038 | 0.883 | no |
| [CIMG_01112](https://fungidb.org/gene/CIMG_01112) | uncharacterized protein CIMG_01112 | 255 | 23.7 | 191 | 2.03 | 2.4e-30 | SP | none | 0 | 0.012 | 0.996 | no |
| [CIMG_02883](https://fungidb.org/gene/CIMG_02883) | uncharacterized protein CIMG_02883 | 257 | 25.3 | 204 | 2.01 | 1.3e-24 | SP | high | 0 | 0.031 | 0.936 | no |

## 6. Cys-rich spherule-up genes (41)

Cys fraction is the share of Cys in the longest protein of the gene. Known *Coccidioides* surface proteins
are in the Cys-rich range: SOWgp 0.065,
PRA3 0.054, PRA2
0.073. The count of 41 against
about 28 expected is a rough comparison. I did not test it. Spherule-up genes were chosen
by expression, and Cys content may correlate with short length. Most of these genes have no signal
peptide. SignalP can miss short Cys-rich proteins, and I did not run a CFEM or other Cys-domain HMM scan
yet.

| Gene | Product | aa | Cys | Cys frac | Pro frac | TPM 48 h | log2FC 48 h | padj | SignalP | PredGPI | TM | Prevalence | Specific |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [CIMG_03342](https://fungidb.org/gene/CIMG_03342) | uncharacterized protein CIMG_0 | 136 | 18 | 0.13 | 0.066 | 429 | 2.48 | 7e-34 | OTHER | none | 0 | 0.522 | no |
| [CIMG_00509](https://fungidb.org/gene/CIMG_00509) | uncharacterized protein CIMG_0 | 99 | 8 | 0.081 | 0.061 | 108 | 3.15 | 7.6e-29 | SP | none | 0 | 0.764 | yes |
| [CIMG_10862](https://fungidb.org/gene/CIMG_10862) | uncharacterized protein CIMG_1 | 58 | 4 | 0.069 | 0.086 | 46.1 | 3.08 | 3.9e-13 | OTHER | none | 0 |  | NA |
| [CIMG_12130](https://fungidb.org/gene/CIMG_12130) | uncharacterized protein CIMG_1 | 118 | 8 | 0.068 | 0.076 | 56.4 | 3.84 | 2.6e-36 | OTHER | none | 0 |  | NA |
| [CIMG_04613](https://fungidb.org/gene/CIMG_04613) | uncharacterized protein CIMG_0 | 324 | 21 | 0.065 | 0.16 | 1.5e+04 | 9.09 | 0 | SP | high | 0 | 0.92 | yes |
| [CIMG_09338](https://fungidb.org/gene/CIMG_09338) | uncharacterized protein CIMG_0 | 143 | 9 | 0.063 | 0.014 | 13.3 | 3.45 | 0.0074 | OTHER | none | 0 | 0.859 | yes |
| [CIMG_12945](https://fungidb.org/gene/CIMG_12945) | uncharacterized protein CIMG_1 | 130 | 8 | 0.061 | 0.1 | 57.9 | 2.91 | 1.2e-19 | OTHER | none | 0 | 0.0164 | yes |
| [CIMG_13539](https://fungidb.org/gene/CIMG_13539) | uncharacterized protein CIMG_1 | 116 | 7 | 0.06 | 0 | 6.17 | 4.16 | 0.026 | OTHER | none | 0 | 0.002 | yes |
| [CIMG_13156](https://fungidb.org/gene/CIMG_13156) | uncharacterized protein CIMG_1 | 259 | 15 | 0.058 | 0.077 | 2.65 | 4.76 | 0.013 | OTHER | none | 0 | 0.266 | yes |
| [CIMG_07502](https://fungidb.org/gene/CIMG_07502) | uncharacterized protein CIMG_0 | 191 | 11 | 0.058 | 0.052 | 459 | 2.98 | 1.4e-63 | SP | high | 0 | 0.996 | no |
| [CIMG_13151](https://fungidb.org/gene/CIMG_13151) | uncharacterized protein CIMG_1 | 280 | 16 | 0.057 | 0.029 | 2.49 | 3.41 | 0.044 | OTHER | none | 0 | 0.002 | yes |
| [CIMG_13220](https://fungidb.org/gene/CIMG_13220) | uncharacterized protein CIMG_1 | 198 | 11 | 0.056 | 0.025 | 1.67 | 3.79 | 0.022 | OTHER | none | 0 | 0.295 | yes |
| [CIMG_07389](https://fungidb.org/gene/CIMG_07389) | uncharacterized protein CIMG_0 | 150 | 8 | 0.053 | 0.033 | 7.63 | 5.07 | 0.0056 | OTHER | none | 0 | 0.199 | yes |
| [CIMG_03617](https://fungidb.org/gene/CIMG_03617) | uncharacterized protein CIMG_0 | 76 | 4 | 0.053 | 0.013 | 14.6 | 2.52 | 3.2e-05 | SP | none | 1 |  | NA |
| [CIMG_11958](https://fungidb.org/gene/CIMG_11958) | uncharacterized protein CIMG_1 | 173 | 9 | 0.052 | 0.029 | 13.7 | 4.55 | 0.004 | OTHER | none | 0 | 0.314 | yes |
| [CIMG_09105](https://fungidb.org/gene/CIMG_09105) | uncharacterized protein CIMG_0 | 139 | 7 | 0.05 | 0.079 | 46.1 | 2.58 | 2.9e-06 | OTHER | none | 0 |  | NA |
| [CIMG_04176](https://fungidb.org/gene/CIMG_04176) | uncharacterized protein CIMG_0 | 159 | 8 | 0.05 | 0.11 | 109 | 2.07 | 3e-10 | OTHER | none | 0 | 0.0041 | yes |
| [CIMG_12595](https://fungidb.org/gene/CIMG_12595) | uncharacterized protein CIMG_1 | 100 | 5 | 0.05 | 0.05 | 113 | 6.55 | 4.8e-06 | OTHER | none | 0 | 0.002 | yes |
| [CIMG_03730](https://fungidb.org/gene/CIMG_03730) | uncharacterized protein CIMG_0 | 41 | 2 | 0.049 | 0.073 | 6.27e+03 | 6.34 | 3.3e-19 | OTHER | none | 0 | 0.002 | yes |
| [CIMG_03347](https://fungidb.org/gene/CIMG_03347) | uncharacterized protein CIMG_0 | 104 | 5 | 0.048 | 0.2 | 940 | 3 | 1e-114 | OTHER | none | 0 | 0.982 | no |
| [CIMG_11663](https://fungidb.org/gene/CIMG_11663) | uncharacterized protein CIMG_1 | 84 | 4 | 0.048 | 0.036 | 221 | 2.84 | 1.4e-09 | OTHER | none | 0 |  | NA |
| [CIMG_11522](https://fungidb.org/gene/CIMG_11522) | uncharacterized protein CIMG_1 | 85 | 4 | 0.047 | 0.012 | 346 | 6.62 | 6.7e-113 | OTHER | none | 1 | 0.0492 | yes |
| [CIMG_11509](https://fungidb.org/gene/CIMG_11509) | uncharacterized protein CIMG_1 | 132 | 6 | 0.045 | 0.076 | 130 | 2.82 | 1.4e-11 | OTHER | none | 0 |  | NA |
| [CIMG_13214](https://fungidb.org/gene/CIMG_13214) | uncharacterized protein CIMG_1 | 133 | 6 | 0.045 | 0.06 | 21.3 | 6.26 | 0.00011 | OTHER | none | 0 | 0.844 | yes |
| [CIMG_03444](https://fungidb.org/gene/CIMG_03444) | uncharacterized protein CIMG_0 | 133 | 6 | 0.045 | 0.083 | 6.8 | 2.27 | 0.0047 | OTHER | none | 0 | 0.15 | yes |
| [CIMG_13082](https://fungidb.org/gene/CIMG_13082) | uncharacterized protein CIMG_1 | 111 | 5 | 0.045 | 0.099 | 17.4 | 2.02 | 0.00016 | OTHER | none | 0 | 0.0963 | yes |
| [CIMG_11847](https://fungidb.org/gene/CIMG_11847) | uncharacterized protein CIMG_1 | 178 | 8 | 0.045 | 0.034 | 24.7 | 2.12 | 1.2e-11 | OTHER | none | 0 | 0.189 | yes |
| [CIMG_11763](https://fungidb.org/gene/CIMG_11763) | uncharacterized protein CIMG_1 | 357 | 16 | 0.045 | 0.059 | 4.81 | 2.45 | 0.016 | OTHER | none | 2 | 0.709 | yes |
| [CIMG_01257](https://fungidb.org/gene/CIMG_01257) | methionine-R-sulfoxide reducta | 158 | 7 | 0.044 | 0.063 | 189 | 3.66 | 4.7e-92 | OTHER | none | 0 | 0.992 | no |
| [CIMG_10013](https://fungidb.org/gene/CIMG_10013) | uncharacterized protein CIMG_1 | 114 | 5 | 0.044 | 0.12 | 133 | 6.18 | 7.1e-37 | OTHER | none | 0 | 0.0369 | yes |
| [CIMG_12358](https://fungidb.org/gene/CIMG_12358) | uncharacterized protein CIMG_1 | 114 | 5 | 0.044 | 0.088 | 36.5 | 2.87 | 1.4e-09 | OTHER | none | 0 |  | NA |
| [CIMG_12995](https://fungidb.org/gene/CIMG_12995) | uncharacterized protein CIMG_1 | 116 | 5 | 0.043 | 0.043 | 14 | 5.31 | 0.0018 | OTHER | none | 0 | 0.266 | yes |
| [CIMG_13457](https://fungidb.org/gene/CIMG_13457) | uncharacterized protein CIMG_1 | 95 | 4 | 0.042 | 0.032 | 15.9 | 2.25 | 2.6e-05 | OTHER | none | 0 | 0.221 | yes |
| [CIMG_12808](https://fungidb.org/gene/CIMG_12808) | uncharacterized protein CIMG_1 | 143 | 6 | 0.042 | 0.1 | 7.21 | 6.63 | 2.6e-05 | OTHER | none | 0 | 0.0143 | yes |
| [CIMG_07724](https://fungidb.org/gene/CIMG_07724) | uncharacterized protein CIMG_0 | 97 | 4 | 0.041 | 0.083 | 34.3 | 2.48 | 0.0041 | OTHER | none | 0 |  | NA |
| [CIMG_12971](https://fungidb.org/gene/CIMG_12971) | uncharacterized protein CIMG_1 | 272 | 11 | 0.04 | 0.081 | 5.67 | 3.57 | 0.0053 | OTHER | none | 0 | 0.266 | yes |
| [CIMG_06955](https://fungidb.org/gene/CIMG_06955) | histone-lysine N-methyltransfe | 446 | 18 | 0.04 | 0.04 | 102 | 2.82 | 9.4e-45 | OTHER | none | 0 | 0.996 | no |
| [CIMG_12913](https://fungidb.org/gene/CIMG_12913) | uncharacterized protein CIMG_1 | 75 | 3 | 0.04 | 0.067 | 18.3 | 4.54 | 0.011 | OTHER | none | 0 | 0.0061 | yes |
| [CIMG_13084](https://fungidb.org/gene/CIMG_13084) | uncharacterized protein CIMG_1 | 200 | 8 | 0.04 | 0.07 | 18.4 | 3.6 | 1.4e-05 | OTHER | none | 0 | 0.65 | yes |
| [CIMG_03569](https://fungidb.org/gene/CIMG_03569) | uncharacterized protein CIMG_0 | 100 | 4 | 0.04 | 0.1 | 7.96 | 2.48 | 0.0027 | OTHER | none | 1 | 0.111 | yes |
| [CIMG_06169](https://fungidb.org/gene/CIMG_06169) | uncharacterized protein CIMG_0 | 127 | 5 | 0.039 | 0.087 | 31.5 | 2.62 | 0.0032 | OTHER | none | 0 |  | NA |

## 7. Extreme genes: log2FC at least 6 and padj below 0.05 (36)

| Gene | Product | aa | TPM myc | TPM 48 h | log2FC 48 h | padj | SignalP | PredGPI | TM | Cys frac | Prevalence | Specific |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [CIMG_04613](https://fungidb.org/gene/CIMG_04613) | uncharacterized protein CIMG_04613 | 324 | 13.1 | 1.5e+04 | 9.09 | 0 | SP | high | 0 | 0.065 | 0.92 | yes |
| [CIMG_04662](https://fungidb.org/gene/CIMG_04662) | magnesium-translocating P-type ATPase | 929 | 3.33 | 528 | 8.85 | 1.9e-176 | OTHER | none | 8 | 0.018 | 0.99 | yes |
| [CIMG_10264](https://fungidb.org/gene/CIMG_10264) | aldehyde reductase | 314 | 3.17 | 1.95e+03 | 8.33 | 7.4e-234 | OTHER | none | 0 | 0.0032 | 0.984 | no |
| [CIMG_13657](https://fungidb.org/gene/CIMG_13657) | uncharacterized protein CIMG_13657 | 162 | 0 | 61.3 | 8.26 | 8.9e-08 | OTHER | none | 0 | 0.012 | 0.0615 | yes |
| [CIMG_09753](https://fungidb.org/gene/CIMG_09753) | ABC multidrug transporter | 1501 | 2.07 | 942 | 7.86 | 0 | OTHER | none | 16 | 0.014 | 0.996 | no |
| [CIMG_10037](https://fungidb.org/gene/CIMG_10037) | copper transporter | 193 | 3.57 | 1.08e+03 | 7.23 | 6.7e-226 | OTHER | none | 2 | 0.0052 | 0.996 | no |
| [CIMG_10877](https://fungidb.org/gene/CIMG_10877) | uncharacterized protein CIMG_10877 | 163 | 0.239 | 87.5 | 7.22 | 3.8e-07 | OTHER | none | 0 | 0.0061 | 0.002 | yes |
| [CIMG_12615](https://fungidb.org/gene/CIMG_12615) | uncharacterized protein CIMG_12615 | 301 | 0 | 14.9 | 7.19 | 4.7e-06 | OTHER | none | 0 | 0.013 | 0.0615 | yes |
| [CIMG_10032](https://fungidb.org/gene/CIMG_10032) | expression library immunization antige | 224 | 18 | 5.14e+03 | 7.16 | 0 | SP | high | 0 | 0.018 | 0.996 | no |
| [CIMG_12484](https://fungidb.org/gene/CIMG_12484) | uncharacterized protein CIMG_12484 | 124 | 1.56 | 458 | 7.16 | 7.2e-21 | OTHER | none | 0 | 0.016 | 0.0779 | yes |
| [CIMG_08364](https://fungidb.org/gene/CIMG_08364) | uncharacterized protein CIMG_08364 | 179 | 0.217 | 74.9 | 7.12 | 5.4e-07 | OTHER | none | 0 | 0.022 | 0.352 | yes |
| [CIMG_13332](https://fungidb.org/gene/CIMG_13332) | uncharacterized protein CIMG_13332 | 263 | 0 | 16.8 | 7.1 | 7.4e-06 | OTHER | none | 0 | 0.0076 | 0.0615 | yes |
| [CIMG_01115](https://fungidb.org/gene/CIMG_01115) | thiamine thiazole synthase | 328 | 47.5 | 2.96e+03 | 6.98 | 0 | OTHER | none | 0 | 0.015 | 0.402 | no |
| [CIMG_08103](https://fungidb.org/gene/CIMG_08103) | opsin 1 | 289 | 5.35 | 1.32e+03 | 6.98 | 5.2e-278 | OTHER | none | 7 | 0.0035 | 0.998 | no |
| [CIMG_02628](https://fungidb.org/gene/CIMG_02628) | Arp2/3 complex subunit Arc16 | 312 | 0.632 | 165 | 6.96 | 2.1e-33 | OTHER | none | 0 | 0.0096 | 0.994 | no |
| [CIMG_12573](https://fungidb.org/gene/CIMG_12573) | uncharacterized protein CIMG_12573 | 238 | 0 | 16.4 | 6.86 | 1.5e-05 | OTHER | none | 0 | 0.0084 | 0.0615 | yes |
| [CIMG_12808](https://fungidb.org/gene/CIMG_12808) | uncharacterized protein CIMG_12808 | 143 | 1.61e-09 | 7.21 | 6.63 | 2.6e-05 | OTHER | none | 0 | 0.042 | 0.0143 | yes |
| [CIMG_11522](https://fungidb.org/gene/CIMG_11522) | uncharacterized protein CIMG_11522 | 85 | 1.76 | 346 | 6.62 | 6.7e-113 | OTHER | none | 1 | 0.047 | 0.0492 | yes |
| [CIMG_13230](https://fungidb.org/gene/CIMG_13230) | uncharacterized protein CIMG_13230 | 104 | 0.959 | 215 | 6.62 | 1.1e-09 | OTHER | none | 0 | 0.019 | 0.201 | yes |
| [CIMG_10007](https://fungidb.org/gene/CIMG_10007) | uncharacterized protein CIMG_10007 | 125 | 1.35 | 255 | 6.59 | 3.5e-107 | OTHER | none | 0 | 0.008 | 0.17 | yes |
| [CIMG_12595](https://fungidb.org/gene/CIMG_12595) | uncharacterized protein CIMG_12595 | 100 | 0.515 | 113 | 6.55 | 4.8e-06 | OTHER | none | 0 | 0.05 | 0.002 | yes |
| [CIMG_06250](https://fungidb.org/gene/CIMG_06250) | uncharacterized protein CIMG_06250 | 320 | 3.78 | 692 | 6.53 | 4.9e-241 | OTHER | none | 1 | 0.019 | 1 | yes |
| [CIMG_13004](https://fungidb.org/gene/CIMG_13004) | uncharacterized protein CIMG_13004 | 289 | 0.32 | 65.8 | 6.48 | 4.9e-13 | OTHER | none | 0 | 0.01 | 0.0615 | yes |
| [CIMG_12994](https://fungidb.org/gene/CIMG_12994) | uncharacterized protein CIMG_12994 | 130 | 0 | 25.5 | 6.47 | 6.7e-05 | OTHER | none | 0 | 0 | 0.002 | yes |
| [CIMG_12953](https://fungidb.org/gene/CIMG_12953) | uncharacterized protein CIMG_12953 | 170 | 0 | 17.3 | 6.4 | 7.6e-05 | OTHER | none | 0 | 0.012 | 0.002 | yes |
| [CIMG_12716](https://fungidb.org/gene/CIMG_12716) | uncharacterized protein CIMG_12716 | 77 | 0 | 62.9 | 6.38 | 0.00011 | OTHER | none | 0 | 0 | 0.002 | yes |
| [CIMG_13499](https://fungidb.org/gene/CIMG_13499) | uncharacterized protein CIMG_13499 | 205 | 0.168 | 33.4 | 6.37 | 9.4e-06 | OTHER | none | 0 | 0.019 | 0.0061 | yes |
| [CIMG_03730](https://fungidb.org/gene/CIMG_03730) | uncharacterized protein CIMG_03730 | 41 | 510 | 6.27e+03 | 6.34 | 3.3e-19 | OTHER | none | 0 | 0.049 | 0.002 | yes |
| [CIMG_13125](https://fungidb.org/gene/CIMG_13125) | uncharacterized protein CIMG_13125 | 139 | 0 | 19.9 | 6.28 | 0.00015 | OTHER | none | 0 | 0.029 | 0.002 | yes |
| [CIMG_13214](https://fungidb.org/gene/CIMG_13214) | uncharacterized protein CIMG_13214 | 133 | 0 | 21.3 | 6.26 | 0.00011 | OTHER | none | 0 | 0.045 | 0.844 | yes |
| [CIMG_10013](https://fungidb.org/gene/CIMG_10013) | uncharacterized protein CIMG_10013 | 114 | 0.919 | 133 | 6.18 | 7.1e-37 | OTHER | none | 0 | 0.044 | 0.0369 | yes |
| [CIMG_07662](https://fungidb.org/gene/CIMG_07662) | uncharacterized protein CIMG_07662 | 91 | 0 | 36.9 | 6.07 | 0.00028 | OTHER | none | 0 | 0.033 |  | NA |
| [CIMG_12654](https://fungidb.org/gene/CIMG_12654) | uncharacterized protein CIMG_12654 | 92 | 0.601 | 93.8 | 6.05 | 3.3e-05 | OTHER | none | 0 | 0.033 | 0.002 | yes |
| [CIMG_09539](https://fungidb.org/gene/CIMG_09539) | uncharacterized protein CIMG_09539 | 387 | 22.3 | 2.78e+03 | 6.02 | 0 | OTHER | none | 0 | 0.01 | 0.996 | no |
| [CIMG_12854](https://fungidb.org/gene/CIMG_12854) | uncharacterized protein CIMG_12854 | 253 | 0 | 8.61 | 6.02 | 0.00026 | OTHER | none | 0 | 0.028 | 0.998 | yes |
| [CIMG_12625](https://fungidb.org/gene/CIMG_12625) | uncharacterized protein CIMG_12625 | 121 | 1.26 | 161 | 6 | 2.5e-11 | OTHER | none | 0 | 0.017 | 0.0615 | yes |

The full list of 564 spherule-up genes is `analysis/cocci_spherule/shortlist_spherule_up.tsv`. The
full table with 9757 genes is `analysis/cocci_spherule/spherule_surface_table.tsv.gz`.

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
  for them. 1 extreme gene is not in the ranking.
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

Check the 26 specific extreme genes without a signal peptide and the 25
Cys-rich spherule-up genes without a signal peptide that pass the specificity flag: gene-model support
(RS2 and pangenome models), orthologs in *C. posadasii* and other Onygenales with a more sensitive
search, Cys-domain HMM scans (CFEM and others), and sequence features. This is the next work item on this
branch.
