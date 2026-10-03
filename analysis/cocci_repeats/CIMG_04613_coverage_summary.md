# CIMG_04613 read depth across 559 C. immitis strains

Date: 2026-10-01

## What was examined
Mean read depth over the gene CIMG_04613 (GG704914:970,094..971,690, 1,597 bp, 1-based inclusive) in each of 559 strains. The reads were aligned to the C. immitis RS reference.

## Tools and versions
- mosdepth 0.3.12 (module `mosdepth/0.3.12`; it loads samtools 1.22.1 and htslib 1.22.1)
- GNU parallel (one mosdepth process per CRAM, 2 threads each)
- Python 3.12 (`/usr/bin/python3.12`) for table assembly
- SLURM job 29335355 (partition short, 64 CPUs, 2 min 35 s, exit 0)

## Full paths
Project root: `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS`

- Alignments (559 files): `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/aln/*.cram`
- Reference: `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/genome/FungiDB-68_CimmitisRS_Genome.fasta`
- Config: `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/config.txt`
- Gene BED (0-based): `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/coverage/mosdepth_gene/CIMG_04613.bed` (GG704914, 970093, 971690)
- Run script: `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/pipeline/11b_mosdepth_gene.sh`
- Table builder: `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/scripts/gene_coverage_table.py`
- Per-strain mosdepth output: `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/coverage/mosdepth_gene/<strain>.{all,Q20}.regions.bed.gz`
- Genome-wide baseline: `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/coverage/mosdepth/<strain>.5000bp.mosdepth.summary.txt` (the `total` row, mean column). These come from the earlier run of `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/pipeline/11_mosdepth.sh`.
- Result table (original): `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/summary_stats/CIMG_04613.coverage.tsv`
- Result table (copy): `/rhome/jstajich/projects/adhesionPred/analysis/cocci_repeats/CIMG_04613.coverage.tsv`

## Method
For each CRAM, mosdepth ran with `-n --by CIMG_04613.bed` twice:
1. Default filters. Duplicates and secondary reads are excluded. No MAPQ filter. Column: `mean_depth`.
2. With `-Q 20` (MAPQ >= 20). Column: `mean_depth_Q20`.

`genome_mean` is the genome-wide mean depth from the existing 5000 bp-window summary. It has no MAPQ filter. `ratio` = `mean_depth` / `genome_mean`. I did not compute a Q20 ratio, because the baseline has no MAPQ filter.

Test: the pipeline ran on 1M0 and 21SD before the full run. The full run produced 559 rows.

## Results
| stat | value |
|---|---|
| strains | 559 |
| mean of per-strain gene depth | 88.1 |
| min / max gene depth | 0.00 / 1036.04 |
| strains with zero depth | 2 |
| ratio < 0.5 | 11 strains |
| ratio > 1.5 | 6 strains |

### Bottom 11 by ratio
| strain | mean_depth | mean_depth_Q20 | genome_mean | ratio |
|---|---|---|---|---|
| CA12 | 0.00 | 0.00 | 0.27 | 0.000 |
| NM_7898 | 0.00 | 0.00 | 2.62 | 0.000 |
| NM_4297 | 0.39 | 0.34 | 10.54 | 0.037 |
| TX13 | 0.18 | 0.13 | 2.76 | 0.065 |
| CA10 | 0.12 | 0.05 | 1.34 | 0.090 |
| NM_3894 | 1.94 | 1.38 | 21.52 | 0.090 |
| NM_9861 | 6.10 | 4.45 | 48.18 | 0.127 |
| NM_0317 | 4.35 | 3.43 | 20.84 | 0.209 |
| NM_3957 | 3.23 | 2.37 | 15.33 | 0.211 |
| NM_459 | 17.69 | 10.04 | 64.24 | 0.275 |
| COCPO_103717 | 6.95 | 4.23 | 13.92 | 0.499 |

### Top 8 by ratio
| strain | mean_depth | mean_depth_Q20 | genome_mean | ratio |
|---|---|---|---|---|
| TX8 | 61.69 | 39.16 | 37.47 | 1.646 |
| GT-161 | 71.86 | 44.50 | 45.30 | 1.586 |
| SOIL_407-0_L_OLD_CPA0001 | 47.02 | 28.81 | 30.21 | 1.556 |
| CA14 | 62.04 | 47.12 | 40.08 | 1.548 |
| TX9 | 88.45 | 57.61 | 57.22 | 1.546 |
| UTAH_23742X189 | 328.06 | 209.12 | 218.54 | 1.501 |
| COCPO_175712 | 111.59 | 63.62 | 76.63 | 1.456 |
| COCPO_175714 | 106.91 | 66.91 | 74.55 | 1.434 |

### Highest gene depth
CA15 1036.04 (ratio 1.185), GT-102 424.11 (1.138), UTAH_23742X205 414.59 (1.117), COCPO_175716 336.44 (1.349), UTAH_23742X189 328.06 (1.501). These are deeply sequenced samples. They are not gene amplifications.

## Interpretation and limits
- CA12, CA10, NM_7898 and TX13 have genome means of 0.27 to 2.76. Their low gene depth reflects too few reads. You cannot call a deletion from them.
- NM_4297, NM_3894, NM_9861, NM_0317, NM_3957 and NM_459 have genome means of 15 to 64 and gene depth of 0.4 to 18. This could be a deletion or reads from a divergent allele that do not map to the reference. I did not test either explanation.
- In TX8, GT-161, CA14 and TX9, Q20 depth is much lower than unfiltered depth. This suggests multi-mapped reads in the gene.
- A ratio of 1.5 or more is consistent with duplication, but I did not test for it. Coverage in flanking windows is not examined here.
