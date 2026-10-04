# Audit of the SOWgp depth report (2026-10-01)

**Working report, 2026-10-04.** Scripts 50, 51, 60 to 67 in `analysis/cocci_repeats/`.
Audits `REPORT_2026-10-01_sowgp_depth_vs_repeats.md`.

> **Status: computational results only. No experimental validation.** All depth is measured
> against the *C. immitis* RS reference. Numbers come from the logs of the scripts listed in
> section 8. Logs are git-ignored (`*.log`). Rerun the scripts to regenerate them.

## 1. Questions

1. Do the numbers in the 2026-10-01 report reproduce from the raw files?
2. Do its conclusions hold after adjusting for read length, sequencing center and relatedness?
3. Does the depth method recover a known repeat-unit count? (A read simulation answers this.)
4. What do the short-read assemblies contain at the SOWgp locus?

## 2. Summary

| # | Claim in the 2026-10-01 report | Audit result |
|---|---|---|
| 1 | The 39 strains with no SOWgp gene model have normal depth. They are not deletions. | **Holds.** Ratios are 0.92 to 1.43 (lowest CA23, 0.92). In the assemblies, 37 of the 39 have an N gap in the array. |
| 2 | These strains have a higher depth ratio than full-length strains (1.19 against 1.04, p = 8e-8). The report reads this as repeat collapse. | **Explained by read length.** After adjustment the fold difference is 1.024 (95% CI 0.974 to 1.076, p = 0.36). The repeat-collapse reading is not needed. |
| 3 | Seven strains are deletion candidates. | **Not tested at read level.** The set is stable for ratio below 0.3 to 0.5 at 10x (6 to 7 strains). Two stay at 30x. |
| 4 | Array depth does not scale with unit count as predicted (slope 0.25 per unit). | **Holds, and the pipeline is not the cause.** The simulation recovers a slope of 0.245. Real slopes are 0.110 (*C. immitis*) and 0.016 (*C. posadasii*). |
| 5 | Pooled correlation of Q20 fraction with units is negative. | **Species effect.** Within species, the partial correlation is +0.123 (p = 0.08). |
| 6 | No multiple-testing correction. | **Done.** 23 tests ran, 14 are in the report. Section 7. |

## 3. Reproduction (script 60, section R0 to R1)

- The coverage and array tables rebuilt from the raw mosdepth files match the committed tables.
  Largest difference: 0.0005 (ratio and array/flank).
- Tip matching reproduces: 484 of 488 tips matched, 75 depth strains with no tip.
- Status counts reproduce: full-length 202, fragment only 243, no gene model 39, not in tree 75
  (559 strains).
- Script 40 hard-codes the flank lengths 249 and 165. These equal the BED lengths.

## 4. Strains with no gene model and the deletion candidates

### 4.1 Read length explains the higher ratio (script 61, Q1)

Read-length class is strongly linked to the ratio in full-length strains
(Kruskal-Wallis p = 3.4e-8). Median ratio falls from 1.130 (reads up to 101 nt) to 0.925 (301 nt).
The no-model strains are mostly short-read strains: 24 of 38 have reads up to 101 nt and 13 have
126 nt (strains with genome depth of at least 10x).

| Model (n = 482) | No-model against full-length, fold in ratio | 95% CI | p |
|---|---|---|---|
| Status only | 1.147 | 1.086 to 1.210 | 8.4e-7 |
| + read-length class | 1.022 | 0.971 to 1.076 | 0.40 |
| + sequencing center | 1.109 | 1.051 to 1.170 | 1.8e-4 |
| + log insert size | 1.138 | 1.081 to 1.198 | 1.2e-6 |
| + all covariates | 1.024 | 0.974 to 1.076 | 0.36 |
| Phylogenetic GLS (CDS tree) + all covariates | 1.009 | 0.961 to 1.059 | 0.72 |

Read length alone removes the effect. Center and insert size do not. Log genome depth has no
effect (coefficient -0.000, p = 0.97). Sequencing center is unrecorded for 43 full-length, 76
fragment-only and 8 no-model strains, so the center adjustment is partial.

### 4.2 The seven deletion candidates (script 60, R3)

Candidates at ratio below 0.5 and genome depth of at least 10x: COCPO_103717, NM_0317, NM_3894,
NM_3957, NM_4297, NM_459, NM_9861.

| Ratio below | Strains at genome depth of at least 10x | at least 30x |
|---|---|---|
| 0.3 | 6 | 2 (NM_459, NM_9861) |
| 0.5 | 7 | 2 |
| 0.7 | 13 | 6 |

`audit_deletion_sensitivity.tsv` has the full grid. A divergent allele that does not map to RS
gives the same signal as a deletion. The reads at the locus have not been inspected.

### 4.3 CA10 and CA12 (unresolved)

- In the RS CRAMs, CA12 has genome depth 0.27x and CA10 has 1.34x.
- Their assemblies are normal in size (32.8 Mb) and `asm_stats` lists 17.3x and 134.6x coverage.
- Each assembly aligns to its own species' long-read reference at 0.07 to 0.08% (about 22 kb of
  28 Mb). I aligned each assembly to the reference of the other species as well. Neither gave any
  alignment block of 500 bp or more (`minimap2 -x asm20`, one-off command, not in a script).
- So these are not a species-label swap. The CRAM reads and the assembly probably do not come
  from the same sample, or the assemblies are not *Coccidioides* sequence. I did not find the cause.
  The strains have no tree tip, so they are in no depth model.

## 5. Depth against repeat-unit count

### 5.1 The simulation (scripts 65, 66)

Paired reads from alleles with a known unit count, through bwa-mem2, fixmate and sort, Picard
MarkDuplicates and mosdepth. These are the tool versions in the CRAM headers. Grid: read length
100 to 300, coverage 30x and 75x, error 0.2%, 336 runs. Allele sources: RS windows with 197 nt
genomic periods added or removed (each edit checked by translation), and four real long-read
alleles.

| Subset | Slope of array/flank per unit | 95% CI | Prediction |
|---|---|---|---|
| RS edits, u = 2 to 5 (n = 192), all reads | 0.245 | 0.232 to 0.257 | 0.25 |
| same, Q20 reads only | 0.207 | 0.186 to 0.229 | 0.25 |
| RS edits, u = 2 to 6 (n = 240), all reads | 0.217 | 0.206 to 0.227 | 0.25 |

- Read length changes the Q20 slope. At 100 nt it is 0.126. At 300 nt it is 0.201 to 0.215.
- Real long-read alleles, median array/flank against u/4: CiB10637 (u = 5) 1.185 against 1.25.
  Silveira 2022 (u = 4) 0.955 against 1.00. VFC140 (u = 2) 0.515 against 0.50. Cpos1038 (u = 5)
  1.214 against 1.25.
- The Q20 fraction of array reads is 0.92 to 0.94 for the RS-like and VFC140 alleles. It is 0.47
  for Cpos1038 and 0.70 for Silveira 2022. Divergence from RS lowers Q20 mapping in the
  *C. posadasii* alleles. Q20 analyses are biased for *C. posadasii*.
- Limits of the simulation: uniform coverage, no GC or PCR bias, and one substitution model.

### 5.2 The real data (scripts 60, 61)

Full-length strains with genome depth of at least 10x. Slope of array/flank per gene-model unit.
The prediction is 0.25.

| Species, model | Slope | 95% CI | p (slope = 0.25) |
|---|---|---|---|
| *C. immitis* (n = 73), OLS | 0.110 | 0.039 to 0.182 | 2.2e-4 |
| *C. immitis*, OLS + covariates | 0.076 | 0.012 to 0.140 | 9.0e-7 |
| *C. immitis*, phylogenetic GLS + covariates (lambda 0.98) | 0.077 | 0.017 to 0.136 | 2.1e-7 |
| *C. posadasii* (n = 129), OLS | 0.016 | -0.032 to 0.064 | 7.3e-17 |
| *C. posadasii*, OLS + covariates | 0.066 | 0.018 to 0.115 | 9.7e-12 |
| *C. posadasii*, phylogenetic GLS + covariates (lambda 0.00) | 0.064 | 0.017 to 0.111 | 2.3e-12 |

Covariates: log genome depth, sequencing center, read-length class, log insert size.
Every slope is far below 0.25. The *C. posadasii* slope is near zero in the raw data and about
0.065 after adjustment. Depth of *C. immitis* shows strong phylogenetic signal (lambda 0.98).
*C. posadasii* shows none (lambda 0.00).

Array/flank also depends on read-length class (Kruskal-Wallis p = 5.7e-5 for *C. immitis*,
p = 0.004 for *C. posadasii*). The simulation shows little read-length effect, so this is not a
mapping effect of the pipeline.

### 5.3 What this means

- The simulation recovers the predicted slope. So the model "array/flank = u/4" is not the
  reason for the low real slopes. Explanation 2 of section 5 of the 2026-10-01 report is not
  supported.
- Two explanations stay open. (a) The gene-model unit counts are wrong in many strains.
  (b) Depth carries a bias that the covariates do not capture. I did not test either.
- Evidence on (a). In contiguous, N-free arrays of the annotated scaffolds (script 63), the unit
  count from the assembly equals the gene-model count in 164 of 176 full-length strains. In the 12
  others, the assembly gives 1 more unit (11 strains) or 2 more (1 strain). This covers only the 87% of full-length strains whose array is
  contiguous. Strains with a gap in the array (55% of fragment-only strains) are not covered.
- Evidence from depth. Depth-implied units (script 67, calibrated on the simulation, residual SD
  0.41 units) exceed the gene-model count by at least 1 in 42% of strains in both species. This
  uses the same depth data, so it does not independently test (a).
- In strains with a contiguous, N-free array, the slope is 0.114 (95% CI 0.036 to 0.192) for
  *C. immitis* and 0.056 (0.006 to 0.107) for *C. posadasii*. Restricting to clean arrays does
  not bring the slope to 0.25.

### 5.4 Species and Q20 fraction (script 60, R6)

The pooled Spearman correlation of Q20 fraction in the array with units is -0.168 (p = 0.017).
Median Q20 fraction in the array is 0.896 in *C. immitis* and 0.465 in *C. posadasii*. Units
differ by species (mean 3.34 and 3.79). The partial Spearman correlation given species is +0.123
(95% CI -0.030 to 0.264, p = 0.08). A within-species permutation test gives p = 0.094. The negative
pooled value is a species effect.

## 6. Assemblies

### 6.1 The SOWgp locus in the assemblies (scripts 62, 63)

Query: RS GG704914:968,094-973,690 (CIMG_04613 plus 2 kb on each side). Classes use minimap2
and blastn.

| Pangenome status | Strains | Contiguous array | Gap in array | Absent |
|---|---|---|---|---|
| Full-length | 203 | 176 (87%) | 27 | 0 |
| Fragment only | 243 | 110 (45%) | 133 | 0 |
| No gene model | 39 | 0 | 37 | 2 |

- In the 36 no-model strains with N bases in the array, the N count is 90 to 101 (median 99).
- The two "absent" no-model strains have the locus in the contigs that the AAFTF mitochondrial
  screen removed. The mitochondrial bin is large in many strains (median 146 kb, 206 strains
  over 200 kb, maximum 1.87 Mb).
- The share of arrays with a gap depends on read length: 97 of 132 strains with reads up to
  101 nt, against 2 of 31 with 301 nt.

### 6.2 Assemblies against long-read references (scripts 50, 51)

Per 10 kb bin of the reference, the fraction covered by the assembly. 486 strains.

| Species | Reference | Strains | Median fraction covered | Median missing |
|---|---|---|---|---|
| *C. immitis* | CiB10637 | 150 | 0.956 | 1.24 Mb |
| *C. posadasii* | Silveira 2022 | 336 | 0.939 | 1.73 Mb |

These values mean "absent from the assembly". They do not mean "deleted in the strain". The
medians include all strains. Only CA10 and CA12 (section 4.3) are below 0.5.
The per-bin table (7.7 MB) is not committed. Rerun script 50 to regenerate it.

## 7. Multiple testing (script 60, R8)

Scripts 40 to 42 ran 23 hypothesis tests. The 2026-10-01 report shows 14. The table
`audit_tests.tsv` has raw p, Benjamini-Hochberg q and Holm p within test families.

- Survive correction: the three status-against-ratio tests (q at most 0.0006), and the
  *C. immitis* array/flank correlation (rho +0.37, q 0.005).
- Do not survive at q below 0.05 within the family: pooled Q20 fraction against units
  (q 0.062), and the *C. posadasii* Q20 fraction (q 0.127).
- The status tests are confounded by read length (section 4.1). The correction does not
  address that.

## 8. Limits

1. The seven deletion candidates have no read-level check.
2. Causes (a) and (b) of the low slope are open (section 5.3).
3. CA10 and CA12 are unresolved (section 4.3).
4. Sequencing center is unrecorded for 127 strains (full-length 43, fragment 76, no-model 8).
5. The phylogenetic models assume Brownian motion within each species. The species split is a
   fixed effect. The CDS tree and the SNP tree give different lambda values in *C. posadasii*
   (0.00 against 0.95), and the SNP tree has fewer strains (88).
6. The simulation has uniform coverage. It does not model GC or PCR bias.
7. The scripts have no automated tests. R0 of script 60 is the only check against committed output.
8. The 2026-10-01 report is not edited. Its sections 4.1 (repeat-collapse reading), 5 (explanation 2)
   and 6 (limits 1 and 7) are superseded by this page.

## 9. Files

| File | Contents |
|---|---|
| `50_asm_vs_ref.sh`, `51_paf_bins.py` | assembly against long-read reference, per 10 kb bin (SLURM) |
| `60_depth_audit_recompute.py` | recompute every number of the 2026-10-01 report; tests with CIs |
| `61_depth_adjusted_models.R` | covariate and phylogenetic models (ape, phytools, sandwich) |
| `62_asm_locus_check.sh`, `63_asm_locus_parse.py` | SOWgp locus in each assembly (SLURM, parse) |
| `64_cram_readlen.sh` | read length and insert size per strain from the CRAMs |
| `65_sim_array_depth.sh`, `65_sim_array_reads.py`, `66_sim_array_summary.py` | read simulation and summary |
| `67_unit_depth_calibration.py` | depth-implied unit count, calibrated on the simulation |
| `audit_tests.tsv`, `audit_deletion_sensitivity.tsv`, `audit_adjusted_models.tsv` | audit tables |
| `asm_locus_summary.tsv`, `asm_locus_query.fa` | per-strain locus class; query sequence |
| `asm_vs_ref_strains.tsv`, `asm_vs_ref_strain_summary.tsv` | strain list; per-strain coverage |
| `cram_readlen.tsv`, `sim_alleles.tsv`, `sim_array_depth.tsv`, `sim_array_depth.regions.tsv.gz` | read length; simulation inputs and output |
| `plots/sim_array_depth.png` | simulation figure |
| Not committed | `asm_vs_ref_bins.tsv.gz` (7.7 MB), `asm_locus_raw.tar.gz` (raw PAF and BLAST output), logs |

Run order: 64, 65 and 66, 62 and 63, 50, then 60, 61 and 67. Run Python scripts with
`/usr/bin/python3.12`. SLURM scripts use `$SLURM_SUBMIT_DIR` and `$SCRATCH`.
