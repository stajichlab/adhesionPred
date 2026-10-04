# SOWgp read depth at CIMG_04613: absence and repeat length

**Working report, 2026-10-01.** Scripts 40 and 41 in `analysis/cocci_repeats/`.
Follows `REPORT_2026-09-29_sowgp_repeat_structure.md` and
`docs/reports/2026-09-27-coccidioides-antigen-findings.md` (section 4.4).

> **Audit 2026-10-04:** `REPORT_2026-10-04_sowgp_depth_audit.md` supersedes sections 4.1 (repeat-collapse
> reading), 5 (explanation 2) and 6 (limits 1 and 7) of this page.

> **Status: computational results only. No experimental validation.** Depth is measured against
> the *C. immitis* RS reference. It does not show gene content in strains whose sequence differs
> from RS.

---

## 1. Questions

1. Do the strains with no SOWgp gene model in the pangenome (the "8% absent") have low read depth
   at the RS locus, as a deletion would?
2. Is read depth related to SOWgp repeat allele length (number of 47 aa units)?

## 2. Data and tools

| Item | Detail |
|---|---|
| Alignments | 559 CRAMs: `/bigdata/stajichlab/shared/projects/Population_Genomics/Coccidioides/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/aln/*.cram`, BWA alignments to the RS reference |
| Reference | `.../all_C_immitis_ref_RS/genome/FungiDB-68_CimmitisRS_Genome.fasta` |
| Gene model | `.../genome/FungiDB-68_CimmitisRS.gff`, `CIMG_04613-t26_1`, `GG704914:970,094..971,690`, + strand, 7 exons |
| Depth tool | mosdepth 0.3.12 (`module load mosdepth`), `-n --by <BED> -t 2`; second run adds `-Q 20` |
| Whole-gene depth | `pipeline/11b_mosdepth_gene.sh`, SLURM job 29335355. Table: `CIMG_04613.coverage.tsv`. Details: `CIMG_04613_coverage_summary.md` |
| Array depth | `pipeline/11c_mosdepth_regions.sh`, SLURM job 29336045 (1 min 32 s). BED: `.../coverage/mosdepth_gene/CIMG_04613.array_regions.bed`. Output: `.../coverage/mosdepth_regions/<strain>.CIMG_04613_array.{all,Q20}.regions.bed.gz` |
| Genome-wide depth | `.../coverage/mosdepth/<strain>.5000bp.mosdepth.summary.txt`, `total` row (from `pipeline/11_mosdepth.sh`) |
| Allele data | `sowgp_tree_units.tsv` (script 13): status, unit count, clade species |
| Statistics | Python 3.12, SciPy (Spearman, Mann-Whitney, linear regression) |

Ratio = gene (or region) depth / genome-wide mean depth of the same strain. Genome-wide depth
has no MAPQ filter, so no Q20 ratio to genome depth is reported.

Strain names: the tree tips carry a `Coccidioides_<species>_` prefix. I removed it to join. 484 of
488 tips matched a depth strain. The four unmatched tips are the references `CimmitisRS`,
`CimmitisWA211`, `CposadasiiSilveira2022` and `Coccidioides_posadasii_Coahuilla_2`. 75 depth
strains have no tree tip.

## 3. Gene structure at the RS locus

From the GFF and the translated CDS (975 nt, 324 aa):

- The anchor `PTDCYGDC` occurs at aa 84, 130, 177 and 224. This is 4 units.
- The CDS has 7 exons. Exons 4, 5 and 6 are 141 nt each (47 aa), separated by 56 nt introns. The
  introns are not at the unit boundaries. Each unit has one intron, between its aa 31 and aa 32,
  so each 141 nt exon holds the end of one unit (aa 32-47) and the start of the next (aa 1-31).
  This gives a 197 nt distance between repeat starts in the genome. It fits the 188-197 nt spacing
  seen by tblastn in section 5 of `REPORT_2026-09-29_sowgp_repeat_structure.md`. The same 194-197
  nt spacing and intron are present in all five long-read assemblies
  (`REPORT_2026-10-01_sowgp_longread_annotation.md`).
- Repeat array (aa 84-270): `GG704914:970,511-971,295`. Spliced CDS in the array: 561 nt.
- Flank (control, within the same gene): CDS before the array (249 nt) and after it (165 nt).
  Total 414 nt.

## 4. Results

### 4.1 Strains with no gene model are not low-depth strains

Script 41. Strain status comes from the pangenome (section 4.6 of the 2026-09-29 report).

| Pangenome SOWgp status | Strains with depth | Median ratio | Ratio < 0.5 | Ratio < 0.5 and genome depth >= 10x | Median genome depth |
|---|---|---|---|---|---|
| full-length copy | 202 | 1.04 | 0 | 0 | 74.0 |
| fragment only | 243 | 1.11 | 1 | 1 | 63.0 |
| **no SOWgp gene model** | **39** | **1.19** | **0** | **0** | 46.7 |
| not in tree | 75 | 1.04 | 10 | 6 | 162.4 |

- All 39 strains with no gene model have ratios of 0.92 to 1.43. Reads cover the gene in each of
  them. Table: `sowgp_missing_model_depth.tsv`.
- The lowest of the 39 is CA23 (ratio 0.92). A deletion would give a ratio near 0.
- Strains with no gene model have a higher median ratio than full-length strains (1.19 against
  1.04, Mann-Whitney p = 8e-8). Fragment-only strains are also higher than full-length ones
  (1.11, p = 4e-4). Extra depth fits repeat collapse in the assembly, or reads from more than one
  locus mapping here. This was not tested.
- The 7 strains with ratio < 0.5 and genome depth >= 10x are `NM_4297`, `NM_3894`, `NM_9861`,
  `NM_0317`, `NM_3957`, `NM_459` and `COCPO_103717` (ratio 0.04 to 0.50). Six have no tree tip.
  `NM_459` is fragment-only. None is among the 39 strains with no gene model. These 7 are the
  deletion candidates. Low depth could also come from a divergent allele that does not map to RS.
  This was not tested.
- Four other low-ratio strains (`CA12`, `CA10`, `NM_7898`, `TX13`) have genome depth below 3x.
  Their low gene depth reflects too few reads.

**Conclusion.** The pangenome's 8% "absent" (39 of 488) is a gene-model gap. Read depth shows the
locus is present in all 39. The deletion candidates are a different set of strains and are
mostly not in the pangenome.

### 4.2 Whole-gene depth against unit count

Done inline in an earlier session. It is not yet in a script. 202 full-length strains.

| Subset | Spearman rho (ratio vs units) | p |
|---|---|---|
| both species | +0.18 | 0.01 |
| *immitis* (n = 73) | +0.30 | 0.011 |
| *posadasii* (n = 129) | +0.04 | 0.69 |

The gene is 1,597 bp. One unit is 141 bp. A one-unit change has little effect on whole-gene
depth. Section 4.3 uses the array only.

### 4.3 Array depth against unit count

**Prediction, written before the data were seen.** RS has 4 units. If reads from a strain with
*u* units map inside the array, then array depth / flank depth is about *u*/4. The slope per
unit is then 0.25.

Script 40. Strains: full-length copy and genome depth >= 10x (n = 202). Table:
`sowgp_array_depth.tsv`. Figure: `plots/sowgp_array_depth.png`.

| Species | Units | n | Median array/flank | IQR | Predicted |
|---|---|---|---|---|---|
| *immitis* | 3 | 48 | 1.00 | 0.93-1.10 | 0.75 |
| *immitis* | 4 | 22 | 1.11 | 1.01-1.31 | 1.00 |
| *immitis* | 5 | 2 | 1.14 | 1.11-1.17 | 1.25 |
| *posadasii* | 3 | 45 | 1.11 | 1.00-1.28 | 0.75 |
| *posadasii* | 4 | 60 | 1.13 | 1.01-1.30 | 1.00 |
| *posadasii* | 5 | 22 | 1.11 | 1.05-1.37 | 1.25 |

(Strains with 2 units: 1 *immitis*, 2 *posadasii*. Not shown.)

| Species | Spearman rho | p | Slope per unit (95% CI) | Predicted slope |
|---|---|---|---|---|
| *immitis* | +0.37 | 0.0012 | 0.110 (0.040 to 0.181) | 0.250 |
| *posadasii* | +0.07 | 0.46 | 0.016 (-0.032 to 0.063) | 0.250 |

- **The prediction does not hold.** The observed slope is below 0.25 in both species. The
  *immitis* confidence interval excludes 0.25. The *posadasii* interval is close to 0.
- *immitis* shows a small positive trend in the predicted direction. *posadasii* shows none.
- Array depth is about equal to flank depth at every unit count. A 3-unit strain has about the
  same array depth as a 4-unit strain.

Q20 results (MAPQ >= 20), same strains:

| Measure | *immitis* rho (p) | *posadasii* rho (p) |
|---|---|---|
| array/flank, Q20 | +0.25 (0.033) | +0.29 (0.001) |
| fraction of array reads kept at Q20 | -0.06 (0.63) | +0.17 (0.052) |

Across both species, the fraction of reads kept at Q20 falls with unit count (rho = -0.17,
p = 0.017). Species explains this. *posadasii* strains keep a lower fraction at Q20 in the flank
too. Within each species the effect is not significant.

I do not have an explanation for the positive Q20 trend. More units would be expected to give
more multi-mapped reads, which Q20 removes.

## 5. What the unit-count result does and does not show

Three explanations fit "array depth does not scale with unit count". I have not tested any of them.

1. The gene-model unit counts in short-read assemblies are wrong for many strains (collapsed
   arrays), so unit count is not a good measure of true length.
2. The *u*/4 model is too simple. Reads from near-identical units may not spread evenly over the
   4 RS units. Reads at the array edge map uniquely to the flanks. Unit classes differ at a few
   sites. A simulation with known unit counts would test the model.
3. Reads from other loci with similar sequence map into the array.

A test of explanation 1 uses the 5 long-read strains (UArizona), where unit count is known
from the assembly. Done in `REPORT_2026-10-01_sowgp_longread_annotation.md`: 3 of 6 loci had split
gene models, which supports explanation 1. Those strains have no short-read depth, so the depth
model itself is still untested. A test of explanation 2 would simulate reads from RS-like arrays with 3, 4 and
5 units and map them to RS.

## 6. Limits

1. The strains are not independent. Strains from one population are related. The p-values are
   too small for this reason. The 2026-09-29 report found no phylogenetic signal in unit count.
2. All depth is against RS. For *posadasii* strains, flank sequence differs from RS (base length
   281 against 277 aa). Mapping rates to RS may differ for this reason. This was not tested.
3. The gene-model unit counts come from mostly short-read assemblies (see the 2026-09-29 report,
   section 6).
4. Only strains with a full-length copy and genome depth >= 10x are in section 4.3. This removes
   the strains most likely to have collapsed or fragmented arrays.
5. The array region is 561 nt of spliced CDS. Reads that span introns are not counted over the
   intron bases, and the BED excludes introns.
6. Section 4.2 was not run from a script.
7. No multiple-testing correction was applied.

## 7. Next steps

1. Test unit-count accuracy on the 5 long-read strains (assembly unit count against short-read
   array depth if short reads exist for them).
2. Simulate reads from arrays with 3-5 units and check the depth model.
3. Join the 7 deletion candidates to a read-level check: look at the split reads and discordant
   pairs at the locus, or assemble the reads that map to the locus.
4. Move section 4.2 into script 41.
5. Update `docs/reports/2026-09-27-coccidioides-antigen-findings.md` (done for the 8% figure; see
   below) and section 4.2 of `docs/reports/2026-09-27-cocci-repeat-surface-proteins.md`.

## 8. Files

| File | Contents |
|---|---|
| `40_sowgp_array_depth.py` | array depth vs unit count; writes `sowgp_array_depth.tsv`, `plots/sowgp_array_depth.png` |
| `41_sowgp_status_depth.py` | status vs depth join; writes `sowgp_status_depth.tsv`, `sowgp_missing_model_depth.tsv` |
| `CIMG_04613.coverage.tsv` | whole-gene depth, 559 strains |
| `CIMG_04613_coverage_summary.md` | whole-gene depth method and outliers |
| `.../Genotyping/all_C_immitis_ref_RS/pipeline/11b_mosdepth_gene.sh`, `11c_mosdepth_regions.sh` | mosdepth SLURM scripts |
| `.../Genotyping/all_C_immitis_ref_RS/scripts/gene_coverage_table.py` | whole-gene table builder |

Run: `sbatch pipeline/11c_mosdepth_regions.sh <BED> <TAG>` in the Genotyping directory, then
`/usr/bin/python3.12 40_sowgp_array_depth.py` and `41_sowgp_status_depth.py` here. The scripts
run in seconds.

**Not experimentally validated.** Every statement here is from sequence alignment and read depth.
