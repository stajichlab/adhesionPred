# SOWgp in the five UArizona long-read assemblies: unit counts and gene models

**Working report, 2026-10-01.** Scripts 43, 44 and 45 in `analysis/cocci_repeats/`.
Follows `REPORT_2026-10-01_sowgp_depth_vs_repeats.md` (section 5, explanation 1) and
`REPORT_2026-09-29_sowgp_repeat_structure.md` (section 5 and the 2026-09-30 correction).

> **Status: computational results only. No experimental validation.** No RNA evidence was used.
> No short reads for these five strains were found in the depth set (see section 6).

---

## 1. Question

Does the SOWgp unit count from the gene models match the unit count in the genome sequence?
The five long-read assemblies have the best chance of representing the repeat array. If the
gene models are wrong here, short-read gene models are probably worse.

## 2. Data and tools

| Item | Detail |
|---|---|
| Assemblies | `/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc/<strain>/*.scaffolds.fa` for CiB10637, CiB10992, VFC140, Cpos1038, Cpos3700 |
| Gene models | `.../For_Marc/<strain>/*.gff3` (funannotate) |
| Unit count in DNA (44) | tblastn, NCBI BLAST 2.13.0+ (`ncbi-blast/2.13.0+`), SEG off, E <= 1e-5. Query: one RS unit, CIMG_04613 aa 130-176 (`PTDCYGDCEDGYDYSPPPPPKKYGDCDYDDGYCDGPSKTSMKPEPPK`) |
| Independent gene model (43) | miniprot 0.2-r116 (`module load miniprot`), `--gff -G 1000`, queries from `sowgp_seed.fa` (8 proteins) |
| Comparison (45) | `45_longread_locus_summary.py` writes `longread_locus_summary.tsv`. Anchor counts come from `sowgp_anchored_longread.tsv` (origin not recorded; see README) |
| Outputs | `longread_unit_tblastn/<strain>.tsv`, `longread_miniprot/<strain>.sowgp.gff` |

**How units are counted in DNA.** Each unit gives two tblastn hits. One covers unit aa 1-31 and
one covers aa 32-47. A 56 nt intron lies between them. I count the hits that start at query
position 1 (one per unit). Spacing between units is 194-197 nt (141 nt of coding sequence plus
56 nt of intron).

## 3. Results

| Strain | Scaffold | Units in DNA | Spacing (nt) | funannotate models at the locus | Anchors in those models | miniprot full-length model |
|---|---|---|---|---|---|---|
| CiB10637 | scaffold_2 | **5** | 194, 197, 197, 197 | 1 (`CIB10637_003943`, 371 aa) | 5 | 100% identity to CiB10637_003943, aa 1-371 |
| CiB10992 | scaffold_2 | **5** | 197 x 4 | 1 (`CIB10992_003451`, 371 aa) | 5 | 100% identity, aa 1-371 |
| VFC140 | scaffold_2 | **2** | 197 | 1 (`VFC140_004201`, 234 aa) | 2 | not tested (no 2-unit query) |
| Cpos1038 | scaffold_2 | **5** | 197, 197, 197, 196 | **2** (`003233` + `003234`) | 2 + 3 = 5 | **98.7% identity to SOWgp66, aa 1-375** |
| Cpos3700 | scaffold_2 | **5** | 197, 197, 194, 196 | **2** (`003932` + `003933`) | 3 (in `003933`) | **100% identity to SOWgp66, aa 1-375** |
| Cpos3700 | scaffold_84 | **5** | 194, 197, 197, 197 | **2** (`010247` + `010248`) | 1 + 1 = 2 | **100% identity to CiB10637_003943, aa 1-371** |

Details:

- **CiB10637, CiB10992 and VFC140: gene model and genome agree.** Unit counts in the model equal
  those in the DNA (5, 5 and 2).
- **Cpos1038: one gene, called as two.** The DNA has 5 units. Two funannotate models
  (`CPOS1038_003233` at 2,242,608-2,243,292 and `CPOS1038_003234` at 2,243,373-2,244,147) are 80 bp
  apart. Together they hold the 5 anchors. miniprot with the 375 aa SOWgp66 protein gives one
  model over 2,242,608-2,244,147, 98.7% identical, with no frameshift flag.
- **Cpos3700 scaffold_2: one gene, called as two, and one unit group is missing.** The DNA has
  5 units. `CPOS3700_003932` (4,405,199-4,405,506) and `CPOS3700_003933` (4,405,569-4,406,734)
  are 63 bp apart. Only `003933` has anchors (3 of 5). miniprot with SOWgp66 gives one model over
  4,405,199-4,406,734, **100% identical to published SOWgp66** (375 aa). Both Cpos3700 loci
  share `003933`, which is why the earlier report listed it as truncated (aa 127-328).
- **Cpos3700 scaffold_84: a second SOWgp copy.** The DNA has 5 units. The two funannotate models
  (`CPOS3700_010247`, 37,069-37,543, and `CPOS3700_010248`, 37,624-38,597, 80 bp apart) hold
  2 anchors in total. miniprot with the CiB10637 protein gives one model, 37,069-38,597, **100%
  identical to the *C. immitis* CiB10637 allele** (371 aa). scaffold_84 is 48,103 bp long and has
  15 annotated genes.
- **VFC140_004201:** the miniprot alignment of the 8 seed proteins to the VFC140 locus has low
  identity (0.62 to 0.74) and frameshift flags. This is expected for a 2-unit target aligned to
  3-6 unit queries. It is not evidence of a problem in the VFC140 model.

**Totals.** Six loci in five strains. In 3 of the 6 loci (Cpos1038, both Cpos3700 loci) the gene
model is split and the unit count in the models (5, 3 and 2) does not match the genome (5, 5, 5)
unless the pieces are joined. In the other 3 (CiB10637, CiB10992, VFC140) the models match the
genome.

## 4. What this changes

1. **The 2026-09-29 report, section 5, is now resolved.**
   - Cpos1038 and Cpos3700 have full SOWgp genes. Each has 5 units at the DNA level.
   - The 188-197 nt repeat spacing there is the intron inside each unit. The earlier text said
     this "has not been checked". It is now checked in all five assemblies.
   - Cpos3700 has two SOWgp loci, not one. The earlier text said this was "not established".
     The two loci are on different scaffolds (scaffold_2 and scaffold_84), and each has 5 units.
2. **Split gene models are common at this locus.** Funannotate split 3 of 6 loci at an 80 bp (or
   63 bp) gap. The two ends were not joined. I did not test why. A short intron inside a
   repeated gene is one candidate.
3. **The unit counts in the pangenome are probably too low for strains with split models.** Of
   the 3 split loci, 2 were counted as short alleles in the model (anchors 2+3 or 3). This
   supports explanation 1 in section 5 of `REPORT_2026-10-01_sowgp_depth_vs_repeats.md`
   (gene-model unit counts are wrong for some strains). It is evidence from 5 long-read assemblies
   only. It does not show the size of the effect in the 205 short-read strains.
4. **Unit count and depth.** The long-read strains give no depth test, because these five strains
   have no CRAMs in the depth set (section 6).

## 5. Open question: Cpos3700 scaffold_84

The *C. posadasii* Cpos3700 assembly holds a SOWgp protein with 100% identity over 371 aa to the
*C. immitis* CiB10637 allele. I do not know why. Possible explanations, none tested:

- an *immitis*-derived sequence in a *posadasii* strain (introgression or hybrid ancestry),
- a contaminating or mis-assigned scaffold in the assembly,
- a real second locus that is conserved between the species.

A check would compare the scaffold_84 flanking genes with the CiB10637 and Cpos3700 scaffold_2
neighbourhoods, and compare the codon-level sequence of the two alleles.

## 6. Limits

1. Unit counts in DNA come from tblastn hits of one RS unit. Hits with a diverged unit may be
   missed. Cpos1038 has an extra weak hit (66.7% identity) at 2,243,501. I did not count it.
2. The miniprot models use protein alignment only. They have no RNA or splice-site training.
   The intron at each unit has a canonical size (56 nt) but I did not check donor and acceptor sites.
3. VFC140 has no clean query of matching length, so I did not build a miniprot model for it.
4. Five assemblies. No statement about frequency in the population.
5. **No depth data for these five strains.** I looked for them in `samples.csv`. Matches on the
   strain names `CiB10637`, `CiB10992`, `VFC140`, `Cpos1038` did not give a CRAM. One row mentions
   `NACVFR_Cocci_3700_S62` (sample `SA2`). I did not confirm that this is Cpos3700, so I did not
   use it.
6. No RNA-seq or protein evidence. Gene calls were not experimentally checked.

## 7. Next steps

1. Fix the three split models (a combined model per locus) and write them to a GFF for the
   pangenome.
2. Find out why funannotate split them (intron length filter, repeat masking, or gene finder
   boundaries).
3. Check the strain identity of `SA2` and the five long-read strains' short reads, if they exist,
   and test depth against the known unit counts. This is the clean test of the depth model.
4. Check scaffold_84 in Cpos3700 (section 5).
5. Re-count the pangenome: search the short-read assemblies for split models next to each other
   at 60-80 bp gaps with anchors on both pieces.

## 8. Files

| File | Contents |
|---|---|
| `43_longread_sowgp_miniprot.sh` | miniprot of the seed proteins to the five assemblies |
| `44_longread_unit_tblastn.sh` | tblastn of one RS unit to the five assemblies |
| `45_longread_locus_summary.py` | joins DNA unit count, funannotate models and miniprot models |
| `longread_locus_summary.tsv` | the table in section 3 |
| `longread_miniprot/`, `longread_unit_tblastn/` | raw outputs |

Run 43 and 44 on a compute node (`module load` is needed for miniprot and BLAST; 44 uses
`$SCRATCH`). Then `/usr/bin/python3.12 45_longread_locus_summary.py`.

**Not experimentally validated.** Every statement here is from sequence alignment.
