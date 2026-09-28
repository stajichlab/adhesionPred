# Coccidioides antigen ranking — rebuilt 2026-09-27

Rewrite of `data/curated/antigens/coccidioides_candidates.tsv` after an independent review
found four defects. Scored on the **488-proteome Coccidioides pangenome** on HPCC
(`shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome`), not on two reference
proteomes.

## What was wrong with v1, and what changed

| defect | fix |
|---|---|
| **SOWgp and Ag2/PRA absent from the ranking.** The *C. posadasii* reference used was strain **C735 ΔSOWgp — a SOWgp deletion strain** | scored on the pangenome; all 9 anchors map and are present |
| **Uncalibrated**: CF antigen scored 4/8 tied with 678 others (63% of rows) | continuous score; explicit acceptance test |
| **No fungal cross-reactivity term** — only human. But cross-reaction with *Histoplasma*/*Blastomyces*/*Paracoccidioides*/*Aspergillus* is the real clinical failure mode of cocci serology | added, and it is now the dominant negative term |
| 63% of rows tied | continuous; plus one entry per orthogroup, since OG0000001's 14 members had flooded 14 of the top 50 |

Also removed: copy-number variability from the score. It is interesting biology and it answers
the separate pangenome question, but for a serodiagnostic a variable multi-copy family is a
liability. Now reported, not scored.

## The score

Two axes, kept separate because v1 wrongly merged them:
- **antigenicity** = 2.5·prevalence + 2.0·(homology to curated IEDB fungal antigens)
- **specificity** = −3.0·(max % identity to a confounder fungus) − 0.5·(breadth of confounders) − 1.5·(human homology)

## Acceptance test — the ranking is only partly calibrated

| anchor | combined | antigenicity | specificity | max cross-react | verdict |
|---|---|---|---|---|---|
| PRA3 | 6.1% | 75.0% | **6.0%** | **0%** | PASS |
| Ag2/PRA | 7.2% | **0.0%** | 50.0% | 60% | PASS |
| **SOWgp** (all 3 alleles) | 7.3% | 79.7% | **7.3%** | **0%** | PASS |
| PRA2 | 10.7% | 0.5% | 71.7% | 69% | FAIL (just outside) |
| CF antigen (CiX1) | 56.1% | 17.9% | 63.4% | 64% | PASS as a *negative* specificity control |

**3/4 Coccidioides-specific anchors in the top decile.** The script prints
"NOT CALIBRATED" whenever that is not 4/4, and that warning should be respected: this is a
ranking with 4 usable controls.

## The finding worth acting on

**SOWgp and PRA3 have no detectable homolog in *Histoplasma*, *Blastomyces*,
*Paracoccidioides* or *A. fumigatus*. Ag2/PRA and PRA2 do.**

Ag2/PRA and PRA2 both hit the *same* genes in all three dimorphic confounders
(*Histoplasma* FBAD1291_004826, *Blastomyces* F00FD2C2_005998, *Paracoccidioides*
FD6225C6_004143) at 55–70% identity with **full query coverage** — these are genuine
orthologs, not low-complexity artifacts. That is a concrete, testable prediction of
serological cross-reactivity for the Ag2/PRA family, and a reason to prefer **SOWgp and PRA3**
as species-specific markers.

## A real limitation of the antigenicity axis

SOWgp and PRA3 score *badly* on antigenicity (75th, 80th percentile) precisely because they
are Coccidioides-specific: the axis is homology to IEDB antigens, which are dominated by other
fungi, so **species-specific antigens are systematically penalised**. The specificity axis is
doing all the useful work. Fixing this needs direct evidence, not homology:
- **spherule-phase expression** (249 *Coccidioides* RNA-seq runs are in SRA) — SOWgp and
  Ag2/PRA matter *because* they are spherule-phase; this is the most discriminating filter
  available and costs only CPU
- B-cell epitope surface accessibility, scored separately from T-cell evidence

## Outputs

- `cocci_antigen_ranking.tsv` — 9,139 proteins / 8,542 orthogroup representatives, with
  prevalence across 488 proteomes, copy number, cross-reactivity and both score axes
- `shortlist_specific_universal.tsv` — **595 orthogroups** present in ≥95% of 488 proteomes
  with no detectable confounder or human homolog (557 of them single-copy)
- `anchors.m8` — anchor → reference gene mapping

## For goal 3 (pangenome variability), already visible here

- 1,580 orthogroups are **accessory** (<95% of 488 proteomes); 286 have copy-number CV > 0.3
- **SOWgp is at 92% prevalence — it is not universal**, which is itself a finding
- The three SOWgp alleles hit one locus (`CIMG_04613`) at 96 / 85 / 74% identity, i.e. the
  repeat-number variation is visible in the sequence data

Caveat carried forward: repeat-copy-number, SOWgp's defining feature, is collapsed by
short-read assemblies. The prevalence and CN numbers here come from an assembly-based
orthogroup table and should be confirmed from read depth before being relied on.

## Reproduce

```bash
srun -p epyc -c 8 --mem 16G -t 90 ./01_build_inputs.sh   # mmseqs maps + cross-reactivity
srun -p epyc -c 4 --mem 12G -t 30 python3 02_score_antigens.py \
    --work cocci_antigens --antigens /dev/null --iedb-hits cocci_antigens/iedb.m8 \
    --out cocci_antigens/cocci_antigen_ranking.tsv
```
