# Task 06. Controls for the specificity ranking (`cocci_specificity_rank_top15`, `serodiagnostic_marker_candidate`)

*Read `COMMON-RULES.md` first. Written 2026-10-06.*

## Goal

Find independent positive and negative controls for the *Coccidioides* genus-specificity ranking.
Check whether a leave-species-out test is possible (decision C3).

## Why this matters

- The ranking scores 9,139 *C. immitis* RS proteins. It combines **antigenicity** (prevalence in 488
  proteomes and homology to curated IEDB fungal antigens) and **specificity** (low identity to
  *Histoplasma*, *Blastomyces*, *Paracoccidioides*, *A. fumigatus* and human). It is not epitope
  prediction. No antibody or T-cell measurement supports it.
- The only truth today is four anchors (PRA3, Ag2/PRA, SOWgp, PRA2) and the CF antigen as a
  cross-reactive control. In the pre-set top-decile test, 3 of 4 anchors pass (PRA2 fails). The 15% cut
  was set after the anchors were seen. Status: `smoke`, `leakage=tuned_on_truth`.
- The anchors **tuned** the ranking. They cannot test it.
- The testable prediction in the report: sera from histoplasmosis, blastomycosis and
  paracoccidioidomycosis patients should react with Ag2/PRA and the CF antigen and **not** with SOWgp
  or PRA3. No such test exists.

## Known today

- `data/curated/antigens/antigens.tsv`: 86 rows from IEDB (classes `antigen_tcell`, `antigen_bcell`,
  `antigen_both`, `antigen_elution_only`; 5 rows with `taxon_role=focus`, 81 `comparison`). 6 rows are
  *Coccidioides* (4 T-cell, 2 both).
- `analysis/cocci_antigens/` (scripts 00 to 06 and the ranking), `docs/reports/2026-09-27-coccidioides-antigen-findings.md`.
- The ranking applies to *C. immitis* RS only (`--applicable-taxa 246410`).

## Two label axes (do not mix them)

| Axis | Positive | Negative |
|---|---|---|
| **Antigenic** (reacts in an assay) | Protein shown to react with patient or immunised-animal sera, or to give a T-cell response (cite assay, host, PMID) | Protein tested in the same assay and **not** reactive |
| **Specific** (does not cross-react) | Protein tested against sera from other mycoses and **not** reactive | Protein shown to react with sera from other mycoses (cross-reactive) |

A cross-reactive antigen is antigenic positive and specific negative. The ranking is meant to find
antigenic **and** specific proteins.

## Tasks

1. **Count first.** In IEDB and the literature, how many *Coccidioides* proteins have a recorded
   assay with a **negative** outcome? How many have cross-reactivity data? How many *Coccidioides*
   proteins have any assay at all? Report the numbers before you build anything else. If the
   negatives are fewer than 20 clusters, say that the ranking cannot get a specificity here.
2. **Feasibility of leave-species-out (C3).** For *Histoplasma*, *Blastomyces*, *Paracoccidioides*
   and *Coccidioides posadasii*: how many proteins with antigen assays exist, per species? Which
   species are in the confounder set of the ranking (they cannot be test species without changing the
   score)? List a confounder-free design, or say that none exists.
3. **Positive controls** (antigenic): published *Coccidioides* antigens beyond the four anchors (for
   example proteins named in serology or vaccine papers). Check each in the literature. Names to
   check, not to assume: the CF antigen (CTS1/CiX1), SOWgp, Ag2/PRA, PRA2, PRA3, and any others you
   find. Each needs PMID, assay, and a quote.
4. **Mark leakage.** SOWgp, PRA3, Ag2/PRA, PRA2 and the CF antigen are `tuned_on_truth`. Proteins
   with at least 30% identity to a curated IEDB antigen that the ranking uses are
   `in_reference` for the antigenicity axis.

## Size target and caution

Report clusters of antigenic positives, antigenic negatives, specific negatives and specific
positives, separately. The 20-cluster floor will probably not be met for specificity. A result of
"cannot be calibrated with the available data" is a valid result.

## Do not

- Do not call the ranking an epitope predictor or a diagnostic test.
- Do not use the spherule expression table as a label (a protein can be down in spherules and still
  antigenic; PRA3 is).
- Do not use tier lists (`TIER1_candidates.tsv`) as truth; they came from the ranking.
- Do not contact anyone or submit any request.
