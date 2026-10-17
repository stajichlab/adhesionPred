# Cell wall remodeling and cell fusion family report (no call)

Date: 2026-10-17. Pattern: option c-i of the lectin report. Plan: docs/reports/cell-wall-fusion/04-repo-coverage-and-gaps.md section 4.

This report makes no call. It writes no status file. It has no call or status column. A hit is a hit to a Pfam family model. It is NOT evidence of adhesion, of cell fusion or of a cell wall remodeling function.

Cell wall enzymes are known false-positive sources for surface-glycoprotein detection. In the S288C check of docs/reports/cell-wall-fusion/07-model-hit-check.md the ESM classifier called 6 of 13 wall enzymes above 0.5 (CRH1, UTR2, SCW11, CTS1, SUN4, KRE9), and none of the 13 hit an active family-table row. Do not read a family hit here as an adhesin.

## Method

- Search: `hmmsearch --cut_ga` (HMMER 3.4 (Aug 2023); http://hmmer.org/) with the models of `wall_fusion_models.tsv` against Pfam 38.2 (`/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`). Set: all.
- The cutoff is the Pfam GA cutoff of each model. A hit is a protein with a sequence score at or above the sequence GA. Scores in the tables are sequence scores.
- Column `domain_ge_ga` shows whether the best domain also reaches the domain GA.
- Proteomes: 1. No signal peptide, TM or other condition is applied to the hits.
- `n_proteins` counts unique protein IDs (first word of the header). Records that repeat an ID are merged: one row per protein and model, with the best sequence score. Each repeat is reported as a warning below.

## Models (all are family-only; none is a call)

| set | model | name | family | verification | hit proteins | proteomes with hit | best score |
|---|---|---|---|---|---|---|---|
| wall_remodeling | PF03198 | Glyco_hydro_72 | GH72 | 05_interpro | 1 | 1 | 504.4 |
| wall_remodeling | PF00722 | Glyco_hydro_16 | GH16 | 05_interpro | 2 | 1 | 170.8 |
| wall_remodeling | PF00332 | Glyco_hydro_17 | GH17 | 05_interpro | 2 | 1 | 353.2 |
| wall_remodeling | PF03663 | Glyco_hydro_76 | GH76 | 05_interpro | 2 | 1 | 547.4 |
| wall_remodeling | PF00704 | Glyco_hydro_18 | GH18 | 05_interpro | 1 | 1 | 81.6 |
| wall_remodeling | PF00150 | Cellulase | GH5 | 05_interpro | 1 | 1 | 43.7 |
| wall_remodeling | PF01522 | Polysacc_deac_1 | chitin_deacetylase | 05_interpro | 0 | 0 | NA |
| wall_remodeling | PF29664 | BtgE-BtgC_glyc-hydro | BtgE-BtgC | 07_hmmscan_only | 3 | 1 | 283.2 |
| wall_synthesis | PF01644 | Chitin_synth_1 | chitin_synthase_1 | 05_interpro | 0 | 0 | NA |
| wall_synthesis | PF03142 | Chitin_synth_2 | chitin_synthase_2 | 05_interpro | 0 | 0 | NA |
| wall_synthesis | PF08407 | Chitin_synth_1N | chitin_synthase_1N | 05_interpro | 0 | 0 | NA |
| wall_synthesis | PF02364 | Glucan_synthase | glucan_synthase | 05_interpro | 0 | 0 | NA |
| wall_structural_other | PF10342 | Kre9_KNH | Kre9_Knh1 | 05_interpro | 1 | 1 | 46.7 |
| wall_structural_other | PF05390 | Kre9_KNH1_C | Kre9_Knh1_C | 07_hmmscan_only | 1 | 1 | 166.2 |
| wall_structural_other | PF03856 | SUN | SUN | 05_interpro | 1 | 1 | 330.1 |
| wall_structural_other | PF12454 | Ecm33 | Ecm33 | 05_interpro | 0 | 0 | NA |
| wall_structural_other | PF10287 | YJL171C_Tos1_C | Tos1_C | 05_interpro | 0 | 0 | NA |
| wall_structural_other | PF10290 | YJL171C_Tos1_N | Tos1_N | 05_interpro | 0 | 0 | NA |
| wall_structural_other | PF07983 | X8 | X8 | 05_interpro | 1 | 1 | 102.4 |
| cell_fusion_candidate_models | PF27758 | PRM1 | Prm1 | 07_hmmscan_only | 1 | 1 | 475.9 |
| cell_fusion_candidate_models | PF12351 | Fig1 | Fig1 | 07_hmmscan_only | 1 | 1 | 192.3 |

Supporting report and verification status of every model are in `wall_fusion_models.tsv`.

## Per proteome and per set

Distinct proteins with a hit to at least one model of the set. The sets are not added up: enzymes (wall_remodeling, wall_synthesis) and non-enzyme models (wall_structural_other, cell_fusion_candidate_models) are never in one total.

| proteome | proteins | wall_remodeling | wall_synthesis | wall_structural_other | cell_fusion_candidate_models |
|---|---|---|---|---|---|
| test_set | 23 | 10 | 0 | 3 | 2 |

## Notes

- Fusion proteins such as Prm1 and Fig1 are membrane proteins. They lie outside the signal-peptide population of step 1 of the sorting hat. A hit to the Prm1 or Fig1 model does not show a role in cell fusion. The Fus1 SH3 hit (SH3_1) is not a fusion family and is not scanned.
- Ecm33 (PF12454): S. cerevisiae Ecm33 itself has no hit to this model (07 row 22). A hit is not an Ecm33 call.
- wall_structural_other holds models that are not an enzyme family by themselves (Kre9/Knh1, Ecm33, Tos1, X8). SUN is kept here because its role is unclear: the InterPro name says beta-glucosidase, the reports show no enzyme role.
- Large families (GH16, GH17, GH18, GH5) have intracellular and non-wall members. A hit says only that the protein carries the domain.
- Inclusion rule: a model is listed if it is verified in docs/reports/cell-wall-fusion/05-pfam-id-verification.md (05_interpro) or measured in 07-model-hit-check.md (07_hmmscan_only, name from Pfam-A, scores are hmmscan best-domain scores there). The script stops if a Pfam model NAME differs from the name in the model TSV.
- Models left out on purpose: PF00026 (generic aspartyl protease), SH3_1 (PF00018, the Fus1 hit, generic domain), the GPI-anchor biosynthesis models (ER enzymes), WSC (a sensor), PF10333 (Pga1, an ER GPI enzyme), and the chitin-binding models (PF00187, PF01607, LysM; in the lectin report).
- Scores: 07 lists the hmmscan best-domain score; this report stores the sequence score (`seq_score`) and the best-domain score (`domain_score`). For PRM1 the values are 475.9 and 475.5. Same HMMER 3.4 and Pfam 38.2.
