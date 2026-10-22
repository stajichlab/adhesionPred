# Counts: wall, fusion and adhesion controls from Kiss 2019

*Written by `analysis/wall_fusion_truth/build_control_set.py` (DRAFT, Task 5.2). Every row is `label=positive` for the call named in `stratum` and `serves`. A wall row is never a negative and never a truth row of the adhesion call (D20). No number is pooled across slices.*

## Summary (`controls.tsv`: rows with a public UniProt sequence; checked by the L1 validator)

| Item | Value |
|---|---|
| Rows (proteins) | 112 |
| Sequences in `sequences.faa` | 112 |
| Independent clusters (Kiss cluster per slice; a proxy until Task 5.3) | 57 |
| Rows positive | 112 |
| Sequences positive | 112 |
| Clusters positive | 57 |
| Rows negative | 0 |
| Sequences negative | 0 |
| Clusters negative | 0 |

This set has no negatives on purpose: hard negatives come from the shared set (`data/controls/hard-negatives-shared`) and are chosen by the owner in Task 5.3.

## Rows by slice

| Slice | In `controls.tsv` | In `controls_local_sequence.tsv` | Both |
|---|---|---|---|
| adhesion | 31 | 5 | 36 |
| wall | 76 | 15 | 91 |
| fusion | 5 | 0 | 5 |

## Rows kept out of `controls.tsv` because the sequence is not public

`controls_local_sequence.tsv` holds 20 rows in 16 Kiss clusters. Their sequence is in a JGI file that stays on the owner's disk (key `jgi:<portal>:<proteinId>`; `sequence_sha256` is the check). The L1 validator does not read this file: it requires every accession to have a sequence in `sequences.faa`. `build_kiss_fasta.py` joins them back for Task 5.3.
Both files together: 132 rows in 66 clusters.

## Evidence level (both files)

| evidence_level | rows |
|---|---|
| row_verified | 79 |
| locus_tag | 28 |
| locus_confirmed | 5 |
| portal_validated | 20 |

## Other files

| File | Rows |
|---|---|
| `family_reps.tsv` (never test controls) | 35 rows, 30 sequences |
| `rejected.tsv` | 13 |

## Size targets (COMMON-RULES rule 7)

The rule gives 20 clusters per class for `estimated`. The L1 validator reports the distance only. A clade floor of 5 clusters per stratum applies to the hold-out (owner decision 2026-10-22); see `docs/reports/2026-10-18-kiss-feasibility.md` section 10 and `analysis/wall_fusion_truth/results/clusters_per_clade.tsv`.
