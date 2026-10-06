# Task 02. More step 1 controls for the species that already have a status

*Read `COMMON-RULES.md` first. Written 2026-10-06.*

## Goal

Find out why four of the six Phase C species stay `smoke`, and list the extra direct-evidence
controls that would let them reach `estimated` for rule R0.

## Known today (Phase C, rule R0; re-derived 2026-10-05/06)

| Species set | Positives / negatives | Clusters (pos / neg) | Sensitivity half-width | Specificity half-width | Status |
|---|---|---|---|---|---|
| *S. cerevisiae* | 79 / 3,785 | 58 / 3,156 | 0.104 | 0.007 | `smoke` (label) |
| *C. albicans* | 153 / 459 | 113 / 410 | 0.118 | 0.028 | `smoke` (label) |
| *A. fumigatus* | 19 / 45 | 17 / 38 | 0.131 | 0.046 | `smoke` (19 positives) |
| *A. nidulans* | 109 / 164 | 100 / 151 | 0.094 | 0.023 | `estimated` |
| *C. neoformans* | 7 / 32 | 6 / 31 | 0.271 | 0.102 | `smoke` |
| *U. maydis* | 9 / 28 | 9 / 24 | 0.150 | 0.144 | `smoke` |

Sources: `_workdir/step1_compare/phasec/metrics.json`, `clusters.tsv.gz`, `eval_table.tsv.gz`;
`docs/paper/02-training-and-testing-ledger.md` table B7.

## Questions to answer first (no new data needed)

1. Why does Phase C label the *S. cerevisiae* and *C. albicans* per-species sets `smoke test` when
   they have 79 and 153 positives? Read `docs/superpowers/specs/2026-10-01-step1-phase-c-evaluation-design.md`
   (section 4, ruling C-8) and `analysis/step1_compare/phasec/`. The label may depend on the ML
   candidates and not on R0. State the exact rule and which candidate fails it.
2. The truth set counts 23 direct-evidence surface genes for *A. fumigatus* but the Phase C set has
   19 positives. Where do the 4 go (no sequence, dedupe, ambiguous, other)? Use `sequence_counts.tsv`,
   `unmatched_ids.tsv`, `eval_dedupe_log.tsv`.
3. What is the effect on the status if the label rule uses only R0?

Report these answers in `findings.md` **before** you add controls.

## Positive and negative controls to add

Only for species with a gap: *S. cerevisiae* (sensitivity half-width 0.104, just above the 0.10
limit), *C. albicans* (0.118), *A. fumigatus* (19 positives, one short of the floor, and a half-width
of 0.131). Do not estimate how many rows are needed; report what you find.
- Positives: GO cell wall / extracellular / cell surface terms with direct evidence codes (IDA, HDA,
  IMP, IGI, EXP) that are **newer than** the annotation files used in Phase A, or that Phase A
  rejected for a stated reason that you can now resolve. Compare with
  `analysis/step1_compare/01_extract_go_truth.py` and `labels.py` (rule D1).
- Negatives: only if a stratum is thin (N-sec, PM-TM). Same evidence rules.
- Literature rows: PMID and quote.

## Size target

Report, per species, how many new clusters of each class you found. Do not exceed what the evidence
supports.

## Deliverables

`data/controls/step1-existing-species/` as in `COMMON-RULES.md`, plus `findings.md` with the
answers to the three questions. Do not change anything in `_workdir/step1_compare/phasec/`. A new
Phase C run is a separate task for the owner.

## Do not

- Do not move a protein from negative to positive because R0 misses it.
- Do not use UniProt keywords as truth (they feed the same predictors).
- Do not remove a hard control because it makes the numbers worse.
