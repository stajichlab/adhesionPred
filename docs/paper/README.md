# Paper notes: logic of the training and testing

*Written 2026-10-06 for the owner and co-authors. These notes summarise work that is already
documented in this repository. They add no new result. Where a note and a source file differ, the
source file wins.*

## Purpose

The paper must explain **why** each training and testing step was done, **what** each step
measured, and **what it did not measure**. These notes give that logic in one place. They do not
replace the specs, reports or plans.

## Files

| File | Content |
|---|---|
| `01-logic-of-the-pipeline.md` | The question, the steps, and the design rules, in the order of the argument. |
| `02-training-and-testing-ledger.md` | One row per training, test or validation activity: data, split, metric, result, status, source. This is the table to build the Methods and Results from. |
| `03-status-and-validation-rules.md` | How a result becomes a status (`estimated`, `smoke`, `unvalidated`). The leakage, species and interval rules. The decisions of 2026-10-06. |
| `04-limits-and-open-questions.md` | What is not measured, what the paper must not claim, and the open decisions. |

## Rules used in these notes

1. A number in these notes has a source file next to it. A number that I re-derived in the
   2026-10-05 and 2026-10-06 sessions says so. A number that I copied from a report and did not
   re-derive says "copied".
2. "Measured" means a run produced the number. "Planned" means a plan or spec describes it and no
   run exists. "Not measured" means nobody has tested it.
3. No claim of function, mechanism or clinical use is made here. All candidate lists are
   hypotheses. No wet-lab validation exists for any candidate.
4. Names follow the 2026-10-05 renames (`signal_peptide_protein`, `tandem_repeat_protein`,
   `wall_family_domain`, `cocci_specificity_rank_top15`, `serodiagnostic_marker_candidate`,
   `iuis_allergen_similarity`, `iuis_allergen_homolog`). Older reports use the older names.

## Suggested outline of the paper

| Section | Built from |
|---|---|
| 1. Problem: "adhesin" is an outcome, not one sequence feature | `01`, section 1 |
| 2. Review of the first classifier and why it was stopped | `02`, rows A1 to A6 |
| 3. Step 1: surface protein prediction, truth set, held-out tests | `02`, rows B1 to B6; `03` |
| 4. Step 2: mechanism evidence modules (repeats, domains) | `02`, rows C1 to C4 |
| 5. Step 3: purpose layers (serodiagnostic marker ranking, allergen similarity) | `02`, rows D1 to D5 |
| 6. The sorting tool: evidence modules, three-valued logic, status per species | `01`, section 4; `03` |
| 7. Software tests and run-level checks | `02`, rows E1 to E4 |
| 8. Limits | `04` |

## Status of the work at the time of writing

- No trained model ships. The step 1 choice (rule, ML or hybrid gate) was decided as a hybrid on
  2026-10-02 in principle. No gate values are set.
- The sorting tool (`cellsurface_sorting_hat`) has its core engine (PR #63, merged) and its module
  wrappers, calibration commands and job scripts (PR #64, open at the time of writing).
- The HPCC runs that would calibrate the modules (Tasks 11 to 15 of Plan 2) have **not** been run.
  Every module status except the Phase C entries for rule R0 is `unvalidated`. No status file has
  been written in a real work directory.
