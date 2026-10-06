# Agent task directives: control datasets for calibration

*Written 2026-10-06. Each file is a task that an agent can do alone, later. Each task builds
positive and negative control datasets for one line of evidence of the sorting tool
(`cellsurface_sorting_hat`). The controls let us measure a module and give it a status
(`estimated`, `smoke`, `unvalidated`). Background: `docs/paper/03-status-and-validation-rules.md`.*

## How to hand over a task

1. Give the agent `COMMON-RULES.md` and one task file. Nothing else is needed to start.
2. Tell it to work on a new branch and to open a PR. It must not merge.
3. The owner reads `counts.md` and `unverified_notes.md` in the PR before any calibration uses the set.

A one-line hand-over:

> Read `docs/agent-tasks/COMMON-RULES.md` and `docs/agent-tasks/<task file>`. Do the task on a new
> branch. Write the files named there. Open a PR. Do not merge.

## The tasks

| File | Evidence line (module, call) | What exists today | Main gap |
|---|---|---|---|
| `01-onygenales-step1-truth.md` | R0 in *Coccidioides* and relatives (`step1_rule@R0`, `signal_peptide_protein`) | Six species with Phase C truth. None is Onygenales | No Onygenales truth. R0 is `unvalidated` in *Coccidioides* |
| `02-step1-more-positives-existing-species.md` | R0 in the six Phase C species | Only *A. nidulans* is `estimated` | Reasons for `smoke` in the other five; extra direct-evidence controls |
| `03-repeat-mechanism-controls.md` | `tandem_repeat_protein` (`repeat02`, `repeat14`) | 98 E1/E2 adhesin rows in 50 clusters; 42 hard negatives in 33 clusters (count of 2026-10-06) | No mechanism label; 5 clusters flagged repeat-mediated from text; two detectors share one call |
| `04-wall-family-domain-controls.md` | `wall_family_domain` (`pfam_adhesion`) | 15 families, all inactive | No reviewed member and non-member lists |
| `05-allergen-controls.md` | `iuis_allergen_similarity`, `iuis_allergen_homolog` | 111 usable IUIS sequences | No negatives; count of IgE-negative records unknown |
| `06-serodiagnostic-antigen-controls.md` | `cocci_specificity_rank_top15`, `serodiagnostic_marker_candidate` | 4 anchors (tuned); 86 IEDB rows | No independent truth; feasibility of leave-species-out unknown |
| `07-hard-negatives-shared.md` | Shared by step 1, repeat and family-domain calls | 31 seed rows; 42 rows in `adhesins.tsv` | Mostly *S. cerevisiae*; few non-yeast look-alikes |

No control set is needed for the evidence-only modules (Cys-rich finder, expression, TMHMM). They
make no accuracy claim.

Not control datasets, so not here: the *A. fumigatus* strain stability comparison and the module
agreement (kappa) measurement. They are Plan 2, Task 14, and need HPCC runs.

## Suggested order and why

| Order | Task | Reason |
|---|---|---|
| 1 | 07 hard negatives | Cheap. Tasks 03 and 04 use its output. |
| 2 | 03 repeat mechanism controls | The 20-cluster floor is the nearest to being met (50 adhesin clusters, 5 labelled repeat). The label is the gap. |
| 3 | 01 Onygenales step 1 truth | The paper applies the tool to *Coccidioides*. R0 has no status there. |
| 4 | 02 more positives, existing species | Needs no new data to start. It gives the reason that four species stay `smoke`. |
| 5 | 05 allergen controls | Start with the count of IgE-negative records. The answer may close the question. |
| 6 | 04 family domain controls | Needs an HPCC job and a human review of hits. The owner signs off each family. |
| 7 | 06 serodiagnostic controls | The least likely to reach a calibration. Start with the counts. |

This order is a recommendation. The owner decides.

## What these tasks do not do

- They do not change module code, thresholds, `categories.yaml` or any status file.
- They do not run calibration. After the owner accepts a set, the calibration commands
  (`cellsurface_sorting_hat_calibrate truth`) are run as a separate step. `truth` takes one species
  and one call that reads one module.
- They do not claim that a protein is an adhesin, antigen or allergen.
