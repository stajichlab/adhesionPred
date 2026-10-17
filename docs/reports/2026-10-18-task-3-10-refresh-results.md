# Task 3.10: results of the refresh of the sorting hat runs

*2026-10-18. Author: Claude Code (claude-sonnet-5-5). Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`, Task 3.10. Scripts: `scripts/sorting_hat/refresh_3_10/`. Run records: `docs/reports/data/sorting_hat/runs_2026-10-17_pr3/`. Research use only.*

## 1. What was run and why

PR 3 changed the engine rules, `categories.yaml` and the family table. Old status files and module tables no longer fit. So these were run again on HPCC (2026-10-17):

- the four Pfam modules for 14 proteomes (they must match the new family table);
- the R0 status file for four work directories (it must carry the leakage exception and the negative basis);
- the antigen module for the eight proteomes that were not *C. immitis* (see section 4);
- the classifier on nine proteomes, twice: `out_pre` (to measure the repeat call) and `out_final` (the source of all counts here);
- the status of `tandem_repeat_protein` in S288C and *C. albicans*.

The run pin is commit `51252bb`; the clone commit is `fdb75a5`. The `categories.yaml` sha256 is `68a35bad9dedcb8c34e0368153817edeb1ef489da0e12ed9d9243d4e52cc3721`. All jobs ended with state COMPLETED and exit code 0.

## 2. Jobs

| job ID | script | state | elapsed | peak memory |
|---|---|---|---|---|
| not stated | job 0, `00_stage_workdirs.sbatch` | COMPLETED | not recorded here | not recorded here |
| not stated | job 1, `01_phasec_r0.sbatch` | COMPLETED | 3 min 30 s | not recorded here |
| not stated | job 2, `02_pfam_modules.sbatch` | COMPLETED | not recorded here | not recorded here |
| not stated | job 2b, `02b_antigen_modules.sbatch` | COMPLETED | not recorded here | about 0.19 GB |
| 30043811 | job 3a, `03_core_runs.sbatch`, `out_pre`, 8 proteomes | COMPLETED | not recorded here | about 1.17 GB |
| 30044277 | job 4, `04_truth_repeat_call.sbatch` | COMPLETED | not recorded here | about 0.96 GB |
| 30044416 | job 3b, `03_core_runs.sbatch`, `out_final` | COMPLETED | not recorded here | about 1.17 GB |

The earlier job IDs 30042875, 30043011, 30043012, 30043515 and 30043066 belong to jobs 0 to 3a. The task record does not say which ID is which job, so they are not assigned here. Peak memory of about 1.17 GB is the value for the core runs.

## 3. R0 result

The stop point of S1 was: A. nidulans (taxon 162425) must stay `estimated`, and taxa 4932, 5476, 746128, 5207 and 5270 must stay `smoke`. The new R0 file shows exactly this. A. nidulans has 163 negatives in 150 clusters (and 109 positives in 100 clusters). The five other species are `smoke`.

The file is identical in the four work directories. One copy is committed at `docs/reports/data/sorting_hat/call_status/R0/step1_rule@R0.json` (sha256 `285818ddb8d9466d761e56ecb5509f61d157521e01408e186518f83c8513bf7d`). It carries a `leakage_exception` for every entry; without it the status would be capped at `smoke`. See the README in that directory. The repo loader accepts the file.

Note: the plan does not contain a table of the S1 numbers. I compared the file with the numbers given in the task (163 and 150 for A. nidulans, status per taxon). They match.

## 4. The antigen module finding

The `antigen` module was absent in 13 of 14 work directories (only `Cimm_RS` had it). With decision D5, an absent input makes the `other_*` calls and `primary_label` `not_assessable`. For S288C this affected 6,688 of 6,722 proteins. The fix for this refresh was to run the antigen module for the eight other proteomes. It writes `not_applicable` outside *Coccidioides*, so the call `cocci_specificity_rank_top15` is not an obstacle any more. The structural fix (the engine should not depend on a module that does not apply) is issue #93 and is not done here.

## 5. The repeat call

The two new status files are `tandem_repeat_protein.owner_definition.config-of-2026-10-17-pr3.json` in `docs/reports/data/sorting_hat/call_status/Scer_S288C/` and `.../Calb_SC5314/`.

| species | sensitivity [95% interval] | specificity [95% interval] | positives | negatives | status |
|---|---|---|---|---|---|
| S288C (taxon 4932) | 0.435 [0.067, 0.667] | 0.987 [0.962, 1.000] | 23 | 297 | smoke |
| *C. albicans* (taxon 5476) | 0.462 [0.000, 0.788] | 0.976 [0.938, 0.996] | 13 | 211 | smoke |

The values are unchanged from the committed values. The negative basis is `assumed` (negatives come from absent annotation). The leakage is `tuned_on_truth`, so the status stays `smoke`. The other ten `tandem_repeat_protein.*` files stay as history under older truth definitions or configs (owner decision Q5).

## 6. The narrowing effect (new rule against old rule, same run)

`compare_old_new_rule.py` read the nine `out_final` runs. Old rule: `signal_peptide_protein` AND (`tandem_repeat_protein` OR `wall_family_domain`). New rule: the narrow `cell_wall_adhesion_candidate` as written. Numbers are proteins called.

| proteome | proteins | old rule | new rule | removed |
|---|---|---|---|---|
| Afum_A1163 | 9942 | 11 | 6 | 5 |
| Afum_Af293_Fungi5k | 9161 | 10 | 8 | 2 |
| Afum_Af293_UniProt | 9647 | 13 | 8 | 5 |
| Afum_W72310 | 10556 | 15 | 9 | 6 |
| Bder_ER3 | 8107 | 12 | 8 | 4 |
| Calb_SC5314 | 6212 | 50 | 29 | 21 |
| Cimm_RS | 9910 | 14 | 7 | 7 |
| Scer_S288C | 6722 | 20 | 17 | 3 |
| Scer_S288C_nodubious | 6039 | 16 | 15 | 1 |
| total | 76,296 | 161 | 107 | 54 |

The 54 removed proteins are split by family: CFEM (PF05730) 34, PIR1-like_C (PF22799) 14, Hyr1 (PF15789) 6. `new_only` is 0 and `rule_mismatch` is 0 in every proteome. So the new call never adds a protein, and the written call equals the rule recomputed from its inputs.

## 7. `primary_label` per proteome

Counted with `python3 -I` from `calls.wide.tsv.gz` of `out_final` (the tables are not committed).

| proteome | other_not_surface | other_surface_no_mechanism | mechanism_without_class | surface_attachment_candidate | cell_wall_adhesion_candidate | not_assessable |
|---|---|---|---|---|---|---|
| Afum_A1163 | 9170 | 728 | 21 | 17 | 6 | 0 |
| Afum_Af293_Fungi5k | 8452 | 678 | 15 | 8 | 8 | 0 |
| Afum_Af293_UniProt | 8877 | 727 | 17 | 18 | 8 | 0 |
| Afum_W72310 | 9758 | 754 | 15 | 20 | 9 | 0 |
| Bder_ER3 | 7728 | 350 | 13 | 8 | 8 | 0 |
| Calb_SC5314 | 5834 | 317 | 11 | 21 | 29 | 0 |
| Cimm_RS | 7449 | 300 | 1375 | 9 | 7 | 770 |
| Scer_S288C | 6390 | 290 | 14 | 3 | 17 | 8 |
| Scer_S288C_nodubious | 5727 | 280 | 8 | 1 | 15 | 8 |
| total | 69,385 | 4,424 | 1,489 | 105 | 107 | 786 |

The `cell_wall_adhesion_candidate` column equals the new rule column of section 6. The two *C. immitis* values (1,375 and 770) are large against the other proteomes. I did not look for the cause; the `inconsistent` record (section 9) is for the same proteome.

## 8. What these numbers are not

- They are counts of labels. They are not precision and they are not recall. No label here was checked against truth.
- The 2026-10-07 numbers (`docs/reports/2026-10-07-classifier-runs-and-scale.md`, commit `9a3202c`, config `f97cb085`) span several changes: CFEM and fifteen other families were activated on 10-08 and 10-09, hydrophobin calls were added, and the D4/D12 narrow sign-off was given. A difference between 10-07 and now is not the effect of one change. The only isolating comparison is old rule against new rule on the same run (section 6).

## 9. Caveats

- The runs reuse the 2026-10-08 Pfam domain tables and the raw tool tables. No tool was run again. The job checked that each domain table is complete; it cannot prove that all 24 models were searched.
- Two proteomes (`Afum_Af293_Fungi5k` and `Scer_S288C_nodubious`) have no `cys8_pattern` and no relaxed-hydrophobin table. Their comparison with the others mixes module set and rule.
- `inconsistent` for Cimm_RS is `{antigen_lookup: 6, expression: 1}`. It is unchanged from the committed 2026-10-07 run.
- `run.json` shows version `0+unknown`: the jobs used `PYTHONPATH` and the package was not installed.
- The ten history repeat files are stale by design. They never raise a status.
- The elapsed time and memory of several jobs are not recorded in this report; use `sacct` on HPCC.
