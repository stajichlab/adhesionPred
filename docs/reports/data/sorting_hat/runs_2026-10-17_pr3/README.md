# Classifier runs of 2026-10-17 (Task 3.10, PR 3 config)

Run records of the nine proteomes that were run again on HPCC with the PR 3 engine and `categories.yaml`. The calls tables and evidence tables are not committed.

## Files

| file | what it is |
|---|---|
| `<proteome>/report.md` | the report the engine wrote for `out_final` |
| `<proteome>/run.json` | the run record of `out_final` (config sha256, module identities, state counts, call status sources) |
| `old_new_counts.tsv` | output of `scripts/sorting_hat/refresh_3_10/compare_old_new_rule.py`: per proteome, the old rule against the new rule on the same `out_final` run |
| `old_new_counts.tsv.by_family.tsv` | the proteins that only the old rule calls, split by Pfam family |

The nine proteomes are Afum_A1163, Afum_Af293_Fungi5k, Afum_Af293_UniProt, Afum_W72310, Bder_ER3, Calb_SC5314, Cimm_RS, Scer_S288C and Scer_S288C_nodubious. There is no `jobs.txt`: the tool jobs of the 2026-10-07 runs (`../runs/<proteome>/jobs.txt`) were not run again (the raw Pfam domain tables and tool tables of 2026-10-08 were reused).

## What the runs are tied to

| item | value |
|---|---|
| code and data commit pin | `51252bbdd11ec2fe5d81e2e6ec205b7b034f0752` |
| clone commit on HPCC | `fdb75a5` (branch `sorting-hat-engine-changes`, PR #92) |
| `data/sorting_hat/family_table.tsv` sha256 | `f15c6f72dc79a79a2f4a128db311551353260efb80f8d63bd843e56cb75031af` |
| `src/cellsurface_sorting_hat/categories.yaml` sha256 | `68a35bad9dedcb8c34e0368153817edeb1ef489da0e12ed9d9243d4e52cc3721` (every `run.json` shows this as `config_sha256`) |

## HPCC jobs

| job ID | what |
|---|---|
| 30043811 | `out_pre`, 8 proteomes |
| 30044277 | the repeat call status (job 4; the two files in `../call_status/{Scer_S288C,Calb_SC5314}/tandem_repeat_protein.owner_definition.config-of-2026-10-17-pr3.json`) |
| 30044416 | `out_final`, the source of every file in this directory |
| 30042875, 30043011, 30043012, 30043515, 30043066 | earlier jobs of the same refresh (before the repeat call); the task did not state which ID is which, so the roles are not assigned here |

The roles of the scripts are in `scripts/sorting_hat/refresh_3_10/README.md`. The R0 status file made by this refresh is in `../call_status/R0/`.

## Reading these files

- `run.json` shows `version` `0+unknown`: the jobs ran with `PYTHONPATH` and the package was not installed.
- The counts are counts of labels. They are not precision or recall.
- The 2026-10-07 numbers (`../runs/`, commit `9a3202c`, config `f97cb085`) span several changes. Use `old_new_counts.tsv` for the effect of the narrow rule.
