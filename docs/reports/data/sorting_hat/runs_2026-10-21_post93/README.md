# Classifier runs of 2026-10-21 (Task 3.10 re-run after issue #93)

Run records of the nine proteomes run again on HPCC after PR #96 (species-limited modules, `ENGINE_SEMANTICS` 3). They replace `../runs_2026-10-17_pr3/` as the record that fits the current `categories.yaml`. The older directory is kept as the record of the first run. Calls tables and evidence tables are not committed.

| file | what it is |
|---|---|
| `<proteome>/report.md` | the report the engine wrote for `out_final` |
| `<proteome>/run.json` | the run record of `out_final` (adds `absent_species_limited` and `applicable_taxa_mismatch`, both empty) |

## What the runs are tied to

| item | value |
|---|---|
| code and data commit pin | `7ed206553ed141e77205a557afbc6a5a8a4eb44a` (HPCC clone `adhesionPred-3_10`, clean tree) |
| `src/cellsurface_sorting_hat/categories.yaml` sha256 | `f50e16844ab62bd590f10cdbe70ddb5217c67ff8dd702e7afbf17ba1b538da37` |
| `data/sorting_hat/family_table.tsv` sha256 | `f15c6f72dc79a79a2f4a128db311551353260efb80f8d63bd843e56cb75031af` (unchanged) |

## Jobs

| job ID | what |
|---|---|
| 30057968 | job 3a, `out_pre`, 9 proteomes, all COMPLETED |
| 30058018 | job 4, repeat call status for S288C and *C. albicans*, both COMPLETED |
| 30219753 | job 3b, `out_final`, 9 proteomes, array 0-8, all COMPLETED (30 to 57 s per task) |

## Result

Every proteome has an `antigen_lookup` table (job 2b), so the species-limited rule does nothing here. Compared cell by cell with the first run (`out_final_before_93`), `calls.wide`, `calls.long`, `evidence` and `proteins` have 0 differing lines in all nine proteomes. `config_sha256` changed from `68a35bad` to `f50e1684`; the number of proteins is the same in each. In `run.json` the `config_sha256` and the two new empty fields are the only differences from the first run.

The two repeat call status files are in `../call_status/{Scer_S288C,Calb_SC5314}/tandem_repeat_protein.owner_definition.config-of-2026-10-21-post93.json`. They differ from the `config-of-2026-10-17-pr3` files only in `config_sha256` (`f50e1684`) and `call_hash` (`8644ca87`). Status `smoke`, `negative_basis` assumed, the measured sensitivity and specificity are the same (S288C 0.435 and 0.987; *C. albicans* 0.462 and 0.976). The `config-of-2026-10-17-pr3` files are stale for the current config and stay as the record.
