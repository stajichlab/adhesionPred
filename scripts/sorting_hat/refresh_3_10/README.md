# Task 3.10: HPCC job scripts (DRAFT, nothing submitted)

*2026-10-17. Author: Claude Code (claude-sonnet-5-5). Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`, Task 3.10, step 1. Decisions: D21 and D22 of the spec `docs/superpowers/specs/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing-design.md` (section 3.14).*

**Status (updated 2026-10-18).** These scripts were run on HPCC on 2026-10-17, each with the owner's go-ahead, and all jobs COMPLETED with exit code 0: jobs 0, 1, 2, 2b, 3a, 4 and 3b. Job IDs known from the task record: 30043811 (job 3a, `out_pre`, 8 proteomes), 30044277 (job 4, repeat call status), 30044416 (job 3b, `out_final`); the earlier jobs were 30042875, 30043011, 30043012, 30043515 and 30043066. Task 3.10 steps 3 to 5 are done in the branch `task-3-10-results`: the new status files and run records are committed, and the results are in `docs/reports/2026-10-18-task-3-10-refresh-results.md`. The text below is the draft as written before the runs; it is kept unchanged.

## Revisions after the 2026-10-17 run (issue #95)

The scripts changed after the run, from the review of #92. The behaviour of the jobs that ran is unchanged, except these four points. The sections below are still the draft as written before the runs.

1. **Pfam model record and check (finding 1).** `scripts/sorting_hat/pfam_hmmsearch.sbatch` now writes `models` (the searched models, each with `acc` and `name`) and `family_table_sha256` into `raw/pfam/provenance.json`, in addition to the old keys (`pfam_hmm`, `resolved`, `release`, `sha256`, `hmmer`). The readers of the file (job 0 and `analysis/sorting_hat_run/run_converters.sh`) use only the old keys and still work. `02_pfam_modules.sbatch` runs `refresh_checks.py models-searched` for each proteome and FAILS if a row of `data/sorting_hat/family_table.tsv` (accession and name) is not in the record. It also FAILS if the provenance has no model list, which is the case for all 14 existing work directories, unless the owner exports **`ALLOW_UNRECORDED_MODELS=1`** (for example `sbatch --export=ALL,PROJ_ROOT=...,ALLOW_UNRECORDED_MODELS=1 ...`). With the override the job logs one `WARNING` line per proteome and goes on; a table row that is missing from a recorded list still fails. Why: a model that was never searched gives no hits, and that looks like a clean zero.
2. **Job 4 on a mismatch (finding 2).** `04_truth_repeat_call.sbatch`: if `refresh_checks.py measure` fails, the result is written as `status/calls/tandem_repeat_protein.json.rejected` and `.json.sha256.rejected`, not under the final name. The job says so and exits with the same code. On success it behaves as before. Why: a re-run after a fix was blocked by the unexpected file.
3. **Headers (finding 3).** `lib.sh`, `00`, `01`, `02`, `02b`, `03` and `04` (`.sbatch`): the header said DRAFT and "not submitted". It now says that the jobs ran on 2026-10-17 and points to the results report. Comments only.
4. **Pinned fallback (finding 4).** `lib.sh`, `make_venv`: without `WHEEL_DIR` it now installs `pyyaml==6.0.3 numpy==2.5.4`, the same pins as with `WHEEL_DIR`.

Tests: `tests/cellsurface_sorting_hat/test_refresh_checks.py`.

**Update 2026-10-17.** The owner answered the six questions of the first draft. The answers are recorded below as `Owner decision 2026-10-17` paragraphs and in section 6. The scripts were not run after this update either.

## 1. What is tied to what

| item | value |
|---|---|
| branch | `sorting-hat-engine-changes` |
| HEAD when I wrote the scripts | `51252bbdd11ec2fe5d81e2e6ec205b7b034f0752` (the code and data commit the jobs are pinned to) |
| `data/sorting_hat/family_table.tsv` sha256 | `f15c6f72dc79a79a2f4a128db311551353260efb80f8d63bd843e56cb75031af` |
| `src/cellsurface_sorting_hat/categories.yaml` sha256 | `68a35bad9dedcb8c34e0368153817edeb1ef489da0e12ed9d9243d4e52cc3721` |
| `data/sorting_hat/phasec_set_species.tsv` sha256 | `b715e19464156dba1ff8d93ca3afb81684fc3ceb5a489bc418284da97346e2a0` |
| truth files (sha256) | `truth_v3.Scer_S288C.tsv` `a83c61e374896943db3b83073bd946b9e965d1c38b0e50a0267a07f8b310acce`; `truth_v3.Calb_SC5314.tsv` `883abcca2bb5861ebc58b41dba0a894ffa2c55d720a23cb6c3576bbd5b10d4ba` |
| detectors (sha256) | `02_repeat_profile.py` `22e6c5a799751858a6e1c8de6ccd200f8dc900042175d25c907afd4283bb8fec`; `14_repeat_detect_general.py` `9d11c7dfa4405484e428e356a7d5a3bd7a433077abb551a0552ca49557eaa62e` (the repeat tables are reused, not re-run) |
| Pfam | release 38.2, `Pfam-A.hmm` sha256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`, HMMER 3.4 (from the old `provenance.json` files) |
| baseline for "before" (D12) | the committed report and tables of 2026-10-07: commit `9a3202c`, config sha256 `f97cb085...`, `docs/reports/data/sorting_hat/runs/*/run.json` |

Every job calls `check_pin` (in `lib.sh`). It stops unless: commit `51252bb` exists in `PROJ_ROOT`; `src`, `data/sorting_hat`, the two detectors and the truth directory are unchanged between `51252bb` and `HEAD`; the working tree is clean; and the two sha256 values above match. If later PR 3 tasks change `src/` or the table, the jobs stop on purpose. The owner then sets `CODE_COMMIT`, `FAMILY_SHA256` and `CATEGORIES_SHA256` in the `--export` list. The report must say so.

## 2. Run order and dependencies

Variables used below: `PROJ_ROOT` = `/bigdata/stajichlab/jstajich/projects/adhesionPred-3_10` (a clone of this branch; section 5); `NEW=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat_refresh_3_10` (all outputs); `OLD=.../_workdir/sorting_hat` (read only).

| order | script | needs | submit as | task count |
|---|---|---|---|---|
| 0 | `00_stage_workdirs.sbatch` | the owner has made `$NEW/logs` and a clone `PROJ_ROOT` | one job | 1 |
| 1 | `01_phasec_r0.sbatch` | job 0 | one job | 1 (4 proteomes inside) |
| 2 | `02_pfam_modules.sbatch` | job 0 (not job 1) | one job | 1 (14 proteomes inside) |
| 2b | `02b_antigen_modules.sbatch` | job 0; the owner's go-ahead (and `$NEW/Cimm_RS/modules/antigen_lookup.json`) | one job | 1 (8 proteomes inside) |
| 3a | `03_core_runs.sbatch` with `OUTNAME=out_pre` | jobs 1, 2 and 2b (a fresh `out_pre`, see Job 2b) | array 0-8 | 9 |
| 4 | `04_truth_repeat_call.sbatch` | job 3a for tasks 5 (Calb) and 7 (S288C); job 2 | array 0-1 | 2 |
| 3b | `03_core_runs.sbatch` with `OUTNAME=out_final` | job 4 | array 0-8 | 9 |

Jobs 1 and 2 are independent of each other and can run together. `sbatch --dependency=afterok:` is not used on purpose: the owner approves each submission, so each one is a separate decision.

Submit commands (the owner copies and runs them):

```
NEW=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat_refresh_3_10
PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred-3_10
WHEEL_DIR=/bigdata/stajichlab/jstajich/projects/adhesionPred-3_10_wheels
mkdir -p $NEW/logs        # the only manual write: sbatch needs the log directory
S=$PROJ_ROOT/scripts/sorting_hat/refresh_3_10
E=ALL,PROJ_ROOT=$PROJ_ROOT,WHEEL_DIR=$WHEEL_DIR
sbatch --export=$E -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log $S/00_stage_workdirs.sbatch
sbatch --export=$E -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log $S/01_phasec_r0.sbatch
sbatch --export=$E -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log $S/02_pfam_modules.sbatch
sbatch --array=0-8%3 --export=$E,OUTNAME=out_pre -o $NEW/logs/%x.%A_%a.log -e $NEW/logs/%x.%A_%a.log $S/03_core_runs.sbatch
sbatch --array=0-1   --export=$E -o $NEW/logs/%x.%A_%a.log -e $NEW/logs/%x.%A_%a.log $S/04_truth_repeat_call.sbatch
sbatch --array=0-8%3 --export=$E,OUTNAME=out_final -o $NEW/logs/%x.%A_%a.log -e $NEW/logs/%x.%A_%a.log $S/03_core_runs.sbatch
```

Common to all jobs: `#!/bin/bash -l`, `set -euo pipefail`, partition `short` (no account line: the existing scripts name none; `sacctmgr` shows account `stajichlab` for `short`, which is the default), `${SCRATCH:?}` with a `trap` that removes it, outputs copied to `/bigdata` under a temporary name and renamed, no `/tmp`, no `BASH_SOURCE`, absolute paths. `module load hmmer/3.4` appears only in job 2 (as in `pfam_hmmsearch.sbatch`), only to read the version.

**Python.** `/usr/bin/python3.12` (the rule of `docs/agent-tasks/COMMON-RULES.md`), in a virtual environment made on `$SCRATCH` by `make_venv` in `lib.sh`. The packages come from the wheelhouse `WHEEL_DIR`, with `pip install --no-index --find-links "$WHEEL_DIR"` (no network from a compute node). Reason: the package needs `pyyaml` and `numpy`; `/usr/bin/python3.12` has neither (its `site-packages` is empty), and the job venv does not see the owner's `~/.local/lib/python3.12/site-packages` (PyYAML 6.0.3, numpy 2.5.2 there). The job environment pins the wheelhouse versions: numpy 2.5.4 and PyYAML 6.0.3 (section 5). The package is not installed: the jobs set `PYTHONPATH=$PROJ_ROOT/src` and run `python -m cellsurface_sorting_hat.cli`, `...modules.cli` and `...calibration.cli`. A consequence: `run.json` will show `version 0+unknown` instead of `0.2.0.dev421`. I assume this is not part of any identity hash (module identities use `tools`, `params` and artefact hashes); the owner should confirm it.

## 3. The jobs

Estimates: run times and memory of jobs 1 to 4 were **not measured** (I ran nothing). They come from the input sizes in appendix A and from the measured times in `docs/reports/2026-10-07-classifier-runs-and-scale.md` (Pfam 2.4 min, repeats 2.2 min for a 8,107-protein proteome). Memory has margin because `nodes.dmp` (232 MB) and `names.dmp` (316 MB) become Python dicts. Record the real values with `sacct` after each job and correct this table.

### Job 0: `00_stage_workdirs.sbatch`
- **Purpose.** Make 14 new work directories under `$NEW` from copies. The old directories are the record of the committed runs and are never written.
- **Inputs.** `$OLD/<p>/modules/*` except the four Pfam modules (5.6 to 7.1 MB per proteome), `$OLD/<p>/raw/pfam/{domtbl.txt,provenance.json}` (domtbl 63 KB for S288C, under 0.5 MB each), `$OLD/<p>/status/pfam_hydrophobin.json` for 8 proteomes; the committed hydrophobin status files in the repo.
- **Command.** In the script: `cp` of the files after `sidecar_ok` (each module file against its `.sha256` file); `refresh_checks.py json-equal` between the HPCC hydrophobin status file and the committed copy; Pfam-A.hmm sha256 and release read from `provenance.json`; the FASTA protein count against `modules/tm.json` `n_rows` and `inputs.tsv`.
- **Resources.** 1 CPU, 2 GB, 30 min. Under 150 MB is copied. Estimate 5 to 10 minutes (copy plus about 1 minute for the venv).
- **Outputs.** `$NEW/<p>/modules/`, `$NEW/<p>/raw/pfam/`, `$NEW/<p>/status/pfam_hydrophobin.json` (+ `.sha256`) for `hyd=yes`.
- **Success.** The log ends with `staged 14 work directories`; `ls $NEW` shows 14 directories and `logs`; no hidden `.<p>.tmp.*` directory.
- **Failure leaves.** Nothing in `$OLD`. A directory is renamed into place only after its copy finished; a stop leaves earlier proteomes complete and removes this job's temporary directory (the `trap`). A re-run stops at the first existing `$NEW/<p>`: the owner removes the directories of a failed job 0 (or sets `NEW` to a new path).
- **Owner decides.** (a) The location `NEW`. (b) Which proteomes: I staged the 9 run proteomes and 5 hydrophobin-only proteomes (Bbas_ARSEF2860, Fful_Race5, Fgra_PH-1, Pexp_MD-8, Post_PC9). I left out `Cneo_H99` (4 modules, a 512-byte domain table, no committed status) and `Calb_SC5314_withdups`. 
  **Owner decision 2026-10-17 (Q1, Q2).** (c) The venv source is the wheelhouse `WHEEL_DIR=/bigdata/stajichlab/jstajich/projects/adhesionPred-3_10_wheels`. (d) `PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred-3_10` exists: a clone of the pushed branch (single branch `sorting-hat-engine-changes`, head `35900de`, clean tree; the sha256 of `family_table.tsv`, `categories.yaml` and `phasec_set_species.tsv` equal section 1). `35900de` holds no change of `src`, `data/sorting_hat` or the pinned analysis files against `51252bb`, so `check_pin` passes. (a) and (b) are not answered in the 2026-10-17 note: the location `NEW` and the proteome list stay as drafted until the owner says otherwise.
- **Owner decision 2026-10-17 (item 7).** `NEW=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat_refresh_3_10`; job 0 stages all 14 proteome work directories of `inputs.tsv` (the 9 run proteomes and the 5 hydrophobin-only ones, which are needed for the hydrophobin identity check). `OLD` stays read-only.

### Job 1: `01_phasec_r0.sbatch`
- **Purpose.** Write `status/step1_rule@R0.json` again with `calibrate phasec`, so it carries `leakage_exception`, the reworded note, `negative_basis` per stratum (N-int and N-sec `curated`, PM-TM `assumed`) and the per-stratum rates. Without it R0 drops below `estimated` (spec D21, `docs/CLASSES.md:11`). **There is no `--negative-basis` option on `calibrate phasec`:** `phasec.py` sets the curated basis itself (D16); `--negative-basis curated` exists only on `calibrate truth`.
- **Inputs.** `$DATA_ROOT/_workdir/step1_compare/phasec/metrics.json` (10 MB), `eval_table.tsv.gz` (1.5 MB), `clusters.tsv.gz` (1.5 MB); `$OLD/taxdump/names.dmp` (316 MB) and `nodes.dmp` (232 MB); `data/sorting_hat/phasec_set_species.tsv`; `$NEW/<p>/modules/step1_rule@R0.json`.
- **Command (per proteome in `r0=yes`: Afum_Af293_UniProt, Calb_SC5314, Scer_S288C, Scer_S288C_nodubious; these are the four old work directories that had an R0 status file).**
  `python -m cellsurface_sorting_hat.calibration.cli phasec --workdir <scratch copy> --metrics .../metrics.json --clusters .../clusters.tsv.gz --eval-table .../eval_table.tsv.gz --phasec-signalp-module signalp/6-gpu --phasec-signalp-mode fast --set-species data/sorting_hat/phasec_set_species.tsv --names-dmp .../names.dmp --nodes-dmp .../nodes.dmp` (the options of `analysis/sorting_hat_run/run_converters.sh`).
- **Resources.** 1 CPU, 6 GB, 45 min. Estimate 3 to 10 minutes per proteome, 15 to 30 minutes in all (it parses a 10 MB JSON, two 1.5 MB tables and 316 MB of names four times).
- **Outputs.** `$NEW/<p>/status/step1_rule@R0.json` (+ `.sha256`). The log shows `diff` old against new.
- **Success.** `r0-check` prints one line per set with `leakage_exception True negative_basis True`; the log ends with four `wrote` lines. After the run the owner reads the diff: only the new keys and the reworded note should change; a change in the numbers of an old key needs an explanation.
- **Failure leaves.** Nothing in `$NEW` (all four results are built in scratch first; copied only if all four worked). The command refuses on a mismatch of the SignalP module or mode, or of the counts of `metrics.json` and the eval table (`_curated_headline`): that is a real finding to report, not to work around.
- **Owner decision 2026-10-17 (item 8).** The settings are `--phasec-signalp-module signalp/6-gpu --phasec-signalp-mode fast`; the proteomes are the four work directories with an R0 status file (Afum_Af293_UniProt, Calb_SC5314, Scer_S288C, Scer_S288C_nodubious). `check_signalp` compares the module and mode with the work directory's R0 run record and stops on a mismatch (fail-safe). **Added success check:** the new R0 status files must still show A. nidulans (taxon 162425) as `estimated` and taxa 4932, 5476, 746128, 5207 and 5270 as `smoke` (the claim that stop point S1 protected); a different result is a finding to report, not to work around.

### Job 2: `02_pfam_modules.sbatch`
- **Purpose.** Rebuild `pfam_adhesion`, `pfam_hydrophobin`, `pfam_hsba` and `pfam_allergen` so their identity matches the new family table (the digest of `pfam_adhesion` now includes `tier` and `narrow_ok`; its params have `narrow_families`; it has the new column `hit_narrow`). Old tables stop the engine with `missing column`.
- **Inputs.** `$NEW/<p>/raw/pfam/domtbl.txt` and `provenance.json`, the FASTA `$OLD/<p>.faa` (3 to 6.5 MB), the `step1_rule@R0` and `tm` tables of the new directory (SP and no-TM conditions), `data/sorting_hat/family_table.tsv`.
- **Command (per proteome, 14).** `python -m cellsurface_sorting_hat.modules.cli pfam --fasta $OLD/<p>.faa --workdir <scratch copy> --domtbl .../domtbl.txt --family-table data/sorting_hat/family_table.tsv --pfam-release 38.2 --pfam-sha256 4b0da639... --hmmer-version 3.4 --sp-module step1_rule@R0 --tm-module tm`.
- **It does not run hmmsearch.** It converts the domain tables of 2026-10-08. Basis: the family table has held the same 24 models since commit `3819f00` (2026-10-08 09:00); every domain table is dated 2026-10-08 12:43 or later; the job checks that each table ends with `# [ok]`, shows `--cut_ga`, and has only hits of accessions in the table. It cannot prove that all 24 models were searched (a model with no hit leaves no trace). If the owner wants proof, run `scripts/sorting_hat/pfam_hmmsearch.sbatch` per proteome first (about 2.4 minutes each, 8 CPUs, 8 GB; a separate decision).
- **Hydrophobin check.** For the 8 proteomes with a committed file (`hyd=yes`), `refresh_checks.py identity-equal` compares `module`, `version`, `params_hash` and `artefact_hash` of the rebuilt `pfam_hydrophobin.json` with the same four keys of `docs/reports/data/sorting_hat/call_status/hydrophobin/pfam_hydrophobin.<p>.json`. The family digest of `pfam_hydrophobin` uses the short tuple (accession, module, second_condition, active), so it should be unchanged by the new columns. If any of the 8 differs the job copies nothing and says which key differed (the file would be stale, `call_status.py:179,194`).
- **Resources.** 1 CPU, 4 GB, 30 min. Estimate 3 to 6 minutes (14 conversions of under 15 seconds each, plus the venv).
- **Outputs.** `$NEW/<p>/modules/pfam_{adhesion,allergen,hydrophobin,hsba}.{json,tsv.gz}` with `.sha256` files. The log lists old and new identities of each module.
- **Success.** The log ends with 14 `wrote four Pfam modules` lines and 8 `same identity` lines. `pfam_adhesion` must show a **new** artefact hash (expected); the other three must show an unchanged artefact hash for `pfam_hydrophobin`, and for `pfam_allergen` and `pfam_hsba` an unchanged one too if their rows did not change (a difference there is informative, not an error).
- **Failure leaves.** Nothing in `$NEW` (copy happens last, after every check); nothing in `$OLD`.
- **Owner decision 2026-10-17 (Q3).** The job converts the existing 2026-10-08 domain tables; no new hmmsearch. If its model check fails for a proteome, only that proteome falls back to a fresh search, and the failure is reported to the owner before anything is run. **As drafted, the script does not do the per-proteome fallback:** any failed check stops the whole job with nothing copied (`fail`), and no search is started. That fits the decision (stop and report); the fresh search per proteome (`scripts/sorting_hat/pfam_hmmsearch.sbatch`, a separate job) would be written and submitted only after the owner has read the failure. Still open: that `pfam_allergen` stays `unavailable` (no active allergen family, as in the committed runs).
- **Owner decision 2026-10-17 (item 9).** `pfam_allergen` is rebuilt as `unavailable` (both allergen families, AltA1 and Allergen_Asp_f_4, are inactive in `family_table.tsv`); its identity must be unchanged and `iuis_allergen_homolog` stays `not_assessable`, as today. Activating those families needs a specificity test and the owner's sign-off, and is not part of this refresh (issue #19).

### Job 2b: `02b_antigen_modules.sbatch` (DRAFT 2026-10-17, not submitted)
- **Owner decision 2026-10-17, option (a).** Only `Cimm_RS` has an `antigen_lookup` record. For the other proteomes the engine sees the call `cocci_specificity_rank_top15` as `not_assessable` (module absent), and D5 turns the other_* calls and `primary_label` into `not_assessable` (S288C: 6,688 of 6,722 proteins). The owner decided to run the `antigen` module on every proteome so it writes `not_applicable` outside *Coccidioides* (`modules/lookups.py` `antigen_rows`). No code or config change. The structural fix is issue #93 and is not done here.
- **Inputs, all found (read-only on HPCC and in the clone, 2026-10-17).** The record `$NEW/Cimm_RS/modules/antigen_lookup.json` holds no input paths, only `artefact_hash` `040ba813...` (the hash of the sorted `name:sha256` lines of the two files). I identified the files by recomputing that hash from them (it matches):

| input | path | size | sha256 |
|---|---|---|---|
| ranking | `$PROJ_ROOT/analysis/cocci_antigens/cocci_antigen_ranking.tsv` (in git; the same bytes are in `/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/cocci_antigens/` and in the `adhesionPred-3_10` clone) | 1,124,804 B, 9,139 rows | `cf04b9d4d26f25b9c78cd073276aebf4bd19966693232a1e3c5903b5369318bc` |
| protein map | `$DATA_ROOT/_workdir/cocci_spherule/ref/protein_map.tsv` (**not in git and not in the clone**; made 2026-10-03 by `analysis/cocci_spherule/00_prepare_refseq_proteins.py`) | 633,921 B, 9,910 rows | `ca0615211e865559207f98f00e3c676769552331871625ef3dbcea7483c83806` |
| proteome FASTA | `$OLD/<p>.faa` | 4.3 to 6.5 MB each (Afum_A1163 6.47, Afum_Af293_Fungi5k 4.94, Afum_Af293_UniProt 6.52, Afum_W72310 5.81, Bder_ER3 4.28, Calb_SC5314 3.13, Scer_S288C 5.51, Scer_S288C_nodubious 5.23) | n/a |
| taxon | `--taxon` from `inputs.tsv` (451804, 330879, 330879, 746128, 559297, 237561, 559292, 559292) | | |
| applicable taxa | `--applicable-taxa 246410` (the record's `params.applicable_taxa`; also the code default) | | |

  The taxon is one integer per proteome, so no taxon file is needed. (`_taxa_args` also accepts `--taxon-map`, a TSV `protein ID<TAB>taxon ID`; it is not used.) The reference record's `params.genes_with_several_ranking_rows` is 131, which the same ranking reproduces; I ran the command on a toy FASTA on a laptop and got `params_hash` `0584b8a2...`, the one of the Cimm record, so the new records should carry the same `params_hash` and `artefact_hash` as the Cimm_RS record.
- **Command (one per proteome, inside one job).** `python -m cellsurface_sorting_hat.modules.cli antigen --fasta $OLD/<p>.faa --workdir <scratch>/<p> --taxon <taxon> --applicable-taxa 246410 --ranking <ranking> --protein-map <protein_map>`. The command that made the Cimm_RS record on 2026-10-07 01:24 is not logged anywhere I could find (`submit_modules.sh` and `run_converters.sh` do not call `antigen`; `Cimm_RS/jobs.txt` lists only the five tool jobs); this is the same command by its options, confirmed by the identity hashes, not by a log.
- **Resources.** `short`, 1 CPU, 2 GB, 20 min. Inputs are 1 to 7 MB each, the module is a table lookup: seconds per proteome (estimate, not measured).
- **Expected outputs.** `$NEW/<p>/modules/antigen_lookup.{json,tsv.gz}` and the two `.sha256` files, for the 8 proteomes. Built in scratch, copied under a temporary name, renamed. Existing files stop the job; nothing is overwritten.
- **Success check (the job fails if any item is false).** Both input hashes equal the table above; each new record has the `params_hash`, `artefact_hash`, `applicable_taxa` and `genes_with_several_ranking_rows` of the Cimm_RS record; `n_rows` equals the FASTA count of `inputs.tsv`; every row is `not_applicable` (a protein that is not a valid sequence is `na_invalid`; these are counted and shown in the log, and a table with any `ok`, `not_in_reference` or `error` row fails, so no protein is called); the sha256 of every other module file in the work directory is the same before and after. After the job the owner can check: `ls $NEW/*/modules/antigen_lookup.json | wc -l` is 9.
- **A failure leaves.** Nothing in `$NEW` if a check fails in the build step (the scratch is removed). A failure in the copy step can leave some proteomes with the module and others without; the stray temporary files (`.antigen_lookup.*.tmp.<pid>`) are removed by the trap, and the finished files are complete. Re-running then stops on the first existing file, so the owner moves the finished files aside, or runs the missing proteomes by hand; I did not draft a partial mode.
- **Not covered / open.** (1) `Afum_A1163`, `Afum_W72310`, `Bder_ER3`: nothing assumed beyond the table above; their taxon IDs are those of `inputs.tsv`. (2) The module note will read "0 applicable proteins: every protein is not_applicable for this table" (the Cimm record, made earlier, has an empty note and no `condition_snapshots`; neither is part of the identity hashes). (3) Whether the engine and `_check_run_identity` accept a module that is `not_applicable` for every protein: this is the code path of the owner's decision and I read `antigen_rows`, but I did not run the engine (no run on HPCC); the first `out_pre` shows it.
- **Steps after job 2b (not part of the script, each needs the owner's approval).** (1) Move each existing `out_pre` aside, a reversible `mv`: `for p in <the 9 run proteomes>; do mv $NEW/$p/out_pre $NEW/$p/out_pre_noantigen; done` (the 8 above plus `Cimm_RS`, which keeps its own record and is re-run only to keep the nine `out_pre` runs on one footing; the owner may leave Cimm_RS out). (2) Submit job 3a again with `OUTNAME=out_pre` (job 3 refuses an existing directory, so the move is needed); the repeat call status of job 4 is measured on this `out_pre`. (3) Check that the other_* calls and `primary_label` are no longer `not_assessable` because of `cocci_specificity_rank_top15` (S288C: the 6,688 proteins). (4) Run job 4 and `out_final` as in the run order.

### Job 3 (array, twice): `03_core_runs.sbatch`
- **Purpose.** The classifier run on the nine proteomes of `docs/reports/2026-10-07-classifier-runs-and-scale.md`, with the new config, for the before/after counts of D12. Run `out_pre` first (the repeat call status is measured on it), then `out_final` after job 4 (same inputs plus all status files). The counts come from `out_final`.
- **Inputs.** `$OLD/<p>.faa`, `$NEW/<p>` (modules, status), `$OLD/taxdump/nodes.dmp` (232 MB), `categories.yaml`; taxon and protein count from `inputs.tsv` (the taxa of the committed `run.json` files: 451804, 330879, 330879, 746128, 559297, 237561, 246410, 559292, 559292).
- **Command.** `python -m cellsurface_sorting_hat.cli --fasta $OLD/<p>.faa --taxon <taxon> --taxdump $OLD/taxdump/nodes.dmp --workdir $NEW/<p> --out <scratch>/out --config <clone>/src/cellsurface_sorting_hat/categories.yaml` (the call of `run_converters.sh`, with the config named).
- **Resources.** 1 CPU, 8 GB, 1 h per task, `%3` at a time. Estimate 2 to 10 minutes per task (not measured; the engine reads about 7 MB of tables and `nodes.dmp`).
- **Outputs.** `$NEW/<p>/out_pre` and later `$NEW/<p>/out_final`: `calls.long.tsv.gz`, `calls.wide.tsv.gz`, `evidence.tsv.gz`, `proteins.tsv.gz`, `report.md`, `run.json`.
- **Success.** Each log prints `config_sha256` (must be the sha256 of the pinned `categories.yaml` as the engine computes it), `taxonomy_sha256` (the committed runs show `78392a6926b01fadda65ab819b75e0685e62593fffc0efe09c70b04d26cb7f6b`; a difference means another taxdump, say so), and `wrote .../out_*`. `run.json` `n_proteins` equals `inputs.tsv`.
- **Failure leaves.** An existing `$NEW/<p>/$OUTNAME` stops the task before any work. A failed run leaves no output directory (copy and rename are the last step).
- **Owner decision 2026-10-17 (Q6).** (a) The committed 2026-10-07 counts (commit `9a3202c`, config `f97cb085`) are the historical "before". The dated note must say that 10-07 to now spans several changes: CFEM and fifteen other families activated on 10-08 and 10-09, the hydrophobin calls, and the D4/D12 narrow sign-off. The comparison that isolates the D4 change is computed from `out_final` with no extra job (section 6). No third run with the old config is drafted or needed. Still open: (b) **Module sets differ.** `Afum_Af293_Fungi5k` and `Scer_S288C_nodubious` have no `cys8_pattern` and no `hydrophobin_relaxed` table (no relaxed-hydrophobin raw table exists for them), so those two runs will show these modules as absent; the other seven have them. The before/after for those two mixes the module set with the config change; the dated note says so. (c) The proteome list.
- **Owner decision 2026-10-17 (item 10).** The nine proteomes of `inputs.tsv` with `run=yes` (the same nine as the committed run directories under `docs/reports/data/sorting_hat/runs/`), submitted twice (`out_pre`, then `out_final`) as `--array=0-8%3` (three tasks at a time, 8G each), each submission with its own go-ahead.

### Job 4 (array): `04_truth_repeat_call.sbatch`
- **Purpose.** Measure the status of `tandem_repeat_protein` again in S288C (taxon 4932) and *C. albicans* (taxon 5476) with the new config (M2d, M5b2). It writes `status/calls/tandem_repeat_protein.json` into the new work directories.
- **Inputs (the same as the committed files `.../call_status/<p>/tandem_repeat_protein.owner_definition.json`).** `analysis/calibration_truth/repeat_call_truth/truth_v3.<p>.tsv` (321 and 225 lines, committed); the run `$NEW/<p>/out_pre/{calls.long.tsv.gz,run.json,proteins.tsv.gz}` (about 0.7 MB); module records of `repeat02` and `repeat14`; `nodes.dmp`.
- **Command.** `python -m cellsurface_sorting_hat.calibration.cli truth --workdir <scratch> --call-status --config categories.yaml --calls-long <scratch>/out/calls.long.tsv.gz --call tandem_repeat_protein --truth truth_v3.<p>.tsv --calibration-set repeat_motif_definition_<S288C|Calb> --taxa <4932|5476> --nodes-dmp nodes.dmp --notes "<text read from the committed file>" --leakage tuned_on_truth --negative-basis assumed` (`--n-boot 2000 --seed 1` are the defaults).
- **Parameters I reused and one I did not change.** `--leakage tuned_on_truth` and the calibration-set names are those of the committed entries. `--negative-basis assumed`: the committed note says "negatives assumed from absent annotation" (D7, D16). Writing `curated` would be false. The notes text is read from the committed file, so the wording is identical.
- **Resources.** 1 CPU, 6 GB, 30 min. Estimate 2 to 5 minutes (reads `nodes.dmp`; the bootstrap of 2,000 draws on about 320 rows takes seconds).
- **Outputs.** `$NEW/<p>/status/calls/tandem_repeat_protein.json` (+ `.sha256`).
- **Success.** The log shows `expected ... MATCH`: S288C sensitivity 0.435 [0.067, 0.667], specificity 0.987, 23 positives (16 clusters), 297 negatives (227 clusters); *C. albicans* 0.462 [0.000, 0.788], 0.976, 13 and 211. If the numbers differ the job still writes the file (a new path), exits 3, and the owner stops and reads it (the plan of 2026-10-08 says: if they differ, stop and report). The status stays `smoke` (the leakage cap).
- **Failure leaves.** Nothing, unless only the comparison failed (exit 3): then the new file exists in `$NEW` and nothing is committed. The command refuses a run whose `config_sha256`, module states or module identities do not fit (`_check_call_run`).
- **Owner decision 2026-10-17 (Q4, Q5).** (Q4) `--negative-basis assumed` is final: the repeat negatives are assumed from missing annotation (D7, D16); the status stays `smoke`. Expected values: S288C sensitivity 0.435, specificity 0.987; *C. albicans* 0.462 and 0.976. (Q5) Only `owner_definition` is refreshed per species. The other ten `tandem_repeat_protein.*` files (per species `broad_truth`, `curated_only` and three `owner_definition.config-*`) stay in place as dated history and are labelled so (section 6). They are not re-measured. The CLASSES.md interim sentence names the current set (the two refreshed repeat files and the eight hydrophobin files) instead of waiting for all 20 files.

## 4. After the jobs (not part of the scripts)
Task 3.10 steps 3 to 5: copy the new status files to `docs/reports/data/sorting_hat/call_status/` as new files (old ones stay), check `stale_reason` is empty for each with the final config (Task 3.1 test 6), take the counts from `out_final/report.md`, write the CLASSES.md sentence (naming the current set, Q5) and the dated note in the classifier-runs report (Q6; section 6). The R0 status files are not in git (`docs/CLASSES.md:11`); decide whether to commit them (they are about 9.6 KB each).

## 5. Found, assumed, open
- **Clone exists (Q1).** `/bigdata/stajichlab/jstajich/projects/adhesionPred-3_10`: single branch `sorting-hat-engine-changes`, head `35900de`, clean tree, the three hashes equal section 1.
- **Wheelhouse exists (Q2).** `WHEEL_DIR=/bigdata/stajichlab/jstajich/projects/adhesionPred-3_10_wheels`, beside the clone and not inside it (so the clean-tree pin rule holds), holding:
  `numpy-2.5.4-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl` (sha256 `fbde6962867ee75b48b0ee29b2b9372ec5d617799dbaf38e82dc0596f2f7738a`) and
  `pyyaml-6.0.3-cp312-cp312-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl` (sha256 `ba1cc08a7ccde2d2ec775841541641e4548226580ab850948cbfda66a1befcdc`).
  A clean-venv install with `--no-index --find-links` from it was tested on the login node (imports yaml 6.0.3, numpy 2.5.4).
- **PyYAML on HPCC (corrects the first draft).** The first draft said HPCC has no PyYAML anywhere. The owner's `~/.local/lib/python3.12/site-packages` has PyYAML 6.0.3 and numpy 2.5.2; the job venv does not see them. The job environment pins the wheelhouse versions (PyYAML 6.0.3, numpy 2.5.4).
- **Pfam domain tables are reused** (assumption above). Raw tool tables (SignalP, TMHMM, repeats, BLAST) are reused unchanged: no change of config or table affects them.
- **Module identity of reused modules.** `step1_rule@R0`, `tm`, `repeat02`, `repeat14`, `allergen_homology`, `cys8_pattern`, `hydrophobin_relaxed` and the Cimm lookups are copied unchanged. I assume their identities do not depend on the config or family table; the engine and `_check_run_identity` will refuse a mismatch.
- **`run.json` version** will read `0+unknown` (package not installed).
- **Truth files and clusters** are the committed `truth_v3` files. Their provenance (UniProt release 2026_03, MMseqs2 30%) is in `docs/reports/2026-10-08-repeat-call-calibration.md`. The HPCC copy of the old S288C status file was written under config `e4f3d63d` (a later config than the committed `f97cb085`); the committed files cover both.
- **Memory and run times** are estimates (section 3).
- **Calb input.** The Calb run uses `$OLD/Calb_SC5314.faa` (6,212 proteins, the deduplicated proteome, as in the committed provenance); the truth ids (for example `C1_12510W_A`) are from that proteome.
- **Not drafted:** a baseline re-run with the old config; a fresh `hmmsearch`; `Cneo_H99`; the other ten repeat status files; HPCC work of PR 5 (clustering and similar).

## 6. The isolating comparison and the history files (Owner decision 2026-10-17, Q5 and Q6)

**Before and after.** The committed 2026-10-07 counts (commit `9a3202c`, config `f97cb085`) stay the historical "before". From 10-07 to now the runs span several changes, not only D4: CFEM and fifteen other families were activated on 10-08 and 10-09, the hydrophobin calls were added, and the D4/D12 narrow sign-off was given. Two proteomes also differ in module set (section 3, job 3). So the 10-07 to now difference must not be read as the effect of D4 alone.

**The isolating comparison, with no extra job.** `compare_old_new_rule.py` reads one or more `out_final` directories (read only) and compares, on the same run, the old rule (`signal_peptide_protein` AND (`tandem_repeat_protein` OR `wall_family_domain`); the broad flag reads all active families) with the new rule (the narrow `cell_wall_adhesion_candidate` as written). It counts `called`, `not_called`, `not_assessable` and `not_applicable` apart, the proteins only the old rule calls, and the same counts for `surface_attachment_candidate`. The old-only proteins are split by the accession in the `pfam_adhesion` `families` evidence field (the families of the broad flag). It makes no call. It runs on a laptop or the login node in seconds; it needs `calls.long.tsv.gz` and/or `calls.wide.tsv.gz` and `evidence.tsv.gz`, nothing else:

```
python $PROJ_ROOT/scripts/sorting_hat/refresh_3_10/compare_old_new_rule.py --out old_vs_new_rule.tsv $NEW/*/out_final
```

A non-zero `rule_mismatch` means the written narrow call differs from the narrow rule recomputed from its inputs: stop and read it.

**History files.** Only `owner_definition` is refreshed per species. These ten stay in `docs/reports/data/sorting_hat/call_status/` unchanged, as dated history under older truth definitions or configs (they are stale under the new config and never raise a status):
`tandem_repeat_protein.broad_truth`, `tandem_repeat_protein.curated_only` and the three `tandem_repeat_protein.owner_definition.config-of-2026-10-08-hydrophobin`, `...config-with-family-evidence` and `...config-with-surface-attachment` files, in each of the `Scer_S288C` and `Calb_SC5314` directories. They are named here and in the CLASSES.md sentence as history. The current set is the two refreshed `owner_definition` repeat files and the eight hydrophobin files.

## 7. What remains before the first submission
The owner's go-ahead for each job (section 2, now including job 2b, which must run before the first `out_pre`; the 2026-10-17 option (a) is recorded in Job 2b); `mkdir -p $NEW/logs`; the run times and memory of the first job checked with `sacct`. The four items that the first six answers did not cover (job 0 location and proteome list, job 1 SignalP settings, job 2 `pfam_allergen`, job 3 proteome list) were decided by the owner on 2026-10-17 (items 7 to 10 above).

## Appendix A: read-only HPCC listings I ran (2026-10-17)

All through `ssh hpcc` (SLURM tools need `bash -lc`). No job was submitted and nothing was written.

| command | result (abridged) |
|---|---|
| `ls -a $P` (`P=/bigdata/stajichlab/jstajich/projects/adhesionPred`); `cat $P/.git/HEAD` | a git clone; `ref: refs/heads/allergen-predictor-wip`; no `.venv` |
| `ls $P/_workdir` | `sorting_hat`, `step1_compare`, ... |
| `ls $P/_workdir/sorting_hat` | 14 proteome work directories plus `Calb_SC5314_withdups`, `calibration`, `taxdump`, `proteomes`, `rbh`, `iuis_fungal_allergens.faa`, `<name>.faa` files |
| `du -sh $S/*` | work directories 4 to 34 MB; `.faa` 3 to 6.5 MB; `taxdump` 603 MB |
| `ls -la $S/taxdump` | `nodes.dmp` 232,276,813 B; `names.dmp` 316,654,187 B; `merged.dmp`; `taxdump.tar.gz` (dated 2026-10-06) |
| `ls`/`cat` of `Scer_S288C/{modules,raw,status,out,logs}`, `Calb_SC5314/...` | 11 module records; `raw/pfam/{domtbl.txt (63 KB),provenance.json}`; `status/{step1_rule@R0.json, calls/tandem_repeat_protein.json}`; `out*` has 6 run directories |
| `cat Scer_S288C/modules/pfam_hydrophobin.json`, `raw/pfam/provenance.json`, `modules/tm.json`, `*.sha256`, `raw/signalp/version.txt` | release 38.2, Pfam-A.hmm sha256 `4b0da639...`, HMMER 3.4; `.sha256` files hold the bare hash; `module=signalp/6-gpu` |
| per work directory: modules, raw, status, out | R0 status in 4 (Afum_Af293_UniProt, Calb_SC5314, Scer_S288C, Scer_S288C_nodubious); `pfam_hydrophobin` status in 8; `Afum_Af293_Fungi5k` and `Scer_S288C_nodubious` have 7 modules (no `cys8_pattern`, `hydrophobin_relaxed`, `pfam_hsba`, `pfam_hydrophobin`) |
| `ls -la $P/_workdir/step1_compare/phasec` and `du -sh` | `metrics.json` 10 MB, `eval_table.tsv.gz` 1.5 MB, `clusters.tsv.gz` 1.5 MB present |
| `grep -c '>' Scer_S288C.faa Calb_SC5314.faa` | 6,722 and 6,212 |
| `ls ~`, `ls ~/.conda/envs`, `ls /usr/lib*/python3.12/site-packages`, `ls $E/lib/python3.14/site-packages \| grep yaml` | conda env `adhesionPred` is Python 3.14, no yaml listed; `/usr/bin/python3.12` site-packages empty; no csh venv |
| `bash -lc 'sinfo -s; squeue -u $USER; scontrol show partition short; sacctmgr -n -P show assoc user=$USER ...; module avail hmmer'` | `short` 2 h limit, QoS short; account `stajichlab` has `short`; `hmmer/3.4` (default) exists; queue holds other jobs of the owner, one pending with `AssocGrpMemLimit` (so request little memory) |
| `df -h /bigdata/stajichlab` | 18 TB free (95% used) |
| job 2b: `ls -la`, `cat` (head) of `$NEW/Cimm_RS/modules/antigen_lookup.json`; `ls -la`, `sha256sum`, `head -3`, `wc -l` of `_workdir/cocci_spherule/ref/protein_map.tsv`; `sha256sum` of the ranking in the old clone and in `adhesionPred-3_10`; `ls -la` of the nine `$OLD/<p>.faa`; `ls` of `$OLD/Cimm_RS` and its `logs`, `jobs.txt`, `rerun_Cimm_RS.log` | record holds `artefact_hash 040ba813...`, `applicable_taxa [246410]`, `genes_with_several_ranking_rows 131`, written 14:00 (the `$OLD` record of 2026-10-07 01:24 has the same size, 398 B; not compared byte by byte); protein map and ranking sha256 in Job 2b reproduce the hash; no log of the `antigen` command exists; Scer_S288C in `$NEW` has no `antigen_lookup` |

## Re-run after issue #93

Issue #93 added the `species_limited` block to `categories.yaml` and raised `ENGINE_SEMANTICS` to 3. The config hash and every call hash changed. The two `tandem_repeat_protein.owner_definition.config-of-2026-10-17-pr3.json` files (S288C and C. albicans) are stale. The call is `unvalidated` until they are made again. The eight `pfam_hydrophobin` files and the R0 file are module statuses and stay valid. The refresh scripts are not changed.

Run again in the order 3a, 4, 3b. The scripts refuse to overwrite, so move the old results aside first. Every step is reversible (`mv` back). Each submission needs the owner's go-ahead. Do not edit the `.sbatch` or `lib.sh` files.

`NEW` is the new work directory root (see `lib.sh`), `PROJ_ROOT` is the HPCC clone. The nine run proteomes are the ones with `run=yes` in `inputs.tsv`; the two with a repeat call status are `Scer_S288C` and `Calb_SC5314`.

1. **Move the old outputs aside.** Job 3 refuses an existing `out_pre` or `out_final` (`03_core_runs.sbatch`, "exists; nothing is overwritten"). Job 4 refuses an existing `status/calls/tandem_repeat_protein.json` (`04_truth_repeat_call.sbatch`), and job 3b needs the new one. In each of the nine work directories:
   ```
   mv $NEW/$p/out_pre   $NEW/$p/out_pre_before_93
   mv $NEW/$p/out_final $NEW/$p/out_final_before_93
   ```
   and in `Scer_S288C` and `Calb_SC5314`:
   ```
   mv $NEW/$p/status/calls/tandem_repeat_protein.json        $NEW/$p/status/calls/tandem_repeat_protein.json.before_93
   mv $NEW/$p/status/calls/tandem_repeat_protein.json.sha256 $NEW/$p/status/calls/tandem_repeat_protein.json.sha256.before_93
   ```
   Keep the `.before_93` suffix: the status loader reads every `*.json` in `status/calls/`, and a moved file with that suffix is not read.
2. **Update the clone to the merge commit of PR #96.** `check_pin` (`lib.sh`) needs `CODE_COMMIT` in the clone, `src`, `data/sorting_hat` and the detector and truth files unchanged from `CODE_COMMIT` to `HEAD`, and a clean tree. In `$PROJ_ROOT`: `git fetch origin`, `git checkout main`, `git pull --ff-only`. Then `git status --short` must print nothing.
3. **Export the pins.**
   ```
   export CODE_COMMIT=$(git -C $PROJ_ROOT rev-parse HEAD)   # the merge commit of PR #96
   export CATEGORIES_SHA256=f50e16844ab62bd590f10cdbe70ddb5217c67ff8dd702e7afbf17ba1b538da37
   ```
   `CATEGORIES_SHA256` is the sha256 of `src/cellsurface_sorting_hat/categories.yaml`. The value above is the one of the branch `species-limited-modules` before the merge (it is also in `tests/cellsurface_sorting_hat/golden/report.expected.md`). Merging does not change the file, but check the value after the merge: `sha256sum $PROJ_ROOT/src/cellsurface_sorting_hat/categories.yaml` must print the same hash. If it does not, stop; use the printed value only after the owner has seen the diff. `FAMILY_SHA256` does not change; keep the value of the first run.
4. **Run, each with its own go-ahead:** job 3a (`OUTNAME=out_pre`, `--array=0-8`), then job 4 (`--array=0-1`), then job 3b (`OUTNAME=out_final`, `--array=0-8`), as in the run order.
5. **Check the result.** All nine work directories already have an `antigen_lookup` table (job 2b). With a table present the new rule does nothing, so the counts of the new `out_final` must equal the counts of `out_final_before_93`. Only the hashes (`config_sha256`, call hashes) change. A difference in a count is a finding: report it, do not explain it away. The `*_before_93` items stay until the owner removes them.

Nothing is submitted without the owner's go-ahead.
