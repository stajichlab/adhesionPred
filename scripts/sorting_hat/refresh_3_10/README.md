# Task 3.10: HPCC job scripts (DRAFT, nothing submitted)

*2026-10-17. Author: Claude Code (claude-sonnet-5-5). Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`, Task 3.10, step 1. Decisions: D21 and D22 of the spec `docs/superpowers/specs/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing-design.md` (section 3.14).*

**Status.** These scripts were written and syntax-checked (`bash -n`) on a laptop. They have **not** been run and **not** submitted. I used `ssh hpcc` only for read-only listings (appendix A). Nothing was written on HPCC. Each submission needs the owner's explicit go-ahead.

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

Variables used below: `PROJ_ROOT` = a clone of this branch on HPCC; `NEW=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat_refresh_3_10` (all outputs); `OLD=.../_workdir/sorting_hat` (read only).

| order | script | needs | submit as | task count |
|---|---|---|---|---|
| 0 | `00_stage_workdirs.sbatch` | the owner has made `$NEW/logs` and a clone `PROJ_ROOT` | one job | 1 |
| 1 | `01_phasec_r0.sbatch` | job 0 | one job | 1 (4 proteomes inside) |
| 2 | `02_pfam_modules.sbatch` | job 0 (not job 1) | one job | 1 (14 proteomes inside) |
| 3a | `03_core_runs.sbatch` with `OUTNAME=out_pre` | jobs 1 and 2 | array 0-8 | 9 |
| 4 | `04_truth_repeat_call.sbatch` | job 3a for tasks 5 (Calb) and 7 (S288C); job 2 | array 0-1 | 2 |
| 3b | `03_core_runs.sbatch` with `OUTNAME=out_final` | job 4 | array 0-8 | 9 |

Jobs 1 and 2 are independent of each other and can run together. `sbatch --dependency=afterok:` is not used on purpose: the owner approves each submission, so each one is a separate decision.

Submit commands (the owner copies and runs them; `<clone>` is the clone path):

```
mkdir -p $NEW/logs        # the only manual write: sbatch needs the log directory
S=<clone>/scripts/sorting_hat/refresh_3_10
sbatch --export=ALL,PROJ_ROOT=<clone> -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log $S/00_stage_workdirs.sbatch
sbatch --export=ALL,PROJ_ROOT=<clone> -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log $S/01_phasec_r0.sbatch
sbatch --export=ALL,PROJ_ROOT=<clone> -o $NEW/logs/%x.%j.log -e $NEW/logs/%x.%j.log $S/02_pfam_modules.sbatch
sbatch --array=0-8%3 --export=ALL,PROJ_ROOT=<clone>,OUTNAME=out_pre -o $NEW/logs/%x.%A_%a.log -e $NEW/logs/%x.%A_%a.log $S/03_core_runs.sbatch
sbatch --array=0-1   --export=ALL,PROJ_ROOT=<clone> -o $NEW/logs/%x.%A_%a.log -e $NEW/logs/%x.%A_%a.log $S/04_truth_repeat_call.sbatch
sbatch --array=0-8%3 --export=ALL,PROJ_ROOT=<clone>,OUTNAME=out_final -o $NEW/logs/%x.%A_%a.log -e $NEW/logs/%x.%A_%a.log $S/03_core_runs.sbatch
```

Common to all jobs: `#!/bin/bash -l`, `set -euo pipefail`, partition `short` (no account line: the existing scripts name none; `sacctmgr` shows account `stajichlab` for `short`, which is the default), `${SCRATCH:?}` with a `trap` that removes it, outputs copied to `/bigdata` under a temporary name and renamed, no `/tmp`, no `BASH_SOURCE`, absolute paths. `module load hmmer/3.4` appears only in job 2 (as in `pfam_hmmsearch.sbatch`), only to read the version.

**Python.** `/usr/bin/python3.12` (the rule of `docs/agent-tasks/COMMON-RULES.md`), in a virtual environment made on `$SCRATCH` with `pip install pyyaml numpy` (as Plan 2, Task 12 step 1). Reason: the package needs `pyyaml` and `numpy`; `/usr/bin/python3.12` has neither (its `site-packages` is empty), and the conda env named in `analysis/step1_compare/README.md` (`/rhome/jstajich/.conda/envs/adhesionPred`) is Python 3.14 and shows no `yaml` in its `site-packages`. The package is not installed: the jobs set `PYTHONPATH=$PROJ_ROOT/src` and run `python -m cellsurface_sorting_hat.cli`, `...modules.cli` and `...calibration.cli`. A consequence: `run.json` will show `version 0+unknown` instead of `0.2.0.dev421`. I assume this is not part of any identity hash (module identities use `tools`, `params` and artefact hashes); the owner should confirm it. If compute nodes cannot reach PyPI, set `WHEEL_DIR=<dir with pyyaml and numpy wheels>` in `--export` (the scripts then use `--no-index`).

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
- **Owner decides.** (a) The location `NEW`. (b) Which proteomes: I staged the 9 run proteomes and 5 hydrophobin-only proteomes (Bbas_ARSEF2860, Fful_Race5, Fgra_PH-1, Pexp_MD-8, Post_PC9). I left out `Cneo_H99` (4 modules, a 512-byte domain table, no committed status) and `Calb_SC5314_withdups`. (c) The venv source (PyPI or `WHEEL_DIR`). (d) That `PROJ_ROOT` is a clone with `51252bb` in it (the branch is not pushed; the HPCC clone `.../projects/adhesionPred` is on branch `allergen-predictor-wip`, so it cannot be used as is).

### Job 1: `01_phasec_r0.sbatch`
- **Purpose.** Write `status/step1_rule@R0.json` again with `calibrate phasec`, so it carries `leakage_exception`, the reworded note, `negative_basis` per stratum (N-int and N-sec `curated`, PM-TM `assumed`) and the per-stratum rates. Without it R0 drops below `estimated` (spec D21, `docs/CLASSES.md:11`). **There is no `--negative-basis` option on `calibrate phasec`:** `phasec.py` sets the curated basis itself (D16); `--negative-basis curated` exists only on `calibrate truth`.
- **Inputs.** `$DATA_ROOT/_workdir/step1_compare/phasec/metrics.json` (10 MB), `eval_table.tsv.gz` (1.5 MB), `clusters.tsv.gz` (1.5 MB); `$OLD/taxdump/names.dmp` (316 MB) and `nodes.dmp` (232 MB); `data/sorting_hat/phasec_set_species.tsv`; `$NEW/<p>/modules/step1_rule@R0.json`.
- **Command (per proteome in `r0=yes`: Afum_Af293_UniProt, Calb_SC5314, Scer_S288C, Scer_S288C_nodubious; these are the four old work directories that had an R0 status file).**
  `python -m cellsurface_sorting_hat.calibration.cli phasec --workdir <scratch copy> --metrics .../metrics.json --clusters .../clusters.tsv.gz --eval-table .../eval_table.tsv.gz --phasec-signalp-module signalp/6-gpu --phasec-signalp-mode fast --set-species data/sorting_hat/phasec_set_species.tsv --names-dmp .../names.dmp --nodes-dmp .../nodes.dmp` (the options of `analysis/sorting_hat_run/run_converters.sh`).
- **Resources.** 1 CPU, 6 GB, 45 min. Estimate 3 to 10 minutes per proteome, 15 to 30 minutes in all (it parses a 10 MB JSON, two 1.5 MB tables and 316 MB of names four times).
- **Outputs.** `$NEW/<p>/status/step1_rule@R0.json` (+ `.sha256`). The log shows `diff` old against new.
- **Success.** `r0-check` prints one line per set with `leakage_exception True negative_basis True`; the log ends with four `wrote` lines. After the run the owner reads the diff: only the new keys and the reworded note should change; a change in the numbers of an old key needs an explanation.
- **Failure leaves.** Nothing in `$NEW` (all four results are built in scratch first; copied only if all four worked). The command refuses on a mismatch of the SignalP module or mode, or of the counts of `metrics.json` and the eval table (`_curated_headline`): that is a real finding to report, not to work around.
- **Owner decides.** That `signalp/6-gpu` and mode `fast` are what Phase C used (taken from `run_converters.sh`); that the four proteomes are the right set.

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
- **Owner decides.** Convert the existing domain tables, or search again (see above). That `pfam_allergen` stays `unavailable` (no active allergen family, as in the committed runs).

### Job 3 (array, twice): `03_core_runs.sbatch`
- **Purpose.** The classifier run on the nine proteomes of `docs/reports/2026-10-07-classifier-runs-and-scale.md`, with the new config, for the before/after counts of D12. Run `out_pre` first (the repeat call status is measured on it), then `out_final` after job 4 (same inputs plus all status files). The counts come from `out_final`.
- **Inputs.** `$OLD/<p>.faa`, `$NEW/<p>` (modules, status), `$OLD/taxdump/nodes.dmp` (232 MB), `categories.yaml`; taxon and protein count from `inputs.tsv` (the taxa of the committed `run.json` files: 451804, 330879, 330879, 746128, 559297, 237561, 246410, 559292, 559292).
- **Command.** `python -m cellsurface_sorting_hat.cli --fasta $OLD/<p>.faa --taxon <taxon> --taxdump $OLD/taxdump/nodes.dmp --workdir $NEW/<p> --out <scratch>/out --config <clone>/src/cellsurface_sorting_hat/categories.yaml` (the call of `run_converters.sh`, with the config named).
- **Resources.** 1 CPU, 8 GB, 1 h per task, `%3` at a time. Estimate 2 to 10 minutes per task (not measured; the engine reads about 7 MB of tables and `nodes.dmp`).
- **Outputs.** `$NEW/<p>/out_pre` and later `$NEW/<p>/out_final`: `calls.long.tsv.gz`, `calls.wide.tsv.gz`, `evidence.tsv.gz`, `proteins.tsv.gz`, `report.md`, `run.json`.
- **Success.** Each log prints `config_sha256` (must be the sha256 of the pinned `categories.yaml` as the engine computes it), `taxonomy_sha256` (the committed runs show `78392a6926b01fadda65ab819b75e0685e62593fffc0efe09c70b04d26cb7f6b`; a difference means another taxdump, say so), and `wrote .../out_*`. `run.json` `n_proteins` equals `inputs.tsv`.
- **Failure leaves.** An existing `$NEW/<p>/$OUTNAME` stops the task before any work. A failed run leaves no output directory (copy and rename are the last step).
- **Owner decides.** (a) That the committed 2026-10-07 numbers are the "before" (they came from commit `9a3202c` and config `f97cb085`, not from a re-run). If the owner wants a measured baseline on the same modules, a third run with the old config is needed; I did not draft it. (b) **Module sets differ.** `Afum_Af293_Fungi5k` and `Scer_S288C_nodubious` have no `cys8_pattern` and no `hydrophobin_relaxed` table (no relaxed-hydrophobin raw table exists for them), so those two runs will show these modules as absent; the other seven have them. The before/after for those two mixes the module set with the config change. (c) The proteome list.

### Job 4 (array): `04_truth_repeat_call.sbatch`
- **Purpose.** Measure the status of `tandem_repeat_protein` again in S288C (taxon 4932) and *C. albicans* (taxon 5476) with the new config (M2d, M5b2). It writes `status/calls/tandem_repeat_protein.json` into the new work directories.
- **Inputs (the same as the committed files `.../call_status/<p>/tandem_repeat_protein.owner_definition.json`).** `analysis/calibration_truth/repeat_call_truth/truth_v3.<p>.tsv` (321 and 225 lines, committed); the run `$NEW/<p>/out_pre/{calls.long.tsv.gz,run.json,proteins.tsv.gz}` (about 0.7 MB); module records of `repeat02` and `repeat14`; `nodes.dmp`.
- **Command.** `python -m cellsurface_sorting_hat.calibration.cli truth --workdir <scratch> --call-status --config categories.yaml --calls-long <scratch>/out/calls.long.tsv.gz --call tandem_repeat_protein --truth truth_v3.<p>.tsv --calibration-set repeat_motif_definition_<S288C|Calb> --taxa <4932|5476> --nodes-dmp nodes.dmp --notes "<text read from the committed file>" --leakage tuned_on_truth --negative-basis assumed` (`--n-boot 2000 --seed 1` are the defaults).
- **Parameters I reused and one I did not change.** `--leakage tuned_on_truth` and the calibration-set names are those of the committed entries. `--negative-basis assumed`: the committed note says "negatives assumed from absent annotation". Writing `curated` would be false. The notes text is read from the committed file, so the wording is identical.
- **Resources.** 1 CPU, 6 GB, 30 min. Estimate 2 to 5 minutes (reads `nodes.dmp`; the bootstrap of 2,000 draws on about 320 rows takes seconds).
- **Outputs.** `$NEW/<p>/status/calls/tandem_repeat_protein.json` (+ `.sha256`).
- **Success.** The log shows `expected ... MATCH`: S288C sensitivity 0.435 [0.067, 0.667], specificity 0.987, 23 positives (16 clusters), 297 negatives (227 clusters); *C. albicans* 0.462 [0.000, 0.788], 0.976, 13 and 211. If the numbers differ the job still writes the file (a new path), exits 3, and the owner stops and reads it (the plan of 2026-10-08 says: if they differ, stop and report). The status stays `smoke` (the leakage cap).
- **Failure leaves.** Nothing, unless only the comparison failed (exit 3): then the new file exists in `$NEW` and nothing is committed. The command refuses a run whose `config_sha256`, module states or module identities do not fit (`_check_call_run`).
- **Owner decides.** (a) Only `owner_definition` is refreshed (it is the current measurement). The ten other `tandem_repeat_protein.*` files (per species: `broad_truth`, `curated_only`, and three `owner_definition.config-*` files) are history under older truth definitions or configs; they cannot be refreshed without re-deciding what they mean. The plan says "all 20 files" must be current before the interim sentence goes. With 8 hydrophobin files kept valid and 2 refreshed, **10 of 20 are current**; the owner decides whether the other 10 are marked as history and the sentence goes, or are re-measured (broad_truth would need `truth_v2.*` and its own calibration-set name: `repeat_curated_plus_uniprot_secreted_S288C` and `_Calb`). (b) Whether `--negative-basis assumed` is final.

## 4. After the jobs (not part of the scripts)
Task 3.10 steps 3 to 5: copy the new status files to `docs/reports/data/sorting_hat/call_status/` as new files (old ones stay), check `stale_reason` is empty for each with the final config (Task 3.1 test 6), take the counts from `out_final/report.md`, write the CLASSES.md sentence and the dated note in the classifier-runs report. The R0 status files are not in git (`docs/CLASSES.md:11`); decide whether to commit them (they are about 9.6 KB each).

## 5. Not found, assumed, open
- **No persistent Python environment with PyYAML on HPCC.** The plan's venv recipe needs network from a compute node (not tested).
- **No clone of this branch on HPCC.** The branch is local (not pushed). The owner pushes it (or copies it) and makes `PROJ_ROOT`.
- **Pfam domain tables are reused** (assumption above). Raw tool tables (SignalP, TMHMM, repeats, BLAST) are reused unchanged: no change of config or table affects them.
- **Module identity of reused modules.** `step1_rule@R0`, `tm`, `repeat02`, `repeat14`, `allergen_homology`, `cys8_pattern`, `hydrophobin_relaxed` and the Cimm lookups are copied unchanged. I assume their identities do not depend on the config or family table; the engine and `_check_run_identity` will refuse a mismatch.
- **`run.json` version** will read `0+unknown` (package not installed).
- **Truth files and clusters** are the committed `truth_v3` files. Their provenance (UniProt release 2026_03, MMseqs2 30%) is in `docs/reports/2026-10-08-repeat-call-calibration.md`. The HPCC copy of the old S288C status file was written under config `e4f3d63d` (a later config than the committed `f97cb085`); the committed files cover both.
- **Memory and run times** are estimates (section 3).
- **Calb input.** The Calb run uses `$OLD/Calb_SC5314.faa` (6,212 proteins, the deduplicated proteome, as in the committed provenance); the truth ids (for example `C1_12510W_A`) are from that proteome.
- **Not drafted:** a baseline re-run with the old config; a fresh `hmmsearch`; `Cneo_H99`; the other ten repeat status files; HPCC work of PR 5 (clustering and similar).

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
