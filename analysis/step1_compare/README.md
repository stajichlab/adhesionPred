# step1_compare: the step 1 GO truth set (Phase A)

This folder implements deliverables D1, D8, D10 and the extraction half of D9 of
`docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md`.
The plan is `docs/superpowers/plans/2026-10-01-step1-truth-set.md`.
`COLUMNS.md` lists the columns of every output file.

## Requirements

- Use Python 3.12 or later: `PY=/usr/bin/python3.12`.
- The code uses `dataclass(slots=True)`. The default `python3` on the HPCC (3.9) fails.
- Every module in this folder uses the Python standard library only. A test checks this with
  the `ast` module (`test_every_module_imports_only_stdlib_or_local` in
  `tests/step1_compare/test_paths.py`). The subfolder `jobs/` is outside this rule: its scripts
  need numpy, torch and fair-esm (conda env `adhesionPred`) or the PredGPI module Python.
- No script uses `BASH_SOURCE`. A test checks this (`test_no_shell_script_uses_bash_source`).
- Scripts 00 to 07 need no SLURM job. Run them on a login or interactive node. The Phase B
  jobs J0, J1 and J2 (`jobs/*.sh`) need a GPU. They write to `${SCRATCH:?}` and copy the
  results to /bigdata (see "Phase B" below).
- Data never go into the repository.

## Work directory and variables

| Variable | Meaning |
|---|---|
| `PROJ_ROOT` | Repository root. If it is not set, the code uses the path two levels above `paths.py`. |
| `STEP1_WORKDIR` | Output root. If it is not set, it is `<workdir in config/site.yaml>/step1_compare`. |
| `STEP1_GO_DIR` | Test only. A folder that holds the real GO files (for example `/tmp/glyco_spec`). |

Downloads go to `$STEP1_WORKDIR/downloads`. All other outputs go to `$STEP1_WORKDIR`.
`paths.workdir()` always returns an absolute path. A relative `STEP1_WORKDIR` is resolved
against the current directory. A relative `workdir` in `config/site.yaml` is resolved against
the repository root, not against the current directory.
- With `STEP1_WORKDIR` and `PROJ_ROOT` unset, the default is `<workdir in config/site.yaml>/step1_compare`.
  `config/site.yaml` gives the absolute path of the main checkout's `_workdir`
  (`/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir`, git-ignored). A run from a git
  worktree therefore writes into the main checkout's `_workdir`, not into the worktree.

## Re-run

```bash
export PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred
export STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
PY=/usr/bin/python3.12
cd "$PROJ_ROOT/analysis/step1_compare"
$PY 00_fetch_inputs.py                 # download; check SHA-256 against manifest.tsv
$PY 01_extract_go_truth.py             # truth_set.tsv.gz, counts.tsv, extract_log.json
$PY 02_attach_sequences.py             # truth_sequences.tsv.gz, unmatched_ids.tsv,
                                       #   sequence_counts.tsv, sequence_run.json
$PY 03_triage_pm.py                    # d8_triage.tsv, d8_gpi_outside_pext.tsv, d8_counts.tsv,
                                       #   truth_set_triaged.tsv.gz, d8_uniprot/, d8_run.json
$PY 04_build_keyword_tier.py --fetch   # keyword_tier.tsv.gz, keyword_tier_removed.tsv,
                                       #   keyword_tier_run.json, keyword_sequences.fasta.gz,
                                       #   keyword_sequences.json
```

Run the scripts in this order. Run `<script> --help` for all options.
The download sizes in `manifest.tsv` add up to 83,309,908 bytes (19 files).

### Dependencies between the scripts

| Script | Reads from the work directory | Reads from the repository |
|---|---|---|
| 00 | nothing (it downloads) | `manifest.tsv` |
| 01 | the files that 00 downloaded | `species.tsv`, `manifest.tsv` |
| 02 | the output of 01 (`truth_set.tsv.gz`, `extract_log.json`) and the FASTA files from 00 | `species.tsv`, `manifest.tsv` |
| 03 | the output of 01 (`truth_set.tsv.gz`, `extract_log.json`). It does not read the output of 02. | `species.tsv`, `curated_gpi.tsv` |
| 04 | the output of 01 (`extract_log.json`, and the SHA-256 of `truth_set.tsv.gz`) and of 02 (`truth_sequences.tsv.gz`, `sequence_run.json`). It does not read the output of 03 or `d8_run.json`. | `species.tsv`, `data/curated/surface/surface.tsv`, `data/curated/adhesins/eurotiomycetes_seeds.tsv` |

Scripts 03 and 04 are independent of each other.

### Network access

- 00 needs the network.
- 01 and 02 need no network. Script 02 reads the FASTA files that 00 downloaded.
- 03 always queries UniProtKB REST. It has no option to work offline.
- 03 skips a source with no PM candidate: it sends no UniProt query for it (not even the
  organism GPI query), and the source has no row in `d8_counts.tsv`.
- 04 needs the network only with `--fetch` when `keyword_sequences.fasta.gz` is missing.
  If the file exists, 04 uses it and does not fetch. 04 refuses a cached file that has no
  matching sidecar `keyword_sequences.json`.
- 04 downloads the cache before it checks `all_sources` and the truth table hash. A STOP from
  those checks leaves the downloaded cache in place.

### Stale outputs

- Script 02 records `truth_set_sha256` (the SHA-256 of the `truth_set.tsv.gz` it read) in
  `sequence_run.json`.
- Script 04 compares that value with the SHA-256 of the current `truth_set.tsv.gz`. On a
  difference, 04 prints `STOP:`, tells you to re-run 02, exits with code 2 and writes no output.
  A missing `sequence_run.json` or a missing value also stops 04. `--allow-partial-truth-set`
  does not disable this check.
- Script 04 records `truth_set_sha256` in `keyword_tier_run.json`.
- Script 03 records `truth_set_sha256` in `d8_run.json`, so its outputs show which truth table
  they came from. Nothing compares this value or refuses stale 03 outputs.
- After any re-run of 01, or a run of 01 with `--sources`, re-run 02, 03 and 04.

## The STOP contract

- A script that cannot continue prints `STOP: <reason>` to stderr and exits with code 2.
- Script 00 differs. A strict hash mismatch prints `FAILED ... rerun with --update-manifest`.
  Script 00 then continues with the other files and exits with code 2 at the end. It prints no
  `STOP:` line for this case. Only a network error or an unknown `--only` name prints `STOP:`.
- Script 00 skips a download only for a strict file that is already in `downloads/` with the
  right hash. It fetches record-mode files every time.
- Causes include a SHA-256 mismatch in strict mode, a missing input, a malformed GAF row, an
  unknown `--sources` value, an HTTP error, and a missing sequence.
- A STOP in 01, 03 or 04 writes no output tables and leaves the previous outputs in place.
- A STOP in 02 differs. It still writes `unmatched_ids.tsv`, `sequence_counts.tsv` and
  `sequence_run.json`, so you can read the cause. It deletes an old `truth_sequences.tsv.gz`.
  It is the only script that deletes an output on a STOP.
- Read the `--help` text of each script for the exact STOP conditions.

## Atomic outputs

- Scripts write outputs to temporary names. They move the files into place with `os.replace`
  after all outputs are complete.
- A failed run leaves the earlier outputs untouched, except in script 02 (see the STOP contract).
- The shared code is `runinfo.atomic_write_all`. `runinfo.py` also holds `git_commit`,
  `python_version` and the `all_sources` guard (`require_full`).

## Partial runs and the `all_sources` markers

- `--sources` selects a subset of the sources in `species.tsv`. A partial run replaces the full
  outputs. Check the marker before you use a file.
- Each run log records `all_sources`. It is `true` only when the run processed every source.
  Script 04 has no `--sources` option. Its `all_sources` is `true` only when both
  `extract_log.json` and `sequence_run.json` say `all_sources: true`.

| Marker file | Written by |
|---|---|
| `extract_log.json` | `01_extract_go_truth.py` |
| `sequence_run.json` | `02_attach_sequences.py` |
| `d8_run.json` | `03_triage_pm.py` |
| `keyword_tier_run.json` | `04_build_keyword_tier.py` |

- Scripts 02, 03 and 04 refuse a truth set whose `extract_log.json` does not say
  `all_sources: true`. Script 04 also checks `sequence_run.json`.
- The flag `--allow-partial-truth-set` overrides this check. Use it only on purpose.
- Script 04 has two more flags: `--allow-missing-sequences` and
  `--allow-missing-seed-sequences`. The second flag disables the hash rule for the seeds
  that lack a sequence.
- Run logs hold no time stamps. They hold input hashes, the git commit, the Python version
  and the arguments.

## Hash-pinned inputs

- `manifest.tsv` pins every input by SHA-256. A row has mode `strict` or `record`.
- In `strict` mode a hash mismatch stops the run. This applies to the GO OBO file, the GAF
  files and the two Cryptococcus `.goa` files (`20846.C_neoformans_JEC21.goa` and
  `313589.C_neoformans_var_grubii_H99.goa`). In `record` mode the hash is logged only. The
  FASTA files use `record` mode.
- Script 00 appends one row per download to `fetch_log.tsv` in `downloads/`.
- `00_fetch_inputs.py --update-manifest` rewrites `manifest.tsv` inside the repository.
- `01_extract_go_truth.py` checks each input against the manifest before it reads the input.
- Accept a new upstream file only on purpose:
  `$PY 00_fetch_inputs.py --update-manifest`. Then rerun the tests and commit the new
  `manifest.tsv`.

## Tests

```bash
cd "$PROJ_ROOT"
/usr/bin/python3.12 -m pytest tests/step1_compare -q                         # fixtures only
STEP1_GO_DIR=/tmp/glyco_spec /usr/bin/python3.12 -m pytest tests/step1_compare -q   # real files
```

- `PROJ_ROOT` need not be set to run the tests. Observed on 2026-10-01: 220 tests collected, and 215 passed
  with 5 skipped, with `PROJ_ROOT` and `STEP1_WORKDIR` unset. With
  `STEP1_GO_DIR=/tmp/glyco_spec`: 219 passed, 1 skipped.
- Without the real files, the tests that need them are skipped.
- Those tests run when the GO files are in `$STEP1_GO_DIR` or in `$STEP1_WORKDIR/downloads`.
- The real-data tests include `test_reproduces_spec_counts_on_real_files`
  and `test_direct_evidence_counts_on_real_files`.
- The CI job for the analysis tests has `continue-on-error: true`. This suite cannot fail the build.

## Truth subsets are filters on stored columns

Every owner choice is a filter on columns in `truth_set.tsv.gz`. A new choice needs no new
extraction.

| Question | Filter |
|---|---|
| Direct-evidence truth (headline metrics and gates) | `homology_only == "no" and label == X`. The `direct_*` columns of `counts.tsv` count it. |
| All non-IEA truth (reported beside it) | `label == X` |
| Hard negatives (N-sec) that also have internal evidence, for example cytosolic enzymes with one high-throughput plasma-membrane annotation | `label == "N-sec" and internal_evidence != ""`. Example: CDC19 (*S. cerevisiae*) is N-sec by rule (plasma membrane HDA; cytosol RCA). The filter is broad: 658 Scer_SGD genes match it (`truth_set.tsv.gz`, 2026-10-01). |
| Ambiguous genes with only high-throughput internal evidence (R-A) | `label == "ambiguous" and internal_evidence_htp_only == "yes"`. `counts.tsv` column `ambiguous_htp_only` counts it. |
| Not used as truth | `label_no_homology == X`. It is stored for traceability. The `nohom_*` columns of `counts.tsv` reproduce `d1_count.py`. |
| *A. nidulans* as second Eurotiomycetes test species | `source_id == "Anid_EMENI"` and `role == "test_clade"` |
| Which *C. neoformans* file is the reference | `Cneo_H99_GOA` has `role == "test_clade"`. `Cneo_JEC21_GOA` and `Cneo_CRYD1` have `role == "alternate_file"`. |
| Basidiomycota truth option (still open) | rows with `in_clade == "Basidiomycota"`. `Umay_MYCMD` has `role == "undecided"`. |

- P-gpi is a list, not a scored stratum. The lists are `d8_triage.tsv` and
  `d8_gpi_outside_pext.tsv`.
- `curated_gpi.tsv` has only a header row. This is by design. Literature curation fills it
  later as separate work.
- Genes with GPI evidence whose label is not P-ext go to `d8_gpi_outside_pext.tsv`. The script
  does not relabel them. The D8 test fixture shows this for GAS1. No real-data run of
  script 03 was checked for GAS1 or TIP1.

## Counts reproduced from the real files

Command: `01_extract_go_truth.py --input-dir /tmp/glyco_spec` with `STEP1_WORKDIR=/tmp/s1w`.
The table is valid for the pinned `manifest.tsv` (verified 2026-10-01) and for the GO files in the
directory given by `STEP1_GO_DIR`. These are the P-ext counts in `counts.tsv`, with the
direct-evidence count (`direct_p_ext`) in brackets.

| source_id | p_ext (direct_p_ext) |
|---|---|
| Scer_SGD | 125 (88) |
| Calb_CGD | 259 (211) |
| Spom_PomBase | 59 (42) |
| Spom_SCHPO-mod | 57 (41) |
| Afum_ASPFU | 132 (23) |
| Anid_EMENI | 211 (113) |
| Cneo_H99_GOA | 11 (9) |
| Cneo_JEC21_GOA | 32 (0) |
| Cneo_CRYD1 | 32 (0) |
| Umay_MYCMD | 62 (10) |

## Phase B: features and embeddings (D2, D3)

The plan is `docs/superpowers/plans/2026-10-01-step1-features-embeddings.md`. `COLUMNS.md`
lists every Phase B output. All Phase B outputs go to `$STEP1_WORKDIR/phaseb/`.

| Step | Where | Command | Output |
|---|---|---|---|
| 05 | login node | `$PY 05_prepare_sequences.py` | members, unique sequences, FASTA |
| J0 | exfab GPU | `sbatch ... jobs/j0_pilot.sh` | `j0/throughput.json`, `j0/gpu_cpu_diff.json` |
| 06 | login node | `$PY 06_plan_embedding.py` | `chunk_plan.tsv`, `job_plan.json` |
| J1 | exfab GPU | `sbatch ... jobs/j1_features.sh` | `signalp/`, `predgpi/` |
| J2 | exfab GPU | `sbatch --array ... jobs/j2_embed.sh` | `emb/<model>/<chunk_id>.npy` |
| 07 | login node | `$PY 07_build_features.py` | `features.tsv.gz`, `feature_coverage.tsv` |
| assemble | login node | `$ENV_PY jobs/assemble_embeddings.py` | `emb/<model>.nterm.npy`, `.cterm.npy` |

Variables used in the commands below (set them in the shell that runs the commands):

```bash
export PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred-feat   # the checkout
export STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
PY=/usr/bin/python3.12
ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python
S1=$PROJ_ROOT/analysis/step1_compare
LOG=$STEP1_WORKDIR/logs
mkdir -p "$LOG"
cd "$S1"
```

Full `sbatch` lines for J0 and J1:

```bash
sbatch --export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR" \
  -o "$LOG/j0.%j.log" -e "$LOG/j0.%j.log" "$S1/jobs/j0_pilot.sh"
sbatch --export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR" \
  -o "$LOG/j1.%j.log" -e "$LOG/j1.%j.log" "$S1/jobs/j1_features.sh"
```

J2 needs two more arguments (see "Rules for the jobs").

- The `jobs/` Python scripts run with the conda env Python
  (`/rhome/jstajich/.conda/envs/adhesionPred/bin/python`) and
  `PYTHONPATH=$PROJ_ROOT/src:$PROJ_ROOT/analysis/step1_compare:$PROJ_ROOT/analysis/step1_compare/jobs`.
  `jobs/predgpi_scores.py` runs with the Python of `module load predgpi/202001` (3.9).
- Every sequence is embedded and scored once: members with the same `seq_sha256` share one row.
- The chunk size and the J2 job count come from the measured J0 throughput (formula in
  `chunk_plan.py`). `06_plan_embedding.py --rate` makes a dry plan from an assumed rate and
  records `rate_source: assumed`.
- J1 and J2 resume: a finished part or chunk is skipped after its hash is checked.
- 07 does not need J2 or assemble.
- Script 07 reads `sequence_run.json`, `d8_run.json` and `phaseb/prepare_run.json`. It stops
  if their `truth_set_sha256` values differ. It also stops if the SHA-256 of the truth sequences
  file recorded in `prepare_run.json` differs from the file now. Re-run 02, 03 and 05 on the
  same truth set.
- Embedding tests need torch and fair-esm:
  `PYTHONPATH=src /rhome/jstajich/.conda/envs/adhesionPred/bin/python -m pytest tests/step1_compare -q`.
  With `/usr/bin/python3.12` those test files are skipped.

### Rules for the jobs

- J0 must run on a GPU (device `cuda`). `jobs/j0_pilot.sh` sets `--device cuda`. Script 06
  stops if the J0 record has a device other than `cuda`. It also stops if the best rate of a
  model is below 1,000 residues/s (`MIN_RATE` in `06_plan_embedding.py`): such a rate means a
  CPU or failed run.
- J0 writes `phaseb/j0/throughput.json`, `phaseb/j0/gpu_cpu_diff.json` and
  `phaseb/j0/nvidia_smi.csv` (GPU name, driver and memory of the J0 node; see `COLUMNS.md`).
  It copies each result to `phaseb/j0/` as soon as it exists. If the job fails or times out
  after the throughput step, `throughput.json` is kept.
- The readers in `feature_parsers.py` stop with `file:line` on a value that is not a number, on
  a duplicate id (SignalP table, SignalP GFF, PredGPI table) and on a PredGPI row whose
  `gpi_prob`, `omega` and `fpr` do not fit its `gpi_call`. `jobs/predgpi_scores.py` removes its
  temporary file and prints `STOP:` (exit 2) on any scoring error.
- J1 splits the FASTA into `J1_PARTS` parts (default 8). If you change `J1_PARTS`, delete
  `phaseb/signalp` and `phaseb/predgpi` first. Old parts have other members. Script 07 stops
  when one sequence is in two parts.
- J2 is an array job. The `#SBATCH` lines cannot read variables. Give `--array=0-<n_jobs-1>`
  and `--time=<time_minutes>` on the `sbatch` command line. Read both values from
  `phaseb/job_plan.json`. Set `J2_JOB_COUNT` to `n_jobs` in `--export`. `embed_chunks.py` stops
  (exit 2) if `J2_JOB_COUNT` differs from `n_jobs`.
- J2 stops (exit 2) if `job_plan.json` has `rate_source` other than `J0` (a dry plan from
  `06 --rate`). `--allow-assumed-plan` overrides this for tests. `j2_embed.sh` never sets it.
- J2 counts the "Warning: batch failed" lines that `get_esm_embeddings` prints for each chunk
  and stores the count as `batch_failures` in the sidecar. It stops (exit 2) when the count is
  above 0, unless `--allow-batch-failures` is given. A failed batch is retried one sequence at
  a time and is slow. Use a smaller `--batch-size`.
- J2 skips a finished chunk only if its sidecar names the same model, `repr_layer` 6 and window
  and the array width equals the model width (320 for `esm2_t6_8M_UR50D`, 480 for
  `esm2_t12_35M_UR50D`; `jobs/embed_constants.py`). The assembly applies the same checks.
- Evaluation scripts must check that `features_run.json` `input_sha256["unique_sequences.tsv.gz"]` equals `embedding_run.json` `unique_sequences_sha256`.

## Phase C: rule versus ML evaluation (E1 to E4)

The spec is `docs/superpowers/specs/2026-10-01-step1-phase-c-evaluation-design.md`; the plan is
`docs/superpowers/plans/2026-10-01-step1-phase-c-evaluation.md`. `COLUMNS.md` lists every
Phase C output. All Phase C outputs go to `$STEP1_WORKDIR/phasec/`.
No Phase C job has been run on the real data. The tests run the scripts on a small fixture with
a stub MMseqs2.

| Step | Where | Command | Output |
|---|---|---|---|
| 08 | login node | `$ENV_PY phasec/08_build_eval_tables.py` | `eval_table.tsv.gz`, `eval_sequences.fasta.gz` |
| C1 | epyc CPU | `sbatch ... phasec/c1_evaluate.sh` | 09: `split_members.tsv.gz`; 10: `scores.tsv.gz`; 11: `metrics.json`, `findings.json`, `proteome_calls.tsv.gz` |
| 12 | login node | `$ENV_PY phasec/12_report.py` | `report.md` |

```bash
export PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred
export STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python
S1=$PROJ_ROOT/analysis/step1_compare
export PYTHONPATH=$PROJ_ROOT/src:$S1:$S1/phasec
$ENV_PY "$S1/phasec/08_build_eval_tables.py"
sbatch --export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR" \
  -o "$STEP1_WORKDIR/logs/c1.%j.log" -e "$STEP1_WORKDIR/logs/c1.%j.log" "$S1/phasec/c1_evaluate.sh"
$ENV_PY "$S1/phasec/12_report.py"
```

- The `phasec/` scripts run with the conda env Python and
  `PYTHONPATH=$PROJ_ROOT/src:$PROJ_ROOT/analysis/step1_compare:$PROJ_ROOT/analysis/step1_compare/phasec`.
  They may import numpy, scipy and scikit-learn; `test_phasec_imports_only_allowed` checks
  this. The top-level modules stay standard library only.
- 08 checks the Phase B hash chain: `features_run.json`
  `input_sha256["unique_sequences.tsv.gz"]` must equal `embedding_run.json`
  `unique_sequences_sha256` and the current file.
- Every script checks the hashes of the files that the earlier scripts wrote (08, then 09, then
  10, then 11, then 12). A changed file stops the next script. Re-run the scripts in order.
- Every T-c `taxon_id` needs a row in `phasec/tc_taxon_clades.tsv`. Add the clade (or
  `other`) when `keyword_tier.tsv.gz` gains a taxon.
- `module load MMseqs2/17-b804f` puts an AVX2 build on PATH. It stops with "Illegal
  instruction" (exit 132) on CPUs without AVX2, for example the abu_dhabi nodes (c01). C1 runs
  on partition epyc. On another node set `STEP1_MMSEQS` to
  `/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs` (no AVX2 needed).
- A job that runs an AVX2 tool must request a node feature that has AVX2 (--constraint=ryzen on
  partition epyc). The rule covers the scripts that run MMseqs2: `phasec/c1_evaluate.sh` and
  `phasec/09_cluster_and_split.sh`. The Phase B jobs (`jobs/j0_pilot.sh`, `jobs/j1_features.sh`,
  `jobs/j2_embed.sh`) run on exfab GPU nodes and use no MMseqs2, so the rule does not cover
  them. `test_avx2_tools_have_a_cpu_constraint` checks every `*.sh` file.
- `phasec/c1_evaluate.sh` requests partition epyc, `--constraint=ryzen`, `-c 16`, `--mem=32G`
  and `--time=4:00:00`. These values are review estimates, not measurements. The first run is a
  pilot. Size later runs from its `wall_seconds` lines. The memory estimate is (workers + 1)
  times the size of one universe (reviewers measured 0.40 GB per worker and 0.78 GB peak for
  one full fit on a synthetic unit; not measured on C1).
- C1 prints `C1 step <step> wall_seconds=<seconds>` per step and copies these lines to
  `phasec/logs/` on every exit, also after a failed step. `C1_STEPS` (default `09 10 11`) runs
  a subset of the steps. Step 11 runs one process with one BLAS thread.
- 09 writes the MMseqs2 output to `phasec/logs/mmseqs.log`. It never deletes a directory that you
  give with `--tmp-dir`; it works in a new subdirectory and removes only that.
- Phase C tests need numpy and scikit-learn:
  `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare -q`. A Python without these libraries
  skips those test files.
