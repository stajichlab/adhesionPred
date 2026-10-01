# step1_compare: the step 1 GO truth set (Phase A)

This folder implements deliverables D1, D8, D10 and the extraction half of D9 of
`docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md`.
The plan is `docs/superpowers/plans/2026-10-01-step1-truth-set.md`.
`COLUMNS.md` lists the columns of every output file.

## Requirements

- Use Python 3.12 or later: `PY=/usr/bin/python3.12`.
- The code uses `dataclass(slots=True)`. The default `python3` on the HPCC (3.9) fails.
- Every module uses the Python standard library only. A test checks this with the `ast` module
  (`test_every_module_imports_only_stdlib_or_local` in `tests/step1_compare/test_paths.py`).
- No script uses `BASH_SOURCE`. The folder has no shell script. A test checks this
  (`test_no_shell_script_uses_bash_source`).
- No SLURM job is needed. Run every script on a login or interactive node.
  If you use a job later, write to `${SCRATCH:?}` and copy the results to /bigdata.
- Data never go into the repository.

## Work directory and variables

| Variable | Meaning |
|---|---|
| `PROJ_ROOT` | Repository root. If it is not set, the code uses the path two levels above `paths.py`. |
| `STEP1_WORKDIR` | Output root. If it is not set, it is `<workdir in config/site.yaml>/step1_compare`. |
| `STEP1_GO_DIR` | Test only. A folder that holds the real GO files (for example `/tmp/glyco_spec`). |

Downloads go to `$STEP1_WORKDIR/downloads`. All other outputs go to `$STEP1_WORKDIR`.

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

| Script | Reads from the work directory |
|---|---|
| 00 | nothing (it downloads) |
| 01 | the files that 00 downloaded |
| 02 | the output of 01 (`truth_set.tsv.gz`, `extract_log.json`) and the FASTA files from 00 |
| 03 | the output of 01 (`truth_set.tsv.gz`, `extract_log.json`), plus `species.tsv` and `curated_gpi.tsv`. It does not read the output of 02. |
| 04 | the output of 01 (`extract_log.json`) and of 02 (`truth_sequences.tsv.gz`, `sequence_run.json`). It does not read the output of 03 or `d8_run.json`. |

Scripts 03 and 04 are independent of each other.

### Network access

- 00 needs the network.
- 01 and 02 need no network. Script 02 reads the FASTA files that 00 downloaded.
- 03 always queries UniProtKB REST. It has no option to work offline.
- 04 needs the network only with `--fetch` when `keyword_sequences.fasta.gz` is missing.
  If the file exists, 04 uses it and does not fetch. 04 refuses a cached file that has no
  matching sidecar `keyword_sequences.json`.

### Stale outputs

- Nothing detects stale downstream outputs.
- Script 02 records `truth_set_sha256` in `sequence_run.json`. Scripts 03 and 04 do not compare
  it, or their own inputs, with the current `truth_set.tsv.gz`.
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

## Partial runs and the `all_sources` markers

- `--sources` selects a subset of the sources in `species.tsv`. A partial run replaces the full
  outputs. Check the marker before you use a file.
- Each run log records `all_sources`. It is `true` only when the run processed every source.

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

- `PROJ_ROOT` need not be set to run the tests. Observed on 2026-10-01: 180 tests collected, and 175 passed
  with 5 skipped, with `PROJ_ROOT` and `STEP1_WORKDIR` unset.
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
