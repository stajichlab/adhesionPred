# Step 1 Phase C: Rule-versus-ML Evaluation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the Phase C evaluation (deliverables E1 to E4): the evaluation table with the class mapping, dedupe and T-c precedence; MMseqs2 clusters and the S1, S2, S3 splits with the T-c leakage controls; the nested fit of B0, B1, R0 to R2, M8, M35, M8-C, M35-C and H under V-go and V-kw; cluster-bootstrap metrics, strata, agreement, calibration and findings in `metrics.json` and `findings.json`; and `report.md`, so that a controller can submit one CPU job after an independent review.

**Architecture:** A new folder `analysis/step1_compare/phasec/` holds seven library modules from the spec (`labelmap`, `dedupe`, `splits`, `rules`, `models`, `metrics`, `bootstrap`), four helper modules (`evalio`, `universe`, `findings`, `agreement`) and the scripts 08 to 12. 08 (login node) builds one row per unique sequence from the Phase A and Phase B outputs. One SLURM CPU job C1 (`phasec/c1_evaluate.sh`) runs 09 (MMseqs2, splits), 10 (fits, in parallel processes) and 11 (metrics with 2,000 cluster-bootstrap resamples as a weight matrix shared by all candidates). 12 (login node) writes `report.md` and refuses any number that is not in `metrics.json` or `findings.json`. Every script follows the Phase B STOP contract and records the SHA-256 of its outputs, which the next script checks.

**Tech Stack:** conda env `adhesionPred` Python 3.14.2 (`ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python`) with numpy 2.4.2, scipy 1.17.0, scikit-learn 1.8.0 for `phasec/`; the Phase A standard-library modules (`truth_table`, `runinfo`, `manifest`, `paths`) through `PYTHONPATH`; MMseqs2 17-b804f (`module load MMseqs2/17-b804f`, AVX2 build; non-AVX2 binary at `/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs`); SLURM partition `epyc`; bash; pytest 9.1.1; ruff 0.3.5 through pre-commit.

**Spec:** `docs/superpowers/specs/2026-10-01-step1-phase-c-evaluation-design.md` (binding; rulings C-1 to C-15; owner answers of 2026-10-01 on the `t` grid, the estimate floor and the `hard_negative` rows). Parent spec: `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md` (sections 2 to 6, 10, 11; in the worktree `/bigdata/stajichlab/jstajich/projects/adhesionPred-spec`). Phase B plan and code (conventions): `docs/superpowers/plans/2026-10-01-step1-features-embeddings.md`, `analysis/step1_compare/`.

**Scope.** E1 (library, scripts 08 to 12, the tests of spec section 6), the code that produces E2 and E3, and E4 (COLUMNS.md and README.md sections). The plan submits no job. The section "Run (controller, after review)" gives the commands and acceptance checks R0 to R6 for the real run. Gate values, a model card and Basidiomycota curation are out of scope (spec 1).

## Facts measured for this plan (2026-10-01)

All numbers come from read-only runs on `$STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare` (the prototype wrote only into a scratch copy that links to these files). "Not run" marks what was not executed.

| Fact | Value | Source |
|---|---|---|
| Phase B hash chain (M3) | `features_run.json` `input_sha256["unique_sequences.tsv.gz"]` = `embedding_run.json` `unique_sequences_sha256` = `e8dadcb0…d7d8` = the current file; `truth_set_sha256` `516daebe…b4d9` in `d8_run.json`, `features_run.json`, `keyword_tier_run.json`; all three say `all_sources: true` | the run JSONs |
| `truth_set_triaged.tsv.gz` | 51,106 genes, 27 columns (`TRUTH_COLUMNS` + `d8_class`, `d8_reason`); `d8_class` values: empty 50,988, `pm-unresolved` 77, `PM-TM` 39, `P-gpi` 2 | prototype count |
| `features.tsv.gz` | 123,567 members, 24 columns; `truth` members 50,844; `uniprot_kw` 3,573; proteome sets: Anid 10,561, Cimm_RS 9,910, Afum 9,647, H99 7,427, Umay 6,805, JEC21 6,740, Scer 6,722, Calb 6,212, Spom 5,126. It has no `internal_evidence_htp_only`, `taxon_id`, `species` or `in_clade` column, so 08 joins `truth_set_triaged.tsv.gz` by `(source_id, gene_id)` | header |
| `keyword_tier.tsv.gz` | 3,092 T-c rows, 3,034 unique hashes, 7 taxa: 330879 *A. fumigatus* 666, 246410 *C. immitis* 555, 443226 *C. posadasii* 511, 237561 *C. albicans* 413, 559292 S288C 350, 498019 *Candidozyma auris* 327, 284593 *Nakaseomyces glabratus* 270 | prototype count |
| GO members (non-alternate, labelled) | 20,873; by source pos / neg / excluded: Scer 113 / 4,193 / 57; Calb 200 / 2,683 / 137; Spom 49 / 4,178 / 20; Afum 123 / 3,081 / 9; Anid 206 / 3,122 / 41; H99 7 / 48 / 4; Umay 59 / 2,540 / 3 | 08 prototype, `build_run.json` |
| Dedupe | 0 hashes with pos and neg; 2 hashes with pos and an ambiguous gene (S000005921 / S000006203; CAL0000180850 / CAL0000182697); 2 multi-source hashes (Afum + Anid); 12 merged hashes whose members differ in `homology_only` | prototype count |
| Ruling C-6 on the real data | 497 T-c rows dropped (`go_label_wins`); by the class of the GO members of their hash: neg 224 rows, excluded 65, pos 208. GO members that share a hash with a T-c row: N-sec 207, N-int 8, PM-TM 9, pm-unresolved 27, ambiguous 38, P-ext 206, P-gpi 1. The spec's reviewer counts (224 negatives, 65 excluded) are confirmed | 08 prototype |
| Labelled genes without a feature row (real) | `build_run.json` `labelled_genes_without_sequence`: `Calb_CGD` 72, `Scer_SGD` 3, `Spom_PomBase` 2. The 75 of the two training sources are negatives (they are the difference between the spec's 6,951 negative members and the 6,876 negative members that 08 counts) | 08 prototype, `build_run.json` |
| 08 on the real data | 23,338 table rows (go pos 752, go neg 19,776, go excluded 266, tc pos 2,544); 23,351 sequences in `eval_sequences.fasta.gz` (13 literature sequences not in the table); 20 literature rows with an accession, 19 positives (HSP60 not; the 9 `hard_negative` rows count, owner decision 2026-10-01, ruling C-9 amended), YPS3 without accession; 17.5 s wall time on c01. The spec's "about 23,400 (reviewer 23,413)" is 23,351 | 08 prototype |
| MMseqs2 module on c01 | `module load MMseqs2/17-b804f` puts `…/17-b804f-avx2/bin/mmseqs` on PATH; `mmseqs version` dies with "Illegal instruction" (exit 132). c01 is an AMD Opteron 6376 (abu_dhabi, no AVX2). The non-AVX2 binary `…/17-b804f/bin/mmseqs` runs (version `b804fbe384e6f6c9fe96322ec0e92d48bccd0a42`). Partition `epyc` (default, AMD Milan) has AVX2 | `module show`, `/proc/cpuinfo`, `sinfo` |
| 09 on the real data (non-AVX2 binary, 2 threads, c01) | 9,737 clusters for 23,351 sequences; 18 min 48 s wall (clustering 2 min 14 s; five `easy-search -s 7.5` runs about 3 min each) | 09 prototype |
| Split sizes (real) | V-go training pool (FULL): 308 pos, 6,842 neg (7,150 unique; spec: 313 and 6,951 members before dedupe; ruling C-14: positives 313 - 2 made `excluded` (share a hash with an ambiguous gene) - 3 merged = 308; negatives 6,951 - 75 without a feature row - 34 merged = 6,842); S1 fold test sizes 1,454 to 1,483; T-c rows removed: S1 0; S2-Calb_CGD 506 (b) + 164 (c); S2-Scer_SGD 633 (b) + 69 (c); S2-Spom_PomBase 464 (b); S3-Basidiomycota 293 (b); S3-Eurotiomycetes 674 (b) + 1,235 (c) | `splits_run.json` of the prototype |
| Direct-evidence test sets (real, pos / neg) | S1 pooled 232 / 4,244; Calb 153 / 459; Scer 79 / 3,785; Spom 33 / 3,668; Afum 19 / 45; Anid 109 / 164; H99 7 / 32; Umay 9 / 28; literature 19 / 0. Below the count floor of 20 direct positives (ruling C-8): Afum, H99, Umay, literature; these are smoke tests whatever their interval widths | prototype count |
| Fit time, one real unit (S1 fold 0, V-go, 5,720 training rows, one BLAS thread, c01) | B0 1.3 s, B1 1.5 s, each rule 1.2 s, M8 6.0 s, M35 7.0 s, M8-C 6.0 s, M35-C 6.3 s, H 23.1 s; each call includes the 3-fold split (about 1.2 s); chosen C = 0.01 for all LR candidates, H variant M8-C, R2 g = highly_probable, t = 0.10. These timings and choices used the old grids (C in 0.01 to 10, `t` from 0.10). The new grids (C-12: six C values; `t` 0.20 to 0.40) were not timed on the real data (not run) | `unit_timing.py` prototype run |
| Bootstrap time (real S1 pooled direct rows: 4,476 rows, 3,477 clusters) | weight matrix 2,000 x 4,476 in 0.4 s; stratum `all` for 20 candidate-variant pairs (binary metrics for all, curve metrics for 14) in 28.0 s on c01 | `boot_timing.py` prototype run |
| Universe load (69,941 sequences, 4 matrices) | 9.2 s with SHA-256 checks, 6.6 s without | prototype |
| Logistic regression on 7,152 x 480 (35M) | 1.4 s (C = 0.01) to 7.3 s (C = 10) per fit, no convergence warning at `max_iter=5000` | prototype |

## Global Constraints

- `phasec/*.py` imports only the standard library, numpy, scipy, scikit-learn and local modules (`phasec/` and the top-level `analysis/step1_compare/*.py`); `test_phasec_imports_only_allowed` pins this. The top-level modules stay standard library only (`test_every_module_imports_only_stdlib_or_local` scans only `analysis/step1_compare/*.py`).
- Run `phasec/` code with `ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python` and `PYTHONPATH=$PROJ_ROOT/src:$PROJ_ROOT/analysis/step1_compare:$PROJ_ROOT/analysis/step1_compare/phasec`.
- STOP contract: a script that cannot continue prints `STOP: <reason>` to stderr, exits with code 2 and writes no output; outputs go to temp names and then `os.replace` (`runinfo.atomic_write_all` through `evalio.write_outputs`).
- Run JSONs record `all_sources`, `truth_set_sha256`, input hashes, `git_commit`, library versions (Python, numpy, scipy, scikit-learn), MMseqs2 version (09), seeds and arguments; no time stamp. Each run JSON records `outputs_sha256`; the next script stops when a file differs (`evalio.require_current`).
- 08 stops unless `features_run.json` `input_sha256["unique_sequences.tsv.gz"]` equals `embedding_run.json` `unique_sequences_sha256` (Phase B review item M3).
- Seed 20261001 (`evalio.SEED`) for every fold split and every bootstrap; 2,000 resamples (`bootstrap.N_RESAMPLES`).
- No file uses `BASH_SOURCE`, not even in a comment.
- SLURM scripts: `#!/bin/bash -l`, `set -euo pipefail`, `${SCRATCH:?}` node-local work, `PROJ_ROOT` and `STEP1_WORKDIR` required (`: "${PROJ_ROOT:?...}"`), explicit `#SBATCH --time`; results on `/bigdata` before the job ends.
- A job that runs an AVX2 tool must request a node feature that has AVX2 (`-p epyc` with `#SBATCH --constraint=ryzen`; ruling C-15). `test_avx2_tools_have_a_cpu_constraint` scans every `*.sh` file of `analysis/step1_compare/` (Task 11).
- Job size: about 1 to 1.5 h of real run time per job. No run time is measured on the real data for C1; the first submission is the pilot and later submissions are sized from its `wall_seconds` lines.
- Compress large text outputs with gzip (`truth_table.write_tsv` with a `.gz` path, byte-stable).
- ruff 0.3.5 through pre-commit, line length 100: `pre-commit run --all-files` before every commit; re-stage if a hook reformats.
- Commit trailer for every implementation commit: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
- Prose in code, docs, logs and the report: Simplified Technical English; no claim without a measured number or a test.

## Review Focus

The spec does not name these five inputs. Each one gets a test in its owning task.

1. **A stratum or test set without positives or without negatives** (for example the N-sec stratum, or the 19 literature rows that have no negatives). Expected: the metric is `null` in `metrics.json` (never NaN, never an exception) and the label is `smoke test` when a recall interval is not defined. A truth block without negatives (the literature set) also stores precision, PR-AUC and precision at recall as `null` (spec 3.2, review I-1), and `report.md` prints recall only for it. Pinned by `test_undefined_metrics_are_nan` (Task 3), `test_estimate_label_rule` and `test_undefined_metrics_are_null` (Task 9) and `test_literature_section_reports_recall_only` (Task 10).
2. **A bootstrap resample without positives** (all positives in one cluster that a resample misses). Expected: that resample does not count; the interval uses the defined resamples and records `n_defined`. Pinned by `test_interval_skips_undefined_resamples` and `test_all_positives_in_one_cluster_gives_some_empty_resamples` (Task 4).
3. **A missing or non-finite feature value** (an empty `sp_prob` or a NaN score). Expected: 08 stops and names the hash; the metric functions refuse non-finite scores. Pinned by `test_missing_feature_value_stops` (Task 2) and `test_non_finite_scores_are_refused` (Task 3).
4. **One sequence in two sources** (2 Afum/Anid hashes exist; a hash in a training and a test source would leak). Expected: one table row that lists both sources and counts in each source's test set; a hash in a training source and a test source stops 09. Pinned by `test_identical_sequence_in_two_sources_is_one_row` (Task 2) and `test_hash_in_a_training_and_a_test_source_stops` (Task 7).
5. **An outer training set with fewer positive clusters than inner folds** (S2-Calb_CGD trains on 109 positives of S. cerevisiae only). Expected: 10 stops and names the unit and the cluster counts, instead of a fold without positives and an undefined Youden threshold. Pinned by `test_inner_folds_stop_without_enough_positive_clusters` (Task 6) and `test_too_few_positive_clusters_stops` (Task 8).

A sixth case found while prototyping: a recall interval of zero width (for example 2 of 2 positives found in every resample) meets the half-width rule. The owner decision of 2026-10-01 adds a count floor (ruling C-8): a test set is an estimate only with at least 20 direct-evidence positives. `metrics.json` records `n_direct_positives`, `floor_met` and `zero_width_recall_interval` beside the label (Tasks 9 and 10).

---

## Design

### Data model

- **Table row** (`eval_table.tsv.gz`): one row per unique `seq_sha256` of the GO truth members (non-alternate sources, label not `unlabelled`) and of the T-c rows that survive ruling C-6. GO members with one hash are merged (`dedupe.merge_group`): list columns hold the sorted unique values; `homology_only` is `no` when any member has direct evidence; `internal_evidence_htp_only` is `yes` when any member says so.
- **Literature row** (`eval_literature.tsv`): one row per seed with an accession; `literature_positive` per spec 3.2. It is a separate test set (S3-Eurotiomycetes), so a literature hash may also be a table row (7 *A. fumigatus* P-ext genes: CspA, rodA, rodB, cfmA, cfmC, gel1, ecm33).
- **Split member** (`split_members.tsv.gz`): one row per split, fold and hash with its `part` (`train`, `train_tc`, `test`, `test_tc`, `test_lit`) and cluster.
- **Unit**: one outer training set, key `<split>|<fold>|<variant>`. 11 units x 2 variants = 22 units: S1 folds 0 to 4, S2-Calb_CGD, S2-Scer_SGD, S2-Spom_PomBase, S3-Basidiomycota, S3-Eurotiomycetes, FULL.
- **Score** (`scores.tsv.gz`): one row per unit, candidate and scored hash: decision value, Platt probability, call.
- **Universe** (`universe.Universe`): one row per Phase B unique sequence (69,941): features, amino-acid composition, embedding rows. FULL units score all of it.

### Splits (spec 3.4, parent 4 step 8)

| Split | Train (GO) | Test | T-c removals |
|---|---|---|---|
| S1 fold k | GO pos/neg rows of the `train` sources (S288C, *C. albicans*) outside fold k | GO rows of fold k (all classes); T-c rows of fold k as `test_tc` (scored, never truth) | a, b (b removes 0 rows by construction: folds are by cluster) |
| S2-Calb_CGD | Scer_SGD | Calb_CGD | a, b, c (taxon 237561) |
| S2-Scer_SGD | Calb_CGD | Scer_SGD | a, b, c (taxon 559292) |
| S2-Spom_PomBase | both yeasts | Spom_PomBase | a, b, c (taxon 284812) |
| S3-Eurotiomycetes | both yeasts | Afum_ASPFU, Anid_EMENI, literature rows | a, b, c (clade Eurotiomycetes: *A. fumigatus*, *C. immitis*, *C. posadasii* T-c rows) |
| S3-Basidiomycota | both yeasts | Cneo_H99_GOA, Umay_MYCMD | a, b, c (clade Basidiomycota: no T-c taxon) |
| FULL | both yeasts | none | none; all T-c rows train under V-kw |

S1 folds come from `StratifiedGroupKFold(5, shuffle=True, random_state=20261001)` over the GO pos/neg training rows with clusters as groups. Excluded rows of the training sources and T-c rows take the fold of their cluster; a cluster without a fold gets `int(sha256(cluster_id)[:16], 16) mod 5`. So every training-source GO row and every T-c row gets an out-of-fold score (spec 4, `score_source`).

In S2 and S3 a GO training protein may share a cluster with a test protein (homology across species). This is by design; the maximum-identity stratum measures it (ruling C-4). `splits.check_no_cluster_spans` therefore checks S1 folds (all training parts) and the T-c training rows of S2 and S3 (ruling C-5), not the GO training rows of S2 and S3.

### Nested protocol (spec 3.3; `models.fit_unit`)

For each unit: inner `StratifiedGroupKFold(3, shuffle=True, random_state=20261001)` over the training rows (groups = clusters). Logistic-regression candidates: for each C in {0.001, 0.003, 0.01, 0.1, 1, 10} (ruling C-12) the inner out-of-fold decision values; C = the highest inner PR-AUC (ties: smaller C). H: the same over (ESM variant, C). With that setting: Youden threshold and Platt (a, b) from the inner out-of-fold values; the final pipeline (StandardScaler + LogisticRegression(class_weight="balanced")) is refitted on all training rows. Rules: g and t (`t` in 0.20 to 0.40, owner decision 2026-10-01) by Youden's J on the training rows (a rule has no fitted model, so this equals J on the pooled inner out-of-fold rows). The fit receives only the training rows' labels (`10.tasks` builds them); `test_threshold_fit_uses_train_only` permutes every outer test label and asserts identical settings.

### Test sets and truth (spec 3.2, 4)

S1:all (pooled out-of-fold), S1:Scer_SGD, S1:Calb_CGD; S2-<source>:<source>; S3-<clade>:clade, S3-<clade>:<source>, S3-Eurotiomycetes:literature. Truth `direct` (`homology_only == no`; literature rows count as direct) is the headline; truth `all` (all non-IEA labels) is beside it. Strata: all, wall, extracellular-only, N-int, N-sec, PM-TM, long (> 1,022 aa), identity_below_0.3 (S2 and S3). Ambiguous genes: score distribution only, with the `internal_evidence_htp_only` sub-stratum. pm-unresolved and P-gpi: lists.

### Bootstrap (spec 4)

`bootstrap.cluster_weights` draws clusters with replacement and returns W (2,000 x n); every metric reads W. One W per test set and truth (seed from `bootstrap.seed_for(20261001, "<test set>|<truth>")`) serves every candidate, variant and stratum, so differences (variant effect, ML minus comparator) are paired. Row 0 of the evaluated matrix is all ones (point value).

### score_source (spec 4, `agreement.score_source`)

`oof` when the hash is in an S1 test part (`test` or `test_tc`); `in_sample` when it is in a FULL training table but in no S1 fold (no such hash exists by construction; the rule is a guard); `final` otherwise (FULL model).

### Storage

`scores.tsv.gz` has about 2.0 million rows on the real data: 20 candidate-variant pairs x (69,941 FULL rows + 9,884 S1 rows + 11,553 S2 rows + 9,261 S3 rows = 100,639), computed from the prototype split sizes; the file was not written on the real data. All tables are gzip. `metrics.json` was 7.9 MB on the test fixture (structure size; the real file is of the same order, not run).

## File structure

| File | Responsibility | Task |
|---|---|---|
| `analysis/step1_compare/phasec/evalio.py` | STOP error, `write_outputs` (run JSON with `outputs_sha256`), `require_current`, seed, variants | 1 |
| `analysis/step1_compare/phasec/labelmap.py` | `class_of`, `stratum_of` (spec 3.1) | 1 |
| `tests/step1_compare/conftest.py` (modify) | `phasec/` on `sys.path`, `load_phasec`; Task 8 adds the session fixture `phasec_chain` | 1, 8 |
| `analysis/step1_compare/phasec/dedupe.py` | dedupe within the GO truth, precedence over T-c (C-6) | 2 |
| `analysis/step1_compare/phasec/tc_taxon_clades.tsv` | clade per T-c `taxon_id` (C-7) | 2 |
| `analysis/step1_compare/phasec/08_build_eval_tables.py` | input checks, table, literature, FASTA | 2 |
| `tests/step1_compare/phasec_fixture.py` | small realistic work directory; Task 7 appends the chain helpers | 2, 7 |
| `analysis/step1_compare/phasec/metrics.py` | weighted metrics on W | 3 |
| `analysis/step1_compare/phasec/bootstrap.py` | cluster weights, seeds, intervals | 4 |
| `analysis/step1_compare/phasec/rules.py` | R0, R1, R2, Youden fit of g and t | 5 |
| `analysis/step1_compare/phasec/universe.py` | per-sequence inputs, C-terminal row choice | 6 |
| `analysis/step1_compare/phasec/models.py` | nested protocol, Platt, work units | 6 |
| `analysis/step1_compare/phasec/splits.py` | clusters, S1-S3, T-c removals, checks, identity parsing | 7 |
| `analysis/step1_compare/phasec/09_make_splits.py` | MMseqs2 calls and split tables | 7 |
| `tests/step1_compare/fixtures/phasec/stub_mmseqs/mmseqs` | test stub of MMseqs2 | 7 |
| `analysis/step1_compare/phasec/10_fit_and_score.py` | units, parallel fits, `scores.tsv.gz` | 8 |
| `analysis/step1_compare/phasec/findings.py` | estimate label, findings (a) to (c) | 9 |
| `analysis/step1_compare/phasec/agreement.py` | score_source, named panel, agreement counts | 9 |
| `analysis/step1_compare/phasec/11_evaluate.py` | metrics.json, findings.json, proteome_calls.tsv.gz | 9 |
| `tests/step1_compare/fixtures/phasec/golden_metrics.json.gz` | frozen metrics.json of the golden fixture | 9 |
| `analysis/step1_compare/phasec/12_report.py` | report.md and the number check | 10 |
| `analysis/step1_compare/phasec/09_cluster_and_split.sh`, `phasec/c1_evaluate.sh` | MMseqs2 wrapper; SLURM job C1 | 11 |
| `analysis/step1_compare/COLUMNS.md`, `README.md` (modify) | Phase C sections | 12 |
| `tests/step1_compare/test_phasec_*.py` (12 files) | tests | 1 to 12 |

## How to run the tests

```bash
cd /bigdata/stajichlab/jstajich/projects/adhesionPred-eval    # the worktree of branch step1-eval
PY=/usr/bin/python3.12
ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python
PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare -q       # full run
$PY -m pytest tests/step1_compare -q                          # stdlib run: numpy test files skip
STEP1_MMSEQS=/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs \
  PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_splits.py -q   # real MMseqs2
```

Baseline before Task 1 (measured on a copy of the worktree at `7b735ec`): full run `409 passed, 1 skipped` in 306 s; stdlib run `360 passed, 7 skipped, 1 failed`. The stdlib failure is old: `test_phaseb_docs.py::test_chunk_manifest_columns_equal_the_code_constant` imports `surface_glyco` without `PYTHONPATH=src`. It is not part of this plan.

The real-MMseqs2 test skips where no working `mmseqs` exists. On c01 the module build stops with exit 132, so set `STEP1_MMSEQS` to the non-AVX2 binary as above.

---

### Task 1: Phase C package scaffold, import rule, class mapping

**Files:**
- Create: `analysis/step1_compare/phasec/evalio.py`, `analysis/step1_compare/phasec/labelmap.py`
- Modify: `tests/step1_compare/conftest.py` (put `phasec/` on `sys.path`; add `load_phasec`)
- Test: `tests/step1_compare/test_phasec_scaffold.py`

**Interfaces:**
- Consumes: Phase A `manifest.sha256_file(path) -> str`, `runinfo.atomic_write_all(out_dir, writers)`, `runinfo.write_json(path, obj)`.
- Produces: `evalio.PHASEC_DIR: Path`, `evalio.SEED = 20261001`, `evalio.VARIANTS = ("V-go", "V-kw")`, `evalio.StopError(ValueError)`, `out_dir(work) -> Path` (`work/"phasec"`), `library_versions() -> dict`, `read_json(path) -> dict` (StopError when unreadable), `write_outputs(out, writers: dict[str, Callable[[Path], None]], run_name: str, log: dict) -> dict` (adds `outputs_sha256`), `require_current(out, run_name, names, producer) -> dict` (the run JSON), `float_or_stop(text, what) -> float`. `labelmap.POS, NEG, EXCLUDED`, `CLASSES`, `D8_CLASSES = ("", "P-gpi", "PM-TM", "pm-unresolved")`, `LABELS`, `SUBSETS`, `class_of(label, d8_class) -> str`, `stratum_of(label, subset, d8_class) -> str`. conftest: `load_phasec(name)` imports `phasec/<name>.py`.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_scaffold.py`:

```python
"""Phase C package rules: allowed imports, no BASH_SOURCE, the class mapping (spec 3.1)."""

import ast
import sys

import labelmap
import paths
import pytest

PHASEC = paths.STEP1_DIR / "phasec"
ALLOWED_THIRD_PARTY = {"numpy", "scipy", "sklearn"}


def imported_names(source: str) -> list[str]:
    names = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            names += [a.name.split(".")[0] for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names.append((node.module or "").split(".")[0])
    return names


def forbidden_imports(source: str, local: set[str]) -> list[str]:
    return [
        n
        for n in imported_names(source)
        if n not in sys.stdlib_module_names and n not in local and n not in ALLOWED_THIRD_PARTY
    ]


def _local() -> set[str]:
    return {p.stem for p in PHASEC.glob("*.py")} | {p.stem for p in paths.STEP1_DIR.glob("*.py")}


def test_phasec_imports_only_allowed():
    modules = sorted(PHASEC.glob("*.py"))
    assert len(modules) >= 2, modules  # evalio and labelmap from Task 1 on
    for path in modules:
        assert forbidden_imports(path.read_text(), _local()) == [], path.name


def test_phasec_import_check_catches_a_forbidden_import():
    path = PHASEC / "labelmap.py"
    copy = path.read_text() + "\nimport pandas\nfrom torch import nn\n"
    assert forbidden_imports(copy, _local()) == ["pandas", "torch"]


def test_no_phasec_file_uses_bash_source():
    files = sorted(PHASEC.rglob("*.py")) + sorted(PHASEC.rglob("*.sh"))
    assert files
    for f in files:
        assert "BASH_SOURCE" not in f.read_text(), f


CASES = [
    ("P-ext", "", "pos", "wall"),
    ("P-ext", "P-gpi", "pos", "wall"),
    ("P-ext", "PM-TM", "neg", "PM-TM"),
    ("P-ext", "pm-unresolved", "excluded", "pm-unresolved"),
]


@pytest.mark.parametrize("label,d8,cls,stratum", CASES)
def test_class_of_truth_table_p_ext(label, d8, cls, stratum):
    assert labelmap.class_of(label, d8) == cls
    assert labelmap.stratum_of(label, "wall", d8) == stratum


def test_class_of_truth_table():
    # every label x d8_class pair of spec 3.1
    expected = {
        ("P-ext", ""): "pos",
        ("P-ext", "P-gpi"): "pos",
        ("P-ext", "PM-TM"): "neg",
        ("P-ext", "pm-unresolved"): "excluded",
    }
    for label in labelmap.LABELS:
        for d8 in labelmap.D8_CLASSES:
            if label == "P-ext":
                want = expected[(label, d8)]
            else:
                want = {"N-int": "neg", "N-sec": "neg"}.get(label, "excluded")
            assert labelmap.class_of(label, d8) == want, (label, d8)
    assert labelmap.stratum_of("P-ext", "extracellular-only", "") == "extracellular-only"
    assert labelmap.stratum_of("N-sec", "", "") == "N-sec"
    assert labelmap.stratum_of("ambiguous", "", "") == "ambiguous"


def test_class_of_refuses_unknown_values():
    with pytest.raises(ValueError, match="unknown label"):
        labelmap.class_of("P-wall", "")
    with pytest.raises(ValueError, match="unknown d8_class"):
        labelmap.class_of("P-ext", "TM")
    with pytest.raises(ValueError, match="subset"):
        labelmap.stratum_of("P-ext", "", "")
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_scaffold.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'labelmap'` (prototype: `1 error`).

- [ ] **Step 3: Modify `tests/step1_compare/conftest.py`**

```diff
@@ -10,6 +10,7 @@ STEP1_DIR = Path(__file__).resolve().parents[2] / "analysis" / "step1_compare"
 TESTS_DIR = Path(__file__).resolve().parent
 sys.path.insert(0, str(STEP1_DIR))
 sys.path.insert(0, str(STEP1_DIR / "jobs"))
+sys.path.insert(0, str(STEP1_DIR / "phasec"))
 sys.path.insert(0, str(TESTS_DIR))
 
 
@@ -21,6 +22,15 @@ def load_script(name: str):
     return module
 
 
+def load_phasec(name: str):
+    """Import a numbered Phase C script such as 08_build_eval_tables.py as a module."""
+    path = STEP1_DIR / "phasec" / f"{name}.py"
+    spec = importlib.util.spec_from_file_location(f"phasec_{name}", path)
+    module = importlib.util.module_from_spec(spec)
+    spec.loader.exec_module(module)
+    return module
+
+
 @pytest.fixture
 def fixtures_dir() -> Path:
     return TESTS_DIR / "fixtures"
```

- [ ] **Step 4: Write `analysis/step1_compare/phasec/evalio.py`**

```python
"""Shared helpers of the Phase C scripts: STOP error, run JSON, output hashes.

Every Phase C script writes its outputs with `write_outputs`: each output goes to a temp name,
the run JSON is written last and records the SHA-256 of every other output
(`outputs_sha256`), then all files are moved into place (runinfo.atomic_write_all). The next
script calls `require_current` and stops when a file differs from its recorded SHA-256.
"""

import importlib
import json
import platform
from pathlib import Path

import manifest
import runinfo

PHASEC_DIR = Path(__file__).resolve().parent
SEED = 20261001
VARIANTS = ("V-go", "V-kw")


class StopError(ValueError):
    """An input is missing, stale or inconsistent. The script prints STOP and exits 2."""


def out_dir(work) -> Path:
    """Phase C outputs go to $STEP1_WORKDIR/phasec/."""
    return Path(work) / "phasec"


def library_versions() -> dict:
    """Python and library versions for the run JSON (None when a library is missing)."""
    out = {"python": platform.python_version()}
    for name in ("numpy", "scipy", "sklearn"):
        try:
            out[name] = importlib.import_module(name).__version__
        except ImportError:
            out[name] = None
    return out


def read_json(path) -> dict:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise StopError(f"cannot read {path}: {exc}") from exc


def write_outputs(out, writers: dict, run_name: str, log: dict) -> dict:
    """Write the outputs and then the run JSON, which records `outputs_sha256`.

    The run JSON writer reads the temp files `.tmp.<name>` that runinfo.atomic_write_all has
    written before it (writers run in dict order). No file gets its final name before all
    writers have succeeded."""
    out = Path(out)
    names = list(writers)

    def write_run(path):
        log["outputs_sha256"] = {n: manifest.sha256_file(out / f".tmp.{n}") for n in names}
        runinfo.write_json(path, log)

    runinfo.atomic_write_all(out, {**writers, run_name: write_run})
    return log


def require_current(out, run_name: str, names, producer: str) -> dict:
    """Stop unless every file in `names` has the SHA-256 that `run_name` recorded."""
    out = Path(out)
    log = read_json(out / run_name)
    recorded = log.get("outputs_sha256", {})
    for name in names:
        path = out / name
        if not path.exists():
            raise StopError(f"{path} is missing; run {producer}")
        if recorded.get(name) != manifest.sha256_file(path):
            raise StopError(
                f"{name} differs from the SHA-256 in {run_name} (stale or edited input); "
                f"re-run {producer}"
            )
    return log


def float_or_stop(text: str, what: str) -> float:
    """float(text) for a finite number; StopError for empty, NaN or inf."""
    try:
        value = float(text)
    except (TypeError, ValueError) as exc:
        raise StopError(f"{what}: {text!r} is not a number") from exc
    if value != value or value in (float("inf"), float("-inf")):
        raise StopError(f"{what}: {text!r} is not a finite number")
    return value
```

- [ ] **Step 5: Write `analysis/step1_compare/phasec/labelmap.py`**

```python
"""Phase C class mapping (Phase C spec 3.1): one function maps a truth row to its class.

| label     | d8_class        | class    | reporting stratum                   |
|-----------|-----------------|----------|-------------------------------------|
| P-ext     | "" or P-gpi     | pos      | subset (wall, extracellular-only)   |
| P-ext     | PM-TM           | neg      | PM-TM                               |
| P-ext     | pm-unresolved   | excluded | pm-unresolved (ruling C-1)          |
| N-int     | any             | neg      | N-int                               |
| N-sec     | any             | neg      | N-sec                               |
| ambiguous | any             | excluded | ambiguous                           |
| unlabelled| any             | excluded | unlabelled (never a negative)       |

Read `subset`, not `stratum`, for positives: 03_triage_pm.py writes d8_class over stratum.
"""

POS = "pos"
NEG = "neg"
EXCLUDED = "excluded"
CLASSES = (POS, NEG, EXCLUDED)
D8_CLASSES = ("", "P-gpi", "PM-TM", "pm-unresolved")
LABELS = ("P-ext", "N-int", "N-sec", "ambiguous", "unlabelled")
SUBSETS = ("wall", "extracellular-only")


def class_of(label: str, d8_class: str) -> str:
    """Return pos, neg or excluded. ValueError for a label or P-ext d8_class not in the table."""
    if label == "P-ext":
        if d8_class in ("", "P-gpi"):
            return POS
        if d8_class == "PM-TM":
            return NEG
        if d8_class == "pm-unresolved":
            return EXCLUDED
        raise ValueError(f"unknown d8_class {d8_class!r} for a P-ext gene")
    if label in ("N-int", "N-sec"):
        return NEG
    if label in ("ambiguous", "unlabelled"):
        return EXCLUDED
    raise ValueError(f"unknown label {label!r}")


def stratum_of(label: str, subset: str, d8_class: str) -> str:
    """Reporting stratum of one truth row (the table above)."""
    cls = class_of(label, d8_class)
    if label == "P-ext" and cls == POS:
        if subset not in SUBSETS:
            raise ValueError(f"P-ext gene with subset {subset!r}; expected one of {SUBSETS}")
        return subset
    if label == "P-ext":
        return d8_class
    return label
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_scaffold.py tests/step1_compare/test_paths.py -q`
Expected: `18 passed`. Also with `$PY` (no numpy needed): the same count.

- [ ] **Step 7: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/evalio.py analysis/step1_compare/phasec/labelmap.py \
  tests/step1_compare/conftest.py tests/step1_compare/test_phasec_scaffold.py
git commit -m "phasec: package scaffold, import rule, class mapping (spec 3.1)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 2: Dedupe, precedence and the evaluation table (08)

**Files:**
- Create: `analysis/step1_compare/phasec/dedupe.py`, `analysis/step1_compare/phasec/tc_taxon_clades.tsv`, `analysis/step1_compare/phasec/08_build_eval_tables.py`, `tests/step1_compare/phasec_fixture.py` (first part; Task 7 appends)
- Test: `tests/step1_compare/test_phasec_build.py`

**Interfaces:**
- Consumes: `evalio.*`, `labelmap.class_of`, `labelmap.stratum_of` (Task 1); Phase A `truth_table.read_tsv`, `truth_table.write_tsv`, `manifest.sha256_file`, `runinfo.git_commit`, `paths.STEP1_DIR`, `paths.workdir`, `paths.repo_root`; Phase B files `features.tsv.gz`, `features_unique.tsv.gz`, `unique_sequences.tsv.gz`, `features_run.json`, `emb/embedding_run.json`.
- Produces: `dedupe.TABLE_COLUMNS`, `dedupe.LOG_COLUMNS`, `merge_group(h, group, cls) -> dict`, `merge_go(members) -> (rows_by_hash, log)`, `apply_precedence(go_hashes: dict[str, str], tc_rows) -> (tc_rows_by_hash, log)`. Script `08_build_eval_tables.py`: `OUTPUT_NAMES`, `LITERATURE_COLUMNS`, `check_chain(work) -> dict`, `iea_problems(truth_rows) -> list[str]`, `read_tc_clades(path) -> dict[str, str]`, `check_tc_clades(kw_rows, clades)`, `go_members(...)`, `literature_rows(seeds, kw_by_acc)`, `check_features(hashes, features_by_hash)`, `run(work, species_path, seeds_path, clades_path, arguments=())`, `main(argv=None) -> int`. Outputs in `$STEP1_WORKDIR/phasec/`: `eval_table.tsv.gz`, `eval_literature.tsv`, `eval_dedupe_log.tsv`, `eval_sequences.fasta.gz`, `build_run.json`. Test helper `phasec_fixture.make_work(root, seed=7) -> dict` (keys `work`, `species`, `seeds`, `clades`, `sets`, `named`, `seqs`), `build_argv(fx)`, `rewrite_json(path, **changes)`, `gz_lines(path)`, `_sha(path)`, the column constants and `TRUTH_SHA`.

Decisions in this task (state them in review):
- A hash with pos (or neg) members and an excluded member becomes one `excluded` row (`class_and_excluded`; ruling C-14). The real data have 2 such hashes (S288C and *C. albicans*, each a P-ext and an ambiguous gene with one sequence). Cost if wrong: 2 positives fewer in training.
- `homology_only` of a merged row is `no` when any member has direct evidence (12 real hashes have members that differ). The sequence is identical, so the direct evidence of one gene applies to it.
- 08 also writes the literature rows (spec 3.2): positives are rows with an accession, a `uniprot_kw` sequence and `moonlighting` not `YES`; this includes the 9 `hard_negative` rows (owner decision 2026-10-01, ruling C-9 amended; `test_literature_rows` pins LIT2, a `hard_negative` row, as a positive).

- [ ] **Step 1: Write the test fixture `tests/step1_compare/phasec_fixture.py`**

```python
"""A small, realistic Phase A + Phase B work directory for the Phase C tests.

`make_work(root)` writes every file that 08_build_eval_tables.py reads, with the real column
sets, consistent run JSONs (the M3 hash chain holds) and small embedding matrices (8M: 6
columns, 35M: 8 columns). Sources: Scer_SGD and Calb_CGD (train), Spom_PomBase (test_species),
Spom_SCHPO-mod (alternate_file), Afum_ASPFU (test_clade, Eurotiomycetes), Cneo_H99_GOA
(test_clade, Basidiomycota), Umay_MYCMD (undecided, Basidiomycota). Features and embeddings
carry a class signal plus noise, so the models can learn.

Planted cases (the tests rely on them):
- `SHARED_NSEC`, `SHARED_UNRES`, `SHARED_AMBIG`, `SHARED_POS`: T-c rows with the hash of an
  N-sec, a pm-unresolved, an ambiguous and a P-ext gene of Scer_SGD (ruling C-6).
- `TC_DUP`: two T-c accessions with one sequence.
- the alternate source Spom_SCHPO-mod has one gene (`ALT_ONLY`) whose sequence is in no other
  source; it must not enter the table.
- `LONG_1022` and `LONG_1023`: Calb_CGD proteins of 1,022 and 1,023 aa.
- families: groups of genes that share the first 6 residues (the stub MMseqs2 clusters them).
- proteome sets Scer_proteome and Cimm_RS_proteome (members of truth and new sequences).
"""

import gzip
import hashlib
import json
import random
from pathlib import Path

import numpy as np
import seqhash
import truth_table

AA = "ACDEFGHIKLMNPQRSTVWY"
MODELS = {"esm2_t6_8M_UR50D": 6, "esm2_t12_35M_UR50D": 8}
TRUTH_SHA = "5" * 64
SPECIES_COLUMNS = (
    "source_id", "species", "taxon_id", "taxon_filter", "in_clade", "role", "role_note",
    "gaf_file", "fasta_file", "id_mapping",
)  # fmt: skip
SOURCES = [
    # source_id, species, taxon, clade, role, (pos, nint, nsec, pmtm, unres, ambig)
    ("Scer_SGD", "Saccharomyces cerevisiae S288C", "559292", "Saccharomycotina", "train",
     (16, 16, 12, 2, 1, 2)),
    ("Calb_CGD", "Candida albicans SC5314", "237561", "Saccharomycotina", "train",
     (16, 16, 12, 2, 1, 2)),
    ("Spom_PomBase", "Schizosaccharomyces pombe 972h-", "284812", "Taphrinomycotina",
     "test_species", (6, 8, 6, 0, 0, 1)),
    ("Spom_SCHPO-mod", "Schizosaccharomyces pombe 972h-", "284812", "Taphrinomycotina",
     "alternate_file", (0, 0, 0, 0, 0, 0)),
    ("Afum_ASPFU", "Aspergillus fumigatus Af293", "330879", "Eurotiomycetes", "test_clade",
     (6, 6, 6, 0, 1, 0)),
    ("Cneo_H99_GOA", "Cryptococcus neoformans H99", "235443", "Basidiomycota", "test_clade",
     (4, 4, 4, 0, 0, 0)),
    ("Umay_MYCMD", "Ustilago maydis 521", "5270", "Basidiomycota", "undecided",
     (4, 4, 4, 0, 0, 0)),
]  # fmt: skip
TC_TAXA = [("330879", "Aspergillus fumigatus Af293", 5), ("246410", "Coccidioides immitis RS", 3),
           ("237561", "Candida albicans SC5314", 3), ("559292", "Saccharomyces cerevisiae S288C", 3),
           ("498019", "Candidozyma auris B8441", 2)]  # fmt: skip
TC_CLADES = {"330879": "Eurotiomycetes", "246410": "Eurotiomycetes", "237561": "Saccharomycotina",
             "559292": "Saccharomycotina", "498019": "Saccharomycotina"}  # fmt: skip


def _seq(rng: random.Random, n: int, st: float, prefix: str = "") -> str:
    body = []
    for _ in range(n - len(prefix)):
        body.append(rng.choice("ST") if rng.random() < st else rng.choice(AA))
    return prefix + "".join(body)


def _truth_row(src, gene_id, label, subset, d8, homology_only, htp, symbol=""):
    source_id, species, taxon, clade, role, _ = src
    stratum = d8 or subset or label
    return {
        "source_id": source_id, "species": species, "taxon_id": taxon, "in_clade": clade,
        "role": role, "gene_id": gene_id, "symbol": symbol, "synonym1": "", "label": label,
        "subset": subset, "stratum": stratum, "tier": "T-a", "label_no_homology": label,
        "label_experimental": label, "homology_only": homology_only,
        "pm_candidate": "yes" if d8 else "no", "evidence_codes": "IDA,IEA",
        "surface_evidence": "IDA" if label in ("P-ext", "ambiguous") else "",
        "internal_evidence": "HDA" if htp else ("IDA" if label in ("N-int", "ambiguous") else ""),
        "internal_evidence_htp_only": "yes" if htp else "no",
        "secretory_evidence": "IDA" if label == "N-sec" else "", "source_file": "f.gaf",
        "source_sha256": "0" * 64, "source_date": "2026-05-21", "obo_sha256": "1" * 64,
        "d8_class": d8, "d8_reason": "",
    }  # fmt: skip


def _features(rng: random.Random, cls: str, seq: str) -> dict:
    if cls == "pos":
        sp = rng.random() < 0.85
        gpi = rng.choices(["highly_probable", "probable", "weakly", "none"], [4, 1, 1, 4])[0]
    elif cls == "nsec":
        sp = rng.random() < 0.5
        gpi = rng.choices(["highly_probable", "none"], [1, 9])[0]
    else:
        sp = rng.random() < 0.05
        gpi = "none"
    if len(seq) <= 40:
        gpi = "too_short"
    sp_prob = rng.uniform(0.7, 1.0) if sp else rng.uniform(0.0, 0.2)
    gpi_prob = {"highly_probable": 1.0, "probable": 0.7, "weakly": 0.55}.get(gpi, 0.0)
    st = (seq.count("S") + seq.count("T")) / len(seq)
    return {
        "ser_thr_frac": f"{st:.6f}", "sp_prediction": "SP" if sp else "OTHER",
        "sp_prob": f"{sp_prob:.4f}", "sp_other_prob": f"{1 - sp_prob:.4f}",
        "sp_cs_end": "22" if sp else "", "sp_cs_prob": "0.9" if sp else "", "gpi_call": gpi,
        "gpi_prob": f"{gpi_prob}", "gpi_omega": "", "gpi_fpr": "0.5", "gpi_svm": "-1.0",
    }  # fmt: skip


def _gz_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_work(root: Path, seed: int = 7) -> dict:
    """Write the work directory under root/work and the side files under root. Return names."""
    rng = random.Random(seed)
    root = Path(root)
    work = root / "work"
    (work / "phaseb" / "emb").mkdir(parents=True)
    seqs: dict[str, str] = {}  # hash -> sequence
    cls_of: dict[str, str] = {}  # hash -> signal class for features and embeddings
    truth, members = [], []
    named: dict[str, str] = {}
    families = [_seq(rng, 6, 0.0) for _ in range(6)]

    def add_member(set_id, source_id, gene_id, seq, cls):
        h = seqhash.seq_sha256(seq)
        seqs[h] = seq
        cls_of.setdefault(h, cls)
        members.append({"set_id": set_id, "source_id": source_id, "gene_id": gene_id,
                        "seq_sha256": h, "length": str(len(seq))})  # fmt: skip
        return h

    for src in SOURCES:
        source_id = src[0]
        npos, nint, nsec, pmtm, unres, ambig = src[5]
        k = 0
        plan = (
            [("P-ext", "wall" if i % 3 else "extracellular-only", "", "pos") for i in range(npos)]
            + [("N-int", "", "", "nint")] * nint
            + [("N-sec", "", "", "nsec")] * nsec
            + [("P-ext", "wall", "PM-TM", "nsec")] * pmtm
            + [("P-ext", "wall", "pm-unresolved", "pos")] * unres
            + [("ambiguous", "", "", "pos")] * ambig
        )
        for label, subset, d8, cls in plan:
            k += 1
            gene_id = f"{source_id[:4].upper()}{k:04d}"
            st = 0.30 if cls == "pos" else 0.12
            n = rng.randint(80, 400)
            prefix = families[k % 6] if k % 7 == 0 else ""
            if source_id == "Calb_CGD" and k == 1:
                n = 1022
            if source_id == "Calb_CGD" and k == 2:
                n = 1023
            seq = _seq(rng, n, st, prefix)
            hom = "yes" if k % 5 == 0 else "no"
            htp = label == "ambiguous" and k % 2 == 0
            truth.append(_truth_row(src, gene_id, label, subset, d8, hom, htp))
            h = add_member("truth", source_id, gene_id, seq, cls)
            if source_id == "Calb_CGD" and k in (1, 2):
                named[f"LONG_{n}"] = h
            if source_id == "Scer_SGD":
                key = {
                    "N-sec": "SHARED_NSEC",
                    "pm-unresolved": "SHARED_UNRES",
                    "ambiguous": "SHARED_AMBIG",
                    "P-ext": "SHARED_POS",
                }.get(d8 or label)
                if key and key not in named:
                    named[key] = h
            if source_id == "Spom_PomBase":
                truth.append(_truth_row(SOURCES[3], gene_id, label, subset, d8, hom, htp))
                add_member("truth", "Spom_SCHPO-mod", gene_id, seq, cls)
        # one unlabelled gene per source: never in the table
        k += 1
        seq = _seq(rng, 120, 0.1)
        truth.append(_truth_row(src, f"{source_id[:4].upper()}{k:04d}", "unlabelled", "", "",
                                "no", False))  # fmt: skip
        add_member("truth", source_id, truth[-1]["gene_id"], seq, "nint")
    alt_seq = _seq(rng, 150, 0.1)
    truth.append(_truth_row(SOURCES[3], "ALT0001", "N-int", "", "", "no", False))
    named["ALT_ONLY"] = add_member("truth", "Spom_SCHPO-mod", "ALT0001", alt_seq, "nint")

    # T-c rows (keyword tier) and the uniprot_kw members
    kw = []
    j = 0
    for taxon, genome, n in TC_TAXA:
        for _ in range(n):
            j += 1
            seq = _seq(rng, rng.randint(100, 300), 0.30)
            acc = f"Q{j:05d}"
            kw.append({"accession": acc, "gene": "", "genome": genome, "taxon_id": taxon,
                       "length": str(len(seq)), "seq_sha256": add_member(
                           "uniprot_kw", "uniprot_kw", acc, seq, "pos"), "tier": "T-c"})  # fmt: skip
    for key in ("SHARED_NSEC", "SHARED_UNRES", "SHARED_AMBIG", "SHARED_POS"):
        j += 1
        acc = f"Q{j:05d}"
        h = named[key]
        kw.append({"accession": acc, "gene": "", "genome": "Saccharomyces cerevisiae S288C",
                   "taxon_id": "559292", "length": str(len(seqs[h])), "seq_sha256": h,
                   "tier": "T-c"})  # fmt: skip
        add_member("uniprot_kw", "uniprot_kw", acc, seqs[h], cls_of[h])
    dup = _seq(rng, 210, 0.30)
    for acc in ("QDUP01", "QDUP02"):
        h = add_member("uniprot_kw", "uniprot_kw", acc, dup, "pos")
        kw.append({"accession": acc, "gene": "", "genome": "Candidozyma auris B8441",
                   "taxon_id": "498019", "length": str(len(dup)), "seq_sha256": h,
                   "tier": "T-c"})  # fmt: skip
    named["TC_DUP"] = h

    # literature seeds: an adhesin, a hard_negative, a moonlighting row, a row without accession
    lit_seqs = {"P90001": _seq(rng, 300, 0.35), "P90002": _seq(rng, 250, 0.30),
                "P90003": _seq(rng, 500, 0.08)}  # fmt: skip
    for acc, seq in lit_seqs.items():
        add_member("uniprot_kw", "uniprot_kw", acc, seq, "pos" if acc != "P90003" else "nint")
    seeds = root / "seeds.tsv"
    seeds.write_text(
        "# test seeds\n"
        "gene\tuniprot_query\tspecies\torder\tfamily\tclass\tevidence_level\tmoonlighting\t"
        "pmids\tevidence_summary\n"
        "LIT1\taccession:P90001\tCoccidioides immitis\tOnygenales\tf\tadhesin\tE1\tno\t1\ts\n"
        "LIT2\taccession:P90002\tAspergillus fumigatus\tEurotiales\tf\thard_negative\tN1\tno\t-\ts\n"
        "LIT3\taccession:P90003\tHistoplasma capsulatum\tOnygenales\tf\tadhesin\tE1\tYES\t1\ts\n"
        "LIT4\t\tHistoplasma capsulatum\tOnygenales\tf\tadhesin\tE2\tno\t1\ts\n"
    )

    # proteome sets: two truth sequences, one T-c sequence and new sequences
    scer_truth = [m for m in members if m["source_id"] == "Scer_SGD"][:2]
    for i, m in enumerate(scer_truth):
        add_member("Scer_proteome", "Scer_proteome", f"YPR{i:03d}W", seqs[m["seq_sha256"]], "pos")
    add_member("Scer_proteome", "Scer_proteome", "YPR900W", _seq(rng, 200, 0.1), "nint")
    for i in range(3):
        add_member("Cimm_RS_proteome", "Cimm_RS_proteome", f"CIMG_{i:05d}-t26_1-p1",
                   _seq(rng, 220, 0.25 if i == 0 else 0.1), "pos" if i == 0 else "nint")  # fmt: skip
    add_member(
        "Cimm_RS_proteome",
        "Cimm_RS_proteome",
        "CIMG_09999-t26_1-p1",
        seqs[kw[5]["seq_sha256"]],
        "pos",
    )

    # unique sequences, features, embeddings
    hashes = sorted(seqs)
    unique, fu, crow_of = [], [], {}
    crow = 0
    for i, h in enumerate(hashes):
        seq = seqs[h]
        long = len(seq) > 1022
        crow_of[h] = str(crow) if long else ""
        unique.append({"row": str(i), "seq_sha256": h, "length": str(len(seq)),
                       "cterm_row": crow_of[h], "sequence": seq})  # fmt: skip
        fu.append({"row": str(i), "seq_sha256": h, "length": str(len(seq)),
                   "cterm_row": crow_of[h], **_features(rng, cls_of[h], seq)})  # fmt: skip
        crow += long
    row_of = {h: str(i) for i, h in enumerate(hashes)}
    fu_by = {r["seq_sha256"]: r for r in fu}
    truth_by = {(t["source_id"], t["gene_id"]): t for t in truth}
    frows = []
    for m in members:
        t = truth_by.get((m["source_id"], m["gene_id"])) if m["set_id"] == "truth" else None
        tcols = {c: (t[c] if t else "") for c in
                 ("label", "subset", "stratum", "d8_class", "homology_only", "role")}  # fmt: skip
        f = fu_by[m["seq_sha256"]]
        frows.append({**m, **tcols, **{c: f[c] for c in FEATURE_COLUMNS},
                      "emb_row": row_of[m["seq_sha256"]],
                      "emb_cterm_row": crow_of[m["seq_sha256"]]})  # fmt: skip
    pb = work / "phaseb"
    truth_table.write_tsv(pb / "unique_sequences.tsv.gz", UNIQUE_COLUMNS, unique)
    truth_table.write_tsv(pb / "features_unique.tsv.gz", UNIQUE_FEATURE_COLUMNS, fu)
    truth_table.write_tsv(pb / "features.tsv.gz", MEMBER_FEATURE_COLUMNS, frows)
    truth_table.write_tsv(pb / "sequence_members.tsv.gz", MEMBER_COLUMNS, members)
    truth_table.write_tsv(work / "truth_set_triaged.tsv.gz", list(truth[0]), truth)
    truth_table.write_tsv(work / "keyword_tier.tsv.gz", KEYWORD_COLUMNS, kw)
    nprng = np.random.default_rng(seed)
    signal = {"pos": 1.5, "nsec": -0.5, "nint": -1.5}
    emb_models = {}
    n_long = sum(1 for h in hashes if crow_of[h] != "")
    for model, dim in MODELS.items():
        nterm = nprng.normal(size=(len(hashes), dim)).astype(np.float32)
        nterm[:, 0] += np.array([signal[cls_of[h]] for h in hashes], dtype=np.float32)
        cterm = nprng.normal(size=(n_long, dim)).astype(np.float32)
        np.save(pb / "emb" / f"{model}.nterm.npy", nterm)
        np.save(pb / "emb" / f"{model}.cterm.npy", cterm)
        emb_models[model] = {
            w: {"array_sha256": hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest(),
                "dtype": "float32", "shape": list(a.shape)}
            for w, a in (("nterm", nterm), ("cterm", cterm))
        }  # fmt: skip
    unique_sha = _sha(pb / "unique_sequences.tsv.gz")
    _gz_json(pb / "emb" / "embedding_run.json", {"models": emb_models,
             "unique_sequences": len(hashes), "unique_sequences_sha256": unique_sha})  # fmt: skip
    _gz_json(pb / "features_run.json", {
        "all_sources": True, "truth_set_sha256": TRUTH_SHA,
        "input_sha256": {"truth_set_triaged.tsv.gz": _sha(work / "truth_set_triaged.tsv.gz"),
                         "sequence_members.tsv.gz": _sha(pb / "sequence_members.tsv.gz"),
                         "unique_sequences.tsv.gz": unique_sha}})  # fmt: skip
    _gz_json(work / "d8_run.json", {"all_sources": True, "truth_set_sha256": TRUTH_SHA})
    _gz_json(work / "keyword_tier_run.json", {"all_sources": True, "truth_set_sha256": TRUTH_SHA})
    species = root / "species.tsv"
    truth_table.write_tsv(species, SPECIES_COLUMNS, [
        {"source_id": s[0], "species": s[1], "taxon_id": s[2], "taxon_filter": "",
         "in_clade": s[3], "role": s[4], "role_note": "", "gaf_file": "", "fasta_file": "",
         "id_mapping": ""} for s in SOURCES])  # fmt: skip
    clades = root / "tc_taxon_clades.tsv"
    truth_table.write_tsv(clades, ("taxon_id", "organism", "clade"), [
        {"taxon_id": t, "organism": g, "clade": TC_CLADES[t]} for t, g, _ in TC_TAXA])  # fmt: skip
    sets = root / "sequence_sets.tsv"
    sets.write_text("set_id\tkind\tlocation\tnote\ntruth\ttruth\ttruth_sequences.tsv.gz\t\n"
                    "uniprot_kw\tkeyword\tkeyword_sequences.fasta.gz\t\n"
                    "Scer_proteome\tdownload\tx.fasta.gz\t\n"
                    "Cimm_RS_proteome\tsite\tcocci:x.fasta\t\n")  # fmt: skip
    return {"work": work, "species": species, "seeds": seeds, "clades": clades, "sets": sets,
            "named": named, "seqs": seqs}  # fmt: skip


def build_argv(fx: dict) -> list[str]:
    return ["--work-dir", str(fx["work"]), "--species", str(fx["species"]),
            "--seeds", str(fx["seeds"]), "--tc-clades", str(fx["clades"])]  # fmt: skip


def rewrite_json(path: Path, **changes) -> None:
    obj = json.loads(Path(path).read_text())
    for key, value in changes.items():
        obj[key] = value
    Path(path).write_text(json.dumps(obj))


def gz_lines(path: Path) -> list[str]:
    return gzip.decompress(Path(path).read_bytes()).decode().splitlines()


FEATURE_COLUMNS = (
    "ser_thr_frac", "sp_prediction", "sp_prob", "sp_other_prob", "sp_cs_end", "sp_cs_prob",
    "gpi_call", "gpi_prob", "gpi_omega", "gpi_fpr", "gpi_svm",
)  # fmt: skip
UNIQUE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", "sequence")
UNIQUE_FEATURE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", *FEATURE_COLUMNS)
MEMBER_COLUMNS = ("set_id", "source_id", "gene_id", "seq_sha256", "length")
MEMBER_FEATURE_COLUMNS = (
    "set_id", "source_id", "gene_id", "seq_sha256", "length", "label", "subset", "stratum",
    "d8_class", "homology_only", "role", *FEATURE_COLUMNS, "emb_row", "emb_cterm_row",
)  # fmt: skip
KEYWORD_COLUMNS = ("accession", "gene", "genome", "taxon_id", "length", "seq_sha256", "tier")
```

- [ ] **Step 2: Write the failing test**

`tests/step1_compare/test_phasec_build.py`:

```python
"""08_build_eval_tables.py: class mapping on real columns, dedupe, precedence (C-6), checks."""

import json

import pytest

np = pytest.importorskip("numpy")

import dedupe  # noqa: E402
import phasec_fixture as pf  # noqa: E402
import truth_table  # noqa: E402
from conftest import load_phasec, load_script  # noqa: E402


@pytest.fixture
def fx(tmp_path):
    return pf.make_work(tmp_path)


def _run(fx):
    return load_phasec("08_build_eval_tables").main(pf.build_argv(fx))


def _table(fx):
    return {
        r["seq_sha256"]: r
        for r in truth_table.read_tsv(fx["work"] / "phasec" / "eval_table.tsv.gz")
    }


def test_fixture_columns_equal_the_phase_b_constants():
    import keyword_tier
    import seqsets

    f07 = load_script("07_build_features")
    assert pf.MEMBER_FEATURE_COLUMNS == f07.MEMBER_FEATURE_COLUMNS
    assert pf.UNIQUE_FEATURE_COLUMNS == f07.UNIQUE_FEATURE_COLUMNS
    assert pf.UNIQUE_COLUMNS == seqsets.UNIQUE_COLUMNS
    assert pf.MEMBER_COLUMNS == seqsets.MEMBER_COLUMNS
    assert pf.KEYWORD_COLUMNS == keyword_tier.KEYWORD_COLUMNS


def test_build_writes_all_outputs_with_recorded_hashes(fx):
    import manifest

    assert _run(fx) == 0
    out = fx["work"] / "phasec"
    log = json.loads((out / "build_run.json").read_text())
    names = load_phasec("08_build_eval_tables").OUTPUT_NAMES
    assert sorted(p.name for p in out.iterdir()) == sorted(names)
    for name in names[:-1]:
        assert log["outputs_sha256"][name] == manifest.sha256_file(out / name)
    assert log["all_sources"] is True and log["truth_set_sha256"] == pf.TRUTH_SHA


def test_alternate_files_excluded(fx):
    assert _run(fx) == 0
    table = _table(fx)
    assert fx["named"]["ALT_ONLY"] not in table
    assert not any("Spom_SCHPO-mod" in r["source_ids"] for r in table.values())
    # the Spom genes are in the table once, from Spom_PomBase only
    spom = [r for r in table.values() if r["source_ids"] == "Spom_PomBase"]
    assert len(spom) == 6 + 8 + 6 + 1


def test_unlabelled_genes_are_not_in_the_table(fx):
    assert _run(fx) == 0
    assert not any("unlabelled" in r["label"] for r in _table(fx).values())


def test_tc_row_yields_to_go_label(fx):
    assert _run(fx) == 0
    table = _table(fx)
    n = fx["named"]
    assert table[n["SHARED_NSEC"]]["class"] == "neg" and table[n["SHARED_NSEC"]]["origin"] == "go"
    assert table[n["SHARED_UNRES"]]["class"] == "excluded"
    assert table[n["SHARED_UNRES"]]["stratum"] == "pm-unresolved"
    assert table[n["SHARED_AMBIG"]]["class"] == "excluded"
    assert table[n["SHARED_AMBIG"]]["stratum"] == "ambiguous"
    assert table[n["SHARED_POS"]]["origin"] == "go"
    log = truth_table.read_tsv(fx["work"] / "phasec" / "eval_dedupe_log.tsv")
    wins = [r for r in log if r["reason"] == "go_label_wins"]
    assert {r["seq_sha256"] for r in wins} == {
        n[k] for k in ("SHARED_NSEC", "SHARED_UNRES", "SHARED_AMBIG", "SHARED_POS")
    }
    run = json.loads((fx["work"] / "phasec" / "build_run.json").read_text())
    assert run["tc_rows_dropped_by_go_class"] == {"excluded": 2, "neg": 1, "pos": 1}
    tc = [r for r in table.values() if r["origin"] == "tc"]
    assert len(tc) == 16 + 1  # 16 T-c rows with own sequences, one row for the duplicate pair
    assert table[n["TC_DUP"]]["gene_ids"] == "QDUP01,QDUP02"
    assert all(r["class"] == "pos" for r in tc)


def test_dedupe_drops_a_hash_present_in_both_classes():
    base = {"label": "P-ext", "subset": "wall", "stratum": "wall", "d8_class": "",
            "homology_only": "no", "internal_evidence_htp_only": "no", "species": "s",
            "role": "train", "clade": "c", "taxon_id": "1", "length": "100", "emb_row": "0",
            "emb_cterm_row": ""}  # fmt: skip
    pos = {**base, "seq_sha256": "a", "class": "pos", "source_id": "S", "gene_id": "g1"}
    neg = {**base, "seq_sha256": "a", "class": "neg", "source_id": "S", "gene_id": "g2",
           "label": "N-sec", "subset": "", "stratum": "N-sec"}  # fmt: skip
    keep = {**base, "seq_sha256": "b", "class": "pos", "source_id": "S", "gene_id": "g3"}
    keep2 = {**keep, "source_id": "T", "gene_id": "g4", "homology_only": "yes"}
    rows, log = dedupe.merge_go([pos, neg, keep, keep2])
    assert set(rows) == {"b"}
    assert rows["b"]["source_ids"] == "S,T" and rows["b"]["homology_only"] == "no"
    assert [(r["gene_id"], r["reason"]) for r in log] == [
        ("g1", "both_classes"),
        ("g2", "both_classes"),
    ]
    excl = {
        **keep,
        "gene_id": "g5",
        "class": "excluded",
        "label": "ambiguous",
        "stratum": "ambiguous",
    }
    rows, log = dedupe.merge_go([keep, excl])
    assert rows["b"]["class"] == "excluded" and {r["reason"] for r in log} == {"class_and_excluded"}


def test_truth_set_has_no_iea(fx):
    path = fx["work"] / "truth_set_triaged.tsv.gz"
    rows = truth_table.read_tsv(path)
    row = next(r for r in rows if r["label"] == "N-sec")
    row["evidence_codes"] = "IEA"
    truth_table.write_tsv(path, list(rows[0]), rows)
    pf.rewrite_json(fx["work"] / "phaseb" / "features_run.json", input_sha256={
        **json.loads((fx["work"] / "phaseb" / "features_run.json").read_text())["input_sha256"],
        "truth_set_triaged.tsv.gz": pf._sha(path)})  # fmt: skip
    m = load_phasec("08_build_eval_tables")
    assert m.iea_problems(rows) == [
        f"{row['source_id']}:{row['gene_id']} has label N-sec with IEA evidence only"
    ]
    assert _run(fx) == 2
    assert not (fx["work"] / "phasec").exists()


def test_iea_in_an_evidence_column_is_found():
    m = load_phasec("08_build_eval_tables")
    row = {"source_id": "S", "gene_id": "g", "label": "P-ext", "evidence_codes": "IDA,IEA",
           "surface_evidence": "IDA,IEA", "internal_evidence": "", "secretory_evidence": ""}  # fmt: skip
    assert m.iea_problems([row]) == ["S:g has IEA in surface_evidence"]


def test_tc_taxon_table_is_complete(fx, capsys):
    rows = truth_table.read_tsv(fx["clades"])
    truth_table.write_tsv(
        fx["clades"], list(rows[0]), [r for r in rows if r["taxon_id"] != "246410"]
    )
    assert _run(fx) == 2
    assert "tc_taxon_clades.tsv has no row for taxon_id 246410" in capsys.readouterr().err


def test_committed_tc_taxon_table_covers_the_phase_a_taxa():
    # The seven taxon_ids of the Phase A keyword_tier.tsv.gz (3,092 rows, counted 2026-10-01).
    m = load_phasec("08_build_eval_tables")
    import evalio

    clades = m.read_tc_clades(evalio.PHASEC_DIR / "tc_taxon_clades.tsv")
    assert clades == {
        "237561": "Saccharomycotina", "246410": "Eurotiomycetes", "284593": "Saccharomycotina",
        "330879": "Eurotiomycetes", "443226": "Eurotiomycetes", "498019": "Saccharomycotina",
        "559292": "Saccharomycotina",
    }  # fmt: skip
    species = truth_table.read_tsv(m.paths.STEP1_DIR / "species.tsv")
    assert set(clades.values()) <= {s["in_clade"] for s in species} | {"other"}


def test_stale_input_stops(fx, capsys):
    # M3: features and embeddings from different unique sequence sets
    pf.rewrite_json(
        fx["work"] / "phaseb" / "emb" / "embedding_run.json", unique_sequences_sha256="0" * 64
    )
    assert _run(fx) == 2
    assert "differs from embedding_run.json unique_sequences_sha256" in capsys.readouterr().err
    assert not (fx["work"] / "phasec").exists()


def test_changed_unique_sequences_file_stops(fx, capsys):
    path = fx["work"] / "phaseb" / "unique_sequences.tsv.gz"
    path.write_bytes(path.read_bytes() + b"\n")
    assert _run(fx) == 2
    assert "differs from the file 07 read" in capsys.readouterr().err


def test_truth_hash_mismatch_stops(fx, capsys):
    pf.rewrite_json(fx["work"] / "keyword_tier_run.json", truth_set_sha256="9" * 64)
    assert _run(fx) == 2
    assert "truth_set_sha256 differs" in capsys.readouterr().err


def test_literature_rows(fx):
    assert _run(fx) == 0
    lit = {
        r["gene"]: r for r in truth_table.read_tsv(fx["work"] / "phasec" / "eval_literature.tsv")
    }
    assert set(lit) == {"LIT1", "LIT2", "LIT3"}
    assert [lit[g]["literature_positive"] for g in ("LIT1", "LIT2", "LIT3")] == ["yes", "yes", "no"]
    run = json.loads((fx["work"] / "phasec" / "build_run.json").read_text())
    assert run["literature_no_accession"] == ["LIT4"] and run["literature_positives"] == 2


def test_fasta_holds_table_and_literature_sequences(fx):
    assert _run(fx) == 0
    lines = pf.gz_lines(fx["work"] / "phasec" / "eval_sequences.fasta.gz")
    heads = [x[1:] for x in lines if x.startswith(">")]
    table = _table(fx)
    lit = truth_table.read_tsv(fx["work"] / "phasec" / "eval_literature.tsv")
    assert heads == sorted(set(table) | {r["seq_sha256"] for r in lit})
    assert lines[lines.index(">" + heads[0]) + 1] == fx["seqs"][heads[0]]


def test_missing_feature_value_stops(fx, capsys):
    # Review Focus 3: a table hash with an empty sp_prob must stop, not become NaN
    path = fx["work"] / "phaseb" / "features_unique.tsv.gz"
    rows = truth_table.read_tsv(path)
    target = fx["named"]["SHARED_POS"]
    for r in rows:
        if r["seq_sha256"] == target:
            r["sp_prob"] = ""
    truth_table.write_tsv(path, pf.UNIQUE_FEATURE_COLUMNS, rows)
    assert _run(fx) == 2
    assert f"{target} sp_prob: '' is not a number" in capsys.readouterr().err


def test_identical_sequence_in_two_sources_is_one_row(fx):
    # Review Focus 4: one row; both sources listed
    path = fx["work"] / "phaseb" / "features.tsv.gz"
    rows = truth_table.read_tsv(path)
    afum = next(r for r in rows if r["source_id"] == "Afum_ASPFU" and r["label"] == "N-sec")
    h99 = next(r for r in rows if r["source_id"] == "Cneo_H99_GOA" and r["label"] == "N-sec")
    for c in ("seq_sha256", "length", "emb_row", "emb_cterm_row"):
        h99[c] = afum[c]
    truth_table.write_tsv(path, pf.MEMBER_FEATURE_COLUMNS, rows)
    assert _run(fx) == 0
    row = _table(fx)[afum["seq_sha256"]]
    assert row["source_ids"] == "Afum_ASPFU,Cneo_H99_GOA"
    assert row["clades"] == "Basidiomycota,Eurotiomycetes"
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_build.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'dedupe'`.

- [ ] **Step 4: Write `analysis/step1_compare/phasec/dedupe.py`**

```python
"""Dedupe and precedence for the Phase C table (Phase C spec 3.4 item 2, ruling C-6).

This module does not import the surface_glyco dedupe function; it has its own rules and tests.

GO members (truth genes) are grouped by seq_sha256:
- one class in the group: one table row; the members' values are merged (see merge_group).
- pos and neg in the group: the hash is dropped (parent spec 4 step 4); every member is logged
  with reason `both_classes`.
- pos or neg together with excluded: the row becomes excluded and every member is logged with
  reason `class_and_excluded` (plan decision: an excluded label on the same sequence makes the
  label uncertain; 2 hashes in the Phase A data).

T-c rows: a T-c row whose hash has any GO member (pos, neg or excluded, and also a hash dropped
as `both_classes`) is dropped with reason `go_label_wins`: the GO label wins (ruling C-6). The
other T-c rows are grouped by hash into one positive row each.
"""

from labelmap import EXCLUDED, NEG, POS

TABLE_COLUMNS = (
    "seq_sha256",
    "origin",
    "class",
    "label",
    "subset",
    "stratum",
    "d8_class",
    "homology_only",
    "internal_evidence_htp_only",
    "source_ids",
    "gene_ids",
    "species",
    "roles",
    "clades",
    "taxon_ids",
    "length",
    "emb_row",
    "emb_cterm_row",
)
LOG_COLUMNS = ("origin", "source_id", "gene_id", "seq_sha256", "class", "reason", "detail")


def _join(values) -> str:
    return ",".join(sorted({v for v in values if v != ""}))


def _check_same(group: list[dict], column: str, h: str) -> str:
    values = {m[column] for m in group}
    if len(values) != 1:
        raise ValueError(f"hash {h}: members disagree on {column} ({sorted(values)})")
    return values.pop()


def merge_group(h: str, group: list[dict], cls: str) -> dict:
    """One table row from GO members that share seq_sha256 `h`.

    homology_only is `no` when any member has direct evidence (the same sequence carries the
    label without homology codes); internal_evidence_htp_only is `yes` when any member says so;
    list columns hold the sorted unique values, comma separated."""
    return {
        "seq_sha256": h,
        "origin": "go",
        "class": cls,
        "label": _join(m["label"] for m in group),
        "subset": _join(m["subset"] for m in group),
        "stratum": _join(m["stratum"] for m in group),
        "d8_class": _join(m["d8_class"] for m in group),
        "homology_only": "no" if any(m["homology_only"] == "no" for m in group) else "yes",
        "internal_evidence_htp_only": "yes"
        if any(m["internal_evidence_htp_only"] == "yes" for m in group)
        else "no",
        "source_ids": _join(m["source_id"] for m in group),
        "gene_ids": _join(m["gene_id"] for m in group),
        "species": _join(m["species"] for m in group),
        "roles": _join(m["role"] for m in group),
        "clades": _join(m["clade"] for m in group),
        "taxon_ids": _join(m["taxon_id"] for m in group),
        "length": _check_same(group, "length", h),
        "emb_row": _check_same(group, "emb_row", h),
        "emb_cterm_row": _check_same(group, "emb_cterm_row", h),
    }


def merge_go(members: list[dict]) -> tuple[dict[str, dict], list[dict]]:
    """Group GO members by hash. Return (table rows by hash, log rows)."""
    groups: dict[str, list[dict]] = {}
    for m in members:
        groups.setdefault(m["seq_sha256"], []).append(m)
    rows, log = {}, []
    for h in sorted(groups):
        group = groups[h]
        classes = {m["class"] for m in group}
        if {POS, NEG} <= classes:
            for m in group:
                log.append(_log("go", m, m["class"], "both_classes", _join(classes)))
            continue
        if EXCLUDED in classes and len(classes) > 1:
            for m in group:
                log.append(_log("go", m, m["class"], "class_and_excluded", _join(classes)))
            rows[h] = merge_group(h, group, EXCLUDED)
            continue
        rows[h] = merge_group(h, group, classes.pop())
    return rows, log


def _log(origin, m, cls, reason, detail) -> dict:
    return {
        "origin": origin,
        "source_id": m.get("source_id", "T-c"),
        "gene_id": m.get("gene_id", m.get("accession", "")),
        "seq_sha256": m["seq_sha256"],
        "class": cls,
        "reason": reason,
        "detail": detail,
    }


def apply_precedence(
    go_hashes: dict[str, str], tc_rows: list[dict]
) -> tuple[dict[str, dict], list[dict]]:
    """Drop T-c rows whose hash has any GO member; merge the rest by hash.

    `go_hashes` maps every hash of a GO member (kept or dropped) to its class detail (for the
    log). Each T-c row needs seq_sha256, accession, taxon_id, clade, length, emb_row,
    emb_cterm_row. Return (T-c table rows by hash, log rows)."""
    kept: dict[str, list[dict]] = {}
    log = []
    for r in tc_rows:
        h = r["seq_sha256"]
        if h in go_hashes:
            log.append(_log("tc", r, POS, "go_label_wins", go_hashes[h]))
            continue
        kept.setdefault(h, []).append(r)
    rows = {}
    for h in sorted(kept):
        group = kept[h]
        rows[h] = {
            "seq_sha256": h,
            "origin": "tc",
            "class": POS,
            "label": "T-c",
            "subset": "",
            "stratum": "T-c",
            "d8_class": "",
            "homology_only": "",
            "internal_evidence_htp_only": "",
            "source_ids": "T-c",
            "gene_ids": _join(r["accession"] for r in group),
            "species": _join(r["genome"] for r in group),
            "roles": "tc",
            "clades": _join(r["clade"] for r in group),
            "taxon_ids": _join(r["taxon_id"] for r in group),
            "length": _check_same(group, "length", h),
            "emb_row": _check_same(group, "emb_row", h),
            "emb_cterm_row": _check_same(group, "emb_cterm_row", h),
        }
    return rows, log
```

- [ ] **Step 5: Write `analysis/step1_compare/phasec/tc_taxon_clades.tsv`**

TAB separated; one row per `taxon_id` of the Phase A `keyword_tier.tsv.gz` (7 taxa, counted 2026-10-01). Clade names are the `in_clade` values of `species.tsv`.

```
taxon_id	organism	clade
237561	Candida albicans SC5314	Saccharomycotina
246410	Coccidioides immitis RS	Eurotiomycetes
284593	Nakaseomyces glabratus CBS138	Saccharomycotina
330879	Aspergillus fumigatus Af293	Eurotiomycetes
443226	Coccidioides posadasii C735 delta SOWgp	Eurotiomycetes
498019	Candidozyma auris B8441	Saccharomycotina
559292	Saccharomyces cerevisiae S288C	Saccharomycotina
```

- [ ] **Step 6: Write `analysis/step1_compare/phasec/08_build_eval_tables.py`**

```python
#!/usr/bin/env python3
"""Phase C step 1: the evaluation table (Phase C spec 3.1, 3.2, 3.4 items 1, 2 and 5c).

Reads from $STEP1_WORKDIR: truth_set_triaged.tsv.gz (03), keyword_tier.tsv.gz (04),
phaseb/features.tsv.gz, phaseb/features_unique.tsv.gz (07), phaseb/unique_sequences.tsv.gz
(05), and the run JSONs d8_run.json, keyword_tier_run.json, phaseb/features_run.json,
phaseb/emb/embedding_run.json. Reads species.tsv, data/curated/adhesins/eurotiomycetes_seeds.tsv
and phasec/tc_taxon_clades.tsv. Writes to $STEP1_WORKDIR/phasec/:

  eval_table.tsv.gz        one row per unique sequence: GO truth genes of non-alternate sources
                           with a label other than `unlabelled`, and the T-c rows that survive
                           the precedence rule (ruling C-6)
  eval_literature.tsv      one row per literature seed with an accession (spec 3.2)
  eval_dedupe_log.tsv      one row per dropped or reclassified member
  eval_sequences.fasta.gz  the sequences of the table and the literature rows (MMseqs2 input)
  build_run.json           counts, input hashes, all_sources, git commit, library versions

STOP (exit 2, no output): a run JSON without all_sources: true; features_run.json
input_sha256["unique_sequences.tsv.gz"] differs from embedding_run.json
unique_sequences_sha256 or from the current file (Phase B review item M3); the truth set hash
differs between d8_run.json, features_run.json and keyword_tier_run.json; the current
truth_set_triaged.tsv.gz differs from the file 07 read; a labelled truth row whose only evidence
is IEA, or IEA in an evidence column; a T-c taxon_id without a row in tc_taxon_clades.tsv; a
source without a row in species.tsv or with another role there; a table or literature hash
without a valid feature row (empty or non-finite sp_prob, gpi_prob, ser_thr_frac; unknown
gpi_call) or without a sequence.
"""

import argparse
import gzip
import io
import sys
from collections import Counter
from pathlib import Path

import dedupe
import evalio
import labelmap
import manifest
import paths
import runinfo
import truth_table

OUTPUT_NAMES = (
    "eval_table.tsv.gz",
    "eval_literature.tsv",
    "eval_dedupe_log.tsv",
    "eval_sequences.fasta.gz",
    "build_run.json",
)
LITERATURE_COLUMNS = (
    "accession",
    "gene",
    "lit_class",
    "moonlighting",
    "species",
    "order",
    "seq_sha256",
    "length",
    "emb_row",
    "emb_cterm_row",
    "literature_positive",
)
GPI_CALLS = ("highly_probable", "probable", "weakly", "none", "too_short")
EVIDENCE_COLUMNS = ("surface_evidence", "internal_evidence", "secretory_evidence")
DEFAULT_SEEDS = Path("data/curated/adhesins/eurotiomycetes_seeds.tsv")


def check_chain(work: Path) -> dict:
    """Run JSON checks. Return the hashes for build_run.json."""
    work = Path(work)
    fr = evalio.read_json(work / "phaseb" / "features_run.json")
    er = evalio.read_json(work / "phaseb" / "emb" / "embedding_run.json")
    d8 = evalio.read_json(work / "d8_run.json")
    kt = evalio.read_json(work / "keyword_tier_run.json")
    for name, log in (
        ("features_run.json", fr),
        ("d8_run.json", d8),
        ("keyword_tier_run.json", kt),
    ):
        if log.get("all_sources") is not True:
            raise evalio.StopError(f"{name} does not say all_sources: true")
    unique = work / "phaseb" / "unique_sequences.tsv.gz"
    recorded = fr.get("input_sha256", {}).get("unique_sequences.tsv.gz")
    if recorded != er.get("unique_sequences_sha256"):
        raise evalio.StopError(
            f"features_run.json input_sha256['unique_sequences.tsv.gz'] ({recorded}) differs "
            f"from embedding_run.json unique_sequences_sha256 "
            f"({er.get('unique_sequences_sha256')}); features and embeddings come from "
            "different unique sequence sets; re-run 07 and the assembly"
        )
    current = manifest.sha256_file(unique)
    if recorded != current:
        raise evalio.StopError(
            f"phaseb/unique_sequences.tsv.gz ({current}) differs from the file 07 read "
            f"({recorded}); re-run 07 and the assembly"
        )
    truth_hashes = {
        "d8_run.json": d8.get("truth_set_sha256"),
        "features_run.json": fr.get("truth_set_sha256"),
        "keyword_tier_run.json": kt.get("truth_set_sha256"),
    }
    if len(set(truth_hashes.values())) != 1 or None in truth_hashes.values():
        raise evalio.StopError(f"truth_set_sha256 differs between run JSONs: {truth_hashes}")
    triaged = work / "truth_set_triaged.tsv.gz"
    triaged_sha = manifest.sha256_file(triaged)
    if fr.get("input_sha256", {}).get("truth_set_triaged.tsv.gz") != triaged_sha:
        raise evalio.StopError(
            "truth_set_triaged.tsv.gz differs from the file that 07 read; re-run 07"
        )
    return {
        "truth_set_sha256": truth_hashes["d8_run.json"],
        "unique_sequences_sha256": current,
        "truth_set_triaged_sha256": triaged_sha,
    }


def iea_problems(truth_rows) -> list[str]:
    """Labelled rows whose evidence is IEA only, and IEA in a label evidence column."""
    bad = []
    for r in truth_rows:
        key = f"{r['source_id']}:{r['gene_id']}"
        if r["label"] == "unlabelled":
            continue
        codes = set(r["evidence_codes"].split(",")) - {""}
        if not codes - {"IEA"}:
            bad.append(f"{key} has label {r['label']} with IEA evidence only")
        for col in EVIDENCE_COLUMNS:
            if "IEA" in r[col].split(","):
                bad.append(f"{key} has IEA in {col}")
    return bad


def read_tc_clades(path: Path) -> dict[str, str]:
    rows = truth_table.read_tsv(path)
    out = {}
    for r in rows:
        if r["taxon_id"] in out:
            raise evalio.StopError(f"{path}: taxon_id {r['taxon_id']} occurs twice")
        if not r["clade"]:
            raise evalio.StopError(f"{path}: taxon_id {r['taxon_id']} has an empty clade")
        out[r["taxon_id"]] = r["clade"]
    return out


def check_tc_clades(kw_rows, clades: dict[str, str]) -> None:
    missing = sorted({r["taxon_id"] for r in kw_rows} - set(clades))
    if missing:
        raise evalio.StopError(
            f"tc_taxon_clades.tsv has no row for taxon_id {', '.join(missing)}; add the clade "
            "(or `other`) for every taxon_id of keyword_tier.tsv.gz"
        )


def go_members(truth_rows, feature_rows, species_rows) -> tuple[list[dict], Counter]:
    """GO members of non-alternate sources with a label other than unlabelled.

    Return (members, count of labelled genes without a sequence per source)."""
    species = {r["source_id"]: r for r in species_rows}
    truth = {(r["source_id"], r["gene_id"]): r for r in truth_rows}
    for r in truth_rows:
        s = species.get(r["source_id"])
        if s is None:
            raise evalio.StopError(f"source {r['source_id']} has no row in species.tsv")
        if s["role"] != r["role"]:
            raise evalio.StopError(
                f"source {r['source_id']}: role {r['role']!r} in the truth set, "
                f"{s['role']!r} in species.tsv"
            )
    members, seen = [], set()
    for f in feature_rows:
        if f["set_id"] != "truth":
            continue
        key = (f["source_id"], f["gene_id"])
        t = truth.get(key)
        if t is None:
            raise evalio.StopError(f"truth member {key[0]}:{key[1]} has no truth row")
        seen.add(key)
        if t["role"] == "alternate_file" or t["label"] == "unlabelled":
            continue
        members.append(
            {
                "seq_sha256": f["seq_sha256"],
                "class": labelmap.class_of(t["label"], t["d8_class"]),
                "label": t["label"],
                "subset": t["subset"],
                "stratum": labelmap.stratum_of(t["label"], t["subset"], t["d8_class"]),
                "d8_class": t["d8_class"],
                "homology_only": t["homology_only"],
                "internal_evidence_htp_only": t["internal_evidence_htp_only"],
                "source_id": t["source_id"],
                "gene_id": t["gene_id"],
                "species": t["species"],
                "role": t["role"],
                "clade": species[t["source_id"]]["in_clade"],
                "taxon_id": t["taxon_id"],
                "length": f["length"],
                "emb_row": f["emb_row"],
                "emb_cterm_row": f["emb_cterm_row"],
            }
        )
    unmatched = Counter(
        r["source_id"]
        for key, r in truth.items()
        if key not in seen and r["role"] != "alternate_file" and r["label"] != "unlabelled"
    )
    return members, unmatched


def kw_features(feature_rows) -> dict[str, dict]:
    """uniprot_kw members by accession (gene_id)."""
    return {f["gene_id"]: f for f in feature_rows if f["set_id"] == "uniprot_kw"}


def tc_members(kw_rows, kw_by_acc: dict, clades: dict) -> list[dict]:
    out = []
    for r in kw_rows:
        f = kw_by_acc.get(r["accession"])
        if f is None or f["seq_sha256"] != r["seq_sha256"]:
            raise evalio.StopError(
                f"T-c row {r['accession']} has no uniprot_kw feature row with its seq_sha256; "
                "re-run 05 and 07 after 04"
            )
        out.append(
            {
                **{k: r[k] for k in ("accession", "genome", "taxon_id", "seq_sha256")},
                "clade": clades[r["taxon_id"]],
                "length": f["length"],
                "emb_row": f["emb_row"],
                "emb_cterm_row": f["emb_cterm_row"],
            }
        )
    return out


def read_seeds(path: Path) -> list[dict]:
    import csv

    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def literature_rows(seeds, kw_by_acc: dict) -> tuple[list[dict], list[str]]:
    """Rows with an accession. Positive: a sequence in the keyword set and moonlighting not YES.

    Return (rows, genes without an accession)."""
    rows, no_accession = [], []
    for s in seeds:
        query = (s.get("uniprot_query") or "").strip()
        if not query.startswith("accession:"):
            no_accession.append(s["gene"])
            continue
        acc = query.split(":", 1)[1].strip()
        f = kw_by_acc.get(acc)
        moon = (s.get("moonlighting") or "").strip().upper() == "YES"
        rows.append(
            {
                "accession": acc,
                "gene": s["gene"],
                "lit_class": s["class"],
                "moonlighting": s["moonlighting"],
                "species": s["species"],
                "order": s["order"],
                "seq_sha256": f["seq_sha256"] if f else "",
                "length": f["length"] if f else "",
                "emb_row": f["emb_row"] if f else "",
                "emb_cterm_row": f["emb_cterm_row"] if f else "",
                "literature_positive": "yes" if f and not moon else "no",
            }
        )
    return rows, no_accession


def check_features(hashes, features_by_hash: dict) -> None:
    for h in sorted(hashes):
        f = features_by_hash.get(h)
        if f is None:
            raise evalio.StopError(f"hash {h} has no row in features_unique.tsv.gz")
        for col in ("sp_prob", "gpi_prob", "ser_thr_frac"):
            evalio.float_or_stop(f[col], f"features_unique.tsv.gz {h} {col}")
        if not f["sp_prediction"]:
            raise evalio.StopError(f"features_unique.tsv.gz {h}: empty sp_prediction")
        if f["gpi_call"] not in GPI_CALLS:
            raise evalio.StopError(f"features_unique.tsv.gz {h}: gpi_call {f['gpi_call']!r}")


def fasta_bytes(hashes, seq_by_hash: dict) -> bytes:
    buf = io.StringIO()
    for h in sorted(hashes):
        seq = seq_by_hash.get(h)
        if seq is None:
            raise evalio.StopError(f"hash {h} has no row in unique_sequences.tsv.gz")
        buf.write(f">{h}\n{seq}\n")
    raw = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(buf.getvalue().encode())
    return raw.getvalue()


def run(work: Path, species_path: Path, seeds_path: Path, clades_path: Path, arguments=()):
    work = Path(work)
    chain = check_chain(work)
    truth_rows = truth_table.read_tsv(work / "truth_set_triaged.tsv.gz")
    problems = iea_problems(truth_rows)
    if problems:
        raise evalio.StopError(f"{len(problems)} IEA problems, first: {problems[0]}")
    features = truth_table.read_tsv(work / "phaseb" / "features.tsv.gz")
    fu = {
        r["seq_sha256"]: r for r in truth_table.read_tsv(work / "phaseb" / "features_unique.tsv.gz")
    }
    kw_rows = truth_table.read_tsv(work / "keyword_tier.tsv.gz")
    clades = read_tc_clades(clades_path)
    check_tc_clades(kw_rows, clades)
    species_rows = truth_table.read_tsv(species_path)
    members, unmatched = go_members(truth_rows, features, species_rows)
    go_rows, go_log = dedupe.merge_go(members)
    go_hashes: dict[str, set] = {}
    for m in members:
        go_hashes.setdefault(m["seq_sha256"], set()).add(m["class"])
    go_detail = {h: ",".join(sorted(v)) for h, v in go_hashes.items()}
    kw_by_acc = kw_features(features)
    tc_rows, tc_log = dedupe.apply_precedence(go_detail, tc_members(kw_rows, kw_by_acc, clades))
    table = sorted([*go_rows.values(), *tc_rows.values()], key=lambda r: r["seq_sha256"])
    lit, no_accession = literature_rows(read_seeds(seeds_path), kw_by_acc)
    hashes = {r["seq_sha256"] for r in table} | {r["seq_sha256"] for r in lit if r["seq_sha256"]}
    check_features(hashes, fu)
    seqs = {
        r["seq_sha256"]: r["sequence"]
        for r in truth_table.read_tsv(work / "phaseb" / "unique_sequences.tsv.gz")
        if r["seq_sha256"] in hashes
    }
    fasta = fasta_bytes(hashes, seqs)
    log_rows = go_log + tc_log
    shared = Counter(
        (m["class"], m["label"], m["d8_class"])
        for m in members
        if m["seq_sha256"] in {r["seq_sha256"] for r in kw_rows}
    )
    log = {
        "all_sources": True,
        "truth_set_sha256": chain["truth_set_sha256"],
        "input_sha256": {
            "truth_set_triaged.tsv.gz": chain["truth_set_triaged_sha256"],
            "features.tsv.gz": manifest.sha256_file(work / "phaseb" / "features.tsv.gz"),
            "features_unique.tsv.gz": manifest.sha256_file(
                work / "phaseb" / "features_unique.tsv.gz"
            ),
            "unique_sequences.tsv.gz": chain["unique_sequences_sha256"],
            "keyword_tier.tsv.gz": manifest.sha256_file(work / "keyword_tier.tsv.gz"),
            "species.tsv": manifest.sha256_file(species_path),
            "eurotiomycetes_seeds.tsv": manifest.sha256_file(seeds_path),
            "tc_taxon_clades.tsv": manifest.sha256_file(clades_path),
        },
        "go_members_by_source_class": {
            s: dict(sorted(Counter(m["class"] for m in members if m["source_id"] == s).items()))
            for s in sorted({m["source_id"] for m in members})
        },
        "labelled_genes_without_sequence": dict(sorted(unmatched.items())),
        "table_rows_by_origin_class": {
            f"{o}:{c}": n
            for (o, c), n in sorted(Counter((r["origin"], r["class"]) for r in table).items())
        },
        "log_rows_by_reason": dict(sorted(Counter(r["reason"] for r in log_rows).items())),
        "tc_rows": len(kw_rows),
        "tc_rows_dropped_by_go_class": dict(sorted(Counter(r["detail"] for r in tc_log).items())),
        "go_members_sharing_a_tc_hash": {
            f"{c}:{lab}:{d8}": n for (c, lab, d8), n in sorted(shared.items())
        },
        "literature_rows": len(lit),
        "literature_positives": sum(r["literature_positive"] == "yes" for r in lit),
        "literature_no_accession": sorted(no_accession),
        "sequences": len(hashes),
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    out = evalio.out_dir(work)
    evalio.write_outputs(
        out,
        {
            "eval_table.tsv.gz": lambda p: truth_table.write_tsv(p, dedupe.TABLE_COLUMNS, table),
            "eval_literature.tsv": lambda p: truth_table.write_tsv(p, LITERATURE_COLUMNS, lit),
            "eval_dedupe_log.tsv": lambda p: truth_table.write_tsv(p, dedupe.LOG_COLUMNS, log_rows),
            "eval_sequences.fasta.gz": lambda p: Path(p).write_bytes(fasta),
        },
        "build_run.json",
        log,
    )
    return table, lit, log_rows, log


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--seeds", default=None, help=f"default: $PROJ_ROOT/{DEFAULT_SEEDS}")
    parser.add_argument("--tc-clades", default=str(evalio.PHASEC_DIR / "tc_taxon_clades.tsv"))
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    seeds = Path(args.seeds) if args.seeds else paths.repo_root() / DEFAULT_SEEDS
    try:
        table, lit, log_rows, log = run(
            work,
            Path(args.species),
            seeds,
            Path(args.tc_clades),
            list(argv) if argv is not None else sys.argv[1:],
        )
    except (evalio.StopError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(f"table_rows={len(table)} sequences={log['sequences']} literature={len(lit)}")
    for key, n in log["table_rows_by_origin_class"].items():
        print(f"  {key}={n}")
    for key, n in log["tc_rows_dropped_by_go_class"].items():
        print(f"  tc_dropped_go_class_{key}={n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_build.py -q`
Expected: `17 passed`.

- [ ] **Step 8: Run 08 on a scratch link of the real inputs (read-only check)**

The real work directory must not receive Phase C files before review. Link the inputs into a scratch work directory:

```bash
RW=${SCRATCH:?run on an interactive node, where SCRATCH is set}/phasec_realwork; W=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
mkdir -p $RW/phaseb/emb
for f in truth_set_triaged.tsv.gz keyword_tier.tsv.gz d8_run.json keyword_tier_run.json; do ln -sf $W/$f $RW/$f; done
for f in features.tsv.gz features_unique.tsv.gz unique_sequences.tsv.gz features_run.json sequence_members.tsv.gz; do ln -sf $W/phaseb/$f $RW/phaseb/$f; done
ln -sf $W/phaseb/emb/*.npy $W/phaseb/emb/embedding_run.json $RW/phaseb/emb/
PYTHONPATH=src:analysis/step1_compare:analysis/step1_compare/phasec $ENV_PY analysis/step1_compare/phasec/08_build_eval_tables.py --work-dir $RW
```

Expected (prototype, 17.5 s): `table_rows=23338 sequences=23351 literature=20`, `go:excluded=266`, `go:neg=19776`, `go:pos=752`, `tc:pos=2544`, `tc_dropped_go_class_excluded=65`, `tc_dropped_go_class_neg=224`, `tc_dropped_go_class_pos=208`.

- [ ] **Step 9: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/dedupe.py analysis/step1_compare/phasec/tc_taxon_clades.tsv \
  analysis/step1_compare/phasec/08_build_eval_tables.py tests/step1_compare/phasec_fixture.py \
  tests/step1_compare/test_phasec_build.py
git commit -m "phasec: evaluation table, dedupe, GO-over-T-c precedence (C-6), input checks (08)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 3: Weighted metrics

**Files:**
- Create: `analysis/step1_compare/phasec/metrics.py`
- Test: `tests/step1_compare/test_phasec_metrics.py`

**Interfaces:**
- Consumes: numpy.
- Produces: every function takes W of shape (B, n) or (n,) and returns an array of length B: `recall(W, y, call)`, `precision(W, y, call)`, `fpr(W, y, call)`, `roc_auc(W, y, score)`, `pr_auc(W, y, score)`, `precision_at_recall(W, y, score, level)`, `recall_at_fpr(W, y, score, level)` (level: scalar or one value per row of W), `fpr_at_recall(W, y, score, level, neg_mask)`, `brier(W, y, prob)`, `scored_summary(W, y, score, recall_levels=(0.8, 0.9), fpr_level=0.01) -> dict` with keys `roc_auc`, `pr_auc`, `precision_at_recall_0.8`, `precision_at_recall_0.9`, `recall_at_fpr_0.01`; `curve(W, y, score, extra=None)`; `precision_at_prevalence(rec, rate, prevalence)`; `as_weights(W, n)`. Undefined values are NaN.

The 10-row case of `test_metrics_hand_computed` (y, score): (1, 0.9) (1, 0.8) (0, 0.8) (1, 0.7) (0, 0.6) (0, 0.5) (1, 0.5) (0, 0.4) (0, 0.3) (0, 0.1); call = score >= 0.6. By hand: TP 3, FP 2, P 4, N 6, so recall 3/4, precision 3/5, FPR 2/6. ROC-AUC by pairs: the positive at 0.9 beats 6 negatives; at 0.8, 5 and a tie (0.5); at 0.7, 5; at 0.5, 3 and a tie (0.5): 20 of 24 pairs. Thresholds (TP, FP): 0.9 (1, 0), 0.8 (2, 1), 0.7 (3, 1), 0.6 (3, 2), 0.5 (4, 3), so AP = 0.25 x (1 + 2/3 + 3/4 + 4/7) = 251/336. Precision at recall 0.5: recall first reaches 0.5 at 0.8 (precision 2/3), but 0.7 gives 3/4; the step-function rule takes the highest, 3/4.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_metrics.py`:

```python
"""metrics.py: hand-computed values, ties, weights, undefined cases (spec 4, section 6)."""

import math

import pytest

np = pytest.importorskip("numpy")
sk = pytest.importorskip("sklearn.metrics")

import metrics as m  # noqa: E402

# 10 proteins; scores with two tie groups (0.8: pos + neg; 0.5: neg + pos).
Y = np.array([1, 1, 0, 1, 0, 0, 1, 0, 0, 0], dtype=bool)
S = np.array([0.9, 0.8, 0.8, 0.7, 0.6, 0.5, 0.5, 0.4, 0.3, 0.1])
CALL = S >= 0.6  # rows 0-4
ONE = np.ones(10)


def test_metrics_hand_computed():
    # P = 4, N = 6. Called: rows 0-4 -> TP 3 (rows 0, 1, 3), FP 2 (rows 2, 4).
    assert m.recall(ONE, Y, CALL)[0] == 0.75
    assert m.precision(ONE, Y, CALL)[0] == 0.6
    assert m.fpr(ONE, Y, CALL)[0] == pytest.approx(2 / 6, abs=1e-15)
    # ROC-AUC by pairs: 6 + 5.5 + 5 + 3.5 = 20 of 24 pairs (ties count 0.5)
    assert m.roc_auc(ONE, Y, S)[0] == pytest.approx(20 / 24, abs=1e-15)
    # AP: 0.25 * (1 + 2/3 + 3/4 + 4/7) = 251/336
    assert m.pr_auc(ONE, Y, S)[0] == pytest.approx(251 / 336, abs=1e-15)


def test_precision_at_recall_reads_the_step_curve_with_ties():
    # thresholds (TP, FP): 0.9 (1,0) 0.8 (2,1) 0.7 (3,1) 0.6 (3,2) 0.5 (4,3) ...
    # recall >= 0.5 first at 0.8 (precision 2/3), but 0.7 gives 3/4: the highest wins
    assert m.precision_at_recall(ONE, Y, S, 0.5)[0] == pytest.approx(0.75, abs=1e-15)
    assert m.precision_at_recall(ONE, Y, S, 0.8)[0] == pytest.approx(4 / 7, abs=1e-15)
    assert m.precision_at_recall(ONE, Y, S, 0.9)[0] == pytest.approx(4 / 7, abs=1e-15)


def test_recall_at_fpr():
    assert m.recall_at_fpr(ONE, Y, S, 0.01)[0] == 0.25  # only the 0.9 threshold has FPR 0
    assert m.recall_at_fpr(ONE, Y, S, 0.17)[0] == 0.75  # FPR 1/6 at 0.8 and at 0.7
    assert m.recall_at_fpr(ONE, Y, S, 1.0)[0] == 1.0
    # review I-3: the top-scored protein is a negative, so every non-empty call set has FPR > 0;
    # the empty call set (recall 0, FPR 0) is the only one within FPR 0.01: the result is 0.0
    y = np.array([0, 1, 1, 0, 0], dtype=bool)
    s = np.array([0.9, 0.8, 0.7, 0.2, 0.1])
    got = m.recall_at_fpr(np.ones(5), y, s, 0.01)[0]
    assert got == 0.0 and not math.isnan(got)
    assert m.scored_summary(np.ones(5), y, s)["recall_at_fpr_0.01"][0] == 0.0


def test_fpr_at_recall_inside_a_stratum():
    nsec = np.zeros(10, dtype=bool)
    nsec[[4, 5]] = True  # two negatives with scores 0.6 and 0.5
    all_neg = ~Y
    assert m.fpr_at_recall(ONE, Y, S, 0.75, all_neg)[0] == pytest.approx(1 / 6, abs=1e-15)
    assert m.fpr_at_recall(ONE, Y, S, 0.75, nsec)[0] == 0.0  # threshold 0.7 calls neither
    assert m.fpr_at_recall(ONE, Y, S, 1.0, nsec)[0] == 1.0  # threshold 0.5 calls both


def test_brier_and_prevalence():
    assert m.brier(np.ones(2), [1, 0], [0.9, 0.1])[0] == pytest.approx(0.01, abs=1e-15)
    # recall 0.8, FPR 0.1 at prevalence 0.05: 0.04 / (0.04 + 0.095)
    got = m.precision_at_prevalence(0.8, 0.1, 0.05)
    assert got == pytest.approx(0.04 / 0.135, abs=1e-15)


def test_weights_equal_duplicated_rows():
    w = np.array([2, 1, 1, 3, 1, 0, 1, 2, 1, 1], dtype=float)
    idx = np.repeat(np.arange(10), w.astype(int))
    for f in (m.roc_auc, m.pr_auc):
        assert f(w, Y, S)[0] == pytest.approx(f(np.ones(len(idx)), Y[idx], S[idx])[0], abs=1e-12)
    assert m.recall(w, Y, CALL)[0] == pytest.approx(
        m.recall(np.ones(len(idx)), Y[idx], CALL[idx])[0]
    )


def test_weighted_values_match_sklearn():
    rng = np.random.default_rng(3)
    for _ in range(20):
        n = 60
        y = rng.random(n) < 0.3
        s = np.round(rng.random(n), 1)  # many ties
        w = rng.integers(0, 4, size=n).astype(float)
        if w[y].sum() == 0 or w[~y].sum() == 0:
            continue
        assert m.roc_auc(w, y, s)[0] == pytest.approx(sk.roc_auc_score(y, s, sample_weight=w))
        assert m.pr_auc(w, y, s)[0] == pytest.approx(
            sk.average_precision_score(y, s, sample_weight=w)
        )


def test_rows_of_w_are_independent():
    W = np.stack([ONE, np.r_[np.ones(5), np.zeros(5)]])
    got = m.recall(W, Y, CALL)
    assert got[0] == 0.75 and got[1] == 1.0  # second row: positives 0, 1, 3, all called


def test_undefined_metrics_are_nan():
    # Review Focus 1: a stratum with no positives or no negatives gives NaN, not an error
    neg_only = np.zeros(4, dtype=bool)
    s = np.array([0.1, 0.2, 0.3, 0.4])
    assert math.isnan(m.recall(np.ones(4), neg_only, s > 0.2)[0])
    assert m.fpr(np.ones(4), neg_only, s > 0.2)[0] == 0.5
    for f in (m.roc_auc, m.pr_auc):
        assert math.isnan(f(np.ones(4), neg_only, s)[0])
    assert math.isnan(m.precision_at_recall(np.ones(4), neg_only, s, 0.8)[0])
    pos_only = np.ones(4, dtype=bool)
    assert math.isnan(m.fpr(np.ones(4), pos_only, s > 0.2)[0])
    assert math.isnan(m.roc_auc(np.ones(4), pos_only, s)[0])
    assert m.pr_auc(np.ones(4), pos_only, s)[0] == 1.0
    assert math.isnan(m.precision(np.ones(4), pos_only, np.zeros(4, dtype=bool))[0])
    assert math.isnan(m.recall(np.zeros((1, 0)), np.zeros(0, bool), np.zeros(0, bool))[0])


def test_non_finite_scores_are_refused():
    # Review Focus 3: a NaN score must not be sorted silently
    with pytest.raises(ValueError, match="finite"):
        m.roc_auc(ONE, Y, np.r_[S[:9], np.nan])


def test_level_can_differ_per_resample():
    W = np.stack([ONE, ONE])
    got = m.precision_at_recall(W, Y, S, np.array([0.5, 0.8]))
    assert got == pytest.approx([0.75, 4 / 7])


def test_scored_summary_equals_the_single_functions():
    W = np.stack([ONE, np.r_[np.ones(5), 2 * np.ones(5)]])
    got = m.scored_summary(W, Y, S, (0.8, 0.9), 0.01)
    assert set(got) == {"roc_auc", "pr_auc", "precision_at_recall_0.8",
                        "precision_at_recall_0.9", "recall_at_fpr_0.01"}  # fmt: skip
    assert got["roc_auc"] == pytest.approx(m.roc_auc(W, Y, S))
    assert got["pr_auc"] == pytest.approx(m.pr_auc(W, Y, S))
    assert got["precision_at_recall_0.9"] == pytest.approx(m.precision_at_recall(W, Y, S, 0.9))
    assert got["recall_at_fpr_0.01"] == pytest.approx(m.recall_at_fpr(W, Y, S, 0.01))
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_metrics.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'metrics'`.

- [ ] **Step 3: Write `analysis/step1_compare/phasec/metrics.py`**

```python
"""Weighted classification metrics for Phase C (Phase C spec 4). Needs numpy.

Every function takes a weight matrix W of shape (B, n) (or a vector of length n, read as one
row) and returns one value per row of W. A row of W holds the number of times each protein is
in one bootstrap resample; a row of ones gives the point estimate. A metric that is not defined
for a row (no positives, no negatives, no positive calls) is NaN for that row.

Definitions (thresholds are the distinct scores; a protein is called when score >= threshold):
- recall = TP / P; precision = TP / (TP + FP); FPR = FP / N (binary calls).
- ROC-AUC: area under the ROC step curve with ties counted as half (trapezoid per tie group).
- PR-AUC: average precision, sum over thresholds of (R_k - R_k-1) * P_k (sklearn definition).
- precision at recall r: the highest precision at any threshold with recall >= r.
- recall at FPR f: the highest recall at any threshold with FPR <= f (the empty call set, with
  recall 0 and FPR 0, counts as a threshold).
- FPR of a subset at recall r: the FPR inside `neg_mask` at the first (highest) threshold whose
  recall reaches r.
"""

import numpy as np

TOL = 1e-12


def as_weights(W, n: int) -> np.ndarray:
    W = np.asarray(W, dtype=np.float64)
    if W.ndim == 1:
        W = W[None, :]
    if W.ndim != 2 or W.shape[1] != n:
        raise ValueError(f"weights must have shape (B, {n}), got {W.shape}")
    return W


def _div(a, b) -> np.ndarray:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    out = np.full(np.broadcast(a, b).shape, np.nan)
    ok = b > 0
    np.divide(a, b, out=out, where=ok)
    return out


def _bool(x) -> np.ndarray:
    return np.asarray(x, dtype=bool)


def recall(W, y, call) -> np.ndarray:
    y, call = _bool(y), _bool(call)
    W = as_weights(W, len(y))
    return _div(W @ (y & call), W @ y)


def precision(W, y, call) -> np.ndarray:
    y, call = _bool(y), _bool(call)
    W = as_weights(W, len(y))
    return _div(W @ (y & call), W @ call)


def fpr(W, y, call) -> np.ndarray:
    y, call = _bool(y), _bool(call)
    W = as_weights(W, len(y))
    return _div(W @ (~y & call), W @ ~y)


def curve(W, y, score, extra=None):
    """Cumulative weighted TP and FP at each distinct threshold, from the highest score down.

    Return (tp, fp, p, n, extra_cum): tp and fp have shape (B, G) for G distinct scores; p and n
    have shape (B,); extra_cum is the cumulative weight of the rows in the boolean mask `extra`
    (shape (B, G)) or None."""
    y = _bool(y)
    score = np.asarray(score, dtype=np.float64)
    if not np.isfinite(score).all():
        raise ValueError("scores must be finite")
    W = as_weights(W, len(y))
    order = np.argsort(-score, kind="stable")
    s = score[order]
    ends = np.r_[np.nonzero(s[1:] != s[:-1])[0], len(s) - 1] if len(s) else np.array([], int)
    Ws = W[:, order]
    tp = np.cumsum(Ws * y[order], axis=1)[:, ends]
    fp = np.cumsum(Ws * ~y[order], axis=1)[:, ends]
    ex = None
    if extra is not None:
        ex = np.cumsum(Ws * _bool(extra)[order], axis=1)[:, ends]
    return tp, fp, W @ y, W @ ~y, ex


def _roc(tp, fp, p, n) -> np.ndarray:
    tp0 = np.c_[np.zeros(len(p)), tp]
    fp0 = np.c_[np.zeros(len(p)), fp]
    area = np.sum((fp0[:, 1:] - fp0[:, :-1]) * (tp0[:, 1:] + tp0[:, :-1]) / 2.0, axis=1)
    return _div(area, p * n)


def _prec(tp, fp, empty) -> np.ndarray:
    called = tp + fp
    return np.where(called > 0, tp / np.where(called > 0, called, 1.0), empty)


def _ap(tp, fp, p) -> np.ndarray:
    dtp = np.diff(np.c_[np.zeros(len(p)), tp], axis=1)
    return _div(np.sum(dtp * _prec(tp, fp, 0.0), axis=1), p)


def _prec_at_recall(tp, fp, p, level) -> np.ndarray:
    level = np.broadcast_to(np.asarray(level, dtype=np.float64), p.shape)[:, None]
    prec = _prec(tp, fp, np.nan)
    ok = (_div(tp, p[:, None]) >= level - TOL) & ~np.isnan(prec)
    best = np.where(ok, prec, -np.inf).max(axis=1, initial=-np.inf)
    best[~np.isfinite(best) | ~(p > 0)] = np.nan
    return best


def _recall_at_fpr(tp, fp, p, n, level) -> np.ndarray:
    level = np.broadcast_to(np.asarray(level, dtype=np.float64), p.shape)[:, None]
    rec = np.c_[np.zeros(len(p)), _div(tp, p[:, None])]
    rate = np.c_[np.zeros(len(p)), _div(fp, n[:, None])]
    best = np.where(rate <= level + TOL, rec, -np.inf).max(axis=1)
    best[~((p > 0) & (n > 0))] = np.nan
    return best


def roc_auc(W, y, score) -> np.ndarray:
    tp, fp, p, n, _ = curve(W, y, score)
    return _roc(tp, fp, p, n)


def pr_auc(W, y, score) -> np.ndarray:
    tp, fp, p, _, _ = curve(W, y, score)
    return _ap(tp, fp, p)


def precision_at_recall(W, y, score, level) -> np.ndarray:
    tp, fp, p, _, _ = curve(W, y, score)
    return _prec_at_recall(tp, fp, p, level)


def recall_at_fpr(W, y, score, level) -> np.ndarray:
    tp, fp, p, n, _ = curve(W, y, score)
    return _recall_at_fpr(tp, fp, p, n, level)


def scored_summary(W, y, score, recall_levels=(0.8, 0.9), fpr_level=0.01) -> dict:
    """ROC-AUC, PR-AUC, precision at each recall level, recall at fpr_level; one curve."""
    tp, fp, p, n, _ = curve(W, y, score)
    out = {"roc_auc": _roc(tp, fp, p, n), "pr_auc": _ap(tp, fp, p)}
    for r in recall_levels:
        out[f"precision_at_recall_{r}"] = _prec_at_recall(tp, fp, p, r)
    out[f"recall_at_fpr_{fpr_level}"] = _recall_at_fpr(tp, fp, p, n, fpr_level)
    return out


def fpr_at_recall(W, y, score, level, neg_mask) -> np.ndarray:
    """FPR inside the rows of `neg_mask` at the first threshold where recall >= level."""
    neg_mask = _bool(neg_mask) & ~_bool(y)
    tp, _, p, _, sub = curve(W, y, score, extra=neg_mask)
    W = as_weights(W, len(neg_mask))
    total = W @ neg_mask
    level = np.broadcast_to(np.asarray(level, dtype=np.float64), p.shape)[:, None]
    reached = _div(tp, p[:, None]) >= level - TOL
    out = np.full(len(p), np.nan)
    has = reached.any(axis=1) & (p > 0) & (total > 0)
    first = reached.argmax(axis=1)
    rows = np.nonzero(has)[0]
    out[rows] = sub[rows, first[rows]] / total[rows]
    return out


def brier(W, y, prob) -> np.ndarray:
    y = _bool(y).astype(np.float64)
    prob = np.asarray(prob, dtype=np.float64)
    W = as_weights(W, len(y))
    return _div(W @ (prob - y) ** 2, W.sum(axis=1))


def precision_at_prevalence(rec, rate, prevalence: float) -> np.ndarray:
    """Precision at an assumed prevalence (spec 4, ruling C-11):
    recall * pi / (recall * pi + FPR * (1 - pi))."""
    rec = np.asarray(rec, dtype=np.float64)
    rate = np.asarray(rate, dtype=np.float64)
    return _div(rec * prevalence, rec * prevalence + rate * (1.0 - prevalence))
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_metrics.py -q`
Expected: `12 passed`.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/metrics.py tests/step1_compare/test_phasec_metrics.py
git commit -m "phasec: weighted metrics for the cluster bootstrap (spec 4)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 4: Cluster bootstrap

**Files:**
- Create: `analysis/step1_compare/phasec/bootstrap.py`
- Test: `tests/step1_compare/test_phasec_bootstrap.py`

**Interfaces:**
- Consumes: `metrics.recall`, `metrics.roc_auc` (Task 3, tests only).
- Produces: `bootstrap.N_RESAMPLES = 2000`, `CI = (2.5, 97.5)`, `seed_for(base_seed, name) -> int`, `cluster_weights(cluster_ids, n_resamples, seed) -> np.ndarray` (int32, shape (n_resamples, n)), `interval(values) -> {"lo", "hi", "n_defined"}` (None when no resample is defined), `half_width(ci) -> float | None`.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_bootstrap.py`:

```python
"""bootstrap.py: whole clusters move together, seeds, pairing, empty resamples (spec 4, 6)."""

import pytest

np = pytest.importorskip("numpy")

import bootstrap as bs  # noqa: E402
import metrics as m  # noqa: E402

CLUSTERS = ["c1", "c1", "c1", "c2", "c3", "c3", "c4", "c5", "c5", "c5"]


def test_bootstrap_resamples_whole_clusters():
    W = bs.cluster_weights(CLUSTERS, 500, seed=1)
    assert W.shape == (500, 10) and W.dtype == np.int32
    for members in ([0, 1, 2], [4, 5], [7, 8, 9]):
        assert (W[:, members] == W[:, members[:1]]).all()
    # each resample draws as many clusters as there are (5)
    first = [0, 3, 4, 6, 7]  # one member per cluster
    assert (W[:, first].sum(axis=1) == 5).all()


def test_bootstrap_seeded():
    a = bs.cluster_weights(CLUSTERS, 50, seed=11)
    b = bs.cluster_weights(CLUSTERS, 50, seed=11)
    c = bs.cluster_weights(CLUSTERS, 50, seed=12)
    assert (a == b).all() and not (a == c).all()
    assert bs.seed_for(20261001, "S1:all") == bs.seed_for(20261001, "S1:all")
    assert bs.seed_for(20261001, "S1:all") != bs.seed_for(20261001, "S1:Scer_SGD")


def test_row_order_does_not_change_the_cluster_draws():
    perm = [9, 0, 5, 3, 1, 8, 2, 7, 4, 6]
    a = bs.cluster_weights(CLUSTERS, 30, seed=5)
    b = bs.cluster_weights([CLUSTERS[i] for i in perm], 30, seed=5)
    assert (a[:, perm] == b).all()


def test_paired_bootstrap_uses_same_indices():
    y = np.array([1, 1, 0, 1, 0, 0, 1, 0, 0, 0], dtype=bool)
    s = np.array([0.9, 0.8, 0.8, 0.7, 0.6, 0.5, 0.5, 0.4, 0.3, 0.1])
    W = bs.cluster_weights(CLUSTERS, 400, seed=2)
    a = m.roc_auc(W, y, s)
    b = m.roc_auc(W, y, s)  # a second candidate with the same scores
    diff = bs.interval(a - b)
    assert diff["lo"] == 0.0 and diff["hi"] == 0.0  # paired: identical candidates differ by 0
    unpaired = m.roc_auc(bs.cluster_weights(CLUSTERS, 400, seed=3), y, s)
    assert bs.interval(a - unpaired)["hi"] > 0  # different indices do not cancel


def test_interval_skips_undefined_resamples():
    # Review Focus 2: a resample without positives has NaN recall; the interval uses the rest
    vals = np.array([np.nan, 0.5, 0.6, 0.7, np.nan])
    ci = bs.interval(vals)
    assert (
        ci["n_defined"] == 3
        and ci["lo"] == pytest.approx(0.505)
        and ci["hi"] == pytest.approx(0.695)
    )
    assert bs.interval([np.nan, np.nan]) == {"lo": None, "hi": None, "n_defined": 0}
    assert bs.half_width({"lo": None, "hi": None, "n_defined": 0}) is None
    assert bs.half_width({"lo": 0.6, "hi": 0.8, "n_defined": 9}) == pytest.approx(0.1)


def test_all_positives_in_one_cluster_gives_some_empty_resamples():
    y = np.array([1, 1, 0, 0, 0, 0], dtype=bool)
    clusters = ["p", "p", "a", "b", "c", "d"]
    W = bs.cluster_weights(clusters, 300, seed=4)
    rec = m.recall(W, y, np.array([1, 0, 0, 0, 0, 0], dtype=bool))
    assert np.isnan(rec).any() and not np.isnan(rec).all()
    ci = bs.interval(rec)
    assert 0 < ci["n_defined"] < 300 and ci["lo"] == ci["hi"] == 0.5


def test_empty_test_set_gives_an_empty_matrix():
    assert bs.cluster_weights([], 7, seed=1).shape == (7, 0)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_bootstrap.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'bootstrap'`.

- [ ] **Step 3: Write `analysis/step1_compare/phasec/bootstrap.py`**

```python
"""Cluster bootstrap for Phase C (Phase C spec 4). Needs numpy.

The unit is the MMseqs2 cluster: a resample draws clusters with replacement, and every member
of a drawn cluster enters the resample as often as its cluster was drawn. The result is a
weight matrix W (B x n): W[b, i] = number of times protein i is in resample b. Metrics read W
(metrics.py), so one W serves every candidate, variant and stratum of a test set: the
differences are paired. The models are not refitted inside the bootstrap; the intervals
describe sampling of the test set only.
"""

import hashlib

import numpy as np

N_RESAMPLES = 2000
CI = (2.5, 97.5)


def seed_for(base_seed: int, name: str) -> int:
    """A fixed seed per test set: the same name and base seed give the same resamples."""
    digest = hashlib.sha256(f"{base_seed}:{name}".encode()).hexdigest()
    return int(digest[:16], 16)


def cluster_weights(cluster_ids, n_resamples: int, seed: int) -> np.ndarray:
    """Weight matrix (n_resamples x len(cluster_ids)), int32.

    Clusters are sorted by id before drawing, so the row order of the proteins does not change
    which clusters a resample holds."""
    ids = np.asarray(cluster_ids, dtype=object)
    if len(ids) == 0:
        return np.zeros((n_resamples, 0), dtype=np.int32)
    uniq, inverse = np.unique(ids.astype(str), return_inverse=True)
    k = len(uniq)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, k, size=(n_resamples, k))
    offset = (np.arange(n_resamples) * k)[:, None]
    counts = np.bincount((draws + offset).ravel(), minlength=n_resamples * k)
    counts = counts.reshape(n_resamples, k).astype(np.int32)
    return counts[:, inverse]


def interval(values) -> dict:
    """Percentile interval over the resamples where the metric is defined (not NaN)."""
    v = np.asarray(values, dtype=np.float64)
    ok = v[~np.isnan(v)]
    if len(ok) == 0:
        return {"lo": None, "hi": None, "n_defined": 0}
    lo, hi = np.percentile(ok, CI)
    return {"lo": float(lo), "hi": float(hi), "n_defined": int(len(ok))}


def half_width(ci: dict):
    """(hi - lo) / 2, or None when the interval is not defined."""
    if ci["lo"] is None or ci["hi"] is None:
        return None
    return (ci["hi"] - ci["lo"]) / 2.0
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_bootstrap.py -q`
Expected: `7 passed`.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/bootstrap.py tests/step1_compare/test_phasec_bootstrap.py
git commit -m "phasec: cluster bootstrap weights, seeds, intervals (spec 4)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 5: Rules R0, R1, R2

**Files:**
- Create: `analysis/step1_compare/phasec/rules.py`
- Test: `tests/step1_compare/test_phasec_rules.py`

**Interfaces:**
- Consumes: numpy.
- Produces: `rules.GPI_RANK`, `G_VALUES = ("highly_probable", "probable", "weakly")`, `T_VALUES = (0.20, 0.25, 0.30, 0.35, 0.40)` (owner decision 2026-10-01; review I-2: at 0.10, R2 equals R0 on the real data), `RULES = ("R0", "R1", "R2")`, `gpi_rank(calls) -> np.ndarray` (ValueError on an unknown call), `rule_call(rule, sp, rank, st, g=None, t=None) -> np.ndarray[bool]`, `youden(y, call) -> float`, `fit_rule(rule, y, sp, rank, st) -> {"g", "t", "j"}`.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_rules.py`:

```python
"""rules.py: R0-R2 truth table, cut-off values, Youden fit (spec 3.3, section 6)."""

import itertools
import math

import pytest

np = pytest.importorskip("numpy")

import rules  # noqa: E402

GPI = ("highly_probable", "probable", "weakly", "none", "too_short")


def test_rule_truth_table():
    # every SP x GPI class x Ser+Thr combination, with Ser+Thr exactly at the cut-off t = 0.25
    st_values = (0.2499, 0.25, 0.30)
    for sp, gpi, st in itertools.product((True, False), GPI, st_values):
        rank = rules.gpi_rank([gpi])
        for g in rules.G_VALUES:
            gpi_ok = rules.GPI_RANK[gpi] >= rules.GPI_RANK[g]
            assert rules.rule_call("R0", [sp], rank, [st])[0] == sp
            assert rules.rule_call("R1", [sp], rank, [st], g)[0] == (sp and gpi_ok)
            want = sp and (gpi_ok or st >= 0.25)
            assert rules.rule_call("R2", [sp], rank, [st], g, 0.25)[0] == want, (sp, gpi, st, g)


def test_cut_off_values_parse_equal_to_the_grid():
    # ser_thr_frac is stored with 6 decimals; "0.300000" must reach t = 0.30
    for t in rules.T_VALUES:
        st = float(f"{t:.6f}")
        assert rules.rule_call("R2", [True], rules.gpi_rank(["none"]), [st], "weakly", t)[0]
        below = float(f"{t - 0.000001:.6f}")
        assert not rules.rule_call("R2", [True], rules.gpi_rank(["none"]), [below], "weakly", t)[0]
    assert rules.T_VALUES[2] == 0.30 and 0.2 + 0.05 * 2 != 0.30  # why the grid is literal


def test_t_grid_is_0_20_to_0_40():
    # owner decision 2026-10-01: the grid starts at 0.20 (at 0.10, R2 equals R0 on the real data)
    assert rules.T_VALUES == (0.20, 0.25, 0.30, 0.35, 0.40)
    # positives at Ser+Thr 0.15, negatives at 0.05: a grid with 0.10 or 0.15 would give J = 1
    y = np.array([1] * 4 + [0] * 4, dtype=bool)
    rank = rules.gpi_rank(["none"] * 8)
    st = np.array([0.15] * 4 + [0.05] * 4)
    fit = rules.fit_rule("R2", y, np.ones(8, dtype=bool), rank, st)
    assert fit == {"g": "highly_probable", "t": 0.40, "j": 0.0}  # no t of the grid calls them


def test_too_short_is_never_gpi():
    rank = rules.gpi_rank(["too_short"])
    assert not rules.rule_call("R1", [True], rank, [0.0], "weakly")[0]


def test_unknown_gpi_call_is_refused():
    with pytest.raises(ValueError, match="unknown gpi_call"):
        rules.gpi_rank(["maybe"])


def test_youden():
    y = [1, 1, 0, 0]
    assert rules.youden(y, [1, 0, 1, 0]) == 0.0
    assert rules.youden(y, [1, 1, 0, 0]) == 1.0
    assert math.isnan(rules.youden([0, 0], [1, 0]))


def test_fit_rule_finds_the_planted_cut_and_breaks_ties_to_the_stricter_rule():
    # positives: SP, no GPI, Ser+Thr 0.32; negatives: SP, no GPI, Ser+Thr 0.20
    y = np.array([1] * 5 + [0] * 5, dtype=bool)
    sp = np.ones(10, dtype=bool)
    rank = rules.gpi_rank(["none"] * 10)
    st = np.array([0.32] * 5 + [0.20] * 5)
    fit = rules.fit_rule("R2", y, sp, rank, st)
    # t = 0.25 and t = 0.30 both give J = 1; the grid runs from 0.40 down, so 0.30 wins
    assert fit == {"g": "highly_probable", "t": 0.30, "j": 1.0}
    assert rules.fit_rule("R0", y, sp, rank, st)["j"] == 0.0


def test_fit_rule_g_for_r1():
    y = np.array([1, 1, 0, 0], dtype=bool)
    rank = rules.gpi_rank(["probable", "weakly", "weakly", "none"])
    fit = rules.fit_rule("R1", y, np.ones(4, dtype=bool), rank, np.zeros(4))
    # g = probable: J = 0.5; g = weakly: J = 1 - 0.5 = 0.5; tie -> probable (stricter)
    assert fit["g"] == "probable" and fit["j"] == 0.5


def test_fit_rule_refuses_rows_without_a_class():
    with pytest.raises(ValueError, match="lack a class"):
        rules.fit_rule(
            "R2", np.ones(3, bool), np.ones(3, bool), rules.gpi_rank(["none"] * 3), np.zeros(3)
        )
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_rules.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'rules'`.

- [ ] **Step 3: Write `analysis/step1_compare/phasec/rules.py`**

```python
"""The rule candidates R0, R1, R2 (Phase C spec 3.3, rulings C-2, C-3). Needs numpy.

R0: SignalP call is SP.
R1: SP and PredGPI class at or above g.
R2: SP and (PredGPI class at or above g, or ser_thr_frac >= t). R2 is "the rule".

g takes highly_probable, probable, weakly (PredGPI scores 1.0, 0.70, 0.55); t takes 0.20 to 0.40
in steps of 0.05 (owner decision 2026-10-01: at t = 0.10 R2 equals R0 on the real data, because
96.7% of the SP proteins have ser_thr_frac >= 0.10). g and t are fitted by Youden's J (recall minus FPR) on the training rows of
one outer training set. A rule has no fitted model, so J on the training rows equals J on the
pooled inner out-of-fold rows. Ties in J: the first setting in grid order wins (g from
highly_probable down, then t from 0.40 down), that is the stricter rule.
"""

import numpy as np

GPI_RANK = {"highly_probable": 3, "probable": 2, "weakly": 1, "none": 0, "too_short": 0}
G_VALUES = ("highly_probable", "probable", "weakly")
T_VALUES = (0.20, 0.25, 0.30, 0.35, 0.40)  # literals: 0.2 + 0.05 * 2 != 0.30
RULES = ("R0", "R1", "R2")


def gpi_rank(calls) -> np.ndarray:
    try:
        return np.array([GPI_RANK[c] for c in calls], dtype=np.int8)
    except KeyError as exc:
        raise ValueError(f"unknown gpi_call {exc.args[0]!r}") from exc


def rule_call(rule: str, sp, rank, st, g: str | None = None, t: float | None = None):
    """Boolean call of `rule` for arrays sp (bool), rank (gpi_rank), st (ser_thr_frac)."""
    sp = np.asarray(sp, dtype=bool)
    if rule == "R0":
        return sp.copy()
    gpi = np.asarray(rank) >= GPI_RANK[g]
    if rule == "R1":
        return sp & gpi
    if rule == "R2":
        return sp & (gpi | (np.asarray(st, dtype=np.float64) >= t))
    raise ValueError(f"unknown rule {rule!r}")


def youden(y, call) -> float:
    """Recall minus FPR; NaN when y has no positives or no negatives."""
    y = np.asarray(y, dtype=bool)
    call = np.asarray(call, dtype=bool)
    p, n = y.sum(), (~y).sum()
    if p == 0 or n == 0:
        return float("nan")
    return float((y & call).sum() / p - (~y & call).sum() / n)


def fit_rule(rule: str, y, sp, rank, st) -> dict:
    """g and t with the highest Youden's J on these rows. R0 has no parameter."""
    if rule == "R0":
        return {"g": None, "t": None, "j": youden(y, rule_call("R0", sp, rank, st))}
    best = None
    t_grid = T_VALUES[::-1] if rule == "R2" else (None,)
    for g in G_VALUES:
        for t in t_grid:
            j = youden(y, rule_call(rule, sp, rank, st, g, t))
            if np.isnan(j):
                raise ValueError("Youden's J is not defined: the rows lack a class")
            if best is None or j > best["j"]:
                best = {"g": g, "t": t, "j": j}
    return best
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_rules.py -q`
Expected: `9 passed`.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/rules.py tests/step1_compare/test_phasec_rules.py
git commit -m "phasec: rules R0-R2 and the Youden fit of g and t (spec 3.3, C-2, C-3)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 6: Universe and models (nested protocol)

**Files:**
- Create: `analysis/step1_compare/phasec/universe.py`, `analysis/step1_compare/phasec/models.py`
- Test: `tests/step1_compare/test_phasec_models.py`

**Interfaces:**
- Consumes: `evalio.StopError`, `evalio.SEED`, `evalio.float_or_stop`, `evalio.read_json` (Task 1); `metrics.pr_auc`, `metrics.roc_auc` (Task 3); `rules.*` (Task 5); `truth_table.read_tsv`.
- Produces: `universe.Universe` (fields `hashes`, `length`, `sp`, `sp_prob`, `rank`, `gpi_prob`, `st`, `comp`, `cterm_row`, `nterm`, `cterm`, `index`; method `rows(hashes) -> np.ndarray`), `universe.EMBEDDINGS`, `MODELS`, `MAX_RESIDUES = 1022`, `check_cterm_rows(length, cterm_row)`, `composition(sequence)`, `embedding(u, name, idx)`, `features(u, candidate, idx, h_variant=None)`, `load(phaseb, verify=True) -> Universe`. `models.CANDIDATES = ("B0", "B1", "R0", "R1", "R2", "M8", "M35", "M8-C", "M35-C", "H")`, `LR_CANDIDATES`, `ML_CANDIDATES`, `H_VARIANTS`, `C_GRID = (0.001, 0.003, 0.01, 0.1, 1.0, 10.0)` (ruling C-12), `INNER_FOLDS = 3`, `ModelError(ValueError)`, `inner_folds(y, groups, seed)`, `make_lr(C)`, `oof_scores(X, y, folds, C, notes)`, `youden_threshold(y, score) -> float`, `fit_platt(y, score) -> (a, b)`, `platt(score, a, b)`, `fit_unit(u, train_idx, y, groups, score_idx, seed, candidates=CANDIDATES) -> (params, scores)`; work-unit functions for Task 8: `init_worker(phaseb: str)`, `run_task(task: dict, u=None) -> (key, params, scores)` (task keys `key`, `train`, `y`, `groups`, `score`, `seed`, `candidates`).

Measured on a real unit (S1 fold 0, V-go, 5,720 rows, 1 BLAS thread, c01): all candidates together about 45 s with the old grid of four C values (see Facts). The grid of ruling C-12 has six values, so each LR candidate makes 18 instead of 12 inner fits and H 72 instead of 48; the added values are small C, which fit fastest (1.4 s at C = 0.01 against 7.3 s at C = 10 on 7,152 x 480). The new grid was not timed on the real data (not run). The chance-level test is a real contamination check: in the prototype, a mutation that trained every outer fold on all rows gave an out-of-fold AUC of 0.970 on the same shuffled labels (the test requires 0.45 to 0.55; five seeds gave 0.493 to 0.532 without the mutation).

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_models.py`:

```python
"""models.py and universe.py: rows, nested protocol, learnability, chance level (spec 3.3, 6)."""

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import evalio  # noqa: E402
import models  # noqa: E402
import universe  # noqa: E402
from metrics import roc_auc  # noqa: E402
from sklearn.model_selection import StratifiedGroupKFold  # noqa: E402

M8, M35 = universe.MODELS


def toy_universe(n, dim, seed, y=None, signal=0.0, lengths=None, cterm=None):
    """Universe with random features; embedding column 0 carries `signal` for y == 1."""
    rng = np.random.default_rng(seed)
    y = np.zeros(n, dtype=bool) if y is None else np.asarray(y, dtype=bool)
    lengths = np.full(n, 300) if lengths is None else np.asarray(lengths)
    crow = np.full(n, -1)
    long = np.nonzero(lengths > 1022)[0]
    crow[long] = np.arange(len(long))
    nterm = rng.normal(size=(n, dim))
    nterm[:, 0] += signal * y
    u = universe.Universe(
        hashes=[f"h{i:05d}" for i in range(n)],
        length=lengths,
        sp=rng.random(n) < 0.5,
        sp_prob=rng.random(n),
        rank=rng.integers(0, 4, size=n).astype(np.int8),
        gpi_prob=rng.random(n),
        st=rng.random(n) * 0.4,
        comp=rng.dirichlet(np.ones(20), size=n),
        cterm_row=crow,
    )
    for model in universe.MODELS:
        u.nterm[model] = nterm
        u.cterm[model] = rng.normal(size=(len(long), dim)) if cterm is None else cterm
    return u


def test_cterm_variant_row_selection():
    # proteins of 1,022 aa (no C-terminal row) and 1,023 aa (C-terminal row 0)
    cterm = np.full((1, 4), 7.0)
    u = toy_universe(2, 4, seed=1, lengths=[1022, 1023], cterm=cterm)
    assert list(u.cterm_row) == [-1, 0]
    for name in ("M8-C", "M35-C"):
        X = universe.features(u, name, [0, 1])
        assert (X[0] == u.nterm[M8][0]).all()
        assert (X[1] == 7.0).all()
    for name in ("M8", "M35"):
        assert (universe.features(u, name, [1]) == u.nterm[M8][1]).all()


def test_cterm_rows_must_match_lengths():
    with pytest.raises(evalio.StopError, match="disagree"):
        universe.check_cterm_rows(np.array([1023]), np.array([-1]))
    with pytest.raises(evalio.StopError, match="disagree"):
        universe.check_cterm_rows(np.array([1022]), np.array([0]))


def test_feature_matrices_of_the_baselines_and_h():
    u = toy_universe(3, 4, seed=2, lengths=[100, 200, 400])
    assert universe.features(u, "B0", [0, 2])[:, 0] == pytest.approx(np.log([100, 400]))
    assert universe.features(u, "B1", [1]).shape == (1, 21)
    h = universe.features(u, "H", [0], "M35")
    assert h.shape == (1, 7) and h[0, 4] == u.sp_prob[0] and h[0, 6] == u.st[0]


def test_youden_threshold_and_platt():
    y = np.array([1, 1, 0, 1, 0, 0], dtype=bool)
    s = np.array([3.0, 2.0, 2.0, 1.0, 0.0, -1.0])
    # thresholds 3: J=1/3; 2: 2/3-1/3=1/3; 1: 1-1/3=2/3; 0: 1-2/3 -> threshold 1.0
    assert models.youden_threshold(y, s) == 1.0
    a, b = models.fit_platt(np.r_[y, y], np.r_[s, s])
    p = models.platt(s, a, b)
    assert a > 0 and (np.diff(p[[5, 4, 3, 0]]) > 0).all()


def test_inner_folds_stop_without_enough_positive_clusters():
    # Review Focus 5: 2 positive clusters cannot fill 3 inner folds
    y = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0], dtype=bool)
    groups = ["a", "a", "b", "c", "d", "e", "f", "g", "h"]
    with pytest.raises(models.ModelError, match="2 positive and 6 negative clusters"):
        models.inner_folds(y, groups, seed=1)


def test_inner_folds_keep_clusters_whole():
    rng = np.random.default_rng(0)
    y = rng.random(90) < 0.3
    groups = np.array([f"c{i // 3}" for i in range(90)])
    for train, test in models.inner_folds(y, groups, seed=5):
        assert not set(groups[train]) & set(groups[test])
        assert y[test].any() and not y[test].all()


def test_planted_signal_is_recovered():
    n = 400
    rng = np.random.default_rng(11)
    y = rng.random(n) < 0.25
    u = toy_universe(n, 10, seed=12, y=y, signal=2.5)
    train, test = np.arange(0, 300), np.arange(300, 400)
    params, scores = models.fit_unit(
        u, train, y[train], [f"c{i}" for i in train], test, 3, ("M8", "H")
    )
    assert roc_auc(np.ones(100), y[test], scores["M8"]["score"])[0] > 0.9
    assert params["H"]["h_variant"] in models.H_VARIANTS
    assert params["M8"]["C"] in models.C_GRID and params["M8"]["n_train"] == 300


def _outer_oof(u, y, groups, candidates, seed):
    oof = np.empty(len(y))
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    for train, test in cv.split(np.zeros(len(y)), y, groups):
        _, s = models.fit_unit(u, train, y[train], groups[train], test, seed, candidates)
        oof[test] = s[candidates[0]]["score"]
    return oof


def test_shuffled_labels_give_chance_auc():
    # Labels drawn at random, independent of the features. 600 proteins x 300 dimensions: the
    # model can memorise its training rows, so a test fold that leaked into training would give
    # an AUC far above 0.5 (0.97 in the prototype when every fold trained on all rows).
    n = 600
    rng = np.random.default_rng(21)
    y = rng.random(n) < 0.3
    u = toy_universe(n, 300, seed=121)
    groups = np.array([f"c{i}" for i in range(n)])
    oof = _outer_oof(u, y, groups, ("M8",), seed=evalio.SEED)
    auc = roc_auc(np.ones(n), y, oof)[0]
    assert 0.45 <= auc <= 0.55, auc


def test_rules_in_fit_unit_use_training_rows_only():
    n = 60
    rng = np.random.default_rng(4)
    y = rng.random(n) < 0.4
    u = toy_universe(n, 3, seed=5, y=y)
    train = np.arange(40)
    params, scores = models.fit_unit(
        u, train, y[train], [f"c{i}" for i in train], np.arange(40, 60), 1, ("R0", "R2")
    )
    assert params["R0"]["g"] is None and params["R2"]["g"] in (
        "highly_probable",
        "probable",
        "weakly",
    )
    assert scores["R2"]["score"] is None and scores["R2"]["call"].dtype == bool
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_models.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'universe'`.

- [ ] **Step 3: Write `analysis/step1_compare/phasec/universe.py`**

```python
"""Per-sequence inputs of the Phase C models: features, composition, embeddings. Needs numpy.

One row per unique sequence of Phase B (`row` of features_unique.tsv.gz and of
emb/<model>.nterm.npy). The C-terminal candidates (M8-C, M35-C) read emb/<model>.cterm.npy row
`cterm_row` for proteins longer than 1,022 aa and the N-terminal row for shorter proteins
(Phase C spec 3.3).
"""

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import evalio
import numpy as np
import rules
import truth_table

AA20 = "ACDEFGHIKLMNPQRSTVWY"
MAX_RESIDUES = 1022
EMBEDDINGS = {
    "M8": ("esm2_t6_8M_UR50D", False),
    "M35": ("esm2_t12_35M_UR50D", False),
    "M8-C": ("esm2_t6_8M_UR50D", True),
    "M35-C": ("esm2_t12_35M_UR50D", True),
}
MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")


@dataclass
class Universe:
    hashes: list[str]
    length: np.ndarray
    sp: np.ndarray
    sp_prob: np.ndarray
    rank: np.ndarray
    gpi_prob: np.ndarray
    st: np.ndarray
    comp: np.ndarray
    cterm_row: np.ndarray
    nterm: dict = field(default_factory=dict)
    cterm: dict = field(default_factory=dict)
    index: dict = field(default_factory=dict)

    def __post_init__(self):
        self.index = {h: i for i, h in enumerate(self.hashes)}
        check_cterm_rows(self.length, self.cterm_row)

    def rows(self, hashes) -> np.ndarray:
        try:
            return np.array([self.index[h] for h in hashes], dtype=np.int64)
        except KeyError as exc:
            raise evalio.StopError(f"hash {exc.args[0]} is not a Phase B unique sequence") from exc


def check_cterm_rows(length, cterm_row) -> None:
    long = np.asarray(length) > MAX_RESIDUES
    has = np.asarray(cterm_row) >= 0
    bad = np.nonzero(long != has)[0]
    if len(bad):
        raise evalio.StopError(
            f"row {bad[0]}: length {int(length[bad[0]])} and cterm_row {int(cterm_row[bad[0]])} "
            f"disagree (a C-terminal row exists exactly for proteins over {MAX_RESIDUES} aa)"
        )


def composition(sequence: str) -> list[float]:
    n = len(sequence)
    return [sequence.count(a) / n for a in AA20]


def embedding(u: Universe, name: str, idx) -> np.ndarray:
    """Embedding rows of candidate `name` (M8, M35, M8-C, M35-C) for universe rows idx."""
    model, use_cterm = EMBEDDINGS[name]
    X = np.asarray(u.nterm[model][idx], dtype=np.float64)
    if use_cterm:
        crow = u.cterm_row[idx]
        long = crow >= 0
        if long.any():
            X[long] = np.asarray(u.cterm[model][crow[long]], dtype=np.float64)
    return X


def features(u: Universe, candidate: str, idx, h_variant: str | None = None) -> np.ndarray:
    idx = np.asarray(idx, dtype=np.int64)
    loglen = np.log(u.length[idx].astype(np.float64))[:, None]
    if candidate == "B0":
        return loglen
    if candidate == "B1":
        return np.hstack([u.comp[idx], loglen])
    if candidate in EMBEDDINGS:
        return embedding(u, candidate, idx)
    if candidate == "H":
        extra = np.c_[u.sp_prob[idx], u.gpi_prob[idx], u.st[idx]]
        return np.hstack([embedding(u, h_variant, idx), extra])
    raise ValueError(f"candidate {candidate!r} has no feature matrix")


def array_sha256(arr) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def load(phaseb: Path, verify: bool = True) -> Universe:
    """Read features_unique, unique_sequences and the four embedding matrices of Phase B."""
    phaseb = Path(phaseb)
    fu = truth_table.read_tsv(phaseb / "features_unique.tsv.gz")
    seqs = truth_table.read_tsv(phaseb / "unique_sequences.tsv.gz")
    if [r["seq_sha256"] for r in fu] != [r["seq_sha256"] for r in seqs]:
        raise evalio.StopError("features_unique and unique_sequences list other hashes")
    for i, r in enumerate(fu):
        if r["row"] != str(i):
            raise evalio.StopError(f"features_unique.tsv.gz row {i} has row value {r['row']}")
    num = {}
    for col in ("sp_prob", "gpi_prob", "ser_thr_frac"):
        num[col] = np.array([evalio.float_or_stop(r[col], f"{r['seq_sha256']} {col}") for r in fu])
    run = evalio.read_json(phaseb / "emb" / "embedding_run.json")
    u = Universe(
        hashes=[r["seq_sha256"] for r in fu],
        length=np.array([int(r["length"]) for r in fu], dtype=np.int64),
        sp=np.array([r["sp_prediction"] == "SP" for r in fu], dtype=bool),
        sp_prob=num["sp_prob"],
        rank=rules.gpi_rank([r["gpi_call"] for r in fu]),
        gpi_prob=num["gpi_prob"],
        st=num["ser_thr_frac"],
        comp=np.array([composition(r["sequence"]) for r in seqs], dtype=np.float64),
        cterm_row=np.array([int(r["cterm_row"]) if r["cterm_row"] else -1 for r in fu]),
    )
    n_long = int((u.cterm_row >= 0).sum())
    for model in MODELS:
        for window, store, rows in (("nterm", u.nterm, len(fu)), ("cterm", u.cterm, n_long)):
            arr = np.load(phaseb / "emb" / f"{model}.{window}.npy", mmap_mode="r")
            want = run.get("models", {}).get(model, {}).get(window, {})
            if list(arr.shape) != want.get("shape") or arr.shape[0] != rows:
                raise evalio.StopError(
                    f"emb/{model}.{window}.npy has shape {arr.shape}; embedding_run.json says "
                    f"{want.get('shape')} and the features need {rows} rows"
                )
            if verify and array_sha256(arr) != want.get("array_sha256"):
                raise evalio.StopError(f"emb/{model}.{window}.npy differs from embedding_run.json")
            store[model] = arr
    return u
```

- [ ] **Step 4: Write `analysis/step1_compare/phasec/models.py`**

```python
"""Fit and score the Phase C candidates on one outer training set (Phase C spec 3.3).

Nested protocol (one rule for every fitted quantity). For one outer training set:
1. Inner folds: StratifiedGroupKFold(3, shuffle, seed) over the training rows, groups =
   MMseqs2 clusters. Each inner test fold must hold both classes, else ModelError.
2. Logistic-regression candidates (B0, B1, M8, M35, M8-C, M35-C; pipeline StandardScaler +
   LogisticRegression(class_weight="balanced")): for each C in C_GRID, the inner out-of-fold
   decision values; C = the value with the highest inner out-of-fold PR-AUC (ties: the smaller
   C). H: the same over (ESM variant, C); ties: the first variant in M8, M35, M8-C, M35-C order.
3. With the chosen setting: the decision threshold maximises Youden's J on the inner
   out-of-fold values (ties: the higher threshold); Platt scaling (a, b) is fitted on the same
   values (ruling C-10). The final pipeline is refitted on all training rows.
4. Rules: g and t by Youden's J on the training rows (rules.fit_rule).
Only the training rows' labels are read. Scores are decision values; `prob` is the Platt
probability; `call` is score >= threshold.
"""

import warnings

import numpy as np
import rules
import universe
from metrics import pr_auc
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

CANDIDATES = ("B0", "B1", "R0", "R1", "R2", "M8", "M35", "M8-C", "M35-C", "H")
LR_CANDIDATES = ("B0", "B1", "M8", "M35", "M8-C", "M35-C")
ML_CANDIDATES = ("M8", "M35", "M8-C", "M35-C", "H")
H_VARIANTS = ("M8", "M35", "M8-C", "M35-C")
C_GRID = (0.001, 0.003, 0.01, 0.1, 1.0, 10.0)  # ruling C-12: S1 fold 0 chose 0.01, the old edge
INNER_FOLDS = 3
MAX_ITER = 5000
SCORE_CHUNK = 20000


class ModelError(ValueError):
    """An outer training set cannot be fitted (for example an inner fold without positives)."""


def inner_folds(y, groups, seed: int, n_folds: int = INNER_FOLDS):
    y = np.asarray(y, dtype=bool)
    groups = np.asarray(groups).astype(str)
    n_pos, n_neg = len(set(groups[y])), len(set(groups[~y]))
    if n_pos < n_folds or n_neg < n_folds:
        raise ModelError(
            f"{n_pos} positive and {n_neg} negative clusters; the inner {n_folds}-fold split "
            "needs at least one cluster of each class per fold"
        )
    cv = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    folds = list(cv.split(np.zeros(len(y)), y, groups))
    for k, (_, test) in enumerate(folds):
        if y[test].all() or not y[test].any():
            raise ModelError(
                f"inner fold {k} has no {'negatives' if y[test].all() else 'positives'}"
            )
    return folds


def make_lr(C: float):
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(C=C, class_weight="balanced", max_iter=MAX_ITER),
    )


def _fit(model, X, y, notes: dict):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(X, y)
    notes["convergence_warnings"] = notes.get("convergence_warnings", 0) + sum(
        issubclass(w.category, ConvergenceWarning) for w in caught
    )
    return model


def oof_scores(X, y, folds, C: float, notes: dict) -> np.ndarray:
    out = np.empty(len(y))
    for train, test in folds:
        out[test] = _fit(make_lr(C), X[train], y[train], notes).decision_function(X[test])
    return out


def youden_threshold(y, score) -> float:
    """The score with the highest recall - FPR when calls are score >= threshold."""
    y = np.asarray(y, dtype=bool)
    score = np.asarray(score, dtype=np.float64)
    order = np.argsort(-score, kind="stable")
    s = score[order]
    ends = np.r_[np.nonzero(s[1:] != s[:-1])[0], len(s) - 1]
    tp = np.cumsum(y[order])[ends] / y.sum()
    fp = np.cumsum(~y[order])[ends] / (~y).sum()
    return float(s[ends][int(np.argmax(tp - fp))])


def fit_platt(y, score) -> tuple[float, float]:
    """Platt scaling: p = 1 / (1 + exp(-(a * score + b))), fitted without class weights."""
    lr = LogisticRegression(C=1e6, max_iter=MAX_ITER)
    lr.fit(np.asarray(score, dtype=np.float64)[:, None], np.asarray(y, dtype=bool))
    return float(lr.coef_[0, 0]), float(lr.intercept_[0])


def platt(score, a: float, b: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-(a * np.asarray(score, dtype=np.float64) + b)))


def _score(model, u, candidate, idx, h_variant) -> np.ndarray:
    out = np.empty(len(idx))
    for start in range(0, len(idx), SCORE_CHUNK):
        part = idx[start : start + SCORE_CHUNK]
        out[start : start + len(part)] = model.decision_function(
            universe.features(u, candidate, part, h_variant)
        )
    return out


def fit_unit(u, train_idx, y, groups, score_idx, seed: int, candidates=CANDIDATES):
    """Fit every candidate on the training rows; score the rows score_idx.

    Return (params, scores): params[candidate] is a JSON-ready dict; scores[candidate] is a
    dict with `score` and `prob` (arrays, or None for rules) and `call` (bool array)."""
    train_idx = np.asarray(train_idx, dtype=np.int64)
    score_idx = np.asarray(score_idx, dtype=np.int64)
    y = np.asarray(y, dtype=bool)
    folds = inner_folds(y, groups, seed)
    params, scores = {}, {}
    base = {"n_train": int(len(y)), "n_train_pos": int(y.sum())}
    for cand in candidates:
        if cand in rules.RULES:
            fit = rules.fit_rule(cand, y, u.sp[train_idx], u.rank[train_idx], u.st[train_idx])
            call = rules.rule_call(
                cand, u.sp[score_idx], u.rank[score_idx], u.st[score_idx], fit["g"], fit["t"]
            )
            params[cand] = {**base, **fit}
            scores[cand] = {"score": None, "prob": None, "call": call}
            continue
        notes: dict = {}
        settings = (
            [(v, c) for v in H_VARIANTS for c in C_GRID]
            if cand == "H"
            else [(None, c) for c in C_GRID]
        )
        best, inner = None, {}
        cache = {}
        for variant, C in settings:
            if variant not in cache:
                cache[variant] = universe.features(u, cand, train_idx, variant)
            oof = oof_scores(cache[variant], y, folds, C, notes)
            ap = float(pr_auc(np.ones(len(y)), y, oof)[0])
            inner[f"{variant}:{C}" if variant else f"{C}"] = ap
            if best is None or ap > best[0]:
                best = (ap, variant, C, oof)
        ap, variant, C, oof = best
        threshold = youden_threshold(y, oof)
        a, b = fit_platt(y, oof)
        model = _fit(make_lr(C), cache[variant], y, notes)
        s = _score(model, u, cand, score_idx, variant)
        params[cand] = {
            **base,
            "C": C,
            "h_variant": variant,
            "threshold": threshold,
            "platt_a": a,
            "platt_b": b,
            "inner_pr_auc": inner,
            "convergence_warnings": notes.get("convergence_warnings", 0),
        }
        scores[cand] = {"score": s, "prob": platt(s, a, b), "call": s >= threshold}
    return params, scores


# --- work units for 10_fit_and_score.py (module level, so worker processes can import them) ---

_WORKER: dict = {}


def init_worker(phaseb: str) -> None:
    """Process-pool initializer: load the universe once per worker (10 verified the files)."""
    _WORKER["u"] = universe.load(phaseb, verify=False)


def run_task(task: dict, u=None):
    """Fit one outer training set. task: key, train, y, groups, score, seed, candidates."""
    u = u if u is not None else _WORKER["u"]
    try:
        params, scores = fit_unit(
            u,
            u.rows(task["train"]),
            task["y"],
            task["groups"],
            u.rows(task["score"]),
            task["seed"],
            task["candidates"],
        )
    except ModelError as exc:
        raise ModelError(f"{task['key']}: {exc}") from exc
    return task["key"], params, scores
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_models.py -q`
Expected: `9 passed`.

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/universe.py analysis/step1_compare/phasec/models.py \
  tests/step1_compare/test_phasec_models.py
git commit -m "phasec: universe and nested protocol for B0, B1, M8, M35, M8-C, M35-C, H (spec 3.3)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 7: Clusters, splits and T-c removals (09)

**Files:**
- Create: `analysis/step1_compare/phasec/splits.py`, `analysis/step1_compare/phasec/09_make_splits.py`, `tests/step1_compare/fixtures/phasec/stub_mmseqs/mmseqs` (executable)
- Modify: `tests/step1_compare/phasec_fixture.py` (append the chain helpers)
- Test: `tests/step1_compare/test_phasec_splits.py`

**Interfaces:**
- Consumes: `evalio.*` (Task 1); 08 outputs and `build_run.json` (Task 2); sklearn `StratifiedGroupKFold`.
- Produces: `splits.S1_FOLDS = 5`, `LITERATURE_CLADE = "Eurotiomycetes"`, `PARTS`, `TRAIN_PARTS`, `TEST_PARTS`, `MEMBER_COLUMNS`, `REMOVED_COLUMNS`, `IDENTITY_COLUMNS`, `CLUSTER_COLUMNS`, `IDENTITY_CUT = 0.3`, `SplitError(ValueError)`, `SplitDef` (frozen dataclass: `split_id`, `train_sources`, `test_sources`, `test_taxa`, `test_clades`, `literature`, `folds`), `read_cluster_tsv(lines) -> dict`, `check_clusters(cluster_of, hashes)`, `split_definitions(species_rows) -> list[SplitDef]`, `s1_folds(table, cluster_of, train_sources, seed)`, `tc_removals(tc_rows, test_rows, cluster_of, taxa, clades) -> (kept, [(row, rule)])`, `build(table, literature, cluster_of, defs, seed) -> (members, removed)`, `check_test_truth(members)`, `check_no_shared_hash(members)`, `check_no_cluster_spans(members)`, `parse_hits(lines) -> dict[str, float]`, `identity_rows(split_id, test_hashes, best)`. Script `09_make_splits.py`: `OUTPUT_NAMES`, `CLUSTER_ARGS`, `SEARCH_ARGS`, `mmseqs_version(mmseqs)`, `run(work, species_path, mmseqs, tmp, threads, arguments=())`, `main(argv=None)`; options `--work-dir`, `--species`, `--mmseqs` (default `$STEP1_MMSEQS` or `mmseqs`), `--tmp-dir` (required), `--threads`. Fixture helpers: `phasec_fixture.STUB_MMSEQS`, `run_chain(root, upto, load_phasec, workers=1, candidates=None) -> dict`, `copy_work(chain, dest) -> dict`, `eval_argv(fx, n_resamples=50)`, `fix_recorded_hash(run_json, name, path, key="outputs_sha256")`.

Notes found while prototyping (keep them):
- `module load MMseqs2/17-b804f` gives the AVX2 build; on c01 (Opteron 6376) `mmseqs version` exits with 132. `mmseqs_version` names this case in its STOP message. The job runs on partition `epyc` (Task 11).
- On the real data the prototype run (non-AVX2 binary, 2 threads of c01) took 18 min 48 s and gave 9,737 clusters and the removal counts in the Facts table.
- MMseqs2 `easy-search` keeps its default `--max-seqs 300` and E-value 1e-3; `max_identity` is the highest `fident` among the reported hits. Not verified: whether a better hit is lost behind the prefilter limit for some query.

- [ ] **Step 1: Write the MMseqs2 stub (test only)**

```bash
mkdir -p tests/step1_compare/fixtures/phasec/stub_mmseqs
```

`tests/step1_compare/fixtures/phasec/stub_mmseqs/mmseqs`:

```python
#!/usr/bin/env python3
"""Test stub of MMseqs2 (runs with Python 3.6 or later): version, easy-cluster, easy-search.

easy-cluster: sequences with the same first 6 residues form one cluster; the representative is
the smallest id. easy-search: a query hits every target with the same first 6 residues;
fident = identical positions / longer length. STUB_MMSEQS_FAIL=<command> makes it exit 1;
STUB_MMSEQS_LOG=<file> appends each argument list to that file.
"""
import os
import sys


def read_fasta(path):
    seqs, head = {}, None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                head = line[1:].split()[0]
                seqs[head] = ""
            elif head:
                seqs[head] += line
    return seqs


def main(argv):
    if os.environ.get("STUB_MMSEQS_LOG"):
        with open(os.environ["STUB_MMSEQS_LOG"], "a") as log:
            log.write("\t".join(argv) + "\n")
    cmd = argv[0]
    if cmd == "version":
        print("stub-17")
        return 0
    if os.environ.get("STUB_MMSEQS_FAIL") == cmd:
        print("stub failure", file=sys.stderr)
        return 1
    if cmd == "easy-cluster":
        seqs = read_fasta(argv[1])
        groups = {}
        for h, s in seqs.items():
            groups.setdefault(s[:6], []).append(h)
        with open(argv[2] + "_cluster.tsv", "w") as out:
            for members in groups.values():
                rep = min(members)
                for m in sorted(members):
                    out.write(rep + "\t" + m + "\n")
        return 0
    if cmd == "easy-search":
        query, target = read_fasta(argv[1]), read_fasta(argv[2])
        with open(argv[3], "w") as out:
            for q, qs in sorted(query.items()):
                for t, ts in sorted(target.items()):
                    if qs[:6] == ts[:6]:
                        same = sum(a == b for a, b in zip(qs, ts))
                        out.write("%s\t%s\t%.3f\n" % (q, t, same / max(len(qs), len(ts))))
        return 0
    print("unknown command " + cmd, file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
```

Run: `chmod +x tests/step1_compare/fixtures/phasec/stub_mmseqs/mmseqs`

- [ ] **Step 2: Append the chain helpers to `tests/step1_compare/phasec_fixture.py`**

Append this text at the end of the file:

```python
STUB_MMSEQS = Path(__file__).resolve().parent / "fixtures" / "phasec" / "stub_mmseqs" / "mmseqs"


def run_chain(root: Path, upto: int, load_phasec, workers: int = 1, candidates=None) -> dict:
    """Build the fixture and run 08 (and 09 with the stub MMseqs2, 10, 11) up to step `upto`."""
    fx = make_work(root)
    assert load_phasec("08_build_eval_tables").main(build_argv(fx)) == 0
    if upto >= 9:
        argv = ["--work-dir", str(fx["work"]), "--species", str(fx["species"]),
                "--mmseqs", str(STUB_MMSEQS), "--tmp-dir", str(Path(root) / "scratch")]  # fmt: skip
        assert load_phasec("09_make_splits").main(argv) == 0
    if upto >= 10:
        argv = ["--work-dir", str(fx["work"]), "--workers", str(workers)]
        if candidates:
            argv += ["--candidates", candidates]
        assert load_phasec("10_fit_and_score").main(argv) == 0
    if upto >= 11:
        argv = ["--work-dir", str(fx["work"]), "--sets", str(fx["sets"]),
                "--species", str(fx["species"]), "--n-resamples", "50"]  # fmt: skip
        assert load_phasec("11_evaluate").main(argv) == 0
    return fx


def copy_work(chain: dict, dest: Path) -> dict:
    """A writable copy of a chain's work directory (and its side files)."""
    import shutil

    dest = Path(dest)
    shutil.copytree(chain["work"], dest / "work", symlinks=True)
    return {**chain, "work": dest / "work"}


def eval_argv(fx: dict, n_resamples: int = 50) -> list[str]:
    return ["--work-dir", str(fx["work"]), "--sets", str(fx["sets"]), "--species",
            str(fx["species"]), "--n-resamples", str(n_resamples)]  # fmt: skip


def fix_recorded_hash(run_json: Path, name: str, path: Path, key: str = "outputs_sha256") -> None:
    """Record the current SHA-256 of `path` under run_json[key][name] (for tamper tests)."""
    obj = json.loads(Path(run_json).read_text())
    obj[key][name] = _sha(path)
    Path(run_json).write_text(json.dumps(obj))
```

`run_chain` calls 10 and 11 only when `upto` is 10 or more; Tasks 8 and 9 add those scripts.

- [ ] **Step 3: Write the failing test**

`tests/step1_compare/test_phasec_splits.py`:

```python
"""splits.py and 09_make_splits.py: clusters, S1-S3, T-c removals, max identity (spec 3.4)."""

import json
import os
import shutil
import subprocess

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import paths  # noqa: E402
import phasec_fixture as pf  # noqa: E402
import splits  # noqa: E402
import truth_table  # noqa: E402
from conftest import load_phasec  # noqa: E402

STUB = paths.STEP1_DIR.parents[1] / "tests/step1_compare/fixtures/phasec/stub_mmseqs/mmseqs"
SPECIES = truth_table.read_tsv(paths.STEP1_DIR / "species.tsv")


def go(h, cls, source, gene=None, **kw):
    return {"seq_sha256": h, "origin": "go", "class": cls, "source_ids": source,
            "gene_ids": gene or f"g_{h}", "taxon_ids": "", "clades": "", **kw}  # fmt: skip


def tc(h, acc, taxon="1", clade="Saccharomycotina"):
    return {"seq_sha256": h, "origin": "tc", "class": "pos", "source_ids": "T-c",
            "gene_ids": acc, "taxon_ids": taxon, "clades": clade}  # fmt: skip


def test_read_cluster_tsv():
    got = splits.read_cluster_tsv(["a\ta\n", "a\tb\n", "c\tc\n", "\n"])
    assert got == {"a": "a", "b": "a", "c": "c"}
    with pytest.raises(splits.SplitError, match="in two clusters"):
        splits.read_cluster_tsv(["a\ta", "c\ta"])
    with pytest.raises(splits.SplitError, match="expected 2 fields"):
        splits.read_cluster_tsv(["a b"])
    with pytest.raises(splits.SplitError, match="1 sequences without a cluster"):
        splits.check_clusters({"a": "a"}, ["a", "z"])


def test_split_definitions_from_species_tsv():
    defs = {d.split_id: d for d in splits.split_definitions(SPECIES)}
    assert list(defs) == ["S1", "S2-Calb_CGD", "S2-Scer_SGD", "S2-Spom_PomBase",
                          "S3-Basidiomycota", "S3-Eurotiomycetes", "FULL"]  # fmt: skip
    assert defs["S2-Calb_CGD"].train_sources == ("Scer_SGD",)
    assert defs["S2-Calb_CGD"].test_taxa == frozenset({"237561"})
    assert defs["S3-Eurotiomycetes"].test_sources == ("Afum_ASPFU", "Anid_EMENI")
    assert defs["S3-Eurotiomycetes"].literature and not defs["S3-Basidiomycota"].literature
    assert defs["S3-Basidiomycota"].test_sources == ("Cneo_H99_GOA", "Umay_MYCMD")
    tested = {s for d in defs.values() for s in d.test_sources}
    assert not tested & {"Spom_SCHPO-mod", "Cneo_JEC21_GOA", "Cneo_CRYD1"}  # alternate files


def _toy(n_per_class=12):
    table, cluster_of = [], {}
    for src in ("Scer_SGD", "Calb_CGD", "Spom_PomBase", "Afum_ASPFU"):
        for i in range(n_per_class):
            for cls in ("pos", "neg"):
                h = f"{src[:4]}{cls}{i:02d}"
                table.append(go(h, cls, src))
                cluster_of[h] = f"cl_{h}" if i % 4 else f"fam{i % 3}_{cls}"
        h = f"{src[:4]}excl"
        table.append(go(h, "excluded", src))
        cluster_of[h] = f"cl_{h}"
    for i in range(6):
        h = f"tc{i:02d}"
        table.append(tc(h, f"Q{i}"))
        cluster_of[h] = f"cl_{h}"
    return table, cluster_of


def _species():
    return [
        {"source_id": "Scer_SGD", "taxon_id": "559292", "in_clade": "Saccharomycotina", "role": "train"},
        {"source_id": "Calb_CGD", "taxon_id": "237561", "in_clade": "Saccharomycotina", "role": "train"},
        {"source_id": "Spom_PomBase", "taxon_id": "284812", "in_clade": "Taphrinomycotina", "role": "test_species"},
        {"source_id": "Afum_ASPFU", "taxon_id": "330879", "in_clade": "Eurotiomycetes", "role": "test_clade"},
    ]  # fmt: skip


def _build(table, cluster_of, lit=()):
    return splits.build(table, list(lit), cluster_of, splits.split_definitions(_species()), 1)


def test_s1_folds_keep_clusters_whole_and_cover_tc_rows():
    table, cluster_of = _toy()
    members, _ = _build(table, cluster_of)
    s1 = [m for m in members if m["split_id"] == "S1"]
    fold_of = {}
    for m in s1:
        if m["part"] in ("test", "test_tc"):
            assert m["seq_sha256"] not in fold_of
            fold_of[m["seq_sha256"]] = m["fold"]
    by_cluster = {}
    for h, f in fold_of.items():
        by_cluster.setdefault(cluster_of[h], set()).add(f)
    assert all(len(v) == 1 for v in by_cluster.values())
    assert {f"tc{i:02d}" for i in range(6)} <= set(fold_of)  # every T-c row has an oof fold
    assert "Scerexcl" in fold_of and "Spompos00" not in fold_of
    splits.check_no_cluster_spans(members)
    splits.check_no_shared_hash(members)


def test_no_cluster_spans_train_and_test():
    table, cluster_of = _toy()
    members, _ = _build(table, cluster_of)
    s1 = [dict(m) for m in members if m["split_id"] == "S1" and m["fold"] == "0"]
    test = next(m for m in s1 if m["part"] == "test")
    train = next(m for m in s1 if m["part"] == "train")
    train["cluster_id"] = test["cluster_id"]  # plant one spanning cluster
    with pytest.raises(splits.SplitError, match="spans train and test"):
        splits.check_no_cluster_spans(s1)


def test_s2_s3_allow_go_homologs_but_not_tc_cluster_mates():
    table, cluster_of = _toy()
    members, _ = _build(table, cluster_of)
    s3 = [m for m in members if m["split_id"] == "S3-Eurotiomycetes"]
    # fam clusters hold Scer, Calb and Afum rows: GO train and test share them by design
    train_c = {m["cluster_id"] for m in s3 if m["part"] == "train"}
    test_c = {m["cluster_id"] for m in s3 if m["part"] == "test"}
    assert train_c & test_c
    splits.check_no_cluster_spans(members)  # does not raise
    bad = [dict(m) for m in s3]
    next(m for m in bad if m["part"] == "train_tc")["cluster_id"] = next(iter(test_c))
    with pytest.raises(splits.SplitError, match="spans train and test"):
        splits.check_no_cluster_spans(bad)


def test_test_truth_sources():
    table, cluster_of = _toy()
    members, _ = _build(table, cluster_of)
    splits.check_test_truth(members)
    assert not [m for m in members if m["origin"] == "tc" and m["part"] in ("test", "test_lit")]
    routed = [dict(m) for m in members]
    next(m for m in routed if m["origin"] == "tc")["part"] = "test"
    with pytest.raises(splits.SplitError, match="T-c row .* is in a test set"):
        splits.check_test_truth(routed)


def test_tc_excludes_test_proteins():
    table, cluster_of = _toy()
    table.append(tc("Afumpos00", "Q99"))  # same hash as an A. fumigatus test protein
    table.append(tc("tcacc", "g_Afumneg03"))  # accession of a test protein
    cluster_of["tcacc"] = "cl_tcacc"
    # the hash row would be dropped by 08 (GO label wins); here it tests rule a directly
    kept, removed = splits.tc_removals(
        [r for r in table if r["origin"] == "tc"],
        [r for r in table if r["origin"] == "go" and r["source_ids"] == "Afum_ASPFU"],
        cluster_of,
    )
    rules = {r["gene_ids"]: rule for r, rule in removed}
    assert rules == {"Q99": "a_test_protein", "g_Afumneg03": "a_test_protein"}
    assert len(kept) == 6


def test_tc_excludes_test_cluster_mates():
    table, cluster_of = _toy()
    cluster_of["tc00"] = cluster_of["Spompos05"]  # T-c row in the cluster of a S. pombe protein
    members, removed = _build(table, cluster_of)
    got = [(r["split_id"], r["rule"]) for r in removed if r["seq_sha256"] == "tc00"]
    assert got == [("S2-Spom_PomBase", "b_cluster_mate")]
    assert any(m["seq_sha256"] == "tc00" and m["part"] == "train_tc" for m in members
               if m["split_id"] == "S2-Calb_CGD")  # fmt: skip


def test_tc_excludes_test_taxa():
    table, cluster_of = _toy()
    for h, taxon, clade in (
        ("tcalb", "237561", "Saccharomycotina"),
        ("tcfum", "330879", "Eurotiomycetes"),
    ):
        table.append(tc(h, f"A_{h}", taxon, clade))
        cluster_of[h] = f"cl_{h}"
    members, removed = _build(table, cluster_of)
    got = {(r["split_id"], r["seq_sha256"]): r["rule"] for r in removed}
    assert got[("S2-Calb_CGD", "tcalb")] == "c_test_taxon"  # S2: taxon of the test species
    assert got[("S3-Eurotiomycetes", "tcfum")] == "c_test_taxon"  # S3: taxon in the test clade
    assert ("S2-Scer_SGD", "tcalb") not in got and ("S2-Calb_CGD", "tcfum") not in got
    s1_tc = {m["seq_sha256"] for m in members if m["split_id"] == "S1" and m["part"] == "train_tc"}
    assert {"tcalb", "tcfum"} <= s1_tc  # S1 keeps them (spec 4, context note)


def test_hash_in_a_training_and_a_test_source_stops():
    # Review Focus 4: one sequence labelled in S. cerevisiae and in S. pombe
    table, cluster_of = _toy()
    table.append(go("both", "neg", "Scer_SGD,Spom_PomBase"))
    cluster_of["both"] = "cl_both"
    with pytest.raises(splits.SplitError, match="belongs to a training source and a test source"):
        _build(table, cluster_of)


def test_literature_rows_are_tested_in_s3_eurotiomycetes_only():
    table, cluster_of = _toy()
    lit = [{"accession": "P1", "seq_sha256": "lit1", "literature_positive": "yes"},
           {"accession": "P2", "seq_sha256": "lit2", "literature_positive": "no"},
           {"accession": "P3", "seq_sha256": "", "literature_positive": "no"}]  # fmt: skip
    cluster_of.update({"lit1": "cl_lit1", "lit2": "cl_lit2"})
    members, _ = _build(table, cluster_of, lit)
    got = {(m["split_id"], m["seq_sha256"], m["class"]) for m in members if m["part"] == "test_lit"}
    assert got == {("S3-Eurotiomycetes", "lit1", "pos"), ("S3-Eurotiomycetes", "lit2", "excluded")}


def test_max_identity_excludes_self_hits():
    lines = ["q1\tq1\t1.000", "q1\tt1\t0.420", "q1\tt2\t0.610", "q2\tt1\t0.250", "q3\tq3\t1.0"]
    best = splits.parse_hits(lines)
    assert best == {"q1": 0.61, "q2": 0.25}
    rows = {r["seq_sha256"]: r for r in splits.identity_rows("S2-x", ["q1", "q2", "q3"], best)}
    assert rows["q1"]["max_identity"] == "0.6100" and rows["q1"]["below_0.3"] == "no"
    assert rows["q2"]["below_0.3"] == "yes"
    assert rows["q3"]["max_identity"] == "" and rows["q3"]["below_0.3"] == "yes"  # self hit only


# --- 09_make_splits.py with the stub MMseqs2 ---


def _built(tmp_path):
    fx = pf.make_work(tmp_path)
    assert load_phasec("08_build_eval_tables").main(pf.build_argv(fx)) == 0
    return fx


def _run09(fx, tmp_path, mmseqs=STUB):
    argv = ["--work-dir", str(fx["work"]), "--species", str(fx["species"]),
            "--mmseqs", str(mmseqs), "--tmp-dir", str(tmp_path / "scratch")]  # fmt: skip
    return load_phasec("09_make_splits").main(argv)


def test_make_splits_with_the_stub(tmp_path, monkeypatch):
    fx = _built(tmp_path)
    log_file = tmp_path / "stub.log"
    monkeypatch.setenv("STUB_MMSEQS_LOG", str(log_file))
    assert _run09(fx, tmp_path) == 0
    out = fx["work"] / "phasec"
    run = json.loads((out / "splits_run.json").read_text())
    assert run["mmseqs_version"] == "stub-17" and run["seed"] == 20261001
    calls = [line.split("\t") for line in log_file.read_text().splitlines()]
    cluster = next(c for c in calls if c[0] == "easy-cluster")
    assert cluster[4:10] == ["--min-seq-id", "0.3", "-c", "0.5", "--cov-mode", "0"]
    search = [c for c in calls if c[0] == "easy-search"]
    assert len(search) == 5  # S2-Calb_CGD, S2-Scer_SGD, S2-Spom_PomBase, S3-Basidio, S3-Euro
    assert search[0][5:13] == ["-s", "7.5", "-c", "0.5", "--cov-mode", "0", "--format-output",
                               "query,target,fident"]  # fmt: skip
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    assert {m["split_id"] for m in members} == set(run["splits"])
    s1_test = {m["seq_sha256"] for m in members if m["split_id"] == "S1" and m["part"] == "test"}
    table = truth_table.read_tsv(out / "eval_table.tsv.gz")
    train_go = {r["seq_sha256"] for r in table if r["origin"] == "go" and r["roles"] == "train"}
    assert s1_test == train_go  # every training-source GO row is scored out of fold once
    ident = truth_table.read_tsv(out / "max_identity.tsv.gz")
    assert {r["split_id"] for r in ident} == {
        "S2-Calb_CGD",
        "S2-Scer_SGD",
        "S2-Spom_PomBase",
        "S3-Basidiomycota",
        "S3-Eurotiomycetes",
    }


def test_make_splits_is_deterministic(tmp_path):
    fx = _built(tmp_path)
    assert _run09(fx, tmp_path) == 0
    out = fx["work"] / "phasec"
    first = {n: (out / n).read_bytes() for n in ("clusters.tsv.gz", "split_members.tsv.gz")}
    assert _run09(fx, tmp_path) == 0
    assert first == {n: (out / n).read_bytes() for n in first}


def test_stale_input_stops(tmp_path, capsys):
    fx = _built(tmp_path)
    table = fx["work"] / "phasec" / "eval_table.tsv.gz"
    table.write_bytes(table.read_bytes() + b"\n")
    assert _run09(fx, tmp_path) == 2
    assert "eval_table.tsv.gz differs from the SHA-256 in build_run.json" in capsys.readouterr().err
    assert not (fx["work"] / "phasec" / "splits_run.json").exists()


def test_mmseqs_failure_stops_without_output(tmp_path, monkeypatch, capsys):
    fx = _built(tmp_path)
    monkeypatch.setenv("STUB_MMSEQS_FAIL", "easy-search")
    assert _run09(fx, tmp_path) == 2
    assert "easy-search exited with 1" in capsys.readouterr().err
    assert not (fx["work"] / "phasec" / "split_members.tsv.gz").exists()


def test_mmseqs_version_failure_names_avx2(tmp_path, capsys):
    fx = _built(tmp_path)
    bad = tmp_path / "mmseqs"
    bad.write_text("#!/bin/bash\nexit 132\n")
    bad.chmod(0o755)
    assert _run09(fx, tmp_path, bad) == 2
    assert "AVX2" in capsys.readouterr().err


def _real_mmseqs():
    exe = os.environ.get("STEP1_MMSEQS") or shutil.which("mmseqs")
    if not exe:
        return None
    try:
        ok = subprocess.run([exe, "version"], capture_output=True).returncode == 0
    except OSError:
        return None
    return exe if ok else None


@pytest.mark.skipif(
    _real_mmseqs() is None, reason="needs a working MMseqs2 (module or STEP1_MMSEQS)"
)
def test_make_splits_with_real_mmseqs(tmp_path):
    fx = _built(tmp_path)
    assert _run09(fx, tmp_path, _real_mmseqs()) == 0
    out = fx["work"] / "phasec"
    clusters = truth_table.read_tsv(out / "clusters.tsv.gz")
    seqs = pf.gz_lines(out / "eval_sequences.fasta.gz")
    assert len(clusters) == sum(x.startswith(">") for x in seqs)
    run = json.loads((out / "splits_run.json").read_text())
    assert run["mmseqs_version"] and run["clusters"] >= 1
```

- [ ] **Step 4: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_splits.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'splits'`.

- [ ] **Step 5: Write `analysis/step1_compare/phasec/splits.py`**

```python
"""Clusters, split schemes S1, S2, S3, the T-c leakage controls (Phase C spec 3.4 items 3-6).

Split ids:
- S1: homology-grouped 5-fold CV over the GO rows of the training sources (role `train`).
  StratifiedGroupKFold(5, shuffle, seed) on the pos and neg rows, groups = clusters. Every
  other row of the training sources and every T-c row takes the fold of its cluster; a cluster
  without a fold gets int(sha256(cluster_id)[:16], 16) mod 5. So T-c rows get out-of-fold
  scores too (spec 4, score_source).
- S2-<source>: leave-species-out. For each training source B: train on the other training
  sources, test B. For each `test_species` source: train on all training sources.
- S3-<clade>: leave-clade-out. Train on all training sources; test every `test_clade` or
  `undecided` source of the clade. S3-Eurotiomycetes also tests the literature rows.
- FULL: train on all training sources and all T-c rows; no test. Its models score the
  proteins that are in no training table (score_source `final`).

Parts: train (GO pos/neg training rows), train_tc (T-c rows, V-kw only), test (GO rows of the
test set, all classes), test_tc (S1: T-c rows of the test fold; scored, never truth), test_lit
(literature rows). T-c rows are removed from train_tc, first match wins:
a_test_protein (accession or hash of a test protein), b_cluster_mate (cluster holds a test
protein, ruling C-5), c_test_taxon (S2: taxon_id of the test species; S3: taxon in the test
clade, ruling C-7).

In S2 and S3 a GO training protein may share a cluster with a test protein: that is homology
across species, which the maximum-identity stratum measures (ruling C-4). Only S1 folds and
the T-c rows must not share clusters with the test set.
"""

import hashlib
from dataclasses import dataclass, field

import numpy as np
from sklearn.model_selection import StratifiedGroupKFold

S1_FOLDS = 5
LITERATURE_CLADE = "Eurotiomycetes"
PARTS = ("train", "train_tc", "test", "test_tc", "test_lit")
TRAIN_PARTS = ("train", "train_tc")
TEST_PARTS = ("test", "test_tc", "test_lit")
MEMBER_COLUMNS = ("split_id", "fold", "seq_sha256", "part", "origin", "class", "cluster_id")
REMOVED_COLUMNS = ("split_id", "fold", "seq_sha256", "gene_ids", "taxon_ids", "rule")
IDENTITY_COLUMNS = ("split_id", "seq_sha256", "max_identity", "below_0.3")
CLUSTER_COLUMNS = ("seq_sha256", "cluster_id")
IDENTITY_CUT = 0.3


class SplitError(ValueError):
    """A split breaks a leakage rule or the cluster table does not fit the sequences."""


@dataclass(frozen=True)
class SplitDef:
    split_id: str
    train_sources: tuple
    test_sources: tuple = ()
    test_taxa: frozenset = field(default_factory=frozenset)
    test_clades: frozenset = field(default_factory=frozenset)
    literature: bool = False
    folds: int = 1


def read_cluster_tsv(lines) -> dict[str, str]:
    """MMseqs2 `<prefix>_cluster.tsv` lines (representative TAB member) -> member: cluster."""
    out = {}
    for n, line in enumerate(lines, 1):
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            raise SplitError(f"cluster table line {n}: expected 2 fields, got {len(parts)}")
        rep, member = parts
        if member in out:
            raise SplitError(f"cluster table: {member} is in two clusters")
        out[member] = rep
    return out


def check_clusters(cluster_of: dict, hashes) -> None:
    hashes = set(hashes)
    missing = sorted(hashes - set(cluster_of))
    extra = sorted(set(cluster_of) - hashes)
    if missing or extra:
        raise SplitError(
            f"cluster table: {len(missing)} sequences without a cluster, {len(extra)} unknown "
            f"members (for example {(missing or extra)[0]})"
        )


def split_definitions(species_rows) -> list[SplitDef]:
    train = tuple(sorted(r["source_id"] for r in species_rows if r["role"] == "train"))
    taxon = {r["source_id"]: r["taxon_id"] for r in species_rows}
    defs = [SplitDef("S1", train, folds=S1_FOLDS)]
    if len(train) > 1:
        for b in train:
            rest = tuple(s for s in train if s != b)
            defs.append(SplitDef(f"S2-{b}", rest, (b,), frozenset({taxon[b]})))
    for r in sorted(species_rows, key=lambda r: r["source_id"]):
        if r["role"] == "test_species":
            s = r["source_id"]
            defs.append(SplitDef(f"S2-{s}", train, (s,), frozenset({taxon[s]})))
    clades: dict[str, list[str]] = {}
    for r in species_rows:
        if r["role"] in ("test_clade", "undecided"):
            clades.setdefault(r["in_clade"], []).append(r["source_id"])
    for clade in sorted(clades):
        defs.append(
            SplitDef(
                f"S3-{clade}",
                train,
                tuple(sorted(clades[clade])),
                test_clades=frozenset({clade}),
                literature=clade == LITERATURE_CLADE,
            )
        )
    defs.append(SplitDef("FULL", train))
    return defs


def _sources(row) -> set[str]:
    return set(row["source_ids"].split(","))


def _hash_fold(cluster_id: str, n_folds: int) -> int:
    return int(hashlib.sha256(cluster_id.encode()).hexdigest()[:16], 16) % n_folds


def s1_folds(table, cluster_of: dict, train_sources, seed: int, n_folds=S1_FOLDS):
    """Fold of every S1 row (GO rows of the training sources, T-c rows), by cluster."""
    train_sources = set(train_sources)
    fit = sorted(
        (r for r in table if r["origin"] == "go" and r["class"] in ("pos", "neg")
         and _sources(r) <= train_sources),
        key=lambda r: r["seq_sha256"],
    )  # fmt: skip
    y = np.array([r["class"] == "pos" for r in fit])
    groups = np.array([cluster_of[r["seq_sha256"]] for r in fit])
    cv = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    fold_of_cluster = {}
    for k, (_, test) in enumerate(cv.split(np.zeros(len(fit)), y, groups)):
        for c in groups[test]:
            fold_of_cluster[c] = k
    out = {}
    for r in table:
        in_s1 = r["origin"] == "tc" or (r["origin"] == "go" and _sources(r) <= train_sources)
        if not in_s1:
            continue
        c = cluster_of[r["seq_sha256"]]
        out[r["seq_sha256"]] = fold_of_cluster.get(c, _hash_fold(c, n_folds))
    return out


def tc_removals(tc_rows, test_rows, cluster_of: dict, taxa=frozenset(), clades=frozenset()):
    """Split T-c rows into (kept, [(row, rule)]) for one test set (rules a, b, c)."""
    test_acc = {a for r in test_rows for a in r["gene_ids"].split(",") if a}
    test_hash = {r["seq_sha256"] for r in test_rows}
    test_clusters = {cluster_of[r["seq_sha256"]] for r in test_rows}
    kept, removed = [], []
    for r in tc_rows:
        accs = set(r["gene_ids"].split(","))
        if accs & test_acc or r["seq_sha256"] in test_hash:
            removed.append((r, "a_test_protein"))
        elif cluster_of[r["seq_sha256"]] in test_clusters:
            removed.append((r, "b_cluster_mate"))
        elif set(r["taxon_ids"].split(",")) & taxa or set(r["clades"].split(",")) & clades:
            removed.append((r, "c_test_taxon"))
        else:
            kept.append(r)
    return kept, removed


def _member(split_id, fold, row, part, cls, cluster_of):
    return {
        "split_id": split_id,
        "fold": str(fold),
        "seq_sha256": row["seq_sha256"],
        "part": part,
        "origin": row.get("origin", "lit"),
        "class": cls,
        "cluster_id": cluster_of[row["seq_sha256"]],
    }


def _lit_as_rows(literature):
    rows = []
    for r in literature:
        if r["seq_sha256"]:
            rows.append(
                {
                    "seq_sha256": r["seq_sha256"],
                    "origin": "lit",
                    "class": "pos" if r["literature_positive"] == "yes" else "excluded",
                    "gene_ids": r["accession"],
                }
            )
    return rows


def build(table, literature, cluster_of: dict, defs, seed: int):
    """Return (members, removed) for every split definition."""
    go = [r for r in table if r["origin"] == "go"]
    tc = [r for r in table if r["origin"] == "tc"]
    lit = _lit_as_rows(literature)
    members, removed = [], []
    for d in defs:
        train_sources = set(d.train_sources)
        if d.split_id == "S1":
            fold = s1_folds(table, cluster_of, d.train_sources, seed, d.folds)
            s1_go = [r for r in go if _sources(r) <= train_sources]
            for f in range(d.folds):
                test = [r for r in s1_go if fold[r["seq_sha256"]] == f]
                train = [
                    r for r in s1_go if fold[r["seq_sha256"]] != f and r["class"] != "excluded"
                ]
                tc_test = [r for r in tc if fold[r["seq_sha256"]] == f]
                tc_kept, tc_out = tc_removals(
                    [r for r in tc if fold[r["seq_sha256"]] != f], test, cluster_of
                )
                members += _parts(d.split_id, f, train, tc_kept, test, tc_test, [], cluster_of)
                removed += _removed(d.split_id, f, tc_out)
            continue
        if d.split_id == "FULL":
            train = [r for r in go if _sources(r) <= train_sources and r["class"] != "excluded"]
            members += _parts("FULL", 0, train, tc, [], [], [], cluster_of)
            continue
        test_sources = set(d.test_sources)
        for r in go:
            if _sources(r) & test_sources and _sources(r) & train_sources:
                raise SplitError(
                    f"{d.split_id}: sequence {r['seq_sha256']} belongs to a training source and "
                    f"a test source ({r['source_ids']})"
                )
        train = [r for r in go if _sources(r) <= train_sources and r["class"] != "excluded"]
        test = [r for r in go if _sources(r) & test_sources]
        test_lit = lit if d.literature else []
        tc_kept, tc_out = tc_removals(tc, test + test_lit, cluster_of, d.test_taxa, d.test_clades)
        members += _parts(d.split_id, 0, train, tc_kept, test, [], test_lit, cluster_of)
        removed += _removed(d.split_id, 0, tc_out)
    return members, removed


def _parts(split_id, fold, train, train_tc, test, test_tc, test_lit, cluster_of):
    out = [_member(split_id, fold, r, "train", r["class"], cluster_of) for r in train]
    out += [_member(split_id, fold, r, "train_tc", "pos", cluster_of) for r in train_tc]
    out += [_member(split_id, fold, r, "test", r["class"], cluster_of) for r in test]
    out += [_member(split_id, fold, r, "test_tc", "pos", cluster_of) for r in test_tc]
    out += [_member(split_id, fold, r, "test_lit", r["class"], cluster_of) for r in test_lit]
    return out


def _removed(split_id, fold, tc_out):
    return [
        {
            "split_id": split_id,
            "fold": str(fold),
            "seq_sha256": r["seq_sha256"],
            "gene_ids": r["gene_ids"],
            "taxon_ids": r["taxon_ids"],
            "rule": rule,
        }
        for r, rule in tc_out
    ]


def _groups(members):
    out: dict[tuple, list] = {}
    for m in members:
        out.setdefault((m["split_id"], m["fold"]), []).append(m)
    return out


def check_test_truth(members) -> None:
    """T-c rows are never test truth (spec 2.3, Q4)."""
    bad = [m for m in members if m["origin"] == "tc" and m["part"] in ("test", "test_lit")]
    if bad:
        raise SplitError(f"{bad[0]['split_id']}: T-c row {bad[0]['seq_sha256']} is in a test set")


def check_no_shared_hash(members) -> None:
    for (split_id, fold), rows in _groups(members).items():
        train = {m["seq_sha256"] for m in rows if m["part"] in TRAIN_PARTS}
        test = {m["seq_sha256"] for m in rows if m["part"] in TEST_PARTS}
        both = train & test
        if both:
            raise SplitError(f"{split_id} fold {fold}: {sorted(both)[0]} is in train and test")


def check_no_cluster_spans(members) -> None:
    """S1: no cluster in both train and test of a fold. S2, S3: no T-c training row shares a
    cluster with the test set (ruling C-5)."""
    for (split_id, fold), rows in _groups(members).items():
        if split_id == "S1":
            train = {m["cluster_id"] for m in rows if m["part"] in TRAIN_PARTS}
        else:
            train = {m["cluster_id"] for m in rows if m["part"] == "train_tc"}
        test = {m["cluster_id"] for m in rows if m["part"] in TEST_PARTS}
        both = train & test
        if both:
            raise SplitError(
                f"{split_id} fold {fold}: cluster {sorted(both)[0]} spans train and test"
            )


def parse_hits(lines) -> dict[str, float]:
    """Highest fident per query from `query TAB target TAB fident` lines; self hits skipped."""
    best: dict[str, float] = {}
    for n, line in enumerate(lines, 1):
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 3:
            raise SplitError(f"search result line {n}: expected 3 fields, got {len(parts)}")
        query, target, fident = parts
        if query == target:
            continue
        value = float(fident)
        if value > best.get(query, -1.0):
            best[query] = value
    return best


def identity_rows(split_id: str, test_hashes, best: dict) -> list[dict]:
    rows = []
    for h in sorted(set(test_hashes)):
        value = best.get(h)
        rows.append(
            {
                "split_id": split_id,
                "seq_sha256": h,
                "max_identity": "" if value is None else f"{value:.4f}",
                "below_0.3": "yes" if value is None or value < IDENTITY_CUT else "no",
            }
        )
    return rows
```

- [ ] **Step 6: Write `analysis/step1_compare/phasec/09_make_splits.py`**

```python
#!/usr/bin/env python3
"""Phase C step 2: MMseqs2 clusters, split tables, T-c removals, maximum identity.

Reads $STEP1_WORKDIR/phasec/eval_table.tsv.gz, eval_literature.tsv and eval_sequences.fasta.gz
(08; their SHA-256 must equal build_run.json) and species.tsv (its SHA-256 must equal the one
08 used). Runs MMseqs2 (spec 3.4 items 3 and 6):

  easy-cluster  --min-seq-id 0.3 -c 0.5 --cov-mode 0           over all table and literature
                                                               sequences
  easy-search   -s 7.5 -c 0.5 --cov-mode 0                     per S2 and S3 split: query =
                --format-output query,target,fident            test proteins, target = GO
                                                               training proteins (same for
                                                               V-go and V-kw)

Writes to $STEP1_WORKDIR/phasec/:

  clusters.tsv.gz       seq_sha256 -> cluster_id (the representative's hash)
  split_members.tsv.gz  one row per split, fold and sequence: part, origin, class, cluster
  tc_removed.tsv        one row per T-c row removed from a split (rule a, b or c)
  max_identity.tsv.gz   one row per S2/S3 split and test sequence; no hit = below 0.3
  splits_run.json       MMseqs2 version and commands, seed, counts, input hashes

STOP (exit 2, no output): stale 08 outputs; another species.tsv than 08 used; `mmseqs version`
fails (exit 132 means an AVX2 build on a CPU without AVX2: run on partition epyc or set
STEP1_MMSEQS to a non-AVX2 binary); an MMseqs2 command fails; a sequence without a cluster; a
split that breaks a leakage rule (splits.check_*).
"""

import argparse
import gzip
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import evalio
import manifest
import paths
import runinfo
import splits
import truth_table

OUTPUT_NAMES = (
    "clusters.tsv.gz",
    "split_members.tsv.gz",
    "tc_removed.tsv",
    "max_identity.tsv.gz",
    "splits_run.json",
)
CLUSTER_ARGS = ("--min-seq-id", "0.3", "-c", "0.5", "--cov-mode", "0")
SEARCH_ARGS = (
    "-s",
    "7.5",
    "-c",
    "0.5",
    "--cov-mode",
    "0",
    "--format-output",
    "query,target,fident",
)


def mmseqs_version(mmseqs: str) -> str:
    try:
        done = subprocess.run([mmseqs, "version"], capture_output=True, text=True)
    except OSError as exc:
        raise evalio.StopError(
            f"cannot run {mmseqs} ({exc}); module load MMseqs2/17-b804f"
        ) from exc
    if done.returncode != 0:
        raise evalio.StopError(
            f"`{mmseqs} version` exited with {done.returncode} (132 = illegal instruction: an "
            "AVX2 build on a CPU without AVX2; run on partition epyc or set STEP1_MMSEQS to "
            "the non-AVX2 binary)"
        )
    return done.stdout.strip()


def run_mmseqs(argv: list[str], log_path: Path) -> None:
    with open(log_path, "a") as log:
        done = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
    if done.returncode != 0:
        raise evalio.StopError(f"mmseqs {argv[1]} exited with {done.returncode}; see {log_path}")


def write_fasta(path: Path, hashes, seqs: dict) -> None:
    with open(path, "w") as handle:
        for h in sorted(set(hashes)):
            handle.write(f">{h}\n{seqs[h]}\n")


def read_fasta_gz(path: Path) -> dict[str, str]:
    seqs, head = {}, None
    for line in gzip.decompress(Path(path).read_bytes()).decode().splitlines():
        if line.startswith(">"):
            head = line[1:].strip()
            seqs[head] = ""
        elif head is not None:
            seqs[head] += line.strip()
    return seqs


def run(work: Path, species_path: Path, mmseqs: str, tmp: Path, threads: int, arguments=()):
    out = evalio.out_dir(work)
    build = evalio.require_current(
        out,
        "build_run.json",
        ("eval_table.tsv.gz", "eval_literature.tsv", "eval_sequences.fasta.gz"),
        "08_build_eval_tables.py",
    )
    if manifest.sha256_file(species_path) != build["input_sha256"]["species.tsv"]:
        raise evalio.StopError(f"{species_path} differs from the species.tsv that 08 used")
    version = mmseqs_version(mmseqs)
    tmp = Path(tmp)
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    log_path = tmp / "mmseqs.log"
    seqs = read_fasta_gz(out / "eval_sequences.fasta.gz")
    fasta = tmp / "all.fasta"
    write_fasta(fasta, seqs, seqs)
    commands = []
    argv = [mmseqs, "easy-cluster", str(fasta), str(tmp / "clu"), str(tmp / "w_clu"), *CLUSTER_ARGS]
    argv += ["--threads", str(threads)]
    commands.append(["easy-cluster", *CLUSTER_ARGS])
    run_mmseqs(argv, log_path)
    raw = tmp / "clu_cluster.tsv"
    cluster_of = splits.read_cluster_tsv(raw.read_text().splitlines())
    splits.check_clusters(cluster_of, seqs)
    table = truth_table.read_tsv(out / "eval_table.tsv.gz")
    lit = truth_table.read_tsv(out / "eval_literature.tsv")
    defs = splits.split_definitions(truth_table.read_tsv(species_path))
    members, removed = splits.build(table, lit, cluster_of, defs, evalio.SEED)
    splits.check_test_truth(members)
    splits.check_no_shared_hash(members)
    splits.check_no_cluster_spans(members)
    identity = []
    for d in defs:
        if not d.test_sources:
            continue
        rows = [m for m in members if m["split_id"] == d.split_id]
        query = [m["seq_sha256"] for m in rows if m["part"] in ("test", "test_lit")]
        target = [m["seq_sha256"] for m in rows if m["part"] == "train"]
        q, t, hits = (
            tmp / f"{d.split_id}.q.fasta",
            tmp / f"{d.split_id}.t.fasta",
            tmp / f"{d.split_id}.m8",
        )
        write_fasta(q, query, seqs)
        write_fasta(t, target, seqs)
        argv = [mmseqs, "easy-search", str(q), str(t), str(hits), str(tmp / f"w_{d.split_id}")]
        argv += [*SEARCH_ARGS, "--threads", str(threads)]
        run_mmseqs(argv, log_path)
        best = splits.parse_hits(hits.read_text().splitlines())
        identity += splits.identity_rows(d.split_id, query, best)
    commands.append(["easy-search", *SEARCH_ARGS])
    clusters = [{"seq_sha256": h, "cluster_id": cluster_of[h]} for h in sorted(cluster_of)]
    log = {
        "all_sources": build["all_sources"],
        "truth_set_sha256": build["truth_set_sha256"],
        "input_sha256": {
            n: build["outputs_sha256"][n]
            for n in ("eval_table.tsv.gz", "eval_literature.tsv", "eval_sequences.fasta.gz")
        }
        | {"species.tsv": manifest.sha256_file(species_path)},
        "mmseqs_version": version,
        "mmseqs_commands": commands,
        "cluster_tsv_sha256": manifest.sha256_file(raw),
        "seed": evalio.SEED,
        "sequences": len(seqs),
        "clusters": len(set(cluster_of.values())),
        "splits": [d.split_id for d in defs],
        "members_by_split_fold_part": {
            f"{s}|{f}|{p}": n
            for (s, f, p), n in sorted(
                Counter((m["split_id"], m["fold"], m["part"]) for m in members).items()
            )
        },
        "tc_removed_by_split_fold_rule": {
            f"{s}|{f}|{r}": n
            for (s, f, r), n in sorted(
                Counter((x["split_id"], x["fold"], x["rule"]) for x in removed).items()
            )
        },
        "identity_below_0.3_by_split": dict(
            sorted(Counter(r["split_id"] for r in identity if r["below_0.3"] == "yes").items())
        ),
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    evalio.write_outputs(
        out,
        {
            "clusters.tsv.gz": lambda p: truth_table.write_tsv(p, splits.CLUSTER_COLUMNS, clusters),
            "split_members.tsv.gz": lambda p: truth_table.write_tsv(
                p, splits.MEMBER_COLUMNS, members
            ),
            "tc_removed.tsv": lambda p: truth_table.write_tsv(p, splits.REMOVED_COLUMNS, removed),
            "max_identity.tsv.gz": lambda p: truth_table.write_tsv(
                p, splits.IDENTITY_COLUMNS, identity
            ),
        },
        "splits_run.json",
        log,
    )
    return log


def main(argv=None) -> int:
    import os

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--mmseqs", default=os.environ.get("STEP1_MMSEQS", "mmseqs"))
    parser.add_argument("--tmp-dir", required=True, help="node-local work dir, e.g. $SCRATCH/x")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        log = run(
            work,
            Path(args.species),
            args.mmseqs,
            Path(args.tmp_dir),
            args.threads,
            list(argv) if argv is not None else sys.argv[1:],
        )
    except (evalio.StopError, splits.SplitError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(f"sequences={log['sequences']} clusters={log['clusters']} mmseqs={log['mmseqs_version']}")
    for key, n in log["tc_removed_by_split_fold_rule"].items():
        print(f"  tc_removed {key}={n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_splits.py -q`
Expected: `17 passed, 1 skipped` (the real-MMseqs2 test skips without a working `mmseqs`).
Run: `STEP1_MMSEQS=/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_splits.py -q`
Expected: `18 passed` (prototype: 46 s on c01; each MMseqs2 call costs a few seconds).

- [ ] **Step 8: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/splits.py analysis/step1_compare/phasec/09_make_splits.py \
  tests/step1_compare/fixtures/phasec/stub_mmseqs/mmseqs tests/step1_compare/phasec_fixture.py \
  tests/step1_compare/test_phasec_splits.py
git commit -m "phasec: MMseqs2 clusters, S1-S3 splits, T-c removals a-c, maximum identity (09)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 8: Fit and score every unit (10)

**Files:**
- Create: `analysis/step1_compare/phasec/10_fit_and_score.py`
- Modify: `tests/step1_compare/conftest.py` (add the session fixture `phasec_chain`)
- Test: `tests/step1_compare/test_phasec_fit.py`

**Interfaces:**
- Consumes: `evalio.*` (Task 1); `universe.load`, `universe.MODELS`, `models.run_task`, `models.init_worker`, `models.CANDIDATES`, `models.C_GRID`, `models.INNER_FOLDS`, `models.ModelError` (Task 6); `split_members.tsv.gz`, `clusters.tsv.gz`, `splits_run.json` (Task 7); `build_run.json` (Task 2).
- Produces: `10_fit_and_score.py` with `OUTPUT_NAMES = ("scores.tsv.gz", "scores_run.json")`, `SCORE_COLUMNS`, `SCORED_PARTS`, `unit_order(members, splits_order)`, `tasks(members, splits_order, all_hashes, candidates) -> list[dict]`, `run_tasks(task_list, u, workers, phaseb)`, `write_scores(path, task_list, results)`, `run(work, workers, candidates=models.CANDIDATES, arguments=())`, `main(argv=None)`; options `--work-dir`, `--workers`, `--candidates` (comma list; tests and the golden file use subsets). conftest fixture `phasec_chain` (session scope): the fixture chain 08, 09 (stub), 10 and, once `phasec/11_evaluate.py` exists, 11.

Notes:
- Worker processes use the `spawn` start method (Python 3.14 defaults to `forkserver` on Linux; `spawn` is explicit and does not copy the parent's BLAS state). Each worker loads the universe once (6.6 s on the real data, without the SHA-256 checks that the parent already made). `test_workers_give_the_same_scores` pins byte-identical output for 1 and 2 workers.
- `score` is written with 10 significant digits (`format(x, ".10g")`).

- [ ] **Step 1: Add the session fixture to `tests/step1_compare/conftest.py`**

Append:

```python
@pytest.fixture(scope="session")
def phasec_chain(tmp_path_factory):
    """Phase C fixture chain 08 -> 09 (stub MMseqs2) -> 10 -> 11, built once per session.

    Step 11 runs only when phasec/11_evaluate.py exists (it is written after 10). Tests that
    change a file must copy chain["work"] first (phasec_fixture.copy_work)."""
    pytest.importorskip("sklearn")
    import phasec_fixture

    upto = 11 if (STEP1_DIR / "phasec" / "11_evaluate.py").exists() else 10
    return phasec_fixture.run_chain(tmp_path_factory.mktemp("phasec_chain"), upto, load_phasec)
```

- [ ] **Step 2: Write the failing test**

`tests/step1_compare/test_phasec_fit.py`:

```python
"""10_fit_and_score.py on the fixture chain (08 -> 09 with stub MMseqs2 -> 10)."""

import json

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import phasec_fixture as pf  # noqa: E402
import splits  # noqa: E402
import truth_table  # noqa: E402
from conftest import load_phasec  # noqa: E402

SUBSET = "R2,M8,H"  # re-runs fit only these candidates (time); the rest is fitted in the chain


def _run10(work, *extra):
    return load_phasec("10_fit_and_score").main(["--work-dir", str(work), *extra])


def _units(work):
    return json.loads((work / "phasec" / "scores_run.json").read_text())["units"]


def _score_rows(work):
    rows = truth_table.read_tsv(work / "phasec" / "scores.tsv.gz")
    return [r for r in rows if r["candidate"] in SUBSET.split(",")]


def _subset(units):
    return {k: {c: v[c] for c in SUBSET.split(",")} for k, v in units.items()}


def test_scores_cover_every_unit_candidate_and_scored_row(phasec_chain):
    out = phasec_chain["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    units = {(m["split_id"], m["fold"]) for m in members}
    assert set(_units(phasec_chain["work"])) == {
        f"{s}|{f}|{v}" for s, f in units for v in ("V-go", "V-kw")
    }
    rows = truth_table.read_tsv(out / "scores.tsv.gz")
    n_unique = len(truth_table.read_tsv(phasec_chain["work"] / "phaseb" / "features_unique.tsv.gz"))
    full = [
        r
        for r in rows
        if r["split_id"] == "FULL" and r["variant"] == "V-go" and r["candidate"] == "M8"
    ]
    assert len(full) == n_unique and {r["part"] for r in full} == {"all"}
    s1 = {
        (r["seq_sha256"], r["part"])
        for r in rows
        if r["split_id"] == "S1" and r["candidate"] == "R2"
    }
    want = {(m["seq_sha256"], m["part"]) for m in members
            if m["split_id"] == "S1" and m["part"] in ("test", "test_tc")}  # fmt: skip
    assert s1 == want
    assert all(r["score"] == "" and r["prob"] == "" and r["call"] in ("0", "1")
               for r in rows if r["candidate"] == "R2")  # fmt: skip
    assert all(r["score"] != "" and 0.0 <= float(r["prob"]) <= 1.0
               for r in rows if r["candidate"] == "H")  # fmt: skip


def test_fitted_settings_are_recorded(phasec_chain):
    units = _units(phasec_chain["work"])
    unit = units["S1|0|V-go"]
    assert unit["M8"]["C"] in (0.001, 0.003, 0.01, 0.1, 1.0, 10.0)  # ruling C-12
    assert unit["H"]["h_variant"] in ("M8", "M35", "M8-C", "M35-C")
    assert unit["R2"]["g"] in ("highly_probable", "probable", "weakly")
    assert unit["R2"]["t"] in (0.20, 0.25, 0.30, 0.35, 0.40)  # owner decision 2026-10-01
    assert set(unit["M35"]["inner_pr_auc"]) == {"0.001", "0.003", "0.01", "0.1", "1.0", "10.0"}
    assert len(unit["H"]["inner_pr_auc"]) == 4 * 6  # (ESM variant, C) pairs
    assert units["S1|0|V-kw"]["M8"]["n_train_pos"] > unit["M8"]["n_train_pos"]  # T-c positives


def test_training_rows_keep_homology_only_rows(phasec_chain):
    # spec 3.2: the direct-evidence filter applies to test rows only (11 applies it)
    out = phasec_chain["work"] / "phasec"
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    train = [m for m in members if m["split_id"] == "FULL" and m["part"] == "train"]
    assert any(table[m["seq_sha256"]]["homology_only"] == "yes" for m in train)
    assert _units(phasec_chain["work"])["FULL|0|V-go"]["B0"]["n_train"] == len(train)


def test_threshold_fit_uses_train_only(phasec_chain, tmp_path):
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    test = [
        m for m in members if m["part"] in ("test", "test_lit") and m["class"] in ("pos", "neg")
    ]
    before = [m["class"] for m in test]
    for m, c in zip(test, np.random.default_rng(0).permutation(before), strict=True):
        m["class"] = str(c)
    assert [m["class"] for m in test] != before  # the outer test labels really changed
    truth_table.write_tsv(out / "split_members.tsv.gz", splits.MEMBER_COLUMNS, members)
    pf.fix_recorded_hash(
        out / "splits_run.json", "split_members.tsv.gz", out / "split_members.tsv.gz"
    )
    assert _run10(fx["work"], "--candidates", SUBSET) == 0
    # t, g, C, the ML threshold, Platt a and b and the H variant: all unchanged
    assert _units(fx["work"]) == _subset(_units(phasec_chain["work"]))
    # the scaler is not stored; identical scores show that it did not change either
    assert _score_rows(fx["work"]) == _score_rows(phasec_chain["work"])


def test_workers_give_the_same_scores(phasec_chain, tmp_path):
    one = pf.copy_work(phasec_chain, tmp_path / "one")
    two = pf.copy_work(phasec_chain, tmp_path / "two")
    assert _run10(one["work"], "--candidates", SUBSET) == 0
    assert _run10(two["work"], "--candidates", SUBSET, "--workers", "2") == 0
    a = (one["work"] / "phasec" / "scores.tsv.gz").read_bytes()
    assert (two["work"] / "phasec" / "scores.tsv.gz").read_bytes() == a
    assert _units(one["work"]) == _subset(_units(phasec_chain["work"]))


def test_unknown_candidate_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    assert _run10(fx["work"], "--candidates", "M8,XGB") == 2
    assert "unknown candidates ['XGB']" in capsys.readouterr().err


def test_stale_input_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    path = fx["work"] / "phasec" / "clusters.tsv.gz"
    path.write_bytes(path.read_bytes() + b"\n")
    before = (fx["work"] / "phasec" / "scores.tsv.gz").read_bytes()
    assert _run10(fx["work"]) == 2
    assert "clusters.tsv.gz differs from the SHA-256 in splits_run.json" in capsys.readouterr().err
    assert (fx["work"] / "phasec" / "scores.tsv.gz").read_bytes() == before  # not replaced


def test_changed_embedding_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    emb = fx["work"] / "phaseb" / "emb" / "esm2_t6_8M_UR50D.nterm.npy"
    arr = np.load(emb)
    arr[0, 0] += 1.0
    np.save(emb, arr)
    assert _run10(fx["work"]) == 2
    assert "differs from embedding_run.json" in capsys.readouterr().err


def test_too_few_positive_clusters_stops(phasec_chain, tmp_path, capsys):
    # Review Focus 5: an S2 training set whose positives fall into 2 clusters
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    for m in members:
        if m["split_id"] == "S2-Calb_CGD" and m["part"] == "train" and m["class"] == "pos":
            m["cluster_id"] = "one_cluster"
    truth_table.write_tsv(out / "split_members.tsv.gz", splits.MEMBER_COLUMNS, members)
    pf.fix_recorded_hash(
        out / "splits_run.json", "split_members.tsv.gz", out / "split_members.tsv.gz"
    )
    assert _run10(fx["work"], "--candidates", "R2,M8") == 2
    err = capsys.readouterr().err
    assert "S2-Calb_CGD|0|V-go: 1 positive and" in err and "clusters" in err
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_fit.py -q`
Expected: `9 errors`: the session fixture `phasec_chain` stops with `FileNotFoundError: ... phasec/10_fit_and_score.py`.

- [ ] **Step 4: Write `analysis/step1_compare/phasec/10_fit_and_score.py`**

```python
#!/usr/bin/env python3
"""Phase C step 3: fit every candidate on every outer training set and score (spec 3.3).

Reads $STEP1_WORKDIR/phasec/split_members.tsv.gz and clusters.tsv.gz (09; SHA-256 checked
against splits_run.json), and the Phase B universe: phaseb/features_unique.tsv.gz and
phaseb/unique_sequences.tsv.gz (SHA-256 checked against build_run.json) and the four
embedding matrices (SHA-256 checked against embedding_run.json).

One work unit = one (split, fold, variant). Training rows: part `train` (V-go) or `train` and
`train_tc` (V-kw); y = class pos. Scored rows: parts test, test_tc, test_lit; for FULL every
Phase B unique sequence. The nested protocol is in models.py. Units run in --workers processes.

Writes to $STEP1_WORKDIR/phasec/:

  scores.tsv.gz    one row per unit, candidate and scored sequence: score (decision value;
                   empty for rules), prob (Platt), call (1/0)
  scores_run.json  fitted settings per unit and candidate (C, g, t, threshold, Platt a and b,
                   H variant, inner PR-AUC), seed, input hashes, library versions

STOP (exit 2, no output): stale 08 or 09 outputs; changed Phase B files; an outer training set
that cannot be fitted (for example fewer than 3 positive clusters).
"""

import argparse
import gzip
import io
import multiprocessing
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import evalio
import manifest
import models
import paths
import runinfo
import truth_table
import universe

OUTPUT_NAMES = ("scores.tsv.gz", "scores_run.json")
SCORE_COLUMNS = (
    "split_id",
    "fold",
    "variant",
    "candidate",
    "seq_sha256",
    "part",
    "score",
    "prob",
    "call",
)
SCORED_PARTS = ("test", "test_tc", "test_lit")


def unit_order(members, splits_order) -> list[tuple[str, str]]:
    keys = {(m["split_id"], m["fold"]) for m in members}
    rank = {s: i for i, s in enumerate(splits_order)}
    return sorted(keys, key=lambda k: (rank[k[0]], int(k[1])))


def tasks(members, splits_order, all_hashes, candidates) -> list[dict]:
    """One task per (split, fold, variant). Only training rows carry labels into a task."""
    by_unit: dict[tuple, list] = {}
    for m in members:
        by_unit.setdefault((m["split_id"], m["fold"]), []).append(m)
    out = []
    for key in unit_order(members, splits_order):
        rows = by_unit[key]
        if key[0] == "FULL":
            score, parts = list(all_hashes), {}
        else:
            scored = [m for m in rows if m["part"] in SCORED_PARTS]
            score = [m["seq_sha256"] for m in scored]
            parts = {m["seq_sha256"]: m["part"] for m in scored}
        for variant in evalio.VARIANTS:
            use = ("train", "train_tc") if variant == "V-kw" else ("train",)
            train = [m for m in rows if m["part"] in use]
            bad = [m for m in train if m["class"] not in ("pos", "neg")]
            if bad:
                raise evalio.StopError(
                    f"{key}: training row {bad[0]['seq_sha256']} has class {bad[0]['class']}"
                )
            out.append(
                {
                    "key": f"{key[0]}|{key[1]}|{variant}",
                    "split_id": key[0],
                    "fold": key[1],
                    "variant": variant,
                    "train": [m["seq_sha256"] for m in train],
                    "y": [m["class"] == "pos" for m in train],
                    "groups": [m["cluster_id"] for m in train],
                    "score": score,
                    "parts": parts,
                    "seed": evalio.SEED,
                    "candidates": candidates,
                }
            )
    return out


def run_tasks(task_list, u, workers: int, phaseb: Path):
    if workers <= 1:
        return [models.run_task(t, u) for t in task_list]
    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=ctx, initializer=models.init_worker, initargs=(str(phaseb),)
    ) as pool:
        return list(pool.map(models.run_task, task_list))


def _fmt(x) -> str:
    return format(float(x), ".10g")


def write_scores(path: Path, task_list, results) -> None:
    raw = open(path, "wb")
    with raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
        text.write("\t".join(SCORE_COLUMNS) + "\n")
        for task, (key, _, scores) in zip(task_list, results, strict=True):
            assert key == task["key"]
            head = f"{task['split_id']}\t{task['fold']}\t{task['variant']}\t"
            for cand in task["candidates"]:
                s = scores[cand]
                for i, h in enumerate(task["score"]):
                    part = task["parts"].get(h, "all")
                    score = "" if s["score"] is None else _fmt(s["score"][i])
                    prob = "" if s["prob"] is None else _fmt(s["prob"][i])
                    call = "1" if s["call"][i] else "0"
                    text.write(f"{head}{cand}\t{h}\t{part}\t{score}\t{prob}\t{call}\n")
        text.flush()
        text.detach()


def run(work: Path, workers: int, candidates=models.CANDIDATES, arguments=()):
    work = Path(work)
    out = evalio.out_dir(work)
    build = evalio.require_current(out, "build_run.json", ("eval_table.tsv.gz",), "08")
    split_log = evalio.require_current(
        out, "splits_run.json", ("clusters.tsv.gz", "split_members.tsv.gz"), "09_make_splits.py"
    )
    phaseb = work / "phaseb"
    for name in ("features_unique.tsv.gz", "unique_sequences.tsv.gz"):
        if manifest.sha256_file(phaseb / name) != build["input_sha256"][name]:
            raise evalio.StopError(f"phaseb/{name} changed after 08 read it; re-run 08 to 10")
    emb_run = evalio.read_json(phaseb / "emb" / "embedding_run.json")
    if emb_run.get("unique_sequences_sha256") != build["input_sha256"]["unique_sequences.tsv.gz"]:
        raise evalio.StopError(
            "embedding_run.json names another unique_sequences.tsv.gz than 08 read"
        )
    u = universe.load(phaseb, verify=True)
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    task_list = tasks(members, split_log["splits"], u.hashes, tuple(candidates))
    results = run_tasks(task_list, u, workers, phaseb)
    log = {
        "all_sources": build["all_sources"],
        "truth_set_sha256": build["truth_set_sha256"],
        "input_sha256": {
            "split_members.tsv.gz": split_log["outputs_sha256"]["split_members.tsv.gz"],
            "clusters.tsv.gz": split_log["outputs_sha256"]["clusters.tsv.gz"],
            "features_unique.tsv.gz": build["input_sha256"]["features_unique.tsv.gz"],
            "unique_sequences.tsv.gz": build["input_sha256"]["unique_sequences.tsv.gz"],
        },
        "embedding_array_sha256": {
            f"{m}.{w}": emb_run["models"][m][w]["array_sha256"]
            for m in universe.MODELS
            for w in ("nterm", "cterm")
        },
        "seed": evalio.SEED,
        "candidates": list(candidates),
        "c_grid": list(models.C_GRID),
        "inner_folds": models.INNER_FOLDS,
        "units": {key: params for key, params, _ in results},
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    evalio.write_outputs(
        out,
        {"scores.tsv.gz": lambda p: write_scores(p, task_list, results)},
        "scores_run.json",
        log,
    )
    return log


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--candidates", default=",".join(models.CANDIDATES), help="comma list (tests, golden file)"
    )
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        cands = tuple(args.candidates.split(","))
        unknown = sorted(set(cands) - set(models.CANDIDATES))
        if unknown:
            raise evalio.StopError(f"unknown candidates {unknown}")
        cands = tuple(c for c in models.CANDIDATES if c in cands)
        log = run(work, args.workers, cands, list(argv) if argv is not None else sys.argv[1:])
    except (evalio.StopError, models.ModelError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(f"units={len(log['units'])} candidates={len(log['candidates'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_fit.py -q`
Expected: `9 passed`.

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/10_fit_and_score.py tests/step1_compare/conftest.py \
  tests/step1_compare/test_phasec_fit.py
git commit -m "phasec: fit and score all candidates, both variants, all splits (10)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 9: Evaluation, findings and proteome calls (11)

**Files:**
- Create: `analysis/step1_compare/phasec/findings.py`, `analysis/step1_compare/phasec/agreement.py`, `analysis/step1_compare/phasec/11_evaluate.py`, `tests/step1_compare/fixtures/phasec/golden_metrics.json.gz` (generated in Step 6)
- Test: `tests/step1_compare/test_phasec_evaluate.py`

**Interfaces:**
- Consumes: `evalio.*` (Task 1); `metrics.*` (Task 3); `bootstrap.*` (Task 4); `splits.check_test_truth`, `splits.IDENTITY_CUT` (Task 7); `scores.tsv.gz`, `scores_run.json` (Task 8); `phaseb/sequence_members.tsv.gz`, `phaseb/features_run.json`, `sequence_sets.tsv`.
- Produces: `findings.ESTIMATE_HALF_WIDTH = 0.10`, `SATURATION_AUC = 0.99`, `LABEL_CANDIDATES = ("R2", "M8", "M35", "M8-C", "M35-C", "H")`, `ML`, `COMPARATORS = ("B1", "R2")`, `MIN_DIRECT_POSITIVES = 20`, `floor_met(n_direct_positives) -> bool`, `estimate_label(recall_half_widths, n_direct_positives) -> str`, `finding_a(b1_auc) -> dict`, `beats(diff) -> bool`, `finding_b(diffs) -> dict`, `finding_c(per_test_set) -> dict`. `agreement.PANEL` (22 entries of parent spec 6), `score_source(h, s1_hashes, full_train) -> str`, `find_panel_hashes(panel, members)`, `agreement_counts(rule_calls, ml_calls) -> dict`. `11_evaluate.py`: `OUTPUT_NAMES`, `METRICS_SCHEMA = "step1-phasec-metrics/1"`, `FINDINGS_SCHEMA = "step1-phasec-findings/1"`, `PRECISION_METRICS`, `STRATA`, `TRUTHS`, `PREVALENCES`, `PROTEOME_BASE`, `proteome_columns(cands) -> tuple`, `test_sets(...)`, `evaluate_truth(...)`, `build_findings(test_blocks) -> dict`, `run(work, sets_path, species_path, n_resamples, arguments=())`, `main(argv=None)`; options `--work-dir`, `--sets`, `--species`, `--n-resamples` (default 2,000).

Decisions in this task (state them in review):
- C for every LR candidate is chosen by inner out-of-fold PR-AUC (ties: the smaller C), the criterion the spec gives for H's variant; the spec now states it (spec 3.3, ruling C-12).
- Findings (b) and (c) are evaluated per ML candidate; `holds_for` lists the candidates for which the statement holds; (c) holds for a candidate only when it holds on all three S2 test sets. No headline ML candidate is named (ruling C-13); the report says "any".
- The literature test set counts as direct truth (it has no homology flag). It has positives only: FPR and AUC are `null` because they are not defined, and `evaluate_truth` sets precision, PR-AUC, precision at recall 0.8 and 0.9 (`PRECISION_METRICS`) and the precision at the rule's recall to `null` for a truth block without negatives (review I-1). Without this step they would be 1.0 by construction. The metric functions themselves are unchanged: `test_undefined_metrics_are_nan` (Task 3) still asserts that `pr_auc` of a positives-only input returns 1.0, as a property of the function; that value is no longer reported.
- The estimate label (ruling C-8, owner decision 2026-10-01): every label candidate has a recall half-width of at most 0.10 and the test set has at least 20 direct-evidence positives (`floor_met`). `n_direct_positives`, `floor_met` and `zero_width_recall_interval` (the label candidates whose recall interval has lo equal to hi) are stored beside the label.
- Calibration bins are 10 equal-width bins of the Platt probability on the S2 test sets (direct truth). The Brier score has no interval in `calibration` (the spec asks for the score and the reliability data).

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_evaluate.py`:

```python
"""11_evaluate.py: metrics.json, labels, prevalence, proteome calls, findings, golden file."""

import gzip
import json
import math
import os

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import findings  # noqa: E402
import metrics as mt  # noqa: E402
import phasec_fixture as pf  # noqa: E402
import splits  # noqa: E402
import truth_table  # noqa: E402
from conftest import TESTS_DIR, load_phasec  # noqa: E402

GOLDEN = TESTS_DIR / "fixtures" / "phasec" / "golden_metrics.json.gz"


def _metrics(work):
    return json.loads((work / "phasec" / "metrics.json").read_text())


def _run11(fx, n=50):
    return load_phasec("11_evaluate").main(pf.eval_argv(fx, n))


def test_metrics_json_has_the_test_sets_and_no_nan(phasec_chain):
    text = (phasec_chain["work"] / "phasec" / "metrics.json").read_text()

    def refuse(name):
        raise AssertionError(f"metrics.json holds {name}")

    m = json.loads(text, parse_constant=refuse)
    assert set(m["test_sets"]) == {
        "S1:all", "S1:Calb_CGD", "S1:Scer_SGD", "S2-Calb_CGD:Calb_CGD", "S2-Scer_SGD:Scer_SGD",
        "S2-Spom_PomBase:Spom_PomBase", "S3-Basidiomycota:clade",
        "S3-Basidiomycota:Cneo_H99_GOA", "S3-Basidiomycota:Umay_MYCMD",
        "S3-Eurotiomycetes:clade", "S3-Eurotiomycetes:Afum_ASPFU",
        "S3-Eurotiomycetes:literature",
    }  # fmt: skip
    s1 = m["test_sets"]["S1:all"]["truth"]["direct"]["n"]
    assert "identity_below_0.3" not in s1  # S1 has no maximum-identity stratum
    s2 = m["test_sets"]["S2-Spom_PomBase:Spom_PomBase"]["truth"]["direct"]["n"]
    assert "identity_below_0.3" in s2
    assert all(ts["label"] in ("estimate", "smoke test") for ts in m["test_sets"].values())


def test_undefined_metrics_are_null(phasec_chain):
    # Review Focus 1: a stratum without positives, a test set without negatives
    m = _metrics(phasec_chain["work"])
    nsec = m["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]["N-sec"]["V-go"]["M8"]
    assert nsec["recall"]["value"] is None and nsec["fpr"]["value"] is not None
    lit = m["test_sets"]["S3-Eurotiomycetes:literature"]["truth"]["direct"]["metrics"]["all"][
        "V-go"
    ]
    assert lit["M8"]["roc_auc"]["value"] is None and lit["M8"]["fpr"]["value"] is None
    assert lit["M8"]["recall"]["value"] is not None
    # review I-1: without negatives precision is 1 by construction; spec 3.2 says recall only
    for c in ("M8", "R2", "B1"):
        assert lit[c]["precision"]["value"] is None, c
    assert lit["M8"]["pr_auc"]["value"] is None
    assert lit["M8"]["precision_at_recall_0.8"]["value"] is None
    assert lit["M8"]["precision_at_recall_0.9"]["value"] is None
    lit_truth = m["test_sets"]["S3-Eurotiomycetes:literature"]["truth"]["direct"]
    vs_rule = lit_truth["vs_rule"]["V-go"]["M8"]["precision_at_rule_recall"]
    assert vs_rule["value"] is None
    # a test set with negatives keeps its precision
    s1 = m["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]["all"]["V-go"]["M8"]
    assert s1["precision"]["value"] is not None and s1["pr_auc"]["value"] is not None


def test_direct_filter_applies_to_test_rows_only(phasec_chain):
    out = phasec_chain["work"] / "phasec"
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    test = [
        m
        for m in members
        if m["split_id"] == "S1" and m["part"] == "test" and m["class"] in ("pos", "neg")
    ]
    direct = [m for m in test if table[m["seq_sha256"]]["homology_only"] == "no"]
    assert len(direct) < len(test)
    n = _metrics(phasec_chain["work"])["test_sets"]["S1:all"]["truth"]
    assert n["all"]["n"]["all"]["pos"] + n["all"]["n"]["all"]["neg"] == len(test)
    assert n["direct"]["n"]["all"]["pos"] + n["direct"]["n"]["all"]["neg"] == len(direct)
    train = [
        m for m in members if m["split_id"] == "S1" and m["fold"] == "0" and m["part"] == "train"
    ]
    assert any(table[m["seq_sha256"]]["homology_only"] == "yes" for m in train)


def test_estimate_label_rule():
    ok = {c: 0.10 for c in findings.LABEL_CANDIDATES}
    assert findings.estimate_label(ok, 20) == "estimate"
    assert findings.estimate_label({**ok, "H": 0.1001}, 20) == "smoke test"
    assert findings.estimate_label({**ok, "R2": 0.0999}, 20) == "estimate"
    assert findings.estimate_label({**ok, "M35": None}, 20) == "smoke test"  # no positives
    assert findings.estimate_label({c: 0.05 for c in ("R2", "M8")}, 50) == "smoke test"  # missing
    assert "B1" not in findings.LABEL_CANDIDATES and "R0" not in findings.LABEL_CANDIDATES
    # count floor (owner decision 2026-10-01): a zero-width interval of 7 positives is a smoke test
    zero = {c: 0.0 for c in findings.LABEL_CANDIDATES}
    assert findings.estimate_label(zero, 7) == "smoke test"
    assert findings.estimate_label(zero, 19) == "smoke test"
    small = {c: 0.04 for c in findings.LABEL_CANDIDATES}
    assert findings.estimate_label(small, 20) == "estimate"
    assert findings.floor_met(20) and not findings.floor_met(19)


def test_label_in_metrics_follows_the_half_widths(phasec_chain):
    for name, ts in _metrics(phasec_chain["work"])["test_sets"].items():
        n_pos = ts["truth"]["direct"]["n"]["all"]["pos"]
        assert ts["n_direct_positives"] == n_pos, name
        assert ts["floor_met"] is (n_pos >= findings.MIN_DIRECT_POSITIVES), name
        assert ts["label"] == findings.estimate_label(ts["recall_half_width"], n_pos), name
        direct = ts["truth"]["direct"]["metrics"]["all"]["V-go"]
        r2 = direct["R2"]["recall"]
        if r2["lo"] is not None:
            assert ts["recall_half_width"]["R2"] == pytest.approx((r2["hi"] - r2["lo"]) / 2)


def test_prevalence_table_formula(phasec_chain):
    # hand-computed: recall 0.8, FPR 0.1, prevalence 0.05 -> 0.04 / (0.04 + 0.095)
    assert mt.precision_at_prevalence(0.8, 0.1, 0.05) == pytest.approx(0.2962962962962963)
    block = _metrics(phasec_chain["work"])["test_sets"]["S1:all"]["truth"]["direct"]
    r2 = block["metrics"]["all"]["V-go"]["R2"]
    for key, pi in (("0.01", 0.01), ("0.05", 0.05), ("0.1", 0.10)):
        want = mt.precision_at_prevalence(r2["recall"]["value"], r2["fpr"]["value"], pi)
        assert block["prevalence"]["V-go"]["R2"][key]["value"] == pytest.approx(float(want))


def test_proteome_calls_use_oof_for_training_hashes(phasec_chain):
    out = phasec_chain["work"] / "phasec"
    calls = truth_table.read_tsv(out / "proteome_calls.tsv.gz")
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    s1 = {
        m["seq_sha256"]: m["fold"]
        for m in members
        if m["split_id"] == "S1" and m["part"] in ("test", "test_tc")
    }
    train = {m["seq_sha256"] for m in members if m["part"] in ("train", "train_tc")}
    scores = {(r["split_id"], r["fold"], r["variant"], r["candidate"], r["seq_sha256"]): r
              for r in truth_table.read_tsv(out / "scores.tsv.gz")}  # fmt: skip
    seen_oof = 0
    for r in calls:
        h = r["seq_sha256"]
        if h in train:
            assert r["score_source"] != "final", h  # a training hash never gets `final`
        if h in s1:
            seen_oof += 1
            assert r["score_source"] == "oof"
            want = scores[("S1", s1[h], "V-go", "M8", h)]
            assert r["score_M8_V-go"] == want["score"] and r["call_M8_V-go"] == want["call"]
        else:
            assert r["score_source"] == "final"
            assert r["score_M8_V-go"] == scores[("FULL", "0", "V-go", "M8", h)]["score"]
    assert seen_oof >= 2  # the fixture puts truth and T-c sequences into the proteome sets
    assert {r["set_id"] for r in calls} == {"Scer_proteome", "Cimm_RS_proteome"}


def test_findings_booleans():
    assert findings.finding_a({"S1:all": 0.995, "S2-x:x": 0.95})["holds"] is False  # saturated
    assert findings.finding_a({"S1:all": 0.97, "S2-x:x": 0.95})["holds"] is True
    assert findings.finding_a({"S1:all": None})["holds"] is False
    assert findings.finding_a({})["holds"] is False
    up = {"value": 0.2, "lo": 0.05, "hi": 0.3}
    zero = {"value": 0.1, "lo": -0.01, "hi": 0.3}
    b = findings.finding_b({"M8": {"B1": up, "R2": up}, "M35": {"B1": up, "R2": zero}})
    assert b["holds_for"] == ["M8"] and b["candidates"]["M35"]["beats_R2"] is False
    c = findings.finding_c(
        {"S2-a:a": b, "S2-b:b": findings.finding_b({"M8": {"B1": zero, "R2": up}})}
    )
    assert c["holds_for"] == []
    c = findings.finding_c({"S2-a:a": b, "S2-b:b": b})
    assert c["holds_for"] == ["M8"]


def test_findings_json_reads_metrics_json(phasec_chain):
    m = _metrics(phasec_chain["work"])
    f = json.loads((phasec_chain["work"] / "phasec" / "findings.json").read_text())
    b1 = m["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]["all"]["V-go"]["B1"]["roc_auc"]
    assert f["a_b1_not_saturated"]["b1_roc_auc"]["S1:all"] == b1["value"]
    diff = m["test_sets"]["S1:all"]["truth"]["direct"]["nsec_fpr_at_rule_recall"]["V-go"]["diff"]
    assert f["b_ml_beats_b1_and_r2_on_nsec_s1"]["candidates"]["M8"]["B1"] == diff["M8"]["B1"]
    assert set(f["c_same_under_s2"]["test_sets"]) == {
        "S2-Calb_CGD:Calb_CGD",
        "S2-Scer_SGD:Scer_SGD",
        "S2-Spom_PomBase:Spom_PomBase",
    }


def test_named_panel_and_context(phasec_chain):
    m = _metrics(phasec_chain["work"])
    names = [p["name"] for p in m["named_panel"]]
    assert names[:2] == ["FLO1", "SAG1"] and "LIT2" in names  # hard_negative row appended
    lit2 = next(p for p in m["named_panel"] if p["name"] == "LIT2")
    assert lit2["found"] and lit2["literature_class"] == "hard_negative"
    assert not next(p for p in m["named_panel"] if p["name"] == "FLO1")["found"]
    ctx = m["context"]["onygenales_tc_rows_in_vkw_training"]
    assert ctx.get("S1|0", 0) + ctx.get("S1|1", 0) + ctx.get("S1|2", 0) >= 1
    assert "S3-Eurotiomycetes|0" not in ctx  # removed by rule c (ruling C-7)


def test_stale_input_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    path = fx["work"] / "phasec" / "scores.tsv.gz"
    path.write_bytes(path.read_bytes() + b"\n")
    assert _run11(fx) == 2
    assert "scores.tsv.gz differs from the SHA-256 in scores_run.json" in capsys.readouterr().err


def test_tc_row_in_a_test_set_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    next(m for m in members if m["origin"] == "tc" and m["split_id"] == "S2-Calb_CGD")["part"] = (
        "test"
    )
    truth_table.write_tsv(out / "split_members.tsv.gz", splits.MEMBER_COLUMNS, members)
    pf.fix_recorded_hash(
        out / "splits_run.json", "split_members.tsv.gz", out / "split_members.tsv.gz"
    )
    pf.fix_recorded_hash(out / "scores_run.json", "split_members.tsv.gz",
                         out / "split_members.tsv.gz", key="input_sha256")  # fmt: skip
    assert _run11(fx) == 2
    assert "is in a test set" in capsys.readouterr().err


def _close(a, b, path="$"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b), path
        for k in a:
            _close(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _close(x, y, f"{path}[{i}]")
    elif isinstance(a, float) or isinstance(b, float):
        assert a is not None and b is not None and math.isclose(a, b, rel_tol=0, abs_tol=1e-9), path
    else:
        assert a == b, path


def test_golden_metrics(tmp_path):
    fx = pf.run_chain(tmp_path, 10, load_phasec, candidates="B1,R2,M8")
    assert _run11(fx, 50) == 0
    got = _metrics(fx["work"])
    # hand check: R2 recall on S1:all (direct truth, V-go) from the scores and the tables
    out = fx["work"] / "phasec"
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    calls = {r["seq_sha256"]: r["call"] for r in truth_table.read_tsv(out / "scores.tsv.gz")
             if r["split_id"] == "S1" and r["variant"] == "V-go" and r["candidate"] == "R2"}  # fmt: skip
    pos = [h for h, r in table.items() if r["origin"] == "go" and r["roles"] == "train"
           and r["class"] == "pos" and r["homology_only"] == "no"]  # fmt: skip
    r2 = got["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]["all"]["V-go"]["R2"]["recall"]
    assert r2["value"] == pytest.approx(sum(calls[h] == "1" for h in pos) / len(pos), abs=1e-12)
    if os.environ.get("PHASEC_WRITE_GOLDEN") == "1":
        GOLDEN.write_bytes(
            gzip.compress(json.dumps(got, indent=1, sort_keys=True).encode(), mtime=0)
        )
        pytest.skip("golden file written")
    _close(got, json.loads(gzip.decompress(GOLDEN.read_bytes())))
    run = json.loads((out / "evaluate_run.json").read_text())
    assert run["library_versions"]["numpy"] and run["library_versions"]["sklearn"]
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_evaluate.py -q`
Expected: collection error `ModuleNotFoundError: No module named 'findings'`.

- [ ] **Step 3: Write `analysis/step1_compare/phasec/findings.py`**

```python
"""Estimate label and machine-checked findings (Phase C spec 4, rulings C-8, C-11).

All functions read values that 11_evaluate.py has already put into metrics.json, so every
boolean in findings.json can be traced to a number there.
"""

ESTIMATE_HALF_WIDTH = 0.10
MIN_DIRECT_POSITIVES = 20  # owner decision 2026-10-01 (ruling C-8)
SATURATION_AUC = 0.99
LABEL_CANDIDATES = ("R2", "M8", "M35", "M8-C", "M35-C", "H")  # ruling C-8
ML = ("M8", "M35", "M8-C", "M35-C", "H")
COMPARATORS = ("B1", "R2")


def floor_met(n_direct_positives: int) -> bool:
    """The test set has at least MIN_DIRECT_POSITIVES direct-evidence positives."""
    return n_direct_positives >= MIN_DIRECT_POSITIVES


def estimate_label(recall_half_widths: dict, n_direct_positives: int) -> str:
    """`estimate` when the recall half-width is <= 0.10 for every listed candidate and the test
    set has at least 20 direct-evidence positives, else `smoke test` (also when a half-width is
    not defined, for example without positives). The count floor stops a zero-width interval
    of a small set (for example 7 of 7 positives found in every resample) from passing."""
    values = [recall_half_widths.get(c) for c in LABEL_CANDIDATES]
    widths_ok = all(v is not None and v <= ESTIMATE_HALF_WIDTH + 1e-12 for v in values)
    if widths_ok and floor_met(n_direct_positives):
        return "estimate"
    return "smoke test"


def finding_a(b1_auc: dict) -> dict:
    """(a) B1 does not saturate: ROC-AUC below 0.99 on every listed test set."""
    holds = bool(b1_auc) and all(v is not None and v < SATURATION_AUC for v in b1_auc.values())
    return {"threshold": SATURATION_AUC, "b1_roc_auc": b1_auc, "holds": holds}


def beats(diff: dict) -> bool:
    """The comparator's N-sec FPR minus the ML candidate's: the interval lies above 0."""
    return diff.get("lo") is not None and diff["lo"] > 0


def finding_b(diffs: dict) -> dict:
    """(b) for one test set: diffs[ml][comparator] = {value, lo, hi}. ML beats B1 and R2 when
    both intervals lie above 0."""
    out = {}
    for ml in ML:
        if ml not in diffs:
            continue
        flags = {f"beats_{c}": beats(diffs[ml].get(c, {})) for c in COMPARATORS}
        out[ml] = {
            **{c: diffs[ml].get(c) for c in COMPARATORS},
            **flags,
            "holds": all(flags.values()),
        }
    return {"candidates": out, "holds_for": [m for m in out if out[m]["holds"]]}


def finding_c(per_test_set: dict) -> dict:
    """(c) the statement of (b) on every S2 test set. per_test_set[ts] = finding_b(...)."""
    names = sorted(per_test_set)
    holds_for = [m for m in ML if names and all(m in per_test_set[ts]["holds_for"] for ts in names)]
    return {"test_sets": per_test_set, "holds_for": holds_for}
```

- [ ] **Step 4: Write `analysis/step1_compare/phasec/agreement.py`**

```python
"""Score lookup for whole proteomes, agreement counts, named panel (Phase C spec 4).

score_source of a sequence (one value per hash):
- oof: the hash is in an S1 test fold (GO rows of the training sources and T-c rows); its
  score comes from the S1 model that did not train on it or its cluster.
- in_sample: the hash is in a FULL training table (V-go or V-kw) but in no S1 fold.
- final: the hash is in no training table; its score comes from the FULL model.
"""

PANEL = (
    # (name, lookup id, why) from the parent spec section 6; the lookup id is a UniProt
    # accession of the keyword set, a truth gene_id, or a proteome gene_id
    ("FLO1", "P32768", "1,537 aa GPI adhesin (C-terminal variant)"),
    ("SAG1", "P20840", "GPI wall protein"),
    ("CWP1", "P28319", "GPI wall protein"),
    ("CCW12", "Q12127", "GPI wall protein, 133 aa"),
    ("GAS1", "P22146", "GPI enzyme; ambiguous under D1"),
    ("PIR1", "Q03178", "wall protein without GPI; ambiguous under D1"),
    ("MSB2", "P32334", "PM-TM negative"),
    ("HKR1", "P41809", "PM-TM, 1,802 aa; surface evidence IBA only"),
    ("SUC2", "P00724", "secreted enzyme; ambiguous"),
    ("PHO5", "P00635", "secreted enzyme; P-ext wall"),
    ("ALS3", "Q59L12", "C. albicans GPI protein"),
    ("HWP1", "P46593", "C. albicans GPI protein"),
    ("SAP9", "Q59SU1", "C. albicans GPI protein"),
    ("ECM33 (C. albicans)", "A0A1D8PCY4", "P-ext wall; gpi_anchor=no in keywords"),
    ("ENO1 (C. albicans)", "CAL0000185645", "ambiguous stratum"),
    ("TDH3 (C. albicans)", "CAL0000197744", "ambiguous; internal term IBA only"),
    ("SOWgp58", "Q8NK60", "Onygenales Pro-rich surface protein"),
    ("SOWgp (RS)", "CIMG_04613-t26_1-p1", "C. immitis RS proteome gene CIMG_04613"),
    ("CTS1", "Q1E3R8", "Eurotiomycetes literature row"),
    ("CspA", "Q4WXC4", "literature row; P-ext wall in ASPFU"),
    ("cfmA", "Q4WLB9", "literature row; P-ext wall in ASPFU"),
    ("HSP60", "P50142", "moonlighting; never positive"),
)
SOURCE_ORDER = ("oof", "in_sample", "final")


def score_source(h: str, s1_hashes, full_train) -> str:
    """s1_hashes and full_train: containers of hashes (a dict of S1 folds works too)."""
    if h in s1_hashes:
        return "oof"
    if h in full_train:
        return "in_sample"
    return "final"


def find_panel_hashes(panel, members) -> dict[str, dict]:
    """lookup id -> {seq_sha256, set_ids}. A keyword-set accession wins over other sets."""
    found: dict[str, dict] = {}
    wanted = {p[1] for p in panel}
    for m in members:
        gid = m["gene_id"]
        if gid not in wanted:
            continue
        hit = found.setdefault(gid, {"seq_sha256": m["seq_sha256"], "set_ids": []})
        hit["set_ids"].append(m["set_id"])
        if m["set_id"] == "uniprot_kw":
            hit["seq_sha256"] = m["seq_sha256"]
    return found


def agreement_counts(rule_calls, ml_calls) -> dict:
    """Counts of (rule call, ML call) pairs."""
    out = {"rule1_ml1": 0, "rule1_ml0": 0, "rule0_ml1": 0, "rule0_ml0": 0}
    for r, m in zip(rule_calls, ml_calls, strict=True):
        out[f"rule{int(bool(r))}_ml{int(bool(m))}"] += 1
    return out
```

- [ ] **Step 5: Write `analysis/step1_compare/phasec/11_evaluate.py`**

```python
#!/usr/bin/env python3
"""Phase C step 4: metrics with cluster-bootstrap intervals, strata, agreement, findings.

Reads from $STEP1_WORKDIR/phasec/ (SHA-256 checked against the run JSONs): eval_table.tsv.gz,
eval_literature.tsv (08), split_members.tsv.gz, clusters.tsv.gz, max_identity.tsv.gz (09),
scores.tsv.gz (10). Reads phaseb/sequence_members.tsv.gz (its SHA-256 must equal
features_run.json) and sequence_sets.tsv (--sets) for the proteome sets.

Test sets: S1:all (out-of-fold, all folds pooled) and S1:<source>; S2-<source>:<source>;
S3-<clade>:clade, S3-<clade>:<source>, S3-Eurotiomycetes:literature. Truth `direct` =
homology_only no (headline); `all` = all non-IEA labels (beside it). Literature rows have
positives only: a truth block without negatives stores precision, PR-AUC, precision at recall
and precision at the rule's recall as null (spec 3.2: recall only). The bootstrap weights (bootstrap.py) are drawn once per test set and truth and
serve every candidate, variant and stratum (paired).

Writes to $STEP1_WORKDIR/phasec/:

  metrics.json           values with 95% intervals per test set, truth, stratum, variant,
                         candidate; estimate or smoke-test label; variant effect; comparison
                         with the rule; prevalence table; calibration; agreement; named panel
  findings.json          findings (a), (b), (c) of spec 4 as booleans with their numbers
  proteome_calls.tsv.gz  one row per protein of every proteome set: call and score of every
                         candidate and variant, score_source
  evaluate_run.json      input hashes, seed, resamples, git commit, library versions

STOP (exit 2, no output): a stale input; scores that do not cover a test row; a T-c row in a
test set.
"""

import argparse
import csv
import gzip
import math
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import agreement
import bootstrap
import evalio
import findings
import manifest
import metrics
import numpy as np
import paths
import runinfo
import splits
import truth_table

OUTPUT_NAMES = ("metrics.json", "findings.json", "proteome_calls.tsv.gz", "evaluate_run.json")
METRICS_SCHEMA = "step1-phasec-metrics/1"
FINDINGS_SCHEMA = "step1-phasec-findings/1"
BINARY = ("recall", "precision", "fpr")
RECALL_LEVELS = (0.8, 0.9)
FPR_LEVEL = 0.01
SCORED = (
    "roc_auc",
    "pr_auc",
    "precision_at_recall_0.8",
    "precision_at_recall_0.9",
    "recall_at_fpr_0.01",
)
# undefined without negatives: precision is 1 (or n/a) by construction (spec 3.2, review I-1)
PRECISION_METRICS = (
    "precision",
    "pr_auc",
    "precision_at_recall_0.8",
    "precision_at_recall_0.9",
)
STRATA = (
    "all",
    "wall",
    "extracellular-only",
    "N-int",
    "N-sec",
    "PM-TM",
    "long",
    "identity_below_0.3",
)
TRUTHS = ("direct", "all")
PREVALENCES = (0.01, 0.03, 0.05, 0.10)
CALIBRATION_BINS = 10
QUANTILES = (0.1, 0.25, 0.5, 0.75, 0.9)
PROTEOME_KINDS = ("download", "site")
ONYGENALES_TC_TAXA = ("246410", "443226")
PROTEOME_BASE = (
    "set_id",
    "source_id",
    "gene_id",
    "seq_sha256",
    "class",
    "label",
    "origin",
    "score_source",
)
SCORE_COLUMNS = (
    "split_id", "fold", "variant", "candidate", "seq_sha256", "part", "score", "prob", "call",
)  # fmt: skip
LONG_CUTOFF = 1022


def num(x):
    """JSON number or None (NaN and inf become None)."""
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def summary(arr) -> dict:
    """Point value (row 0) and the bootstrap interval (rows 1..B)."""
    return {"value": num(arr[0]), **bootstrap.interval(arr[1:])}


@dataclass
class TestSet:
    name: str
    split: str
    kind: str
    rows: list


def read_scores(path: Path) -> dict:
    raw: dict = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if tuple(next(reader)) != SCORE_COLUMNS:
            raise evalio.StopError(f"{path}: unexpected header")
        for split, fold, variant, cand, h, _part, score, prob, call in reader:
            d = raw.setdefault((split, fold, variant, cand), ([], [], [], []))
            d[0].append(h)
            d[1].append(float(score) if score else math.nan)
            d[2].append(float(prob) if prob else math.nan)
            d[3].append(call == "1")
    out = {}
    for key, (hs, s, p, c) in raw.items():
        out[key] = {
            "index": {h: i for i, h in enumerate(hs)},
            "score": np.array(s),
            "prob": np.array(p),
            "call": np.array(c, dtype=bool),
        }
    return out


def test_sets(members, table, literature, identity, species_rows) -> list[TestSet]:
    lit = {r["seq_sha256"]: r for r in literature if r["seq_sha256"]}
    ident = {(r["split_id"], r["seq_sha256"]): r["below_0.3"] for r in identity}
    train_sources = sorted(r["source_id"] for r in species_rows if r["role"] == "train")
    rows_by_split: dict[str, list] = {}
    for m in members:
        if m["part"] not in ("test", "test_lit"):
            continue
        if m["origin"] == "tc":
            raise evalio.StopError(f"T-c row {m['seq_sha256']} is in test set {m['split_id']}")
        h = m["seq_sha256"]
        if m["part"] == "test_lit":
            info = {"label": "literature", "subset": "", "stratum": "literature", "d8_class": "",
                    "homology_only": "no", "internal_evidence_htp_only": "no",
                    "source_ids": "literature", "gene_ids": lit[h]["accession"],
                    "length": lit[h]["length"]}  # fmt: skip
        else:
            info = table[h]
        rows_by_split.setdefault(m["split_id"], []).append(
            {
                **{
                    k: info[k]
                    for k in (
                        "label",
                        "subset",
                        "stratum",
                        "d8_class",
                        "homology_only",
                        "internal_evidence_htp_only",
                        "source_ids",
                        "gene_ids",
                    )
                },
                "length": int(info["length"]),
                "seq_sha256": h,
                "split": m["split_id"],
                "fold": m["fold"],
                "class": m["class"],
                "cluster": m["cluster_id"],
                "part": m["part"],
                "below": ident.get((m["split_id"], h), ""),
            }  # fmt: skip
        )
    out = []
    for split in sorted(rows_by_split):
        rows = rows_by_split[split]
        go_rows = [r for r in rows if r["part"] == "test"]
        lit_rows = [r for r in rows if r["part"] == "test_lit"]
        sources = sorted({s for r in go_rows for s in r["source_ids"].split(",")})
        if split == "S1":
            out.append(TestSet("S1:all", split, "pooled", go_rows))
            for s in train_sources:
                out.append(
                    TestSet(
                        f"S1:{s}",
                        split,
                        "species",
                        [r for r in go_rows if s in r["source_ids"].split(",")],
                    )
                )
        elif split.startswith("S2-"):
            out.append(TestSet(f"{split}:{split[3:]}", split, "species", go_rows))
        else:
            out.append(TestSet(f"{split}:clade", split, "clade", go_rows))
            for s in sources:
                out.append(
                    TestSet(
                        f"{split}:{s}",
                        split,
                        "species",
                        [r for r in go_rows if s in r["source_ids"].split(",")],
                    )
                )
            if lit_rows:
                out.append(TestSet(f"{split}:literature", split, "literature", lit_rows))
    return out


def gather(scores, rows, variant: str, cand: str):
    """score, prob, call arrays for the rows (each row knows its split and fold)."""
    n = len(rows)
    s, p, c = np.full(n, math.nan), np.full(n, math.nan), np.zeros(n, dtype=bool)
    for i, r in enumerate(rows):
        d = scores.get((r["split"], r["fold"], variant, cand))
        j = None if d is None else d["index"].get(r["seq_sha256"])
        if j is None:
            raise evalio.StopError(
                f"no {cand} {variant} score for {r['seq_sha256']} in {r['split']} fold {r['fold']}"
            )
        s[i], p[i], c[i] = d["score"][j], d["prob"][j], d["call"][j]
    return s, p, c


def strata_masks(rows, y) -> dict:
    def has(col, value):
        return np.array([value in r[col].split(",") for r in rows], dtype=bool)

    out = {
        "all": np.ones(len(rows), dtype=bool),
        "wall": y & has("subset", "wall"),
        "extracellular-only": y & has("subset", "extracellular-only"),
        "N-int": ~y & has("stratum", "N-int"),
        "N-sec": ~y & has("stratum", "N-sec"),
        "PM-TM": ~y & has("stratum", "PM-TM"),
        "long": np.array([r["length"] > LONG_CUTOFF for r in rows], dtype=bool),
    }
    if any(r["below"] for r in rows):
        out["identity_below_0.3"] = np.array([r["below"] == "yes" for r in rows], dtype=bool)
    return out


def evaluate_truth(ts, truth, scores, cands, n_resamples, seed):
    rows = [
        r for r in ts.rows
        if r["class"] in ("pos", "neg") and (truth == "all" or r["homology_only"] == "no")
    ]  # fmt: skip
    y = np.array([r["class"] == "pos" for r in rows], dtype=bool)
    W = bootstrap.cluster_weights([r["cluster"] for r in rows], n_resamples,
                                  bootstrap.seed_for(seed, f"{ts.name}|{truth}"))  # fmt: skip
    Wx = np.vstack([np.ones((1, len(rows))), W])
    no_negatives = not bool((~y).any())
    masks = strata_masks(rows, y)
    data = {(v, c): gather(scores, rows, v, c) for v in evalio.VARIANTS for c in cands}
    arrays: dict = {}
    block = {"n": {}, "metrics": {}, "variant_effect": {}}
    for stratum, mask in masks.items():
        Wm, ym = Wx[:, mask], y[mask]
        block["n"][stratum] = {"pos": int(ym.sum()), "neg": int((~ym).sum())}
        per = block["metrics"][stratum] = {}
        for v in evalio.VARIANTS:
            per[v] = {}
            for c in cands:
                s, _, call = data[(v, c)]
                vals = {"recall": metrics.recall(Wm, ym, call[mask]),
                        "precision": metrics.precision(Wm, ym, call[mask]),
                        "fpr": metrics.fpr(Wm, ym, call[mask])}  # fmt: skip
                if not np.isnan(s).any():
                    vals.update(metrics.scored_summary(Wm, ym, s[mask], RECALL_LEVELS, FPR_LEVEL))
                if no_negatives:
                    for name in PRECISION_METRICS:
                        if name in vals:
                            vals[name] = np.full(len(Wm), np.nan)
                for name, arr in vals.items():
                    arrays[(stratum, v, c, name)] = arr
                per[v][c] = {name: summary(arr) for name, arr in vals.items()}
        block["variant_effect"][stratum] = {
            c: {
                name: summary(
                    arrays[(stratum, "V-kw", c, name)] - arrays[(stratum, "V-go", c, name)]
                )
                for name in per["V-go"][c]
            }
            for c in cands
        }
    block["vs_rule"], block["prevalence"], block["nsec_fpr_at_rule_recall"] = {}, {}, {}
    for v in evalio.VARIANTS:
        rule_rec = arrays.get(("all", v, "R2", "recall"))
        rule_fpr = arrays.get(("all", v, "R2", "fpr"))
        block["vs_rule"][v] = {}
        block["prevalence"][v] = {
            c: {
                f"{pi}": summary(metrics.precision_at_prevalence(
                    arrays[("all", v, c, "recall")], arrays[("all", v, c, "fpr")], pi))
                for pi in PREVALENCES
            }
            for c in cands
        }  # fmt: skip
        if rule_rec is None:
            continue
        nsec = masks["N-sec"]
        comp = {"R2": arrays[("N-sec", v, "R2", "fpr")]}
        for c in cands:
            s = data[(v, c)][0]
            if c == "R2" or np.isnan(s).any():
                continue
            if c in findings.ML:
                prec = metrics.precision_at_recall(Wx, y, s, rule_rec)
                if no_negatives:
                    prec = np.full(len(Wx), np.nan)
                block["vs_rule"][v][c] = {
                    "precision_at_rule_recall": summary(prec),
                    "recall_at_rule_fpr": summary(metrics.recall_at_fpr(Wx, y, s, rule_fpr)),
                }  # fmt: skip
            comp[c] = metrics.fpr_at_recall(Wx, y, s, rule_rec, nsec)
        block["nsec_fpr_at_rule_recall"][v] = {
            "fpr": {c: summary(a) for c, a in comp.items()},
            "diff": {
                ml: {k: summary(comp[k] - comp[ml]) for k in findings.COMPARATORS if k in comp}
                for ml in findings.ML
                if ml in comp
            },
        }
    return block, rows, y, data


def calibration(rows, y, data, cands) -> dict:
    out = {}
    for (v, c), (_, prob, _) in data.items():
        if c not in findings.ML or np.isnan(prob).any():
            continue
        bins = np.clip((prob * CALIBRATION_BINS).astype(int), 0, CALIBRATION_BINS - 1)
        out.setdefault(v, {})[c] = {
            "brier": num(metrics.brier(np.ones(len(y)), y, prob)[0]),
            "bins": [
                {"lo": k / CALIBRATION_BINS, "hi": (k + 1) / CALIBRATION_BINS,
                 "n": int((bins == k).sum()),
                 "mean_prob": num(prob[bins == k].mean()) if (bins == k).any() else None,
                 "frac_pos": num(y[bins == k].mean()) if (bins == k).any() else None}
                for k in range(CALIBRATION_BINS)
            ],
        }  # fmt: skip
    return out


def excluded_views(ts, scores, cands) -> dict:
    amb = [r for r in ts.rows if r["class"] == "excluded" and "ambiguous" in r["label"].split(",")]
    out = {}
    for name, rows in (("ambiguous", amb),
                       ("ambiguous_htp_only", [r for r in amb if r["internal_evidence_htp_only"] == "yes"])):  # fmt: skip
        view = {"n": len(rows)}
        for v in evalio.VARIANTS:
            for c in cands:
                if not rows:
                    continue
                s, _, call = gather(scores, rows, v, c)
                entry = {"call_fraction": num(call.mean())}
                if not np.isnan(s).any():
                    entry["score_quantiles"] = {f"{q}": num(np.quantile(s, q)) for q in QUANTILES}
                view.setdefault(v, {})[c] = entry
        out[name] = view
    lists = {}
    for name, test in (("pm-unresolved", lambda r: "pm-unresolved" in r["stratum"].split(",")),
                       ("P-gpi", lambda r: "P-gpi" in r["d8_class"].split(","))):  # fmt: skip
        picked = [r for r in ts.rows if test(r)]
        entries = []
        for r in picked:
            calls = {}
            for v in evalio.VARIANTS:
                for c in cands:
                    calls[f"{c}|{v}"] = int(gather(scores, [r], v, c)[2][0])
            entries.append(
                {"gene_ids": r["gene_ids"], "seq_sha256": r["seq_sha256"], "calls": calls}
            )
        lists[name] = entries
    out["lists"] = lists
    return out


def lookup_tables(members):
    s1 = {m["seq_sha256"]: m["fold"] for m in members
          if m["split_id"] == "S1" and m["part"] in ("test", "test_tc")}  # fmt: skip
    full_train = {m["seq_sha256"] for m in members
                  if m["split_id"] == "FULL" and m["part"] in ("train", "train_tc")}  # fmt: skip
    return s1, full_train


def hash_score(scores, s1, h, v, c):
    key = ("S1", s1[h], v, c) if h in s1 else ("FULL", "0", v, c)
    d = scores[key]
    j = d["index"][h]
    return d["score"][j], d["call"][j]


def proteome_rows(seq_members, set_ids, table, scores, s1, full_train, cands):
    rows = []
    for m in seq_members:
        if m["set_id"] not in set_ids:
            continue
        h = m["seq_sha256"]
        t = table.get(h, {})
        row = {"set_id": m["set_id"], "source_id": m["source_id"], "gene_id": m["gene_id"],
               "seq_sha256": h, "class": t.get("class", ""), "label": t.get("label", ""),
               "origin": t.get("origin", ""),
               "score_source": agreement.score_source(h, s1, full_train)}  # fmt: skip
        for c in cands:
            for v in evalio.VARIANTS:
                s, call = hash_score(scores, s1, h, v, c)
                row[f"call_{c}_{v}"] = "1" if call else "0"
                row[f"score_{c}_{v}"] = "" if math.isnan(s) else format(float(s), ".10g")
        rows.append(row)
    return rows


def proteome_columns(cands) -> tuple:
    return PROTEOME_BASE + tuple(
        f"{kind}_{c}_{v}" for c in cands for v in evalio.VARIANTS for kind in ("call", "score")
    )


def agreement_block(prows, table, scores, s1, cands) -> dict:
    out = {"proteomes": {}, "truth": {}}
    ml = [c for c in cands if c in findings.ML]
    for set_id in sorted({r["set_id"] for r in prows}):
        rows = [r for r in prows if r["set_id"] == set_id]
        out["proteomes"][set_id] = {
            v: {c: agreement.agreement_counts([r[f"call_R2_{v}"] == "1" for r in rows],
                                              [r[f"call_{c}_{v}"] == "1" for r in rows])
                for c in ml}
            for v in evalio.VARIANTS
        }  # fmt: skip
    for cls in ("pos", "neg"):
        hashes = [h for h, t in table.items() if t["origin"] == "go" and t["class"] == cls]
        out["truth"][cls] = {}
        for v in evalio.VARIANTS:
            rule = [hash_score(scores, s1, h, v, "R2")[1] for h in hashes]
            out["truth"][cls][v] = {
                c: agreement.agreement_counts(rule, [hash_score(scores, s1, h, v, c)[1] for h in hashes])
                for c in ml
            }  # fmt: skip
    return out


def named_panel(seq_members, table, literature, scores, s1, full_train, cands) -> list:
    panel = list(agreement.PANEL)
    have = {p[1] for p in panel}
    for r in literature:
        if r["lit_class"] == "hard_negative" and r["accession"] not in have:
            panel.append((r["gene"], r["accession"], "literature hard_negative row (ruling C-9)"))
    found = agreement.find_panel_hashes(panel, seq_members)
    lit = {r["accession"]: r for r in literature}
    out = []
    for name, pid, why in panel:
        hit = found.get(pid)
        entry = {"name": name, "id": pid, "why": why, "found": hit is not None}
        if hit:
            h = hit["seq_sha256"]
            t = table.get(h, {})
            entry.update({"seq_sha256": h, "set_ids": sorted(set(hit["set_ids"])),
                          "class": t.get("class", ""), "label": t.get("label", ""),
                          "stratum": t.get("stratum", ""),
                          "literature_class": lit.get(pid, {}).get("lit_class", ""),
                          "score_source": agreement.score_source(h, s1, full_train),
                          "calls": {}, "scores": {}})  # fmt: skip
            for c in cands:
                for v in evalio.VARIANTS:
                    s, call = hash_score(scores, s1, h, v, c)
                    entry["calls"][f"{c}|{v}"] = int(call)
                    entry["scores"][f"{c}|{v}"] = num(s)
        out.append(entry)
    return out


def build_findings(test_blocks: dict) -> dict:
    s2 = sorted(n for n, b in test_blocks.items() if b["split"].startswith("S2-"))
    b1 = {}
    for name in ["S1:all", *s2]:
        if name in test_blocks:
            m = test_blocks[name]["truth"]["direct"]["metrics"]["all"]["V-go"].get("B1", {})
            b1[name] = m.get("roc_auc", {}).get("value")

    def diffs(name):
        block = test_blocks[name]["truth"]["direct"]["nsec_fpr_at_rule_recall"].get("V-go", {})
        return block.get("diff", {})

    return {
        "schema": FINDINGS_SCHEMA,
        "a_b1_not_saturated": findings.finding_a(b1),
        "b_ml_beats_b1_and_r2_on_nsec_s1": {
            "test_set": "S1:all",
            **findings.finding_b(diffs("S1:all") if "S1:all" in test_blocks else {}),
        },
        "c_same_under_s2": findings.finding_c({n: findings.finding_b(diffs(n)) for n in s2}),
    }


def run(work: Path, sets_path: Path, species_path: Path, n_resamples: int, arguments=()):
    work = Path(work)
    out = evalio.out_dir(work)
    build = evalio.require_current(
        out, "build_run.json", ("eval_table.tsv.gz", "eval_literature.tsv"), "08"
    )
    split_log = evalio.require_current(
        out,
        "splits_run.json",
        ("split_members.tsv.gz", "clusters.tsv.gz", "max_identity.tsv.gz"),
        "09",
    )
    score_log = evalio.require_current(
        out, "scores_run.json", ("scores.tsv.gz",), "10_fit_and_score.py"
    )
    if (
        score_log["input_sha256"]["split_members.tsv.gz"]
        != split_log["outputs_sha256"]["split_members.tsv.gz"]
    ):
        raise evalio.StopError(
            "scores.tsv.gz was made from another split_members.tsv.gz; re-run 10"
        )
    fr = evalio.read_json(work / "phaseb" / "features_run.json")
    seq_members_path = work / "phaseb" / "sequence_members.tsv.gz"
    if manifest.sha256_file(seq_members_path) != fr["input_sha256"]["sequence_members.tsv.gz"]:
        raise evalio.StopError("phaseb/sequence_members.tsv.gz differs from the file 07 read")
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    literature = truth_table.read_tsv(out / "eval_literature.tsv")
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    splits.check_test_truth(members)
    identity = truth_table.read_tsv(out / "max_identity.tsv.gz")
    if manifest.sha256_file(species_path) != build["input_sha256"]["species.tsv"]:
        raise evalio.StopError(f"{species_path} differs from the species.tsv that 08 used")
    species = truth_table.read_tsv(species_path)
    scores = read_scores(out / "scores.tsv.gz")
    cands = tuple(score_log["candidates"])
    blocks = {}
    for ts in test_sets(members, table, literature, identity, species):
        res = {"split": ts.split, "kind": ts.kind, "truth": {}}
        for truth in TRUTHS:
            block, rows, y, data = evaluate_truth(
                ts, truth, scores, cands, n_resamples, evalio.SEED
            )
            res["truth"][truth] = block
            if truth == "direct" and ts.split.startswith("S2-"):
                res["calibration"] = calibration(rows, y, data, cands)
        direct = res["truth"]["direct"]["metrics"]["all"]["V-go"]
        hw = {
            c: bootstrap.half_width(direct[c]["recall"])
            for c in findings.LABEL_CANDIDATES
            if c in direct
        }
        res["recall_half_width"] = hw
        res["fpr_half_width"] = {c: bootstrap.half_width(direct[c]["fpr"]) for c in hw}
        n_pos = res["truth"]["direct"]["n"]["all"]["pos"]
        res["n_direct_positives"] = n_pos
        res["floor_met"] = findings.floor_met(n_pos)
        res["label"] = findings.estimate_label(hw, n_pos)
        for key, src in (
            ("max_recall_half_width", hw),
            ("max_fpr_half_width", res["fpr_half_width"]),
        ):
            vals = list(src.values())
            res[key] = None if not vals or None in vals else max(vals)
        res["zero_width_recall_interval"] = sorted(
            c for c in hw if direct[c]["recall"]["lo"] is not None
            and direct[c]["recall"]["lo"] == direct[c]["recall"]["hi"]
        )  # fmt: skip
        res.update(excluded_views(ts, scores, cands))
        blocks[ts.name] = res
    s1, full_train = lookup_tables(members)
    set_rows = truth_table.read_tsv(sets_path)
    proteome_sets = {r["set_id"] for r in set_rows if r["kind"] in PROTEOME_KINDS}
    seq_members = truth_table.read_tsv(seq_members_path)
    prows = proteome_rows(seq_members, proteome_sets, table, scores, s1, full_train, cands)
    onygenales = Counter(
        f"{m['split_id']}|{m['fold']}" for m in members
        if m["part"] == "train_tc"
        and set(table[m["seq_sha256"]]["taxon_ids"].split(",")) & set(ONYGENALES_TC_TAXA)
    )  # fmt: skip
    metrics_json = {
        "schema": METRICS_SCHEMA,
        "settings": {
            "n_resamples": n_resamples,
            "seed": evalio.SEED,
            "ci_level": 95,
            "estimate_half_width": findings.ESTIMATE_HALF_WIDTH,
            "estimate_min_direct_positives": findings.MIN_DIRECT_POSITIVES,
            "saturation_auc": findings.SATURATION_AUC,
            "prevalences": list(PREVALENCES),
            "recall_levels": list(RECALL_LEVELS),
            "fpr_level": FPR_LEVEL,
            "long_cutoff": LONG_CUTOFF,
            "identity_cutoff": splits.IDENTITY_CUT,
            "calibration_bins": CALIBRATION_BINS,
            "candidates": list(cands),
            "variants": list(evalio.VARIANTS),
        },  # fmt: skip
        "test_sets": blocks,
        "agreement": agreement_block(prows, table, scores, s1, cands),
        "score_sources": {
            s: dict(sorted(Counter(r["score_source"] for r in prows if r["set_id"] == s).items()))
            for s in sorted(proteome_sets)
        },
        "named_panel": named_panel(seq_members, table, literature, scores, s1, full_train, cands),
        "context": {"onygenales_tc_rows_in_vkw_training": dict(sorted(onygenales.items()))},
        "tc_removed": split_log["tc_removed_by_split_fold_rule"],
    }
    findings_json = build_findings(blocks)
    log = {
        "all_sources": build["all_sources"],
        "truth_set_sha256": build["truth_set_sha256"],
        "input_sha256": {
            "scores.tsv.gz": score_log["outputs_sha256"]["scores.tsv.gz"],
            "split_members.tsv.gz": split_log["outputs_sha256"]["split_members.tsv.gz"],
            "eval_table.tsv.gz": build["outputs_sha256"]["eval_table.tsv.gz"],
            "sequence_members.tsv.gz": fr["input_sha256"]["sequence_members.tsv.gz"],
            "sequence_sets.tsv": manifest.sha256_file(sets_path),
        },
        "seed": evalio.SEED,
        "n_resamples": n_resamples,
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    columns = proteome_columns(cands)

    evalio.write_outputs(
        out,
        {
            "metrics.json": lambda p: runinfo.write_json(p, metrics_json),
            "findings.json": lambda p: runinfo.write_json(p, findings_json),
            "proteome_calls.tsv.gz": lambda p: truth_table.write_tsv(p, columns, prows),
        },
        "evaluate_run.json",
        log,
    )
    return metrics_json, findings_json


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sets", default=str(paths.STEP1_DIR / "sequence_sets.tsv"))
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--n-resamples", type=int, default=bootstrap.N_RESAMPLES)
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        m, f = run(work, Path(args.sets), Path(args.species), args.n_resamples,
                   list(argv) if argv is not None else sys.argv[1:])  # fmt: skip
    except (evalio.StopError, splits.SplitError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for name, block in m["test_sets"].items():
        print(f"{name}\t{block['label']}")
    print(f"finding_a_holds={f['a_b1_not_saturated']['holds']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 6: Generate the golden file once, then run the tests**

Run: `PHASEC_WRITE_GOLDEN=1 PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_evaluate.py::test_golden_metrics -q -rs`
Expected: `1 skipped` with reason `golden file written`; the file `tests/step1_compare/fixtures/phasec/golden_metrics.json.gz` exists (prototype: 54,129 bytes, SHA-256 `233a0333518551229725e5e7da97a717a56f0e8ec46665da5f9811ea7b9883c8` with numpy 2.4.2 and scikit-learn 1.8.0 on c01; regenerated after the plan-review fixes, which change the C and `t` grids, the label rule and the literature precision). The same SHA-256 means the code was transcribed exactly; another value is not an error by itself, because the test compares numbers with a tolerance of 1e-9. Before you generate the file, the hand check in `test_golden_metrics` (R2 recall on S1:all from the raw score rows) must pass.

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_evaluate.py -q`
Expected: `13 passed`.

- [ ] **Step 7: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/findings.py analysis/step1_compare/phasec/agreement.py \
  analysis/step1_compare/phasec/11_evaluate.py tests/step1_compare/test_phasec_evaluate.py \
  tests/step1_compare/fixtures/phasec/golden_metrics.json.gz
git commit -m "phasec: metrics.json, findings.json, proteome calls, named panel (11)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 10: Report (12)

**Files:**
- Create: `analysis/step1_compare/phasec/12_report.py`
- Test: `tests/step1_compare/test_phasec_report.py`

**Interfaces:**
- Consumes: `evalio.*` (Task 1); `metrics.json`, `findings.json`, `evaluate_run.json` (Task 9).
- Produces: `12_report.py` with `OUTPUT_NAMES = ("report.md", "report_run.json")`, `NUMBER` (regex of a number token), `HEADLINE`, `LITERATURE_HEADLINE = ("recall",)`, `CAVEATS` (three fixed lines: reviews M-3, M-4, M-7), `CANNOT_SHOW`, `fmt(x) -> str` (int: thousands separator; float: 3 decimals; bool: yes/no; None: n/a), `ci(m) -> str`, `allowed_numbers(*objs) -> set[str]`, `unsupported_numbers(text, allowed) -> list[str]`, `render(m, f) -> str`, `run(work, arguments=())`, `main(argv=None)`.

The number check: a token is `-?digits[,digits][.digits]` not preceded by a letter, digit, `_`, `.` or `-` and not followed by a letter or digit, so identifiers such as `M35`, `S2-Calb_CGD`, `C-11` and `precision_at_recall_0.8` hold no number. The allowed set is every numeric leaf of `metrics.json` and `findings.json` printed with `fmt` (and integers also without separators) plus the number tokens inside their keys and strings (for example the fold in `S1|0`). The report has 2,390 lines on the fixture. Prose of the report must not hold a number token that is not in the JSON files: the caveat lines and the literature note name the spec in words, not by section number.

The literature section (positives only) prints the recall column only and no comparison or prevalence table (review I-1). The report names no headline ML candidate: findings (b) and (c) read "Any ML candidate" and list `holds_for` (ruling C-13).

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_report.py`:

```python
"""12_report.py: every number comes from metrics.json or findings.json; findings are quoted."""

import json
import re

import pytest

pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import phasec_fixture as pf  # noqa: E402
from conftest import load_phasec  # noqa: E402


@pytest.fixture
def reported(phasec_chain, tmp_path):
    fx = pf.copy_work(phasec_chain, tmp_path)
    assert load_phasec("12_report").main(["--work-dir", str(fx["work"])]) == 0
    out = fx["work"] / "phasec"
    return (
        fx,
        (out / "report.md").read_text(),
        json.loads((out / "metrics.json").read_text()),
        json.loads((out / "findings.json").read_text()),
    )


def test_report_numbers_come_from_metrics(reported):
    _, text, m, f = reported
    rep = load_phasec("12_report")
    allowed = rep.allowed_numbers(m, f)
    assert rep.NUMBER.findall(text)  # the check sees numbers at all
    assert rep.unsupported_numbers(text, allowed) == []
    first = re.search(r"\| (\d\.\d{3}) \[", text).group(1)
    altered = text.replace(f"| {first} [", "| 0.4242 [", 1)
    assert rep.unsupported_numbers(altered, allowed) == ["0.4242"]
    assert rep.unsupported_numbers("negative -0.4242 here", allowed) == ["-0.4242"]
    assert rep.unsupported_numbers("ruling C-11, model M35, set S2-x, R2", allowed) == []


def test_report_stops_on_a_number_without_source(phasec_chain, tmp_path, monkeypatch, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    rep = load_phasec("12_report")
    real = rep.render
    monkeypatch.setattr(rep, "render", lambda m, f: real(m, f) + "\nRecall is 0.4242.\n")
    assert rep.main(["--work-dir", str(fx["work"])]) == 2
    assert "0.4242" in capsys.readouterr().err
    assert not (fx["work"] / "phasec" / "report.md").exists()


def test_report_quotes_findings_and_labels(reported):
    _, text, m, f = reported
    holds = "yes" if f["a_b1_not_saturated"]["holds"] else "no"
    assert f"holds = {holds}." in text
    # ruling C-13: no headline ML candidate; (b) and (c) say "any" and list holds_for
    b = ", ".join(f["b_ml_beats_b1_and_r2_on_nsec_s1"]["holds_for"]) or "no ML candidate"
    want = "- (b) Any ML candidate beats B1 and R2 on N-sec FPR at the recall of R2 (S1:all): "
    assert f"{want}holds for {b}." in text
    c = ", ".join(f["c_same_under_s2"]["holds_for"]) or "no ML candidate"
    assert f"- (c) Any ML candidate, the same on every S2 test set: holds for {c}." in text
    assert "What the data cannot show" in text and "assumed, not measured" in text
    for name, ts in m["test_sets"].items():
        assert f"## {name} ({ts['label']})" in text
    # every metric row of a test set carries the label of that test set
    section = text.split("## S1:all (")[1].split("\n## ")[0]
    label = m["test_sets"]["S1:all"]["label"]
    rows = [r for r in section.splitlines() if r.startswith("| ") and "[" in r]
    assert rows and all(f"| {label} |" in r for r in rows)


def test_report_states_the_caveats(reported):
    # review M-3, M-4, M-7: fixed caveat lines, not free text about accuracy
    _, text, m, _ = reported
    rep = load_phasec("12_report")
    caveats = text.split("## Caveats\n")[1].split("\n## ")[0]
    assert len(rep.CAVEATS) == 3
    for line in rep.CAVEATS:
        assert f"- {line}" in caveats
    for review in ("M-3", "M-4", "M-7"):
        assert review in caveats
    floor = m["settings"]["estimate_min_direct_positives"]
    assert f"at least {floor} direct-evidence positives" in text


def test_literature_section_reports_recall_only(reported):
    # review I-1: the literature set has no negatives; spec 3.2 asks for recall only
    _, text, m, _ = reported
    lit = m["test_sets"]["S3-Eurotiomycetes:literature"]
    section = text.split(f"## S3-Eurotiomycetes:literature ({lit['label']})")[1].split("\n## ")[0]
    assert "| Candidate | Variant | Label | recall |" in section
    for word in ("precision", "pr_auc", "roc_auc", "fpr", "prevalence"):
        assert word not in section, word
    # a set with negatives still shows the precision columns
    s1 = text.split("## S1:all (")[1].split("\n## ")[0]
    assert "| precision |" in s1 and "Precision at R2 recall" in s1


def test_stale_input_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    path = fx["work"] / "phasec" / "metrics.json"
    path.write_text(path.read_text().replace('"ci_level": 95', '"ci_level": 90'))
    assert load_phasec("12_report").main(["--work-dir", str(fx["work"])]) == 2
    assert "metrics.json differs from the SHA-256 in evaluate_run.json" in capsys.readouterr().err
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_report.py -q`
Expected: `2 failed, 4 errors`, each with `FileNotFoundError: ... phasec/12_report.py` (the four tests that use the fixture `reported` error).

- [ ] **Step 3: Write `analysis/step1_compare/phasec/12_report.py`**

```python
#!/usr/bin/env python3
"""Phase C step 5: report.md from metrics.json and findings.json (Phase C spec 4, E3).

Reads $STEP1_WORKDIR/phasec/metrics.json and findings.json (SHA-256 checked against
evaluate_run.json). Writes report.md and report_run.json to $STEP1_WORKDIR/phasec/.

Every number in report.md is a value of metrics.json or findings.json, printed by `fmt`
(integers with thousands separators, other numbers with 3 decimals). Before it writes, the
script checks this with `unsupported_numbers` and stops if a number has no source. The report
quotes the findings fields; it adds no free-text claim about accuracy. The literature test set
has no negatives, so its tables give recall only (spec 3.2). The fixed caveat lines (CAVEATS)
name known limits of the method.

STOP (exit 2, no output): stale inputs; a number in the report without a source.
"""

import argparse
import re
import sys
from pathlib import Path

import evalio
import paths
import runinfo

OUTPUT_NAMES = ("report.md", "report_run.json")
NUMBER = re.compile(r"(?<![\w.\-])-?\d[\d,]*(?:\.\d+)?(?![\w])")
HEADLINE = (
    "recall",
    "precision",
    "fpr",
    "roc_auc",
    "pr_auc",
    "precision_at_recall_0.8",
    "precision_at_recall_0.9",
    "recall_at_fpr_0.01",
)
LITERATURE_HEADLINE = ("recall",)  # positives only: recall is the one defined metric
CAVEATS = (
    "The pooled S1:all ROC-AUC and PR-AUC rank decision values from five fold models. Each "
    "fold model has its own score scale, so the pooled values mix scales (review M-3; the "
    "per-fold diagnostic is not computed).",
    "The ML decision threshold and the Platt scaling are fitted on the inner out-of-fold "
    "predictions and then applied to the model refitted on all training rows. The refitted "
    "model can have another score scale, so binary calls and calibrated probabilities can "
    "shift (review M-4; not measured).",
    "The report-number check is set membership: each printed number equals some value in "
    "metrics.json or findings.json. The check does not show that a number is in the correct "
    "cell (review M-7).",
)
CANNOT_SHOW = (
    "The whole-proteome prevalence of surface proteins. It is not measured; the prevalence "
    "table uses assumed values (ruling C-11).",
    "Basidiomycota performance. Cneo_H99_GOA and Umay_MYCMD are smoke tests, not the Q8 "
    "validation (owner decision on the Basidiomycota question).",
    "Whether SignalP under-calls C. immitis RS proteins. No truth set shows it; the literature "
    "rows give recall only.",
    "Performance on P-gpi proteins. P-gpi is a list until curated_gpi.tsv has literature rows "
    "(R-B).",
    "Variance from refitting the models. The intervals describe sampling of the test set only.",
)


def fmt(x) -> str:
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, int):
        return f"{x:,}"
    return f"{x:.3f}"


def ci(m: dict | None) -> str:
    if not m or m.get("value") is None:
        return "n/a"
    return f"{fmt(m['value'])} [{fmt(m['lo'])}, {fmt(m['hi'])}]"


def _leaves(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from _leaves(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _leaves(v)
    else:
        yield obj


def allowed_numbers(*objs) -> set[str]:
    out = set()
    for obj in objs:
        for x in _leaves(obj):
            if isinstance(x, bool) or x is None:
                continue
            if isinstance(x, str):
                out |= set(NUMBER.findall(x))
                continue
            if isinstance(x, int):
                out |= {f"{x:,}", str(x)}
            else:
                out.add(f"{x:.3f}")
    return out


def unsupported_numbers(text: str, allowed: set[str]) -> list[str]:
    return [t for t in NUMBER.findall(text) if t not in allowed]


def table(header, rows) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out + [""]


def render(m: dict, f: dict) -> str:
    st = m["settings"]
    cands, variants = st["candidates"], st["variants"]
    lines = [
        "# Step 1 Phase C: rule versus ML",
        "",
        f"Intervals: {fmt(st['ci_level'])}% cluster bootstrap, {fmt(st['n_resamples'])} "
        f"resamples, seed {st['seed']}. The intervals describe sampling of the test set only. "
        "The models are fitted once per split and are not refitted inside the bootstrap.",
        "Headline truth: direct evidence (`homology_only == no`). All non-IEA truth is beside it.",
        f"A test set is an estimate when the recall half-width is at most "
        f"{fmt(st['estimate_half_width'])} for R2 and every ML candidate (V-go, direct truth) "
        f"and the set has at least {fmt(st['estimate_min_direct_positives'])} direct-evidence "
        "positives; otherwise it is a smoke test.",
        "",
        "## Caveats",
        "",
        *[f"- {x}" for x in CAVEATS],
        "",
        "## Estimate or smoke test",
        "",
    ]
    rows = []
    for name, ts in m["test_sets"].items():
        n = ts["truth"]["direct"]["n"]["all"]
        rows.append([name, ts["label"], fmt(n["pos"]), fmt(n["neg"]), fmt(ts["floor_met"]),
                     fmt(ts["max_recall_half_width"]), fmt(ts["max_fpr_half_width"]),
                     ", ".join(ts["zero_width_recall_interval"]) or "none"])  # fmt: skip
    lines += table(["Test set", "Label", "Positives", "Negatives", "Count floor met",
                    "Max recall half-width", "Max FPR half-width", "Zero-width recall interval"],
                   rows)  # fmt: skip
    lines += ["A zero-width interval (every resample gives the same recall, for example 1 of 1) "
              "meets the half-width rule but carries no information about precision of the "
              "estimate. The count floor keeps such a small set a smoke test.", ""]  # fmt: skip
    for name, ts in m["test_sets"].items():
        lit = ts["kind"] == "literature"
        heads = LITERATURE_HEADLINE if lit else HEADLINE
        lines += [f"## {name} ({ts['label']})", ""]
        if lit:
            lines += [
                "Literature rows have positives only: the tables give recall only (Phase C spec). "
                "Context: the V-kw training of S1 and S2 keeps the C. immitis and "
                "C. posadasii T-c rows (counts in the context section).",
                "",
            ]
        for truth in ("direct", "all"):
            block = ts["truth"][truth]
            n = block["n"]["all"]
            lines += [f"Truth `{truth}`: {fmt(n['pos'])} positives, {fmt(n['neg'])} negatives.", ""]
            rows = []
            for v in variants:
                for c in cands:
                    mm = block["metrics"]["all"][v][c]
                    rows.append([c, v, ts["label"], *[ci(mm.get(k)) for k in heads]])
            lines += table(["Candidate", "Variant", "Label", *heads], rows)
        block = ts["truth"]["direct"]
        cols = LITERATURE_HEADLINE if lit else ("recall", "fpr", "roc_auc")
        rows = []
        for stratum, n in block["n"].items():
            for c in cands:
                mm = block["metrics"][stratum]["V-go"][c]
                rows.append([stratum, c, ts["label"], fmt(n["pos"]), fmt(n["neg"]),
                             *[ci(mm.get(k)) for k in cols]])  # fmt: skip
        lines += ["Strata (direct truth, V-go):", ""]
        lines += table(["Stratum", "Candidate", "Label", "Positives", "Negatives", *cols], rows)
        cols = LITERATURE_HEADLINE if lit else ("recall", "fpr", "roc_auc", "pr_auc")
        rows = []
        for c in cands:
            eff = block["variant_effect"]["all"][c]
            rows.append([c, ts["label"], *[ci(eff.get(k)) for k in cols]])
        lines += ["Variant effect, V-kw minus V-go (paired intervals; test-set sampling only):", ""]
        lines += table(["Candidate", "Label", *cols], rows)
        amb = ts["ambiguous"]
        lines += [f"Ambiguous genes: {fmt(amb['n'])} (score distribution only); with "
                  f"high-throughput-only internal evidence: {fmt(ts['ambiguous_htp_only']['n'])}.", ""]  # fmt: skip
        if lit:
            continue  # no negatives: no comparison at the rule's FPR, no prevalence table
        rows = []
        for v in variants:
            for c, vr in block["vs_rule"].get(v, {}).items():
                rows.append(
                    [
                        c,
                        v,
                        ts["label"],
                        ci(vr["precision_at_rule_recall"]),
                        ci(vr["recall_at_rule_fpr"]),
                    ]
                )
        if rows:
            lines += [
                "ML against the rule R2 (precision at the rule's recall, recall at the rule's FPR):",
                "",
            ]
            lines += table(
                ["Candidate", "Variant", "Label", "Precision at R2 recall", "Recall at R2 FPR"],
                rows,
            )
        rows = []
        for c in cands:
            pv = block["prevalence"]["V-go"][c]
            rows.append([c, ts["label"], *[ci(pv[k]) for k in pv]])
        prev = [fmt(p) for p in st["prevalences"]]
        lines += ["Precision at assumed prevalence (assumed, not measured), V-go:", ""]
        lines += table(["Candidate", "Label", *prev], rows)
        if "calibration" in ts:
            rows = []
            for v, per in ts["calibration"].items():
                for c, cal in per.items():
                    rows.append([c, v, ts["label"], fmt(cal["brier"])])
            lines += [
                "Calibration after Platt scaling on inner out-of-fold predictions (Brier score):",
                "",
            ]
            lines += table(["Candidate", "Variant", "Label", "Brier"], rows)
            rows = []
            for v, per in ts["calibration"].items():
                for c, cal in per.items():
                    cells = [f"{fmt(b['frac_pos'])} ({fmt(b['n'])})" for b in cal["bins"]]
                    rows.append([c, v, ts["label"], *cells])
            bins = next(iter(next(iter(ts["calibration"].values())).values()))["bins"]
            head = [f"{fmt(b['lo'])} to {fmt(b['hi'])}" for b in bins]
            lines += ["Reliability: observed positive fraction (count) per bin of the calibrated "
                      "probability:", ""]  # fmt: skip
            lines += table(["Candidate", "Variant", "Label", *head], rows)
    lines += ["## Agreement of R2 with ML (proteome_calls.tsv.gz holds the protein IDs)", ""]
    rows = []
    for set_id, per in m["agreement"]["proteomes"].items():
        for v, cc in per.items():
            for c, k in cc.items():
                rows.append([set_id, c, v, fmt(k["rule1_ml1"]), fmt(k["rule1_ml0"]),
                             fmt(k["rule0_ml1"]), fmt(k["rule0_ml0"])])  # fmt: skip
    lines += table(["Set", "ML", "Variant", "both", "rule only", "ML only", "neither"], rows)
    lines += ["The C. immitis RS row is read as held-out, but V-kw training of S1 and S2 keeps "
              "the C. immitis and C. posadasii T-c rows (context section).", ""]  # fmt: skip
    rows = []
    for cls, per in m["agreement"]["truth"].items():
        for v, cc in per.items():
            for c, k in cc.items():
                rows.append([cls, c, v, fmt(k["rule1_ml1"]), fmt(k["rule1_ml0"]),
                             fmt(k["rule0_ml1"]), fmt(k["rule0_ml0"])])  # fmt: skip
    lines += ["Agreement on the GO truth rows (all sources), by class:", ""]
    lines += table(["Class", "ML", "Variant", "both", "rule only", "ML only", "neither"], rows)
    lines += ["## Named panel", ""]
    rows = []
    for p in m["named_panel"]:
        calls = [fmt(p["calls"].get(f"{c}|V-go")) if p["found"] else "n/a" for c in cands]
        rows.append([p["name"], p["id"], fmt(p["found"]), p.get("label", ""), p.get("class", ""),
                     p.get("score_source", ""), *calls])  # fmt: skip
    lines += ["Calls under V-go (1 = called); scores and V-kw calls are in metrics.json.", ""]
    lines += table(["Protein", "ID", "Found", "Label", "Class", "Score source", *cands], rows)
    a, b, c3 = f["a_b1_not_saturated"], f["b_ml_beats_b1_and_r2_on_nsec_s1"], f["c_same_under_s2"]
    lines += ["## Findings (quoted from findings.json)", ""]
    lines += [f"- (a) B1 ROC-AUC below {fmt(a['threshold'])} on S1 and S2 (direct truth, V-go): "
              f"holds = {fmt(a['holds'])}."]  # fmt: skip
    lines += [f"  - {k}: {fmt(v)}" for k, v in a["b1_roc_auc"].items()]
    lines += ["Findings (b) and (c) are read per ML candidate; no headline ML candidate is named "
              "(ruling C-13). \"Any\" means at least one listed candidate.", ""]  # fmt: skip
    lines += [f"- (b) Any ML candidate beats B1 and R2 on N-sec FPR at the recall of R2 "
              f"({b['test_set']}): holds for {', '.join(b['holds_for']) or 'no ML candidate'}."]  # fmt: skip
    for cand, entry in b["candidates"].items():
        lines += [f"  - {cand}: B1 minus ML {ci(entry['B1'])}; R2 minus ML {ci(entry['R2'])}"]
    lines += [
        f"- (c) Any ML candidate, the same on every S2 test set: holds for "
        f"{', '.join(c3['holds_for']) or 'no ML candidate'}.",
        "",
    ]
    lines += ["## Context", ""]
    for key, n in m["context"]["onygenales_tc_rows_in_vkw_training"].items():
        lines += [f"- {key}: {fmt(n)} C. immitis and C. posadasii T-c rows in V-kw training."]
    lines += ["", "## What the data cannot show", ""]
    lines += [f"- {x}" for x in CANNOT_SHOW]
    return "\n".join(lines) + "\n"


def run(work: Path, arguments=()):
    out = evalio.out_dir(work)
    ev = evalio.require_current(out, "evaluate_run.json", ("metrics.json", "findings.json"), "11")
    m = evalio.read_json(out / "metrics.json")
    f = evalio.read_json(out / "findings.json")
    text = render(m, f)
    bad = unsupported_numbers(text, allowed_numbers(m, f))
    if bad:
        raise evalio.StopError(f"report numbers without a source in metrics.json: {bad[:5]}")
    log = {
        "all_sources": ev["all_sources"],
        "truth_set_sha256": ev["truth_set_sha256"],
        "input_sha256": {n: ev["outputs_sha256"][n] for n in ("metrics.json", "findings.json")},
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    evalio.write_outputs(
        out, {"report.md": lambda p: Path(p).write_text(text)}, "report_run.json", log
    )
    return text


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        text = run(work, list(argv) if argv is not None else sys.argv[1:])
    except (evalio.StopError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(f"report.md: {len(text.splitlines())} lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_report.py -q`
Expected: `6 passed`.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/12_report.py tests/step1_compare/test_phasec_report.py
git commit -m "phasec: report.md from metrics.json and findings.json with a number check (12)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 11: SLURM job C1 and the MMseqs2 wrapper

**Files:**
- Create: `analysis/step1_compare/phasec/09_cluster_and_split.sh`, `analysis/step1_compare/phasec/c1_evaluate.sh`
- Test: `tests/step1_compare/test_phasec_jobs.py`

**Interfaces:**
- Consumes: 09, 10, 11 (Tasks 7 to 9).
- Produces: environment contract: `PROJ_ROOT`, `STEP1_WORKDIR` (required); `STEP1_ENV_PY`, `STEP1_MMSEQS` (a binary; skips `module load`), `PHASEC_MMSEQS_MODULE` (default `MMseqs2/17-b804f`), `PHASEC_SPECIES`, `PHASEC_SETS`, `PHASEC_N_RESAMPLES` (default 2000), `PHASEC_CANDIDATES` (default all ten), `C1_STEPS` (default `09 10 11`). C1 prints `C1 step <step> wall_seconds=<s>` per step and copies these lines to `$STEP1_WORKDIR/phasec/logs/wall.<job id>.txt` from an EXIT trap, so the copy also happens after a failed step (review M-1). Both scripts request `-p epyc` and `--constraint=ryzen` (ruling C-15; owner rule: a job that runs an AVX2 tool must request a node feature that has AVX2). Test helpers `needs_avx2(text)`, `has_avx2_constraint(text)`, `uses_non_avx2_binary(text)`, `AVX2_FEATURES = ("ryzen", "milan", "genoa", "rome")`.

Sizing (global rule; no run time measured on epyc yet). The three steps depend on each other, so one job is the smallest job count. Estimate from c01 measurements (Facts): 09 about 19 min with 2 threads (less with 16); 10 about 22 units x 45 to 90 s single-threaded = 17 to 33 CPU-minutes with the old C grid, at most 1.5 x that with the six-value grid of ruling C-12 (estimate from the grid size, not measured), spread over 16 worker processes; 11 about 25 min single-threaded (28 s per 4,476 rows and 20 candidate-variant pairs for stratum `all`, scaled to all test sets, strata and truths). Total estimate: under 1 h on 16 epyc cores. This is an estimate, not a measurement: the first submission is the pilot with `--time=4:00:00`, and the run section sizes any further submission from its `wall_seconds` lines. Memory: 48 GB for 16 workers that each hold the universe (about 0.7 GB, estimate) and 11's weight matrices (2,001 x 7,150 float64, 114 MB per test set).

Step 11 runs one process with `OMP_NUM_THREADS=1` (review M-2). Its pilot `wall_seconds` is therefore a one-core number, not a 16-core number; 15 of the 16 requested cores are idle during step 11. R3 sizes from the measured value as it is.

The Phase B jobs (`jobs/j0_pilot.sh`, `jobs/j1_features.sh`, `jobs/j2_embed.sh`) run on exfab GPU nodes and use no MMseqs2. `test_avx2_tools_have_a_cpu_constraint` must not flag them, and this plan does not change them.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_jobs.py`:

```python
"""Static checks of the Phase C shell scripts and a C1 smoke run with the stub MMseqs2."""

import json
import os
import re
import subprocess

import paths
import pytest

PHASEC = paths.STEP1_DIR / "phasec"
SCRIPTS = sorted(PHASEC.glob("*.sh"))
# node features of AVX2 CPUs that sbatch accepts on UCR HPCC (epyc nodes: ryzen, amd, milan)
AVX2_FEATURES = ("ryzen", "milan", "genoa", "rome")
NON_AVX2_MMSEQS = "/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs"


def needs_avx2(text: str) -> bool:
    """The script loads the MMseqs2 module, names the AVX2 build, or is a SLURM script that
    runs mmseqs."""
    loads = re.search(r"module\s+load\b[^\n]*MMseqs2", text) or "17-b804f-avx2" in text
    return bool(loads) or ("#SBATCH" in text and "mmseqs" in text.lower())


def has_avx2_constraint(text: str) -> bool:
    found = re.findall(r"^#SBATCH\s+--constraint=(\S+)", text, flags=re.M)
    return any(f in AVX2_FEATURES for c in found for f in re.split(r"[&|,]", c))


def uses_non_avx2_binary(text: str) -> bool:
    return re.search(rf'STEP1_MMSEQS="?{re.escape(NON_AVX2_MMSEQS)}"?', text) is not None


def test_two_shell_scripts_exist():
    assert [s.name for s in SCRIPTS] == ["09_cluster_and_split.sh", "c1_evaluate.sh"]


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_shell_script_conventions(script):
    text = script.read_text()
    assert text.startswith("#!/bin/bash -l\n")
    assert "set -euo pipefail" in text
    assert "BASH_SOURCE" not in text
    assert '"${SCRATCH:?' in text
    assert ': "${PROJ_ROOT:?' in text and ': "${STEP1_WORKDIR:?' in text
    subprocess.run(["bash", "-n", str(script)], check=True)


def test_c1_requests_cpu_on_an_avx2_partition_with_a_time_limit():
    text = (PHASEC / "c1_evaluate.sh").read_text()
    assert "#SBATCH -p epyc" in text and "#SBATCH --time=" in text
    assert "#SBATCH --constraint=ryzen" in text  # ruling C-15
    assert """trap 'cp "$TMP/wall.txt" "$WALL_OUT"' EXIT""" in text  # review M-1
    assert "--gres" not in text  # no GPU (spec 5)
    assert "OMP_NUM_THREADS=1" in text
    assert 'STEPS="${C1_STEPS:-09 10 11}"' in text
    assert "wall_seconds=" in text


def test_avx2_tools_have_a_cpu_constraint():
    # owner rule: a job that runs an AVX2 tool must request a node feature that has AVX2
    scripts = sorted(paths.STEP1_DIR.rglob("*.sh"))
    flagged = {s.relative_to(paths.STEP1_DIR).as_posix(): s.read_text() for s in scripts
               if needs_avx2(s.read_text())}  # fmt: skip
    assert {"phasec/c1_evaluate.sh", "phasec/09_cluster_and_split.sh"} <= set(flagged)
    # the Phase B jobs run on exfab GPU nodes and use no MMseqs2
    assert not set(flagged) & {"jobs/j0_pilot.sh", "jobs/j1_features.sh", "jobs/j2_embed.sh"}
    bad = [
        n for n, t in flagged.items() if not has_avx2_constraint(t) and not uses_non_avx2_binary(t)
    ]
    assert bad == [], f"no #SBATCH --constraint with an AVX2 feature: {bad}"


def test_avx2_check_flags_a_planted_script():
    head = "#!/bin/bash -l\n#SBATCH -p epyc\n"
    plain = head + 'module load "${M:-MMseqs2/17-b804f}"\nmmseqs easy-cluster in out tmp\n'
    assert needs_avx2(plain) and not has_avx2_constraint(plain)
    assert has_avx2_constraint(plain.replace(head, head + "#SBATCH --constraint=ryzen\n"))
    assert not has_avx2_constraint(plain.replace(head, head + "#SBATCH --constraint=intel\n"))
    assert uses_non_avx2_binary(f"STEP1_MMSEQS={NON_AVX2_MMSEQS}\n")
    assert not needs_avx2("#!/bin/bash -l\n#SBATCH -p exfab\nsignalp6 --help\n")


def _env(tmp_path, fx, **extra):
    import phasec_fixture as pf

    (tmp_path / "scratch").mkdir(exist_ok=True)
    return {
        **os.environ,
        "SCRATCH": str(tmp_path / "scratch"),
        "PROJ_ROOT": str(paths.STEP1_DIR.parents[1]),
        "STEP1_WORKDIR": str(fx["work"]),
        "STEP1_MMSEQS": str(pf.STUB_MMSEQS),
        "STEP1_ENV_PY": os.environ.get("STEP1_ENV_PY", __import__("sys").executable),
        "SLURM_CPUS_PER_TASK": "2",
        "PHASEC_SPECIES": str(fx["species"]),
        "PHASEC_SETS": str(fx["sets"]),
        "PHASEC_N_RESAMPLES": "20",
        "PHASEC_CANDIDATES": "B1,R2,M8",
        **extra,
    }


def test_c1_smoke_run_with_the_stub(tmp_path):
    pytest.importorskip("sklearn")
    import phasec_fixture as pf
    from conftest import load_phasec

    fx = pf.make_work(tmp_path)
    assert load_phasec("08_build_eval_tables").main(pf.build_argv(fx)) == 0
    run = subprocess.run(["bash", str(PHASEC / "c1_evaluate.sh")], env=_env(tmp_path, fx),
                         capture_output=True, text=True)  # fmt: skip
    assert run.returncode == 0, run.stdout[-2000:] + run.stderr[-2000:]
    for step in ("09", "10", "11"):
        assert f"C1 step {step} wall_seconds=" in run.stdout
    out = fx["work"] / "phasec"
    m = json.loads((out / "metrics.json").read_text())
    assert m["settings"]["n_resamples"] == 20 and m["settings"]["candidates"] == ["B1", "R2", "M8"]
    assert list((out / "logs").glob("wall.*.txt"))


def test_c1_stops_when_mmseqs_fails(tmp_path):
    pytest.importorskip("sklearn")
    import phasec_fixture as pf
    from conftest import load_phasec

    fx = pf.make_work(tmp_path)
    assert load_phasec("08_build_eval_tables").main(pf.build_argv(fx)) == 0
    env = _env(tmp_path, fx, STUB_MMSEQS_FAIL="easy-cluster")
    run = subprocess.run(["bash", str(PHASEC / "c1_evaluate.sh")], env=env,
                         capture_output=True, text=True)  # fmt: skip
    assert run.returncode == 2 and "STOP: mmseqs easy-cluster exited with 1" in run.stderr
    assert "C1 step 09 wall_seconds" not in run.stdout
    assert not (fx["work"] / "phasec" / "split_members.tsv.gz").exists()


def test_c1_keeps_the_wall_times_after_a_failed_step(tmp_path):
    # review M-1: step 09 passes, step 10 stops; wall.<job>.txt must still hold the 09 line
    pytest.importorskip("sklearn")
    import phasec_fixture as pf
    from conftest import load_phasec

    fx = pf.make_work(tmp_path)
    assert load_phasec("08_build_eval_tables").main(pf.build_argv(fx)) == 0
    env = _env(tmp_path, fx, SLURM_JOB_ID="m1test", PHASEC_CANDIDATES="M8,XGB")
    run = subprocess.run(["bash", str(PHASEC / "c1_evaluate.sh")], env=env,
                         capture_output=True, text=True)  # fmt: skip
    assert run.returncode == 2 and "unknown candidates ['XGB']" in run.stderr
    assert "C1 step 09 wall_seconds=" in run.stdout and "C1 done" not in run.stdout
    wall = (fx["work"] / "phasec" / "logs" / "wall.m1test.txt").read_text()
    assert wall.startswith("C1 step 09 wall_seconds=") and "step 10" not in wall
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_jobs.py -q`
Expected: `6 failed, 1 passed, 1 skipped`: `test_two_shell_scripts_exist` (`[] == ['09_cluster_and_split.sh', 'c1_evaluate.sh']`), `test_c1_requests_cpu_on_an_avx2_partition_with_a_time_limit` (no file), `test_avx2_tools_have_a_cpu_constraint` (the two scripts are not flagged because they do not exist) and the three C1 runs (bash exit 127); `test_shell_script_conventions` is skipped (empty parameter set). `test_avx2_check_flags_a_planted_script` passes: it tests the checker on planted text, not the scripts.

- [ ] **Step 3: Write `analysis/step1_compare/phasec/09_cluster_and_split.sh`**

```bash
#!/bin/bash -l
# Phase C step 2 (spec 5): MMseqs2 clustering, maximum-identity search and the split tables.
# c1_evaluate.sh runs this script inside its SLURM job; it is not submitted alone.
# PROJ_ROOT and STEP1_WORKDIR come from the environment, never from the script location.
# MMseqs2 temp files go to node-local $SCRATCH; 09_make_splits.py writes its outputs to
# $STEP1_WORKDIR/phasec/ (on /bigdata) with a temp name and os.replace.
# The module MMseqs2/17-b804f puts the AVX2 build on PATH. It stops with "Illegal instruction"
# (exit 132) on CPUs without AVX2 (for example the abu_dhabi nodes of partition batch), so
# c1_evaluate.sh requests partition epyc with --constraint=ryzen. STEP1_MMSEQS (a binary path)
# skips the module (tests pass the stub there). A job that runs an AVX2 tool must request a node
# feature that has AVX2; the two #SBATCH lines below keep this script on such nodes if it is
# ever submitted alone (c1_evaluate.sh requests the same).
#SBATCH -p epyc
#SBATCH --constraint=ryzen
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root)}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR}"
TMP="${SCRATCH:?SCRATCH is not set; run this inside a SLURM job}/phasec_09"
ENV_PY="${STEP1_ENV_PY:-/rhome/jstajich/.conda/envs/adhesionPred/bin/python}"
S1="$PROJ_ROOT/analysis/step1_compare"
export PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/phasec"
if [ -z "${STEP1_MMSEQS:-}" ]; then
  module load "${PHASEC_MMSEQS_MODULE:-MMseqs2/17-b804f}"
  STEP1_MMSEQS=$(command -v mmseqs)
fi
"$ENV_PY" "$S1/phasec/09_make_splits.py" --work-dir "$STEP1_WORKDIR" \
  --species "${PHASEC_SPECIES:-$S1/species.tsv}" --mmseqs "$STEP1_MMSEQS" \
  --tmp-dir "$TMP" --threads "${SLURM_CPUS_PER_TASK:-4}"
```

- [ ] **Step 4: Write `analysis/step1_compare/phasec/c1_evaluate.sh`**

```bash
#!/bin/bash -l
# Phase C job C1 (spec 5): 09 (clusters and splits), 10 (fit and score), 11 (metrics) on CPU.
# Submit (see the plan's run section):
#   sbatch --export=ALL,PROJ_ROOT=...,STEP1_WORKDIR=... -o <log> -e <log> c1_evaluate.sh
# One job, not one job per split: the three steps depend on each other, and the global rule
# says not to split work into jobs of minutes. No run time is measured on the real data yet;
# the first submission is the pilot. It prints `C1 step <n> wall_seconds=<s>` for every step,
# and the run section sizes later submissions from these lines (C1_STEPS selects steps).
# The time limit of 4 h is a bound for the pilot, not a measurement (estimate in the plan).
# A job that runs an AVX2 tool must request a node feature that has AVX2 (--constraint=ryzen on
# partition epyc). Partition epyc and constraint ryzen (ruling C-15): the epyc nodes carry the
# features ryzen, amd, milan; their CPUs have AVX2, which the MMseqs2 module build needs (it
# stops with exit 132 on the abu_dhabi Opterons of partition batch). Step 11 runs one process
# with one BLAS thread, so its wall time is a one-core number, not a 16-core number.
# PROJ_ROOT comes from the environment, never from the script location. Temp files go to
# node-local $SCRATCH; the scripts write their outputs to $STEP1_WORKDIR/phasec/ (on /bigdata)
# with temp names and os.replace, so the results are on /bigdata when the job ends.
#SBATCH -p epyc
#SBATCH --constraint=ryzen
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH -J step1_c1
set -euo pipefail

: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${STEP1_WORKDIR:?export STEP1_WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/step1_c1"
ENV_PY="${STEP1_ENV_PY:-/rhome/jstajich/.conda/envs/adhesionPred/bin/python}"
S1="$PROJ_ROOT/analysis/step1_compare"
export PYTHONPATH="$PROJ_ROOT/src:$S1:$S1/phasec"
# one BLAS thread per process: 10 runs one process per CPU
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CPUS="${SLURM_CPUS_PER_TASK:-4}"
STEPS="${C1_STEPS:-09 10 11}"
LOGDIR="$STEP1_WORKDIR/phasec/logs"
mkdir -p "$TMP" "$LOGDIR"
: > "$TMP/wall.txt"
WALL_OUT="$LOGDIR/wall.${SLURM_JOB_ID:-local}.txt"
# copy the wall times on every exit path, also after a failed step (review M-1)
trap 'cp "$TMP/wall.txt" "$WALL_OUT"' EXIT

timed() {  # timed NAME CMD...: run CMD and record its wall time
  local name=$1 t0=$SECONDS
  shift
  "$@"
  echo "C1 step $name wall_seconds=$((SECONDS - t0))" | tee -a "$TMP/wall.txt"
}

for step in $STEPS; do
  case "$step" in
    09) timed 09 bash "$S1/phasec/09_cluster_and_split.sh" ;;
    10) timed 10 "$ENV_PY" "$S1/phasec/10_fit_and_score.py" --work-dir "$STEP1_WORKDIR" \
          --workers "$CPUS" --candidates "${PHASEC_CANDIDATES:-B0,B1,R0,R1,R2,M8,M35,M8-C,M35-C,H}" ;;
    11) timed 11 "$ENV_PY" "$S1/phasec/11_evaluate.py" --work-dir "$STEP1_WORKDIR" \
          --sets "${PHASEC_SETS:-$S1/sequence_sets.tsv}" \
          --species "${PHASEC_SPECIES:-$S1/species.tsv}" \
          --n-resamples "${PHASEC_N_RESAMPLES:-2000}" ;;
    *) echo "C1: unknown step $step" >&2; exit 2 ;;
  esac
done
echo "C1 done: $STEPS"
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_jobs.py tests/step1_compare/test_paths.py -q`
Expected: `18 passed` (9 + 9; `test_no_shell_script_uses_bash_source` in `test_paths.py` scans `phasec/*.sh` too).

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/phasec/09_cluster_and_split.sh analysis/step1_compare/phasec/c1_evaluate.sh \
  tests/step1_compare/test_phasec_jobs.py
git commit -m "phasec: SLURM job C1 (09, 10, 11 on epyc) and the MMseqs2 wrapper

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 12: Documentation of the Phase C outputs

**Files:**
- Modify: `analysis/step1_compare/COLUMNS.md` (append), `analysis/step1_compare/README.md` (append)
- Test: `tests/step1_compare/test_phasec_docs.py`

**Interfaces:**
- Consumes: `OUTPUT_NAMES` and column tuples of 08 to 12, `dedupe`, `splits` (Tasks 2 to 10); `section`, `table_names`, `key_names`, `COLUMNS`, `README` of `test_phaseb_docs.py`.
- Produces: one `## phasec/<name>` section in `COLUMNS.md` per Phase C output; a README section "Phase C".

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_phasec_docs.py`:

```python
"""COLUMNS.md and README.md document the Phase C outputs exactly as the code writes them.

Same rules as test_phaseb_docs.py: the set of names in a table section must equal the code
tuple; the bullet list after a `**Keys...:**` marker must equal the keys of a real output.
"""

import json

import pytest

pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import dedupe  # noqa: E402
import paths  # noqa: E402
import splits  # noqa: E402
from conftest import load_phasec  # noqa: E402
from test_phaseb_docs import COLUMNS, README, key_names, section, table_names  # noqa: E402

SCRIPTS = ("08_build_eval_tables", "09_make_splits", "10_fit_and_score", "11_evaluate", "12_report")


def _json(work, name):
    return json.loads((work / "phasec" / name).read_text())


def test_phasec_outputs_have_columns_md_headings():
    names = [f"phasec/{n}" for s in SCRIPTS for n in load_phasec(s).OUTPUT_NAMES]
    names.append("phasec/logs/wall.<job>.txt")
    headings = [line[3:].split()[0] for line in COLUMNS.splitlines() if line.startswith("## ")]
    missing = [n for n in names if n not in headings]
    assert not missing, f"COLUMNS.md has no heading for {missing}"


def test_table_columns_equal_the_code_constants():
    s08 = load_phasec("08_build_eval_tables")
    s10 = load_phasec("10_fit_and_score")
    expected = {
        "phasec/eval_table.tsv.gz": dedupe.TABLE_COLUMNS,
        "phasec/eval_literature.tsv": s08.LITERATURE_COLUMNS,
        "phasec/eval_dedupe_log.tsv": dedupe.LOG_COLUMNS,
        "phasec/clusters.tsv.gz": splits.CLUSTER_COLUMNS,
        "phasec/split_members.tsv.gz": splits.MEMBER_COLUMNS,
        "phasec/tc_removed.tsv": splits.REMOVED_COLUMNS,
        "phasec/max_identity.tsv.gz": splits.IDENTITY_COLUMNS,
        "phasec/scores.tsv.gz": s10.SCORE_COLUMNS,
    }
    for heading, columns in expected.items():
        assert table_names(heading) == set(columns), heading


def test_proteome_calls_columns_equal_the_code(phasec_chain):
    s11 = load_phasec("11_evaluate")
    patterns = {"call_<candidate>_<variant>", "score_<candidate>_<variant>"}
    got = table_names("phasec/proteome_calls.tsv.gz")
    assert got == set(s11.PROTEOME_BASE) | patterns
    m = _json(phasec_chain["work"], "metrics.json")
    cols = s11.proteome_columns(tuple(m["settings"]["candidates"]))
    header = (
        __import__("gzip")
        .decompress((phasec_chain["work"] / "phasec" / "proteome_calls.tsv.gz").read_bytes())
        .decode()
        .splitlines()[0]
        .split("\t")
    )
    assert header == list(cols)
    dynamic = set(cols) - set(s11.PROTEOME_BASE)
    assert {c.split("_")[0] + "_<candidate>_<variant>" for c in dynamic} == patterns


def test_run_json_keys_equal_the_code(phasec_chain, tmp_path):
    import phasec_fixture as pf

    fx = pf.copy_work(phasec_chain, tmp_path)
    assert load_phasec("12_report").main(["--work-dir", str(fx["work"])]) == 0
    w = fx["work"]
    build = _json(w, "build_run.json")
    assert key_names("phasec/build_run.json") == set(build)
    assert key_names("phasec/build_run.json", "**Keys of `input_sha256`:**") == set(
        build["input_sha256"]
    )
    assert key_names("phasec/splits_run.json") == set(_json(w, "splits_run.json"))
    scores = _json(w, "scores_run.json")
    assert key_names("phasec/scores_run.json") == set(scores)
    assert key_names("phasec/scores_run.json", "**Keys of `input_sha256`:**") == set(
        scores["input_sha256"]
    )
    unit = scores["units"]["S1|0|V-go"]
    assert key_names("phasec/scores_run.json", "**Keys of a rule object in `units`:**") == set(
        unit["R2"]
    )
    lr = "**Keys of a logistic-regression object in `units`:**"
    assert key_names("phasec/scores_run.json", lr) == set(unit["M8"]) == set(unit["H"])
    assert key_names("phasec/evaluate_run.json") == set(_json(w, "evaluate_run.json"))
    assert key_names("phasec/report_run.json") == set(_json(w, "report_run.json"))


def test_metrics_and_findings_keys_equal_the_code(phasec_chain):
    m = _json(phasec_chain["work"], "metrics.json")
    f = _json(phasec_chain["work"], "findings.json")
    assert key_names("phasec/metrics.json") == set(m)
    assert key_names("phasec/findings.json") == set(f)
    s2 = m["test_sets"]["S2-Spom_PomBase:Spom_PomBase"]
    s1 = m["test_sets"]["S1:all"]
    assert key_names("phasec/metrics.json", "**Keys of a test set object:**") == set(s2)
    assert set(s1) == set(s2) - {"calibration"}  # calibration only for S2 test sets
    assert key_names("phasec/metrics.json", "**Keys of a truth object:**") == set(
        s1["truth"]["direct"]
    )
    assert f"`{m['schema']}`" in section("phasec/metrics.json")
    assert f"`{f['schema']}`" in section("phasec/findings.json")


def test_readme_names_the_phase_c_commands_and_rules():
    for text in ("phasec/08_build_eval_tables.py", '"$S1/phasec/c1_evaluate.sh"', "phasec/12_report.py",
                 "C1_STEPS", "STEP1_MMSEQS", "Illegal", "tc_taxon_clades.tsv", "wall_seconds"):  # fmt: skip
        assert text in README, text
    job = (paths.STEP1_DIR / "phasec" / "c1_evaluate.sh").read_text()
    assert "C1_STEPS" in job and "wall_seconds" in job
    assert '`input_sha256["unique_sequences.tsv.gz"]` must equal `embedding_run.json`' in README
    rule = "A job that runs an AVX2 tool must request a node feature that has AVX2"
    assert rule in " ".join(README.split()) and rule in " ".join(job.split())


def test_grids_in_columns_md_equal_the_code():
    import findings
    import models
    import rules

    def bullet(heading, key):
        lines = section(heading).splitlines()
        return next(x for x in lines if x.startswith(f"- `{key}` "))

    c_line = bullet("phasec/scores_run.json", "c_grid")
    values = c_line.split(":", 1)[1].split("(")[0].split(",")
    assert tuple(float(x) for x in values) == models.C_GRID  # ruling C-12
    t_line = bullet("phasec/scores_run.json", "t")
    lo, hi = (float(x) for x in t_line.split("cut, ")[1].split(" in steps")[0].split(" to "))
    assert (lo, hi) == (rules.T_VALUES[0], rules.T_VALUES[-1])
    floor = bullet("phasec/metrics.json", "floor_met")
    assert f"at least {findings.MIN_DIRECT_POSITIVES} " in floor
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare/test_phasec_docs.py -q`
Expected: `7 failed`, for example `COLUMNS.md has no heading for ['phasec/eval_table.tsv.gz', ...]` and `'phasec/08_build_eval_tables.py' in README`.

- [ ] **Step 3: Append the Phase C section to `analysis/step1_compare/COLUMNS.md`**

Append this text at the end of the file:

````markdown

# Phase C outputs (evaluation)

Scripts 08 to 12 in `phasec/` write to `$STEP1_WORKDIR/phasec/`. Each run JSON records
`outputs_sha256` for its other outputs; the next script stops when a file differs from it.
Tables join on `seq_sha256`. List columns hold sorted unique values, comma separated.

## phasec/eval_table.tsv.gz (08_build_eval_tables.py, one row per unique sequence)

| Column | Meaning |
|---|---|
| seq_sha256 | Hash of the sequence (Phase A `seqhash.seq_sha256`). |
| origin | `go` (truth genes) or `tc` (keyword tier T-c rows that survive ruling C-6). |
| class | `pos`, `neg` or `excluded` (`labelmap.class_of`; `excluded` also for a hash with a class and an excluded gene). |
| label | Truth labels of the members (list); `T-c` for T-c rows. |
| subset | `wall`, `extracellular-only` (list) for P-ext members. |
| stratum | Reporting strata of the members (list): `wall`, `extracellular-only`, `PM-TM`, `pm-unresolved`, `N-int`, `N-sec`, `ambiguous`; `T-c` for T-c rows. |
| d8_class | D8 classes of the members (list). |
| homology_only | `no` when any member has the label without homology codes (direct evidence); `yes` otherwise; empty for T-c rows. |
| internal_evidence_htp_only | `yes` when any member says so; empty for T-c rows. |
| source_ids | Truth sources of the members (list); `T-c` for T-c rows. |
| gene_ids | Gene IDs of the members, or the UniProt accessions of the T-c rows (list). |
| species | Species (truth) or genome (T-c) names (list). |
| roles | Source roles from `species.tsv` (list); `tc` for T-c rows. |
| clades | `in_clade` of the sources, or the clade from `phasec/tc_taxon_clades.tsv` (list). |
| taxon_ids | NCBI taxon IDs (list). |
| length | Sequence length. |
| emb_row, emb_cterm_row | As in `features.tsv.gz`. |

## phasec/eval_literature.tsv (08_build_eval_tables.py, one row per literature seed with an accession)

| Column | Meaning |
|---|---|
| accession | UniProt accession from `uniprot_query`. |
| gene, species, order, moonlighting | From `eurotiomycetes_seeds.tsv`. |
| lit_class | The seed file's `class` (`adhesin`, `hard_negative`). It describes adhesion, not location. |
| seq_sha256, length, emb_row, emb_cterm_row | From the `uniprot_kw` member; empty without a sequence. |
| literature_positive | `yes` when the row has a sequence and `moonlighting` is not `YES` (spec 3.2). |

## phasec/eval_dedupe_log.tsv (08_build_eval_tables.py, one row per dropped or reclassified member)

| Column | Meaning |
|---|---|
| origin | `go` or `tc`. |
| source_id, gene_id | Truth source and gene, or `T-c` and the accession. |
| seq_sha256 | Hash of the member. |
| class | Class of the member. |
| reason | `both_classes` (hash dropped), `class_and_excluded` (row becomes excluded), `go_label_wins` (T-c row dropped, ruling C-6). |
| detail | The classes of the GO members with this hash. |

## phasec/eval_sequences.fasta.gz (08_build_eval_tables.py)

One record per hash of `eval_table.tsv.gz` and `eval_literature.tsv`, sorted by hash. The header
is the `seq_sha256`. The MMseqs2 input of 09.

## phasec/build_run.json (08_build_eval_tables.py)

**Keys:**

- `all_sources` `true` (08 stops otherwise)
- `truth_set_sha256` shared by `d8_run.json`, `features_run.json`, `keyword_tier_run.json`
- `input_sha256` see below
- `go_members_by_source_class` GO members per source and class, before dedupe
- `labelled_genes_without_sequence` per source
- `table_rows_by_origin_class` key `<origin>:<class>`
- `log_rows_by_reason`
- `tc_rows` rows of `keyword_tier.tsv.gz`
- `tc_rows_dropped_by_go_class` T-c rows dropped by ruling C-6, by the classes of the GO members
- `go_members_sharing_a_tc_hash` key `<class>:<label>:<d8_class>`
- `literature_rows`
- `literature_positives`
- `literature_no_accession` seed genes without an accession
- `sequences` records in `eval_sequences.fasta.gz`
- `git_commit`
- `library_versions` Python, numpy, scipy, scikit-learn
- `arguments`
- `outputs_sha256` the other four outputs

**Keys of `input_sha256`:**

- `truth_set_triaged.tsv.gz`
- `features.tsv.gz`
- `features_unique.tsv.gz`
- `unique_sequences.tsv.gz`
- `keyword_tier.tsv.gz`
- `species.tsv`
- `eurotiomycetes_seeds.tsv`
- `tc_taxon_clades.tsv`

## phasec/clusters.tsv.gz (09_make_splits.py, one row per sequence)

| Column | Meaning |
|---|---|
| seq_sha256 | Hash of a sequence of `eval_sequences.fasta.gz`. |
| cluster_id | Hash of the MMseqs2 cluster representative (`easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0`). |

## phasec/split_members.tsv.gz (09_make_splits.py, one row per split, fold and sequence)

| Column | Meaning |
|---|---|
| split_id | `S1`, `S2-<source>`, `S3-<clade>` or `FULL` (`splits.py`). |
| fold | S1 fold 0 to 4; `0` for the other splits. |
| seq_sha256 | Hash of the sequence. |
| part | `train`, `train_tc` (V-kw only), `test`, `test_tc` (S1, scored, never truth), `test_lit`. |
| origin | `go`, `tc` or `lit`. |
| class | `pos`, `neg` or `excluded`. |
| cluster_id | As in `clusters.tsv.gz`. |

## phasec/tc_removed.tsv (09_make_splits.py, one row per T-c row removed from a split)

| Column | Meaning |
|---|---|
| split_id, fold | As in `split_members.tsv.gz`. |
| seq_sha256 | Hash of the T-c row. |
| gene_ids, taxon_ids | As in `eval_table.tsv.gz`. |
| rule | `a_test_protein`, `b_cluster_mate` (ruling C-5), `c_test_taxon` (ruling C-7). |

## phasec/max_identity.tsv.gz (09_make_splits.py, one row per S2 or S3 split and test sequence)

| Column | Meaning |
|---|---|
| split_id | An S2 or S3 split. |
| seq_sha256 | A test or literature sequence. |
| max_identity | Highest `fident` against the GO training proteins of the split (`easy-search -s 7.5 -c 0.5 --cov-mode 0`); empty when there is no hit. Self hits are skipped. |
| below_0.3 | `yes` when there is no hit or the identity is below 0.3 (ruling C-4). |

## phasec/splits_run.json (09_make_splits.py)

**Keys:**

- `all_sources`
- `truth_set_sha256`
- `input_sha256` the three 08 inputs and `species.tsv`
- `mmseqs_version` output of `mmseqs version`
- `mmseqs_commands` the two MMseqs2 commands without paths and threads
- `cluster_tsv_sha256` the raw MMseqs2 cluster table
- `seed` 20261001
- `sequences`
- `clusters`
- `splits` split ids in order
- `members_by_split_fold_part` key `<split>|<fold>|<part>`
- `tc_removed_by_split_fold_rule` key `<split>|<fold>|<rule>`
- `identity_below_0.3_by_split`
- `git_commit`
- `library_versions`
- `arguments`
- `outputs_sha256`

## phasec/scores.tsv.gz (10_fit_and_score.py, one row per unit, candidate and scored sequence)

| Column | Meaning |
|---|---|
| split_id, fold | The outer training set (unit). |
| variant | `V-go` or `V-kw`. |
| candidate | `B0`, `B1`, `R0`, `R1`, `R2`, `M8`, `M35`, `M8-C`, `M35-C`, `H`. |
| seq_sha256 | The scored sequence. |
| part | Its part in the split; `all` for FULL (every Phase B unique sequence). |
| score | Logistic-regression decision value; empty for rules. |
| prob | Platt probability (ruling C-10); empty for rules. |
| call | `1` when score >= the unit's threshold (Youden's J), or the rule call; else `0`. |

## phasec/scores_run.json (10_fit_and_score.py)

**Keys:**

- `all_sources`
- `truth_set_sha256`
- `input_sha256` see below
- `embedding_array_sha256` key `<model>.<window>`
- `seed`
- `candidates`
- `c_grid` values of C tried: 0.001, 0.003, 0.01, 0.1, 1, 10 (ruling C-12)
- `inner_folds`
- `units` key `<split>|<fold>|<variant>`, one object per candidate
- `git_commit`
- `library_versions`
- `arguments`
- `outputs_sha256`

**Keys of `input_sha256`:**

- `split_members.tsv.gz`
- `clusters.tsv.gz`
- `features_unique.tsv.gz`
- `unique_sequences.tsv.gz`

**Keys of a rule object in `units`:**

- `n_train`
- `n_train_pos`
- `g` PredGPI class cut (null for R0)
- `t` Ser+Thr cut, 0.20 to 0.40 in steps of 0.05 (null for R0 and R1)
- `j` Youden's J on the training rows

**Keys of a logistic-regression object in `units`:**

- `n_train`
- `n_train_pos`
- `C`
- `h_variant` the ESM variant of H; null for the other candidates
- `threshold` decision value cut (Youden's J on the inner out-of-fold values)
- `platt_a`
- `platt_b`
- `inner_pr_auc` inner out-of-fold PR-AUC per setting
- `convergence_warnings`

## phasec/metrics.json (11_evaluate.py)

Values are objects `{value, lo, hi, n_defined}`: the point value, the 2.5th and 97.5th
percentiles over the resamples where the metric is defined, and their number. `null` means
not defined (for example recall without positives). A truth block without negatives (the
literature set) has `null` precision, PR-AUC, precision at recall and precision at the rule's
recall (spec 3.2: recall only).

**Keys:**

- `schema` `step1-phasec-metrics/1`
- `settings` resamples, seed, cut-offs, candidates, variants
- `test_sets` key = test set name, see below
- `agreement` R2 against each ML candidate per proteome set and per truth class
- `score_sources` per proteome set: counts of `oof`, `in_sample`, `final`
- `named_panel` one object per panel protein (parent spec section 6)
- `context` T-c rows of C. immitis and C. posadasii in V-kw training, per split and fold
- `tc_removed` as `tc_removed_by_split_fold_rule` in `splits_run.json`

**Keys of a test set object:**

- `split`
- `kind` `pooled`, `species`, `clade` or `literature`
- `truth` key `direct` or `all`, see below
- `recall_half_width` per label candidate (R2 and ML; V-go; direct truth)
- `fpr_half_width`
- `n_direct_positives` positives of the direct truth (stratum `all`)
- `floor_met` `true` when `n_direct_positives` is at least 20 (ruling C-8)
- `label` `estimate` (every recall half-width at most 0.10 and `floor_met`) or `smoke test` (ruling C-8)
- `max_recall_half_width`
- `max_fpr_half_width`
- `zero_width_recall_interval` label candidates whose recall interval has lo equal to hi
- `ambiguous` score distribution of the ambiguous genes
- `ambiguous_htp_only` the same for the high-throughput-only sub-stratum
- `lists` `pm-unresolved` and `P-gpi` genes with their calls
- `calibration` S2 test sets only: Brier score and 10 reliability bins per ML candidate

**Keys of a truth object:**

- `n` positives and negatives per stratum
- `metrics` per stratum, variant and candidate
- `variant_effect` V-kw minus V-go per stratum, candidate and metric (paired)
- `vs_rule` precision at the recall of R2 and recall at the FPR of R2 per ML candidate
- `prevalence` precision at assumed prevalence 0.01, 0.03, 0.05, 0.1 (assumed, not measured)
- `nsec_fpr_at_rule_recall` N-sec FPR at the recall of R2 and the paired differences

## phasec/findings.json (11_evaluate.py)

**Keys:**

- `schema` `step1-phasec-findings/1`
- `a_b1_not_saturated` B1 ROC-AUC below 0.99 on S1:all and every S2 test set (direct, V-go)
- `b_ml_beats_b1_and_r2_on_nsec_s1` per ML candidate: B1 and R2 N-sec FPR minus the ML N-sec FPR at the recall of R2, interval above 0
- `c_same_under_s2` the statement of (b) on every S2 test set

## phasec/proteome_calls.tsv.gz (11_evaluate.py, one row per protein of every proteome set)

| Column | Meaning |
|---|---|
| set_id, source_id, gene_id, seq_sha256 | As in `sequence_members.tsv.gz`; sets with kind `download` or `site` in `sequence_sets.tsv`. |
| class, label, origin | From `eval_table.tsv.gz` when the hash is there; else empty. |
| score_source | `oof` (hash in an S1 fold), `in_sample` (in a FULL training table but in no S1 fold), `final` (in no training table). |
| call_<candidate>_<variant> | `1` or `0`. |
| score_<candidate>_<variant> | Decision value; empty for rules. |

## phasec/evaluate_run.json (11_evaluate.py)

**Keys:**

- `all_sources`
- `truth_set_sha256`
- `input_sha256` scores, split members, table, sequence members, sequence sets
- `seed`
- `n_resamples`
- `git_commit`
- `library_versions`
- `arguments`
- `outputs_sha256`

## phasec/report.md (12_report.py)

The owner's report. Every number in it is a value of `metrics.json` or `findings.json`;
12 checks this before it writes.

## phasec/report_run.json (12_report.py)

**Keys:**

- `all_sources`
- `truth_set_sha256`
- `input_sha256` `metrics.json` and `findings.json`
- `git_commit`
- `library_versions`
- `arguments`
- `outputs_sha256`

## phasec/logs/wall.<job>.txt (phasec/c1_evaluate.sh)

One line `C1 step <step> wall_seconds=<seconds>` per finished step of job C1.
````

- [ ] **Step 4: Append the Phase C section to `analysis/step1_compare/README.md`**

Append this text at the end of the file:

````markdown

## Phase C: rule versus ML evaluation (E1 to E4)

The spec is `docs/superpowers/specs/2026-10-01-step1-phase-c-evaluation-design.md`; the plan is
`docs/superpowers/plans/2026-10-01-step1-phase-c-evaluation.md`. `COLUMNS.md` lists every
Phase C output. All Phase C outputs go to `$STEP1_WORKDIR/phasec/`.

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
- C1 prints `C1 step <step> wall_seconds=<seconds>` per step and copies these lines to
  `phasec/logs/` on every exit, also after a failed step. `C1_STEPS` (default `09 10 11`) runs
  a subset of the steps. Step 11 runs one process with one BLAS thread.
- Phase C tests need numpy and scikit-learn:
  `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare -q`. A Python without these libraries
  skips those test files.
````

- [ ] **Step 5: Run the whole suite**

Run: `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare -q`
Expected: `533 passed, 2 skipped` (prototype after the plan-review fixes, Python 3.14.2, 663 s; 409 old + 124 new passed; skipped: one old test and the real-MMseqs2 test)
Run: `$PY -m pytest tests/step1_compare -q`
Expected: `484 passed, 1 failed, 8 skipped` (prototype after the plan-review fixes, Python 3.12.14, 209 s; the failure is the old surface_glyco import). `~/.local` gives Python 3.12 numpy 2.5.2, scipy 1.18.1 and scikit-learn 1.9.1, so the Phase C tests also ran and passed there, the golden file included (tolerance 1e-9, other library versions); without these libraries they skip.

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add analysis/step1_compare/COLUMNS.md analysis/step1_compare/README.md \
  tests/step1_compare/test_phasec_docs.py
git commit -m "step1_compare: document the Phase C outputs and the run order

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

## Run (controller, after review)

Run these steps in this order. The plan author ran none of them on the real work directory and submitted no job. Numbers marked "prototype" come from the scratch run described in "Verification".

```bash
export PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred-eval   # checkout of the reviewed commit
export STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
ENV_PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python
S1=$PROJ_ROOT/analysis/step1_compare
export PYTHONPATH=$PROJ_ROOT/src:$S1:$S1/phasec
LOG=$STEP1_WORKDIR/logs
cd "$S1"
```

**R0. Preflight (login node).**

```bash
git -C "$PROJ_ROOT" status --short                     # expect no output
$ENV_PY -m pytest "$PROJ_ROOT/tests/step1_compare" -q
grep -h '"all_sources"' "$STEP1_WORKDIR"/{d8_run,keyword_tier_run}.json "$STEP1_WORKDIR"/phaseb/features_run.json
ls "$STEP1_WORKDIR/phasec" 2>/dev/null                 # expect: no such directory
```
Accept: the suite passes (counts in Task 12 Step 5); three lines `"all_sources": true`; no `phasec/` directory yet.

**R1. 08 (login node, about 20 s).**

```bash
$ENV_PY phasec/08_build_eval_tables.py
$ENV_PY -c 'import json, os; b = json.load(open(os.environ["STEP1_WORKDIR"] + "/phasec/build_run.json")); print("labelled_genes_without_sequence", b["labelled_genes_without_sequence"], "literature_positives", b["literature_positives"])'
```
Accept: exit 0 and `table_rows=23338 sequences=23351 literature=20`, `go:excluded=266 go:neg=19776 go:pos=752 tc:pos=2544`, `tc_dropped_go_class_excluded=65 tc_dropped_go_class_neg=224 tc_dropped_go_class_pos=208` (prototype values on the same inputs); `literature_positives 19` (the 9 `hard_negative` rows included, ruling C-9 amended). Another value means an input changed: read `phasec/build_run.json` before you go on. A STOP names the failed check. Record in the run notes `labelled_genes_without_sequence` (prototype: `Calb_CGD` 72, `Scer_SGD` 3, `Spom_PomBase` 2; the 75 of the two training sources are labelled negatives without a feature row, review M-5). They have no sequence, so they are in no table; this explains part of the difference between the spec's 6,951 negative members and the 6,842 negative training sequences (ruling C-14).

**R2. C1 pilot (epyc CPU; the first submission is the pilot).**

```bash
sbatch --export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR" \
  -o "$LOG/c1.%j.log" -e "$LOG/c1.%j.log" "$S1/phasec/c1_evaluate.sh"
```
Accept: the log ends with `C1 done: 09 10 11` and holds three lines `C1 step <step> wall_seconds=<s>`; `phasec/logs/wall.<job id>.txt` exists (it also exists after a failed step, with the lines of the finished steps). Then check the outputs:

```bash
$ENV_PY - <<'EOF'
import json, os, gzip
d = os.environ["STEP1_WORKDIR"] + "/phasec/"
s = json.load(open(d + "splits_run.json"))
print("sequences", s["sequences"], "clusters", s["clusters"], "mmseqs", s["mmseqs_version"])
print(json.dumps(s["tc_removed_by_split_fold_rule"], indent=1))
r = json.load(open(d + "scores_run.json"))
assert len(r["units"]) == 22, len(r["units"])
warn = {k: c for k, u in r["units"].items() for c, p in u.items() if p.get("convergence_warnings")}
print("convergence warnings", warn)
m = json.load(open(d + "metrics.json"))
assert m["settings"]["n_resamples"] == 2000 and len(m["test_sets"]) == 13, len(m["test_sets"])
for name, ts in m["test_sets"].items():
    assert ts["floor_met"] == (ts["n_direct_positives"] >= 20), name
    assert ts["label"] == "smoke test" or ts["floor_met"], name  # ruling C-8 floor
    print(name, ts["label"], ts["n_direct_positives"], ts["floor_met"],
          ts["max_recall_half_width"], ts["zero_width_recall_interval"])
print("score sources", m["score_sources"])
n = sum(1 for _ in gzip.open(d + "proteome_calls.tsv.gz", "rt")) - 1
print("proteome_calls rows", n)
EOF
```
Accept: `sequences 23351`; 22 units; 13 test sets; `S3-Eurotiomycetes:Afum_ASPFU`, `S3-Basidiomycota:Cneo_H99_GOA`, `S3-Basidiomycota:Umay_MYCMD` and `S3-Eurotiomycetes:literature` are `smoke test` (19, 7, 9 and 19 direct positives in the prototype counts, below the floor); `proteome_calls rows 69150` (members of the nine proteome sets in `sequence_sets.tsv`, prototype count); no `in_sample` in `score_sources`. Record in the run notes: `clusters` (prototype: 9,737 with the non-AVX2 binary; the AVX2 build of the job was not run, so another count is possible and is not a STOP), the removal counts (prototype: S2-Calb_CGD b 506, c 164; S2-Scer_SGD b 633, c 69; S2-Spom_PomBase b 464; S3-Basidiomycota b 293; S3-Eurotiomycetes b 674, c 1,235; S1 none), any convergence warnings, and every label.

**R3. Size from the pilot (global job-size rule).**

Read the three `wall_seconds` values. Step 11 is a one-process number (review M-2): it does not shrink with more cores. The steps depend on each other, so one job is the smallest job count; there is no unit of work to merge into it.
- Total at most 5,400 s (1.5 h): keep one job. A later full re-run uses `--time` = ceil(1.5 x total / 60) minutes, at least 30 minutes.
- Total above 5,400 s: submit the steps as separate jobs (`--export=...,C1_STEPS=09`, then `C1_STEPS=10`, then `C1_STEPS=11`, each with `--dependency=afterok:<previous job id>` and `--time` = ceil(1.5 x step seconds / 60) minutes). If one step alone is above 5,400 s, record it in the run notes; the plan has no split of a step.
- The pilot timed out: read the finished steps in the log, resubmit the unfinished steps with `C1_STEPS` and `--time=8:00:00`, and size again from that run.

**R4. 12 (login node).**

```bash
$ENV_PY phasec/12_report.py
```
Accept: exit 0 and `report.md: <n> lines`. 12 itself stops if a number has no source in `metrics.json` or `findings.json`.

**R5. Determinism of 11.**

```bash
cp "$STEP1_WORKDIR/phasec/metrics.json" "$STEP1_WORKDIR/phasec/logs/metrics.pilot.json"
sbatch --export=ALL,PROJ_ROOT="$PROJ_ROOT",STEP1_WORKDIR="$STEP1_WORKDIR",C1_STEPS=11 \
  -o "$LOG/c1.%j.log" -e "$LOG/c1.%j.log" "$S1/phasec/c1_evaluate.sh"
cmp "$STEP1_WORKDIR/phasec/metrics.json" "$STEP1_WORKDIR/phasec/logs/metrics.pilot.json" && echo identical
```
Accept: `identical` (same seed, same inputs). The 1-versus-2-worker identity of 10 is covered by `test_workers_give_the_same_scores`.

**R6. Hand-over.** Give the owner `phasec/report.md` with `metrics.json`, `findings.json` and the run notes (wall times, cluster count, removal counts, labels). The independent reviewer checks the code and this run (spec 7, definition of done). No gate value, card or README accuracy statement follows from this run (spec 1).

## Verification done while writing this plan

A prototype of every file was built in a scratch copy of the worktree (all tracked files of `step1-eval` at `7b735ec`, outside the repository). The code blocks of this plan are those files after `ruff check` and `ruff format` (ruff 0.3.5): `All checks passed!`.

**Plan-review fixes (2026-10-01).** The fixes of the plan review and the owner answers (rulings C-8, C-9 amended, C-12 to C-15; reviews I-1, I-3, M-1, M-2, M-5, M-3/M-4/M-7 caveats; the AVX2 rule) were prototyped in a second scratch copy (tracked files of `step1-eval` at `0b05bdd`, plan code extracted, changed, then copied back into this plan). A fresh extraction of every code block of this plan into another copy gives files identical to the prototype (`diff -r`), and the golden command there gives the same SHA-256 (`233a0333...83c8`). Before the fixes the same extraction reproduced the old golden SHA-256 `5f55c304...26a8` and `526 passed, 2 skipped`.

| Check | Result |
|---|---|
| Full suite, `PYTHONPATH=src $ENV_PY -m pytest tests/step1_compare -q` (Python 3.14.2) | `533 passed, 2 skipped` in 663 s (after the plan-review fixes; before them `526 passed, 2 skipped` in 449 s); skipped: one old test and the real-MMseqs2 test |
| Stdlib suite, `$PY -m pytest tests/step1_compare -q` | `484 passed, 1 failed, 8 skipped` in 209 s (Python 3.12.14 with `~/.local` numpy 2.5.2, scipy 1.18.1, scikit-learn 1.9.1, so the Phase C tests ran, the golden file included); the failure is the old `test_chunk_manifest_columns_equal_the_code_constant` (see "How to run the tests") |
| New tests per file | scaffold 9, build 17, metrics 12, bootstrap 7, rules 9, models 9, splits 17 (+1 skipped), fit 9, evaluate 13, report 6, jobs 9, docs 7; 124 new tests pass (117 before the fixes; added: `test_t_grid_is_0_20_to_0_40`, `test_report_states_the_caveats`, `test_literature_section_reports_recall_only`, `test_avx2_tools_have_a_cpu_constraint`, `test_avx2_check_flags_a_planted_script`, `test_c1_keeps_the_wall_times_after_a_failed_step`, `test_grids_in_columns_md_equal_the_code`; extended: `test_recall_at_fpr`, `test_cut_off_values_parse_equal_to_the_grid`, `test_fitted_settings_are_recorded`, `test_undefined_metrics_are_null`, `test_estimate_label_rule`, `test_label_in_metrics_follows_the_half_widths`, `test_report_quotes_findings_and_labels`, `test_c1_requests_cpu_on_an_avx2_partition_with_a_time_limit`, `test_readme_names_the_phase_c_commands_and_rules`) |
| Each task's tests before its implementation (scratch copies with the files of that task and later tasks removed) | the "Expected: FAIL" lines of the tasks are these runs; after the fixes re-run for Tasks 10 (`2 failed, 4 errors`), 11 (`6 failed, 1 passed, 1 skipped`) and 12 (`7 failed`); Tasks 1 to 9 not re-run (their red step is a missing module or script, which the fixes do not change) |
| Real MMseqs2 (non-AVX2 binary) on the fixture | `test_phasec_splits.py`: 18 passed in 46 s (before the fixes; not re-run, Task 7 did not change) |
| Mutation checks before the fixes (each change made in a copy, then the named test run) | 16 of 16 mutations make their test fail: labelmap pm-unresolved -> neg -> `test_class_of_truth_table`; no GO-over-T-c precedence -> `test_tc_row_yields_to_go_label`; no M3 check -> `test_stale_input_stops`; metrics: precision at first threshold reaching recall -> `test_precision_at_recall_reads_the_step_curve_with_ties`; bootstrap per protein -> `test_bootstrap_resamples_whole_clusters`; R2 uses > t -> `test_rule_truth_table`; M8-C ignores the C-terminal row -> `test_cterm_variant_row_selection`; no T-c rule b -> `test_tc_excludes_test_cluster_mates`; no T-c rule c -> `test_tc_excludes_test_taxa`; self hits kept -> `test_max_identity_excludes_self_hits`; 10 trains on test labels -> `test_threshold_fit_uses_train_only`; 11 no direct filter -> `test_direct_filter_applies_to_test_rows_only`; score_source always final -> `test_proteome_calls_use_oof_for_training_hashes`; estimate rule uses < -> `test_estimate_label_rule`; beats reads hi -> `test_findings_booleans`; report regex without '-' guard -> `test_report_numbers_come_from_metrics`. Not re-run after the fixes |
| Mutation checks of the plan-review fixes (each change made in a copy, then the named test run; a kill counts only when the test itself fails, not its fixture) | 15 of 15 make their test fail: `T_VALUES` keeps 0.10 and 0.15 -> `test_t_grid_is_0_20_to_0_40`; old C grid -> `test_fitted_settings_are_recorded`; no count floor -> `test_estimate_label_rule`; floor 19 -> `test_estimate_label_rule`; literature precision kept -> `test_undefined_metrics_are_null` (`1.0 is None`); literature precision at the rule's recall kept -> `test_undefined_metrics_are_null`; literature section prints all columns -> `test_literature_section_reports_recall_only`; no caveat lines -> `test_report_states_the_caveats`; findings without "Any" -> `test_report_quotes_findings_and_labels`; `recall_at_fpr` without the empty call set -> `test_recall_at_fpr` (`-inf == 0.0`); no EXIT trap in C1 -> `test_c1_keeps_the_wall_times_after_a_failed_step` (no `wall.m1test.txt`); no `--constraint` in `c1_evaluate.sh` -> `test_avx2_tools_have_a_cpu_constraint`; no `--constraint` in `09_cluster_and_split.sh` -> the same test; `hard_negative` not positive -> `test_literature_rows`; old C grid in COLUMNS.md -> `test_grids_in_columns_md_equal_the_code` |
| Leakage probe | shuffled labels, 600 x 300: out-of-fold AUC 0.5143 (seeds 21 to 25: 0.493 to 0.532); every fold trained on all rows: 0.970 (before the fixes; not re-run) |
| 08 on the real inputs (scratch links) | Facts table; 17.5 s |
| 09 on the real inputs (non-AVX2 binary, 2 threads, c01) | 9,737 clusters, removal counts in Facts; 18 min 48 s |
| One real unit of 10 (S1 fold 0, V-go, all candidates, 1 BLAS thread) | about 45 s fit time with the old grids (Facts); not run with the new grids |
| Bootstrap at real size (S1 pooled direct rows) | 28.0 s for stratum `all` and 20 candidate-variant pairs |
| Fixture chain (08, 09 stub, 10 all candidates, 11 with 50 resamples, 12) | 08 + 09 13 s; 10 78 s (under load); 11 13 to 25 s; report 2,390 lines after the fixes (2,412 before) |
| Not run | 10 and 11 on the real data (old or new grids); C1 on epyc; the AVX2 MMseqs2 build (it cannot run on c01); any SLURM submission; `sbatch` parsing of `--constraint=ryzen` (no submission). Checked read-only: `sinfo -p epyc -o '%N %f'` lists `r[21-23,25-40] ryzen,amd,milan` (19 nodes) |

## Self-review

1. **Spec coverage.** Spec 3.1: `labelmap` (Task 1). 3.2: direct filter on test rows only (Tasks 8, 9), literature rows and HSP60 (Task 2), C-9 (Tasks 2, 9). 3.3: candidates (Tasks 5, 6), g and t grids (Task 5), C grid and `class_weight` (Task 6), nested protocol for scaler, C, g, t, threshold and H variant (Tasks 6, 8), Youden (C-3; Tasks 5, 6), both variants (Task 8). 3.4: table and dedupe and precedence C-6 (Task 2), clusters (Task 7), S1 to S3 (Task 7), T-c removals a, b (C-5), c (C-7) with `tc_taxon_clades.tsv` and its completeness test (Tasks 2, 7), maximum identity with `fident`, `-s 7.5`, coverage and self-hit exclusion (C-4; Task 7). 4: bootstrap 2,000 by cluster, paired, seeded (Tasks 4, 9), metrics list (Task 3), estimate rule C-8 with FPR half-width (Task 9), prevalence table C-11 (Tasks 3, 9), calibration C-10 (Tasks 6, 9), variant effect (Task 9), agreement and `proteome_calls.tsv.gz` with `score_source` (Task 9), named panel (Task 9), findings (a) to (c) (Task 9), context note (Tasks 9, 10). 5: layout, allowed imports, STOP contract, run JSONs, M3 check, compute (Tasks 1, 2, 11). 6: every named test exists: `test_class_of_truth_table` (1), `test_alternate_files_excluded`, `test_tc_row_yields_to_go_label`, `test_truth_set_has_no_iea`, `test_tc_taxon_table_is_complete` (2), `test_metrics_hand_computed` (3), `test_bootstrap_resamples_whole_clusters`, `test_bootstrap_seeded`, `test_paired_bootstrap_uses_same_indices` (4), `test_rule_truth_table` (5), `test_shuffled_labels_give_chance_auc`, `test_planted_signal_is_recovered`, `test_cterm_variant_row_selection` (6), `test_no_cluster_spans_train_and_test`, `test_test_truth_sources`, `test_tc_excludes_test_proteins`, `test_tc_excludes_test_cluster_mates`, `test_tc_excludes_test_taxa`, `test_max_identity_excludes_self_hits` (7), `test_threshold_fit_uses_train_only` (8), `test_direct_filter_applies_to_test_rows_only`, `test_estimate_label_rule`, `test_prevalence_table_formula`, `test_proteome_calls_use_oof_for_training_hashes`, `test_findings_booleans`, `test_golden_metrics` (9), `test_report_numbers_come_from_metrics` (10), `test_stale_input_stops` (in Tasks 2, 7, 8, 9, 10, one per script), `test_phasec_imports_only_allowed` (1). 7: E1 to E4 (Tasks 1 to 12). Not covered by design: gate values, card, Basidiomycota curation (spec 1).
2. **Placeholder scan.** No TBD or "similar to" steps; every code step holds a whole file, an exact append or an exact diff. The golden file is generated by a command (Task 9 Step 6) because it is binary; its prototype SHA-256 is given.
3. **Type consistency.** `evalio.write_outputs` and `require_current` are used with the same names in 08 to 12; `split_members.tsv.gz` parts (`train`, `train_tc`, `test`, `test_tc`, `test_lit`) are read with these names in 10 and 11; unit keys `<split>|<fold>|<variant>` in `scores_run.json` and in the tests; `models.run_task` task keys equal those that `10.tasks` builds; `metrics.scored_summary` keys equal `11.SCORED`; candidate names are `models.CANDIDATES` everywhere (`M8-C`, `M35-C` with a hyphen).
4. **Review Focus.** Each of the five lines has its test in the owning task (named in the section).

## Spec issues found while planning

1. **The estimate rule accepts zero-width intervals (defect).** The percentile interval of recall has width 0 when every resample gives the same recall, for example when all positives of a small test set are found (H99 has 7 direct positives). Such a test set meets "half-width 0.10 or less" and gets the label `estimate`, although 7 positives cannot support it. In the fixture the literature set (2 positives, all found) was labelled `estimate` this way. The plan keeps ruling C-8 as written and adds `zero_width_recall_interval` beside the label in `metrics.json` and `report.md`. Owner decision needed: add a minimum count of direct positives, or use an interval that does not collapse (for example Wilson), or accept the rule. **Resolved (owner, 2026-10-01):** a count floor of 20 direct-evidence positives (ruling C-8 amended); `n_direct_positives` and `floor_met` are stored beside the label, and `zero_width_recall_interval` stays.
2. **Literature positives versus ruling C-9 (ambiguity).** Spec 3.2 defines literature positives as the rows with an accession, a sequence and `moonlighting` not `YES` (19 rows; this includes the 9 `hard_negative` rows). Ruling C-9 says the `hard_negative` rows are not negatives and "their scores appear in the named panel only", which can also be read as "not positives either". The plan follows the explicit 3.2 definition (19 positives, as in the parent's "about 19 usable") and also lists the `hard_negative` rows in the named panel. If C-9 means "panel only", the literature set has 10 positives; the change is one line in `08.literature_rows`. **Resolved (owner, 2026-10-01):** the `hard_negative` rows count as positives (19 literature positives; HSP60 never); ruling C-9 is amended to say so, and the named panel still lists every `hard_negative` row.
3. **Findings (b) and (c) do not name the ML candidate.** The plan evaluates each ML candidate (M8, M35, M8-C, M35-C, H) against B1 and R2 and reports `holds_for`; (c) holds for a candidate only on all three S2 test sets. Comparator FPRs are read at the threshold where the comparator reaches the recall of R2, per resample (paired). **Resolved (ruling C-13):** no headline ML candidate; the report says "any" and prints the per-candidate table.
4. **No criterion for C.** Spec 3.3 says C is "chosen on the inner folds" without a criterion; the plan uses the inner out-of-fold PR-AUC, the criterion that the spec gives for H's variant (ties: the smaller C). On the real S1 fold 0 every LR candidate chose C = 0.01, the smallest value of the grid; the grid may be too narrow at the low end (not tested further). **Resolved (ruling C-12):** the grid is {0.001, 0.003, 0.01, 0.1, 1, 10}; the spec now names the criterion (inner out-of-fold PR-AUC, ties: the smaller C).
5. **S2 test sets.** Spec 3.4 item 4 says S2 and S3 test sets are fixed by role (`test_species`, `test_clade`, `undecided`); the parent spec 4 step 8 also has "train S288C, test *C. albicans*; the reverse". The plan has both swaps (S2-Calb_CGD, S2-Scer_SGD) and S2-Spom_PomBase, and C-7 removes the T-c rows of the test species in each.
6. **Cluster check for S2 and S3.** `test_no_cluster_spans_train_and_test` cannot apply to the GO training rows of S2 and S3: orthologs across species share clusters (by design of leave-species-out; ruling C-4 measures it). The plan applies the check to S1 folds and to the T-c training rows of S2 and S3 (C-5), and says so in `splits.py`.
7. **Dedupe with an excluded gene.** Spec 3.4 item 2 covers a hash in pos and neg; it does not cover pos or neg with an excluded gene (2 real hashes). The plan makes such a row `excluded` (Task 2). **Resolved (ruling C-14):** the spec now states this rule.
8. **`in_sample` cannot occur.** Every hash of a FULL training table is in an S1 fold (T-c rows are assigned by cluster), so `score_source` is `oof` or `final`. The rule stays as a guard; R2 checks that `in_sample` does not appear.
9. **Rule b in S1 removes nothing.** S1 folds are formed by cluster, so no T-c training row can share a cluster with the test fold. The code still applies rule b and logs 0.
10. **Counts.** "About 23,400 (reviewer 23,413)" sequences to cluster: measured 23,351 (23,338 table rows and 13 literature-only sequences). V-go training pool: 308 pos and 6,842 neg unique sequences (spec: 313 and 6,951 members before dedupe). Direct-evidence *C. albicans* positives in the test set: 153 (spec table: 154; one positive merged into an excluded row, item 7). **Resolved (ruling C-14):** the spec explains 308 / 6,842 against 313 / 6,951: positives 313 - 2 made excluded - 3 merged; negatives 6,951 - 75 without a feature row (`Calb_CGD` 72, `Scer_SGD` 3) - 34 merged (counts from the prototype `build_run.json` and `eval_dedupe_log.tsv`).
11. **Environment (not a spec defect).** Spec 3.4 item 3 names `module MMseqs2/17-b804f`. Its build needs AVX2 and dies on the abu_dhabi nodes; C1 requests partition `epyc`. **Resolved (ruling C-15):** both phasec shell scripts request `-p epyc --constraint=ryzen`; `test_avx2_tools_have_a_cpu_constraint` scans every shell script of `analysis/step1_compare/`.
12. **Output location (deviation).** Parent spec 6 says `results/step1_compare/metrics.json`. `results/*` is git-ignored and the Phase A README says data never go into the repository, so the plan writes every Phase C output to `$STEP1_WORKDIR/phasec/`, as Phase B did for its outputs. The owner can copy `metrics.json` and `report.md` where D7 needs them.
13. **Observation for the owner (one unit only).** On S1 fold 0 the Youden fit gave R2 t = 0.10, the lowest grid value; then R2 is close to "SP and Ser+Thr at least 10%", which most SP proteins meet. This is a property of Youden's J (ruling C-3), not a code error; the R0-to-R2 rows of the report show the effect. **Resolved (owner, 2026-10-01):** the `t` grid is 0.20 to 0.40; the plan review measured that 96.7% of the SP proteins have `ser_thr_frac` >= 0.10, so R2 at t = 0.10 equals R0 on the real data.

## Deviations from the suggested decomposition (each with a reason)

- Four helper modules beyond the seven of spec 5: `evalio` (STOP error, hashes, run JSON; every script needs it), `universe` (per-sequence inputs; `models` would otherwise mix file reading and fitting), `findings` and `agreement` (they keep `11_evaluate.py` readable and are unit-tested alone).
- `test_threshold_fit_uses_train_only` is in Task 8, not Task 6: outer test labels exist only once 10 builds the units from `split_members.tsv.gz`. Task 6 has the model-level checks (chance level, planted signal, inner folds).
- `test_direct_filter_applies_to_test_rows_only` checks the test side in Task 9 and the training side in Task 9 and Task 8 (`test_training_rows_keep_homology_only_rows`).
- Task 11 has one SLURM job (C1) for 09 to 11 and a wrapper `09_cluster_and_split.sh` that the job calls; spec 5 names the wrapper, and the global rule asks for as few jobs as the dependencies allow.
- `10_fit_and_score.py --candidates` and `11_evaluate.py --n-resamples` exist for tests, the golden file and the smoke run; the job uses all candidates and 2,000 resamples by default.
- The golden file is a gzip JSON produced by a command in Task 9 (binary data cannot be pasted); the plan records its prototype SHA-256.

## Assumptions

- The AVX2 MMseqs2 build on epyc gives the same clusters as the non-AVX2 build on c01 (not verified; R2 records the count).
- The C1 run time is under 1 h on 16 epyc cores (estimate from c01 measurements; the pilot measures it).
- 48 GB is enough for C1 (estimate; 16 workers with mmap-shared embeddings).
- `StratifiedGroupKFold(shuffle=True, random_state=20261001)` gives the same folds with scikit-learn 1.8.0 on every node (the library version is in every run JSON).
- Literature rows count as direct-evidence truth (they have no homology flag).
