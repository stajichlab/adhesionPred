# AGA1-class modules for the sorting hat and a cross-proteome aggregator: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax. Tests come first (red, then green); each task ends in a commit; each work package is one PR.

**Goal:** Make the AGA1-class predictors (PredGPI, AIUPred, the broad and narrow rule flags, the feature classifier) modules and calls of `cellsurface_sorting_hat`, with statuses, and add `cellsurface_sorting_hat_aggregate`, which turns the outputs of many runs plus a grouping of taxa into tables, a Markdown report and static figures.

**Architecture:** Three wrappers in `src/cellsurface_sorting_hat/modules/` (`gpi.py`, `disorder.py`, `aga1_model.py`) registered in `modules/cli.py`; a signed thresholds table `data/sorting_hat/aga1_model_thresholds.tsv`; six calls in `categories.yaml`; status files under `docs/reports/data/sorting_hat/call_status/`; a new `src/cellsurface_sorting_hat/aggregate/` package with its own CLI. The engine does not run tools and gets no clade logic. No scikit-learn or matplotlib import in the engine or the wrappers; matplotlib is imported only by the aggregator figure code.

**Tech Stack:** Python 3.11+, PyYAML, standard library for the wrappers; pandas or csv for the aggregator; matplotlib for figures (already a dependency); pytest; ruff 0.3.5 (CI pin).

**Spec:** `docs/superpowers/specs/2026-10-26-aga1-sorting-hat-modules-and-aggregator-design.md` (revision 2, after the adversarial review). Engine spec: `docs/superpowers/specs/2026-10-04-orchestrator-design.md`. Owner decisions: memory note `aga1-sorting-hat-integration-decisions-2026-10-26` (all nine accepted, including the revised status plan of 2.4, accepted 2026-10-26).

## Global Constraints

- The engine never runs tools. A wrapper reads one tool's output and writes `<workdir>/modules/<module>.tsv.gz` + `.json` through `modules/base.write_module`. Required columns `id`, `state`; a protein the tool did not return gets `state: error`; an invalid protein `invalid_row`. (Engine spec, module contract.)
- No silent zeros: a value that does not apply is empty or `not_assessable`, never 0. A protein of 40 residues or fewer has no GPI call: `state: not_assessable`.
- Protein ids are the first token of the FASTA header. Each wrapper stops on an id that does not match its input.
- Streaming: the AIUPred per-residue file (13.76 million rows for 27,811 proteins) is read as a stream. Residues are never copied into the work directory or the repository.
- AIUPred scores up to 0.05 outside [0, 1] are accepted and used as written; further out stops the wrapper.
- Thresholds: `disorder_run_min` = 200 in `thresholds:` of `categories.yaml`. Model thresholds come only from `data/sorting_hat/aga1_model_thresholds.tsv`; the engine reads no model threshold.
- Module identity (`ModuleSpec`): `aga1_model` carries `artefact_digest` = hash of the threshold rows used + sha256 of the score table, so a changed threshold or score table makes status files stale.
- Call rules of the engine (verified): a call with `ref` is derived, cannot have a status file; an ungated call cannot `ref` a per_variant call; no literal step 1 module in a call; `calibrate truth --module` takes a one-module call, `--call-status` a call that reads two or more modules; leakage values are `none, partial, tuned_on_truth, in_reference, unknown`.
- The narrow flag reads `repeat14` only (as measured). The broad and narrow flags and the model call are gated, per_variant, derived.
- Statuses (owner-accepted): `gpi_long_disorder_protein` `smoke` for *S. cerevisiae* (4932) and *C. albicans* (5476) against the class truth, leakage `tuned_on_truth`, notes say it measures agreement with the AGA1 class labels; `gpi_anchored_protein` against UniProt GPI features if the 850 table has enough (decided in Task 7 after counting), else `unvalidated`; `long_disordered_protein` and `model_surface_candidate` `unvalidated`.
- Aggregator: calls are read for the default gate only; per 1,000 **assessable** proteins; `assessable == 0` gives `NA` with a warning; status columns always present; `--measured-only` restricts group summaries; mixed `config_sha256` stops unless `--allow-mixed-config`.
- Python files: ruff `line-length = 100`, CI pins ruff 0.3.5 (`uv run --no-project --with ruff==0.3.5 ruff format . && ... ruff check .`). Report and docs text: Simplified Technical English, no claim that is not measured.
- Golden files and `CLASSES.md`: any change of `categories.yaml` changes the config hash and `tests/cellsurface_sorting_hat/golden/`; regenerate them in their own commit.
- Nothing in this plan is submitted to SLURM. The Pfam and TMHMM arrays are separate and need the owner's go-ahead per submission.

## Work packages (one PR each)

| WP | Tasks | Content |
|---|---|---|
| 1 | 1 | Baseline of the test suite |
| 2 | 2, 3, 4 | The three wrappers, the thresholds table, their tests |
| 3 | 5, 6 | Calls, descriptions, `CLASSES.md`, label rules, golden regeneration |
| 4 | 7 | Local engine run on the 10 proteomes (acceptance), status files |
| 5 | 8, 9 | The aggregator and its tests |
| 6 | 10 | First aggregate report, README and technical notes |

## Task 1: Baseline

- [ ] Run `PYTHONPATH=src pytest tests/cellsurface_sorting_hat tests/class_table tests/aga1_class -q` on a clean checkout of `main`; record the number of tests and the ruff result in the PR text.
- [ ] Confirm `docs/CLASSES.md` is current (`tests/class_table`), and note the config sha256 (`run.json` of the golden run) for comparison after Task 5.

## Task 2: `gpi` wrapper

Files: `src/cellsurface_sorting_hat/modules/gpi.py`, `modules/cli.py` (new subcommand `gpi`), `tests/cellsurface_sorting_hat/modules/test_gpi.py`, a fixture gff3.

- [ ] Tests first (synthetic gff3, bare accessions): `Chain` row gives `call=not_called`; `GPI-anchor` row gives `called` with `omega_site`, `confidence_class` (score 1.0 high, 0.70 probable, 0.55 weak) and `fpr_upper_bound` (0.0015, 0.005, 0.01); a protein of 40 residues or fewer (FASTA length) gives `state=not_assessable`; a protein missing from the file gives `state=error` and run state `partial`; invalid FASTA protein gives `invalid_row`.
- [ ] Stop conditions with named tests: unknown row type; `GPI-anchor` with start != end; `Chain` that does not span the protein; wrong column count; unknown confidence score; duplicate id; id in the file but not in the FASTA; empty file.
- [ ] Implement `parse_gpi(path, lengths)` and `gpi_rows(proteins, parsed)` after `modules/signalp.py`; reuse the validation of `analysis/aga1_class/convert_predictions.py` (same rules, no import from `analysis/`).
- [ ] Register `gpi` in `build_parser`/`run` with `--fasta --workdir --gff3 --predgpi-version`; the run record has the PredGPI script sha256 given by the caller.
- [ ] Run the new tests, then the whole `tests/cellsurface_sorting_hat`; ruff; commit "gpi module wrapper".

## Task 3: `disorder` wrapper

Files: `modules/disorder.py`, `modules/cli.py` (`disorder`), `tests/cellsurface_sorting_hat/modules/test_disorder.py`.

- [ ] Tests first: the AIUPred layout (banner, `#>ID` lines, `pos aa score` rows, blank lines) gives `longest_run`, `frac`, `mean`, `n_res` per protein (threshold 0.5 for a residue; hand-computed values); an IUPred3-style layout (`>ID` lines) is also read; the function streams (a test feeds a generator of lines and asserts memory is not required: it never builds the list of all rows).
- [ ] Stop conditions: score more than 0.05 outside [0, 1]; positions out of order; residue count different from the FASTA length after the trailing `*` strip; duplicate record; record not in the FASTA; FASTA protein without a record gives `error`.
- [ ] Implement `parse_disorder_stream(lines, threshold)` returning per-protein summaries without keeping scores; module columns `longest_run, frac, mean, n_res`; no `call` column.
- [ ] Commit "disorder module wrapper".

## Task 4: Thresholds table and `aga1_model` wrapper

Files: `data/sorting_hat/aga1_model_thresholds.tsv`, `src/cellsurface_sorting_hat/modules/aga1_thresholds.py` (loader and validator), `modules/aga1_model.py`, `modules/cli.py` (`aga1_model`), tests under `tests/cellsurface_sorting_hat/modules/`.

- [ ] Table (columns in the spec 2.3): rows Saccharomycotina (taxon 147537; 4891 is the class Saccharomycetes and does not contain C. albicans, C. auris or C. parapsilosis) logreg 0.56834, rf 0.25448, gboost 0.11414, `calibrated=yes` with the exact meaning of the spec, signed by the owner, date 2026-10-26, evidence `docs/reports/data/aga1_proteome_test_2026-10-25/thresholds.tsv`; fallback rows for taxon 1 with the same values and `calibrated=no`.
- [ ] Validator tests: header exact; one row per (clade, model); the three models present for every clade; thresholds numeric in (0, 1); `calibrated` is yes or no; a signed row needs `signed_by` and `signed_date`; duplicate (clade, model) stops; the Saccharomycotina values equal the committed record `thresholds.tsv` (a test reads both).
- [ ] Wrapper tests (synthetic score table): the deepest ancestor row is chosen with `Lineage.ancestors`; a taxon outside Saccharomycotina gets the fallback row and `calibrated=no`; a taxon not in `nodes.dmp` stops (`TaxonError`); 2-of-3 logic at each combination; a score exactly equal to the threshold counts as at or above; a missing score gives `error`; module rows carry `thr_*`, `threshold_clade`, `calibrated`.
- [ ] Identity tests: changing a threshold in the table, or one byte of the score table, changes `artefact_hash` in the run record; `params` lists the model names and the rule.
- [ ] The wrapper imports no scikit-learn (a test asserts `sklearn` is not in `sys.modules` after import).
- [ ] Commit "aga1_model wrapper and thresholds table". Open WP 2 PR.

## Task 5: Calls, descriptions, `CLASSES.md`

Files: `src/cellsurface_sorting_hat/categories.yaml`, `docs/class_descriptions.yaml`, `docs/CLASSES.md` (generated), `src/cellsurface_sorting_hat/labels.py`, `tests/cellsurface_sorting_hat/test_aga1_calls.py`.

- [ ] Tests first: a synthetic work directory (`engine_helpers.py`) with every combination of `called / not_called / not_assessable` for `gpi`, `disorder`, `step1_rule@R0` and `repeat14`; expected values of the six calls by hand (Kleene logic); `gpi_long_disorder_protein` reads two modules and has no `ref`; the three composites are per_variant, gated and derived; `flocculin_like_candidate` uses `repeat14` and not `tandem_repeat_protein`; `thresholds.disorder_run_min` is read.
- [ ] Add to `categories.yaml`: `thresholds.disorder_run_min: 200`; the six calls of the spec 2.2 in the category "cell wall and adhesion"; the `evidence:` entries `aga1_model.threshold_clade`, `aga1_model.calibrated`, `gpi.confidence_class`, `disorder.longest_run`.
- [x] Primary label: NOT changed (the mechanism lists accept only earlier ungated calls, and adding the new calls would make the existing `other_*` calls `not_assessable` in runs without the new tables); a test pins that the existing labels and calls are unchanged.
- [ ] `docs/class_descriptions.yaml`: description and "known limits" sentence per call (spec 2.2); regenerate `docs/CLASSES.md` with `analysis/class_table/make_classes_md.py`; the CLASSES tests pass.
- [ ] Commit "AGA1 calls in categories.yaml" (code, yaml and docs, without golden files).

## Task 6: Golden files

- [ ] Run the golden test; it fails on the config hash and the new absent-module warnings. Regenerate the golden files and `report.expected.md` with the repo's documented command (see `tests/cellsurface_sorting_hat/golden/README` or `test_cli.py`), review the diff: only the config hash, the new call columns and the "modules without a table" warnings may differ.
- [ ] Commit "regenerate golden files for the AGA1 calls" (alone). Open WP 3 PR.

## Task 7: Local engine run on the 10 proteomes (acceptance) and statuses

Files: `analysis/sorting_hat_run/run_aga1_modules_local.sh` (documented, no SLURM), `docs/reports/data/sorting_hat/call_status/<species>/...json`, `docs/reports/data/aga1_sorting_hat_acceptance_2026-10-26/` (record + README with hashes), tests.

- [ ] Script: for each proteome, build module tables with the new wrappers from the committed or locally held tool outputs (SignalP via the existing `signalp` wrapper, repeat14, PredGPI, AIUPred, the classifier score tables of PR 149, 151) and run the engine; modules without a table (Pfam, TM, allergen) stay absent and give `not_assessable`. If the full Pfam and TMHMM arrays have finished, add the `pfam` and `tm` modules for the family-table calls (optional, noted in the record).
- [ ] Acceptance test (committed small expected counts): the engine's call counts per proteome equal the flag counts of `analysis/aga1_class/` (broad rule 11, 25, 25, 10, 14, 7, 8, 2, 2, 28; narrow 5, 6, 4, 2, 5, 3, 2, 0, 1, 9 for the 10 proteomes of the records; the exact list is read from the committed `summary_by_species.tsv` files), and `model_surface_candidate` equals "2 of 3 models and a signal peptide" of `proteome_summary.py`. Any difference is a bug in a wrapper or a call.
- [ ] Count UniProt GPI feature coverage in the 850 table (`uniprot_gpi_feature`); decide `gpi_anchored_protein` measured or `unvalidated` and state the numbers in the record.
- [ ] Truth tables for *S. cerevisiae* and *C. albicans* from `docs/reports/data/aga1_class_scoping/proteins.tsv` (id, label 1 for tier A, 0 for hard negatives and background, cluster); run the engine on the truth proteins; `calibrate truth --call-status` for `gpi_long_disorder_protein` with `--taxa 4932` and `--taxa 5476`, `--leakage tuned_on_truth`, `--negative-basis assumed`, notes as in the spec; commit the status files; a test loads them and checks status `smoke` and that they are valid for the committed config hash.
- [ ] Re-run Task 6's golden check and the CLASSES test (the cited status files appear in `CLASSES.md`; regenerate it).
- [ ] Commit; open WP 4 PR.

## Task 8: Aggregator core

Files: `src/cellsurface_sorting_hat/aggregate/__init__.py`, `inputs.py`, `groups.py`, `tables.py`, `cli.py`; `pyproject.toml` entry point `cellsurface_sorting_hat_aggregate`; tests under `tests/cellsurface_sorting_hat/aggregate/`.

- [ ] Tests first with synthetic runs (written by the engine test helpers): reading a run directory (the three files); the default gate only (a per_variant call with three variants is counted once); `--variant` selects another; duplicate proteome labels stop; runs with different `config_sha256` stop unless `--allow-mixed-config`; a missing file stops with the file name.
- [ ] Grouping tests: lineage rank grouping with `calibration/cli.read_nodes` ranks; a taxon not in `nodes.dmp` stops; a group-table row overrides the lineage group for that taxon and for descendants without a more specific row; malformed group table stops; optional `names.dmp` gives labels, else taxon ids.
- [ ] Table tests: `counts_by_proteome.tsv` columns and values by hand (called, assessable, per 1,000 assessable, status, status species, measured, unvalidated marker); `assessable == 0` gives `NA` and a warning; `counts_by_group.tsv` (n_proteomes, n_measured, mean/median/min/max over all and over measured proteomes, total called); `--measured-only`; `proteomes.tsv`.
- [ ] Status and evidence: the status column comes from the engine's resolved status in `calls.long.tsv.gz`; `threshold_clade` and `calibrated` come from `evidence.tsv.gz`.
- [ ] Commit "aggregator tables".

## Task 9: Aggregator report and figures

Files: `aggregate/report.py`, `aggregate/figures.py`, tests.

- [ ] Report tests: sections (what was aggregated, group summaries per call category, status coverage, warnings, known limits), warnings for mixed configs, proteomes with many invalid proteins, groups of one or two proteomes, absent module tables, proteomes with many `not_assessable` GPI rows; counts printed beside every group mean.
- [ ] Figure tests: the three PNGs (heatmap with unvalidated cells hatched, group bars with n_proteomes and error bars only for `n_proteomes >= 3`, status-coverage panel) are written, non-empty, referenced from the report; with `metadata={"Software": None}` and fixed fonts two runs give identical bytes (the test compares sha256). Matplotlib is imported only here, with the Agg backend.
- [ ] CLI test (end to end on three synthetic runs); `--help` text; refuses an existing output directory.
- [ ] Commit "aggregator report and figures"; open WP 5 PR.

## Task 10: First aggregate report and docs

- [ ] Run `cellsurface_sorting_hat_aggregate` on the 10 proteomes with `--rank class` and a small group table (Saccharomycotina yeasts, *Candida* clade, Aspergilli, Coccidioides, Basidiomycota) as an owner-overridable example; commit the tables, report and figures with a README (command, inputs, hashes) under `docs/reports/data/aga1_aggregate_2026-10-26/`.
- [ ] Read the report for what it can and cannot say (unvalidated counts, uncalibrated model call outside yeast, absent Pfam and TM modules unless the arrays finished) and write those lines into the report's known limits if the tool did not already.
- [ ] Update `README.md` (new calls, the aggregator command, general text only) and `docs/TECHNICAL-NOTES.md` (conventions, the thresholds table, the acceptance counts); update `AGENTS.md` with the contribution record.
- [ ] Commit; open WP 6 PR.

## Review Focus

- Every `not_assessable` path: absent module table, protein of 40 residues or fewer, wrapper `error` rows; none may become 0 in a count.
- Status staleness: a changed threshold, score table, `categories.yaml` or call definition must make the status files stale (tests in Tasks 4, 5, 7).
- The acceptance equality of Task 7 is the guard against a silent difference between the analysis scripts and the engine.
- Aggregator denominators and the default-gate rule (Task 8) are the likeliest source of double counting.
