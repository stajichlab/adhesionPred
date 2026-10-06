# Changelog

All notable changes to this project are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- `cellsurface_sorting_hat` core engine: reads module result tables, applies three-valued category rules and writes calls, evidence, a report and run records.
- `cellsurface_sorting_hat_module`: wrappers that turn tool output (SignalP, Pfam `hmmsearch`, repeat detectors, BLASTP against the IUIS allergens, TMHMM) and lookup tables (antigen ranking, Cys-rich tiers, spherule expression) into module tables. Each input is checked against the FASTA before any output is written. A run that ends with no usable result is refused, and some `error` rows make the run `partial`.
- `cellsurface_sorting_hat_calibrate`: `phasec` (status entries for R0 from the Phase C metrics, one per species, with cluster counts from the Phase C cluster files and a check of the SignalP module and mode), `truth` (sensitivity and specificity of one call that reads one module against a truth table, with a leakage cap, a taxon check and a check of the `run.json` module identity), `pfam-specificity`, `allergen-lso` and `panel`.
- Job scripts in `scripts/sorting_hat/` and the data tables `data/sorting_hat/family_table.tsv` (all families inactive) and `phasec_set_species.tsv`.
- No calibration has been run on HPCC data yet (Tasks 11 to 15 of the plan). Every module status is `unvalidated`. `calibrate phasec` can write R0 status entries from the Phase C numbers, but it has not been run in a real work directory.
- Continuous integration: lint (ruff) and unit tests on every push and pull request.
- Model card (JSON) written next to each trained model; `predict` refuses a `--model-name` that differs from the card.
- Curated label tables, a stage-2 classifier prototype and model-review analyses (PR #18).
- `predict` reports how many sequences are truncated at 1022 residues. Cards of new models record library versions.

### Changed
- **BREAKING:** the package `adhesion_predict` is now `surface_glyco`. Entry points: `surface_glyco_predict`, `surface_glyco_train`, `surface_glyco_evaluate`. Output column `probability_adhesion` is now `surface_glycoprotein_score`; labels `Adhesion` / `Non-adhesion` are now `surface_glycoprotein` / `other`; default output files end in `.surface_glyco.csv`; model files are `surface_glyco_model_<esm>.pkl`. No aliases, and result files written by the old package are no longer read by the analysis code. Reason: the model scored FLO/ALS-like surface glycoproteins, not adhesins.
- **BREAKING:** the two models that shipped with 0.1.0 are removed. They were trained on FLO/ALS homologs versus random proteins with padding-inclusive pooling (issue #25). No model ships until a surface-glycoprotein model is trained and validated.
- `predict` and `evaluate` take the ESM model, layer and pooling from the model's card and stop with a message naming the field on a mismatch, a missing card or an unsupported pooling.
- A `./models` directory in the working directory no longer replaces the packaged models. Set `SURFACE_GLYCO_MODELS_DIR` to use another directory. `surface_glyco_train` writes to `./models`. `surface_glyco_predict` and `surface_glyco_evaluate` need `--model` or `SURFACE_GLYCO_MODELS_DIR` to find that model.
- `requires-python` is now `>=3.11`.
- `kingdom_survey` counts called proteins instead of all rows in a result file, and finds `*.surface_glyco.csv`.
- `predict` writes the score of every sequence to the output CSV. `--show-all` now only controls what is printed (PR #23).
- Embeddings average over residue tokens only, so a sequence's embedding no longer depends on the other sequences in its batch (PR #18). The models trained before this change (`models/*.pkl`) are removed (see the BREAKING entry above and issue #25).

### Fixed
- Training drops sequences that occur in both classes and empty sequences.
- Batches that failed to embed no longer shift labels silently (PR #18).
- Output ids that contain commas are quoted correctly (PR #18).

### Removed
- The conda build step in the release workflow. The repository has no conda recipe.

## [0.1.0] - 2026-02-16

Initial tagged version (`v0.1.0`).
