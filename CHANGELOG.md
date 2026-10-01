# Changelog

All notable changes to this project are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- Continuous integration: lint (ruff) and unit tests on every push and pull request.
- Model card (JSON) written next to each trained model; `predict` refuses a `--model-name` that differs from the card.
- Curated label tables, a stage-2 classifier prototype and model-review analyses (PR #18).

### Changed
- `predict` writes the score of every sequence to the output CSV. `--show-all` now only controls what is printed (PR #23).
- Embeddings average over residue tokens only, so a sequence's embedding no longer depends on the other sequences in its batch (PR #18). Models trained before this change (`models/*.pkl`) were trained on the old pooling; see issue #25.

### Fixed
- Batches that failed to embed no longer shift labels silently (PR #18).
- Output ids that contain commas are quoted correctly (PR #18).

### Removed
- The conda build step in the release workflow. The repository has no conda recipe.
