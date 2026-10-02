# Independent design and test review (Fable, 2026-09-30)

*Requested before any new tool or model code. Reviewer read docs, code and GitHub issues and ran
`pytest tests/adhesion_predict` (11 passed). Items marked "inferred" were not measured. I checked
items 30 and 31 myself (missing workflows, missing CHANGELOG, cwd-based model path); the rest is
the reviewer's report. Step 1 job 29301182 had not finished when this was written.*

Status codes: FV fixed and verified, FU fixed but unverified, OP open with a plan, ON open with no solution.

## 1. Register

| # | Problem | Where | Status |
|---|---|---|---|
| 1 | Padding in mean-pool | `embeddings.py:108-113`; test is CPU, 2 sequences | FV (CPU only) |
| 2 | Skipped batch misaligns labels | `embeddings.py:177-195`, `train.py:103` | FU: no test of the skip/retry path |
| 3 | Layer 6 hard-coded; card not read back | `embeddings.py:12`, `predict.py:47-54` | OP (#9 plan 3.3) |
| 4 | Truncation at 1022; 57% of positives cut; no truncation count | `embeddings.py:14,163` | OP (#10, plan 3.3) |
| 5 | Bare pickle, no scaler, C untuned | `model.py:32,45,59` | OP (#9) |
| 6 | Shipped models have no card | `models/*.json` absent | OP (#9, #25) |
| 7 | Only positive calls written | PR #23 | FV on main |
| 8 | Random-split CV on a redundant set; `evaluate.py` scores its own training dirs | `model.py:24-41`, `evaluate.py:80-90` | ON (#12 names it, no design) |
| 9 | `features.py` no-op | test `treats_j_as_l` | FV |
| 10 | File handles / case-duplicate globs | `io.py:50-53,121`; test does not cover the case-duplicate case; `load_sequences_from_dir:148` lower-case only | FU |
| 11 | README command; FungiDB 401 | README:16 | FV (401 not retested) |
| 12 | Negative-sampling seed unknown | PR #23 logs future seeds; the original is unrecoverable | OP partial |
| 13 | Batch size ignores sequence length | `embeddings.py:50-81` | ON, no issue |
| 14 | 167 exact duplicates | PR #24; shipped pickles trained with them | FV for new training |
| 15 | Legacy pickles scored with new pooling (#25) | job 29301182 | OP, measurement pending |
| 16 | Evaluated model (scaled) differs from shipped (unscaled) | plan fact 1 | OP (step 1 job) |
| 17 | Threshold 0.5; training prior 414:1500 vs ~0.3% genomic | `model.py:113` | ON until T5 |
| 18 | Positives are FLO/ALS homologs; ~12% adhesin precision; scope is step 1 | review 4.3, ARCH 2.0 | OP (R1), proposal only |
| 19 | Rename CLI, column, label | ARCH 2.0 | OP, no migration spec |
| 20 | Fungi_5k survey used pre-fix code, 35M model, p>0.5 rows only | STATUS 5.6 | OP (#11, #15, #16) |
| 21 | Adhesin families are clade-specific | review 4.5 | design decided, no code |
| 22 | Stage-2 prototype unpackaged | STATUS 5.4 | OP (#15) |
| 23 | ESM+LR vs SignalP+GPI agreement not measured (R3) | plan 8 | ON |
| 24 | Stage 1 SignalP under-calls (4%; SOWgp none) | ARCH Stage 1 | ON |
| 25 | Repeat-detector floor; `14` loses Pro/Cys candidates | ARCH "What the repeat detector measures" | ON |
| 26 | `src/.../models/*_t6_8M.pkl` untracked though packaged | git status, #26 | ON |
| 27 | `config.ESM2_MODELS` lists non-existent models | `config.py:41-46` | OP (delete) |
| 28 | Curated labels are a draft (`needs_review=yes` nearly all); `adhesins.tsv` row 1 has no accession | `data/curated/README.md` | ON (#13, #14) |
| 29 | Repo hygiene; 15 MB tracked file | #26 | ON |
| **30** | **No CI test workflow.** `publish_release.yml:10,13` call `build_and_test.yml` and `conda_build_and_test.yml`, which do not exist. `version_bump.yml` needs `CHANGELOG.md`, which does not exist | `.github/` | ON, untracked |
| **31** | `MODELS_DIR` and `DATA_DIR` resolve from the working directory; a `./models` in cwd silently replaces the packaged model | `config.py:16,31` | ON, untracked |
| **32** | Dedupe keys on (label, sequence): a sequence in both classes stays in both, and the test codifies this | `train.py:47` | ON, untracked |
| **33** | AGENTS.md says the model cache is thread-safe; it is a plain dict | `embeddings.py:9` | doc error, untracked |

## 2. What the 11 tests prove, and what they miss
- Real but narrow: batch independence (CPU, 2-3 sequences), `sanitize_sequence`, layer range, CSV quoting, J-to-L, gzip read.
- Card test proves JSON write and read only. It does not prove `predict` refuses a mismatch, because `main()` is untested end to end.
- `find_fasta_files` test uses distinct names, so the case-duplicate bug is not exercised.
- No test for: card enforcement, legacy pickle load, `train_classifier`, `evaluate`, CLI arguments, truncation, retry/skip alignment, `sequences_sha256`, parallel file processing, or any golden score on a fixed fixture.
- Nothing in CI runs any tests (item 30).

## 3. Tier testability (review section 6)
Only T0 has a binary criterion. No tier has a numeric acceptance value, and none is wired to the model card.
T1 is saturated at 0.999 and cannot fail. T3 (only Saccharomycotina plus 22 Eurotiomycetes labels) and T6 (needs #16) are not testable now.

## 4. Where "adhesin" now misleads
`pyproject.toml` (name, description, 3 entry points); `predict.py` banner, `Adhesion` labels, `probability_adhesion`, output suffix `.adhesion_predict.csv`; `evaluate.py:62`; `train.py` banner and model filename `adhesion_model_<esm>.pkl`; README title and lines 3, 16, 33, 40, 55; AGENTS.md; `analysis/kingdom_survey` (`adhesion_fraction`, 25 references in its tests); `analysis/adhesion_properties`; `docs/superpowers/*`; `run_fungi5k*.sh`. `surface.tsv` `adhesion_status` is a stage-2 label and is legitimate.
Migration is one sentence, not a plan: no dual-column period, no downstream migration for `kingdom_survey`, no filename plan, no CHANGELOG, no decision on relabelling existing Fungi_5k CSVs.

## 5. Risks not yet tracked
- The card check never fires for the default model (no card), so the "assume defaults" path is the normal path.
- `torch`, `fair-esm`, `scikit-learn` are unpinned; ESM weights are fetched by name without a hash. sklearn 1.8.0 loads the pickles without a warning, so the training-time version is unknown.
- Tests have run under Python 3.9 and 3.14 with no canonical interpreter recorded.
- Legacy-pooling emulation assumes batch 4 in file order and cannot be verified against the Feb 2026 training.
- Inferred, not measured: GPU fp32 could differ from CPU by more than 1e-5.
- Checkout note: branch `issue-9-plan` was cut before PR #23, so code work there must start from current `main`.

## 6. Recommended order
1. Sync with main (#23), merge #24. Add `.github/workflows/test.yml` (ruff + `pytest tests/adhesion_predict`, CPU torch, cached 8M weights). Fix or delete the broken release workflows. Add `CHANGELOG.md`.
2. Finish the step 1 measurement (job 29301182) and report it, including the pooling mismatch.
3. #9/#25: bundle and card v2, tests first: card-driven refusal (layer, pooling, max_residues), legacy pickle warns, end-to-end `main()` on a 5-sequence fixture, truncation count.
4. Rename (R1) in the same release: both column names for one release, alias entry points, model filename with a `legacy/` copy. Tests: CSV header, alias import, legacy load.
5. Retrain the step-1 model (residue pooling, dedupe). Acceptance script on S288C: recall on the 9 known adhesins, calls on the 20 hard-negative panel, total calls.
6. `analysis/evaluation/` harness (#12): T0 in CI; T1/T2/T4/T5 as HPCC scripts that emit JSON for the card; each with a numeric gate.
7. Then #10, #15, #16.

Layout: `tests/adhesion_predict/` fast unit tests (CI); `tests/integration/` marked `slow` with a golden-score file; `analysis/evaluation/` for the tier harness (HPCC, `mmseqs`, GPU). Results enter a card only through the harness.

## 7. Decisions needed from the owner before coding
1. Final names: package, entry points, column, label, model filename, alias duration.
2. Threshold: keep 0.5, or choose on T5 at a stated target precision.
3. Numeric acceptance values: T0 tolerance, T4 FPR, T5 recall and precision on S288C.
4. What step 1 is validated against: SignalP+GPI, curated `surface.tsv`, or both (R3).
5. Dependency pins, canonical Python version, whether CI ever uses a GPU.
6. Whether curated labels may serve as test fixtures before expert review.
7. Fungi_5k outputs from the legacy model: freeze as "legacy pooling, 35M", or re-run.
8. Repo hygiene (#26): track the 8M pickle or stop shipping two copies; the 15 MB file; the cross-label duplicate rule.
