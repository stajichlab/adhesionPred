# Plan for #9: model bundle, model card, explicit embedding config

*Drafted 2026-09-30. Revised 2026-09-30 after an independent Opus 5.5 review (section 7). Status: proposal, not implemented.*

## 0. Scope and wording (clarified 2026-09-30)

The model this plan covers is **step 1: a surface glycoprotein predictor**, not an adhesin
predictor (`docs/TOOL-ARCHITECTURE.md` section 2.0). The bundle, card and `predict` output in
this plan therefore describe a *surface glycoprotein score*. Where this plan or the code says
"adhesion", read "surface glycoprotein" until the rename is done. Renaming the CLI, the
`probability_adhesion` column and the `Adhesion` label is a schema change. It is **not** part of
steps 1-6 below. It is listed as open item R1 in section 8 and needs the design review.
Mechanism-specific adhesin classes are separate tools with their own models and are out of scope here.

## 1. Problem, with what is already fixed

Issue #9 has three parts. Two are partly done in `main`.

| Item in #9 | State in `main` |
|---|---|
| Record the embedding config in a card | **Done** (`save_model_card`): ESM model, layer, pooling, max residues, sklearn version, data sha256, counts, holdout accuracy |
| `predict` refuses a mismatched embedding config | **Partly.** It checks `--model-name` only. It does not check layer, pooling or truncation |
| Stop hard-coding repr layer 6 | **Not done.** `DEFAULT_REPR_LAYER = 6` is a module constant. The card records it, but nothing reads it back |
| `Pipeline(StandardScaler, clf)` | **Not done.** `train_classifier` fits a bare `LogisticRegression(max_iter=1000)` |
| Card has CV tiers, results, threshold | **Not done.** The card has one random 80/20 holdout accuracy |
| Tune `C` (currently 1, untuned) | **Not done** |
| joblib instead of pickle | **Not done** |

## 2. Facts that shape the plan

1. **The evaluated model is not the shipped model.** `analysis/model_review/02_cv_and_proteome_eval.py`
   uses `make_pipeline(StandardScaler(), LogisticRegression(C=1.0))`. `model.py` fits an unscaled
   `LogisticRegression`. Every CV figure in the 2026-09-27 review describes the first. The shipped
   `.pkl` files are the second. **Nobody has measured the shipped model with the review's CV.**
2. **Accuracy is unlikely to move.** Reported CV is already ROC-AUC 0.999 for ESM-2 8M, and 0.995
   for amino-acid composition plus length (leave-family-out). Scaling or tuning `C` cannot add much
   there. The value of #9 is reproducibility and safety, not score. This is a prediction, not a
   measurement. The plan includes the measurement (step 1).
3. **Layer 6 means different things.** For `esm2_t6_8M` it is the final layer. For
   `esm2_t12_35M` (used for the Fungi_5k screen) it is the middle layer. Changing the default
   would silently change every existing embedding. The default must stay 6 for existing models.
4. **The shipped pickles were trained on a different pooling than `predict` now computes.**
   The pickles were committed 2026-02-15 (8M) and 2026-09-22 (35M). The pooling fix is commit
   `9dfe519` (2026-09-27). Before it, pooling was `token_representations.mean(dim=1)`, which
   averages BOS, EOS and padding. Now it is `residue_mean`. So `predict` today scores the legacy
   pickles on embeddings they were not trained on. This is a live error, separate from #9's
   goals, and the plan must not preserve it. (Verified in git history; the effect on scores is
   not measured.)
5. **Shipped models have no card.** `models/*.pkl` and `src/adhesion_predict/models/*.pkl` have
   no JSON. They must keep working.
6. **The shipped models were probably trained with 167 duplicate sequences** (PR #24 removes them for
   new training). Not confirmed which data the pickles saw. A retrained bundle is therefore not comparable to the old pickles on the
   training data, only on an external set.
7. `config.ESM2_MODELS` lists `esm2_t6_35M_UR50D` and `esm2_t6_150M_UR50D`. These do not exist in
   fair-esm (real names: `esm2_t12_35M_UR50D`, `esm2_t30_150M_UR50D`). The dict is unused in `src/`
   and `tests/`. Delete it.
8. Only 0.5 thresholding exists. The review (§5.4, §6 T5) says the threshold must be chosen on
   proteome-level data, which needs the T5 inventory (#12, #13, #14). #9 cannot choose it.

## 3. Design

### 3.1 Bundle format

One directory or one file per model, `adhesion_model_<esm>.joblib`, plus `.json` card beside it
(same layout as now, so `model_card_path` is unchanged). The joblib object is a
`sklearn.pipeline.Pipeline([("scale", StandardScaler()), ("clf", LogisticRegression(...))])`.

Loading rules:
- Legacy pickles (no card) get an implied card with `pooling="legacy_all_tokens"`. Silent "assume
  defaults" is not acceptable, because the defaults are wrong for these files.
  **Exact reproduction is impossible.** Legacy pooling averaged padding tokens, so each training
  embedding depended on the other sequences in its batch. The Feb 2026 batch size and order were
  not recorded (the 2026-09-27 review script assumes batch 4 on CPU in file order). An emulation
  can approximate the legacy input but cannot be shown equal to it. So: emulate behind an
  explicit `--legacy-pooling` flag with a printed warning, never as a default, and retrain
  (step 6) so the legacy path can be retired.
- `load_model` accepts `.joblib` and legacy `.pkl`. A bare classifier is returned as is.
- A legacy model with no card prints one warning naming what is unverified.
- Pickle and joblib both execute code on load. State in the README that model files must come
  from a trusted source. Do not add a format that claims to fix this.

### 3.2 Card schema (version 2)

Existing fields stay. New fields:

| field | content |
|---|---|
| `card_version` | `2` |
| `embedding` | `{esm_model, repr_layer, pooling, max_residues, truncation}`; the existing flat fields are kept for version-1 readers |
| `estimator` | `{scaler: bool, C, class_weight, max_iter, random_state}` |
| `selection` | how `C` was chosen: `"fixed"` or `"nested_cv"`, with the grid |
| `validation` | list of `{scheme, n_folds, group_source, roc_auc, pr_auc, recall, fpr}`; empty until a scheme is run |
| `threshold` | `{value, chosen_on}`; `0.5` and `"default (uncalibrated)"` until T5 exists |
| `duplicates_removed` | already added in PR #24 |
| `code_version` | package version, git hash if available |
| `environment` | numpy, scikit-learn, torch, fair-esm versions. Legacy pickles need numpy >= 2 and warn under sklearn 1.9.1 |

The card must not claim a tier that was not run. An empty `validation` list is the honest value.

### 3.3 `predict` compatibility check

`predict` reads `esm_model`, `repr_layer` and `pooling` **from the card** and uses them. The CLI
value is only an override, and a disagreement with the card is an error that names the field.
Today the CLI is the authority and the card only vetoes it (`predict.py:47-54`).
Truncation: the card records the policy (`truncate_at_1022` now; versioned, because #10 changes
it). `predict` counts truncated sequences and reports the count. Model path lookup: try
`.joblib`, then `.pkl`. No card: use the legacy pooling above and warn.

### 3.4 Training

- `train_classifier` builds the Pipeline. `random_state` is a CLI argument (default 42), stored in the card.
- Hyperparameter choice: **default is fixed `C=1.0`**, matching the review. A `--tune-c` option runs
  nested `StratifiedGroupKFold` over a small grid, and records the grid and chosen value.
  Tuning is opt-in because fact 2 predicts no gain.
- Grouping for any CV uses MMseqs2 30% clusters, as in the review (`run.sh`). If `mmseqs` is not
  on PATH, `--tune-c` and validation are unavailable and the card says so.

### 3.5 Out of scope

Threshold calibration (needs T5, #12). Sliding-window embedding (#10). New features (#15).
ESM C / ESM-2 150M support. Retraining and replacing the shipped pickles is a separate,
explicit step (section 4, step 6).

## 4. Steps

Each step is one commit. Steps 1-3 need no GPU.

1. **Measure, on the right axis.** This needs a small code change: give `02_cv_and_proteome_eval.py`
   a `--no-scaler` flag (the scaler is hard-coded at line 43-44). Run on `train_legacy.npy` (the
   features the shipped pickles were trained on) with and without the scaler, and on
   `train_masked.npy` for the current pooling. Nothing is cached, so `run.sh` runs first
   (about 1.5 h CPU, needs `mmseqs`). The embeddings contain the 167 duplicates, which inflate
   random-fold CV. The report must say so. Also score S288C with the old pickle under both
   poolings to measure the live mismatch in fact 4.
2. **Tests first.** Card v2 round trip. Legacy pickle loads with implied `legacy_all_tokens`
   pooling. Mismatched layer or pooling is refused. Pipeline bundle round-trips through joblib.
   Truncated-sequence count is reported.
3. **Legacy pooling in `get_esm_embeddings`, card v2, card-driven `predict`.** Delete
   `config.ESM2_MODELS`.
4. **Pipeline in `train_classifier`; `--repr-layer`, `--seed`.**
5. **Optional `--tune-c`.**
6. **Retrain both bundles; compare with old pickles on S288C before replacing.** Keep old files
   under `legacy/` for one release. Test loading in the target conda env.
7. README, `docs/TOOL-ARCHITECTURE.md`, AGENTS.md.

## 5. Acceptance

- `pytest tests/adhesion_predict` passes, with new tests for each step-2 case.
- `predict` on a v2 bundle and on each legacy pickle both run; a wrong-layer call exits non-zero.
- The card of a newly trained bundle lists every field in 3.2 and sets no value it did not compute.
- Step 1 and step 6 comparisons are in a report with numbers, including the case where nothing changed.
- T0: for a fixed 100-sequence set and a fixed batch size, the `--legacy-pooling` path gives the
  same scores as the pre-`9dfe519` code. Identity is not claimed for other batchings.

## 6. Risks

| risk | handling |
|---|---|
| Retraining changes calls that other work depends on (Fungi_5k survey used the old model) | keep `legacy/` models; record model hash in any survey output |
| Scaler changes score scale, so the 0.5 threshold means something different | compare call counts on S288C before and after; do not ship without that |
| joblib/sklearn/numpy drift breaks loading (legacy pickles fail under numpy < 2) | store all versions; test load in the target env |
| Card claims more than was measured | schema allows empty `validation`; no default metrics |
| Grouped CV needs MMseqs2 | make it optional and recorded |

## 7. Independent review (Opus 5.5, 2026-09-30)

Seven findings; all are applied above.

| # | Severity | Finding | Applied in |
|---|---|---|---|
| 1 | High | Pooling changed after the pickles were made; legacy models get wrong input today | fact 4, 3.1, 3.3, step 3 |
| 2 | High | Step 1 was not "no code change" and used the wrong features | step 1 |
| 3 | Medium | Predict must take layer and model from the card; truncation and model path were unhandled | 3.3 |
| 4 | Medium | Version drift includes numpy and torch | 3.2, risks |
| 5 | Low | Fixed `C=1` is defensible but untested | step 1, 3.4 |
| 6 | Low | T0 acceptance was ambiguous | section 5 |
| 7 | Low | Step order encoded the wrong legacy default | steps 2-3 |

The reviewer confirmed the card, bare-LR and pickle claims against `train.py`, `model.py` and
`predict.py`, and loaded both shipped pickles (C=1.0, shapes (1,320) and (1,480)). It corrected
one statement: the 0.995 figure is for `aa_comp+length`, not composition alone.

## 8. Open items raised after review

| id | item | state |
|---|---|---|
| R1 | Rename CLI, column and label from adhesion to surface glycoprotein (section 0) | proposed in TOOL-ARCHITECTURE 2.0; waits for design review |
| R2 | Pooling mismatch for shipped pickles: measure and document | job 29301182 (step 1); issue filed |
| R3 | Agreement between the ESM + LR step-1 CLI and a SignalP + GPI call | not measured |
| R4 | `src/adhesion_predict/models/adhesion_model_esm2_t6_8M_UR50D.pkl` is untracked although `pyproject.toml` packages `models/*` | open; identical to `models/` copy (md5 checked) |
