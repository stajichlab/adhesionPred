# Plan for #9: model bundle, model card, explicit embedding config

*Drafted 2026-09-30. Status: proposal, not implemented. Review requested before any code.*

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
   for amino-acid composition alone (leave-family-out). Scaling or tuning `C` cannot add much
   there. The value of #9 is reproducibility and safety, not score. This is a prediction, not a
   measurement. The plan includes the measurement (step 2).
3. **Layer 6 means different things.** For `esm2_t6_8M` it is the final layer. For
   `esm2_t12_35M` (used for the Fungi_5k screen) it is the middle layer. Changing the default
   would silently change every existing embedding. The default must stay 6 for existing models.
4. **Shipped models have no card.** `models/*.pkl` and `src/adhesion_predict/models/*.pkl` have
   no JSON. They must keep working.
5. **The shipped models were trained with 167 duplicate sequences** (PR #24 removes them for
   new training). A retrained bundle is therefore not comparable to the old pickles on the
   training data, only on an external set.
6. `config.ESM2_MODELS` lists `esm2_t6_35M_UR50D` and `esm2_t6_150M_UR50D`. I do not know that
   these ESM-2 checkpoints exist. Only two names are in `ESM2_MODEL_CHOICES`. Verify or delete
   before writing them into any card or test.
7. Only 0.5 thresholding exists. The review (§5.4, §6 T5) says the threshold must be chosen on
   proteome-level data, which needs the T5 inventory (#12, #13, #14). #9 cannot choose it.

## 3. Design

### 3.1 Bundle format

One directory or one file per model, `adhesion_model_<esm>.joblib`, plus `.json` card beside it
(same layout as now, so `model_card_path` is unchanged). The joblib object is a
`sklearn.pipeline.Pipeline([("scale", StandardScaler()), ("clf", LogisticRegression(...))])`.

Loading rules:
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

The card must not claim a tier that was not run. An empty `validation` list is the honest value.

### 3.3 `predict` compatibility check

Compare `(esm_model, repr_layer, pooling, max_residues)` in the card against what this
installation will compute. On mismatch, exit with an error that names the field. Add
`--repr-layer` to `predict` and `train`, defaulting to the card value for `predict`. No card:
warn, assume defaults, continue.

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

Each step is one commit; steps 1-3 need no GPU.

1. **Measure the gap (no code change).** Re-run `02_cv_and_proteome_eval.py` twice on the same
   embeddings: scaled (as reviewed) and unscaled (as shipped). Report both. If the unscaled model
   is materially worse on leave-family-out, the plan matters more than fact 2 predicts. If not,
   say so and keep the Pipeline for robustness only.
   *Needs the cached `train_masked.npy`; `run.sh` makes it in about 1.5 h on CPU.*
2. **Tests first (TDD).** Card round trip with v2 fields; legacy model loads with a warning;
   mismatched layer is refused; a Pipeline bundle round-trips through joblib.
3. **Card v2 + `load_model` for both formats + predict checks.** No change to training output yet.
4. **Pipeline in `train_classifier`; `--repr-layer`, `--seed`.**
5. **Optional `--tune-c` with grouped nested CV.**
6. **Retrain both bundles from deduplicated data; compare with the old pickles on a fixed
   external set (S288C, `analysis/model_review/s288c_scores.tsv`) before replacing them.**
   Keep the old files under a `legacy/` name for one release.
7. Update README, `docs/TOOL-ARCHITECTURE.md` status, and AGENTS.md.

## 5. Acceptance

- `pytest tests/adhesion_predict` passes, with new tests for each step-2 case.
- `predict` on a v2 bundle and on each legacy pickle both run; a wrong-layer call exits non-zero.
- The card of a newly trained bundle lists every field in 3.2 and sets no value it did not compute.
- Step 1 and step 6 comparisons are in a report with numbers, including the case where nothing changed.
- Scores for a fixed 100-sequence set are identical between the old pickle path and the new
  bundle path when both use the same scaler setting (T0 determinism).

## 6. Risks

| risk | handling |
|---|---|
| Retraining changes calls that other work depends on (Fungi_5k survey used the old model) | keep `legacy/` models; record model hash in any survey output |
| Scaler changes score scale, so the 0.5 threshold means something different | compare call counts on S288C before and after; do not ship without that |
| joblib/sklearn version drift breaks loading | store `sklearn_version`; warn on mismatch at load |
| Card claims more than was measured | schema allows empty `validation`; no default metrics |
| Grouped CV needs MMseqs2 | make it optional and recorded |
