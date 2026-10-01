# surface_glyco rename and legacy safety fixes: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rename the frozen legacy ESM-2 + logistic regression package from `adhesion_predict` to `surface_glyco` (breaking, no aliases), and add the card, pooling and refusal checks that stop it from scoring with the wrong embedding settings.

**Architecture:** One mechanical rename task first, so every later task works in the new namespace. Then small test-first changes: output schema, model location, shared card/settings logic (`card.py`, no torch import), a legacy pooling mode in `embeddings.py`, card-driven `predict`/`evaluate`/`train`, and a shared reader that lets downstream analysis code read both old and new result files. Release is a minor bump to 0.2.0 through the existing `version_bump` workflow.

**Tech Stack:** Python 3.12, setuptools + setuptools-scm, fair-esm (ESM-2 8M for tests), scikit-learn, pytest, ruff 0.3.5 (pre-commit), GitHub Actions.

**Spec:** `docs/PLAN-2026-09-30-pipeline-and-decisions.md` (decisions 2, 4, 5 and defaults 1-6) and `docs/plans/2026-09-30-design-review-fable.md` (register items 1-4, 6, 14, 26, 27, 30-33).

## Global Constraints

- **Names, exactly:** package `surface_glyco`; entry points `surface_glyco_predict`, `surface_glyco_train`, `surface_glyco_evaluate`; output column `surface_glycoprotein_score`; labels `surface_glycoprotein` and `other`; model files `surface_glyco_model_<esm>.pkl`. Repo name `adhesionPred` stays.
- **Breaking rename, no aliases** for the package, entry points, column and labels (decision 4). Frozen legacy result files are the one exception: downstream readers must still read them (Task 4).
- **No retraining, no tuning, no scaler pipeline** (decision 2). The two shipped pickles are renamed and given cards, nothing else.
- **Default model** stays `esm2_t6_8M_UR50D`. Default `repr_layer` stays 6. `MAX_RESIDUES` stays 1022. Threshold stays 0.5 and is documented as uncalibrated.
- **Legacy pooling cannot be reproduced exactly** (batch composition was not recorded). It is emulated only behind an explicit flag, never by default.
- **Python 3.12 is canonical; CI is CPU only.** Lint is ruff 0.3.5 via pre-commit, line length 100. Run `pre-commit run --all-files` before each commit.
- **Version:** the version is derived from git tags by setuptools-scm. 0.2.0 is produced by the `version_bump` workflow with `bump_type=minor`, not by editing a file.
- **Do not edit historical documents:** `docs/superpowers/plans/2026-09-04-*`, `docs/superpowers/specs/2026-09-04-*`, `Changes.md`, `docs/model-review/*`, and `analysis/chytrid_batrach/*` (one-off analyses of frozen legacy result files) keep the old names.
- Each commit message ends with: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
- Run tests with `PYTHONPATH=$PWD/src python -m pytest ...` only until Task 1 step 7 reinstalls the package; after that use plain `pytest` in an environment where `pip install -e .` was run on this checkout.

## Review Focus

1. **Frozen legacy result files** (calls only, header `id,prediction,probability_adhesion`, label `Adhesion`) must still be read correctly by `kingdom_survey` and `adhesion_properties` after the rename. Pinned in Task 4.
2. **New all-rows result files** (since PR #23 every protein is written, not only calls): `kingdom_survey.join.count_result_rows` currently counts every data row as a call, so it overcounts. Pinned in Task 4.
3. **A `./models` directory in the current working directory** must not silently replace the packaged model. Pinned in Task 3.
4. **A model card that is missing, has an unknown `pooling`, or names a different ESM model** must stop `predict`/`evaluate` with an error that names the field. Pinned in Task 6.
5. **Sequences that are empty, longer than 1022 residues, or present in both classes** must not corrupt counts or labels. Pinned in Tasks 5 and 7.
6. **Installed-name mismatch:** importing the package from a checkout where it is not installed under the new name must not raise `PackageNotFoundError`; a built wheel must contain the model files. Pinned in Tasks 1 and 9.

---

## File structure

| File | Responsibility | Task |
|---|---|---|
| `src/surface_glyco/` (renamed from `src/adhesion_predict/`) | the whole legacy package | 1 |
| `src/surface_glyco/card.py` (new) | card constants, `EmbeddingSettings`, `resolve_embedding_settings`, `new_card`; no torch import | 5 |
| `src/surface_glyco/results.py` (new) | read result CSVs of either schema | 4 |
| `src/surface_glyco/embeddings.py` | add `pooling` argument, `count_truncated` | 6 |
| `src/surface_glyco/model.py` | card load/save (unchanged API) | 5 |
| `src/surface_glyco/scripts/{predict,train,evaluate}.py` | schema, card use, flags | 2, 7 |
| `src/surface_glyco/config.py` | model directory, delete `ESM2_MODELS` | 3 |
| `src/surface_glyco/models/surface_glyco_model_*.pkl` + `.json` | shipped legacy models and their cards | 3, 5 |
| `tests/surface_glyco/` (renamed from `tests/adhesion_predict/`) | package tests | 1-7 |
| `analysis/kingdom_survey/join.py`, `analysis/adhesion_properties/universe.py` | use `results.read_called` | 4 |
| `pyproject.toml`, `.coveragerc`, `CITATION.cff`, `.github/workflows/build_and_test.yml` | names | 1 |
| `CHANGELOG.md`, `README.md`, `AGENTS.md`, `docs/TOOL-ARCHITECTURE.md` | release and docs | 8 |

---

### Task 1: Rename the package, entry points and test folder

**Files:**
- Rename: `src/adhesion_predict/` to `src/surface_glyco/`; `tests/adhesion_predict/` to `tests/surface_glyco/`
- Modify: `pyproject.toml:10,37-39,73`; `.coveragerc:3`; `CITATION.cff:5`; `.github/workflows/build_and_test.yml` (2 lines); `src/surface_glyco/__init__.py`; every `adhesion_predict` import in `src/` and `tests/surface_glyco/`; `analysis/embedding_clustering/{00_feasibility_check,embed,recover_missing_esm2_embeddings}.py`; `tests/embedding_clustering/*` if they import it
- Test: `tests/surface_glyco/test_package.py` (new)

**Interfaces:**
- Produces: importable package `surface_glyco`; `surface_glyco.__version__` (str) that never raises; console scripts `surface_glyco_predict|train|evaluate` pointing at `surface_glyco.scripts.<name>:cli`.

- [ ] **Step 1: Write the failing test**

Create `tests/surface_glyco/test_package.py` (the folder does not exist until step 3, so create it in step 3 first; the test text is final):

```python
"""Package identity: names, version fallback, entry points."""

import importlib
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_new_package_imports_and_has_a_version_string():
    import surface_glyco

    assert isinstance(surface_glyco.__version__, str) and surface_glyco.__version__


def test_old_package_name_is_gone():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("adhesion_predict")


def test_pyproject_names_match_decision_4():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert cfg["project"]["name"] == "surface_glyco"
    assert cfg["project"]["scripts"] == {
        "surface_glyco_predict": "surface_glyco.scripts.predict:cli",
        "surface_glyco_train": "surface_glyco.scripts.train:cli",
        "surface_glyco_evaluate": "surface_glyco.scripts.evaluate:cli",
    }
    assert cfg["tool"]["setuptools"]["package-data"] == {"surface_glyco": ["models/*"]}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd <repo> && PYTHONPATH=$PWD/src python -m pytest tests/adhesion_predict -q` first to record the baseline (12 passed), then the new test after step 3 creates the folder.
Expected after step 3 and before step 4: FAIL `ModuleNotFoundError: No module named 'surface_glyco'`.

- [ ] **Step 3: Move the folders**

```bash
git mv src/adhesion_predict src/surface_glyco
git mv tests/adhesion_predict tests/surface_glyco
```
Then create `tests/surface_glyco/test_package.py` with the text from step 1.

- [ ] **Step 4: Rewrite names**

```bash
grep -rl "adhesion_predict" src tests/surface_glyco analysis/embedding_clustering tests/embedding_clustering .coveragerc CITATION.cff .github/workflows/build_and_test.yml \
  | xargs sed -i 's/adhesion_predict/surface_glyco/g'
```
Then edit by hand, because the sed above does not know the entry-point names:

`pyproject.toml`:
```toml
[project]
name = "surface_glyco"
description = "Legacy ESM-2 + logistic regression scorer for FLO/ALS-like surface glycoproteins"

[project.scripts]
surface_glyco_predict = "surface_glyco.scripts.predict:cli"
surface_glyco_train = "surface_glyco.scripts.train:cli"
surface_glyco_evaluate = "surface_glyco.scripts.evaluate:cli"

[tool.setuptools.package-data]
"surface_glyco" = ["models/*"]
```
Also update the comment at the top of the file that mentions the `adhesion_` script prefix, and in `.github/workflows/build_and_test.yml` change `pytest tests/adhesion_predict` to `pytest tests/surface_glyco` and the `--ignore=` path if present.

Replace the body of `src/surface_glyco/__init__.py` header so the version never raises:

```python
"""surface_glyco - legacy ESM-2 + logistic regression scorer for surface glycoproteins."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("surface_glyco")
except PackageNotFoundError:  # running from a checkout that was not pip-installed
    __version__ = "0+unknown"
```
Keep the existing `from surface_glyco.io import ...` and `__all__` lines below it.

- [ ] **Step 5: Run tests**

Run: `PYTHONPATH=$PWD/src python -m pytest tests/surface_glyco -q -p no:cacheprovider`
Expected: 15 passed (12 existing + 3 new). `test_old_package_name_is_gone` passes only if no stale `adhesion_predict` is on `PYTHONPATH`; unset any editable install of the old name in the shell used for the run.

- [ ] **Step 6: Verify no stale references in code**

Run: `git grep -n "adhesion_predict\|adhesion_train\|adhesion_evaluate" -- src tests analysis/embedding_clustering analysis/kingdom_survey analysis/adhesion_properties pyproject.toml .coveragerc CITATION.cff .github`
Expected: only `analysis/model_review/02_cv_and_proteome_eval.py` (handled in Task 3) and comment/doc strings you decide to keep. Fix any other hit.

- [ ] **Step 7: Reinstall and run the entry points**

Run: `pip install -e . --no-deps && surface_glyco_predict --help | head -3 && python -c "import surface_glyco; print(surface_glyco.__version__)"`
Expected: help text prints; a version string prints. If `adhesion_predict` is still installed in this environment run `pip uninstall -y adhesion_predict` first.

- [ ] **Step 8: Commit**

```bash
pre-commit run --all-files
git add -A
git commit -m "rename: package adhesion_predict to surface_glyco, with new entry points

Breaking, no aliases. __version__ falls back to 0+unknown when the package
is not installed under the new name.

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Output schema and wording

**Files:**
- Modify: `src/surface_glyco/scripts/predict.py` (lines for `write_results`, `build_results`, banner, help text, default output name), `src/surface_glyco/scripts/evaluate.py` (`target_names`, description), `src/surface_glyco/scripts/train.py` (banner, description)
- Test: `tests/surface_glyco/test_io_model_predict.py`

**Interfaces:**
- Produces: `build_results(seq_ids, predictions, probabilities) -> list[dict]` with keys `id`, `prediction` (`"surface_glycoprotein"` or `"other"`), `surface_glycoprotein_score` (float). `write_results(results, path)` writes header `id,prediction,surface_glycoprotein_score`. Constants `LABEL_POSITIVE = "surface_glycoprotein"`, `LABEL_NEGATIVE = "other"`, `SCORE_COLUMN = "surface_glycoprotein_score"` exported from `surface_glyco.scripts.predict`. Default output file name `<stem>.surface_glyco.csv`.

- [ ] **Step 1: Write the failing tests** (replace the two existing tests that use the old names, and add the third)

```python
def test_write_results_header_and_quoting(tmp_path):
    out = tmp_path / "r.csv"
    write_results(
        [{"id": "sp|P1,x", "prediction": "surface_glycoprotein", "surface_glycoprotein_score": 0.91234}],
        out,
    )
    assert out.read_text().splitlines()[0] == "id,prediction,surface_glycoprotein_score"
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows == [
        {"id": "sp|P1,x", "prediction": "surface_glycoprotein", "surface_glycoprotein_score": "0.9123"}
    ]


def test_build_results_keeps_every_row_with_new_labels():
    rows = build_results(["a", "b"], [1, 0], [[0.1, 0.9], [0.8, 0.2]])
    assert [r["prediction"] for r in rows] == ["surface_glycoprotein", "other"]
    assert rows[1]["surface_glycoprotein_score"] == 0.2
    assert "probability_adhesion" not in rows[0]


def test_default_output_name_uses_new_suffix(tmp_path, monkeypatch):
    from surface_glyco.scripts.predict import default_output_path

    assert default_output_path(tmp_path / "Scer.pep.fa").name == "Scer.pep.surface_glyco.csv"
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_io_model_predict.py -q`
Expected: FAIL (old header, old labels, `default_output_path` missing).

- [ ] **Step 3: Implement**

In `predict.py` add near the top:

```python
LABEL_POSITIVE = "surface_glycoprotein"
LABEL_NEGATIVE = "other"
SCORE_COLUMN = "surface_glycoprotein_score"


def default_output_path(input_path):
    """Output CSV in the current directory, named after the input."""
    input_path = Path(input_path)
    stem = input_path.name if input_path.is_dir() else input_path.stem
    return Path.cwd() / f"{stem}.surface_glyco.csv"
```
Change `write_results` to write `["id", "prediction", SCORE_COLUMN]` and read `result[SCORE_COLUMN]`. Change `build_results` to emit `LABEL_POSITIVE if predictions[i] == 1 else LABEL_NEGATIVE` and key `SCORE_COLUMN`. In `main`, replace the `if output_file is None:` block with `output_file = default_output_path(input_path)`; replace the stdout filter `result["prediction"] == "Adhesion"` with `== LABEL_POSITIVE`, and the printed `p=` field with `result[SCORE_COLUMN]`. Banner: `"Surface glycoprotein scoring (legacy ESM-2 + logistic regression)"`. Parser description: `"Score proteins as FLO/ALS-like surface glycoproteins (legacy model)"`. Error hint: `"Run training first: surface_glyco_train"`. In `evaluate.py` set `target_names = ["other", "surface_glycoprotein"]` and description to `"Evaluate the legacy surface glycoprotein model"`. In `train.py` change the banner to `"Surface glycoprotein classifier - legacy training"` and the description to match.

- [ ] **Step 4: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "schema: surface_glycoprotein_score column and surface_glycoprotein/other labels

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Model files and model directory

**Files:**
- Move: `models/adhesion_model_esm2_t6_8M_UR50D.pkl` and `models/adhesion_model_esm2_t12_35M_UR50D.pkl` to `src/surface_glyco/models/surface_glyco_model_<esm>.pkl` (the `src/` copies are byte-identical; keep one set)
- Modify: `src/surface_glyco/config.py` (directory lookup, delete `ESM2_MODELS`), the three scripts (default model path), `analysis/model_review/02_cv_and_proteome_eval.py:95`
- Delete: root `models/` after the move
- Test: `tests/surface_glyco/test_config_paths.py` (new)

**Interfaces:**
- Produces: `surface_glyco.config.get_models_dir() -> Path` (packaged directory first; the current-directory `models/` is used only when `SURFACE_GLYCO_MODELS_DIR` is set); `surface_glyco.config.model_filename(esm_model: str) -> str` returning `f"surface_glyco_model_{esm_model}.pkl"`.

- [ ] **Step 1: Write the failing tests**

```python
"""Model directory lookup and file names."""

import hashlib
from pathlib import Path

from surface_glyco import config

PKG_MODELS = Path(config.__file__).parent / "models"


def test_packaged_models_exist_under_new_names():
    for esm in ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D"):
        assert (PKG_MODELS / config.model_filename(esm)).is_file()


def test_stray_models_dir_in_cwd_does_not_override(tmp_path, monkeypatch):
    (tmp_path / "models").mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("SURFACE_GLYCO_MODELS_DIR", raising=False)
    assert config.get_models_dir() == PKG_MODELS


def test_env_var_overrides(tmp_path, monkeypatch):
    monkeypatch.setenv("SURFACE_GLYCO_MODELS_DIR", str(tmp_path))
    assert config.get_models_dir() == tmp_path


def test_shipped_pickles_are_unchanged_bytes():
    # md5 of the pickles as committed on 2026-02-14; the rename must not alter them.
    expected = {
        "esm2_t6_8M_UR50D": "ce6dc8a40809d6afb8269ee3280c9910",
        "esm2_t12_35M_UR50D": "5f7a08cc0122e3e43bc18c1539f77d7f",
    }
    for esm, md5 in expected.items():
        data = (PKG_MODELS / config.model_filename(esm)).read_bytes()
        assert hashlib.md5(data).hexdigest() == md5
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_config_paths.py -q`
Expected: FAIL (`model_filename` missing, files not moved).

- [ ] **Step 3: Move the pickles**

```bash
git mv src/surface_glyco/models/adhesion_model_esm2_t12_35M_UR50D.pkl src/surface_glyco/models/surface_glyco_model_esm2_t12_35M_UR50D.pkl
git add src/surface_glyco/models/ # the t6_8M copy is untracked on some checkouts; add it from models/
cp models/adhesion_model_esm2_t6_8M_UR50D.pkl src/surface_glyco/models/surface_glyco_model_esm2_t6_8M_UR50D.pkl
git add src/surface_glyco/models/surface_glyco_model_esm2_t6_8M_UR50D.pkl
git rm -r models/
```

- [ ] **Step 4: Implement `config.py`**

Replace `get_models_dir` and delete `ESM2_MODELS`:

```python
import os

PACKAGED_MODELS_DIR = Path(__file__).parent / "models"


def get_models_dir():
    """Directory holding trained models.

    The packaged directory is used unless SURFACE_GLYCO_MODELS_DIR is set. A ./models
    directory in the current working directory is deliberately ignored, so a stray
    folder cannot replace the shipped model.
    """
    override = os.environ.get("SURFACE_GLYCO_MODELS_DIR")
    return Path(override) if override else PACKAGED_MODELS_DIR


def model_filename(esm_model):
    """File name of the model trained on embeddings from esm_model."""
    return f"surface_glyco_model_{esm_model}.pkl"
```
Keep `MODELS_DIR = get_models_dir()`. In `predict.py`, `train.py` and `evaluate.py` replace `f"adhesion_model_{...}.pkl"` with `model_filename(...)` (import it from `surface_glyco.config`). In `analysis/model_review/02_cv_and_proteome_eval.py` change the path to `str(REPO) + "/src/surface_glyco/models/surface_glyco_model_esm2_t6_8M_UR50D.pkl"`.

- [ ] **Step 5: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS. If `test_shipped_pickles_are_unchanged_bytes` fails, stop: the pickle bytes changed, which violates the "no retraining" constraint.

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "models: ship both legacy pickles under surface_glyco names; ignore ./models in cwd

Packaged directory is used unless SURFACE_GLYCO_MODELS_DIR is set. Root
models/ removed (byte-identical copy). The t6_8M pickle is now tracked in
the package folder. ESM2_MODELS removed (listed models that do not exist).

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Read both result-file schemas downstream

**Files:**
- Create: `src/surface_glyco/results.py`
- Modify: `analysis/kingdom_survey/join.py` (`count_result_rows`, `probability_stats`), `analysis/adhesion_properties/universe.py` (`read_result_ids`)
- Test: `tests/surface_glyco/test_results.py` (new); add cases to `tests/kingdom_survey/test_join.py` and `tests/adhesion_properties/test_universe.py`

**Interfaces:**
- Produces: `surface_glyco.results.read_called(path) -> list[tuple[str, float]]`: `(id, score)` for rows whose label is called. Accepts the legacy header (`probability_adhesion`, label `Adhesion`) and the new header (`surface_glycoprotein_score`, label `surface_glycoprotein`). Rows labelled `Non-adhesion` or `other` are skipped. A file with neither score column raises `ValueError` naming the file. `read_all(path) -> list[tuple[str, str, float]]` returns every row as `(id, label, score)`.
- Consumes: nothing from earlier tasks except the package name.

- [ ] **Step 1: Write the failing tests**

```python
"""Reading result CSVs of the legacy and the current schema."""

import pytest

from surface_glyco.results import read_all, read_called

LEGACY = "id,prediction,probability_adhesion\nA,Adhesion,0.9\nB,Adhesion,0.7\n"
NEW = (
    "id,prediction,surface_glycoprotein_score\n"
    "A,surface_glycoprotein,0.9\nB,other,0.1\nC,surface_glycoprotein,0.6\n"
)
LEGACY_ALL_ROWS = (  # a legacy-header file written after PR #23, with non-calls
    "id,prediction,probability_adhesion\nA,Adhesion,0.9\nB,Non-adhesion,0.1\n"
)


def _write(tmp_path, text):
    p = tmp_path / "r.csv"
    p.write_text(text)
    return p


def test_legacy_calls_only_file(tmp_path):
    assert read_called(_write(tmp_path, LEGACY)) == [("A", 0.9), ("B", 0.7)]


def test_new_schema_returns_only_called_rows(tmp_path):
    assert read_called(_write(tmp_path, NEW)) == [("A", 0.9), ("C", 0.6)]


def test_legacy_header_with_non_calls_skips_them(tmp_path):
    assert read_called(_write(tmp_path, LEGACY_ALL_ROWS)) == [("A", 0.9)]


def test_header_only_file_is_empty_not_an_error(tmp_path):
    assert read_called(_write(tmp_path, "id,prediction,surface_glycoprotein_score\n")) == []


def test_ids_with_commas_survive(tmp_path):
    text = 'id,prediction,surface_glycoprotein_score\n"sp|P1,x",surface_glycoprotein,0.5\n'
    assert read_called(_write(tmp_path, text)) == [("sp|P1,x", 0.5)]


def test_unknown_header_raises_naming_the_file(tmp_path):
    p = _write(tmp_path, "id,prediction,score\nA,x,0.5\n")
    with pytest.raises(ValueError, match="r.csv"):
        read_called(p)


def test_read_all_returns_every_row(tmp_path):
    assert read_all(_write(tmp_path, NEW))[1] == ("B", "other", 0.1)
```

And add to `tests/kingdom_survey/test_join.py`:

```python
def test_count_result_rows_counts_calls_not_all_rows(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text(
        "id,prediction,surface_glycoprotein_score\n"
        "A,surface_glycoprotein,0.9\nB,other,0.1\nC,other,0.2\n"
    )
    assert join.count_result_rows(p) == 1
    assert join.probability_stats(p) == (0.9, 0.9)
```
(Use the module alias that file already uses for `join`.) Add one case to `tests/adhesion_properties/test_universe.py` calling `universe.read_result_ids` on the new-schema file above and asserting `[("A", 0.9), ("C", 0.6)]` for a three-row version.

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_results.py -q`
Expected: FAIL (`ModuleNotFoundError: surface_glyco.results`).

- [ ] **Step 3: Implement `results.py`**

```python
"""Read prediction CSVs written by surface_glyco_predict, current or legacy schema.

Legacy files (adhesion_predict, before 0.2.0) use the column probability_adhesion and the
labels Adhesion / Non-adhesion. They are frozen outputs on disk and must stay readable.
"""

import csv
from pathlib import Path

SCORE_COLUMNS = ("surface_glycoprotein_score", "probability_adhesion")
CALLED_LABELS = ("surface_glycoprotein", "Adhesion")


def _score_column(fieldnames, path):
    for name in SCORE_COLUMNS:
        if fieldnames and name in fieldnames:
            return name
    raise ValueError(f"{path}: no score column; expected one of {SCORE_COLUMNS}")


def read_all(path):
    """Return (id, label, score) for every row."""
    path = Path(path)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        column = _score_column(reader.fieldnames, path)
        return [(row["id"], row["prediction"], float(row[column])) for row in reader]


def read_called(path):
    """Return (id, score) for rows called as surface glycoproteins."""
    return [(i, score) for i, label, score in read_all(path) if label in CALLED_LABELS]
```
Edge case: a header-only file has `fieldnames` set, so it returns `[]` (test 4). A completely empty file has `fieldnames is None` and raises `ValueError`; that is intended.

- [ ] **Step 4: Use it downstream**

`analysis/kingdom_survey/join.py`: add `from surface_glyco.results import read_called`; make `count_result_rows` return `len(read_called(result_csv_path))` and `probability_stats` compute mean and median of `[s for _, s in read_called(result_csv_path)]` (keep the `nan, nan` return for empty). Update both docstrings to say "called proteins". `analysis/adhesion_properties/universe.py`: `read_result_ids` returns `read_called(result_csv_path)` and its docstring says "(protein_id, score) for called proteins".

- [ ] **Step 5: Run all affected suites separately** (each folder's `conftest.py` edits `sys.path`)

```bash
pytest tests/surface_glyco -q
pytest tests/kingdom_survey -q
pytest tests/adhesion_properties -q
```
Expected: all PASS. The existing legacy-header fixtures in those two folders must pass unchanged; that is the proof that frozen files still read.

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "results: one reader for legacy and current result CSVs

kingdom_survey counted every data row as a call, which overcounts since all
scores are written (PR #23). It now counts called rows only, for both schemas.

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Card v2, settings resolution, and cards for the shipped pickles

**Files:**
- Create: `src/surface_glyco/card.py`; `src/surface_glyco/models/surface_glyco_model_esm2_t6_8M_UR50D.json`, `...t12_35M_UR50D.json`
- Modify: `src/surface_glyco/embeddings.py` (import constants from `card.py`)
- Test: `tests/surface_glyco/test_card.py` (new)

**Interfaces:**
- Produces, all in `surface_glyco.card` (no torch import):
  - Constants `CARD_VERSION = 2`, `POOLING_RESIDUE_MEAN = "residue_mean"`, `POOLING_LEGACY = "all_tokens_legacy"`, `KNOWN_POOLINGS`, `MAX_RESIDUES = 1022`, `DEFAULT_REPR_LAYER = 6`.
  - `class ModelCardError(Exception)`.
  - `@dataclass(frozen=True) class EmbeddingSettings: esm_model: str; repr_layer: int; pooling: str; max_residues: int; implied: bool`.
  - `resolve_embedding_settings(card, cli_model_name, default_model, allow_legacy_pooling=False) -> EmbeddingSettings`.
  - `new_card(esm_model, repr_layer, pooling, n_positive, n_negative, **extra) -> dict` (v2 layout).

- [ ] **Step 1: Write the failing tests**

```python
"""Card contents and settings resolution."""

import pytest

from surface_glyco.card import (
    MAX_RESIDUES,
    POOLING_LEGACY,
    POOLING_RESIDUE_MEAN,
    ModelCardError,
    new_card,
    resolve_embedding_settings,
)

DEFAULT = "esm2_t6_8M_UR50D"


def _card(**over):
    card = new_card(DEFAULT, 6, POOLING_RESIDUE_MEAN, n_positive=10, n_negative=20)
    card.update(over)
    return card


def test_new_card_has_v2_fields():
    card = _card()
    assert card["card_version"] == 2
    assert card["pooling"] == POOLING_RESIDUE_MEAN
    assert card["max_residues"] == MAX_RESIDUES
    assert card["truncation"] == "truncate_at_max_residues"
    assert card["validation"] == []  # nothing claimed that was not run
    assert card["threshold"] == {"value": 0.5, "chosen_on": "default (uncalibrated)"}


def test_card_settings_are_used_when_cli_model_is_omitted():
    s = resolve_embedding_settings(_card(esm_model="esm2_t12_35M_UR50D"), None, DEFAULT)
    assert (s.esm_model, s.repr_layer, s.pooling, s.implied) == (
        "esm2_t12_35M_UR50D", 6, POOLING_RESIDUE_MEAN, False)


def test_cli_model_that_disagrees_with_card_names_the_field():
    with pytest.raises(ModelCardError, match="esm_model"):
        resolve_embedding_settings(_card(), "esm2_t12_35M_UR50D", DEFAULT)


def test_unknown_pooling_is_refused_naming_the_field():
    with pytest.raises(ModelCardError, match="pooling"):
        resolve_embedding_settings(_card(pooling="max"), None, DEFAULT)


def test_unsupported_max_residues_is_refused():
    with pytest.raises(ModelCardError, match="max_residues"):
        resolve_embedding_settings(_card(max_residues=2000), None, DEFAULT)


def test_legacy_pooling_needs_the_flag():
    card = _card(pooling=POOLING_LEGACY)
    with pytest.raises(ModelCardError, match="allow-legacy-pooling"):
        resolve_embedding_settings(card, None, DEFAULT)
    s = resolve_embedding_settings(card, None, DEFAULT, allow_legacy_pooling=True)
    assert s.pooling == POOLING_LEGACY


def test_missing_card_is_implied_legacy_and_needs_the_flag():
    with pytest.raises(ModelCardError, match="allow-legacy-pooling"):
        resolve_embedding_settings(None, None, DEFAULT)
    s = resolve_embedding_settings(None, None, DEFAULT, allow_legacy_pooling=True)
    assert (s.esm_model, s.repr_layer, s.pooling, s.implied) == (DEFAULT, 6, POOLING_LEGACY, True)


def test_shipped_cards_describe_the_legacy_training(tmp_path):
    import json
    from pathlib import Path

    from surface_glyco import config

    for esm in ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D"):
        card = json.loads(
            (Path(config.__file__).parent / "models" / f"surface_glyco_model_{esm}.json").read_text()
        )
        assert card["esm_model"] == esm
        assert card["pooling"] == POOLING_LEGACY
        assert card["legacy"] is True
        assert "batch" in card["pooling_note"]
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_card.py -q`
Expected: FAIL (`ModuleNotFoundError: surface_glyco.card`).

- [ ] **Step 3: Implement `card.py`**

```python
"""Model card contents and the rules for turning a card into embedding settings.

This module must not import torch or esm: predict/evaluate validate a model before
loading any weights, and tests of these rules run without a GPU stack.
"""

from dataclasses import dataclass

CARD_VERSION = 2
POOLING_RESIDUE_MEAN = "residue_mean"
POOLING_LEGACY = "all_tokens_legacy"
KNOWN_POOLINGS = (POOLING_RESIDUE_MEAN, POOLING_LEGACY)
MAX_RESIDUES = 1022  # ESM-2 context is 1024 tokens including BOS and EOS
DEFAULT_REPR_LAYER = 6


class ModelCardError(Exception):
    """The model card does not allow scoring with the requested settings."""


@dataclass(frozen=True)
class EmbeddingSettings:
    esm_model: str
    repr_layer: int
    pooling: str
    max_residues: int
    implied: bool  # True when there was no card and the settings are assumed


def new_card(esm_model, repr_layer, pooling, n_positive, n_negative, **extra):
    """Card for a newly trained model. Fields not measured are left empty, not guessed."""
    card = {
        "card_version": CARD_VERSION,
        "esm_model": esm_model,
        "repr_layer": repr_layer,
        "pooling": pooling,
        "max_residues": MAX_RESIDUES,
        "truncation": "truncate_at_max_residues",
        "n_positive": n_positive,
        "n_negative": n_negative,
        "validation": [],
        "threshold": {"value": 0.5, "chosen_on": "default (uncalibrated)"},
    }
    card.update(extra)
    return card


def _refuse_legacy(why):
    raise ModelCardError(
        f"{why} Legacy pooling averages padding tokens, so scores depend on batch composition "
        "and cannot be reproduced exactly. Re-run with --allow-legacy-pooling to emulate it "
        "(batch size 4, input order), or use a model trained with residue_mean pooling."
    )


def resolve_embedding_settings(card, cli_model_name, default_model, allow_legacy_pooling=False):
    """Return the embedding settings a model must be scored with, or raise ModelCardError."""
    if card is None:
        if not allow_legacy_pooling:
            _refuse_legacy("The model has no card, so its embedding settings are unknown.")
        return EmbeddingSettings(
            cli_model_name or default_model, DEFAULT_REPR_LAYER, POOLING_LEGACY, MAX_RESIDUES, True
        )

    esm_model = card.get("esm_model")
    if cli_model_name is not None and cli_model_name != esm_model:
        raise ModelCardError(
            f"esm_model: the model was trained on {esm_model} embeddings, "
            f"but --model-name is {cli_model_name}"
        )
    pooling = card.get("pooling", POOLING_LEGACY)
    if pooling not in KNOWN_POOLINGS:
        raise ModelCardError(f"pooling: unknown value {pooling!r}; known: {KNOWN_POOLINGS}")
    max_residues = card.get("max_residues", MAX_RESIDUES)
    if max_residues != MAX_RESIDUES:
        raise ModelCardError(
            f"max_residues: the model used {max_residues}, this build supports {MAX_RESIDUES}"
        )
    if pooling == POOLING_LEGACY and not allow_legacy_pooling:
        _refuse_legacy("pooling: the model was trained with all_tokens_legacy pooling.")
    return EmbeddingSettings(
        esm_model or default_model,
        int(card.get("repr_layer", DEFAULT_REPR_LAYER)),
        pooling,
        max_residues,
        False,
    )
```
In `embeddings.py` replace the local `DEFAULT_REPR_LAYER`, `MAX_RESIDUES` and `POOLING` definitions with `from surface_glyco.card import DEFAULT_REPR_LAYER, MAX_RESIDUES, POOLING_LEGACY, POOLING_RESIDUE_MEAN, KNOWN_POOLINGS` and keep `POOLING = POOLING_RESIDUE_MEAN` for `train.py`'s import until Task 7 removes it.

- [ ] **Step 4: Write the two shipped cards**

`surface_glyco_model_esm2_t6_8M_UR50D.json` (same with `t12_35M` and its name for the second file):

```json
{
  "card_version": 2,
  "legacy": true,
  "esm_model": "esm2_t6_8M_UR50D",
  "repr_layer": 6,
  "pooling": "all_tokens_legacy",
  "pooling_note": "Trained 2026-02-14 (commit 679b667) with mean over all token positions including BOS, EOS and padding, so each training embedding depended on its batch. Batch size and order were not recorded. Emulation (batch 4, input order) is approximate.",
  "max_residues": 1022,
  "truncation": "truncate_at_max_residues",
  "classifier": "LogisticRegression, C=1.0, unscaled",
  "training_data": "FLO/ALS homologs (positive) vs random proteins (negative); trained before duplicate removal (167 exact duplicates)",
  "validation": [],
  "threshold": {"value": 0.5, "chosen_on": "default (uncalibrated)"},
  "scope": "FLO/ALS-like surface glycoprotein score. Not an adhesin predictor: about 12% of S288C calls are known adhesins."
}
```

- [ ] **Step 5: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "card: v2 fields, settings resolution, and cards for the shipped pickles

Settings come from the card; a missing card or legacy pooling is refused
unless --allow-legacy-pooling is given (resolve_embedding_settings).

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Legacy pooling mode and truncation count in `embeddings.py`

**Files:**
- Modify: `src/surface_glyco/embeddings.py` (`_embed_batch`, `get_esm_embeddings`, new `count_truncated`)
- Test: `tests/surface_glyco/test_embeddings.py`

**Interfaces:**
- Consumes: `POOLING_RESIDUE_MEAN`, `POOLING_LEGACY`, `KNOWN_POOLINGS`, `MAX_RESIDUES` from `surface_glyco.card` (Task 5).
- Produces: `get_esm_embeddings(sequences, model_name=..., batch_size=None, device=None, repr_layer=DEFAULT_REPR_LAYER, return_indices=False, pooling=POOLING_RESIDUE_MEAN)`. With `pooling=POOLING_LEGACY`: input order (no length sorting), default batch size 4, mean over all positions. `count_truncated(sequences) -> int`: how many sequences exceed `MAX_RESIDUES` after sanitizing.

- [ ] **Step 1: Write the failing tests** (the file already has `_seqs` and `_cpu` helpers; reuse them)

```python
def test_legacy_pooling_depends_on_batch_mates_residue_mean_does_not():
    short = ("s", "MKT")
    long = ("l", "MKTAYIAKQRQISFVKSHFSRQ")
    solo = get_esm_embeddings(_seqs(short), device=_cpu(), batch_size=4)[0][0]
    pair_res = get_esm_embeddings(_seqs(short, long), device=_cpu(), batch_size=4)[0][0]
    assert np.allclose(solo, pair_res, atol=1e-5)

    solo_l = get_esm_embeddings(_seqs(short), device=_cpu(), pooling=POOLING_LEGACY)[0][0]
    pair_l = get_esm_embeddings(_seqs(short, long), device=_cpu(), pooling=POOLING_LEGACY)[0][0]
    assert not np.allclose(solo_l, pair_l, atol=1e-5)


def test_legacy_pooling_matches_a_direct_mean_over_all_tokens():
    model, alphabet = get_cached_model("esm2_t6_8M_UR50D", _cpu())
    _, _, toks = alphabet.get_batch_converter()([("a", "MKT"), ("b", "MKTAYIAKQRQ")])
    with torch.no_grad():
        reps = model(toks, repr_layers=[6])["representations"][6]
    expected = reps.mean(dim=1).numpy()
    got = get_esm_embeddings(
        _seqs(("a", "MKT"), ("b", "MKTAYIAKQRQ")), device=_cpu(), pooling=POOLING_LEGACY,
        batch_size=2,
    )[0]
    assert np.allclose(got, expected, atol=1e-5)


def test_unknown_pooling_raises():
    with pytest.raises(ValueError, match="pooling"):
        get_esm_embeddings(_seqs(("a", "MKT")), device=_cpu(), pooling="max")


def test_count_truncated_counts_sequences_over_the_limit():
    seqs = _seqs(("a", "M" * 1022), ("b", "M" * 1023), ("c", ""), ("d", "M*" * 600))
    assert count_truncated(seqs) == 2  # b, and d (1200 residues after '*' is stripped)
```
Add the imports at the top of the file: `import numpy as np`, `import torch`, and `from surface_glyco.embeddings import count_truncated, get_cached_model` and `from surface_glyco.card import POOLING_LEGACY` (merge with existing import lines).

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_embeddings.py -q`
Expected: FAIL (`pooling` argument and `count_truncated` do not exist). Requires the 8M weights (downloaded on first run, about 30 MB).

- [ ] **Step 3: Implement**

`_embed_batch` gets a `pooling` parameter:

```python
def _embed_batch(model, alphabet, batch_converter, batch, repr_layer, device, pooling=POOLING_RESIDUE_MEAN):
    """Embed one batch of (id, sequence) pairs.

    residue_mean: mean over residue tokens only, so an embedding does not depend on batch mates.
    all_tokens_legacy: mean over every position including BOS, EOS and padding, as the shipped
    models were trained (before 2026-09-27). Batch-dependent by construction.
    """
    _, _, tokens = batch_converter(batch)
    tokens = tokens.to(device)
    with torch.no_grad():
        results = model(tokens, repr_layers=[repr_layer], return_contacts=False)
    reps = results["representations"][repr_layer]
    if pooling == POOLING_LEGACY:
        return reps.mean(dim=1).float().cpu().numpy()
    mask = (
        (tokens != alphabet.padding_idx)
        & (tokens != alphabet.cls_idx)
        & (tokens != alphabet.eos_idx)
    ).unsqueeze(-1)
    pooled = (reps * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)
    return pooled.float().cpu().numpy()
```
In `get_esm_embeddings`: add `pooling=POOLING_RESIDUE_MEAN` to the signature; after the device block add

```python
    if pooling not in KNOWN_POOLINGS:
        raise ValueError(f"pooling={pooling!r}; known: {KNOWN_POOLINGS}")
```
replace the batch-size block with

```python
    if batch_size is None:
        batch_size = 4 if pooling == POOLING_LEGACY else get_optimal_batch_size(device, model_name)
    if pooling == POOLING_LEGACY:
        print("Warning: emulating legacy all-token pooling (input order, batch size "
              f"{batch_size}); scores are approximate, see the model card.")
```
replace the `order = sorted(...)` line with

```python
    if pooling == POOLING_LEGACY:
        order = list(range(len(sequences)))  # input order, as when the models were trained
    else:
        order = sorted(range(len(sequences)), key=lambda k: len(cleaned[k]))
```
and pass `pooling` to both `_embed_batch` calls (the batch call and the one-at-a-time retry). Add:

```python
def count_truncated(sequences):
    """Number of sequences longer than MAX_RESIDUES after sanitizing (they get truncated)."""
    return sum(len(sanitize_sequence(s["sequence"])) > MAX_RESIDUES for s in sequences)
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS. If `test_legacy_pooling_matches_a_direct_mean_over_all_tokens` fails, the emulation is wrong: do not loosen the tolerance, fix the pooling.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "embeddings: legacy all-token pooling mode and truncation count

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Card-driven `predict`, `evaluate`, `train`

**Files:**
- Modify: `src/surface_glyco/scripts/predict.py`, `evaluate.py`, `train.py`
- Test: `tests/surface_glyco/test_scripts_cards.py` (new); update `test_io_model_predict.py` dedupe tests

**Interfaces:**
- Consumes: `resolve_embedding_settings`, `ModelCardError`, `new_card`, `POOLING_RESIDUE_MEAN` (Task 5); `get_esm_embeddings(..., pooling=)` and `count_truncated` (Task 6); `model_filename` (Task 3); `LABEL_*`/`SCORE_COLUMN` (Task 2).
- Produces:
  - `predict.main(input_path, model_path, output_file, model_name=None, silent=False, show_all=False, max_workers=None, allow_legacy_pooling=False)`; exits with status 1 and a message on stderr when the card refuses.
  - CLI flags `--allow-legacy-pooling` on `predict` and `evaluate`; `--model-name` default `None` on both.
  - `train.dedupe_sequences(sequences) -> list` now also drops a sequence that occurs with both labels (all copies, both classes); `train.prepare_data` returns `(sequences, n_duplicates_removed, n_conflicting_removed)`.
  - Trained model card built with `new_card(...)` plus `n_duplicates_removed`, `n_conflicting_removed`, `training_sequences_sha256`, and an `environment` dict (`numpy`, `scikit_learn`, `torch`, `fair_esm` version strings).

- [ ] **Step 1: Write the failing tests**

```python
"""predict/evaluate/train honour the model card."""

import json

import numpy as np
import pytest
from sklearn.linear_model import LogisticRegression

from surface_glyco.card import POOLING_LEGACY, POOLING_RESIDUE_MEAN, new_card
from surface_glyco.model import save_model, save_model_card
from surface_glyco.scripts import predict as predict_mod
from surface_glyco.scripts import train as train_mod


def _model(tmp_path, card=None, n_features=4):
    rng = np.random.default_rng(0)
    clf = LogisticRegression().fit(rng.normal(size=(20, n_features)), [0, 1] * 10)
    path = tmp_path / "m.pkl"
    save_model(clf, path)
    if card is not None:
        save_model_card(path, card)
    return path


@pytest.fixture
def fasta(tmp_path):
    p = tmp_path / "in.fa"
    p.write_text(">a\nMKT\n>b\n" + "M" * 1100 + "\n")
    return p


@pytest.fixture
def fake_embed(monkeypatch):
    calls = {}

    def fake(sequences, **kw):
        calls.update(kw)
        n = len(sequences)
        return np.ones((n, 4)), [s["id"] for s in sequences]

    monkeypatch.setattr(predict_mod, "get_esm_embeddings", fake)
    return calls


def test_card_settings_reach_the_embedder(tmp_path, fasta, fake_embed):
    card = new_card("esm2_t6_8M_UR50D", 6, POOLING_RESIDUE_MEAN, 1, 1)
    out = tmp_path / "o.csv"
    predict_mod.main(fasta, _model(tmp_path, card), out, None, silent=True)
    assert fake_embed["pooling"] == POOLING_RESIDUE_MEAN
    assert fake_embed["repr_layer"] == 6
    assert fake_embed["model_name"] == "esm2_t6_8M_UR50D"
    assert out.read_text().splitlines()[0] == "id,prediction,surface_glycoprotein_score"


def test_mismatched_model_name_exits_naming_the_field(tmp_path, fasta, fake_embed, capsys):
    card = new_card("esm2_t6_8M_UR50D", 6, POOLING_RESIDUE_MEAN, 1, 1)
    with pytest.raises(SystemExit) as e:
        predict_mod.main(fasta, _model(tmp_path, card), tmp_path / "o.csv", "esm2_t12_35M_UR50D")
    assert e.value.code == 1
    assert "esm_model" in capsys.readouterr().err
    assert fake_embed == {}  # refused before any embedding was requested


def test_legacy_model_is_refused_without_the_flag(tmp_path, fasta, fake_embed, capsys):
    card = new_card("esm2_t6_8M_UR50D", 6, POOLING_LEGACY, 1, 1)
    with pytest.raises(SystemExit):
        predict_mod.main(fasta, _model(tmp_path, card), tmp_path / "o.csv", None)
    assert "--allow-legacy-pooling" in capsys.readouterr().err
    assert fake_embed == {}


def test_legacy_model_with_flag_uses_legacy_pooling(tmp_path, fasta, fake_embed):
    card = new_card("esm2_t6_8M_UR50D", 6, POOLING_LEGACY, 1, 1)
    predict_mod.main(
        fasta, _model(tmp_path, card), tmp_path / "o.csv", None, silent=True,
        allow_legacy_pooling=True,
    )
    assert fake_embed["pooling"] == POOLING_LEGACY


def test_model_without_card_is_refused_without_the_flag(tmp_path, fasta, fake_embed):
    with pytest.raises(SystemExit):
        predict_mod.main(fasta, _model(tmp_path), tmp_path / "o.csv", None)


def test_truncated_sequences_are_reported(tmp_path, fasta, fake_embed, capsys):
    card = new_card("esm2_t6_8M_UR50D", 6, POOLING_RESIDUE_MEAN, 1, 1)
    predict_mod.main(fasta, _model(tmp_path, card), tmp_path / "o.csv", None, silent=True)
    assert "1 of 2 sequences" in capsys.readouterr().out


def test_dedupe_drops_sequences_present_in_both_classes():
    seqs = [
        {"id": "a", "sequence": "MK", "label": 1},
        {"id": "b", "sequence": "MK", "label": 0},
        {"id": "c", "sequence": "MKT", "label": 1},
        {"id": "d", "sequence": "MKT", "label": 1},
    ]
    kept = train_mod.dedupe_sequences(seqs)
    assert [s["id"] for s in kept] == ["c"]


def test_empty_sequences_are_not_trained_on():
    seqs = [{"id": "a", "sequence": "", "label": 1}, {"id": "b", "sequence": "MK", "label": 0}]
    assert [s["id"] for s in train_mod.dedupe_sequences(seqs)] == ["b"]
```
Replace the existing `test_dedupe_sequences_keeps_first_id_per_label` in `test_io_model_predict.py` with the two dedupe tests above (its expectation, that a sequence stays in both classes, is the behaviour being removed). Add one end-to-end test with the real 8M model:

```python
def test_end_to_end_five_sequences_with_real_embeddings(tmp_path):
    seqs = [
        {"id": f"s{i}", "sequence": "MKTAYIAKQRQISFVKSHFSRQ"[: 8 + i], "label": i % 2}
        for i in range(6)
    ]
    from surface_glyco.embeddings import get_esm_embeddings

    emb, ids = get_esm_embeddings(seqs, model_name="esm2_t6_8M_UR50D")
    clf = LogisticRegression(max_iter=200).fit(emb, [s["label"] for s in seqs])
    card = new_card("esm2_t6_8M_UR50D", 6, POOLING_RESIDUE_MEAN, 3, 3)
    model = tmp_path / "m.pkl"
    save_model(clf, model)
    save_model_card(model, card)
    fa = tmp_path / "in.fa"
    fa.write_text("".join(f">{s['id']}\n{s['sequence']}\n" for s in seqs))
    out = tmp_path / "o.csv"
    predict_mod.main(fa, model, out, None, silent=True)
    rows = out.read_text().splitlines()
    assert len(rows) == 7 and rows[0] == "id,prediction,surface_glycoprotein_score"
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_scripts_cards.py -q`
Expected: FAIL (`main` has no `allow_legacy_pooling`, no card-driven settings).

- [ ] **Step 3: Implement `predict.py`**

Change the signature to the one in **Interfaces**. Replace the card block with:

```python
    classifier = load_model(model_path)
    card = load_model_card(model_path)
    try:
        settings = resolve_embedding_settings(
            card, model_name, DEFAULT_MODEL, allow_legacy_pooling=allow_legacy_pooling
        )
    except ModelCardError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    if settings.implied:
        print("Warning: no model card; assuming legacy embedding settings", file=sys.stderr)
```
(Order: load the model and card, resolve settings, only then read the FASTA and embed.) After the sequences are loaded:

```python
    n_trunc = count_truncated(all_sequences)
    if n_trunc:
        print(f"{n_trunc} of {len(all_sequences)} sequences are longer than "
              f"{settings.max_residues} residues and will be truncated")
```
Call the embedder as `get_esm_embeddings(all_sequences, model_name=settings.esm_model, repr_layer=settings.repr_layer, pooling=settings.pooling)`. In `cli()`: `--model-name` gets `default=None`; add `--allow-legacy-pooling` (`store_true`); when `--model` is omitted use `MODELS_DIR / model_filename(args.model_name or DEFAULT_MODEL)`; pass the flag to `main`.

- [ ] **Step 4: Implement `evaluate.py` and `train.py`**

`evaluate.py`: same `resolve_embedding_settings` block (error to stderr, `sys.exit(1)`), same `--model-name default=None` and `--allow-legacy-pooling` flag, and pass `repr_layer` and `pooling` to `get_esm_embeddings`. `train.py`: remove the `POOLING` import; `dedupe_sequences` becomes

```python
def dedupe_sequences(sequences):
    """Drop empty sequences, exact duplicates within a class, and sequences in both classes.

    A sequence present with both labels cannot be learned and leaks across any split, so every
    copy is dropped. Duplicates within one class keep their first id.
    """
    labels_by_seq = {}
    for seq in sequences:
        labels_by_seq.setdefault(seq["sequence"], set()).add(seq["label"])
    seen = set()
    unique = []
    for seq in sequences:
        key = (seq["label"], seq["sequence"])
        if not seq["sequence"] or len(labels_by_seq[seq["sequence"]]) > 1 or key in seen:
            continue
        seen.add(key)
        unique.append(seq)
    return unique
```
`prepare_data` returns `(unique, n_duplicates_removed, n_conflicting_removed)` where `n_conflicting_removed` is the number of input records whose sequence has both labels; `main` unpacks all three. Build the card with `new_card(model_name, DEFAULT_REPR_LAYER, POOLING_RESIDUE_MEAN, int(labels.sum()), int(len(labels) - labels.sum()), classifier=type(classifier).__name__, n_duplicates_removed=..., n_conflicting_removed=..., n_not_embedded=..., training_sequences_sha256=sequences_sha256(sequences), positive_dir=str(positive_dir), negative_dir=str(negative_dir), holdout_accuracy=float(test_acc), environment={...})` where `environment` maps `numpy`, `scikit_learn`, `torch`, `fair_esm` to `importlib.metadata.version(<dist>)` (`"fair-esm"` for the last), using `"unknown"` on `PackageNotFoundError`.

- [ ] **Step 5: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS (about 35 tests).

- [ ] **Step 6: Check the real shipped model end to end**

Run (needs the 8M weights; CPU is fine):
```bash
surface_glyco_predict --input tests/input_tests/Saccharomyces.pep --output /tmp/x.csv
```
Expected: exits 1 with a message containing `--allow-legacy-pooling` (the shipped card is legacy). Then
```bash
surface_glyco_predict --input tests/input_tests/Saccharomyces.pep --output /tmp/x.csv --allow-legacy-pooling
```
Expected: exits 0, prints the emulation warning, `/tmp/x.csv` has the new header. (If the input file does not exist in your checkout, use any small FASTA.)

- [ ] **Step 7: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "predict/evaluate/train: use the model card, refuse legacy settings, report truncation

Models with legacy pooling or no card need --allow-legacy-pooling. Training
now drops sequences present in both classes and empty sequences, and writes
a v2 card with environment versions.

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Docs, CHANGELOG, release prep

**Files:**
- Modify: `CHANGELOG.md`, `README.md`, `AGENTS.md`, `docs/TOOL-ARCHITECTURE.md` (section 2.0), `analysis/chytrid_batrach/NOTES.md` (one note)
- Test: `tests/surface_glyco/test_docs.py` (new)

**Interfaces:**
- Produces: `CHANGELOG.md` with a `## [0.1.0] - 2026-02-16` section so that `bump minor` yields 0.2.0, and a `### Changed` entry marked `**BREAKING**` under Unreleased.

- [ ] **Step 1: Write the failing test**

```python
"""Docs and changelog agree with the code."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_changelog_lists_the_breaking_rename_and_a_0_1_0_baseline():
    text = (ROOT / "CHANGELOG.md").read_text()
    assert "## [0.1.0] - 2026-02-16" in text
    unreleased = text.split("## [Unreleased]")[1].split("## [0.1.0]")[0]
    assert "BREAKING" in unreleased and "surface_glyco" in unreleased


def test_readme_uses_current_names_only():
    text = (ROOT / "README.md").read_text()
    assert "adhesion_train" not in text and "adhesion_predict " not in text
    assert "surface_glyco_train" in text


def test_agents_does_not_claim_a_thread_safe_cache():
    assert "thread-safe" not in (ROOT / "AGENTS.md").read_text().lower()
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/surface_glyco/test_docs.py -q`
Expected: FAIL (3).

- [ ] **Step 3: Edit the docs**

`CHANGELOG.md`: under `## [Unreleased]` add

```markdown
### Changed
- **BREAKING:** the package `adhesion_predict` is now `surface_glyco`. Entry points are `surface_glyco_predict`, `surface_glyco_train`, `surface_glyco_evaluate`. The output column `probability_adhesion` is now `surface_glycoprotein_score`; labels `Adhesion` / `Non-adhesion` are now `surface_glycoprotein` / `other`; model files are `surface_glyco_model_<esm>.pkl`; default output files end in `.surface_glyco.csv`. No aliases. Reason: the model scores FLO/ALS-like surface glycoproteins, not adhesins (about 12% of S288C calls are known adhesins).
- Models are read with their card: `predict` and `evaluate` take the ESM model, layer and pooling from it and stop with a message naming the field on a mismatch.
- A `./models` directory in the working directory no longer replaces the packaged models. Set `SURFACE_GLYCO_MODELS_DIR` to use another directory.
- `kingdom_survey` counts called proteins instead of all rows in a result file.

### Added
- `--allow-legacy-pooling` on `predict` and `evaluate`. The two shipped models were trained with padding-inclusive pooling (issue #25) and are refused without it.
- Cards for the shipped models; truncation counts in `predict`; environment versions in new cards.

### Fixed
- Training drops sequences that occur in both classes and empty sequences.
```
and append at the end of the file, below the existing sections:

```markdown
## [0.1.0] - 2026-02-16

Initial tagged version (`v0.1.0`).
```
`README.md`: replace the old names (`adhesion_train`, `adhesion_predict`, `adhesion_evaluate`, `probability_adhesion`, `Adhesion`) with the new ones, show `--allow-legacy-pooling` in the usage example for the shipped models, and keep the scope paragraph added in PR #28. `AGENTS.md`: delete the "thread-safe" claim; change the sentence to "model cache is a module-level dict and is not thread-safe". `docs/TOOL-ARCHITECTURE.md` section 2.0: change "Proposed names (proposal only; no code ... renamed)" to "Names (implemented in 0.2.0)" and mark the rows done. `analysis/chytrid_batrach/NOTES.md`: add "These scripts read frozen result files written by `adhesion_predict` < 0.2.0 (`*.adhesion_predict.csv`, column `probability_adhesion`). They were not updated."

- [ ] **Step 4: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "docs: changelog (breaking), README, AGENTS, architecture names

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Whole-branch verification and release

**Files:** none modified unless a check fails.

- [ ] **Step 1: Stale-name scan**

Run: `git grep -n "adhesion_predict\|adhesion_train\|adhesion_evaluate\|adhesion_model_" -- . ':!docs/superpowers/plans/2026-09-04-*' ':!docs/superpowers/specs/2026-09-04-*' ':!Changes.md' ':!docs/model-review' ':!analysis/chytrid_batrach' ':!CHANGELOG.md' ':!docs/PLAN-2026-09-30-pipeline-and-decisions.md' ':!docs/plans' ':!docs/superpowers/plans/2026-09-30-*'`
Expected: no output. Anything listed is either fixed or added to the allowed list with a reason in the PR description.

- [ ] **Step 2: All test folders, one process each**

```bash
for d in surface_glyco kingdom_survey adhesion_properties embedding_clustering; do pytest tests/$d -q || exit 1; done
```
Expected: all PASS.

- [ ] **Step 3: Build a wheel and check its contents** (pins the packaging gap in Review Focus 6)

```bash
pip wheel . --no-deps -w /tmp/wheel_check -q
python - <<'E'
import glob, zipfile
w = glob.glob("/tmp/wheel_check/surface_glyco-*.whl")[0]
names = zipfile.ZipFile(w).namelist()
need = [f"surface_glyco/models/surface_glyco_model_{e}.{x}"
        for e in ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D") for x in ("pkl", "json")]
missing = [n for n in need if n not in names]
assert not missing, missing
assert not any(n.startswith("adhesion_predict/") for n in names)
print("wheel ok:", w)
E
```
Expected: `wheel ok`. A `0+unknown` or `0.1.dev...` version in the wheel name is acceptable on an untagged branch.

- [ ] **Step 4: Push and open the PR; wait for CI**

```bash
git push -u origin <branch>
gh pr create --base main --title "Rename adhesion_predict to surface_glyco; card-driven legacy safety (0.2.0)" --body "<summary, breaking list from CHANGELOG, Review Focus results>"
gh pr checks --watch
```
Expected: lint, unit tests and the three analysis jobs pass.

- [ ] **Step 5: After merge, cut 0.2.0 (owner action or explicit go-ahead)**

Run the `CD/version bump` workflow with `bump_type = minor` (Actions tab, or `gh workflow run version_bump.yml -f bump_type=minor`). It bumps the changelog to 0.2.0 and dispatches the release. Verify: `gh release view v0.2.0` shows the notes, and `git describe --tags` on `main` after the release commit reports `v0.2.0`. This is the first run of this flow; if it fails, report the log and do not retry blindly.

- [ ] **Step 6: Update the other docs**

Edit `docs/PLAN-2026-09-30-pipeline-and-decisions.md` section 9 (paths to the merged `docs/plans/*` files) and mark the item 1 work done in its section 7. Commit and open a small docs PR.

---

## Self-review

**Spec coverage.** Decision 2 (freeze, safety only): Tasks 5-7; no retrain (constraint, byte-hash test in Task 3). Decision 4 names: Tasks 1-3. Decision 5 new tools get their own packages: nothing here adds a tool. Default 3 (one model location, drop root copies): Task 3. Default 4 (both-class sequences): Task 7. Default 1 (threshold 0.5 documented): card `threshold` field, Task 5. Default 6 (Fungi_5k outputs legacy): Task 4 readers. Fable items: 3/4 (card read back, truncation count) Tasks 5-7; 5/6 (cards for pickles) Task 5; 26 (t6 pickle untracked) Task 3 and the wheel check; 27 (`ESM2_MODELS`) Task 3; 31 (cwd `./models`) Task 3; 32 (cross-label duplicates) Task 7; 33 (thread-safe claim) Task 8; 30 (CI) already merged, path updated in Task 1. Not covered here, by design: Fable items 8, 13, 17-25, 28 (separate projects).

**Placeholders.** None intended. Two places depend on the reader's checkout: the existence of `tests/input_tests/Saccharomyces.pep` (Task 7 step 6 says what to do if absent) and the branch name in Task 9 step 4.

**Type consistency.** `EmbeddingSettings` fields (`esm_model`, `repr_layer`, `pooling`, `max_residues`, `implied`) are used with those names in Tasks 5 and 7. `pooling` constants are defined in `card.py` and imported, not redefined. `SCORE_COLUMN`/`LABEL_*` (Task 2) match the strings in `results.CALLED_LABELS` (Task 4) and the `new_card`/test strings. `prepare_data` returns 3 values in Task 7; no earlier task unpacks it.

**Risks and what is not verified.**
- The legacy emulation (batch 4, input order) is approximate. Task 6 pins that it equals a direct all-token mean for a given batch, not that it equals the Feb 2026 training embeddings.
- `bump minor` on a changelog whose only release is the hand-written 0.1.0 stub has never been run; Task 9 step 5 is its first test.
- Renaming the distribution means anyone with `adhesion_predict` installed must reinstall; Task 1 step 7 covers this checkout only.
- GPU-versus-CPU embedding equality is not tested (CI is CPU only).
