# surface_glyco rename and model-card framework: Implementation Plan (revision 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rename the package `adhesion_predict` to `surface_glyco` (breaking, no aliases), remove the old models and the legacy-pooling problem, and build the model-card framework that the new surface-glycoprotein model will be trained and scored with.

**Architecture:** One mechanical rename task first. Then small test-first changes: output schema, model location, a result reader and file discovery for the new schema, a torch-free `card.py` that decides what embedding settings a model needs, card-driven `predict`/`evaluate`/`train`. The two shipped pickles (FLO/ALS vs random labels, padding-inclusive pooling) are deleted; the package ships no model until the new model exists. No release is cut by this plan.

**Tech Stack:** Python 3.12, setuptools + setuptools-scm, fair-esm (ESM-2 8M in tests), scikit-learn, pytest, ruff 0.3.5 via pre-commit, GitHub Actions.

**Spec:** `docs/PLAN-2026-09-30-pipeline-and-decisions.md` and the changes of 2026-09-30 (later) recorded in Task 8; `docs/plans/2026-09-30-design-review-fable.md`. This revision incorporates the independent review of revision 1 (listed in the last section).

## Global Constraints

- **Names, exactly:** package `surface_glyco`; entry points `surface_glyco_predict`, `surface_glyco_train`, `surface_glyco_evaluate`; output column `surface_glycoprotein_score`; labels `surface_glycoprotein` and `other`; model files `surface_glyco_model_<esm>.pkl`; result files `<name>.surface_glyco.csv`. Repo name `adhesionPred` stays.
- **Breaking rename, no aliases, no reader for the old result schema.** The tool is not in circulation. Frozen files written by `adhesion_predict` are not read by the analysis code after this change.
- **No retraining in this plan.** The glycoprotein-label model is a separate project with its own spec, plan and review. This plan must not train or ship a model.
- **Old pickles are deleted.** `models/adhesion_model_*.pkl` and `src/adhesion_predict/models/adhesion_model_*.pkl` (identical copies, md5 `ce6dc8a40809d6afb8269ee3280c9910` for t6_8M and `5f7a08cc0122e3e43bc18c1539f77d7f` for t12_35M) stay available in git history only.
- **Pooling:** only `residue_mean`. A card with any other `pooling` value is refused. `repr_layer` default 6, `MAX_RESIDUES` 1022, default ESM model `esm2_t6_8M_UR50D`, threshold 0.5 documented as uncalibrated.
- **Python 3.12 is canonical; CI is CPU only.** `requires-python` becomes `>=3.11`. Lint is ruff 0.3.5 via pre-commit, line length 100; run `pre-commit run --all-files` before each commit. Each commit must pass `ruff check .` (import only what the task uses).
- **No release is cut here.** The 0.2.0 tag is made later, when a validated model ships, through the `version_bump` workflow (`bump_type=minor`). The version comes from git tags (setuptools-scm).
- **Do not edit historical documents:** `docs/superpowers/plans/2026-09-04-*`, `docs/superpowers/specs/2026-09-04-*`, `Changes.md`, `docs/model-review/*`, `analysis/embedding_clustering/REPORT.md` (generated), and `analysis/chytrid_batrach/*` (one-off analyses of frozen old result files).
- Each commit message ends with: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
- Until Task 1 step 7 reinstalls the package, run tests with `PYTHONPATH=$PWD/src`. Another checkout in the same conda environment may have an editable `adhesion_predict` install; do not uninstall it (Task 1 handles this in the test).

## Review Focus

1. **Result discovery:** after the rename, `kingdom_survey` and `adhesion_properties` must find `*.surface_glyco.csv` files, not only read them. Pinned in Task 4.
2. **All-rows result files** (every protein is written since PR #23): called proteins are counted, not rows. Pinned in Task 4.
3. **A `./models` directory in the working directory** must not replace the packaged model; **no model present** must give a clear message, not a traceback. Pinned in Task 3.
4. **A card that is malformed JSON, lacks `esm_model`, has a non-integer `repr_layer`, an unknown `pooling`, or names a different ESM model** must stop `predict`/`evaluate` with an error that names the field. Pinned in Tasks 5 and 7.
5. **Empty sequences, sequences over 1022 residues, sequences present in both classes, and a batch that fails to embed** must not corrupt counts or labels. Pinned in Tasks 6 and 7.
6. **Import without install:** importing the package from a checkout where it is not installed under the new name must not raise; a built wheel must contain the package and no pickles. Pinned in Tasks 1 and 9.

---

## File structure

| File | Responsibility | Task |
|---|---|---|
| `src/surface_glyco/` (from `src/adhesion_predict/`) | the package | 1 |
| `src/surface_glyco/models/README.md` | explains that no model ships yet | 3 |
| `src/surface_glyco/card.py` (new) | card constants, `EmbeddingSettings`, `resolve_embedding_settings`, `new_card`; no torch import | 5 |
| `src/surface_glyco/results.py` (new) | read result CSVs | 4 |
| `src/surface_glyco/embeddings.py` | `count_truncated`; constants from `card.py` | 6 |
| `src/surface_glyco/model.py` | card load/save (malformed JSON becomes `ModelCardError`) | 5 |
| `src/surface_glyco/scripts/{predict,train,evaluate}.py` | schema, card use | 2, 7 |
| `src/surface_glyco/config.py` | model directory, `model_filename`; delete `ESM2_MODELS` | 3 |
| `tests/surface_glyco/` (from `tests/adhesion_predict/`) | package tests | 1-7 |
| `analysis/kingdom_survey/{join.py,01_build_species_table.py}`, `analysis/adhesion_properties/universe.py` | result suffix, `read_called` | 4 |
| `analysis/model_review/02_cv_and_proteome_eval.py` | `features.py` path | 1 |
| `pyproject.toml`, `.coveragerc`, `CITATION.cff`, `.github/workflows/build_and_test.yml` | names, Python floor | 1 |
| `CHANGELOG.md`, `README.md`, `AGENTS.md`, `docs/TOOL-ARCHITECTURE.md`, `docs/PLAN-2026-09-30-pipeline-and-decisions.md` | docs | 8 |

---

### Task 1: Rename the package, entry points and test folder

**Files:**
- Rename: `src/adhesion_predict/` to `src/surface_glyco/`; `tests/adhesion_predict/` to `tests/surface_glyco/`
- Modify: `pyproject.toml` (name, scripts, package-data, `requires-python`, drop the dead `importlib-metadata` line), `.coveragerc:3`, `CITATION.cff` (title; add `version:` and `date-released:` lines), `.github/workflows/build_and_test.yml` (2 lines), `src/surface_glyco/__init__.py`, imports in `src/`, `tests/surface_glyco/`, `analysis/embedding_clustering/{00_feasibility_check,embed,recover_missing_esm2_embeddings}.py`, `analysis/model_review/02_cv_and_proteome_eval.py:21-23`
- Test: `tests/surface_glyco/test_package.py` (new)

**Interfaces:**
- Produces: importable package `surface_glyco`; `surface_glyco.__version__` (str, never raises); console scripts `surface_glyco_predict|train|evaluate` at `surface_glyco.scripts.<name>:cli`.

- [ ] **Step 1: Write the failing test** (create the file in step 3; the text is final)

```python
"""Package identity: names, version fallback, entry points."""

import importlib
import sys
from pathlib import Path

import pytest

tomllib = pytest.importorskip("tomllib")  # Python >= 3.11

ROOT = Path(__file__).resolve().parents[2]


def test_new_package_imports_and_has_a_version_string():
    import surface_glyco

    assert isinstance(surface_glyco.__version__, str) and surface_glyco.__version__


def test_old_package_is_gone_from_this_checkout():
    assert not (ROOT / "src" / "adhesion_predict").exists()
    try:
        old = importlib.import_module("adhesion_predict")
    except ModuleNotFoundError:
        return
    # An editable install from another checkout may exist in this environment; ours must not.
    assert not Path(old.__file__).resolve().is_relative_to(ROOT)
    sys.modules.pop("adhesion_predict", None)


def test_pyproject_names_match_the_decision():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert cfg["project"]["name"] == "surface_glyco"
    assert cfg["project"]["scripts"] == {
        "surface_glyco_predict": "surface_glyco.scripts.predict:cli",
        "surface_glyco_train": "surface_glyco.scripts.train:cli",
        "surface_glyco_evaluate": "surface_glyco.scripts.evaluate:cli",
    }
    assert cfg["tool"]["setuptools"]["package-data"] == {"surface_glyco": ["models/*"]}
    assert cfg["project"]["requires-python"] == ">=3.11"
```

- [ ] **Step 2: Record the baseline**

Run: `PYTHONPATH=$PWD/src python -m pytest tests/adhesion_predict -q -p no:cacheprovider`
Expected: 12 passed. After step 3 and before step 4, `import surface_glyco` raises `PackageNotFoundError` (version lookup) or `ModuleNotFoundError`; either is the expected failing state.

- [ ] **Step 3: Move the folders and add the test**

```bash
git mv src/adhesion_predict src/surface_glyco
git mv tests/adhesion_predict tests/surface_glyco
```
Create `tests/surface_glyco/test_package.py` with the text from step 1.

- [ ] **Step 4: Rewrite names**

```bash
grep -rl "adhesion_predict" src tests/surface_glyco analysis/embedding_clustering/*.py tests/embedding_clustering \
  analysis/model_review/02_cv_and_proteome_eval.py .coveragerc CITATION.cff .github/workflows/build_and_test.yml \
  | xargs sed -i 's/adhesion_predict/surface_glyco/g'
```
Do not run sed over `analysis/embedding_clustering/REPORT.md` (the glob above is `*.py`). Then edit by hand:

`pyproject.toml`:
```toml
requires-python = ">=3.11"

[project]
name = "surface_glyco"
description = "ESM-2 based scorer for fungal cell-surface glycoproteins"

[project.scripts]
surface_glyco_predict = "surface_glyco.scripts.predict:cli"
surface_glyco_train = "surface_glyco.scripts.train:cli"
surface_glyco_evaluate = "surface_glyco.scripts.evaluate:cli"

[tool.setuptools.package-data]
"surface_glyco" = ["models/*"]
```
(`requires-python` stays inside `[project]` where it already is; delete the `importlib-metadata; python_version<'3.8'` dependency; update the header comment about the `adhesion_` script prefix.) In `CITATION.cff` set `title: surface_glyco` and add the two lines the release workflow edits with sed:
```yaml
version: 0.1.0
date-released: '2026-02-16'
```
In `.github/workflows/build_and_test.yml` change `pytest tests/adhesion_predict` to `pytest tests/surface_glyco` and the `--ignore=tests/adhesion_predict` path to `--ignore=tests/surface_glyco`.

`src/surface_glyco/__init__.py`:
```python
"""surface_glyco - ESM-2 based scorer for fungal cell-surface glycoproteins."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("surface_glyco")
except PackageNotFoundError:  # a checkout that was not pip-installed under this name
    __version__ = "0+unknown"
```
and keep the existing `from surface_glyco.io import ...` / `from surface_glyco.model import ...` and `__all__`.

In `analysis/model_review/02_cv_and_proteome_eval.py` the `features.py` path (lines 21-23) becomes `str(REPO) + "/src/surface_glyco/features.py"` (sed in step 4 already does this; confirm). The `shipped = pickle.load(...)` block near line 95 reads a deleted file after Task 3; wrap it as `if (REPO / "src/surface_glyco/models/<file>").exists():` in Task 3.

- [ ] **Step 5: Run tests**

Run: `PYTHONPATH=$PWD/src python -m pytest tests/surface_glyco -q -p no:cacheprovider`
Expected: 15 passed (12 existing + 3 new).

- [ ] **Step 6: Stale-reference scan of code**

Run: `git grep -n "adhesion_predict\|adhesion_train\|adhesion_evaluate" -- src tests analysis/embedding_clustering analysis/kingdom_survey analysis/adhesion_properties analysis/model_review pyproject.toml .coveragerc CITATION.cff .github`
Expected: only the strings Task 4 and Task 3 change later (result suffix, pickle path). Note each remaining hit; none may be an `import`.

- [ ] **Step 7: Reinstall in a throwaway environment and run the entry points**

```bash
python -m venv /tmp/sg_venv --system-site-packages && /tmp/sg_venv/bin/pip install -e . --no-deps -q
/tmp/sg_venv/bin/surface_glyco_predict --help | head -3
/tmp/sg_venv/bin/python -c "import surface_glyco; print(surface_glyco.__version__)"
```
Expected: help prints; a version string prints (not `0+unknown`). Using a throwaway environment leaves other checkouts' editable installs alone.

- [ ] **Step 8: Commit**

```bash
pre-commit run --all-files
git add -A
git commit -m "rename: package adhesion_predict to surface_glyco, with new entry points

Breaking, no aliases. __version__ falls back to 0+unknown when the package
is not installed under the new name. requires-python is now >=3.11.

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Output schema and wording

**Files:**
- Modify: `src/surface_glyco/scripts/predict.py`, `evaluate.py`, `train.py`
- Test: `tests/surface_glyco/test_io_model_predict.py`

**Interfaces:**
- Produces, in `surface_glyco.scripts.predict`: `LABEL_POSITIVE = "surface_glycoprotein"`, `LABEL_NEGATIVE = "other"`, `SCORE_COLUMN = "surface_glycoprotein_score"`; `build_results(seq_ids, predictions, probabilities) -> list[dict]` with keys `id`, `prediction`, `surface_glycoprotein_score`; `write_results(results, path)` writing header `id,prediction,surface_glycoprotein_score`; `default_output_path(input_path) -> Path` ending `.surface_glyco.csv` in the current directory.

- [ ] **Step 1: Write the failing tests** (replace the two existing tests that use the old names)

```python
def test_write_results_header_and_quoting(tmp_path):
    out = tmp_path / "r.csv"
    write_results(
        [{"id": "sp|P1,x", "prediction": "surface_glycoprotein",
          "surface_glycoprotein_score": 0.91234}],
        out,
    )
    assert out.read_text().splitlines()[0] == "id,prediction,surface_glycoprotein_score"
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert rows == [{"id": "sp|P1,x", "prediction": "surface_glycoprotein",
                     "surface_glycoprotein_score": "0.9123"}]


def test_build_results_keeps_every_row_with_new_labels():
    rows = build_results(["a", "b"], [1, 0], [[0.1, 0.9], [0.8, 0.2]])
    assert [r["prediction"] for r in rows] == ["surface_glycoprotein", "other"]
    assert rows[1]["surface_glycoprotein_score"] == 0.2
    assert "probability_adhesion" not in rows[0]


def test_default_output_name_uses_new_suffix(tmp_path, monkeypatch):
    from surface_glyco.scripts.predict import default_output_path

    monkeypatch.chdir(tmp_path)
    assert default_output_path(Path("x/Scer.pep.fa")).name == "Scer.pep.surface_glyco.csv"
    (tmp_path / "proteomes").mkdir()
    assert default_output_path(tmp_path / "proteomes").name == "proteomes.surface_glyco.csv"
```
(Add `from pathlib import Path` to the test imports.)

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_io_model_predict.py -q`
Expected: FAIL (old header and labels; `default_output_path` missing).

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
`write_results` writes `["id", "prediction", SCORE_COLUMN]` and reads `result[SCORE_COLUMN]`. `build_results` emits `LABEL_POSITIVE if predictions[i] == 1 else LABEL_NEGATIVE` under key `SCORE_COLUMN`. In `main`, replace the `if output_file is None:` block with `if output_file is None: output_file = default_output_path(input_path)`; the stdout filter compares to `LABEL_POSITIVE` and prints `result[SCORE_COLUMN]`. Banner `"Surface glycoprotein scoring"`; parser description `"Score proteins as fungal cell-surface glycoproteins"`; error hint `"Run training first: surface_glyco_train"`. `evaluate.py`: `target_names = ["other", "surface_glycoprotein"]`, description `"Evaluate a surface glycoprotein model"`. `train.py`: banner `"Surface glycoprotein classifier training"`, description to match.

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

### Task 3: Model directory, no shipped model, delete the old pickles

**Files:**
- Delete: `models/adhesion_model_esm2_t6_8M_UR50D.pkl`, `models/adhesion_model_esm2_t12_35M_UR50D.pkl`, `src/surface_glyco/models/adhesion_model_esm2_t12_35M_UR50D.pkl` (and the untracked `src/.../adhesion_model_esm2_t6_8M_UR50D.pkl` if present on your checkout)
- Create: `src/surface_glyco/models/README.md`
- Modify: `src/surface_glyco/config.py`, the three scripts (default model path and missing-model message), `analysis/model_review/02_cv_and_proteome_eval.py` (~line 95)
- Test: `tests/surface_glyco/test_config_paths.py` (new)

**Interfaces:**
- Produces: `surface_glyco.config.get_models_dir() -> Path` (the packaged directory unless env `SURFACE_GLYCO_MODELS_DIR` is set; a `./models` in the working directory is ignored); `model_filename(esm_model: str) -> str` returning `f"surface_glyco_model_{esm_model}.pkl"`; `PACKAGED_MODELS_DIR`.

- [ ] **Step 1: Write the failing tests**

```python
"""Model directory lookup, file names, and the no-model case."""

import sys
from pathlib import Path

import pytest

from surface_glyco import config
from surface_glyco.scripts import predict as predict_mod


def test_model_filename():
    assert config.model_filename("esm2_t6_8M_UR50D") == "surface_glyco_model_esm2_t6_8M_UR50D.pkl"


def test_packaged_dir_has_readme_and_no_pickles():
    pkg = config.PACKAGED_MODELS_DIR
    assert (pkg / "README.md").is_file()
    assert list(pkg.glob("*.pkl")) == []


def test_stray_models_dir_in_cwd_does_not_override(tmp_path, monkeypatch):
    (tmp_path / "models").mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("SURFACE_GLYCO_MODELS_DIR", raising=False)
    assert config.get_models_dir() == config.PACKAGED_MODELS_DIR


def test_env_var_overrides(tmp_path, monkeypatch):
    monkeypatch.setenv("SURFACE_GLYCO_MODELS_DIR", str(tmp_path))
    assert config.get_models_dir() == tmp_path


def test_predict_without_a_model_explains_what_to_do(tmp_path, monkeypatch, capsys):
    fa = tmp_path / "in.fa"
    fa.write_text(">a\nMKT\n")
    monkeypatch.setenv("SURFACE_GLYCO_MODELS_DIR", str(tmp_path / "empty"))
    monkeypatch.setattr(sys, "argv", ["surface_glyco_predict", "--input", str(fa)])
    with pytest.raises(SystemExit) as e:
        predict_mod.cli()
    assert e.value.code == 1
    err = capsys.readouterr().out + capsys.readouterr().err
    assert "surface_glyco_train" in err or "--model" in err
```
(`config.MODELS_DIR` is read at import time; `cli()` must call `get_models_dir()` at run time, not use the import-time constant, for the last test to work. Step 3 does this.)

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_config_paths.py -q`
Expected: FAIL (`model_filename`, `PACKAGED_MODELS_DIR` missing; pickles present).

- [ ] **Step 3: Implement**

```bash
git rm -q models/adhesion_model_esm2_t6_8M_UR50D.pkl models/adhesion_model_esm2_t12_35M_UR50D.pkl \
          src/surface_glyco/models/adhesion_model_esm2_t12_35M_UR50D.pkl
rm -f src/surface_glyco/models/adhesion_model_esm2_t6_8M_UR50D.pkl   # untracked copy, if present
```
`src/surface_glyco/models/README.md`:
```markdown
# Models

No trained model ships with this package yet. The earlier models (trained on FLO/ALS homologs
versus random proteins, with padding-inclusive pooling) were removed; they remain in git history.
Train one with `surface_glyco_train`, or pass `--model path/to/model.pkl`. A model is accompanied
by a JSON card (`model.json`) that `surface_glyco_predict` reads. Set `SURFACE_GLYCO_MODELS_DIR`
to use a different directory.
```
`config.py`: replace `get_models_dir`, delete `ESM2_MODELS`, keep `DEFAULT_MODEL`:
```python
import os

PACKAGED_MODELS_DIR = Path(__file__).parent / "models"


def get_models_dir():
    """Directory holding trained models.

    The packaged directory is used unless SURFACE_GLYCO_MODELS_DIR is set. A ./models directory
    in the working directory is deliberately ignored so a stray folder cannot replace a model.
    """
    override = os.environ.get("SURFACE_GLYCO_MODELS_DIR")
    return Path(override) if override else PACKAGED_MODELS_DIR


def model_filename(esm_model):
    """File name of the model trained on embeddings from esm_model."""
    return f"surface_glyco_model_{esm_model}.pkl"
```
In `predict.py`, `train.py`, `evaluate.py`: import `get_models_dir`, `model_filename`; compute the default path inside `cli()` as `get_models_dir() / model_filename(args.model_name or DEFAULT_MODEL)` (train uses `args.model`). When the file does not exist, `predict` and `evaluate` print `Error: model file not found at <path>. Train one with surface_glyco_train, or pass --model.` and `sys.exit(1)`. In `analysis/model_review/02_cv_and_proteome_eval.py` wrap the `shipped = pickle.load(...)` block in `if SHIPPED.exists():` where `SHIPPED = REPO / "src/surface_glyco/models/surface_glyco_model_esm2_t6_8M_UR50D.pkl"`, with a comment that the model was removed on 2026-09-30 and can be restored from git history (`git show 7f97c9a:models/adhesion_model_esm2_t6_8M_UR50D.pkl`).

- [ ] **Step 4: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "models: remove the old pickles; ignore ./models in cwd; clear no-model message

The old models (FLO/ALS vs random labels, padding-inclusive pooling) stay in
git history. ESM2_MODELS removed (listed models that do not exist in fair-esm).

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Result files: discovery, reader, and call counting

**Files:**
- Create: `src/surface_glyco/results.py`
- Modify: `analysis/kingdom_survey/join.py` (`RESULT_SUFFIX` line 14; `count_result_rows`, `probability_stats`), `analysis/kingdom_survey/01_build_species_table.py:59`, `analysis/adhesion_properties/universe.py` (`RESULT_SUFFIX` line 18; `read_result_ids`)
- Test: `tests/surface_glyco/test_results.py` (new); update fixtures in `tests/kingdom_survey/test_join.py` and `tests/adhesion_properties/test_universe.py`

**Interfaces:**
- Produces: `surface_glyco.results.read_all(path) -> list[tuple[str, str, float]]` as `(id, label, score)` and `read_called(path) -> list[tuple[str, float]]` for rows labelled `surface_glycoprotein`. A file without the `surface_glycoprotein_score` column raises `ValueError` naming the file (this includes old-schema files). `RESULT_SUFFIX = ".surface_glyco.csv"` in both analysis modules.

- [ ] **Step 1: Write the failing tests**

```python
"""Reading result CSVs."""

import pytest

from surface_glyco.results import read_all, read_called

NEW = (
    "id,prediction,surface_glycoprotein_score\n"
    "A,surface_glycoprotein,0.9\nB,other,0.1\nC,surface_glycoprotein,0.6\n"
)


def _write(tmp_path, text):
    p = tmp_path / "r.csv"
    p.write_text(text)
    return p


def test_only_called_rows_are_returned(tmp_path):
    assert read_called(_write(tmp_path, NEW)) == [("A", 0.9), ("C", 0.6)]


def test_header_only_file_is_empty_not_an_error(tmp_path):
    assert read_called(_write(tmp_path, "id,prediction,surface_glycoprotein_score\n")) == []


def test_ids_with_commas_survive(tmp_path):
    text = 'id,prediction,surface_glycoprotein_score\n"sp|P1,x",surface_glycoprotein,0.5\n'
    assert read_called(_write(tmp_path, text)) == [("sp|P1,x", 0.5)]


def test_old_schema_is_rejected_with_the_file_name(tmp_path):
    p = _write(tmp_path, "id,prediction,probability_adhesion\nA,Adhesion,0.9\n")
    with pytest.raises(ValueError, match="r.csv"):
        read_called(p)


def test_read_all_returns_every_row(tmp_path):
    assert read_all(_write(tmp_path, NEW))[1] == ("B", "other", 0.1)
```
In `tests/kingdom_survey/test_join.py` and `tests/adhesion_properties/test_universe.py`: change every result file name from `.adhesion_predict.csv` to `.surface_glyco.csv`, every header `id,prediction,probability_adhesion` to `id,prediction,surface_glycoprotein_score`, and every label `Adhesion` to `surface_glycoprotein` in those fixtures. Add to `test_join.py` (use that file's existing import alias for `join`):
```python
def test_count_result_rows_counts_calls_not_all_rows(tmp_path):
    p = tmp_path / "r.surface_glyco.csv"
    p.write_text(
        "id,prediction,surface_glycoprotein_score\n"
        "A,surface_glycoprotein,0.9\nB,other,0.1\nC,other,0.2\n"
    )
    assert join.count_result_rows(p) == 1
    assert join.probability_stats(p) == (0.9, 0.9)


def test_results_are_discovered_by_the_new_suffix(tmp_path):
    (tmp_path / "Sp_one.surface_glyco.csv").write_text(
        "id,prediction,surface_glycoprotein_score\n"
    )
    (tmp_path / "Sp_old.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\n"
    )
    assert join.result_stems(tmp_path) == {"Sp_one"}
```
`join.result_stems(results_dir)` does not exist yet: extract line 162's set comprehension into it (see step 3). Add one case to `test_universe.py` calling `universe.read_result_ids` on a three-row file (two called, one `other`) and asserting `[("A", 0.9), ("C", 0.6)]`.

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_results.py -q`
Expected: FAIL (`ModuleNotFoundError: surface_glyco.results`).

- [ ] **Step 3: Implement**

`results.py`:
```python
"""Read prediction CSVs written by surface_glyco_predict."""

import csv
from pathlib import Path

SCORE_COLUMN = "surface_glycoprotein_score"
CALLED_LABEL = "surface_glycoprotein"


def read_all(path):
    """Return (id, label, score) for every row."""
    path = Path(path)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or SCORE_COLUMN not in reader.fieldnames:
            raise ValueError(f"{path}: expected a {SCORE_COLUMN!r} column (header: {reader.fieldnames})")
        return [(row["id"], row["prediction"], float(row[SCORE_COLUMN])) for row in reader]


def read_called(path):
    """Return (id, score) for rows called as surface glycoproteins."""
    return [(i, score) for i, label, score in read_all(path) if label == CALLED_LABEL]
```
`analysis/kingdom_survey/join.py`: `RESULT_SUFFIX = ".surface_glyco.csv"`; `from surface_glyco.results import read_called`; `count_result_rows` returns `len(read_called(result_csv_path))`; `probability_stats` computes mean and median over `[s for _, s in read_called(path)]` (keep the `nan, nan` return when empty); docstrings say "called proteins". Add
```python
def result_stems(results_dir: Path) -> set[str]:
    """Species stems that have a result file."""
    return {p.name[: -len(RESULT_SUFFIX)] for p in results_dir.glob(f"*{RESULT_SUFFIX}")}
```
and use it at the former line 162. `01_build_species_table.py:59`: `glob("*.surface_glyco.csv")`. `universe.py`: `RESULT_SUFFIX = ".surface_glyco.csv"`, `read_result_ids` returns `read_called(result_csv_path)`, docstring "(protein_id, score) for called proteins".

- [ ] **Step 4: Run all affected suites, one process each**

```bash
pytest tests/surface_glyco -q
pytest tests/kingdom_survey -q
pytest tests/adhesion_properties -q
```
Expected: PASS. These analysis modules now import `surface_glyco`; the package must be installed in the environment (CI does `pip install -e .`). Add one line to `analysis/kingdom_survey/README.md` and `analysis/adhesion_properties/README.md` (create if absent): "Requires `pip install -e .` from the repository root (imports `surface_glyco.results`)."

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "results: new-schema reader and discovery for the survey and property analyses

kingdom_survey counted every data row as a call, which overcounts since all
scores are written (PR #23). Result files are now found by .surface_glyco.csv.
Old-schema files are rejected with a clear error (not in circulation).

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Card framework (`card.py`) and card loading

**Files:**
- Create: `src/surface_glyco/card.py`
- Modify: `src/surface_glyco/model.py` (`load_model_card`), `src/surface_glyco/embeddings.py` (constants only)
- Test: `tests/surface_glyco/test_card.py` (new)

**Interfaces:**
- Produces, in `surface_glyco.card` (no torch import): `CARD_VERSION = 2`, `POOLING_RESIDUE_MEAN = "residue_mean"`, `KNOWN_POOLINGS = (POOLING_RESIDUE_MEAN,)`, `MAX_RESIDUES = 1022`, `DEFAULT_REPR_LAYER = 6`; `class ModelCardError(Exception)`; `@dataclass(frozen=True) class EmbeddingSettings: esm_model: str; repr_layer: int; pooling: str; max_residues: int`; `new_card(esm_model, repr_layer, pooling, n_positive, n_negative, **extra) -> dict`; `resolve_embedding_settings(card, cli_model_name, default_model) -> EmbeddingSettings`.
- `surface_glyco.model.load_model_card(path)` returns `None` when no card file exists, and raises `ModelCardError` (naming the file) when the JSON is malformed.

- [ ] **Step 1: Write the failing tests**

```python
"""Card contents and settings resolution."""

import pytest

from surface_glyco.card import (
    MAX_RESIDUES,
    POOLING_RESIDUE_MEAN,
    ModelCardError,
    new_card,
    resolve_embedding_settings,
)
from surface_glyco.model import load_model_card

DEFAULT = "esm2_t6_8M_UR50D"


def _card(**over):
    card = new_card(DEFAULT, 6, POOLING_RESIDUE_MEAN, n_positive=10, n_negative=20)
    card.update(over)
    return card


def test_new_card_has_v2_fields_and_claims_nothing_unmeasured():
    card = _card()
    assert card["card_version"] == 2
    assert (card["pooling"], card["max_residues"]) == (POOLING_RESIDUE_MEAN, MAX_RESIDUES)
    assert card["truncation"] == "truncate_at_max_residues"
    assert card["validation"] == []
    assert card["threshold"] == {"value": 0.5, "chosen_on": "default (uncalibrated)"}


def test_card_settings_are_used_when_cli_model_is_omitted():
    s = resolve_embedding_settings(_card(esm_model="esm2_t12_35M_UR50D"), None, DEFAULT)
    assert (s.esm_model, s.repr_layer, s.pooling) == ("esm2_t12_35M_UR50D", 6, POOLING_RESIDUE_MEAN)


def test_matching_cli_model_is_accepted():
    assert resolve_embedding_settings(_card(), DEFAULT, DEFAULT).esm_model == DEFAULT


@pytest.mark.parametrize(
    "card, field",
    [
        (_card(esm_model="esm2_t12_35M_UR50D"), "esm_model"),  # with cli DEFAULT below
        (_card(pooling="all_tokens_legacy"), "pooling"),
        (_card(pooling="max"), "pooling"),
        (_card(max_residues=2000), "max_residues"),
        (_card(repr_layer="six"), "repr_layer"),
    ],
)
def test_bad_cards_are_refused_naming_the_field(card, field):
    with pytest.raises(ModelCardError, match=field):
        resolve_embedding_settings(card, DEFAULT, DEFAULT)


def test_card_without_esm_model_is_refused():
    card = _card()
    del card["esm_model"]
    with pytest.raises(ModelCardError, match="esm_model"):
        resolve_embedding_settings(card, None, DEFAULT)


def test_missing_card_is_refused():
    with pytest.raises(ModelCardError, match="card"):
        resolve_embedding_settings(None, None, DEFAULT)


def test_malformed_card_json_names_the_file(tmp_path):
    model = tmp_path / "m.pkl"
    model.write_bytes(b"x")
    (tmp_path / "m.json").write_text("{not json")
    with pytest.raises(ModelCardError, match="m.json"):
        load_model_card(model)
    (tmp_path / "m.json").unlink()
    assert load_model_card(model) is None
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_card.py -q`
Expected: FAIL (`ModuleNotFoundError: surface_glyco.card`).

- [ ] **Step 3: Implement `card.py`**

```python
"""Model card contents and the rules for turning a card into embedding settings.

This module must not import torch or esm: predict and evaluate validate a model before
loading any weights, and these rules are tested without a GPU stack.
"""

from dataclasses import dataclass

CARD_VERSION = 2
POOLING_RESIDUE_MEAN = "residue_mean"
KNOWN_POOLINGS = (POOLING_RESIDUE_MEAN,)
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


def new_card(esm_model, repr_layer, pooling, n_positive, n_negative, **extra):
    """Card for a newly trained model. Fields that were not measured are left empty."""
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


def resolve_embedding_settings(card, cli_model_name, default_model):
    """Return the embedding settings a model must be scored with, or raise ModelCardError."""
    if card is None:
        raise ModelCardError(
            "card: this model has no card, so its embedding settings are unknown. "
            "Train a model with surface_glyco_train, which writes one."
        )
    esm_model = card.get("esm_model")
    if not esm_model:
        raise ModelCardError("esm_model: missing from the model card")
    if cli_model_name is not None and cli_model_name != esm_model:
        raise ModelCardError(
            f"esm_model: the model was trained on {esm_model} embeddings, "
            f"but --model-name is {cli_model_name}"
        )
    pooling = card.get("pooling")
    if pooling not in KNOWN_POOLINGS:
        raise ModelCardError(f"pooling: {pooling!r} is not supported; supported: {KNOWN_POOLINGS}")
    max_residues = card.get("max_residues", MAX_RESIDUES)
    if max_residues != MAX_RESIDUES:
        raise ModelCardError(
            f"max_residues: the model used {max_residues}, this build supports {MAX_RESIDUES}"
        )
    try:
        repr_layer = int(card.get("repr_layer", DEFAULT_REPR_LAYER))
    except (TypeError, ValueError):
        raise ModelCardError(f"repr_layer: {card.get('repr_layer')!r} is not an integer") from None
    return EmbeddingSettings(esm_model, repr_layer, pooling, max_residues)
```
`model.py`: in `load_model_card`, wrap `json.load` in `try/except json.JSONDecodeError as e: raise ModelCardError(f"{path}: not valid JSON ({e})") from e` and import `ModelCardError` from `surface_glyco.card`. `embeddings.py`: replace the local `DEFAULT_REPR_LAYER`, `MAX_RESIDUES` definitions with `from surface_glyco.card import DEFAULT_REPR_LAYER, MAX_RESIDUES` and keep `POOLING = "residue_mean"` as `from surface_glyco.card import POOLING_RESIDUE_MEAN as POOLING` only if `train.py` still imports `POOLING` (Task 7 removes that import). Import nothing that this task does not use.

- [ ] **Step 4: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS. (`test_bad_cards_are_refused_naming_the_field`'s first case relies on the CLI model being `DEFAULT` while the card says t12; confirm the match text contains `esm_model`.)

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "card: v2 framework and settings resolution; malformed cards are card errors

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Truncation count and skip/retry alignment in `embeddings.py`

**Files:**
- Modify: `src/surface_glyco/embeddings.py`
- Test: `tests/surface_glyco/test_embeddings.py`

**Interfaces:**
- Produces: `count_truncated(sequences) -> int`: how many sequences exceed `MAX_RESIDUES` after `sanitize_sequence`.
- Pins existing behaviour: when a sequence cannot be embedded, `get_esm_embeddings(..., return_indices=True)` returns the indices of the sequences that were embedded, so labels can be aligned.

- [ ] **Step 1: Write the failing tests** (the file has `_seqs` and `_cpu` helpers; add imports `import pytest` and `from surface_glyco import embeddings as emb`)

```python
def test_count_truncated_counts_sequences_over_the_limit():
    seqs = _seqs(("a", "M" * 1022), ("b", "M" * 1023), ("c", ""), ("d", "M*" * 1200))
    assert emb.count_truncated(seqs) == 2  # b, and d (1200 residues once '*' is stripped)


def test_failed_sequence_is_skipped_and_indices_stay_aligned(monkeypatch):
    real = emb._embed_batch

    def flaky(model, alphabet, bc, batch, layer, device):
        if any(seq == "BAD" for _, seq in batch):
            raise RuntimeError("boom")
        return real(model, alphabet, bc, batch, layer, device)

    monkeypatch.setattr(emb, "_embed_batch", flaky)
    seqs = _seqs(("a", "MKTAYIAK"), ("bad", "BAD"), ("c", "MKTAYIAKQR"))
    vecs, ids, kept = emb.get_esm_embeddings(seqs, device=_cpu(), batch_size=3, return_indices=True)
    assert ids == ["a", "c"] and kept == [0, 2] and len(vecs) == 2
```
(`"BAD"` is sanitized to `BAD`, a valid residue string; the monkeypatch makes any batch containing it raise, including its one-at-a-time retry.)

- [ ] **Step 2: Run to verify the first fails and the second passes or fails for the right reason**

Run: `pytest tests/surface_glyco/test_embeddings.py -q`
Expected: first FAIL (`count_truncated` missing). The second documents current behaviour; if it fails, stop and report: register item 2 would then be a real bug.

- [ ] **Step 3: Implement**

```python
def count_truncated(sequences):
    """Number of sequences longer than MAX_RESIDUES after sanitizing (they get truncated)."""
    return sum(len(sanitize_sequence(s["sequence"])) > MAX_RESIDUES for s in sequences)
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "embeddings: count_truncated; test that a failed sequence keeps indices aligned

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Card-driven `predict`, `evaluate`, `train`

**Files:**
- Modify: `src/surface_glyco/scripts/predict.py`, `evaluate.py`, `train.py`
- Test: `tests/surface_glyco/test_scripts_cards.py` (new); replace the dedupe test in `test_io_model_predict.py`

**Interfaces:**
- Consumes: `resolve_embedding_settings`, `ModelCardError`, `new_card`, `POOLING_RESIDUE_MEAN`, `DEFAULT_REPR_LAYER` (Task 5); `count_truncated` (Task 6); `model_filename`, `get_models_dir` (Task 3); `SCORE_COLUMN` (Task 2).
- Produces:
  - `predict.main(input_path, model_path, output_file, model_name=None, silent=False, show_all=False, max_workers=None)`; on `ModelCardError` prints `Error: <message>` to stderr and exits 1 before reading FASTA or embedding.
  - `--model-name` defaults to `None` on `predict` and `evaluate`.
  - `train.dedupe_sequences(sequences)` drops empty sequences, exact duplicates within a class (first id kept) and every copy of a sequence present in both classes.
  - `train.prepare_data` returns `(sequences, n_duplicates_removed, n_conflicting_removed)`.
  - A trained model's card is `new_card(...)` plus `n_duplicates_removed`, `n_conflicting_removed`, `n_not_embedded`, `training_sequences_sha256`, `positive_dir`, `negative_dir`, `holdout_accuracy`, `classifier`, and `environment` (`numpy`, `scikit_learn`, `torch`, `fair_esm` version strings, `"unknown"` if not installed).

- [ ] **Step 1: Write the failing tests**

```python
"""predict/evaluate/train honour the model card."""

import sys

import numpy as np
import pytest
from sklearn.linear_model import LogisticRegression

from surface_glyco.card import POOLING_RESIDUE_MEAN, new_card
from surface_glyco.model import save_model, save_model_card
from surface_glyco.scripts import predict as predict_mod
from surface_glyco.scripts import train as train_mod

ESM = "esm2_t6_8M_UR50D"


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
        return np.ones((len(sequences), 4)), [s["id"] for s in sequences]

    monkeypatch.setattr(predict_mod, "get_esm_embeddings", fake)
    return calls


def test_card_settings_reach_the_embedder(tmp_path, fasta, fake_embed):
    out = tmp_path / "o.csv"
    card = new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1)
    predict_mod.main(fasta, _model(tmp_path, card), out, None, silent=True)
    assert (fake_embed["model_name"], fake_embed["repr_layer"]) == (ESM, 6)
    assert out.read_text().splitlines()[0] == "id,prediction,surface_glycoprotein_score"


@pytest.mark.parametrize(
    "card, cli_model, field",
    [
        (new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1), "esm2_t12_35M_UR50D", "esm_model"),
        (new_card(ESM, 6, "all_tokens_legacy", 1, 1), None, "pooling"),
        (None, None, "card"),
    ],
)
def test_refusals_exit_1_naming_the_field_before_embedding(
    tmp_path, fasta, fake_embed, capsys, card, cli_model, field
):
    with pytest.raises(SystemExit) as e:
        predict_mod.main(fasta, _model(tmp_path, card), tmp_path / "o.csv", cli_model)
    assert e.value.code == 1
    assert field in capsys.readouterr().err
    assert fake_embed == {}  # nothing was embedded


def test_truncated_sequences_are_reported(tmp_path, fasta, fake_embed, capsys):
    card = new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1)
    predict_mod.main(fasta, _model(tmp_path, card), tmp_path / "o.csv", None, silent=True)
    assert "1 of 2 sequences" in capsys.readouterr().out


def test_cli_passes_model_and_default_model_name(tmp_path, fasta, monkeypatch):
    got = {}
    monkeypatch.setattr(predict_mod, "main", lambda *a, **k: got.update(args=a, kwargs=k))
    model = _model(tmp_path, new_card(ESM, 6, POOLING_RESIDUE_MEAN, 1, 1))
    monkeypatch.setattr(
        sys, "argv", ["surface_glyco_predict", "--input", str(fasta), "--model", str(model)]
    )
    predict_mod.cli()
    assert got["args"][3] is None  # model_name left to the card


def test_dedupe_drops_sequences_present_in_both_classes():
    seqs = [
        {"id": "a", "sequence": "MK", "label": 1},
        {"id": "b", "sequence": "MK", "label": 0},
        {"id": "c", "sequence": "MKT", "label": 1},
        {"id": "d", "sequence": "MKT", "label": 1},
    ]
    assert [s["id"] for s in train_mod.dedupe_sequences(seqs)] == ["c"]


def test_empty_sequences_are_not_trained_on():
    seqs = [{"id": "a", "sequence": "", "label": 1}, {"id": "b", "sequence": "MK", "label": 0}]
    assert [s["id"] for s in train_mod.dedupe_sequences(seqs)] == ["b"]


def test_end_to_end_with_real_embeddings(tmp_path):
    from surface_glyco.embeddings import get_esm_embeddings

    seqs = [
        {"id": f"s{i}", "sequence": "MKTAYIAKQRQISFVKSHFSRQ"[: 8 + i], "label": i % 2}
        for i in range(6)
    ]
    emb, _ = get_esm_embeddings(seqs, model_name=ESM)
    clf = LogisticRegression(max_iter=200).fit(emb, [s["label"] for s in seqs])
    model = tmp_path / "m.pkl"
    save_model(clf, model)
    save_model_card(model, new_card(ESM, 6, POOLING_RESIDUE_MEAN, 3, 3))
    fa = tmp_path / "in.fa"
    fa.write_text("".join(f">{s['id']}\n{s['sequence']}\n" for s in seqs))
    out = tmp_path / "o.csv"
    predict_mod.main(fa, model, out, None, silent=True)
    rows = out.read_text().splitlines()
    assert len(rows) == 7 and rows[0] == "id,prediction,surface_glycoprotein_score"
```
Delete the old `test_dedupe_sequences_keeps_first_id_per_label` from `test_io_model_predict.py` (its expectation, that a sequence stays in both classes, is the behaviour being removed). Add a train-side card test:
```python
def test_trained_model_card_has_environment_and_counts(tmp_path, monkeypatch):
    # prepare_data and embeddings are stubbed (40 records: train_classifier runs 5-fold CV);
    # this pins what train.main writes to the card
    from surface_glyco.model import load_model_card

    monkeypatch.setattr(train_mod, "prepare_data", lambda p, n: (
        [{"id": f"s{i}", "sequence": "MK" * (i + 2), "label": i % 2} for i in range(40)], 3, 2))
    monkeypatch.setattr(train_mod, "get_esm_embeddings", lambda seqs, **kw: (
        np.random.default_rng(1).normal(size=(len(seqs), 4)), [s["id"] for s in seqs],
        list(range(len(seqs)))))
    out = tmp_path / "m.pkl"
    train_mod.main(tmp_path, tmp_path, out, ESM, 0.2)
    card = load_model_card(out)
    assert card["card_version"] == 2 and card["pooling"] == POOLING_RESIDUE_MEAN
    assert (card["n_duplicates_removed"], card["n_conflicting_removed"]) == (3, 2)
    assert set(card["environment"]) == {"numpy", "scikit_learn", "torch", "fair_esm"}
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/surface_glyco/test_scripts_cards.py -q`
Expected: FAIL (`main` does not resolve a card; `dedupe_sequences` keeps both-class sequences; `prepare_data` returns two values).

- [ ] **Step 3: Implement `predict.py`**

New `main` head (after `load_model`):
```python
    classifier = load_model(model_path)
    try:
        settings = resolve_embedding_settings(load_model_card(model_path), model_name, DEFAULT_MODEL)
    except ModelCardError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
```
Resolve settings before reading any FASTA. After the sequences are loaded:
```python
    n_trunc = count_truncated(all_sequences)
    if n_trunc:
        print(f"{n_trunc} of {len(all_sequences)} sequences are longer than "
              f"{settings.max_residues} residues and will be truncated")
```
Call `get_esm_embeddings(all_sequences, model_name=settings.esm_model, repr_layer=settings.repr_layer)`. In `cli()`: `--model-name` default `None` (keep `choices`); default model path `get_models_dir() / model_filename(args.model_name or DEFAULT_MODEL)`; `main(args.input, args.model, args.output, args.model_name, silent=..., show_all=..., max_workers=...)`. `load_model_card` raising `ModelCardError` for malformed JSON must also be caught: wrap both calls in the same `try`.

- [ ] **Step 4: Implement `evaluate.py` and `train.py`**

`evaluate.py`: same `try/except ModelCardError` block, `--model-name default=None`, pass `repr_layer=settings.repr_layer` and `model_name=settings.esm_model` to `get_esm_embeddings`. `train.py`: remove the `POOLING` import (and the alias in `embeddings.py` if nothing else uses it); `dedupe_sequences`:
```python
def dedupe_sequences(sequences):
    """Drop empty sequences, exact duplicates within a class, and sequences in both classes.

    A sequence with both labels cannot be learned and leaks across any split, so every copy is
    dropped. Duplicates within one class keep their first id.
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
`prepare_data` computes `n_conflicting_removed = sum(len(labels_by_seq[s["sequence"]]) > 1 for s in all_sequences)` (records, not sequences) and `n_duplicates_removed = len(all_sequences) - len(unique) - n_conflicting_removed_empty_adjusted`; simplest and testable: `n_duplicates_removed = len(all_sequences) - len(unique) - n_conflicting_removed`, so that the two counts and the kept count sum to the input. It returns `(unique, n_duplicates_removed, n_conflicting_removed)`. `main` unpacks the three and builds the card with `new_card(model_name, DEFAULT_REPR_LAYER, POOLING_RESIDUE_MEAN, int(labels.sum()), int(len(labels) - labels.sum()), classifier=type(classifier).__name__, n_duplicates_removed=..., n_conflicting_removed=..., n_not_embedded=len(sequences) - len(kept), training_sequences_sha256=sequences_sha256(sequences), positive_dir=str(positive_dir), negative_dir=str(negative_dir), holdout_accuracy=float(test_acc), environment=_environment())` with
```python
def _environment():
    """Versions of the libraries a model's scores depend on."""
    from importlib.metadata import PackageNotFoundError, version

    out = {}
    for key, dist in (("numpy", "numpy"), ("scikit_learn", "scikit-learn"),
                      ("torch", "torch"), ("fair_esm", "fair-esm")):
        try:
            out[key] = version(dist)
        except PackageNotFoundError:
            out[key] = "unknown"
    return out
```
Add `a one-line assertion test` that `n_duplicates_removed + n_conflicting_removed + len(unique) == len(all_sequences)` using a small fixture in `test_scripts_cards.py`.

- [ ] **Step 5: Run tests**

Run: `pytest tests/surface_glyco -q`
Expected: PASS (about 48 tests).

- [ ] **Step 6: Manual check of the refusal path with a real, tiny model**

```bash
python - <<'E'
import numpy as np
from sklearn.linear_model import LogisticRegression
from surface_glyco.model import save_model
rng = np.random.default_rng(0)
save_model(LogisticRegression().fit(rng.normal(size=(20, 320)), [0, 1] * 10), "/tmp/nocard.pkl")
E
surface_glyco_predict --input tests/input_tests/Saccharomyces.pep --model /tmp/nocard.pkl --output /tmp/x.csv; echo "exit=$?"
```
Expected: `Error: card: this model has no card ...` on stderr, `exit=1`, no `/tmp/x.csv`.

- [ ] **Step 7: Commit**

```bash
pre-commit run --all-files
git add -A && git commit -m "predict/evaluate/train: use the model card; report truncation; stricter dedupe

Settings come from the card; a missing card or an unsupported pooling is refused
before any embedding. Training drops sequences present in both classes and empty
sequences and writes a v2 card with environment versions.

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Docs and decision records

**Files:**
- Modify: `CHANGELOG.md`, `README.md`, `AGENTS.md`, `docs/TOOL-ARCHITECTURE.md` (section 2.0 and line ~47 "keep old entry points as deprecated aliases"), `docs/PLAN-2026-09-30-pipeline-and-decisions.md`, `analysis/chytrid_batrach/NOTES.md`, `docs/plans/2026-09-30-model-bundle-issue-9.md` (banner)
- Test: `tests/surface_glyco/test_docs.py` (new)

- [ ] **Step 1: Write the failing test**

```python
"""Docs and changelog agree with the code."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_changelog_lists_the_breaking_rename_under_unreleased():
    text = (ROOT / "CHANGELOG.md").read_text()
    unreleased = text.split("## [Unreleased]")[1].split("\n## [")[0]
    assert "BREAKING" in unreleased and "surface_glyco" in unreleased


def test_readme_uses_current_names_only():
    text = (ROOT / "README.md").read_text()
    assert "adhesion_train" not in text and "adhesion_predict " not in text
    assert "surface_glyco_train" in text


def test_agents_does_not_claim_a_thread_safe_cache():
    assert "thread-safe" not in (ROOT / "AGENTS.md").read_text().lower()


def test_decisions_doc_records_the_2026_09_30_changes():
    text = (ROOT / "docs" / "PLAN-2026-09-30-pipeline-and-decisions.md").read_text()
    assert "Changes made later on 2026-09-30" in text
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/surface_glyco/test_docs.py -q`
Expected: FAIL (4).

- [ ] **Step 3: Edit the docs**

`CHANGELOG.md`, under `## [Unreleased]` add:
```markdown
### Changed
- **BREAKING:** the package `adhesion_predict` is now `surface_glyco`. Entry points: `surface_glyco_predict`, `surface_glyco_train`, `surface_glyco_evaluate`. Output column `probability_adhesion` is now `surface_glycoprotein_score`; labels `Adhesion` / `Non-adhesion` are now `surface_glycoprotein` / `other`; default output files end in `.surface_glyco.csv`; model files are `surface_glyco_model_<esm>.pkl`. No aliases, and result files written by the old package are no longer read by the analysis code. Reason: the model scored FLO/ALS-like surface glycoproteins, not adhesins.
- **BREAKING:** the two models that shipped with 0.1.0 are removed. They were trained on FLO/ALS homologs versus random proteins with padding-inclusive pooling (issue #25). No model ships until a surface-glycoprotein model is trained and validated.
- `predict` and `evaluate` take the ESM model, layer and pooling from the model's card and stop with a message naming the field on a mismatch, a missing card or an unsupported pooling.
- A `./models` directory in the working directory no longer replaces the packaged models. Set `SURFACE_GLYCO_MODELS_DIR` to use another directory.
- `requires-python` is now `>=3.11`.
- `kingdom_survey` counts called proteins instead of all rows in a result file, and finds `*.surface_glyco.csv`.

### Added
- `predict` reports how many sequences are truncated at 1022 residues. Cards of new models record library versions.

### Fixed
- Training drops sequences that occur in both classes and empty sequences.
```
and append at the end of the file (so a later `bump minor` has a base release):
```markdown
## [0.1.0] - 2026-02-16

Initial tagged version (`v0.1.0`).
```
`README.md`: replace old names with the new ones, state that no model ships yet and how to train one, and keep the scope paragraph from PR #28. `AGENTS.md`: replace the thread-safety claim with "the model cache is a module-level dict and is not thread-safe". `docs/TOOL-ARCHITECTURE.md` section 2.0: rename the table heading to "Names (implemented, unreleased)", mark rows done, and change the "keep old entry points as deprecated aliases" row to "no aliases (decision 4)". `analysis/chytrid_batrach/NOTES.md`: add "These scripts read frozen result files written by `adhesion_predict` 0.1.0 (`*.adhesion_predict.csv`, column `probability_adhesion`) and were not updated." `docs/plans/2026-09-30-model-bundle-issue-9.md`: extend the banner with "The legacy-pooling emulation, `--allow-legacy-pooling` and the cards for the old pickles were dropped on 2026-09-30 (later): the old models were deleted instead."

`docs/PLAN-2026-09-30-pipeline-and-decisions.md`: add a section **"Changes made later on 2026-09-30"** after section 4:
```markdown
## 4a. Changes made later on 2026-09-30

| Was | Now | Reason |
|---|---|---|
| ESM + LR model frozen as legacy; pickles kept, legacy pooling emulated behind a flag (decision 2) | The old pickles are **deleted**. No legacy mode, no `--allow-legacy-pooling`. | The tool is not in circulation. The emulation of the original padding-inclusive pooling could only be approximate and could not be verified. |
| Old result files stay readable by downstream code | Old-schema files are **not** read; result discovery uses `*.surface_glyco.csv`. | Not in circulation; the Fungi_5k outputs came from the pre-fix model and are already due for a re-run (#16). |
| Step 1 is a SignalP + GPI + Ser/Thr rule (decision 1); ESM model only a possible step 2 base | Step 1 is **two candidates**: the rule, and an ESM-based model **retrained on a surface-glycoprotein label** (surface vs non-surface). Both are tested on the same GO truth set (decision 3); rule, ML or hybrid is decided from that comparison. | Owner is open to retraining on the correct target. This also measures the rule-vs-ESM agreement (R3). |
| 0.2.0 ships with the rename | **0.2.0 is cut when a validated model ships.** The rename and card framework merge to `main` earlier, unreleased. | The package has no useful default model until the new one exists. |
| Simple retrain on the current labels in 0.2.0 | **Skipped.** | Replaced by the surface-glycoprotein-label model. |

The surface-glycoprotein-label model needs its own spec (positive and negative definitions,
homology-grouped CV, GO truth set, head-to-head with the rule) and independent review before
any plan or code.
```
Also change section 7 item 1 of that file to "Done in the rename PR" only once merged; in this PR write "In progress (PR for `rename-plan`)".

- [ ] **Step 4: Run tests and commit**

```bash
pytest tests/surface_glyco -q
pre-commit run --all-files
git add -A && git commit -m "docs: changelog (breaking), README, AGENTS, architecture, decision record

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```
Expected: PASS.

---

### Task 9: Whole-branch verification

**Files:** none modified unless a check fails.

- [ ] **Step 1: Stale-name scan**

Run: `git grep -n "adhesion_predict\|adhesion_train\|adhesion_evaluate\|adhesion_model_\|probability_adhesion" -- . ':!docs/superpowers/plans/2026-09-04-*' ':!docs/superpowers/specs/2026-09-04-*' ':!Changes.md' ':!docs/model-review' ':!analysis/chytrid_batrach' ':!analysis/embedding_clustering/REPORT.md' ':!CHANGELOG.md' ':!docs/PLAN-2026-09-30-pipeline-and-decisions.md' ':!docs/plans' ':!docs/superpowers/plans/2026-09-30-*' ':!docs/TOOL-ARCHITECTURE.md'`
Expected: no output. A hit is fixed, or added to the exclusion list with a reason in the PR description.

- [ ] **Step 2: All test folders, one process each**

```bash
for d in surface_glyco kingdom_survey adhesion_properties embedding_clustering; do pytest tests/$d -q || exit 1; done
```
Expected: all PASS.

- [ ] **Step 3: Build a wheel and check its contents**

```bash
pip wheel . --no-deps -w /tmp/wheel_check -q
python - <<'E'
import glob, zipfile
w = glob.glob("/tmp/wheel_check/surface_glyco-*.whl")[0]
names = zipfile.ZipFile(w).namelist()
assert "surface_glyco/models/README.md" in names, "models/README.md missing"
assert not any(n.endswith(".pkl") for n in names), "a pickle is in the wheel"
assert not any(n.startswith("adhesion_predict/") for n in names)
assert "surface_glyco/card.py" in names and "surface_glyco/results.py" in names
print("wheel ok:", w)
E
```
Expected: `wheel ok`. A `0.1.dev...` or `0+unknown` version in the file name is acceptable on an untagged branch.

- [ ] **Step 4: Push, open the PR, wait for CI**

```bash
git push -u origin <branch>
gh pr create --base main --title "Rename adhesion_predict to surface_glyco; model-card framework (unreleased)" --body "<summary, breaking list from CHANGELOG, Review Focus results>"
gh pr checks --watch
```
Expected: lint, unit tests and the three analysis jobs pass. Do **not** run the `version_bump` workflow.

- [ ] **Step 5: Close out the pooling issue**

Comment on issue #25: the old pickles were deleted instead of emulated, with the PR link; leave the issue open until the new model replaces them, or close it with the owner's agreement.

---

## Self-review

**Spec coverage.** Rename (decision 4, names): Tasks 1-3. Old models removed, no legacy mode (changes of 2026-09-30 later): Tasks 3, 5, 7. Card-driven settings with refusal naming the field: Tasks 5, 7. Result reading, discovery and call counting: Task 4. Both-class and empty sequences: Task 7. Truncation count: Tasks 6-7. Failed-batch alignment (register item 2): Task 6. `ESM2_MODELS` (27): Task 3. `./models` in cwd (31): Task 3. Thread-safe claim (33): Task 8. CI path (30): Task 1. Python floor / `tomllib`: Task 1. Decision record: Task 8. Release: deliberately not in this plan. Not covered, by design: the glycoprotein-label model, the step 1 rule, the orchestrator, the GO truth set (separate projects, each needs a spec).

**Placeholders.** None intended. Two items depend on the reader's checkout: the existence of `src/.../models/adhesion_model_esm2_t6_8M_UR50D.pkl` (untracked on some checkouts; `rm -f` handles both) and `<branch>` in Task 9 step 4.

**Type consistency.** `EmbeddingSettings` has four fields (`esm_model`, `repr_layer`, `pooling`, `max_residues`) in Tasks 5 and 7; the `implied` field of revision 1 is gone. `SCORE_COLUMN`/`LABEL_*` (Task 2) match `results.SCORE_COLUMN`/`CALLED_LABEL` (Task 4) and the strings in tests. `prepare_data` returns three values in Task 7 only. `result_stems` is introduced in Task 4 and used only there.

**Not verified.** The code blocks have been read against the real code by a reviewer who ran the revision-1 versions; revision 2 has not been run. The GitHub Actions run is the first full check. `bump minor` has never been run; this plan avoids it.

## What changed from revision 1

Removed (owner decisions of 2026-09-30, later): legacy pooling mode, `--allow-legacy-pooling`, cards for the old pickles, old-schema result reader, the 0.2.0 release steps, the pickle-checksum test. Added from the independent review of revision 1: B1 `count_truncated` test fixed; B4 result discovery and suffix (three files); S1 environment-independent old-package test; S2 `requires-python >=3.11`; S3 `02_cv_and_proteome_eval.py` features path; S4 card edge cases (missing `esm_model`, bad `repr_layer`, malformed JSON); S5 `cli()` test; S6 README note on `pip install -e .`; S7 `CITATION.cff` version lines; S8 import only what each task uses; register item 2 retry test; Task 1 expected-failure text corrected.
