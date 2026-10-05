# cellsurface_sorting_hat core engine: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the part of `cellsurface_sorting_hat` that needs no GPU and no HPCC: FASTA checks, taxonomy-aware module status, the three-valued category engine, cache and provenance helpers, outputs and a command line driver that reads module tables from a work directory.

**Architecture:** A new package `src/cellsurface_sorting_hat/`, separate from `surface_glyco`. The driver reads per-module result tables and status sources from `<workdir>`, evaluates the rules in a packaged `categories.yaml` with Kleene logic, and writes `calls.long.tsv.gz`, `calls.wide.tsv.gz`, `report.md` and `run.json`. Wrappers that produce the module tables (SignalP, PredGPI, Pfam, repeat detectors, allergen homology, antigen lookup) are a second plan (section "Plan 2").

**Tech Stack:** Python 3.11+ (tests run with `/usr/bin/python3.12`), PyYAML, the standard library (`csv`, `gzip`, `hashlib`, `json`), pytest, ruff. No pandas, torch or GPU is imported by this package.

**Spec:** `docs/superpowers/specs/2026-10-04-orchestrator-design.md` (revision 5). Reviews: `...-review-1.md`, `...-review-2.md`.

## Global Constraints

Values are copied from the spec.

- `--taxon` or `--taxon-map` is required. `--taxon-map` overrides `--taxon` per protein. A protein with no taxon is an error. (D7)
- Applicability and status are separate. A module that is not applicable gives `not_assessable`. An applicable module with status `unvalidated` still gives `called` or `not_called`. (D7, revised)
- Status values: `estimated`, `smoke`, `unvalidated`. A status applies to a taxon only if the taxon is a tested taxon or a descendant of one. A shared broad label is not enough. Most specific tested taxon wins.
- A `status_source` is refused (status `unvalidated`) if its module name, version, `params_hash` or `artefact_hash` differs from the running module.
- Call values: `called`, `not_called`, `not_assessable`. Kleene three-valued logic. The value column never holds the status. Status is in a separate `_status`/`status` field.
- Default gate: `step1_rule@R0` (D10). It is a config item. Ungated calls are always written.
- Antigen call: combined `percentile` (all proteins) at most 15 (D12). Allergen call: `identity` at least 70 and `coverage` at least 80, or an allergen-specific Pfam hit; no surface gate (D6).
- `other_*` is computed over the categories that are not unknown for that protein, and `other_basis` lists those left out. It is unknown only when `surface_glycoprotein[v]` is unknown.
- Only available step 1 variants get records. Unavailable variants are listed in the report.
- Output tables are `.tsv.gz`. Writes are atomic (temporary file, then rename) with a `.sha256` sidecar. Cache tables are locked with `fcntl` for parallel runs; a damaged cache table is a cache miss.
- Outputs: `calls.long.tsv.gz`, `calls.wide.tsv.gz`, `evidence.tsv.gz`, `proteins.tsv.gz`, `report.md`, `run.json`.
- Input errors name the file (and the line for the taxon map) and exit with code 2. Program errors are not caught.
- Cache key covers module name, version, `params_hash`, `artefact_hash` and tool versions. Results are cached per sequence sha256.
- Python files: ruff `line-length = 100`, rules from `pyproject.toml`. CI pins ruff 0.3.5. Run the repo's pinned ruff (pre-commit) before committing.
- Report and docs text: Simplified Technical English. No claims that are not measured.
- Do not use `$(dirname "${BASH_SOURCE[0]}")` in any script that can run under SLURM. (No SLURM script is in this plan.)

## Module output contract

Plan 2 wrappers and the tests in this plan both rely on this contract.

- `<workdir>/modules/<module>.tsv.gz` (always gzip): tab-separated, header row. Required columns: `id` (the FASTA ID) and `state`; a table without them, or with a duplicate `id`, stops the run. Other columns are module fields. `state` is `ok` or a non-ok value (`not_applicable`, `not_in_reference`, `na_too_short`, `na_window`, `na_invalid`, `error`). Only `state == ok` rows are used.
- A call module has a `call` field with `called` or `not_called`. A flag module has a `hit` field with `1` or `0`. A number module has numeric fields (for example `percentile`, `identity`, `coverage`, `aligned_length`). An empty or non-numeric number is unknown. A `call` value that is not `called`/`not_called`, or a `hit` that is not `0`/`1`, makes the row state `bad_value` (unknown, counted in the report).
- `<workdir>/modules/<module>.json`: run record with `module`, `version`, `params_hash`, `artefact_hash` and optional `run_state` (`ok`, `partial`, `unavailable`, `not_run`, `error`; `unavailable`, `not_run` and `error` make the module absent). A run record without a table is listed in the report with its state (`error` if it gives none). A table in which no ID matches a FASTA ID makes the module `error`.
- `<workdir>/status/<module>.json`: status source with `module`, `version`, `params_hash`, `artefact_hash` and `entries`: a list of `{"taxa": [taxon IDs], "status": "estimated"|"smoke"|"unvalidated", "source": "path", "measure": {...}}`. `measure` is optional. It is where Plan 2 records the calibration: `calibration_set` (text), `truth_source` (path), `n_pos`, `n_neg` (non-negative integers), `sensitivity` and `specificity` (each `{"value", "lo", "hi"}` with 0 <= lo <= value <= hi <= 1) and `notes`. An unknown key in an entry or in `measure` is an error. The report prints a "Module calibration" table (one row per module and run taxon; "not measured" where there is no `measure`) and `run.json` carries the same rows.
- All text inputs may start with a UTF-8 byte order mark. A module row whose field count differs from the header, a truncated gzip file, invalid JSON and a run record that is not a JSON object stop the run with exit code 2 and a message that names the file (and the line).
- A protein that is in the FASTA but not in a module table gets a row with `state: error`, and the module's run state becomes `partial`.
- Modules used by the packaged `categories.yaml`: `step1_rule@R0`, `step1_rule@R1`, `step1_rule@R2`, `step1_ml@card` (call); `repeat02`, `repeat14` (call); `pfam_adhesion`, `pfam_allergen` (field `hit`); `antigen_lookup` (field `percentile`); `allergen_homology` (fields `identity`, `coverage`, `aligned_length`, `allergen_name`; a protein with no hit has `0`, `0`, `0`, state `ok`).
- Evidence fields copied to `evidence.tsv.gz` (list in `categories.yaml`): `antigen_lookup.{percentile,antigenicity,specificity,prevalence,max_crossreact}`, `allergen_homology.{allergen_name,identity,coverage,aligned_length}`, `cys_rich.tier`, `expression.log2fc`. A module that is absent contributes none.

## Clarifications of the spec made by this plan

1. Call names. The ungated form of `cell_wall_adhesion_candidate` is `cell_wall_adhesion_ungated`. The gated forms keep the spec names with the variant label, for example `cell_wall_adhesion_candidate[R0]`.
2. Status of a call. It is the weakest status among the modules that **decided** the result. For an AND that is false: the false inputs. For an AND that is true: all inputs. For an OR that is true: the true inputs. For an OR that is false: all inputs. Unknown inputs never contribute. The spec text for `not_called` ("all required modules") is replaced by this rule (the plan's Task 7 edits the spec).
3. `--taxdump` takes the path of `nodes.dmp` from the NCBI taxonomy dump. Its sha256 is written to the report.
4. The long table column is `other_basis` (the spec text said `basis`). The `categories.yaml` hash is written to `report.md` and `run.json`, not to the calls files.
5. The golden test checks hand-derived values for named proteins and calls. It does not compare a stored `calls` file. Review of the plan flagged this difference from spec 4 item 2; a stored file can be added in Plan 2 when real module outputs exist.
6. `allergen_homolog_hit` (identity >= 35%, aligned length >= 80 aa) is a call in `categories.yaml`. It is evidence and is not a mechanism category.
7. The module tables are always `.tsv.gz`. Plain and `.zst` readers (spec 3.6) apply to the FASTA input only.

## Review Focus

Input classes the spec implies and the tests pin. Each has a test in the task named.

1. A Windows FASTA (CRLF line ends, lowercase residues, trailing `*`) must be read, not rejected. (Task 7, `test_crlf_lowercase_trailing_stop_fasta_is_read`; Task 3)
2. A `--taxon` that is not in the taxonomy stops the run with a message that names the taxon, exit code 2, also when no status file exists. (Task 7, `test_unknown_taxon_stops_the_run`, `test_an_unknown_taxon_stops_the_run_even_with_no_status_file`)
3. A module table that lacks some proteins must make those proteins unknown (not false) and the module `partial` in the report. (Task 7, `test_partial_module_output_is_reported_and_missing_proteins_are_unknown`)
4. Two proteins with identical sequences and different IDs must get identical calls, and the report must warn when a module gives them different rows. (Task 7, `test_identical_sequences_with_different_ids_get_identical_calls`, `test_identical_sequences_with_different_module_rows_are_flagged`; Task 3)
5. A run with no module output at all must finish, write a report with a warning about the default gate, and write no per-variant records. (Task 7, `test_no_module_at_all_still_writes_a_report_with_a_warning`)
6. A module table with an unusable shape or content must be refused or counted, never read silently: no `state` column, duplicate rows, values that are not call values, IDs that match no protein. (Task 7, `test_module_table_without_a_state_column_stops_the_run`, `test_duplicate_rows_in_a_module_table_stop_the_run`, `test_a_bad_call_value_is_unknown_and_counted`, `test_a_module_whose_ids_match_no_protein_is_an_error`)
7. Two runs that write the same cache table at once must not lose rows or leave an unreadable table. (Task 4, `test_parallel_updates_lose_no_rows_and_leave_a_readable_table`)

---

### Task 0: Branch

**Files:** none.

- [ ] **Step 1: Create the working branch**

```bash
git fetch origin
git checkout -b sorting-hat-core origin/main   # after the spec PR is merged; before that, branch from the spec PR branch tip
git branch --show-current                        # expected: sorting-hat-core
mkdir -p src/cellsurface_sorting_hat tests/cellsurface_sorting_hat
```

- [ ] **Step 2: Tell isort the new package is first-party**

Without this, the pinned ruff reports `I001` for every file of this plan. In `pyproject.toml`, change the line in `[tool.ruff.lint.isort]` to:

```toml
known-first-party = ["surface_glyco", "cellsurface_sorting_hat"]
```

```bash
git add pyproject.toml
git commit -m "build(sorting-hat): ruff isort knows the new package

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

The package directory is found by `[tool.setuptools.packages.find] where = ["src"]` in `pyproject.toml`. Task 7 adds the package data and the script entry.

### Task 1: Kleene logic

**Files:**
- Create: `src/cellsurface_sorting_hat/__init__.py`, `src/cellsurface_sorting_hat/logic.py`
- Test: `tests/cellsurface_sorting_hat/test_logic.py`

**Interfaces:**
- Produces: `CALLED`, `NOT_CALLED`, `NOT_ASSESSABLE` (strings `called`, `not_called`, `not_assessable`); `k_not(a)`, `k_and(*values)`, `k_or(*values)`; each raises `ValueError` for a value that is not a call value.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_logic.py` with exactly this content:

```python
"""The Kleene tables of spec section 3.3, every case."""

import itertools

import pytest

from cellsurface_sorting_hat.logic import (
    CALLED as T,
)
from cellsurface_sorting_hat.logic import (
    NOT_ASSESSABLE as U,
)
from cellsurface_sorting_hat.logic import (
    NOT_CALLED as F,
)
from cellsurface_sorting_hat.logic import (
    k_and,
    k_not,
    k_or,
)

AND = {
    (T, T): T,
    (T, F): F,
    (T, U): U,
    (F, T): F,
    (F, F): F,
    (F, U): F,
    (U, T): U,
    (U, F): F,
    (U, U): U,
}
OR = {
    (T, T): T,
    (T, F): T,
    (T, U): T,
    (F, T): T,
    (F, F): F,
    (F, U): U,
    (U, T): T,
    (U, F): U,
    (U, U): U,
}


@pytest.mark.parametrize("a,b", list(itertools.product([T, F, U], repeat=2)))
def test_and_or_tables(a, b):
    assert k_and(a, b) == AND[(a, b)]
    assert k_or(a, b) == OR[(a, b)]


@pytest.mark.parametrize("a,expected", [(T, F), (F, T), (U, U)])
def test_not(a, expected):
    assert k_not(a) == expected


def test_a_true_input_keeps_an_or_true_when_another_input_is_unknown():
    assert k_or(U, T) == T


def test_a_false_input_makes_an_and_false_when_another_input_is_unknown():
    assert k_and(U, F) == F


def test_empty_and_or():
    assert k_and() == T
    assert k_or() == F


def test_rejects_values_that_are_not_call_values():
    with pytest.raises(ValueError):
        k_and("yes", T)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_logic.py -q`
Expected: collection error, `ModuleNotFoundError: No module named 'cellsurface_sorting_hat.logic'` (Task 0 made the package directory; `__init__.py` is created in Step 3)

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/__init__.py` with exactly this content:

```python
"""cellsurface_sorting_hat - sort fungal proteins into cell surface categories."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("surface_glyco")
except PackageNotFoundError:  # a checkout that was not pip-installed under this name
    __version__ = "0+unknown"
```

Create `src/cellsurface_sorting_hat/logic.py` with exactly this content:

```python
"""Kleene three-valued logic for category calls.

Values are the strings used in every output file: ``called`` (true), ``not_called`` (false) and
``not_assessable`` (unknown).
"""

CALLED = "called"
NOT_CALLED = "not_called"
NOT_ASSESSABLE = "not_assessable"
VALUES = (CALLED, NOT_CALLED, NOT_ASSESSABLE)


def _check(*values):
    for v in values:
        if v not in VALUES:
            raise ValueError(f"not a call value: {v!r}")


def k_not(a):
    _check(a)
    if a == CALLED:
        return NOT_CALLED
    if a == NOT_CALLED:
        return CALLED
    return NOT_ASSESSABLE


def k_and(*values):
    """False if any input is false; else unknown if any is unknown; else true."""
    _check(*values)
    if NOT_CALLED in values:
        return NOT_CALLED
    if NOT_ASSESSABLE in values:
        return NOT_ASSESSABLE
    return CALLED


def k_or(*values):
    """True if any input is true; else unknown if any is unknown; else false."""
    _check(*values)
    if CALLED in values:
        return CALLED
    if NOT_ASSESSABLE in values:
        return NOT_ASSESSABLE
    return NOT_CALLED
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_logic.py -q`
Expected: `16 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-core
git add src/cellsurface_sorting_hat/__init__.py src/cellsurface_sorting_hat/logic.py tests/cellsurface_sorting_hat/test_logic.py
git commit -m "feat(sorting-hat): Kleene three-valued logic

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 2: Taxonomy and module status

**Files:**
- Create: `src/cellsurface_sorting_hat/taxonomy.py`, `src/cellsurface_sorting_hat/status.py`, `tests/cellsurface_sorting_hat/conftest.py`
- Test: `tests/cellsurface_sorting_hat/test_taxonomy_status.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `Lineage.from_nodes_dmp(path)`, `Lineage.ancestors(taxon) -> list[int]`, `Lineage.is_descendant_or_self(taxon, ancestor) -> bool`, `Lineage.depth(taxon) -> int`, `TaxonError`.
- Produces: `ModuleIdentity(name, version, params_hash, artefact_hash)`, `StatusEntry(taxa, status, source, measure)`, `StatusRecord(identity, entries)`, `load_status_source(path) -> StatusRecord` (validates the optional `measure` object), `resolve_entry(record_or_None, running_identity, taxon, lineage) -> (entry_or_None, tested_taxon, reason)`, `resolve_status(...) -> (status, basis)`, `weakest(statuses) -> str`, constants `ESTIMATED`, `SMOKE`, `UNVALIDATED`.
- Produces fixtures: `nodes_dmp` (toy taxonomy file path), `write_module(workdir, name, rows, meta=None, status=None)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/conftest.py` with exactly this content:

```python
"""Shared fixtures: a toy taxonomy and a helper that writes module tables."""

import csv
import gzip
import json

import pytest

# taxid, parent. 1 root; 10 Fungi; 20 Eurotiomycetes; 30 Eurotiales; 31 Onygenales;
# 40 Aspergillus fumigatus (in 30); 41 Coccidioides immitis (in 31); 42 Aspergillus nidulans (in 30)
TOY_NODES = {1: 1, 10: 1, 20: 10, 30: 20, 31: 20, 40: 30, 41: 31, 42: 30}


@pytest.fixture
def nodes_dmp(tmp_path):
    path = tmp_path / "nodes.dmp"
    path.write_text("".join(f"{t}\t|\t{p}\t|\tno rank\t|\n" for t, p in TOY_NODES.items()))
    return path


@pytest.fixture
def write_module():
    def _write(workdir, name, rows, meta=None, status=None):
        folder = workdir / "modules"
        folder.mkdir(parents=True, exist_ok=True)
        columns = sorted({k for r in rows for k in r} - {"id"})
        with gzip.open(folder / f"{name}.tsv.gz", "wt") as fh:
            writer = csv.DictWriter(fh, ["id"] + columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
        record = {"module": name, "version": "1", "params_hash": "p", "artefact_hash": "a"}
        record.update(meta or {})
        (folder / f"{name}.json").write_text(json.dumps(record))
        if status is not None:
            sfolder = workdir / "status"
            sfolder.mkdir(exist_ok=True)
            (sfolder / f"{name}.json").write_text(
                json.dumps(
                    {
                        "module": name,
                        "version": "1",
                        "params_hash": "p",
                        "artefact_hash": "a",
                        "entries": status,
                    }
                )
            )

    return _write
```

Create `tests/cellsurface_sorting_hat/test_taxonomy_status.py` with exactly this content:

```python
import json

import pytest

from cellsurface_sorting_hat.status import (
    ModuleIdentity,
    StatusEntry,
    StatusRecord,
    load_status_source,
    resolve_entry,
    resolve_status,
    weakest,
)
from cellsurface_sorting_hat.taxonomy import Lineage, TaxonError

IDENT = ModuleIdentity("step1_rule@R0", "1", "p", "a")


@pytest.fixture
def lineage(nodes_dmp):
    return Lineage.from_nodes_dmp(nodes_dmp)


def test_ancestors_and_depth(lineage):
    assert lineage.ancestors(41) == [41, 31, 20, 10, 1]
    assert lineage.depth(1) == 0
    assert lineage.depth(41) == 4


def test_unknown_taxon_raises(lineage):
    with pytest.raises(TaxonError):
        lineage.ancestors(999)


def test_descendant_or_self(lineage):
    assert lineage.is_descendant_or_self(40, 30)
    assert lineage.is_descendant_or_self(40, 40)
    assert not lineage.is_descendant_or_self(30, 40)


def test_weakest_status():
    assert weakest(["estimated", "smoke"]) == "smoke"
    assert weakest(["estimated", "unvalidated", "smoke"]) == "unvalidated"
    assert weakest([]) == "unvalidated"
    with pytest.raises(ValueError):
        weakest(["good"])


def _record(entries, ident=IDENT):
    return StatusRecord(ident, tuple(StatusEntry(taxa, status) for taxa, status in entries))


def test_status_applies_to_a_tested_taxon(lineage):
    rec = _record([((40,), "estimated")])
    assert resolve_status(rec, IDENT, 40, lineage) == ("estimated", "taxon:40")


def test_status_does_not_pass_to_a_sibling_clade_with_the_same_broad_label(lineage):
    # tested on A. fumigatus (Eurotiales); Coccidioides is Onygenales, also in Eurotiomycetes
    rec = _record([((40,), "estimated")])
    assert resolve_status(rec, IDENT, 41, lineage) == ("unvalidated", "taxon not tested")


def test_status_passes_to_descendants_of_a_tested_taxon(lineage):
    rec = _record([((30,), "smoke")])
    assert resolve_status(rec, IDENT, 42, lineage) == ("smoke", "taxon:30")


def test_most_specific_tested_taxon_wins(lineage):
    rec = _record([((20,), "smoke"), ((40,), "estimated")])
    assert resolve_status(rec, IDENT, 40, lineage)[0] == "estimated"
    assert resolve_status(rec, IDENT, 42, lineage) == ("smoke", "taxon:20")


@pytest.mark.parametrize("field", ["name", "version", "params_hash", "artefact_hash"])
def test_stale_status_source_is_refused(lineage, field):
    rec = _record([((40,), "estimated")])
    running = ModuleIdentity(**{**IDENT.__dict__, field: "changed"})
    status, basis = resolve_status(rec, running, 40, lineage)
    assert status == "unvalidated"
    assert field in basis


def test_no_status_source(lineage):
    assert resolve_status(None, IDENT, 40, lineage) == ("unvalidated", "no status_source")


def test_load_status_source(tmp_path):
    path = tmp_path / "s.json"
    path.write_text(
        json.dumps(
            {
                "module": "m",
                "version": "2",
                "params_hash": "p",
                "artefact_hash": "a",
                "entries": [{"taxa": [40, 42], "status": "estimated", "source": "x"}],
            }
        )
    )
    rec = load_status_source(path)
    assert rec.identity == ModuleIdentity("m", "2", "p", "a")
    assert rec.entries == (StatusEntry((40, 42), "estimated", "x", None),)


def test_load_status_source_rejects_unknown_status(tmp_path):
    path = tmp_path / "s.json"
    path.write_text(
        json.dumps(
            {
                "module": "m",
                "version": "2",
                "params_hash": "p",
                "artefact_hash": "a",
                "entries": [{"taxa": [40], "status": "validated"}],
            }
        )
    )
    with pytest.raises(ValueError):
        load_status_source(path)


MEASURE = {
    "calibration_set": "S1:all",
    "truth_source": "phasec/metrics.json",
    "n_pos": 232,
    "n_neg": 4244,
    "sensitivity": {"value": 0.603, "lo": 0.55, "hi": 0.65},
    "specificity": {"value": 0.963, "lo": 0.95, "hi": 0.97},
    "notes": "R0",
}


def _status_file(tmp_path, measure, extra_entry=None):
    entry = {"taxa": [40], "status": "estimated", "source": "phasec"}
    if measure is not None:
        entry["measure"] = measure
    entry.update(extra_entry or {})
    path = tmp_path / "s.json"
    path.write_text(
        json.dumps(
            {
                "module": "m",
                "version": "1",
                "params_hash": "p",
                "artefact_hash": "a",
                "entries": [entry],
            }
        )
    )
    return path


def test_measure_with_sensitivity_and_specificity_is_loaded_and_returned(tmp_path, lineage):
    rec = load_status_source(_status_file(tmp_path, MEASURE))
    ident = ModuleIdentity("m", "1", "p", "a")
    entry, tested, reason = resolve_entry(rec, ident, 40, lineage)
    assert (tested, reason) == (40, "")
    assert entry.measure["sensitivity"]["value"] == 0.603 and entry.measure["n_pos"] == 232
    assert resolve_entry(rec, ident, 41, lineage)[0] is None  # taxon 41 is not covered


def test_an_entry_without_a_measure_is_allowed(tmp_path):
    assert load_status_source(_status_file(tmp_path, None)).entries[0].measure is None


@pytest.mark.parametrize(
    "bad",
    [
        {**MEASURE, "sensitivity": {"value": 1.2, "lo": 1.0, "hi": 1.3}},
        {**MEASURE, "specificity": {"value": 0.5, "lo": 0.6, "hi": 0.7}},
        {**MEASURE, "sensitivity": {"value": 0.5}},
        {**MEASURE, "n_pos": -1},
        {**MEASURE, "surprise": 1},
    ],
)
def test_a_bad_measure_is_refused(tmp_path, bad):
    with pytest.raises(ValueError):
        load_status_source(_status_file(tmp_path, bad))


def test_an_unknown_entry_key_is_refused(tmp_path):
    with pytest.raises(ValueError, match="unknown key"):
        load_status_source(_status_file(tmp_path, None, {"surprise": 1}))


def test_a_non_integer_parent_in_nodes_dmp_names_the_line(tmp_path):
    path = tmp_path / "nodes.dmp"
    path.write_text("1\t|\t1\t|\tno rank\t|\n2\t|\tx\t|\tno rank\t|\n")
    with pytest.raises(TaxonError, match="nodes.dmp:2"):
        Lineage.from_nodes_dmp(path)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_taxonomy_status.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.status`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/taxonomy.py` with exactly this content:

```python
"""NCBI taxonomy lineage lookups (from nodes.dmp) for status matching."""

from pathlib import Path


class TaxonError(ValueError):
    """A taxon ID is not in the taxonomy."""


class Lineage:
    def __init__(self, parent):
        self._parent = dict(parent)

    @classmethod
    def from_nodes_dmp(cls, path):
        """Read ``nodes.dmp`` (fields separated by ``\\t|\\t``; field 0 = taxid, field 1 = parent)."""
        parent = {}
        text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
        for n, line in enumerate(text.splitlines(), 1):
            fields = line.split("\t|\t")
            if len(fields) < 2:
                continue
            try:
                parent[int(fields[0])] = int(fields[1].replace("\t|", "").strip())
            except ValueError:
                raise TaxonError(f"{path}:{n}: taxon IDs are not integers") from None
        return cls(parent)

    def ancestors(self, taxon):
        """The taxon first, then its parents up to the root."""
        taxon = int(taxon)
        if taxon not in self._parent:
            raise TaxonError(f"taxon {taxon} is not in the taxonomy")
        chain = [taxon]
        while True:
            up = self._parent[chain[-1]]
            if up == chain[-1] or up in chain:  # root, or a cycle: stop
                return chain
            if up not in self._parent:
                raise TaxonError(f"parent {up} of taxon {chain[-1]} is not in the taxonomy")
            chain.append(up)

    def is_descendant_or_self(self, taxon, ancestor):
        return int(ancestor) in self.ancestors(taxon)

    def depth(self, taxon):
        return len(self.ancestors(taxon)) - 1
```

Create `src/cellsurface_sorting_hat/status.py` with exactly this content:

```python
"""Validation status of a module, as a function of (module, version, taxon).

A status applies to a protein's taxon only if that taxon is a tested taxon or a descendant of one.
A shared broad label (for example "Eurotiomycetes") is not enough.

Each entry of a status source can carry an optional ``measure`` object with the calibration set and
the sensitivity and specificity measured on it, so reports can show them.
"""

import json
from dataclasses import dataclass
from pathlib import Path

ESTIMATED, SMOKE, UNVALIDATED = "estimated", "smoke", "unvalidated"
_STRENGTH = {ESTIMATED: 2, SMOKE: 1, UNVALIDATED: 0}
ENTRY_KEYS = {"taxa", "status", "source", "measure"}
MEASURE_KEYS = {
    "calibration_set",
    "truth_source",
    "n_pos",
    "n_neg",
    "sensitivity",
    "specificity",
    "notes",
}


def weakest(statuses):
    """The weakest of the given statuses; ``unvalidated`` if there are none."""
    statuses = list(statuses)
    for s in statuses:
        if s not in _STRENGTH:
            raise ValueError(f"unknown status: {s!r}")
    if not statuses:
        return UNVALIDATED
    return min(statuses, key=_STRENGTH.__getitem__)


@dataclass(frozen=True)
class ModuleIdentity:
    name: str
    version: str
    params_hash: str
    artefact_hash: str


@dataclass(frozen=True)
class StatusEntry:
    taxa: tuple
    status: str
    source: str = ""
    measure: dict | None = None


@dataclass(frozen=True)
class StatusRecord:
    identity: ModuleIdentity
    entries: tuple  # of StatusEntry


def _check_rate(rate, where):
    if (
        not isinstance(rate, dict)
        or set(rate) != {"value", "lo", "hi"}
        or not all(isinstance(rate[k], int | float) for k in rate)
        or not 0 <= rate["lo"] <= rate["value"] <= rate["hi"] <= 1
    ):
        raise ValueError(f"{where}: needs value, lo, hi with 0 <= lo <= value <= hi <= 1")


def _check_measure(measure, where):
    if not isinstance(measure, dict):
        raise ValueError(f"{where}: measure must be an object")
    unknown = set(measure) - MEASURE_KEYS
    if unknown:
        raise ValueError(f"{where}: unknown measure key(s): {sorted(unknown)}")
    for key in ("sensitivity", "specificity"):
        if key in measure:
            _check_rate(measure[key], f"{where}.{key}")
    for key in ("n_pos", "n_neg"):
        if key in measure and not (isinstance(measure[key], int) and measure[key] >= 0):
            raise ValueError(f"{where}.{key}: must be a non-negative integer")


def load_status_source(path):
    data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    identity = ModuleIdentity(
        data["module"], data["version"], data["params_hash"], data["artefact_hash"]
    )
    entries = []
    for i, e in enumerate(data["entries"]):
        where = f"{path} entries[{i}]"
        unknown = set(e) - ENTRY_KEYS
        if unknown:
            raise ValueError(f"{where}: unknown key(s): {sorted(unknown)}")
        if e["status"] not in _STRENGTH:
            raise ValueError(f"unknown status {e['status']!r} in {path}")
        measure = e.get("measure")
        if measure is not None:
            _check_measure(measure, where)
        entries.append(
            StatusEntry(
                tuple(int(t) for t in e["taxa"]), e["status"], str(e.get("source", "")), measure
            )
        )
    return StatusRecord(identity, tuple(entries))


def resolve_entry(record, running, taxon, lineage):
    """Return (entry, tested_taxon, reason). ``entry`` is None when nothing applies.

    A stale record is refused: name, version, params hash and artefact hash must match.
    """
    if record is None:
        return None, None, "no status_source"
    for field in ("name", "version", "params_hash", "artefact_hash"):
        if getattr(record.identity, field) != getattr(running, field):
            return None, None, f"status_source stale: {field} differs"
    best = None  # (depth, tested taxon, entry)
    for entry in record.entries:
        for tested in entry.taxa:
            if lineage.is_descendant_or_self(taxon, tested):
                depth = lineage.depth(tested)
                if best is None or depth > best[0]:
                    best = (depth, tested, entry)
    if best is None:
        return None, None, "taxon not tested"
    return best[2], best[1], ""


def resolve_status(record, running, taxon, lineage):
    """Return (status, basis)."""
    entry, tested, reason = resolve_entry(record, running, taxon, lineage)
    if entry is None:
        return UNVALIDATED, reason
    return entry.status, f"taxon:{tested}"
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_taxonomy_status.py -q`
Expected: `24 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-core
git add src/cellsurface_sorting_hat/taxonomy.py src/cellsurface_sorting_hat/status.py tests/cellsurface_sorting_hat/conftest.py tests/cellsurface_sorting_hat/test_taxonomy_status.py
git commit -m "feat(sorting-hat): taxonomy lineage and module status resolution

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 3: FASTA reading and checks

**Files:**
- Create: `src/cellsurface_sorting_hat/fasta.py`
- Test: `tests/cellsurface_sorting_hat/test_fasta.py`

**Interfaces:**
- Produces: `read_fasta(path) -> list[Protein]` (plain, `.gz`, `.zst`); `Protein(id, sequence, sha256, state, note, ambiguous_fraction)` with `state` `ok` or `na_invalid`; constants `OK`, `NA_INVALID`; `FastaError` for duplicate IDs, an empty file and a file with no valid protein.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_fasta.py` with exactly this content:

```python
import gzip
import hashlib
import shutil
import subprocess

import pytest

from cellsurface_sorting_hat.fasta import NA_INVALID, OK, FastaError, read_fasta


def _write(tmp_path, text, name="p.faa"):
    path = tmp_path / name
    path.write_bytes(text.encode())
    return path


def test_reads_proteins_and_hashes_the_sequence(tmp_path):
    p = read_fasta(_write(tmp_path, ">A first protein\nMKT\nAYI\n>B\nMKTAYI\n"))
    assert [x.id for x in p] == ["A", "B"]
    assert p[0].sequence == "MKTAYI"
    assert p[0].sha256 == hashlib.sha256(b"MKTAYI").hexdigest()
    assert p[0].sha256 == p[1].sha256  # identical sequences, different IDs, are allowed


def test_trailing_stop_is_stripped_silently(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKTAYI*\n>B\nMKT\n"))
    assert (p[0].sequence, p[0].state) == ("MKTAYI", OK)


def test_internal_stop_is_invalid(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKT*AYI\n>B\nMKT\n"))
    assert p[0].state == NA_INVALID and "stop" in p[0].note


def test_non_residue_characters_are_invalid(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKT-AY1\n>B\nMKT\n"))
    assert p[0].state == NA_INVALID and "-" in p[0].note


def test_ambiguous_residues_are_allowed_and_counted(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKXXAYIB\n>B\nMKT\n"))
    assert p[0].state == OK
    assert p[0].ambiguous_fraction == pytest.approx(3 / 8)


def test_empty_sequence_is_invalid(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\n>B\nMKT\n"))
    assert p[0].state == NA_INVALID


def test_crlf_and_lowercase(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\r\nmktayi\r\n>B\r\nMKT\r\n"))
    assert p[0].sequence == "MKTAYI"


def test_duplicate_ids_stop_the_run(tmp_path):
    with pytest.raises(FastaError, match="duplicate ID: A"):
        read_fasta(_write(tmp_path, ">A\nMKT\n>A\nMKS\n"))


def test_empty_file_and_no_valid_protein_stop_the_run(tmp_path):
    with pytest.raises(FastaError, match="no records"):
        read_fasta(_write(tmp_path, ""))
    with pytest.raises(FastaError, match="no valid proteins"):
        read_fasta(_write(tmp_path, ">A\nMK*T\n"))


def test_gzip(tmp_path):
    path = tmp_path / "p.faa.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(">A\nMKT\n")
    assert read_fasta(path)[0].sequence == "MKT"


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd is not installed")
def test_zstd(tmp_path):
    plain = _write(tmp_path, ">A\nMKT\n")
    subprocess.run(["zstd", "-q", str(plain), "-o", str(tmp_path / "p.faa.zst")], check=True)
    assert read_fasta(tmp_path / "p.faa.zst")[0].sequence == "MKT"


def test_trailing_stop_is_flagged(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKT*\n>B\nMKT\n"))
    assert (p[0].trailing_stop, p[1].trailing_stop) == (True, False)


def test_sequence_lines_before_the_first_header_stop_the_run(tmp_path):
    with pytest.raises(FastaError, match="before the first header"):
        read_fasta(_write(tmp_path, "MKT\n>A\nMKT\n"))


def test_a_byte_order_mark_is_ignored(tmp_path):
    path = tmp_path / "p.faa"
    path.write_bytes(b"\xef\xbb\xbf>A\nMKT\n")
    assert read_fasta(path)[0].id == "A"


def test_bytes_that_are_not_utf8_stop_the_run_with_the_file_name(tmp_path):
    path = tmp_path / "p.faa"
    path.write_bytes(b">A \xe9\nMKT\n")
    with pytest.raises(FastaError, match="p.faa: cannot be read as text"):
        read_fasta(path)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_fasta.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.fasta`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/fasta.py` with exactly this content:

```python
"""Read and check a protein FASTA (plain, .gz or .zst).

Rules (spec section 3.8): a trailing ``*`` is stripped; an internal ``*`` or a character that is not
a residue makes the protein ``na_invalid``; ``X B Z U J O`` are allowed and counted; duplicate IDs,
an empty file and a file with no valid protein stop the run.
"""

import gzip
import hashlib
import subprocess
import zlib
from dataclasses import dataclass
from pathlib import Path

OK = "ok"
NA_INVALID = "na_invalid"
STANDARD = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS = frozenset("XBZUJO")


class FastaError(ValueError):
    """The FASTA cannot be used; the run must stop."""


@dataclass(frozen=True)
class Protein:
    id: str
    sequence: str
    sha256: str
    state: str
    note: str
    ambiguous_fraction: float
    trailing_stop: bool = False


def _lines(path):
    path = Path(path)
    try:
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="utf-8-sig") as fh:
                yield from fh
        elif path.suffix == ".zst":
            out = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, check=True)
            yield from out.stdout.decode("utf-8-sig").splitlines()
        else:
            with open(path, encoding="utf-8-sig") as fh:
                yield from fh
    except (UnicodeDecodeError, EOFError, zlib.error) as err:
        raise FastaError(f"{path}: cannot be read as text ({err.__class__.__name__})") from err


def _records(path):
    header, chunks = None, []
    for raw in _lines(path):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:].strip(), []
        elif header is not None:
            chunks.append(line)
        else:
            raise FastaError("sequence lines before the first header")
    if header is not None:
        yield header, "".join(chunks)


def _check(sequence):
    if not sequence:
        return NA_INVALID, "empty sequence", 0.0
    if "*" in sequence:
        return NA_INVALID, "internal stop codon", 0.0
    bad = sorted(set(sequence) - STANDARD - AMBIGUOUS)
    if bad:
        return NA_INVALID, "not residues: " + "".join(bad), 0.0
    ambiguous = sum(1 for c in sequence if c in AMBIGUOUS) / len(sequence)
    return OK, "", ambiguous


def read_fasta(path):
    proteins, seen = [], set()
    for header, raw_seq in _records(path):
        if not header:
            raise FastaError("a record has an empty header")
        pid = header.split()[0]
        if pid in seen:
            raise FastaError(f"duplicate ID: {pid}")
        seen.add(pid)
        sequence = raw_seq.upper()
        trailing = sequence.endswith("*")
        if trailing:
            sequence = sequence[:-1]
        state, note, ambiguous = _check(sequence)
        sha = hashlib.sha256(sequence.encode()).hexdigest()
        proteins.append(Protein(pid, sequence, sha, state, note, ambiguous, trailing))
    if not proteins:
        raise FastaError("the FASTA has no records")
    if all(p.state != OK for p in proteins):
        raise FastaError("the FASTA has no valid proteins")
    return proteins
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_fasta.py -q`
Expected: `15 passed` (the zstd test is skipped when `zstd` is not installed, then `14 passed, 1 skipped`)

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-core
git add src/cellsurface_sorting_hat/fasta.py tests/cellsurface_sorting_hat/test_fasta.py
git commit -m "feat(sorting-hat): FASTA reader and checks

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 4: Atomic writes and the module cache

**Files:**
- Create: `src/cellsurface_sorting_hat/cache.py`
- Test: `tests/cellsurface_sorting_hat/test_cache.py`

**Interfaces:**
- Consumes: `ModuleIdentity` from Task 2.
- Produces: `identity_key(identity, tool_versions=None) -> str`; `write_atomic(path, data: bytes)`; `read_verified(path) -> bytes` (raises `CacheError` for a missing or wrong `.sha256` sidecar); `ModuleCache(directory, key)` with `.load() -> dict[sha256, row]` and `.update(rows)`. `update` holds an exclusive `fcntl` lock and `load` a shared one; a table that fails its checksum loads as `{}` (a cache miss).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_cache.py` with exactly this content:

```python
import os
from concurrent.futures import ProcessPoolExecutor

import pytest

from cellsurface_sorting_hat.cache import (
    CacheError,
    ModuleCache,
    identity_key,
    read_verified,
    write_atomic,
)
from cellsurface_sorting_hat.status import ModuleIdentity

IDENT = ModuleIdentity("pfam", "1", "p", "a")


def test_identity_key_changes_with_every_field():
    base = identity_key(IDENT, {"hmmer": "3.4"})
    assert base == identity_key(IDENT, {"hmmer": "3.4"})
    for changed in (
        ModuleIdentity("pfam", "2", "p", "a"),
        ModuleIdentity("pfam", "1", "q", "a"),
        ModuleIdentity("pfam", "1", "p", "b"),
        ModuleIdentity("other", "1", "p", "a"),
    ):
        assert identity_key(changed, {"hmmer": "3.4"}) != base
    assert identity_key(IDENT, {"hmmer": "3.5"}) != base


def test_write_atomic_leaves_data_and_checksum_and_no_temporary_file(tmp_path):
    path = tmp_path / "d" / "f.bin"
    write_atomic(path, b"abc")
    assert read_verified(path) == b"abc"
    assert sorted(p.name for p in path.parent.iterdir()) == ["f.bin", "f.bin.sha256"]


def test_read_verified_refuses_a_missing_or_wrong_checksum(tmp_path):
    path = tmp_path / "f.bin"
    path.write_bytes(b"abc")
    with pytest.raises(CacheError, match="no checksum"):
        read_verified(path)
    write_atomic(path, b"abc")
    path.write_bytes(b"abd")
    with pytest.raises(CacheError, match="does not match"):
        read_verified(path)


def test_module_cache_merges_rows_by_sha256(tmp_path):
    cache = ModuleCache(tmp_path, "k1")
    assert cache.load() == {}
    cache.update([{"sha256": "aa", "call": "called"}, {"sha256": "bb", "call": "not_called"}])
    cache.update([{"sha256": "bb", "call": "called"}, {"sha256": "cc", "call": "called", "x": "1"}])
    table = ModuleCache(tmp_path, "k1").load()
    assert set(table) == {"aa", "bb", "cc"}
    assert table["bb"]["call"] == "called"
    assert table["cc"]["x"] == "1"
    assert ModuleCache(tmp_path, "k2").load() == {}  # another identity key: another table


def test_a_damaged_table_is_a_cache_miss_and_the_next_update_rebuilds_it(tmp_path):
    cache = ModuleCache(tmp_path, "k1")
    cache.update([{"sha256": "aa", "call": "called"}])
    cache.path.write_bytes(b"not a gzip file")  # the sidecar no longer matches
    assert cache.load() == {}
    cache.update([{"sha256": "bb", "call": "called"}])
    assert set(cache.load()) == {"bb"}


def _worker(args):
    directory, worker = args
    cache = ModuleCache(directory, "shared")
    for i in range(15):
        cache.update([{"sha256": f"w{worker}-{i}", "call": "called"}])
        cache.load()  # readers run while other processes write
    return worker


def test_parallel_updates_lose_no_rows_and_leave_a_readable_table(tmp_path):
    with ProcessPoolExecutor(max_workers=4) as pool:
        list(pool.map(_worker, [(str(tmp_path), w) for w in range(4)]))
    table = ModuleCache(tmp_path, "shared").load()
    assert len(table) == 60
    assert read_verified(ModuleCache(tmp_path, "shared").path)  # data and checksum agree


def test_a_reader_can_read_when_the_lock_file_cannot_be_opened(tmp_path):
    cache = ModuleCache(tmp_path, "k1")
    cache.update([{"sha256": "aa", "call": "called"}])
    lock = cache._lock
    lock.chmod(0o444)
    try:
        if os.access(lock, os.W_OK):
            pytest.skip("this user can write to a 0444 file (for example root)")
        assert set(cache.load()) == {"aa"}
    finally:
        lock.chmod(0o644)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_cache.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.cache`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/cache.py` with exactly this content:

```python
"""Module result cache and atomic writes.

Results are stored per sequence sha256 in one table per module identity. The identity key covers
module name, version, params hash, artefact hash and tool versions, so a change in any of them
starts a new table.

Concurrency: ``ModuleCache`` takes an exclusive ``fcntl`` lock on ``<table>.lock`` for an update and
a shared lock for a read, so a reader never sees a data file with a sidecar that belongs to another
version. A table that fails its checksum (for example after a crash between the two renames) is
treated as a cache miss and is rebuilt by the next update. ``fcntl`` locks may not be reliable on
every network file system; use a node-local or otherwise tested directory for parallel runs.
"""

import csv
import fcntl
import gzip
import hashlib
import io
import json
import os
import uuid
from contextlib import contextmanager
from pathlib import Path


class CacheError(RuntimeError):
    """A cached file is missing its checksum or does not match it."""


def identity_key(identity, tool_versions=None):
    payload = {
        "name": identity.name,
        "version": identity.version,
        "params_hash": identity.params_hash,
        "artefact_hash": identity.artefact_hash,
        "tools": dict(sorted((tool_versions or {}).items())),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _tmp_name(path):
    return path.with_name(f"{path.name}.tmp{os.getpid()}-{uuid.uuid4().hex[:8]}")


def write_atomic(path, data):
    """Write bytes to ``path`` through a temporary file, then write ``path + '.sha256'``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp_name(path)
    tmp.write_bytes(data)
    os.replace(tmp, path)
    side = path.with_name(path.name + ".sha256")
    side_tmp = _tmp_name(side)
    side_tmp.write_text(hashlib.sha256(data).hexdigest() + "\n")
    os.replace(side_tmp, side)


def read_verified(path):
    path = Path(path)
    side = path.with_name(path.name + ".sha256")
    if not side.exists():
        raise CacheError(f"{path} has no checksum file")
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != side.read_text().strip():
        raise CacheError(f"{path} does not match its checksum")
    return data


@contextmanager
def _locked(lock_path, exclusive):
    """Hold a lock on ``lock_path``. A reader that cannot open the lock file (a read-only
    directory) reads without a lock. A steady stream of readers can delay a writer."""
    if exclusive:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fh = open(lock_path, "a")
    except OSError:
        if exclusive:
            raise
        yield
        return
    with fh:
        fcntl.flock(fh, fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


class ModuleCache:
    """A table of module rows keyed by sequence sha256, for one module identity."""

    def __init__(self, directory, key):
        self.path = Path(directory) / f"{key}.tsv.gz"
        self._lock = self.path.with_name(self.path.name + ".lock")

    def _load_unlocked(self):
        if not self.path.exists():
            return {}
        try:
            text = gzip.decompress(read_verified(self.path)).decode()
        except (CacheError, OSError, EOFError):
            return {}  # a damaged table is a cache miss
        return {row["sha256"]: row for row in csv.DictReader(io.StringIO(text), delimiter="\t")}

    def load(self):
        with _locked(self._lock, exclusive=False):
            return self._load_unlocked()

    def update(self, rows):
        """Merge ``rows`` (dicts with a ``sha256`` key) into the table and rewrite it atomically."""
        with _locked(self._lock, exclusive=True):
            table = self._load_unlocked()
            for row in rows:
                table[row["sha256"]] = row
            columns = ["sha256"] + sorted({k for r in table.values() for k in r} - {"sha256"})
            buf = io.StringIO()
            writer = csv.DictWriter(buf, columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for sha in sorted(table):
                writer.writerow(table[sha])
            write_atomic(self.path, gzip.compress(buf.getvalue().encode(), mtime=0))
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_cache.py -q`
Expected: `7 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-core
git add src/cellsurface_sorting_hat/cache.py tests/cellsurface_sorting_hat/test_cache.py
git commit -m "feat(sorting-hat): atomic writes and per-sequence module cache

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 5: Rules file and the category engine

**Files:**
- Create: `src/cellsurface_sorting_hat/categories.yaml`, `src/cellsurface_sorting_hat/engine.py`
- Test: `tests/cellsurface_sorting_hat/test_engine.py`

**Interfaces:**
- Consumes: `k_and`, `k_or`, `k_not`, the call constants (Task 1); `weakest`, `UNVALIDATED` (Task 2).
- Produces: `load_config(path=None) -> Config` (packaged file by default; raises `ConfigError`); `ModuleTable(name, rows)` with `.get(protein_id)`; `evaluate(cfg, protein_ids, taxa, modules, status_of) -> list[CallRecord]`; `CallRecord(protein, call, variant, value, status, status_basis, other_basis)`; `available_variants(cfg, modules)`; `variant_label(module_name)`.
- `status_of(module, taxon) -> (status, basis)` is supplied by the caller.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_engine.py` with exactly this content:

```python
"""Rules, gating, status and the `other_*` formulas."""

from pathlib import Path

import pytest
import yaml

from cellsurface_sorting_hat import engine
from cellsurface_sorting_hat.engine import (
    ConfigError,
    ModuleTable,
    collect_evidence,
    evaluate,
    load_config,
    referenced_modules,
)

R0 = "step1_rule@R0"
R2 = "step1_rule@R2"


def table(name, rows):
    return ModuleTable(name, {k: dict(v) for k, v in rows.items()})


def run(modules, status_of=None, ids=("P",), taxon=40):
    cfg = load_config()
    status_of = status_of or (lambda module, t: ("unvalidated", "x"))
    records = evaluate(cfg, list(ids), dict.fromkeys(ids, taxon), modules, status_of)
    return {(r.protein, r.call, r.variant): r for r in records}


def ok(**fields):
    return {"state": "ok", **fields}


def base_modules(step1="not_called", **over):
    mods = {
        R0: table(R0, {"P": ok(call=step1)}),
        "repeat02": table("repeat02", {"P": ok(call="not_called")}),
        "repeat14": table("repeat14", {"P": ok(call="not_called")}),
        "pfam_adhesion": table("pfam_adhesion", {"P": ok(hit="0")}),
        "pfam_allergen": table("pfam_allergen", {"P": ok(hit="0")}),
        "antigen_lookup": table("antigen_lookup", {"P": ok(percentile="50")}),
        "allergen_homology": table("allergen_homology", {"P": ok(identity="0", coverage="0")}),
    }
    mods.update(over)
    return mods


def test_packaged_config_loads_and_has_the_spec_calls():
    names = [c["name"] for c in load_config().calls]
    assert names == [
        "surface_glycoprotein",
        "adhesion_repeat",
        "adhesion_domain",
        "cell_wall_adhesion_ungated",
        "cell_wall_adhesion_candidate",
        "antigen_candidate",
        "antigen_candidate_surface",
        "allergen_homolog_hit",
        "allergen_candidate",
        "other_not_surface",
        "other_surface_no_mechanism",
    ]


def test_only_available_step1_variants_get_per_variant_records():
    res = run(base_modules())
    variants = {v for (_, call, v) in res if call == "surface_glycoprotein"}
    assert variants == {"R0"}


def test_a_true_or_input_wins_over_an_unknown_input():
    mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})
    r = run(mods)[("P", "adhesion_repeat", "")]
    assert r.value == "called"


def test_unknown_or_input_with_a_false_other_is_unknown():
    mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
    assert run(mods)[("P", "adhesion_repeat", "")].value == "not_assessable"


def test_false_surface_makes_gated_antigen_false_even_when_antigen_is_unknown():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
    res = run(mods)
    assert res[("P", "antigen_candidate", "")].value == "not_assessable"
    assert res[("P", "antigen_candidate_surface", "R0")].value == "not_called"


def test_antigen_threshold_is_inclusive():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="15")}))
    assert run(mods)[("P", "antigen_candidate", "")].value == "called"
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="15.01")}))
    assert run(mods)[("P", "antigen_candidate", "")].value == "not_called"


def test_empty_or_nan_number_is_unknown():
    for bad in ("", "nan", "abc"):
        mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile=bad)}))
        assert run(mods)[("P", "antigen_candidate", "")].value == "not_assessable", bad


@pytest.mark.parametrize(
    "identity,coverage,pfam,expected",
    [
        ("82", "90", "0", "called"),
        ("69.9", "90", "0", "not_called"),
        ("82", "79", "0", "not_called"),
        ("40", "85", "1", "called"),
        ("0", "0", "0", "not_called"),
    ],
)
def test_allergen_two_tier_call(identity, coverage, pfam, expected):
    mods = base_modules(
        allergen_homology=table(
            "allergen_homology", {"P": ok(identity=identity, coverage=coverage)}
        ),
        pfam_allergen=table("pfam_allergen", {"P": ok(hit=pfam)}),
    )
    assert run(mods)[("P", "allergen_candidate", "")].value == expected


def test_allergen_needs_no_surface_call():
    mods = base_modules(step1="not_called")
    mods["allergen_homology"] = table("allergen_homology", {"P": ok(identity="90", coverage="95")})
    res = run(mods)
    assert res[("P", "surface_glycoprotein", "R0")].value == "not_called"
    assert res[("P", "allergen_candidate", "")].value == "called"


def test_other_not_surface_leaves_out_unknown_calls_and_names_them():
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": {"state": "not_applicable"}}))
    r = run(mods)[("P", "other_not_surface", "R0")]
    assert (r.value, r.other_basis) == ("called", "antigen_candidate")


def test_other_not_surface_is_false_when_a_mechanism_call_is_true():
    mods = base_modules(
        allergen_homology=table("allergen_homology", {"P": ok(identity="90", coverage="95")})
    )
    assert run(mods)[("P", "other_not_surface", "R0")].value == "not_called"


def test_other_surface_no_mechanism():
    res = run(base_modules(step1="called"))
    assert res[("P", "other_surface_no_mechanism", "R0")].value == "called"
    assert res[("P", "other_not_surface", "R0")].value == "not_called"


def test_other_is_unknown_when_the_surface_call_is_unknown():
    mods = base_modules()
    mods[R0] = table(R0, {"P": {"state": "error"}})
    res = run(mods)
    assert res[("P", "other_not_surface", "R0")].value == "not_assessable"
    assert res[("P", "other_surface_no_mechanism", "R0")].value == "not_assessable"


def test_other_can_be_true_in_one_variant_and_unknown_in_another():
    mods = base_modules()
    mods[R2] = table(R2, {"P": {"state": "error"}})
    res = run(mods)
    assert res[("P", "other_not_surface", "R0")].value == "called"
    assert res[("P", "other_not_surface", "R2")].value == "not_assessable"


def test_status_is_the_weakest_of_the_modules_that_decided_the_result():
    statuses = {
        R0: ("estimated", "taxon:40"),
        "repeat02": ("smoke", "taxon:30"),
        "repeat14": ("unvalidated", "no status_source"),
    }
    mods = base_modules(step1="called")
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})  # repeat14 stays not_called
    res = run(mods, status_of=lambda m, t: statuses.get(m, ("unvalidated", "")))
    # adhesion_repeat is true because of repeat02 only: status smoke, not unvalidated
    assert res[("P", "adhesion_repeat", "")].status == "smoke"
    # the gated call is true because of repeat02 and step 1: weakest is smoke
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert (r.value, r.status) == ("called", "smoke")
    assert "step1_rule@R0:taxon:40" in r.status_basis


def test_invalid_protein_rows_make_every_call_unknown():
    mods = base_modules()
    for name in mods:
        mods[name] = table(name, {"P": {"state": "na_invalid"}})
    res = run(mods)
    assert {r.value for r in res.values()} == {"not_assessable"}


def _config_with(tmp_path, mutate):
    packaged = Path(engine.__file__).with_name("categories.yaml")
    data = yaml.safe_load(packaged.read_text())
    mutate(data)
    path = tmp_path / "c.yaml"
    path.write_text(yaml.safe_dump(data))
    return path


@pytest.mark.parametrize(
    "mutate,message",
    [
        (lambda d: d.update(default_gate="nope"), "default_gate"),
        (lambda d: d["calls"].append(dict(d["calls"][0])), "duplicate call"),
        (lambda d: d["calls"][1].update(expr={"ref": "later_call"}), "unknown or later"),
        (lambda d: d["calls"][1].update(expr={"bogus": 1}), "bad node"),
        (lambda d: d["calls"][5]["expr"]["test"].update(op="~"), "bad operator"),
        (lambda d: d["calls"][5]["expr"]["test"].update(value="$nope"), "unknown threshold"),
        (
            lambda d: d["calls"][1].update(expr={"ref": "surface_glycoprotein"}),
            "ungated call cannot use",
        ),
    ],
)
def test_bad_config_is_refused(tmp_path, mutate, message):
    with pytest.raises(ConfigError, match=message):
        load_config(_config_with(tmp_path, mutate))


def test_an_unknown_result_has_status_unvalidated_and_no_basis():
    # R0 is in state error; repeat02 is called with status smoke. The gated call is unknown.
    mods = base_modules(step1="called")
    mods[R0] = table(R0, {"P": {"state": "error"}})
    mods["repeat02"] = table("repeat02", {"P": ok(call="called")})
    res = run(mods, status_of=lambda m, t: ("smoke", "t"))
    r = res[("P", "cell_wall_adhesion_candidate", "R0")]
    assert (r.value, r.status, r.status_basis) == ("not_assessable", "unvalidated", "")
    # an OR with a false input and an unknown input is unknown and also has no deciding module
    mods = base_modules(repeat14=table("repeat14", {"P": {"state": "error"}}))
    r = run(mods, status_of=lambda m, t: ("estimated", "t"))[("P", "adhesion_repeat", "")]
    assert (r.value, r.status) == ("not_assessable", "unvalidated")


@pytest.mark.parametrize(
    "identity,length,expected",
    [
        ("35", "80", "called"),
        ("34.9", "120", "not_called"),
        ("60", "79", "not_called"),
        ("", "90", "not_assessable"),
    ],
)
def test_allergen_homolog_hit_is_the_35_percent_80_aa_evidence_tier(identity, length, expected):
    mods = base_modules(
        allergen_homology=table(
            "allergen_homology", {"P": ok(identity=identity, aligned_length=length, coverage="0")}
        )
    )
    assert run(mods)[("P", "allergen_homolog_hit", "")].value == expected


def test_evidence_rows_are_collected_only_for_ok_rows_with_a_value():
    cfg = load_config()
    mods = {
        "antigen_lookup": table(
            "antigen_lookup",
            {
                "A": ok(percentile="7.1", antigenicity=""),
                "B": {"state": "not_applicable", "percentile": "9"},
            },
        )
    }
    rows = collect_evidence(cfg, ["A", "B"], mods)
    assert rows == [("A", "antigen_lookup", "percentile", "7.1")]


def test_referenced_modules_expand_the_step1_variants():
    names = referenced_modules(load_config())
    assert "step1_rule@R0" in names and "step1_ml@card" in names
    assert {
        "repeat02",
        "repeat14",
        "pfam_adhesion",
        "pfam_allergen",
        "antigen_lookup",
        "allergen_homology",
    } <= set(names)
    assert not any("{step1}" in n for n in names)


@pytest.mark.parametrize(
    "mutate,message",
    [
        (
            lambda d: d["calls"][-2].update(surface="adhesion_repeat"),
            "needs an earlier per_variant surface",
        ),
        (
            lambda d: d["calls"][1].update(expr={"call": "{step1}"}),
            "only allowed in a per_variant call",
        ),
        (lambda d: d.update(evidence=["no_dot"]), "evidence entry"),
    ],
)
def test_more_bad_config_is_refused(tmp_path, mutate, message):
    with pytest.raises(ConfigError, match=message):
        load_config(_config_with(tmp_path, mutate))


@pytest.mark.parametrize("bad", ["inf", "-inf", "Infinity"])
def test_a_non_finite_number_is_unknown(bad):
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile=bad)}))
    assert run(mods)[("P", "antigen_candidate", "")].value == "not_assessable"


def test_a_false_and_takes_its_status_from_the_false_inputs_only():
    # surface (R0) is false with status smoke; antigen is true with status estimated
    statuses = {R0: ("smoke", "t"), "antigen_lookup": ("estimated", "t")}
    mods = base_modules(antigen_lookup=table("antigen_lookup", {"P": ok(percentile="5")}))
    res = run(mods, status_of=lambda m, t: statuses.get(m, ("unvalidated", "")))
    r = res[("P", "antigen_candidate_surface", "R0")]
    assert (r.value, r.status) == ("not_called", "smoke")
    assert r.status_basis == "step1_rule@R0:t"
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_engine.py -q`
Expected: collection error, `ImportError: cannot import name 'engine'`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/categories.yaml` with exactly this content:

```yaml
# cellsurface_sorting_hat: category rules (version 1). See
# docs/superpowers/specs/2026-10-04-orchestrator-design.md section 3.4.
#
# Node kinds:
#   {call: MODULE}                  the module's `call` column (called / not_called)
#   {flag: MODULE.FIELD}            a 1/0 field
#   {test: {module, field, op, value}}   a number compared with value; "$name" reads `thresholds`
#   {ref: CALL}                     another call (same step 1 variant if that call is per_variant)
#   {and: [...]}, {or: [...]}, {not: NODE}
# `{step1}` in a module name is replaced by the step 1 variant module for per_variant calls.
# A field that is empty, or a module whose state is not `ok`, gives not_assessable.

version: 1
default_gate: "step1_rule@R0"
step1_variants: ["step1_rule@R0", "step1_rule@R1", "step1_rule@R2", "step1_ml@card"]

thresholds:
  antigen_percentile_max: 15
  allergen_hit_identity_min: 35
  allergen_hit_length_min: 80
  allergen_identity_min: 70
  allergen_coverage_min: 80

# Fields copied to evidence.tsv.gz when the module row is `ok` and the field is not empty.
evidence:
  - antigen_lookup.percentile
  - antigen_lookup.antigenicity
  - antigen_lookup.specificity
  - antigen_lookup.prevalence
  - antigen_lookup.max_crossreact
  - allergen_homology.allergen_name
  - allergen_homology.identity
  - allergen_homology.coverage
  - allergen_homology.aligned_length
  - cys_rich.tier
  - expression.log2fc

calls:
  - name: surface_glycoprotein
    per_variant: true
    expr: {call: "{step1}"}
  - name: adhesion_repeat
    expr:
      or: [{call: repeat02}, {call: repeat14}]
  - name: adhesion_domain
    expr: {flag: pfam_adhesion.hit}
  - name: cell_wall_adhesion_ungated
    expr:
      or: [{ref: adhesion_repeat}, {ref: adhesion_domain}]
  - name: cell_wall_adhesion_candidate
    per_variant: true
    expr:
      and:
        - or: [{ref: adhesion_repeat}, {ref: adhesion_domain}]
        - {ref: surface_glycoprotein}
  - name: antigen_candidate
    expr:
      test: {module: antigen_lookup, field: percentile, op: "<=", value: "$antigen_percentile_max"}
  - name: antigen_candidate_surface
    per_variant: true
    expr:
      and: [{ref: antigen_candidate}, {ref: surface_glycoprotein}]
  - name: allergen_homolog_hit
    expr:
      and:
        - test: {module: allergen_homology, field: identity, op: ">=", value: "$allergen_hit_identity_min"}
        - test: {module: allergen_homology, field: aligned_length, op: ">=", value: "$allergen_hit_length_min"}
  - name: allergen_candidate
    expr:
      or:
        - and:
            - test: {module: allergen_homology, field: identity, op: ">=", value: "$allergen_identity_min"}
            - test: {module: allergen_homology, field: coverage, op: ">=", value: "$allergen_coverage_min"}
        - {flag: pfam_allergen.hit}
  - name: other_not_surface
    per_variant: true
    kind: other
    surface: surface_glycoprotein
    surface_is: not_called
    mechanism: [adhesion_repeat, adhesion_domain, antigen_candidate, allergen_candidate]
  - name: other_surface_no_mechanism
    per_variant: true
    kind: other
    surface: surface_glycoprotein
    surface_is: called
    mechanism: [adhesion_repeat, adhesion_domain, antigen_candidate, allergen_candidate]
```

Create `src/cellsurface_sorting_hat/engine.py` with exactly this content:

```python
"""Category engine: evaluate the rules in categories.yaml for each protein.

Inputs are module tables (``ModuleTable``) and a function ``status_of(module, taxon)``. The engine
does no I/O except reading the config.
"""

import hashlib
import math
import operator
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from cellsurface_sorting_hat.logic import (
    CALLED,
    NOT_ASSESSABLE,
    NOT_CALLED,
    k_and,
    k_not,
    k_or,
)
from cellsurface_sorting_hat.status import UNVALIDATED, weakest

OK_STATE = "ok"
_OPS = {">=": operator.ge, "<=": operator.le, ">": operator.gt, "<": operator.lt, "==": operator.eq}
_NODE_KEYS = {"call", "flag", "test", "ref", "and", "or", "not"}


class ConfigError(ValueError):
    """categories.yaml is not valid."""


@dataclass
class ModuleTable:
    """Per-protein rows of one module: id -> dict of strings (``state``, ``call``, other fields)."""

    name: str
    rows: dict = field(default_factory=dict)

    def get(self, protein_id):
        return self.rows.get(protein_id)


@dataclass
class Config:
    default_gate: str
    step1_variants: list
    thresholds: dict
    calls: list
    evidence: list
    sha256: str

    def call_by_name(self, name):
        for c in self.calls:
            if c["name"] == name:
                return c
        raise ConfigError(f"unknown call: {name}")


@dataclass(frozen=True)
class CallRecord:
    protein: str
    call: str
    variant: str
    value: str
    status: str
    status_basis: str
    other_basis: str


def variant_label(step1_module):
    """``step1_rule@R0`` -> ``R0``."""
    return step1_module.split("@", 1)[1] if "@" in step1_module else step1_module


def load_config(path=None):
    path = Path(path) if path else Path(__file__).with_name("categories.yaml")
    raw = path.read_bytes()
    data = yaml.safe_load(raw)
    cfg = Config(
        default_gate=data["default_gate"],
        step1_variants=list(data["step1_variants"]),
        thresholds=dict(data["thresholds"]),
        calls=list(data["calls"]),
        evidence=list(data.get("evidence", [])),
        sha256=hashlib.sha256(raw).hexdigest(),
    )
    _validate(cfg)
    return cfg


def _validate(cfg):
    if cfg.default_gate not in cfg.step1_variants:
        raise ConfigError("default_gate must be one of step1_variants")
    for item in cfg.evidence:
        if "." not in item:
            raise ConfigError(f"evidence entry needs MODULE.FIELD: {item!r}")
    seen = {}
    for call in cfg.calls:
        name = call["name"]
        if name in seen:
            raise ConfigError(f"duplicate call name: {name}")
        per_variant = bool(call.get("per_variant"))
        if call.get("kind") == "other":
            if not per_variant or not seen.get(call.get("surface")):
                raise ConfigError(
                    f"{name}: an `other` call needs an earlier per_variant surface call"
                )
            if call.get("surface_is") not in (CALLED, NOT_CALLED):
                raise ConfigError(f"{name}: surface_is must be called or not_called")
            for m in call["mechanism"]:
                if m not in seen or seen[m]:
                    raise ConfigError(f"{name}: mechanism {m} must be an earlier ungated call")
        else:
            _validate_node(call["expr"], seen, per_variant, cfg.thresholds, name)
        seen[name] = per_variant


def _check_module_name(module, per_variant, where):
    if "{step1}" in module and not per_variant:
        raise ConfigError(f"{where}: {{step1}} is only allowed in a per_variant call")


def _validate_node(node, seen, per_variant, thresholds, where):
    if not isinstance(node, dict) or len(node) != 1 or next(iter(node)) not in _NODE_KEYS:
        raise ConfigError(f"{where}: bad node {node!r}")
    kind, arg = next(iter(node.items()))
    if kind in ("and", "or"):
        if not isinstance(arg, list) or not arg:
            raise ConfigError(f"{where}: {kind} needs a non-empty list")
        for child in arg:
            _validate_node(child, seen, per_variant, thresholds, where)
    elif kind == "not":
        _validate_node(arg, seen, per_variant, thresholds, where)
    elif kind == "ref":
        if arg not in seen:
            raise ConfigError(f"{where}: ref to unknown or later call {arg}")
        if seen[arg] and not per_variant:
            raise ConfigError(f"{where}: an ungated call cannot use per_variant call {arg}")
    elif kind == "test":
        _check_module_name(arg["module"], per_variant, where)
        if arg["op"] not in _OPS:
            raise ConfigError(f"{where}: bad operator {arg['op']!r}")
        value = arg["value"]
        if isinstance(value, str) and (not value.startswith("$") or value[1:] not in thresholds):
            raise ConfigError(f"{where}: unknown threshold {value!r}")
    elif kind == "flag":
        if "." not in arg:
            raise ConfigError(f"{where}: flag needs MODULE.FIELD")
        _check_module_name(arg.split(".", 1)[0], per_variant, where)
    elif kind == "call":
        _check_module_name(arg, per_variant, where)


@dataclass
class _Result:
    value: str
    contributors: frozenset


def _leaf(value, module):
    return _Result(value, frozenset([module]) if value != NOT_ASSESSABLE else frozenset())


def _row(modules, module, protein_id):
    table = modules.get(module)
    if table is None:
        return None
    row = table.get(protein_id)
    if row is None or row.get("state", OK_STATE) != OK_STATE:
        return None
    return row


def _eval(node, ctx):
    kind, arg = next(iter(node.items()))
    if kind == "call":
        module = arg.replace("{step1}", ctx["step1"] or "")
        row = _row(ctx["modules"], module, ctx["protein"])
        value = row.get("call") if row else None
        return _leaf(value if value in (CALLED, NOT_CALLED) else NOT_ASSESSABLE, module)
    if kind == "flag":
        module, fld = arg.split(".", 1)
        row = _row(ctx["modules"], module, ctx["protein"])
        raw = row.get(fld, "") if row else ""
        value = {"1": CALLED, "0": NOT_CALLED}.get(str(raw).strip(), NOT_ASSESSABLE)
        return _leaf(value, module)
    if kind == "test":
        row = _row(ctx["modules"], arg["module"], ctx["protein"])
        limit = arg["value"]
        if isinstance(limit, str):
            limit = ctx["thresholds"][limit[1:]]
        try:
            number = float(row[arg["field"]]) if row else None
        except (KeyError, TypeError, ValueError):
            number = None
        if number is None or not math.isfinite(number):  # missing, empty, NaN or infinite
            return _leaf(NOT_ASSESSABLE, arg["module"])
        return _leaf(CALLED if _OPS[arg["op"]](number, limit) else NOT_CALLED, arg["module"])
    if kind == "ref":
        callee = ctx["cfg"].call_by_name(arg)
        label = ctx["label"] if callee.get("per_variant") else ""
        return ctx["results"][(arg, label)]
    if kind == "not":
        child = _eval(arg, ctx)
        return _Result(k_not(child.value), child.contributors)
    children = [_eval(c, ctx) for c in arg]
    values = [c.value for c in children]
    if kind == "and":
        out = k_and(*values)
        decisive = [c for c in children if out != NOT_CALLED or c.value == NOT_CALLED]
    else:
        out = k_or(*values)
        decisive = [c for c in children if out != CALLED or c.value == CALLED]
    if out == NOT_ASSESSABLE:  # an unknown result has no deciding module
        return _Result(out, frozenset())
    return _Result(out, frozenset().union(*(c.contributors for c in decisive)))


def _eval_other(call, ctx):
    surface = ctx["results"][(call["surface"], ctx["label"])]
    if surface.value == NOT_ASSESSABLE:
        return _Result(NOT_ASSESSABLE, frozenset()), ""
    if surface.value != call["surface_is"]:
        return _Result(NOT_CALLED, surface.contributors), ""
    contributors, left_out = set(surface.contributors), []
    for name in call["mechanism"]:
        res = ctx["results"][(name, "")]
        if res.value == CALLED:
            return _Result(NOT_CALLED, frozenset(res.contributors)), ""
        if res.value == NOT_ASSESSABLE:
            left_out.append(name)
        else:
            contributors |= res.contributors
    return _Result(CALLED, frozenset(contributors)), ",".join(left_out)


def available_variants(cfg, modules):
    return [v for v in cfg.step1_variants if v in modules]


def evaluate(cfg, protein_ids, taxa, modules, status_of):
    """Evaluate every call for every protein; return a list of ``CallRecord``.

    ``taxa`` maps protein ID to taxon ID. ``modules`` maps module name to ``ModuleTable``.
    ``status_of(module, taxon)`` returns ``(status, basis)``.
    """
    variants = available_variants(cfg, modules)
    status_cache, records = {}, []

    def status_for(module, taxon):
        key = (module, taxon)
        if key not in status_cache:
            status_cache[key] = status_of(module, taxon)
        return status_cache[key]

    for pid in protein_ids:
        taxon = taxa[pid]
        results = {}
        for call in cfg.calls:
            if call.get("per_variant"):
                labels = [(v, variant_label(v)) for v in variants]
            else:
                labels = [(None, "")]
            for step1, label in labels:
                ctx = {
                    "protein": pid,
                    "modules": modules,
                    "thresholds": cfg.thresholds,
                    "cfg": cfg,
                    "results": results,
                    "step1": step1,
                    "label": label,
                }
                other_basis = ""
                if call.get("kind") == "other":
                    res, other_basis = _eval_other(call, ctx)
                else:
                    res = _eval(call["expr"], ctx)
                results[(call["name"], label)] = res
                if res.contributors:
                    names = sorted(res.contributors)
                    pairs = [status_for(m, taxon) for m in names]
                    status = weakest(s for s, _ in pairs)
                    basis = ";".join(f"{m}:{b}" for m, (_, b) in zip(names, pairs, strict=True))
                else:
                    status, basis = UNVALIDATED, ""
                records.append(
                    CallRecord(pid, call["name"], label, res.value, status, basis, other_basis)
                )
    return records


def referenced_modules(cfg):
    """Names of all modules the rules use (``{step1}`` expanded to the step 1 variants)."""
    names = set()

    def walk(node):
        kind, arg = next(iter(node.items()))
        if kind in ("and", "or"):
            for child in arg:
                walk(child)
        elif kind == "not":
            walk(arg)
        elif kind == "call":
            names.add(arg)
        elif kind == "flag":
            names.add(arg.split(".", 1)[0])
        elif kind == "test":
            names.add(arg["module"])

    for call in cfg.calls:
        if call.get("kind") != "other":
            walk(call["expr"])
    expanded = set()
    for n in names:
        if "{step1}" in n:
            expanded.update(n.replace("{step1}", v) for v in cfg.step1_variants)
        else:
            expanded.add(n)
    return sorted(expanded)


def collect_evidence(cfg, protein_ids, modules):
    """Rows (protein, module, field, value) for the configured evidence fields."""
    rows = []
    for item in cfg.evidence:
        module, fld = item.split(".", 1)
        table = modules.get(module)
        if table is None:
            continue
        for pid in protein_ids:
            row = table.get(pid)
            if row and row.get("state", OK_STATE) == OK_STATE and str(row.get(fld, "")).strip():
                rows.append((pid, module, fld, str(row[fld]).strip()))
    return rows
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_engine.py -q`
Expected: `41 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-core
git add src/cellsurface_sorting_hat/categories.yaml src/cellsurface_sorting_hat/engine.py tests/cellsurface_sorting_hat/test_engine.py
git commit -m "feat(sorting-hat): category engine and packaged rules

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 6: Output writers

**Files:**
- Create: `src/cellsurface_sorting_hat/outputs.py`
- Test: `tests/cellsurface_sorting_hat/test_outputs.py`

**Interfaces:**
- Consumes: `write_atomic` (Task 4); `CallRecord` (Task 5); the call constants (Task 1).
- Produces: `RunInfo` (dataclass, fields as in the code, including `calibration`), `LONG_COLUMNS`, `write_atomic_text(path, text)`, `write_long(path, records)`, `write_wide(path, records, protein_ids)`, `write_evidence(path, rows)`, `write_proteins(path, proteins, taxa)`, `write_run_json(path, info)`, `render_report(info, records) -> str`, `write_report(path, info, records)`, `KNOWN_LIMITS`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_outputs.py` with exactly this content:

```python
"""Unit tests of the writers, without the command line."""

import csv
import gzip

from cellsurface_sorting_hat.engine import CallRecord
from cellsurface_sorting_hat.outputs import (
    LONG_COLUMNS,
    RunInfo,
    render_report,
    write_long,
    write_wide,
)


def info(**over):
    base = {
        "version": "x",
        "taxa": {40: 2},
        "taxonomy_sha256": "t",
        "config_sha256": "c",
        "default_gate": "step1_rule@R0",
        "thresholds": {"antigen_percentile_max": 15},
        "n_proteins": 2,
        "n_invalid": 0,
        "invalid": [],
        "n_trailing_stop": 0,
        "module_states": {"step1_rule@R0": "ok"},
    }
    base.update(over)
    return RunInfo(**base)


def rec(protein, call, variant, value, status="unvalidated", basis="", other=""):
    return CallRecord(protein, call, variant, value, status, basis, other)


def read(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def test_long_table_has_the_documented_columns(tmp_path):
    write_long(tmp_path / "l.tsv.gz", [rec("A", "adhesion_repeat", "", "called", "smoke", "m:t")])
    rows = read(tmp_path / "l.tsv.gz")
    assert rows[0] == LONG_COLUMNS
    assert rows[1] == ["A", "adhesion_repeat", "", "called", "smoke", "m:t", ""]


def test_wide_table_names_gated_columns_with_the_variant(tmp_path):
    records = [
        rec("A", "surface_glycoprotein", "R0", "called", "estimated"),
        rec("A", "other_not_surface", "R0", "called", "estimated", other="antigen_candidate"),
        rec("B", "surface_glycoprotein", "R0", "not_called"),
        rec("B", "other_not_surface", "R0", "not_assessable"),
    ]
    write_wide(tmp_path / "w.tsv.gz", records, ["A", "B"])
    rows = read(tmp_path / "w.tsv.gz")
    assert rows[0] == [
        "protein",
        "surface_glycoprotein[R0]",
        "surface_glycoprotein[R0]_status",
        "other_not_surface[R0]",
        "other_not_surface[R0]_status",
        "other_not_surface[R0]_basis",
    ]
    assert rows[1][1:3] == ["called", "estimated"] and rows[1][-1] == "antigen_candidate"


def test_report_warns_about_error_partial_and_inconsistent_modules():
    text = render_report(
        info(
            module_states={"a": "error", "b": "partial", "c": "ok"},
            module_notes={"a": "no result table"},
            inconsistent={"c": 2},
        ),
        [rec("A", "adhesion_repeat", "", "called")],
    )
    assert "WARNING: Module a is in state error. no result table" in text
    assert "WARNING: Module b is in state partial." in text
    assert "WARNING: Module c gives different rows to 2 group(s)" in text


def test_report_truncates_a_long_invalid_list():
    invalid = [[f"P{i}", "internal stop codon"] for i in range(25)]
    text = render_report(info(n_invalid=25, invalid=invalid), [rec("A", "x", "", "called")])
    assert "- P19: internal stop codon" in text and "- P20:" not in text
    assert "... and 5 more (see proteins.tsv.gz)" in text


def test_report_counts_values_per_call_and_variant():
    records = [
        rec("A", "surface_glycoprotein", "R0", v) for v in ("called", "called", "not_called")
    ]
    text = render_report(info(), records)
    assert "| surface_glycoprotein | R0 | 2 | 1 | 0 |" in text
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_outputs.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.outputs`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/outputs.py` with exactly this content:

```python
"""Write calls.long.tsv.gz, calls.wide.tsv.gz, evidence.tsv.gz, proteins.tsv.gz, report.md, run.json."""

import csv
import gzip
import io
import json
from collections import Counter
from dataclasses import asdict, dataclass, field

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.logic import CALLED, NOT_ASSESSABLE, NOT_CALLED

LONG_COLUMNS = ["protein", "call", "variant", "value", "status", "status_basis", "other_basis"]
MAX_LISTED_INVALID = 20

KNOWN_LIMITS = [
    "GPI-anchored and secreted enzymes, and non-adhesive structural wall proteins, get only "
    "`surface_glycoprotein`. There is no `cell_wall_protein` call in version 1.",
    "`surface_glycoprotein` is defined by GO cell wall and extracellular region evidence. It is not "
    'evidence of glycosylation. With the step 1 rule R0 the call means "SignalP calls a signal '
    'peptide" and nothing more.',
    "CFEM is filed under adhesion because class 2b-i is. Its confirmed fold is a hemophore. Binding "
    "to a host receptor is not shown.",
    "The repeat detectors have no clade truth set. Their calls are hypotheses.",
    "The antigen call is the top 15% of a fixed Coccidioides ranking. It is a weak label, and the "
    "ranking prints NOT CALIBRATED (3 of 4 anchors pass the top-decile test).",
    "Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, "
    "moonlighting proteins and biofilm are not categories.",
    "The taxon you give is recorded as given. It is not checked against the sequences.",
]


@dataclass
class RunInfo:
    version: str
    taxa: dict  # taxon ID -> number of proteins
    taxonomy_sha256: str
    config_sha256: str
    default_gate: str
    thresholds: dict
    n_proteins: int
    n_invalid: int
    invalid: list  # [id, note] pairs
    n_trailing_stop: int
    module_states: dict  # module name -> run state
    module_notes: dict = field(default_factory=dict)  # module name -> text
    state_counts: dict = field(default_factory=dict)  # module -> {protein state: count}
    absent_modules: list = field(default_factory=list)  # needed by the rules, not found
    unavailable_variants: list = field(default_factory=list)
    inconsistent: dict = field(default_factory=dict)  # module -> groups of identical sequences
    map_ids_not_in_fasta: int = 0
    calibration: list = field(default_factory=list)  # rows from cli.calibration_rows


def _gz(text):
    return gzip.compress(text.encode(), mtime=0)


def _tsv(header, rows):
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return _gz(buf.getvalue())


def column_name(call, variant):
    return f"{call}[{variant}]" if variant else call


def write_long(path, records):
    rows = [
        [r.protein, r.call, r.variant, r.value, r.status, r.status_basis, r.other_basis]
        for r in records
    ]
    write_atomic(path, _tsv(LONG_COLUMNS, rows))


def write_wide(path, records, protein_ids):
    columns, table = [], {pid: {} for pid in protein_ids}
    for r in records:
        name = column_name(r.call, r.variant)
        if name not in columns:
            columns.append(name)
        row = table[r.protein]
        row[name] = r.value
        row[name + "_status"] = r.status
        if r.call.startswith("other_"):
            row[name + "_basis"] = r.other_basis
    header = ["protein"]
    for name in columns:
        header += [name, name + "_status"]
        if name.startswith("other_"):
            header.append(name + "_basis")
    rows = [[pid] + [table[pid].get(h, "") for h in header[1:]] for pid in protein_ids]
    write_atomic(path, _tsv(header, rows))


def write_evidence(path, rows):
    write_atomic(path, _tsv(["protein", "module", "field", "value"], rows))


def write_proteins(path, proteins, taxa):
    header = ["id", "sha256", "taxon", "state", "note", "trailing_stop", "ambiguous_fraction"]
    rows = [
        [
            p.id,
            p.sha256,
            taxa[p.id],
            p.state,
            p.note,
            int(p.trailing_stop),
            f"{p.ambiguous_fraction:.4f}",
        ]
        for p in proteins
    ]
    write_atomic(path, _tsv(header, rows))


def write_run_json(path, info):
    write_atomic(path, (json.dumps(asdict(info), indent=2, sort_keys=True) + "\n").encode())


def _warnings(info):
    out = []
    if info.default_gate in info.unavailable_variants:
        out.append(
            f"The default gate {info.default_gate} is not available. Gated calls are not written for it."
        )
    for module, state in sorted(info.module_states.items()):
        if state in ("error", "partial"):
            note = info.module_notes.get(module, "")
            out.append(f"Module {module} is in state {state}. {note}".strip())
    for module, n in sorted(info.inconsistent.items()):
        out.append(
            f"Module {module} gives different rows to {n} group(s) of identical sequences. "
            "Check the module."
        )
    if info.absent_modules:
        out.append(
            "Modules the rules need and that were not found: " + ", ".join(info.absent_modules)
        )
    if info.map_ids_not_in_fasta:
        out.append(f"{info.map_ids_not_in_fasta} ID(s) in --taxon-map are not in the FASTA.")
    return out


def render_report(info, records):
    counts = Counter((r.call, r.variant, r.value) for r in records)
    keys = sorted({(r.call, r.variant) for r in records})
    lines = [
        "# cellsurface_sorting_hat report",
        "",
        f"- version: {info.version}",
        f"- proteins: {info.n_proteins} ({info.n_invalid} invalid, excluded from all modules; "
        f"{info.n_trailing_stop} had a trailing `*`, removed)",
        f"- taxa (taxon ID: proteins): {', '.join(f'{t}: {n}' for t, n in sorted(info.taxa.items()))}",
        f"- taxonomy file sha256: {info.taxonomy_sha256}",
        f"- categories.yaml sha256: {info.config_sha256}",
        f"- default gate: {info.default_gate}",
        "- thresholds: " + ", ".join(f"{k} = {v}" for k, v in sorted(info.thresholds.items())),
        "",
    ]
    warnings = _warnings(info)
    if warnings:
        lines += ["## Warnings", ""] + [f"- **WARNING: {w}**" for w in warnings] + [""]
    if info.invalid:
        shown = info.invalid[:MAX_LISTED_INVALID]
        lines += ["## Invalid proteins", ""] + [f"- {i}: {note}" for i, note in shown]
        if len(info.invalid) > len(shown):
            lines.append(f"- ... and {len(info.invalid) - len(shown)} more (see proteins.tsv.gz)")
        lines.append("")
    lines += ["## Module run states", "", "| module | run state |", "|---|---|"]
    lines += [f"| {m} | {s} |" for m, s in sorted(info.module_states.items())]
    if info.unavailable_variants:
        lines += ["", "Unavailable step 1 variants: " + ", ".join(info.unavailable_variants)]
    if any(info.state_counts.values()):
        lines += [
            "",
            "## Proteins that are not `ok` in a module",
            "",
            "| module | state | proteins |",
            "|---|---|---|",
        ]
        for m, states in sorted(info.state_counts.items()):
            lines += [f"| {m} | {s} | {n} |" for s, n in sorted(states.items())]
    if info.calibration:
        lines += [
            "",
            "## Module calibration",
            "",
            "| module | taxon | status | calibration set | positives | negatives "
            "| sensitivity [95% CI] | specificity [95% CI] |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for c in info.calibration:
            n_pos = c["n_pos"] if c["n_pos"] != "" else "-"
            n_neg = c["n_neg"] if c["n_neg"] != "" else "-"
            lines.append(
                f"| {c['module']} | {c['taxon']} | {c['status']} | {c['calibration_set'] or '-'} "
                f"| {n_pos} | {n_neg} | {c['sensitivity']} | {c['specificity']} |"
            )
    lines += [
        "",
        "## Calls",
        "",
        "| call | variant | called | not_called | not_assessable |",
        "|---|---|---|---|---|",
    ]
    for call, variant in keys:
        c = [counts[(call, variant, v)] for v in (CALLED, NOT_CALLED, NOT_ASSESSABLE)]
        lines.append(f"| {call} | {variant or '-'} | {c[0]} | {c[1]} | {c[2]} |")
    basis = Counter(r.other_basis for r in records if r.call.startswith("other_") and r.other_basis)
    if basis:
        lines += ["", "## `other_basis` (categories left out because they were not assessable)", ""]
        lines += [f"- {b}: {n}" for b, n in sorted(basis.items())]
    lines += ["", "## Known limits", ""] + [f"{i}. {t}" for i, t in enumerate(KNOWN_LIMITS, 1)]
    return "\n".join(lines) + "\n"


def write_atomic_text(path, text):
    write_atomic(path, text.encode())


def write_report(path, info, records):
    write_atomic_text(path, render_report(info, records))
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_outputs.py -q`
Expected: `5 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-core
git add src/cellsurface_sorting_hat/outputs.py tests/cellsurface_sorting_hat/test_outputs.py
git commit -m "feat(sorting-hat): table, report and run-record writers

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 7: Command line driver

**Files:**
- Create: `src/cellsurface_sorting_hat/cli.py`
- Test: `tests/cellsurface_sorting_hat/test_cli.py`

**Interfaces:**
- Consumes: `read_fasta`, `NA_INVALID`, `FastaError` (Task 3); `Lineage`, `TaxonError` (Task 2); `ModuleIdentity`, `load_status_source`, `resolve_entry`, `UNVALIDATED` (Task 2); `ConfigError`, `ModuleTable`, `evaluate`, `load_config`, `collect_evidence`, `referenced_modules` (Task 5); the writers and `RunInfo` (Task 6).
- Produces: `main(argv=None) -> int` (0 on success, 2 on a user or input error with the message on stderr; program errors are not caught); `parse_args`, `run(args)`, `read_taxon_map` (names file and line, refuses a repeated ID), `assign_taxa`, `load_modules`, `check_identical_sequences`, `StatusResolver(workdir, identities, lineage)` (callable `(module, taxon) -> (status, basis)`; `.entry(module, taxon)`), `calibration_rows(resolver, modules, taxa)`, `RunError`, `InputError`.
- Output files in `--out`: `calls.long.tsv.gz` (columns `protein, call, variant, value, status, status_basis, other_basis`), `calls.wide.tsv.gz`, `evidence.tsv.gz`, `proteins.tsv.gz`, `report.md`, `run.json`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_cli.py` with exactly this content:

```python
"""End to end: FASTA + module tables + status sources -> calls, wide table, report."""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.cli import main

FASTA = (
    ">SOW1 SOWgp-like\nMKTAYIAKQRQISFVKSHFSRQ\n"
    ">ENZ1 enzyme\nMNLLPQWERTYIPASDFG\n"
    ">DUP1 same sequence as SOW1\nMKTAYIAKQRQISFVKSHFSRQ\n"
    ">BAD1 internal stop\nMKT*AYI\n"
    ">STAR1 CRLF, lowercase, trailing stop\nmktayhhhh*\n"
)
TAXA = {"SOW1": 41, "DUP1": 41, "ENZ1": 40, "BAD1": 40, "STAR1": 40}
VALID = ["SOW1", "ENZ1", "DUP1", "STAR1"]


def build(tmp_path, write_module):
    fasta = tmp_path / "p.faa"
    fasta.write_bytes(FASTA.encode())
    taxon_map = tmp_path / "taxa.tsv"
    taxon_map.write_text("".join(f"{k}\t{v}\n" for k, v in TAXA.items()))
    wd = tmp_path / "wd"

    def call(values):
        return [{"id": i, "state": "ok", "call": v} for i, v in zip(VALID, values, strict=True)]

    write_module(
        wd,
        "step1_rule@R0",
        call(["called", "not_called", "called", "called"]),
        status=[{"taxa": [40], "status": "estimated", "source": "phasec"}],
    )
    write_module(wd, "repeat02", call(["called", "not_called", "called", "not_called"]))
    write_module(wd, "repeat14", call(["not_called"] * 4))
    write_module(wd, "pfam_adhesion", [{"id": i, "state": "ok", "hit": "0"} for i in VALID])
    write_module(wd, "pfam_allergen", [{"id": i, "state": "ok", "hit": "0"} for i in VALID])
    write_module(
        wd,
        "antigen_lookup",
        [
            {"id": "SOW1", "state": "ok", "percentile": "7.1"},
            {"id": "DUP1", "state": "ok", "percentile": "7.1"},
            {"id": "ENZ1", "state": "not_applicable", "percentile": ""},
            {"id": "STAR1", "state": "not_applicable", "percentile": ""},
        ],
    )
    write_module(
        wd,
        "allergen_homology",
        [
            {"id": "SOW1", "state": "ok", "identity": "0", "coverage": "0"},
            {"id": "ENZ1", "state": "ok", "identity": "82", "coverage": "90"},
            {"id": "DUP1", "state": "ok", "identity": "0", "coverage": "0"},
            {"id": "STAR1", "state": "ok", "identity": "40", "coverage": "85"},
        ],
    )
    return fasta, taxon_map, wd


def run_cli(tmp_path, write_module, nodes_dmp, extra=None, **over):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    out = tmp_path / "out"
    argv = [
        "--fasta",
        str(fasta),
        "--taxon-map",
        str(taxon_map),
        "--taxdump",
        str(nodes_dmp),
        "--workdir",
        str(wd),
        "--out",
        str(out),
    ] + (extra or [])
    return main(argv), out, wd


def read_long(out):
    with gzip.open(out / "calls.long.tsv.gz", "rt") as fh:
        return {
            (r["protein"], r["call"], r["variant"]): r for r in csv.DictReader(fh, delimiter="\t")
        }


def test_golden_calls(tmp_path, write_module, nodes_dmp):
    code, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    assert code == 0
    got = read_long(out)

    def v(p, c, var=""):
        r = got[(p, c, var)]
        return r["value"], r["status"]

    # SOW1: Coccidioides (taxon 41), not covered by the step 1 estimate (tested on taxon 40)
    assert v("SOW1", "surface_glycoprotein", "R0") == ("called", "unvalidated")
    assert v("SOW1", "adhesion_repeat") == ("called", "unvalidated")
    assert v("SOW1", "cell_wall_adhesion_candidate", "R0") == ("called", "unvalidated")
    assert v("SOW1", "antigen_candidate") == ("called", "unvalidated")
    assert v("SOW1", "antigen_candidate_surface", "R0") == ("called", "unvalidated")
    assert v("SOW1", "allergen_candidate")[0] == "not_called"
    assert v("SOW1", "other_not_surface", "R0")[0] == "not_called"
    assert v("SOW1", "other_surface_no_mechanism", "R0")[0] == "not_called"
    assert got[("SOW1", "surface_glycoprotein", "R0")]["status_basis"] == (
        "step1_rule@R0:taxon not tested"
    )
    # ENZ1: A. fumigatus (40), step 1 estimated; antigen not applicable; allergen homology called
    assert v("ENZ1", "surface_glycoprotein", "R0") == ("not_called", "estimated")
    assert v("ENZ1", "antigen_candidate")[0] == "not_assessable"
    assert v("ENZ1", "antigen_candidate_surface", "R0") == ("not_called", "estimated")
    assert v("ENZ1", "allergen_candidate")[0] == "called"
    assert v("ENZ1", "other_not_surface", "R0")[0] == "not_called"
    # STAR1: surface, no mechanism called, antigen left out of the basis
    r = got[("STAR1", "other_surface_no_mechanism", "R0")]
    assert (r["value"], r["other_basis"]) == ("called", "antigen_candidate")
    # BAD1: invalid protein, every call unknown
    assert {r["value"] for k, r in got.items() if k[0] == "BAD1"} == {"not_assessable"}


def test_identical_sequences_with_different_ids_get_identical_calls(
    tmp_path, write_module, nodes_dmp
):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    got = read_long(out)
    keys = {(c, var) for (p, c, var) in got if p == "SOW1"}
    assert keys == {(c, var) for (p, c, var) in got if p == "DUP1"}
    for c, var in keys:
        a, b = got[("SOW1", c, var)], got[("DUP1", c, var)]
        assert (a["value"], a["status"], a["other_basis"]) == (
            b["value"],
            b["status"],
            b["other_basis"],
        )


def test_crlf_lowercase_trailing_stop_fasta_is_read(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    fasta.write_bytes(FASTA.replace("\n", "\r\n").encode())
    out = tmp_path / "out"
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon-map",
            str(taxon_map),
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(out),
        ]
    )
    assert code == 0
    assert read_long(out)[("STAR1", "surface_glycoprotein", "R0")]["value"] == "called"


def test_unknown_taxon_stops_the_run(tmp_path, write_module, nodes_dmp, capsys):
    fasta, _, wd = build(tmp_path, write_module)
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon",
            "999",
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(tmp_path / "out"),
        ]
    )
    assert code == 2
    assert "taxon 999" in capsys.readouterr().err


def test_a_protein_with_no_taxon_stops_the_run(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text("SOW1\t41\n")  # the others have no taxon and no --taxon
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon-map",
            str(taxon_map),
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(tmp_path / "out"),
        ]
    )
    assert code == 2
    assert "no taxon for protein" in capsys.readouterr().err


def test_taxon_map_overrides_taxon(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text("SOW1\t41\n")
    out = tmp_path / "out"
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon",
            "40",
            "--taxon-map",
            str(taxon_map),
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(out),
        ]
    )
    assert code == 0
    got = read_long(out)
    assert got[("SOW1", "surface_glycoprotein", "R0")]["status"] == "unvalidated"  # taxon 41
    assert got[("STAR1", "surface_glycoprotein", "R0")]["status"] == "estimated"  # --taxon 40


def test_partial_module_output_is_reported_and_missing_proteins_are_unknown(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(wd, "repeat14", [{"id": "SOW1", "state": "ok", "call": "not_called"}])
    out = tmp_path / "out"
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon-map",
            str(taxon_map),
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(out),
        ]
    )
    assert code == 0
    assert "| repeat14 | partial |" in (out / "report.md").read_text()
    got = read_long(out)
    # ENZ1: repeat02 false and repeat14 missing -> unknown, not false
    assert got[("ENZ1", "adhesion_repeat", "")]["value"] == "not_assessable"
    # SOW1: repeat02 true wins over anything
    assert got[("SOW1", "adhesion_repeat", "")]["value"] == "called"


def test_module_with_run_state_unavailable_is_treated_as_absent(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(wd, "step1_rule@R2", [], meta={"run_state": "unavailable"})
    out = tmp_path / "out"
    assert (
        main(
            [
                "--fasta",
                str(fasta),
                "--taxon-map",
                str(taxon_map),
                "--taxdump",
                str(nodes_dmp),
                "--workdir",
                str(wd),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    report = (out / "report.md").read_text()
    assert "| step1_rule@R2 | unavailable |" in report
    assert "Unavailable step 1 variants:" in report and "step1_rule@R2" in report
    assert not any(k[2] == "R2" for k in read_long(out))


def test_no_module_at_all_still_writes_a_report_with_a_warning(tmp_path, nodes_dmp):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKTAYI\n")
    (tmp_path / "wd").mkdir()
    out = tmp_path / "out"
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon",
            "40",
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(tmp_path / "wd"),
            "--out",
            str(out),
        ]
    )
    assert code == 0
    report = (out / "report.md").read_text()
    assert "WARNING: The default gate step1_rule@R0 is not available" in report
    got = read_long(out)
    assert {r["value"] for r in got.values()} == {"not_assessable"}
    assert not any(k[2] for k in got)  # no per-variant records without a step 1 variant


def test_wide_table_and_run_json(tmp_path, write_module, nodes_dmp):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    with gzip.open(out / "calls.wide.tsv.gz", "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["protein"] for r in rows] == ["SOW1", "ENZ1", "DUP1", "BAD1", "STAR1"]
    star = rows[-1]
    assert star["other_surface_no_mechanism[R0]"] == "called"
    assert star["other_surface_no_mechanism[R0]_basis"] == "antigen_candidate"
    assert star["surface_glycoprotein[R0]_status"] == "estimated"
    run = json.loads((out / "run.json").read_text())
    assert run["n_proteins"] == 5 and run["n_invalid"] == 1
    assert run["taxa"] == {"40": 3, "41": 2}


def test_report_has_the_known_limits_and_counts(tmp_path, write_module, nodes_dmp):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    report = (out / "report.md").read_text()
    assert "## Known limits" in report
    assert 'means "SignalP calls a signal peptide"' in report
    assert "| surface_glycoprotein | R0 | 3 | 1 | 1 |" in report
    assert "- antigen_candidate: " in report


def _run(tmp_path, nodes_dmp, fasta, taxon_map, wd, out="out"):
    argv = [
        "--fasta",
        str(fasta),
        "--taxdump",
        str(nodes_dmp),
        "--workdir",
        str(wd),
        "--out",
        str(tmp_path / out),
    ]
    if taxon_map is not None:
        argv += ["--taxon-map", str(taxon_map)]
    return main(argv), tmp_path / out


def test_evidence_and_protein_tables(tmp_path, write_module, nodes_dmp):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    with gzip.open(out / "evidence.tsv.gz", "rt") as fh:
        ev = {
            (r["protein"], r["module"], r["field"]): r["value"]
            for r in csv.DictReader(fh, delimiter="\t")
        }
    assert ev[("SOW1", "antigen_lookup", "percentile")] == "7.1"
    assert ev[("ENZ1", "allergen_homology", "identity")] == "82"
    assert ("ENZ1", "antigen_lookup", "percentile") not in ev  # state not_applicable
    with gzip.open(out / "proteins.tsv.gz", "rt") as fh:
        prot = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert prot["SOW1"]["sha256"] == prot["DUP1"]["sha256"]
    assert prot["STAR1"]["trailing_stop"] == "1" and prot["SOW1"]["trailing_stop"] == "0"
    assert prot["BAD1"]["state"] == "na_invalid" and prot["SOW1"]["taxon"] == "41"


def test_report_shows_thresholds_invalid_proteins_and_trailing_stops(
    tmp_path, write_module, nodes_dmp
):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    report = (out / "report.md").read_text()
    assert "antigen_percentile_max = 15" in report and "allergen_identity_min = 70" in report
    assert "- BAD1: internal stop codon" in report
    assert "1 had a trailing `*`" in report
    assert "| antigen_lookup | not_applicable | 2 |" in report


def test_identical_sequences_with_different_module_rows_are_flagged(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(
        wd,
        "repeat02",
        [
            {"id": i, "state": "ok", "call": "called" if i == "SOW1" else "not_called"}
            for i in VALID
        ],
    )
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert (
        "Module repeat02 gives different rows to 1 group(s) of identical sequences"
        in (out / "report.md").read_text()
    )


def test_a_bad_call_value_is_unknown_and_counted(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(
        wd,
        "repeat14",
        [
            {"id": i, "state": "ok", "call": "Called" if i == "SOW1" else "not_called"}
            for i in VALID
        ],
    )
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert "| repeat14 | bad_value | 1 |" in (out / "report.md").read_text()
    got = read_long(out)
    assert got[("SOW1", "adhesion_repeat", "")]["value"] == "called"  # repeat02 is true
    assert got[("STAR1", "adhesion_repeat", "")]["value"] == "not_called"


def test_a_module_whose_ids_match_no_protein_is_an_error(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(
        wd, "repeat14", [{"id": f"sp|{i}|X", "state": "ok", "call": "not_called"} for i in VALID]
    )
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    report = (out / "report.md").read_text()
    assert "| repeat14 | error |" in report
    assert "no FASTA ID matches the 4 row(s)" in report
    assert "WARNING: Module repeat14 is in state error" in report


def test_module_table_without_a_state_column_stops_the_run(
    tmp_path, write_module, nodes_dmp, capsys
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    with gzip.open(wd / "modules" / "repeat14.tsv.gz", "wt") as fh:
        fh.write("id\tcall\nSOW1\tcalled\n")
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert "repeat14.tsv.gz: missing column 'state'" in capsys.readouterr().err


def test_duplicate_rows_in_a_module_table_stop_the_run(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    with gzip.open(wd / "modules" / "repeat14.tsv.gz", "wt") as fh:
        fh.write("id\tstate\tcall\nSOW1\tok\tcalled\nSOW1\tok\tnot_called\n")
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert "duplicate row for ID 'SOW1'" in capsys.readouterr().err


def test_a_run_record_without_a_table_is_listed_with_its_state(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "modules" / "step1_rule@R1.json").write_text(json.dumps({"run_state": "unavailable"}))
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert "| step1_rule@R1 | unavailable |" in (out / "report.md").read_text()


def test_absent_required_modules_are_listed(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "modules" / "pfam_adhesion.tsv.gz").unlink()
    (wd / "modules" / "pfam_adhesion.json").unlink()
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert (
        "Modules the rules need and that were not found: pfam_adhesion"
        in (out / "report.md").read_text()
    )
    assert read_long(out)[("SOW1", "adhesion_domain", "")]["value"] == "not_assessable"


@pytest.mark.parametrize(
    "text,message",
    [
        ("SOW1 41\n", "taxa.tsv:1: expected"),
        ("protein\ttaxon\n", "taxa.tsv:1: taxon ID is not an integer"),
    ],
)
def test_bad_taxon_map_lines_name_the_line(
    tmp_path, write_module, nodes_dmp, capsys, text, message
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text(text)
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert message in capsys.readouterr().err


def test_taxon_map_ids_not_in_the_fasta_are_reported(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    with open(taxon_map, "a") as fh:
        fh.write("GHOST\t40\n")
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert "1 ID(s) in --taxon-map are not in the FASTA" in (out / "report.md").read_text()


def test_a_broken_status_source_stops_the_run_and_names_the_file(
    tmp_path, write_module, nodes_dmp, capsys
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "status" / "step1_rule@R0.json").write_text(json.dumps({"module": "x"}))
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert "step1_rule@R0.json: not a valid status source" in capsys.readouterr().err


def test_a_program_error_is_not_hidden_as_a_user_error(
    tmp_path, write_module, nodes_dmp, monkeypatch
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    import cellsurface_sorting_hat.cli as cli

    def boom(*a, **k):
        raise KeyError("bug")

    monkeypatch.setattr(cli, "evaluate", boom)
    with pytest.raises(KeyError):
        _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)


def _module_bytes(wd, name, data):
    with gzip.open(wd / "modules" / f"{name}.tsv.gz", "wb") as fh:
        fh.write(data)


@pytest.mark.parametrize(
    "data,message",
    [
        (
            b"id\tstate\tcall\nSOW1\tok\tcalled\textra\n",
            "repeat14.tsv.gz:2: 4 fields, the header has 3",
        ),
        (b"id\tstate\tcall\nSOW1\tok\n", "repeat14.tsv.gz:2: 2 fields"),
    ],
)
def test_a_module_row_with_the_wrong_number_of_fields_names_file_and_line(
    tmp_path, write_module, nodes_dmp, capsys, data, message
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    _module_bytes(wd, "repeat14", data)
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and message in capsys.readouterr().err


def test_a_truncated_module_table_names_the_file(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    path = wd / "modules" / "repeat14.tsv.gz"
    path.write_bytes(path.read_bytes()[:20])
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and "repeat14.tsv.gz: cannot be read" in capsys.readouterr().err


@pytest.mark.parametrize("text", ["{not json", "[]"])
def test_a_bad_run_record_names_the_file(tmp_path, write_module, nodes_dmp, capsys, text):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "modules" / "repeat14.json").write_text(text)
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and "repeat14.json:" in capsys.readouterr().err


def test_a_failed_run_leaves_no_partial_output(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    _module_bytes(wd, "repeat14", b"id\tcall\nSOW1\tcalled\n")  # no state column
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and not out.exists()


def test_byte_order_marks_in_the_taxon_map_and_module_table_are_ignored(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_bytes(b"\xef\xbb\xbf" + taxon_map.read_bytes())
    path = wd / "modules" / "repeat14.tsv.gz"
    with gzip.open(path, "rb") as fh:
        raw = fh.read()
    _module_bytes(wd, "repeat14", b"\xef\xbb\xbf" + raw)
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert read_long(out)[("SOW1", "surface_glycoprotein", "R0")]["value"] == "called"


def test_a_repeated_id_in_the_taxon_map_names_both_lines(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text("SOW1\t41\nENZ1\t40\nSOW1\t40\n")
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert "taxa.tsv:3: ID 'SOW1' already on line 1" in capsys.readouterr().err


def test_an_unknown_taxon_stops_the_run_even_with_no_status_file(tmp_path, nodes_dmp, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKTAYI\n")
    (tmp_path / "wd").mkdir()
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon",
            "999",
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(tmp_path / "wd"),
            "--out",
            str(tmp_path / "out"),
        ]
    )
    assert code == 2 and "taxon 999" in capsys.readouterr().err


def test_calibration_is_shown_when_the_status_source_has_a_measure(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    status = json.loads((wd / "status" / "step1_rule@R0.json").read_text())
    status["entries"][0]["measure"] = {
        "calibration_set": "S1:all",
        "n_pos": 232,
        "n_neg": 4244,
        "sensitivity": {"value": 0.603, "lo": 0.55, "hi": 0.65},
        "specificity": {"value": 0.963, "lo": 0.95, "hi": 0.97},
    }
    (wd / "status" / "step1_rule@R0.json").write_text(json.dumps(status))
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    report = (out / "report.md").read_text()
    assert "## Module calibration" in report
    assert (
        "| step1_rule@R0 | 40 | estimated | S1:all | 232 | 4244 | 0.603 [0.550, 0.650] | 0.963 [0.950, 0.970] |"
        in report
    )
    # taxon 41 is not covered by the entry; a module with no status file says so
    assert (
        "| step1_rule@R0 | 41 | unvalidated | - | - | - | not measured | not measured |" in report
    )
    assert "| repeat02 | 40 | unvalidated | - | - | - | not measured | not measured |" in report
    run = json.loads((out / "run.json").read_text())
    row = next(c for c in run["calibration"] if c["module"] == "step1_rule@R0" and c["taxon"] == 40)
    assert row["n_pos"] == 232 and row["matched_taxon"] == 40
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_cli.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.cli`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/cli.py` with exactly this content:

```python
"""Command line driver: read module tables from a workdir, evaluate the rules, write outputs.

Version 1 of the driver does not run modules. It reads their outputs from
``<workdir>/modules/<module>.tsv.gz`` (columns ``id``, ``state``, ``call`` and module fields) and
optional ``<module>.json`` (run record) and ``<workdir>/status/<module>.json`` (status source).
"""

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import zlib
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from cellsurface_sorting_hat import __version__
from cellsurface_sorting_hat.engine import (
    ConfigError,
    ModuleTable,
    collect_evidence,
    evaluate,
    load_config,
    referenced_modules,
)
from cellsurface_sorting_hat.fasta import NA_INVALID, FastaError, read_fasta
from cellsurface_sorting_hat.outputs import (
    RunInfo,
    render_report,
    write_atomic_text,
    write_evidence,
    write_long,
    write_proteins,
    write_run_json,
    write_wide,
)
from cellsurface_sorting_hat.status import (
    UNVALIDATED,
    ModuleIdentity,
    load_status_source,
    resolve_entry,
)
from cellsurface_sorting_hat.taxonomy import Lineage, TaxonError

UNUSABLE_RUN_STATES = {"unavailable", "not_run", "error"}
CALL_VALUES = {"called", "not_called"}
FLAG_VALUES = {"0", "1"}


class RunError(RuntimeError):
    """The run cannot continue; the message is shown to the user."""


class InputError(RunError):
    """An input file is malformed; the message names the file."""


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="cellsurface_sorting_hat", description=__doc__)
    p.add_argument("--fasta", required=True)
    p.add_argument("--taxon", type=int, help="NCBI taxon ID for every protein not in --taxon-map")
    p.add_argument("--taxon-map", help="TSV: protein ID, taxon ID (overrides --taxon per protein)")
    p.add_argument("--taxdump", required=True, help="path to nodes.dmp from the NCBI taxonomy")
    p.add_argument("--workdir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--config", help="categories.yaml (default: the packaged file)")
    return p.parse_args(argv)


def read_taxon_map(path):
    out, first_line = {}, {}
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as err:
        raise InputError(f"{path}: cannot be read as text") from err
    for n, line in enumerate(text.splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = [f.strip() for f in line.split("\t")]
        if len(fields) < 2 or not fields[0]:
            raise InputError(f"{path}:{n}: expected 'protein ID<TAB>taxon ID'")
        try:
            taxon = int(fields[1])
        except ValueError:
            raise InputError(f"{path}:{n}: taxon ID is not an integer: {fields[1]!r}") from None
        if fields[0] in out:
            raise InputError(
                f"{path}:{n}: ID {fields[0]!r} already on line {first_line[fields[0]]}"
            )
        out[fields[0]] = taxon
        first_line[fields[0]] = n
    return out


def assign_taxa(proteins, default_taxon, taxon_map):
    taxa = {}
    for p in proteins:
        taxon = taxon_map.get(p.id, default_taxon)
        if taxon is None:
            raise RunError(f"no taxon for protein {p.id}: give --taxon or list it in --taxon-map")
        taxa[p.id] = taxon
    return taxa


@dataclass
class LoadedModules:
    tables: dict = field(default_factory=dict)
    states: dict = field(default_factory=dict)
    identities: dict = field(default_factory=dict)
    notes: dict = field(default_factory=dict)
    state_counts: dict = field(default_factory=dict)


def _read_json(path):
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except (json.JSONDecodeError, UnicodeDecodeError) as err:
        raise InputError(f"{path}: not valid JSON ({err.__class__.__name__})") from err
    if not isinstance(data, dict):
        raise InputError(f"{path}: must be a JSON object")
    return data


def _read_module_table(path):
    rows = {}
    try:
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            columns = next(reader, [])
            for needed in ("id", "state"):
                if needed not in columns:
                    raise InputError(f"{path.name}: missing column {needed!r}")
            for n, values in enumerate(reader, 2):
                if len(values) != len(columns):
                    raise InputError(
                        f"{path.name}:{n}: {len(values)} fields, the header has {len(columns)}"
                    )
                row = dict(zip(columns, values, strict=True))
                if row["id"] in rows:
                    raise InputError(f"{path.name}: duplicate row for ID {row['id']!r}")
                rows[row["id"]] = row
    except (OSError, EOFError, zlib.error, UnicodeDecodeError) as err:
        raise InputError(f"{path.name}: cannot be read ({err.__class__.__name__})") from err
    return columns, rows


def load_modules(workdir, protein_ids, invalid_ids):
    """Read module tables. Missing IDs, invalid proteins and bad values become unknown rows."""
    loaded = LoadedModules()
    folder = Path(workdir) / "modules"
    if not folder.exists():
        return loaded
    names = sorted(
        {p.name[: -len(".tsv.gz")] for p in folder.glob("*.tsv.gz")}
        | {p.stem for p in folder.glob("*.json")}
    )
    wanted = set(protein_ids)
    for name in names:
        meta_path = folder / f"{name}.json"
        meta = _read_json(meta_path) if meta_path.exists() else {}
        table_path = folder / f"{name}.tsv.gz"
        loaded.identities[name] = ModuleIdentity(
            name,
            str(meta.get("version", "")),
            str(meta.get("params_hash", "")),
            str(meta.get("artefact_hash", "")),
        )
        state = meta.get("run_state", "ok" if table_path.exists() else "error")
        loaded.states[name] = state
        if not table_path.exists():
            loaded.notes[name] = "no result table"
            continue
        if state in UNUSABLE_RUN_STATES:
            continue
        columns, rows = _read_module_table(table_path)
        matched = wanted & set(rows)
        if not matched:
            loaded.states[name] = "error"
            loaded.notes[name] = (
                f"no FASTA ID matches the {len(rows)} row(s); check the ID format of the module"
            )
            continue
        table, counts = ModuleTable(name), Counter()
        for pid in protein_ids:
            if pid in invalid_ids:
                row = {"state": NA_INVALID}
            elif pid not in rows:
                row = {"state": "error"}
            else:
                row = dict(rows[pid])
                if row["state"] == "ok" and _bad_value(columns, row):
                    row["state"] = "bad_value"
            table.rows[pid] = row
            if row["state"] != "ok":
                counts[row["state"]] += 1
        missing = sum(1 for pid in protein_ids if pid not in rows and pid not in invalid_ids)
        if missing:
            loaded.states[name] = "partial"
            loaded.notes[name] = f"{missing} protein(s) have no row"
        loaded.tables[name] = table
        loaded.state_counts[name] = dict(counts)
    return loaded


def _bad_value(columns, row):
    if "call" in columns and row.get("call") not in CALL_VALUES:
        return True
    if "hit" in columns and row.get("hit") not in FLAG_VALUES:
        return True
    return False


def check_identical_sequences(proteins, tables):
    """Per module, count groups of identical sequences whose module rows differ."""
    groups = defaultdict(list)
    for p in proteins:
        if p.state != NA_INVALID:
            groups[p.sha256].append(p.id)
    groups = [ids for ids in groups.values() if len(ids) > 1]
    out = {}
    for name, table in tables.items():
        bad = 0
        for ids in groups:
            seen = {
                tuple(sorted((k, v) for k, v in table.rows[i].items() if k != "id")) for i in ids
            }
            bad += len(seen) > 1
        if bad:
            out[name] = bad
    return out


class StatusResolver:
    """Resolve the status of a module for a taxon, and the matched calibration entry."""

    def __init__(self, workdir, identities, lineage):
        self.folder = Path(workdir) / "status"
        self.identities = identities
        self.lineage = lineage
        self._records = {}

    def _record(self, module):
        if module not in self._records:
            path = self.folder / f"{module}.json"
            try:
                self._records[module] = load_status_source(path) if path.exists() else None
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as err:
                raise InputError(f"{path}: not a valid status source ({err})") from err
        return self._records[module]

    def entry(self, module, taxon):
        """(entry, tested taxon, reason); reason is empty when an entry applies."""
        identity = self.identities.get(module)
        if identity is None:
            return None, None, "no module run record"
        return resolve_entry(self._record(module), identity, taxon, self.lineage)

    def __call__(self, module, taxon):
        entry, tested, reason = self.entry(module, taxon)
        if entry is None:
            return UNVALIDATED, reason
        return entry.status, f"taxon:{tested}"


def _rate_text(rate):
    if not rate:
        return "not measured"
    return f"{rate['value']:.3f} [{rate['lo']:.3f}, {rate['hi']:.3f}]"


def calibration_rows(resolver, modules, taxa):
    """One row per (module, taxon): status, calibration set, counts, sensitivity, specificity."""
    rows = []
    for module in sorted(modules):
        for taxon in sorted(set(taxa)):
            entry, tested, reason = resolver.entry(module, taxon)
            measure = (entry.measure if entry else None) or {}
            rows.append(
                {
                    "module": module,
                    "taxon": taxon,
                    "status": entry.status if entry else UNVALIDATED,
                    "matched_taxon": tested if entry else "",
                    "reason": reason,
                    "calibration_set": measure.get("calibration_set", ""),
                    "n_pos": measure.get("n_pos", ""),
                    "n_neg": measure.get("n_neg", ""),
                    "sensitivity": _rate_text(measure.get("sensitivity")),
                    "specificity": _rate_text(measure.get("specificity")),
                }
            )
    return rows


def run(args):
    cfg = load_config(args.config)
    lineage = Lineage.from_nodes_dmp(args.taxdump)
    proteins = read_fasta(args.fasta)
    taxon_map = read_taxon_map(args.taxon_map) if args.taxon_map else {}
    if args.taxon is None and not taxon_map:
        raise RunError("give --taxon or --taxon-map")
    taxa = assign_taxa(proteins, args.taxon, taxon_map)
    for taxon in set(taxa.values()):
        lineage.ancestors(taxon)  # raises TaxonError for an unknown taxon
    ids = [p.id for p in proteins]
    invalid = {p.id for p in proteins if p.state == NA_INVALID}
    loaded = load_modules(args.workdir, ids, invalid)
    tables = loaded.tables
    resolver = StatusResolver(args.workdir, loaded.identities, lineage)
    records = evaluate(cfg, ids, taxa, tables, resolver)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    counts = Counter(taxa.values())
    info = RunInfo(
        version=__version__,
        taxa=dict(counts),
        taxonomy_sha256=hashlib.sha256(Path(args.taxdump).read_bytes()).hexdigest(),
        config_sha256=cfg.sha256,
        default_gate=cfg.default_gate,
        thresholds=cfg.thresholds,
        n_proteins=len(proteins),
        n_invalid=len(invalid),
        invalid=[[p.id, p.note] for p in proteins if p.state == NA_INVALID],
        n_trailing_stop=sum(p.trailing_stop for p in proteins),
        module_states=loaded.states,
        module_notes=loaded.notes,
        state_counts=loaded.state_counts,
        absent_modules=[
            m
            for m in referenced_modules(cfg)
            if m not in tables and m not in cfg.step1_variants and m not in loaded.states
        ],
        unavailable_variants=[v for v in cfg.step1_variants if v not in tables],
        inconsistent=check_identical_sequences(proteins, tables),
        map_ids_not_in_fasta=len(set(taxon_map) - set(ids)),
        calibration=calibration_rows(resolver, tables, taxa.values()),
    )
    report = render_report(info, records)  # render first: a failure leaves no partial output
    write_long(out / "calls.long.tsv.gz", records)
    write_wide(out / "calls.wide.tsv.gz", records, ids)
    write_evidence(out / "evidence.tsv.gz", collect_evidence(cfg, ids, tables))
    write_proteins(out / "proteins.tsv.gz", proteins, taxa)
    write_atomic_text(out / "report.md", report)
    write_run_json(out / "run.json", info)
    return records


def main(argv=None):
    try:
        run(parse_args(argv))
    except (
        RunError,
        FastaError,
        TaxonError,
        ConfigError,
        OSError,
        subprocess.CalledProcessError,
    ) as err:
        print(f"cellsurface_sorting_hat: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_cli.py -q`
Expected: `35 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-core
git add src/cellsurface_sorting_hat/cli.py tests/cellsurface_sorting_hat/test_cli.py
git commit -m "feat(sorting-hat): command line driver

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 8: Packaging, CI, docs, whole-suite check

**Files:**
- Modify: `pyproject.toml`, `tests/surface_glyco/test_package.py`, `.github/workflows/build_and_test.yml`, `README.md`, `CHANGELOG.md`, `AGENTS.md`

**Interfaces:**
- Consumes: every earlier task.
- Produces: the command `cellsurface_sorting_hat` and the CI job that runs the tests.

- [ ] **Step 1: Add the package data and the script entry to `pyproject.toml`**

In `[project.scripts]`, add below the existing three lines:

```toml
cellsurface_sorting_hat = "cellsurface_sorting_hat.cli:main"
```

In `[tool.setuptools.package-data]`, add below the `surface_glyco` line:

```toml
"cellsurface_sorting_hat" = ["categories.yaml"]
```

(`known-first-party` was changed in Task 0.)

- [ ] **Step 2: Update the existing test that pins `pyproject.toml`**

`tests/surface_glyco/test_package.py::test_pyproject_names_match_the_decision` requires exactly three scripts and one package-data entry. Without this step the CI job fails. Replace the two assertions (lines 35 to 40) with:

```python
    assert cfg["project"]["scripts"] == {
        "surface_glyco_predict": "surface_glyco.scripts.predict:cli",
        "surface_glyco_train": "surface_glyco.scripts.train:cli",
        "surface_glyco_evaluate": "surface_glyco.scripts.evaluate:cli",
        "cellsurface_sorting_hat": "cellsurface_sorting_hat.cli:main",
    }
    assert cfg["tool"]["setuptools"]["package-data"] == {
        "surface_glyco": ["models/*"],
        "cellsurface_sorting_hat": ["categories.yaml"],
    }
```

Run it (it needs only the standard library):

```bash
PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/surface_glyco/test_package.py::test_pyproject_names_match_the_decision -q
```

Expected: `1 passed`. (Other tests in that file import `surface_glyco` and may need torch; run the whole `tests/surface_glyco` folder only in an environment that has the project dependencies.)

- [ ] **Step 3: Add the tests to the CI unit-test job**

In `.github/workflows/build_and_test.yml`, change the job name to `Unit tests (tests/surface_glyco and tests/cellsurface_sorting_hat, CPU)` and the run line to:

```yaml
        run: pytest tests/surface_glyco tests/cellsurface_sorting_hat
```

`tests/surface_glyco` has no `conftest.py`, so the one in `tests/cellsurface_sorting_hat` is the only one in that process.

- [ ] **Step 4: Add a short README section, a CHANGELOG line and an AGENTS.md entry**

Insert this into `README.md` before the final `# Author` heading (the file ends with that heading; text after it would nest under it). Use the heading level that matches the file:

```markdown
## cellsurface_sorting_hat (in development)

`cellsurface_sorting_hat` sorts the proteins of a proteome into cell surface categories
(surface glycoprotein, cell wall and adhesion candidate, antigen candidate, allergen candidate,
other). Design: `docs/superpowers/specs/2026-10-04-orchestrator-design.md`.

The core engine reads module result tables from a work directory:

    cellsurface_sorting_hat --fasta proteome.faa --taxon TAXON_ID --taxdump nodes.dmp \
        --workdir work --out out

It writes `calls.long.tsv.gz`, `calls.wide.tsv.gz`, `evidence.tsv.gz`, `proteins.tsv.gz`,
`report.md` and `run.json`. The modules that make the result tables (SignalP, Pfam, repeat
detectors, allergen homology, antigen lookup) are not part of this release. Calls from a module
that has no measurement are marked `unvalidated`.
```

(`TAXON_ID` is the NCBI taxon ID of the proteome. Look it up in the NCBI taxonomy.)

In `CHANGELOG.md`, under `## [Unreleased]` and `### Added`, add:

```markdown
- `cellsurface_sorting_hat` core engine: reads module result tables, applies three-valued category rules and writes calls, evidence, a report and run records. No module wrappers yet.
```

In `AGENTS.md`, add a short section after the "Model Review and Framework Plan" section. Replace `MODEL` with the model that carries out this plan:

```markdown
## cellsurface_sorting_hat core engine (October 2026)

### Agent: Claude Code (MODEL)
**Task:** Implement the core engine of the orchestrator from `docs/superpowers/plans/2026-10-04-cellsurface-sorting-hat-core.md`.

- Code: `src/cellsurface_sorting_hat/`; tests: `tests/cellsurface_sorting_hat/`.
- The engine reads module result tables. Wrappers that make them are a second plan.
```

- [ ] **Step 5: Run the whole suite and lint with the pinned ruff**

```bash
PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q
pre-commit run --files src/cellsurface_sorting_hat/*.py tests/cellsurface_sorting_hat/*.py pyproject.toml
```

Expected: `143 passed` (or `142 passed, 1 skipped` without `zstd`); pre-commit passes. If the pinned ruff 0.3.5 reformats a file, stage the result and commit again; do not change behavior.

- [ ] **Step 6: Mutation checks (confirm the tests can fail)**

Make each change in turn, run `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q`, confirm the named tests fail, then undo the change with `git checkout -- <file>`.

| Change | File | Tests that fail (measured) |
|---|---|---|
| In `k_and`, return `NOT_CALLED` when any input is `NOT_ASSESSABLE` | `logic.py` | 3 cases of `test_logic.py::test_and_or_tables`, and 5 tests elsewhere (8 in all) |
| In `weakest`, use `max` instead of `min` | `status.py` | `test_weakest_status`, `test_status_is_the_weakest_of_the_modules_that_decided_the_result` |
| In `resolve_status`, replace `if lineage.is_descendant_or_self(taxon, tested):` by `if True:` | `status.py` | 4 tests, among them `test_status_does_not_pass_to_a_sibling_clade_with_the_same_broad_label` and `test_golden_calls` |
| In `load_modules`, replace the line `loaded.states[name] = "partial"` by `pass` | `cli.py` | `test_partial_module_output_is_reported_and_missing_proteins_are_unknown` |
| In `_eval`, delete the two lines under `if out == NOT_ASSESSABLE:` | `engine.py` | `test_an_unknown_result_has_status_unvalidated_and_no_basis` |
| In `_locked`, replace the `fcntl.flock(...)` line by `pass` | `cache.py` | `test_parallel_updates_lose_no_rows_and_leave_a_readable_table` (a race test; run it three times) |
| In `run()`, replace `lineage.ancestors(taxon)` by `pass` | `cli.py` | `test_an_unknown_taxon_stops_the_run_even_with_no_status_file` |
| In `_eval`, replace the AND `decisive` list by `list(children)` | `engine.py` | `test_a_false_and_takes_its_status_from_the_false_inputs_only` |
| In `_eval`, replace `not math.isfinite(number)` by `number != number` | `engine.py` | `test_a_non_finite_number_is_unknown` (3 cases) |

- [ ] **Step 7: Install check and commit**

`pip` and `python` on the HPCC login node resolve to Python 3.9, which is below the project's `>=3.11`. Use the system Python 3.12 in a throwaway virtual environment:

```bash
/usr/bin/python3.12 -m venv "$SCRATCH/csh-venv"
"$SCRATCH/csh-venv/bin/pip" install -e . --no-deps   # needs network access for the build tools
"$SCRATCH/csh-venv/bin/pip" install pyyaml
"$SCRATCH/csh-venv/bin/python" -c "from cellsurface_sorting_hat.engine import load_config; print(len(load_config().calls))"   # expected: 11
"$SCRATCH/csh-venv/bin/cellsurface_sorting_hat" --help | head -3
git branch --show-current   # must print sorting-hat-core
git add pyproject.toml tests/surface_glyco/test_package.py .github/workflows/build_and_test.yml README.md CHANGELOG.md AGENTS.md
git commit -m "build(sorting-hat): script entry, package data, CI job, README, changelog

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


---

## Self-review against the spec

| Spec section | Where it is implemented |
|---|---|
| 3.2 modules: applicability separate from status; status by (module, version, taxon); stale `status_source` refused | Task 2 (`resolve_entry`, `resolve_status`), Task 7 (`StatusResolver`) |
| 3.2 taxa: lineage, tested taxa only, most specific wins, `--taxon-map` overrides | Task 2, Task 7 |
| 3.2 run states, `partial`, `unavailable`, JSON-only run records | Task 7 (`load_modules`), Task 6 (report) |
| 3.2 kind K (`not_in_reference`, not applicable) | Contract only: the engine treats any non-`ok` state as unknown (Task 5). The lookup module and its ID mapping (95% identity and 90% mutual coverage as proposed defaults) are Plan 2. |
| 3.3 Kleene tables, status columns, deciding modules, unknown has no basis | Task 1, Task 5, Task 6, Task 7 |
| 3.4 categories, gated and ungated calls, `other_*` with `other_basis`, thresholds, `allergen_homolog_hit` | Task 5 (`categories.yaml`, `_eval_other`) |
| 3.4 known limits in the report header; thresholds printed with every run | Task 6 (`KNOWN_LIMITS`, `render_report`) |
| 3.1 evidence output (antigen axes, allergen hit fields, Cys-rich, expression) | Task 5 (`collect_evidence`), Task 6 (`write_evidence`), Task 7 |
| 3.6 cache key, per-sequence cache, atomic writes, lock, protein key (ID, sha256), output schema, unavailable variants omitted | Task 3, Task 4, Task 6, Task 7 (`proteins.tsv.gz`, identical-sequence check) |
| Owner comment: where each module is calibrated and its Sn/Sp | Task 2 (`measure` in the status source), Task 7 (`calibration_rows`), Task 6 (report section "Module calibration" and `run.json`). The numbers themselves come from Plan 2. |
| 3.8 failure modes: trailing `*` (counted), internal `*`, non-residues, ambiguous residues, empty file, duplicate IDs, partial output, no valid proteins, bad module values, no ID match, sequence before header | Task 3, Task 7 |
| 3.8 GPU out of memory retry, module timeouts, PredGPI `too_short`, ESM window | Plan 2 (they belong to module wrappers) |
| 4 acceptance: software correctness on fixtures, Kleene cases, lineage tests, identical sequences | Tasks 1, 2, 4, 5, 6, 7 |
| 4 run-level check on Af293 on HPCC | Plan 2 |
| 4.1 five proteomes | Plan 2 (downloads and the run) |
| 5 statistics (bootstrap, kappa with 2x2) | Plan 2 (report-only, needs truth tables) |

Placeholder scan: none. All names used in later tasks are defined in earlier ones (`ModuleIdentity`, `Lineage`, `weakest`, `read_fasta`, `write_atomic`, `ModuleTable`, `evaluate`, `RunInfo`, `collect_evidence`, `referenced_modules`).

Known gaps in this plan, stated so the reviewer can check them:
1. `cache.ModuleCache` and `identity_key` are tested but not used by the driver. The driver reads module tables; the cache is used by the Plan 2 wrappers. The plan keeps it here so the key rules have one owner. The driver cannot yet expand per-sha256 results to IDs; it checks instead that identical sequences have identical module rows.
2. `ambiguous_fraction` is computed and written to `proteins.tsv.gz`; module wrappers decide their own limit (spec 3.8).
3. The packaged `categories.yaml` lists the step 1 variants `R1`, `R2` and `ml@card` as names. They are unavailable until their modules exist, so no records are written for them.
4. `merged.dmp` of the NCBI taxonomy is not read. A retired taxon ID stops the run with the taxon named. A status entry that lists a taxon not in the taxonomy never matches and gives no warning.
5. `fcntl` locks may not be reliable on every network file system. The lock test was run on node-local storage only.
6. The `--taxon` value is not checked against the sequences.
7. A steady stream of readers can delay a cache writer (shared locks are served at once). The driver does not use the cache yet.
8. The report gives counts of `not_assessable` per call. It does not give the reason or the share per call (spec section 5, "applicability hides most of the genome"). This is a Plan 2 item.

## Plan 2 (not written)

A second plan, after this one is merged, covers the pieces that need HPCC or external data:

1. Wrappers that write module tables: SignalP 6 and PredGPI with the R0 rule (`step1_rule@R0`), Pfam scan with the family table (`pfam_adhesion`, `pfam_allergen`), repeat detectors (`repeat02`, `repeat14`), allergen homology against the WHO/IUIS fungal set (`allergen_homology`, with `aligned_length` and `allergen_name`), antigen lookup with ID mapping (`antigen_lookup`, with the antigen axes), the Cys-rich finder tier (`cys_rich.tier`), expression (`expression.log2fc`), TM evidence (Phobius, TMHMM).
2. Measurement JSONs for `status_source` from the Phase C harness, filled with the `measure` object (calibration set, truth source, counts, sensitivity and specificity with intervals) for every module and taxon where truth exists. For each module, Plan 2 states the calibration set, the operating point (threshold) and which taxa it covers. Where no truth exists the entry stays `unvalidated`.
3. Download of the *A. fumigatus* A1163 (UniProt UP000001699) and W72310 (NCBI GCA_040167795.1) proteomes; check of the S288C protein count.
4. The run-level check on *A. fumigatus* Af293 on HPCC, with wall time and resources.
5. Report items that need truth tables: recall and false-positive rate per category and clade, kappa with 2x2 counts.
6. SLURM job scripts (`$SCRATCH`, no `BASH_SOURCE`), polling, resubmission, GPU out-of-memory retry, resource table.
7. A stored `calls` golden file made from real module outputs.
