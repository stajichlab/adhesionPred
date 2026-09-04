# Kingdom-wide Fungal Adhesion Survey Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Aggregate the already-computed `results/*.adhesion_predict.csv` files against Fungi_5k taxonomy (`samples.csv`) into a per-species table, run taxonomic-rank comparisons (naive + genus-averaged + mixed-effects) with figures, and produce a citable Markdown report showing whether specific fungal clades carry more adhesion-predicted proteins than others.

**Architecture:** A standalone `analysis/kingdom_survey/` directory, independent of the `src/adhesion_predict` package. Three importable modules (`join.py`, `stats.py`, `figures.py`) hold all logic and are unit-tested with small synthetic fixtures; three numbered CLI scripts drive them over the real 5,803-species dataset; a final report-assembly script plus a short manual write-up produces `REPORT.md`.

**Tech Stack:** Python 3.9, pandas, numpy, scipy, seaborn/matplotlib, statsmodels (already installed), scikit-posthocs (new dependency), pytest.

**Spec:** `docs/superpowers/specs/2026-09-04-kingdom-adhesion-survey-design.md` (as revised after Opus design review — read it before starting; this plan implements it section by section).

## Global Constraints

- Python `>=3.9` (per `pyproject.toml`) — no PEP 604 `X | None` union syntax; use `typing.Optional`/`typing.Union`.
- Ruff: line-length 100, rules `E,F,W,I,N,UP,B,C4`, ignoring `E501,N803,N806,N999` — run `ruff check analysis/` before each commit.
- Never re-run the model or modify `src/adhesion_predict/` — this plan only reads `results/` and external Fungi_5k paths.
- Real, non-negotiable counts from the spec/design-review that every task's output must reconcile against: **5,803** result files; **10** missing-result species, all genus *Colletotrichum*; **exactly 1** mismatched-annotation species, `Neohortaea_acidophila_CBS_113389`; SUBCLASS dropped from all aggregation (>50% blank).
- Paths: `RESULTS_DIR = <repo>/results`, `INPUT_DIR = /bigdata/stajichlab/shared/projects/Fungi_5k/input`, `SAMPLES_CSV = /bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv`.

---

## File Structure

```
analysis/kingdom_survey/
  join.py                     # join/reconciliation logic (Task 2)
  stats.py                    # aggregation + KW/Dunn's + genus-avg + mixedlm (Task 4)
  figures.py                  # all plotting functions (Task 6)
  01_build_species_table.py   # CLI: writes tables/species_adhesion_summary.csv + QC tables (Task 3)
  02_taxonomic_stats.py       # CLI: writes stats_by_*.csv, *_genus_avg.csv, mixedlm_by_*.csv (Task 5)
  03_make_figures.py          # CLI: writes figures/*.png (Task 7)
  04_generate_report.py       # CLI: assembles REPORT.md skeleton from real table/figure output (Task 8)
  tables/                     # generated, gitignored except a .gitkeep
  figures/                    # generated, gitignored except a .gitkeep
  REPORT.md                   # final report (committed)
tests/kingdom_survey/
  conftest.py                 # adds analysis/kingdom_survey to sys.path for import
  test_join.py                # Task 2
  test_stats.py                # Task 4
  test_figures.py              # Task 6
```

`tables/` and `figures/` outputs are large/regeneratable — gitignored (with a `.gitkeep` each) except `REPORT.md`, which is the citable deliverable and is committed along with the final table/figure set referenced by it (see Task 8).

---

### Task 1: Add `scikit-posthocs` dependency

**Files:**
- Modify: `pyproject.toml`

**Interfaces:**
- Produces: `scikit_posthocs` importable in the environment, for Task 4.

- [ ] **Step 1: Add the dependency**

Edit `pyproject.toml`'s `[project]` `dependencies` list (currently ends `"importlib-metadata; python_version<'3.8'",`) to add a new line:

```toml
    "scikit-posthocs>=0.7",
```

- [ ] **Step 2: Install it in the active environment**

Run: `pip install scikit-posthocs>=0.7`

- [ ] **Step 3: Verify import**

Run: `python -c "import scikit_posthocs; print(scikit_posthocs.__version__)"`
Expected: prints a version string, no `ModuleNotFoundError`.

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml
git commit -m "Add scikit-posthocs dependency for Dunn's post-hoc test"
```

---

### Task 2: `join.py` — reconciliation and join logic

**Files:**
- Create: `analysis/kingdom_survey/join.py`
- Create: `tests/kingdom_survey/conftest.py`
- Test: `tests/kingdom_survey/test_join.py`

**Interfaces:**
- Produces (used by Task 3):
  - `locustag_from_id(seq_id: str) -> str`
  - `first_id_from_fai(fai_path: Path) -> str`
  - `count_fai_lines(fai_path: Path) -> int`
  - `first_id_from_result_csv(result_csv_path: Path) -> Optional[str]`
  - `count_result_rows(result_csv_path: Path) -> int`
  - `probability_stats(result_csv_path: Path) -> Tuple[float, float]` (mean, median)
  - `load_samples_taxonomy(samples_csv_path: Path) -> Dict[str, dict]`
  - `build_species_table(results_dir: Path, input_dir: Path, taxonomy: Dict[str, dict]) -> Tuple[List[dict], List[dict], List[dict]]` (matched, unmatched, mismatched)
  - `find_missing_results(results_dir: Path, input_dir: Path) -> List[dict]`

- [ ] **Step 1: Create the test scaffolding**

Create `tests/kingdom_survey/conftest.py`:

```python
"""Make analysis/kingdom_survey importable as plain modules for tests."""
import sys
from pathlib import Path

KINGDOM_SURVEY_DIR = Path(__file__).resolve().parents[2] / "analysis" / "kingdom_survey"
sys.path.insert(0, str(KINGDOM_SURVEY_DIR))
```

- [ ] **Step 2: Write the failing tests**

Create `tests/kingdom_survey/test_join.py`:

```python
"""Tests for analysis/kingdom_survey/join.py."""
from join import (
    build_species_table,
    count_fai_lines,
    count_result_rows,
    find_missing_results,
    first_id_from_fai,
    first_id_from_result_csv,
    load_samples_taxonomy,
    locustag_from_id,
    probability_stats,
)


def test_locustag_from_id():
    assert locustag_from_id("F07B100A_000481-T1") == "F07B100A"


def test_first_id_from_fai(tmp_path):
    fai = tmp_path / "sp.proteins.fa.fai"
    fai.write_text("LOC1_000001-T1\t100\t10\t60\t61\nLOC1_000002-T1\t100\t150\t60\t61\n")
    assert first_id_from_fai(fai) == "LOC1_000001-T1"


def test_count_fai_lines(tmp_path):
    fai = tmp_path / "sp.proteins.fa.fai"
    fai.write_text("a\tb\tc\td\te\n" * 5)
    assert count_fai_lines(fai) == 5


def test_first_id_from_result_csv(tmp_path):
    result = tmp_path / "sp.adhesion_predict.csv"
    result.write_text("id,prediction,probability_adhesion\nLOC1_000002-T1,Adhesion,0.9\n")
    assert first_id_from_result_csv(result) == "LOC1_000002-T1"


def test_count_result_rows(tmp_path):
    result = tmp_path / "sp.adhesion_predict.csv"
    result.write_text(
        "id,prediction,probability_adhesion\n"
        "LOC1_000002-T1,Adhesion,0.9\n"
        "LOC1_000003-T1,Adhesion,0.6\n"
    )
    assert count_result_rows(result) == 2


def test_probability_stats(tmp_path):
    result = tmp_path / "sp.adhesion_predict.csv"
    result.write_text(
        "id,prediction,probability_adhesion\n"
        "LOC1_000002-T1,Adhesion,0.9\n"
        "LOC1_000003-T1,Adhesion,0.7\n"
    )
    mean, median = probability_stats(result)
    assert mean == 0.8
    assert median == 0.8


def test_load_samples_taxonomy(tmp_path):
    samples = tmp_path / "samples.csv"
    samples.write_text(
        "ASMID,SPECIESIN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,"
        "CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,LOCUSTAG\n"
        "ASM1,Sp one,strain1,PRJ1,1,dikarya,Ascomycota,Pezizomycotina,"
        "Sordariomycetes,,Hypocreales,Nectriaceae,Fusarium,Fusarium one,LOC1\n"
    )
    taxonomy = load_samples_taxonomy(samples)
    assert taxonomy["LOC1"]["phylum"] == "Ascomycota"
    assert taxonomy["LOC1"]["genus"] == "Fusarium"


def test_build_species_table_matched(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_one.proteins.fa.fai").write_text(
        "LOC1_000001-T1\t100\t10\t60\t61\n"
        "LOC1_000002-T1\t100\t10\t60\t61\n"
        "LOC1_000003-T1\t100\t10\t60\t61\n"
        "LOC1_000004-T1\t100\t10\t60\t61\n"
    )
    (results_dir / "Sp_one.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nLOC1_000002-T1,Adhesion,0.9\n"
    )
    taxonomy = {
        "LOC1": {
            "asmid": "ASM1", "phylum": "Ascomycota", "subphylum": "Pezizomycotina",
            "class": "Sordariomycetes", "order": "Hypocreales", "family": "Nectriaceae",
            "genus": "Fusarium", "species_name": "Fusarium one",
        }
    }
    matched, unmatched, mismatched = build_species_table(results_dir, input_dir, taxonomy)
    assert len(matched) == 1
    assert unmatched == []
    assert mismatched == []
    row = matched[0]
    assert row["locustag"] == "LOC1"
    assert row["total_proteins"] == 4
    assert row["adhesion_count"] == 1
    assert row["adhesion_fraction"] == 0.25


def test_build_species_table_mismatched_locustag(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_two.proteins.fa.fai").write_text("LOC2_000001-T1\t100\t10\t60\t61\n")
    (results_dir / "Sp_two.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nXP_999999.1,Adhesion,0.9\n"
    )
    taxonomy = {
        "LOC2": {
            "asmid": "ASM2", "phylum": "Ascomycota", "subphylum": None, "class": None,
            "order": None, "family": None, "genus": None, "species_name": "Sp two",
        }
    }
    matched, unmatched, mismatched = build_species_table(results_dir, input_dir, taxonomy)
    assert matched == []
    assert len(mismatched) == 1
    assert mismatched[0]["fai_locustag"] == "LOC2"
    assert mismatched[0]["result_locustag"] == "XP"


def test_build_species_table_unmatched_locustag(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_three.proteins.fa.fai").write_text("LOC3_000001-T1\t100\t10\t60\t61\n")
    (results_dir / "Sp_three.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nLOC3_000001-T1,Adhesion,0.9\n"
    )
    matched, unmatched, mismatched = build_species_table(results_dir, input_dir, {})
    assert matched == []
    assert len(unmatched) == 1
    assert "LOC3" in unmatched[0]["reason"]


def test_find_missing_results(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_one.proteins.fa.fai").write_text("LOC1_000001-T1\t100\t10\t60\t61\n")
    (input_dir / "Sp_two.proteins.fa.fai").write_text("LOC2_000001-T1\t100\t10\t60\t61\n")
    (results_dir / "Sp_one.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nLOC1_000001-T1,Adhesion,0.9\n"
    )
    missing = find_missing_results(results_dir, input_dir)
    assert len(missing) == 1
    assert missing[0]["stem"] == "Sp_two"
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `pytest tests/kingdom_survey/test_join.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'join'` (the module doesn't exist yet).

- [ ] **Step 4: Implement `join.py`**

Create `analysis/kingdom_survey/join.py`:

```python
"""Join logic for building the kingdom-wide adhesion species table.

Join key: the LOCUSTAG prefix embedded in protein ids (e.g.
"F07B100A_000481-T1" -> "F07B100A"). Anchored on each species' .fai index
(verified 100%-reliable during design review), with the result CSV's own
first-row LOCUSTAG used as a cross-check that catches annotation mismatches.
"""

import csv
import statistics
from pathlib import Path
from typing import Dict, List, Optional, Tuple

RESULT_SUFFIX = ".adhesion_predict.csv"
FAI_SUFFIX = ".proteins.fa.fai"


def locustag_from_id(seq_id: str) -> str:
    """Extract the LOCUSTAG prefix from a protein id like 'F07B100A_000481-T1'."""
    return seq_id.split("_", 1)[0]


def first_id_from_fai(fai_path: Path) -> str:
    """Return the id in the first column of the first line of a .fai index."""
    with open(fai_path) as fh:
        first_line = fh.readline()
    return first_line.split("\t", 1)[0]


def count_fai_lines(fai_path: Path) -> int:
    """Count total records (= total proteins) in a .fai index."""
    with open(fai_path) as fh:
        return sum(1 for _ in fh)


def first_id_from_result_csv(result_csv_path: Path) -> Optional[str]:
    """Return the id in the first data row of a result CSV, or None if header-only."""
    with open(result_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            return row["id"]
    return None


def count_result_rows(result_csv_path: Path) -> int:
    """Count adhesion-called proteins (data rows) in a result CSV."""
    with open(result_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        return sum(1 for _ in reader)


def probability_stats(result_csv_path: Path) -> Tuple[float, float]:
    """Return (mean, median) of probability_adhesion over adhesion-called proteins."""
    probs = []
    with open(result_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            probs.append(float(row["probability_adhesion"]))
    if not probs:
        return (float("nan"), float("nan"))
    return (statistics.mean(probs), statistics.median(probs))


def load_samples_taxonomy(samples_csv_path: Path) -> Dict[str, dict]:
    """Load samples.csv into a dict keyed by LOCUSTAG.

    Blank taxonomy fields are stored as None (not "") so downstream
    per-rank grouping can exclude a species from a rank rather than
    pooling it into a fake blank group. SUBCLASS is intentionally not
    loaded (>50% blank in samples.csv; dropped from analysis entirely
    per design review).
    """
    lookup: Dict[str, dict] = {}
    with open(samples_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            locustag = row["LOCUSTAG"].strip()
            if not locustag:
                continue
            lookup[locustag] = {
                "asmid": row["ASMID"],
                "phylum": row["PHYLUM"].strip() or None,
                "subphylum": row["SUBPHYLUM"].strip() or None,
                "class": row["CLASS"].strip() or None,
                "order": row["ORDER"].strip() or None,
                "family": row["FAMILY"].strip() or None,
                "genus": row["GENUS"].strip() or None,
                "species_name": row["SPECIES"].strip() or None,
            }
    return lookup


def build_species_table(
    results_dir: Path, input_dir: Path, taxonomy: Dict[str, dict]
) -> Tuple[List[dict], List[dict], List[dict]]:
    """Build the master species table plus QC rows.

    Returns (matched, unmatched, mismatched), where each is a list of
    dict rows. matched + unmatched + mismatched always sums to the
    number of result files found in results_dir.
    """
    matched: List[dict] = []
    unmatched: List[dict] = []
    mismatched: List[dict] = []

    for result_path in sorted(results_dir.glob(f"*{RESULT_SUFFIX}")):
        stem = result_path.name[: -len(RESULT_SUFFIX)]
        fai_path = input_dir / f"{stem}{FAI_SUFFIX}"

        if not fai_path.exists():
            unmatched.append({"stem": stem, "reason": "missing_fai"})
            continue

        fai_locustag = locustag_from_id(first_id_from_fai(fai_path))
        if fai_locustag not in taxonomy:
            unmatched.append(
                {"stem": stem, "reason": f"locustag_{fai_locustag}_not_in_samples_csv"}
            )
            continue

        result_first_id = first_id_from_result_csv(result_path)
        result_locustag = locustag_from_id(result_first_id) if result_first_id else None
        if result_locustag is not None and result_locustag != fai_locustag:
            mismatched.append(
                {
                    "stem": stem,
                    "fai_locustag": fai_locustag,
                    "result_locustag": result_locustag,
                }
            )
            continue

        tax = taxonomy[fai_locustag]
        total_proteins = count_fai_lines(fai_path)
        adhesion_count = count_result_rows(result_path)
        mean_prob, median_prob = probability_stats(result_path)

        matched.append(
            {
                "locustag": fai_locustag,
                "asmid": tax["asmid"],
                "species_name": tax["species_name"],
                "phylum": tax["phylum"],
                "subphylum": tax["subphylum"],
                "class": tax["class"],
                "order": tax["order"],
                "family": tax["family"],
                "genus": tax["genus"],
                "total_proteins": total_proteins,
                "adhesion_count": adhesion_count,
                "adhesion_fraction": adhesion_count / total_proteins,
                "mean_adhesion_prob": mean_prob,
                "median_adhesion_prob": median_prob,
            }
        )

    return matched, unmatched, mismatched


def find_missing_results(results_dir: Path, input_dir: Path) -> List[dict]:
    """Species with a .fai in the input dir but no corresponding result file."""
    result_stems = {p.name[: -len(RESULT_SUFFIX)] for p in results_dir.glob(f"*{RESULT_SUFFIX}")}
    missing = []
    for fai_path in sorted(input_dir.glob(f"*{FAI_SUFFIX}")):
        stem = fai_path.name[: -len(FAI_SUFFIX)]
        if stem not in result_stems:
            missing.append({"stem": stem, "fai_path": str(fai_path)})
    return missing
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/kingdom_survey/test_join.py -v`
Expected: all PASS.

- [ ] **Step 6: Lint**

Run: `ruff check analysis/kingdom_survey/join.py tests/kingdom_survey/`
Expected: no errors.

- [ ] **Step 7: Commit**

```bash
git add analysis/kingdom_survey/join.py tests/kingdom_survey/conftest.py tests/kingdom_survey/test_join.py
git commit -m "Add join/reconciliation logic for kingdom adhesion survey"
```

---

### Task 3: `01_build_species_table.py` — run the join on real data

**Files:**
- Create: `analysis/kingdom_survey/01_build_species_table.py`
- Create: `analysis/kingdom_survey/tables/.gitkeep`
- Modify: `.gitignore` (add `analysis/kingdom_survey/tables/*.csv`, `analysis/kingdom_survey/figures/*.png`, `analysis/kingdom_survey/figures/*.svg`)

**Interfaces:**
- Consumes: everything from Task 2 (`join.py`).
- Produces: `analysis/kingdom_survey/tables/species_adhesion_summary.csv` with columns `locustag, asmid, species_name, phylum, subphylum, class, order, family, genus, total_proteins, adhesion_count, adhesion_fraction, mean_adhesion_prob, median_adhesion_prob` — consumed by Task 5/6/7.
- Produces: `tables/unmatched_species.csv`, `tables/mismatched_locustag.csv`, `tables/missing_results.csv`.

- [ ] **Step 1: Update `.gitignore`**

Append to `.gitignore`:

```
analysis/kingdom_survey/tables/*.csv
analysis/kingdom_survey/figures/*.png
analysis/kingdom_survey/figures/*.svg
!analysis/kingdom_survey/tables/.gitkeep
!analysis/kingdom_survey/figures/.gitkeep
```

- [ ] **Step 2: Create placeholder dirs**

```bash
mkdir -p analysis/kingdom_survey/tables analysis/kingdom_survey/figures
touch analysis/kingdom_survey/tables/.gitkeep analysis/kingdom_survey/figures/.gitkeep
```

- [ ] **Step 3: Write the CLI script**

Create `analysis/kingdom_survey/01_build_species_table.py`:

```python
#!/usr/bin/env python
"""Build the master per-species adhesion summary table for the kingdom survey."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from join import build_species_table, find_missing_results, load_samples_taxonomy  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results"
INPUT_DIR = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/input")
SAMPLES_CSV = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv")
OUT_DIR = Path(__file__).resolve().parent / "tables"

MATCHED_FIELDS = [
    "locustag", "asmid", "species_name", "phylum", "subphylum", "class",
    "order", "family", "genus", "total_proteins", "adhesion_count",
    "adhesion_fraction", "mean_adhesion_prob", "median_adhesion_prob",
]
UNMATCHED_FIELDS = ["stem", "reason"]
MISMATCHED_FIELDS = ["stem", "fai_locustag", "result_locustag"]
MISSING_FIELDS = ["stem", "fai_path"]


def _write_csv(path: Path, rows: list, fieldnames: list) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    taxonomy = load_samples_taxonomy(SAMPLES_CSV)
    matched, unmatched, mismatched = build_species_table(RESULTS_DIR, INPUT_DIR, taxonomy)
    missing = find_missing_results(RESULTS_DIR, INPUT_DIR)

    _write_csv(OUT_DIR / "species_adhesion_summary.csv", matched, MATCHED_FIELDS)
    _write_csv(OUT_DIR / "unmatched_species.csv", unmatched, UNMATCHED_FIELDS)
    _write_csv(OUT_DIR / "mismatched_locustag.csv", mismatched, MISMATCHED_FIELDS)
    _write_csv(OUT_DIR / "missing_results.csv", missing, MISSING_FIELDS)

    n_results = sum(1 for _ in RESULTS_DIR.glob("*.adhesion_predict.csv"))
    print(f"Result files found: {n_results}")
    print(f"Matched: {len(matched)}  Unmatched: {len(unmatched)}  Mismatched: {len(mismatched)}")
    print(f"Missing results (fai present, no result file): {len(missing)}")

    assert len(matched) + len(unmatched) + len(mismatched) == n_results, (
        "Reconciliation failed: matched+unmatched+mismatched != result file count"
    )


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run it on real data**

Run: `python analysis/kingdom_survey/01_build_species_table.py`

Expected output includes `Result files found: 5803`, and the assertion passes (no `AssertionError`).

- [ ] **Step 5: Verify against the design-review numbers**

Run:
```bash
wc -l analysis/kingdom_survey/tables/missing_results.csv   # expect 11 (header + 10)
grep -c "Colletotrichum" <(cut -d, -f1 analysis/kingdom_survey/tables/missing_results.csv) || true
cat analysis/kingdom_survey/tables/mismatched_locustag.csv  # expect exactly 1 data row: Neohortaea_acidophila_CBS_113389
```
If `missing_results.csv` isn't 10 rows or all *Colletotrichum* (cross-check the stem names against `samples.csv`'s GENUS column, since `missing_results.csv` itself doesn't carry taxonomy — look up each stem's species name), or if `mismatched_locustag.csv` doesn't contain exactly the `Neohortaea_acidophila_CBS_113389` row, stop and re-examine `join.py` before proceeding — these are the exact numbers the design review verified by hand.

- [ ] **Step 6: Lint**

Run: `ruff check analysis/kingdom_survey/01_build_species_table.py`

- [ ] **Step 7: Commit**

```bash
git add analysis/kingdom_survey/01_build_species_table.py analysis/kingdom_survey/tables/.gitkeep analysis/kingdom_survey/figures/.gitkeep .gitignore
git commit -m "Add CLI to build master species adhesion table from real data"
```

---

### Task 4: `stats.py` — taxonomic aggregation and statistics

**Files:**
- Create: `analysis/kingdom_survey/stats.py`
- Test: `tests/kingdom_survey/test_stats.py`

**Interfaces:**
- Consumes: a `pandas.DataFrame` matching `species_adhesion_summary.csv`'s columns (Task 3's output, loaded via `pandas.read_csv`).
- Produces (used by Task 5 and Task 6):
  - `MIN_GROUP_N: int` (= 5)
  - `aggregate_summary(df, rank_col: str) -> pandas.DataFrame` with columns `[rank_col, "n_species", "median_fraction", "iqr_fraction", "median_count", "iqr_count", "median_prob", "small_n"]`
  - `kruskal_wallis_by_rank(df, rank_col: str, value_col: str) -> dict` with keys `rank, value_col, h_stat, p_value, epsilon_squared, n_groups, n_total`
  - `dunn_posthoc(df, rank_col: str, value_col: str) -> pandas.DataFrame` (square matrix, index/columns = group names)
  - `genus_average(df) -> pandas.DataFrame` with columns `[genus, phylum, class, order, family, n_species, adhesion_fraction, adhesion_count, mean_adhesion_prob]`
  - `mixedlm_by_rank(df, rank_col: str, value_col: str = "adhesion_fraction") -> dict` with keys `rank, value_col, converged, genus_variance, residual_variance, llf, summary`

- [ ] **Step 1: Write the failing tests**

Create `tests/kingdom_survey/test_stats.py`:

```python
"""Tests for analysis/kingdom_survey/stats.py."""
import numpy as np
import pandas as pd

from stats import (
    aggregate_summary,
    dunn_posthoc,
    genus_average,
    kruskal_wallis_by_rank,
    mixedlm_by_rank,
)


def _sample_df():
    rows = []
    rng = np.random.default_rng(42)
    for _ in range(6):
        rows.append(
            {
                "genus": "GenusA", "phylum": "PhylumX", "class": "ClassX", "order": "OrderX",
                "family": "FamilyX", "adhesion_fraction": 0.10 + rng.normal(0, 0.01),
                "adhesion_count": 10, "mean_adhesion_prob": 0.7,
            }
        )
    for _ in range(6):
        rows.append(
            {
                "genus": "GenusB", "phylum": "PhylumY", "class": "ClassY", "order": "OrderY",
                "family": "FamilyY", "adhesion_fraction": 0.30 + rng.normal(0, 0.01),
                "adhesion_count": 30, "mean_adhesion_prob": 0.8,
            }
        )
    return pd.DataFrame(rows)


def test_aggregate_summary():
    df = _sample_df()
    summary = aggregate_summary(df, "phylum")
    assert set(summary["phylum"]) == {"PhylumX", "PhylumY"}
    y_frac = summary.loc[summary["phylum"] == "PhylumY", "median_fraction"].iloc[0]
    x_frac = summary.loc[summary["phylum"] == "PhylumX", "median_fraction"].iloc[0]
    assert y_frac > x_frac


def test_kruskal_wallis_by_rank_detects_difference():
    df = _sample_df()
    result = kruskal_wallis_by_rank(df, "phylum", "adhesion_fraction")
    assert result["n_groups"] == 2
    assert result["p_value"] < 0.05


def test_kruskal_wallis_excludes_small_groups():
    df = _sample_df()
    small = pd.DataFrame(
        [
            {
                "genus": "GenusC", "phylum": "PhylumZ", "class": "ClassZ", "order": "OrderZ",
                "family": "FamilyZ", "adhesion_fraction": 0.99, "adhesion_count": 1,
                "mean_adhesion_prob": 0.99,
            }
        ]
    )
    df = pd.concat([df, small], ignore_index=True)
    result = kruskal_wallis_by_rank(df, "phylum", "adhesion_fraction")
    assert result["n_groups"] == 2  # PhylumZ (n=1) excluded, not counted


def test_dunn_posthoc_shape():
    df = _sample_df()
    posthoc = dunn_posthoc(df, "phylum", "adhesion_fraction")
    assert "PhylumX" in posthoc.index
    assert "PhylumY" in posthoc.columns


def test_genus_average_collapses_species():
    df = _sample_df()
    averaged = genus_average(df)
    assert len(averaged) == 2
    assert set(averaged["n_species"]) == {6}


def test_mixedlm_by_rank_runs():
    df = _sample_df()
    result = mixedlm_by_rank(df, "phylum", "adhesion_fraction")
    assert result["converged"] is True
    assert "genus_variance" in result
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/kingdom_survey/test_stats.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'stats'`.

- [ ] **Step 3: Implement `stats.py`**

Create `analysis/kingdom_survey/stats.py`:

```python
"""Taxonomic aggregation and statistics for the kingdom-wide adhesion survey.

Runs three parallel views at each rank, per the design's taxonomy-based
pseudo-phylogenetic correction: naive per-species tests (kruskal_wallis_by_rank/
dunn_posthoc directly on the input df), genus-averaged tests (same functions
called on genus_average(df)'s output), and a genus-random-intercept mixed
model (mixedlm_by_rank). None of this uses real branch lengths.
"""

import statsmodels.formula.api as smf
from scipy.stats import kruskal
import pandas as pd
import scikit_posthocs as sp

MIN_GROUP_N = 5


def aggregate_summary(df: pd.DataFrame, rank_col: str) -> pd.DataFrame:
    """Per-group N, median/IQR of adhesion_fraction and adhesion_count."""
    sub = df[df[rank_col].notna()]
    rows = []
    for group, gdf in sub.groupby(rank_col):
        rows.append(
            {
                rank_col: group,
                "n_species": len(gdf),
                "median_fraction": gdf["adhesion_fraction"].median(),
                "iqr_fraction": (
                    gdf["adhesion_fraction"].quantile(0.75)
                    - gdf["adhesion_fraction"].quantile(0.25)
                ),
                "median_count": gdf["adhesion_count"].median(),
                "iqr_count": (
                    gdf["adhesion_count"].quantile(0.75) - gdf["adhesion_count"].quantile(0.25)
                ),
                "median_prob": gdf["mean_adhesion_prob"].median(),
                "small_n": len(gdf) < MIN_GROUP_N,
            }
        )
    return pd.DataFrame(rows).sort_values("median_fraction", ascending=False).reset_index(
        drop=True
    )


def kruskal_wallis_by_rank(df: pd.DataFrame, rank_col: str, value_col: str) -> dict:
    """Kruskal-Wallis omnibus test across groups with N >= MIN_GROUP_N.

    Reports epsilon-squared effect size alongside p-value always — with
    ~5,800 species and wildly unbalanced phylum sample sizes, p-value alone
    is not a reliable signal of a biologically meaningful difference.
    """
    sub = df[df[rank_col].notna()]
    counts = sub[rank_col].value_counts()
    eligible_groups = counts[counts >= MIN_GROUP_N].index
    samples = [sub.loc[sub[rank_col] == g, value_col].to_numpy() for g in eligible_groups]

    if len(samples) < 2:
        return {
            "rank": rank_col, "value_col": value_col, "h_stat": float("nan"),
            "p_value": float("nan"), "epsilon_squared": float("nan"),
            "n_groups": len(samples), "n_total": sum(len(s) for s in samples),
        }

    h_stat, p_value = kruskal(*samples)
    n_total = sum(len(s) for s in samples)
    k = len(samples)
    epsilon_sq = (h_stat - k + 1) / (n_total - k) if n_total > k else float("nan")
    return {
        "rank": rank_col, "value_col": value_col, "h_stat": h_stat, "p_value": p_value,
        "epsilon_squared": epsilon_sq, "n_groups": k, "n_total": n_total,
    }


def dunn_posthoc(df: pd.DataFrame, rank_col: str, value_col: str) -> pd.DataFrame:
    """Dunn's post-hoc pairwise comparisons (BH-corrected) among N>=MIN_GROUP_N groups."""
    sub = df[df[rank_col].notna()].copy()
    counts = sub[rank_col].value_counts()
    eligible = counts[counts >= MIN_GROUP_N].index
    sub = sub[sub[rank_col].isin(eligible)]
    return sp.posthoc_dunn(sub, val_col=value_col, group_col=rank_col, p_adjust="fdr_bh")


def genus_average(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse to one row per genus: mean of numeric metrics, first taxonomy seen.

    This is the "genus-averaged" leg of the taxonomy-based pseudo-
    phylogenetic correction: it removes the dominant pseudoreplication
    source (a genus with many sequenced species dominating a coarser
    rank's median).
    """
    sub = df[df["genus"].notna()]
    agg = (
        sub.groupby("genus")
        .agg(
            phylum=("phylum", "first"),
            class_=("class", "first"),
            order=("order", "first"),
            family=("family", "first"),
            n_species=("adhesion_fraction", "size"),
            adhesion_fraction=("adhesion_fraction", "mean"),
            adhesion_count=("adhesion_count", "mean"),
            mean_adhesion_prob=("mean_adhesion_prob", "mean"),
        )
        .reset_index()
    )
    return agg.rename(columns={"class_": "class"})


def mixedlm_by_rank(df: pd.DataFrame, rank_col: str, value_col: str = "adhesion_fraction") -> dict:
    """Fit value_col ~ C(rank_col) with genus as a random intercept.

    Genus (not a fully nested order/family/genus structure) is used as
    the single random-effects level: statsmodels' MixedLM does not
    cleanly support 3-4 levels of nesting without variance-component
    machinery, and genus is the finest per-species grouping, capturing
    the largest share of non-independence. Still an approximation of a
    true phylogenetic comparative method (equal "branch lengths").
    """
    sub = df[df[rank_col].notna() & df["genus"].notna()].copy()
    if sub[rank_col].nunique() < 2 or sub["genus"].nunique() < 2:
        return {"rank": rank_col, "value_col": value_col, "converged": False}

    model = smf.mixedlm(f"{value_col} ~ C({rank_col})", sub, groups=sub["genus"])
    fit = model.fit(reml=False)
    return {
        "rank": rank_col,
        "value_col": value_col,
        "converged": bool(fit.converged),
        "genus_variance": float(fit.cov_re.iloc[0, 0]),
        "residual_variance": float(fit.scale),
        "llf": float(fit.llf),
        "summary": fit.summary().as_text(),
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/kingdom_survey/test_stats.py -v`
Expected: all PASS. (If `test_mixedlm_by_rank_runs` fails to converge with the small synthetic fixture, that's a fixture-realism issue, not a bug in the model call — increase the per-genus sample count in `_sample_df` to 12 and rerun before treating it as a real failure.)

- [ ] **Step 5: Lint**

Run: `ruff check analysis/kingdom_survey/stats.py tests/kingdom_survey/test_stats.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/kingdom_survey/stats.py tests/kingdom_survey/test_stats.py
git commit -m "Add taxonomic aggregation and stats module (naive + genus-avg + mixedlm)"
```

---

### Task 5: `02_taxonomic_stats.py` — run stats on real data

**Files:**
- Create: `analysis/kingdom_survey/02_taxonomic_stats.py`

**Interfaces:**
- Consumes: `stats.py` from Task 4; `tables/species_adhesion_summary.csv` from Task 3.
- Produces: `tables/stats_by_{phylum,class,order}.csv`, `tables/stats_by_{phylum,class,order}_posthoc.csv`, `tables/stats_by_{phylum,class,order}_genus_avg.csv` (+ posthoc), `tables/mixedlm_by_{phylum,class,order}.csv` — all consumed by Task 8's report assembly.

- [ ] **Step 1: Write the CLI script**

Create `analysis/kingdom_survey/02_taxonomic_stats.py`:

```python
#!/usr/bin/env python
"""Run taxonomic-rank statistics (naive, genus-averaged, mixed model) on the
master species table produced by 01_build_species_table.py."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats import (  # noqa: E402
    aggregate_summary,
    dunn_posthoc,
    genus_average,
    kruskal_wallis_by_rank,
    mixedlm_by_rank,
)

TABLES_DIR = Path(__file__).resolve().parent / "tables"
RANKS = ["phylum", "class", "order"]


def _run_rank(df: pd.DataFrame, rank: str, suffix: str) -> None:
    summary = aggregate_summary(df, rank)
    omnibus = kruskal_wallis_by_rank(df, rank, "adhesion_fraction")
    summary.to_csv(TABLES_DIR / f"stats_by_{rank}{suffix}.csv", index=False)
    pd.DataFrame([omnibus]).to_csv(
        TABLES_DIR / f"stats_by_{rank}{suffix}_omnibus.csv", index=False
    )
    if omnibus["p_value"] == omnibus["p_value"] and omnibus["p_value"] < 0.05:  # not NaN
        posthoc = dunn_posthoc(df, rank, "adhesion_fraction")
        posthoc.to_csv(TABLES_DIR / f"stats_by_{rank}{suffix}_posthoc.csv")
    print(f"{rank}{suffix}: n_groups={omnibus['n_groups']} p={omnibus['p_value']:.4g} "
          f"eps2={omnibus['epsilon_squared']:.4g}")


def main() -> None:
    df = pd.read_csv(TABLES_DIR / "species_adhesion_summary.csv")
    genus_df = genus_average(df)

    for rank in RANKS:
        _run_rank(df, rank, "")
        _run_rank(genus_df, rank, "_genus_avg")

        mixed = mixedlm_by_rank(df, rank)
        pd.DataFrame([{k: v for k, v in mixed.items() if k != "summary"}]).to_csv(
            TABLES_DIR / f"mixedlm_by_{rank}.csv", index=False
        )
        (TABLES_DIR / f"mixedlm_by_{rank}_summary.txt").write_text(
            mixed.get("summary", "did not converge")
        )
        print(f"{rank} mixedlm: converged={mixed['converged']}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python analysis/kingdom_survey/02_taxonomic_stats.py`

Expected: prints one line per rank/view with no exceptions; `tables/stats_by_phylum.csv`, `tables/stats_by_phylum_genus_avg.csv`, `tables/mixedlm_by_phylum.csv` (and class/order equivalents) all exist afterward.

- [ ] **Step 3: Sanity-check output**

Run: `python -c "import pandas as pd; df = pd.read_csv('analysis/kingdom_survey/tables/species_adhesion_summary.csv'); assert df['adhesion_fraction'].between(0,1).all(); assert (df['adhesion_count'] >= 0).all(); print('bounds OK, n =', len(df))"`

Expected: `bounds OK, n = <matched count from Task 3>`, no `AssertionError`.

- [ ] **Step 4: Lint**

Run: `ruff check analysis/kingdom_survey/02_taxonomic_stats.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/kingdom_survey/02_taxonomic_stats.py
git commit -m "Add CLI to run taxonomic stats (naive/genus-avg/mixedlm) on real data"
```

---

### Task 6: `figures.py` — plotting functions

Before writing this task's code, invoke the **dataviz** skill for palette/consistency guidance (the design spec requires it) — it governs color choices and mark specs even though these are static PNG/SVG, not an interactive artifact.

**Files:**
- Create: `analysis/kingdom_survey/figures.py`
- Test: `tests/kingdom_survey/test_figures.py`

**Interfaces:**
- Consumes: a `pandas.DataFrame` matching `species_adhesion_summary.csv`'s columns; `aggregate_summary()`'s output from Task 4 (for `ranked_summary_table`).
- Produces (used by Task 7):
  - `MIN_GROUP_N: int` (= 5, same threshold as `stats.py`)
  - `boxplot_by_rank(df, rank_col: str, value_col: str, out_path: Path, top_n: Optional[int] = None, title: str = "") -> None`
  - `scatter_proteome_vs_adhesion(df, out_path: Path) -> None`
  - `probability_boxplot_by_rank(df, rank_col: str, out_path: Path, top_n: Optional[int] = None) -> None`
  - `ranked_summary_table(summary_df, rank_col: str) -> str` (markdown table text)

- [ ] **Step 1: Write the failing tests**

Create `tests/kingdom_survey/test_figures.py`:

```python
"""Smoke tests for analysis/kingdom_survey/figures.py — verifies each
function runs without error and writes a non-empty file; does not
inspect pixel content."""
import pandas as pd

from figures import (
    boxplot_by_rank,
    probability_boxplot_by_rank,
    ranked_summary_table,
    scatter_proteome_vs_adhesion,
)
from stats import aggregate_summary


def _sample_df():
    rows = []
    for i in range(8):
        rows.append(
            {
                "phylum": "PhylumX" if i < 4 else "PhylumY",
                "adhesion_fraction": 0.1 + 0.01 * i,
                "adhesion_count": 10 + i,
                "total_proteins": 1000 + 50 * i,
                "mean_adhesion_prob": 0.6 + 0.01 * i,
            }
        )
    return pd.DataFrame(rows)


def test_boxplot_by_rank_creates_file(tmp_path):
    df = _sample_df()
    out = tmp_path / "box.png"
    boxplot_by_rank(df, "phylum", "adhesion_fraction", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_scatter_creates_file(tmp_path):
    df = _sample_df()
    out = tmp_path / "scatter.png"
    scatter_proteome_vs_adhesion(df, out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_probability_boxplot_creates_file(tmp_path):
    df = _sample_df()
    out = tmp_path / "prob.png"
    probability_boxplot_by_rank(df, "phylum", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_ranked_summary_table_is_markdown():
    df = _sample_df()
    summary = aggregate_summary(df, "phylum")
    md = ranked_summary_table(summary, "phylum")
    assert "PhylumX" in md
    assert "|" in md
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/kingdom_survey/test_figures.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'figures'`.

- [ ] **Step 3: Implement `figures.py`**

Create `analysis/kingdom_survey/figures.py`:

```python
"""Figure generation for the kingdom-wide adhesion survey report.

Every group figure annotates each group with its N and visually
desaturates groups with N < MIN_GROUP_N, since sample-size imbalance
across fungal phyla/orders is the main confound in this analysis and
needs to be visible in the figure itself, not only in a table.
"""
from pathlib import Path
from typing import List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

MIN_GROUP_N = 5
_FULL_COLOR = "#4C72B0"
_SMALL_N_COLOR = "#B0B0B0"


def _order_groups(df: pd.DataFrame, rank_col: str, value_col: str, top_n: Optional[int]) -> List:
    medians = df.groupby(rank_col)[value_col].median().sort_values(ascending=False)
    if top_n is not None:
        counts = df[rank_col].value_counts()
        top_groups = counts.sort_values(ascending=False).head(top_n).index
        medians = medians[medians.index.isin(top_groups)]
    return list(medians.index)


def boxplot_by_rank(
    df: pd.DataFrame,
    rank_col: str,
    value_col: str,
    out_path: Path,
    top_n: Optional[int] = None,
    title: str = "",
) -> None:
    """Ordered box plot of value_col across rank_col groups, N-annotated, small-N greyed."""
    sub = df[df[rank_col].notna()].copy()
    order = _order_groups(sub, rank_col, value_col, top_n)
    sub = sub[sub[rank_col].isin(order)]
    counts = sub[rank_col].value_counts()

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.5), 6))
    palette = {g: (_FULL_COLOR if counts[g] >= MIN_GROUP_N else _SMALL_N_COLOR) for g in order}
    sns.boxplot(
        data=sub, x=rank_col, y=value_col, order=order, hue=rank_col,
        palette=palette, legend=False, ax=ax,
    )
    labels = [f"{g} (n={counts[g]})" for g in order]
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(labels, rotation=90 if len(order) > 8 else 45, ha="right")
    ax.set_title(title)
    ax.set_xlabel(rank_col.capitalize())
    ax.set_ylabel(value_col.replace("_", " "))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def scatter_proteome_vs_adhesion(df: pd.DataFrame, out_path: Path) -> None:
    """Scatter of total_proteins vs adhesion_count, colored by phylum, with a
    reference line for 'expected count under the kingdom-wide median
    fraction' — separates genome-size effects from real enrichment."""
    sub = df[df["phylum"].notna()].copy()
    kingdom_median_fraction = sub["adhesion_fraction"].median()

    fig, ax = plt.subplots(figsize=(9, 7))
    sns.scatterplot(
        data=sub, x="total_proteins", y="adhesion_count", hue="phylum",
        alpha=0.6, s=20, ax=ax, legend="brief",
    )
    x_range = [sub["total_proteins"].min(), sub["total_proteins"].max()]
    ax.plot(
        x_range,
        [x * kingdom_median_fraction for x in x_range],
        color="black",
        linestyle="--",
        linewidth=1,
        label=f"expected @ kingdom median fraction ({kingdom_median_fraction:.4f})",
    )
    ax.set_xlabel("Total proteins in proteome")
    ax.set_ylabel("Adhesion-called proteins")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def probability_boxplot_by_rank(
    df: pd.DataFrame, rank_col: str, out_path: Path, top_n: Optional[int] = None
) -> None:
    """Distribution of mean_adhesion_prob by rank_col — secondary signal of
    'how confidently adhesive', same N-annotation/greying as boxplot_by_rank."""
    boxplot_by_rank(
        df, rank_col, "mean_adhesion_prob", out_path, top_n=top_n,
        title=f"Mean adhesion probability by {rank_col}",
    )


def ranked_summary_table(summary_df: pd.DataFrame, rank_col: str) -> str:
    """Render a sorted markdown table (replaces a min-max-normalized heatmap,
    which would visually inflate trivial differences over ~9 phyla)."""
    cols = [rank_col, "n_species", "median_fraction", "median_count", "median_prob"]
    sorted_df = summary_df.sort_values("median_fraction", ascending=False)[cols]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for _, row in sorted_df.iterrows():
        values = []
        for c in cols:
            v = row[c]
            values.append(f"{v:.4f}" if isinstance(v, float) else str(v))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/kingdom_survey/test_figures.py -v`
Expected: all PASS.

- [ ] **Step 5: Lint**

Run: `ruff check analysis/kingdom_survey/figures.py tests/kingdom_survey/test_figures.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/kingdom_survey/figures.py tests/kingdom_survey/test_figures.py
git commit -m "Add figure-generation module for kingdom adhesion survey"
```

---

### Task 7: `03_make_figures.py` — run figures on real data

**Files:**
- Create: `analysis/kingdom_survey/03_make_figures.py`

**Interfaces:**
- Consumes: `figures.py` (Task 6), `stats.py`'s `aggregate_summary`/`genus_average` (Task 4), `tables/species_adhesion_summary.csv` (Task 3).
- Produces: `figures/box_fraction_by_phylum.png`, `figures/box_fraction_by_order.png`, `figures/scatter_proteome_vs_adhesion.png`, `figures/box_prob_by_phylum.png`, `figures/box_prob_by_order.png` — consumed by Task 8.

- [ ] **Step 1: Write the CLI script**

Create `analysis/kingdom_survey/03_make_figures.py`:

```python
#!/usr/bin/env python
"""Generate all figures for the kingdom-wide adhesion survey report."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figures import (  # noqa: E402
    boxplot_by_rank,
    probability_boxplot_by_rank,
    scatter_proteome_vs_adhesion,
)

TABLES_DIR = Path(__file__).resolve().parent / "tables"
FIGURES_DIR = Path(__file__).resolve().parent / "figures"


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(TABLES_DIR / "species_adhesion_summary.csv")

    boxplot_by_rank(
        df, "phylum", "adhesion_fraction", FIGURES_DIR / "box_fraction_by_phylum.png",
        title="Adhesion fraction by phylum",
    )
    boxplot_by_rank(
        df, "order", "adhesion_fraction", FIGURES_DIR / "box_fraction_by_order.png",
        top_n=20, title="Adhesion fraction by order (top 20 by species count)",
    )
    scatter_proteome_vs_adhesion(df, FIGURES_DIR / "scatter_proteome_vs_adhesion.png")
    probability_boxplot_by_rank(df, "phylum", FIGURES_DIR / "box_prob_by_phylum.png")
    probability_boxplot_by_rank(
        df, "order", FIGURES_DIR / "box_prob_by_order.png", top_n=20,
    )
    print(f"Wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python analysis/kingdom_survey/03_make_figures.py`

- [ ] **Step 3: Visual sanity check**

Open each PNG in `analysis/kingdom_survey/figures/` (or use the Read tool to view them as images) and confirm: axis labels are present and readable, phylum/order names aren't overlapping illegibly, small-N groups are visibly greyed out, and the scatter plot's reference line looks like a plausible diagonal through the bulk of the points.

- [ ] **Step 4: Lint**

Run: `ruff check analysis/kingdom_survey/03_make_figures.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/kingdom_survey/03_make_figures.py
git commit -m "Add CLI to generate figures for kingdom adhesion survey"
```

---

### Task 8: Assemble and write `REPORT.md`

**Files:**
- Create: `analysis/kingdom_survey/04_generate_report.py`
- Create: `analysis/kingdom_survey/REPORT.md` (generated, then hand-edited)

**Interfaces:**
- Consumes: every table under `tables/` and every figure under `figures/` produced by Tasks 3, 5, 7.
- Produces: the final `REPORT.md` deliverable.

- [ ] **Step 1: Write the report-assembly script**

Create `analysis/kingdom_survey/04_generate_report.py` — this pulls real numbers out of the generated tables (no invented content) into a Markdown skeleton with an explicit placeholder section for the headline-findings paragraph, which Step 3 below fills in by hand after inspecting real results:

```python
#!/usr/bin/env python
"""Assemble REPORT.md from the tables/figures already generated by
01/02/03. Pulls all numbers directly from the CSVs — the only prose this
script does NOT generate is the headline-findings interpretation, which
is written by hand in Step 3 after inspecting the real output."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figures import ranked_summary_table  # noqa: E402

BASE = Path(__file__).resolve().parent
TABLES = BASE / "tables"


def _count_rows(csv_path: Path) -> int:
    return len(pd.read_csv(csv_path)) if csv_path.exists() else 0


def main() -> None:
    species = pd.read_csv(TABLES / "species_adhesion_summary.csv")
    n_matched = len(species)
    n_unmatched = _count_rows(TABLES / "unmatched_species.csv")
    n_mismatched = _count_rows(TABLES / "mismatched_locustag.csv")
    n_missing = _count_rows(TABLES / "missing_results.csv")

    phylum_summary = pd.read_csv(TABLES / "stats_by_phylum.csv")
    order_summary = pd.read_csv(TABLES / "stats_by_order.csv")
    phylum_omnibus = pd.read_csv(TABLES / "stats_by_phylum_omnibus.csv").iloc[0]
    order_omnibus = pd.read_csv(TABLES / "stats_by_order_omnibus.csv").iloc[0]
    phylum_omnibus_genus = pd.read_csv(TABLES / "stats_by_phylum_genus_avg_omnibus.csv").iloc[0]

    report = f"""# Kingdom-wide Fungal Adhesion Protein Survey

## Scope

- {n_matched} species successfully joined to taxonomy and proteome-size data.
- {n_unmatched} species excluded (no LOCUSTAG match or missing .fai) — see `tables/unmatched_species.csv`.
- {n_mismatched} species excluded for a mismatched protein-annotation source — see `tables/mismatched_locustag.csv`.
- {n_missing} species in Fungi_5k have no `adhesion_predict` result file at all — see `tables/missing_results.csv`.

## Headline findings

_(Written by hand in Task 8 Step 3, after inspecting the real per-rank
tables and figures below — naming the actual top/bottom phylum and order
by median adhesion_fraction, their epsilon-squared effect sizes, and
whether the naive/genus-averaged/mixed-model views agree.)_

![Adhesion fraction by phylum](figures/box_fraction_by_phylum.png)

![Proteome size vs adhesion count](figures/scatter_proteome_vs_adhesion.png)

## By phylum

Naive per-species test: H={phylum_omnibus['h_stat']:.3g}, p={phylum_omnibus['p_value']:.3g}, \
epsilon-squared={phylum_omnibus['epsilon_squared']:.3g}, n_groups={int(phylum_omnibus['n_groups'])}.

Genus-averaged test: H={phylum_omnibus_genus['h_stat']:.3g}, p={phylum_omnibus_genus['p_value']:.3g}, \
epsilon-squared={phylum_omnibus_genus['epsilon_squared']:.3g}.

{ranked_summary_table(phylum_summary, "phylum")}

![Mean adhesion probability by phylum](figures/box_prob_by_phylum.png)

## By order (top 20 by species count)

Naive per-species test: H={order_omnibus['h_stat']:.3g}, p={order_omnibus['p_value']:.3g}, \
epsilon-squared={order_omnibus['epsilon_squared']:.3g}, n_groups={int(order_omnibus['n_groups'])}.

{ranked_summary_table(order_summary, "order")}

![Adhesion fraction by order](figures/box_fraction_by_order.png)

![Mean adhesion probability by order](figures/box_prob_by_order.png)

## Family/genus drill-down

_(Written by hand: drill into the top 3 orders by epsilon-squared effect
size from the table above; label this section exploratory /
hypothesis-generating, per the pre-registered rule in the design spec.)_

## Caveats

1. **Non-independence beyond genus.** Genus-averaging and a genus-random-
   intercept mixed model correct for oversampled genera, but neither uses
   real branch lengths or divergence times — a true phylogenetic
   comparative method would. Residual non-independence above the genus
   level (shared family/order ancestry) is not corrected for here.
2. **Annotation quality varies by genome** and could confound raw counts;
   mitigated partly by using `adhesion_fraction`, not raw count alone.
3. **Classifier training-set bias.** The model was trained on a small
   positive set (FLO/ALS1-like proteins from *Saccharomyces*/*Candida*);
   it may under-call adhesion-like proteins in very divergent lineages
   whose domains don't resemble the training set. Treat phyla with
   unusually low fractions as a training-bias hypothesis to test, not
   automatically a biological conclusion.
4. **Taxonomic sampling is wildly unbalanced** (phylum N ranges from
   thousands down to 1); Kruskal-Wallis reads as significant on almost
   any real difference at this scale, so epsilon-squared effect size is
   the headline statistic here, not the p-value alone. Groups with
   N < {5} are shown for reference but excluded from formal tests and
   visually greyed out in figures.
5. **`adhesion_fraction` is fixed to the model's internal 0.5 probability
   threshold.** Results files contain only positive calls, so no
   threshold-sensitivity analysis is possible from existing data — a
   different cutoff could shift which clades look "enriched."

## Follow-up

A phylogeny-aware re-analysis (real branch lengths from `nf_phyling`
output or another tree source, e.g. phylogenetic ANOVA / PGLS) is the
natural next step once a tree is available for this species set.
"""
    (BASE / "REPORT.md").write_text(report)
    print(f"Wrote {BASE / 'REPORT.md'}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python analysis/kingdom_survey/04_generate_report.py`

- [ ] **Step 3: Fill in the two hand-written sections**

Open `analysis/kingdom_survey/REPORT.md`, read the actual `tables/stats_by_phylum.csv` and `tables/stats_by_order.csv` (already inlined as markdown tables in the report) plus the omnibus numbers already inlined, and replace:
- The `_(Written by hand ...)_` placeholder under "Headline findings" with a concrete 3-5 sentence paragraph naming the actual top and bottom phylum/order by median `adhesion_fraction`, their epsilon-squared values, and whether the naive vs. genus-averaged results agree or diverge (if they diverge for a given rank, say so explicitly and call it a likely oversampling artifact, per the design spec).
- The `_(Written by hand ...)_` placeholder under "Family/genus drill-down" with the top-3-by-epsilon-squared orders' family/genus breakdown (`aggregate_summary` on a filtered subset, or by reading `species_adhesion_summary.csv` directly for those orders), clearly labeled exploratory.

- [ ] **Step 4: Verify the report renders**

Run: `python -c "import pathlib; print(pathlib.Path('analysis/kingdom_survey/REPORT.md').read_text()[:200])"` to confirm the file is non-empty and starts as expected; open it in a Markdown previewer (or view the embedded PNGs directly) to confirm the images actually display from the relative paths.

- [ ] **Step 5: Commit**

```bash
git add analysis/kingdom_survey/04_generate_report.py analysis/kingdom_survey/REPORT.md
git add -f analysis/kingdom_survey/tables/species_adhesion_summary.csv analysis/kingdom_survey/tables/stats_by_*.csv analysis/kingdom_survey/tables/mixedlm_by_*.csv analysis/kingdom_survey/tables/unmatched_species.csv analysis/kingdom_survey/tables/mismatched_locustag.csv analysis/kingdom_survey/tables/missing_results.csv
git add -f analysis/kingdom_survey/figures/*.png
git commit -m "Add kingdom-wide fungal adhesion survey report, tables, and figures"
```

(The `-f` is needed because `tables/*.csv` and `figures/*.png` are gitignored by default per Task 3 — this is the one deliberate exception: the exact tables/figures the committed report cites are checked in so the report is self-contained and reproducible from this commit, even though the pipeline can regenerate them.)

---

## Self-Review Notes

- **Spec coverage:** every section of the design doc has a task — join logic (Task 2/3), blank-taxonomy/SUBCLASS handling (Task 2's `load_samples_taxonomy`), taxonomic aggregation + naive stats (Task 4/5), taxonomy-based pseudo-phylogenetic correction (Task 4's `genus_average`/`mixedlm_by_rank`), all 5 figures (Task 6/7, heatmap replaced with `ranked_summary_table` per the design's rejection of it), report + caveats (Task 8), verification plan (reconciliation asserts in Task 3, bounds check in Task 5).
- **Placeholder scan:** the only literal placeholder text (`_(Written by hand ...)_`) is inside the generated `REPORT.md`, not this plan, and Task 8 Step 3 gives it concrete, code-referenced instructions for what real content replaces it and where that content comes from — this is not deferred/vague work, it's work that genuinely requires the real numbers to exist first.
- **Type/name consistency:** checked `matched`/`unmatched`/`mismatched` return order is consistent across `join.py`'s docstring, tests, and `01_build_species_table.py`'s unpacking; `aggregate_summary`'s column names (`n_species`, `median_fraction`, `median_count`, `median_prob`) match what `ranked_summary_table` and `04_generate_report.py` read; `genus_average`'s output column name `genus` (not `Genus`) matches `mixedlm_by_rank`'s `groups=sub["genus"]`.
