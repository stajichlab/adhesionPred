# Adhesion Protein Molecular Properties Survey Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract length, amino-acid-composition, and functional-domain (Pfam/CAZy/MEROPS/signalp/tmhmm/targetp) properties for every adhesion-predicted protein and a matched per-species background sample, aggregate by taxonomic clade, and produce figures/tables/report showing how these properties differ from background and vary across clades.

**Architecture:** A standalone `analysis/adhesion_properties/` directory, parallel to and reusing generic pieces of `analysis/kingdom_survey/` (its `join.py` locustag-matching primitives, `stats.py`'s generic `kruskal_wallis_by_rank`, and `figures.py`'s validated color constants — all imported via `sys.path.insert`, not duplicated). A DuckDB read-only connection (`db.py`) queries `functionalDB/function.duckdb` for domain hits and sequence data, filtered against the specific protein-ID subset needed (never a full-table scan).

**Tech Stack:** Python 3.9, duckdb, pandas, numpy, scipy, statsmodels, matplotlib/seaborn (all already dependencies of `analysis/kingdom_survey/`), pytest.

**Spec:** `docs/superpowers/specs/2026-09-04-adhesion-protein-properties-design.md` — read it before starting; this plan implements it task by task.

## Global Constraints

- Python `>=3.9` — no PEP 604 `X | None` union syntax; use `typing.Optional`.
- Ruff: line-length 100, rules `E,F,W,I,N,UP,B,C4`, ignoring `E501,N803,N806,N999`.
- `functionalDB/function.duckdb` is 278GB and **read-only** — every query must filter against a specific protein-ID list (registered as a temp view/DataFrame), never scan a full table unfiltered.
- Data-source facts established during design review, binding on every task that touches them: `gene_proteins.protein_id` carries a `.protein` suffix the other functional tables don't; `pfam` has multiple rows per protein (avg 2.1) so presence/identity tallies need `DISTINCT protein_id`; `targetp` has no row at all for a protein with no predicted signal (absence of a row = "none", not a category value); `cazy_overview.cazyme_fam` embeds an alignment-range suffix (`GH114(59-261)`) that must be stripped before tallying; `merops_id` is a per-sequence accession, not a family code.
- `MIN_GROUP_N = 5` (imported from `kingdom_survey.stats`, not re-defined) governs which clades are eligible for formal statistical tests throughout.
- Reuse, don't duplicate: `kingdom_survey/join.py`'s `first_id_from_fai`/`locustag_from_id`, `kingdom_survey/stats.py`'s `MIN_GROUP_N`/`kruskal_wallis_by_rank`, and `kingdom_survey/figures.py`'s `_FULL_COLOR`/`_OTHER_COLOR` are imported via `sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "kingdom_survey"))` at the top of the module that needs them — never copy their code.

---

## File Structure

```
analysis/adhesion_properties/
  db.py                          # DuckDB connection + filtered-query helpers (Task 1)
  properties.py                  # AA composition / sequence property functions (Task 2)
  universe.py                    # adhesion set + background sample builder (Task 3)
  01_build_protein_universe.py   # CLI (Task 4)
  domains.py                     # domain-flag + top-domain-enrichment functions (Task 5)
  02_query_functional_domains.py # CLI (Task 6)
  03_compute_sequence_properties.py  # CLI (Task 7)
  stats_ext.py                   # per-clade aggregation + Mann-Whitney (Task 8)
  04_aggregate_and_test.py       # CLI (Task 9)
  figures_ext.py                 # plotting functions (Task 10)
  05_make_figures.py             # CLI (Task 11)
  06_generate_report.py          # CLI (Task 12)
  tables/
  figures/
  REPORT.md
tests/adhesion_properties/
  conftest.py                    # adds analysis/adhesion_properties to sys.path
  test_db.py
  test_properties.py
  test_universe.py
  test_domains.py
  test_stats_ext.py
  test_figures_ext.py
```

`tables/` and `figures/` outputs are gitignored except explicitly committed final artifacts (same convention as `kingdom_survey`).

---

### Task 1: `db.py` — DuckDB connection and filtered-query helpers

**Files:**
- Create: `analysis/adhesion_properties/db.py`
- Create: `tests/adhesion_properties/conftest.py`
- Test: `tests/adhesion_properties/test_db.py`

**Interfaces:**
- Produces (used by Tasks 5, 7):
  - `FUNCTION_DB_PATH: Path`
  - `connect(db_path: Path = FUNCTION_DB_PATH) -> duckdb.DuckDBPyConnection`
  - `fetch_lengths(con, protein_ids: Iterable[str]) -> pd.DataFrame` with columns `protein_id, length, peptide`
  - `fetch_domain_hits(con, table: str, protein_ids: Iterable[str], id_col: str = "protein_id") -> pd.DataFrame` (raw rows, no dedup)

- [ ] **Step 1: Write the failing tests**

Create `tests/adhesion_properties/conftest.py`:

```python
"""Make analysis/adhesion_properties importable as plain modules for tests."""
import sys
from pathlib import Path

ADHESION_PROPERTIES_DIR = Path(__file__).resolve().parents[2] / "analysis" / "adhesion_properties"
sys.path.insert(0, str(ADHESION_PROPERTIES_DIR))
```

Create `tests/adhesion_properties/test_db.py`:

```python
"""Tests for analysis/adhesion_properties/db.py — uses an in-memory DuckDB
with tiny fake tables matching the real functionalDB schema, never the
real 278GB database."""
import duckdb
import pandas as pd

from db import fetch_domain_hits, fetch_lengths


def _fake_connection():
    con = duckdb.connect(":memory:")
    con.execute(
        "CREATE TABLE gene_proteins (protein_id VARCHAR, length BIGINT, peptide VARCHAR)"
    )
    con.execute(
        "INSERT INTO gene_proteins VALUES "
        "('LOC1_000001-T1.protein', 4, 'SSTT'), "
        "('LOC1_000002-T1.protein', 3, 'PPP')"
    )
    con.execute("CREATE TABLE pfam (protein_id VARCHAR, pfam_id VARCHAR)")
    con.execute(
        "INSERT INTO pfam VALUES ('LOC1_000001-T1', 'DomainA'), ('LOC1_000001-T1', 'DomainB')"
    )
    return con


def test_fetch_lengths_strips_protein_suffix_and_matches_ids():
    con = _fake_connection()
    result = fetch_lengths(con, ["LOC1_000001-T1", "LOC1_000002-T1"])
    result = result.sort_values("protein_id").reset_index(drop=True)
    assert list(result["protein_id"]) == ["LOC1_000001-T1", "LOC1_000002-T1"]
    assert list(result["length"]) == [4, 3]
    assert list(result["peptide"]) == ["SSTT", "PPP"]


def test_fetch_lengths_only_returns_requested_ids():
    con = _fake_connection()
    result = fetch_lengths(con, ["LOC1_000001-T1"])
    assert len(result) == 1
    assert result.iloc[0]["protein_id"] == "LOC1_000001-T1"


def test_fetch_domain_hits_returns_matching_rows():
    con = _fake_connection()
    result = fetch_domain_hits(con, "pfam", ["LOC1_000001-T1"])
    assert len(result) == 2
    assert set(result["pfam_id"]) == {"DomainA", "DomainB"}


def test_fetch_domain_hits_rejects_unknown_table():
    con = _fake_connection()
    try:
        fetch_domain_hits(con, "not_a_real_table", ["LOC1_000001-T1"])
        assert False, "expected ValueError"
    except ValueError:
        pass
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/adhesion_properties/test_db.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'db'`.

- [ ] **Step 3: Implement `db.py`**

Create `analysis/adhesion_properties/db.py`:

```python
"""DuckDB query helpers for the functional annotation database.

functionalDB/function.duckdb is 278GB and opened read-only; every
function here filters against a caller-supplied protein-id list
(registered as a temp view) rather than scanning a full table.
"""

from pathlib import Path
from typing import Iterable

import duckdb
import pandas as pd

FUNCTION_DB_PATH = Path(
    "/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb"
)

_ALLOWED_DOMAIN_TABLES = {"pfam", "cazy_overview", "merops", "signalp", "tmhmm", "targetp"}


def connect(db_path: Path = FUNCTION_DB_PATH) -> "duckdb.DuckDBPyConnection":
    """Open a read-only connection to the functional annotation database."""
    return duckdb.connect(str(db_path), read_only=True)


def fetch_lengths(con, protein_ids: Iterable[str]) -> pd.DataFrame:
    """Return protein_id (without the .protein suffix) + length + peptide
    for the given ids, joined against gene_proteins.

    gene_proteins.protein_id carries a '.protein' suffix (confirmed during
    design review); this strips it so the returned protein_id matches the
    convention every other table and this analysis's own tables use.
    """
    ids_df = pd.DataFrame({"protein_id": list(protein_ids)})
    con.register("_wanted_ids", ids_df)
    query = """
        SELECT
            replace(gp.protein_id, '.protein', '') AS protein_id,
            gp.length AS length,
            gp.peptide AS peptide
        FROM gene_proteins gp
        JOIN _wanted_ids w ON gp.protein_id = w.protein_id || '.protein'
    """
    result = con.execute(query).fetchdf()
    con.unregister("_wanted_ids")
    return result


def fetch_domain_hits(
    con, table: str, protein_ids: Iterable[str], id_col: str = "protein_id"
) -> pd.DataFrame:
    """Return all rows from a functional-annotation table for the given
    protein ids. No dedup — callers decide how to collapse multi-row-
    per-protein tables (pfam) vs. single-row tables (tmhmm/targetp/signalp).

    `table` is restricted to a fixed allow-list (never interpolated from
    external input) since it's spliced directly into the SQL statement.
    """
    if table not in _ALLOWED_DOMAIN_TABLES:
        raise ValueError(f"table must be one of {_ALLOWED_DOMAIN_TABLES}, got {table!r}")

    ids_df = pd.DataFrame({id_col: list(protein_ids)})
    con.register("_wanted_ids", ids_df)
    query = f"SELECT t.* FROM {table} t JOIN _wanted_ids w ON t.{id_col} = w.{id_col}"  # noqa: S608
    result = con.execute(query).fetchdf()
    con.unregister("_wanted_ids")
    return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/adhesion_properties/test_db.py -v`
Expected: all PASS.

- [ ] **Step 5: Lint**

Run: `ruff check analysis/adhesion_properties/db.py tests/adhesion_properties/`

- [ ] **Step 6: Commit**

```bash
git add analysis/adhesion_properties/db.py tests/adhesion_properties/conftest.py tests/adhesion_properties/test_db.py
git commit -m "Add DuckDB query helpers for adhesion protein properties analysis"
```

---

### Task 2: `properties.py` — sequence property computation

**Files:**
- Create: `analysis/adhesion_properties/properties.py`
- Test: `tests/adhesion_properties/test_properties.py`

**Interfaces:**
- Produces (used by Task 7): `sequence_properties(peptide: str) -> dict` with keys `pct_ser, pct_thr, pct_pro, pct_ser_thr_pro, pct_cys, aromaticity, mean_hydrophobicity, net_charge_ph7`.

- [ ] **Step 1: Write the failing tests**

Create `tests/adhesion_properties/test_properties.py`:

```python
"""Tests for analysis/adhesion_properties/properties.py."""
import math

from properties import sequence_properties


def test_sequence_properties_pure_ser_thr_pro():
    result = sequence_properties("SSSTTTPPP")
    assert math.isclose(result["pct_ser"], 100 / 3, rel_tol=1e-6)
    assert math.isclose(result["pct_thr"], 100 / 3, rel_tol=1e-6)
    assert math.isclose(result["pct_pro"], 100 / 3, rel_tol=1e-6)
    assert math.isclose(result["pct_ser_thr_pro"], 100.0, rel_tol=1e-6)
    assert result["aromaticity"] == 0.0
    assert result["net_charge_ph7"] == 0.0


def test_sequence_properties_aromaticity():
    result = sequence_properties("FWYA")
    assert math.isclose(result["aromaticity"], 75.0, rel_tol=1e-6)


def test_sequence_properties_net_charge():
    result = sequence_properties("KRDE")
    assert result["net_charge_ph7"] == 0.0
    result2 = sequence_properties("KKRR")
    assert result2["net_charge_ph7"] == 4.0


def test_sequence_properties_empty_string_returns_nan():
    result = sequence_properties("")
    assert math.isnan(result["pct_ser"])
    assert math.isnan(result["mean_hydrophobicity"])


def test_sequence_properties_handles_unknown_residue():
    """An unusual/ambiguous residue code (X) should not raise, and should
    still be counted in length for percentage denominators."""
    result = sequence_properties("SSXX")
    assert math.isclose(result["pct_ser"], 50.0, rel_tol=1e-6)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/adhesion_properties/test_properties.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'properties'`.

- [ ] **Step 3: Implement `properties.py`**

Create `analysis/adhesion_properties/properties.py`:

```python
"""Sequence property computation for adhesion-protein analysis."""

from typing import Dict

_KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
_AROMATIC = set("FWY")
_POSITIVE = set("KR")
_NEGATIVE = set("DE")

_NAN_RESULT = {
    "pct_ser": float("nan"), "pct_thr": float("nan"), "pct_pro": float("nan"),
    "pct_ser_thr_pro": float("nan"), "pct_cys": float("nan"),
    "aromaticity": float("nan"), "mean_hydrophobicity": float("nan"),
    "net_charge_ph7": float("nan"),
}


def sequence_properties(peptide: str) -> Dict[str, float]:
    """Compute composition/property features for one peptide sequence.

    All percentage features are in [0, 100]. Unrecognized residue codes
    (X, *, etc.) are counted toward length/percentages but contribute 0 to
    hydrophobicity, not excluded — a coarse approximation, not a crash.
    net_charge_ph7 is a simple Lys+Arg minus Asp+Glu residue count, not a
    full pKa model — stated as an approximation in the report.
    """
    seq = peptide.upper()
    n = len(seq)
    if n == 0:
        return dict(_NAN_RESULT)

    counts: Dict[str, int] = {}
    for aa in seq:
        counts[aa] = counts.get(aa, 0) + 1

    pct_ser = 100.0 * counts.get("S", 0) / n
    pct_thr = 100.0 * counts.get("T", 0) / n
    pct_pro = 100.0 * counts.get("P", 0) / n
    pct_cys = 100.0 * counts.get("C", 0) / n
    aromaticity = 100.0 * sum(counts.get(aa, 0) for aa in _AROMATIC) / n
    hydrophobicity_sum = sum(_KYTE_DOOLITTLE.get(aa, 0.0) * c for aa, c in counts.items())
    mean_hydrophobicity = hydrophobicity_sum / n
    net_charge = sum(counts.get(aa, 0) for aa in _POSITIVE) - sum(
        counts.get(aa, 0) for aa in _NEGATIVE
    )

    return {
        "pct_ser": pct_ser,
        "pct_thr": pct_thr,
        "pct_pro": pct_pro,
        "pct_ser_thr_pro": pct_ser + pct_thr + pct_pro,
        "pct_cys": pct_cys,
        "aromaticity": aromaticity,
        "mean_hydrophobicity": mean_hydrophobicity,
        "net_charge_ph7": float(net_charge),
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/adhesion_properties/test_properties.py -v`
Expected: all PASS.

- [ ] **Step 5: Lint**

Run: `ruff check analysis/adhesion_properties/properties.py tests/adhesion_properties/test_properties.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/adhesion_properties/properties.py tests/adhesion_properties/test_properties.py
git commit -m "Add amino-acid composition property functions"
```

---

### Task 3: `universe.py` — adhesion set + background sample builder

**Files:**
- Create: `analysis/adhesion_properties/universe.py`
- Test: `tests/adhesion_properties/test_universe.py`

**Interfaces:**
- Consumes: `kingdom_survey/join.py`'s `first_id_from_fai`, `locustag_from_id` (imported, not duplicated).
- Produces (used by Task 4):
  - `load_taxonomy_by_locustag(summary_csv_path: Path) -> Dict[str, dict]`
  - `read_result_ids(result_csv_path: Path) -> List[Tuple[str, float]]`
  - `read_fai_ids(fai_path: Path) -> List[str]`
  - `sample_background_ids(all_ids: List[str], adhesion_ids: set, n_wanted: int, rng: random.Random) -> List[str]`
  - `build_protein_universe(results_dir: Path, input_dir: Path, summary_csv_path: Path, seed: int = 42) -> Tuple[List[dict], List[dict]]` (adhesion_rows, background_rows), each row a dict with keys `locustag, phylum, class, order, family, genus, protein_id, probability_adhesion` (background rows have `probability_adhesion=None`)

- [ ] **Step 1: Write the failing tests**

Create `tests/adhesion_properties/test_universe.py`:

```python
"""Tests for analysis/adhesion_properties/universe.py."""
import random

from universe import (
    build_protein_universe,
    load_taxonomy_by_locustag,
    read_fai_ids,
    read_result_ids,
    sample_background_ids,
)

_SUMMARY_HEADER = (
    "locustag,asmid,species_name,phylum,subphylum,class,order,family,genus,"
    "total_proteins,adhesion_count,adhesion_fraction,mean_adhesion_prob,median_adhesion_prob\n"
)


def test_read_result_ids(tmp_path):
    result = tmp_path / "sp.adhesion_predict.csv"
    result.write_text(
        "id,prediction,probability_adhesion\n"
        "LOC1_000001-T1,Adhesion,0.9\n"
        "LOC1_000002-T1,Adhesion,0.7\n"
    )
    assert read_result_ids(result) == [("LOC1_000001-T1", 0.9), ("LOC1_000002-T1", 0.7)]


def test_read_fai_ids(tmp_path):
    fai = tmp_path / "sp.proteins.fa.fai"
    fai.write_text("LOC1_000001-T1\t100\t10\t60\t61\nLOC1_000002-T1\t100\t10\t60\t61\n")
    assert read_fai_ids(fai) == ["LOC1_000001-T1", "LOC1_000002-T1"]


def test_load_taxonomy_by_locustag(tmp_path):
    summary = tmp_path / "species_adhesion_summary.csv"
    summary.write_text(
        _SUMMARY_HEADER
        + "LOC1,ASM1,Sp one,Ascomycota,Pezizomycotina,Sordariomycetes,Hypocreales,"
        "Nectriaceae,Fusarium,100,5,0.05,0.8,0.8\n"
    )
    lookup = load_taxonomy_by_locustag(summary)
    assert lookup["LOC1"]["phylum"] == "Ascomycota"
    assert lookup["LOC1"]["genus"] == "Fusarium"


def test_sample_background_ids_returns_all_when_scarce():
    rng = random.Random(1)
    result = sample_background_ids(["A", "B", "C"], {"A"}, n_wanted=10, rng=rng)
    assert sorted(result) == ["B", "C"]


def test_sample_background_ids_caps_at_n_wanted():
    rng = random.Random(1)
    all_ids = [f"P{i}" for i in range(100)]
    result = sample_background_ids(all_ids, set(), n_wanted=10, rng=rng)
    assert len(result) == 10
    assert set(result).issubset(set(all_ids))


def test_build_protein_universe_matched_species(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_one.proteins.fa.fai").write_text(
        "".join(f"LOC1_{i:06d}-T1\t100\t10\t60\t61\n" for i in range(1, 11))
    )
    (results_dir / "Sp_one.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\n"
        "LOC1_000001-T1,Adhesion,0.9\n"
        "LOC1_000002-T1,Adhesion,0.8\n"
    )
    summary = tmp_path / "species_adhesion_summary.csv"
    summary.write_text(
        _SUMMARY_HEADER
        + "LOC1,ASM1,Sp one,Ascomycota,Pezizomycotina,Sordariomycetes,Hypocreales,"
        "Nectriaceae,Fusarium,10,2,0.2,0.85,0.85\n"
    )

    adhesion_rows, background_rows = build_protein_universe(results_dir, input_dir, summary, seed=1)

    assert len(adhesion_rows) == 2
    assert {r["protein_id"] for r in adhesion_rows} == {"LOC1_000001-T1", "LOC1_000002-T1"}
    assert adhesion_rows[0]["phylum"] == "Ascomycota"
    assert len(background_rows) == 2
    assert all(
        r["protein_id"] not in {"LOC1_000001-T1", "LOC1_000002-T1"} for r in background_rows
    )


def test_build_protein_universe_skips_unmatched_species(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_two.proteins.fa.fai").write_text("LOC2_000001-T1\t100\t10\t60\t61\n")
    (results_dir / "Sp_two.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nLOC2_000001-T1,Adhesion,0.9\n"
    )
    summary = tmp_path / "species_adhesion_summary.csv"
    summary.write_text(_SUMMARY_HEADER)  # no rows — LOC2 not present

    adhesion_rows, background_rows = build_protein_universe(results_dir, input_dir, summary, seed=1)
    assert adhesion_rows == []
    assert background_rows == []
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/adhesion_properties/test_universe.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'universe'`.

- [ ] **Step 3: Implement `universe.py`**

Create `analysis/adhesion_properties/universe.py`:

```python
"""Build the adhesion-protein set and a matched per-species background
sample of non-adhesion proteins, for the protein-properties analysis.

Reuses kingdom_survey's join.py locustag-matching primitives (the .fai-
anchored join key, verified during that project's design review) rather
than re-deriving them.
"""

import csv
import random
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

_KINGDOM_SURVEY_DIR = Path(__file__).resolve().parent.parent / "kingdom_survey"
sys.path.insert(0, str(_KINGDOM_SURVEY_DIR))
from join import first_id_from_fai, locustag_from_id  # noqa: E402

RESULT_SUFFIX = ".adhesion_predict.csv"
FAI_SUFFIX = ".proteins.fa.fai"
_TAXONOMY_FIELDS = ["phylum", "class", "order", "family", "genus"]


def load_taxonomy_by_locustag(summary_csv_path: Path) -> Dict[str, dict]:
    """Load kingdom_survey's species_adhesion_summary.csv keyed by locustag,
    for attaching taxonomy to individual protein rows."""
    lookup = {}
    with open(summary_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lookup[row["locustag"]] = row
    return lookup


def read_result_ids(result_csv_path: Path) -> List[Tuple[str, float]]:
    """Return [(protein_id, probability_adhesion), ...] for one result CSV."""
    ids = []
    with open(result_csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            ids.append((row["id"], float(row["probability_adhesion"])))
    return ids


def read_fai_ids(fai_path: Path) -> List[str]:
    """Return every protein id in a .fai index (first column)."""
    ids = []
    with open(fai_path) as fh:
        for line in fh:
            ids.append(line.split("\t", 1)[0])
    return ids


def sample_background_ids(
    all_ids: List[str], adhesion_ids: Set[str], n_wanted: int, rng: random.Random
) -> List[str]:
    """Sample up to n_wanted protein ids from all_ids that are NOT in
    adhesion_ids. Returns all available non-adhesion ids if fewer than
    n_wanted exist."""
    candidates = [pid for pid in all_ids if pid not in adhesion_ids]
    if len(candidates) <= n_wanted:
        return candidates
    return rng.sample(candidates, n_wanted)


def build_protein_universe(
    results_dir: Path, input_dir: Path, summary_csv_path: Path, seed: int = 42
) -> Tuple[List[dict], List[dict]]:
    """Build (adhesion_rows, background_rows) for every species with a
    result file whose locustag is in the taxonomy lookup. Species not
    found in the lookup (unmatched/mismatched per kingdom_survey) are
    skipped — this analysis only covers species kingdom_survey already
    successfully joined.
    """
    rng = random.Random(seed)
    taxonomy = load_taxonomy_by_locustag(summary_csv_path)

    adhesion_rows: List[dict] = []
    background_rows: List[dict] = []

    for result_path in sorted(results_dir.glob(f"*{RESULT_SUFFIX}")):
        stem = result_path.name[: -len(RESULT_SUFFIX)]
        fai_path = input_dir / f"{stem}{FAI_SUFFIX}"
        if not fai_path.exists():
            continue
        locustag = locustag_from_id(first_id_from_fai(fai_path))
        if locustag not in taxonomy:
            continue
        tax = taxonomy[locustag]
        tax_fields = {field: tax[field] for field in _TAXONOMY_FIELDS}

        result_ids = read_result_ids(result_path)
        adhesion_id_set = {pid for pid, _ in result_ids}
        for pid, prob in result_ids:
            adhesion_rows.append(
                {"locustag": locustag, **tax_fields, "protein_id": pid, "probability_adhesion": prob}
            )

        all_ids = read_fai_ids(fai_path)
        background_ids = sample_background_ids(all_ids, adhesion_id_set, len(result_ids), rng)
        for pid in background_ids:
            background_rows.append(
                {"locustag": locustag, **tax_fields, "protein_id": pid, "probability_adhesion": None}
            )

    return adhesion_rows, background_rows
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/adhesion_properties/test_universe.py -v`
Expected: all PASS.

- [ ] **Step 5: Lint**

Run: `ruff check analysis/adhesion_properties/universe.py tests/adhesion_properties/test_universe.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/adhesion_properties/universe.py tests/adhesion_properties/test_universe.py
git commit -m "Add adhesion-protein universe and background-sample builder"
```

---

### Task 4: `01_build_protein_universe.py` — run on real data

**Files:**
- Create: `analysis/adhesion_properties/01_build_protein_universe.py`
- Create: `analysis/adhesion_properties/tables/.gitkeep`, `analysis/adhesion_properties/figures/.gitkeep`
- Modify: `.gitignore` (add `analysis/adhesion_properties/tables/*.csv`, `analysis/adhesion_properties/figures/*.png`, with `.gitkeep` negations, mirroring the `kingdom_survey` entries)

**Interfaces:**
- Consumes: `universe.py` (Task 3).
- Produces: `analysis/adhesion_properties/tables/protein_universe.csv` with columns `protein_id, group, locustag, phylum, class, order, family, genus, probability_adhesion` — consumed by Tasks 6, 7, 9.

- [ ] **Step 1: Update `.gitignore`**

Append to `.gitignore`:

```
analysis/adhesion_properties/tables/*.csv
analysis/adhesion_properties/figures/*.png
analysis/adhesion_properties/figures/*.svg
!analysis/adhesion_properties/tables/.gitkeep
!analysis/adhesion_properties/figures/.gitkeep
```

- [ ] **Step 2: Create placeholder dirs**

```bash
mkdir -p analysis/adhesion_properties/tables analysis/adhesion_properties/figures
touch analysis/adhesion_properties/tables/.gitkeep analysis/adhesion_properties/figures/.gitkeep
```

- [ ] **Step 3: Write the CLI script**

Create `analysis/adhesion_properties/01_build_protein_universe.py`:

```python
#!/usr/bin/env python
"""Build the adhesion-protein set and matched background sample, writing
tables/protein_universe.csv (one row per protein, group=adhesion|background)."""
import csv
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from universe import build_protein_universe  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results"
INPUT_DIR = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/input")
SUMMARY_CSV = REPO_ROOT / "analysis" / "kingdom_survey" / "tables" / "species_adhesion_summary.csv"
OUT_DIR = Path(__file__).resolve().parent / "tables"

FIELDS = [
    "protein_id", "group", "locustag", "phylum", "class", "order", "family",
    "genus", "probability_adhesion",
]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adhesion_rows, background_rows = build_protein_universe(
        RESULTS_DIR, INPUT_DIR, SUMMARY_CSV, seed=42
    )

    with open(OUT_DIR / "protein_universe.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        for row in adhesion_rows:
            writer.writerow({**row, "group": "adhesion"})
        for row in background_rows:
            writer.writerow({**row, "group": "background", "probability_adhesion": ""})

    print(f"Adhesion proteins: {len(adhesion_rows)}")
    print(f"Background proteins: {len(background_rows)}")

    summary = pd.read_csv(SUMMARY_CSV)
    expected = int(summary["adhesion_count"].sum())
    assert len(adhesion_rows) == expected, (
        f"Reconciliation failed: {len(adhesion_rows)} adhesion rows built, "
        f"expected {expected} from species_adhesion_summary.csv"
    )


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run it on real data**

Run: `python analysis/adhesion_properties/01_build_protein_universe.py`

This reads every matched species' result file and `.fai` (thousands of small files, similar to `kingdom_survey`'s Task 3 — expect it to take a few minutes). Expected output: `Adhesion proteins: <N>` matching `species_adhesion_summary.csv`'s summed `adhesion_count` (no `AssertionError`), and `Background proteins:` roughly similar in magnitude.

- [ ] **Step 5: Lint**

Run: `ruff check analysis/adhesion_properties/01_build_protein_universe.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/adhesion_properties/01_build_protein_universe.py analysis/adhesion_properties/tables/.gitkeep analysis/adhesion_properties/figures/.gitkeep .gitignore
git commit -m "Add CLI to build the adhesion protein universe from real data"
```

---

### Task 5: `domains.py` — domain-flag and top-domain-enrichment functions

**Files:**
- Create: `analysis/adhesion_properties/domains.py`
- Test: `tests/adhesion_properties/test_domains.py`

**Interfaces:**
- Consumes: `db.fetch_domain_hits` (Task 1).
- Produces (used by Task 6):
  - `clean_cazy_family(cazyme_fam: str) -> str`
  - `compute_domain_flags(con, protein_ids: Iterable[str]) -> pd.DataFrame` with columns `protein_id, has_pfam, has_cazy, has_merops, has_signal_peptide, has_tm_helix, targetp_category`
  - `top_domain_table(con, table: str, id_field: str, adhesion_ids: Iterable[str], background_ids: Iterable[str], top_n: int = 20) -> pd.DataFrame` with columns `domain_id, adhesion_count, adhesion_rate, background_count, background_rate, enrichment_ratio`

- [ ] **Step 1: Write the failing tests**

Create `tests/adhesion_properties/test_domains.py`:

```python
"""Tests for analysis/adhesion_properties/domains.py — uses an in-memory
DuckDB with tiny fake tables, never the real 278GB database."""
import duckdb

from domains import clean_cazy_family, compute_domain_flags, top_domain_table


def _fake_connection():
    con = duckdb.connect(":memory:")
    con.execute("CREATE TABLE pfam (protein_id VARCHAR, pfam_id VARCHAR)")
    con.execute(
        "INSERT INTO pfam VALUES "
        "('A1', 'DomainX'), ('A1', 'DomainX'), ('A2', 'DomainY'), ('B1', 'DomainX')"
    )
    con.execute("CREATE TABLE cazy_overview (protein_id VARCHAR, cazyme_fam VARCHAR)")
    con.execute("INSERT INTO cazy_overview VALUES ('A1', 'GH114(59-261)')")
    con.execute("CREATE TABLE merops (protein_id VARCHAR, merops_id VARCHAR)")
    con.execute("INSERT INTO merops VALUES ('A2', 'MER0000001')")
    con.execute("CREATE TABLE signalp (protein_id VARCHAR, probability DOUBLE)")
    con.execute("INSERT INTO signalp VALUES ('A1', 0.99)")
    con.execute("CREATE TABLE tmhmm (protein_id VARCHAR, PredHel BIGINT)")
    con.execute("INSERT INTO tmhmm VALUES ('A1', 0), ('A2', 2)")
    con.execute("CREATE TABLE targetp (protein_id VARCHAR, prediction VARCHAR)")
    con.execute("INSERT INTO targetp VALUES ('A1', 'SP')")
    return con


def test_clean_cazy_family_strips_range_suffix():
    assert clean_cazy_family("GH114(59-261)") == "GH114"
    assert clean_cazy_family("GT2_Chitin_synth_2(661-1184)") == "GT2_Chitin_synth_2"
    assert clean_cazy_family("GH20") == "GH20"


def test_compute_domain_flags():
    con = _fake_connection()
    flags = compute_domain_flags(con, ["A1", "A2", "B2"]).set_index("protein_id")
    assert bool(flags.loc["A1", "has_pfam"]) is True
    assert bool(flags.loc["A1", "has_cazy"]) is True
    assert bool(flags.loc["A1", "has_merops"]) is False
    assert bool(flags.loc["A1", "has_signal_peptide"]) is True
    assert bool(flags.loc["A1", "has_tm_helix"]) is False  # PredHel=0
    assert flags.loc["A1", "targetp_category"] == "SP"
    assert bool(flags.loc["A2", "has_tm_helix"]) is True  # PredHel=2
    assert flags.loc["A2", "targetp_category"] == "none"
    assert bool(flags.loc["B2", "has_pfam"]) is False


def test_top_domain_table_counts_distinct_proteins_not_raw_rows():
    con = _fake_connection()
    # A1 has 2 pfam rows for DomainX but should count once as a protein
    top = top_domain_table(
        con, "pfam", "pfam_id", adhesion_ids=["A1", "A2"], background_ids=["B1"], top_n=5
    )
    domain_x = top.set_index("domain_id").loc["DomainX"]
    assert domain_x["adhesion_count"] == 1  # A1 only, not counted twice
    assert domain_x["background_count"] == 1  # B1


def test_top_domain_table_cleans_cazy_family_names():
    con = _fake_connection()
    top = top_domain_table(
        con, "cazy_overview", "cazyme_fam", adhesion_ids=["A1"], background_ids=[], top_n=5
    )
    assert "GH114" in set(top["domain_id"])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/adhesion_properties/test_domains.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'domains'`.

- [ ] **Step 3: Implement `domains.py`**

Create `analysis/adhesion_properties/domains.py`:

```python
"""Functional domain flag computation for the protein-properties analysis."""

import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from db import fetch_domain_hits  # noqa: E402

_CAZY_FAM_RE = re.compile(r"^([^(]+)")
MIN_ADHESION_COUNT_FOR_ENRICHMENT = 5


def clean_cazy_family(cazyme_fam: str) -> str:
    """Strip the alignment-range suffix from a CAZy family string,
    e.g. 'GH114(59-261)' -> 'GH114'."""
    match = _CAZY_FAM_RE.match(cazyme_fam)
    return match.group(1) if match else cazyme_fam


def compute_domain_flags(con, protein_ids: Iterable[str]) -> pd.DataFrame:
    """Return one row per protein_id with boolean/categorical domain flags."""
    protein_ids = list(protein_ids)
    base = pd.DataFrame({"protein_id": protein_ids})

    pfam_hits = set(fetch_domain_hits(con, "pfam", protein_ids)["protein_id"])
    cazy_hits = set(fetch_domain_hits(con, "cazy_overview", protein_ids)["protein_id"])
    merops_hits = set(fetch_domain_hits(con, "merops", protein_ids)["protein_id"])
    signalp_hits = set(fetch_domain_hits(con, "signalp", protein_ids)["protein_id"])

    tmhmm = fetch_domain_hits(con, "tmhmm", protein_ids)
    tm_hits = set(tmhmm.loc[tmhmm["PredHel"] >= 1, "protein_id"]) if not tmhmm.empty else set()

    targetp = fetch_domain_hits(con, "targetp", protein_ids)
    targetp_map = dict(zip(targetp["protein_id"], targetp["prediction"])) if not targetp.empty else {}

    base["has_pfam"] = base["protein_id"].isin(pfam_hits)
    base["has_cazy"] = base["protein_id"].isin(cazy_hits)
    base["has_merops"] = base["protein_id"].isin(merops_hits)
    base["has_signal_peptide"] = base["protein_id"].isin(signalp_hits)
    base["has_tm_helix"] = base["protein_id"].isin(tm_hits)
    base["targetp_category"] = base["protein_id"].map(targetp_map).fillna("none")

    return base


def top_domain_table(
    con,
    table: str,
    id_field: str,
    adhesion_ids: Iterable[str],
    background_ids: Iterable[str],
    top_n: int = 20,
) -> pd.DataFrame:
    """Distinct-protein frequency of each domain identity among adhesion vs.
    background proteins, restricted to domains with at least
    MIN_ADHESION_COUNT_FOR_ENRICHMENT adhesion hits (avoids a single rare
    hit with zero background hits producing a meaningless infinite ratio),
    sorted by enrichment ratio (adhesion_rate / background_rate).
    """
    adhesion_ids = list(adhesion_ids)
    background_ids = list(background_ids)
    n_adhesion = len(adhesion_ids)
    n_background = len(background_ids)

    adhesion_hits = fetch_domain_hits(con, table, adhesion_ids)
    background_hits = fetch_domain_hits(con, table, background_ids)

    if table == "cazy_overview":
        if not adhesion_hits.empty:
            adhesion_hits = adhesion_hits.assign(
                **{id_field: adhesion_hits[id_field].map(clean_cazy_family)}
            )
        if not background_hits.empty:
            background_hits = background_hits.assign(
                **{id_field: background_hits[id_field].map(clean_cazy_family)}
            )

    def _distinct_counts(hits: pd.DataFrame) -> pd.Series:
        if hits.empty:
            return pd.Series(dtype=int)
        return hits.drop_duplicates(["protein_id", id_field]).groupby(id_field)["protein_id"].nunique()

    adhesion_counts = _distinct_counts(adhesion_hits)
    background_counts = _distinct_counts(background_hits)

    all_ids = sorted(set(adhesion_counts.index) | set(background_counts.index))
    rows = []
    for domain_id in all_ids:
        a_count = int(adhesion_counts.get(domain_id, 0))
        b_count = int(background_counts.get(domain_id, 0))
        a_rate = a_count / n_adhesion if n_adhesion else 0.0
        b_rate = b_count / n_background if n_background else 0.0
        enrichment = (a_rate / b_rate) if b_rate > 0 else float("inf")
        rows.append(
            {
                "domain_id": domain_id,
                "adhesion_count": a_count,
                "adhesion_rate": a_rate,
                "background_count": b_count,
                "background_rate": b_rate,
                "enrichment_ratio": enrichment,
            }
        )

    result = pd.DataFrame(rows)
    if result.empty:
        return result
    frequent = result[result["adhesion_count"] >= MIN_ADHESION_COUNT_FOR_ENRICHMENT]
    return frequent.sort_values("enrichment_ratio", ascending=False).head(top_n)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/adhesion_properties/test_domains.py -v`
Expected: all PASS. (`test_top_domain_table_counts_distinct_proteins_not_raw_rows` and `test_top_domain_table_cleans_cazy_family_names` use fewer than `MIN_ADHESION_COUNT_FOR_ENRICHMENT` hits per domain in their fixtures — if the frequency filter empties the result, lower the module constant is NOT the fix; instead the test fixtures must produce at least `MIN_ADHESION_COUNT_FOR_ENRICHMENT` adhesion hits for the domain under test, or call the underlying counting logic directly. Adjust the fixture's row count for the specific domain being asserted on, not the production threshold.)

- [ ] **Step 5: Lint**

Run: `ruff check analysis/adhesion_properties/domains.py tests/adhesion_properties/test_domains.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/adhesion_properties/domains.py tests/adhesion_properties/test_domains.py
git commit -m "Add functional domain flag and top-domain enrichment functions"
```

---

### Task 6: `02_query_functional_domains.py` — run on real data

**Files:**
- Create: `analysis/adhesion_properties/02_query_functional_domains.py`

**Interfaces:**
- Consumes: `db.py` (Task 1), `domains.py` (Task 5), `tables/protein_universe.csv` (Task 4).
- Produces: `tables/protein_domain_flags.csv`, `tables/top_domains_pfam.csv`, `tables/top_domains_cazy.csv`, `tables/top_domains_merops.csv` — consumed by Tasks 9, 11, 12.

- [ ] **Step 1: Write the CLI script**

Create `analysis/adhesion_properties/02_query_functional_domains.py`:

```python
#!/usr/bin/env python
"""Query functional domain annotations for the adhesion/background protein
universe, writing per-protein domain flags and top-domain enrichment tables."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from db import connect  # noqa: E402
from domains import compute_domain_flags, top_domain_table  # noqa: E402

TABLES_DIR = Path(__file__).resolve().parent / "tables"

DOMAIN_TABLES = [
    ("pfam", "pfam_id", "top_domains_pfam.csv"),
    ("cazy_overview", "cazyme_fam", "top_domains_cazy.csv"),
    ("merops", "merops_id", "top_domains_merops.csv"),
]


def main() -> None:
    universe = pd.read_csv(TABLES_DIR / "protein_universe.csv")
    adhesion_ids = universe.loc[universe["group"] == "adhesion", "protein_id"].tolist()
    background_ids = universe.loc[universe["group"] == "background", "protein_id"].tolist()

    con = connect()

    flags = compute_domain_flags(con, universe["protein_id"].tolist())
    flags.to_csv(TABLES_DIR / "protein_domain_flags.csv", index=False)
    print(f"Domain flags computed for {len(flags)} proteins")

    for table, id_field, out_name in DOMAIN_TABLES:
        top = top_domain_table(con, table, id_field, adhesion_ids, background_ids)
        top.to_csv(TABLES_DIR / out_name, index=False)
        print(f"{out_name}: {len(top)} rows")

    con.close()


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python analysis/adhesion_properties/02_query_functional_domains.py`

This queries the 278GB DuckDB filtered against ~1.5M protein ids (adhesion + background combined) across 6 tables — allow it time to complete. Expected: 4 CSVs written, non-empty flag table, each top-domain table with up to 20 rows.

- [ ] **Step 3: Spot-check against a known protein**

Run:
```bash
python3 -c "
import pandas as pd
flags = pd.read_csv('analysis/adhesion_properties/tables/protein_domain_flags.csv')
print(flags[flags['protein_id'] == 'F07B100A_000481-T1'])
"
```
This protein was confirmed during design review to have no Pfam/CAZy hit — expect `has_pfam` and `has_cazy` both `False` for it (if it appears in the flags table at all; it may or may not be in the adhesion/background universe depending on sampling — if absent, spot-check a different protein_id actually present in `protein_universe.csv` instead).

- [ ] **Step 4: Lint**

Run: `ruff check analysis/adhesion_properties/02_query_functional_domains.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/adhesion_properties/02_query_functional_domains.py
git commit -m "Add CLI to query functional domain annotations on real data"
```

---

### Task 7: `03_compute_sequence_properties.py` — run on real data

**Files:**
- Create: `analysis/adhesion_properties/03_compute_sequence_properties.py`

**Interfaces:**
- Consumes: `db.py` (Task 1), `properties.py` (Task 2), `tables/protein_universe.csv` (Task 4).
- Produces: `tables/protein_sequence_properties.csv` (columns: everything from `protein_universe.csv` plus `length, pct_ser, pct_thr, pct_pro, pct_ser_thr_pro, pct_cys, aromaticity, mean_hydrophobicity, net_charge_ph7`) — consumed by Task 9.

- [ ] **Step 1: Write the CLI script**

Create `analysis/adhesion_properties/03_compute_sequence_properties.py`:

```python
#!/usr/bin/env python
"""Compute length and amino-acid-composition properties for every protein
in the adhesion/background universe."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from db import connect, fetch_lengths  # noqa: E402
from properties import sequence_properties  # noqa: E402

TABLES_DIR = Path(__file__).resolve().parent / "tables"


def main() -> None:
    universe = pd.read_csv(TABLES_DIR / "protein_universe.csv")
    con = connect()
    seq_data = fetch_lengths(con, universe["protein_id"].tolist())
    con.close()

    prop_rows = []
    for _, row in seq_data.iterrows():
        props = sequence_properties(row["peptide"] or "")
        prop_rows.append({"protein_id": row["protein_id"], "length": row["length"], **props})

    props_df = pd.DataFrame(prop_rows)
    merged = universe.merge(props_df, on="protein_id", how="inner")
    merged.to_csv(TABLES_DIR / "protein_sequence_properties.csv", index=False)

    n_missing = len(universe) - len(merged)
    print(f"Computed properties for {len(merged)} proteins ({n_missing} not found in gene_proteins)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python analysis/adhesion_properties/03_compute_sequence_properties.py`

This pulls peptide sequences for ~1.5M proteins from `gene_proteins` — allow it time. If `n_missing` is large (more than a handful), investigate before proceeding — it likely means the `.protein` suffix join in `db.fetch_lengths` isn't matching as expected at real-data scale, not that the data is genuinely absent (design review confirmed the convention holds for at least one real ID).

- [ ] **Step 3: Sanity-check bounds**

Run:
```bash
python3 -c "
import pandas as pd
df = pd.read_csv('analysis/adhesion_properties/tables/protein_sequence_properties.csv')
assert df['pct_ser_thr_pro'].dropna().between(0, 100).all()
assert (df['length'] > 0).all()
print('bounds OK, n =', len(df))
"
```

- [ ] **Step 4: Lint**

Run: `ruff check analysis/adhesion_properties/03_compute_sequence_properties.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/adhesion_properties/03_compute_sequence_properties.py
git commit -m "Add CLI to compute sequence properties on real data"
```

---

### Task 8: `stats_ext.py` — per-clade aggregation and Mann-Whitney comparison

**Files:**
- Create: `analysis/adhesion_properties/stats_ext.py`
- Test: `tests/adhesion_properties/test_stats_ext.py`

**Interfaces:**
- Consumes: `kingdom_survey/stats.py`'s `MIN_GROUP_N`, `kruskal_wallis_by_rank` (imported, not duplicated).
- Produces (used by Task 9):
  - `aggregate_by_clade(df: pd.DataFrame, rank_col: str, value_col: str) -> pd.DataFrame` with columns `rank_col, n, median, iqr, small_n`
  - `mannwhitney_within_clade(adhesion_df, background_df, rank_col: str, value_col: str) -> pd.DataFrame` with columns `rank_col, n_adhesion, n_background, median_adhesion, median_background, u_stat, p_value, p_value_bh`
  - re-exports `MIN_GROUP_N`, `kruskal_wallis_by_rank` from `kingdom_survey.stats`

- [ ] **Step 1: Write the failing tests**

Create `tests/adhesion_properties/test_stats_ext.py`:

```python
"""Tests for analysis/adhesion_properties/stats_ext.py."""
import numpy as np
import pandas as pd

from stats_ext import MIN_GROUP_N, aggregate_by_clade, mannwhitney_within_clade


def _clade_df(clade_values, feature_values):
    return pd.DataFrame({"phylum": clade_values, "length": feature_values})


def test_aggregate_by_clade_basic():
    df = _clade_df(
        ["PhylumX"] * 6 + ["PhylumY"] * 6,
        [100, 110, 120, 130, 140, 150, 200, 210, 220, 230, 240, 250],
    )
    summary = aggregate_by_clade(df, "phylum", "length")
    assert set(summary["phylum"]) == {"PhylumX", "PhylumY"}
    y_median = summary.loc[summary["phylum"] == "PhylumY", "median"].iloc[0]
    x_median = summary.loc[summary["phylum"] == "PhylumX", "median"].iloc[0]
    assert y_median > x_median
    assert not summary["small_n"].any()  # both groups have 6 >= MIN_GROUP_N


def test_aggregate_by_clade_flags_small_n():
    df = _clade_df(["PhylumZ"] * 2, [100, 110])
    summary = aggregate_by_clade(df, "phylum", "length")
    assert summary.iloc[0]["small_n"] is np.True_ or summary.iloc[0]["small_n"] is True
    assert 2 < MIN_GROUP_N


def test_mannwhitney_within_clade_detects_difference():
    rng = np.random.default_rng(1)
    adhesion = _clade_df(["PhylumX"] * 10, list(rng.normal(200, 5, 10)))
    background = _clade_df(["PhylumX"] * 10, list(rng.normal(100, 5, 10)))
    result = mannwhitney_within_clade(adhesion, background, "phylum", "length")
    assert len(result) == 1
    assert result.iloc[0]["p_value"] < 0.05
    assert "p_value_bh" in result.columns


def test_mannwhitney_within_clade_excludes_small_n_clades():
    adhesion = _clade_df(["PhylumSmall"] * 2, [100, 110])
    background = _clade_df(["PhylumSmall"] * 10, list(range(100, 110)))
    result = mannwhitney_within_clade(adhesion, background, "phylum", "length")
    assert result.empty
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/adhesion_properties/test_stats_ext.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'stats_ext'`.

- [ ] **Step 3: Implement `stats_ext.py`**

Create `analysis/adhesion_properties/stats_ext.py`:

```python
"""Per-clade statistics for the adhesion-protein properties analysis.

Reuses kingdom_survey's MIN_GROUP_N/kruskal_wallis_by_rank directly
(already generic over an arbitrary value_col) rather than duplicating
them; adds a generic per-clade summary (kingdom_survey's aggregate_summary
is hardcoded to its own adhesion_fraction/adhesion_count column names, so
is not reusable here) and a within-clade adhesion-vs-background test,
which kingdom_survey has no equivalent of.
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

_KINGDOM_SURVEY_DIR = Path(__file__).resolve().parent.parent / "kingdom_survey"
sys.path.insert(0, str(_KINGDOM_SURVEY_DIR))
from stats import MIN_GROUP_N, kruskal_wallis_by_rank  # noqa: E402,F401


def aggregate_by_clade(df: pd.DataFrame, rank_col: str, value_col: str) -> pd.DataFrame:
    """Per-clade N, median, IQR of an arbitrary value_col."""
    sub = df[df[rank_col].notna()]
    rows = []
    for group, gdf in sub.groupby(rank_col):
        rows.append(
            {
                rank_col: group,
                "n": len(gdf),
                "median": gdf[value_col].median(),
                "iqr": gdf[value_col].quantile(0.75) - gdf[value_col].quantile(0.25),
                "small_n": len(gdf) < MIN_GROUP_N,
            }
        )
    return pd.DataFrame(rows).sort_values("median", ascending=False).reset_index(drop=True)


def mannwhitney_within_clade(
    adhesion_df: pd.DataFrame, background_df: pd.DataFrame, rank_col: str, value_col: str
) -> pd.DataFrame:
    """For each clade with both groups meeting MIN_GROUP_N, Mann-Whitney U
    comparing adhesion vs. background value_col; BH-corrected across clades."""
    a_counts = adhesion_df[rank_col].value_counts()
    b_counts = background_df[rank_col].value_counts()
    eligible = sorted(
        set(a_counts[a_counts >= MIN_GROUP_N].index) & set(b_counts[b_counts >= MIN_GROUP_N].index)
    )

    rows = []
    for group in eligible:
        a_vals = adhesion_df.loc[adhesion_df[rank_col] == group, value_col].dropna()
        b_vals = background_df.loc[background_df[rank_col] == group, value_col].dropna()
        if len(a_vals) < 2 or len(b_vals) < 2:
            continue
        u_stat, p_value = mannwhitneyu(a_vals, b_vals, alternative="two-sided")
        rows.append(
            {
                rank_col: group,
                "n_adhesion": len(a_vals),
                "n_background": len(b_vals),
                "median_adhesion": a_vals.median(),
                "median_background": b_vals.median(),
                "u_stat": u_stat,
                "p_value": p_value,
            }
        )

    result = pd.DataFrame(rows)
    if not result.empty:
        _, corrected, _, _ = multipletests(result["p_value"], method="fdr_bh")
        result["p_value_bh"] = corrected
    return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/adhesion_properties/test_stats_ext.py -v`
Expected: all PASS.

- [ ] **Step 5: Lint**

Run: `ruff check analysis/adhesion_properties/stats_ext.py tests/adhesion_properties/test_stats_ext.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/adhesion_properties/stats_ext.py tests/adhesion_properties/test_stats_ext.py
git commit -m "Add per-clade aggregation and Mann-Whitney comparison functions"
```

---

### Task 9: `04_aggregate_and_test.py` — run on real data

**Files:**
- Create: `analysis/adhesion_properties/04_aggregate_and_test.py`

**Interfaces:**
- Consumes: `stats_ext.py` (Task 8), `tables/protein_sequence_properties.csv` (Task 7), `tables/protein_domain_flags.csv` (Task 6).
- Produces: `tables/adhesion_protein_properties.csv` (master merged table), `tables/clade_{rank}_{feature}_{adhesion|background}.csv`, `tables/clade_{rank}_{feature}_omnibus_{adhesion|background}.csv`, `tables/clade_{rank}_{feature}_mannwhitney.csv` for each of 3 ranks × 5 continuous features, `tables/clade_{rank}_domain_fractions.csv` for each of 3 ranks — consumed by Tasks 11, 12.

- [ ] **Step 1: Write the CLI script**

Create `analysis/adhesion_properties/04_aggregate_and_test.py`:

```python
#!/usr/bin/env python
"""Merge sequence properties + domain flags into one master table, and
run per-clade aggregation/statistics on adhesion vs. background proteins."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats_ext import aggregate_by_clade, kruskal_wallis_by_rank, mannwhitney_within_clade  # noqa: E402

TABLES_DIR = Path(__file__).resolve().parent / "tables"
RANKS = ["phylum", "class", "order"]
CONTINUOUS_FEATURES = [
    "length", "pct_ser_thr_pro", "aromaticity", "mean_hydrophobicity", "net_charge_ph7",
]
DOMAIN_FLAGS = ["has_pfam", "has_cazy", "has_merops", "has_signal_peptide", "has_tm_helix"]


def main() -> None:
    props = pd.read_csv(TABLES_DIR / "protein_sequence_properties.csv")
    flags = pd.read_csv(TABLES_DIR / "protein_domain_flags.csv")
    master = props.merge(flags, on="protein_id", how="inner")
    master.to_csv(TABLES_DIR / "adhesion_protein_properties.csv", index=False)
    print(f"Master table: {len(master)} proteins")

    adhesion = master[master["group"] == "adhesion"]
    background = master[master["group"] == "background"]

    for rank in RANKS:
        for feature in CONTINUOUS_FEATURES:
            aggregate_by_clade(adhesion, rank, feature).to_csv(
                TABLES_DIR / f"clade_{rank}_{feature}_adhesion.csv", index=False
            )
            aggregate_by_clade(background, rank, feature).to_csv(
                TABLES_DIR / f"clade_{rank}_{feature}_background.csv", index=False
            )
            omnibus_adhesion = kruskal_wallis_by_rank(adhesion, rank, feature)
            omnibus_background = kruskal_wallis_by_rank(background, rank, feature)
            pd.DataFrame([omnibus_adhesion]).to_csv(
                TABLES_DIR / f"clade_{rank}_{feature}_omnibus_adhesion.csv", index=False
            )
            pd.DataFrame([omnibus_background]).to_csv(
                TABLES_DIR / f"clade_{rank}_{feature}_omnibus_background.csv", index=False
            )
            mw = mannwhitney_within_clade(adhesion, background, rank, feature)
            mw.to_csv(TABLES_DIR / f"clade_{rank}_{feature}_mannwhitney.csv", index=False)

        rows = []
        for group_name, gdf in [("adhesion", adhesion), ("background", background)]:
            sub = gdf[gdf[rank].notna()]
            for clade, cdf in sub.groupby(rank):
                row = {rank: clade, "group": group_name, "n": len(cdf)}
                for flag in DOMAIN_FLAGS:
                    row[flag] = cdf[flag].mean()
                rows.append(row)
        pd.DataFrame(rows).to_csv(TABLES_DIR / f"clade_{rank}_domain_fractions.csv", index=False)

    print("Aggregation complete.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python analysis/adhesion_properties/04_aggregate_and_test.py`

Expected: `Master table: <N> proteins` (should equal or be close to the sum of Task 7's output), then `Aggregation complete.` with no exceptions.

- [ ] **Step 3: Sanity-check output**

Run:
```bash
python3 -c "
import pandas as pd
mw = pd.read_csv('analysis/adhesion_properties/tables/clade_phylum_pct_ser_thr_pro_mannwhitney.csv')
print(mw)
"
```
Confirm it has rows with sane `p_value`/`p_value_bh` values in [0, 1].

- [ ] **Step 4: Lint**

Run: `ruff check analysis/adhesion_properties/04_aggregate_and_test.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/adhesion_properties/04_aggregate_and_test.py
git commit -m "Add CLI to aggregate and test protein properties by clade on real data"
```

---

### Task 10: `figures_ext.py` — plotting functions

Before writing this task's code, invoke the **dataviz** skill for palette/consistency guidance, same as `kingdom_survey`'s figures task — reuse `kingdom_survey/figures.py`'s already-validated `_FULL_COLOR`/`_OTHER_COLOR` constants (imported, not re-derived) so both reports share one visual language: "adhesion" = full color, "background" = muted gray, matching the semantic kingdom_survey already established for "the reference/baseline group."

**Files:**
- Create: `analysis/adhesion_properties/figures_ext.py`
- Test: `tests/adhesion_properties/test_figures_ext.py`

**Interfaces:**
- Consumes: `kingdom_survey/figures.py`'s `_FULL_COLOR`, `_OTHER_COLOR` (imported, not duplicated).
- Produces (used by Task 11):
  - `paired_boxplot_by_rank(master_df, rank_col: str, value_col: str, out_path: Path, top_n: Optional[int] = None) -> None`
  - `domain_fraction_barplot(fractions_df, rank_col: str, flag_col: str, out_path: Path) -> None`
  - `top_domain_enrichment_barplot(top_df, out_path: Path, title: str, ratio_cap: float = 50.0) -> None`
  - `scatter_feature_vs_probability(adhesion_df, feature_col: str, out_path: Path) -> None`
  - `scatter_species_mean_length_vs_fraction(adhesion_df, kingdom_summary_df, out_path: Path) -> None`

- [ ] **Step 1: Write the failing tests**

Create `tests/adhesion_properties/test_figures_ext.py`:

```python
"""Smoke tests for analysis/adhesion_properties/figures_ext.py — verifies
each function runs without error and writes a non-empty file."""
import pandas as pd

from figures_ext import (
    domain_fraction_barplot,
    paired_boxplot_by_rank,
    scatter_feature_vs_probability,
    scatter_species_mean_length_vs_fraction,
    top_domain_enrichment_barplot,
)


def _master_df():
    rows = []
    for i in range(12):
        rows.append(
            {
                "phylum": "PhylumX" if i < 6 else "PhylumY",
                "group": "adhesion" if i % 2 == 0 else "background",
                "length": 200 + 10 * i,
                "probability_adhesion": 0.5 + 0.02 * i,
                "locustag": f"LOC{i % 3}",
            }
        )
    return pd.DataFrame(rows)


def test_paired_boxplot_by_rank_creates_file(tmp_path):
    df = _master_df()
    out = tmp_path / "box.png"
    paired_boxplot_by_rank(df, "phylum", "length", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_domain_fraction_barplot_creates_file(tmp_path):
    fractions = pd.DataFrame(
        [
            {"phylum": "PhylumX", "group": "adhesion", "n": 6, "has_pfam": 0.5},
            {"phylum": "PhylumX", "group": "background", "n": 6, "has_pfam": 0.3},
            {"phylum": "PhylumY", "group": "adhesion", "n": 6, "has_pfam": 0.6},
            {"phylum": "PhylumY", "group": "background", "n": 6, "has_pfam": 0.4},
        ]
    )
    out = tmp_path / "bar.png"
    domain_fraction_barplot(fractions, "phylum", "has_pfam", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_top_domain_enrichment_barplot_creates_file(tmp_path):
    top = pd.DataFrame(
        [
            {"domain_id": "DomainA", "enrichment_ratio": 3.5},
            {"domain_id": "DomainB", "enrichment_ratio": float("inf")},
        ]
    )
    out = tmp_path / "top.png"
    top_domain_enrichment_barplot(top, out, title="Test")
    assert out.exists()
    assert out.stat().st_size > 0


def test_scatter_feature_vs_probability_creates_file(tmp_path):
    df = _master_df()
    adhesion = df[df["group"] == "adhesion"]
    out = tmp_path / "scatter.png"
    scatter_feature_vs_probability(adhesion, "length", out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_scatter_species_mean_length_vs_fraction_creates_file(tmp_path):
    df = _master_df()
    adhesion = df[df["group"] == "adhesion"]
    kingdom_summary = pd.DataFrame(
        {"locustag": ["LOC0", "LOC1", "LOC2"], "adhesion_fraction": [0.01, 0.02, 0.03]}
    )
    out = tmp_path / "species_scatter.png"
    scatter_species_mean_length_vs_fraction(adhesion, kingdom_summary, out)
    assert out.exists()
    assert out.stat().st_size > 0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/adhesion_properties/test_figures_ext.py -v`
Expected: FAIL/ERROR — `ModuleNotFoundError: No module named 'figures_ext'`.

- [ ] **Step 3: Implement `figures_ext.py`**

Create `analysis/adhesion_properties/figures_ext.py`:

```python
"""Figure generation for the adhesion-protein properties analysis.

Reuses kingdom_survey/figures.py's validated color constants (dataviz
skill palette) for consistency across both reports: "adhesion" = full
color, "background" = the same muted gray kingdom_survey uses for its
de-emphasized/reference states.
"""

import sys
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

_KINGDOM_SURVEY_DIR = Path(__file__).resolve().parent.parent / "kingdom_survey"
sys.path.insert(0, str(_KINGDOM_SURVEY_DIR))
from figures import _FULL_COLOR, _OTHER_COLOR  # noqa: E402

_GROUP_PALETTE = {"adhesion": _FULL_COLOR, "background": _OTHER_COLOR}


def paired_boxplot_by_rank(
    master_df: pd.DataFrame,
    rank_col: str,
    value_col: str,
    out_path: Path,
    top_n: Optional[int] = None,
) -> None:
    """Boxplot of value_col by rank_col, adhesion vs. background as paired
    boxes per clade (hue='group'), ordered by adhesion-group median."""
    sub = master_df[master_df[rank_col].notna()].copy()
    adhesion_medians = (
        sub[sub["group"] == "adhesion"]
        .groupby(rank_col)[value_col]
        .median()
        .sort_values(ascending=False)
    )
    if top_n is not None:
        counts = sub[sub["group"] == "adhesion"][rank_col].value_counts()
        top_groups = counts.sort_values(ascending=False).head(top_n).index
        adhesion_medians = adhesion_medians[adhesion_medians.index.isin(top_groups)]
    order = list(adhesion_medians.index)
    sub = sub[sub[rank_col].isin(order)]

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.6), 6))
    sns.boxplot(
        data=sub, x=rank_col, y=value_col, hue="group", order=order,
        palette=_GROUP_PALETTE, ax=ax,
    )
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90 if len(order) > 8 else 45, ha="right")
    ax.set_xlabel(rank_col.capitalize())
    ax.set_ylabel(value_col.replace("_", " "))
    ax.legend(title="")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def domain_fraction_barplot(
    fractions_df: pd.DataFrame, rank_col: str, flag_col: str, out_path: Path
) -> None:
    """Grouped bar chart of domain-presence fraction by clade, adhesion vs
    background."""
    order = (
        fractions_df[fractions_df["group"] == "adhesion"]
        .sort_values(flag_col, ascending=False)[rank_col]
        .tolist()
    )
    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.6), 6))
    sns.barplot(
        data=fractions_df, x=rank_col, y=flag_col, hue="group", order=order,
        palette=_GROUP_PALETTE, ax=ax,
    )
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90 if len(order) > 8 else 45, ha="right")
    ax.set_xlabel(rank_col.capitalize())
    ax.set_ylabel(f"Fraction with {flag_col.replace('has_', '').replace('_', ' ')}")
    ax.legend(title="")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def top_domain_enrichment_barplot(
    top_df: pd.DataFrame, out_path: Path, title: str, ratio_cap: float = 50.0
) -> None:
    """Horizontal bar chart of enrichment_ratio for the top domains,
    capping infinite/very large ratios at ratio_cap for display (noted in
    the axis label, not silently)."""
    sub = top_df.copy()
    sub["display_ratio"] = sub["enrichment_ratio"].clip(upper=ratio_cap)
    sub = sub.sort_values("display_ratio")

    fig, ax = plt.subplots(figsize=(8, max(4, len(sub) * 0.35)))
    ax.barh(sub["domain_id"], sub["display_ratio"], color=_FULL_COLOR)
    ax.set_xlabel(f"Adhesion enrichment ratio (capped at {ratio_cap:.0f})")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def scatter_feature_vs_probability(adhesion_df: pd.DataFrame, feature_col: str, out_path: Path) -> None:
    """Scatter of feature_col vs. probability_adhesion, adhesion proteins only."""
    sub = adhesion_df.dropna(subset=[feature_col, "probability_adhesion"])
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(sub[feature_col], sub["probability_adhesion"], alpha=0.3, s=10, color=_FULL_COLOR)
    ax.set_xlabel(feature_col.replace("_", " "))
    ax.set_ylabel("Adhesion probability")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def scatter_species_mean_length_vs_fraction(
    adhesion_df: pd.DataFrame, kingdom_summary_df: pd.DataFrame, out_path: Path
) -> None:
    """Scatter of per-species mean adhesion-protein length vs. that
    species' adhesion_fraction (from kingdom_survey's summary table)."""
    mean_length = adhesion_df.groupby("locustag")["length"].mean().rename("mean_length")
    merged = kingdom_summary_df.set_index("locustag").join(mean_length, how="inner")

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(merged["mean_length"], merged["adhesion_fraction"], alpha=0.4, s=15, color=_FULL_COLOR)
    ax.set_xlabel("Mean adhesion-protein length (species)")
    ax.set_ylabel("Species adhesion fraction")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/adhesion_properties/test_figures_ext.py -v`
Expected: all PASS.

- [ ] **Step 5: Lint**

Run: `ruff check analysis/adhesion_properties/figures_ext.py tests/adhesion_properties/test_figures_ext.py`

- [ ] **Step 6: Commit**

```bash
git add analysis/adhesion_properties/figures_ext.py tests/adhesion_properties/test_figures_ext.py
git commit -m "Add figure-generation module for adhesion protein properties"
```

---

### Task 11: `05_make_figures.py` — run on real data

**Files:**
- Create: `analysis/adhesion_properties/05_make_figures.py`

**Interfaces:**
- Consumes: `figures_ext.py` (Task 10), tables from Tasks 6, 9.
- Produces: `figures/box_length_by_phylum.png`, `figures/box_pct_ser_thr_pro_by_phylum.png`, `figures/bar_{flag}_by_phylum.png` (5 flags), `figures/top_domains_{pfam,cazy,merops}.png`, `figures/scatter_length_vs_probability.png`, `figures/scatter_pct_ser_thr_pro_vs_probability.png`, `figures/scatter_species_mean_length_vs_fraction.png` — consumed by Task 12.

- [ ] **Step 1: Write the CLI script**

Create `analysis/adhesion_properties/05_make_figures.py`:

```python
#!/usr/bin/env python
"""Generate all figures for the adhesion-protein properties report."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figures_ext import (  # noqa: E402
    domain_fraction_barplot,
    paired_boxplot_by_rank,
    scatter_feature_vs_probability,
    scatter_species_mean_length_vs_fraction,
    top_domain_enrichment_barplot,
)

TABLES_DIR = Path(__file__).resolve().parent / "tables"
FIGURES_DIR = Path(__file__).resolve().parent / "figures"
KINGDOM_SUMMARY_CSV = (
    Path(__file__).resolve().parents[1] / "kingdom_survey" / "tables" / "species_adhesion_summary.csv"
)

DOMAIN_FLAGS = ["has_pfam", "has_cazy", "has_merops", "has_signal_peptide", "has_tm_helix"]


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(TABLES_DIR / "adhesion_protein_properties.csv")
    adhesion = master[master["group"] == "adhesion"]

    paired_boxplot_by_rank(master, "phylum", "length", FIGURES_DIR / "box_length_by_phylum.png")
    paired_boxplot_by_rank(
        master, "phylum", "pct_ser_thr_pro", FIGURES_DIR / "box_pct_ser_thr_pro_by_phylum.png"
    )

    fractions = pd.read_csv(TABLES_DIR / "clade_phylum_domain_fractions.csv")
    for flag in DOMAIN_FLAGS:
        domain_fraction_barplot(fractions, "phylum", flag, FIGURES_DIR / f"bar_{flag}_by_phylum.png")

    for domain_type in ["pfam", "cazy", "merops"]:
        top = pd.read_csv(TABLES_DIR / f"top_domains_{domain_type}.csv")
        if top.empty:
            print(f"Skipping top_domains_{domain_type}.png — no enriched domains found")
            continue
        top_domain_enrichment_barplot(
            top,
            FIGURES_DIR / f"top_domains_{domain_type}.png",
            title=f"Top enriched {domain_type.upper()} domains among adhesion proteins",
        )

    scatter_feature_vs_probability(adhesion, "length", FIGURES_DIR / "scatter_length_vs_probability.png")
    scatter_feature_vs_probability(
        adhesion, "pct_ser_thr_pro", FIGURES_DIR / "scatter_pct_ser_thr_pro_vs_probability.png"
    )

    kingdom_summary = pd.read_csv(KINGDOM_SUMMARY_CSV)
    scatter_species_mean_length_vs_fraction(
        adhesion, kingdom_summary, FIGURES_DIR / "scatter_species_mean_length_vs_fraction.png"
    )

    print(f"Wrote figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python analysis/adhesion_properties/05_make_figures.py`

- [ ] **Step 3: Visual sanity check**

View each PNG in `analysis/adhesion_properties/figures/` (via the Read tool). Confirm: axis labels present and legible, boxplot phylum labels aren't overlapping illegibly, the top-domain bar charts show real domain names (not empty), and the scatter plots show real point clouds (not blank/garbled).

- [ ] **Step 4: Lint**

Run: `ruff check analysis/adhesion_properties/05_make_figures.py`

- [ ] **Step 5: Commit**

```bash
git add analysis/adhesion_properties/05_make_figures.py
git commit -m "Add CLI to generate figures for adhesion protein properties"
```

---

### Task 12: `06_generate_report.py` — assemble and write `REPORT.md`

**Files:**
- Create: `analysis/adhesion_properties/06_generate_report.py`
- Create: `analysis/adhesion_properties/REPORT.md` (generated, then hand-edited)

**Interfaces:**
- Consumes: every table under `tables/` and every figure under `figures/` produced by Tasks 4, 6, 7, 9, 11.
- Produces: the final `REPORT.md` deliverable.

- [ ] **Step 1: Write the report-assembly script**

Create `analysis/adhesion_properties/06_generate_report.py`:

```python
#!/usr/bin/env python
"""Assemble REPORT.md from the tables/figures already generated by 01-05.
Pulls all numbers directly from the CSVs — the headline-findings prose is
written by hand in Step 3 after inspecting the real output, same pattern
as kingdom_survey's report."""
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent
TABLES = BASE / "tables"


def main() -> None:
    universe = pd.read_csv(TABLES / "protein_universe.csv")
    n_adhesion = int((universe["group"] == "adhesion").sum())
    n_background = int((universe["group"] == "background").sum())

    flags = pd.read_csv(TABLES / "protein_domain_flags.csv")
    merged = universe.merge(flags, on="protein_id", how="inner")
    adhesion_flags = merged[merged["group"] == "adhesion"]
    background_flags = merged[merged["group"] == "background"]

    domain_summary_rows = []
    for flag in ["has_pfam", "has_cazy", "has_merops", "has_signal_peptide", "has_tm_helix"]:
        domain_summary_rows.append(
            {
                "feature": flag,
                "adhesion_fraction": adhesion_flags[flag].mean(),
                "background_fraction": background_flags[flag].mean(),
            }
        )
    domain_summary_md = pd.DataFrame(domain_summary_rows).to_string(index=False)

    top_pfam = pd.read_csv(TABLES / "top_domains_pfam.csv").head(10)
    top_pfam_md = top_pfam.to_string(index=False) if not top_pfam.empty else "(none found)"

    mw_length_phylum = pd.read_csv(TABLES / "clade_phylum_length_mannwhitney.csv")
    mw_stp_phylum = pd.read_csv(TABLES / "clade_phylum_pct_ser_thr_pro_mannwhitney.csv")

    report = f"""# Adhesion Protein Molecular Properties & Functional Domain Survey

## Scope

- {n_adhesion} adhesion-predicted proteins across the matched species set.
- {n_background} background (non-adhesion) proteins sampled 1:1 per species
  for comparison (seed=42, reproducible).

## Domain/topology presence: adhesion vs. background

```
{domain_summary_md}
```

## Length and Ser/Thr/Pro content: adhesion vs. background, by phylum

Mann-Whitney U test (adhesion vs. background, within each well-powered
phylum), BH-corrected across phyla:

**Length:**
```
{mw_length_phylum.to_string(index=False)}
```

**Ser+Thr+Pro %:**
```
{mw_stp_phylum.to_string(index=False)}
```

![Length by phylum](figures/box_length_by_phylum.png)

![Ser/Thr/Pro composition by phylum](figures/box_pct_ser_thr_pro_by_phylum.png)

## Domain presence by phylum

![Pfam presence by phylum](figures/bar_has_pfam_by_phylum.png)

## Top enriched Pfam domains among adhesion proteins

```
{top_pfam_md}
```

![Top enriched Pfam domains](figures/top_domains_pfam.png)

![Top enriched CAZy families](figures/top_domains_cazy.png)

![Top enriched MEROPS hits](figures/top_domains_merops.png)

## Correlation with model confidence

![Length vs. adhesion probability](figures/scatter_length_vs_probability.png)

![Ser/Thr/Pro vs. adhesion probability](figures/scatter_pct_ser_thr_pro_vs_probability.png)

![Species mean length vs. adhesion fraction](figures/scatter_species_mean_length_vs_fraction.png)

## Headline findings

_(Written by hand after inspecting the real tables/figures above: name the
actual direction and magnitude of the Ser/Thr/Pro-content difference
between adhesion and background proteins and in which phyla it's
statistically significant after BH correction, whether the classifier's
probability correlates with length or Ser/Thr/Pro content, and which
specific Pfam/CAZy/MEROPS hits are most enriched and what they suggest
about function — e.g. GPI-anchor-related or glycosyl hydrolase domains
would be a biologically meaningful pattern worth calling out by name.)_

## Caveats

1. **MEROPS results are accession-level, not family-level** — this
   database has no peptidase family/clan lookup, so MEROPS hit identities
   are individual sequence accessions (e.g. `MER0080922`), not the
   family codes (e.g. "S08") domain experts would recognize; treat the
   MEROPS enrichment table as "has a hit in this specific reference
   sequence's neighborhood," not a family-level functional claim.
2. **`net_charge_ph7` is a coarse approximation** (simple Lys+Arg minus
   Asp+Glu residue count), not a full pKa-based charge model.
3. **Background sampling is random** (seeded, reproducible) — a different
   seed would produce a slightly different background set; the 1:1
   per-species matching mitigates but doesn't eliminate this.
4. **Non-independence across clades** (species within a genus/family
   share ancestry) applies here exactly as it does in the kingdom-wide
   survey — this pass does not repeat that survey's genus-averaging/
   mixed-model correction; treat clade-level patterns here as
   descriptive/hypothesis-generating, not confirmatory.
5. **Multiple testing across features and clades**: the Mann-Whitney
   tables are BH-corrected *within* each (rank, feature) combination, not
   jointly across every feature/rank combination in this report — a
   stricter joint correction would raise the significance bar further.

## Follow-up

Whether an ESM-2 embedding-space clustering reveals distinct sub-types
within the adhesion-predicted set (e.g. FLO11-like vs. Als-like vs. other
architectures) is deferred to a follow-up phase, informed by which
domain/property patterns found here are worth validating against.
"""
    (BASE / "REPORT.md").write_text(report)
    print(f"Wrote {BASE / 'REPORT.md'}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python analysis/adhesion_properties/06_generate_report.py`

- [ ] **Step 3: Fill in the hand-written section**

Open `analysis/adhesion_properties/REPORT.md`, read the actual Mann-Whitney tables and top-Pfam-domain table already inlined, and replace the `_(Written by hand ...)_` placeholder under "Headline findings" with a concrete paragraph naming: the actual direction/magnitude of the Ser/Thr/Pro difference and which phyla it's significant in (after BH correction), whether length/Ser-Thr-Pro correlate visually with probability in the scatter plots (view them), and the 2-3 most enriched specific Pfam/CAZy domains by name with a one-sentence note on what they suggest about adhesin function if recognizable (e.g. a GPI-anchor, lectin, or cell-wall-related domain).

- [ ] **Step 4: Verify the report renders**

Confirm the file is non-empty and every referenced image path actually exists under `figures/`.

- [ ] **Step 5: Commit**

```bash
git add analysis/adhesion_properties/06_generate_report.py analysis/adhesion_properties/REPORT.md
git add -f analysis/adhesion_properties/tables/*.csv analysis/adhesion_properties/figures/*.png
git commit -m "Add adhesion protein properties report, tables, and figures"
```

(The `-f` is needed because `tables/*.csv` and `figures/*.png` are gitignored by default — this is the one deliberate exception, mirroring `kingdom_survey`'s Task 8: the exact tables/figures the committed report cites are checked in so the report is self-contained and reproducible from this commit.)

---

## Self-Review Notes

- **Spec coverage:** every section of the design doc has a task — data-source facts (Task 1's `db.py` docstring/tests), protein universe + background sampling (Tasks 3-4), sequence properties (Tasks 2, 7), domain flags + top-domain enrichment (Tasks 5-6), per-clade aggregation + all 3 statistical framings (Task 8-9: `aggregate_by_clade` per-group, `kruskal_wallis_by_rank` reused for the within-adhesion/within-background clade-variation check, `mannwhitney_within_clade` for the adhesion-vs-background-within-clade check), all 7 figure types (Tasks 10-11), report + caveats (Task 12).
- **Placeholder scan:** the only literal placeholder text (`_(Written by hand ...)_`) is inside the generated `REPORT.md`, not this plan, and Task 12 Step 3 gives concrete instructions for what replaces it, sourced from real numbers already inlined in the same generated file — consistent with how `kingdom_survey`'s Task 8 handled the identical situation.
- **Type/name consistency:** checked `universe.build_protein_universe`'s row dict keys match `01_build_protein_universe.py`'s `FIELDS` list; `domains.compute_domain_flags`'s output columns match `04_aggregate_and_test.py`'s `DOMAIN_FLAGS` list and `figures_ext.domain_fraction_barplot`'s expected `flag_col`; `properties.sequence_properties`'s returned keys match `04_aggregate_and_test.py`'s `CONTINUOUS_FEATURES` list; `stats_ext.py` re-exports `kruskal_wallis_by_rank` so `04_aggregate_and_test.py` can import it from `stats_ext` alone rather than needing two separate `sys.path` insertions.
