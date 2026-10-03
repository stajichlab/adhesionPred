# D8 curated GPI rows: TM rule, conflict file, input checks (plan 1 of 2 for issue #50)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make step 03 (D8 triage) treat curated GPI literature rows as the owner decided: a TM feature blocks a literature row unless the owner sets `override_tm=yes`; blocked rows go to a review file; the input file is checked.

**Architecture:** Add a `read_curated_gpi` reader, a `check_curated_gpi` match check and a `tm_conflicts` helper to `d8_triage.py` (standard library only). Change `classify_pm` to take an `override_tm` flag. Step 03 reads the validated rows, passes a dictionary instead of a set, and writes two new files.

**Tech Stack:** Python 3.12, pytest, ruff 0.3.5. No network use in the tests (the existing tests pass a fake `fetch`).

**Spec:** `docs/superpowers/specs/2026-10-02-step1-basidiomycota-and-gpi-truth-design.md`, sections 3.2, 4, 5, 9 and 11 (E3, and the step 03 part of E4). Owner decisions 4, 7, 8 and 12 apply.

**Scope.** This plan covers step 03 and `curated_gpi.tsv` only. It does **not** cover the curated Basidiomycota file, the step 01 merge, the new truth columns, or the Phase C tier strata. Those are plan 2. Plan 2 is written after the first curated rows exist, so that its code is tested on real rows. This plan changes no result of the 2026-10-02 run, because `curated_gpi.tsv` has no rows.

## Global Constraints

- Run Python with `/usr/bin/python3.12` (or `python3.12`). Tests run as `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py` from the repo root.
- `d8_triage.py` uses the standard library only (module docstring).
- Reviewed UniProt entries with ECO:0000269 GPI evidence keep winning over a TM feature. This behaviour does not change.
- `override_tm` is `yes` or `no`. The default for a new row is `no`.
- `gene_id` in `curated_gpi.tsv` is the native identifier of the source (for example SGD `S000004924`). D8 matches on `(source_id, gene_id)`.
- A `curated_gpi.tsv` without the new columns stops step 03 with exit code 2 and a `STOP:` message.
- Predictor output never counts as evidence for a curated row (decision Q2). Nothing in this plan reads predictor output.
- New output files need a heading in `analysis/step1_compare/COLUMNS.md` (`tests/step1_compare/test_paths.py::test_every_output_file_has_a_columns_md_heading` checks this).
- Documents use Simplified Technical English: short sentences, active voice, no idioms.
- Lint with `uvx ruff@0.3.5 check analysis tests` (CI pins ruff 0.3.5). The pre-commit hook may reformat a file and fail the commit once. Stage and commit again.
- Every commit message ends with: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.

## Review Focus

Each line has a test in the task named in brackets.

1. A `curated_gpi.tsv` with the old five-column header stops step 03 with a clear message, and writes no output. [Task 2, Task 3]
2. `override_tm` with any value other than `yes` or `no` (for example `Yes` or empty) stops the run. [Task 2]
3. The same `(source_id, gene_id)` twice stops the run. [Task 2]
4. A row whose gene is missing from the truth set, outside P-ext, or not a plasma-membrane candidate is listed in `curated_gpi_unmatched.tsv`, and the run still succeeds. [Task 2, Task 3]
5. A literature row on a gene that has a TM feature **and** a reviewed UniProt ECO:0000269 GPI entry gives P-gpi and no conflict row. [Task 1, Task 3]
6. A literature row on a TM gene with `override_tm=no` gives PM-TM and one conflict row. With `override_tm=yes` it gives P-gpi and no conflict row. [Task 1, Task 3]

## File Structure

| File | Change |
|---|---|
| `analysis/step1_compare/d8_triage.py` | add `CURATED_GPI_COLUMNS`, `CuratedGpiError`, `read_curated_gpi`, `check_curated_gpi`, `tm_conflicts`; change `classify_pm` |
| `analysis/step1_compare/03_triage_pm.py` | read validated rows; dictionary instead of set; write `d8_curated_conflicts.tsv` and `curated_gpi_unmatched.tsv` |
| `analysis/step1_compare/curated_gpi.tsv` | new header (no rows) |
| `analysis/step1_compare/COLUMNS.md`, `README.md` | document the two new files and the new columns |
| `tests/step1_compare/test_d8.py` | new tests; five changes to existing tests, found by content (see Task 3) |

---

### Task 1: `classify_pm` with `override_tm`

**Files:**
- Modify: `analysis/step1_compare/d8_triage.py:108-120`
- Test: `tests/step1_compare/test_d8.py` (add after `test_classify_pm`, line 68)

**Interfaces:**
- Consumes: `UniprotEvidence` (existing), constants `P_GPI`, `PM_TM`, `PM_UNRESOLVED`.
- Produces: `classify_pm(entries: list[UniprotEvidence], literature: bool, override_tm: bool = False) -> tuple[str, str]`. Task 3 calls it.

- [ ] **Step 1: Write the failing test.** Add to `tests/step1_compare/test_d8.py` after `test_classify_pm`:

```python
def test_classify_pm_literature_row_is_blocked_by_tm_unless_override():
    p = d8_triage.parse_entry
    cls, reason = d8_triage.classify_pm([p(MSB2)], True)
    assert cls == "PM-TM"
    assert "blocked" in reason and "P32334 1 TM ECO:0000255" in reason
    assert d8_triage.classify_pm([p(MSB2)], True, override_tm=True) == (
        "P-gpi",
        "literature row in curated_gpi.tsv",
    )
    assert d8_triage.classify_pm([p(YPS1)], True)[0] == "P-gpi"  # no TM feature
    assert d8_triage.classify_pm([], True)[0] == "P-gpi"  # no UniProt entry at all


def test_classify_pm_uniprot_experimental_gpi_still_wins_over_tm():
    both = d8_triage.parse_entry(entry("P5", True, gpi_eco=["ECO:0000269"], tm=["ECO:0000255"]))
    assert d8_triage.classify_pm([both], False)[0] == "P-gpi"
    assert d8_triage.classify_pm([both], True)[0] == "P-gpi"
    assert d8_triage.classify_pm([both], True, override_tm=False)[1] == (
        "literature row in curated_gpi.tsv"
    )
```

- [ ] **Step 2: Run the test to verify it fails.**

Run: `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py::test_classify_pm_literature_row_is_blocked_by_tm_unless_override -q`
Expected: FAIL with `AssertionError` (`assert 'P-gpi' == 'PM-TM'`), because the old `classify_pm` returns P-gpi for any literature row. The keyword argument is not reached.

- [ ] **Step 3: Write the implementation.** Replace `classify_pm` in `analysis/step1_compare/d8_triage.py`:

```python
def classify_pm(
    entries: list[UniprotEvidence], literature: bool, override_tm: bool = False
) -> tuple[str, str]:
    """Class and reason for one P-ext gene with a plasma-membrane term.

    A literature row gives P-gpi, except when a UniProt TM feature blocks it: the gene has a TM
    feature, `override_tm` is false, and no reviewed UniProt entry has experimental GPI evidence.
    A blocked gene stays PM-TM (owner decision 8, 2026-10-02). The reason says so."""
    with_tm = [e for e in entries if e.tm_count > 0]
    uniprot_curated = any(e.curated_gpi for e in entries)
    blocked = literature and bool(with_tm) and not override_tm and not uniprot_curated
    if literature and not blocked:
        return P_GPI, "literature row in curated_gpi.tsv"
    for e in entries:
        if e.curated_gpi:
            return P_GPI, f"{e.accession} reviewed GPI-anchor {','.join(e.gpi_eco)}"
    if with_tm:
        e = with_tm[0]
        reason = f"{e.accession} {e.tm_count} TM {','.join(e.tm_eco) or 'no ECO'}"
        if blocked:
            reason += "; literature row blocked by the TM feature (override_tm=no)"
        return PM_TM, reason
    if not entries:
        return PM_UNRESOLVED, "no UniProt entry found"
    return PM_UNRESOLVED, "no curated GPI evidence and no TM feature"
```

- [ ] **Step 4: Run the tests to verify they pass.**

Run: `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py -q -k "classify_pm or pm_tm_needs"`
Expected: PASS (the existing `test_classify_pm` and `test_pm_tm_needs_a_tm_feature` still pass, because they use `literature=False` or an entry without TM).

- [ ] **Step 5: Commit.**

```bash
git add analysis/step1_compare/d8_triage.py tests/step1_compare/test_d8.py
git commit -m "feat(d8): a TM feature blocks a literature GPI row unless override_tm is set

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 2: reader, match check and conflict helper in `d8_triage.py`

**Files:**
- Modify: `analysis/step1_compare/d8_triage.py` (imports at lines 16-17; add code at the end of the file)
- Test: `tests/step1_compare/test_d8.py` (add after the tests of Task 1)

**Interfaces:**
- Consumes: `PM_TM` constant.
- Produces (Task 3 calls these):
  - `CURATED_GPI_COLUMNS: tuple[str, ...]` (12 names, order below)
  - `class CuratedGpiError(ValueError)`
  - `read_curated_gpi(path) -> list[dict[str, str]]`
  - `check_curated_gpi(curated_rows, truth_rows) -> list[dict[str, str]]` with columns `UNMATCHED_COLUMNS`
  - `tm_conflicts(triage_rows, literature) -> list[dict[str, str]]` with columns `CONFLICT_COLUMNS`

- [ ] **Step 1: Write the failing tests.** Add to `tests/step1_compare/test_d8.py`. First the helpers near the top of the file, after the line `BARE = entry("Q00001", True, sgd="S000000010")`:

```python
CURATED_HEADER = "\t".join(d8_triage.CURATED_GPI_COLUMNS) + "\n"


def curated_row(gene_id, symbol, override="no", source_id="Scer", **changes):
    values = {
        "source_id": source_id,
        "gene_id": gene_id,
        "symbol": symbol,
        "pmid": "12345678",
        "note": "n",
        "species": "Saccharomyces cerevisiae",
        "uniprot_accession": "",
        "evidence_level": "direct",
        "evidence_note": "quoted sentence; retrieved 2026-10-02",
        "reviewer": "test",
        "review_date": "2026-10-02",
        "override_tm": override,
    }
    values.update(changes)
    return "\t".join(values[c] for c in d8_triage.CURATED_GPI_COLUMNS) + "\n"
```

Then the tests, after the Task 1 tests:

```python
def test_read_curated_gpi_accepts_a_valid_file_and_a_header_only_file(tmp_path):
    path = tmp_path / "curated_gpi.tsv"
    path.write_text(CURATED_HEADER)
    assert d8_triage.read_curated_gpi(path) == []
    path.write_text(CURATED_HEADER + curated_row("S000003246", "MSB2", "yes"))
    rows = d8_triage.read_curated_gpi(path)
    assert [(r["gene_id"], r["override_tm"]) for r in rows] == [("S000003246", "yes")]


@pytest.mark.parametrize(
    "text, message",
    [
        ("source_id\tgene_id\tsymbol\tpmid\tnote\n", "missing columns"),  # old header
        (CURATED_HEADER + curated_row("G1", "A", override="Yes"), "override_tm"),
        (CURATED_HEADER + curated_row("G1", "A", override=""), "override_tm"),
        (CURATED_HEADER + curated_row("G1", "A", evidence_level="maybe"), "evidence_level"),
        (CURATED_HEADER + curated_row("G1", "A", pmid=""), "pmid"),
        (CURATED_HEADER + curated_row("G1", "A", reviewer=""), "reviewer"),
        (CURATED_HEADER + curated_row("", "A"), "source_id and gene_id"),
        (CURATED_HEADER + curated_row("G1", "A") + curated_row("G1", "A"), "twice"),
    ],
)
def test_read_curated_gpi_rejects_bad_input(tmp_path, text, message):
    path = tmp_path / "curated_gpi.tsv"
    path.write_text(text)
    with pytest.raises(d8_triage.CuratedGpiError, match=message):
        d8_triage.read_curated_gpi(path)


def test_check_curated_gpi_lists_rows_that_cannot_act():
    truth = [
        {"source_id": "Scer", "gene_id": "G1", "label": "P-ext", "pm_candidate": "yes"},
        {"source_id": "Scer", "gene_id": "G2", "label": "P-ext", "pm_candidate": "no"},
        {"source_id": "Scer", "gene_id": "G3", "label": "ambiguous", "pm_candidate": "no"},
    ]
    curated = [
        {"source_id": "Scer", "gene_id": g, "symbol": g.lower()} for g in ("G1", "G2", "G3", "G4")
    ]
    got = d8_triage.check_curated_gpi(curated, truth)
    assert [(r["gene_id"], r["reason"], r["label"]) for r in got] == [
        ("G2", "not_pm_candidate", "P-ext"),
        ("G3", "outside_p_ext", "ambiguous"),
        ("G4", "no_truth_gene", ""),
    ]


def test_tm_conflicts_lists_only_blocked_literature_genes():
    def triage_row(gene_id, symbol, d8_class, accession, tm_count, tm_eco):
        return {
            "source_id": "Scer",
            "gene_id": gene_id,
            "symbol": symbol,
            "d8_class": d8_class,
            "uniprot_accessions": accession,
            "tm_count": tm_count,
            "tm_eco": tm_eco,
            "d8_reason": "r",
        }

    triage = [
        triage_row("G1", "A", "PM-TM", "P1", "1", "ECO:0000255"),
        triage_row("G2", "B", "PM-TM", "P2", "2", "ECO:0000255"),
        triage_row("G3", "C", "P-gpi", "P3", "0", ""),
    ]
    literature = {("Scer", "G1"): False, ("Scer", "G3"): True}
    got = d8_triage.tm_conflicts(triage, literature)
    assert [(r["gene_id"], r["override_tm"]) for r in got] == [("G1", "no")]
```

- [ ] **Step 2: Run the tests to verify they fail.**

Run: `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py -q -k "curated_gpi or tm_conflicts"`
Expected: FAIL with `AttributeError: module 'd8_triage' has no attribute 'CURATED_GPI_COLUMNS'` (the whole file fails at import of the helper line).

- [ ] **Step 3: Write the implementation.** In `analysis/step1_compare/d8_triage.py`, change the imports and add the code. Replace lines 16-17:

```python
import csv
import urllib.parse
from dataclasses import dataclass, field
```

Add at the end of the file:

```python
CURATED_GPI_COLUMNS = (
    "source_id",
    "gene_id",
    "symbol",
    "pmid",
    "note",
    "species",
    "uniprot_accession",
    "evidence_level",
    "evidence_note",
    "reviewer",
    "review_date",
    "override_tm",
)
CURATED_GPI_REQUIRED = ("pmid", "evidence_note", "reviewer", "review_date")
EVIDENCE_LEVELS = frozenset({"direct", "transfer"})
UNMATCHED_COLUMNS = ("source_id", "gene_id", "symbol", "reason", "label")
CONFLICT_COLUMNS = (
    "source_id",
    "gene_id",
    "symbol",
    "uniprot_accessions",
    "tm_count",
    "tm_eco",
    "d8_reason",
    "override_tm",
)


class CuratedGpiError(ValueError):
    """curated_gpi.tsv has a missing column, a bad value or a repeated gene."""


def read_curated_gpi(path) -> list[dict[str, str]]:
    """Read and check curated_gpi.tsv. A header-only file gives an empty list.

    The PMID is checked for presence only. The reviewer opens each PMID (spec section 4)."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in CURATED_GPI_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise CuratedGpiError(f"{path}: missing columns {missing}")
        rows = list(reader)
    seen = set()
    for number, row in enumerate(rows, start=2):
        where = f"{path} line {number}"
        key = (row["source_id"], row["gene_id"])
        if not all(key):
            raise CuratedGpiError(f"{where}: source_id and gene_id are required")
        if key in seen:
            raise CuratedGpiError(f"{where}: {key[0]} {key[1]} appears twice")
        seen.add(key)
        if row["override_tm"] not in ("yes", "no"):
            raise CuratedGpiError(f"{where}: override_tm must be yes or no")
        if row["evidence_level"] not in EVIDENCE_LEVELS:
            raise CuratedGpiError(f"{where}: evidence_level must be direct or transfer")
        for column in CURATED_GPI_REQUIRED:
            if not row[column]:
                raise CuratedGpiError(f"{where}: {column} is required")
    return rows


def check_curated_gpi(curated_rows, truth_rows) -> list[dict[str, str]]:
    """Rows of curated_gpi.tsv that cannot change a D8 class, with the reason.

    D8 reads a row only for a P-ext gene that is a plasma-membrane candidate."""
    truth = {(r["source_id"], r["gene_id"]): r for r in truth_rows}
    problems = []
    for row in curated_rows:
        gene = truth.get((row["source_id"], row["gene_id"]))
        if gene is None:
            reason, label = "no_truth_gene", ""
        elif gene["label"] != "P-ext":
            reason, label = "outside_p_ext", gene["label"]
        elif gene["pm_candidate"] != "yes":
            reason, label = "not_pm_candidate", gene["label"]
        else:
            continue
        problems.append(
            {
                "source_id": row["source_id"],
                "gene_id": row["gene_id"],
                "symbol": row["symbol"],
                "reason": reason,
                "label": label,
            }
        )
    return problems


def tm_conflicts(triage_rows, literature) -> list[dict[str, str]]:
    """Literature genes that a TM feature blocked. They stay PM-TM until the owner reviews them.

    `literature` maps (source_id, gene_id) to the row's override_tm flag. A PM-TM gene with a
    literature row is always a blocked gene, because override_tm=yes would have given P-gpi."""
    return [
        {
            "source_id": t["source_id"],
            "gene_id": t["gene_id"],
            "symbol": t["symbol"],
            "uniprot_accessions": t["uniprot_accessions"],
            "tm_count": t["tm_count"],
            "tm_eco": t["tm_eco"],
            "d8_reason": t["d8_reason"],
            "override_tm": "no",
        }
        for t in triage_rows
        if t["d8_class"] == PM_TM and (t["source_id"], t["gene_id"]) in literature
    ]
```

- [ ] **Step 4: Run the tests to verify they pass.**

Run: `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py -q -k "curated_gpi or tm_conflicts or classify_pm"`
Expected: PASS. Other tests in the file may fail until Task 3 changes `_work` (they still use the old header). Do not run the whole file yet.

- [ ] **Step 5: Commit.**

```bash
git add analysis/step1_compare/d8_triage.py tests/step1_compare/test_d8.py
git commit -m "feat(d8): read, check and report curated GPI rows

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 3: step 03 uses the validated rows and writes two new files

**Files:**
- Modify: `analysis/step1_compare/03_triage_pm.py` (lines 41-48, 77-98, 197-234)
- Modify: `tests/step1_compare/test_d8.py` (five existing places, found by their content, not by line number: Tasks 1 and 2 add about 127 lines first)
- Test: `tests/step1_compare/test_d8.py` (new tests at the end)

**Interfaces:**
- Consumes: `d8_triage.read_curated_gpi`, `check_curated_gpi`, `tm_conflicts`, `CONFLICT_COLUMNS`, `UNMATCHED_COLUMNS`, `classify_pm(entries, literature, override_tm)` from Tasks 1 and 2.
- Produces: `triage_source(sp, rows, literature, fetch)` where `literature` is `dict[tuple[str, str], bool]` (the value is `override_tm == "yes"`). Output files `d8_curated_conflicts.tsv` and `curated_gpi_unmatched.tsv`.

- [ ] **Step 1: Update the existing tests and write the new failing tests.** In `tests/step1_compare/test_d8.py`:

Find each existing edit point by its content (the line numbers of the original file are about 127 lines lower than in the file after Tasks 1 and 2: the original lines were 118, 142, 181-187, 334-344 and 429).

1. In `test_triage_source_and_apply`: replace `triage.triage_source(sp, rows, set(), fake_fetch)` with `triage.triage_source(sp, rows, {}, fake_fetch)`.
2. In `test_d8_counts_report_no_uniprot_entry_and_gpi_feature_no_evidence`: replace `triage.triage_source(sp, rows, set(), fetch)` with `triage.triage_source(sp, rows, {}, fetch)`.
3. In `_work`: replace `(tmp_path / "curated_gpi.tsv").write_text("source_id\tgene_id\tsymbol\tpmid\tnote\n")` with `(tmp_path / "curated_gpi.tsv").write_text(CURATED_HEADER)`.
4. Replace the `OUTPUTS` tuple (after the comment `# ---- hardening tests`) with:

```python
OUTPUTS = (
    "d8_triage.tsv",
    "d8_gpi_outside_pext.tsv",
    "d8_curated_conflicts.tsv",
    "curated_gpi_unmatched.tsv",
    "d8_counts.tsv",
    "truth_set_triaged.tsv.gz",
    "d8_run.json",
)
```

5. Replace the whole function `test_literature_row_makes_p_gpi_and_missing_curated_file_stops` (it runs from its `def` line to the `assert "STOP:" in capsys.readouterr().err` line; 11 lines in the original file) with:

```python
def test_literature_row_with_override_makes_p_gpi_and_missing_curated_file_stops(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    (tmp_path / "curated_gpi.tsv").write_text(
        CURATED_HEADER + curated_row("S000003246", "MSB2", override="yes")
    )
    assert triage.main(argv, fetch=good_fetch) == 0
    assert truth_table.read_tsv(tmp_path / "d8_triage.tsv")[0]["d8_class"] == "P-gpi"
    assert truth_table.read_tsv(tmp_path / "d8_curated_conflicts.tsv") == []
    (tmp_path / "curated_gpi.tsv").unlink()
    assert triage.main(argv, fetch=good_fetch) == 2
    assert "STOP:" in capsys.readouterr().err
```

Add the new tests at the end of the file:

```python
def test_literature_row_on_a_tm_gene_stays_pm_tm_and_is_reported(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    (tmp_path / "curated_gpi.tsv").write_text(
        CURATED_HEADER + curated_row("S000003246", "MSB2", override="no")
    )
    assert triage.main(argv, fetch=good_fetch) == 0
    row = truth_table.read_tsv(tmp_path / "d8_triage.tsv")[0]
    assert row["d8_class"] == "PM-TM" and "blocked" in row["d8_reason"]
    conflicts = truth_table.read_tsv(tmp_path / "d8_curated_conflicts.tsv")
    assert [(c["gene_id"], c["override_tm"], c["tm_count"]) for c in conflicts] == [
        ("S000003246", "no", "1")
    ]
    assert "d8_curated_conflicts.tsv" in capsys.readouterr().err


def test_literature_row_on_a_tm_gene_with_uniprot_gpi_evidence_gives_p_gpi_without_conflict(
    tmp_path,
):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    (tmp_path / "curated_gpi.tsv").write_text(
        CURATED_HEADER + curated_row("S000003246", "MSB2", override="no")
    )
    both = entry("P32334", True, gpi_eco=["ECO:0000269"], tm=["ECO:0000255"], sgd="S000003246")

    def fetch(url, tag):
        return [both], "2026_03"

    assert triage.main(argv, fetch=fetch) == 0
    assert truth_table.read_tsv(tmp_path / "d8_triage.tsv")[0]["d8_class"] == "P-gpi"
    assert truth_table.read_tsv(tmp_path / "d8_curated_conflicts.tsv") == []


def test_unmatched_curated_rows_are_listed_and_the_run_succeeds(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    (tmp_path / "curated_gpi.tsv").write_text(CURATED_HEADER + curated_row("S999999999", "NOPE"))
    assert triage.main(argv, fetch=good_fetch) == 0
    got = truth_table.read_tsv(tmp_path / "curated_gpi_unmatched.tsv")
    assert [(r["gene_id"], r["reason"]) for r in got] == [("S999999999", "no_truth_gene")]
    assert "curated_gpi_unmatched.tsv" in capsys.readouterr().err


def test_old_five_column_curated_file_stops_and_writes_nothing(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    (tmp_path / "curated_gpi.tsv").write_text("source_id\tgene_id\tsymbol\tpmid\tnote\n")
    assert triage.main(argv, fetch=good_fetch) == 2
    assert "missing columns" in capsys.readouterr().err
    assert not [n for n in OUTPUTS if (tmp_path / n).exists()]


def test_header_only_curated_file_gives_empty_review_files(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    assert triage.main(_work(tmp_path), fetch=good_fetch) == 0
    assert truth_table.read_tsv(tmp_path / "d8_curated_conflicts.tsv") == []
    assert truth_table.read_tsv(tmp_path / "curated_gpi_unmatched.tsv") == []
    assert "need review" not in capsys.readouterr().err
```

- [ ] **Step 2: Run the tests to verify they fail.**

Run: `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py -q`
Expected: FAIL. Several tests fail because step 03 does not yet write the two files or read the new columns (for example `AssertionError` on `(tmp_path / name).exists()`, and a `KeyError` or wrong class for the MSB2 literature row).

- [ ] **Step 3: Write the implementation.** In `analysis/step1_compare/03_triage_pm.py`:

(a) Replace the `OUTPUT_NAMES` tuple (lines 42-48):

```python
OUTPUT_NAMES = (
    "d8_triage.tsv",
    "d8_gpi_outside_pext.tsv",
    "d8_curated_conflicts.tsv",
    "curated_gpi_unmatched.tsv",
    "d8_counts.tsv",
    "truth_set_triaged.tsv.gz",
    "d8_run.json",
)
```

(b) In `triage_source`, rename the parameter and pass the flag. Replace the `def` line (line 77) and the `classify_pm` call (lines 96-98):

```python
def triage_source(sp, rows, literature, fetch):
```

```python
        key = (sp["source_id"], r["gene_id"])
        d8_class, reason = d8_triage.classify_pm(
            mine, key in literature, literature.get(key, False)
        )
```

Also update the docstring of the script (lines 4-6). Replace these three lines:

```
Reads $STEP1_WORKDIR/truth_set.tsv.gz, species.tsv and curated_gpi.tsv. Queries UniProtKB REST.
Writes d8_triage.tsv, d8_gpi_outside_pext.tsv, d8_counts.tsv, truth_set_triaged.tsv.gz and the
raw UniProt JSON pages (d8_uniprot/; pages of older runs are not deleted and can remain) to the
```

with:

```
Reads $STEP1_WORKDIR/truth_set.tsv.gz, species.tsv and curated_gpi.tsv (checked by
d8_triage.read_curated_gpi). Queries UniProtKB REST.
Writes d8_triage.tsv, d8_gpi_outside_pext.tsv, d8_curated_conflicts.tsv,
curated_gpi_unmatched.tsv, d8_counts.tsv, truth_set_triaged.tsv.gz and the
raw UniProt JSON pages (d8_uniprot/; pages of older runs are not deleted and can remain) to the
```

(c) In `main`, replace the block that builds `literature_ids` (lines 199-201):

```python
        curated_rows = d8_triage.read_curated_gpi(args.curated_gpi)
        literature = {
            (r["source_id"], r["gene_id"]): r["override_tm"] == "yes" for r in curated_rows
        }
```

(d) Change the call (line 207) to `t, o, c = triage_source(sp, rows, literature, fetch)`.

(e) After `triaged = apply_triage(truth_rows, triage)` (line 211), add:

```python
        conflicts = d8_triage.tm_conflicts(triage, literature)
        unmatched = d8_triage.check_curated_gpi(curated_rows, truth_rows)
```

(f) In the `tables` dictionary (line 223), add two entries after `"d8_gpi_outside_pext.tsv"`:

```python
            "d8_curated_conflicts.tsv": (d8_triage.CONFLICT_COLUMNS, conflicts),
            "curated_gpi_unmatched.tsv": (d8_triage.UNMATCHED_COLUMNS, unmatched),
```

(g) After the `except` block and before the `for c in counts:` print loop (line 245), add:

```python
    for name, rows_for_review in (
        ("d8_curated_conflicts.tsv", conflicts),
        ("curated_gpi_unmatched.tsv", unmatched),
    ):
        if rows_for_review:
            print(
                f"NOTE: {len(rows_for_review)} rows need review; see {name}",
                file=sys.stderr,
            )
```

- [ ] **Step 4: Run the tests to verify they pass.**

Run: `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py -q`
Expected: PASS for the whole file (the dry run counted 48 passed at this point).

Do not run `tests/step1_compare/test_paths.py` yet. Its heading test fails until Task 4 adds the two `COLUMNS.md` headings. If CI runs on every commit, commit Tasks 3 and 4 together.

- [ ] **Step 5: Commit.**

```bash
git add analysis/step1_compare/03_triage_pm.py tests/step1_compare/test_d8.py
git commit -m "feat(d8): step 03 checks curated_gpi.tsv, honours override_tm, writes review files

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 4: tracked header, documentation, full check

**Files:**
- Modify: `analysis/step1_compare/curated_gpi.tsv`
- Modify: `analysis/step1_compare/COLUMNS.md` (insert after the `d8_gpi_outside_pext.tsv` section, line 127)
- Modify: `analysis/step1_compare/README.md` (lines 50-51 and 195-196)
- Test: `tests/step1_compare/test_paths.py` (existing test checks the headings)

**Interfaces:**
- Consumes: the file names and columns from Tasks 2 and 3.
- Produces: the tracked header that step 03 accepts; documentation that the heading test checks.

- [ ] **Step 1: Write the failing check.** Run the heading test and a header test now:

Run: `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_paths.py::test_every_output_file_has_a_columns_md_heading -q`
Expected: FAIL with `COLUMNS.md has no heading for ['d8_curated_conflicts.tsv', 'curated_gpi_unmatched.tsv']`.

Add this test to `tests/step1_compare/test_d8.py` (it pins the tracked file to the reader):

```python
def test_tracked_curated_gpi_file_passes_the_reader():
    from conftest import STEP1_DIR

    assert d8_triage.read_curated_gpi(STEP1_DIR / "curated_gpi.tsv") == []
```

Run: `PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py::test_tracked_curated_gpi_file_passes_the_reader -q`
Expected: FAIL with `CuratedGpiError: ... missing columns`.

- [ ] **Step 2: Write the tracked header.**

Run:

```bash
python3.12 - <<'EOF'
import sys
sys.path.insert(0, "analysis/step1_compare")
import d8_triage
open("analysis/step1_compare/curated_gpi.tsv", "w").write("\t".join(d8_triage.CURATED_GPI_COLUMNS) + "\n")
EOF
```

- [ ] **Step 3: Write the documentation.** In `analysis/step1_compare/COLUMNS.md`, insert this text after the table of the `d8_gpi_outside_pext.tsv` section (after the line that starts `| gpi_eco | ECO codes of the GPI-anchor features of that entry`):

```markdown

## curated_gpi.tsv (input of 03_triage_pm.py, one row per gene)

Literature rows with experimental proof of a GPI anchor. A row acts only for a P-ext gene that is a
plasma-membrane candidate. Step 03 stops if a column is missing, if a value is invalid, or if a
`(source_id, gene_id)` pair appears twice. The file in the repo has a header and no rows.

| Column | Meaning |
|---|---|
| source_id | Source in `species.tsv`. |
| gene_id | The native identifier of the source (for example SGD `S000004924`). D8 matches on `(source_id, gene_id)`. |
| symbol | Gene symbol. |
| pmid | PubMed identifier. Required. The reviewer opens it and checks that it resolves. |
| note | Free text. |
| species | Species name, for readability. |
| uniprot_accession | UniProt accession. |
| evidence_level | `direct` (the paper measures the anchor in this protein) or `transfer` (in an ortholog). |
| evidence_note | The sentence of the paper that states the evidence, and the retrieval date. Required. |
| reviewer | Initials or model name of the second check. Required. |
| review_date | `YYYY-MM-DD`. Required. |
| override_tm | `no` (default) or `yes`. With `no`, a UniProt TM feature blocks the row and the gene stays `PM-TM`. The owner sets `yes` after review. |

## d8_curated_conflicts.tsv (03_triage_pm.py, one row per gene)

Genes that have a row in `curated_gpi.tsv` and a UniProt TM feature, with `override_tm` equal to `no`
and no reviewed UniProt entry with experimental GPI evidence. The TM feature blocks the row, so the
gene stays `PM-TM`. The owner reviews this file. If the owner agrees that the gene is GPI-anchored,
the owner sets `override_tm` to `yes` in `curated_gpi.tsv` and re-runs step 03. The file is empty
when no row is blocked.

| Column | Meaning |
|---|---|
| source_id, gene_id, symbol | As in `d8_triage.tsv`. |
| uniprot_accessions, tm_count, tm_eco, d8_reason | As in `d8_triage.tsv`. |
| override_tm | Always `no` in this file. |

## curated_gpi_unmatched.tsv (03_triage_pm.py, one row per curated row that cannot act)

D8 reads a `curated_gpi.tsv` row only for a P-ext gene that is a plasma-membrane candidate. This file
lists the rows that fail that condition. The run still succeeds. The file is empty when every row
can act.

| Column | Meaning |
|---|---|
| source_id, gene_id, symbol | From the `curated_gpi.tsv` row. |
| reason | `no_truth_gene` (no gene with this `source_id` and `gene_id` in `truth_set.tsv.gz`), `outside_p_ext` (the gene is not P-ext) or `not_pm_candidate` (the gene is P-ext without a non-IEA plasma-membrane term). |
| label | The label of the gene in `truth_set.tsv.gz`. Empty for `no_truth_gene`. |
```

In `analysis/step1_compare/README.md`, replace lines 50-51:

```
$PY 03_triage_pm.py                    # d8_triage.tsv, d8_gpi_outside_pext.tsv, d8_counts.tsv,
                                       #   truth_set_triaged.tsv.gz, d8_uniprot/, d8_run.json
```

with:

```
$PY 03_triage_pm.py                    # d8_triage.tsv, d8_gpi_outside_pext.tsv, d8_counts.tsv,
                                       #   d8_curated_conflicts.tsv, curated_gpi_unmatched.tsv,
                                       #   truth_set_triaged.tsv.gz, d8_uniprot/, d8_run.json
```

and replace the bullet at lines 195-196:

```
- `curated_gpi.tsv` has only a header row. This is by design. Literature curation fills it
  later as separate work.
```

with:

```
- `curated_gpi.tsv` has only a header row. This is by design. Literature curation fills it
  later as separate work (issue #50). Step 03 stops if a column is missing or a value is invalid.
  A row gives P-gpi unless a UniProt TM feature blocks it. The owner reviews blocked rows in
  `d8_curated_conflicts.tsv` and sets `override_tm=yes` to allow a row.
```

- [ ] **Step 4: Run all checks.**

Run:

```bash
PYTHONPATH=src python3.12 -m pytest tests/step1_compare/test_d8.py tests/step1_compare/test_paths.py -q
uvx ruff@0.3.5 check analysis tests
uvx ruff@0.3.5 format --check analysis tests
```

Expected: all tests pass; ruff reports no findings. If `ruff format --check` lists a file that this plan changed, run `uvx ruff@0.3.5 format <file>` and stage it.

Then run the Phase C and extract tests that read step 03 outputs, to check nothing else depends on the old header (about 5 minutes):

```bash
PYTHONPATH=src python3.12 -m pytest tests/step1_compare -q -x
```

Expected: PASS. A failure in a test that does not mention `curated_gpi` or `d8` means the plan broke something else. Stop and report it.

- [ ] **Step 5: Commit.**

```bash
git add analysis/step1_compare/curated_gpi.tsv analysis/step1_compare/COLUMNS.md analysis/step1_compare/README.md tests/step1_compare/test_d8.py
git commit -m "docs(d8): new curated_gpi.tsv header and the two review files

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage (section 3.2, 5, 9, 11 E3 and the step 03 part of E4):**
- TM blocks a literature row unless `override_tm=yes`: Task 1, Task 3.
- UniProt ECO:0000269 behaviour unchanged: Task 1 (second test), Task 3.
- `d8_curated_conflicts.tsv` written by step 03: Task 3.
- New columns, `evidence_level`, `override_tm` default, header check: Task 2, Task 4.
- Match check and `curated_gpi_unmatched.tsv` with the three reasons: Task 2, Task 3.
- Changes to the five existing test places in `test_d8.py` (original lines 118, 142, 181-187, 334-344, 429): Task 3.
- COLUMNS.md headings required by `test_paths.py`: Task 4.
- Not covered here, by design: step 01 merge, `curated_basidiomycota.tsv`, new truth columns, Phase C tiers, `dedupe.merge_group` (plan 2).

**Placeholder scan:** no TBD, no "add validation" without code. Every code step shows the code.

**Type consistency:** `literature` is `dict[tuple[str, str], bool]` in `triage_source`, `main` and `tm_conflicts`. `classify_pm(entries, literature, override_tm)` is called with `(mine, key in literature, literature.get(key, False))`. `CURATED_GPI_COLUMNS` has 12 names and `curated_row` in the tests uses the same keys. `CONFLICT_COLUMNS` and `UNMATCHED_COLUMNS` are used in Task 2 and Task 3 with the same names.

**Known limits.** `read_curated_gpi` checks that `pmid` is present, not that it resolves in PubMed. The reviewer does that (spec section 4). The new required fields (`evidence_note`, `reviewer`, `review_date`) make a draft row fail the reader. Draft rows belong in a separate file until a reviewer checks them.
