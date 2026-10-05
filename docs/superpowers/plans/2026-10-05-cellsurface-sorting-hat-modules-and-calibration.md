# cellsurface_sorting_hat module wrappers and calibration: Implementation Plan (Plan 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Tasks 1 to 10 are code with tests that run anywhere. Tasks 11 to 15 run on the HPCC (not in CI) and write measured numbers to a report.

**Goal:** Produce the module tables that the core engine reads (SignalP rule R0, Pfam families, two repeat detectors, allergen homology, antigen ranking lookup, Cys-rich tiers, spherule expression, TM helices), and record for each module where it was calibrated and what sensitivity and specificity were measured there.

**Architecture:** One small wrapper per tool turns the tool's output into a module table plus a run record (the contract of Plan 1). Heavy tools run as SLURM jobs from scripts in `scripts/sorting_hat/`; parsing and decisions run in Python and are tested on small fixtures. A calibration package computes sensitivity and specificity with cluster-bootstrap intervals, writes them into the `measure` object of the status sources, and refuses to call a module `estimated` unless both rates were measured on enough independent proteins.

**Tech Stack:** Python 3.11+ (`/usr/bin/python3.12`), numpy, PyYAML, pytest, ruff. Tools on HPCC: SignalP 6 GPU, hmmer 3.4, ncbi-blast 2.14.0+, TMHMM 2.0c. No new Python dependency.

**Spec:** `docs/superpowers/specs/2026-10-04-orchestrator-design.md` (revision 5). **Plan 1:** `docs/superpowers/plans/2026-10-04-cellsurface-sorting-hat-core.md` (must be implemented first; its package, `ModuleTable` contract and status sources are used here).

## Global Constraints

- The module table, run record and status source follow the "Module output contract" of Plan 1. Module tables are `.tsv.gz`; writes are atomic.
- A module's identity (`version`, `params_hash`, `artefact_hash`) belongs to the **tool, database and parameters**, never to one input file. A measurement of a rule must stay valid for another proteome run with the same tool version. (A SignalP result file is an input, not an artefact.)
- Status values: `estimated` needs at least 20 positives, at least 20 negatives and a 95% interval half-width of at most 0.10 for **both** sensitivity and specificity. A module that was never tested on negatives is at most `smoke`. A truth set that helped to set the rule or its cutoffs caps the status at `smoke` (`--leakage` is required and recorded).
- Intervals are cluster bootstrap (homology clusters), 2,000 resamples, fixed seed. Models are not refitted.
- Specificity is never inferred from absence in a database. A negative needs an independent reason to be a negative; if there is none, the module gets no specificity and the report says so.
- Tested taxa in a status entry are species level (resolved from `names.dmp`; a name with zero or several IDs stops the run). A protein's taxon may be a strain below the species.
- Antigen, Cys-rich and expression lookups apply only to the taxon IDs listed in `--applicable-taxa` (default: *C. immitis* RS, 246410, exact match, because the tables belong to one annotation).
- SLURM scripts: use `$SCRATCH` (`${SCRATCH:?}`), copy results back with a temporary name then `mv`, never `BASH_SOURCE`, never `/tmp`. Wait for jobs; do not poll in a tight loop.
- Use `/usr/bin/python3.12` for repository scripts. `pip`/`python` on the login node are Python 3.9.
- Compress large outputs (`.gz`); raw tool outputs under `$WORKDIR/raw/` are kept uncompressed only while they are small.
- Report and docs text: Simplified Technical English. Say when a number was not measured.

## What is calibrated, where, and what can be measured

This table is the calibration design. "Sn" is sensitivity, "Sp" is specificity. A cell that says *none* means no truth exists today. Numbers marked (measured) were produced from files in this repository on 2026-10-05.

| Module (call) | Calibration set and truth | Sn | Sp | Status ceiling now | What limits it |
|---|---|---|---|---|---|
| `step1_rule@R0` (`surface_glycoprotein`) | Phase C GO-evidence truth, `metrics.json`, pooled test sets. S1:all (*S. cerevisiae* + *C. albicans*): 232 positives, 4,244 negatives. Eurotiomycetes (*A. fumigatus* + *A. nidulans*): 128 and 208. Basidiomycota (*C. neoformans*, *U. maydis*): 16 and 60. | S1:all 0.603 [0.510, 0.683]; Euro 0.727 [0.640, 0.802]; Basid 0.938 [0.786, 1.000] (measured) | 0.963 [0.956, 0.971]; 0.990 [0.975, 1.000]; 0.917 [0.824, 0.984] (measured) | `estimated`, `estimated`, `smoke` | R0 has no fitted parameter, so no leakage. Negatives are GO-annotated non-surface proteins; secreted proteins that are not wall proteins are called by R0 by design (N-sec FPR 0.08 in S1:all). *Coccidioides* (Onygenales) is not covered: `unvalidated` there. |
| `adhesion_repeat` (`repeat02`, `repeat14`) | None with clade truth. Detectors were measured on a synthetic divergence series and SOWgp. A truth table of curated repeat adhesins and non-repeat surface proteins must be built (owner decision C1). | none | none | `unvalidated` | Calling known repeat adhesins is circular for detector settings that were set on SOWgp. |
| `adhesion_domain` (`pfam_adhesion`) | Per family: hits on the proteomes of 4 clades, non-member hits read by a person (Task 13). Members come from curated tables, not from Pfam itself. | per family, from curated members | per family, from reviewed non-member hits | `smoke` at best (Bys1 has 1 curated member) | Domain presence is not function: Pth11-like GPCRs carry CFEM (handled by `no_tm`). |
| `allergen_homology` (`allergen_homolog_hit`, `allergen_candidate`) | WHO/IUIS fungal sequences (116), leave-cluster-out against themselves (clusters at 40% identity, 70% coverage). | recovered share per cutoff (measured in Task 13) | **none**: no negatives exist | `smoke` (Sn only) | Absence from IUIS means "never tested", not "not an allergen". An identity cutoff cannot separate an allergen from a housekeeping paralog (Af293: 105 proteins hit at 35%/80 aa, 41 at 70% identity and 80% coverage, 30 at 95% or more; measured). |
| `antigen_lookup` (`antigen_candidate`) | Four *Coccidioides* anchors (PRA3, Ag2/PRA, SOWgp, PRA2) and the CF antigen as a cross-reactive control. | 4 of 4 inside the top 15%; PRA2 is at percentile 10.7 and the cut of 15 was chosen after the anchors were seen | none independent | `smoke`, `leakage = tuned_on_truth` | Four proteins. The ranking measures similarity to IEDB antigens and absence of orthologs in confounder fungi; it is not epitope prediction and not serology. An independent test needs IEDB assay data or serology (owner decision C3). |
| `cys_rich` | No accuracy claim (finder README). | n/a | n/a | evidence only | Not a category. |
| `expression` | n/a (RNA-seq log2 fold change, 2 replicates per condition). | n/a | n/a | evidence only | Not a category. |
| `tm` (TMHMM) | n/a here. | n/a | n/a | evidence only | TMHMM can read a signal peptide as a helix. |

Owner decisions that block more calibration (not part of the code tasks):
- **C1** Which proteins are truth for `adhesion_repeat` and `adhesion_domain`, and are the detector settings counted as tuned on them?
- **C2** Which proteins are negatives for the allergen module, if any? (Option: fungal proteins from families with a documented lack of IgE binding. None is listed today.)
- **C3** Which external data counts as independent truth for `antigen_candidate` (IEDB assay positives, serology, a held-out species)?
- **C4** Who signs off that a Pfam family is made `active` in `data/sorting_hat/family_table.tsv` after its specificity test.

## Review Focus

Input classes and failure modes this plan pins with tests.

1. A result file whose IDs do not match the FASTA (SignalP keeps the whole header; other tools cut at a space) must stop the module with exit code 2 and write nothing. (Task 9, `test_a_result_file_whose_ids_do_not_match_the_fasta_is_refused`)
2. A measurement of rule R0 must still apply to a different proteome run with the same SignalP version. (Task 9, `test_the_signalp_module_identity_does_not_depend_on_the_proteome`)
3. A rule never tested on negatives, or a truth set that overlaps the tuning data, must not be called `estimated`. (Task 6, `test_status_rule_follows_phase_c_and_needs_negatives_for_an_estimate`; Task 9, `test_leakage_other_than_none_caps_an_estimate_at_smoke`)
4. A CFEM domain in a transmembrane receptor must not count as an adhesion family hit. (Task 2, `test_no_tm_condition_drops_a_domain_in_a_protein_with_transmembrane_helices`; Task 9, `test_pfam_command_applies_the_no_tm_condition_from_the_tm_table`)
5. A protein with no BLAST hit has identity 0 (state `ok`), not unknown; a protein of another taxon has `not_applicable` for the *Coccidioides* lookups, never a made-up value. (Task 4, `test_blast_best_hit_is_by_bit_score_and_fields_are_derived`; Task 5, `test_antigen_states_and_fields`)
6. An allergen set cannot be calibrated against itself without removing close relatives; recall must be computed leave-cluster-out. (Task 4, `test_leave_cluster_out_recall_ignores_hits_inside_the_cluster`)
7. Two strains of one species: calls for identical sequences must be compared by sequence hash, not by ID. (Task 7, `test_only_identical_sequences_are_compared_by_sha256_not_by_id`)

---

### Task 0: Branch

**Files:** none.

- [ ] **Step 1: Create the working branch from the merged Plan 1 implementation**

```bash
git fetch origin
git checkout -b sorting-hat-modules origin/main      # Plan 1 (sorting-hat-core) must be merged first
git branch --show-current                            # expected: sorting-hat-modules
ls src/cellsurface_sorting_hat/engine.py             # expected: the file exists (Plan 1)
mkdir -p src/cellsurface_sorting_hat/modules src/cellsurface_sorting_hat/calibration
mkdir -p tests/cellsurface_sorting_hat/modules tests/cellsurface_sorting_hat/calibration
mkdir -p data/sorting_hat scripts/sorting_hat
```

If Plan 1 is not merged, branch from the Plan 1 branch tip instead and say so in the pull request.

### Task 1: Module writer and the SignalP wrapper

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/__init__.py`, `src/cellsurface_sorting_hat/modules/base.py`, `src/cellsurface_sorting_hat/modules/signalp.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_base_signalp.py`

**Interfaces:**
- Consumes: `write_atomic` (Plan 1 cache), `NA_INVALID`, `Protein` (Plan 1 fasta).
- Produces: `sha256_file`, `params_hash`, `artefact_hash`, `ModuleSpec(name, version, params, artefacts, tools, artefact_digest)`, `invalid_row(protein)`, `write_module(workdir, spec, columns, rows, run_state='ok', note='') -> run record dict`.
- Produces: `parse_signalp(path) -> {id: {prediction, sp_prob}}`, `signalp_rows(proteins, parsed)`, `COLUMNS`, `SignalPFormatError`. Only prediction `SP` is a call (rule R0); `LIPO`, `TAT` and the others are `not_called`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_base_signalp.py` with exactly this content:

```python
import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.base import (
    ModuleSpec,
    artefact_hash,
    params_hash,
    sha256_file,
    write_module,
)
from cellsurface_sorting_hat.modules.signalp import (
    COLUMNS,
    SignalPFormatError,
    parse_signalp,
    signalp_rows,
)

SIGNALP = (
    "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 20260930100126\n"
    "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
    "CIMG_00013-t26_1-p1 | transcript=CIMG_00013-t26_1 | gene=CIMG_00013 | organism=Coccidioides_immitis_RS"
    "\tOTHER\t1.000000\t0.000000\t\n"
    "CIMG_04613-t26_1-p1 | gene=CIMG_04613 | protein_length=324\tSP\t0.000100\t0.999800\tCS pos: 19-20. Pr: 0.9\n"
    "XP_3\tLIPO\t0.1\t0.2\tCS pos: 20-21. Pr: 0.8\n"
)


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


def test_parse_signalp_uses_the_first_token_of_the_header_as_id(tmp_path):
    path = tmp_path / "prediction_results.txt"
    path.write_text(SIGNALP)
    got = parse_signalp(path)
    assert set(got) == {"CIMG_00013-t26_1-p1", "CIMG_04613-t26_1-p1", "XP_3"}
    assert got["CIMG_04613-t26_1-p1"] == {"prediction": "SP", "sp_prob": 0.9998}


def test_only_sp_is_a_call_lipo_is_not(tmp_path):
    path = tmp_path / "p.txt"
    path.write_text(SIGNALP)
    rows = signalp_rows(
        [prot("CIMG_04613-t26_1-p1"), prot("XP_3"), prot("CIMG_00013-t26_1-p1")],
        parse_signalp(path),
    )
    assert [r["call"] for r in rows] == ["called", "not_called", "not_called"]
    assert rows[1]["prediction"] == "LIPO"


def test_missing_and_invalid_proteins_are_not_ok(tmp_path):
    path = tmp_path / "p.txt"
    path.write_text(SIGNALP)
    rows = signalp_rows([prot("ABSENT"), prot("BAD", state="na_invalid")], parse_signalp(path))
    assert [(r["id"], r["state"]) for r in rows] == [("ABSENT", "error"), ("BAD", "na_invalid")]


@pytest.mark.parametrize(
    "text,message",
    [
        ("", "no prediction rows"),
        ("A\tSP\t0.1\n", "at least 4"),
        ("A\tSP\t0.1\tx\n", "not a number"),
        ("A\tSP\t0.1\t0.9\nA\tSP\t0.1\t0.9\n", "duplicate ID"),
    ],
)
def test_bad_signalp_files_are_refused(tmp_path, text, message):
    path = tmp_path / "p.txt"
    path.write_text(text)
    with pytest.raises(SignalPFormatError, match=message):
        parse_signalp(path)


def test_hashes(tmp_path):
    a, b = tmp_path / "a.bin", tmp_path / "b.bin"
    a.write_bytes(b"1")
    b.write_bytes(b"2")
    assert sha256_file(a) != sha256_file(b)
    assert params_hash({"x": 1, "y": 2}) == params_hash({"y": 2, "x": 1})
    assert params_hash({"x": 1}) != params_hash({"x": 2})
    assert artefact_hash([a, b]) == artefact_hash([b, a])  # order does not matter
    assert artefact_hash([a]) != artefact_hash([b])


def test_write_module_writes_the_table_and_the_run_record(tmp_path):
    db = tmp_path / "db.hmm"
    db.write_bytes(b"hmm")
    spec = ModuleSpec("step1_rule@R0", "1", {"mode": "fast"}, (db,), {"signalp": "6.0h"})
    rows = [
        {"id": "A", "state": "ok", "call": "called", "sp_prob": "0.9", "prediction": "SP"},
        {"id": "B", "state": "na_invalid"},
    ]
    rec = write_module(tmp_path, spec, COLUMNS, rows)
    with gzip.open(tmp_path / "modules" / "step1_rule@R0.tsv.gz", "rt") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    assert table[0]["call"] == "called" and table[1]["call"] == ""
    on_disk = json.loads((tmp_path / "modules" / "step1_rule@R0.json").read_text())
    assert on_disk == rec
    assert rec["params_hash"] == params_hash({"mode": "fast"})
    assert rec["artefact_hash"] == artefact_hash([db]) and rec["n_rows"] == 2


def test_write_module_refuses_a_duplicate_id(tmp_path):
    with pytest.raises(ValueError, match="duplicate row"):
        write_module(
            tmp_path,
            ModuleSpec("m", "1"),
            [],
            [{"id": "A", "state": "ok"}, {"id": "A", "state": "ok"}],
        )
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_base_signalp.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.modules`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/__init__.py` with exactly this content:

```python
"""Module wrappers: turn tool outputs into the module tables that the engine reads."""
```

Create `src/cellsurface_sorting_hat/modules/base.py` with exactly this content:

```python
"""Shared helpers: hashes and the writer for a module table plus its run record."""

import csv
import gzip
import hashlib
import io
import json
from dataclasses import dataclass, field
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.fasta import NA_INVALID


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def params_hash(params):
    """Hash of a parameter dict (canonical JSON)."""
    return hashlib.sha256(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()


def artefact_hash(paths):
    """Hash of the databases or models a module uses: sorted ``name:sha256`` lines."""
    lines = sorted(f"{Path(p).name}:{sha256_file(p)}" for p in paths)
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


@dataclass(frozen=True)
class ModuleSpec:
    name: str
    version: str
    params: dict = field(default_factory=dict)
    artefacts: tuple = ()
    tools: dict = field(default_factory=dict)
    artefact_digest: str = (
        ""  # use instead of hashing large files (for example a recorded Pfam sha256)
    )


def invalid_row(protein):
    return {"id": protein.id, "state": NA_INVALID}


def write_module(workdir, spec, columns, rows, run_state="ok", note=""):
    """Write ``<workdir>/modules/<name>.tsv.gz`` and ``<name>.json`` (see the module contract).

    ``columns`` are the fields after ``id`` and ``state``. A row that lacks a column gets an empty
    value. Every row needs ``id`` and ``state``.
    """
    folder = Path(workdir) / "modules"
    header = ["id", "state"] + list(columns)
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    seen = set()
    for r in rows:
        if r["id"] in seen:
            raise ValueError(f"{spec.name}: duplicate row for ID {r['id']!r}")
        seen.add(r["id"])
        writer.writerow([r.get(h, "") for h in header])
    write_atomic(folder / f"{spec.name}.tsv.gz", gzip.compress(buf.getvalue().encode(), mtime=0))
    record = {
        "module": spec.name,
        "version": spec.version,
        "params_hash": params_hash(spec.params),
        "artefact_hash": spec.artefact_digest or artefact_hash(spec.artefacts),
        "run_state": run_state,
        "params": spec.params,
        "tools": spec.tools,
        "n_rows": len(seen),
        "note": note,
    }
    write_atomic(
        folder / f"{spec.name}.json", (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
    )
    return record
```

Create `src/cellsurface_sorting_hat/modules/signalp.py` with exactly this content:

```python
"""SignalP 6 results -> module ``step1_rule@R0`` (rule R0: SignalP calls a signal peptide, ``SP``)."""

from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

MODULE = "step1_rule@R0"


class SignalPFormatError(ValueError):
    """The SignalP file does not have the expected shape."""


def parse_signalp(path):
    """Return ``{protein ID: {"prediction": str, "sp_prob": float}}``.

    ``prediction_results.txt`` has two ``#`` header lines, then one TAB separated row per protein:
    ID (the FASTA header, the ID is its first token), prediction, OTHER probability, SP(Sec/SPI)
    probability, cleavage site. Prediction ``SP`` is the Sec/SPI signal peptide that rule R0 uses.
    """
    out = {}
    for n, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise SignalPFormatError(f"{path}:{n}: expected at least 4 TAB separated fields")
        pid = fields[0].split()[0]
        try:
            sp_prob = float(fields[3])
        except ValueError:
            raise SignalPFormatError(f"{path}:{n}: SP probability is not a number") from None
        if pid in out:
            raise SignalPFormatError(f"{path}:{n}: duplicate ID {pid!r}")
        out[pid] = {"prediction": fields[1], "sp_prob": sp_prob}
    if not out:
        raise SignalPFormatError(f"{path}: no prediction rows")
    return out


def signalp_rows(proteins, parsed):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed:
            rows.append({"id": p.id, "state": "error"})
        else:
            hit = parsed[p.id]
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "called" if hit["prediction"] == "SP" else "not_called",
                    "sp_prob": f"{hit['sp_prob']:.6f}",
                    "prediction": hit["prediction"],
                }
            )
    return rows


COLUMNS = ["call", "sp_prob", "prediction"]
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_base_signalp.py -q`
Expected: `10 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules tests/cellsurface_sorting_hat/modules/test_base_signalp.py
git commit -m "feat(sorting-hat): module writer and SignalP R0 wrapper

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 2: Pfam family table, domain table parser and the specificity report

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/pfam.py`, `data/sorting_hat/family_table.tsv`
- Test: `tests/cellsurface_sorting_hat/modules/test_pfam.py`

**Interfaces:**
- Consumes: `invalid_row` (Task 1).
- Produces: `load_family_table(path) -> list[Family]`, `parse_domtblout(path)`, `pfam_rows(proteins, hits, families, module, sp_calls=None, tm_counts=None)`, `specificity_report(hit_ids, member_ids, universe_ids)`, `FamilyTableError`, `FAMILY_COLUMNS`, `MODULES`, `COLUMNS`.
- A family counts only if `active = yes` (with `active_by` and `active_date`). `second_condition` is empty, `signal_peptide` or `no_tm`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_pfam.py` with exactly this content:

```python
import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.pfam import (
    FAMILY_COLUMNS,
    FamilyTableError,
    load_family_table,
    parse_domtblout,
    pfam_rows,
    specificity_report,
)

DOMTBL = (
    "# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description\n"
    "P1 - 289 CFEM PF05730.17 70 1e-20 60.1 8.9 1 1 1e-21 2e-20 59.0 8.9 1 70 20 90 20 91 0.9 -\n"
    "P2 - 400 Hydrophobin PF01185.24 60 1e-12 40.0 0.0 1 1 1e-13 3e-12 39.0 0.0 1 60 5 65 5 66 0.9 -\n"
    "P3 - 300 AltA1 PF16541.11 150 1e-30 90.0 0.0 1 1 1e-31 1e-30 89.0 0.0 1 150 10 160 10 161 0.9 -\n"
    "P4 - 500 Asp PF00026.29 300 1e-40 120.0 0.0 1 1 1e-41 1e-40 119.0 0.0 1 300 10 310 10 311 0.9 -\n"
)


def write_table(path, rows):
    lines = ["\t".join(FAMILY_COLUMNS)]
    for r in rows:
        lines.append("\t".join(r.get(c, "") for c in FAMILY_COLUMNS))
    path.write_text("\n".join(lines) + "\n")


def fam(acc, module="pfam_adhesion", active="yes", second="", **extra):
    base = {
        "pfam_acc": acc,
        "name": acc,
        "class": "x",
        "module": module,
        "source_pmid": "1",
        "pfam_release": "38.2",
        "specificity_note": "n",
        "second_condition": second,
        "active": active,
        "active_by": "owner" if active == "yes" else "",
        "active_date": "2026-10-05" if active == "yes" else "",
    }
    base.update(extra)
    return base


def prot(pid):
    return Protein(pid, "MKT", "x" * 64, "ok", "", 0.0)


def test_domtblout_is_parsed_without_the_accession_version(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text(DOMTBL)
    hits = parse_domtblout(path)
    assert [(h["target"], h["acc"]) for h in hits] == [
        ("P1", "PF05730"),
        ("P2", "PF01185"),
        ("P3", "PF16541"),
        ("P4", "PF00026"),
    ]
    assert hits[0]["ievalue"] == 2e-20


def test_only_active_families_of_the_module_count(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(
        path, [fam("PF05730"), fam("PF01185", active="no"), fam("PF16541", module="pfam_allergen")]
    )
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    prots = [prot(p) for p in ("P1", "P2", "P3", "P4")]
    adh = {r["id"]: r["hit"] for r in pfam_rows(prots, hits, fams, "pfam_adhesion")}
    all_ = {r["id"]: r["hit"] for r in pfam_rows(prots, hits, fams, "pfam_allergen")}
    assert adh == {"P1": "1", "P2": "0", "P3": "0", "P4": "0"}  # PF01185 is inactive
    assert all_ == {"P1": "0", "P2": "0", "P3": "1", "P4": "0"}


def test_second_condition_needs_a_signal_peptide_call(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF00026", second="signal_peptide")])
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    p4 = [prot("P4")]
    assert pfam_rows(p4, hits, fams, "pfam_adhesion")[0]["hit"] == "0"  # no SP information
    assert pfam_rows(p4, hits, fams, "pfam_adhesion", {"P4": "not_called"})[0]["hit"] == "0"
    assert pfam_rows(p4, hits, fams, "pfam_adhesion", {"P4": "called"})[0]["hit"] == "1"


@pytest.mark.parametrize(
    "mutate,message",
    [
        (lambda r: r.update(pfam_acc="PF05730.17"), "no version"),
        (lambda r: r.update(module="other"), "module must be"),
        (lambda r: r.update(active="maybe"), "active must be"),
        (lambda r: r.update(second_condition="tm"), "second_condition"),
        (lambda r: r.update(active="yes", active_by=""), "needs active_by"),
    ],
)
def test_bad_family_rows_are_refused(tmp_path, mutate, message):
    row = fam("PF05730")
    mutate(row)
    path = tmp_path / "f.tsv"
    write_table(path, [row])
    with pytest.raises(FamilyTableError, match=message):
        load_family_table(path)


def test_duplicate_families_are_refused(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF05730"), fam("PF05730")])
    with pytest.raises(FamilyTableError, match="duplicate"):
        load_family_table(path)


def test_specificity_report_counts_and_lists_the_non_member_hits():
    rep = specificity_report(
        hit_ids={"a", "b", "x"}, member_ids={"a", "b", "c"}, universe_ids=set("abcxyz")
    )
    assert (rep["tp"], rep["fp"], rep["fn"], rep["tn"]) == (2, 1, 1, 2)
    assert rep["nonmember_hits"] == ["x"] and rep["missed_members"] == ["c"]
    assert rep["sensitivity"] == pytest.approx(2 / 3) and rep["specificity"] == pytest.approx(2 / 3)


def test_specificity_report_needs_members_inside_the_universe():
    with pytest.raises(ValueError):
        specificity_report({"a"}, {"zz"}, {"a"})


def test_no_tm_condition_drops_a_domain_in_a_protein_with_transmembrane_helices(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF05730", second="no_tm")])
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    p1 = [prot("P1")]
    assert pfam_rows(p1, hits, fams, "pfam_adhesion")[0]["hit"] == "0"  # no TM information
    assert (
        pfam_rows(p1, hits, fams, "pfam_adhesion", tm_counts={"P1": 7})[0]["hit"] == "0"
    )  # a receptor
    assert pfam_rows(p1, hits, fams, "pfam_adhesion", tm_counts={"P1": 0})[0]["hit"] == "1"
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_pfam.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.modules.pfam`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/pfam.py` with exactly this content:

```python
"""Pfam family table and ``hmmsearch --domtblout`` results -> modules ``pfam_adhesion``, ``pfam_allergen``."""

import csv
from dataclasses import dataclass
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

FAMILY_COLUMNS = [
    "pfam_acc",
    "name",
    "class",
    "module",
    "source_pmid",
    "pfam_release",
    "specificity_note",
    "second_condition",
    "active",
    "active_by",
    "active_date",
]
MODULES = ("pfam_adhesion", "pfam_allergen")
SECOND_CONDITIONS = ("", "signal_peptide", "no_tm")
COLUMNS = ["hit", "families"]


class FamilyTableError(ValueError):
    """The family table is not valid."""


@dataclass(frozen=True)
class Family:
    pfam_acc: str
    name: str
    cls: str
    module: str
    second_condition: str
    active: bool


def load_family_table(path):
    rows, seen = [], set()
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(FAMILY_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise FamilyTableError(f"{path}: missing column(s) {sorted(missing)}")
        for n, r in enumerate(reader, 2):
            where = f"{path}:{n}"
            acc = r["pfam_acc"].strip()
            if not acc.startswith("PF") or "." in acc:
                raise FamilyTableError(f"{where}: pfam_acc must look like PF05730 (no version)")
            if acc in seen:
                raise FamilyTableError(f"{where}: duplicate {acc}")
            seen.add(acc)
            if r["module"] not in MODULES:
                raise FamilyTableError(f"{where}: module must be one of {MODULES}")
            if r["second_condition"] not in SECOND_CONDITIONS:
                raise FamilyTableError(
                    f"{where}: second_condition must be one of {SECOND_CONDITIONS}"
                )
            if r["active"] not in ("yes", "no"):
                raise FamilyTableError(f"{where}: active must be yes or no")
            if r["active"] == "yes" and not (r["active_by"].strip() and r["active_date"].strip()):
                raise FamilyTableError(f"{where}: an active family needs active_by and active_date")
            rows.append(
                Family(
                    acc,
                    r["name"],
                    r["class"],
                    r["module"],
                    r["second_condition"],
                    r["active"] == "yes",
                )
            )
    return rows


def parse_domtblout(path):
    """Return a list of ``{"target", "acc", "ievalue", "score"}`` (accession without version)."""
    hits = []
    for n, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 22:
            raise ValueError(f"{path}:{n}: expected at least 22 fields, found {len(f)}")
        hits.append(
            {
                "target": f[0],
                "acc": f[4].split(".")[0],
                "ievalue": float(f[12]),
                "score": float(f[13]),
            }
        )
    return hits


def pfam_rows(proteins, hits, families, module, sp_calls=None, tm_counts=None):
    """Rows for ``module``. A protein is a hit if it has a domain of an active family of this module.

    ``second_condition = signal_peptide`` counts a hit only when ``sp_calls[id]`` is ``called``.
    ``second_condition = no_tm`` counts a hit only when ``tm_counts[id]`` is 0 (for example a CFEM
    domain in a receptor with transmembrane helices is not a cell wall CFEM protein). If the needed
    table is None, such a family never counts.
    """
    wanted = {f.pfam_acc: f for f in families if f.module == module and f.active}
    found = {}
    for h in hits:
        fam = wanted.get(h["acc"])
        if fam is None:
            continue
        if (
            fam.second_condition == "signal_peptide"
            and (sp_calls or {}).get(h["target"]) != "called"
        ):
            continue
        if fam.second_condition == "no_tm" and (tm_counts or {}).get(h["target"]) != 0:
            continue
        found.setdefault(h["target"], set()).add(fam.pfam_acc)
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        accs = sorted(found.get(p.id, ()))
        rows.append(
            {"id": p.id, "state": "ok", "hit": "1" if accs else "0", "families": ",".join(accs)}
        )
    return rows


def specificity_report(hit_ids, member_ids, universe_ids):
    """Compare the proteins hit by a family with the known members (a specificity test).

    ``universe_ids`` are all proteins searched. Returns counts and the non-member hits, which a
    person must read before the family is made active.
    """
    hit, members, universe = set(hit_ids), set(member_ids), set(universe_ids)
    if not members <= universe:
        raise ValueError("member_ids must be a subset of universe_ids")
    tp, fp, fn = hit & members, hit - members, members - hit
    tn = universe - hit - members
    return {
        "tp": len(tp),
        "fp": len(fp),
        "fn": len(fn),
        "tn": len(tn),
        "sensitivity": len(tp) / len(members) if members else None,
        "specificity": len(tn) / (len(tn) + len(fp)) if (len(tn) + len(fp)) else None,
        "nonmember_hits": sorted(fp),
        "missed_members": sorted(fn),
    }
```

The family table is TAB separated (the first line is the header). All ten families start inactive. A family is made active only after its specificity test and a sign-off (Task 13).

Create `data/sorting_hat/family_table.tsv` with exactly this content:

```text
pfam_acc	name	class	module	source_pmid	pfam_release	specificity_note	second_condition	active	active_by	active_date
PF05730	CFEM	2b-i CFEM-domain surface protein	pfam_adhesion	docs/reports/2026-09-29-class2b-structure.md; PMID 28513415 names CFEM proteins	38.2	hemophore fold, not shown to bind a host receptor; Pth11-like GPCRs also carry a CFEM domain, so no_tm keeps only proteins without TMHMM helices; specificity test not run	no_tm	no		
PF04681	Bys1	2b-iii Bys1 invasin (CalA)	pfam_adhesion	docs/reports/2026-09-29-class2b-structure.md	38.2	Bys1 is not a thaumatin Pfam; calB and calC paralogs are the control; specificity test not run		no		
PF01185	Hydrophobin	2c hydrophobin	pfam_adhesion	PMID 28513415 names RodA and RodB	38.2	specificity test not run		no		
PF06766	Hydrophobin_2	2c hydrophobin	pfam_adhesion	PMID 28513415 names RodA and RodB	38.2	specificity test not run		no		
PF28987	DewD	2c hydrophobin (DewD)	pfam_adhesion	analysis/cys_candidates/01_known_family_hmm.sh	38.2	specificity test not run		no		
PF22354	Eas	2c hydrophobin (Eas)	pfam_adhesion	Pfam 38.2 name Eas, Hydrophobin	38.2	not in the earlier HMM script; specificity test not run		no		
PF11766	Candida_ALS_N	2a Als adhesin N-terminal domain	pfam_adhesion	PMID 28513415 names the Als family	38.2	Candida-specific; specificity test not run		no		
PF05792	Candida_ALS	2a Als adhesin repeat region	pfam_adhesion	PMID 28513415 names the Als family	38.2	Candida-specific; specificity test not run		no		
PF16541	AltA1	Alt a 1 allergen family	pfam_allergen	WHO/IUIS Alt a 1; docs/reports/2026-10-04-fungal-allergen-scoping.md	38.2	allergen-specific; specificity test not run		no		
PF25312	Allergen_Asp_f_4	Asp f 4 allergen family	pfam_allergen	WHO/IUIS Asp f 4; docs/reports/2026-10-04-fungal-allergen-scoping.md	38.2	allergen-specific; specificity test not run		no
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_pfam.py -q`
Expected: `12 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/pfam.py data/sorting_hat/family_table.tsv tests/cellsurface_sorting_hat/modules/test_pfam.py
git commit -m "feat(sorting-hat): Pfam family table and domain table parser

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 3: Repeat detector wrappers

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/repeats.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_repeats.py`

**Interfaces:**
- Consumes: `invalid_row` (Task 1).
- Produces: `parse_repeat_table(path) -> {id: {period, copies, coverage}}`, `repeat_rows(proteins, parsed, min_coverage=0.25, min_copies=2.5)`, `run_detector(repo_root, module, fasta, out_tsv, python=sys.executable, extra=())`, `COLUMNS`.
- The call is `period > 0 and coverage >= min_coverage and copies >= min_copies`, the same rule as `analysis/cocci_repeats/03_repeat_surface_candidates.py`. The detector scripts are not changed.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_repeats.py` with exactly this content:

```python
import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.repeats import parse_repeat_table, repeat_rows, run_detector


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


REPEATS = (
    "strain\tprotein\tlength\trep_period\trep_score\trep_start\trep_end\trep_n_copies\trep_coverage\n"
    "S\tA\t300\t47\t0.8\t80\t270\t4.0\t0.62\n"
    "S\tB\t300\t0\t0.1\t0\t0\t0\t0.0\n"
    "S\tC\t300\t12\t0.5\t10\t60\t4.1\t0.17\n"
    "S\tD\t300\t12\t0.5\t10\t60\t2.4\t0.40\n"
)


def test_repeat_call_needs_period_coverage_and_copies(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    rows = repeat_rows([prot(x) for x in "ABCD"] + [prot("E")], parse_repeat_table(path))
    assert [r.get("call") for r in rows] == [
        "called",
        "not_called",
        "not_called",
        "not_called",
        None,
    ]
    assert rows[-1]["state"] == "error"  # not in the detector output


def test_repeat_thresholds_are_parameters(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(REPEATS)
    rows = repeat_rows([prot("C")], parse_repeat_table(path), min_coverage=0.10)
    assert rows[0]["call"] == "called"


def test_repeat_table_needs_its_columns_and_unique_proteins(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text("protein\trep_period\nA\t3\n")
    with pytest.raises(ValueError, match="missing column"):
        parse_repeat_table(path)
    path.write_text("protein\trep_period\trep_n_copies\trep_coverage\nA\t3\t3\t0.5\nA\t3\t3\t0.5\n")
    with pytest.raises(ValueError, match="duplicate protein"):
        parse_repeat_table(path)


def test_run_detector_builds_the_command(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr("subprocess.run", lambda cmd, check: calls.append((cmd, check)))
    cmd = run_detector(
        tmp_path, "repeat14", "in.faa", "out.tsv", python="py", extra=["--mode", "exact"]
    )
    assert calls == [(cmd, True)]
    assert cmd[0] == "py" and cmd[1].endswith("analysis/cocci_repeats/14_repeat_detect_general.py")
    assert cmd[2:] == ["in.faa", "--out", "out.tsv", "--mode", "exact"]
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_repeats.py -q`
Expected: collection error, `ModuleNotFoundError: ... modules.repeats`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/repeats.py` with exactly this content:

```python
"""Repeat detector tables -> modules ``repeat02`` and ``repeat14`` (call = repeat per 03_repeat_surface_candidates.py)."""

import csv
import subprocess
import sys
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

DEFAULT_MIN_COVERAGE = 0.25
DEFAULT_MIN_COPIES = 2.5
SCRIPTS = {"repeat02": "02_repeat_profile.py", "repeat14": "14_repeat_detect_general.py"}


def parse_repeat_table(path):
    """Return ``{protein ID: {"period": int, "copies": float, "coverage": float}}``."""
    out = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {"protein", "rep_period", "rep_n_copies", "rep_coverage"}
        missing = need - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
        for n, r in enumerate(reader, 2):
            if r["protein"] in out:
                raise ValueError(f"{path}:{n}: duplicate protein {r['protein']!r}")
            out[r["protein"]] = {
                "period": int(float(r["rep_period"] or 0)),
                "copies": float(r["rep_n_copies"] or 0),
                "coverage": float(r["rep_coverage"] or 0),
            }
    return out


def repeat_rows(proteins, parsed, min_coverage=DEFAULT_MIN_COVERAGE, min_copies=DEFAULT_MIN_COPIES):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed:
            rows.append({"id": p.id, "state": "error"})
        else:
            r = parsed[p.id]
            called = r["period"] > 0 and r["coverage"] >= min_coverage and r["copies"] >= min_copies
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "called" if called else "not_called",
                    "period": r["period"],
                    "copies": f"{r['copies']:.2f}",
                    "coverage": f"{r['coverage']:.3f}",
                }
            )
    return rows


COLUMNS = ["call", "period", "copies", "coverage"]


def run_detector(repo_root, module, fasta, out_tsv, python=sys.executable, extra=()):
    """Run an existing detector script from ``analysis/cocci_repeats`` on one FASTA."""
    script = Path(repo_root) / "analysis" / "cocci_repeats" / SCRIPTS[module]
    cmd = [python, str(script), str(fasta), "--out", str(out_tsv), *extra]
    subprocess.run(cmd, check=True)
    return cmd
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_repeats.py -q`
Expected: `4 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/repeats.py tests/cellsurface_sorting_hat/modules/test_repeats.py
git commit -m "feat(sorting-hat): repeat detector wrappers

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 4: Allergen homology and the leave-cluster-out recall

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/allergen.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_allergen.py`

**Interfaces:**
- Consumes: `invalid_row` (Task 1).
- Produces: `build_allergen_fasta(isoallergen_tsv, out_fasta) -> int`, `parse_blast(path) -> {query: best hit}`, `allergen_rows(proteins, best)` (fields `identity`, `aligned_length`, `coverage`, `allergen_name`; no hit gives `0`, `0`, `0`, `''` with state `ok`), `parse_blast_hits`, `cluster_by_hits`, `best_hit_outside_cluster`, `recall_by_cutoff`, `lco_report(path, cutoffs, cluster_identity=40.0, cluster_coverage=70.0)`, `BLAST_FIELDS`, `COLUMNS`.
- The best hit is the highest bit score. Coverage is the alignment length as a percent of the allergen sequence (capped at 100).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_allergen.py` with exactly this content:

```python
import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.allergen import (
    allergen_name,
    allergen_rows,
    best_hit_outside_cluster,
    build_allergen_fasta,
    cluster_by_hits,
    parse_blast,
    recall_by_cutoff,
)


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


IUIS = (
    "AllergenID\tIsoName\tName\tSequence\n"
    "11\tAsp f 1.0101\tAsp f 1\tMKTAYIAKQR\n"
    "12\tAsp f 2.0101\tAsp f 2\tMNLLPQWERT\n"
    "13\tNoSeq\tNo seq\t\n"
)


def test_allergen_fasta_has_clean_ids_and_skips_rows_without_sequence(tmp_path):
    src, out = tmp_path / "i.tsv", tmp_path / "a.faa"
    src.write_text(IUIS)
    assert build_allergen_fasta(src, out) == 2
    assert out.read_text() == ">Asp_f_1.0101|11\nMKTAYIAKQR\n>Asp_f_2.0101|12\nMNLLPQWERT\n"


BLAST = (
    "P1\tAsp_f_1.0101|11\t99.0\t100\t120\t125\t200\t1e-50\n"
    "P1\tAsp_f_2.0101|12\t45.0\t90\t120\t300\t60\t1e-5\n"
    "P2\tAsp_f_2.0101|12\t38.0\t85\t200\t300\t50\t1e-3\n"
)


def test_blast_best_hit_is_by_bit_score_and_fields_are_derived(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text(BLAST)
    best = parse_blast(path)
    rows = {
        r["id"]: r
        for r in allergen_rows(
            [prot("P1"), prot("P2"), prot("P3"), prot("BAD", "na_invalid")], best
        )
    }
    assert rows["P1"]["identity"] == "99.00" and rows["P1"]["allergen_name"] == "Asp_f_1.0101"
    assert rows["P1"]["coverage"] == "80.0"  # 100 of 125 residues
    assert (rows["P2"]["identity"], rows["P2"]["aligned_length"], rows["P2"]["coverage"]) == (
        "38.00",
        "85",
        "28.3",
    )
    assert (rows["P3"]["identity"], rows["P3"]["aligned_length"], rows["P3"]["allergen_name"]) == (
        "0",
        "0",
        "",
    )
    assert rows["BAD"]["state"] == "na_invalid"


def test_blast_line_with_the_wrong_field_count_is_refused(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text("P1\tS\t99.0\n")
    with pytest.raises(ValueError, match="expected 8 fields"):
        parse_blast(path)


def test_allergen_name():
    assert allergen_name("Asp_f_1.0101|11") == "Asp_f_1.0101"


def h(q, s, ident, length, qlen=100, slen=100, bits=100.0):
    return {
        "query": q,
        "subject": s,
        "identity": ident,
        "length": length,
        "qlen": qlen,
        "slen": slen,
        "bitscore": bits,
    }


def test_clusters_join_close_relatives_only():
    hits = [
        h("a", "b", 90, 95),
        h("b", "c", 30, 95),
        h("c", "d", 50, 40),
    ]  # b-c too distant; c-d covers 40% only
    cl = cluster_by_hits(["a", "b", "c", "d"], hits)
    assert cl["a"] == cl["b"] and len({cl["a"], cl["c"], cl["d"]}) == 3


def test_leave_cluster_out_recall_ignores_hits_inside_the_cluster():
    ids = ["a", "b", "c"]
    cl = {"a": "C1", "b": "C1", "c": "C2"}
    hits = [
        h("a", "b", 99, 100, bits=300),
        h("a", "c", 45, 80, bits=60),
        h("c", "a", 45, 80, bits=60),
        h("b", "a", 99, 100, bits=300),
    ]
    best = best_hit_outside_cluster(hits, cl)
    assert (
        set(best) == {"a", "c"} and best["a"]["subject"] == "c"
    )  # b has no hit outside its cluster
    rec = recall_by_cutoff(ids, best, [(35, 50), (70, 80)])
    assert [(r["min_identity"], r["recovered"], r["n"]) for r in rec] == [(35, 2, 3), (70, 0, 3)]
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_allergen.py -q`
Expected: collection error, `ModuleNotFoundError: ... modules.allergen`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/allergen.py` with exactly this content:

```python
"""WHO/IUIS fungal allergen sequences, BLASTP results and the leave-cluster-out recall.

``allergen_homology`` module fields: ``identity`` (percent identity of the best local alignment),
``aligned_length`` (alignment length in residues), ``coverage`` (percent of the allergen sequence
covered by the alignment) and ``allergen_name``. A protein with no hit has ``0``, ``0``, ``0`` and an
empty name (state ``ok``). The 35%/80 aa rule (``allergen_homolog_hit``) and the stricter candidate
call (>= 70% identity and >= 80% coverage) are made by the engine from these fields.
"""

import csv
import re
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

BLAST_FIELDS = "qseqid sseqid pident length qlen slen bitscore evalue"
COLUMNS = ["identity", "aligned_length", "coverage", "allergen_name"]
_ID_BAD = re.compile(r"[^A-Za-z0-9_.\-]+")


def build_allergen_fasta(isoallergen_tsv, out_fasta):
    """Write the IUIS fungal sequences as FASTA. ID = ``<IsoAllergenName>|<AllergenID>``.

    Rows without a sequence are skipped. Returns the number of sequences written.
    """
    n = 0
    with open(isoallergen_tsv, encoding="utf-8-sig", newline="") as fh, open(out_fasta, "w") as out:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq = (r.get("Sequence") or "").strip().upper()
            if not seq or seq == "NAN":
                continue
            name = _ID_BAD.sub("_", (r.get("IsoName") or r.get("Name") or "").strip())
            out.write(f">{name}|{r['AllergenID']}\n{seq}\n")
            n += 1
    return n


def parse_blast(path):
    """Return the best local alignment per query: ``{query: {...}}`` (highest bit score)."""
    best = {}
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) != len(BLAST_FIELDS.split()):
            raise ValueError(
                f"{path}:{n}: expected {len(BLAST_FIELDS.split())} fields, found {len(f)}"
            )
        q, s, pident, length, qlen, slen, bits, ev = f
        hit = {
            "subject": s,
            "identity": float(pident),
            "length": int(length),
            "qlen": int(qlen),
            "slen": int(slen),
            "bitscore": float(bits),
            "evalue": float(ev),
        }
        cur = best.get(q)
        if cur is None or (hit["bitscore"], hit["identity"]) > (cur["bitscore"], cur["identity"]):
            best[q] = hit
    return best


def allergen_name(subject):
    """``Asp_f_1.0101|1234`` -> ``Asp_f_1.0101``."""
    return subject.split("|", 1)[0]


def allergen_rows(proteins, best):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        h = best.get(p.id)
        if h is None:
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "identity": "0",
                    "aligned_length": "0",
                    "coverage": "0",
                    "allergen_name": "",
                }
            )
        else:
            cov = min(100.0, 100.0 * h["length"] / h["slen"])
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "identity": f"{h['identity']:.2f}",
                    "aligned_length": str(h["length"]),
                    "coverage": f"{cov:.1f}",
                    "allergen_name": allergen_name(h["subject"]),
                }
            )
    return rows


def cluster_by_hits(ids, hits, min_identity=40.0, min_coverage=70.0):
    """Single-linkage clusters of allergens from an all-against-all BLAST table.

    ``hits`` is a list of dicts with ``query``, ``subject``, ``identity``, ``length``, ``qlen``,
    ``slen``. Two sequences join when the alignment reaches ``min_identity`` percent identity and
    covers ``min_coverage`` percent of the shorter sequence. Returns ``{id: cluster id}``.
    """
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for h in hits:
        q, s = h["query"], h["subject"]
        if q == s or q not in parent or s not in parent:
            continue
        short = min(h["qlen"], h["slen"])
        if h["identity"] >= min_identity and 100.0 * h["length"] / short >= min_coverage:
            parent[find(q)] = find(s)
    roots = {}
    return {i: roots.setdefault(find(i), f"C{len(roots) + 1:04d}") for i in sorted(ids)}


def best_hit_outside_cluster(hits, cluster_of):
    """For every query, the best hit (bit score) to a subject in a different cluster."""
    best = {}
    for h in hits:
        q, s = h["query"], h["subject"]
        if q == s or cluster_of.get(q) == cluster_of.get(s):
            continue
        cur = best.get(q)
        if cur is None or (h["bitscore"], h["identity"]) > (cur["bitscore"], cur["identity"]):
            best[q] = h
    return best


def recall_by_cutoff(all_ids, best_outside, cutoffs):
    """Share of allergens whose best hit in another cluster reaches each (identity, coverage) cutoff.

    ``cutoffs`` is a list of ``(min_identity, min_coverage)`` pairs. This estimates how many
    allergens a homology rule finds when no close relative is in the reference set (sensitivity
    for new families). It says nothing about false calls.
    """
    out = []
    for min_id, min_cov in cutoffs:
        k = 0
        for i in all_ids:
            h = best_outside.get(i)
            if h and h["identity"] >= min_id and 100.0 * h["length"] / h["slen"] >= min_cov:
                k += 1
        out.append(
            {"min_identity": min_id, "min_coverage": min_cov, "recovered": k, "n": len(all_ids)}
        )
    return out


def parse_blast_hits(path):
    """All rows of a BLAST table as dicts with the keys that ``cluster_by_hits`` and
    ``best_hit_outside_cluster`` use."""
    hits = []
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) != len(BLAST_FIELDS.split()):
            raise ValueError(
                f"{path}:{n}: expected {len(BLAST_FIELDS.split())} fields, found {len(f)}"
            )
        hits.append(
            {
                "query": f[0],
                "subject": f[1],
                "identity": float(f[2]),
                "length": int(f[3]),
                "qlen": int(f[4]),
                "slen": int(f[5]),
                "bitscore": float(f[6]),
                "evalue": float(f[7]),
            }
        )
    return hits


def lco_report(path, cutoffs, cluster_identity=40.0, cluster_coverage=70.0):
    """Leave-cluster-out recall of the allergen set against itself (all-against-all BLAST table)."""
    hits = parse_blast_hits(path)
    ids = sorted({h["query"] for h in hits})
    clusters = cluster_by_hits(ids, hits, cluster_identity, cluster_coverage)
    best = best_hit_outside_cluster(hits, clusters)
    return {
        "n_sequences": len(ids),
        "n_clusters": len(set(clusters.values())),
        "recall": recall_by_cutoff(ids, best, cutoffs),
    }
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_allergen.py -q`
Expected: `6 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/allergen.py tests/cellsurface_sorting_hat/modules/test_allergen.py
git commit -m "feat(sorting-hat): allergen homology module and leave-cluster-out recall

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 5: Lookups: antigen ranking, Cys-rich tiers, expression, TM helices

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/lookups.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_lookups.py`

**Interfaces:**
- Consumes: `invalid_row` (Task 1).
- Produces: `read_table(path, key)`, `load_protein_map(path)`, `ranking_by_gene(path)`, `gene_of_ranking_id`, `antigen_rows(proteins, taxa, protein_map, by_gene, applicable_taxa)`, `cys_rows`, `expression_rows`, `tm_rows`, and the `*_COLUMNS` lists.
- States: `not_applicable` (taxon not in `applicable_taxa`), `not_in_reference` (no map entry or table row), `ok`. For *C. immitis* RS the RefSeq protein ID (`XP_...`) maps to a gene (`CIMG_...`) through `protein_map.tsv`; several ranking rows per gene keep the lowest rank.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_lookups.py` with exactly this content:

```python
import gzip

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.lookups import (
    antigen_rows,
    cys_rows,
    expression_rows,
    gene_of_ranking_id,
    load_protein_map,
    ranking_by_gene,
    read_table,
    tm_rows,
)

RANKING = (
    "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
    "CIMG_04613-t26_1-p1\t651\t7.12\t2.5\t0.0\t0.9201\t0.0\n"
    "CIMG_04613-t26_2-p1\t900\t9.9\t2.0\t0.0\t0.9201\t0.0\n"
    "CIMG_09560-t26_1-p1\t986\t10.79\t0.1\t1.0\t0.9877\t69.4\n"
)
PMAP = "protein_id\tgene_id\tproduct\tlength\nXP_1\tCIMG_04613\tp\t324\nXP_2\tCIMG_09560\tp\t100\nXP_9\tCIMG_99999\tp\t10\n"


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


def test_gene_of_ranking_id():
    assert gene_of_ranking_id("CIMG_04613-t26_1-p1") == "CIMG_04613"


def test_ranking_keeps_the_best_row_per_gene_and_counts_genes_with_several(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(RANKING)
    by_gene, several = ranking_by_gene(path)
    assert by_gene["CIMG_04613"]["rank"] == "651" and several == 1


def test_antigen_states_and_fields(tmp_path):
    (tmp_path / "r.tsv").write_text(RANKING)
    (tmp_path / "m.tsv").write_text(PMAP)
    by_gene, _ = ranking_by_gene(tmp_path / "r.tsv")
    pmap = load_protein_map(tmp_path / "m.tsv")
    taxa = {
        "XP_1": 246410,
        "XP_2": 246410,
        "XP_9": 246410,
        "XP_X": 246410,
        "AF1": 746128,
        "BAD": 246410,
    }
    rows = antigen_rows(
        [prot(p, "na_invalid" if p == "BAD" else "ok") for p in taxa], taxa, pmap, by_gene, {246410}
    )
    got = {r["id"]: r for r in rows}
    assert (
        got["XP_1"]["state"] == "ok"
        and got["XP_1"]["percentile"] == "7.12"
        and got["XP_1"]["max_crossreact"] == "0.0"
    )
    assert got["XP_2"]["percentile"] == "10.79"
    assert got["XP_9"]["state"] == "not_in_reference"  # gene not in the ranking
    assert got["XP_X"]["state"] == "not_in_reference"  # protein not in the map
    assert got["AF1"]["state"] == "not_applicable"  # another taxon
    assert got["BAD"]["state"] == "na_invalid"


def test_read_table_counts_repeated_keys(tmp_path):
    path = tmp_path / "t.tsv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("k\tv\na\t1\na\t2\nb\t3\n")
    rows, repeated = read_table(path, "k")
    assert rows["a"]["v"] == "1" and repeated == 1
    with pytest.raises(ValueError, match="missing column"):
        read_table(path, "nope")


def test_cys_expression_and_tm_rows():
    taxa = {"A": 246410, "B": 246410, "C": 746128}
    cys = cys_rows(
        [prot(p) for p in "ABC"],
        taxa,
        {"A": {"tier": "cys_rich_sp_unassigned", "cys_frac": "0.1"}},
        {246410},
    )
    assert [(r["id"], r["state"]) for r in cys] == [
        ("A", "ok"),
        ("B", "not_in_reference"),
        ("C", "not_applicable"),
    ]
    assert cys[0]["tier"] == "cys_rich_sp_unassigned"
    expr = expression_rows(
        [prot("A"), prot("C")],
        taxa,
        {"A": "CIMG_1"},
        {"CIMG_1": {"log2fc_48h": "9.09", "padj_48h": "1e-5", "log2fc_8d": "8.0"}},
        {246410},
    )
    assert expr[0]["log2fc"] == "9.09" and expr[1]["state"] == "not_applicable"
    tm = tm_rows([prot("A"), prot("B")], {"A": {"pred_hel": "7", "topology": "o10-32i"}})
    assert tm[0]["n_tm"] == "7" and tm[1]["state"] == "error"
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_lookups.py -q`
Expected: collection error, `ModuleNotFoundError: ... modules.lookups`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/lookups.py` with exactly this content:

```python
"""Lookups into precomputed tables: antigen ranking, Cys-rich tiers, spherule expression, TM helices.

These modules need an ID map. For *C. immitis* RS the RefSeq protein ID (``XP_...``) maps to a
gene ID (``CIMG_...``) through ``protein_map.tsv``. The antigen ranking and the spherule table are
keyed on that gene ID. A protein with no map entry or no table row is ``not_in_reference``. A protein
from a taxon where the table does not apply is ``not_applicable``.
"""

import csv
import gzip
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

ANTIGEN_COLUMNS = [
    "percentile",
    "antigenicity",
    "specificity",
    "prevalence",
    "max_crossreact",
    "rank",
]
CYS_COLUMNS = ["tier", "cys_frac"]
EXPRESSION_COLUMNS = ["log2fc", "padj", "log2fc_8d"]
TM_COLUMNS = ["n_tm", "topology"]


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def read_table(path, key, delimiter="\t"):
    """Read a TSV into ``{row[key]: row}``. A repeated key keeps the first row and is counted."""
    rows, repeated = {}, 0
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if key not in (reader.fieldnames or []):
            raise ValueError(f"{path}: missing column {key!r}")
        for r in reader:
            if r[key] in rows:
                repeated += 1
            else:
                rows[r[key]] = r
    return rows, repeated


def load_protein_map(path):
    rows, _ = read_table(path, "protein_id")
    return {pid: r["gene_id"] for pid, r in rows.items()}


def gene_of_ranking_id(protein):
    """``CIMG_04613-t26_1-p1`` -> ``CIMG_04613``."""
    return protein.split("-", 1)[0]


def ranking_by_gene(ranking_tsv):
    """Best (lowest rank) ranking row per gene. Returns ``(rows, n_genes_with_several_rows)``."""
    rows, several = {}, set()
    with _open_text(ranking_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene = gene_of_ranking_id(r["protein"])
            cur = rows.get(gene)
            if cur is not None:
                several.add(gene)
            if cur is None or int(r["rank"]) < int(cur["rank"]):
                rows[gene] = r
    return rows, len(several)


def _applicable(taxa, pid, applicable_taxa):
    return taxa[pid] in applicable_taxa


def antigen_rows(proteins, taxa, protein_map, by_gene, applicable_taxa):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        else:
            gene = protein_map.get(p.id) or (
                gene_of_ranking_id(p.id) if p.id.startswith("CIMG_") else None
            )
            r = by_gene.get(gene)
            if r is None:
                rows.append({"id": p.id, "state": "not_in_reference"})
            else:
                rows.append(
                    {
                        "id": p.id,
                        "state": "ok",
                        "percentile": r["percentile"],
                        "antigenicity": r["antigenicity"],
                        "specificity": r["specificity"],
                        "prevalence": r["prevalence"],
                        "max_crossreact": r["max_fungal_crossreact_pid"],
                        "rank": r["rank"],
                    }
                )
    return rows


def cys_rows(proteins, taxa, cys_table, applicable_taxa):
    """``cys_table``: ``{protein_id: row}`` from ``candidates.tsv.gz`` of ``analysis/cys_candidates``."""
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        elif p.id not in cys_table:
            rows.append({"id": p.id, "state": "not_in_reference"})
        else:
            r = cys_table[p.id]
            rows.append({"id": p.id, "state": "ok", "tier": r["tier"], "cys_frac": r["cys_frac"]})
    return rows


def expression_rows(proteins, taxa, protein_map, spherule, applicable_taxa):
    """``spherule``: ``{gene_id: row}`` from ``spherule_surface_table.tsv.gz``."""
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        else:
            r = spherule.get(protein_map.get(p.id, ""))
            if r is None:
                rows.append({"id": p.id, "state": "not_in_reference"})
            else:
                rows.append(
                    {
                        "id": p.id,
                        "state": "ok",
                        "log2fc": r["log2fc_48h"],
                        "padj": r["padj_48h"],
                        "log2fc_8d": r["log2fc_8d"],
                    }
                )
    return rows


def tm_rows(proteins, tmhmm):
    """``tmhmm``: ``{protein_id: row}`` from the TMHMM table (columns ``pred_hel``, ``topology``)."""
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in tmhmm:
            rows.append({"id": p.id, "state": "error"})
        else:
            r = tmhmm[p.id]
            rows.append(
                {"id": p.id, "state": "ok", "n_tm": r["pred_hel"], "topology": r["topology"]}
            )
    return rows
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_lookups.py -q`
Expected: `5 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/lookups.py tests/cellsurface_sorting_hat/modules/test_lookups.py
git commit -m "feat(sorting-hat): antigen, Cys-rich, expression and TM lookups

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 6: Intervals and measures

**Files:**
- Create: `src/cellsurface_sorting_hat/calibration/__init__.py`, `src/cellsurface_sorting_hat/calibration/intervals.py`, `src/cellsurface_sorting_hat/calibration/measure.py`
- Test: `tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py`

**Interfaces:**
- Consumes: `write_atomic` (Plan 1), `load_status_source`, status constants (Plan 1).
- Produces: `wilson(k, n)`, `cluster_bootstrap(y, call, clusters, n_boot=2000, seed=1) -> {sensitivity, specificity}`, `build_measure(...)`, `status_from_measure(measure)`, `make_entry(taxa, measure, source='', cap=None)`, `write_status_source(workdir, module, entries)` (reads `modules/<module>.json` for the identity, refuses two entries that share a tested taxon, validates the result).

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py` with exactly this content:

```python
import json

import numpy as np
import pytest

from cellsurface_sorting_hat.calibration.intervals import cluster_bootstrap, wilson
from cellsurface_sorting_hat.calibration.measure import (
    build_measure,
    make_entry,
    status_from_measure,
    write_status_source,
)
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source


def test_wilson_matches_known_values():
    v, lo, hi = wilson(8, 10)
    assert (
        v == 0.8 and lo == pytest.approx(0.4902, abs=1e-3) and hi == pytest.approx(0.9433, abs=1e-3)
    )
    assert wilson(0, 10)[1] == 0.0 and wilson(10, 10)[2] == pytest.approx(1.0)
    assert wilson(0, 0) is None


def test_bootstrap_point_estimates_and_ordering():
    y = np.array([1] * 30 + [0] * 70)
    call = np.array([True] * 24 + [False] * 6 + [True] * 7 + [False] * 63)
    clusters = [f"c{i}" for i in range(100)]
    r = cluster_bootstrap(y, call, clusters, n_boot=500, seed=3)
    assert r["sensitivity"]["value"] == pytest.approx(0.8) and r["specificity"][
        "value"
    ] == pytest.approx(0.9)
    for k in r.values():
        assert 0 <= k["lo"] <= k["value"] <= k["hi"] <= 1
    assert r == cluster_bootstrap(y, call, clusters, n_boot=500, seed=3)  # same seed, same result


def test_clusters_widen_the_interval_compared_with_single_proteins():
    y = np.array([1] * 40 + [0] * 40)
    call = np.array(([True] * 20 + [False] * 20) * 2)
    single = cluster_bootstrap(y, call, [str(i) for i in range(80)], n_boot=800, seed=1)[
        "sensitivity"
    ]
    # 8 clusters of 5 positives that are all called or all missed: far less independent information
    pos_clusters = [f"p{i // 5}" for i in range(40)]
    call2 = np.array([(i // 5) % 2 == 0 for i in range(40)] + [False] * 40)
    grouped = cluster_bootstrap(
        y, call2, pos_clusters + [f"n{i}" for i in range(40)], n_boot=800, seed=1
    )["sensitivity"]
    assert (grouped["hi"] - grouped["lo"]) > (single["hi"] - single["lo"])


def test_bootstrap_without_negatives_has_no_specificity():
    r = cluster_bootstrap([1, 1, 1], [True, False, True], ["a", "b", "c"], n_boot=100)
    assert r["specificity"] is None and r["sensitivity"] is not None


@pytest.mark.parametrize(
    "args", [([1, 0], [True], ["a", "b"]), ([2, 0], [True, False], ["a", "b"])]
)
def test_bootstrap_refuses_bad_input(args):
    with pytest.raises(ValueError):
        cluster_bootstrap(*args)


def test_status_rule_follows_phase_c_and_needs_negatives_for_an_estimate():
    def m(n_pos, sens, n_neg=None, spec=None):
        out = {"n_pos": n_pos, "sensitivity": sens}
        if spec:
            out.update(n_neg=n_neg, specificity=spec)
        return out

    tight = {"value": 0.6, "lo": 0.55, "hi": 0.65}
    wide = {"value": 0.8, "lo": 0.5, "hi": 0.95}
    sp_ok = {"value": 0.96, "lo": 0.95, "hi": 0.97}
    sp_wide = {"value": 0.9, "lo": 0.6, "hi": 0.99}
    assert status_from_measure(m(232, tight, 4244, sp_ok)) == "estimated"
    assert status_from_measure(m(25, wide, 100, sp_ok)) == "smoke"  # sensitivity half-width 0.225
    assert status_from_measure(m(16, tight, 100, sp_ok)) == "smoke"  # fewer than 20 positives
    assert status_from_measure(m(232, tight, 15, sp_ok)) == "smoke"  # fewer than 20 negatives
    assert status_from_measure(m(232, tight, 4244, sp_wide)) == "smoke"  # specificity too wide
    assert status_from_measure(m(232, tight)) == "smoke"  # never tested on negatives
    # negatives were counted but no specificity was measured: still not an estimate
    assert status_from_measure({"n_pos": 232, "n_neg": 500, "sensitivity": tight}) == "smoke"
    assert status_from_measure({"n_pos": 0}) == "unvalidated"


def test_a_leakage_cap_limits_an_estimate_to_smoke():
    measure = {
        "n_pos": 232,
        "n_neg": 4244,
        "sensitivity": {"value": 0.6, "lo": 0.55, "hi": 0.65},
        "specificity": {"value": 0.96, "lo": 0.95, "hi": 0.97},
    }
    assert make_entry([1], measure)["status"] == "estimated"
    assert make_entry([1], measure, cap="smoke")["status"] == "smoke"


def test_build_measure_counts_and_notes():
    m = build_measure(
        "set",
        "truth.tsv",
        [1, 1, 0, 0],
        [True, False, False, False],
        ["a", "b", "c", "d"],
        notes="n",
        n_boot=50,
    )
    assert (m["n_pos"], m["n_neg"], m["calibration_set"], m["notes"]) == (2, 2, "set", "n")


def test_write_status_source_uses_the_module_identity_and_validates(tmp_path):
    write_module(
        tmp_path, ModuleSpec("pfam_adhesion", "1", {"a": 1}), [], [{"id": "A", "state": "ok"}]
    )
    m = build_measure("set", "t", [1, 0], [True, False], ["a", "b"], n_boot=20)
    path = write_status_source(tmp_path, "pfam_adhesion", [make_entry([40], m)])
    rec = load_status_source(path)
    run = json.loads((tmp_path / "modules" / "pfam_adhesion.json").read_text())
    assert rec.identity.params_hash == run["params_hash"] and rec.entries[0].taxa == (40,)
    with pytest.raises(ValueError, match="appears in two entries"):
        write_status_source(
            tmp_path, "pfam_adhesion", [make_entry([40], m), make_entry([40, 41], m)]
        )
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.calibration`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/calibration/__init__.py` with exactly this content:

```python
"""Calibration: sensitivity and specificity of each module, and the status sources that record them."""
```

Create `src/cellsurface_sorting_hat/calibration/intervals.py` with exactly this content:

```python
"""Confidence intervals for sensitivity and specificity.

``wilson`` is for simple counts. ``cluster_bootstrap`` resamples homology clusters (not single
proteins), so related proteins do not make the interval too narrow. Models are not refitted; the
interval describes sampling of the calibration set only.
"""

import math

import numpy as np


def wilson(k, n, z=1.959964):
    """Wilson score interval for k successes in n trials: ``(value, lo, hi)``; None when n is 0."""
    if n <= 0:
        return None
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return p, max(0.0, centre - half), min(1.0, centre + half)


def _rates(y, call):
    pos, neg = y == 1, y == 0
    sens = (call & pos).sum() / pos.sum() if pos.sum() else np.nan
    spec = (~call & neg).sum() / neg.sum() if neg.sum() else np.nan
    return sens, spec


def cluster_bootstrap(y, call, clusters, n_boot=2000, seed=1):
    """Sensitivity and specificity with a 95% cluster bootstrap interval.

    ``y`` is 1 for a known positive and 0 for a known negative; ``call`` is True when the module
    called the protein. Returns ``{"sensitivity": {...} or None, "specificity": {...} or None}``,
    each ``{"value", "lo", "hi"}`` with ``lo <= value <= hi`` (the percentile interval is widened to
    include the point estimate if the resamples are one-sided).
    """
    y = np.asarray(y, dtype=int)
    call = np.asarray(call, dtype=bool)
    clusters = np.asarray(clusters, dtype=str)
    if not (len(y) == len(call) == len(clusters)):
        raise ValueError("y, call and clusters must have the same length")
    if not set(np.unique(y)) <= {0, 1}:
        raise ValueError("y must hold 0 and 1 only")
    uniq, inverse = np.unique(clusters, return_inverse=True)
    members = [np.flatnonzero(inverse == i) for i in range(len(uniq))]
    point = _rates(y, call)
    rng = np.random.default_rng(seed)
    boot = np.full((n_boot, 2), np.nan)
    for b in range(n_boot):
        draw = rng.integers(0, len(uniq), size=len(uniq))
        idx = np.concatenate([members[i] for i in draw])
        boot[b] = _rates(y[idx], call[idx])
    out = {}
    for j, name in enumerate(("sensitivity", "specificity")):
        value = point[j]
        if np.isnan(value):
            out[name] = None
            continue
        col = boot[:, j][~np.isnan(boot[:, j])]
        lo, hi = np.percentile(col, [2.5, 97.5]) if len(col) else (value, value)
        out[name] = {
            "value": float(value),
            "lo": float(min(lo, value)),
            "hi": float(max(hi, value)),
        }
    return out
```

Create `src/cellsurface_sorting_hat/calibration/measure.py` with exactly this content:

```python
"""Build ``measure`` objects and write status sources (the files that record calibration)."""

import json
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.calibration.intervals import cluster_bootstrap
from cellsurface_sorting_hat.status import ESTIMATED, SMOKE, UNVALIDATED, load_status_source

MIN_POSITIVES = 20
MIN_NEGATIVES = 20
MAX_HALF_WIDTH = 0.10


def build_measure(calibration_set, truth_source, y, call, clusters, notes="", n_boot=2000, seed=1):
    rates = cluster_bootstrap(y, call, clusters, n_boot=n_boot, seed=seed)
    y = list(y)
    measure = {
        "calibration_set": calibration_set,
        "truth_source": truth_source,
        "n_pos": sum(1 for v in y if v == 1),
        "n_neg": sum(1 for v in y if v == 0),
    }
    for key in ("sensitivity", "specificity"):
        if rates[key] is not None:
            measure[key] = rates[key]
    if notes:
        measure["notes"] = notes
    return measure


def status_from_measure(
    measure, min_positives=MIN_POSITIVES, min_negatives=MIN_NEGATIVES, max_half_width=MAX_HALF_WIDTH
):
    """``estimated`` needs at least 20 positives and 20 negatives and a 95% interval half-width of
    at most 0.10 for both sensitivity and specificity (the Phase C rule). A measure with any
    positive but without a specificity is at most ``smoke``: a rule that was never tested on
    negatives is not an estimate. No sensitivity or no positive gives ``unvalidated``."""
    sens, spec = measure.get("sensitivity"), measure.get("specificity")
    if not sens or measure.get("n_pos", 0) < 1:
        return UNVALIDATED
    narrow = (sens["hi"] - sens["lo"]) / 2 <= max_half_width
    if spec:
        narrow = narrow and (spec["hi"] - spec["lo"]) / 2 <= max_half_width
    enough = measure["n_pos"] >= min_positives and measure.get("n_neg", 0) >= min_negatives
    if spec and enough and narrow:
        return ESTIMATED
    return SMOKE


def make_entry(taxa, measure, source="", cap=None):
    """A status entry. ``cap="smoke"`` limits the status, for example when the truth set overlaps
    the data that set the rule."""
    status = status_from_measure(measure)
    if cap == SMOKE and status == ESTIMATED:
        status = SMOKE
    return {"taxa": [int(t) for t in taxa], "status": status, "source": source, "measure": measure}


def write_status_source(workdir, module, entries):
    """Write ``<workdir>/status/<module>.json`` using the identity in ``modules/<module>.json``.

    The file is validated by reading it back with ``load_status_source``. Entries must not share a
    tested taxon (the engine would pick the first of equal depth).
    """
    record = json.loads((Path(workdir) / "modules" / f"{module}.json").read_text())
    seen = {}
    for e in entries:
        for t in e["taxa"]:
            if t in seen:
                raise ValueError(
                    f"taxon {t} appears in two entries ({seen[t]} and {e['measure']['calibration_set']})"
                )
            seen[t] = e["measure"]["calibration_set"]
    data = {k: record[k] for k in ("module", "version", "params_hash", "artefact_hash")}
    data["entries"] = entries
    path = Path(workdir) / "status" / f"{module}.json"
    write_atomic(path, (json.dumps(data, indent=2, sort_keys=True) + "\n").encode())
    load_status_source(path)  # raises if the file is not valid
    return path
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py -q`
Expected: `10 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/calibration tests/cellsurface_sorting_hat/calibration/test_intervals_measure.py
git commit -m "feat(sorting-hat): cluster bootstrap intervals and status measures

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 7: Phase C converter and the strain stability comparison

**Files:**
- Create: `src/cellsurface_sorting_hat/calibration/phasec.py`, `src/cellsurface_sorting_hat/calibration/stability.py`, `data/sorting_hat/phasec_set_species.tsv`
- Test: `tests/cellsurface_sorting_hat/calibration/test_phasec.py`, `tests/cellsurface_sorting_hat/calibration/test_stability.py`

**Interfaces:**
- Consumes: `make_entry` (Task 6); `TaxonError` (Plan 1).
- Produces: `read_names(names_dmp)`, `species_taxid(names, scientific_name)`, `entries_from_phasec(metrics_json, set_taxa, candidate='R0', variant='V-go', truth='direct')`, `SETS`; `compare_runs(out_a, out_b)`.
- `entries_from_phasec` takes the status Phase C gave (`estimate` or `smoke test`) and the sensitivity and false-positive rate of the **all** stratum; specificity is one minus the false-positive rate with the interval ends swapped.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/calibration/test_phasec.py` with exactly this content:

```python
import json

import pytest

from cellsurface_sorting_hat.calibration.measure import (
    write_status_source,
)
from cellsurface_sorting_hat.calibration.phasec import (
    SETS,
    entries_from_phasec,
    read_names,
    species_taxid,
)
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source
from cellsurface_sorting_hat.taxonomy import TaxonError


def _metrics(tmp_path):
    def cell(r, f):
        return {
            "recall": {"value": r[0], "lo": r[1], "hi": r[2]},
            "fpr": {"value": f[0], "lo": f[1], "hi": f[2]},
        }

    def ts(label, pos, neg, c):
        return {
            "label": label,
            "n_direct_positives": pos,
            "truth": {
                "direct": {
                    "n": {"all": {"pos": pos, "neg": neg}},
                    "metrics": {"all": {"V-go": {"R0": c}}},
                }
            },
        }

    data = {
        "test_sets": {
            "S1:all": ts("estimate", 232, 4244, cell((0.603, 0.51, 0.68), (0.037, 0.029, 0.044))),
            "S3-Eurotiomycetes:clade": ts(
                "estimate", 128, 208, cell((0.727, 0.65, 0.80), (0.010, 0.0, 0.03))
            ),
            "S3-Basidiomycota:clade": ts(
                "smoke test", 16, 60, cell((0.938, 0.7, 1.0), (0.083, 0.03, 0.17))
            ),
        }
    }
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps(data))
    return path


def test_phasec_entries_carry_sensitivity_and_specificity(tmp_path):
    set_taxa = {
        "S1:all": [4932, 5476],
        "S3-Eurotiomycetes:clade": [746128, 162425],
        "S3-Basidiomycota:clade": [5207],
    }
    entries = entries_from_phasec(_metrics(tmp_path), set_taxa)
    by = {e["measure"]["calibration_set"]: e for e in entries}
    s1 = by["S1:all"]
    assert s1["status"] == "estimated" and s1["taxa"] == [4932, 5476]
    assert s1["measure"]["sensitivity"] == {"value": 0.603, "lo": 0.51, "hi": 0.68}
    assert s1["measure"]["specificity"]["value"] == pytest.approx(0.963)
    assert s1["measure"]["specificity"]["lo"] == pytest.approx(0.956) and s1["measure"][
        "specificity"
    ]["hi"] == pytest.approx(0.971)
    assert (s1["measure"]["n_pos"], s1["measure"]["n_neg"]) == (232, 4244)
    assert by["S3-Basidiomycota:clade"]["status"] == "smoke"


def test_phasec_entries_make_a_valid_status_source(tmp_path):
    write_module(tmp_path, ModuleSpec("step1_rule@R0", "1"), [], [{"id": "A", "state": "ok"}])
    set_taxa = {
        "S1:all": [4932, 5476],
        "S3-Eurotiomycetes:clade": [746128, 162425],
        "S3-Basidiomycota:clade": [5207],
    }
    path = write_status_source(
        tmp_path, "step1_rule@R0", entries_from_phasec(_metrics(tmp_path), set_taxa)
    )
    assert len(load_status_source(path).entries) == 3


def test_default_sets_do_not_share_a_source():
    sources = [s for group in SETS.values() for s in group]
    assert len(sources) == len(set(sources))


def test_species_taxid_reads_names_dmp_and_refuses_missing_or_ambiguous_names(tmp_path):
    names = tmp_path / "names.dmp"
    names.write_text(
        "4932\t|\tSaccharomyces cerevisiae\t|\t\t|\tscientific name\t|\n"
        "4932\t|\tbaker's yeast\t|\t\t|\tcommon name\t|\n"
        "111\t|\tTwin\t|\t\t|\tscientific name\t|\n"
        "222\t|\tTwin\t|\t\t|\tscientific name\t|\n"
    )
    table = read_names(names)
    assert species_taxid(table, "Saccharomyces cerevisiae") == 4932
    for bad in ("Twin", "Nothing here"):
        with pytest.raises(TaxonError):
            species_taxid(table, bad)
```

Create `tests/cellsurface_sorting_hat/calibration/test_stability.py` with exactly this content:

```python
import csv
import gzip

from cellsurface_sorting_hat.calibration.stability import compare_runs


def write_run(path, proteins, calls):
    path.mkdir()
    with gzip.open(path / "proteins.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["id", "sha256", "taxon", "state", "note", "trailing_stop", "ambiguous_fraction"]
        )
        for pid, sha in proteins:
            w.writerow([pid, sha, 1, "ok", "", 0, 0])
    with gzip.open(path / "calls.long.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein", "call", "variant", "value", "status", "status_basis", "other_basis"])
        for pid, call, variant, value in calls:
            w.writerow([pid, call, variant, value, "unvalidated", "", ""])


def test_only_identical_sequences_are_compared_by_sha256_not_by_id(tmp_path):
    write_run(
        tmp_path / "a",
        [("a1", "s1"), ("a2", "s2"), ("a3", "s3")],
        [("a1", "x", "", "called"), ("a2", "x", "", "called"), ("a3", "x", "", "called")],
    )
    write_run(
        tmp_path / "b",
        [("b9", "s1"), ("b8", "s2"), ("b7", "s4")],
        [("b9", "x", "", "called"), ("b8", "x", "", "not_called"), ("b7", "x", "", "called")],
    )
    r = compare_runs(tmp_path / "a", tmp_path / "b")
    assert (r["n_a"], r["n_b"], r["n_identical_sequences"]) == (3, 3, 2)
    assert dict(r["per_call"][("x", "")]) == {"agree": 1, "differ": 1}
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/calibration/test_phasec.py tests/cellsurface_sorting_hat/calibration/test_stability.py -q`
Expected: collection error, `ModuleNotFoundError: ... calibration.phasec`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/calibration/phasec.py` with exactly this content:

```python
"""Turn the Phase C ``metrics.json`` into status entries for the parameter-free rule R0.

R0 ("SignalP calls a signal peptide") has no fitted parameter, so its sensitivity and false-positive
rate on a test set are measurements of the rule itself. Only the pooled test sets in ``SETS`` are
used, so no tested taxon appears in two entries. Specificity is ``1 - FPR`` over all negative
classes of the set, with the interval ends swapped.
"""

import json
from pathlib import Path

from cellsurface_sorting_hat.calibration.measure import make_entry
from cellsurface_sorting_hat.taxonomy import TaxonError

# Phase C test set -> the source IDs (rows of species.tsv) it pools.
SETS = {
    "S1:all": ("Scer_SGD", "Calb_CGD"),
    "S3-Eurotiomycetes:clade": ("Afum_ASPFU", "Anid_EMENI"),
    "S3-Basidiomycota:clade": ("Cneo_H99_GOA", "Umay_MYCMD"),
}


def read_names(names_dmp):
    """``names.dmp`` -> ``{scientific name: [taxon IDs]}``."""
    out = {}
    for line in Path(names_dmp).read_text(encoding="utf-8-sig").splitlines():
        f = [x.strip() for x in line.split("|")]
        if len(f) >= 4 and f[3] == "scientific name":
            out.setdefault(f[1], []).append(int(f[0]))
    return out


def species_taxid(names, scientific_name):
    """The one taxon ID for a scientific name; refuses a missing or ambiguous name."""
    ids = names.get(scientific_name, [])
    if len(ids) != 1:
        raise TaxonError(f"{scientific_name!r} has {len(ids)} taxon IDs in names.dmp; expected 1")
    return ids[0]


def rate_from_phasec(cell):
    if not cell or cell.get("value") is None:
        return None
    return {"value": cell["value"], "lo": cell["lo"], "hi": cell["hi"]}


def entries_from_phasec(metrics_json, set_taxa, candidate="R0", variant="V-go", truth="direct"):
    """Return a list of entries (see ``make_entry``) for the pooled sets in ``set_taxa``.

    ``set_taxa`` maps a Phase C test set key to its species-level taxon IDs. The status is the one
    Phase C assigned (``label``: estimate or smoke test); no new rule is applied here.
    """
    m = json.loads(Path(metrics_json).read_text())
    entries = []
    for key, taxa in set_taxa.items():
        ts = m["test_sets"][key]
        cell = ts["truth"][truth]["metrics"]["all"][variant][candidate]
        sens = rate_from_phasec(cell["recall"])
        fpr = rate_from_phasec(cell["fpr"])
        measure = {
            "calibration_set": key,
            "truth_source": f"{Path(metrics_json).name}:test_sets/{key}/truth/{truth}",
            "n_pos": int(ts["truth"][truth]["n"]["all"]["pos"]),
            "n_neg": int(ts["truth"][truth]["n"]["all"]["neg"]),
            "notes": f"Phase C {ts['label']}; rule {candidate}, variant {variant}",
        }
        if sens:
            measure["sensitivity"] = sens
        if fpr:
            measure["specificity"] = {
                "value": 1 - fpr["value"],
                "lo": 1 - fpr["hi"],
                "hi": 1 - fpr["lo"],
            }
        entry = make_entry(taxa, measure, source=measure["truth_source"])
        entry["status"] = "estimated" if ts["label"] == "estimate" else "smoke"
        entries.append(entry)
    return entries
```

Create `src/cellsurface_sorting_hat/calibration/stability.py` with exactly this content:

```python
"""Stability of calls between two runs (for example two strains of one species).

Proteins are matched by sequence sha256 (``proteins.tsv.gz``), not by ID. Only proteins whose
sequence is identical in both runs are compared; the others are counted. A difference in calls for
an identical sequence comes from the taxon, the status or a module, not from the sequence.
"""

import csv
import gzip
from collections import Counter, defaultdict


def _read(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _calls(out_dir):
    calls = defaultdict(dict)
    for r in _read(f"{out_dir}/calls.long.tsv.gz"):
        calls[r["protein"]][(r["call"], r["variant"])] = r["value"]
    return calls


def compare_runs(out_a, out_b):
    """Return ``{"n_a", "n_b", "n_identical_sequences", "per_call": {(call, variant): Counter}}``."""
    pa = {r["id"]: r["sha256"] for r in _read(f"{out_a}/proteins.tsv.gz")}
    pb = {r["id"]: r["sha256"] for r in _read(f"{out_b}/proteins.tsv.gz")}
    by_sha_b = defaultdict(list)
    for pid, sha in pb.items():
        by_sha_b[sha].append(pid)
    ca, cb = _calls(out_a), _calls(out_b)
    per_call, n_identical = defaultdict(Counter), 0
    for pid, sha in pa.items():
        if sha not in by_sha_b:
            continue
        n_identical += 1
        other = sorted(by_sha_b[sha])[0]  # identical sequences have identical calls inside one run
        for key, value in ca[pid].items():
            per_call[key]["agree" if cb[other].get(key) == value else "differ"] += 1
    return {
        "n_a": len(pa),
        "n_b": len(pb),
        "n_identical_sequences": n_identical,
        "per_call": dict(per_call),
    }
```

Create `data/sorting_hat/phasec_set_species.tsv` with exactly this content:

```text
set_key	scientific_name
S1:all	Saccharomyces cerevisiae
S1:all	Candida albicans
S3-Eurotiomycetes:clade	Aspergillus fumigatus
S3-Eurotiomycetes:clade	Aspergillus nidulans
S3-Basidiomycota:clade	Cryptococcus neoformans
S3-Basidiomycota:clade	Ustilago maydis
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/calibration/test_phasec.py tests/cellsurface_sorting_hat/calibration/test_stability.py -q`
Expected: `5 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/calibration/phasec.py src/cellsurface_sorting_hat/calibration/stability.py data/sorting_hat/phasec_set_species.tsv tests/cellsurface_sorting_hat/calibration/test_phasec.py tests/cellsurface_sorting_hat/calibration/test_stability.py
git commit -m "feat(sorting-hat): Phase C converter and strain stability comparison

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 8: Proteome provenance

**Files:**
- Create: `src/cellsurface_sorting_hat/proteomes.py`
- Test: `tests/cellsurface_sorting_hat/test_proteomes.py`

**Interfaces:**
- Consumes: `read_fasta`, `OK` (Plan 1); `sha256_file` (Task 1); `write_atomic` (Plan 1).
- Produces: `write_provenance(path, name, source, fasta, retrieved, expected_min, expected_max) -> dict`, `ProteomeError`.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_proteomes.py` with exactly this content:

```python
import json

import pytest

from cellsurface_sorting_hat.proteomes import ProteomeError, write_provenance


def test_provenance_checks_the_count_and_records_the_digest(tmp_path):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n>B\nMKT\n>C\nMKS\n")
    rec = write_provenance(
        tmp_path / "prov.json", "toy", "https://example.org/x", fasta, "2026-10-05", 2, 5
    )
    assert (rec["n_proteins"], rec["n_unique_sequences"], rec["n_invalid"]) == (3, 2, 0)
    assert json.loads((tmp_path / "prov.json").read_text())["sha256"] == rec["sha256"]
    with pytest.raises(ProteomeError, match="expected 10 to 20"):
        write_provenance(tmp_path / "p2.json", "toy", "x", fasta, "2026-10-05", 10, 20)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_proteomes.py -q`
Expected: collection error, `ModuleNotFoundError: ... cellsurface_sorting_hat.proteomes`

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/proteomes.py` with exactly this content:

```python
"""Provenance record for a proteome FASTA that a run uses."""

import json
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.fasta import OK, read_fasta
from cellsurface_sorting_hat.modules.base import sha256_file


class ProteomeError(ValueError):
    """The proteome does not look like the one that was expected."""


def write_provenance(path, name, source, fasta, retrieved, expected_min, expected_max):
    """Check the FASTA and write ``path`` (JSON). ``source`` is a URL or a file path; ``retrieved``
    is the date as text. The protein count must lie in ``[expected_min, expected_max]``."""
    proteins = read_fasta(fasta)
    n = len(proteins)
    if not expected_min <= n <= expected_max:
        raise ProteomeError(f"{name}: {n} proteins, expected {expected_min} to {expected_max}")
    record = {
        "name": name,
        "source": source,
        "retrieved": retrieved,
        "fasta": str(Path(fasta).resolve()),
        "sha256": sha256_file(fasta),
        "n_proteins": n,
        "n_invalid": sum(p.state != OK for p in proteins),
        "n_unique_sequences": len({p.sha256 for p in proteins}),
        "first_ids": [p.id for p in proteins[:3]],
    }
    write_atomic(path, (json.dumps(record, indent=2, sort_keys=True) + "\n").encode())
    return record
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_proteomes.py -q`
Expected: `1 passed`

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/proteomes.py tests/cellsurface_sorting_hat/test_proteomes.py
git commit -m "feat(sorting-hat): proteome provenance record

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 9: Module command and calibration command

**Files:**
- Create: `src/cellsurface_sorting_hat/modules/cli.py`, `src/cellsurface_sorting_hat/calibration/panel.py`, `src/cellsurface_sorting_hat/calibration/cli.py`
- Modify: `src/cellsurface_sorting_hat/categories.yaml`, `pyproject.toml`, `tests/surface_glyco/test_package.py`
- Test: `tests/cellsurface_sorting_hat/modules/test_module_cli.py`, `tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py`

**Interfaces:**
- Consumes: everything from Tasks 1 to 8; `read_taxon_map`, `assign_taxa`, `RunError`, `InputError` (Plan 1 cli).
- Produces the command `cellsurface_sorting_hat_module {signalp,pfam,repeat02,repeat14,allergen,antigen,cys,expression,tm}` (exit 0 or 2; refuses a table in which every protein is `error`) and the command `cellsurface_sorting_hat_calibrate {phasec,truth,allergen-lco,panel}`. `truth` requires `--leakage {none,partial,tuned_on_truth,unknown}`; anything but `none` caps the status at `smoke`. `panel` is report-only.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/modules/test_module_cli.py` with exactly this content:

```python
"""The module command end to end on tiny inputs."""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.modules.cli import main
from cellsurface_sorting_hat.modules.pfam import FAMILY_COLUMNS

FASTA = ">XP_1 a\nMKTAYIAKQRQ\n>XP_2 b\nMNLLPQWERT\n>BAD\nMK*T\n"


def read(workdir, name):
    with gzip.open(workdir / "modules" / f"{name}.tsv.gz", "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "p.faa"
    path.write_text(FASTA)
    return path


def test_signalp_command_writes_the_r0_module(tmp_path, fasta, capsys):
    res = tmp_path / "prediction_results.txt"
    res.write_text(
        "# h\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nXP_1 x\tSP\t0.1\t0.9\t\nXP_2 y\tOTHER\t0.9\t0.1\t\n"
    )
    wd = tmp_path / "wd"
    code = main(
        [
            "signalp",
            "--fasta",
            str(fasta),
            "--workdir",
            str(wd),
            "--results",
            str(res),
            "--signalp-version",
            "6.0h-gpu",
        ]
    )
    assert code == 0 and json.loads(capsys.readouterr().out)["module"] == "step1_rule@R0"
    rows = read(wd, "step1_rule@R0")
    assert (rows["XP_1"]["call"], rows["XP_2"]["call"], rows["BAD"]["state"]) == (
        "called",
        "not_called",
        "na_invalid",
    )
    run = json.loads((wd / "modules" / "step1_rule@R0.json").read_text())
    assert run["tools"] == {"signalp": "6.0h-gpu"} and run["params"]["rule"] == "R0"


def test_pfam_command_writes_both_modules_and_uses_the_recorded_pfam_digest(tmp_path, fasta):
    table = tmp_path / "f.tsv"
    row = dict.fromkeys(FAMILY_COLUMNS, "")
    row.update(
        pfam_acc="PF05730",
        name="CFEM",
        module="pfam_adhesion",
        **{"class": "2b-i"},
        active="yes",
        active_by="owner",
        active_date="2026-10-05",
    )
    row2 = dict(row, pfam_acc="PF16541", name="AltA1", module="pfam_allergen")
    table.write_text(
        "\n".join(
            "\t".join(r[c] for c in FAMILY_COLUMNS)
            for r in [dict(zip(FAMILY_COLUMNS, FAMILY_COLUMNS, strict=True)), row, row2]
        )
        + "\n"
    )
    dom = tmp_path / "d.domtbl"
    dom.write_text(
        "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    wd = tmp_path / "wd"
    assert (
        main(
            [
                "pfam",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--domtbl",
                str(dom),
                "--family-table",
                str(table),
                "--pfam-release",
                "38.2",
                "--pfam-sha256",
                "abc",
            ]
        )
        == 0
    )
    assert read(wd, "pfam_adhesion")["XP_1"]["hit"] == "1"
    assert read(wd, "pfam_allergen")["XP_1"]["hit"] == "0"
    rec = json.loads((wd / "modules" / "pfam_adhesion.json").read_text())
    assert rec["artefact_hash"].startswith("abc:") and rec["params"]["families"] == ["PF05730"]


def test_repeat_command(tmp_path, fasta):
    tab = tmp_path / "r.tsv"
    tab.write_text(
        "protein\trep_period\trep_n_copies\trep_coverage\nXP_1\t5\t3.0\t0.5\nXP_2\t0\t0\t0\n"
    )
    wd = tmp_path / "wd"
    assert main(["repeat02", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tab)]) == 0
    assert read(wd, "repeat02")["XP_1"]["call"] == "called"


def test_allergen_command(tmp_path, fasta):
    blast = tmp_path / "b.tsv"
    blast.write_text("XP_2\tAsp_f_1.0101|11\t82.0\t90\t120\t100\t150\t1e-20\n")
    ref = tmp_path / "a.faa"
    ref.write_text(">Asp_f_1.0101|11\nMKT\n")
    wd = tmp_path / "wd"
    assert (
        main(
            [
                "allergen",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--blast",
                str(blast),
                "--allergen-fasta",
                str(ref),
                "--blast-version",
                "2.16.0+",
            ]
        )
        == 0
    )
    r = read(wd, "allergen_homology")["XP_2"]
    assert (r["identity"], r["coverage"], r["allergen_name"]) == ("82.00", "90.0", "Asp_f_1.0101")


def test_antigen_command_needs_a_taxon_and_marks_other_taxa_not_applicable(tmp_path, fasta, capsys):
    rank = tmp_path / "r.tsv"
    rank.write_text(
        "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
        "CIMG_1-t26_1-p1\t5\t0.05\t2\t1\t1\t0\n"
    )
    pmap = tmp_path / "m.tsv"
    pmap.write_text("protein_id\tgene_id\tproduct\tlength\nXP_1\tCIMG_1\tp\t11\n")
    wd = tmp_path / "wd"
    base = [
        "antigen",
        "--fasta",
        str(fasta),
        "--workdir",
        str(wd),
        "--ranking",
        str(rank),
        "--protein-map",
        str(pmap),
    ]
    assert main(base) == 2 and "give --taxon or --taxon-map" in capsys.readouterr().err
    assert main(base + ["--taxon", "746128"]) == 0
    assert read(wd, "antigen_lookup")["XP_1"]["state"] == "not_applicable"
    assert main(base + ["--taxon", "246410"]) == 0
    assert read(wd, "antigen_lookup")["XP_1"]["percentile"] == "0.05"


def test_tm_command(tmp_path, fasta):
    tab = tmp_path / "t.tsv"
    tab.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nXP_1\t11\t0\t0\t0\to\nXP_2\t10\t20\t5\t1\ti5-27o\n"
    )
    wd = tmp_path / "wd"
    assert main(["tm", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tab)]) == 0
    assert read(wd, "tm")["XP_2"]["n_tm"] == "1"


def test_input_errors_exit_2(tmp_path, fasta, capsys):
    assert (
        main(
            [
                "signalp",
                "--fasta",
                str(fasta),
                "--workdir",
                str(tmp_path),
                "--results",
                str(tmp_path / "none"),
                "--signalp-version",
                "6",
            ]
        )
        == 2
    )
    assert "cellsurface_sorting_hat_module: error" in capsys.readouterr().err


def test_a_result_file_whose_ids_do_not_match_the_fasta_is_refused(tmp_path, fasta, capsys):
    res = tmp_path / "prediction_results.txt"
    res.write_text(
        "# h\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nsp|Q1|X\tSP\t0.1\t0.9\t\n"
    )
    code = main(
        [
            "signalp",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--results",
            str(res),
            "--signalp-version",
            "6.0h-gpu",
        ]
    )
    assert code == 2 and "no protein of the FASTA has a result" in capsys.readouterr().err
    assert not (tmp_path / "wd" / "modules").exists()  # nothing was written


def test_the_signalp_module_identity_does_not_depend_on_the_proteome(tmp_path, fasta):
    """A measurement of rule R0 belongs to the tool version, so it must stay valid for another proteome."""
    other = tmp_path / "o.faa"
    other.write_text(">Z1\nMKTAYI\n")
    ids = {}
    for name, f, pid in (("a", fasta, "XP_1"), ("b", other, "Z1")):
        res = tmp_path / f"{name}.txt"
        res.write_text(
            f"# h\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n{pid}\tSP\t0.1\t0.9\t\n"
        )
        wd = tmp_path / name
        assert (
            main(
                [
                    "signalp",
                    "--fasta",
                    str(f),
                    "--workdir",
                    str(wd),
                    "--results",
                    str(res),
                    "--signalp-version",
                    "6.0h-gpu",
                ]
            )
            == 0
        )
        ids[name] = json.loads((wd / "modules" / "step1_rule@R0.json").read_text())
    assert (
        ids["a"]["artefact_hash"] == ids["b"]["artefact_hash"]
        and ids["a"]["params_hash"] == ids["b"]["params_hash"]
    )
    changed = tmp_path / "c"
    assert (
        main(
            [
                "signalp",
                "--fasta",
                str(fasta),
                "--workdir",
                str(changed),
                "--results",
                str(tmp_path / "a.txt"),
                "--signalp-version",
                "6.0i-gpu",
            ]
        )
        == 0
    )
    assert (
        json.loads((changed / "modules" / "step1_rule@R0.json").read_text())["artefact_hash"]
        != ids["a"]["artefact_hash"]
    )


def test_pfam_command_applies_the_no_tm_condition_from_the_tm_table(tmp_path, fasta):
    table = tmp_path / "f.tsv"
    row = dict.fromkeys(FAMILY_COLUMNS, "")
    row.update(
        pfam_acc="PF05730",
        name="CFEM",
        module="pfam_adhesion",
        **{"class": "2b-i"},
        second_condition="no_tm",
        active="yes",
        active_by="owner",
        active_date="2026-10-05",
    )
    header = dict(zip(FAMILY_COLUMNS, FAMILY_COLUMNS, strict=True))
    table.write_text(
        "\n".join("\t".join(r[c] for c in FAMILY_COLUMNS) for r in [header, row]) + "\n"
    )
    dom = tmp_path / "d.domtbl"
    dom.write_text(
        "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
        "XP_2 - 10 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    tm = tmp_path / "t.tsv"
    tm.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nXP_1\t11\t0\t0\t0\to\nXP_2\t10\t150\t20\t7\to5-27i\n"
    )
    wd = tmp_path / "wd"
    assert main(["tm", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tm)]) == 0
    assert (
        main(
            [
                "pfam",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--domtbl",
                str(dom),
                "--family-table",
                str(table),
                "--pfam-release",
                "38.2",
                "--pfam-sha256",
                "abc",
                "--tm-module",
                "tm",
            ]
        )
        == 0
    )
    rows = read(wd, "pfam_adhesion")
    assert (rows["XP_1"]["hit"], rows["XP_2"]["hit"]) == (
        "1",
        "0",
    )  # XP_2 has 7 helices: a receptor-like protein
```

Create `tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py` with exactly this content:

```python
import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.calibration.cli import main
from cellsurface_sorting_hat.calibration.panel import panel_check
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source


def write_calls(path, rows):
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein", "call", "variant", "value", "status", "status_basis", "other_basis"])
        for p, c, v, val in rows:
            w.writerow([p, c, v, val, "unvalidated", "", ""])


def write_truth(path, rows):
    path.write_text("id\tlabel\tcluster\n" + "".join(f"{i}\t{y}\t{c}\n" for i, y, c in rows))


def test_truth_command_writes_an_entry_with_sensitivity_and_specificity(tmp_path, capsys):
    write_module(tmp_path, ModuleSpec("allergen_homology", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(40)]
    calls = [
        (i, "allergen_candidate", "", "called" if k < 12 or 20 <= k < 24 else "not_called")
        for k, i in enumerate(ids)
    ]
    calls.append(("U1", "allergen_candidate", "", "not_assessable"))
    write_calls(tmp_path / "c.tsv.gz", calls)
    write_truth(
        tmp_path / "t.tsv",
        [(i, 1 if k < 20 else 0, f"c{k}") for k, i in enumerate(ids)]
        + [("U1", 1, "cu"), ("GONE", 0, "cg")],
    )
    code = main(
        [
            "truth",
            "--workdir",
            str(tmp_path),
            "--module",
            "allergen_homology",
            "--calls-long",
            str(tmp_path / "c.tsv.gz"),
            "--call",
            "allergen_candidate",
            "--truth",
            str(tmp_path / "t.tsv"),
            "--calibration-set",
            "toy",
            "--leakage",
            "none",
            "--taxa",
            "746128",
            "--n-boot",
            "200",
        ]
    )
    assert code == 0
    entry = load_status_source(tmp_path / "status" / "allergen_homology.json").entries[0]
    m = entry.measure
    assert (m["n_pos"], m["n_neg"]) == (20, 20)
    assert m["sensitivity"]["value"] == pytest.approx(12 / 20) and m["specificity"][
        "value"
    ] == pytest.approx(16 / 20)
    assert "truth rows without a call: 1" in m["notes"] and "not assessable: 1" in m["notes"]
    assert entry.status == "smoke"  # 20 positives but a wide interval


def test_a_second_set_is_added_and_the_same_set_is_replaced(tmp_path):
    write_module(tmp_path, ModuleSpec("repeat02", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(10)]
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            (i, "adhesion_repeat", "", "called" if k < 5 else "not_called")
            for k, i in enumerate(ids)
        ],
    )
    write_truth(tmp_path / "t.tsv", [(i, 1 if k < 5 else 0, f"c{k}") for k, i in enumerate(ids)])

    def run(name, taxon):
        return main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "repeat02",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "adhesion_repeat",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                name,
                "--leakage",
                "none",
                "--taxa",
                str(taxon),
                "--n-boot",
                "50",
            ]
        )

    assert run("setA", 4932) == 0 and run("setB", 5476) == 0 and run("setA", 4932) == 0
    sets = [
        e.measure["calibration_set"]
        for e in load_status_source(tmp_path / "status" / "repeat02.json").entries
    ]
    assert sorted(sets) == ["setA", "setB"]


def test_truth_with_no_matching_protein_is_an_error(tmp_path, capsys):
    write_module(tmp_path, ModuleSpec("repeat02", "1"), [], [{"id": "A", "state": "ok"}])
    write_calls(tmp_path / "c.tsv.gz", [("X", "adhesion_repeat", "", "called")])
    write_truth(tmp_path / "t.tsv", [("Y", 1, "c")])
    assert (
        main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "repeat02",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "adhesion_repeat",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                "s",
                "--leakage",
                "none",
                "--taxa",
                "1",
            ]
        )
        == 2
    )
    assert "no truth protein has a call" in capsys.readouterr().err


def test_phasec_command_resolves_species_names_to_taxa(tmp_path):
    names = tmp_path / "names.dmp"
    names.write_text(
        "".join(
            f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n"
            for t, n in [
                (4932, "Saccharomyces cerevisiae"),
                (5476, "Candida albicans"),
                (746128, "Aspergillus fumigatus"),
                (162425, "Aspergillus nidulans"),
                (5207, "Cryptococcus neoformans"),
                (5270, "Ustilago maydis"),
            ]
        )
    )
    sp = tmp_path / "sets.tsv"
    sp.write_text(
        "set_key\tscientific_name\nS1:all\tSaccharomyces cerevisiae\nS1:all\tCandida albicans\n"
        "S3-Eurotiomycetes:clade\tAspergillus fumigatus\nS3-Eurotiomycetes:clade\tAspergillus nidulans\n"
        "S3-Basidiomycota:clade\tCryptococcus neoformans\nS3-Basidiomycota:clade\tUstilago maydis\n"
    )

    def cell(r, f):
        return {
            "recall": {"value": r, "lo": r - 0.05, "hi": r + 0.05},
            "fpr": {"value": f, "lo": f / 2, "hi": f * 2},
        }

    def ts(label, pos, c):
        return {
            "label": label,
            "n_direct_positives": pos,
            "truth": {
                "direct": {
                    "n": {"all": {"pos": pos, "neg": 100}},
                    "metrics": {"all": {"V-go": {"R0": c}}},
                }
            },
        }

    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        json.dumps(
            {
                "test_sets": {
                    "S1:all": ts("estimate", 232, cell(0.6, 0.04)),
                    "S3-Eurotiomycetes:clade": ts("estimate", 128, cell(0.7, 0.01)),
                    "S3-Basidiomycota:clade": ts("smoke test", 16, cell(0.9, 0.08)),
                }
            }
        )
    )
    write_module(
        tmp_path / "wd", ModuleSpec("step1_rule@R0", "1"), [], [{"id": "A", "state": "ok"}]
    )
    assert (
        main(
            [
                "phasec",
                "--workdir",
                str(tmp_path / "wd"),
                "--metrics",
                str(metrics),
                "--set-species",
                str(sp),
                "--names-dmp",
                str(names),
            ]
        )
        == 0
    )
    entries = load_status_source(tmp_path / "wd" / "status" / "step1_rule@R0.json").entries
    assert {frozenset(e.taxa) for e in entries} == {
        frozenset({4932, 5476}),
        frozenset({746128, 162425}),
        frozenset({5207, 5270}),
    }
    assert {e.measure["calibration_set"]: e.status for e in entries}[
        "S3-Basidiomycota:clade"
    ] == "smoke"


def test_panel_check_counts_agreement_and_leaves_known_misses_out(tmp_path):
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            ("A", "adhesion_repeat", "", "called"),
            ("B", "adhesion_repeat", "", "not_called"),
            ("C", "adhesion_repeat", "", "called"),
        ],
    )
    panel = tmp_path / "p.tsv"
    panel.write_text(
        "protein\tcall\tvariant\texpected\tsource\nA\tadhesion_repeat\t\tcalled\tx\nB\tadhesion_repeat\t\tcalled\tx\n"
        "C\tadhesion_repeat\t\tknown_miss\tx\nZ\tadhesion_repeat\t\tcalled\tx\n"
    )
    rows, summary = panel_check(tmp_path / "c.tsv.gz", panel)
    assert summary == {"agree": 1, "disagree": 2, "known_miss": 1}
    assert [r["observed"] for r in rows] == ["called", "not_called", "called", "missing"]


def test_allergen_lco_command_prints_recall_per_cutoff(tmp_path, capsys):
    blast = tmp_path / "b.tsv"
    blast.write_text(
        "a\tb\t99\t100\t100\t100\t300\t1e-50\nb\ta\t99\t100\t100\t100\t300\t1e-50\n"
        "a\tc\t38\t80\t100\t100\t60\t1e-5\nc\ta\t38\t80\t100\t100\t60\t1e-5\n"
        "b\tc\t40\t60\t100\t100\t50\t1e-3\nc\tb\t40\t60\t100\t100\t50\t1e-3\n"
    )
    assert main(["allergen-lco", "--blast", str(blast), "--cutoffs", "35:50,70:80"]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "sequences\t3\tclusters\t2"
    assert out[1].endswith("3/3") and out[2].endswith("0/3")


def test_leakage_other_than_none_caps_an_estimate_at_smoke(tmp_path):
    write_module(tmp_path, ModuleSpec("antigen_lookup", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(60)]
    # 30 positives, all called; 30 negatives, none called: a tight interval that would be an estimate
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            (i, "antigen_candidate", "", "called" if k < 30 else "not_called")
            for k, i in enumerate(ids)
        ],
    )
    write_truth(tmp_path / "t.tsv", [(i, 1 if k < 30 else 0, f"c{k}") for k, i in enumerate(ids)])

    def run(leakage, name):
        return main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "antigen_lookup",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "antigen_candidate",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                name,
                "--taxa",
                "246410" if name == "a" else "5476",
                "--leakage",
                leakage,
                "--n-boot",
                "100",
            ]
        )

    assert run("none", "a") == 0 and run("tuned_on_truth", "b") == 0
    status = {
        e.measure["calibration_set"]: (e.status, e.measure["notes"])
        for e in load_status_source(tmp_path / "status" / "antigen_lookup.json").entries
    }
    assert status["a"][0] == "estimated" and status["b"][0] == "smoke"
    assert "leakage: tuned_on_truth" in status["b"][1]


def test_pfam_specificity_command_lists_non_member_hits_with_their_helices(tmp_path, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n>B\nMKT\n>C\nMKT\n>D\nMKT\n")
    dom = tmp_path / "d.domtbl"
    row = "{} - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    dom.write_text(row.format("A") + row.format("B") + row.format("D"))
    members = tmp_path / "m.tsv"
    members.write_text("pfam_acc\tprotein_id\nPF05730\tA\nPF05730\tC\nPF01185\tD\n")
    tm = tmp_path / "t.tsv"
    tm.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nB\t11\t150\t20\t7\to5-27i\n"
    )
    assert (
        main(
            [
                "pfam-specificity",
                "--family",
                "PF05730",
                "--domtbl",
                str(dom),
                "--members",
                str(members),
                "--universe-fasta",
                str(fasta),
                "--tm-table",
                str(tm),
            ]
        )
        == 0
    )
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "family\tPF05730\thits\t3\tmembers\t2"
    assert "tp\t1" in out and "fp\t2" in out and "fn\t1" in out
    assert (
        "nonmember_hit\tB\tn_tm=7" in out
        and "nonmember_hit\tD\tn_tm=NA" in out
        and "missed_member\tC" in out
    )


def test_pfam_specificity_refuses_members_outside_the_proteome(tmp_path, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n")
    dom = tmp_path / "d.domtbl"
    dom.write_text("")
    members = tmp_path / "m.tsv"
    members.write_text("pfam_acc\tprotein_id\nPF05730\tZ\n")
    assert (
        main(
            [
                "pfam-specificity",
                "--family",
                "PF05730",
                "--domtbl",
                str(dom),
                "--members",
                str(members),
                "--universe-fasta",
                str(fasta),
            ]
        )
        == 2
    )
    assert "not in the proteome" in capsys.readouterr().err
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_module_cli.py tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py -q`
Expected: collection error, `ModuleNotFoundError: ... modules.cli`.

- [ ] **Step 3: Write the implementation**

Create `src/cellsurface_sorting_hat/modules/cli.py` with exactly this content:

```python
"""Command ``cellsurface_sorting_hat_module``: turn one tool's output into a module table."""

import argparse
import json
import sys
from pathlib import Path

from cellsurface_sorting_hat.cli import InputError, RunError, assign_taxa, read_taxon_map
from cellsurface_sorting_hat.fasta import FastaError, read_fasta
from cellsurface_sorting_hat.modules import allergen, lookups, pfam, repeats, signalp
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module

APPLICABLE_RS = {
    246410
}  # C. immitis RS (taxon_id in analysis/step1_compare/phasec/tc_taxon_clades.tsv)


def _common(p):
    p.add_argument("--fasta", required=True)
    p.add_argument("--workdir", required=True)


def _taxa_args(p):
    p.add_argument("--taxon", type=int)
    p.add_argument("--taxon-map")
    p.add_argument("--applicable-taxa", type=int, nargs="+", default=sorted(APPLICABLE_RS))


def build_parser():
    ap = argparse.ArgumentParser(prog="cellsurface_sorting_hat_module")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("signalp", help="SignalP 6 prediction_results.txt -> step1_rule@R0")
    _common(p)
    p.add_argument("--results", required=True)
    p.add_argument("--signalp-version", required=True, help="for example 6.0h-gpu (recorded)")
    p.add_argument("--signalp-mode", default="fast")

    p = sub.add_parser("pfam", help="hmmsearch --domtblout -> pfam_adhesion and pfam_allergen")
    _common(p)
    p.add_argument("--domtbl", required=True)
    p.add_argument("--family-table", required=True)
    p.add_argument("--pfam-release", required=True)
    p.add_argument(
        "--pfam-sha256", required=True, help="sha256 of Pfam-A.hmm, from pfam_provenance.json"
    )
    p.add_argument(
        "--sp-module", help="step1_rule@R0 table in --workdir, for second_condition=signal_peptide"
    )
    p.add_argument("--tm-module", help="tm table in --workdir, for second_condition=no_tm")

    for name in ("repeat02", "repeat14"):
        p = sub.add_parser(name, help="repeat detector table -> " + name)
        _common(p)
        p.add_argument("--table", required=True)
        p.add_argument("--min-coverage", type=float, default=repeats.DEFAULT_MIN_COVERAGE)
        p.add_argument("--min-copies", type=float, default=repeats.DEFAULT_MIN_COPIES)

    p = sub.add_parser("allergen", help="BLASTP results against the IUIS fungal allergens")
    _common(p)
    p.add_argument("--blast", required=True, help="outfmt 6 with: " + allergen.BLAST_FIELDS)
    p.add_argument("--allergen-fasta", required=True)
    p.add_argument("--blast-version", required=True)
    p.add_argument("--evalue", default="1")

    p = sub.add_parser("antigen", help="antigen ranking lookup (Coccidioides)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--ranking", required=True)
    p.add_argument("--protein-map", required=True)

    p = sub.add_parser("cys", help="Cys-rich tiers (analysis/cys_candidates candidates.tsv.gz)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--candidates", required=True)

    p = sub.add_parser("expression", help="spherule table lookup (Coccidioides)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--table", required=True)
    p.add_argument("--protein-map", required=True)

    p = sub.add_parser("tm", help="TMHMM table -> tm")
    _common(p)
    p.add_argument("--table", required=True)
    p.add_argument("--tmhmm-version", default="2.0c")
    return ap


def _taxa(args, proteins):
    tmap = read_taxon_map(args.taxon_map) if args.taxon_map else {}
    if args.taxon is None and not tmap:
        raise RunError("give --taxon or --taxon-map")
    return assign_taxa(proteins, args.taxon, tmap)


def _write(workdir, spec, columns, rows):
    """Write the module, but refuse a table in which every protein is ``error`` (an ID mismatch)."""
    real = [r for r in rows if r["state"] != "na_invalid"]
    if real and all(r["state"] == "error" for r in real):
        raise RunError(
            f"{spec.name}: no protein of the FASTA has a result; check that the IDs agree"
        )
    return write_module(workdir, spec, columns, rows)


def run(args):
    proteins = read_fasta(args.fasta)
    w = args.workdir
    if args.cmd == "signalp":
        spec = ModuleSpec(
            "step1_rule@R0",
            "1",
            {"rule": "R0", "mode": args.signalp_mode},
            (),
            {"signalp": args.signalp_version},
            artefact_digest=_tool_digest("signalp", args.signalp_version),
        )
        rows = signalp.signalp_rows(proteins, signalp.parse_signalp(args.results))
        return _write(w, spec, signalp.COLUMNS, rows)
    if args.cmd == "pfam":
        families = pfam.load_family_table(args.family_table)
        hits = pfam.parse_domtblout(args.domtbl)
        sp_calls = _module_column(w, args.sp_module, "call") if args.sp_module else None
        tm_counts = None
        if args.tm_module:
            raw = _module_column(w, args.tm_module, "n_tm")
            tm_counts = {k: int(v) for k, v in raw.items() if v.isdigit()}
        out = None
        for module in pfam.MODULES:
            params = {
                "pfam_release": args.pfam_release,
                "cut": "ga",
                "families": sorted(f.pfam_acc for f in families if f.module == module and f.active),
            }
            spec = ModuleSpec(
                module,
                "1",
                params,
                (args.family_table,),
                {"hmmer": "3.4"},
                artefact_digest=args.pfam_sha256 + ":" + _file_digest(args.family_table),
            )
            out = write_module(
                w,
                spec,
                pfam.COLUMNS,
                pfam.pfam_rows(proteins, hits, families, module, sp_calls, tm_counts),
            )
        return out
    if args.cmd in ("repeat02", "repeat14"):
        spec = ModuleSpec(
            args.cmd, "1", {"min_coverage": args.min_coverage, "min_copies": args.min_copies}
        )
        rows = repeats.repeat_rows(
            proteins, repeats.parse_repeat_table(args.table), args.min_coverage, args.min_copies
        )
        return _write(w, spec, repeats.COLUMNS, rows)
    if args.cmd == "allergen":
        spec = ModuleSpec(
            "allergen_homology",
            "1",
            {"evalue": args.evalue, "program": "blastp"},
            (args.allergen_fasta,),
            {"blast": args.blast_version},
        )
        rows = allergen.allergen_rows(proteins, allergen.parse_blast(args.blast))
        return _write(w, spec, allergen.COLUMNS, rows)
    if args.cmd == "tm":
        table, _ = lookups.read_table(args.table, "protein_id")
        spec = ModuleSpec(
            "tm",
            "1",
            {},
            (),
            {"tmhmm": args.tmhmm_version},
            artefact_digest=_tool_digest("tmhmm", args.tmhmm_version),
        )
        return _write(w, spec, lookups.TM_COLUMNS, lookups.tm_rows(proteins, table))
    taxa = _taxa(args, proteins)
    applicable = set(args.applicable_taxa)
    if args.cmd == "antigen":
        by_gene, _ = lookups.ranking_by_gene(args.ranking)
        pmap = lookups.load_protein_map(args.protein_map)
        spec = ModuleSpec(
            "antigen_lookup",
            "1",
            {"applicable_taxa": sorted(applicable)},
            (args.ranking, args.protein_map),
        )
        return write_module(
            w,
            spec,
            lookups.ANTIGEN_COLUMNS,
            lookups.antigen_rows(proteins, taxa, pmap, by_gene, applicable),
        )
    if args.cmd == "cys":
        table, _ = lookups.read_table(args.candidates, "protein_id")
        spec = ModuleSpec(
            "cys_rich", "1", {"applicable_taxa": sorted(applicable)}, (args.candidates,)
        )
        return write_module(
            w, spec, lookups.CYS_COLUMNS, lookups.cys_rows(proteins, taxa, table, applicable)
        )
    if args.cmd == "expression":
        table, _ = lookups.read_table(args.table, "gene_id")
        pmap = lookups.load_protein_map(args.protein_map)
        spec = ModuleSpec(
            "expression",
            "1",
            {"applicable_taxa": sorted(applicable)},
            (args.table, args.protein_map),
        )
        return write_module(
            w,
            spec,
            lookups.EXPRESSION_COLUMNS,
            lookups.expression_rows(proteins, taxa, pmap, table, applicable),
        )
    raise AssertionError(args.cmd)


def _tool_digest(tool, version):
    """Identity of a tool: the measurement of a rule belongs to the tool version, not to one input."""
    import hashlib

    return hashlib.sha256(f"{tool}:{version}".encode()).hexdigest()


def _module_column(workdir, module, column):
    import csv
    import gzip

    with gzip.open(Path(workdir) / "modules" / f"{module}.tsv.gz", "rt", newline="") as fh:
        return {r["id"]: r.get(column, "") for r in csv.DictReader(fh, delimiter="\t")}


def _file_digest(path):
    from cellsurface_sorting_hat.modules.base import sha256_file

    return sha256_file(path)


def main(argv=None):
    try:
        record = run(build_parser().parse_args(argv))
    except (RunError, InputError, FastaError, ValueError, OSError) as err:
        print(f"cellsurface_sorting_hat_module: error: {err}", file=sys.stderr)
        return 2
    print(json.dumps({k: record[k] for k in ("module", "run_state", "n_rows")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

Create `src/cellsurface_sorting_hat/calibration/panel.py` with exactly this content:

```python
"""Report-only check of calls against a panel of proteins with expected values."""

import csv
import gzip
from collections import Counter

EXPECTED = ("called", "not_called", "known_miss", "excluded")


def read_calls(path):
    with gzip.open(path, "rt", newline="") as fh:
        return {
            (r["protein"], r["call"], r["variant"]): r["value"]
            for r in csv.DictReader(fh, delimiter="\t")
        }


def panel_check(calls_long, panel_tsv):
    """Compare ``calls_long`` with ``panel_tsv`` (columns protein, call, variant, expected, source).

    ``known_miss`` and ``excluded`` rows are listed and not counted as agreement or disagreement.
    Returns ``(rows, summary)``. Nothing here is an acceptance gate.
    """
    calls = read_calls(calls_long)
    rows, summary = [], Counter()
    with open(panel_tsv, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["expected"] not in EXPECTED:
                raise ValueError(f"{r['protein']}: expected must be one of {EXPECTED}")
            observed = calls.get((r["protein"], r["call"], r["variant"]), "missing")
            if r["expected"] in ("known_miss", "excluded"):
                verdict = r["expected"]
            else:
                verdict = "agree" if observed == r["expected"] else "disagree"
            summary[verdict] += 1
            rows.append({**r, "observed": observed, "verdict": verdict})
    return rows, dict(summary)
```

Create `src/cellsurface_sorting_hat/calibration/cli.py` with exactly this content:

```python
"""Command ``cellsurface_sorting_hat_calibrate``: write status sources and check panels."""

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

from cellsurface_sorting_hat.calibration import phasec
from cellsurface_sorting_hat.calibration.measure import (
    build_measure,
    make_entry,
    write_status_source,
)
from cellsurface_sorting_hat.calibration.panel import panel_check
from cellsurface_sorting_hat.fasta import read_fasta
from cellsurface_sorting_hat.modules import allergen, pfam
from cellsurface_sorting_hat.taxonomy import TaxonError


def build_parser():
    ap = argparse.ArgumentParser(prog="cellsurface_sorting_hat_calibrate")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser(
        "phasec", help="status source for step1_rule@R0 from the Phase C metrics.json"
    )
    p.add_argument("--workdir", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--set-species", required=True, help="TSV: set_key, scientific_name")
    p.add_argument("--names-dmp", required=True)

    p = sub.add_parser(
        "truth", help="sensitivity and specificity of one call against a truth table"
    )
    p.add_argument("--workdir", required=True)
    p.add_argument("--module", required=True, help="module that receives the status entry")
    p.add_argument(
        "--calls-long", required=True, help="calls.long.tsv.gz of a run on the truth proteins"
    )
    p.add_argument("--call", required=True)
    p.add_argument("--variant", default="")
    p.add_argument("--truth", required=True, help="TSV: id, label (1 or 0), cluster")
    p.add_argument("--calibration-set", required=True)
    p.add_argument("--taxa", type=int, nargs="+", required=True, help="tested taxa (species level)")
    p.add_argument("--notes", default="")
    p.add_argument(
        "--leakage",
        required=True,
        choices=["none", "partial", "tuned_on_truth", "unknown"],
        help="did the truth proteins help to set the rule or its cutoffs? anything but 'none' caps the status at smoke",
    )
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=1)

    p = sub.add_parser(
        "allergen-lco", help="leave-cluster-out recall of the allergen set (sensitivity only)"
    )
    p.add_argument(
        "--blast",
        required=True,
        help="allergens against allergens, outfmt 6 (see allergen.BLAST_FIELDS)",
    )
    p.add_argument(
        "--cutoffs", default="35:50,50:70,70:80", help="identity:coverage pairs, comma separated"
    )
    p.add_argument("--cluster-identity", type=float, default=40.0)
    p.add_argument("--cluster-coverage", type=float, default=70.0)

    p = sub.add_parser(
        "pfam-specificity", help="specificity test of one Pfam family on one proteome"
    )
    p.add_argument("--family", required=True, help="for example PF05730")
    p.add_argument("--domtbl", required=True, help="hmmsearch --domtblout of the proteome")
    p.add_argument(
        "--members",
        required=True,
        help="TSV: pfam_acc, protein_id (members known from curation, not from Pfam)",
    )
    p.add_argument("--universe-fasta", required=True, help="the proteome that was searched")
    p.add_argument(
        "--tm-table",
        help="TMHMM table (protein_id, pred_hel) to show helices beside each non-member hit",
    )

    p = sub.add_parser("panel", help="report-only check of calls against a panel")
    p.add_argument("--calls-long", required=True)
    p.add_argument("--panel", required=True)
    return ap


def _set_taxa(path, names_dmp):
    names = phasec.read_names(names_dmp)
    by_set = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by_set.setdefault(r["set_key"], []).append(
                phasec.species_taxid(names, r["scientific_name"])
            )
    return by_set


def _merge_entries(path, new_entries):
    """Entries already in the status source stay, except those of the same calibration set."""
    if not Path(path).exists():
        return list(new_entries)
    old = json.loads(Path(path).read_text())["entries"]
    names = {e["measure"]["calibration_set"] for e in new_entries}
    return [e for e in old if e["measure"]["calibration_set"] not in names] + list(new_entries)


def run(args):
    if args.cmd == "phasec":
        entries = phasec.entries_from_phasec(
            args.metrics, _set_taxa(args.set_species, args.names_dmp)
        )
        path = Path(args.workdir) / "status" / "step1_rule@R0.json"
        return write_status_source(args.workdir, "step1_rule@R0", _merge_entries(path, entries))
    if args.cmd == "truth":
        with gzip.open(args.calls_long, "rt", newline="") as fh:
            calls = {
                r["protein"]: r["value"]
                for r in csv.DictReader(fh, delimiter="\t")
                if r["call"] == args.call and r["variant"] == args.variant
            }
        y, called, clusters, unmatched, unknown = [], [], [], 0, 0
        with open(args.truth, encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                value = calls.get(r["id"])
                if value is None:
                    unmatched += 1
                elif value == "not_assessable":
                    unknown += 1
                else:
                    y.append(int(r["label"]))
                    called.append(value == "called")
                    clusters.append(r["cluster"])
        if not y:
            raise ValueError("no truth protein has a call in --calls-long")
        notes = f"{args.notes} truth rows without a call: {unmatched}; not assessable: {unknown}".strip()
        measure = build_measure(
            args.calibration_set,
            str(args.truth),
            y,
            called,
            clusters,
            notes,
            args.n_boot,
            args.seed,
        )
        path = Path(args.workdir) / "status" / f"{args.module}.json"
        measure["notes"] = f"{measure['notes']} leakage: {args.leakage}".strip()
        cap = None if args.leakage == "none" else "smoke"
        entry = make_entry(args.taxa, measure, source=str(args.truth), cap=cap)
        return write_status_source(args.workdir, args.module, _merge_entries(path, [entry]))
    if args.cmd == "pfam-specificity":
        hits = {h["target"] for h in pfam.parse_domtblout(args.domtbl) if h["acc"] == args.family}
        with open(args.members, encoding="utf-8-sig", newline="") as fh:
            members = {
                r["protein_id"]
                for r in csv.DictReader(fh, delimiter="\t")
                if r["pfam_acc"] == args.family
            }
        universe = {p.id for p in read_fasta(args.universe_fasta)}
        outside = members - universe
        if outside:
            raise ValueError(
                f"{len(outside)} member(s) are not in the proteome, for example {sorted(outside)[0]}"
            )
        rep = pfam.specificity_report(hits & universe, members, universe)
        print(f"family\t{args.family}\thits\t{len(hits & universe)}\tmembers\t{len(members)}")
        for key in ("tp", "fp", "fn", "tn", "sensitivity", "specificity"):
            print(f"{key}\t{rep[key]}")
        tm = {}
        if args.tm_table:
            with open(args.tm_table, encoding="utf-8-sig", newline="") as fh:
                tm = {r["protein_id"]: r["pred_hel"] for r in csv.DictReader(fh, delimiter="\t")}
        for pid in rep["nonmember_hits"]:
            print(f"nonmember_hit\t{pid}\tn_tm={tm.get(pid, 'NA')}")
        for pid in rep["missed_members"]:
            print(f"missed_member\t{pid}")
        return None
    if args.cmd == "allergen-lco":
        cutoffs = [tuple(float(x) for x in c.split(":")) for c in args.cutoffs.split(",")]
        rep = allergen.lco_report(args.blast, cutoffs, args.cluster_identity, args.cluster_coverage)
        print(f"sequences\t{rep['n_sequences']}\tclusters\t{rep['n_clusters']}")
        for r in rep["recall"]:
            print(
                f"identity>={r['min_identity']:g}\tcoverage>={r['min_coverage']:g}\t{r['recovered']}/{r['n']}"
            )
        return None
    rows, summary = panel_check(args.calls_long, args.panel)
    for r in rows:
        print(
            "\t".join(
                [r["protein"], r["call"], r["variant"], r["expected"], r["observed"], r["verdict"]]
            )
        )
    print(json.dumps(summary), file=sys.stderr)
    return None


def main(argv=None):
    try:
        run(build_parser().parse_args(argv))
    except (ValueError, KeyError, TaxonError, OSError) as err:
        print(f"cellsurface_sorting_hat_calibrate: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

In `src/cellsurface_sorting_hat/categories.yaml`, under `evidence:`, add one line after `- expression.log2fc`:

```yaml
  - tm.n_tm
```

In `pyproject.toml`, `[project.scripts]`, add:

```toml
cellsurface_sorting_hat_module = "cellsurface_sorting_hat.modules.cli:main"
cellsurface_sorting_hat_calibrate = "cellsurface_sorting_hat.calibration.cli:main"
```

In `tests/surface_glyco/test_package.py`, add the same two entries to the `scripts` dictionary that `test_pyproject_names_match_the_decision` compares (Plan 1 added `cellsurface_sorting_hat` there).

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/modules/test_module_cli.py tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py -q`
Expected: `19 passed`. Then the whole suite: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q` Expected: `226 passed, 7 skipped` (the 7 skips are the shellcheck tests; `shellcheck` is not installed on the login node).

- [ ] **Step 5: Lint and commit**

```bash
ruff check src tests/cellsurface_sorting_hat && ruff format --check src tests/cellsurface_sorting_hat
git branch --show-current   # must print sorting-hat-modules
git add src/cellsurface_sorting_hat/modules/cli.py src/cellsurface_sorting_hat/calibration/panel.py src/cellsurface_sorting_hat/calibration/cli.py src/cellsurface_sorting_hat/categories.yaml pyproject.toml tests/surface_glyco/test_package.py tests/cellsurface_sorting_hat/modules/test_module_cli.py tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py
git commit -m "feat(sorting-hat): module and calibration commands

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 10: Job scripts

**Files:**
- Create: `scripts/sorting_hat/signalp_gpu.sbatch`, `pfam_hmmsearch.sbatch`, `repeats.sbatch`, `blast_allergen.sbatch`, `tmhmm.sbatch`, `fetch_proteomes.sh`, `submit_modules.sh`
- Test: `tests/cellsurface_sorting_hat/test_data_and_scripts.py`

**Interfaces:**
- Consumes: `load_family_table` (Task 2) and the data files of Tasks 2 and 7.
- Produces: raw tool outputs under `$WORKDIR/raw/{signalp,pfam,repeats,allergen,tmhmm}/`. `pfam_hmmsearch.sbatch` stops if the model names in the database differ from the family table or if `PFAM_RELEASE` is not part of the resolved database path.

- [ ] **Step 1: Write the failing tests**

Create `tests/cellsurface_sorting_hat/test_data_and_scripts.py` with exactly this content:

```python
"""The shipped family table, the Phase C species table and the job scripts."""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest

from cellsurface_sorting_hat.modules.pfam import load_family_table

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = sorted((ROOT / "scripts" / "sorting_hat").glob("*"))


def test_shipped_family_table_loads_and_starts_with_every_family_inactive():
    families = load_family_table(ROOT / "data" / "sorting_hat" / "family_table.tsv")
    assert len(families) == 10
    assert not any(
        f.active for f in families
    )  # a family is made active only after its specificity test
    assert {f.module for f in families} == {"pfam_adhesion", "pfam_allergen"}


def test_phasec_species_table_has_one_row_per_species():
    with open(ROOT / "data" / "sorting_hat" / "phasec_set_species.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len({r["scientific_name"] for r in rows}) == len(rows) == 6


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_scripts_are_valid_bash_and_follow_the_site_rules(script):
    text = script.read_text()
    assert subprocess.run(["bash", "-n", str(script)], capture_output=True).returncode == 0
    assert "BASH_SOURCE" not in text  # breaks under sbatch
    if script.suffix == ".sbatch":
        assert "SCRATCH:?" in text and "#SBATCH" in text and "set -euo pipefail" in text
        assert "/tmp" not in text


def test_there_is_a_script_for_every_job_that_submit_modules_starts():
    names = {p.stem for p in SCRIPTS if p.suffix == ".sbatch"}
    text = (ROOT / "scripts" / "sorting_hat" / "submit_modules.sh").read_text()
    for job in text.split("for job in ", 1)[1].split(";", 1)[0].split():
        assert job in names


@pytest.mark.skipif(shutil.which("shellcheck") is None, reason="shellcheck is not installed")
@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_shellcheck(script):
    assert (
        subprocess.run(["shellcheck", "-S", "warning", str(script)], capture_output=True).returncode
        == 0
    )


def test_the_pfam_job_checks_model_names_against_the_family_table():
    text = (ROOT / "scripts" / "sorting_hat" / "pfam_hmmsearch.sbatch").read_text()
    assert "differ from the family table" in text and "hmmfetch -f" in text and "--cut_ga" in text
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_data_and_scripts.py -q`
Expected: failures: the script glob finds no files, and `test_there_is_a_script_for_every_job_that_submit_modules_starts` cannot read `submit_modules.sh`.

- [ ] **Step 3: Write the scripts**

Make the files executable (`chmod +x scripts/sorting_hat/*`). Each job script takes its inputs from exported variables, uses `$SCRATCH`, and copies results back with a temporary name.

Create `scripts/sorting_hat/signalp_gpu.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# SignalP 6 (GPU build, fast mode) on one proteome FASTA. Output: $WORKDIR/raw/signalp/prediction_results.txt
# Submit: sbatch --export=ALL,FASTA=/path/p.faa,WORKDIR=/path/work signalp_gpu.sbatch
# Measured rate on gpu12: about 192 proteins/s (analysis/step1_compare/jobs/j1_features.sh, job 29280458).
# The GPU build refuses to run without a CUDA device. Do not use the CPU build (about 2 s/protein).
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=1:00:00
#SBATCH -J csh_signalp
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_signalp.$$"
mkdir -p "$TMP" "$WORKDIR/raw/signalp"
module load signalp/6-gpu
command -v signalp6 >/dev/null || { echo "FATAL: signalp6 not on PATH" >&2; exit 1; }
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
signalp6 --fastafile "$TMP/in.fasta" --organism eukarya --output_dir "$TMP/out" \
  --format none --mode fast --write_procs 4 --torch_num_threads 8
cp "$TMP/out/prediction_results.txt" "$WORKDIR/raw/signalp/.tmp.prediction_results.txt"
mv "$WORKDIR/raw/signalp/.tmp.prediction_results.txt" "$WORKDIR/raw/signalp/prediction_results.txt"
signalp6 --version > "$WORKDIR/raw/signalp/version.txt" 2>&1 || true
rm -rf "$TMP"
```

Create `scripts/sorting_hat/pfam_hmmsearch.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# Search the Pfam models of the family table with hmmsearch --cut_ga (as analysis/cys_candidates/01_known_family_hmm.sh).
# Output: $WORKDIR/raw/pfam/domtbl.txt and provenance.json (path, resolved path, sha256 of Pfam-A.hmm).
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=...,FAMILY_TABLE=...,PFAM_RELEASE=38.2 pfam_hmmsearch.sbatch
# PFAM_RELEASE is stated by the person who submits; it is checked against the file name of the resolved path.
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH -J csh_pfam
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
: "${FAMILY_TABLE:?export FAMILY_TABLE (data/sorting_hat/family_table.tsv) before sbatch}"
: "${PFAM_RELEASE:?export PFAM_RELEASE, for example 38.2}"
PFAM_HMM="${PFAM_HMM:-/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_pfam.$$"
OUT="$WORKDIR/raw/pfam"
mkdir -p "$TMP" "$OUT"
module load hmmer/3.4
RESOLVED="$(readlink -f "$PFAM_HMM")"
case "$RESOLVED" in *"$PFAM_RELEASE"*) ;; *) echo "FATAL: $RESOLVED does not contain release $PFAM_RELEASE" >&2; exit 1 ;; esac
tail -n +2 "$FAMILY_TABLE" | cut -f1 | sort -u > "$TMP/keys.txt"
hmmfetch -f "$PFAM_HMM" "$TMP/keys.txt" > "$TMP/families.hmm"
N=$(grep -c '^NAME' "$TMP/families.hmm" || true)
[ "$N" -eq "$(wc -l < "$TMP/keys.txt")" ] || { echo "FATAL: hmmfetch found $N of $(wc -l < "$TMP/keys.txt") models" >&2; exit 1; }
awk '/^NAME/{n=$2} /^ACC/{split($2,a,"."); print a[1]"\t"n}' "$TMP/families.hmm" | sort > "$TMP/fetched.txt"
tail -n +2 "$FAMILY_TABLE" | cut -f1,2 | sort > "$TMP/table.txt"
diff -q "$TMP/fetched.txt" "$TMP/table.txt" > /dev/null || { echo "FATAL: model names in the database differ from the family table" >&2; diff "$TMP/fetched.txt" "$TMP/table.txt" >&2; exit 1; }
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
hmmsearch --cut_ga --cpu "${SLURM_CPUS_ON_NODE:-2}" --noali -o /dev/null --domtblout "$TMP/domtbl.txt" "$TMP/families.hmm" "$TMP/in.fasta"
cp "$TMP/domtbl.txt" "$OUT/.tmp.domtbl.txt" && mv "$OUT/.tmp.domtbl.txt" "$OUT/domtbl.txt"
printf '{"pfam_hmm": "%s", "resolved": "%s", "release": "%s", "sha256": "%s", "hmmer": "%s"}\n' \
  "$PFAM_HMM" "$RESOLVED" "$PFAM_RELEASE" "$(sha256sum "$RESOLVED" | cut -d' ' -f1)" "$(hmmsearch -h | sed -n 2p)" > "$OUT/provenance.json"
rm -rf "$TMP"
```

Create `scripts/sorting_hat/repeats.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# Run both repeat detectors on one proteome FASTA. Output: $WORKDIR/raw/repeats/{repeat02,repeat14}.tsv
# Submit: sbatch --export=ALL,PROJ_ROOT=...,FASTA=...,WORKDIR=... repeats.sbatch
# Run time is not measured for a whole proteome; record it from the job (sacct) in the run report.
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -J csh_repeats
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT (repository root) before sbatch}"
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_repeats.$$"
OUT="$WORKDIR/raw/repeats"
mkdir -p "$TMP" "$OUT"
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
/usr/bin/python3.12 "$PROJ_ROOT/analysis/cocci_repeats/02_repeat_profile.py" "$TMP/in.fasta" --out "$TMP/repeat02.tsv"
/usr/bin/python3.12 "$PROJ_ROOT/analysis/cocci_repeats/14_repeat_detect_general.py" "$TMP/in.fasta" --out "$TMP/repeat14.tsv"
for m in repeat02 repeat14; do cp "$TMP/$m.tsv" "$OUT/.tmp.$m.tsv" && mv "$OUT/.tmp.$m.tsv" "$OUT/$m.tsv"; done
rm -rf "$TMP"
```

Create `scripts/sorting_hat/blast_allergen.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# BLASTP of one proteome against the WHO/IUIS fungal allergen sequences.
# Output: $WORKDIR/raw/allergen/blast.tsv (columns: qseqid sseqid pident length qlen slen bitscore evalue)
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=...,ALLERGEN_FASTA=... blast_allergen.sbatch
# -max_target_seqs 200 is larger than the database (116 sequences), so no hit is lost to the
# early-stop behaviour of -max_target_seqs.
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH -J csh_blast
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
: "${ALLERGEN_FASTA:?export ALLERGEN_FASTA (made by build_allergen_fasta) before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_blast.$$"
OUT="$WORKDIR/raw/allergen"
mkdir -p "$TMP" "$OUT"
module load ncbi-blast/2.14.0+
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
makeblastdb -in "$ALLERGEN_FASTA" -dbtype prot -out "$TMP/allergens" -logfile /dev/null
blastp -query "$TMP/in.fasta" -db "$TMP/allergens" -evalue 1 -max_target_seqs 200 \
  -outfmt "6 qseqid sseqid pident length qlen slen bitscore evalue" -num_threads "${SLURM_CPUS_ON_NODE:-2}" \
  -out "$TMP/blast.tsv"
cp "$TMP/blast.tsv" "$OUT/.tmp.blast.tsv" && mv "$OUT/.tmp.blast.tsv" "$OUT/blast.tsv"
blastp -version | head -1 > "$OUT/version.txt"
rm -rf "$TMP"
```

Create `scripts/sorting_hat/tmhmm.sbatch` with exactly this content:

```bash
#!/bin/bash -l
# TMHMM 2.0c on one proteome FASTA. Output: $WORKDIR/raw/tmhmm/tmhmm.tsv (columns as analysis/cocci_spherule/01b_tmhmm.sh)
# Submit: sbatch --export=ALL,FASTA=...,WORKDIR=... tmhmm.sbatch
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH -J csh_tmhmm
set -euo pipefail
: "${FASTA:?export FASTA before sbatch}"
: "${WORKDIR:?export WORKDIR before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/csh_tmhmm.$$"
OUT="$WORKDIR/raw/tmhmm"
mkdir -p "$TMP" "$OUT"
module load tmhmm/2.0c
case "$FASTA" in *.gz) zcat "$FASTA" > "$TMP/in.fasta" ;; *) cp "$FASTA" "$TMP/in.fasta" ;; esac
tmhmm -short "$TMP/in.fasta" > "$TMP/tmhmm.out"
awk 'BEGIN{OFS="\t"; print "protein_id","len","exp_aa","first60","pred_hel","topology"}
  {for(i=2;i<=NF;i++){sub(/^[^=]*=/,"",$i)} print $1,$2,$3,$4,$5,$6}' "$TMP/tmhmm.out" > "$TMP/tmhmm.tsv"
cp "$TMP/tmhmm.tsv" "$OUT/.tmp.tmhmm.tsv" && mv "$OUT/.tmp.tmhmm.tsv" "$OUT/tmhmm.tsv"
rm -rf "$TMP"
```

Create `scripts/sorting_hat/fetch_proteomes.sh` with exactly this content:

```bash
#!/bin/bash
# Download the two extra A. fumigatus proteomes (decision D13) and write provenance records.
# Run on the login node: WORKDIR=/path/work bash fetch_proteomes.sh
# A1163 (CEA10, FGSC A1163, CBS 144.89): UniProt proteome UP000001699 (9,942 proteins on 2026-10-04).
# W72310: NCBI GCA_040167795.1 (UCR_Afum_W72310_1.0).
set -euo pipefail
: "${WORKDIR:?export WORKDIR}"
DEST="$WORKDIR/proteomes"
mkdir -p "$DEST"
curl -fsSL "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001699&format=fasta" -o "$DEST/.tmp.A1163.faa"
mv "$DEST/.tmp.A1163.faa" "$DEST/Afum_A1163_UP000001699.faa"
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/167/795/GCA_040167795.1_UCR_Afum_W72310_1.0"
curl -fsSL "$BASE/GCA_040167795.1_UCR_Afum_W72310_1.0_protein.faa.gz" -o "$DEST/.tmp.W72310.faa.gz"
mv "$DEST/.tmp.W72310.faa.gz" "$DEST/Afum_W72310_GCA_040167795.1.faa.gz"
echo "downloaded; now write provenance (see Task 12 of the plan)"
```

Create `scripts/sorting_hat/submit_modules.sh` with exactly this content:

```bash
#!/bin/bash
# Submit the module jobs for one proteome. Usage:
#   PROJ_ROOT=... FASTA=... WORKDIR=... FAMILY_TABLE=... PFAM_RELEASE=38.2 ALLERGEN_FASTA=... bash submit_modules.sh
# Prints one job ID per line. Wait for the jobs (squeue / sacct), then run the converters
# (module commands of Task 13) and the core command. All paths come from PROJ_ROOT.
set -euo pipefail
: "${PROJ_ROOT:?export PROJ_ROOT}"
: "${FASTA:?export FASTA}"
: "${WORKDIR:?export WORKDIR}"
: "${FAMILY_TABLE:?export FAMILY_TABLE}"
: "${PFAM_RELEASE:?export PFAM_RELEASE}"
: "${ALLERGEN_FASTA:?export ALLERGEN_FASTA}"
S="$PROJ_ROOT/scripts/sorting_hat"
mkdir -p "$WORKDIR/logs"
for job in signalp_gpu pfam_hmmsearch repeats blast_allergen tmhmm; do
  sbatch --parsable --export=ALL -o "$WORKDIR/logs/$job.%j.log" -e "$WORKDIR/logs/$job.%j.log" "$S/$job.sbatch"
done
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat/test_data_and_scripts.py -q`
Expected: `11 passed, 7 skipped` (skips are `shellcheck` tests; run `shellcheck -S warning scripts/sorting_hat/*` on a machine that has it).

- [ ] **Step 5: Mutation checks for Tasks 1 to 9 (confirm the tests can fail)**

Make each change in turn, run `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q`, confirm that the named tests fail, then undo the change with `git checkout -- <file>`.

| Change | File | Tests that fail (measured) |
|---|---|---|
| In `status_from_measure`, replace `if spec and enough and narrow:` by `if enough and narrow:` | `calibration/measure.py` | `test_status_rule_follows_phase_c_and_needs_negatives_for_an_estimate` |
| In the `signalp` command, set `artefact_digest=_file_digest(args.results)` | `modules/cli.py` | `test_the_signalp_module_identity_does_not_depend_on_the_proteome` |
| In `pfam_rows`, delete the two lines of the `no_tm` condition | `modules/pfam.py` | `test_no_tm_condition_drops_a_domain_in_a_protein_with_transmembrane_helices`, `test_pfam_command_applies_the_no_tm_condition_from_the_tm_table` |
| In `make_entry`, delete the `if cap == SMOKE and status == ESTIMATED:` block | `calibration/measure.py` | `test_a_leakage_cap_limits_an_estimate_to_smoke`, `test_leakage_other_than_none_caps_an_estimate_at_smoke` |
| In `best_hit_outside_cluster`, delete the `cluster_of.get(q) == cluster_of.get(s)` test | `modules/allergen.py` | `test_leave_cluster_out_recall_ignores_hits_inside_the_cluster`, `test_allergen_lco_command_prints_recall_per_cutoff` |
| In `compare_runs`, set `other = pid` | `calibration/stability.py` | `test_only_identical_sequences_are_compared_by_sha256_not_by_id` |

- [ ] **Step 6: Commit**

```bash
git branch --show-current   # must print sorting-hat-modules
git add scripts/sorting_hat tests/cellsurface_sorting_hat/test_data_and_scripts.py
git commit -m "feat(sorting-hat): SLURM job scripts for the module tools

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```


### Task 11: Proteomes and taxonomy (HPCC)

**Files:**
- Create (not committed): `_workdir/sorting_hat/proteomes/*`, `_workdir/sorting_hat/taxdump/*`, `_workdir/sorting_hat/provenance/*.json`

**Interfaces:**
- Consumes: `write_provenance` (Task 8), `fetch_proteomes.sh` (Task 10).
- Produces: five proteome FASTA files with provenance records, and `nodes.dmp`/`names.dmp`.

- [ ] **Step 1: Download the taxonomy and the two extra proteomes**

```bash
export PROJ_ROOT=$PWD WORKDIR=$PWD/_workdir/sorting_hat
mkdir -p "$WORKDIR/taxdump" && cd "$WORKDIR/taxdump"
curl -fsSL https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -o taxdump.tar.gz && tar xzf taxdump.tar.gz nodes.dmp names.dmp
cd "$PROJ_ROOT" && bash scripts/sorting_hat/fetch_proteomes.sh
ls "$WORKDIR/proteomes"
```

Expected: `Afum_A1163_UP000001699.faa` and `Afum_W72310_GCA_040167795.1.faa.gz`.

- [ ] **Step 2: Check that each tested taxon is above the proteome's taxon**

Status applies to a protein only if its taxon is a descendant of a tested taxon. Check that the strain-level taxa you will pass lie under the species-level taxa that Task 13 resolves:

```bash
/usr/bin/python3.12 - <<'E'
import os
from cellsurface_sorting_hat.taxonomy import Lineage
from cellsurface_sorting_hat.calibration.phasec import read_names, species_taxid
W = os.environ["WORKDIR"] + "/taxdump"
lin, names = Lineage.from_nodes_dmp(W + "/nodes.dmp"), read_names(W + "/names.dmp")
afum = species_taxid(names, "Aspergillus fumigatus")
for label, taxon in {"Af293": 330879, "W72310": 746128}.items():
    print(label, taxon, lin.is_descendant_or_self(taxon, afum))
E
```

Expected: both print `True`. If `False`, stop: a measurement of *A. fumigatus* would not apply, and the run report would show `unvalidated`. Read the strain taxon of A1163 from the UniProt proteome record (`https://rest.uniprot.org/proteomes/UP000001699.json`, field `taxonomy.taxonId`) and check it the same way; use that ID for A1163 in Task 14.

- [ ] **Step 3: Write provenance records and check the counts**

```bash
/usr/bin/python3.12 - <<'E'
import os, datetime
from cellsurface_sorting_hat.proteomes import write_provenance
W = os.environ["WORKDIR"]; P = os.environ["PROJ_ROOT"]
today = datetime.date.today().isoformat()
os.makedirs(W + "/provenance", exist_ok=True)
sets = [
 ("Cimm_RS", "NCBI RefSeq GCF_000149335.2", P + "/_workdir/cocci_spherule/ref/GCF_000149335.2_ASM14933v2_protein.faa.gz", 9910, 9910),
 ("Afum_Af293", "Fungi_5k input", "/bigdata/stajichlab/shared/projects/Fungi_5k/input/Aspergillus_fumigatus_Af293.proteins.fa", 9161, 9161),
 ("Scer_S288C", "SGD orf_trans_all.fasta.gz", P + "/_workdir/step1_compare/downloads/orf_trans_all.fasta.gz", 6000, 7500),
 ("Afum_A1163", "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001699", W + "/proteomes/Afum_A1163_UP000001699.faa", 9000, 10500),
 ("Afum_W72310", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/167/795/GCA_040167795.1_UCR_Afum_W72310_1.0/", W + "/proteomes/Afum_W72310_GCA_040167795.1.faa.gz", 9000, 11000),
]
for name, source, fasta, lo, hi in sets:
    r = write_provenance(f"{W}/provenance/{name}.json", name, source, fasta, today, lo, hi)
    print(name, r["n_proteins"], r["n_unique_sequences"], r["n_invalid"])
E
```

Expected: a line per proteome. RS prints 9,910 and Af293 prints 9,161 (counted 2026-10-04 and 2026-10-05). The ranges for S288C, A1163 and W72310 are plausibility ranges chosen from the *A. fumigatus* proteomes already in the repository (9,161 and 9,942) and the S288C step 1 set (about 6,700); they are not measured. Record the three counts you get. A `ProteomeError` means the file is not the expected proteome: stop and find out why.

- [ ] **Step 4: Commit nothing (the files are under `_workdir`, which is git-ignored)**

Copy the five provenance JSON files to `docs/reports/data/sorting_hat/provenance/` and commit them with the report in Task 15.

### Task 12: First end-to-end run on A. fumigatus Af293 (HPCC)

**Files:**
- Create (not committed): `_workdir/sorting_hat/Afum_Af293/{raw,modules,status,out,logs}`

**Interfaces:**
- Consumes: Tasks 1 to 11 and Plan 1 (`cellsurface_sorting_hat`).
- Produces: `out/calls.long.tsv.gz`, `out/report.md`, and the job resource table (Step 4).

- [ ] **Step 1: Install the package in a virtual environment (Python 3.12)**

```bash
/usr/bin/python3.12 -m venv "$SCRATCH/csh-venv"
"$SCRATCH/csh-venv/bin/pip" install -e . --no-deps && "$SCRATCH/csh-venv/bin/pip" install pyyaml numpy
export PATH="$SCRATCH/csh-venv/bin:$PATH"
cellsurface_sorting_hat --help | head -2 && cellsurface_sorting_hat_module --help | head -2 && cellsurface_sorting_hat_calibrate --help | head -2
```

Expected: three usage lines. (`$SCRATCH` is node-local on a login node only if the shell has it; use a job shell or `/scratch/$USER` if it is empty.)

- [ ] **Step 2: Build the allergen database FASTA and submit the module jobs**

```bash
export PROJ_ROOT=$PWD WORKDIR=$PWD/_workdir/sorting_hat/Afum_Af293
export FASTA=/bigdata/stajichlab/shared/projects/Fungi_5k/input/Aspergillus_fumigatus_Af293.proteins.fa
export FAMILY_TABLE=$PROJ_ROOT/data/sorting_hat/family_table.tsv PFAM_RELEASE=38.2
export ALLERGEN_FASTA=$PROJ_ROOT/_workdir/sorting_hat/iuis_fungal_allergens.faa
/usr/bin/python3.12 -c "from cellsurface_sorting_hat.modules.allergen import build_allergen_fasta as b; print(b('analysis/allergen_scoping/iuis_fungal_isoallergens.tsv', '$ALLERGEN_FASTA'))"
bash scripts/sorting_hat/submit_modules.sh | tee "$WORKDIR/jobs.txt"
```

Expected: the allergen FASTA has 116 sequences; five job IDs print. Wait until `squeue -u $USER` shows none of them. If a job fails, read its log in `$WORKDIR/logs/` before resubmitting; do not resubmit in a loop.

- [ ] **Step 3: Convert the tool outputs to module tables and run the core command**

Use one shell for Steps 2 to 4 and for Task 13 (the variables `T`, `M`, `FASTA`, `WORKDIR` are reused).

```bash
T=$PROJ_ROOT/_workdir/sorting_hat/taxdump
M="cellsurface_sorting_hat_module"
$M signalp --fasta $FASTA --workdir $WORKDIR --results $WORKDIR/raw/signalp/prediction_results.txt --signalp-version "$(head -1 $WORKDIR/raw/signalp/version.txt)"
$M tm --fasta $FASTA --workdir $WORKDIR --table $WORKDIR/raw/tmhmm/tmhmm.tsv
$M pfam --fasta $FASTA --workdir $WORKDIR --domtbl $WORKDIR/raw/pfam/domtbl.txt --family-table $FAMILY_TABLE --pfam-release 38.2 \
  --pfam-sha256 "$(python3.12 -c "import json;print(json.load(open('$WORKDIR/raw/pfam/provenance.json'))['sha256'])")" --sp-module step1_rule@R0 --tm-module tm
$M repeat02 --fasta $FASTA --workdir $WORKDIR --table $WORKDIR/raw/repeats/repeat02.tsv
$M repeat14 --fasta $FASTA --workdir $WORKDIR --table $WORKDIR/raw/repeats/repeat14.tsv
$M allergen --fasta $FASTA --workdir $WORKDIR --blast $WORKDIR/raw/allergen/blast.tsv --allergen-fasta $ALLERGEN_FASTA --blast-version "$(head -1 $WORKDIR/raw/allergen/version.txt)"
cellsurface_sorting_hat --fasta $FASTA --taxon 330879 --taxdump $T/nodes.dmp --workdir $WORKDIR --out $WORKDIR/out
head -40 $WORKDIR/out/report.md
```

Expected: exit codes 0. The report lists the modules `repeat02`, `repeat14`, `step1_rule@R0`, `tm`, `pfam_adhesion`, `pfam_allergen`, `allergen_homology` as `ok`, and `antigen_lookup` as absent from the run states (the antigen lookup is run for RS only; its absence is listed under "Modules the rules need and that were not found"). Every status is `unvalidated` until Task 13.

- [ ] **Step 4: Compare with the numbers already measured, and record the job resources**

```bash
zcat $WORKDIR/out/calls.long.tsv.gz | awk -F'\t' '$2=="allergen_homolog_hit" && $4=="called"' | wc -l
zcat $WORKDIR/out/calls.long.tsv.gz | awk -F'\t' '$2=="allergen_candidate" && $4=="called"' | wc -l
sacct -j "$(paste -sd, $WORKDIR/jobs.txt)" --format=JobID,JobName%20,Partition,Elapsed,MaxRSS,AllocCPUS,State
```

Expected: 105 and 41 (measured 2026-10-05 with BLAST 2.14.0+ on this FASTA; the engine's candidate call also needs the Pfam rule, which has no active allergen family yet). A different number is a finding: find the cause (BLAST version, FASTA, filter) before going on. Write the `sacct` table into the report (Task 15); these are the first measured times for the whole-proteome Pfam, repeat and TMHMM jobs.

- [ ] **Step 5: No commit** (outputs are under `_workdir`). Keep `jobs.txt` and the `sacct` output for the report.

### Task 13: Calibrations that can be measured today (HPCC)

**Files:**
- Create (not committed): `$WORKDIR/status/*.json`, `_workdir/sorting_hat/calibration/*`
- Modify: `data/sorting_hat/family_table.tsv` (only for a family that passed review and has a sign-off)

**Interfaces:**
- Consumes: Tasks 6, 7, 9 and the Phase C `metrics.json`.
- Produces: the status source of `step1_rule@R0`, the allergen leave-cluster-out table, one specificity table per Pfam family.

- [ ] **Step 1: Status source for rule R0 from the Phase C metrics**

```bash
cellsurface_sorting_hat_calibrate phasec --workdir $WORKDIR --metrics _workdir/step1_compare/phasec/metrics.json \
  --set-species data/sorting_hat/phasec_set_species.tsv --names-dmp $T/names.dmp
cellsurface_sorting_hat --fasta $FASTA --taxon 330879 --taxdump $T/nodes.dmp --workdir $WORKDIR --out $WORKDIR/out
grep -A12 "Module calibration" $WORKDIR/out/report.md
```

Expected: `step1_rule@R0` for taxon 330879 reports status `estimated`, calibration set `S3-Eurotiomycetes:clade`, 128 positives, 208 negatives, sensitivity 0.727 [0.640, 0.802] and specificity 0.990 [0.975, 1.000] (values in `metrics.json`, read 2026-10-05). Other modules show `not measured`. The `surface_glycoprotein[R0]` rows now carry status `estimated` with basis `step1_rule@R0:taxon:<A. fumigatus species ID>`.

- [ ] **Step 2: Allergen leave-cluster-out recall (sensitivity only)**

```bash
module load ncbi-blast/2.14.0+
mkdir -p _workdir/sorting_hat/calibration && cd _workdir/sorting_hat/calibration
makeblastdb -in $ALLERGEN_FASTA -dbtype prot -out alg -logfile /dev/null
blastp -query $ALLERGEN_FASTA -db alg -evalue 1 -max_target_seqs 200 -outfmt "6 qseqid sseqid pident length qlen slen bitscore evalue" -out iuis_vs_iuis.tsv
cd $PROJ_ROOT && cellsurface_sorting_hat_calibrate allergen-lco --blast _workdir/sorting_hat/calibration/iuis_vs_iuis.tsv --cutoffs 35:50,50:70,70:80
```

Expected: first line `sequences 116 clusters <n>`; then one line per cutoff with `recovered/116`. Record all four numbers. This says how many known allergens a homology rule finds when no relative within 40% identity and 70% coverage is in the reference set. It is not a specificity and writes no status entry.

- [ ] **Step 3: Pfam family specificity tests (one per family)**

For each family in `data/sorting_hat/family_table.tsv`, run `pfam_hmmsearch.sbatch` (with the family table) on the proteomes of four clades (*A. fumigatus* Af293, *S. cerevisiae* S288C, *C. immitis* RS, *C. neoformans* H99; their FASTA paths are in `analysis/step1_compare/species.tsv` and Task 11) and run TMHMM (`tmhmm.sbatch`) on the same files. Make `members.tsv` (TAB separated, header `pfam_acc`, `protein_id`) from the curated tables in `data/curated/` and the literature (the proteins that are known members of the family by function, not because Pfam found a domain in them). Then, for each family and proteome:

```bash
cellsurface_sorting_hat_calibrate pfam-specificity --family PF05730 --domtbl <proteome>/raw/pfam/domtbl.txt \
  --members members.tsv --universe-fasta <proteome FASTA> --tm-table <proteome>/raw/tmhmm/tmhmm.tsv > <proteome>.PF05730.specificity.tsv
```

A person reads every `nonmember_hit` line (annotation, length, domain architecture, TMHMM helices) and writes one line per hit in the `specificity_note` column of the family table. CFEM hits with helices are expected (Pth11-like receptors); the `no_tm` condition removes them.

Expected: for each family a table with `tp`, `fp`, `fn`, `tn` and the list of non-member hits. A family becomes `active` only when the owner signs off (decision C4): set `active = yes`, `active_by`, `active_date` in the family table and commit that change alone with the table of non-member hits as evidence.

- [ ] **Step 4: Truth tables for the other modules (blocked on owner decisions C1 to C3)**

When an owner decision provides a truth table `truth.tsv` (`id`, `label` 0 or 1, `cluster`), the calibration command is:

```bash
cellsurface_sorting_hat_calibrate truth --workdir $WORKDIR --module <module> --calls-long <out>/calls.long.tsv.gz \
  --call <call name> --truth truth.tsv --calibration-set <name> --taxa <species-level taxon IDs> --leakage <none|partial|tuned_on_truth|unknown> --notes "<what the set is>"
```

The run that produced `calls.long.tsv.gz` must contain the truth proteins. Use `--leakage tuned_on_truth` for `antigen_candidate` against the four anchors. Do not write an entry for a module whose truth set has no negatives (the status stays at most `smoke` and the report shows "not measured" for specificity).

- [ ] **Step 5: Commit the family table change and the specificity tables only**

Copy the `*.specificity.tsv` files to `docs/reports/data/sorting_hat/specificity/` first.

```bash
git branch --show-current   # must print sorting-hat-modules
git add data/sorting_hat/family_table.tsv docs/reports/data/sorting_hat/specificity
git commit -m "analysis(sorting-hat): Pfam specificity tables; family table sign-offs

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

### Task 14: C. immitis RS and the stability of A. fumigatus strains (HPCC)

**Files:**
- Create (not committed): `_workdir/sorting_hat/{Cimm_RS,Afum_A1163,Afum_W72310,Scer_S288C}/...`

**Interfaces:**
- Consumes: Tasks 12 and 13; `compare_runs` (Task 7); `panel` command (Task 9).
- Produces: runs for the other four proteomes, a strain comparison, and the report-only panel check.

- [ ] **Step 1: Run each remaining proteome the way Task 12 does**

Set `WORKDIR`, `FASTA` and the taxon for each: *C. immitis* RS (`FASTA` is the RefSeq protein file, taxon 246410), *S. cerevisiae* S288C (`orf_trans_all.fasta.gz`, taxon 559292), A1163 (taxon read in Task 11), W72310 (taxon 746128). For RS also build the lookup modules. Run the Cys-rich pipeline of `analysis/cys_candidates/README.md` on the RefSeq FASTA first (its IDs must be the `XP_...` IDs of the FASTA, or the lookup finds nothing):

```bash
$M antigen --fasta $FASTA --workdir $WORKDIR --taxon 246410 --ranking analysis/cocci_antigens/cocci_antigen_ranking.tsv --protein-map _workdir/cocci_spherule/ref/protein_map.tsv
$M expression --fasta $FASTA --workdir $WORKDIR --taxon 246410 --table analysis/cocci_spherule/spherule_surface_table.tsv.gz --protein-map _workdir/cocci_spherule/ref/protein_map.tsv
$M cys --fasta $FASTA --workdir $WORKDIR --taxon 246410 --candidates _workdir/cys_candidates/candidates.tsv.gz
```

Record, for RS: the number of proteins with `antigen_lookup` state `ok`, `not_in_reference` and `not_applicable` (the match rate of the ID map), and how many of the 14 Tier 1 and 45 Tier 2 candidates have `antigen_candidate` called. Also record how many of the 9,910 RS proteins have a signal peptide under rule R0 (expected 460, counted from `analysis/cocci_repeats/signalp/CimmitisRS_FungiDB/prediction_results.txt` on 2026-10-05).

- [ ] **Step 2: Compare the three *A. fumigatus* runs**

```bash
/usr/bin/python3.12 - <<'E'
import json, os
from cellsurface_sorting_hat.calibration.stability import compare_runs
W = os.environ["PROJ_ROOT"] + "/_workdir/sorting_hat"
for other in ("Afum_A1163", "Afum_W72310"):
    r = compare_runs(f"{W}/Afum_Af293/out", f"{W}/{other}/out")
    print(other, r["n_a"], r["n_b"], r["n_identical_sequences"])
    for key, c in sorted(r["per_call"].items()):
        print("  ", key, dict(c))
E
```

Expected: the number of identical sequences and, per call, how many agree. Differences in calls for identical sequences come from the taxon or status, not the sequence; any `differ` for an unchanged status is a bug to find. Proteins without an identical sequence in the other strain are not compared: report how many.

- [ ] **Step 3: Report-only panel check**

Make `docs/reports/data/sorting_hat/panel.tsv` (columns `protein`, `call`, `variant`, `expected`, `source`) from the panel of spec section 4: SOWgp (`adhesion_repeat`, `antigen_candidate`), Als1, Flo11 (`adhesion_repeat`), Ag2/PRA, Rbt5 (`adhesion_domain`), RodA, CalA, Gel1 (`surface_glycoprotein` only), Asp f 1, Asp f 2, Asp f 34 (`allergen_candidate`), Hsp60 (`known_miss`). Each row needs the protein ID of the proteome it is in and a source (PMID or table row). Then:

```bash
cellsurface_sorting_hat_calibrate panel --calls-long <out>/calls.long.tsv.gz --panel docs/reports/data/sorting_hat/panel.tsv
```

Expected: one line per panel row and a summary of `agree`, `disagree`, `known_miss` on stderr. Nothing is gated on the result; list every disagreement in the report with the reason.

- [ ] **Step 4: No commit** (outputs under `_workdir`); commit `panel.tsv` with the report.

### Task 15: Calibration and run report

**Files:**
- Create: `docs/reports/2026-10-05-sorting-hat-run-and-calibration.md`, `docs/reports/data/sorting_hat/provenance/*.json`, `docs/reports/data/sorting_hat/panel.tsv`

**Interfaces:**
- Consumes: the numbers from Tasks 11 to 14.
- Produces: the report that Plan 3 (acceptance) and the owner read.

- [ ] **Step 1: Write the report**

Sections, each filled from a command above (copy the numbers, do not retype them from memory):
1. Proteomes: the five provenance records (source, date, sha256, counts).
2. Jobs: the `sacct` table (partition, elapsed, MaxRSS) for each module job on Af293.
3. Calibration table: one row per module and calibration set with positives, negatives, Sn and Sp with 95% intervals, status, leakage, and "not measured" where nothing was measured.
4. Allergen leave-cluster-out recall (Task 13 step 2) and the Af293 counts (105, 41, 30 or the numbers you measured).
5. Pfam specificity tables and the non-member hits you read, with the sign-off decision per family.
6. Antigen lookup match rate on RS; Tier 1 and Tier 2 coverage by `antigen_candidate`; a statement that the cut of 15 was chosen after the anchors were seen.
7. Strain stability (Task 14 step 2).
8. Panel check (disagreements and reasons).
9. What is still not measured, by module, and which owner decision (C1 to C4) blocks it.

- [ ] **Step 2: Commit**

```bash
git branch --show-current   # must print sorting-hat-modules
git add docs/reports/2026-10-05-sorting-hat-run-and-calibration.md docs/reports/data/sorting_hat
git commit -m "docs: sorting hat run and calibration report

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

## Self-review against the spec and Plan 1

| Item | Where |
|---|---|
| Module tables and run records per the Plan 1 contract (`id`, `state`, fields; identity of tool, database, parameters) | Tasks 1 to 5, 9 |
| `step1_rule@R0` from SignalP `SP`; `pfam_adhesion`, `pfam_allergen`; `repeat02`, `repeat14`; `allergen_homology`; `antigen_lookup`; `cys_rich`; `expression`; `tm` | Tasks 1 to 5, 9 |
| Family table with `active` sign-off, `signal_peptide` and `no_tm` conditions; specificity test | Task 2, Task 13 |
| Allergen: evidence tier and candidate tier fields; leave-cluster-out recall | Task 4, Task 13 |
| Antigen lookup: exact-taxon applicability, ID map through `protein_map.tsv`, match rate reported | Task 5, Task 14 |
| Status rule needs negatives; leakage cap; Phase C numbers into status sources; Sn/Sp with cluster bootstrap | Tasks 6, 7, 9, 13 |
| `measure` object in status sources shown in the report ("Module calibration" of Plan 1) | Tasks 6, 7, 13 |
| Strain stability by sequence hash | Task 7, Task 14 |
| Proteome downloads with provenance; D13 proteomes | Tasks 8, 10, 11 |
| SLURM job scripts, `$SCRATCH`, no `BASH_SOURCE`, resource table measured | Tasks 10, 12 |
| Run-level check on Af293 on HPCC | Task 12 |
| Spec 3.8 GPU out-of-memory retry, job timeout and preemption handling | **Not implemented.** The jobs have fixed resources and stop on error; the person resubmits after reading the log. Add to Plan 3 if jobs fail in practice. |
| Spec 5 kappa between modules of the same kind (R0 against ML) | Not in this plan: the ML variant has no model. |

Placeholder scan: none. The Pfam specificity test is the tested command `pfam-specificity` (Task 9).

Known gaps, stated so the reviewer can check them:
1. The antigen calibration is circular (cut chosen after the anchors were seen). The plan records this and caps the status; it does not remove it.
2. No negatives exist for the allergen module. The leave-cluster-out run writes no status entry, so `allergen_candidate` stays `unvalidated`; a sensitivity alone could give at most `smoke`.
3. TMHMM, not Phobius, is the TM tool; neither has a calibration here.
4. `hmmsearch --cut_ga` uses the Pfam gathering thresholds; they were not tuned for fungal proteins.
5. Run times of the whole-proteome Pfam, repeat and TMHMM jobs are not measured yet (Task 12 step 4 measures them).
6. The IUIS sequence set is 116 sequences from 31 species; Onygenales is represented by four *Trichophyton* proteases, none from *Coccidioides*.
7. The plan does not add Phobius, an ML step 1 variant, R1/R2 freezing, Fungi_5k batch runs, or a stored golden `calls` file.
