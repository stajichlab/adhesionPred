# Step 1 truth set (Phase A) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the GO truth set for the step 1 surface model: the D1 extractor that reproduces the spec 3.3 counts exactly, sequence attachment, the D8 GPI/TM triage, and the D10 keyword-tier table with test proteins removed.

**Architecture:** A new analysis folder `analysis/step1_compare/` holds small standard-library modules (ontology parser, GAF reader, label rules, truth-table builder, manifest checks, sequence mapping, D8 parser, D10 filter) and five numbered scripts (`00_` to `04_`). Each script reads pinned inputs from a work directory taken from the environment and writes compressed tables there, never to the repository. The truth table stores every evidence code, three label columns (one per evidence policy), the source file hash and date, the clade and the role, so evidence subsets and owner choices are filters on stored columns, not a re-extraction. Direct-evidence truth is an intersection: `label == X and homology_only == "no"`.

**Tech Stack:** Python 3.12 standard library (gzip, csv, hashlib, json, argparse, urllib, dataclasses); pytest; ruff 0.3.5 via pre-commit; GitHub Actions (CPU).

**Spec:** `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md` (sections 2.2, 2.3, 3.1 to 3.4, 4, 7, 11; deliverables D1, D8, D9 extraction half, D10). The counting rule is `/tmp/glyco_spec/d1_count.py` (SHA-256 `25021b6aa9a93849995d3834d9bbb1960eff984657648f920f6b290e458194db`).

## Global Constraints

- Python 3.12 is canonical; CI is CPU only; ruff 0.3.5 via pre-commit, line length 100; every commit passes `pre-commit run --all-files` and `ruff check .`.
- stdlib-only for the extractor and label rules (gzip, csv, hashlib, argparse, json, urllib); no pandas in those modules, so CI needs no extra install; scripts accept plain or .gz input and write compressed output; large intermediates go to the work directory, not the repo.
- Never use `$(dirname "${BASH_SOURCE[0]}")` in any shell script; pass the project root by environment variable. Jobs on SLURM use `$SCRATCH` (node-local) and copy results to /bigdata before ending; this phase needs no SLURM job (downloads are small); say so.
- Run python for any Mycelium script with /usr/bin/python3.12 (global instruction) — only relevant if a Mycelium script is used; you may use /usr/bin/python3.12 for the analysis scripts too.
- Writing style for all prose in the plan: Simplified Technical English (short sentences, one idea per sentence, active voice, no idioms, no flowery language). No inflated claims: never state a count in the plan that you did not compute. Where the real data are needed (SHA-256, counts), compute them now from /tmp/glyco_spec/ and cite the command.
- Commit message trailer: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
- The analysis-tests CI job has `continue-on-error: true`. The `tests/step1_compare` suite is therefore informational in CI: it cannot fail the build. Moving it into a required job is an owner decision; this plan does not change it.
- This plan adds no SLURM job. The downloads total about 83 MB (measured on 2026-10-01). Run every script on a login or interactive node.
- The work directory is `$STEP1_WORKDIR`. If it is not set, it is `<workdir in config/site.yaml>/step1_compare`. The repository root is `$PROJ_ROOT`; if it is not set, `paths.py` uses its own file location (a Python path, not `BASH_SOURCE`).
- Interpreter: `PY=/usr/bin/python3.12`. On this system pytest for that interpreter is in the user site (`~/.local`). Run tests as `$PY -m pytest tests/step1_compare -q` from the repository root.
- Do not edit `/tmp/glyco_spec/d1_count.py`. Do not commit any downloaded file.

### Decisions recorded on 2026-10-01 and how the data keep them reversible

The review fix brief of 2026-10-01 closed four of the five spec section 12 questions and gave three controller rulings (R-A, R-B, R-C). Each choice is a filter on stored columns, so a later change needs no re-extraction.

| Item | Status | Stored columns and filter |
|---|---|---|
| Homology codes as truth | closed: direct-evidence truth is the intersection `label == X and homology_only == "no"`; headline metrics and gates use it; all non-IEA truth (`label == X`) is reported beside it | `label`, `homology_only`; `counts.tsv` `direct_*`. `label_no_homology` is stored for traceability and is not used as truth: recomputing without homology codes moves 83 ambiguous *C. albicans* genes (TDH3, PGK1, ADH1 and others) into P-ext, against Q9 |
| *A. nidulans* in S3 | closed: second Eurotiomycetes test species | `source_id == "Anid_EMENI"`, `role == "test_clade"` |
| H99 or JEC21 | closed: H99 is the reference; a FungiDB H99 GO check and literature curation are separate later work | `Cneo_H99_GOA` `role == "test_clade"`; JEC21 and CRYD1 `role == "alternate_file"` |
| Estimate versus smoke test | closed: "estimate" when the 95% cluster-bootstrap interval for recall has half-width <= 0.10, else "smoke test" | `counts.tsv` per source; applied by the later evaluation plan |
| Basidiomycota truth option | open | `in_clade == "Basidiomycota"`; `Umay_MYCMD` `role == "undecided"` |
| R-A high-throughput internal evidence | D1 rule unchanged; sub-stratum reported | `internal_evidence_htp_only == "yes"` within `label == "ambiguous"`; `counts.tsv` `ambiguous_htp_only` |
| R-B P-gpi | list, not a scored stratum, until `curated_gpi.tsv` has literature rows; filling it is separate curation work | `d8_triage.tsv`, `d8_gpi_outside_pext.tsv` |
| R-C shipping gate | unchanged: no model ships before the Basidiomycota truth is in the evaluation | none in this plan |

## Review Focus

1. **One gene, conflicting evidence rows** (wall IDA, wall IBA, cytosol IBA; or a `NOT` row next to a positive row). Expected: the label follows the spec rule, the table keeps every code, and the gene is excluded from the direct-evidence stratum. Pinned by `test_conflicting_evidence_is_stored_not_collapsed` and `test_direct_stratum_is_an_intersection` (Task 4) and the `NOT` row in `test_gaf_extract_golden` (Task 3).
2. **IDs that differ in case, transcript suffix or isoform suffix between GAF and FASTA** (`spbc21h7.03C` vs `SPBC21H7.03c.1:pep`; `P22146-2` vs `sp|P22146|`), and **one ID twice in a FASTA**. Expected: a match after normalisation; identical duplicates merge; different duplicates stop the run. Pinned by `test_pombase_suffix_and_case_insensitive_match`, `test_uniprot_isoform_and_version_differences`, `test_duplicate_key_with_different_sequence_raises` (Task 7).
3. **A gene with wall and cytosol only at IEA level** (plus a non-IEA plasma-membrane term). Expected: unlabelled, never P-ext, never N-sec. Pinned by `test_truth_set_has_no_iea` and the truth-table row `(set(), {WALL, CYTO})` (Task 4).
4. **A truncated download or an HTML error page saved as `.gaf.gz`.** Expected: the fetch stops with exit 2, keeps no file, and the extractor refuses any input whose SHA-256 differs from the manifest. An HTTP 403 or a broken connection stops with `STOP:` and leaves no `.part` file. Pinned by `test_check_payload_refuses_html_truncated_and_empty`, `test_html_error_page_stops_even_with_update`, `test_strict_hash_mismatch_stops_and_keeps_no_file`, `test_http_403_stops_with_exit_2`, `test_url_error_mid_download_leaves_no_part_file` (Task 5), and `test_extract_stops_on_hash_mismatch`, `test_missing_gaf_stops_with_exit_2` (Task 6).
5. **Label leakage from the keyword tier into held-out species**: the same sequence under a different accession (for example a *C. posadasii* entry identical to SOWgp58 Q8NK60). Expected: removal by exact cleaned-sequence hash, with the matched test protein in the log. Pinned by `test_tc_excludes_test_proteins` (Task 9).

---

## Measured facts this plan relies on (computed on 2026-10-01)

Commands were run with `/usr/bin/python3.12` unless stated.

- `cd /tmp/glyco_spec && /usr/bin/python3.12 d1_count.py` prints exactly the spec 3.3 numbers. A prototype of the code in Tasks 2 to 6 reproduced all 21 count fields that the script prints, for all 10 sources, with zero mismatches (field-by-field comparison against the script output). Its PM-candidate counts (12, 60, 10, 9) equal the spec 2.2 values.
- `sha256sum` of the files in `/tmp/glyco_spec/`: the values are in `manifest.tsv` (Task 5). A fresh download of all 11 GO files on 2026-10-01 gave the same 11 hashes (`diff` of the two `sha256sum` lists: no difference). `curl -sIL` gave HTTP 200 and the `Last-Modified` values in the manifest. The content lengths equal the local file sizes.
- `go-basic.obo` header: `data-version: releases/2026-07-26`. Its HTTP `Last-Modified` is `Sat, 08 Aug 2026 22:27:06 GMT`. The spec date "2026-08-08" is the HTTP date, not the release.
- Protein FASTA (all HTTP 200 on 2026-10-01): SGD `orf_trans_all.fasta.gz` (6,722 records, header field `SGDID:S...`); CGD `C_albicans_SC5314_A22_current_default_protein.fasta.gz` (6,250 records, 6,212 distinct IDs, 38 duplicate records with identical sequences, header = systematic name only, no `CAL` ID; the http URL redirects 301 to https); PomBase `peptide.fa.gz` (5,126 records, IDs with transcript suffix such as `SPAC1002.01.1:pep`).
- UniProt proteomes (queried at `rest.uniprot.org/proteomes/search` on 2026-10-01): UP000002530 Af293 9,647; UP000000560 *A. nidulans* FGSC A4 10,561 (a second proteome, UP000005890, returned 0 entries); UP000010091 H99 7,427; UP000002149 JEC21 6,740; UP000000561 *U. maydis* 6,805. The UniProt `stream` endpoint ended two downloads early: *U. maydis* once, and H99 once (curl exit 92, HTTP/2 stream error). Both files were truncated gzip ("unexpected end of file"). The paginated `search` endpoint with `size=500` gave the complete H99 file in 33 s. The plan uses the paginated endpoint.
- ID mapping on the 2026-10-01 files (prototype of Task 7): SGD 6,052 of 6,056 genes matched; CGD 6,060 of 6,313 (253 unmatched; 72 of them are N-int, with tRNA, rRNA and snoRNA symbols such as `tQ(UUG)6mt`, `RDN58`); PomBase 5,020 of 5,025 (`pombase.gaf.gz`) and 5,008 of 5,008 (`SCHPO-mod.gaf.gz`); all five UniProt-based sources 100%. No unmatched gene is P-ext or ambiguous.
- UniProt 2026_03, reviewed entries with a `GPI-anchor` lipidation feature (query `(organism_id:T) AND (reviewed:true) AND (ft_lipid:GPI-anchor)`): S288C 64 entries, 3 with ECO:0000269 (GAS1 P22146, TIP1 P27654, P53872); *C. albicans* 94, 3 with ECO:0000269; *S. pombe* 17, 0; *A. fumigatus* 22, 0. YPS1, SAG1, CWP2 and CRH1 carry ECO:0000255 (sequence analysis). MSB2 and HKR1 (both reviewed) have one TRANSMEM feature each, evidence ECO:0000255.
- In `sgd.gaf.gz`, GAS1 has fungal-type cell wall (IDA), mitochondrion (HDA) and nuclear periphery (IDA, part of nucleus). The D1 rule therefore labels GAS1 ambiguous, not P-ext; its internal evidence is not high-throughput only.
- Direct-evidence truth (`label == X and homology_only == "no"`) on `truth_set.tsv.gz` built by the Task 6 code, P-ext / N-int / N-sec: S288C 88 / 2,360 / 1,457; *C. albicans* 211 / 179 / 263; *S. pombe* (pombase) 42 / 2,799 / 886; H99 9 / 17 / 15; JEC21 0 / 3 / 1; Af293 23 / 18 / 26; *A. nidulans* 113 / 97 / 66; *U. maydis* 10 / 6 / 21. Ambiguous genes whose label turns P-ext when labels are recomputed without homology codes: S288C 8, *C. albicans* 83, *S. pombe* 2, *A. nidulans* 35.
- Ambiguous genes with only high-throughput internal evidence (`ambiguous_htp_only`): S288C 18 (for example TIP1, CWP2, CCW14, HSP150, ECM33, CTS1), *S. pombe* 4, all other sources 0.
- Extracellular-region genes with an experimental code: S288C 125, *C. albicans* 305. With an internal term at any evidence level: 39 and 102. With an experimental internal code: 32 and 13. Command (run in `analysis/step1_compare/` after Task 6):

  ```python
  import gaf, go_obo, labels
  o = go_obo.parse_obo("/tmp/glyco_spec/go-basic.obo")
  for fn, tx in [("sgd.gaf.gz", None), ("cgd.gaf.gz", "237561")]:
      f = gaf.filter_gaf("/tmp/glyco_spec/" + fn, o, tx)
      ext, any_int, exp_int = set(), set(), set()
      for r in f.cc_rows:
          a = o.ancestors(r.term)
          if r.evidence in gaf.EXPERIMENTAL_CODES and labels.EXTRACELLULAR in a:
              ext.add(r.gene_id)
          if a & labels.INTERNAL:
              any_int.add(r.gene_id)
              if r.evidence in gaf.EXPERIMENTAL_CODES:
                  exp_int.add(r.gene_id)
      print(fn, len(ext), len(ext & any_int), len(ext & exp_int))
  ```
- HKR1, YPS1 (S288C) and chiA (*A. nidulans*) have surface evidence from IBA only, so `homology_only == "yes"` and they are outside the direct-evidence stratum.

## File structure

| File | Responsibility | Task |
|---|---|---|
| `analysis/step1_compare/paths.py` | repository root from `PROJ_ROOT`, work directory from `STEP1_WORKDIR` or `config/site.yaml` | 1 |
| `tests/step1_compare/conftest.py` | puts the analysis folder on `sys.path`; `load_script`; fixtures | 1 |
| `.github/workflows/build_and_test.yml` | adds `step1_compare` to the analysis-tests matrix | 1 |
| `analysis/step1_compare/go_obo.py` | parse `go-basic.obo`: `is_a` + `part_of` ancestors, obsolete terms, alt_ids | 2 |
| `tests/step1_compare/fixtures/mini.obo` | 21-term ontology fixture | 2 |
| `analysis/step1_compare/gaf.py` | GAF reader and the spec 3.1 row filters; evidence-code sets | 3 |
| `tests/step1_compare/gaf_fixture.py` | the golden 50-row GAF and expected labels | 3 |
| `analysis/step1_compare/labels.py` | label rules (P-ext, ambiguous, N-int, N-sec, unlabelled), subset, PM candidate | 4 |
| `analysis/step1_compare/truth_table.py` | truth-table rows (Task 4); counts and TSV I/O (Task 6) | 4, 6 |
| `analysis/step1_compare/manifest.py` | manifest I/O, SHA-256, payload checks, HTTP and UniProt paging | 5 |
| `analysis/step1_compare/manifest.tsv` | pinned inputs: file, URL, SHA-256, size, dates, mode | 5 |
| `analysis/step1_compare/species.tsv` | one row per GO source: species, taxon, clade, role, files, ID mapping | 5 |
| `analysis/step1_compare/00_fetch_inputs.py` | download and verify the pinned inputs | 5 |
| `analysis/step1_compare/01_extract_go_truth.py` | D1: `truth_set.tsv.gz`, `counts.tsv`, `extract_log.json` | 6 |
| `analysis/step1_compare/seqhash.py` | sequence cleaning and SHA-256 | 7 |
| `analysis/step1_compare/sequences.py` | FASTA reader and the four ID mappings | 7 |
| `analysis/step1_compare/02_attach_sequences.py` | D1 step 2: `truth_sequences.tsv.gz`, `unmatched_ids.tsv`, `sequence_counts.tsv` | 7 |
| `analysis/step1_compare/d8_triage.py` | UniProt JSON parser and the P-gpi / PM-TM rule | 8 |
| `analysis/step1_compare/curated_gpi.tsv` | literature GPI evidence (header only at start) | 8 |
| `analysis/step1_compare/03_triage_pm.py` | D8: triage tables, counts, `truth_set_triaged.tsv.gz` | 8 |
| `analysis/step1_compare/keyword_tier.py` | D10 removal rules | 9 |
| `analysis/step1_compare/04_build_keyword_tier.py` | D10: `keyword_tier.tsv.gz`, `keyword_tier_removed.tsv` | 9 |
| `analysis/step1_compare/README.md` | how to re-run; truth subsets and the open question as filters | 10 |

### Output tables (all in `$STEP1_WORKDIR`)

`truth_set.tsv.gz`, one row per gene with at least one aspect-C row that passes the spec 3.1 filters:

| Column | Content |
|---|---|
| `source_id` | `species.tsv` key, for example `Scer_SGD` |
| `species`, `taxon_id` | from `species.tsv` |
| `in_clade` | Saccharomycotina, Taphrinomycotina, Eurotiomycetes or Basidiomycota |
| `role` | `train`, `test_species`, `test_clade`, `alternate_file` or `undecided` |
| `gene_id`, `symbol`, `synonym1` | GAF columns 2, 3 and the first value of column 11 |
| `label` | rule of spec 2.2 on non-IEA evidence: `P-ext`, `ambiguous`, `N-int`, `N-sec`, `unlabelled` |
| `subset` | `wall` or `extracellular-only` for P-ext, else empty |
| `stratum` | `subset` for P-ext, else `label`; D8 replaces it with `P-gpi`, `PM-TM` or `pm-unresolved` in `truth_set_triaged.tsv.gz` |
| `tier` | `T-a` |
| `label_no_homology` | same rule, evidence without IEA and without IBA, IBD, IKR, IRD, ISS, ISO, ISA, ISM, RCA; stored for traceability, not used as truth |
| `label_experimental` | same rule, evidence EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP only |
| `homology_only` | `yes` if `label` differs from `label_no_homology`; the direct-evidence stratum is `label == X and homology_only == "no"` |
| `pm_candidate` | `yes` if P-ext and a non-IEA plasma-membrane (GO:0005886) ancestor |
| `evidence_codes` | sorted codes of all aspect-C rows, IEA included |
| `surface_evidence`, `internal_evidence`, `secretory_evidence` | sorted non-IEA codes on rows whose term reaches wall/extracellular, cytosol/nucleus/mitochondrion, or endomembrane/PM/vacuole |
| `internal_evidence_htp_only` | `yes` if the gene has non-IEA internal evidence and every code is HDA, HMP, HEP, HGI or HTP (R-A) |
| `source_file`, `source_sha256`, `source_date` | GAF file name, its SHA-256, its `!date-generated` header |
| `obo_sha256` | SHA-256 of `go-basic.obo` |

The spec D1 list also names a `cluster` column. Clusters need MMseqs2 on sequences (spec section 4, step 6). That step belongs to the later dataset plan, which appends the column.

`counts.tsv` columns: `source_id, primary_db, genes_cc, genes_noniea_cc, p_ext, p_ext_wall, p_ext_extonly, n_int, n_sec, ambiguous, pm_candidates, cc_iea_triples, cc_triples, cc_iea_frac, all_aspects_iea_frac, obsolete_rows, unknown_term_rows, exp_p_ext, exp_n_int, exp_n_sec, exp_ambiguous, nohom_p_ext, nohom_n_int, nohom_n_sec, nohom_ambiguous, direct_p_ext, direct_n_int, direct_n_sec, direct_ambiguous, ambiguous_htp_only`. The `nohom_*` columns reproduce the `d1_count.py` output; the `direct_*` columns are the spec 3.3 direct-evidence counts.

---

### Task 1: Package skeleton, path helper and CI entry

**Files:**
- Create: `analysis/step1_compare/paths.py`
- Create: `tests/step1_compare/conftest.py`
- Create: `tests/step1_compare/test_paths.py`
- Modify: `.github/workflows/build_and_test.yml:57`

**Interfaces:**
- Consumes: nothing.
- Produces: `paths.STEP1_DIR: Path`; `paths.repo_root() -> Path`; `paths.site_value(key: str, site_yaml: Path | None = None) -> str`; `paths.workdir() -> Path`; `paths.downloads_dir() -> Path`. In tests: `conftest.load_script(name: str) -> module`, fixtures `fixtures_dir -> Path` and `mini_ontology -> go_obo.Ontology` (used from Task 2 on).

- [ ] **Step 1: Write the conftest and the failing test**

`tests/step1_compare/conftest.py`:

```python
"""Make analysis/step1_compare importable as plain modules for tests."""

import importlib.util
import sys
from pathlib import Path

import pytest

STEP1_DIR = Path(__file__).resolve().parents[2] / "analysis" / "step1_compare"
TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(STEP1_DIR))
sys.path.insert(0, str(TESTS_DIR))


def load_script(name: str):
    """Import a numbered script such as 01_extract_go_truth.py as a module."""
    spec = importlib.util.spec_from_file_location(f"step1_{name}", STEP1_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def fixtures_dir() -> Path:
    return TESTS_DIR / "fixtures"


@pytest.fixture
def mini_ontology(fixtures_dir):
    import go_obo

    return go_obo.parse_obo(fixtures_dir / "mini.obo")
```

`tests/step1_compare/test_paths.py`:

```python
import ast
import sys
from pathlib import Path

import paths


def test_workdir_from_env(monkeypatch, tmp_path):
    monkeypatch.setenv("STEP1_WORKDIR", str(tmp_path))
    assert paths.workdir() == tmp_path
    assert paths.downloads_dir() == tmp_path / "downloads"


def test_workdir_from_site_yaml(monkeypatch, tmp_path):
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "site.yaml").write_text("fungi5k_input: /a\nworkdir: /x/y\n")
    monkeypatch.delenv("STEP1_WORKDIR", raising=False)
    monkeypatch.setenv("PROJ_ROOT", str(tmp_path))
    assert paths.repo_root() == tmp_path.resolve()
    assert paths.workdir() == Path("/x/y/step1_compare")


def test_every_module_imports_only_stdlib_or_local():
    modules = sorted(paths.STEP1_DIR.glob("*.py"))
    local = {p.stem for p in modules}
    for path in modules:
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [(node.module or "").split(".")[0]]
            else:
                continue
            for name in names:
                assert name in sys.stdlib_module_names or name in local, f"{path.name}: {name}"


def test_no_shell_script_uses_bash_source():
    for script in paths.STEP1_DIR.rglob("*.sh"):
        assert "BASH_SOURCE" not in script.read_text(), script
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_paths.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'paths'`.

- [ ] **Step 3: Write `analysis/step1_compare/paths.py`**

```python
"""Locate the repository root and the step 1 work directory. Standard library only.

The repository root comes from the PROJ_ROOT environment variable. If PROJ_ROOT is not set,
the root is two levels above this file. The work directory comes from STEP1_WORKDIR. If
STEP1_WORKDIR is not set, it is `<workdir from config/site.yaml>/step1_compare`.
"""

import os
from pathlib import Path

STEP1_DIR = Path(__file__).resolve().parent


def repo_root() -> Path:
    env = os.environ.get("PROJ_ROOT")
    if env:
        return Path(env).resolve()
    return STEP1_DIR.parents[1]


def site_value(key: str, site_yaml: Path | None = None) -> str:
    """Read one top-level `key: value` line from config/site.yaml without PyYAML."""
    path = site_yaml or repo_root() / "config" / "site.yaml"
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(f"{key}:"):
                return line.split(":", 1)[1].strip()
    raise KeyError(f"{key} not found in {path}")


def workdir() -> Path:
    env = os.environ.get("STEP1_WORKDIR")
    if env:
        return Path(env)
    return Path(site_value("workdir")) / "step1_compare"


def downloads_dir() -> Path:
    return workdir() / "downloads"
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `$PY -m pytest tests/step1_compare/test_paths.py -q`
Expected: PASS (4 passed).

- [ ] **Step 5: Add the CI matrix entry**

In `.github/workflows/build_and_test.yml`, change line 57 from

```yaml
        suite: [adhesion_properties, embedding_clustering, kingdom_survey]
```

to

```yaml
        suite: [adhesion_properties, embedding_clustering, kingdom_survey, step1_compare]
```

Each suite runs in its own job and its own pytest process, so the plain module names in `analysis/step1_compare/` do not collide with other analysis folders. The analysis-tests job has `continue-on-error: true`, so this suite is informational in CI. Do not change that setting in this plan; a required job is an owner decision.

- [ ] **Step 6: Lint and commit**

```bash
pre-commit run --all-files   # re-run and re-stage if it rewrites files
ruff check .
git add analysis/step1_compare/paths.py tests/step1_compare/conftest.py \
        tests/step1_compare/test_paths.py .github/workflows/build_and_test.yml
git commit -m "step1_compare: package skeleton, work-directory helper, CI suite

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 2: GO ontology parser

**Files:**
- Create: `analysis/step1_compare/go_obo.py`
- Create: `tests/step1_compare/fixtures/mini.obo`
- Create: `tests/step1_compare/test_go_obo.py`

**Interfaces:**
- Consumes: `conftest.mini_ontology` fixture (Task 1).
- Produces: `go_obo.parse_obo(path: str | Path) -> Ontology`; `Ontology.parents: dict[str, set[str]]`, `Ontology.names: dict[str, str]`, `Ontology.obsolete: set[str]`, `Ontology.alt_ids: dict[str, str]`, `Ontology.data_version: str`, `Ontology.known(term: str) -> bool`, `Ontology.ancestors(term: str) -> frozenset[str]` (includes the term; an unknown term gives `{term}`, as in `d1_count.py`).

- [ ] **Step 1: Write the fixture and the failing test**

`tests/step1_compare/fixtures/mini.obo`:

```text
format-version: 1.2
data-version: releases/2026-07-26

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Term]
id: GO:0110165
name: cellular anatomical structure
namespace: cellular_component
is_a: GO:0005575 ! cellular_component

[Term]
id: GO:0071944
name: cell periphery
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure

[Term]
id: GO:0005618
name: cell wall
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure
relationship: part_of GO:0071944 ! cell periphery

[Term]
id: GO:0009277
name: fungal-type cell wall
namespace: cellular_component
alt_id: GO:0099999
is_a: GO:0005618 ! cell wall

[Term]
id: GO:0005576
name: extracellular region
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure

[Term]
id: GO:0016020
name: membrane
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure

[Term]
id: GO:0005886
name: plasma membrane
namespace: cellular_component
is_a: GO:0016020 ! membrane
relationship: part_of GO:0071944 ! cell periphery

[Term]
id: GO:0005737
name: cytoplasm
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure

[Term]
id: GO:0005829
name: cytosol
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure
relationship: part_of GO:0005737 ! cytoplasm

[Term]
id: GO:0005634
name: nucleus
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure

[Term]
id: GO:0005739
name: mitochondrion
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure

[Term]
id: GO:0005743
name: mitochondrial inner membrane
namespace: cellular_component
is_a: GO:0016020 ! membrane
relationship: part_of GO:0005739 ! mitochondrion

[Term]
id: GO:0012505
name: endomembrane system
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure

[Term]
id: GO:0005783
name: endoplasmic reticulum
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure
relationship: part_of GO:0012505 ! endomembrane system

[Term]
id: GO:0005773
name: vacuole
namespace: cellular_component
is_a: GO:0110165 ! cellular anatomical structure

[Term]
id: GO:0000324
name: fungal-type vacuole
namespace: cellular_component
is_a: GO:0005773 ! vacuole

[Term]
id: GO:0032991
name: protein-containing complex
namespace: cellular_component
is_a: GO:0005575 ! cellular_component
relationship: has_part GO:0005886 ! plasma membrane

[Term]
id: GO:0031225
name: obsolete anchored component of membrane
namespace: cellular_component
is_obsolete: true

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Typedef]
id: part_of
name: part of
is_a: GO:0005575
is_transitive: true
```

`tests/step1_compare/test_go_obo.py`:

```python
def test_ancestors_follow_is_a_and_part_of(mini_ontology):
    assert mini_ontology.ancestors("GO:0009277") == {
        "GO:0009277",
        "GO:0005618",
        "GO:0110165",
        "GO:0005575",
        "GO:0071944",
    }
    assert "GO:0012505" in mini_ontology.ancestors("GO:0005783")  # ER part_of endomembrane
    assert {"GO:0016020", "GO:0005739"} <= mini_ontology.ancestors("GO:0005743")


def test_other_relationships_and_typedefs_are_ignored(mini_ontology):
    assert "GO:0005886" not in mini_ontology.ancestors("GO:0032991")  # has_part ignored
    assert mini_ontology.parents["GO:0008150"] == set()  # Typedef is_a not attached


def test_obsolete_alt_id_and_version(mini_ontology):
    assert mini_ontology.obsolete == {"GO:0031225"}
    assert mini_ontology.alt_ids == {"GO:0099999": "GO:0009277"}
    assert mini_ontology.ancestors("GO:0099999") == {"GO:0099999"}  # as d1_count.py
    assert mini_ontology.data_version == "releases/2026-07-26"
    assert mini_ontology.known("GO:0005576")
    assert not mini_ontology.known("GO:1234567")
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_go_obo.py -q`
Expected: 3 errors at fixture setup with `ModuleNotFoundError: No module named 'go_obo'`.

- [ ] **Step 3: Write `analysis/step1_compare/go_obo.py`**

The parser copies `d1_count.py`: only `[Term]` stanzas; parents from `is_a:` and `relationship: part_of`; other relationships ignored; `alt_id` recorded but not used for ancestors. The ancestor walk is iterative, so deep ontologies do not reach the recursion limit.

```python
"""Parse go-basic.obo: parents over is_a and part_of, obsolete terms (spec 3.1 step 1).

Standard library only. The rules copy /tmp/glyco_spec/d1_count.py: only [Term] stanzas count,
parents are `is_a` and `relationship: part_of`, other relationships are ignored, and alt_id
lines are recorded but not used for ancestors.
"""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Ontology:
    parents: dict[str, set[str]] = field(default_factory=dict)
    names: dict[str, str] = field(default_factory=dict)
    obsolete: set[str] = field(default_factory=set)
    alt_ids: dict[str, str] = field(default_factory=dict)
    data_version: str = ""
    _cache: dict[str, frozenset[str]] = field(default_factory=dict, repr=False)

    def known(self, term: str) -> bool:
        return term in self.names

    def ancestors(self, term: str) -> frozenset[str]:
        """The term itself plus every is_a / part_of ancestor. Unknown terms give {term}."""
        cached = self._cache.get(term)
        if cached is not None:
            return cached
        seen = {term}
        stack = [term]
        while stack:
            for parent in self.parents.get(stack.pop(), ()):
                if parent not in seen:
                    seen.add(parent)
                    stack.append(parent)
        result = frozenset(seen)
        self._cache[term] = result
        return result


def parse_obo(path: str | Path) -> Ontology:
    onto = Ontology()
    current = None
    in_term = False
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("["):
                in_term = line == "[Term]"
                current = None
            elif not in_term and current is None and line.startswith("data-version:"):
                onto.data_version = line.split(":", 1)[1].strip()
            elif in_term and line.startswith("id: GO:"):
                current = line.split()[1]
                onto.names.setdefault(current, "")
                onto.parents.setdefault(current, set())
            elif current is None:
                continue
            elif line.startswith("name: "):
                onto.names[current] = line[6:]
            elif line.startswith("is_a: "):
                onto.parents[current].add(line[6:].split()[0])
            elif line.startswith("relationship: part_of "):
                onto.parents[current].add(line.split()[2])
            elif line == "is_obsolete: true":
                onto.obsolete.add(current)
            elif line.startswith("alt_id: "):
                onto.alt_ids[line[8:].split()[0]] = current
    return onto
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `$PY -m pytest tests/step1_compare/test_go_obo.py -q`
Expected: PASS (5 passed).

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --all-files
ruff check .
git add analysis/step1_compare/go_obo.py tests/step1_compare/fixtures/mini.obo \
        tests/step1_compare/test_go_obo.py
git commit -m "step1_compare: GO obo parser (is_a + part_of ancestors, obsolete terms)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 3: GAF reader with the spec 3.1 filters

**Files:**
- Create: `analysis/step1_compare/gaf.py`
- Create: `tests/step1_compare/gaf_fixture.py`
- Create: `tests/step1_compare/test_gaf.py`

**Interfaces:**
- Consumes: `go_obo.Ontology.obsolete`, `Ontology.known` (Task 2).
- Produces: `gaf.HOMOLOGY_CODES`, `gaf.EXPERIMENTAL_CODES`, `gaf.PROTEIN_TYPES` (frozensets); `gaf.GafFormatError`; `gaf.GafRow` (fields `db, gene_id, symbol, qualifier, term, evidence, aspect, synonyms, object_type, taxon`); `gaf.open_text(path)` (gzip by magic bytes); `gaf.iter_gaf_rows(path)`; `gaf.header_value(path, key) -> str`; `gaf.taxon_matches(taxon_field, taxon_id) -> bool`; `gaf.FilteredGaf` (fields `path, primary_db, cc_rows: list[GafRow], all_aspect_triples: set[tuple[str, str, str, str]], dropped: Counter, unknown_term_rows: int`); `gaf.filter_gaf(path, ontology, taxon_id: str | None = None) -> FilteredGaf`. In tests: `gaf_fixture.ROWS`, `gaf_fixture.EXPECTED_LABELS`, `gaf_fixture.write_golden_gaf(path) -> Path`.

The fixture has 50 data rows. It contains complex rows (ComplexPortal), ncRNA rows, duplicate-database rows (UniProtKB in a MOD file), a `NOT` row, IEA rows, homology-code rows (IBA, ISS), an obsolete term, an other-strain taxon and an `NCBITaxon:` prefix. Hand-computed results: primary DB `SGD`; dropped `other_db` 5, `object_type` 2, `other_aspect` 6, `not_qualifier` 1, `obsolete_term` 1; 35 aspect-C rows kept; 43 all-aspect triples; 20 genes.

- [ ] **Step 1: Write the fixture and the failing test**

`tests/step1_compare/gaf_fixture.py`:

```python
"""Golden 50-row GAF fixture for the D1 extractor (spec section 7, test_gaf_extract_golden).

Each tuple is (db, gene_id, symbol, qualifier, term, evidence, aspect, object_type, taxon).
The comment after each gene block gives its expected label under the non-IEA policy.
"""

from pathlib import Path

TAXON = "taxon:559292"
P, F = "GO:0008150", "GO:0003674"  # biological_process, molecular_function roots
WALL, EXT, PM = "GO:0009277", "GO:0005576", "GO:0005886"
CYTO, NUC, MITO, MITO_IM = "GO:0005829", "GO:0005634", "GO:0005739", "GO:0005743"
ER, VAC, MEM, PERIPH, OBSOLETE = (
    "GO:0005783",
    "GO:0000324",
    "GO:0016020",
    "GO:0071944",
    "GO:0031225",
)

ROWS = [
    # 1-3 G01 wall IDA -> P-ext (wall); F and P aspect rows are dropped from labels
    ("SGD", "G01", "WALL1", "located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G01", "WALL1", "enables", F, "IDA", "F", "protein", TAXON),
    ("SGD", "G01", "WALL1", "involved_in", P, "IEA", "P", "protein", TAXON),
    # 4-5 G02 extracellular IDA (+ IEA copy) -> P-ext (extracellular-only)
    ("SGD", "G02", "SECR1", "located_in", EXT, "IDA", "C", "protein", TAXON),
    ("SGD", "G02", "SECR1", "located_in", EXT, "IEA", "C", "protein", TAXON),
    # 6-7 G03 wall + plasma membrane IDA -> P-ext (wall), pm_candidate
    ("SGD", "G03", "GPIPM", "located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G03", "GPIPM", "located_in", PM, "IDA", "C", "protein", TAXON),
    # 8-9 G04 extracellular IDA + cytosol IDA -> ambiguous
    ("SGD", "G04", "MOON1", "located_in", EXT, "IDA", "C", "protein", TAXON),
    ("SGD", "G04", "MOON1", "located_in", CYTO, "IDA", "C", "protein", TAXON),
    # 10-11 G05 cytosol IDA, nucleus IEA -> N-int
    ("SGD", "G05", "CYTO1", "located_in", CYTO, "IDA", "C", "protein", TAXON),
    ("SGD", "G05", "CYTO1", "located_in", NUC, "IEA", "C", "protein", TAXON),
    # 12-13 G06 nucleus IDA, membrane IEA -> unlabelled (membrane at any level blocks N-int)
    ("SGD", "G06", "NUCMEM", "located_in", NUC, "IDA", "C", "protein", TAXON),
    ("SGD", "G06", "NUCMEM", "located_in", MEM, "IEA", "C", "protein", TAXON),
    # 14 G07 ER IDA -> N-sec (ER is part_of endomembrane system)
    ("SGD", "G07", "ER1", "located_in", ER, "IDA", "C", "protein", TAXON),
    # 15-16 G08 plasma membrane IDA + nucleus IDA -> N-sec (internal terms allowed in N-sec)
    ("SGD", "G08", "PMNUC", "located_in", PM, "IDA", "C", "protein", TAXON),
    ("SGD", "G08", "PMNUC", "located_in", NUC, "IDA", "C", "protein", TAXON),
    # 17-19 G09 wall IEA + cytosol IEA + PM IDA -> unlabelled (IEA wall blocks N-sec)
    ("SGD", "G09", "IEAWALL", "located_in", WALL, "IEA", "C", "protein", TAXON),
    ("SGD", "G09", "IEAWALL", "located_in", CYTO, "IEA", "C", "protein", TAXON),
    ("SGD", "G09", "IEAWALL", "located_in", PM, "IDA", "C", "protein", TAXON),
    # 20-22 G10 wall IDA + wall IBA + cytosol IBA -> ambiguous; P-ext without homology codes
    ("SGD", "G10", "IBAWALL", "located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G10", "IBAWALL", "is_active_in", WALL, "IBA", "C", "protein", TAXON),
    ("SGD", "G10", "IBAWALL", "located_in", CYTO, "IBA", "C", "protein", TAXON),
    # 23 G11 extracellular IBA only -> P-ext (extracellular-only); unlabelled without homology
    ("SGD", "G11", "IBAONLY", "located_in", EXT, "IBA", "C", "protein", TAXON),
    # 24-25 G12 NOT wall IDA + vacuole IDA -> N-sec (NOT row skipped)
    ("SGD", "G12", "NOTWALL", "NOT|located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G12", "NOTWALL", "located_in", VAC, "IDA", "C", "protein", TAXON),
    # 26-27 G13 obsolete term IDA + ER IDA -> N-sec; obsolete row counted
    ("SGD", "G13", "OBSOL", "located_in", OBSOLETE, "IDA", "C", "protein", TAXON),
    ("SGD", "G13", "OBSOL", "located_in", ER, "IDA", "C", "protein", TAXON),
    # 28-29 G14 wall IDA + mitochondrion HDA -> ambiguous (as GAS1 in sgd.gaf.gz)
    ("SGD", "G14", "HDAMITO", "located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G14", "HDAMITO", "located_in", MITO, "HDA", "C", "protein", TAXON),
    # 30 G15 vacuole IEA only -> unlabelled; counted in genes_cc, not in genes_noniea_cc
    ("SGD", "G15", "VAC1", "located_in", VAC, "IEA", "C", "protein", TAXON),
    # 31-32 G16 wall ISS + ER IDA -> P-ext (wall); unlabelled without homology codes
    ("SGD", "G16", "ISSWALL", "located_in", WALL, "ISS", "C", "protein", TAXON),
    ("SGD", "G16", "ISSWALL", "located_in", ER, "IDA", "C", "protein", TAXON),
    # 33 G17 mitochondrial inner membrane IDA -> unlabelled (membrane ancestor blocks N-int)
    ("SGD", "G17", "MITOMEM", "located_in", MITO_IM, "IDA", "C", "protein", TAXON),
    # 34 G18 cell periphery IDA only -> unlabelled
    ("SGD", "G18", "PERIPH", "located_in", PERIPH, "IDA", "C", "protein", TAXON),
    # 35-36 complex rows from ComplexPortal -> dropped (other_db)
    ("ComplexPortal", "CPX-1", "cpx1", "part_of", CYTO, "IDA", "C", "protein_complex", TAXON),
    ("ComplexPortal", "CPX-1", "cpx1", "part_of", NUC, "IDA", "C", "protein_complex", TAXON),
    # 37-38 ncRNA rows from the primary database -> dropped (object_type)
    ("SGD", "R01", "SNR1", "located_in", NUC, "IDA", "C", "ncRNA", TAXON),
    ("SGD", "R01", "SNR1", "located_in", CYTO, "IDA", "C", "ncRNA", TAXON),
    # 39-41 duplicate-database rows (UniProtKB in a MOD file) -> dropped (other_db)
    ("UniProtKB", "P99999", "DUP1", "located_in", WALL, "IBA", "C", "protein", TAXON),
    ("UniProtKB", "P99999", "DUP1", "involved_in", P, "IBA", "P", "protein", TAXON),
    ("UniProtKB", "P99999", "DUP1", "enables", F, "IBA", "F", "protein", TAXON),
    # 42-43 G21 other strain -> P-ext without a taxon filter, dropped with taxon 559292
    ("SGD", "G21", "STRAIN2", "located_in", WALL, "IDA", "C", "protein", "taxon:999999"),
    ("SGD", "G21", "STRAIN2", "involved_in", P, "IDA", "P", "protein", "taxon:999999"),
    # 44 G22 NCBITaxon prefix for the same taxon -> kept, P-ext (extracellular-only)
    ("SGD", "G22", "NCBI1", "located_in", EXT, "IDA", "C", "gene_product", "NCBITaxon:559292"),
    # 45-50 more rows for genes above
    ("SGD", "G05", "CYTO1", "involved_in", P, "IDA", "P", "protein", TAXON),
    ("SGD", "G07", "ER1", "enables", F, "IEA", "F", "protein", TAXON),
    ("SGD", "G03", "GPIPM", "located_in", PM, "IEA", "C", "protein", TAXON),
    ("SGD", "G08", "PMNUC", "located_in", NUC, "HDA", "C", "protein", TAXON),
    ("SGD", "G11", "IBAONLY", "enables", F, "IBA", "F", "protein", TAXON),
    ("SGD", "G04", "MOON1", "located_in", EXT, "EXP", "C", "protein", TAXON),
]

EXPECTED_LABELS = {
    "G01": ("P-ext", "wall"),
    "G02": ("P-ext", "extracellular-only"),
    "G03": ("P-ext", "wall"),
    "G04": ("ambiguous", ""),
    "G05": ("N-int", ""),
    "G06": ("unlabelled", ""),
    "G07": ("N-sec", ""),
    "G08": ("N-sec", ""),
    "G09": ("unlabelled", ""),
    "G10": ("ambiguous", ""),
    "G11": ("P-ext", "extracellular-only"),
    "G12": ("N-sec", ""),
    "G13": ("N-sec", ""),
    "G14": ("ambiguous", ""),
    "G15": ("unlabelled", ""),
    "G16": ("P-ext", "wall"),
    "G17": ("unlabelled", ""),
    "G18": ("unlabelled", ""),
    "G21": ("P-ext", "wall"),
    "G22": ("P-ext", "extracellular-only"),
}


def write_golden_gaf(path: Path) -> Path:
    lines = ["!gaf-version: 2.2", "!date-generated: 2026-05-21T09:00", ""]
    for db, gene, symbol, qual, term, ev, aspect, otype, taxon in ROWS:
        cols = [
            db,
            gene,
            symbol,
            qual,
            term,
            "PMID:1",
            ev,
            "",
            aspect,
            symbol.lower(),
            f"{symbol}|{gene}-syn",
            otype,
            taxon,
            "20260101",
            "SGD",
            "",
            "",
        ]
        lines.append("\t".join(cols))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
```

`tests/step1_compare/test_gaf.py`:

```python
import gzip

import gaf
import pytest
from gaf_fixture import write_golden_gaf


def test_gaf_extract_golden(tmp_path, mini_ontology):
    path = write_golden_gaf(tmp_path / "golden.gaf")
    out = gaf.filter_gaf(path, mini_ontology)
    assert out.primary_db == "SGD"
    assert dict(out.dropped) == {
        "other_db": 5,
        "object_type": 2,
        "other_aspect": 6,
        "not_qualifier": 1,
        "obsolete_term": 1,
    }
    assert len(out.cc_rows) == 35
    assert len(out.all_aspect_triples) == 43
    genes = {r.gene_id for r in out.cc_rows}
    assert {"CPX-1", "R01", "P99999"}.isdisjoint(genes)
    assert len(genes) == 20
    assert out.unknown_term_rows == 0


def test_taxon_filter_accepts_ncbitaxon_prefix(tmp_path, mini_ontology):
    path = write_golden_gaf(tmp_path / "golden.gaf")
    out = gaf.filter_gaf(path, mini_ontology, taxon_id="559292")
    genes = {r.gene_id for r in out.cc_rows}
    assert "G21" not in genes and "G22" in genes
    assert out.dropped["taxon"] == 2


def test_gzip_detected_by_magic_not_suffix(tmp_path, mini_ontology):
    plain = write_golden_gaf(tmp_path / "golden.gaf")
    packed = tmp_path / "golden_no_suffix"
    packed.write_bytes(gzip.compress(plain.read_bytes()))
    assert len(gaf.filter_gaf(packed, mini_ontology).cc_rows) == 35
    assert gaf.header_value(packed, "date-generated") == "2026-05-21T09:00"


def test_short_row_and_empty_file_raise(tmp_path, mini_ontology):
    bad = tmp_path / "bad.gaf"
    bad.write_text("!gaf-version: 2.2\nSGD\tS1\tX\n")
    with pytest.raises(gaf.GafFormatError, match="3 columns"):
        gaf.filter_gaf(bad, mini_ontology)
    empty = tmp_path / "empty.gaf"
    empty.write_text("!gaf-version: 2.2\n")
    with pytest.raises(gaf.GafFormatError, match="no data rows"):
        gaf.filter_gaf(empty, mini_ontology)


def test_evidence_code_sets_match_spec():
    assert gaf.HOMOLOGY_CODES == {"IBA", "IBD", "IKR", "IRD", "ISS", "ISO", "ISA", "ISM", "RCA"}
    assert "IEA" not in gaf.EXPERIMENTAL_CODES | gaf.HOMOLOGY_CODES
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_gaf.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'gaf'`.

- [ ] **Step 3: Write `analysis/step1_compare/gaf.py`**

The filter order copies `d1_count.py` line by line: database and object type, then taxon, then the all-aspect triple, then aspect `C` and `NOT`, then obsolete terms. The taxon check compares the NCBI number, so `taxon:237561` and `NCBITaxon:237561` both match. (`d1_count.py` uses a substring test; on the real files both give the same rows, which the Task 6 regression test confirms.)

```python
"""Read GAF 2.x files and apply the D1 row filters (spec section 3.1, steps 2 to 4).

Standard library only. The filters copy /tmp/glyco_spec/d1_count.py in this order:
1. column 1 must be the most frequent database in the file, and column 12 must be
   `protein` or `gene_product`;
2. if a taxon is given, the first taxon in column 13 must equal it;
3. (gene, term, evidence, aspect) is recorded for the all-aspect IEA fraction;
4. aspect must be `C` and the qualifier must not contain `NOT`;
5. rows with an obsolete term are counted and dropped.
"""

import collections
import gzip
from dataclasses import dataclass, field
from pathlib import Path

HOMOLOGY_CODES = frozenset({"IBA", "IBD", "IKR", "IRD", "ISS", "ISO", "ISA", "ISM", "RCA"})
EXPERIMENTAL_CODES = frozenset(
    {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"}
)
PROTEIN_TYPES = frozenset({"protein", "gene_product"})
GZIP_MAGIC = b"\x1f\x8b"


class GafFormatError(ValueError):
    """A GAF file is empty, truncated or not a GAF file."""


@dataclass(frozen=True, slots=True)
class GafRow:
    db: str
    gene_id: str
    symbol: str
    qualifier: str
    term: str
    evidence: str
    aspect: str
    synonyms: str
    object_type: str
    taxon: str


def open_text(path: str | Path):
    """Open plain or gzip text. Gzip is detected by its magic bytes, not by the suffix."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == GZIP_MAGIC:
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def iter_gaf_rows(path: str | Path):
    with open_text(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            if raw.startswith("!") or not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 13:
                raise GafFormatError(f"{path}:{lineno}: {len(f)} columns, GAF needs at least 13")
            yield GafRow(f[0], f[1], f[2], f[3], f[4], f[6], f[8], f[10], f[11], f[12])


def header_value(path: str | Path, key: str) -> str:
    """Return the value of the first `!key:` header line, or '' if absent."""
    prefix = f"!{key}:"
    with open_text(path) as handle:
        for raw in handle:
            if not raw.startswith("!"):
                break
            if raw.startswith(prefix):
                return raw[len(prefix) :].strip()
    return ""


def taxon_matches(taxon_field: str, taxon_id: str) -> bool:
    """True if the first taxon in GAF column 13 has this NCBI id (`taxon:` or `NCBITaxon:`)."""
    first = taxon_field.split("|")[0]
    return first.split(":")[-1] == str(taxon_id)


@dataclass
class FilteredGaf:
    path: str
    primary_db: str
    cc_rows: list[GafRow]
    all_aspect_triples: set[tuple[str, str, str, str]]
    dropped: collections.Counter = field(default_factory=collections.Counter)
    unknown_term_rows: int = 0


def filter_gaf(path: str | Path, ontology, taxon_id: str | None = None) -> FilteredGaf:
    rows = list(iter_gaf_rows(path))
    if not rows:
        raise GafFormatError(f"{path}: no data rows")
    primary = collections.Counter(r.db for r in rows).most_common(1)[0][0]
    out = FilteredGaf(str(path), primary, [], set())
    for r in rows:
        if r.db != primary:
            out.dropped["other_db"] += 1
            continue
        if r.object_type not in PROTEIN_TYPES:
            out.dropped["object_type"] += 1
            continue
        if taxon_id and not taxon_matches(r.taxon, taxon_id):
            out.dropped["taxon"] += 1
            continue
        out.all_aspect_triples.add((r.gene_id, r.term, r.evidence, r.aspect))
        if r.aspect != "C":
            out.dropped["other_aspect"] += 1
            continue
        if "NOT" in r.qualifier:
            out.dropped["not_qualifier"] += 1
            continue
        if r.term in ontology.obsolete:
            out.dropped["obsolete_term"] += 1
            continue
        if not ontology.known(r.term):
            out.unknown_term_rows += 1
        out.cc_rows.append(r)
    return out
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `$PY -m pytest tests/step1_compare/test_gaf.py -q`
Expected: PASS (5 passed).

- [ ] **Step 5: Lint and commit**

```bash
pre-commit run --all-files
ruff check .
git add analysis/step1_compare/gaf.py tests/step1_compare/gaf_fixture.py \
        tests/step1_compare/test_gaf.py
git commit -m "step1_compare: GAF reader with the D1 row filters and golden fixture

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Label rules and truth-table rows

**Files:**
- Create: `analysis/step1_compare/labels.py`
- Create: `analysis/step1_compare/truth_table.py` (rows only; Task 6 adds counts and I/O)
- Create: `tests/step1_compare/test_labels.py`

**Interfaces:**
- Consumes: `gaf.HOMOLOGY_CODES`, `gaf.EXPERIMENTAL_CODES`, `gaf.FilteredGaf`, `gaf.filter_gaf` (Task 3); `go_obo.Ontology.ancestors` (Task 2); `gaf_fixture` (Task 3).
- Produces: `labels.WALL`, `labels.EXTRACELLULAR`, `labels.PLASMA_MEMBRANE`, `labels.INTERNAL`, `labels.SECRETORY`, `labels.ANY_SECRETORY`, `labels.SURFACE`; label strings `labels.P_EXT = "P-ext"`, `AMBIGUOUS = "ambiguous"`, `N_INT = "N-int"`, `N_SEC = "N-sec"`, `UNLABELLED = "unlabelled"`, `labels.LABELS`; `labels.POLICIES = ("non_iea", "no_homology", "experimental")`; `labels.policy_accepts(policy, evidence) -> bool`; `labels.classify(label_terms, any_terms) -> str`; `labels.subset_of(label, label_terms) -> str`; `labels.is_pm_candidate(label, label_terms) -> bool`; `labels.HIGH_THROUGHPUT_CODES`; `labels.htp_only(internal_codes: set[str]) -> bool`. `truth_table.TRUTH_COLUMNS: tuple[str, ...]`; `truth_table.SourceInfo(source_id, species, taxon_id, in_clade, role, source_file, source_sha256, source_date, obo_sha256)`; `truth_table.GeneRecord`; `truth_table.collect_genes(filtered) -> dict[str, GeneRecord]`; `truth_table.build_truth_rows(filtered, ontology, info) -> list[dict[str, str]]`; `truth_table.in_direct_stratum(row, label) -> bool` (`row["label"] == label and row["homology_only"] == "no"`).

The rules copy `d1_count.py` `sets()`. The labels are mutually exclusive: a surface gene is P-ext or ambiguous; N-int needs no secretory-side term at any level; N-sec needs no surface term at any level and may carry an internal term. Every label column is computed from the same stored rows. Direct-evidence truth is the intersection `in_direct_stratum(row, X)`: the gene keeps the same label with and without homology codes. `label_no_homology` alone is not used as truth, because it moves genes whose only internal term is IBA or ISS (fixture gene G10; 83 genes in *C. albicans*) into P-ext. `internal_evidence_htp_only` marks ambiguous genes whose internal evidence is high-throughput only (R-A); the D1 rule itself does not change.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_labels.py`:

```python
import gaf
import labels
import pytest
import truth_table
from gaf_fixture import EXPECTED_LABELS, write_golden_gaf

WALL, EXT, PM, CYTO, NUC = "GO:0005618", "GO:0005576", "GO:0005886", "GO:0005829", "GO:0005634"
ER, MEM = "GO:0012505", "GO:0016020"


@pytest.mark.parametrize(
    ("label_terms", "any_terms", "expected"),
    [
        ({WALL}, {WALL}, "P-ext"),
        ({EXT}, {EXT, CYTO}, "P-ext"),  # internal term at IEA level only
        ({WALL, CYTO}, {WALL, CYTO}, "ambiguous"),
        ({NUC}, {NUC}, "N-int"),
        ({NUC}, {NUC, MEM}, "unlabelled"),  # membrane at any level blocks N-int
        ({NUC}, {NUC, WALL}, "unlabelled"),
        ({PM}, {PM}, "N-sec"),
        ({PM, NUC}, {PM, NUC}, "N-sec"),
        ({ER}, {ER, EXT}, "unlabelled"),  # surface at any level blocks N-sec
        (set(), {WALL, CYTO}, "unlabelled"),  # wall and cytosol at IEA level only
        (set(), set(), "unlabelled"),
    ],
)
def test_label_rules_truth_table(label_terms, any_terms, expected):
    assert labels.classify(frozenset(label_terms), frozenset(any_terms)) == expected


def test_subset_and_pm_candidate():
    assert labels.subset_of("P-ext", {WALL, EXT}) == "wall"
    assert labels.subset_of("P-ext", {EXT}) == "extracellular-only"
    assert labels.subset_of("N-sec", {PM}) == ""
    assert labels.is_pm_candidate("P-ext", {WALL, PM})
    assert not labels.is_pm_candidate("N-sec", {PM})


def test_policies():
    assert labels.policy_accepts("non_iea", "IBA") and not labels.policy_accepts("non_iea", "IEA")
    assert not labels.policy_accepts("no_homology", "ISS") and labels.policy_accepts(
        "no_homology", "TAS"
    )
    assert labels.policy_accepts("experimental", "HDA") and not labels.policy_accepts(
        "experimental", "TAS"
    )
    with pytest.raises(ValueError):
        labels.policy_accepts("all", "IDA")


def _rows(tmp_path, ontology):
    filtered = gaf.filter_gaf(write_golden_gaf(tmp_path / "g.gaf"), ontology)
    info = truth_table.SourceInfo(
        "Fix",
        "Fixture yeast",
        "559292",
        "Saccharomycotina",
        "train",
        "g.gaf",
        "sha",
        "2026-05-21",
        "obo",
    )
    return {r["gene_id"]: r for r in truth_table.build_truth_rows(filtered, ontology, info)}


def test_golden_labels(tmp_path, mini_ontology):
    rows = _rows(tmp_path, mini_ontology)
    assert {g: (r["label"], r["subset"]) for g, r in rows.items()} == EXPECTED_LABELS
    assert [g for g, r in rows.items() if r["pm_candidate"] == "yes"] == ["G03"]


def test_conflicting_evidence_is_stored_not_collapsed(tmp_path, mini_ontology):
    g10 = _rows(tmp_path, mini_ontology)["G10"]
    assert (g10["label"], g10["label_no_homology"], g10["homology_only"]) == (
        "ambiguous",
        "P-ext",
        "yes",
    )
    assert g10["surface_evidence"] == "IBA,IDA"
    assert g10["internal_evidence"] == "IBA"
    g11 = _rows(tmp_path, mini_ontology)["G11"]
    assert (g11["label"], g11["label_no_homology"], g11["label_experimental"]) == (
        "P-ext",
        "unlabelled",
        "unlabelled",
    )


def test_truth_set_has_no_iea(tmp_path, mini_ontology):
    rows = _rows(tmp_path, mini_ontology)
    support = {
        "P-ext": "surface_evidence",
        "N-int": "internal_evidence",
        "N-sec": "secretory_evidence",
    }
    for r in rows.values():
        for column in ("surface_evidence", "internal_evidence", "secretory_evidence"):
            assert "IEA" not in r[column].split(",")
        if r["label"] in support:
            assert r[support[r["label"]]], f"{r['gene_id']} label has no non-IEA support"
    assert rows["G09"]["label"] == "unlabelled"  # wall + cytosol only at IEA level
    assert rows["G15"]["label"] == "unlabelled"  # IEA-only gene


def test_htp_only_rule():
    assert labels.htp_only({"HDA"}) and labels.htp_only({"HDA", "HTP"})
    assert not labels.htp_only({"HDA", "IDA"}) and not labels.htp_only(set())


def test_direct_stratum_is_an_intersection(tmp_path, mini_ontology):
    rows = _rows(tmp_path, mini_ontology)
    direct = sorted(g for g, r in rows.items() if truth_table.in_direct_stratum(r, "P-ext"))
    assert direct == ["G01", "G02", "G03", "G21", "G22"]
    # G10 (wall IDA + cytosol IBA) turns P-ext when labels are recomputed without homology
    # codes; the intersection rule keeps it out of the direct P-ext stratum.
    assert rows["G10"]["label_no_homology"] == "P-ext"
    assert not truth_table.in_direct_stratum(rows["G10"], "P-ext")
    assert not truth_table.in_direct_stratum(rows["G10"], "ambiguous")


def test_internal_evidence_htp_only_column(tmp_path, mini_ontology):
    rows = _rows(tmp_path, mini_ontology)
    g14, g04, g10 = rows["G14"], rows["G04"], rows["G10"]
    assert (g14["label"], g14["internal_evidence_htp_only"]) == ("ambiguous", "yes")  # mito HDA
    assert (g04["label"], g04["internal_evidence_htp_only"]) == ("ambiguous", "no")  # cytosol IDA
    assert g10["internal_evidence_htp_only"] == "no"  # cytosol IBA only
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_labels.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'labels'`.

- [ ] **Step 3: Write `analysis/step1_compare/labels.py`**

```python
"""Step 1 label rules (spec section 2.2) on GO ancestor sets. Standard library only.

`label_terms` are the ancestors of the rows that an evidence policy accepts. `any_terms` are the
ancestors of all rows of the gene, IEA included. The rules copy /tmp/glyco_spec/d1_count.py.
"""

from gaf import EXPERIMENTAL_CODES, HOMOLOGY_CODES

WALL = "GO:0005618"
EXTRACELLULAR = "GO:0005576"
PLASMA_MEMBRANE = "GO:0005886"
INTERNAL = frozenset({"GO:0005829", "GO:0005634", "GO:0005739"})  # cytosol, nucleus, mito
SECRETORY = frozenset({"GO:0012505", PLASMA_MEMBRANE, "GO:0005773"})  # endomembrane, PM, vacuole
ANY_SECRETORY = SECRETORY | {"GO:0016020", "GO:0071944", WALL, EXTRACELLULAR}
SURFACE = frozenset({WALL, EXTRACELLULAR})
HIGH_THROUGHPUT_CODES = frozenset({"HDA", "HMP", "HEP", "HGI", "HTP"})

P_EXT = "P-ext"
AMBIGUOUS = "ambiguous"
N_INT = "N-int"
N_SEC = "N-sec"
UNLABELLED = "unlabelled"
LABELS = (P_EXT, AMBIGUOUS, N_INT, N_SEC, UNLABELLED)

POLICIES = ("non_iea", "no_homology", "experimental")


def policy_accepts(policy: str, evidence: str) -> bool:
    if policy == "non_iea":
        return evidence != "IEA"
    if policy == "no_homology":
        return evidence != "IEA" and evidence not in HOMOLOGY_CODES
    if policy == "experimental":
        return evidence in EXPERIMENTAL_CODES
    raise ValueError(f"unknown evidence policy {policy!r}")


def classify(label_terms: frozenset[str] | set[str], any_terms: frozenset[str] | set[str]) -> str:
    surface = bool(label_terms & SURFACE)
    internal = bool(label_terms & INTERNAL)
    if surface:
        return AMBIGUOUS if internal else P_EXT
    if internal and not (any_terms & ANY_SECRETORY):
        return N_INT
    if label_terms & SECRETORY and not (any_terms & SURFACE):
        return N_SEC
    return UNLABELLED


def subset_of(label: str, label_terms: frozenset[str] | set[str]) -> str:
    if label != P_EXT:
        return ""
    return "wall" if WALL in label_terms else "extracellular-only"


def htp_only(internal_codes: set[str]) -> bool:
    """True if the gene has non-IEA internal evidence and all of it is high-throughput (R-A)."""
    return bool(internal_codes) and internal_codes <= HIGH_THROUGHPUT_CODES


def is_pm_candidate(label: str, label_terms: frozenset[str] | set[str]) -> bool:
    """P-ext with a non-IEA plasma-membrane term: the input set of D8 triage."""
    return label == P_EXT and PLASMA_MEMBRANE in label_terms
```

- [ ] **Step 4: Write `analysis/step1_compare/truth_table.py` (rows only)**

```python
"""Build the D1 truth table from filtered GAF rows. Standard library only."""

from dataclasses import dataclass, field

import labels

TRUTH_COLUMNS = (
    "source_id",
    "species",
    "taxon_id",
    "in_clade",
    "role",
    "gene_id",
    "symbol",
    "synonym1",
    "label",
    "subset",
    "stratum",
    "tier",
    "label_no_homology",
    "label_experimental",
    "homology_only",
    "pm_candidate",
    "evidence_codes",
    "surface_evidence",
    "internal_evidence",
    "internal_evidence_htp_only",
    "secretory_evidence",
    "source_file",
    "source_sha256",
    "source_date",
    "obo_sha256",
)


@dataclass
class SourceInfo:
    source_id: str
    species: str
    taxon_id: str
    in_clade: str
    role: str
    source_file: str
    source_sha256: str
    source_date: str
    obo_sha256: str


@dataclass
class GeneRecord:
    gene_id: str
    symbol: str
    synonym1: str
    rows: set[tuple[str, str]] = field(default_factory=set)  # (term, evidence)


def collect_genes(filtered) -> dict[str, GeneRecord]:
    genes: dict[str, GeneRecord] = {}
    for r in filtered.cc_rows:
        rec = genes.get(r.gene_id)
        if rec is None:
            rec = GeneRecord(r.gene_id, r.symbol, r.synonyms.split("|")[0])
            genes[r.gene_id] = rec
        rec.rows.add((r.term, r.evidence))
    return genes


def _terms(rec: GeneRecord, ontology, policy: str | None) -> frozenset[str]:
    out: set[str] = set()
    for term, evidence in rec.rows:
        if policy is None or labels.policy_accepts(policy, evidence):
            out |= ontology.ancestors(term)
    return frozenset(out)


def _codes(rec: GeneRecord, ontology, targets: frozenset[str]) -> str:
    codes = {ev for term, ev in rec.rows if ev != "IEA" and ontology.ancestors(term) & targets}
    return ",".join(sorted(codes))


def build_truth_rows(filtered, ontology, info: SourceInfo) -> list[dict[str, str]]:
    rows = []
    for gene_id, rec in sorted(collect_genes(filtered).items()):
        any_terms = _terms(rec, ontology, None)
        by_policy = {p: _terms(rec, ontology, p) for p in labels.POLICIES}
        label = labels.classify(by_policy["non_iea"], any_terms)
        label_nohom = labels.classify(by_policy["no_homology"], any_terms)
        label_exp = labels.classify(by_policy["experimental"], any_terms)
        subset = labels.subset_of(label, by_policy["non_iea"])
        internal = _codes(rec, ontology, labels.INTERNAL)
        rows.append(
            {
                "source_id": info.source_id,
                "species": info.species,
                "taxon_id": info.taxon_id,
                "in_clade": info.in_clade,
                "role": info.role,
                "gene_id": gene_id,
                "symbol": rec.symbol,
                "synonym1": rec.synonym1,
                "label": label,
                "subset": subset,
                "stratum": subset or label,
                "tier": "T-a",
                "label_no_homology": label_nohom,
                "label_experimental": label_exp,
                "homology_only": "yes" if label != label_nohom else "no",
                "pm_candidate": "yes"
                if labels.is_pm_candidate(label, by_policy["non_iea"])
                else "no",
                "evidence_codes": ",".join(sorted({ev for _, ev in rec.rows})),
                "surface_evidence": _codes(rec, ontology, labels.SURFACE),
                "internal_evidence": internal,
                "internal_evidence_htp_only": "yes"
                if labels.htp_only(set(internal.split(",")) - {""})
                else "no",
                "secretory_evidence": _codes(rec, ontology, labels.SECRETORY),
                "source_file": info.source_file,
                "source_sha256": info.source_sha256,
                "source_date": info.source_date,
                "obo_sha256": info.obo_sha256,
            }
        )
    return rows


def in_direct_stratum(row: dict[str, str], label: str) -> bool:
    """Direct-evidence truth (spec 3.3): the label is the same with and without homology codes.

    This is an intersection. Labels recomputed without homology codes are NOT used, because
    that moves genes whose only internal term is IBA or ISS into P-ext (contradicts Q9).
    """
    return row["label"] == label and row["homology_only"] == "no"
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `$PY -m pytest tests/step1_compare/test_labels.py -q`
Expected: PASS (19 passed).

- [ ] **Step 6: Lint and commit**

```bash
pre-commit run --all-files
ruff check .
git add analysis/step1_compare/labels.py analysis/step1_compare/truth_table.py \
        tests/step1_compare/test_labels.py
git commit -m "step1_compare: label rules and truth-table rows (three evidence policies)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Species table, pinned manifest and download script

**Files:**
- Create: `analysis/step1_compare/manifest.py`
- Create: `analysis/step1_compare/manifest.tsv`
- Create: `analysis/step1_compare/species.tsv`
- Create: `analysis/step1_compare/00_fetch_inputs.py`
- Create: `tests/step1_compare/test_fetch.py`

**Interfaces:**
- Consumes: `paths.STEP1_DIR`, `paths.downloads_dir()` (Task 1); `conftest.load_script` (Task 1).
- Produces: `manifest.MANIFEST_COLUMNS`; `manifest.DownloadError`; `manifest.read_manifest(path) -> list[dict[str, str]]`; `manifest.write_manifest(path, rows)`; `manifest.sha256_file(path) -> str`; `manifest.check_payload(path)`; `manifest.verify_against_manifest(path, row) -> str` (returns the SHA-256; raises `DownloadError` on a strict mismatch); `manifest.http_get(url, opener) -> (bytes, dict)`; `manifest.uniprot_pages(url, opener) -> iterator of (bytes, dict)`. `00_fetch_inputs.NetworkError`; `00_fetch_inputs.fetch_one(row, dest_dir, opener) -> dict`; `00_fetch_inputs.main(argv=None, opener=None) -> int` (0 ok, 2 failure; an HTTP or URL error prints `STOP: <url>: <error>`, removes the `.part` file and returns 2 at once). `species.tsv` columns: `source_id, species, taxon_id, taxon_filter, in_clade, role, role_note, gaf_file, fasta_file, id_mapping`.

Manifest modes. `strict`: the 11 GO files. A different SHA-256 stops the run unless `--update-manifest` is given. `record`: the protein FASTA files. PomBase rebuilds `peptide.fa.gz` often (its `Last-Modified` on 2026-10-01 was 02:43 GMT that day), and UniProt changes per release. A different SHA-256 is logged and written to `fetch_log.tsv`, not fatal. Both modes refuse a missing file, an empty file, an HTML page and a truncated gzip. `check_payload` does not catch an empty but valid gzip stream or a truncated plain-text file (for example a cut `go-basic.obo`): for `strict` files the SHA-256 comparison catches both; `record` files have only the gzip end-of-stream check.

The SHA-256 values below come from `sha256sum` on `/tmp/glyco_spec/` (GO files) and from the first run of the script in this task (FASTA files, 2026-10-01). The `last_modified` values are HTTP `Last-Modified` headers read with `curl -sIL` on 2026-10-01. UniProt sends no `Last-Modified`; those rows say `not verified` and give the release (`2026_03`) in `date_generated`.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_fetch.py`:

```python
import gzip
import hashlib
import io
import urllib.error

import manifest
import paths
import pytest
from conftest import load_script

GOOD = gzip.compress(b"!gaf-version: 2.2\nSGD\tS1\n", mtime=0)
GOOD_SHA = hashlib.sha256(GOOD).hexdigest()


class FakeResponse(io.BytesIO):
    def __init__(self, body, headers):
        super().__init__(body)
        self.headers = headers

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def fake_opener(pages):
    """pages: url -> (body, headers). Records every requested URL."""
    calls = []

    def opener(request, timeout=None):
        calls.append(request.full_url)
        body, headers = pages[request.full_url]
        return FakeResponse(body, headers)

    opener.calls = calls
    return opener


def write_manifest(path, rows):
    full = [{c: "" for c in manifest.MANIFEST_COLUMNS} | r for r in rows]
    manifest.write_manifest(path, full)
    return path


def test_check_payload_refuses_html_truncated_and_empty(tmp_path):
    html = tmp_path / "a.gaf.gz"
    html.write_bytes(b"<!DOCTYPE html><html><body>403 Forbidden</body></html>")
    with pytest.raises(manifest.DownloadError, match="HTML"):
        manifest.check_payload(html)
    html_gz = tmp_path / "b.gaf.gz"
    html_gz.write_bytes(gzip.compress(b"<html>error</html>"))
    with pytest.raises(manifest.DownloadError, match="HTML"):
        manifest.check_payload(html_gz)
    cut = tmp_path / "c.gaf.gz"
    noise = b"".join(hashlib.sha256(str(i).encode()).digest() for i in range(3000))
    cut.write_bytes(gzip.compress(noise)[:50000])
    with pytest.raises(manifest.DownloadError, match="truncated"):
        manifest.check_payload(cut)
    empty = tmp_path / "d.gaf.gz"
    empty.write_bytes(b"")
    with pytest.raises(manifest.DownloadError, match="empty"):
        manifest.check_payload(empty)


def test_strict_hash_mismatch_stops_and_keeps_no_file(tmp_path):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/sgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "sgd.gaf.gz", "kind": "http", "url": url, "sha256": "0" * 64, "mode": "strict"}],
    )
    opener = fake_opener({url: (GOOD, {"Last-Modified": "Thu, 28 May 2026 21:18:22 GMT"})})
    dest = tmp_path / "dl"
    assert fetch.main(["--manifest", str(mpath), "--dest", str(dest)], opener=opener) == 2
    assert not (dest / "sgd.gaf.gz").exists()
    assert not (dest / "sgd.gaf.gz.part").exists()


def test_update_manifest_accepts_new_hash(tmp_path):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/sgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "sgd.gaf.gz", "kind": "http", "url": url, "sha256": "0" * 64, "mode": "strict"}],
    )
    opener = fake_opener({url: (GOOD, {"Last-Modified": "Thu, 28 May 2026 21:18:22 GMT"})})
    argv = ["--manifest", str(mpath), "--dest", str(tmp_path / "dl"), "--update-manifest"]
    assert fetch.main(argv, opener=opener) == 0
    row = manifest.read_manifest(mpath)[0]
    assert row["sha256"] == GOOD_SHA
    assert row["last_modified"] == "Thu, 28 May 2026 21:18:22 GMT"
    assert (tmp_path / "dl" / "fetch_log.tsv").read_text().count("sgd.gaf.gz") == 2


def test_html_error_page_stops_even_with_update(tmp_path):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/aspgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "aspgd.gaf.gz", "kind": "http", "url": url, "sha256": "", "mode": "record"}],
    )
    opener = fake_opener({url: (b"<html><body>403</body></html>", {})})
    argv = ["--manifest", str(mpath), "--dest", str(tmp_path / "dl"), "--update-manifest"]
    assert fetch.main(argv, opener=opener) == 2
    assert not (tmp_path / "dl" / "aspgd.gaf.gz").exists()


def test_cached_strict_file_is_not_downloaded_again(tmp_path):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/sgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "sgd.gaf.gz", "kind": "http", "url": url, "sha256": GOOD_SHA, "mode": "strict"}],
    )
    (tmp_path / "dl").mkdir()
    (tmp_path / "dl" / "sgd.gaf.gz").write_bytes(GOOD)
    opener = fake_opener({})
    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 0
    )
    assert opener.calls == []


def test_uniprot_pages_follow_link_header(tmp_path):
    fetch = load_script("00_fetch_inputs")
    first = "https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP1&format=fasta&size=500"
    second = "https://rest.uniprot.org/uniprotkb/search?cursor=abc"
    pages = {
        first: (
            b">sp|P1|A_B\nMK\n",
            {"Link": f'<{second}>; rel="next"', "X-UniProt-Release": "2026_03"},
        ),
        second: (b">tr|Q2|C_D\nMA\n", {}),
    }
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "UP1.fasta.gz", "kind": "uniprot_fasta", "url": first, "mode": "record"}],
    )
    opener = fake_opener(pages)
    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 0
    )
    assert opener.calls == [first, second]
    assert gzip.decompress((tmp_path / "dl" / "UP1.fasta.gz").read_bytes()).count(b">") == 2


def test_http_403_stops_with_exit_2(tmp_path, capsys):
    fetch = load_script("00_fetch_inputs")
    url = "https://example.org/aspgd.gaf.gz"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "aspgd.gaf.gz", "kind": "http", "url": url, "sha256": "", "mode": "strict"}],
    )

    def opener(request, timeout=None):
        raise urllib.error.HTTPError(request.full_url, 403, "Forbidden", {}, None)

    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 2
    )
    assert "STOP:" in capsys.readouterr().err
    assert list((tmp_path / "dl").iterdir()) == []


def test_url_error_mid_download_leaves_no_part_file(tmp_path):
    fetch = load_script("00_fetch_inputs")
    first = "https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP1&format=fasta&size=500"
    second = "https://rest.uniprot.org/uniprotkb/search?cursor=abc"
    mpath = write_manifest(
        tmp_path / "m.tsv",
        [{"file": "UP1.fasta.gz", "kind": "uniprot_fasta", "url": first, "mode": "record"}],
    )

    def opener(request, timeout=None):
        if request.full_url == first:
            return FakeResponse(b">sp|P1|A_B\nMK\n", {"Link": f'<{second}>; rel="next"'})
        raise urllib.error.URLError("connection reset")

    assert (
        fetch.main(["--manifest", str(mpath), "--dest", str(tmp_path / "dl")], opener=opener) == 2
    )
    assert list((tmp_path / "dl").iterdir()) == []


def test_committed_manifest_covers_species_table():
    rows = {r["file"]: r for r in manifest.read_manifest(paths.STEP1_DIR / "manifest.tsv")}
    import csv

    with open(paths.STEP1_DIR / "species.tsv", encoding="utf-8") as handle:
        species = list(csv.DictReader(handle, delimiter="\t"))
    assert "go-basic.obo" in rows
    for sp in species:
        assert rows[sp["gaf_file"]]["mode"] == "strict"
        assert len(rows[sp["gaf_file"]]["sha256"]) == 64
        assert sp["fasta_file"] in rows
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_fetch.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'manifest'`.

- [ ] **Step 3: Write `analysis/step1_compare/manifest.py`**

```python
"""Pinned-input manifest: read it, hash files, validate downloads. Standard library only.

Mode `strict`: a SHA-256 that differs from the manifest stops the run.
Mode `record`: the upstream file changes often (PomBase daily, UniProt per release). A different
SHA-256 is logged, not fatal. The payload checks below still apply to both modes.
"""

import csv
import gzip
import hashlib
import re
import urllib.request
import zlib
from pathlib import Path

MANIFEST_COLUMNS = (
    "file",
    "kind",
    "url",
    "sha256",
    "size",
    "last_modified",
    "date_generated",
    "mode",
    "verified",
)
HTML_PREFIXES = (b"<!doctype html", b"<html", b"<?xml", b"<head")


class DownloadError(RuntimeError):
    """A download is empty, HTML, truncated, or has a SHA-256 that differs from the manifest."""


def read_manifest(path: str | Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        missing = [c for c in MANIFEST_COLUMNS if c not in row]
        if missing:
            raise ValueError(f"{path}: manifest row {row.get('file')} lacks {missing}")
    return rows


def write_manifest(path: str | Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(MANIFEST_COLUMNS), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_payload(path: str | Path) -> None:
    """Refuse a missing file, an empty file, an HTML error page and a truncated gzip file.

    Limits: an empty but valid gzip stream and a truncated plain-text file (for example a cut
    go-basic.obo) pass this check. For `strict` files the SHA-256 comparison catches both. For
    `record` files only the gzip end-of-stream check applies.
    """
    path = Path(path)
    if not path.exists():
        raise DownloadError(f"{path}: file not found")
    if path.stat().st_size == 0:
        raise DownloadError(f"{path.name}: empty file")
    with open(path, "rb") as handle:
        head = handle.read(512)
    if head[:2] == b"\x1f\x8b":
        try:
            with gzip.open(path, "rb") as gz:
                first = gz.read(512)
                while gz.read(1 << 20):
                    pass
        except (EOFError, OSError, zlib.error) as exc:
            raise DownloadError(f"{path.name}: truncated or corrupt gzip ({exc})") from exc
        head = first
    if head.lstrip().lower().startswith(HTML_PREFIXES):
        raise DownloadError(f"{path.name}: HTML page, not data")


def verify_against_manifest(path: str | Path, row: dict[str, str]) -> str:
    """Return the file's SHA-256. Raise DownloadError on a strict-mode mismatch."""
    check_payload(path)
    actual = sha256_file(path)
    if row["mode"] == "strict" and actual != row["sha256"]:
        raise DownloadError(
            f"{Path(path).name}: SHA-256 {actual[:12]} differs from manifest {row['sha256'][:12]}"
        )
    return actual


def http_get(url: str, opener=urllib.request.urlopen):
    """Return (body bytes, headers dict with lower-case keys)."""
    request = urllib.request.Request(url, headers={"User-Agent": "adhesionPred-step1/1"})
    with opener(request, timeout=300) as response:
        body = response.read()
        headers = {k.lower(): v for k, v in response.headers.items()}
    return body, headers


_NEXT = re.compile(r'<([^>]+)>;\s*rel="next"')


def uniprot_pages(url: str, opener=urllib.request.urlopen):
    """Yield (body, headers) for each page of a UniProt REST search, following Link rel=next."""
    while url:
        body, headers = http_get(url, opener)
        yield body, headers
        match = _NEXT.search(headers.get("link", ""))
        url = match.group(1) if match else ""
```

- [ ] **Step 4: Write `analysis/step1_compare/00_fetch_inputs.py`**

```python
#!/usr/bin/env python3
"""Fetch the pinned step 1 inputs into the work directory and check them against manifest.tsv.

Environment: PROJ_ROOT (repository root) and STEP1_WORKDIR (output root). Files go to
$STEP1_WORKDIR/downloads. A strict-mode SHA-256 mismatch stops the run (exit 2) unless
--update-manifest is given. Small downloads only: this needs no SLURM job.
"""

import argparse
import datetime
import gzip
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

import manifest
import paths


class NetworkError(RuntimeError):
    """An HTTP error, a URL error or a timeout during a download."""


def fetch_one(row: dict[str, str], dest_dir: Path, opener) -> dict[str, str]:
    dest = dest_dir / row["file"]
    part = dest.with_name(dest.name + ".part")
    headers: dict[str, str] = {}
    try:
        if row["kind"] == "http":
            body, headers = manifest.http_get(row["url"], opener)
            part.write_bytes(body)
        elif row["kind"] == "uniprot_fasta":
            with (
                open(part, "wb") as raw,
                gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz,
            ):
                for body, page_headers in manifest.uniprot_pages(row["url"], opener):
                    gz.write(body)
                    headers = headers or page_headers
        else:
            raise ValueError(f"{row['file']}: unknown kind {row['kind']!r}")
    except (urllib.error.URLError, TimeoutError) as exc:
        part.unlink(missing_ok=True)
        raise NetworkError(f"{row['url']}: {exc}") from exc
    try:
        manifest.check_payload(part)
    except manifest.DownloadError:
        part.unlink()
        raise
    sha = manifest.sha256_file(part)
    return {
        "part": str(part),
        "dest": str(dest),
        "sha256": sha,
        "size": str(part.stat().st_size),
        "last_modified": headers.get("last-modified", ""),
        "http_date": headers.get("date", ""),
        "uniprot_release": headers.get("x-uniprot-release", ""),
    }


def main(argv=None, opener=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(paths.STEP1_DIR / "manifest.tsv"))
    parser.add_argument("--dest", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument("--only", nargs="*", default=None, help="file names to fetch")
    parser.add_argument("--update-manifest", action="store_true")
    args = parser.parse_args(argv)

    opener = opener or urllib.request.urlopen
    dest_dir = Path(args.dest) if args.dest else paths.downloads_dir()
    dest_dir.mkdir(parents=True, exist_ok=True)
    rows = manifest.read_manifest(args.manifest)
    today = datetime.datetime.now(datetime.UTC).date().isoformat()
    log_rows, failed, changed = [], [], False
    for row in rows:
        if args.only and row["file"] not in args.only:
            continue
        dest = dest_dir / row["file"]
        if (
            dest.exists()
            and row["mode"] == "strict"
            and manifest.sha256_file(dest) == row["sha256"]
        ):
            print(f"ok (cached) {row['file']}")
            continue
        try:
            got = fetch_one(row, dest_dir, opener)
        except NetworkError as exc:
            print(f"STOP: {exc}", file=sys.stderr)
            return 2
        except manifest.DownloadError as exc:
            print(f"FAILED {exc}", file=sys.stderr)
            failed.append(row["file"])
            continue
        mismatch = got["sha256"] != row["sha256"]
        if mismatch and row["mode"] == "strict" and not args.update_manifest:
            Path(got["part"]).unlink()
            print(
                f"FAILED {row['file']}: SHA-256 {got['sha256'][:12]} differs from manifest "
                f"{row['sha256'][:12]}; rerun with --update-manifest to accept it",
                file=sys.stderr,
            )
            failed.append(row["file"])
            continue
        os.replace(got["part"], got["dest"])
        if mismatch:
            print(f"changed {row['file']}: {row['sha256'][:12]} -> {got['sha256'][:12]}")
        if mismatch and args.update_manifest:
            row.update(
                sha256=got["sha256"],
                size=got["size"],
                last_modified=got["last_modified"] or "not verified",
                verified=today,
            )
            changed = True
        log_rows.append({"file": row["file"], **{k: got[k] for k in got if k not in ("part",)}})
    if changed:
        manifest.write_manifest(args.manifest, rows)
    if log_rows:
        log = dest_dir / "fetch_log.tsv"
        write_header = not log.exists()
        with open(log, "a", encoding="utf-8") as handle:
            cols = [
                "file",
                "dest",
                "sha256",
                "size",
                "last_modified",
                "http_date",
                "uniprot_release",
            ]
            if write_header:
                handle.write("\t".join(["fetched_on", *cols]) + "\n")
            for entry in log_rows:
                handle.write("\t".join([today, *(entry.get(c, "") for c in cols)]) + "\n")
    return 2 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Write `analysis/step1_compare/species.tsv` (tab-separated)**

`role` values: `train` (Q6), `test_species` (Q6), `test_clade` (Q7; *A. nidulans* and H99 recorded on 2026-10-01; H99 truth is pending curation, which is separate work), `alternate_file` (another file or strain kept for comparison: SCHPO-mod, JEC21, CRYD1), `undecided` (*U. maydis*: the Basidiomycota option is open).

```text
source_id	species	taxon_id	taxon_filter	in_clade	role	role_note	gaf_file	fasta_file	id_mapping
Scer_SGD	Saccharomyces cerevisiae S288C	559292		Saccharomycotina	train	Q6	sgd.gaf.gz	orf_trans_all.fasta.gz	sgd
Calb_CGD	Candida albicans SC5314	237561	237561	Saccharomycotina	train	Q6	cgd.gaf.gz	C_albicans_SC5314_A22_current_default_protein.fasta.gz	cgd
Spom_PomBase	Schizosaccharomyces pombe 972h-	284812		Taphrinomycotina	test_species	Q6	pombase.gaf.gz	peptide.fa.gz	pombase
Spom_SCHPO-mod	Schizosaccharomyces pombe 972h-	284812		Taphrinomycotina	alternate_file	spec 3.2	SCHPO-mod.gaf.gz	peptide.fa.gz	pombase
Afum_ASPFU	Aspergillus fumigatus Af293	330879		Eurotiomycetes	test_clade	Q7	ASPFU-uniprot.gaf.gz	UP000002530.fasta.gz	uniprot
Anid_EMENI	Aspergillus nidulans FGSC A4	227321		Eurotiomycetes	test_clade	decided 2026-10-01: second Eurotiomycetes test species	EMENI-uniprot.gaf.gz	UP000000560.fasta.gz	uniprot
Cneo_H99_GOA	Cryptococcus neoformans H99	235443		Basidiomycota	test_clade	decided 2026-10-01: Cryptococcus reference; pending curation (separate work)	313589.C_neoformans_var_grubii_H99.goa	UP000010091.fasta.gz	uniprot
Cneo_JEC21_GOA	Cryptococcus deneoformans JEC21	214684		Basidiomycota	alternate_file	H99 is the reference; JEC21 kept for comparison	20846.C_neoformans_JEC21.goa	UP000002149.fasta.gz	uniprot
Cneo_CRYD1	Cryptococcus deneoformans JEC21	214684		Basidiomycota	alternate_file	same genome as Cneo_JEC21_GOA	CRYD1-uniprot.gaf.gz	UP000002149.fasta.gz	uniprot
Umay_MYCMD	Ustilago maydis 521	5270		Basidiomycota	undecided	Basidiomycota truth option is open	MYCMD-uniprot.gaf.gz	UP000000561.fasta.gz	uniprot
```

- [ ] **Step 6: Write `analysis/step1_compare/manifest.tsv` (tab-separated)**

```text
file	kind	url	sha256	size	last_modified	date_generated	mode	verified
go-basic.obo	http	https://current.geneontology.org/ontology/go-basic.obo	b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc	32227785	Sat, 08 Aug 2026 22:27:06 GMT	releases/2026-07-26	strict	2026-10-01
sgd.gaf.gz	http	https://current.geneontology.org/annotations/sgd.gaf.gz	f9c603d5fc358d7a2fa3dc6cd8d354d0723a13e1dfa92e608091d11943d4559e	2602501	Thu, 28 May 2026 21:18:22 GMT	2026-05-21T09:00	strict	2026-10-01
cgd.gaf.gz	http	https://current.geneontology.org/annotations/cgd.gaf.gz	a5d44d9395b2fb7b50822a9175a9488260ddf4eae8600d17d98711bea9adb41b	4601396	Thu, 28 May 2026 21:18:36 GMT	2026-05-21T07:35	strict	2026-10-01
pombase.gaf.gz	http	https://current.geneontology.org/annotations/pombase.gaf.gz	174c236b98191d378bd0f6335625fed6fd33364f9c3095f375b16bedbcd20b1d	1990730	Thu, 28 May 2026 21:18:25 GMT	2026-05-21T08:42	strict	2026-10-01
SCHPO-mod.gaf.gz	http	https://current.geneontology.org/annotations/gaf/SCHPO-mod.gaf.gz	2ac1820b7e62480b97f3de91db4c07cfc033586b2b29a4956f3fdadfb300bc9d	1542262	Sat, 08 Aug 2026 22:15:59 GMT	2026-07-28 16:39	strict	2026-10-01
ASPFU-uniprot.gaf.gz	http	https://current.geneontology.org/annotations/gaf/ASPFU-uniprot.gaf.gz	407a6eef89a6962a28a497924af2152d74bbc7a9304fe558f6fbf0dd5ff096a4	1304170	Sat, 08 Aug 2026 22:16:07 GMT	2026-07-28 19:47	strict	2026-10-01
EMENI-uniprot.gaf.gz	http	https://current.geneontology.org/annotations/gaf/EMENI-uniprot.gaf.gz	e9efa6e27ff4e3784b3aacc0a57076453049577b8cfa604f75c9bd1c05fd95ee	1377233	Sat, 08 Aug 2026 22:15:48 GMT	2026-07-28 19:39	strict	2026-10-01
CRYD1-uniprot.gaf.gz	http	https://current.geneontology.org/annotations/gaf/CRYD1-uniprot.gaf.gz	a7702530d23531375d0c5d8025638eb681e8035f78f2e458d308f76de7447c90	953586	Sat, 08 Aug 2026 22:15:43 GMT	2026-07-28 19:38	strict	2026-10-01
MYCMD-uniprot.gaf.gz	http	https://current.geneontology.org/annotations/gaf/MYCMD-uniprot.gaf.gz	d6ab3b5a0cc94109a398b8e495f51f66a4b2c0f554880685508cbb82fe3c2618	951232	Sat, 08 Aug 2026 22:16:16 GMT	2026-07-28 19:54	strict	2026-10-01
20846.C_neoformans_JEC21.goa	http	https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/20846.C_neoformans_JEC21.goa	851e69c7d16f5c9a4282d4e312eb3e794c739961e3801ce98e6e6ee524290a9a	8984355	Tue, 28 Jul 2026 20:48:14 GMT	2026-07-28	strict	2026-10-01
313589.C_neoformans_var_grubii_H99.goa	http	https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/313589.C_neoformans_var_grubii_H99.goa	528ba85cba1c2569dfae5fcaee120c94dcb580ea2edf606a164b07463937509e	6353529	Tue, 28 Jul 2026 20:48:30 GMT	2026-07-28	strict	2026-10-01
orf_trans_all.fasta.gz	http	https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz	17e8b47e1ae23178c6000fbc4ab548f102d1b250ef9dff5d811feb3f03dd2c5b	2689634	Thu, 13 Jun 2024 22:28:53 GMT	Genome Release 64-5-1	record	2026-10-01
C_albicans_SC5314_A22_current_default_protein.fasta.gz	http	https://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_default_protein.fasta.gz	cd8da74fab4fbc479572a3e2d75c76c661d0317591360f6f5309248b698ee7ad	1809513	Sun, 27 Sep 2026 03:02:58 GMT	Assembly 22	record	2026-10-01
peptide.fa.gz	http	https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/peptide.fa.gz	76c1eb8212c193ed24143ab2693e19ec7ca9a2f2e9127069098b7409cc9ff2d6	1620472	Thu, 01 Oct 2026 02:43:23 GMT		record	2026-10-01
UP000002530.fasta.gz	uniprot_fasta	https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP000002530&format=fasta&size=500	7cfae982f8d29b5b7c51ad85541e04c8d78300265e84390bb657b87577fad84c	3211583	not verified	UniProt 2026_03	record	2026-10-01
UP000000560.fasta.gz	uniprot_fasta	https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP000000560&format=fasta&size=500	ce6b98cb7b74cf3f3da93c175a4aa7858f6d1a33d7056a3efa4182ff80d2ef03	3463564	not verified	UniProt 2026_03	record	2026-10-01
UP000010091.fasta.gz	uniprot_fasta	https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP000010091&format=fasta&size=500	eadcc23f379f8cca60c99617212b8ac2ae0b15e464e4058a549067dfdfd595af	2594162	not verified	UniProt 2026_03	record	2026-10-01
UP000002149.fasta.gz	uniprot_fasta	https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP000002149&format=fasta&size=500	d893b06ed5e98b7c65ea17af8c4935c790a793588e43c5f5fda6b29878b86b0e	2407035	not verified	UniProt 2026_03	record	2026-10-01
UP000000561.fasta.gz	uniprot_fasta	https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP000000561&format=fasta&size=500	b181fadc3185189adb1d2e2cc04c1002ac18b11c2f0aa9df61bfe6934285cf76	2625166	not verified	UniProt 2026_03	record	2026-10-01
```

- [ ] **Step 7: Run the test to verify it passes**

Run: `$PY -m pytest tests/step1_compare/test_fetch.py -q`
Expected: PASS (9 passed).

- [ ] **Step 8: Run the real download (network; about 3 minutes)**

```bash
export PROJ_ROOT=$PWD STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
cd analysis/step1_compare && /usr/bin/python3.12 00_fetch_inputs.py; echo "exit=$?"; cd -
```

Expected: `exit=0`. The 11 GO files pass the strict check (on 2026-10-01 a fresh download gave the pinned hashes). For a `record` file whose upstream copy changed after 2026-10-01, the script prints `changed <file>: <old> -> <new>` and logs it. Read `$STEP1_WORKDIR/downloads/fetch_log.tsv`. If a strict file fails, do not pass `--update-manifest` without the owner's agreement: a new GAF changes the spec 3.3 numbers.

- [ ] **Step 9: Lint and commit**

```bash
pre-commit run --all-files
ruff check .
git add analysis/step1_compare/manifest.py analysis/step1_compare/manifest.tsv \
        analysis/step1_compare/species.tsv analysis/step1_compare/00_fetch_inputs.py \
        tests/step1_compare/test_fetch.py
git commit -m "step1_compare: species table, pinned input manifest, verified downloader

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 6: D1 extractor command and the spec 3.3 regression test

**Files:**
- Modify: `analysis/step1_compare/truth_table.py` (replace with the complete version below)
- Create: `analysis/step1_compare/01_extract_go_truth.py`
- Create: `tests/step1_compare/test_extract.py`

**Interfaces:**
- Consumes: `go_obo.parse_obo` (Task 2); `gaf.filter_gaf`, `gaf.header_value` (Task 3); `truth_table.build_truth_rows`, `truth_table.SourceInfo`, `labels.*` (Task 4); `manifest.verify_against_manifest`, `manifest.read_manifest`, `species.tsv`, `manifest.tsv` (Task 5); `paths.*` (Task 1).
- Produces: `truth_table.COUNT_COLUMNS` (includes `direct_p_ext`, `direct_n_int`, `direct_n_sec`, `direct_ambiguous`, `ambiguous_htp_only`); `truth_table.count_rows(rows, filtered, source_id) -> dict[str, str]`; `truth_table.write_tsv(path, columns, rows)` (a `.gz` path is gzip with `mtime=0` and no file name, so output is byte-stable); `truth_table.read_tsv(path) -> list[dict[str, str]]`. `01_extract_go_truth.read_species(path) -> list[dict[str, str]]`; `01_extract_go_truth.run(species_rows, manifest_rows, input_dir: Path, out_dir: Path, obo_name="go-basic.obo") -> (truth_rows, count_rows)`; `main(argv=None) -> int`. Files: `$STEP1_WORKDIR/truth_set.tsv.gz`, `counts.tsv`, `extract_log.json`.

Fixture counts (hand-computed from the 50 rows): genes 20, non-IEA genes 19, P-ext 7 (wall 4, extracellular-only 3), N-int 1, N-sec 4, ambiguous 3, PM candidates 1, CC triples 35 of which 7 IEA (0.200), all-aspect IEA fraction 9/43 = 0.209, obsolete rows 1; experimental-only 6 / 1 / 4 / 2; recomputed without homology codes 6 / 1 / 4 / 2; direct-evidence intersection 5 / 1 / 4 / 2 (G10, G11, G16 drop out); ambiguous with high-throughput-only internal evidence 1 (G14).

A GAF named in `species.tsv` but absent on disk stops with `STOP: <path>: file not found` and exit 2 (`check_payload` checks existence).

Two regression tests share one real extraction (module-scoped fixture). `test_reproduces_spec_counts_on_real_files` holds the `d1_count.py` numbers, the SCHPO-mod and CRYD1 rows, and `genes_cc` for the five UniProt-based sources. `test_direct_evidence_counts_on_real_files` holds the spec 3.3 direct-evidence counts. They run when the real files are in `$STEP1_GO_DIR` or `$STEP1_WORKDIR/downloads`. CI has no real files, so CI skips it and runs the fixture tests.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_extract.py`:

```python
import os
import shutil
from pathlib import Path

import manifest
import paths
import pytest
import truth_table
from conftest import load_script
from gaf_fixture import write_golden_gaf

extract = load_script("01_extract_go_truth")

FIXTURE_SPECIES = [
    {
        "source_id": "Fix_SGD",
        "species": "Fixture yeast",
        "taxon_id": "559292",
        "taxon_filter": "",
        "in_clade": "Saccharomycotina",
        "role": "train",
        "gaf_file": "golden.gaf",
    }
]


def _inputs(tmp_path, fixtures_dir):
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    shutil.copy(fixtures_dir / "mini.obo", input_dir / "go-basic.obo")
    write_golden_gaf(input_dir / "golden.gaf")
    rows = [
        {"file": f, "mode": "strict", "sha256": manifest.sha256_file(input_dir / f)}
        for f in ("go-basic.obo", "golden.gaf")
    ]
    return input_dir, rows


def test_extract_writes_truth_and_counts(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    out = tmp_path / "out"
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, out)
    truth = truth_table.read_tsv(out / "truth_set.tsv.gz")
    assert tuple(truth[0].keys()) == truth_table.TRUTH_COLUMNS
    assert len(truth) == 20
    counts = truth_table.read_tsv(out / "counts.tsv")[0]
    expected = {
        "primary_db": "SGD",
        "genes_cc": "20",
        "genes_noniea_cc": "19",
        "p_ext": "7",
        "p_ext_wall": "4",
        "p_ext_extonly": "3",
        "n_int": "1",
        "n_sec": "4",
        "ambiguous": "3",
        "pm_candidates": "1",
        "cc_iea_triples": "7",
        "cc_triples": "35",
        "cc_iea_frac": "0.200",
        "all_aspects_iea_frac": "0.209",
        "obsolete_rows": "1",
        "unknown_term_rows": "0",
        "exp_p_ext": "6",
        "exp_n_int": "1",
        "exp_n_sec": "4",
        "exp_ambiguous": "2",
        "nohom_p_ext": "6",
        "nohom_n_int": "1",
        "nohom_n_sec": "4",
        "nohom_ambiguous": "2",
        "direct_p_ext": "5",
        "direct_n_int": "1",
        "direct_n_sec": "4",
        "direct_ambiguous": "2",
        "ambiguous_htp_only": "1",
    }
    assert {k: counts[k] for k in expected} == expected
    assert truth[0]["source_sha256"] == manifest_rows[1]["sha256"]
    assert truth[0]["source_date"] == "2026-05-21T09:00"


def test_extract_output_is_byte_stable(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "a")
    extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "b")
    a = (tmp_path / "a" / "truth_set.tsv.gz").read_bytes()
    assert a == (tmp_path / "b" / "truth_set.tsv.gz").read_bytes()


def test_extract_stops_on_hash_mismatch(tmp_path, fixtures_dir):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    with open(input_dir / "golden.gaf", "a") as handle:
        handle.write(
            "SGD\tG99\tNEW\tlocated_in\tGO:0005576\tPMID:1\tIDA\t\tC\tn\t\tprotein"
            "\ttaxon:559292\t20260101\tSGD\t\t\n"
        )
    with pytest.raises(manifest.DownloadError, match="differs from manifest"):
        extract.run(FIXTURE_SPECIES, manifest_rows, input_dir, tmp_path / "out")
    assert not (tmp_path / "out" / "truth_set.tsv.gz").exists()


def test_missing_gaf_stops_with_exit_2(tmp_path, fixtures_dir, capsys):
    input_dir, manifest_rows = _inputs(tmp_path, fixtures_dir)
    (input_dir / "golden.gaf").unlink()
    species_path = tmp_path / "species.tsv"
    truth_table.write_tsv(species_path, list(FIXTURE_SPECIES[0]), FIXTURE_SPECIES)
    mpath = tmp_path / "manifest.tsv"
    full = [{c: "" for c in manifest.MANIFEST_COLUMNS} | r for r in manifest_rows]
    manifest.write_manifest(mpath, full)
    argv = ["--species", str(species_path), "--manifest", str(mpath)]
    argv += ["--input-dir", str(input_dir), "--out-dir", str(tmp_path / "out")]
    assert extract.main(argv) == 2
    err = capsys.readouterr().err
    assert "STOP:" in err and "golden.gaf" in err and "file not found" in err


# Spec 3.3 numbers, reproduced from /tmp/glyco_spec/d1_count.py on 2026-10-01.
SPEC_COUNTS = {
    "Scer_SGD": {
        "genes_cc": 6056,
        "genes_noniea_cc": 6012,
        "p_ext": 125,
        "p_ext_wall": 104,
        "p_ext_extonly": 21,
        "n_int": 2645,
        "n_sec": 1548,
        "ambiguous": 48,
        "cc_iea_frac": "0.262",
        "nohom_p_ext": 96,
        "nohom_n_int": 2360,
        "nohom_n_sec": 1457,
        "pm_candidates": 12,
    },
    "Calb_CGD": {
        "genes_cc": 6313,
        "genes_noniea_cc": 6160,
        "p_ext": 259,
        "p_ext_wall": 122,
        "p_ext_extonly": 137,
        "n_int": 1772,
        "n_sec": 961,
        "ambiguous": 100,
        "cc_iea_frac": "0.158",
        "nohom_p_ext": 294,
        "nohom_n_int": 179,
        "nohom_n_sec": 263,
        "pm_candidates": 60,
    },
    "Spom_PomBase": {
        "genes_cc": 5025,
        "genes_noniea_cc": 4998,
        "p_ext": 59,
        "p_ext_wall": 40,
        "p_ext_extonly": 19,
        "n_int": 2962,
        "n_sec": 1214,
        "ambiguous": 14,
        "cc_iea_frac": "0.292",
        "nohom_p_ext": 44,
        "nohom_n_int": 2799,
        "nohom_n_sec": 886,
        "pm_candidates": 10,
    },
    "Spom_SCHPO-mod": {"p_ext": 57, "n_int": 2966, "n_sec": 1212, "ambiguous": 14},
    "Cneo_H99_GOA": {
        "genes_cc": 4320,
        "genes_noniea_cc": 77,
        "p_ext": 11,
        "p_ext_wall": 3,
        "p_ext_extonly": 8,
        "n_int": 23,
        "n_sec": 25,
        "ambiguous": 0,
        "cc_iea_frac": "0.986",
        "nohom_p_ext": 9,
        "nohom_n_int": 17,
        "nohom_n_sec": 15,
    },
    "Cneo_JEC21_GOA": {
        "p_ext": 32,
        "p_ext_wall": 10,
        "p_ext_extonly": 22,
        "n_int": 1813,
        "n_sec": 817,
        "ambiguous": 0,
        "cc_iea_frac": "0.520",
        "nohom_p_ext": 0,
        "nohom_n_int": 3,
        "nohom_n_sec": 1,
    },
    "Cneo_CRYD1": {"p_ext": 32, "n_int": 1813, "n_sec": 817, "ambiguous": 0},
    "Afum_ASPFU": {
        "p_ext": 132,
        "p_ext_wall": 35,
        "p_ext_extonly": 97,
        "n_int": 2018,
        "n_sec": 1062,
        "ambiguous": 1,
        "cc_iea_frac": "0.533",
        "nohom_p_ext": 23,
        "nohom_n_int": 18,
        "nohom_n_sec": 26,
        "pm_candidates": 9,
    },
    "Anid_EMENI": {
        "p_ext": 211,
        "p_ext_wall": 43,
        "p_ext_extonly": 168,
        "n_int": 2054,
        "n_sec": 1066,
        "ambiguous": 38,
        "cc_iea_frac": "0.517",
        "nohom_p_ext": 148,
        "nohom_n_int": 97,
        "nohom_n_sec": 66,
    },
    "Umay_MYCMD": {
        "p_ext": 62,
        "p_ext_wall": 14,
        "p_ext_extonly": 48,
        "n_int": 1712,
        "n_sec": 827,
        "ambiguous": 1,
        "cc_iea_frac": "0.519",
        "nohom_p_ext": 10,
        "nohom_n_int": 6,
        "nohom_n_sec": 21,
    },
}


# M2: genes_cc for the five UniProt-based sources.
GENES_CC_UNIPROT = {
    "Cneo_H99_GOA": 4320,
    "Cneo_JEC21_GOA": 4364,
    "Cneo_CRYD1": 4317,
    "Afum_ASPFU": 5613,
    "Anid_EMENI": 5884,
    "Umay_MYCMD": 4206,
}

# Spec 3.3 direct-evidence truth: label == X and homology_only == "no" (P-ext, N-int, N-sec).
DIRECT_COUNTS = {
    "Scer_SGD": (88, 2360, 1457),
    "Calb_CGD": (211, 179, 263),
    "Spom_PomBase": (42, 2799, 886),
    "Cneo_H99_GOA": (9, 17, 15),
    "Cneo_JEC21_GOA": (0, 3, 1),
    "Afum_ASPFU": (23, 18, 26),
    "Anid_EMENI": (113, 97, 66),
    "Umay_MYCMD": (10, 6, 21),
}


def _real_input_dir() -> Path:
    return Path(os.environ.get("STEP1_GO_DIR") or paths.downloads_dir())


needs_real_files = pytest.mark.skipif(
    not (_real_input_dir() / "go-basic.obo").exists(),
    reason="real GO files absent; run 00_fetch_inputs.py or set STEP1_GO_DIR",
)


@pytest.fixture(scope="module")
def real_counts(tmp_path_factory):
    species = extract.read_species(paths.STEP1_DIR / "species.tsv")
    manifest_rows = manifest.read_manifest(paths.STEP1_DIR / "manifest.tsv")
    out = tmp_path_factory.mktemp("real")
    _, counts = extract.run(species, manifest_rows, _real_input_dir(), out)
    return {c["source_id"]: c for c in counts}


@needs_real_files
def test_reproduces_spec_counts_on_real_files(real_counts):
    for source_id, expected in SPEC_COUNTS.items():
        for key, value in expected.items():
            assert real_counts[source_id][key] == str(value), (source_id, key)
    for source_id, value in GENES_CC_UNIPROT.items():
        assert real_counts[source_id]["genes_cc"] == str(value), source_id


@needs_real_files
def test_direct_evidence_counts_on_real_files(real_counts):
    for source_id, (p_ext, n_int, n_sec) in DIRECT_COUNTS.items():
        got = real_counts[source_id]
        assert (got["direct_p_ext"], got["direct_n_int"], got["direct_n_sec"]) == (
            str(p_ext),
            str(n_int),
            str(n_sec),
        ), source_id
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_extract.py -q`
Expected: FAIL (collection error) with `FileNotFoundError` for `01_extract_go_truth.py`.

- [ ] **Step 3: Replace `analysis/step1_compare/truth_table.py` with the complete version**

```python
"""Build the D1 truth table and its counts from filtered GAF rows. Standard library only."""

import csv
import gzip
import io
from dataclasses import dataclass, field
from pathlib import Path

import labels

TRUTH_COLUMNS = (
    "source_id",
    "species",
    "taxon_id",
    "in_clade",
    "role",
    "gene_id",
    "symbol",
    "synonym1",
    "label",
    "subset",
    "stratum",
    "tier",
    "label_no_homology",
    "label_experimental",
    "homology_only",
    "pm_candidate",
    "evidence_codes",
    "surface_evidence",
    "internal_evidence",
    "internal_evidence_htp_only",
    "secretory_evidence",
    "source_file",
    "source_sha256",
    "source_date",
    "obo_sha256",
)

COUNT_COLUMNS = (
    "source_id",
    "primary_db",
    "genes_cc",
    "genes_noniea_cc",
    "p_ext",
    "p_ext_wall",
    "p_ext_extonly",
    "n_int",
    "n_sec",
    "ambiguous",
    "pm_candidates",
    "cc_iea_triples",
    "cc_triples",
    "cc_iea_frac",
    "all_aspects_iea_frac",
    "obsolete_rows",
    "unknown_term_rows",
    "exp_p_ext",
    "exp_n_int",
    "exp_n_sec",
    "exp_ambiguous",
    "nohom_p_ext",
    "nohom_n_int",
    "nohom_n_sec",
    "nohom_ambiguous",
    "direct_p_ext",
    "direct_n_int",
    "direct_n_sec",
    "direct_ambiguous",
    "ambiguous_htp_only",
)


@dataclass
class SourceInfo:
    source_id: str
    species: str
    taxon_id: str
    in_clade: str
    role: str
    source_file: str
    source_sha256: str
    source_date: str
    obo_sha256: str


@dataclass
class GeneRecord:
    gene_id: str
    symbol: str
    synonym1: str
    rows: set[tuple[str, str]] = field(default_factory=set)  # (term, evidence)


def collect_genes(filtered) -> dict[str, GeneRecord]:
    genes: dict[str, GeneRecord] = {}
    for r in filtered.cc_rows:
        rec = genes.get(r.gene_id)
        if rec is None:
            rec = GeneRecord(r.gene_id, r.symbol, r.synonyms.split("|")[0])
            genes[r.gene_id] = rec
        rec.rows.add((r.term, r.evidence))
    return genes


def _terms(rec: GeneRecord, ontology, policy: str | None) -> frozenset[str]:
    out: set[str] = set()
    for term, evidence in rec.rows:
        if policy is None or labels.policy_accepts(policy, evidence):
            out |= ontology.ancestors(term)
    return frozenset(out)


def _codes(rec: GeneRecord, ontology, targets: frozenset[str]) -> str:
    codes = {ev for term, ev in rec.rows if ev != "IEA" and ontology.ancestors(term) & targets}
    return ",".join(sorted(codes))


def build_truth_rows(filtered, ontology, info: SourceInfo) -> list[dict[str, str]]:
    rows = []
    for gene_id, rec in sorted(collect_genes(filtered).items()):
        any_terms = _terms(rec, ontology, None)
        by_policy = {p: _terms(rec, ontology, p) for p in labels.POLICIES}
        label = labels.classify(by_policy["non_iea"], any_terms)
        label_nohom = labels.classify(by_policy["no_homology"], any_terms)
        label_exp = labels.classify(by_policy["experimental"], any_terms)
        subset = labels.subset_of(label, by_policy["non_iea"])
        internal = _codes(rec, ontology, labels.INTERNAL)
        rows.append(
            {
                "source_id": info.source_id,
                "species": info.species,
                "taxon_id": info.taxon_id,
                "in_clade": info.in_clade,
                "role": info.role,
                "gene_id": gene_id,
                "symbol": rec.symbol,
                "synonym1": rec.synonym1,
                "label": label,
                "subset": subset,
                "stratum": subset or label,
                "tier": "T-a",
                "label_no_homology": label_nohom,
                "label_experimental": label_exp,
                "homology_only": "yes" if label != label_nohom else "no",
                "pm_candidate": "yes"
                if labels.is_pm_candidate(label, by_policy["non_iea"])
                else "no",
                "evidence_codes": ",".join(sorted({ev for _, ev in rec.rows})),
                "surface_evidence": _codes(rec, ontology, labels.SURFACE),
                "internal_evidence": internal,
                "internal_evidence_htp_only": "yes"
                if labels.htp_only(set(internal.split(",")) - {""})
                else "no",
                "secretory_evidence": _codes(rec, ontology, labels.SECRETORY),
                "source_file": info.source_file,
                "source_sha256": info.source_sha256,
                "source_date": info.source_date,
                "obo_sha256": info.obo_sha256,
            }
        )
    return rows


def in_direct_stratum(row: dict[str, str], label: str) -> bool:
    """Direct-evidence truth (spec 3.3): the label is the same with and without homology codes.

    This is an intersection. Labels recomputed without homology codes are NOT used, because
    that moves genes whose only internal term is IBA or ISS into P-ext (contradicts Q9).
    """
    return row["label"] == label and row["homology_only"] == "no"


def _n(rows, column, value) -> int:
    return sum(1 for r in rows if r[column] == value)


def count_rows(rows: list[dict[str, str]], filtered, source_id: str) -> dict[str, str]:
    triples = {(r.gene_id, r.term, r.evidence) for r in filtered.cc_rows}
    iea = sum(1 for t in triples if t[2] == "IEA")
    iea_all = sum(1 for t in filtered.all_aspect_triples if t[2] == "IEA")
    counts = {
        "source_id": source_id,
        "primary_db": filtered.primary_db,
        "genes_cc": len(rows),
        "genes_noniea_cc": sum(1 for r in rows if r["evidence_codes"] != "IEA"),
        "p_ext": _n(rows, "label", labels.P_EXT),
        "p_ext_wall": _n(rows, "subset", "wall"),
        "p_ext_extonly": _n(rows, "subset", "extracellular-only"),
        "n_int": _n(rows, "label", labels.N_INT),
        "n_sec": _n(rows, "label", labels.N_SEC),
        "ambiguous": _n(rows, "label", labels.AMBIGUOUS),
        "pm_candidates": _n(rows, "pm_candidate", "yes"),
        "cc_iea_triples": iea,
        "cc_triples": len(triples),
        "cc_iea_frac": f"{iea / len(triples):.3f}",
        "all_aspects_iea_frac": f"{iea_all / len(filtered.all_aspect_triples):.3f}",
        "obsolete_rows": filtered.dropped["obsolete_term"],
        "unknown_term_rows": filtered.unknown_term_rows,
        "exp_p_ext": _n(rows, "label_experimental", labels.P_EXT),
        "exp_n_int": _n(rows, "label_experimental", labels.N_INT),
        "exp_n_sec": _n(rows, "label_experimental", labels.N_SEC),
        "exp_ambiguous": _n(rows, "label_experimental", labels.AMBIGUOUS),
        "nohom_p_ext": _n(rows, "label_no_homology", labels.P_EXT),
        "nohom_n_int": _n(rows, "label_no_homology", labels.N_INT),
        "nohom_n_sec": _n(rows, "label_no_homology", labels.N_SEC),
        "nohom_ambiguous": _n(rows, "label_no_homology", labels.AMBIGUOUS),
        "direct_p_ext": sum(in_direct_stratum(r, labels.P_EXT) for r in rows),
        "direct_n_int": sum(in_direct_stratum(r, labels.N_INT) for r in rows),
        "direct_n_sec": sum(in_direct_stratum(r, labels.N_SEC) for r in rows),
        "direct_ambiguous": sum(in_direct_stratum(r, labels.AMBIGUOUS) for r in rows),
        "ambiguous_htp_only": sum(
            r["label"] == labels.AMBIGUOUS and r["internal_evidence_htp_only"] == "yes"
            for r in rows
        ),
    }
    return {k: str(v) for k, v in counts.items()}


def write_tsv(path: str | Path, columns, rows) -> None:
    """Write a TSV. A path ending in .gz is gzip-compressed with a fixed mtime (byte-stable)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = buffer.getvalue().encode("utf-8")
    if path.suffix == ".gz":
        with (
            open(path, "wb") as raw,
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz,
        ):
            gz.write(data)
    else:
        path.write_bytes(data)


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
```

- [ ] **Step 4: Write `analysis/step1_compare/01_extract_go_truth.py`**

```python
#!/usr/bin/env python3
"""D1: extract the GO truth set (spec 3.1) into truth_set.tsv.gz and counts.tsv.

Inputs are read from --input-dir (default $STEP1_WORKDIR/downloads). Every input is checked
against manifest.tsv before it is read; a strict-mode SHA-256 mismatch stops the run.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import gaf
import go_obo
import manifest
import paths
import truth_table


def read_species(path: str | Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(
    species_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    input_dir: Path,
    out_dir: Path,
    obo_name: str = "go-basic.obo",
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    by_file = {r["file"]: r for r in manifest_rows}
    obo_path = input_dir / obo_name
    obo_sha = manifest.verify_against_manifest(obo_path, by_file[obo_name])
    ontology = go_obo.parse_obo(obo_path)
    truth_rows, count_rows, log = [], [], {"obo": obo_name, "obo_sha256": obo_sha, "sources": []}
    for sp in species_rows:
        gaf_path = input_dir / sp["gaf_file"]
        sha = manifest.verify_against_manifest(gaf_path, by_file[sp["gaf_file"]])
        filtered = gaf.filter_gaf(gaf_path, ontology, sp["taxon_filter"] or None)
        info = truth_table.SourceInfo(
            source_id=sp["source_id"],
            species=sp["species"],
            taxon_id=sp["taxon_id"],
            in_clade=sp["in_clade"],
            role=sp["role"],
            source_file=sp["gaf_file"],
            source_sha256=sha,
            source_date=gaf.header_value(gaf_path, "date-generated"),
            obo_sha256=obo_sha,
        )
        rows = truth_table.build_truth_rows(filtered, ontology, info)
        truth_rows.extend(rows)
        count_rows.append(truth_table.count_rows(rows, filtered, sp["source_id"]))
        log["sources"].append(
            {"source_id": sp["source_id"], "sha256": sha, "dropped": dict(filtered.dropped)}
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    truth_table.write_tsv(out_dir / "truth_set.tsv.gz", truth_table.TRUTH_COLUMNS, truth_rows)
    truth_table.write_tsv(out_dir / "counts.tsv", truth_table.COUNT_COLUMNS, count_rows)
    (out_dir / "extract_log.json").write_text(json.dumps(log, indent=2, sort_keys=True) + "\n")
    return truth_rows, count_rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--manifest", default=str(paths.STEP1_DIR / "manifest.tsv"))
    parser.add_argument("--input-dir", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument("--out-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sources", nargs="*", default=None, help="source_id values to run")
    args = parser.parse_args(argv)
    species_rows = read_species(args.species)
    if args.sources:
        species_rows = [r for r in species_rows if r["source_id"] in args.sources]
    input_dir = Path(args.input_dir) if args.input_dir else paths.downloads_dir()
    out_dir = Path(args.out_dir) if args.out_dir else paths.workdir()
    try:
        _, counts = run(species_rows, manifest.read_manifest(args.manifest), input_dir, out_dir)
    except manifest.DownloadError as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for c in counts:
        print(
            "\t".join(f"{k}={c[k]}" for k in ("source_id", "p_ext", "n_int", "n_sec", "ambiguous"))
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Run the fixture tests to verify they pass**

Run: `$PY -m pytest tests/step1_compare -q`
Expected: PASS (44 passed, 2 skipped: the real-data tests). If the Task 5 downloads are in `$STEP1_WORKDIR/downloads`, the real-data tests run as well: 46 passed.

- [ ] **Step 6: Run the real-data regression test and the real extraction**

```bash
STEP1_GO_DIR=/tmp/glyco_spec /usr/bin/python3.12 -m pytest \
  tests/step1_compare/test_extract.py -q -k real_files
export PROJ_ROOT=$PWD STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
cd analysis/step1_compare && /usr/bin/python3.12 01_extract_go_truth.py; cd -
```

Expected: both real-data tests PASS (about 20 s). The script prints one tab-separated line per source; the first two lines are `source_id=Scer_SGD p_ext=125 n_int=2645 n_sec=1548 ambiguous=48` and `source_id=Calb_CGD p_ext=259 n_int=1772 n_sec=961 ambiguous=100` (tabs shown as spaces). If `/tmp/glyco_spec` is gone, use the Task 5 download directory as `STEP1_GO_DIR`.

- [ ] **Step 7: Lint and commit**

```bash
pre-commit run --all-files
ruff check .
git add analysis/step1_compare/truth_table.py analysis/step1_compare/01_extract_go_truth.py \
        tests/step1_compare/test_extract.py
git commit -m "step1_compare: D1 extractor (truth_set.tsv.gz, counts.tsv) reproducing spec 3.3

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Sequence attachment (D1 step 2)

**Files:**
- Create: `analysis/step1_compare/seqhash.py`
- Create: `analysis/step1_compare/sequences.py`
- Create: `analysis/step1_compare/02_attach_sequences.py`
- Create: `tests/step1_compare/test_sequences.py`

**Interfaces:**
- Consumes: `gaf.open_text` (Task 3); `labels.*` (Task 4); `manifest.verify_against_manifest`, `manifest.read_manifest` (Task 5); `truth_table.read_tsv`, `truth_table.write_tsv` (Task 6); `species.tsv` `fasta_file`, `id_mapping` (Task 5).
- Produces: `seqhash.clean(sequence) -> str`; `seqhash.seq_sha256(sequence) -> str`. `sequences.MAPPINGS`; `sequences.MappingError`; `sequences.read_fasta(path) -> Iterator[tuple[str, str]]`; `sequences.normalize_id(value, mapping) -> str`; `sequences.fasta_key(header, mapping) -> str | None`; `sequences.gene_key(row, mapping) -> str`; `sequences.index_fasta(path, mapping) -> tuple[dict[str, tuple[str, str]], int]`; `sequences.attach(rows, index, mapping) -> (matched, unmatched)`. `02_attach_sequences.run(truth_rows, species_rows, manifest_rows, input_dir, out_dir) -> (seq_rows, unmatched_rows, counts)`; it raises `sequences.MappingError` (main exits 2, no `truth_sequences.tsv.gz`) when a source matches no gene or when a P-ext or ambiguous gene has no sequence; `unmatched_ids.tsv` and `sequence_counts.tsv` are written first so the cause can be read. Files: `truth_sequences.tsv.gz` (`source_id, gene_id, label, fasta_id, length, seq_sha256, sequence`), `unmatched_ids.tsv` (`source_id, gene_id, symbol, synonym1, label`), `sequence_counts.tsv` (`source_id, fasta_file, fasta_sha256, fasta_duplicate_records, genes, matched, unmatched, unmatched_p_ext, unmatched_ambiguous, unmatched_n_int, unmatched_n_sec`).

Mapping table (spec 3.3, corrected for CGD):

| `id_mapping` | GAF side | FASTA side | URL checked (HTTP 200, 2026-10-01) |
|---|---|---|---|
| `sgd` | column 2 `S000...` | header field `SGDID:S000...` | `https://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz` |
| `cgd` | first synonym in column 11 (`C1_00010W_A`) | first header token | `https://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_default_protein.fasta.gz` |
| `pombase` | column 2 `SPAC1002.01` | `SPAC1002.01.1:pep` without `.1:pep` | `https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/peptide.fa.gz` |
| `uniprot` | column 2 accession | `sp|ACC|` or `tr|ACC|` | `https://rest.uniprot.org/uniprotkb/search?query=proteome%3AUP000002530&format=fasta&size=500`; the same for UP000000560, UP000010091, UP000002149, UP000000561 (all five downloaded completely) |

Keys are upper-cased and stripped on both sides. For `uniprot` an isoform suffix (`-2`) is removed. `seqhash.clean` gives the same result as `surface_glyco.io.clean_sequence` for upper-case input (J to L, `*` removed); the test checks this when the package imports.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_sequences.py`:

```python
import gzip

import manifest
import pytest
import seqhash
import sequences
import truth_table
from conftest import load_script

SGD_FASTA = (
    ">YAL001C TFC3 SGDID:S000000001, Chr I from 151006-147594, Genome Release 64-5-1\n"
    "MVLTIYPDELVQ\nIVSDKI*\n"
    ">YAL002W VPS8 SGDID:S000000002, Chr I\nMEQNGLDHDSRS\n"
)
CGD_FASTA = ">C1_00010W_A\nMSTQKA\n>C1_00020C_A\nMKKLLV\n>C1_00020C_A\nMKKLLV\n"
POMBASE_FASTA = ">SPAC1002.01.1:pep mrx11|component\nMSKF\n>SPBC21H7.03c.1:pep x\nMAAA\n"
UNIPROT_FASTA = ">sp|P22146|GAS1_YEAST 1,3-beta-glucanosyltransferase\nMLFKS\n>tr|Q4WXC4|Q4WXC4_ASPFU CspA\nMKVA\n"


def _write(tmp_path, name, text, compress=False):
    path = tmp_path / name
    data = text.encode()
    path.write_bytes(gzip.compress(data) if compress else data)
    return path


def _row(gene_id, synonym1="", label="P-ext", source="X"):
    return {
        "source_id": source,
        "gene_id": gene_id,
        "symbol": gene_id,
        "synonym1": synonym1,
        "label": label,
    }


def test_sgd_mapping_reads_sgdid(tmp_path):
    index, dups = sequences.index_fasta(_write(tmp_path, "s.fa.gz", SGD_FASTA, True), "sgd")
    assert set(index) == {"S000000001", "S000000002"} and dups == 0
    assert index["S000000001"] == ("YAL001C", "MVLTIYPDELVQIVSDKI*")


def test_cgd_mapping_uses_first_synonym_and_merges_identical_duplicates(tmp_path):
    index, dups = sequences.index_fasta(_write(tmp_path, "c.fa", CGD_FASTA), "cgd")
    assert dups == 1
    matched, unmatched = sequences.attach(
        [_row("CAL0000173921", "C1_00010W_A"), _row("CAL0000000001", "snR5a")], index, "cgd"
    )
    assert [m["fasta_id"] for m in matched] == ["C1_00010W_A"]
    assert [u["gene_id"] for u in unmatched] == ["CAL0000000001"]


def test_duplicate_key_with_different_sequence_raises(tmp_path):
    bad = _write(tmp_path, "c.fa", ">C1_00010W_A\nMSTQKA\n>C1_00010W_A\nMSTQKV\n")
    with pytest.raises(sequences.MappingError, match="two different sequences"):
        sequences.index_fasta(bad, "cgd")


def test_pombase_suffix_and_case_insensitive_match(tmp_path):
    index, _ = sequences.index_fasta(_write(tmp_path, "p.fa", POMBASE_FASTA), "pombase")
    matched, unmatched = sequences.attach(
        [_row("SPAC1002.01"), _row("spbc21h7.03C"), _row("SPAC1002.02")], index, "pombase"
    )
    assert [m["fasta_id"] for m in matched] == ["SPAC1002.01.1:pep", "SPBC21H7.03c.1:pep"]
    assert [u["gene_id"] for u in unmatched] == ["SPAC1002.02"]


def test_uniprot_isoform_and_version_differences(tmp_path):
    index, _ = sequences.index_fasta(_write(tmp_path, "u.fa", UNIPROT_FASTA), "uniprot")
    matched, _ = sequences.attach([_row("P22146-2"), _row(" q4wxc4 ")], index, "uniprot")
    assert [m["fasta_id"] for m in matched] == ["sp|P22146|GAS1_YEAST", "tr|Q4WXC4|Q4WXC4_ASPFU"]


def test_clean_and_hash_ignore_case_whitespace_and_stop():
    assert seqhash.clean("mk j*\nA") == "MKLA"
    assert seqhash.seq_sha256("MKLA") == seqhash.seq_sha256("mkja*")


def test_clean_matches_package_clean_sequence():
    io_mod = pytest.importorskip("surface_glyco.io")
    for seq in ("MKJA*", "MSTQKA", "JJ*"):
        assert seqhash.clean(seq) == io_mod.clean_sequence(seq)


def test_attach_script_reports_unmatched_by_label(tmp_path):
    attach = load_script("02_attach_sequences")
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    _write(input_dir, "c.fa", CGD_FASTA)
    species = [{"source_id": "X", "fasta_file": "c.fa", "id_mapping": "cgd"}]
    manifest_rows = [{"file": "c.fa", "mode": "record", "sha256": ""}]
    truth = [
        _row("CAL1", "C1_00010W_A", "P-ext"),
        _row("CAL2", "tRNA1", "N-int"),
        _row("CAL3", "C1_00020C_A", "unlabelled"),
    ]
    seq_rows, unmatched, counts = attach.run(truth, species, manifest_rows, input_dir, tmp_path)
    assert counts[0]["matched"] == "2" and counts[0]["unmatched_n_int"] == "1"
    assert counts[0]["fasta_duplicate_records"] == "1"
    assert counts[0]["fasta_sha256"] == manifest.sha256_file(input_dir / "c.fa")
    assert seq_rows[0]["seq_sha256"] == seqhash.seq_sha256("MSTQKA")
    assert truth_table.read_tsv(tmp_path / "unmatched_ids.tsv")[0]["gene_id"] == "CAL2"


def test_attach_stops_when_a_source_matches_nothing(tmp_path):
    attach = load_script("02_attach_sequences")
    _write(tmp_path, "c.fa", CGD_FASTA)
    species = [{"source_id": "X", "fasta_file": "c.fa", "id_mapping": "cgd"}]
    manifest_rows = [{"file": "c.fa", "mode": "record", "sha256": ""}]
    truth = [_row("CAL9", "C9_99999W_A", "N-sec")]
    with pytest.raises(sequences.MappingError, match="no gene matched"):
        attach.run(truth, species, manifest_rows, tmp_path, tmp_path)
    assert not (tmp_path / "truth_sequences.tsv.gz").exists()


def test_attach_stops_when_a_positive_or_ambiguous_gene_is_unmatched(tmp_path):
    attach = load_script("02_attach_sequences")
    _write(tmp_path, "c.fa", CGD_FASTA)
    species = [{"source_id": "X", "fasta_file": "c.fa", "id_mapping": "cgd"}]
    manifest_rows = [{"file": "c.fa", "mode": "record", "sha256": ""}]
    truth = [_row("CAL1", "C1_00010W_A", "N-sec"), _row("CAL5", "C5_00000W_A", "ambiguous")]
    with pytest.raises(sequences.MappingError, match="ambiguous gene CAL5"):
        attach.run(truth, species, manifest_rows, tmp_path, tmp_path)
    assert truth_table.read_tsv(tmp_path / "unmatched_ids.tsv")[0]["gene_id"] == "CAL5"
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_sequences.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'seqhash'`.

- [ ] **Step 3: Write `analysis/step1_compare/seqhash.py`**

```python
"""Sequence cleaning and hashing shared by D1 step 2 and D10. Standard library only.

`clean` matches surface_glyco.io.clean_sequence (J -> L, '*' removed) and also upper-cases
and removes whitespace, so the same protein gives the same hash from any FASTA source.
"""

import hashlib


def clean(sequence: str) -> str:
    return "".join(sequence.split()).upper().replace("J", "L").replace("*", "")


def seq_sha256(sequence: str) -> str:
    return hashlib.sha256(clean(sequence).encode("ascii")).hexdigest()
```

- [ ] **Step 4: Write `analysis/step1_compare/sequences.py`**

```python
"""Map truth-table gene IDs to protein FASTA records (spec 3.3 table). Standard library only.

Mapping types (column `id_mapping` of species.tsv):
- sgd:     GAF column 2 `S000...`  <-> FASTA header field `SGDID:S000...`
- cgd:     GAF column 11 first synonym (`C1_00010W_A`) <-> FASTA first token. The CGD FASTA
           headers carry no CAL... ID (checked 2026-10-01).
- pombase: GAF column 2 `SPAC1002.01` <-> FASTA `SPAC1002.01.1:pep` (transcript suffix removed)
- uniprot: GAF column 2 accession <-> FASTA `sp|ACC|NAME` or `tr|ACC|NAME`
Keys are compared after `normalize_id`: whitespace stripped, upper-cased, and a UniProt
isoform suffix (`-2`) removed for the uniprot mapping.
"""

import re
from collections.abc import Iterator
from pathlib import Path

from gaf import open_text
from seqhash import clean

_SGD = re.compile(r"SGDID:(S\d+)")
_POMBASE = re.compile(r"^(\S+)\.\d+:pep\b")
_UNIPROT = re.compile(r"^(?:sp|tr)\|([^|]+)\|")
MAPPINGS = ("sgd", "cgd", "pombase", "uniprot")


class MappingError(ValueError):
    """Two FASTA records give the same key, or the mapping type is unknown."""


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    header, chunks = None, []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def normalize_id(value: str, mapping: str) -> str:
    key = value.strip().upper()
    if mapping == "uniprot":
        key = re.sub(r"-\d+$", "", key)
    return key


def fasta_key(header: str, mapping: str) -> str | None:
    if mapping == "sgd":
        match = _SGD.search(header)
        raw = match.group(1) if match else None
    elif mapping == "cgd":
        raw = header.split()[0] if header.split() else None
    elif mapping == "pombase":
        match = _POMBASE.match(header)
        raw = match.group(1) if match else None
    elif mapping == "uniprot":
        match = _UNIPROT.match(header)
        raw = match.group(1) if match else None
    else:
        raise MappingError(f"unknown id_mapping {mapping!r}")
    return normalize_id(raw, mapping) if raw else None


def gene_key(row: dict[str, str], mapping: str) -> str:
    value = row["synonym1"] if mapping == "cgd" else row["gene_id"]
    return normalize_id(value, mapping)


def index_fasta(path: str | Path, mapping: str) -> tuple[dict[str, tuple[str, str]], int]:
    """Return (key -> (fasta_id, sequence), number of identical duplicate records skipped).

    A key that occurs twice with the same cleaned sequence is kept once (the CGD Assembly 22
    file has 38 such records on 2026-10-01). A key that occurs twice with different sequences
    raises MappingError.
    """
    index: dict[str, tuple[str, str]] = {}
    duplicates = 0
    for header, seq in read_fasta(path):
        key = fasta_key(header, mapping)
        if key is None:
            continue
        if key in index:
            if clean(index[key][1]) != clean(seq):
                raise MappingError(f"{Path(path).name}: key {key} has two different sequences")
            duplicates += 1
            continue
        index[key] = (header.split()[0], seq)
    return index, duplicates


def attach(rows: list[dict[str, str]], index: dict[str, tuple[str, str]], mapping: str):
    """Return (matched, unmatched). matched rows gain fasta_id and sequence."""
    matched, unmatched = [], []
    for row in rows:
        hit = index.get(gene_key(row, mapping))
        if hit is None:
            unmatched.append(row)
        else:
            matched.append({**row, "fasta_id": hit[0], "sequence": hit[1]})
    return matched, unmatched
```

- [ ] **Step 5: Write `analysis/step1_compare/02_attach_sequences.py`**

```python
#!/usr/bin/env python3
"""D1 step 2: attach protein sequences to truth-table genes and report unmatched IDs.

Reads $STEP1_WORKDIR/truth_set.tsv.gz and the FASTA files named in species.tsv. Writes
truth_sequences.tsv.gz, unmatched_ids.tsv and sequence_counts.tsv to the work directory.
Stops (exit 2, no truth_sequences.tsv.gz) if a source matches no gene, or if a P-ext or
ambiguous gene has no sequence; unmatched_ids.tsv and sequence_counts.tsv are still written.
"""

import argparse
import sys
from pathlib import Path

import labels
import manifest
import paths
import seqhash
import sequences
import truth_table

SEQUENCE_COLUMNS = (
    "source_id",
    "gene_id",
    "label",
    "fasta_id",
    "length",
    "seq_sha256",
    "sequence",
)
UNMATCHED_COLUMNS = ("source_id", "gene_id", "symbol", "synonym1", "label")
SEQ_COUNT_COLUMNS = (
    "source_id",
    "fasta_file",
    "fasta_sha256",
    "fasta_duplicate_records",
    "genes",
    "matched",
    "unmatched",
    "unmatched_p_ext",
    "unmatched_ambiguous",
    "unmatched_n_int",
    "unmatched_n_sec",
)


def run(truth_rows, species_rows, manifest_rows, input_dir: Path, out_dir: Path):
    by_file = {r["file"]: r for r in manifest_rows}
    indexes: dict[str, tuple[dict, int, str]] = {}
    seq_rows, unmatched_rows, counts = [], [], []
    for sp in species_rows:
        rows = [r for r in truth_rows if r["source_id"] == sp["source_id"]]
        if not rows:
            continue
        fasta = sp["fasta_file"]
        if fasta not in indexes:
            path = input_dir / fasta
            sha = manifest.verify_against_manifest(path, by_file[fasta])
            index, duplicates = sequences.index_fasta(path, sp["id_mapping"])
            indexes[fasta] = (index, duplicates, sha)
        index, duplicates, sha = indexes[fasta]
        matched, unmatched = sequences.attach(rows, index, sp["id_mapping"])
        for m in matched:
            seq = seqhash.clean(m["sequence"])
            seq_rows.append(
                {
                    "source_id": m["source_id"],
                    "gene_id": m["gene_id"],
                    "label": m["label"],
                    "fasta_id": m["fasta_id"],
                    "length": str(len(seq)),
                    "seq_sha256": seqhash.seq_sha256(seq),
                    "sequence": seq,
                }
            )
        unmatched_rows.extend({c: u[c] for c in UNMATCHED_COLUMNS} for u in unmatched)
        counts.append(
            {
                "source_id": sp["source_id"],
                "fasta_file": fasta,
                "fasta_sha256": sha,
                "fasta_duplicate_records": str(duplicates),
                "genes": str(len(rows)),
                "matched": str(len(matched)),
                "unmatched": str(len(unmatched)),
                "unmatched_p_ext": str(sum(u["label"] == labels.P_EXT for u in unmatched)),
                "unmatched_ambiguous": str(sum(u["label"] == labels.AMBIGUOUS for u in unmatched)),
                "unmatched_n_int": str(sum(u["label"] == labels.N_INT for u in unmatched)),
                "unmatched_n_sec": str(sum(u["label"] == labels.N_SEC for u in unmatched)),
            }
        )
    truth_table.write_tsv(out_dir / "unmatched_ids.tsv", UNMATCHED_COLUMNS, unmatched_rows)
    truth_table.write_tsv(out_dir / "sequence_counts.tsv", SEQ_COUNT_COLUMNS, counts)
    problems = [f"{c['source_id']}: no gene matched" for c in counts if c["matched"] == "0"]
    problems += [
        f"{u['source_id']}: {u['label']} gene {u['gene_id']} has no sequence"
        for u in unmatched_rows
        if u["label"] in (labels.P_EXT, labels.AMBIGUOUS)
    ]
    if problems:
        raise sequences.MappingError("; ".join(problems[:10]))
    truth_table.write_tsv(out_dir / "truth_sequences.tsv.gz", SEQUENCE_COLUMNS, seq_rows)
    return seq_rows, unmatched_rows, counts


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--manifest", default=str(paths.STEP1_DIR / "manifest.tsv"))
    parser.add_argument("--input-dir", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    input_dir = Path(args.input_dir) if args.input_dir else paths.downloads_dir()
    species_rows = truth_table.read_tsv(args.species)
    truth_rows = truth_table.read_tsv(work / "truth_set.tsv.gz")
    try:
        _, _, counts = run(
            truth_rows, species_rows, manifest.read_manifest(args.manifest), input_dir, work
        )
    except (manifest.DownloadError, sequences.MappingError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for c in counts:
        print("\t".join(f"{k}={c[k]}" for k in ("source_id", "genes", "matched", "unmatched")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `$PY -m pytest tests/step1_compare/test_sequences.py -q`
Expected: PASS (9 passed, 1 skipped). The skip is `test_clean_matches_package_clean_sequence` when `surface_glyco` (Biopython) does not import under `/usr/bin/python3.12`. With `PYTHONPATH=$PWD/src /rhome/jstajich/.conda/envs/adhesionPred/bin/python -m pytest tests/step1_compare/test_sequences.py -q` all 10 pass.

- [ ] **Step 7: Run on the real files and record the unmatched counts**

```bash
export PROJ_ROOT=$PWD STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
cd analysis/step1_compare && /usr/bin/python3.12 02_attach_sequences.py; cd -
column -t "$STEP1_WORKDIR/sequence_counts.tsv"
```

Expected with the 2026-10-01 FASTA files (prototype run): Scer_SGD 6,052 matched / 4 unmatched (3 N-int: RAF1, REP1, REP2 of the 2-micron plasmid); Calb_CGD 6,060 / 253 (72 N-int; RNA genes typed `gene_product`; 38 duplicate FASTA records); Spom_PomBase 5,020 / 5 (1 N-int, 1 N-sec); Spom_SCHPO-mod 5,008 / 0; all UniProt-based sources 0 unmatched; `unmatched_p_ext` 0 everywhere. A later FASTA release can change these numbers. Step 8 puts the actual table in the commit message body.

- [ ] **Step 8: Lint and commit**

```bash
pre-commit run --all-files
ruff check .
git add analysis/step1_compare/seqhash.py analysis/step1_compare/sequences.py \
        analysis/step1_compare/02_attach_sequences.py tests/step1_compare/test_sequences.py
git commit -m "step1_compare: attach protein sequences, report unmatched IDs" \
  -m "$(column -t "$STEP1_WORKDIR/sequence_counts.tsv" | cut -c1-100)" \
  -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 8: D8 triage of P-ext genes with a plasma-membrane term

**Files:**
- Create: `analysis/step1_compare/d8_triage.py`
- Create: `analysis/step1_compare/curated_gpi.tsv`
- Create: `analysis/step1_compare/03_triage_pm.py`
- Create: `tests/step1_compare/test_d8.py`

**Interfaces:**
- Consumes: `truth_set.tsv.gz` columns `source_id, gene_id, symbol, label, pm_candidate, stratum` (Task 6); `truth_table.read_tsv`, `write_tsv`, `TRUTH_COLUMNS` (Task 6); `manifest.uniprot_pages` (Task 5); `species.tsv` `id_mapping`, `taxon_id` (Task 5).
- Produces: `d8_triage.CURATED_GPI_ECO`, `UNIPROT_FIELDS`, `P_GPI`, `PM_TM`, `PM_UNRESOLVED`; `d8_triage.UniprotError`; `d8_triage.UniprotEvidence` (`accession, reviewed, gpi_eco, tm_count, tm_eco, xrefs`, property `curated_gpi`); `parse_entry(entry: dict) -> UniprotEvidence` (raises `UniprotError` if `primaryAccession` or `entryType` is missing); `query_terms(rows, mapping) -> dict[str, str]`; `search_urls(terms, batch=50) -> list[str]`; `organism_gpi_url(taxon_id) -> str`; `entries_for_gene(gene_id, mapping, entries) -> list[UniprotEvidence]`; `classify_pm(entries, literature: bool) -> tuple[str, str]`. `03_triage_pm.triage_source(sp, rows, literature_ids, fetch) -> (triage, outside, counts)`; `03_triage_pm.apply_triage(truth_rows, triage_rows) -> list[dict]`; `main(argv=None, fetch=None) -> int` (exit 2 with `STOP:` and no output tables on an HTTP or URL error, a response without `results`, an entry without the expected fields, or PM candidates of which none matches a UniProt entry). Files: `d8_triage.tsv`, `d8_gpi_outside_pext.tsv`, `d8_counts.tsv`, `truth_set_triaged.tsv.gz`, raw JSON pages in `d8_uniprot/`.

UniProt REST query. Endpoint `https://rest.uniprot.org/uniprotkb/search`, `format=json`, `size=500`, `fields=accession,reviewed,ft_lipid,ft_transmem,xref_sgd,xref_cgd,xref_pombase`. Candidate genes are queried in batches of 50 terms joined by `OR`: `xref:sgd-<S000...>`, `xref:cgd-<CAL...>`, `xref:pombase-<SP...>` or `accession:<ACC>`. One more query per species lists every reviewed GPI entry: `(organism_id:<taxon>) AND (reviewed:true) AND (ft_lipid:GPI-anchor)`. The parser records, per entry, whether it is reviewed, the ECO codes of every `Lipidation` feature whose description starts with `GPI-anchor`, the number of `Transmembrane` features and their ECO codes.

Rules. P-gpi: a row in `curated_gpi.tsv`, or a reviewed entry with a GPI-anchor ECO code in `CURATED_GPI_ECO = {"ECO:0000269"}`. PM-TM: not P-gpi and at least one Transmembrane feature (the code is recorded; it is often ECO:0000255). pm-unresolved: neither; the spec does not define this case, so the gene stays P-ext with stratum `pm-unresolved` and goes to the owner. `d8_gpi_outside_pext.tsv` lists curated-GPI entries whose gene is not P-ext (for example an ambiguous gene); the plan does not relabel them. Because every ECO code is stored, a wider GPI rule later is a filter on `d8_triage.tsv`. **P-gpi is a list, not a scored stratum** (ruling R-B), until `curated_gpi.tsv` has literature rows. Filling that file is separate curation work, outside this plan.

A prototype run against UniProt 2026_03 on 2026-10-01 gave P-gpi / PM-TM / pm-unresolved: S288C 0 / 3 / 9, *C. albicans* 1 / 22 / 37, *S. pombe* 0 / 4 / 6, *A. fumigatus* 0 / 1 / 8, *A. nidulans* 1 / 2 / 3. The independent reviewer re-ran S288C and *S. pombe* only. Step 7 re-runs the query and records the counts of the day.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_d8.py`:

```python
import urllib.error
import urllib.parse

import d8_triage
import pytest
import truth_table
from conftest import load_script


def entry(acc, reviewed, gpi_eco=(), tm=(), sgd=None):
    features = []
    if gpi_eco:
        features.append(
            {
                "type": "Lipidation",
                "description": "GPI-anchor amidated asparagine",
                "evidences": [{"evidenceCode": c} for c in gpi_eco],
            }
        )
    for code in tm:
        features.append(
            {
                "type": "Transmembrane",
                "description": "Helical",
                "evidences": [{"evidenceCode": code}],
            }
        )
    xrefs = [{"database": "SGD", "id": sgd}] if sgd else []
    kind = "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)"
    return {
        "primaryAccession": acc,
        "entryType": kind,
        "features": features,
        "uniProtKBCrossReferences": xrefs,
    }


GAS1 = entry("P22146", True, gpi_eco=["ECO:0000269"], sgd="S000004924")
YPS1 = entry("P32329", True, gpi_eco=["ECO:0000255"], sgd="S000004110")
MSB2 = entry("P32334", True, tm=["ECO:0000255"], sgd="S000003246")
TREMBL_GPI = entry("A0A000", False, gpi_eco=["ECO:0000269"], sgd="S000000009")
BARE = entry("Q00001", True, sgd="S000000010")


def test_parse_entry_records_eco_and_review_state():
    ev = d8_triage.parse_entry(MSB2)
    assert (ev.accession, ev.reviewed, ev.tm_count, ev.tm_eco) == (
        "P32334",
        True,
        1,
        ["ECO:0000255"],
    )
    assert ev.xrefs == {"SGD": {"S000003246"}}
    assert d8_triage.parse_entry(GAS1).curated_gpi
    assert not d8_triage.parse_entry(YPS1).curated_gpi  # predicted GPI site
    assert not d8_triage.parse_entry(TREMBL_GPI).curated_gpi  # not reviewed


def test_classify_pm():
    p = d8_triage.parse_entry
    assert d8_triage.classify_pm([p(GAS1)], False)[0] == "P-gpi"
    assert d8_triage.classify_pm([p(MSB2)], False) == ("PM-TM", "P32334 1 TM ECO:0000255")
    assert d8_triage.classify_pm([p(YPS1)], False)[0] == "pm-unresolved"
    assert d8_triage.classify_pm([p(YPS1)], True)[0] == "P-gpi"  # literature row
    assert d8_triage.classify_pm([], False) == ("pm-unresolved", "no UniProt entry found")


def test_query_urls():
    terms = d8_triage.query_terms([{"gene_id": "S000003246"}], "sgd")
    assert terms == {"xref:sgd-S000003246": "S000003246"}
    assert d8_triage.query_terms([{"gene_id": "Q4WXC4"}], "uniprot") == {
        "accession:Q4WXC4": "Q4WXC4"
    }
    urls = d8_triage.search_urls([f"accession:P{i:05d}" for i in range(120)])
    assert len(urls) == 3
    query = urllib.parse.parse_qs(urllib.parse.urlparse(urls[0]).query)
    assert query["fields"] == [d8_triage.UNIPROT_FIELDS] and query["format"] == ["json"]
    assert "ft_lipid:GPI-anchor" in urllib.parse.unquote_plus(d8_triage.organism_gpi_url("559292"))


def test_triage_source_and_apply(tmp_path):
    triage = load_script("03_triage_pm")
    rows = [
        {
            "source_id": "Scer",
            "gene_id": "S000003246",
            "symbol": "MSB2",
            "label": "P-ext",
            "pm_candidate": "yes",
            "stratum": "wall",
        },
        {
            "source_id": "Scer",
            "gene_id": "S000004110",
            "symbol": "YPS1",
            "label": "P-ext",
            "pm_candidate": "yes",
            "stratum": "wall",
        },
        {
            "source_id": "Scer",
            "gene_id": "S000004924",
            "symbol": "GAS1",
            "label": "ambiguous",
            "pm_candidate": "no",
            "stratum": "ambiguous",
        },
    ]
    sp = {"source_id": "Scer", "id_mapping": "sgd", "taxon_id": "559292"}

    def fake_fetch(url, tag):
        if "organism_id" in urllib.parse.unquote_plus(url):
            return [GAS1, YPS1], "2026_03"
        return [MSB2, YPS1], "2026_03"

    t, outside, counts = triage.triage_source(sp, rows, set(), fake_fetch)
    assert {r["symbol"]: r["d8_class"] for r in t} == {"MSB2": "PM-TM", "YPS1": "pm-unresolved"}
    assert [o["symbol"] for o in outside] == ["GAS1"]
    assert (counts["pm_tm"], counts["pm_unresolved"], counts["p_gpi"]) == ("1", "1", "0")
    assert counts["uniprot_release"] == "2026_03"
    triaged = triage.apply_triage(rows, t)
    assert [r["stratum"] for r in triaged] == ["PM-TM", "pm-unresolved", "ambiguous"]
    assert triaged[0]["label"] == "P-ext"  # label kept; stratum carries the D8 class


def _work(tmp_path):
    truth = [
        {
            "source_id": "Scer",
            "gene_id": "S000003246",
            "symbol": "MSB2",
            "label": "P-ext",
            "pm_candidate": "yes",
            "stratum": "extracellular-only",
        },
    ]
    truth_table.write_tsv(tmp_path / "truth_set.tsv.gz", list(truth[0]), truth)
    species = [{"source_id": "Scer", "id_mapping": "sgd", "taxon_id": "559292"}]
    truth_table.write_tsv(tmp_path / "species.tsv", list(species[0]), species)
    (tmp_path / "curated_gpi.tsv").write_text("source_id\tgene_id\tsymbol\tpmid\tnote\n")
    return [
        "--species",
        str(tmp_path / "species.tsv"),
        "--work-dir",
        str(tmp_path),
        "--curated-gpi",
        str(tmp_path / "curated_gpi.tsv"),
    ]


def test_no_candidate_matched_stops_and_writes_nothing(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    assert triage.main(argv, fetch=lambda url, tag: ([], "2026_03")) == 2
    assert "none of 1 candidates matched" in capsys.readouterr().err
    assert not (tmp_path / "truth_set_triaged.tsv.gz").exists()
    assert not (tmp_path / "d8_triage.tsv").exists()


def test_http_error_stops(tmp_path):
    triage = load_script("03_triage_pm")

    def broken(url, tag):
        raise urllib.error.HTTPError(url, 500, "Server Error", {}, None)

    assert triage.main(_work(tmp_path), fetch=broken) == 2
    assert not (tmp_path / "truth_set_triaged.tsv.gz").exists()


def test_entry_without_expected_fields_raises():
    with pytest.raises(d8_triage.UniprotError, match="primaryAccession"):
        d8_triage.parse_entry({"features": []})
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_d8.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'd8_triage'`.

- [ ] **Step 3: Write `analysis/step1_compare/d8_triage.py`**

```python
"""D8: split P-ext genes with a non-IEA plasma-membrane term into P-gpi, PM-TM, pm-unresolved.

Standard library only. Evidence comes from UniProtKB REST JSON (fields accession, reviewed,
ft_lipid, ft_transmem, xref_sgd, xref_cgd, xref_pombase) and from curated_gpi.tsv (literature).

Rules (spec 2.2, Q2 and Q3):
- P-gpi: a literature row in curated_gpi.tsv, or a reviewed UniProt entry with a Lipidation
  feature whose description starts with "GPI-anchor" and whose evidence includes a code in
  CURATED_GPI_ECO. Predictor output never counts.
- PM-TM: not P-gpi, and at least one UniProt Transmembrane feature (any evidence code; the
  codes are recorded because they are often ECO:0000255, sequence analysis).
- pm-unresolved: neither. The spec does not define this case; it stays P-ext and is listed.
P-gpi is reported as a list, not a scored stratum, until curated_gpi.tsv has literature rows.
"""

import urllib.parse
from dataclasses import dataclass, field

CURATED_GPI_ECO = frozenset({"ECO:0000269"})  # experimental evidence used in manual assertion
UNIPROT_FIELDS = "accession,reviewed,ft_lipid,ft_transmem,xref_sgd,xref_cgd,xref_pombase"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
XREF_DB = {"sgd": "SGD", "cgd": "CGD", "pombase": "PomBase"}
P_GPI, PM_TM, PM_UNRESOLVED = "P-gpi", "PM-TM", "pm-unresolved"


class UniprotError(RuntimeError):
    """A UniProt response lacks expected fields, or no candidate gene matched any entry."""


@dataclass
class UniprotEvidence:
    accession: str
    reviewed: bool
    gpi_eco: list[str] = field(default_factory=list)
    tm_count: int = 0
    tm_eco: list[str] = field(default_factory=list)
    xrefs: dict[str, set[str]] = field(default_factory=dict)

    @property
    def curated_gpi(self) -> bool:
        return self.reviewed and bool(set(self.gpi_eco) & CURATED_GPI_ECO)


def _eco(feature: dict) -> list[str]:
    return [e["evidenceCode"] for e in feature.get("evidences", []) if "evidenceCode" in e]


def parse_entry(entry: dict) -> UniprotEvidence:
    missing = [k for k in ("primaryAccession", "entryType") if k not in entry]
    if missing:
        raise UniprotError(f"UniProt entry lacks {missing}: {str(entry)[:80]}")
    ev = UniprotEvidence(
        accession=entry["primaryAccession"],
        reviewed=entry.get("entryType", "").startswith("UniProtKB reviewed"),
    )
    gpi, tm = set(), set()
    for feature in entry.get("features", []):
        if feature.get("type") == "Lipidation" and feature.get("description", "").startswith(
            "GPI-anchor"
        ):
            gpi.update(_eco(feature))
        elif feature.get("type") == "Transmembrane":
            ev.tm_count += 1
            tm.update(_eco(feature))
    ev.gpi_eco, ev.tm_eco = sorted(gpi), sorted(tm)
    for xref in entry.get("uniProtKBCrossReferences", []):
        ev.xrefs.setdefault(xref["database"], set()).add(xref["id"])
    return ev


def query_terms(rows: list[dict[str, str]], mapping: str) -> dict[str, str]:
    """Return UniProt query term -> gene_id for each truth row."""
    if mapping == "uniprot":
        return {f"accession:{r['gene_id']}": r["gene_id"] for r in rows}
    return {f"xref:{mapping}-{r['gene_id']}": r["gene_id"] for r in rows}


def search_urls(terms: list[str], batch: int = 50) -> list[str]:
    urls = []
    for start in range(0, len(terms), batch):
        query = " OR ".join(f"({t})" for t in terms[start : start + batch])
        params = {"query": query, "fields": UNIPROT_FIELDS, "format": "json", "size": "500"}
        urls.append(f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}")
    return urls


def organism_gpi_url(taxon_id: str) -> str:
    query = f"(organism_id:{taxon_id}) AND (reviewed:true) AND (ft_lipid:GPI-anchor)"
    params = {"query": query, "fields": UNIPROT_FIELDS, "format": "json", "size": "500"}
    return f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"


def entries_for_gene(gene_id: str, mapping: str, entries: list[UniprotEvidence]):
    if mapping == "uniprot":
        return [e for e in entries if e.accession == gene_id]
    db = XREF_DB[mapping]
    return [e for e in entries if gene_id in e.xrefs.get(db, set())]


def classify_pm(entries: list[UniprotEvidence], literature: bool) -> tuple[str, str]:
    if literature:
        return P_GPI, "literature row in curated_gpi.tsv"
    for e in entries:
        if e.curated_gpi:
            return P_GPI, f"{e.accession} reviewed GPI-anchor {','.join(e.gpi_eco)}"
    with_tm = [e for e in entries if e.tm_count > 0]
    if with_tm:
        e = with_tm[0]
        return PM_TM, f"{e.accession} {e.tm_count} TM {','.join(e.tm_eco) or 'no ECO'}"
    if not entries:
        return PM_UNRESOLVED, "no UniProt entry found"
    return PM_UNRESOLVED, "no curated GPI evidence and no TM feature"
```

- [ ] **Step 4: Write `analysis/step1_compare/curated_gpi.tsv` (header only, tab-separated)**

```text
source_id	gene_id	symbol	pmid	note
```

The owner adds literature rows (one per gene, with a PMID). The file starts empty because no curated list exists in `data/curated/` for GPI anchoring of truth-set genes.

- [ ] **Step 5: Write `analysis/step1_compare/03_triage_pm.py`**

```python
#!/usr/bin/env python3
"""D8: triage P-ext genes with a non-IEA plasma-membrane term (P-gpi, PM-TM, pm-unresolved).

Reads $STEP1_WORKDIR/truth_set.tsv.gz, species.tsv and curated_gpi.tsv. Queries UniProtKB REST.
Writes d8_triage.tsv, d8_gpi_outside_pext.tsv, d8_counts.tsv, truth_set_triaged.tsv.gz and the
raw UniProt JSON pages (d8_uniprot/) to the work directory. Stops (exit 2, no output tables)
on an HTTP error, a response without the expected fields, or when PM candidates exist and none
of them matches a UniProt entry.
"""

import argparse
import json
import sys
import urllib.error
from pathlib import Path

import d8_triage
import manifest
import paths
import truth_table

TRIAGE_COLUMNS = (
    "source_id",
    "gene_id",
    "symbol",
    "d8_class",
    "d8_reason",
    "uniprot_accessions",
    "reviewed",
    "gpi_eco",
    "tm_count",
    "tm_eco",
)
OUTSIDE_COLUMNS = ("source_id", "gene_id", "symbol", "label", "uniprot_accession", "gpi_eco")
D8_COUNT_COLUMNS = (
    "source_id",
    "pm_candidates",
    "p_gpi",
    "pm_tm",
    "pm_unresolved",
    "organism_curated_gpi_entries",
    "curated_gpi_outside_pext",
    "uniprot_release",
)


def uniprot_fetch(url: str, raw_dir: Path, tag: str) -> tuple[list[dict], str]:
    """Return (entries, release). Every page is saved as JSON for provenance."""
    entries, release = [], ""
    raw_dir.mkdir(parents=True, exist_ok=True)
    for n, (body, headers) in enumerate(manifest.uniprot_pages(url)):
        (raw_dir / f"{tag}_{n:03d}.json").write_bytes(body)
        page = json.loads(body)
        if "results" not in page:
            raise d8_triage.UniprotError(f"{url}: response has no 'results' field")
        entries.extend(page["results"])
        release = release or headers.get("x-uniprot-release", "")
    return entries, release


def triage_source(sp, rows, literature_ids, fetch):
    mapping = sp["id_mapping"]
    candidates = [r for r in rows if r["pm_candidate"] == "yes"]
    entries, release = [], ""
    terms = d8_triage.query_terms(candidates, mapping)
    for i, url in enumerate(d8_triage.search_urls(sorted(terms))):
        got, release = fetch(url, f"{sp['source_id']}_candidates_{i}")
        entries.extend(d8_triage.parse_entry(e) for e in got)
    matched_any = any(
        d8_triage.entries_for_gene(r["gene_id"], mapping, entries) for r in candidates
    )
    if candidates and not matched_any:
        raise d8_triage.UniprotError(
            f"{sp['source_id']}: none of {len(candidates)} candidates matched a UniProt entry"
        )
    triage = []
    for r in candidates:
        mine = d8_triage.entries_for_gene(r["gene_id"], mapping, entries)
        d8_class, reason = d8_triage.classify_pm(
            mine, (sp["source_id"], r["gene_id"]) in literature_ids
        )
        triage.append(
            {
                "source_id": r["source_id"],
                "gene_id": r["gene_id"],
                "symbol": r["symbol"],
                "d8_class": d8_class,
                "d8_reason": reason,
                "uniprot_accessions": ",".join(e.accession for e in mine),
                "reviewed": ",".join("yes" if e.reviewed else "no" for e in mine),
                "gpi_eco": ";".join(",".join(e.gpi_eco) for e in mine),
                "tm_count": ";".join(str(e.tm_count) for e in mine),
                "tm_eco": ";".join(",".join(e.tm_eco) for e in mine),
            }
        )
    organism, release2 = fetch(d8_triage.organism_gpi_url(sp["taxon_id"]), f"{sp['source_id']}_gpi")
    curated = [e for e in (d8_triage.parse_entry(x) for x in organism) if e.curated_gpi]
    outside = []
    for r in rows:
        if r["label"] == "P-ext":
            continue
        for e in d8_triage.entries_for_gene(r["gene_id"], mapping, curated):
            outside.append(
                {
                    "source_id": r["source_id"],
                    "gene_id": r["gene_id"],
                    "symbol": r["symbol"],
                    "label": r["label"],
                    "uniprot_accession": e.accession,
                    "gpi_eco": ",".join(e.gpi_eco),
                }
            )
    counts = {
        "source_id": sp["source_id"],
        "pm_candidates": str(len(candidates)),
        "p_gpi": str(sum(t["d8_class"] == d8_triage.P_GPI for t in triage)),
        "pm_tm": str(sum(t["d8_class"] == d8_triage.PM_TM for t in triage)),
        "pm_unresolved": str(sum(t["d8_class"] == d8_triage.PM_UNRESOLVED for t in triage)),
        "organism_curated_gpi_entries": str(len(curated)),
        "curated_gpi_outside_pext": str(len(outside)),
        "uniprot_release": release or release2,
    }
    return triage, outside, counts


def apply_triage(truth_rows, triage_rows):
    """Return truth rows with d8_class and d8_reason; stratum becomes the D8 class."""
    by_gene = {(t["source_id"], t["gene_id"]): t for t in triage_rows}
    out = []
    for r in truth_rows:
        t = by_gene.get((r["source_id"], r["gene_id"]))
        new = {
            **r,
            "d8_class": t["d8_class"] if t else "",
            "d8_reason": t["d8_reason"] if t else "",
        }
        if t:
            new["stratum"] = t["d8_class"]
        out.append(new)
    return out


def main(argv=None, fetch=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--curated-gpi", default=str(paths.STEP1_DIR / "curated_gpi.tsv"))
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sources", nargs="*", default=None)
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    if fetch is None:

        def fetch(url, tag):
            return uniprot_fetch(url, work / "d8_uniprot", tag)

    truth_rows = truth_table.read_tsv(work / "truth_set.tsv.gz")
    literature_ids = {
        (r["source_id"], r["gene_id"]) for r in truth_table.read_tsv(args.curated_gpi)
    }
    triage, outside, counts = [], [], []
    for sp in truth_table.read_tsv(args.species):
        if args.sources and sp["source_id"] not in args.sources:
            continue
        rows = [r for r in truth_rows if r["source_id"] == sp["source_id"]]
        if not any(r["pm_candidate"] == "yes" for r in rows):
            continue
        try:
            t, o, c = triage_source(sp, rows, literature_ids, fetch)
        except (d8_triage.UniprotError, urllib.error.URLError, TimeoutError) as exc:
            print(f"STOP: {exc}", file=sys.stderr)
            return 2
        triage += t
        outside += o
        counts.append(c)
    truth_table.write_tsv(work / "d8_triage.tsv", TRIAGE_COLUMNS, triage)
    truth_table.write_tsv(work / "d8_gpi_outside_pext.tsv", OUTSIDE_COLUMNS, outside)
    truth_table.write_tsv(work / "d8_counts.tsv", D8_COUNT_COLUMNS, counts)
    triaged = apply_triage(truth_rows, triage)
    columns = (*truth_table.TRUTH_COLUMNS, "d8_class", "d8_reason")
    truth_table.write_tsv(work / "truth_set_triaged.tsv.gz", columns, triaged)
    for c in counts:
        print("\t".join(f"{k}={c[k]}" for k in D8_COUNT_COLUMNS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `$PY -m pytest tests/step1_compare/test_d8.py -q`
Expected: PASS (7 passed).

- [ ] **Step 7: Run the real query and record the counts**

```bash
export PROJ_ROOT=$PWD STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
cd analysis/step1_compare && /usr/bin/python3.12 03_triage_pm.py; cd -
column -t -s $'\t' "$STEP1_WORKDIR/d8_counts.tsv"
cat "$STEP1_WORKDIR/d8_gpi_outside_pext.tsv"
```

Expected: one line per source with PM candidates; `pm_candidates` equals `counts.tsv` (`Scer_SGD` 12, `Calb_CGD` 60, `Spom_PomBase` 10, `Afum_ASPFU` 9 on the pinned GAFs). Step 8 puts `d8_counts.tsv` (`p_gpi`, `pm_tm`, `pm_unresolved`, `curated_gpi_outside_pext`, UniProt release) in the commit message body. These are the counts that spec D8 asks for.

- [ ] **Step 8: Lint and commit**

```bash
pre-commit run --all-files
ruff check .
git add analysis/step1_compare/d8_triage.py analysis/step1_compare/curated_gpi.tsv \
        analysis/step1_compare/03_triage_pm.py tests/step1_compare/test_d8.py
git commit -m "step1_compare: D8 GPI/TM triage of P-ext genes with a plasma-membrane term" \
  -m "$(column -t "$STEP1_WORKDIR/d8_counts.tsv")" \
  -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 9: D10 keyword tier with test proteins removed

**Files:**
- Create: `analysis/step1_compare/keyword_tier.py`
- Create: `analysis/step1_compare/04_build_keyword_tier.py`
- Create: `tests/step1_compare/test_keyword_tier.py`

**Interfaces:**
- Consumes: `seqhash.seq_sha256`, `seqhash.clean`, `sequences.read_fasta`, `sequences.fasta_key` (Task 7); `truth_sequences.tsv.gz` (Task 7); `species.tsv` `role`, `id_mapping` (Task 5); `manifest.uniprot_pages`, `manifest.check_payload` (Task 5); `truth_table.read_tsv`, `write_tsv` (Task 6); `data/curated/surface/surface.tsv`; `data/curated/adhesins/eurotiomycetes_seeds.tsv`.
- Produces: `keyword_tier.SPOMBE_TAXA`, `TC_SOURCE`, `KEYWORD_COLUMNS`, `REMOVED_COLUMNS`; `read_seed_accessions(path) -> dict[str, str]`; `accession_fasta_urls(accessions, batch=100) -> list[str]`; `build_keyword_tier(surface_rows, seq_by_acc, seed_accessions, heldout_accessions, heldout_hashes) -> (kept, removed)`. `04_build_keyword_tier.heldout_sets(species_rows, seq_rows) -> (accessions, hashes)`; `main(argv=None) -> int`. Files: `keyword_tier.tsv.gz` (`accession, gene, genome, taxon_id, length, seq_sha256, tier`), `keyword_tier_removed.tsv` (`accession, gene, genome, taxon_id, reason, matched`), `keyword_sequences.fasta.gz`, `keyword_sequences.json` (`sha256`, `uniprot_release`, `accessions`, `fetched_on`). `04_build_keyword_tier.fetch_sequences(accessions, dest) -> dict`; `check_cached(seq_path) -> dict` (raises `manifest.DownloadError` if the sidecar is missing, has no release, or its SHA-256 differs).

`surface.tsv` has no sequence column. The script fetches the sequences of all `surface.tsv` accessions and all seed accessions from UniProt (100 accessions per paginated query) into `keyword_sequences.fasta.gz` in the work directory.

Test proteins are: the Eurotiomycetes literature rows (Q7), every *S. pombe* protein (taxon 284812 or 4896; Q6), and every **labelled** gene (`label` not `unlabelled`) of every source whose `role` is not `train`. The held-out set is wide on purpose: it includes `test_clade` (*A. fumigatus*, *A. nidulans*, H99), `alternate_file` and `undecided` sources, so a later choice on the Basidiomycota option never re-opens a leak. The removal log keeps the reason and the matched protein, so rows can be restored by a filter if the owner narrows the held-out set. The per-fold removal for S1 (spec 2.3, "the test fold in each CV split") belongs to the later training plan.

- [ ] **Step 1: Write the failing test**

`tests/step1_compare/test_keyword_tier.py`:

```python
import gzip
import json

import keyword_tier
import manifest
import pytest
import seqhash
from conftest import load_script

SEEDS_TSV = (
    "# comment line\n"
    "gene\tuniprot_query\tspecies\torder\tfamily\tclass\tevidence_level\tmoonlighting\tpmids"
    "\tevidence_summary\n"
    "SOWgp58\taccession:Q8NK60\tCoccidioides immitis\tOnygenales\tf\tadhesin\tE1\tno\t1\ts\n"
    "YPS3\t\tHistoplasma capsulatum\tOnygenales\tf\tadhesin\tE2\tno\t1\tunresolved\n"
)


def srow(acc, taxon="330879", source="uniprot_surface_kw", gene="g"):
    return {"accession": acc, "gene": gene, "genome": "G", "taxon_id": taxon, "source": source}


def test_read_seed_accessions_skips_comments_and_blank(tmp_path):
    path = tmp_path / "seeds.tsv"
    path.write_text(SEEDS_TSV)
    assert keyword_tier.read_seed_accessions(path) == {"Q8NK60": "SOWgp58"}


def test_tc_excludes_test_proteins():
    surface = [
        srow("Q8NK60", taxon="246410"),  # literature accession
        srow("O14000", taxon="284812"),  # S. pombe
        srow("Q4WAAA"),  # held-out accession (A. fumigatus truth gene)
        srow("Q4WBBB"),  # no sequence found
        srow("A0ACG8", taxon="443226"),  # same sequence as Q8NK60 under another accession
        srow("Q4WCCC"),  # same sequence as a held-out truth protein
        srow("Q4WDDD"),  # kept
        srow("P50142", source="curated_literature"),  # not T-c: ignored
    ]
    seqs = {
        "Q8NK60": "MKSOW",
        "O14000": "MPOM",
        "Q4WAAA": "MAAA",
        "A0ACG8": "mksow*",
        "Q4WCCC": "MCCC",
        "Q4WDDD": "MDDD",
        "P50142": "MHSP",
    }
    kept, removed = keyword_tier.build_keyword_tier(
        surface,
        seqs,
        {"Q8NK60": "SOWgp58"},
        {"Q4WAAA": "Afum_ASPFU:Q4WAAA"},
        {seqhash.seq_sha256("MCCC"): "Afum_ASPFU:Q4WZZZ"},
    )
    assert [r["accession"] for r in kept] == ["Q4WDDD"]
    assert {r["accession"]: r["reason"] for r in removed} == {
        "Q8NK60": "literature_accession",
        "O14000": "spombe_taxon",
        "Q4WAAA": "heldout_accession",
        "Q4WBBB": "no_sequence",
        "A0ACG8": "literature_hash",
        "Q4WCCC": "heldout_hash",
    }
    assert kept[0]["seq_sha256"] == seqhash.seq_sha256("MDDD") and kept[0]["tier"] == "T-c"


def test_heldout_sets_skip_training_and_unlabelled_rows():
    build = load_script("04_build_keyword_tier")
    species = [
        {"source_id": "Scer_SGD", "role": "train", "id_mapping": "sgd"},
        {"source_id": "Afum_ASPFU", "role": "test_clade", "id_mapping": "uniprot"},
        {"source_id": "Spom_PomBase", "role": "test_species", "id_mapping": "pombase"},
    ]
    seq_rows = [
        {"source_id": "Scer_SGD", "gene_id": "S1", "label": "P-ext", "seq_sha256": "h1"},
        {"source_id": "Afum_ASPFU", "gene_id": "Q4W1", "label": "N-sec", "seq_sha256": "h2"},
        {"source_id": "Afum_ASPFU", "gene_id": "Q4W2", "label": "unlabelled", "seq_sha256": "h3"},
        {"source_id": "Spom_PomBase", "gene_id": "SPAC1.01", "label": "P-ext", "seq_sha256": "h4"},
    ]
    accessions, hashes = build.heldout_sets(species, seq_rows)
    assert accessions == {"Q4W1": "Afum_ASPFU:Q4W1"}
    assert hashes == {"h2": "Afum_ASPFU:Q4W1", "h4": "Spom_PomBase:SPAC1.01"}


def test_cached_sequences_without_sidecar_are_refused(tmp_path, capsys):
    build = load_script("04_build_keyword_tier")
    (tmp_path / "keyword_sequences.fasta.gz").write_bytes(gzip.compress(b">sp|P1|X\nMK\n"))
    surface = tmp_path / "surface.tsv"
    surface.write_text("accession\tgene\tgenome\ttaxon_id\tsource\n")
    seeds = tmp_path / "seeds.tsv"
    seeds.write_text(SEEDS_TSV)
    argv = ["--surface", str(surface), "--seeds", str(seeds), "--work-dir", str(tmp_path)]
    assert build.main(argv) == 2
    assert "keyword_sequences.json missing" in capsys.readouterr().err


def test_check_cached_accepts_matching_sidecar_and_rejects_changed_file(tmp_path):
    build = load_script("04_build_keyword_tier")
    seq_path = tmp_path / "keyword_sequences.fasta.gz"
    seq_path.write_bytes(gzip.compress(b">sp|P1|X\nMK\n"))
    meta = {"sha256": manifest.sha256_file(seq_path), "uniprot_release": "2026_03"}
    (tmp_path / "keyword_sequences.json").write_text(json.dumps(meta))
    assert build.check_cached(seq_path)["uniprot_release"] == "2026_03"
    seq_path.write_bytes(gzip.compress(b">sp|P1|X\nMA\n"))
    with pytest.raises(manifest.DownloadError, match="does not match"):
        build.check_cached(seq_path)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `$PY -m pytest tests/step1_compare/test_keyword_tier.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'keyword_tier'`.

- [ ] **Step 3: Write `analysis/step1_compare/keyword_tier.py`**

```python
"""D10: the T-c keyword tier for V-kw with test proteins removed. Standard library only.

A T-c row is removed when it is a test protein (spec 2.3, Q4 and Q7). Reasons, first match wins:
1. literature_accession   accession is a row of eurotiomycetes_seeds.tsv
2. spombe_taxon           taxon is S. pombe (284812 or 4896); S. pombe is test-only (Q6)
3. heldout_accession      accession is a gene of a non-training truth source
4. no_sequence            no sequence was found, so the hash rule cannot be checked
5. literature_hash        exact cleaned sequence equals a literature row's sequence
6. heldout_hash           exact cleaned sequence equals a non-training truth protein
"""

import csv
import urllib.parse

import seqhash

SPOMBE_TAXA = frozenset({"284812", "4896"})
TC_SOURCE = "uniprot_surface_kw"
KEYWORD_COLUMNS = ("accession", "gene", "genome", "taxon_id", "length", "seq_sha256", "tier")
REMOVED_COLUMNS = ("accession", "gene", "genome", "taxon_id", "reason", "matched")


def read_seed_accessions(path) -> dict[str, str]:
    """Return accession -> gene for rows of eurotiomycetes_seeds.tsv that have an accession."""
    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    out = {}
    for row in csv.DictReader(lines, delimiter="\t"):
        query = (row.get("uniprot_query") or "").strip()
        if query.startswith("accession:"):
            out[query.split(":", 1)[1].strip()] = row["gene"]
    return out


def accession_fasta_urls(accessions: list[str], batch: int = 100) -> list[str]:
    urls = []
    for start in range(0, len(accessions), batch):
        query = " OR ".join(f"(accession:{a})" for a in accessions[start : start + batch])
        params = {"query": query, "format": "fasta", "size": "500"}
        urls.append(f"https://rest.uniprot.org/uniprotkb/search?{urllib.parse.urlencode(params)}")
    return urls


def build_keyword_tier(
    surface_rows: list[dict[str, str]],
    seq_by_acc: dict[str, str],
    seed_accessions: dict[str, str],
    heldout_accessions: dict[str, str],
    heldout_hashes: dict[str, str],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Return (kept T-c rows, removal log). heldout_* values name the matched test protein."""
    seed_hashes = {seqhash.seq_sha256(seq_by_acc[a]): a for a in seed_accessions if a in seq_by_acc}
    kept, removed = [], []
    for row in surface_rows:
        if row["source"] != TC_SOURCE:
            continue
        acc = row["accession"]
        seq = seq_by_acc.get(acc)
        digest = seqhash.seq_sha256(seq) if seq else ""
        reason, matched = "", ""
        if acc in seed_accessions:
            reason, matched = "literature_accession", seed_accessions[acc]
        elif row["taxon_id"] in SPOMBE_TAXA:
            reason, matched = "spombe_taxon", row["taxon_id"]
        elif acc in heldout_accessions:
            reason, matched = "heldout_accession", heldout_accessions[acc]
        elif not seq:
            reason = "no_sequence"
        elif digest in seed_hashes:
            reason, matched = "literature_hash", seed_hashes[digest]
        elif digest in heldout_hashes:
            reason, matched = "heldout_hash", heldout_hashes[digest]
        base = {k: row[k] for k in ("accession", "gene", "genome", "taxon_id")}
        if reason:
            removed.append({**base, "reason": reason, "matched": matched})
        else:
            length = str(len(seqhash.clean(seq)))
            kept.append({**base, "length": length, "seq_sha256": digest, "tier": "T-c"})
    return kept, removed
```

- [ ] **Step 4: Write `analysis/step1_compare/04_build_keyword_tier.py`**

```python
#!/usr/bin/env python3
"""D10: build keyword_tier.tsv.gz (T-c rows for V-kw) and keyword_tier_removed.tsv.

Reads data/curated/surface/surface.tsv, data/curated/adhesins/eurotiomycetes_seeds.tsv,
species.tsv and $STEP1_WORKDIR/truth_sequences.tsv.gz. Sequences for surface.tsv and the seed
accessions come from $STEP1_WORKDIR/keyword_sequences.fasta.gz; --fetch downloads that file
from UniProtKB REST if it is missing and records its SHA-256 and UniProt release in
keyword_sequences.json. A cached sequence file without a matching sidecar is refused.
"""

import argparse
import datetime
import gzip
import json
import sys
from pathlib import Path

import keyword_tier
import manifest
import paths
import sequences
import truth_table


def fetch_sequences(accessions: list[str], dest: Path) -> dict[str, str]:
    """Download the sequences and write a JSON sidecar with SHA-256 and UniProt release."""
    part = dest.with_name(dest.name + ".part")
    release = ""
    with open(part, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        for url in keyword_tier.accession_fasta_urls(accessions):
            for body, headers in manifest.uniprot_pages(url):
                gz.write(body)
                release = release or headers.get("x-uniprot-release", "")
    manifest.check_payload(part)
    part.replace(dest)
    meta = {
        "sha256": manifest.sha256_file(dest),
        "uniprot_release": release,
        "accessions": str(len(accessions)),
        "fetched_on": datetime.datetime.now(datetime.UTC).date().isoformat(),
    }
    sidecar(dest).write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return meta


def sidecar(seq_path: Path) -> Path:
    return seq_path.with_name("keyword_sequences.json")


def check_cached(seq_path: Path) -> dict[str, str]:
    """Return the sidecar of a cached sequence file. Refuse a file without a matching sidecar."""
    meta_path = sidecar(seq_path)
    if not meta_path.exists():
        raise manifest.DownloadError(f"{meta_path} missing; delete {seq_path} and rerun --fetch")
    meta = json.loads(meta_path.read_text())
    if not meta.get("uniprot_release") or meta.get("sha256") != manifest.sha256_file(seq_path):
        raise manifest.DownloadError(
            f"{seq_path}: SHA-256 or UniProt release does not match {meta_path.name}"
        )
    return meta


def heldout_sets(species_rows, seq_rows):
    heldout_sources = {r["source_id"]: r for r in species_rows if r["role"] != "train"}
    accessions, hashes = {}, {}
    for r in seq_rows:
        sp = heldout_sources.get(r["source_id"])
        if sp is None or r["label"] == "unlabelled":
            continue
        label = f"{r['source_id']}:{r['gene_id']}"
        if sp["id_mapping"] == "uniprot":
            accessions[r["gene_id"]] = label
        hashes.setdefault(r["seq_sha256"], label)
    return accessions, hashes


def main(argv=None) -> int:
    root = paths.repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", default=str(root / "data/curated/surface/surface.tsv"))
    parser.add_argument(
        "--seeds", default=str(root / "data/curated/adhesins/eurotiomycetes_seeds.tsv")
    )
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()

    surface_rows = truth_table.read_tsv(args.surface)
    seeds = keyword_tier.read_seed_accessions(args.seeds)
    seq_path = work / "keyword_sequences.fasta.gz"
    try:
        if seq_path.exists():
            meta = check_cached(seq_path)
        elif args.fetch:
            wanted = sorted({r["accession"] for r in surface_rows} | set(seeds))
            meta = fetch_sequences(wanted, seq_path)
        else:
            print(f"STOP: {seq_path} missing; rerun with --fetch", file=sys.stderr)
            return 2
    except manifest.DownloadError as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    seq_by_acc = {
        sequences.fasta_key(h, "uniprot"): s
        for h, s in sequences.read_fasta(seq_path)
        if sequences.fasta_key(h, "uniprot")
    }
    accessions, hashes = heldout_sets(
        truth_table.read_tsv(args.species), truth_table.read_tsv(work / "truth_sequences.tsv.gz")
    )
    kept, removed = keyword_tier.build_keyword_tier(
        surface_rows, seq_by_acc, seeds, accessions, hashes
    )
    truth_table.write_tsv(work / "keyword_tier.tsv.gz", keyword_tier.KEYWORD_COLUMNS, kept)
    truth_table.write_tsv(work / "keyword_tier_removed.tsv", keyword_tier.REMOVED_COLUMNS, removed)
    reasons: dict[str, int] = {}
    for r in removed:
        reasons[r["reason"]] = reasons.get(r["reason"], 0) + 1
    print(
        f"kept={len(kept)} removed={len(removed)} {reasons} "
        f"sequences_sha256={meta['sha256'][:12]} uniprot_release={meta['uniprot_release']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `$PY -m pytest tests/step1_compare/test_keyword_tier.py -q`
Expected: PASS (3 passed).

- [ ] **Step 6: Run on the real data and record the removal counts**

```bash
export PROJ_ROOT=$PWD STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
cd analysis/step1_compare && /usr/bin/python3.12 04_build_keyword_tier.py --fetch \
  | tee "$STEP1_WORKDIR/keyword_tier_summary.txt"; cd -
```

Expected (prototype run, 2026-10-01, UniProt 2026_03, about 1 minute): `kept=3092 removed=448 {'literature_accession': 11, 'heldout_accession': 220, 'literature_hash': 1, 'spombe_taxon': 216} sequences_sha256=091887296077 uniprot_release=2026_03`. The role change of *A. nidulans* and H99 to `test_clade` does not change these numbers, because every non-training source was already held out. The sequence hash changes with each fetch. A cached sequence file without its `keyword_sequences.json` sidecar is refused. The `literature_hash` row is A0ACG8D6C9 (*C. posadasii* C735), whose sequence equals SOWgp58 Q8NK60. Step 7 puts the actual line in the commit message body.

- [ ] **Step 7: Lint and commit**

```bash
pre-commit run --all-files
ruff check .
git add analysis/step1_compare/keyword_tier.py analysis/step1_compare/04_build_keyword_tier.py \
        tests/step1_compare/test_keyword_tier.py
git commit -m "step1_compare: D10 keyword tier with test proteins removed by accession and hash" \
  -m "$(cat "$STEP1_WORKDIR/keyword_tier_summary.txt")" \
  -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 10: Verification and re-run documentation

**Files:**
- Create: `analysis/step1_compare/README.md`

**Interfaces:**
- Consumes: every script and test from Tasks 1 to 9.
- Produces: `analysis/step1_compare/README.md`.

- [ ] **Step 1: Write `analysis/step1_compare/README.md`**

````markdown
# step1_compare: the step 1 GO truth set (Phase A)

Implements deliverables D1, D8, D10 and the extraction half of D9 of
`docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md`.
Plan: `docs/superpowers/plans/2026-10-01-step1-truth-set.md`.

All modules use the Python standard library only. Data never go into the repository.
Outputs go to `$STEP1_WORKDIR` (default: `<workdir in config/site.yaml>/step1_compare`).

## Re-run

```bash
export PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred
export STEP1_WORKDIR=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare
PY=/usr/bin/python3.12
cd "$PROJ_ROOT/analysis/step1_compare"
$PY 00_fetch_inputs.py            # downloads, checks SHA-256 against manifest.tsv
$PY 01_extract_go_truth.py        # truth_set.tsv.gz, counts.tsv, extract_log.json
$PY 02_attach_sequences.py        # truth_sequences.tsv.gz, unmatched_ids.tsv, sequence_counts.tsv
$PY 03_triage_pm.py               # d8_triage.tsv, d8_counts.tsv, truth_set_triaged.tsv.gz
$PY 04_build_keyword_tier.py --fetch   # keyword_tier.tsv.gz, keyword_tier_removed.tsv
```

The downloads total about 83 MB (measured 2026-10-01). No SLURM job is needed. Run on a login or
interactive node. If a job is used later, write to `${SCRATCH:?}` and copy results to /bigdata.

A strict-mode SHA-256 mismatch stops `00_fetch_inputs.py` and `01_extract_go_truth.py` (exit 2).
Accept a new upstream file only on purpose: `00_fetch_inputs.py --update-manifest`, then rerun
the regression test and commit the new `manifest.tsv`.

## Tests

```bash
cd "$PROJ_ROOT"
/usr/bin/python3.12 -m pytest tests/step1_compare -q
```

The real-data tests (`test_reproduces_spec_counts_on_real_files`, `test_direct_evidence_counts_on_real_files`) run only when the GO files are present in
`$STEP1_WORKDIR/downloads` or in `$STEP1_GO_DIR`.

## Truth subsets as filters on `truth_set.tsv.gz`

| Subset | Filter |
|---|---|
| Direct-evidence truth (headline metrics and gates) | `label == X and homology_only == "no"`; counts in `counts.tsv` `direct_*` |
| All non-IEA truth (reported beside it) | `label == X` |
| Ambiguous with only high-throughput internal evidence (R-A sub-stratum) | `label == "ambiguous" and internal_evidence_htp_only == "yes"` |
| Not used | `label_no_homology == X`: recomputing without homology codes moves genes whose only internal term is IBA or ISS into P-ext (83 in *C. albicans*), against Q9 |

P-gpi is a list (`d8_triage.tsv`), not a scored stratum, until `curated_gpi.tsv` has literature rows.

## Open owner question and its filter

| Question | Filter |
|---|---|
| Basidiomycota truth option | rows with `in_clade == "Basidiomycota"`; `Cneo_H99_GOA` is the reference (`role == "test_clade"`), `Umay_MYCMD` is `undecided` |
````

- [ ] **Step 2: Run the full test folder (fixtures only)**

Run: `env -u STEP1_WORKDIR -u STEP1_GO_DIR /usr/bin/python3.12 -m pytest tests/step1_compare -q`
Expected: PASS (65 passed, 3 skipped) when `$STEP1_WORKDIR/downloads` under the `config/site.yaml` work directory does not hold the GO files. If Task 5 step 8 filled that directory, the two real-data tests run too: 67 passed, 1 skipped.

- [ ] **Step 3: Run the full test folder with the real files**

Run: `STEP1_GO_DIR=/tmp/glyco_spec /usr/bin/python3.12 -m pytest tests/step1_compare -q`
Expected: PASS (67 passed, 1 skipped; about 20 s).

- [ ] **Step 4: Run the package tests and lint**

```bash
PYTHONPATH=$PWD/src /rhome/jstajich/.conda/envs/adhesionPred/bin/python -m pytest tests/surface_glyco -q
pre-commit run --all-files
ruff check .
ruff format --check .
```

Expected: the package tests pass unchanged (this plan does not touch `src/`); pre-commit and ruff report no change.

- [ ] **Step 5: Confirm no data entered the repository**

Run: `git status --short && git ls-files analysis/step1_compare tests/step1_compare`
Expected: only the files in the File structure table; no `.gaf`, `.gz`, `.fasta` or `.json` file.

- [ ] **Step 6: Commit**

```bash
git add analysis/step1_compare/README.md
git commit -m "step1_compare: README with re-run steps and truth subsets as filters

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

## Self-review

**1. Spec coverage.**

| Spec item | Task |
|---|---|
| 3.1 step 1 (obo, `is_a` + `part_of`) | 2 |
| 3.1 steps 2 to 4 (primary DB, object type, CGD taxon, aspect C, `NOT`, obsolete) | 3 |
| 3.1 step 5 and 2.2 labels P-ext / ambiguous / N-int / N-sec, `subset`, IEA fraction | 4, 6 |
| 3.3 counts reproduced exactly | 6 (`test_reproduces_spec_counts_on_real_files`) |
| 3.2 file pinning, date and hash per file | 5 (`manifest.tsv`), 6 (`source_sha256`, `source_date`) |
| 3.3 ID-to-sequence mapping, unmatched counts (section 4 step 2) | 7 |
| 2.2 P-gpi, PM-TM, D8 | 8 |
| 2.3 and D10 (removal by accession and hash, log) | 9 |
| D9 extraction half: *S. pombe* and the Basidiomycota GO sources | 6 (all ten sources in `species.tsv`) |
| Section 7 tests `test_gaf_extract_golden`, `test_label_rules_truth_table`, `test_truth_set_has_no_iea`, `test_tc_excludes_test_proteins` | 3, 4, 4, 9 |
| 3.3 direct-evidence truth (intersection), R-A sub-stratum, R-B list | 4, 6 (`test_direct_evidence_counts_on_real_files`), 8 |
| Section 12 decisions, rulings R-A to R-C, and the remaining open question as filters | Global Constraints table; columns in Task 6 |

Not covered here, by design (later plans): section 4 steps 4, 6 to 8 (dedupe, MMseqs2 clusters and the `cluster` column, class balance, splits); `test_cterm_window_takes_last_1022`, `test_rule_truth_table`, cluster leakage tests, `test_golden_scores`; D2 to D7; the D9 Basidiomycota curation, which waits for the owner's Q8 option.

**2. Placeholder scan.** Searched for TBD, TODO, "similar to", "add error handling" and "fill in": none. Angle brackets appear only where the text describes a query or path pattern, not as work left to do. Every code step has complete code. Commit bodies that carry real counts take them from the files that the preceding step writes.

**3. Type consistency.** Checked names across tasks: `filter_gaf`, `FilteredGaf.cc_rows`, `all_aspect_triples`, `dropped`; `build_truth_rows`, `count_rows`, `write_tsv`, `read_tsv`, `TRUTH_COLUMNS`, `COUNT_COLUMNS`; `verify_against_manifest`, `uniprot_pages`, `check_payload`; `index_fasta` returns `(index, duplicates)` and `02_attach_sequences.py` unpacks both; `heldout_sets` reads `label` from `truth_sequences.tsv.gz`, which Task 7 writes. After the review fixes of 2026-10-01 the code blocks were extracted from this plan and applied in task order in a clean directory (cumulative: 4, 7, 12, 31, 40, 44 + 2 skipped, 53 + 3, 60 + 3, 65 + 3, 65 + 3); after each task the cumulative tests passed and `ruff check .` and `ruff format --check .` (0.3.5, repository `pyproject.toml`) reported no finding.

**4. Review Focus.** Each of the five items has its test in the owning task (Tasks 3, 4, 5, 6, 7, 9). Further inputs checked by tests: an empty GAF, a row with fewer than 13 columns, a gzip file without a `.gz` suffix, a cached file (no second download), UniProt `Link` paging, an IEA-only gene, a membrane-only internal gene, an unreviewed entry with experimental GPI evidence (does not count), a missing GAF, an HTTP 403, a broken connection during UniProt paging, a source with no matched gene, an unmatched ambiguous gene, a UniProt answer with no entries, an entry without `primaryAccession`, and a cached keyword-sequence file without its sidecar.
