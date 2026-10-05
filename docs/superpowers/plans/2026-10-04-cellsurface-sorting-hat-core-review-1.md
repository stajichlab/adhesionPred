# Review 1 of the core-engine plan, with dispositions

*Reviewer: an independent Opus subagent, 2026-10-04. It replayed the plan, walked the spec,
reviewed the code and ran adversarial inputs in a scratch directory. It changed no repository file.
Verdict on the first draft: **needs rework** (1 blocker, 8 major, 10 minor). The plan was then
revised and checked by me, not by a second reviewer: I extracted all 18 file blocks of the revised
plan into a fresh directory (117 tests pass), ran ruff 0.3.5 and 0.15 (both clean), and ran six
mutation checks (each makes the named tests fail). A second independent review of the revised plan
has not run.*

Codes: **A** applied. **P** left to Plan 2 (stated in the plan). **N** noted, no change.

| # | Sev | Finding (short) | Disposition |
|---|---|---|---|
| B1 | blocker | Task 7 edit to `pyproject.toml` breaks `tests/surface_glyco/test_package.py`; CI goes red. | A. Task 8 step 2 edits that test; step 3 runs it. |
| M1 | major | `I001` lint failures in Tasks 1 to 6 because `known-first-party` is set only in Task 7. | A. Moved to Task 0. Checked with ruff 0.3.5. |
| M2 | major | `pip` and `python` on the login node are Python 3.9. | A. Task 8 step 7 uses `/usr/bin/python3.12 -m venv`. |
| M3 | major | A `not_assessable` result got the status of the true or false inputs. | A. An unknown result has no deciding module: status `unvalidated`, empty basis. Test added; mutation check added. |
| M4 | major | `write_atomic` and `ModuleCache.update` are not safe for parallel runs (1,026 reader errors in 6,000 reads; a lost-update race left a table unreadable). | A. `fcntl` lock (shared for reads, exclusive for updates), unique temporary names, a damaged table is a cache miss. Parallel test with 4 processes added. Not tested on a network file system. |
| M5 | major | Bad module tables accepted silently (empty or wrong `call`, duplicate rows, no ID match, JSON-only modules, missing `state` column). | A. Tables are validated; bad values become `bad_value` and are counted; no ID match makes the module `error`; JSON-only run records are listed. Tests added. |
| M6 | major | The report misses thresholds, per-module non-ok counts, trailing-`*` count, invalid IDs, absent modules, module warnings, taxon-not-checked sentence. | A. Added, with tests. |
| M7 | major | Evidence outputs (`evidence.tsv.gz`, `allergen_homolog_hit`, antigen axes) are missing from the plan and from Plan 2. | A. `evidence:` list in `categories.yaml`, `evidence.tsv.gz`, `allergen_homolog_hit` call. Kind K match thresholds are named as Plan 2 defaults. |
| M8 | major | Review Focus 4 (identical sequences) is not pinned by code; `sha256` is not used or written. | A. `proteins.tsv.gz` carries the sha256; the driver flags identical sequences whose module rows differ. Cache expansion per sha256 stays in Plan 2 (stated). |
| m1 | minor | Task 1 step 2 error text differs. | A. |
| m2 | minor | Two mutation expectations were wrong (the `continue` mutation survives). | A. Table rewritten from measured results; six mutations. |
| m3 | minor | `main()` catches `KeyError` and `ValueError` and hides program bugs. | A. Only project errors, `OSError` and `CalledProcessError` are caught. Test that a `KeyError` from the engine is not hidden. |
| m4 | minor | Config validation gaps (`other` with an ungated surface call; `{step1}` in an ungated call). | A. Checks and tests added. |
| m5 | minor | Taxon-map parse problems. | A. Errors name the file and line; IDs not in the FASTA are reported. |
| m6 | minor | `merged.dmp` not read; status entries for unknown taxa never match silently. | N. Stated as known gaps 4 in the plan. |
| m7 | minor | Sequence lines before the first header dropped; module tables only `.tsv.gz`. | A for the FASTA (error). N for tables: the contract says always gzip. |
| m8 | minor | Spec deviations not stated (`other_basis` column name, config hash location, golden test form). | A. Clarifications 4 to 7 in the plan. |
| m9 | minor | The spec exists only on the local branch; the Step 4 replacement text spans a line break. | A. The spec edits were made directly (no plan step). Task 0 says to branch from the PR tip until the PR is merged. |
| m10 | minor | Tasks 5 and 6 are large; `outputs.py` has no unit tests; trailers; CHANGELOG and AGENTS.md. | A. Task 6 split into outputs (with its own tests) and the driver; CHANGELOG and AGENTS.md entries added. Trailers follow the session's attribution rule. |

Not checked by the reviewer: `pre-commit run`, the full `tests/surface_glyco` suite (needs torch),
`os.replace` behaviour on /bigdata, README rendering.
