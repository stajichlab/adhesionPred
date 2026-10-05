# Review 2 of the core-engine plan, with dispositions

*Reviewer: an independent Fable subagent, 2026-10-05 (a different model from the author and from
reviewer 1). It replayed the revised plan, re-ran the six mutation checks, applied the Task 8
edits to copies of the real repository files, and ran adversarial inputs in a scratch directory.
Verdict: **needs rework (small)**: 0 blocker, 2 major, 8 minor. All review 1 findings were fixed.
Revision 3 of the plan applies the findings below. I checked revision 3 myself (18 file blocks
extracted, 143 tests pass, ruff 0.3.5 clean, nine mutation checks). A third review has not run.*

| # | Sev | Finding (short) | Disposition |
|---|---|---|---|
| 1 | major | Malformed input files crash with a traceback and exit code 1: a module row with an extra field, a truncated gzip, invalid JSON, a JSON run record that is not an object, a non-UTF-8 FASTA header, a non-integer parent in `nodes.dmp`. A row with a bad state could also leave a partial output directory. | A. Module tables are read with `csv.reader` and a field-count check; gzip, JSON, UTF-8 and `nodes.dmp` errors become `InputError`/`FastaError`/`TaxonError` that name the file (and line). The report is rendered before any table is written. Tests added for each shape. |
| 2 | major | The status source cannot hold sensitivity, specificity or the calibration set (owner comment). | A. Optional validated `measure` object per entry; `StatusEntry`, `resolve_entry`; `StatusResolver`, `calibration_rows`; report section "Module calibration" and `run.json` rows ("not measured" without a `measure`). Plan 2 supplies the numbers. |
| 3 | minor | The install check lacks `pyyaml`. | A. |
| 4 | minor | Readers fail when the lock file cannot be opened; writers can starve. | A for the read-only case (read without a lock). N for starvation: stated as known gap 7. |
| 5 | minor | `inf` and `-inf` are read as numbers. | A. `math.isfinite`; test and mutation check. |
| 6 | minor | A repeated ID in `--taxon-map` is accepted. | A. `InputError` naming both lines. |
| 7 | minor | The unknown-taxon test passes only through a lazy lookup. | A. A second test with no status file; mutation check added. |
| 8 | minor | The AND-false status rule has no test. | A. Test and mutation check added. |
| 9 | minor | A byte order mark breaks FASTA, module table and taxon map reading. | A. Text inputs are read as `utf-8-sig`; tests added. |
| 10 | minor | README text lands under `# Author`; AGENTS.md names the implementer in advance; per-call reason and share of `not_assessable` are not reported. | A for README and AGENTS.md. Per-call reason and share are named as Plan 2 item (known gap 8). |

Not checked by the reviewer: `pre-commit run`, the full `tests/surface_glyco` suite (needs torch),
the other analysis test suites, `flock` behaviour on /bigdata.
