# Sorting-hat review fixes and wall/fusion testing: implementation plan

*2026-10-17. Draft 2. Plan for `docs/superpowers/specs/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing-design.md` (owner decisions D1 to D22, confirmed 2026-10-17). Author: Claude Code (claude-sonnet-5-5). No code is written by this plan's PR. Every code task starts with a test that fails on `main`. One commit per task. Nothing is pushed until the owner says so. Items marked PROPOSAL are not owner decisions. Finding ids `FD-xx` are from `docs/reviews/2026-10-17-decisions-d1-d18-review-fable.md`.*

Python 3.12. Run tests with `PYTHONPATH=$PWD/src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q -p no:cacheprovider`. The logic review's baseline is 818 passed, 9 skipped, 2 failed (macOS). The pre-commit hook runs ruff 0.3.5 and ruff-format and aborts a commit that it reformats: re-add and commit again. Line numbers are for `main` at `83d5f60`; re-read them before each edit.

Branches: one per PR, each from the merged state of the one before. Names (PROPOSAL): `sorting-hat-review-spec` (PR 1), `sorting-hat-cheap-fixes` (PR 2, PR #88, open), `sorting-hat-engine-fixes` (PR 3), `sorting-hat-review-tests` (PR 4), `wall-fusion-baseline` (PR 5), `wall-fusion-custom-hmm` (PR 6).

## Rules of this plan (decisions of the spec applied)

P1. **Tests first in PR 3 (D10, D22).** Task 3.0 writes every new test of PR 3 before any engine change. Each test that fails on `main` is marked `xfail(strict=True)` with the task id as the reason. When a task lands its change, the same commit removes the marker. A strict xfail that starts to pass without the marker removed fails the run, so no test stays marked by mistake.
P2. **Golden files in their own commit (D10).** No task regenerates golden files. Task 3.9 does it in one commit that contains only the golden files; the owner reviews that commit line by line. The expected rows of the gap cases are written by hand in Task 3.0 (not copied from engine output).
P3. **One config edit.** `categories.yaml` is edited once (D4) and `ENGINE_SEMANTICS` is bumped once (D5), so call files go stale once. Task 3.10 refreshes them (D21, D22).
P4. **Code location (D15).** The L1 validator and the tier derivation live in `src/cellsurface_sorting_hat/`. Their tests live in a directory that the blocking job runs. No blocking test imports `analysis/`. A new test directory is added to `.github/workflows/build_and_test.yml:46` in the same commit that creates it.
P5. **Defaults that keep old data readable.** A status file without the basis field loads as `assumed`; a measure-less entry is `unvalidated` (D17); a family table without the new columns is refused (the PR rewrites the table).
P6. **Nothing goes to SLURM without the owner (D21).** Job scripts are shown to the owner first; each submission needs an explicit go-ahead; heavy work never runs on the login node.
P7. **Names are PROVISIONAL** (spec Q2). Use the names of the spec; a rename is one search and replace.

## Stop points before PR 3

The owner (or the plan's author, where stated) must finish these before Task 3.0 starts:
- S1. **R0 per-stratum floor gate (D16, Q4).** Compute, from the Phase C `metrics.json`, eval table and cluster files on HPCC (read-only, small job; the owner approves the job), whether N-int plus N-sec alone meet the floors for R0: at least 20 negatives, at least 20 clusters, interval half-width at most 0.10 for specificity. Put the three numbers in this plan (a new line under this list) and in the PR description. If any floor fails, stop: D6 and D13 need the owner again.
  - **S1 result (run 2026-10-17, owner approved the run and, after seeing the numbers, accepted the gate as PASSED).** Read-only script on the HPCC login node over `phasec/eval_table.tsv.gz`, `clusters.tsv.gz` and `metrics.json`; N-int plus N-sec only (PM-TM and proteins in both left out). Specificity interval is Wilson on proteins, widened with Wilson on clusters (the bootstrap interval of `metrics.json` cannot be combined across strata). Floors: at least 20 negatives, at least 20 clusters, half-width at most 0.10.

    | Species (taxon) | Curated negatives | Clusters | Specificity | Half-width | Floors |
    |---|---|---|---|---|---|
    | S. cerevisiae (4932) | 3783 | 3154 | 0.966 | 0.006 | met |
    | C. albicans (5476) | 438 | 399 | 0.952 | 0.021 | met |
    | A. fumigatus (746128) | 45 | 38 | 1.000 | 0.046 | met |
    | A. nidulans (162425) | 163 | 150 | 0.994 | 0.017 | met |
    | C. neoformans H99 (5207) | 32 | 31 | 0.938 | 0.094 | met |
    | U. maydis (5270) | 27 | 23 | 0.926 | 0.116 | not met (half-width) |

    Reading: A. nidulans, the only species where R0 is `estimated`, meets every specificity floor without PM-TM, so D6 and D13 stand (the sensitivity floors are unchanged by D16). U. maydis fails the half-width floor but R0 is already `smoke` there (pooled half-width 0.144). C. neoformans shows 0.094 curated against 0.102 pooled with no PM-TM-only negatives (32 of 32): the difference is the interval method, not the exclusion, and R0 stays `smoke` there. Script and output: `r0_stratum_floor_check.py` and `result.tsv` (HPCC `~/s1_r0_floor/`); the script is to be committed with PR 3 for reproducibility. The strict rule "if any floor fails, stop" was applied by the owner as: the gate concerns species whose status depends on D16, so U. maydis does not stop it.
- S2. The owner confirms the names of the spec's Q2 or accepts them.
- S3. SETTLED by the owner on 2026-10-17: Q7 (the strict cluster rule applies to all measures) and Q8 (`narrow_ok=no` for rows D12 does not name; the column has meaning only on `pfam_adhesion` rows). No action.

## PR 1: this spec and plan

Task 1.1. Files: the three reviews in `docs/reviews/`, the spec, this plan. Steps: copy the reviews with a short header; write spec and plan; commit; no push. Test: none (documents). Done when the owner has read the spec and the stop points are answered or deferred.

## PR 2: cheap fixes (no engine change; PR #88)

### Task 2.1: stale statements about `active` and the narrow call (D2, SH-D2, SH-D3, WF-F1, WF-F2)

Files: `src/cellsurface_sorting_hat/modules/pfam.py` (docstring `:31-32` only), `docs/reports/cell-wall-fusion/04-repo-coverage-and-gaps.md` (`:22`, `:227`), `docs/superpowers/specs/2026-10-10-lectin-class-design.md` (`:266`, `:270`), `docs/reports/2026-10-07-classifier-runs-and-scale.md` (`:30`, `:70`), `tests/cellsurface_sorting_hat/test_docs_consistency.py` (new, in the blocking dir).

Steps:
1. Write the test first. It reads the four files and fails on these phrases: "passed its specificity test", "passes its specificity test", "promises a passed specificity test", "partition the assessable proteins", "a partition of the proteome". The test lists the paths; it does not scan the repo.
2. Reword the docstring: "No family of the module is enabled by the owner." Reword the report and spec lines: `active` means enabled by the owner; the tier column (PR 3) says how well tested.
3. Reword `classifier-runs-and-scale.md:30,70`: the calls are labels, not a partition; some proteins have no label; PR 3 adds `primary_label`.
4. Do not edit `categories.yaml`, `family_table.tsv` or `outputs.py` here. Do not edit the counts of the classifier-runs report: they are stale on main (FD-A9) and get a re-run plus a dated note in Task 3.10.

Tests: the new docs test (fails before, passes after).

### Task 2.2: errata and report fixes from the wall/fusion review (PROPOSAL; spec section 7)

Files: `docs/reports/cell-wall-fusion/01-cell-wall-reviews.md`, `02-cell-fusion.md`, `03-multicellularity-hyphae.md`, `04-repo-coverage-and-gaps.md`, errata logs.
Steps: for each item (WF-F1, WF-F2, WF-F4, WF-F5, WF-F10): re-read the cited line as it is now, apply the marked correction, keep the old text visible, date the marker. The owner confirms the Gel4 wording first. Tests: none; read-through against the review. Stop point: owner reads the diff.

### Task 2.3: macOS test failures and CI coverage (D10 PR 2, D15)

Files: `tests/cellsurface_sorting_hat/test_data_and_scripts.py`, `.github/workflows/build_and_test.yml`.
Steps:
1. Reproduce on macOS: `pytest tests/cellsurface_sorting_hat/test_data_and_scripts.py -q`. Probable names: `test_tmhmm_job_works_with_relative_paths` (`:393`) and `test_submit_modules_passes_absolute_paths_to_sbatch_and_supports_dry_run` (`:408`); confirm them.
2. Fix in the tests, not the HPCC scripts: skip when `realpath -m` is not supported (probe `realpath -m /nonexistent/x`), or use `grealpath`. Do not change `scripts/sorting_hat/*.sh`.
3. CI: add the five directories CI never runs (`cys_candidates`, `lectin_class_report`, `lectin_validation`, `secretome_discovery`, `wall_fusion_report`) to the `analysis-tests` matrix (`build_and_test.yml:57`), one pytest process each (comment at `:55-56`). They stay informational (`continue-on-error`, `:51`); the owner promotes a suite by hand (D15).
Tests: the two macOS tests pass or skip with a reason; the CI run lists five new suites. A suite that fails for a reason that is not the change gets a skip reason in the PR; do not hide it.

## PR 3: engine changes (tests first)

PR 3 is not marked ready until Task 3.10 (HPCC re-runs) is done. PROPOSAL: split into 3a (Tasks 3.2, 3.5, 3.7) and 3b (3.1, 3.3, 3.4, 3.6) if the owner wants smaller reviews; each half still starts with its share of Task 3.0.

### Task 3.0: all new tests first (D10, D22)

Files: the test files named in Tasks 3.1 to 3.7, `tests/cellsurface_sorting_hat/golden/gaps/` (new, with hand-written expected rows).
Steps: write every test listed under "Tests" in Tasks 3.1 to 3.7. Mark each that fails on `main` as strict xfail with the task id. Write the gap-case fixture and its expected `calls.long` and `calls.wide` lines by hand from the spec rules (D14, D18) before the engine is read again: six proteins: (1) no SP, repeat called; (2) SP, hydrophobin only; (3) SP, HsbA only; (4) SP, `pfam_hsba` error; (5) SP, `antigen_lookup` not applicable, all else not called; (6) SP, a hit in an `untested`, `narrow_ok=no` family; plus (7) SP with only the antigen rank called. Expected labels: 1 `mechanism_without_class`; 2, 3 `surface_attachment_candidate`; 4 `not_assessable`; 5 `other_surface_no_mechanism`; 6 `surface_attachment_candidate`; 7 `mechanism_without_class`. The run is `cellsurface_sorting_hat.cli` on the fixture.
Commit: one commit of tests only. Done when the suite is green with the strict xfails in place.

### Task 3.1: loader re-validation, `negative_basis` per stratum, measure-less entries (D7, D16, D17)

Files: `src/cellsurface_sorting_hat/status.py` (`:17-28`, `:74-94`, `:105-133`, `:136-164`), `call_status.py` (`:37-39`), `calibration/measure.py` (`:27-42`, `:56-66`, `:69-93`; measure-less rule at `:79-81`), `calibration/cli.py` (`:92-97`, `:442-529`), `calibration/phasec.py` (`:229`, `:259-269`), `docs/paper/03-status-and-validation-rules.md`, tests `test_taxonomy_status.py`, `calibration/test_calibration_cli.py`, `calibration/test_phasec.py`, `test_call_status.py`.

Produces:
- The basis field per stratum in `MEASURE_KEYS`; `_check_measure` accepts only the allowed values.
- `status_from_measure` uses only `curated` strata for `n_neg`, specificity and interval. An `assumed` stratum (PM-TM) is reported in its own row and does not count toward the floors.
- `MIN_TAXON` moves to `status.py`.
- `load_status_source` computes the derived status and refuses a stronger stored status, a taxon below `MIN_TAXON`, and a leakage value other than `none` above `smoke` (same regex as the call-file loader; skipped for an entry with the D13 key, Task 3.5).
- A measure-less entry is `unvalidated`; any other status is refused (D17). Remove the "cap at smoke" wording from code comments.
- `truth`: option for the basis; `build_measure` stores it; `phasec` sets it per stratum.

Tests:
1. A module status file with taxon 1 and `estimated` raises on load.
2. `estimated` with a measure that derives `smoke` raises; a lower status loads.
3. 25 positives, 25 negatives, 22 and 22 clusters, narrow intervals: all-`curated` gives `estimated`; adding a PM-TM `assumed` stratum does not change the result and is reported in its own row; a missing field gives `smoke`.
4. Leakage `partial` with `estimated` raises; with `smoke` it loads.
5. A measure-less entry with `smoke` raises; with `unvalidated` it loads.
6. Each of the 20 committed files under `docs/reports/data/sorting_hat/call_status/` still loads (all carry a measure; all are `smoke`).
7. Strict cluster rule (D3, S3): a measure without `n_clusters_pos` or `n_clusters_neg` cannot be `estimated`. The 20 committed files still load.
8. The call-file loader gives the same result as before for the existing fixtures.

Mutations: loader ignores the cap; cap ignores a missing key; taxon floor off by one; loader compares status strings instead of strengths; PM-TM counted in `n_neg`; missing cluster key skips the floor.

### Task 3.2: `not_applicable` and the `other_*` rule (D5, D18)

Files: `src/cellsurface_sorting_hat/logic.py` (`:7-45`), `engine.py` (`:30`, `:174-215`, `:236-253`, `:348-349`, `:497-508`), `outputs.py` (`:276-289`), `calibration/cli.py` (`:132`, `:290-293`, `:473-488`), `categories.yaml` (comment `:11`), `docs/class_descriptions.yaml`, tests `test_logic.py`, `test_engine.py`, `calibration/test_calibration_cli.py`, `test_outputs.py`. Also grep `calibration/panel.py` and `analysis/` readers of `calls.long` for the three-value assumption (FD-A2); list each finding in the PR.

Steps:
1. `logic.py`: add `NOT_APPLICABLE`; implement the four-value table of spec 3.5 (NA acts as `not_called` inside `and` and `or`; `and` and `or` never return NA; `not(NA)` is `not_applicable`).
2. `engine._row` returns the state; `_eval` maps a `not_applicable` row to `NOT_APPLICABLE` at a leaf.
3. `_eval_other`: called mechanism gives `not_called`; else unknown gives `not_assessable` with `other_basis` listing it; an NA mechanism is skipped.
4. An NA leaf adds no contributor. The basis lists `not_applicable: <module>`; a `not_called` record caused by an NA leaf never has an empty basis (change at `engine.py:348-349`).
5. Bump `ENGINE_SEMANTICS` (`engine.py:30`) with one comment line.
6. `calibration/cli.py:132,290-293`: accept the value; `:473-488`: count it with the not-assessable rows and report it; never a negative.
7. `outputs.py:276-289`: fourth count column; rename the `other_basis` section title.

Tests:
1. Probe P2 (SP called, `pfam_hsba` row `error`): `other_surface_no_mechanism` is `not_assessable`, `other_basis` is `hsba_domain`.
2. Probe P3 (`antigen_lookup` `not_applicable`): the rank is `not_applicable`; `serodiagnostic_marker_candidate` is `not_called` with basis `not_applicable: antigen_lookup`; `other_surface_no_mechanism` is `called`.
3. An SP protein outside RS is never `serodiagnostic_marker_candidate = called` (FD-A1).
4. A called mechanism with another unknown one: `other_*` is `not_called`.
5. The whole truth table of spec 3.5 as a unit test, all combinations, `not` included (`not(NA)=NA`).
6. `calibrate truth` on a table with `not_applicable` exits 0 and reports the count (fails before: exit 2, FD-A2). The report's Calls table shows the fourth column.
7. `call_hash` differs from the value stored in a fixture with the old semantics.

Mutations: NA treated as unknown; NA dropped from operands; unknown skipped as today; NA leaf contributes a status; semantics not bumped.

### Task 3.3: tier derivation in `src/` and the measurement file (D3, D11, D19)

Files: `src/cellsurface_sorting_hat/family_tier.py` (new, PROPOSAL name), `status.py` (`:128-130`), `data/sorting_hat/family_measure.tsv` (new), `analysis/family_tests/README.md`, `tests/cellsurface_sorting_hat/test_family_tier.py` (new).

Steps:
1. The module imports `status_from_measure`; it does not copy thresholds. Mapping: no row, no positive or no negative gives `untested`; `descriptive_shared_sets` is capped at `smoke`; `estimated` needs `per_family_reviewed` and `status_from_measure == estimated`; a missing `n_clusters_pos` or `n_clusters_neg` does not meet the floor (fix at `status.py:128-130`); positives mostly `unknown` or `in_reference` cap at `smoke` (D19). A command (`python -m cellsurface_sorting_hat.family_tier --check|--write`) compares to or rewrites the table column.
2. Convert `docs/reports/data/sorting_hat/family_tests/family_tests.tsv` into `family_measure.tsv` with `measure_basis=descriptive_shared_sets`, the E1 hit count in a `coverage` column (never `sensitivity`), the sources and date, the Pfam and HMMER versions. Hyr1 has no E1 hit: its row is absent or has zero positives, so it stays `untested`.
3. Write one sentence in the table header, here and in Task 3.8: the family tier describes the family model's own evidence; the call status describes the call's measurement; neither implies the other.
4. Run `--write` and commit the table together with Task 3.4.

Tests:
1. The committed table equals the derived tiers.
2. A copy of the table with one stronger tier fails.
3. Fixtures: 25/25/22/22 with narrow intervals and `per_family_reviewed` derive `estimated`; the same with `descriptive_shared_sets` derives `smoke`; 5/5 derive `smoke`; no row derives `untested`.
4. A row without `n_clusters_pos` derives at most `smoke` (fails before, FD-A3).
5. A family whose positives are mostly `in_reference` derives at most `smoke`.
6. Grep test: no file under `tests/cellsurface_sorting_hat` that is a blocking suite imports from `analysis/`.

Mutations: ignores clusters; ignores `measure_basis`; copies thresholds; `--write` changes other columns.

### Task 3.4: family table columns, narrow flag and narrow call (D2, D4, D12)

Files: `data/sorting_hat/family_table.tsv`, `modules/pfam.py` (`:9-24`, `:35-80`, `:125-171`), `modules/cli.py` (`:225-232`, `:252-300`), `categories.yaml` (`:53-75`), `outputs.py` (`:29-34`, `:292-318`), `docs/class_descriptions.yaml`, `docs/CLASSES.md` (regenerated), `analysis/class_table/make_classes_md.py` if needed, tests `modules/test_pfam.py` (or the file that tests `pfam_rows`), `test_cli.py`, `test_engine.py`.

Steps:
1. Table: add `specificity_tier`, `narrow_ok`, `narrow_ok_by`, `narrow_ok_date`, and an evidence-summary field for the sign-off (members, hits, leakage, near-neighbour negative hits; FD-A15). Tiers come from Task 3.3. `narrow_ok` draft values (D12): `yes` for PA14, Flocculin, Flocculin_t3, Candida_ALS_N, Candida_ALS, Flo11, GLEYA; `no` for CFEM, PIR1-like_C, Hyr1. **The owner signs each row in the PR** (who and date); until signed the column is `no`. Rows D12 does not name: per S3.
2. `pfam.py`: read the columns; `Family` gets `tier` and `narrow_ok`; validation as spec 3.4. `COLUMNS` (`:24`) becomes per-module: `pfam_adhesion` gets `hit_narrow`, the other three do not (FD-A9).
3. `modules/cli.py`: `_family_digest` covers tier and `narrow_ok`; `params` lists the narrow families. The call hash covers `family_table.tsv` through these fields.
4. `categories.yaml`: add `wall_family_narrow_domain`; `cell_wall_adhesion_candidate` reads it; `surface_attachment_candidate` stays on the broad call. This changes `config_sha256`.
5. Report: warn when no family has `narrow_ok=yes` (FD-B8). `outputs.py` and descriptions: say which families feed which call; fix the label at `outputs.py:294`.
6. Record before and after counts of the classifier runs for the PR description (D12); the numbers come from Task 3.10.

Tests:
1. SP protein, hit in an `untested`, `narrow_ok=no` family: broad `called`, narrow `not_called`.
2. Same protein, family `smoke`, `narrow_ok=yes`: both `called`.
3. `load_family_table` refuses: missing new columns; tier outside the vocabulary; `narrow_ok=yes` with tier `untested`; `narrow_ok=yes` without `narrow_ok_by` or date.
4. `hit_narrow == 1` implies `hit == 1`.
5. `pfam_hydrophobin`, `pfam_allergen` and `pfam_hsba` tables have no `hit_narrow` column.
6. `_family_digest` changes when `narrow_ok` or the tier changes.
7. An old `pfam_adhesion` table without `hit_narrow` stops the run with `missing column`.
8. The report warns when every `narrow_ok` is `no`.
9. `tests/class_table/test_class_table.py` passes after regeneration.

Mutations: narrow flag ignores tier; ignores `narrow_ok`; digest ignores `narrow_ok`; narrow call reads the broad flag.

### Task 3.5: `derived` marker, `leakage_exception`, note wording (D6, D13)

Files: `src/cellsurface_sorting_hat/engine.py` (`:334-349`), `status.py` (`:17-28`, `:196-201`), `calibration/phasec.py` (`:41`, `:250-251`, `:269`), `outputs.py` (`:47-48`, `:206`), `docs/CLASSES.md`, tests `test_engine.py`, `calibration/test_phasec.py`, `test_outputs.py`.

Steps:
1. Add the key `leakage_exception` (reason, signer, date) to `MEASURE_KEYS`. `phasec` writes it for R0 entries.
2. Reword `LEAKAGE_NOTE` (`phasec.py:41`) to start `leakage: unknown (` with the prose in brackets (FD-A6). The loader (Task 3.1) skips the leakage cap only for an entry with the key.
3. Where an entry has the key, `status_basis` carries "leakage unmeasured (SignalP training overlap), uncapped by owner decision".
4. Derived marker for composite and `other_*` records (a `derived_from` column naming the leaf calls; format is spec Q5). `calls.long` gains it; `calls.wide` does not.

Tests:
1. A composite record has the marker (fails before).
2. A record whose status comes from the R0 entry has the exception text; a leaf `repeat` record does not.
3. `leakage_of()` of every R0 entry returns a member of `LEAKAGE_VALUES`.
4. The loader accepts an `estimated` R0 entry with the key; it refuses the same entry without the key.
5. `phasec` output for the A. nidulans fixture has the key and the per-stratum basis.

Mutations: marker on leaf calls; marker missing on `other_*`; loader exempts every entry; note not reworded.

### Task 3.6: R0 snapshot check (D6)

Files: `modules/cli.py` (`:428-459`), `cli.py` (`:366-372`, `:380-412`), `outputs.py` (`_warnings`, `RunInfo` at `:53-75`), tests `modules/test_module_cli.py`, `test_cli.py`, `test_outputs.py`.
Steps: record a digest of the `(id, value)` pairs in `_condition_table`; after `load_modules` recompute it from the engine's R0 table and compare with each `pfam_*` run record; a mismatch gives a warning and a `run.json` field; a missing digest in an old record gives a note, not an error.
Tests: `pfam_adhesion` built against R0 table A, run with table B: warning (fails before); same table: no warning; old record without digest: note; `run.json` carries the field.
Mutations: compare only `params_hash`; compare digests of ids only.

### Task 3.7: `primary_label` (D1, D14)

Files: `src/cellsurface_sorting_hat/labels.py` (new), `outputs.py` (`:93-127`, `:276-285`), `cli.py` (call site near `:416-417`), `docs/class_descriptions.yaml`, `docs/CLASSES.md`, tests `test_labels.py` (new), `test_outputs.py`.
Steps: `primary_labels(records, cfg)` returns one label per protein for the default gate with the precedence of spec 3.1: `cell_wall_adhesion_candidate`; `surface_attachment_candidate` (with basis sub-label); `other_surface_no_mechanism`; `other_not_surface`; `mechanism_without_class`; `not_assessable`. Allergen and antigen calls are not label sources. Add the column to `write_wide` and the "Primary label" count table to the report.
Tests: P1 (no SP, repeat02 called) gets `mechanism_without_class`; SP with HsbA only gets `surface_attachment_candidate`; P2 gets `not_assessable`; an SP protein with a mechanism no composite reads (a test config with a `sowgp_ortholog`-style call) gets `mechanism_without_class`; an SP protein with only the antigen rank gets `mechanism_without_class`; every protein in the fixtures and in the committed S288C run gets exactly one label; report counts add up to the protein count.
Mutations: precedence reversed; fallback missing; label computed for a non-default variant.

### Task 3.8: documentation of the new rules (D2, D3, D7, D13, D14, D18)

Files: `docs/paper/03-status-and-validation-rules.md`, `docs/CLASSES.md` (regenerated), `docs/class_descriptions.yaml`, `AGENTS.md` (one entry for this work, by the implementing agent).
Steps: state `active`, tier, `narrow_ok`, `not_applicable`, `derived`, the per-stratum basis, `leakage_exception` and `primary_label` once each; link to the spec. Add the sentence of Task 3.3 step 3 (tier versus call status) to `CLASSES.md`. Add the interim sentence of Task 3.10 until the status files are current. Test: the docs tests of Tasks 2.1 and 4.3 and the class-table test.

### Task 3.9: golden regeneration (own commit; D10)

Files: `tests/cellsurface_sorting_hat/golden/*` only.
Steps: regenerate once after Tasks 3.1 to 3.8 are in. The commit contains no other file. The PR description lists every changed line by kind (new column, new call, changed composite basis, new value). The owner reviews this commit separately. The hand-written gap-case rows of Task 3.0 must match the regenerated output without edits; any difference is a defect to explain, not to overwrite.

### Task 3.10: HPCC re-runs, acceptance gate (D21, D22; FD-A11, FD-B3)

Files: `scripts/sorting_hat/` job scripts (new or existing), `docs/reports/data/sorting_hat/call_status/*` (new files; old files stay as history), `docs/CLASSES.md` (interim sentence), `docs/reports/2026-10-07-classifier-runs-and-scale.md` (dated note).
Steps:
1. Write the job scripts: `calibrate phasec` for R0 (so it carries `leakage_exception`, the reworded note and the per-stratum basis); `calibrate truth --call-status` for `tandem_repeat_protein` in S288C and in C. albicans (as in M2d and M5b2); rebuild of the four Pfam modules; a re-run of the classifier runs for the before/after counts (D12). Use `${SCRATCH:?}`, no `/tmp`, sized per job.
2. **Show the scripts to the owner. Submit nothing to SLURM without the owner's explicit go-ahead, each time. No heavy work on the login node.**
3. Commit the refreshed status tables like the existing ones. Note whether numbers changed. Update the `CLASSES.md` citations.
4. In the same commit as the first refreshed file, add one interim sentence to `CLASSES.md`: call statuses are `unvalidated` where the call file is stale until the refresh is merged. Remove it when all 20 files are current.
5. Add a dated note to the classifier-runs report: counts were measured on a given commit and family table; the PR shows the new counts (FD-A9). This is a re-run, not an edit framed as a D4 consequence.
Tests: every committed new status file loads (Task 3.1 test 6); `stale_reason` (`call_status.py:179,194`) is empty for each new file with the final config.
Done when: the owner has reviewed the new files and the golden commit; then the PR is marked ready.

## PR 4: tests that did not fit PR 3

### Task 4.1: L1 shared control-set validator (D15)

Files: `src/cellsurface_sorting_hat/controls.py` (new, PROPOSAL name), `tests/cellsurface_sorting_hat/test_controls_validator.py` (new, in the blocking directory), `tests/cellsurface_sorting_hat/fixtures/controls/` (tiny good and bad sets), a command entry (`python -m cellsurface_sorting_hat.controls data/controls/<dir>`).
Steps: tests first; the validator checks COMMON-RULES columns (`docs/agent-tasks/COMMON-RULES.md`), `label` in `positive|negative`, `leakage` in the five values, `tuned` in `yes|no`, one row per accession, every accession in `sequences.faa` and `clusters.tsv`, no cluster mixing a positive and a negative unless the row says so (`mixed_cluster_with_positive` in the shared set), counts in `counts.md` equal counts from `controls.tsv`.
Tests: the good fixture passes; each defect (missing column, bad leakage value, missing sequence, row missing from clusters, count mismatch, duplicate accession) fails with a message naming the row. Run the validator on `data/controls/hard-negatives-shared/` and `data/controls/repeat-mechanism/`; report any real defect to the owner, do not fix it silently.
Mutations: ignores `clusters.tsv`; ignores `counts.md`.

### Task 4.2: L3 golden extensions

Files: `tests/cellsurface_sorting_hat/golden/` and `test_golden_gaps.py` (extend the files of Task 3.0).
Steps: add the cases that did not fit PR 3: a fixture with an NA row at every module that writes `not_applicable`; a protein with several unknown mechanisms; a variant run for a `per_variant` call. Expected rows by hand first.
Mutations: swap two precedence entries; skip the narrow flag; read the broad flag in the narrow call.

### Task 4.3: docs-consistency test (L3)

Files: `tests/cellsurface_sorting_hat/test_docs_consistency.py` (extend Task 2.1's file).
Tests: (a) every `pfam_adhesion` family with `narrow_ok=yes` is named in the generated `CLASSES.md` text for the narrow call; (b) every active `pfam_adhesion` family is named for the broad call; (c) `CLASSES.md` mentions each call of `categories.yaml`; (d) the tier and the derived tier agree.
Mutation: remove a family from the description; the test must fail.

## PR 5: feasibility stage, then baseline recall (HPCC)

Stop point before: the owner has read the feasibility report (Task 5.1) and accepted or adjusted the clade floor (spec Q1) before any model is evaluated. The owner confirms the two hard-negative sources (spec section 9, item 12).

### Task 5.1: feasibility stage (D8, D9 step 0; FD-A8)

Files: `analysis/wall_fusion_truth/map_kiss.py` (new), `docs/reports/2026-10-XX-kiss-feasibility.md` (new). No model is run.
Steps:
1. Read `docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` (181 rows; `protein_id_jgi`, `cluster_number`, `node_of_emergence`). Map each JGI ID to a sequence in the repo proteomes or to UniProt.
2. Report: rows with and without an ID (37 have none), rows mapped and lost, rows per organism code, clusters per clade (clade from the organism code, not from `node_of_emergence`, which is the cluster's ancestral node).
3. Reconcile counts: adhesion 51 rows / 26 clusters, wall 122 / 44, fusion 8 / 6, against "23 clusters, 45 curated genes" of report 06.
4. Show how many clades each floor (3, 5, 8 independent clusters) would admit for hold-out; show which slices fall back to descriptive recall and leave-cluster-out. The floor of 5 is provisional.
5. Name the second source for Basidiomycota and non-hyphal clades (UniProt reviewed entries plus literature rows) and count what it adds.
Tests: a unit test on a 3-row fixture checks the mapping and the loss count. Stop point: owner reads the report.

### Task 5.2: build the control set (D8)

Files: `data/controls/wall-fusion-kiss2019/` (COMMON-RULES files), validated by Task 4.1. Do not pad. Rows lost in 5.1 go in `rejected.tsv` and `unverified_notes.md`.

### Task 5.3: cluster and split

Files: `analysis/wall_fusion_truth/cluster_split.sh` (HPCC, MMseqs2 30% identity, 50% coverage, as in `docs/superpowers/plans/2026-10-10-allergen-predictor-plan.md` T0.4), `analysis/wall_fusion_truth/seed_split.py` (new).
Steps: cluster positives and hard negatives; compute identity to the nearest Pfam seed of the models in `analysis/wall_fusion_report/wall_fusion_models.tsv` and of the active family rows (`analysis/lectin_validation/measure.py:328-329,402-409`); label each cluster `seed-like` or `seed-independent`. D21 rules apply: owner approves each submission; `${SCRATCH:?}`; no `/tmp`.
Test: a fixture with identities 35% and 20% lands in the two groups; a cluster with any member at 30% or more is `seed-like`.

### Task 5.4: baseline recall per slice (D8, D9 step 1, D20)

Files: `analysis/wall_fusion_truth/baseline_recall.py`, `analysis/wall_fusion_truth/submit_baseline.sbatch`, `docs/reports/2026-10-XX-wall-fusion-baseline.md`.
Steps: run `hmmsearch --cut_ga`. Report per slice: adhesion rows, recall of the active adhesion Pfam families (family coverage, not the call's sensitivity); wall rows, recall of the 21 wall models; fusion rows, descriptive only. Report recall per cluster and per clade where the floor allows, else leave-cluster-out; report `seed-like` and `seed-independent` apart; label every number with its n; no pooled number. Wall rows are scored as a descriptive false-hit panel for the narrow call and are never negatives and never truth for the adhesion call. No tier is changed by this report. A model with no GA line is reported, not skipped.
Tests: recall on a fixture equals the hand count; the script writes no call and no status; the wall rows never enter a negative count.
Output: the table of clusters Pfam misses. That table, not a guess, names the PR 6 candidates.

### Task 5.5: feed L2 and L5 (PROPOSAL)

Steps: write the family rows that now have family-own positives (D19) into `data/sorting_hat/family_measure.tsv` with the right `measure_basis`; run the tier check. Pin the numbers with a tolerance and the Pfam and HMMER versions (L5, PROPOSAL). Any `narrow_ok` change is an owner decision (D4).

## PR 6: custom HMMs with leave-clade-out

Own spec before work starts, written from the Task 5.4 table. Outline:
- 6.1 For each cluster Pfam misses, align orthologs (MAFFT or the method of `2026-10-08-hydrophobin-custom-hmm-design.md`), `hmmbuild`, set GA from the training clades only.
- 6.2 Leave-clade-out (`analysis/hydrophobin_truth/lco.py` is the precedent), or leave-cluster-out where the clade floor is not met; report per held-out unit.
- 6.3 A custom model gets at most `smoke` until 6.2 is committed (the derivation refuses a higher tier for a row without that measurement). No ESM model (issue #15; gate of 300 clusters).
- Tests: the derivation refuses `estimated` for a custom model without the hold-out field; the held-out unit never appears in the build alignment (a test on the build manifest).

## Out of scope

- ESM or fine-tuned models (D9 step 4, issue #15).
- A residual `other_no_surface_with_mechanism` call (D1 chose a derived label).
- A `secondary_label` for allergen and antigen calls (spec Q6).
- Changes to the `estimated` thresholds (only the missing-cluster rule changes, D3).
- Wall/fusion model-list and label changes (WF-F8, WF-F9) unless the owner decides them.
- The NRP and GPU work of the model-review plan (issues #15, #16).

## Self-review against the spec

| Spec item | Task |
|---|---|
| D1, D14 `primary_label`, wording | 2.1 (wording), 3.0, 3.7, 4.2 |
| D2 tier column, `active` meaning, docs | 2.1, 3.4 |
| D3 derived tier in `src/`, strict cluster rule, tier-versus-status sentence | 3.3, 3.8, 4.3 |
| D4 narrow call, per-module columns, call hash, warning | 3.4 |
| D5, D18 `other_*`, `not_applicable`, truth table, `truth` | 3.2, 4.2 |
| D6 R0 text, `derived`, snapshot check | 3.5, 3.6 |
| D7, D16, D17 loader, per-stratum basis, measure-less | S1, 3.1 |
| D8, D20 truth set, slices | 5.1 to 5.4 |
| D9 stages | 5.1 (step 0), 5.4 (step 1), PR 6 |
| D10 order | PR headings, 3.0, 3.9 |
| D11, D19 family measurement | 3.3 |
| D12 `narrow_ok` | 3.4 |
| D13 key and note | 3.5 |
| D15, D22 code location, tests first | P1 to P4, 3.0, 4.1 |
| D21 HPCC re-runs | 3.10 |
| Testing framework: L0 2.3; L1 4.1; L2 3.3, 5.5; L3 3.0, 4.2, 4.3; L4 and L5 proposals | as listed |
