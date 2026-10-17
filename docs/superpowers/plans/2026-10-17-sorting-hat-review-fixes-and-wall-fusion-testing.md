# Sorting-hat review fixes and wall/fusion testing: implementation plan

*2026-10-17. Draft 1, not reviewed. Plan for `docs/superpowers/specs/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing-design.md` (owner decisions D1 to D10 of 2026-10-17). Author: Claude Code (claude-sonnet-5-5). No code is written by this plan's PR. Every code task starts with a test that fails on `main`. One commit per task. Nothing is pushed until the owner says so. Items marked PROPOSAL are not owner decisions; the owner confirms them before the task starts.*

Python 3.12. Run tests with `PYTHONPATH=$PWD/src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q -p no:cacheprovider`. The logic review's baseline is 818 passed, 9 skipped, 2 failed (macOS). The pre-commit hook runs ruff 0.3.5 and ruff-format and aborts a commit that it reformats: re-add and commit again. Line numbers are for `main` at `83d5f60`; re-read them before each edit.

Branches: one branch per PR, each from the merged state of the one before. Names (PROPOSAL): `sorting-hat-review-spec` (PR 1), `sorting-hat-cheap-fixes` (PR 2, exists and has no commits yet), `sorting-hat-engine-fixes` (PR 3), `sorting-hat-review-tests` (PR 4), `wall-fusion-baseline` (PR 5), `wall-fusion-custom-hmm` (PR 6).

## Decisions taken in this plan (the spec left them open or implicit)

P1. **PR 3 order.** Inside PR 3 the order is: loader and `negative_basis` (D7), `not_applicable` (D5), family-table columns and the narrow flag (D2, D3, D4), `derived` marker and R0 snapshot check (D6), `primary_label` (D1). Reason: D7 and D5 do not touch `categories.yaml`; D4 does and makes every call file stale, so it comes after them. D1 comes last because its tests need the final call set.
P2. **One config edit.** `categories.yaml` is edited once (D4) and `ENGINE_SEMANTICS` is bumped once (D5), both inside PR 3, so call files go stale once. Task 3.9 re-measures them.
P3. **Derived tier, no hand edit.** The family table keeps the tier column (so a reader sees it), but the test of Task 3.5 fails when it differs from the derived value. The derive script has `--write` to update the table. Nobody edits the column by hand.
P4. **Defaults that keep old data readable** (spec section 6): a status file without `negative_basis` loads as `assumed`; a family table without the new columns is refused (the PR rewrites the table).
P5. **Labels in code.** `primary_label` precedence lives in code with a test, not in `categories.yaml` (spec 3.1; OPEN Q1).
P6. **Open questions block tasks.** Q1 (precedence), Q3 (`not_applicable` in logic), Q5 (Phase C negatives and the R0 exception key), Q7 (initial `narrow_ok`) must be answered by the owner before Tasks 3.1, 3.2, 3.5 and 3.6 start. The plan lists the stop points.

## PR 1: this spec and plan

Task 1.1. Files: the two reviews in `docs/reviews/`, the spec, this plan. Steps: copy the reviews with a 3-line header; write spec and plan; one commit; no push. Test: none (documents). Done when the owner has read the spec and answered the stop-point questions or deferred them.

## PR 2: cheap fixes (no engine change)

### Task 2.1: stale statements about `active` and the narrow call (D2, SH-D2, SH-D3, WF-F1, WF-F2)

Files: `src/cellsurface_sorting_hat/modules/pfam.py` (docstring `:31-32` only), `docs/reports/cell-wall-fusion/04-repo-coverage-and-gaps.md` (`:22`, `:227`), `docs/superpowers/specs/2026-10-10-lectin-class-design.md` (`:266`, `:270`), `docs/reports/2026-10-07-classifier-runs-and-scale.md` (`:30`, `:70`, D1 wording), `tests/cellsurface_sorting_hat/test_docs_consistency.py` (new).

Steps:
1. Write the test first. It reads the four files and fails on each phrase: "passed its specificity test", "passes its specificity test", "promises a passed specificity test", "partition the assessable proteins", "a partition of the proteome". The test file lists the paths; it does not scan the whole repo.
2. Reword the docstring: "No family of the module is enabled by the owner." Reword the report and spec lines: `active` means enabled by the owner; the tier column (PR 3) says how well tested.
3. Reword `classifier-runs-and-scale.md:30,70`: the calls are labels, not a partition. Add the sentence that some proteins have no label and that PR 3 adds `primary_label`.
4. Do not edit `categories.yaml`, `family_table.tsv` or `outputs.py` here (they change in PR 3).

Tests: the new docs test (fails before, passes after).

### Task 2.2: errata and report fixes from the wall/fusion review (PROPOSAL; not in D1 to D10)

Files: `docs/reports/cell-wall-fusion/01-cell-wall-reviews.md`, `02-cell-fusion.md`, `03-multicellularity-hyphae.md`, `04-repo-coverage-and-gaps.md`, errata logs.

Steps: for each item of spec section 7 (WF-F1, WF-F2, WF-F4, WF-F5, WF-F10): re-read the cited line in the file as it is now (the review used the errata worktree; line numbers differ), apply the marked correction, keep the old text visible, date the marker. The owner confirms the Gel4 wording first (WF-F4 is an inference).
Tests: none; the check is a read-through against the review. Stop point: owner reads the diff.

### Task 2.3: macOS test failures and CI coverage (D10 PR 2)

Files: `tests/cellsurface_sorting_hat/test_data_and_scripts.py`, `.github/workflows/build_and_test.yml`, possibly `pyproject.toml`.

Steps:
1. Reproduce on macOS: `pytest tests/cellsurface_sorting_hat/test_data_and_scripts.py -q`. The two tests that run job scripts with `realpath -m` are probably `test_tmhmm_job_works_with_relative_paths` (`:393`) and `test_submit_modules_passes_absolute_paths_to_sbatch_and_supports_dry_run` (`:408`); I did not run them. Confirm the names.
2. Fix in the tests, not the HPCC scripts: skip when `realpath -m` is not supported (a `pytest.mark.skipif` on a probe `realpath -m /nonexistent/x`), or call a GNU `grealpath` if present. Do not change `scripts/sorting_hat/*.sh` (they run on Linux); the owner decides if they should become portable.
3. CI: add the five directories that CI never runs, `cys_candidates`, `lectin_class_report`, `lectin_validation`, `secretome_discovery`, `wall_fusion_report`, to the matrix of `analysis-tests` (`build_and_test.yml:57`). Run each locally first; each directory needs its own pytest process because of conftest path rules (comment at `:55-56`).
4. OPEN: the `analysis-tests` job is `continue-on-error: true` (`:51`), so it cannot block a PR. Ask the owner whether the new suites go in a blocking job. Default: informational, as the other seven.

Tests: the two macOS tests pass or skip with a reason on macOS; the CI run lists five new suites. If a suite fails in CI for a reason that is not the change (network, HMMER), mark it with a skip reason and write the reason in the PR; do not hide it.

## PR 3: engine changes

Stop points before this PR: owner answers Q1, Q3, Q5, Q7 and Q8 of the spec.

### Task 3.1: loader re-validation and `negative_basis` (D7)

Files: `src/cellsurface_sorting_hat/status.py` (`:17-28`, `:74-94`, `:105-133`, `:136-164`), `src/cellsurface_sorting_hat/call_status.py` (`:37-39` constants), `calibration/measure.py` (`:27-42`, `:56-66`, `:69-93`), `calibration/cli.py` (`:92-97`, `:442-529`), `calibration/phasec.py` (`:269`), `docs/paper/03-status-and-validation-rules.md`, `tests/cellsurface_sorting_hat/test_taxonomy_status.py`, `tests/cellsurface_sorting_hat/calibration/test_calibration_cli.py`, `.../test_phasec.py`, `.../test_call_status.py`.

Produces:
- `negative_basis` in `MEASURE_KEYS`; `_check_measure` accepts `curated` or `assumed`.
- `status_from_measure` caps at `smoke` unless `negative_basis == "curated"`.
- `MIN_TAXON` moves to `status.py`; `call_status.py` imports it.
- `load_status_source` computes the derived status and refuses a stronger stored status, a taxon below `MIN_TAXON`, and a leakage value other than `none` with a status above `smoke` (from `measure.notes`, same regex as the call-file loader; skipped for an entry with the exception key of Task 3.6).
- `truth`: required `--negatives {curated,assumed}`; `build_measure` stores it. `phasec` sets it per Q5.
- Entries without `measure`: per Q8 (default: capped at `smoke`).

Tests (write first):
1. A module status file with taxon 1 and `estimated` raises on load. (Fails before.)
2. A file with `estimated` and a measure that derives `smoke` raises. A file with a lower status than derived loads.
3. 25 positives, 25 negatives, 22 and 22 clusters, narrow intervals: `curated` gives `estimated`; `assumed` and a missing key give `smoke`.
4. Leakage `partial` with `estimated` raises on load; with `smoke` it loads.
5. `truth` without `--negatives` exits with an argparse error; the written measure carries the value.
6. The call-file loader gives the same result as before for the existing fixtures, except the new cap.
7. Existing committed status files under `docs/reports/data/sorting_hat/call_status/` still load (a test loads each; the hydrophobin and tandem-repeat files are `smoke`).

Mutations: loader ignores the cap; cap ignores a missing key; taxon floor off by one; loader compares status strings instead of strengths.

### Task 3.2: `not_applicable` and the `other_*` rule (D5)

Files: `src/cellsurface_sorting_hat/logic.py`, `engine.py` (`:30`, `:180-215`, `:236-253`, `:497-508`), `outputs.py` (`:276-289`), `calibration/cli.py` (`:473-488`), `categories.yaml` (comment only, `:11`), `docs/class_descriptions.yaml`, tests `test_logic.py`, `test_engine.py`, `calibration/test_calibration_cli.py`, `test_outputs.py`.

Steps:
1. Tests first (below).
2. `logic.py`: add `NOT_APPLICABLE`; `VALUES`; Kleene rules per spec 3.5 (Q3).
3. `engine._row` returns the state too; `_eval` maps state `not_applicable` to `NOT_APPLICABLE`.
4. `_eval_other`: called mechanism gives `not_called`; else unknown mechanism gives `NOT_ASSESSABLE` with `other_basis` listing the unknown mechanisms; `not_applicable` is ignored.
5. Bump `ENGINE_SEMANTICS` to `"2"`. One comment line says why.
6. `truth`: a `not_applicable` call value is counted with the not-assessable rows and reported; it is never a negative.
7. `outputs`: new count column and a changed `other_basis` section title.

Tests (write first):
1. Probe P2 (SP called, `pfam_hsba` row `error`): `other_surface_no_mechanism` is `not_assessable` and `other_basis` is `hsba_domain`. (Fails before: `called`.)
2. Probe P3 (`antigen_lookup` `not_applicable`): rank is `not_applicable`; `other_surface_no_mechanism` is `called`; basis empty.
3. A called mechanism with another unknown one: `other_*` is `not_called`.
4. Truth table with a `not_applicable` value: not in the negatives; counted in the notes.
5. `call_hash` differs from the value stored in a fixture with semantics `"1"`.
6. Logic table for the Kleene rules with the fourth value (all combinations).

Mutations: `not_applicable` treated as unknown; unknown skipped as today; semantics not bumped.

### Task 3.3: family table columns, tier, narrow flag, narrow call (D2, D3, D4)

Files: `data/sorting_hat/family_table.tsv`, `src/cellsurface_sorting_hat/modules/pfam.py` (`:9-24`, `:35-80`, `:125-171`), `modules/cli.py` (`:225-232`, `:252-300`), `categories.yaml` (`:53-75`), `outputs.py` (`:29-34`, `:292-318`), `docs/class_descriptions.yaml`, `docs/CLASSES.md` (regenerated), `analysis/class_table/make_classes_md.py` if needed, tests `modules/test_pfam.py` (or the existing file for `pfam_rows`), `test_cli.py`, `test_engine.py`, `golden/*`.

Steps:
1. Tests first (below).
2. Table: add `specificity_tier`, `narrow_ok`, `narrow_ok_by`, `narrow_ok_date`. Initial tier values come from Task 3.5 (until then `untested`; the test of 3.5 is added in the same PR, so this task ends with all tiers `untested` and `narrow_ok=no`, unless the owner gives Q7 values). Initial `narrow_ok`: owner input (Q7).
3. `pfam.py`: read the columns; `Family` gets `tier` and `narrow_ok`; validation as spec 3.4; `COLUMNS` gets `hit_narrow` for `pfam_adhesion` only; `pfam_rows` writes it.
4. `modules/cli.py`: digest and `params` as spec 3.4.
5. `categories.yaml`: add `wall_family_narrow_domain`; `cell_wall_adhesion_candidate` reads it; `surface_attachment_candidate` stays on the broad call. This changes `config_sha256`.
6. `outputs.py` and descriptions: say which families feed which call. Regenerate `CLASSES.md`.
7. Regenerate golden files; read the diff line by line; the only changes are the new column, the new call and the changed composite.

Tests (write first):
1. SP protein with a hit in an `untested`, `narrow_ok=no` family: `surface_attachment_candidate` `called`, `cell_wall_adhesion_candidate` `not_called`. (Fails before.)
2. Same protein, family `smoke`, `narrow_ok=yes`: both `called`.
3. `load_family_table` refuses: missing new columns; tier not in the vocabulary; `narrow_ok=yes` with tier `untested`; `narrow_ok=yes` without `narrow_ok_by` or `narrow_ok_date`.
4. `hit_narrow == 1` implies `hit == 1` (generated rows).
5. `_family_digest` changes when `narrow_ok` or the tier changes.
6. An old `pfam_adhesion` table without `hit_narrow` stops the run with `missing column`.
7. `tests/class_table/test_class_table.py` still passes after regeneration.

Mutations: narrow flag ignores tier; ignores `narrow_ok`; digest ignores `narrow_ok`; narrow call reads the broad flag.

### Task 3.4: `derived` marker and R0 exception text (D6)

Files: `src/cellsurface_sorting_hat/engine.py` (`:334-349`), `status.py` (`:17-28`, `:196-201`), `calibration/phasec.py` (`:41`, `:250-251`, `:269`), `outputs.py` (`:47-48`, `:206`), `docs/CLASSES.md` (regenerated), tests `test_engine.py`, `calibration/test_phasec.py`, `test_outputs.py`, `golden/*`.

Steps:
1. Tests first.
2. Key `leakage_exception` (Q5) in `MEASURE_KEYS`; `phasec` writes it for R0 entries with the owner decision and date text. The loader (Task 3.1) skips the leakage cap for an entry with this key and no other.
3. Basis text: where an entry has the key, `status_basis` carries "leakage unmeasured (SignalP training overlap), uncapped by owner decision".
4. `derived;` prefix for records of calls that have a `ref` node or are kind `other`.
5. Regenerate golden files.

Tests (write first):
1. A composite record has `status_basis` starting with `derived;`. (Fails before.)
2. A record whose status comes from the R0 entry contains the exception text in the basis; a leaf `repeat` record does not.
3. The loader accepts an `estimated` R0 entry with the key and no `leakage:` value; it refuses the same entry without the key.
4. `phasec` output for the A. nidulans fixture has the key and `negative_basis` per Q5.

Mutations: marker on leaf calls; marker missing on `other_*`; loader exempts every entry.

### Task 3.5: tier derivation and measurement file (D3)

Files: `data/sorting_hat/family_measure.tsv` (new, PROPOSAL name), `analysis/family_tests/derive_tier.py` (new), `analysis/family_tests/README.md`, `tests/cellsurface_sorting_hat/test_family_tier.py` (new).

Steps:
1. Tests first (below).
2. The script imports `status_from_measure` from `cellsurface_sorting_hat.status`. Mapping: no row, no positive or no negative gives `untested`; `status_from_measure == estimated` gives `estimated`; otherwise `smoke`. `--check` compares to the table and exits 1 on a difference; `--write` rewrites the column.
3. Fill `family_measure.tsv` from `docs/reports/data/sorting_hat/family_tests/family_tests.tsv` where it has counts (Q4); rows without lists stay absent. Record the source file and date in each row. This gives counts, not intervals; rows without intervals derive `smoke` at most.
4. Run `--write`; commit the table.

Tests (write first):
1. The committed table equals the derived tiers (fails before: no tier column).
2. A copy of the table with one stronger tier fails the test.
3. Fixture rows: 25/25/22/22 with narrow intervals derive `estimated`; 5/5 derive `smoke`; no row derives `untested`.
4. The script refuses a row whose cluster counts are missing (derives at most `smoke`, with a message).

Mutations: script ignores clusters; script copies the thresholds instead of calling `status_from_measure`; `--write` changes other columns.

### Task 3.6: R0 snapshot check (D6)

Files: `src/cellsurface_sorting_hat/modules/cli.py` (`:428-459`), `cli.py` (`:366-372`, `:380-412`), `outputs.py` (`_warnings`, `RunInfo` at `:53-75`), `docs/TOOL-ARCHITECTURE.md` if it describes run records, tests `modules/test_module_cli.py`, `test_cli.py`, `test_outputs.py`.

Steps: record a digest of the `(id, value)` pairs in `_condition_table`; after `load_modules` recompute the digest from the engine's R0 table and compare with the digest in each `pfam_*` run record's `params.conditions`; mismatch gives a warning and a `run.json` field; a missing digest in an old record gives a note, not an error.
Tests (write first): build `pfam_adhesion` against R0 table A, run with R0 table B: warning present (fails before). Same table: no warning. Old record without digest: note. `run.json` carries the field.
Mutations: compare only `params_hash`; compare digests of ids only.

### Task 3.7: `primary_label` (D1)

Files: `src/cellsurface_sorting_hat/labels.py` (new), `outputs.py` (`:93-127`, `:276-285`), `cli.py` (call site near `:416-417`), `docs/class_descriptions.yaml`, `docs/CLASSES.md` (regenerated), tests `test_labels.py` (new), `test_outputs.py`, `golden/*`.

Steps: tests first; implement `primary_labels(records, cfg)` returning one label per protein for the default gate with the precedence of spec 3.1 (Q1); add the column to `write_wide`; add the "Primary label" count table (including `mechanism_no_class` and `not_assessable`) to the report.
Tests (write first): P1 (no SP, repeat02 called) gets `mechanism_no_class`; SP plus HsbA only gets `surface_attachment_candidate`; P2 gets `not_assessable`; an SP protein with only `sowgp_ortholog`-style mechanism (a test config that adds such a call) gets `mechanism_no_class`; every protein of a generated fixture gets exactly one label; report counts add up to the protein count.
Mutations: precedence reversed; fallback missing; label computed for a non-default variant.

### Task 3.8: documentation of the new rules

Files: `docs/paper/03-status-and-validation-rules.md`, `docs/CLASSES.md` (regenerated), `docs/class_descriptions.yaml`, `AGENTS.md` (one entry for this work, by the implementing agent).
Steps: state `active`, tier, `narrow_ok`, `not_applicable`, `derived`, `negative_basis`, `primary_label` once each; link to the spec. Test: the docs test of Task 2.1 plus the class-table test.

### Task 3.9: re-measure after the config edit

Files: `docs/reports/data/sorting_hat/call_status/*` (new files only; old files stay as history).
Steps: after Tasks 3.2 and 3.3 the existing call files are stale (`call_status.py:179-180`, `engine.py:497-508`). Re-run `calibrate truth --call-status` for `tandem_repeat_protein` in S288C and C. albicans on HPCC (as in M2d and M5b2); note whether the numbers changed; update the `CLASSES.md` citations. Stop point: owner reviews the new files. Not done by this plan's PR.
Tests: the committed new files load (Task 3.1 test 7).

## PR 4: tests

### Task 4.1: L1 shared control-set validator

Files: `analysis/_common/validate_controls.py` (new; `analysis/_common` exists), `tests/controls/test_validate_controls.py` (new), `tests/controls/fixtures/` (tiny good and bad sets), `.github/workflows/build_and_test.yml` (new test path).
Steps: tests first; validator checks COMMON-RULES columns (`docs/agent-tasks/COMMON-RULES.md`), `label` in `positive|negative`, `leakage` in the five values, `tuned` in `yes|no`, one row per accession, every accession in `sequences.faa` and `clusters.tsv`, no cluster mixes a positive and a negative unless the row says so (`mixed_cluster_with_positive` in the shared set), counts in `counts.md` equal counts from `controls.tsv`. A script entry `python analysis/_common/validate_controls.py data/controls/<dir>`.
Tests (write first): good fixture passes; each defect (missing column, bad leakage value, missing sequence, row missing from clusters, count mismatch, duplicate accession) fails with a message that names the row. Run the validator on `data/controls/hard-negatives-shared/` and `data/controls/repeat-mechanism/`; any real defect it finds is reported to the owner, not fixed silently.
Mutations: validator ignores `clusters.tsv`; ignores `counts.md`.

### Task 4.2: L3 golden fixture tests, including the gap cases (SH-D1)

Files: `tests/cellsurface_sorting_hat/golden/` (new fixture set `gaps/`), `tests/cellsurface_sorting_hat/test_golden_gaps.py` (new).
Steps: build module tables for six proteins: (1) no SP, repeat called; (2) SP, hydrophobin only; (3) SP, HsbA only; (4) SP, `pfam_hsba` error; (5) SP, `antigen_lookup` not applicable, all else not called; (6) SP, family hit in an untested family. Store the expected `calls.long` and `calls.wide` lines (including `primary_label`). The test runs the real command (`cellsurface_sorting_hat.cli`) on the fixture and compares.
Tests (write first): the six expected rows, written by hand from the spec rules before the engine is read again. Gap cases: protein 1 has no composite and no `other_*` label and `primary_label` is `mechanism_no_class`; protein 4 is `not_assessable`.
Mutations: swap two precedence entries; skip the narrow flag; read the broad flag in the narrow call.

### Task 4.3: docs-consistency test (L3)

Files: `tests/cellsurface_sorting_hat/test_docs_consistency.py` (extend Task 2.1's file).
Steps: tests that (a) every family of an active `pfam_adhesion` row with `narrow_ok=yes` is named in the generated `CLASSES.md` text for the narrow call, (b) every active `pfam_adhesion` row is named for the broad call, (c) `CLASSES.md` mentions each call of `categories.yaml`, (d) tier and derived tier agree (reuses Task 3.5).
Mutation: remove a family from the description; the test must fail.

## PR 5: baseline recall on the Kiss truth set (HPCC)

Stop point before: owner confirms Q9 and Q10 and the two hard-negative sources (spec section 9, item 6).

### Task 5.1: map the Kiss proteins to sequences (D8)

Files: `analysis/wall_fusion_truth/map_kiss.py` (new), `data/controls/wall-fusion-kiss2019/` (new, COMMON-RULES files).
Steps: read `docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` (`protein_id_jgi`, `cluster_number`, `node_of_emergence`). Map each JGI ID to a protein sequence in the repo proteomes or to UniProt. Count mapped and lost rows per clade and write `counts.md` and `unverified_notes.md`. Do not pad. Test: the L1 validator (Task 4.1) accepts the directory; a unit test on a 3-row fixture checks the mapping and the loss count.

### Task 5.2: cluster and split

Files: `analysis/wall_fusion_truth/cluster_split.sh` (HPCC, MMseqs2 30% identity, 50% coverage, as in `docs/superpowers/plans/2026-10-10-allergen-predictor-plan.md` T0.4), `analysis/wall_fusion_truth/seed_split.py` (new).
Steps: cluster the positives and the hard negatives; compute identity to the nearest Pfam seed of the models in `analysis/wall_fusion_report/wall_fusion_models.tsv` and of the active family rows (method of `analysis/lectin_validation/measure.py:328-329,402-409`); label each cluster `seed-like` or `seed-independent`. Job sized to 1 hour; `${SCRATCH:?}`; no `/tmp`.
Test: fixture with identities 35% and 20% lands in the two groups; a cluster with both kinds is `seed-like` (a cluster is seed-independent only when no member reaches 30%).

### Task 5.3: baseline recall and false-positive rate

Files: `analysis/wall_fusion_truth/baseline_recall.py` (new), `analysis/wall_fusion_truth/submit_baseline.sbatch`, `docs/reports/2026-10-XX-wall-fusion-baseline.md`.
Steps: run `hmmsearch --cut_ga` of the 21 models and the active family models on truth and hard negatives; report recall per cluster and per clade, and separately for `seed-like` and `seed-independent`; report hits on the hard negatives per model with cluster counts; label every number with its n; no pooled number; no tier is changed by this report.
Tests: recall on a fixture equals the hand count; the script writes no call and no status; a model with no GA line is reported, not skipped.
Output: the table of clusters Pfam misses. That table, not a guess, names the PR 6 candidates (D9 stage 2).

### Task 5.4: feed L2 and L5 (PROPOSAL)

Steps: write the family rows that now have positives and negatives into `data/sorting_hat/family_measure.tsv` (Task 3.5); run `derive_tier.py --check`. Pin the baseline numbers with a tolerance and the Pfam and HMMER versions (L5). Owner decides whether any `narrow_ok` follows (D4); the plan does not.

## PR 6: custom HMMs with leave-clade-out

Own spec before work starts, written from the Task 5.3 table. Outline only:
- 6.1 For each cluster Pfam misses, build an alignment of orthologs (MAFFT or the method of `2026-10-08-hydrophobin-custom-hmm-design.md`), `hmmbuild`, and set GA from the training clades only.
- 6.2 Leave-clade-out (`analysis/hydrophobin_truth/lco.py` is the precedent): build from some clades, test on held-out clades, report per held-out clade.
- 6.3 A custom model gets at most `smoke` until 6.2 is committed (the derive script refuses a higher tier for a row without a leave-clade-out measurement). No ESM model (issue #15; gate of 300 clusters).
- Tests: the derive script refuses `estimated` for a custom model without the leave-clade-out field; the held-out clade never appears in the build alignment (a test on the build manifest).

## Out of scope

- ESM or fine-tuned models (D9 stage 4, issue #15).
- A residual `other_no_surface_with_mechanism` call (D1 chose a derived label).
- Changes to `status_from_measure` thresholds.
- Wall/fusion model-list and label changes (WF-F8, WF-F9; spec Q9) unless the owner decides them.
- The NRP and GPU work of the model-review plan (issues #15, #16).

## Self-review against the spec

| Spec item | Task |
|---|---|
| D1 `primary_label`, wording | 2.1 (wording), 3.7, 4.2 |
| D2 tier column, `active` meaning, docs | 2.1, 3.3 |
| D3 derived tier, CI test | 3.5, 4.3 |
| D4 narrow call | 3.3, 4.2 |
| D5 `other_*` and `not_applicable` | 3.2, 4.2 |
| D6 R0 text, `derived`, snapshot check | 3.4, 3.6 |
| D7 loader and `negative_basis` | 3.1 |
| D8 truth set | 5.1, 5.2 |
| D9 model stages | 5.3 (stage 1), PR 6 (stages 2, 3) |
| D10 order | PR headings |
| Testing framework L0 to L5 | L0 2.3; L1 4.1; L2 3.5, 5.4; L3 4.2, 4.3; L4 not planned (PROPOSAL, issue #12); L5 5.4 |
| Risks and defaults | P4, Task 3.9 |
