# Sorting-hat review fixes and wall/fusion testing: design

*2026-10-17. Draft 1. Author: Claude Code (claude-sonnet-5-5). Status: **draft; the owner is still refining it.** Sources: two independent reviews of 2026-10-17 (`docs/reviews/2026-10-17-sorting-hat-logic-review.md`, `docs/reviews/2026-10-17-wall-fusion-products-review.md`; both by an Opus agent, not verified by the owner) and the owner's grilling session of 2026-10-17. Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`. No code, config or existing file is changed by this spec. Line numbers are for `main` at `83d5f60`; the reviews cite `c7eadf0`, and I re-read every line I cite.*

**How to read and edit this document.**
- Each section stands alone. Edit one section without touching the others.
- Three labels are used. **DECIDED** is an owner decision of 2026-10-17 (D1 to D10, section 3). **PROPOSAL** is my suggestion; the owner has not decided it. **OPEN** is a question that section 8 lists.
- Review findings are cited with a prefix so they do not clash with the owner decisions D1 to D10. `SH-D1` is finding D1 of the logic review. `SH-C5` is finding C5. `WF-F2` is finding F2 of the wall/fusion review.
- I added no decision beyond D1 to D10. Anything else is a PROPOSAL or an OPEN question.

## 1. Problem

The reviews read the code and the reports as they stand on `main`. The logic review ran `uv run pytest tests/cellsurface_sorting_hat`: 818 passed, 9 skipped, 2 failed (both failures are macOS-only: `realpath -m` is GNU). I did not re-run the tests. Findings, ranked by what they can do to a reader of a report:

| Rank | Finding | Review ids | What is wrong |
|---|---|---|---|
| 1 | The `other_*` calls do not close the class list | SH-D1 | A protein with no signal peptide and a mechanism hit gets no composite and no `other_*` label. An SP protein with only hydrophobin or HsbA gets `surface_attachment_candidate` but not the narrow call and not `other_surface_no_mechanism`. Two reports call the labels a "partition" (`docs/reports/2026-10-07-classifier-runs-and-scale.md:30,70`). The sowgp-like spec (M5) would widen the gap. |
| 2 | `active` does not mean "tested" | SH-D3, SH-D2, SH-C5, SH-C9, WF-F1 | 8 of the 18 active family rows have no passed specificity test (rows 8, 9, 12, 13, 21, 22, 23, 24 of `data/sorting_hat/family_table.tsv`). All 8 are in `pfam_adhesion`, so they set `wall_family_domain` and the narrow call. PIR1-like_C has 0 of 54 E1 positives and hits 3 N1 negatives (`docs/reports/2026-10-09-family-tests.md:32,38`). Docs still describe only CFEM and PA14. |
| 3 | The narrow call and the broad call read one flag | SH-C5 | Both composites read `pfam_adhesion.hit` (`categories.yaml:54,63,74`). A family cannot be limited to the broad call. |
| 4 | Unknown and not-applicable are mixed in `other_*` | SH-C1 | With an unknown mechanism, `other_surface_no_mechanism` can be `called` while the composite is `not_assessable`. `antigen_lookup` writes `not_applicable` for every non-Coccidioides protein (`modules/lookups.py:156-157`); the engine reads it as unknown. |
| 5 | The status loader trusts files; negatives have no quality field | SH-C7, SH-C4 | `load_status_source` (`status.py:136-164`) does not re-check status against the measure, the leakage cap or the taxon floor. The measure counts assumed negatives like curated ones (`status.py:105-133`). |
| 6 | Status wording hides two facts | SH-C2, SH-C3 | R0 is `estimated` although its leakage is unmeasured (owner exception, `docs/paper/03-status-and-validation-rules.md:91`). Composite negatives inherit `estimated`, and the planned `derived` marker was never added (`2026-10-09-v1-class-list-design.md:114`, O5). |
| 7 | Two R0 snapshots are never compared | SH-C6 | `pfam_adhesion` reads one R0 table when it is built (`modules/cli.py:265-269`). The engine gates on another. Nothing checks that they agree. |
| 8 | Wall/fusion reports carry wrong or unmarked statements | WF-F1, WF-F2, WF-F3, WF-F4, WF-F5 | Errata missed the status of PIR1-like_C and Hyr1; claims that report 07 contradicts are unmarked; "wall enzymes are the main false-positive risk" holds for the ESM classifier, not for the family route. |
| 9 | The wall/fusion script lost two checks of the lectin script | WF-F6, WF-F7 | No duplicate-ID warning; no run-time check of model names; the "verified" test list is hand-copied. |
| 10 | The wall/fusion model list and labels are inconsistent | WF-F8, WF-F9 | Different criteria for similar models; `model_label` names a function where the lectin script names the model. |
| 11 | Test coverage has holes | SH-D1 (probe), CI file | CI runs 3 test dirs plus 7 informational analysis dirs (`.github/workflows/build_and_test.yml:46,57`). Five test dirs never run in CI: `cys_candidates`, `lectin_class_report`, `lectin_validation`, `secretome_discovery`, `wall_fusion_report`. Two job-script tests fail on macOS. |
| 12 | Small items | SH-L1, SH-L2, SH-C8, SH-C9, WF-F10, WF-F11 | Raw `--notes` before `call=` parsing; complete-case sensitivity; thin evidence for some CLASSES.md claims; PRM1 score columns; PA14 negative tests nothing. |

Things the reviews found sound are listed in each review. In short: the three-valued logic, config validation, per-call status design, module wrappers and the c-i report scripts hold.

## 2. Decisions at a glance

| # | Decision (owner, 2026-10-17) | Label | Section |
|---|---|---|---|
| D1 | Calls stay independent flags. One derived column `primary_label`. Drop the word "partition". | DECIDED (precedence order OPEN) | 3.1 |
| D2 | Keep the 8 untested families active. Add `specificity_tier`. `active` means "enabled by owner". | DECIDED | 3.2 |
| D3 | Tier values `untested`, `smoke`, `estimated`, derived by a script from a committed measurement file. CI test. | DECIDED (vocabulary accepted) | 3.3 |
| D4 | Narrow call reads a family only if tier is `smoke` or better and `narrow_ok=yes`. | DECIDED (accepted) | 3.4 |
| D5 | Unknown mechanism makes `other_*` `not_assessable`. New `not_applicable` value. | DECIDED | 3.5 |
| D6 | R0 keeps its owner exception, made visible. `derived` marker. Compare the two R0 snapshots. | DECIDED | 3.6 |
| D7 | Loader re-validates. `negative_basis` field; `assumed` caps at `smoke`. | DECIDED | 3.7 |
| D8 | Truth set for wall remodeling and fusion genes: Kiss 2019 Data 3 and 4. | DECIDED | 3.8 |
| D9 | Models in stages: Pfam baseline, custom HMMs only where Pfam misses, leave-clade-out, no ESM now. | DECIDED | 3.9 |
| D10 | Six PRs in a fixed order. | DECIDED | 5 |

## 3. Decisions in detail

Each block gives the decision, the alternatives the owner considered, the reason, what changes and where, and the acceptance test. "Fails before" means the test fails on `main` and passes after the change.

### 3.1 D1: top-level labels

- **Decision (DECIDED).** The calls stay independent flags. Add one derived column `primary_label` with a documented precedence. Remove the word "partition" from `docs/reports/2026-10-07-classifier-runs-and-scale.md:30` and `:70`.
- **Alternatives.** (a) Add a residual call `other_no_surface_with_mechanism` so the `other_*` calls close the list (SH-D1 option 1). (b) Only reword and print the unlabeled count (SH-D1 option 2). (c) Make the calls a strict partition.
- **Reason.** Independent flags keep every call testable on its own and keep `categories.yaml` and call hashes unchanged. One derived column gives readers one answer without a new call.
- **Changes.**
  - New derived label code (PROPOSAL: a new file `src/cellsurface_sorting_hat/labels.py`; precedence in code, not in `categories.yaml`, so `config_sha256` does not change, see `engine.py:82-95` and `call_status.py:179-180`).
  - `outputs.py:110-127` (`write_wide`): add the column `primary_label` for the default gate (`categories.yaml:14`, R0). `write_long` (`outputs.py:93-98`) is not changed.
  - `outputs.py:199-208` and the "Calls" section (`outputs.py:276-285`): add a table "Primary label" with the count per label, including the fallback labels below.
  - `docs/class_descriptions.yaml` and `analysis/class_table/make_classes_md.py` (regenerate `docs/CLASSES.md`): one paragraph that says the calls are labels, not a partition, and name the fallback labels. `CLASSES.md:22` is accurate and stays.
  - Reports: reword `docs/reports/2026-10-07-classifier-runs-and-scale.md:30,70`.
- **PROPOSAL: precedence (first match wins; OPEN, section 8 Q1).**
  1. `cell_wall_adhesion_candidate` called.
  2. `surface_attachment_candidate` called.
  3. `serodiagnostic_marker_candidate` called.
  4. `other_surface_no_mechanism` called.
  5. `other_not_surface` called.
  6. `mechanism_no_class`: no label above and at least one mechanism call is `called`. This holds the two gap cases of SH-D1: no signal peptide with a mechanism hit, and a signal peptide with only a mechanism that no composite reads (the M5 `sowgp_ortholog` case).
  7. `not_assessable`: no label above and the gate or a deciding input is unknown.
  Allergen calls are outside the mechanism list (`categories.yaml:100,106`) and do not enter `primary_label` (OPEN, Q2).
- **Acceptance tests.**
  - Fails before: a three-protein fixture (the probe of the logic review). P1 = no SP and repeat02 called gets `mechanism_no_class`. There is no `primary_label` column on `main`.
  - Fails before: a docs test that greps `docs/reports/2026-10-07-classifier-runs-and-scale.md` for the word "partition" and expects no match.
  - Property test: every protein in a generated fixture gets exactly one `primary_label`.

### 3.2 D2: family evidence tier; what `active` means

- **Decision (DECIDED).** Keep the 8 untested families active (owner decision of 2026-10-09 stands). Add a column `specificity_tier` to `data/sorting_hat/family_table.tsv`. `active` means "enabled by owner", not "tested". Fix the docs that say otherwise.
- **Alternatives.** (a) Switch the 8 families off until tested. (b) Keep them active and say nothing. (c) Keep them active, but keep them out of the narrow call (this is D4, which the owner also took).
- **Reason.** The owner wants the families to feed the broader call now (decision of 2026-10-09). The label must say how well each family is known.
- **Changes.**
  - `data/sorting_hat/family_table.tsv`: new columns `specificity_tier`, `narrow_ok`, `narrow_ok_by`, `narrow_ok_date` (D3, D4). The file has 24 data rows and 11 columns today (`modules/pfam.py:9-21`).
  - `modules/pfam.py:9-21` (`FAMILY_COLUMNS`), `:35-42` (`Family`), `:45-80` (`load_family_table`): read and validate the new columns.
  - `modules/pfam.py:31-32`: the docstring of `NoActiveFamilyError` says "has passed its specificity test and sign-off". Reword to "no family of the module is enabled by the owner".
  - `docs/reports/cell-wall-fusion/04-repo-coverage-and-gaps.md:22` ("`active` means the family passed a specificity test and an owner sign-off") and `:227` ("`active` already promises a passed specificity test"). Mark both as wrong and say what `active` means now (WF-F1, SH-D3).
  - `docs/superpowers/specs/2026-10-10-lectin-class-design.md:266` (item 4) and `:270` (c-ii cell): same fix.
  - Also stale for SH-D2, to fix with D4: `docs/CLASSES.md:13`, `outputs.py:29-34,294,317`, `2026-10-04-orchestrator-design.md:191`, `2026-10-09-v1-class-list-design.md:26`.
  - Rows 8, 9, 12, 13 say "specificity test not run". The family test of 2026-10-09 exists (`docs/reports/2026-10-09-family-tests.md`). Update the notes from the measurement file (see Q4).
- **Acceptance tests.**
  - Fails before: `load_family_table` on a table without `specificity_tier` raises `FamilyTableError`; today it accepts the old table.
  - Fails before: a docs test that greps `pfam.py` and the three named docs for the phrase "passed its specificity test" (or "passes its specificity test") and expects none.

### 3.3 D3: tier vocabulary and derivation

- **Decision (DECIDED).** The tier reuses the status vocabulary (`status.py:14`): `untested`, `smoke`, `estimated`.
  - `smoke`: reviewed member and non-member lists exist but are below the size floors of `docs/agent-tasks/COMMON-RULES.md` item 7.
  - `estimated`: at least 20 positives and 20 negatives, at least 20 clusters each, 95% half-width at most 0.10 (`status.py:97-133`).
  - A script derives the tier from a committed measurement file. Nobody edits the tier by hand. A CI test fails if `family_table.tsv` disagrees with the derived value.
- **Alternatives.** A new vocabulary (for example `E0` to `E3`); a hand-edited tier column with a review rule.
- **Reason.** One vocabulary for module status, call status and family tier. A derived value cannot drift from its evidence.
- **Changes.**
  - New committed measurement file (PROPOSAL name `data/sorting_hat/family_measure.tsv`, one row per family: counts of positives, negatives, clusters, sensitivity and specificity with intervals, list sources, date, HMMER and Pfam versions). Produced by the L2 job (section 4).
  - New derive script (PROPOSAL `analysis/family_tests/derive_tier.py`). It must call `status_from_measure` (`status.py:105-133`), not copy its rule. Mapping: no row, no positive or no negative gives `untested`; `estimated` when `status_from_measure` says so; otherwise `smoke`.
  - New test `tests/cellsurface_sorting_hat/test_family_tier.py`: derive from the measurement file and compare to the table. CI already runs `tests/cellsurface_sorting_hat` (`build_and_test.yml:46`).
  - The existing descriptive test (`analysis/family_tests/family_tests.py`, `docs/reports/data/sorting_hat/family_tests/family_tests.tsv`) gives E1 and N1 hit counts, not intervals (report header, line 3). It is the first input to the measurement file; it does not fill it by itself (Q4).
- **Acceptance tests.**
  - Fails before: edit one tier in a copy of the table to a stronger value; the test fails. No such test exists.
  - Fails before: a measurement row with 25 positives, 25 negatives, 22 and 22 clusters and narrow intervals must derive `estimated`; a row with 5 and 5 must derive `smoke`; a family with no row must derive `untested`.
  - Mutation: the script ignores the cluster counts; a named test must fail.

### 3.4 D4: the narrow call reads a family only with evidence and sign-off

- **Decision (DECIDED, accepted).** `cell_wall_adhesion_candidate` reads a family only if `specificity_tier` is `smoke` or better and `narrow_ok=yes`. `narrow_ok` is an owner-signed column with who and when, like `active_by` and `active_date`. Untested families still feed the broader `surface_attachment_candidate`.
- **Alternatives.** Leave both composites on one flag (today); switch the families off (D2 option a); split by module name.
- **Reason.** SH-C5: structural wall proteins (Pir-like, CFEM, Hyr1) reach a call named "adhesion" although their own data show non-adhesive hits. The engine cannot restrict a family to the broad call today.
- **Changes.**
  - `modules/pfam.py:24` (`COLUMNS`): add `hit_narrow` (PROPOSAL: written for `pfam_adhesion` only). `pfam_rows` (`:125-171`): a family counts for `hit_narrow` only when `tier` is `smoke` or `estimated` and `narrow_ok` is `yes`. `hit_narrow` implies `hit`.
  - `modules/pfam.py:45-80`: reject `narrow_ok=yes` with tier `untested`, or without `narrow_ok_by` and `narrow_ok_date` (the same pattern as `:68-69`).
  - `modules/cli.py:225-232` (`_family_digest`): add tier and `narrow_ok` to the digest, because they now decide a call. `modules/cli.py:278-283`: add the narrow family list to `params`.
  - `categories.yaml`: keep `wall_family_domain` (`:53-54`, `pfam_adhesion.hit`). Add an evidence call (PROPOSAL name `wall_family_narrow_domain`, `pfam_adhesion.hit_narrow`). `cell_wall_adhesion_candidate` (`:59-64`) reads the narrow call. `surface_attachment_candidate` (`:69-75`) keeps the broad call; add the narrow call to `basis_calls` (`:71`) only if the owner wants it shown (OPEN, Q6). The `mechanism` lists (`:100,106`) keep the broad call: narrow implies broad.
  - `outputs.py:29-34` (known limit 4), `:294,317` (label "wall family domain (PA14, CFEM)"), `docs/class_descriptions.yaml`: say which families feed which call.
  - `cli.py:154-181` (`required_columns`) reads the new flag, so an old `pfam_adhesion` table without `hit_narrow` stops the run with `missing column` (`cli.py:214-216`). Pfam modules must be rebuilt.
- **Acceptance tests.**
  - Fails before: a protein with an SP and a hit in an `untested`, `narrow_ok=no` family: `surface_attachment_candidate` is `called` and `cell_wall_adhesion_candidate` is `not_called`. On `main` both are `called`.
  - A protein with a hit in a `smoke`, `narrow_ok=yes` family: both are `called`.
  - `load_family_table` refuses `narrow_ok=yes` with tier `untested`.
  - Property test: `hit_narrow == 1` implies `hit == 1` for every row.

### 3.5 D5: unknown mechanism input and `not_applicable`

- **Decision (DECIDED).** An unknown mechanism input makes the `other_*` calls `not_assessable`. Today the engine skips the input, makes the call, and lists the input in `left_out` (`engine.py:243-253`). Add a separate value `not_applicable` for species-limited calls such as the Coccidioides antigen rank. It does not count as unknown.
- **Alternatives.** Keep today's rule and report the count of `other_*` with a non-empty basis (SH-C1 suggestion); mark such records `partial`.
- **Reason.** SH-C1: the probe protein P2 gets `other_surface_no_mechanism = called` while `surface_attachment_candidate = not_assessable`.
- **Changes.**
  - `logic.py:7-10`: add `NOT_APPLICABLE` to the value set. PROPOSAL for the Kleene rules (`logic.py:19-45`): `not_applicable` acts as `not_called` in `and`, `or` and `not`, because the call cannot be true for that protein (OPEN, Q3).
  - `engine.py:180-187` (`_row`) and `:190-215` (`_eval`): a row in state `not_applicable` gives `NOT_APPLICABLE`. Today every non-`ok` state gives `None`, which is unknown.
  - `engine.py:236-253` (`_eval_other`): a called mechanism still makes the result `not_called` (`:245-252`). Otherwise, an unknown mechanism gives `NOT_ASSESSABLE`. A `not_applicable` mechanism is ignored. PROPOSAL: `other_basis` then lists the unknown mechanisms on the `not_assessable` record, so the reader sees why.
  - `engine.py:30`: bump `ENGINE_SEMANTICS` (it is in every `call_hash`, `engine.py:497-508`). Every call status file becomes stale.
  - `outputs.py:276-289`: add a `not_applicable` column to the call counts. The section on `other_basis` changes meaning.
  - `calibration/cli.py:473-488`: `truth` handles only `not_assessable` (`:480`). A `not_applicable` call value would fall into `called == False`, a silent negative. Handle it on purpose (count as not assessable, report it).
  - `categories.yaml:11` comment and `docs/CLASSES.md`: describe the new value.
- **Acceptance tests.**
  - Fails before: the probe P2 (SP called, `pfam_hsba` row `error`): `other_surface_no_mechanism` must be `not_assessable`. On `main` it is `called`.
  - Probe P3 (`antigen_lookup` `not_applicable`, all else `not_called`): `cocci_specificity_rank_top15` is `not_applicable`, `other_surface_no_mechanism` is `called`, `other_basis` is empty.
  - A `truth` run on a table with `not_applicable` values does not count them as negatives.

### 3.6 D6: R0 exception made visible; `derived` marker; R0 snapshot check

- **Decision (DECIDED).** R0 (the Phase C SignalP call) keeps its owner exception: no leakage cap, status `estimated` (`docs/paper/03-status-and-validation-rules.md:91`). Make it visible. Wherever R0 `estimated` appears, the basis says "leakage unmeasured (SignalP training overlap), uncapped by owner decision". Composite calls that inherit a status get the planned `derived` marker (`2026-10-09-v1-class-list-design.md:114`, O5). Add a check that compares the two R0 snapshots (`modules/cli.py:265-269,458-459`).
- **Alternatives.** Cap R0 at `smoke` like every other call (SH-C2 option 2); leave the wording as is.
- **Reason.** SH-C2, SH-C3, SH-C6.
- **Changes.**
  - Exception text. Phase C writes it in the measure notes (`calibration/phasec.py:41,250-251`). The status basis is built in `status.py:196-201` (`taxon:N`) and in the engine (`engine.py:341-347`). PROPOSAL: a new optional measure key `leakage_exception` (text with who and when) in `MEASURE_KEYS` (`status.py:17-28`). The basis then carries the text. The loader (D7) uses the same key to skip the leakage cap for this entry only. Without a key the exception cannot be told from an error in the loader (Q5).
  - Derived marker. PROPOSAL: in `evaluate` (`engine.py:334-349`), a record of a call with a `ref` node or of kind `other` gets `status_basis` that starts with `derived;`. `outputs.py:206` keeps its general sentence. `docs/CLASSES.md` status column already says "derived".
  - R0 snapshot check. PROPOSAL: `_condition_table` (`modules/cli.py:428-459`) also records a digest of the `(id, value)` pairs it read, next to `params_hash` and `artefact_hash`. After `load_modules` (`cli.py:366`), the engine recomputes the digest from the R0 table it loaded and compares it to the digest inside each `pfam_*` run record. A mismatch is a warning in `report.md` (`outputs.py` `_warnings`) and a field in `run.json`. The identity hashes alone are not enough: the SignalP identity is tool, version and mode only (`modules/cli.py:240-247`).
  - `outputs.py:47-48` (limit 9): keep; add the same words as the basis.
- **Acceptance tests.**
  - Fails before: build `pfam_adhesion` against R0 table A, then run the engine with R0 table B; the report must show the warning. On `main` nothing fires.
  - Fails before: a record of `cell_wall_adhesion_candidate` has a `status_basis` that starts with `derived;`.
  - Fails before: a record whose status comes from R0 `estimated` has the exception text in its basis.

### 3.7 D7: status loader re-validates; `negative_basis`

- **Decision (DECIDED).** `load_status_source` (`status.py:136-164`) re-validates on load. It recomputes the status from the embedded measure, the leakage cap and the taxon rules. It rejects a file whose status is higher than the derived one. A lower status is allowed. Add `negative_basis` (`curated` or `assumed`) to the measure. An `assumed` measure is capped at `smoke`. A file without the field defaults to `assumed`.
- **Alternatives.** Trust files written by the tool (today); a signed hash in each file.
- **Reason.** SH-C7: hand-edited or old files can claim `estimated` for taxon 1. SH-C4: assumed negatives measure absence of annotation, not specificity. The call-file loader already does the first part (`call_status.py:66-102`).
- **Changes.**
  - `status.py:17-28` (`MEASURE_KEYS`): add `negative_basis`. `:74-94` (`_check_measure`): accept only `curated` or `assumed`.
  - `status.py:105-133` (`status_from_measure`): cap at `smoke` when `negative_basis` is not `curated` (a missing key counts as `assumed`). One place, so `measure.py:63,84-88` and `call_status.py:91` get the rule too.
  - `status.py:136-164`: after `_check_measure` (`:155`) compute the derived status. Apply the leakage cap from `measure.notes` with the regex of `call_status.py:37-38,61-63`, unless the entry has the exception key of D6. Refuse taxa below `MIN_TAXON` (`call_status.py:39`; move the constant to `status.py`). Refuse when the stored status is stronger than the derived one. A rank check needs a taxonomy and stays at write time (`calibration/cli.py:135-180`).
  - `calibration/measure.py:27-42` (`build_measure`) and `calibration/cli.py:92-97,442-529` (`truth`): a required option `--negatives {curated,assumed}`, as `--leakage` is required. `calibration/phasec.py:269` sets it for the R0 entries (Q5).
  - `docs/paper/03-status-and-validation-rules.md`: add the rule.
- **Acceptance tests.**
  - Fails before: a status file with taxon 1 and status `estimated` loads on `main`; it must raise.
  - Fails before: a file that says `estimated` with a measure that derives `smoke` must raise.
  - A measure with 25 positives, 25 negatives, narrow intervals and `negative_basis: assumed` derives `smoke`; the same with `curated` derives `estimated`; the same without the key derives `smoke`.
  - A file with a lower status than the derived one loads.

### 3.8 D8: truth set for recovery of wall remodeling and fusion genes

- **Decision (DECIDED).** The truth set is Kiss et al. 2019 (PMID 31501435) Supplementary Data 3 and 4, as reconciled in `docs/reports/cell-wall-fusion/06-kiss2019-reconciliation.md` and `06-kiss2019-hypha-genes.tsv` (182 lines with the header). Cluster at 30% identity and 50% coverage. Report recall per cluster and per clade. Split into `seed-like` (30% identity or more to a Pfam seed) and `seed-independent`, as in `analysis/lectin_validation/measure.py:328-329,402-409`. Pair with the existing hard-negative sets: `data/curated/adhesins/hard_negative_seeds.tsv` and the shared set of 161 rows (`data/controls/hard-negatives-shared/controls.tsv`, 162 lines with the header).
- **Known weakness.** The set is thin outside hyphal clades. A second source (UniProt reviewed entries, literature rows) may be needed for Basidiomycota and non-hyphal clades.
- **Alternatives.** UniProt reviewed entries alone; the Pfam seed members; the S. cerevisiae set of report 07 (23 proteins).
- **Reason.** It is the only curated, clade-spread list in the repo with reconciled counts (23 clusters, 45 curated genes in 06).
- **Changes.** None in `src/`. New data and a script (PROPOSAL): `data/controls/wall-fusion-kiss2019/` with the files that `docs/agent-tasks/COMMON-RULES.md` asks for (`controls.tsv`, `sequences.faa`, `clusters.tsv`, `counts.md`, `rejected.tsv`, `unverified_notes.md`, `README.md`); a recall script under `analysis/`. The hard-negative rows keep their own `leakage` and `tuned` columns.
- **Acceptance tests.** The L1 validator (section 4) accepts the new directory. A test on a three-row fixture shows recall per cluster and per clade, and puts a protein with 35% identity to a seed in `seed-like` and one with 20% in `seed-independent`.

### 3.9 D9: model approach, in stages

- **Decision (DECIDED).**
  1. Baseline: the current Pfam models (the 21 models of `analysis/wall_fusion_report/wall_fusion_models.tsv` and the models of the active family-table rows) on the D8 truth set.
  2. Custom profile HMMs (HMMER, from aligned ortholog sets) only for clusters that Pfam misses. Candidates: SOFT/HAM-2, HAM-5 and other fusion and signalling genes. The baseline decides; I do not guess.
  3. Any custom model is evaluated leave-clade-out (build from some clades, test on held-out clades). No custom model gets a tier above `smoke` before that.
  4. No embedding (ESM) models now. That is issue #15 with its own gate of 300 clusters.
- **Alternatives.** Embeddings first; custom HMMs for every cluster; fine-tuning a PLM.
- **Reason.** Pfam is cheap and already measured. Custom models and embeddings need evidence that Pfam fails.
- **Precedent.** Custom HMM design: `2026-10-08-hydrophobin-custom-hmm-design.md`. Leave-clade-out code: `analysis/hydrophobin_truth/lco.py` (test: `tests/hydrophobin_truth/test_lco.py`).
- **Changes.** None now. Stage 1 is a job on HPCC (section 5, PR 5). Stage 2 and 3 are PR 6 and have their own spec when the baseline is in.
- **Acceptance tests.** Stage 3: a test that a custom model without a leave-clade-out measurement cannot derive a tier above `smoke` (the D3 script refuses it). Stage 1: the recall table reproduces the counts of `analysis/wall_fusion_report/results/wall_fusion_counts.tsv` for the shared models.

## 4. Testing framework

Six layers. Only the D3 vocabulary and D4 are accepted by the owner. The rest are PROPOSALS awaiting the owner.

| Layer | What it checks | Where it runs | Label |
|---|---|---|---|
| L0 code and CI | Unit tests, lint. Add CI jobs for the 5 test dirs CI skips, and make the 2 macOS job-script tests portable. | GitHub CI | PROPOSAL (also PR 2 of D10) |
| L1 control-set integrity validator | One shared validator for every `data/controls/<task>/` set: COMMON-RULES columns, leakage values, sequences present, cluster file covers every row, counts match `counts.md`. | CI | PROPOSAL (D10 names it in PR 4) |
| L2 family specificity measurement | For each family: positives, negatives, clusters, intervals on the control sets. Writes the committed measurement file. The tier is derived from it (D3). | HPCC | D3 vocabulary and derivation DECIDED; the job is PROPOSAL |
| L3 golden pipeline fixture tests and docs-consistency tests | Small FASTA and module tables through the engine to expected calls, including the gap cases of SH-D1. `CLASSES.md` agrees with `categories.yaml` and the family table. | CI | PROPOSAL (golden tests named in D10 PR 4) |
| L4 T0-T6 benchmark harness | Issue #12, per mechanism class and per clade, never one pooled number (`docs/model-review/2026-09-27-review-and-framework-plan.md`, section 12 and P2). | HPCC | PROPOSAL |
| L5 pinned measurements | Committed measurements with a tolerance. Re-run when Pfam or HMMER changes. | HPCC | PROPOSAL |

- **L1.** Today each control set has its own checks (for example `analysis/lectin_validation/controls/`). No shared validator exists. The D8 directory is the first user.
- **L3 docs-consistency.** `tests/class_table/test_class_table.py` checks that `CLASSES.md` is generated from `class_descriptions.yaml` and `categories.yaml`; it does not read the family table (SH-D2). The new test checks that every family in `wall_family_domain` and `wall_family_narrow_domain` is named in the generated text, and that the tier column agrees with D3.
- **L3 golden.** `tests/cellsurface_sorting_hat/golden/` holds `calls.long.expected.tsv` and `report.expected.md`. D1, D5 and D6 change them. Read the diff line by line before accepting it.
- **L5.** A pinned measurement stores the value, the tolerance, and the versions that made it. A re-run outside the tolerance fails the check and asks for review; it never edits the pin.

## 5. Order of work (D10)

Each item is a separate PR. The plan has the tasks.

| PR | Content | Depends on |
|---|---|---|
| 1 | This spec and plan (with the two reviews). | none |
| 2 | Cheap fixes: stale docs (D2, WF-F1, WF-F2, D1 wording), the 2 macOS test failures, CI jobs for the 5 skipped test dirs. | PR 1 |
| 3 | Engine changes: `primary_label` (D1), `not_applicable` and `other_*` rule (D5), `derived` marker and R0 snapshot check (D6), loader re-validation and `negative_basis` (D7), family-table `specificity_tier` and `narrow_ok` (D2, D3, D4). | PR 2 |
| 4 | Tests: L1 shared validator, L3 golden fixture tests including the gap cases of SH-D1. | PR 3 |
| 5 | Baseline recall on the Kiss truth set (D8, D9 stage 1), on HPCC. | PR 4 (L1) |
| 6 | Custom HMMs with leave-clade-out (D9 stages 2 and 3). | PR 5 |

PR 3 is large. PROPOSAL: it may be split into 3a (D1, D5, D6 marker) and 3b (D2, D3, D4, D7) if the owner wants smaller reviews. The split is not a decision.

## 6. Risk and compatibility

| Change | Who breaks | Default and migration |
|---|---|---|
| `family_table.tsv` gets 4 columns (D2, D4) | Any reader of the table. `load_family_table` requires all of `FAMILY_COLUMNS` (`pfam.py:49-51`). The Pfam command reads the table (`modules/cli.py:253`). | The loader requires the new columns; the PR rewrites the table in the same commit. Initial `narrow_ok` for CFEM, PA14 and the other rows is an owner input (Q7). With all `no`, the narrow call stops using wall families at once. |
| `pfam_adhesion` gets `hit_narrow` (D4) | Old module tables: the run stops with `missing column` (`cli.py:214-216`). | Rebuild the four Pfam modules (HPCC job scripts in `scripts/sorting_hat/`). Module identity changes (`_family_digest`), so module status files go stale. |
| `categories.yaml` edit (D4) | `config_sha256` changes, so every call status file is stale (`call_status.py:179-180`). Today there are call files for `tandem_repeat_protein` in S288C and C. albicans. | Re-measure after PR 3 (as in M5b2 of the sowgp-like spec). A stale file never raises a status (`call_status.py:177-196`). |
| `ENGINE_SEMANTICS` bump (D5) | All call hashes change (`engine.py:497-508`). | Same re-measure. |
| New call value `not_applicable` (D5) | Readers of `calls.long`/`calls.wide`; `calibration/cli.py:480`. | Document in `CLASSES.md`. Handle in `truth` (3.5). |
| Status-file schema (D6, D7): `negative_basis`, `leakage_exception` | Module status files without the field. | Default `assumed`, so `estimated` becomes `smoke` on load. The R0 Phase C files must be rewritten by `calibrate phasec` with `curated` and the exception key (Q5). Without the rewrite, R0 in A. nidulans drops from `estimated` to `smoke` and `CLASSES.md:11` becomes wrong. |
| Loader refusal (D7) | A hand-edited or old file with a status above its measure. | Loader raises with the file name. Fix the file; do not relax the loader. Entries without a `measure` are OPEN (Q8). |
| Report and wide-table columns (D1, D5, D6) | Scripts that read `calls.wide` by position; the golden files. | `primary_label` is appended. `calls.long` gains no column. `status_basis` text changes for composites. Regenerate the golden files and read the diff. |
| Wall/fusion `model_label` wording (WF-F9, PROPOSAL) | Readers of `analysis/wall_fusion_report/results/*.tsv`. | Not part of D1 to D10; Q9. |

## 7. Proposals not decided

These came from the reviews. The owner has not decided them. They are cheap and I list them for PR 2 or for later.

1. WF-F2, WF-F1, WF-F4, WF-F5: mark the contradicted claims in reports 01 to 04, fix the commit attribution for Flocculin, close the Gel4 doubt as "annotation inconsistent with its own cluster", log replacements in the errata.
2. WF-F6: bring the duplicate-ID warning of `lectin_report.py:323-339,741-748,771-778` into `wall_fusion_report.py` (the review cites `:523-537`, `:256-266`).
3. WF-F7: compare the model NAME of Pfam to the TSV `name` at run time; derive the "verified" list from report 05 instead of a hand-copied dict.
4. WF-F10: close the PRM1 score question in the README (`hmmscan` domain score 475.5 versus `hmmsearch` sequence score 475.9).
5. SH-L1: apply `_plain` to `args.notes` (`calibration/cli.py:252-254,499-502`). SH-L2: put the pessimistic sensitivity bound in the measure, not only in the notes.
6. SH-C8: commit the other five hydrophobin status files or soften "smoke in 8 proteomes" in `CLASSES.md:14`.

## 8. Open questions

| # | Question | My default (PROPOSAL) |
|---|---|---|
| Q1 | Precedence of `primary_label` (3.1). Which order? Should it be per step 1 variant or for the default gate only? | The order in 3.1; default gate only; precedence in code, not in `categories.yaml`. |
| Q2 | Do allergen calls enter `primary_label`? | No. They are outside the mechanism list. |
| Q3 | How does `not_applicable` act inside `and`/`or`? As `not_called`, `serodiagnostic_marker_candidate` becomes `not_called` (not `not_assessable`) for every non-Coccidioides protein. Is that wanted? | Yes: the rank cannot be true outside C. immitis RS. |
| Q4 | Does the family test of 2026-10-09 (E1: 54 positives; N1: 76 clusters) count as "reviewed member and non-member lists" for `smoke`? If yes, the Als, Flo11 and GLEYA rows are `smoke` now, and their notes "specificity test not run" are wrong. | Yes for rows that have E1 and N1 counts; confirm per row. PIR1-like_C is `smoke` with a bad result. |
| Q5 | Is the Phase C negative set `curated` or `assumed` (D7)? How do R0 files carry the owner exception (new `leakage_exception` key)? D6 and D7 conflict without an answer: the R0 `estimated` entry has no `leakage:` value in its notes (`phasec.py:41`). | `curated` (GO-based) and the new key; the owner confirms. |
| Q6 | Is the narrow call listed in `basis_calls` of the broad composite? | No; the broad basis lists broad evidence. |
| Q7 | Initial `narrow_ok` values (who signs, which rows). Candidates: CFEM has 12 N1 hits in 9 of 76 clusters and 1 of 54 E1 positives; PA14 passes in 4 S288C proteins. | All `no` until the owner signs row by row. |
| Q8 | Entries without `measure`: the loader (D7) cannot derive a status. Cap at `smoke`, or refuse? | Cap at `smoke`. |
| Q9 | WF-F8, WF-F9: which criterion for the 21-model list, and which `label` text? Not in D1 to D10. | Decide in PR 5, when the baseline uses the list. |
| Q10 | Sequences for the Kiss proteins: the TSV has JGI protein IDs (`protein_id_jgi`). Can they be mapped to proteomes in the repo? How many are lost? | Count first (PR 5, first task); report the loss. |
| Q11 | Is the `leakage` cap on a family measure needed (family models are built from Pfam seeds, which may overlap the truth set)? | Yes: `seed-like` rows carry `in_reference`, as in `analysis/lectin_validation`. |

## 9. Conflicts between the decisions and the code as it is

1. D6 and D7: the R0 entry is `estimated` with unmeasured leakage and no `leakage:` value (`calibration/phasec.py:269`, no `cap=`; `call_status.py:96-100` would refuse it in a call file). A loader that enforces the leakage cap would reject the R0 file. D7's default `assumed` would also cap it at `smoke`. Both need Q5.
2. D5: `calibration/cli.py:480` treats only `not_assessable` as unknown. A new value would be a silent negative unless `truth` is changed in the same PR.
3. D4: the narrow call reads CFEM and PA14 today. With `narrow_ok` empty the narrow call loses them at once. This changes the counts in `docs/reports/2026-10-07-classifier-runs-and-scale.md:31` (6 to 15 candidates per proteome).
4. D2: rows 8, 9, 12, 13 say "specificity test not run", but the family test of 2026-10-09 has counts for them. The notes are older than the test.
5. D10 PR 2 says "CI jobs for the 5 test dirs CI skips". The five are `cys_candidates`, `lectin_class_report`, `lectin_validation`, `secretome_discovery`, `wall_fusion_report`. `tests/input_tests` is data, not tests. The `analysis-tests` matrix is `continue-on-error` (`build_and_test.yml:51`); "informational" jobs cannot fail a PR. Whether the new jobs must block is OPEN (plan, Task 2.3).
6. D8: the hard-negative set of 161 rows is `data/controls/hard-negatives-shared/`. The file `data/curated/adhesins/hard_negative_seeds.tsv` is a smaller seed file (32 lines). Both are named; the owner should confirm that both are meant.

## 10. Not in this spec

- No change to `src/`, `categories.yaml`, `family_table.tsv`, status files or CI by this PR.
- No ESM or fine-tuned model (D9 stage 4, issue #15).
- No new call beyond the evidence call of D4.
- No change to the interval or `estimated` rules (`status.py:97-133`).
- The family-table Pfam lists of the wall/fusion report (21 models) are not changed.

## 11. Sources

Reviews: `docs/reviews/2026-10-17-sorting-hat-logic-review.md`, `docs/reviews/2026-10-17-wall-fusion-products-review.md`. Specs followed for form: `2026-10-07-per-call-status-design.md`, `2026-10-09-sowgp-like-calls-design.md`, `2026-10-10-lectin-class-design.md`. Code read: `engine.py`, `logic.py`, `status.py`, `call_status.py`, `outputs.py`, `cli.py`, `categories.yaml`, `modules/{pfam,cli,lookups,base}.py`, `calibration/{cli,measure,phasec}.py`. Data read: `data/sorting_hat/family_table.tsv`. Docs read: `docs/CLASSES.md`, `docs/paper/03-status-and-validation-rules.md`, `docs/agent-tasks/COMMON-RULES.md`, `docs/reports/2026-10-09-family-tests.md`, `docs/reports/cell-wall-fusion/04`, `06`. CI: `.github/workflows/build_and_test.yml`.
