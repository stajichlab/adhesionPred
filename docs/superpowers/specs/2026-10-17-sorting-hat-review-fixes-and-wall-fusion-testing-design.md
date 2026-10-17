# Sorting-hat review fixes and wall/fusion testing: design

*2026-10-17. Draft 2. Author: Claude Code (claude-sonnet-5-5). Status: **owner decisions D1 to D22 are confirmed (2026-10-17); names and some numbers are still provisional (section 8).** Sources: two independent reviews of 2026-10-17 (`docs/reviews/2026-10-17-sorting-hat-logic-review.md`, `docs/reviews/2026-10-17-wall-fusion-products-review.md`), an independent review of D1 to D18 (`docs/reviews/2026-10-17-decisions-d1-d18-review-fable.md`; all three by AI agents and not verified by the owner), and the owner's grilling session. Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`. No code, config or existing file is changed by this spec. Line numbers are for `main` at `83d5f60`; I re-read every line I cite.*

**How to read and edit this document.**
- Each section stands alone. Edit one section without touching the others.
- Three labels are used. **Owner decision 2026-10-17** is a decision (D1 to D22, section 3). **PROPOSAL** is my suggestion; the owner has not decided it. **OPEN** is a question that section 8 lists.
- Finding prefixes: `SH-D1` is finding D1 of the logic review. `WF-F2` is finding F2 of the wall/fusion review. `FD-A3` is finding A3 of the Fable review (A1 to A15, B1 to B8).
- Where a decision was changed after the Fable review, the block says "Revised 2026-10-17 after the Fable review (FD-xx)" and keeps the old text as "Originally: ...".
- I added no decision beyond D1 to D22. Names of columns and keys are PROVISIONAL (section 8).

## 1. Problem

The logic review ran `uv run pytest tests/cellsurface_sorting_hat`: 818 passed, 9 skipped, 2 failed (both macOS-only: `realpath -m` is GNU). I did not re-run the tests. Findings, ranked by what they can do to a reader of a report:

| Rank | Finding | Review ids | What is wrong |
|---|---|---|---|
| 1 | The `other_*` calls do not close the class list | SH-D1 | A protein with no signal peptide and a mechanism hit gets no composite and no `other_*` label. An SP protein with only hydrophobin or HsbA gets `surface_attachment_candidate` but not the narrow call and not `other_surface_no_mechanism`. Two reports call the labels a "partition" (`docs/reports/2026-10-07-classifier-runs-and-scale.md:30,70`). |
| 2 | `active` does not mean "tested" | SH-D3, SH-D2, SH-C5, SH-C9, WF-F1 | 8 of the 18 active family rows have no passed specificity test. All 8 are in `pfam_adhesion`, so they set `wall_family_domain` and the narrow call. PIR1-like_C has 0 of 54 E1 positives and hits 3 N1 negatives (`docs/reports/2026-10-09-family-tests.md:32,38`). |
| 3 | The narrow call and the broad call read one flag | SH-C5 | Both composites read `pfam_adhesion.hit` (`categories.yaml:54,63,74`). |
| 4 | Unknown and not-applicable are mixed in `other_*` | SH-C1 | With an unknown mechanism, `other_surface_no_mechanism` can be `called` while the composite is `not_assessable`. `antigen_lookup` writes `not_applicable` for every non-Coccidioides protein (`modules/lookups.py:156-157`); the engine reads it as unknown (`engine.py:185-187`). |
| 5 | The status loader trusts files; negatives have no quality field | SH-C7, SH-C4 | `load_status_source` (`status.py:136-164`) does not re-check status against the measure, the leakage cap or the taxon floor. The measure counts assumed negatives like curated ones. |
| 6 | Status wording hides two facts | SH-C2, SH-C3 | R0 is `estimated` although its leakage is unmeasured (owner exception, `docs/paper/03-status-and-validation-rules.md:91`). Composite negatives inherit `estimated`; the planned `derived` marker was never added. |
| 7 | Two R0 snapshots are never compared | SH-C6 | `pfam_adhesion` reads one R0 table when it is built (`modules/cli.py:265-269`). The engine gates on another. |
| 8 | Wall/fusion reports carry wrong or unmarked statements | WF-F1 to WF-F5 | Errata missed the status of PIR1-like_C and Hyr1; contradicted claims are unmarked. |
| 9 | The wall/fusion script lost two checks of the lectin script | WF-F6, WF-F7 | No duplicate-ID warning; no run-time check of model names. |
| 10 | The wall/fusion model list and labels are inconsistent | WF-F8, WF-F9 | Different criteria for similar models; `model_label` names a function. |
| 11 | Test coverage has holes | SH-D1 (probe), CI file | CI blocks 3 test dirs and runs 7 analysis dirs as informational (`.github/workflows/build_and_test.yml:46,51,57`). Five test dirs never run in CI. |
| 12 | Small items | SH-L1, SH-L2, SH-C8, SH-C9, WF-F10, WF-F11 | See section 7. |

The Fable review (FD) found further faults in the first set of decisions; section 3 shows each fix in place.

## 2. Decisions at a glance

All are **Owner decision 2026-10-17**. "Revised" means changed after the Fable review.

| # | Decision | Section |
|---|---|---|
| D1 | Calls stay independent flags. One derived column `primary_label`. Drop "partition". | 3.1 |
| D2 | The 8 untested families stay active. Add `specificity_tier`. `active` = enabled by owner. | 3.2 |
| D3 | Revised. Tier (`untested`/`smoke`/`estimated`) derived by code in `src/` from a committed measurement file with `measure_basis`. | 3.3 |
| D4 | Narrow call reads a family only if tier is `smoke` or better and `narrow_ok=yes`. Implementation notes added. | 3.4 |
| D5 | Revised. Unknown mechanism makes `other_*` `not_assessable`. `not_applicable` value (see D18). | 3.5 |
| D6 | R0 keeps its owner exception, made visible. `derived` marker. Compare the two R0 snapshots. | 3.6 |
| D7 | Loader re-validates. `negative_basis`. Files without it default to `assumed`. | 3.7 |
| D8 | Revised. Truth set Kiss 2019 Data 3 and 4, after a feasibility stage. | 3.8 |
| D9 | Revised. Stages: feasibility (step 0), Pfam baseline, custom HMMs, leave-clade-out. | 3.9 |
| D10 | Revised. Order of six PRs; tests first in PR 3. | 5 |
| D11 | Revised. The 2026-10-09 descriptive test becomes the measurement file, capped at `smoke`. | 3.10 |
| D12 | Draft `narrow_ok` values; the owner signs each in the PR. | 3.11 |
| D13 | Revised. `leakage_exception` key; reword `LEAKAGE_NOTE`. | 3.12 |
| D14 | Revised. `primary_label` precedence with `mechanism_without_class`. | 3.1 |
| D15 | Revised. Tiered gating; guards live in `src/` and blocking suites. | 3.13 |
| D16 | Revised. `negative_basis` per stratum; PM-TM is `assumed` and excluded. | 3.7 |
| D17 | A measure-less entry can only be `unvalidated`. | 3.7 |
| D18 | Revised. `not_applicable` acts as `not_called` in `and`/`or`; `not(not_applicable)` is `not_applicable`; neutral only in the `other_*` list. | 3.5 |
| D19 | New. Per-family sensitivity uses the family's own members; per-positive `leakage`. | 3.10 |
| D20 | New. Each Kiss slice measures its own thing. | 3.8 |
| D21 | New. HPCC re-runs are a numbered acceptance task of PR 3. | 3.14 |
| D22 | New. Code location and test-first order; who refreshes stale call files. | 3.13, 3.14 |

## 3. Decisions in detail

Each block gives the decision, the alternatives, the reason, what changes and where, and the acceptance test. "Fails before" means the test fails on `main` and passes after the change.

### 3.1 D1 and D14: top-level labels

- **Decision (Owner decision 2026-10-17, D1).** The calls stay independent flags. Add one derived column `primary_label` (PROVISIONAL name) with a documented precedence. Remove the word "partition" from `docs/reports/2026-10-07-classifier-runs-and-scale.md:30,70`.
- **Alternatives.** (a) A residual call `other_no_surface_with_mechanism` (SH-D1 option 1). (b) Only reword and print the unlabeled count. (c) A strict partition.
- **Reason.** Independent flags keep every call testable on its own. Precedence in code does not change `config_sha256` (`engine.py:82-95`), so D1 alone makes no call file stale.
- **Decision (Owner decision 2026-10-17, D14).** *Revised 2026-10-17 after the Fable review (FD-A5).* Precedence, first match wins:
  1. `cell_wall_adhesion_candidate` called.
  2. `surface_attachment_candidate` called (with a basis sub-label).
  3. `other_surface_no_mechanism` called.
  4. `other_not_surface` called.
  5. `mechanism_without_class`: any called mechanism not captured by labels 1 and 2, whatever the signal-peptide status. This holds the no-SP case, the SP case with a mechanism that no composite reads (M5 `sowgp_ortholog`) and an SP protein whose only mechanism is the antigen rank (143 proteins in the one proteome where the module applies, FD-A5).
  6. `not_assessable`.
  - *Originally:* `serodiagnostic_marker_candidate` at rank 3 and a bucket `mechanism_no_class` (no label above and at least one mechanism `called`); the later draft renamed it `mechanism_no_signal_peptide`, which held only the no-SP case.
  - Allergen and antigen calls (`iuis_allergen_*`, `serodiagnostic_marker_candidate`, `cocci_specificity_rank_top15` as a label source) stay out of `primary_label`. A `secondary_label` may come later (OPEN, Q6). The antigen rank is still a mechanism for the `other_*` calls (`categories.yaml:100,106`) and so counts for label 5.
  - `surface_attachment_candidate = called` and `other_surface_no_mechanism = called` cannot both hold (`engine.py:242-253`), so labels 2 and 3 do not compete (FD-A5, sound).
- **Changes.**
  - New derived label code (PROPOSAL: `src/cellsurface_sorting_hat/labels.py`; precedence in code, not in `categories.yaml`).
  - `outputs.py:110-127` (`write_wide`): add the column for the default gate. `write_long` (`outputs.py:93-98`) is not changed.
  - `outputs.py:199-208` and the "Calls" section (`outputs.py:276-285`): a table "Primary label" with one count per label.
  - `docs/class_descriptions.yaml`, `analysis/class_table/make_classes_md.py`, `docs/CLASSES.md`: the calls are labels, not a partition.
  - Reports: reword `docs/reports/2026-10-07-classifier-runs-and-scale.md:30,70`.
- **Acceptance tests.**
  - Fails before: a fixture with P1 = no SP and repeat02 called gets `mechanism_without_class`.
  - Fails before: an SP protein with only the antigen rank called gets `mechanism_without_class`.
  - Fails before: a docs test greps the classifier-runs report for "partition" and expects no match.
  - Property test (FD-A5): every protein in the fixtures and in the committed S288C run gets exactly one `primary_label`.

### 3.2 D2: family evidence tier; what `active` means

- **Decision (Owner decision 2026-10-17).** The 8 untested families stay active (decision of 2026-10-09 stands). Add `specificity_tier` to `data/sorting_hat/family_table.tsv`. `active` means "enabled by owner", not "tested". Fix the docs that say otherwise.
- **Alternatives.** (a) Switch the 8 off until tested. (b) Keep them active and say nothing. (c) Keep them active but out of the narrow call (D4, also taken).
- **Reason.** The owner wants the families in the broader call now. The table must say how well each family is known.
- **Changes.**
  - Table: new columns `specificity_tier`, `narrow_ok`, `narrow_ok_by`, `narrow_ok_date` (24 data rows, 11 columns today; `modules/pfam.py:9-21`). Names PROVISIONAL.
  - `modules/pfam.py:9-21`, `:35-42`, `:45-80`: read and validate the new columns.
  - `modules/pfam.py:31-32` (docstring of `NoActiveFamilyError`): reword to "no family of the module is enabled by the owner".
  - `docs/reports/cell-wall-fusion/04-repo-coverage-and-gaps.md:22,227` and `docs/superpowers/specs/2026-10-10-lectin-class-design.md:266,270`: mark as wrong; say what `active` means (WF-F1, SH-D3).
  - Also stale (SH-D2), fixed with D4: `docs/CLASSES.md:13`, `outputs.py:29-34,294,317`, `2026-10-04-orchestrator-design.md:191`, `2026-10-09-v1-class-list-design.md:26`.
  - Notes "specificity test not run" in rows 8, 9, 12, 13 are older than the family test; update them when D11 converts it.
- **Acceptance tests.**
  - Fails before: `load_family_table` on a table without `specificity_tier` must raise `FamilyTableError`; today it accepts it.
  - Fails before: a docs test greps `pfam.py` and the named docs for "passed its specificity test" and expects none.

### 3.3 D3: tier vocabulary and derivation

- **Decision (Owner decision 2026-10-17).** *Revised 2026-10-17 after the Fable review (FD-A3, FD-A4, FD-A10, FD-A12).* The tier reuses the status words (`status.py:14`): `untested`, `smoke`, `estimated`. Code **derives** the tier from a committed measurement file; nobody edits it by hand.
  - The derivation lives in `src/cellsurface_sorting_hat/` (PROPOSAL: a module next to `status.py`), not in `analysis/`. It calls `status_from_measure` (`status.py:105-133`); it does not copy the rule.
  - A measure carries `measure_basis` (PROVISIONAL name): `descriptive_shared_sets` or `per_family_reviewed`. `descriptive_shared_sets` is capped at `smoke`. `estimated` needs `per_family_reviewed`.
  - A missing `n_clusters_pos` or `n_clusters_neg` is treated as NOT meeting the cluster floor. Today the floor is skipped for absent keys (`status.py:128-130`), so a measure without cluster counts can reach `estimated` (FD-A3). All 20 committed call-status entries are `smoke`, so none is lowered; a test proves it.
  - `smoke` means: member and non-member lists exist and are documented, below the floors of `docs/agent-tasks/COMMON-RULES.md` item 7. It does not say "owner-reviewed" (FD-A4).
  - CI fails if `family_table.tsv` disagrees with the derived tier.
  - One sentence goes into this spec, the CLASSES.md plan task and the table header (FD-A12): *the family tier describes the family model's own evidence; the call status describes the call's measurement; neither implies the other.*
  - *Originally:* a script `analysis/family_tests/derive_tier.py`; `smoke` = "reviewed lists exist"; no `measure_basis`; no rule for missing cluster counts.
- **Alternatives.** A new vocabulary (`E0` to `E3`); a hand-edited column with a review rule.
- **Reason.** One vocabulary for module status, call status and family tier. A derived value cannot drift from its evidence. The code location lets the blocking CI job import it without crossing the `analysis/` conftest boundary (FD-A10).
- **Changes.**
  - New committed measurement file (PROPOSAL `data/sorting_hat/family_measure.tsv`): one row per family with counts, clusters, intervals, `measure_basis`, list sources, date, HMMER and Pfam versions, per-positive `leakage` summary (D19).
  - New test `tests/cellsurface_sorting_hat/test_family_tier.py` (blocking suite).
  - `status.py:128-130`: strict cluster rule (Owner decision 2026-10-17, Q7: it applies in `status_from_measure` to ALL module and call measures; the 20 committed call files stay `smoke`, and a test proves they load).
- **Acceptance tests.**
  - Fails before: a copy of the table with one stronger tier fails the test.
  - 25 positives, 25 negatives, 22 and 22 clusters, narrow intervals, `per_family_reviewed`: `estimated`. The same with `descriptive_shared_sets`: `smoke`. 5 and 5: `smoke`. No row: `untested`.
  - A row without `n_clusters_pos` derives at most `smoke` (fails before, FD-A3).
  - Mutations: the derivation ignores cluster counts; ignores `measure_basis`.

### 3.4 D4: the narrow call reads a family only with evidence and sign-off

- **Decision (Owner decision 2026-10-17).** `cell_wall_adhesion_candidate` reads a family only if `specificity_tier` is `smoke` or better AND `narrow_ok=yes`. Untested families still feed `surface_attachment_candidate`. Implementation notes added after the Fable review:
  - `pfam.COLUMNS` (`modules/pfam.py:24`) is one list for all four Pfam modules. It becomes per-module, so only `pfam_adhesion` gets `hit_narrow` (FD-A9).
  - The call hash must cover `family_table.tsv` (tier, `narrow_ok`). `config_sha256` does not, so the narrow call could go stale silently (FD-A9). `_family_digest` (`modules/cli.py:225-232`) adds both fields.
  - The report warns if every `narrow_ok` is `no`: the narrow call is then `not_called` for every protein (FD-B8).
  - The `narrow_ok` sign-off records an evidence summary (members, hits, leakage, near-neighbour negative hits), not a bare `yes` (FD-A15).
  - The published counts in `docs/reports/2026-10-07-classifier-runs-and-scale.md` are **already stale on main** (the wall-family row says "PA14 only"; CFEM and eight more families came later; `outputs.py:294` prints "PA14, CFEM"). They need a re-run and a dated note, not an edit framed as a D4 consequence (FD-A9).
- **Alternatives.** Leave both composites on one flag; switch the families off; split by module name.
- **Reason.** SH-C5: structural wall proteins (Pir-like, CFEM, Hyr1) reach a call named "adhesion" although their own data show non-adhesive hits.
- **Changes.**
  - `modules/pfam.py:24`, `pfam_rows` (`:125-171`): `hit_narrow` only when tier is `smoke` or `estimated` and `narrow_ok` is `yes`; `hit_narrow` implies `hit`.
  - `modules/pfam.py:45-80`: reject `narrow_ok=yes` with tier `untested`, or without `narrow_ok_by` and `narrow_ok_date` (same pattern as `:68-69`).
  - `modules/cli.py:225-232`, `:278-283`: digest and `params`.
  - `categories.yaml`: keep `wall_family_domain` (`:53-54`). Add an evidence call (PROPOSAL `wall_family_narrow_domain`, `pfam_adhesion.hit_narrow`). `cell_wall_adhesion_candidate` (`:59-64`) reads it. `surface_attachment_candidate` (`:69-75`) keeps the broad call. The `mechanism` lists (`:100,106`) keep the broad call.
  - `cli.py:154-181` (`required_columns`): an old `pfam_adhesion` table without `hit_narrow` stops the run (`cli.py:214-216`). Rebuild the Pfam modules.
  - `outputs.py:29-34,294,317`, `docs/class_descriptions.yaml`: say which families feed which call.
- **Acceptance tests.**
  - Fails before: an SP protein with a hit in an `untested`, `narrow_ok=no` family: broad `called`, narrow `not_called`. On `main` both are `called`.
  - A hit in a `smoke`, `narrow_ok=yes` family: both `called`.
  - `load_family_table` refuses `narrow_ok=yes` with tier `untested`.
  - `hit_narrow == 1` implies `hit == 1`.
  - The Pfam modules other than `pfam_adhesion` have no `hit_narrow` column.
  - The family digest changes when tier or `narrow_ok` changes.
  - The report warns when no family has `narrow_ok=yes`.

### 3.5 D5 and D18: unknown mechanism input and `not_applicable`

- **Decision (Owner decision 2026-10-17, D5).** *Revised 2026-10-17 after the Fable review (FD-A2).* An unknown mechanism input makes the `other_*` calls `not_assessable` (today the engine skips it and lists it in `left_out`, `engine.py:243-253`). A new value `not_applicable` marks species-limited calls (D18). The `truth` handling changes in the same PR: `calibration/cli.py` rejects the value at `:132` (`CALL_VALUES`) and `:290-293`, so every `calibrate truth` run on a table with the value would exit 2 before `:480`; `outputs.py:280-285` counts only three values and would drop it silently. The earlier text said "a silent negative"; that was wrong.
- **Decision (Owner decision 2026-10-17, D18).** *Revised 2026-10-17 after the Fable review (FD-A1, FD-A13).* `not_applicable` acts as `not_called` INSIDE `and` and `or`, and `not(not_applicable)` is `not_applicable` (added 2026-10-17 on the owner's confirmation of the Fable suggestion). No existing call changes. It is neutral ONLY in the `other_*` mechanism list (`engine.py:242-253`).
  - *Originally:* "neutral in and/or; dropped from operands". Read literally, `serodiagnostic_marker_candidate = and(rank, signal_peptide_protein)` (`categories.yaml:79-82`) would become `called` for every secreted protein outside C. immitis RS (FD-A1).
  - Truth table (rows: value of the first input; `NA` = `not_applicable`; `C` = called, `N` = not_called, `U` = not_assessable). An NA input is first replaced by `N`, then the Kleene rule applies:

    | input | `not` | `and` with C | `and` with N | `and` with U | `or` with C | `or` with N | `or` with U |
    |---|---|---|---|---|---|---|---|
    | C | N | C | N | U | C | C | C |
    | N | C | N | N | N | C | N | U |
    | U | U | U | N | U | C | U | U |
    | NA | C | N | N | N | C | N | U |

    `and` and `or` never return NA. A leaf `test` or `flag` on an NA row returns NA. `not(NA)` is `not_applicable` (Owner decision 2026-10-17, on the Fable suggestion FD-A1: the neutral value must not turn into `called`). No call in `categories.yaml` uses `not` today (the grammar allows it, `categories.yaml:9`), so this affects only future calls.
  - In `_eval_other`, an NA mechanism is skipped; a called mechanism still gives `not_called`; an unknown mechanism gives `not_assessable`.
  - Status bookkeeping (FD-A13): an NA leaf contributes no status (like `not_assessable` in `engine.py:174-177`). The reason travels in the basis column: the basis lists `not_applicable: <module>`. So a `not_called` record caused by an NA leaf never has an empty basis; the engine currently gives such a record `unvalidated` with an empty basis (`engine.py:348-349`).
- **Alternatives.** Keep today's rule and report the count of `other_*` with a non-empty basis; mark such records `partial`; drop NA from operands (rejected, FD-A1).
- **Reason.** SH-C1: probe P2 gets `other_surface_no_mechanism = called` while `surface_attachment_candidate = not_assessable`. FD-A1 shows the drop rule inverts a call for 283 to 774 proteins per proteome.
- **Changes.**
  - `logic.py:7-10`: add `NOT_APPLICABLE`; `k_not`, `k_and`, `k_or` per the table (`logic.py:19-45`).
  - `engine.py:180-187` (`_row`), `:190-215` (`_eval`): state `not_applicable` gives `NOT_APPLICABLE`.
  - `engine.py:236-253` (`_eval_other`): rule above; `other_basis` lists the unknown mechanisms.
  - `engine.py:30`: bump `ENGINE_SEMANTICS`; it is in every `call_hash` (`engine.py:497-508`). All 20 committed call-status files become stale (`call_status.py:179,194`); D22, D21 say who refreshes them.
  - `outputs.py:276-289`: a fourth count column. `calibration/cli.py:132,290-293,473-488`: `truth` accepts the value, counts it as not assessable and reports the count. Grep `calibration/panel.py` and `analysis/` readers of `calls.long` for the three-value assumption (CANNOT TELL, FD-A2).
  - `categories.yaml:11` comment and `docs/CLASSES.md`.
- **Acceptance tests.**
  - Fails before: probe P2 (SP called, `pfam_hsba` row `error`): `other_surface_no_mechanism` is `not_assessable`.
  - Probe P3 (`antigen_lookup` `not_applicable`, rest `not_called`): the rank is `not_applicable`, `serodiagnostic_marker_candidate` is `not_called` with basis `not_applicable: antigen_lookup`, `other_surface_no_mechanism` is `called`.
  - An SP protein outside RS is never `serodiagnostic_marker_candidate = called` (guards FD-A1).
  - The full table above, all combinations, as a unit test.
  - `calibrate truth` on a table with `not_applicable` succeeds and reports the count.
  - Mutations: NA treated as unknown; NA dropped from `and`; NA contributing a status; semantics not bumped.

### 3.6 D6: R0 exception made visible; `derived` marker; R0 snapshot check

- **Decision (Owner decision 2026-10-17).** R0 (the Phase C SignalP call) keeps its owner exception: no leakage cap, status `estimated` (`docs/paper/03-status-and-validation-rules.md:91`). The basis string says "leakage unmeasured (SignalP training overlap), uncapped by owner decision". Composite calls that inherit a status get a `derived` marker (`2026-10-09-v1-class-list-design.md:114`, O5). Compare the two R0 snapshots (`modules/cli.py:265-269,458-459`).
- **Alternatives.** Cap R0 at `smoke` (SH-C2 option 2); leave the wording.
- **Reason.** SH-C2, SH-C3, SH-C6.
- **Changes.**
  - Exception text via the `leakage_exception` key (D13, 3.12). The status basis is built in `status.py:196-201` and `engine.py:341-347`.
  - Derived marker (PROPOSAL, format OPEN, Q5): a `derived_from` column naming the leaf calls. `evaluate` (`engine.py:334-349`) sets it for composites and `other_*`.
  - R0 snapshot check (PROPOSAL): `_condition_table` (`modules/cli.py:428-459`) records a digest of the `(id, value)` pairs; the engine recomputes it after `load_modules` (`cli.py:366`) and warns in `report.md` and `run.json` on a mismatch. The SignalP identity is only tool, version and mode (`modules/cli.py:240-247`).
  - `outputs.py:47-48` (limit 9): add the same words.
- **Acceptance tests.** Fails before: (a) `pfam_adhesion` built against R0 table A, engine run with table B: the report shows the warning. (b) A composite record carries the derived marker. (c) A record whose status comes from R0 `estimated` has the exception text in its basis.

### 3.7 D7, D16 and D17: status loader, `negative_basis`, measure-less entries

- **Decision (Owner decision 2026-10-17, D7; direction unchanged).** `load_status_source` (`status.py:136-164`) re-validates on load: it recomputes the status from the embedded measure, the leakage cap and the taxon rules, and rejects a status higher than the derived one (lower is allowed). Add `negative_basis` (PROVISIONAL). Files without it default to `assumed`. How R0 passes: see D16 and D13.
- **Decision (Owner decision 2026-10-17, D16).** *Revised 2026-10-17 after the Fable review (FD-A7).* `negative_basis` is recorded PER STRATUM. Status-driving numbers (headline `n_neg`, specificity, interval) use only strata whose basis is `curated`: N-int and N-sec, defined by non-IEA GO location terms (the evidence-code mix is recorded in the notes; a three-value split `curated_experimental` / `curated_inferred` / `assumed` may come later). PM-TM (predicted topology) is `assumed`, is reported in its own row, and does not count toward the floors or the interval.
  - *Originally:* one `negative_basis` per measure, an `assumed` measure capped at `smoke`. Phase C pools all strata into one `n_neg` (`phasec.py:229`; strata at `:39`, `:259-268`), so one value was false either way (FD-A7).
  - **Gate before the engine PR:** compute whether N-int plus N-sec alone meet the floors (20 negatives, 20 clusters, half-width 0.10) for R0. The result goes into the plan. If they do not, R0 falls below `estimated` and D6 and D13 need the owner again (OPEN, Q4). This changes the published R0 specificity number.
- **Decision (Owner decision 2026-10-17, D17; unchanged).** An entry without a `measure` can only be `unvalidated`; a higher status is rejected on load. All 20 committed call-status files carry a measure on every entry, and every entry is `smoke`, so none breaks (FD-A14). Today `measure.py:79-81` allows a measure-less entry below `estimated`; the old "cap at `smoke`" default is removed from the plan and from Q8.
- **Alternatives.** Trust files written by the tool (today); a signed hash in each file; one value per measure.
- **Reason.** SH-C7: hand-edited or old files can claim `estimated` for taxon 1. SH-C4: assumed negatives measure absence of annotation. The call-file loader already does the first part (`call_status.py:66-102`).
- **Changes.**
  - `status.py:17-28` (`MEASURE_KEYS`), `:74-94` (`_check_measure`): accept the per-stratum field.
  - `status.py:105-133`: use only `curated` strata for the status numbers (one place, so `measure.py:63,84-88` and `call_status.py:91` get the rule).
  - `status.py:136-164`: derive the status; apply the leakage cap from `measure.notes` with `call_status.py:37-38,61-63`, unless the entry has the D13 key; refuse taxa below `MIN_TAXON` (`call_status.py:39`, move to `status.py`). Rank checks stay at write time (`calibration/cli.py:135-180`).
  - `calibration/measure.py:27-42` and `calibration/cli.py:92-97,442-529` (`truth`): option for the basis; `calibration/phasec.py:269` sets it per stratum.
  - `docs/paper/03-status-and-validation-rules.md`: add the rules.
- **Acceptance tests.**
  - Fails before: a status file with taxon 1 and `estimated` must raise on load.
  - A file that says `estimated` with a measure that derives `smoke` must raise; a lower status loads.
  - 25 positives, 25 negatives, narrow intervals: all-`curated` gives `estimated`; the same with a PM-TM stratum counted gives the same status as without it (PM-TM does not move the interval); a missing field gives `smoke`.
  - A measure-less entry with `smoke` is refused; with `unvalidated` it loads. All 20 committed files still load.

### 3.8 D8 and D20: truth set for wall remodeling and fusion genes

- **Decision (Owner decision 2026-10-17, D8).** *Revised 2026-10-17 after the Fable review (FD-A8).* The truth set is Kiss et al. 2019 (PMID 31501435) Supplementary Data 3 and 4 (`docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv`, 181 data rows), with a FEASIBILITY STAGE first (D9 step 0):
  - Map the JGI protein IDs to sequences. Report rows per organism and clusters per clade.
  - 37 of 181 rows have no `protein_id_jgi`. `species` is empty for all rows, and `node_of_emergence` (blank for the same 37) is the cluster's ancestral node, not the protein's clade. The clade of a row must come from the organism code in the JGI ID.
  - Reconcile counts with the TSV (checked): adhesion 51 rows in 26 clusters; wall 122 in 44; fusion 8 in 6 (the fusion slice includes 4 pheromone-receptor rows). Report 06 states "23 clusters, 45 curated genes" for the set; the difference must be explained before `counts.md` is written.
  - A clade can be held out only if it has at least 5 independent clusters (PROVISIONAL: the owner was unsure of 5). The feasibility report must show how many clades each floor (for example 3, 5, 8) would admit; the owner may adjust it before any model is evaluated. Below the floor: descriptive recall and leave-cluster-out.
  - Second source for Basidiomycota and non-hyphal clades: UniProt reviewed entries plus literature rows (13 Basidiomycota rows only).
  - Cluster at 30% identity and 50% coverage; split `seed-like` and `seed-independent` (`analysis/lectin_validation/measure.py:328-329,402-409`). Negatives: `data/curated/adhesins/hard_negative_seeds.tsv` and the shared set of 161 rows (`data/controls/hard-negatives-shared/controls.tsv`).
  - *Originally:* the set was taken as given, with a "known weakness" note and the count question as Q10.
- **Decision (Owner decision 2026-10-17, D20).** Each Kiss slice measures its own thing:
  - Adhesion rows (51 / 26 clusters): recall of the active adhesion Pfam families. This is family coverage, not the sensitivity of the call.
  - Wall rows (122 / 44): recall of the `analysis/wall_fusion_report` wall models.
  - Fusion rows (8 / 6): descriptive only; no tier claim.
  - Wall rows are never truth for the adhesion call. They serve only as a DESCRIPTIVE false-hit panel for the narrow call, never as negatives, because the Kiss cell-wall category includes structural proteins (FD-B1).
- **Alternatives.** UniProt reviewed entries alone; Pfam seed members; the S. cerevisiae set of report 07 (23 proteins).
- **Reason.** It is the only curated, clade-spread list in the repo. The Fable count shows 17 organism codes, four of which carry 78% of the rows, and 19 of 26 adhesion clusters hold one organism: per-clade and per-cluster recall are nearly the same number there.
- **Changes.** None in `src/`. New data and scripts under `data/controls/wall-fusion-kiss2019/` and `analysis/` (PROPOSAL), after the feasibility report.
- **Acceptance tests.** The feasibility report lists, per organism and clade, the rows with and without a sequence and the clusters; it states the loss. The L1 validator accepts the new directory. A three-row fixture shows recall per cluster and per clade and puts 35% identity in `seed-like` and 20% in `seed-independent`.

### 3.9 D9: model approach, in stages

- **Decision (Owner decision 2026-10-17).** *Revised 2026-10-17 after the Fable review (FD-A8).* Step 0 is the feasibility stage (3.8).
  1. Baseline: the current Pfam models (21 models of `analysis/wall_fusion_report/wall_fusion_models.tsv` and the active family rows) on the D8 truth set, per the D20 slices.
  2. Custom profile HMMs (HMMER) only where Pfam fails. The baseline decides.
  3. Leave-clade-out (or leave-cluster-out where clades are too thin) is required before any tier above `smoke`.
  4. No ESM models now (issue #15; gate of 300 clusters).
  - *Originally:* no step 0; leave-clade-out only.
- **Alternatives.** Embeddings first; custom HMMs for every cluster; fine-tuning a PLM.
- **Reason.** Pfam is cheap and measured. The data, not the method, limits the design (FD-A8).
- **Precedent.** `2026-10-08-hydrophobin-custom-hmm-design.md`; `analysis/hydrophobin_truth/lco.py` (test `tests/hydrophobin_truth/test_lco.py`).
- **Acceptance tests.** A custom model without a leave-clade-out measurement cannot derive a tier above `smoke`. Stage 1 reproduces `analysis/wall_fusion_report/results/wall_fusion_counts.tsv` for the shared models.

### 3.10 D11 and D19: the family measurement

- **Decision (Owner decision 2026-10-17, D11).** *Revised 2026-10-17 after the Fable review (FD-A3, FD-A4).* The 2026-10-09 descriptive family test (E1: 54 adhesins with sequence; N1: 137 hard negatives in 76 clusters; `docs/reports/2026-10-09-family-tests.md`) is converted into the committed measurement file with `measure_basis=descriptive_shared_sets` (capped at `smoke`). Hyr1 has no E1 hit and stays `untested`. Update the stale "specificity test not run" notes when converting.
  - What the test can and cannot show (FD-A4). Every family is scored against adhesins of all classes, so a family for one clade (Candida_ALS_N, 10 of 54) shows its share of a mixed panel, not its sensitivity. Hyr1 scores 0 of 54 because no E1 protein belongs to its family (`2026-10-09-family-tests.md:39`). The negatives are not relatives of any family, so narrow families look near 1.0 by construction; the two families that hit negatives (CFEM 12, PIR1-like_C 3) are the broad ones. The report itself calls the lists agent-curated and not owner-reviewed (`:61`), the positives unclustered (`:62`), and the result "counts, not estimated sensitivity or specificity" (`:3`).
  - *Originally:* the conversion was proposed without a basis field and `smoke` meant "reviewed lists".
- **Decision (Owner decision 2026-10-17, D19).** A per-family measurement uses the family's own, independently evidenced members as the sensitivity denominator. The all-E1 hit count is reported as coverage only and never called sensitivity. Each positive records `leakage`: 30% identity or more to the family's seed alignment is `in_reference`, else one of the five COMMON-RULES values, with the same cap as other statuses (`call_status.py:37`, `:98-100`). A family measure whose positives are mostly `unknown` or `in_reference` caps at `smoke` (FD-B6, FD-B7).
- **Changes.** Measurement file columns for the basis, the denominator type and the leakage summary; `analysis/family_tests/README.md`: say that E1 gives coverage. `pfam.specificity_report` and `calibrate pfam-specificity` (`pfam.py:174-194`, `calibration/cli.py:530-556`) are the existing path for near-neighbour negatives and may be used by a later per-family measurement (PROPOSAL).
- **Acceptance tests.** A row built from E1 and N1 derives `smoke` at most. A column `sensitivity` for a family uses its own members; the coverage column is named `coverage`. A family whose positives are mostly `in_reference` derives `smoke` at most.

### 3.11 D12: draft `narrow_ok` values

- **Decision (Owner decision 2026-10-17).** Draft `narrow_ok`: `yes` for PA14, Flocculin, Flocculin_t3, Candida_ALS_N, Candida_ALS, Flo11, GLEYA; `no` for CFEM, PIR1-like_C, Hyr1. The owner signs each row in the PR. The PR records before and after counts of the classifier runs.
- **Alternatives.** All `no` until signed row by row (my earlier default).
- **Reason.** The owner reviewed the per-family counts. The criterion behind the order is not yet written down (FD-A15); the sign-off summary of D4 fixes that.
- **Note.** The `yes` rows rest on a descriptive test (D11), so each is `smoke` at most. GLEYA hits one N2 negative; CFEM is `no` on 12 N1 hits (FD-A15). Flocculin (3 of 54) and Flo11 (6 of 54) rest on few observations. Rows not named by D12 (for example Bys1, ALS_M, the allergen and hydrophobin rows) get `narrow_ok=no` (Owner decision 2026-10-17, Q8): the narrow call reads only `pfam_adhesion` rows, so the column has meaning only there. The owner still signs each value.

### 3.12 D13: leakage exception and note wording

- **Decision (Owner decision 2026-10-17).** *Revised 2026-10-17 after the Fable review (FD-A6).* An optional measure key `leakage_exception` (reason, signer, date). The loader accepts an uncapped `estimated` only with it. ALSO reword `LEAKAGE_NOTE` (`calibration/phasec.py:41`) to `leakage: unknown (...)` in the same PR. The regex of `call_status.py:37-38,61-63` takes the first word after `leakage:`; the current note ("leakage: overlap between ...") is read as the value `overlap`, which is not allowed (`call_status.py:37`), so the call-file rules would refuse the R0 entry (`call_status.py:96-98`).
  - *Originally:* only the key, no change of the note.
- **Alternatives.** A flag without text; cap R0 at `smoke`.
- **Reason.** The key lifts the cap on purpose and records who decided; the note then carries a legal token.
- **Acceptance tests.** `leakage_of()` on every R0 entry returns a member of `LEAKAGE_VALUES`. The loader refuses an uncapped `estimated` without the key and accepts it with the key.

### 3.13 D15 and D22: where the guards live; test-first order

- **Decision (Owner decision 2026-10-17, D15).** *Revised 2026-10-17 after the Fable review (FD-A10).* Tiered gating. The L1 control-set validator and the tier derivation live in `src/cellsurface_sorting_hat/`, with their tests in the blocking job. Blocking suites import nothing from `analysis/`. Any new test directory is added to the blocking job in `.github/workflows/build_and_test.yml:46`. Analysis suites stay informational (`:51`, `continue-on-error`) and the owner promotes them by hand.
  - *Originally:* the validator at `analysis/_common/validate_controls.py` with tests in `tests/controls/`, and a derive script under `analysis/family_tests/`; both outside the blocking job.
- **Decision (Owner decision 2026-10-17, D22).** Code location and test-first order: see D15 and D10. All `ENGINE_SEMANTICS` and `categories.yaml` changes make the 20 committed call-status files stale (`call_status.py:179,194`); the plan says who refreshes them and when (D21).
- **Alternatives.** Keep all new checks informational.
- **Reason.** A check that cannot fail a PR does not guard anything.
- **Acceptance tests.** The CI file lists every new test directory in the blocking job. No blocking test imports `analysis/` (a grep test).

### 3.14 D21: HPCC re-runs

- **Decision (Owner decision 2026-10-17).** The HPCC re-runs are a numbered task inside PR 3's acceptance: `calibrate phasec` for R0 (so it carries `leakage_exception` and the reworded note) and the two repeat call-status files. PR 3 is not marked ready until they are done. The refreshed status tables are committed like the existing ones. The same commit adds an interim sentence to CLASSES.md. The job scripts are shown to the owner first, and NOTHING is submitted to SLURM without the owner's explicit go-ahead each time. Heavy work never runs on the login node (FD-A11, FD-B3).
- **Alternatives.** Leave the re-runs outside the PR (the earlier plan, task 3.9).
- **Reason.** R0 and every call file drop in status between the merge and a re-run, and the docs go wrong in the interval (FD-A11). The R0 status files are not in git (`docs/CLASSES.md:11`); only `calibrate phasec` writes them.
- **Acceptance tests.** The committed new status files load. CLASSES.md contains the interim sentence until the files are current.

## 4. Testing framework

Six layers. Accepted layers, and the decision that accepts each: L0 (D15), L1 and L3 (D15, D22), L2 (D3, D11, D19). L4 and L5 stay proposals.

| Layer | What it checks | Where it runs | Label |
|---|---|---|---|
| L0 code and CI | Unit tests, lint. CI jobs for the 5 skipped test dirs; the 2 macOS job-script tests made portable. | GitHub CI | Accepted (D15; PR 2) |
| L1 control-set integrity validator | One validator for every `data/controls/<task>/` set: COMMON-RULES columns, leakage values, sequences present, cluster file covers every row, counts match `counts.md`. Code in `src/cellsurface_sorting_hat/`. | CI, blocking | Accepted (D15, D22) |
| L2 family measurement | Counts, clusters, intervals, `measure_basis`, per-positive leakage on the control sets. Tier derived (D3). | HPCC | Accepted (D3, D11, D19) |
| L3 golden pipeline fixtures and docs-consistency | Small FASTA and module tables through the engine to expected calls, including the gap cases. `CLASSES.md` agrees with `categories.yaml` and the family table. | CI, blocking | Accepted (D15, D22) |
| L4 T0-T6 benchmark harness | Issue #12, per mechanism class and clade, never one pooled number (`docs/model-review/2026-09-27-review-and-framework-plan.md`, section 12 and P2). | HPCC | PROPOSAL |
| L5 pinned measurements | Committed measurements with a tolerance, re-run when Pfam or HMMER changes. A re-run outside tolerance asks for review; it never edits the pin. | HPCC | PROPOSAL |

- **L1.** Each control set has its own checks today. No shared validator exists. The D8 directory is the first user.
- **L3 docs-consistency.** `tests/class_table/test_class_table.py` checks that `CLASSES.md` is generated; it does not read the family table (SH-D2). The new test checks that every family of the broad and narrow calls is named in the generated text and that the tier column agrees with D3.
- **L3 golden.** `tests/cellsurface_sorting_hat/golden/` holds `calls.long.expected.tsv` and `report.expected.md`. The expected rows for the gap cases are written by hand before the engine changes (D10). Golden files are regenerated in their own commit, which the owner reviews separately.

## 5. Order of work (D10)

*Revised 2026-10-17 after the Fable review (FD-A10).* Each item is a separate PR. The plan has the tasks.

| PR | Content | Depends on |
|---|---|---|
| 1 | This spec and plan (with the three reviews). | none |
| 2 | Cheap fixes (PR #88, open): stale docs, the 2 macOS test failures, CI jobs for the 5 skipped test dirs. | PR 1 |
| 3 | Engine changes, STARTING with the new tests (strict xfail, including the gap cases) which are flipped as each change lands: D1/D14, D5/D18, D6/D13, D7/D16/D17, D2/D3/D4/D11/D12. Golden files regenerated in their own commit. The HPCC re-runs (D21) are a numbered acceptance task. | PR 2 |
| 4 | Further tests that did not fit PR 3: L1 control-set validator tests and L3 golden extensions. | PR 3 |
| 5 | Feasibility stage, then baseline recall on HPCC (D8, D9 steps 0 and 1). | PR 4 |
| 6 | Custom HMMs with leave-clade-out (D9 steps 2 and 3). | PR 5 |

*Originally:* PR 3 had the tests after the code and PR 4 held all the tests; PR 5 started with recall. PROPOSAL: PR 3 may be split into 3a (D1, D5, D6) and 3b (D2, D3, D4, D7); the split is not a decision.

## 6. Risk and compatibility

| Change | Who breaks | Default and migration |
|---|---|---|
| `family_table.tsv` gets 4 columns (D2, D4) | Any reader of the table; `load_family_table` requires all `FAMILY_COLUMNS` (`pfam.py:49-51`); the Pfam command (`modules/cli.py:253`). | The loader requires the new columns; the PR rewrites the table in the same commit with the D12 draft values, signed by the owner. |
| `pfam_adhesion` gets `hit_narrow` (D4) | Old module tables stop with `missing column` (`cli.py:214-216`). | Rebuild the four Pfam modules (job scripts in `scripts/sorting_hat/`; D21 rules apply). Module identity changes (`_family_digest`). |
| `categories.yaml` edit (D4) and `ENGINE_SEMANTICS` bump (D5) | `config_sha256` and every `call_hash` change; all 20 call-status files are stale (`call_status.py:179,194`; `engine.py:497-508`). A stale file never raises a status. | D21: refresh in PR 3 acceptance; interim sentence in CLASSES.md. |
| New call value `not_applicable` (D5, D18) | Readers of `calls.long`/`calls.wide`; `calibration/cli.py:132,290-293`; `outputs.py:280-285`. | Document in `CLASSES.md`; handle in `truth`; grep other readers. |
| Strict cluster rule (D3) | Any measure without cluster counts that today reaches `estimated`. | The 20 committed call files are `smoke`; a test proves they load. Scope of the rule: all measures (Q7, settled). |
| Status-file schema (D7, D13, D16) | Module status files without the new fields. | Default `assumed`. R0 files are rewritten by `calibrate phasec` (D21). Without the rewrite R0 in A. nidulans drops below `estimated` and `CLASSES.md:11` is wrong. |
| Loader refusal (D7, D17) | A hand-edited or old file with a status above its measure, or a measure-less entry above `unvalidated`. | Loader raises with the file name. Fix the file; do not relax the loader. |
| Report and wide-table columns (D1, D5, D6) | Scripts that read `calls.wide` by position; the golden files. | New columns are appended. `calls.long` gains the derived marker column (Q5). Regenerate golden files in their own commit. |

## 7. Proposals not decided

These came from the reviews. They are cheap and fit PR 2 or later.

1. WF-F2, WF-F1, WF-F4, WF-F5: mark the contradicted claims in reports 01 to 04, fix the commit attribution for Flocculin, close the Gel4 doubt as "annotation inconsistent with its own cluster", log replacements in the errata.
2. WF-F6: bring the duplicate-ID warning of `lectin_report.py:323-339,741-748,771-778` into `wall_fusion_report.py`.
3. WF-F7: compare the model NAME of Pfam to the TSV `name` at run time; derive the "verified" list from report 05.
4. WF-F10: close the PRM1 score question in the README (`hmmscan` domain score 475.5 versus `hmmsearch` sequence score 475.9).
5. SH-L1: apply `_plain` to `args.notes` (`calibration/cli.py:252-254,499-502`). SH-L2: put the pessimistic sensitivity bound in the measure.
6. SH-C8: commit the other five hydrophobin status files or soften "smoke in 8 proteomes" in `CLASSES.md:14`.
7. WF-F8, WF-F9: the 21-model criterion and the `model_label` text; decide in PR 5.

## 8. Open questions

Questions Q1 to Q11 of the first draft are settled or moved. Settled: old Q1 and Q2 by D1/D14; old Q3 by D18; old Q4 by D11 (descriptive test, capped at `smoke`); old Q5 by D13 and D16 (the key, per-stratum basis); old Q6 (no decision text; the narrow call stays out of `basis_calls`, PROPOSAL unchanged); old Q7 by D12; old Q8 by D17; old Q10 by D8 (feasibility stage); old Q11 by D19. Old Q9 is now item 7 of section 7.

STILL OPEN or PROVISIONAL:

| # | Question | My default (PROPOSAL) |
|---|---|---|
| Q1 | The floor of 5 independent clusters for holding out a clade (D8). The owner was unsure. | The feasibility report shows how many clades each floor admits; the owner adjusts it before any model is evaluated. |
| Q2 | Exact column and key names: `specificity_tier`, `narrow_ok`, `leakage_exception`, `negative_basis`, `measure_basis`, `primary_label`, `primary_basis`, `derived_from`. | Use these names until the owner changes them. |
| Q3 | Whether the Kiss JGI IDs map to sequences, and how many rows are lost (37 of 181 have no ID). | The feasibility stage answers it. |
| Q4 | Per-stratum floor check for R0 (D16): do N-int and N-sec alone meet 20 negatives, 20 clusters and half-width 0.10? | Compute it before the engine PR; if not, return to the owner. |
| Q5 | Format of the `derived` marker. | A `derived_from` column naming the leaf calls. |
| Q6 | Allergen and antigen `secondary_label`. | Later; out of `primary_label` now. |
| Q7 | SETTLED (owner, 2026-10-17): the strict cluster rule (D3) applies in `status_from_measure` to all measures. | The 20 committed call files stay `smoke`; a test proves they load. |
| Q8 | SETTLED (owner, 2026-10-17): `narrow_ok=no` for rows D12 does not name (Bys1, ALS_M, allergen, hydrophobin, HsbA rows); the column has meaning only on `pfam_adhesion` rows. | The owner signs each value in the PR. |

## 9. Conflicts between the decisions and the code, and how each was resolved

| # | Conflict | Resolved by |
|---|---|---|
| 1 | The R0 `estimated` entry has no `leakage:` value; the note is read as `overlap` (`calibration/phasec.py:41`; `call_status.py:37-38,61-63,96-98`). D7's cap would reject or demote it. | D13 (key and reworded note) and D16 (per-stratum basis); Q4 is the remaining gate. |
| 2 | `calibration/cli.py:132,290-293` rejects a fourth value, and `outputs.py:280-285` drops it. The first draft located the fix at `:480` and called the effect a silent negative. | D5 (changed in the same PR). |
| 3 | The "dropped from operands" rule would make `serodiagnostic_marker_candidate` `called` for every secreted protein outside RS (FD-A1). | D18 (acts as `not_called` in `and`/`or`; `not(NA)=NA`). |
| 4 | D4 changes the counts in `2026-10-07-classifier-runs-and-scale.md`; those counts are already stale on main. | D4 note (re-run and dated note), D12 (before/after counts). |
| 5 | Rows 8, 9, 12, 13 say "specificity test not run" but the 2026-10-09 test has counts. | D11 (update notes at conversion). |
| 6 | The status loader would skip the cluster floor for absent keys (`status.py:128-130`), so the E1/N1 test could derive `estimated`. | D3 (missing counts do not meet the floor; `descriptive_shared_sets` capped). |
| 7 | One `negative_basis` cannot describe pooled Phase C negatives (`phasec.py:229,259-268`). | D16 (per stratum). |
| 8 | New tests were planned outside the blocking CI job (`.github/workflows/build_and_test.yml:46,51`). | D15, D22. |
| 9 | No task re-ran `phasec` or the call files. | D21. |
| 10 | Tier and call status use the same words for different objects (FD-A12). | D3 sentence. |
| 11 | The Kiss TSV counts (51/26, 122/44, 8/6) differ from "23 clusters, 45 genes" of report 06. | D8 (reconcile in the feasibility stage). |
| 12 | `data/curated/adhesins/hard_negative_seeds.tsv` is a smaller seed file (32 lines) next to the shared set of 161 rows. | D8 names both; the owner confirms both are meant. |
| 13 | D18 as first worded gave `not(NA)=called`; the Fable review suggested `not(NA)=NA` (FD-A1). No call in `categories.yaml` uses `not`. | RESOLVED (owner, 2026-10-17): `not(NA)` is `not_applicable`; `and` and `or` treat NA as `not_called`. The truth table in 3.5 and its unit test follow this. |

## 10. Not in this spec

- No change to `src/`, `categories.yaml`, `family_table.tsv`, status files or CI by this PR.
- No ESM or fine-tuned model (D9 step 4, issue #15).
- No new call beyond the evidence call of D4.
- The family-table Pfam lists of the wall/fusion report (21 models) are not changed.

## 11. Sources

Reviews: `docs/reviews/2026-10-17-sorting-hat-logic-review.md`, `docs/reviews/2026-10-17-wall-fusion-products-review.md`, `docs/reviews/2026-10-17-decisions-d1-d18-review-fable.md`. Specs followed for form: `2026-10-07-per-call-status-design.md`, `2026-10-09-sowgp-like-calls-design.md`, `2026-10-10-lectin-class-design.md`. Code read: `engine.py`, `logic.py`, `status.py`, `call_status.py`, `outputs.py`, `cli.py`, `categories.yaml`, `modules/{pfam,cli,lookups}.py`, `calibration/{cli,measure,phasec}.py`. Data read: `data/sorting_hat/family_table.tsv`, `docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` (counted). Docs read: `docs/CLASSES.md`, `docs/paper/03-status-and-validation-rules.md`, `docs/agent-tasks/COMMON-RULES.md`, `docs/reports/2026-10-09-family-tests.md`. CI: `.github/workflows/build_and_test.yml`.
