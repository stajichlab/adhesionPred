> What this is: an independent adversarial review of the owner decisions D1 to D18, by a Fable agent. Verbatim copy; the text below the rule is unchanged.
> Date: 2026-10-17. Not verified by the owner. Finding ids are A1 to A15 and B1 to B8; the spec and plan cite them as FD-A1 etc.
> Source: scratchpad copy of fable_decisions_review.md.

---

# Independent review of owner decisions D1-D18 (cellsurface_sorting_hat), 2026-10-17

Read-only review. Reviewer: Claude Code (Fable 5.1). Repo `main` at `83d5f60`; spec and plan read on
branch `sorting-hat-review-spec` in `/Users/jstajich/projects/adhesionPred-spec`.

**Verified (read or ran).** Spec and plan of 2026-10-17 and both input reviews; `engine.py`,
`logic.py`, `status.py`, `call_status.py`, `categories.yaml`, `modules/pfam.py`, `cli.py`
(`required_columns`, `load_modules`), `modules/cli.py` (pfam command, `_family_digest`),
`calibration/{cli,measure,phasec}.py`, `outputs.py:270-300`; `data/sorting_hat/family_table.tsv`;
`docs/reports/2026-10-09-family-tests.md` + `docs/reports/data/sorting_hat/family_tests/family_tests.tsv`
+ `analysis/family_tests/README.md`; all 20 files under `docs/reports/data/sorting_hat/call_status/`
(parsed with a throwaway python3 -I one-liner, no files written); `data/controls/hard-negatives-shared/controls.tsv`
(161 rows); `docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` (counted with python3 -I);
`.github/workflows/build_and_test.yml:46-68`; `docs/CLASSES.md`; `docs/paper/03-status-and-validation-rules.md`;
`docs/agent-tasks/COMMON-RULES.md`; PR #88 file list via `gh`.
**Not run:** the test suite (no claim needed it). **Inferred** items are marked.

Severity order. Verdicts: SOUND / CONCERN / DEFECT / CANNOT TELL.

---

## A. Findings ranked by severity

### A1. DEFECT (D18 vs D5, vs the code): "not_applicable dropped from operands" turns `serodiagnostic_marker_candidate` into a mass false positive

*Verified by reading.*

- `cocci_specificity_rank_top15` is a `test` node on `antigen_lookup.percentile`
  (`categories.yaml:76-78`). `serodiagnostic_marker_candidate` is
  `and: [{ref: cocci_specificity_rank_top15}, {ref: signal_peptide_protein}]` (`categories.yaml:79-82`).
- `antigen_lookup` writes state `not_applicable` for every protein outside *C. immitis* RS
  (`modules/lookups.py:157`). Today `_row` returns `None` for any non-`ok` state (`engine.py:185-187`),
  the test gives `NOT_ASSESSABLE` (`engine.py:212-213`), and `k_and(NA, called) = not_assessable`
  (`logic.py:28-35`). The call is therefore never `called` outside RS.
- D18 says `not_applicable` is "dropped from operands". Applied literally to this `and`, the operand
  list becomes `[signal_peptide_protein]`, so **every secreted protein in every non-Coccidioides
  proteome becomes `serodiagnostic_marker_candidate = called`** — 283 to 774 proteins per proteome
  (`docs/reports/2026-10-07-classifier-runs-and-scale.md`, "surface, no mechanism evidence" row is the
  SP count minus candidates; the report's own serodiagnostic row is 0 everywhere except RS: 143).
- D18 also contradicts D5 as written: spec 3.5 (`...design.md:135`) proposes `not_applicable` acts as
  `not_called` in `and`/`or`/`not`, and Q3 (`:265`) records that default with the owner's reason ("the
  rank cannot be true outside C. immitis RS"). D18 replaces "acts as false" with "neutral" without
  saying so. Only one of the two can be implemented.
- D18 says nothing about `not`. `k_not` exists and `_eval` supports it (`engine.py:221-222`,
  `logic.py:19-25`); the value set must cover it or the engine raises `ValueError: not a call value`
  from `_check` (`logic.py:13-16`).

**Fix.** Keep D5's rule inside `and`/`or`/`not` (NA behaves as `not_called`; `k_not(NA) = NA`), and
apply D18's "neutral/dropped" rule **only** to the `other_*` mechanism list (`engine.py:242-253`),
which is where it is needed and where the "all operands not_applicable -> not_applicable" clause makes
sense. If the owner wants D18's clause in `and`/`or` too, it must read: *all* operands
`not_applicable` -> `not_applicable`; otherwise NA counts as `not_called`. Write the full 4x4 truth
table into the spec before Task 3.2 starts.

### A2. DEFECT (D18 + D5 vs the code): the `truth` command fails earlier and in a different place than D18 says

*Verified.*

- D18 and spec section 9 item 2 both say the fix is `calibration/cli.py` "around line 480", where
  `not_assessable` is handled, and that an unhandled `not_applicable` would be "a silent negative".
- It would not. `_read_call_values` validates every value against
  `CALL_VALUES = ("called", "not_called", "not_assessable")` (`calibration/cli.py:132, 290-293`) and
  raises `value must be one of ...`. A run output containing `not_applicable` makes **every**
  `calibrate truth` invocation exit 2 before line 480 is reached.
- Same blind spot in the report: `outputs.py:280-285` counts only the three values, so a
  `not_applicable` record disappears from the "Calls" table without a warning.
- `calibration/panel.py` and `analysis/` readers of `calls.long` were not checked for the same
  assumption (CANNOT TELL; they should be grepped in Task 3.2).

**Fix.** Task 3.2 must touch `calibration/cli.py:132` (value vocabulary) and `outputs.py:280-285`
(fourth count column), not only `:473-488`. Add a test that `truth` on a table containing
`not_applicable` succeeds and reports the count.

### A3. DEFECT (D3 + D11 + D12): the derived tier can reach `estimated` on the E1/N1 test, because the cluster floor is skipped when cluster counts are absent

*Verified.*

- D3 delegates the tier to `status_from_measure` (`status.py:105-133`). That function applies the
  20-cluster floor **only for keys that are present**: `for key in ("n_clusters_pos","n_clusters_neg"):
  if key in measure and measure[key] < min_clusters` (`status.py:128-130`). A measure that omits
  `n_clusters_pos` keeps `enough = True`.
- The E1 positives are **not clustered** (`analysis/family_tests/README.md`, Limits; report
  `docs/reports/2026-10-09-family-tests.md:62`). So the natural measurement row has `n_pos = 54`,
  `n_neg = 137`, `n_clusters_neg = 76`, and **no** `n_clusters_pos`.
- With 54 positives and 137 negatives, a family such as GLEYA (12/54, 0 N1 hits) has specificity 1.0
  on a 137-row set (interval very narrow after Wilson on 76 clusters) and a sensitivity half-width near
  0.10. `status_from_measure` can therefore return `estimated` — a tier the owner intends to mean
  "measured" — from a descriptive test that the report itself says "gives counts, not estimated
  sensitivity or specificity" (`2026-10-09-family-tests.md:3`).
- The plan already senses the hole: Task 3.5 test 4 (`plan:162`) says a row with missing cluster
  counts "derives at most `smoke`". That rule is **not** in D3 and **not** in `status_from_measure`.
  Spec 3.3 (`:103`) explicitly says the script must call `status_from_measure` and not copy the rule,
  so the plan's test would fail against the decision as written.

**Fix.** Make it part of D3: a family measure without both cluster counts is capped at `smoke`.
Implement as a wrapper in the derive script *and* (better) as an extra argument or a sibling function
in `status.py`, so the same hole does not stay open for module and call measures.

### A4. DEFECT (D3 + D11 vs the data): the E1/N1 test is not a per-family specificity measurement, and "reviewed lists" is contradicted by the report it comes from

*Verified.*

- **Denominator.** `pos_E1_hit` is "of 54" for every family (`family_tests.tsv` header and column).
  The 54 E1 adhesins are adhesins of all classes. A Candida-specific family (Candida_ALS_N, 10/54)
  cannot hit the non-Candida positives; its 0.185 is not its sensitivity but the share of a mixed
  panel it covers. Hyr1 scores 0/54 because *no E1 protein is a Hyr1/Iff protein*
  (`2026-10-09-family-tests.md:39`), not because the model fails.
- **Negatives.** The 161-row shared set is cross-species (S. cerevisiae 27, A. fumigatus 24,
  V. dahliae 12, ...), chosen as adhesin-architecture hard negatives, not as relatives of each family.
  A narrow family hits 0 of them by construction, so every narrow family gets specificity ~1.0 with a
  tight interval. The two families that DO hit negatives (CFEM 12 in 9 clusters; PIR1-like_C 3 in 1
  cluster) are the broad/structural ones — exactly inverted from what a per-family specificity test
  should reward.
- **"Reviewed".** D3 defines `smoke` as "reviewed member and non-member lists exist". The source
  report says the opposite: "The hard-negative set and the positives were curated by an agent and have
  not been owner-reviewed" (`2026-10-09-family-tests.md:61`). D11 converts that same test into the
  committed measurement file. As written, D11 asserts a review that D3 requires and the repo denies.
- A defensible per-family denominator: a curated member list **for that family** (the positives the
  family claims to detect, from literature, not from Pfam's own seed/full alignment), clustered at 30/50,
  plus negatives drawn from the family's **near neighbours** (for PA14: GH3 beta-glucosidases, already
  known to be the false-hit class, `family_table.tsv:16`; for PIR/PIR1-like_C: structural wall
  proteins; for CFEM: Pth11-like GPCRs, `family_table.tsv:2`). The existing
  `pfam.specificity_report` + `calibrate pfam-specificity` path (`pfam.py:174-194`,
  `calibration/cli.py:530-556`) is built for exactly this and is unused by D11.

**Fix.** Either (a) rename the tier derived from E1/N1 to something that does not claim specificity
(for example `panel_screen`), and forbid it from feeding D4's gate; or (b) keep D11 as the *initial*
tier only, cap every E1/N1-derived row at `smoke` by rule (A3's fix does this), and record in each
measurement row that the denominator is a mixed panel. Also fix D3's `smoke` wording: "lists exist and
are documented", not "reviewed", unless the owner reviews them in the PR.

### A5. DEFECT (D14): the precedence list is not exhaustive; it drops the serodiagnostic case and mislabels the SP gap case

*Verified.*

- D14 drops `serodiagnostic_marker_candidate` from the precedence (spec 3.1 had it at rank 3,
  `...design.md:65`) and says antigen calls stay out of `primary_label`. Consider an RS protein with an
  SP and a top-15 rank and no other mechanism: it is not a `*_candidate` under ranks 1-2; it is
  **not** `other_surface_no_mechanism`, because the rank is in the mechanism list
  (`categories.yaml:106`); it is not `other_not_surface`; it has an SP, so
  `mechanism_no_signal_peptide` is factually false; it is assessable, so `not_assessable` is false.
  **It gets no label.** In the one proteome where the module applies this is 143 proteins
  (`2026-10-07-classifier-runs-and-scale.md`, serodiagnostic row).
- The same hole swallows the M5 `sowgp_ortholog` case that D1 was created to close: SP + a mechanism
  that no composite reads. The spec's name for this bucket was `mechanism_no_class`
  (`...design.md:68`), which covered both gap cases. D14 renamed it `mechanism_no_signal_peptide`,
  which by its own words can only hold the no-SP case.
- **D1+D14 cross-check (asked):** `surface_attachment_candidate = called` and
  `other_surface_no_mechanism = called` cannot both happen — the first needs a mechanism `called`, the
  second needs every assessable mechanism not `called` (`engine.py:242-253`). SOUND; precedence is not
  load-bearing there.

**Fix.** Restore a general residual at the end of the mechanism branch (`mechanism_no_composite`, or
keep the spec's `mechanism_no_class`) covering "at least one mechanism called, no composite called",
and decide explicitly whether the antigen rank counts as a mechanism for the label (it does for the
`other_*` calls, so consistency says yes). Add the property test the spec already proposes
("exactly one label per protein") and a fixture for the RS rank-only protein.

### A6. DEFECT (D6 + D7 + D13): the leakage regex reads the R0 note as `leakage: overlap`

*Verified.*

- D7 says the loader applies the leakage cap "with the regex of `call_status.py:37-38,61-63`". That
  regex is `(?:^|[\s;])leakage:\s*([a-z_]+)` and takes the last match.
- The Phase C note is `LEAKAGE_NOTE = "leakage: overlap between the Phase C positives and the SignalP 6
  training data was not measured"` (`calibration/phasec.py:41`), appended last to the measure notes
  (`phasec.py:247-252`). The regex captures **`overlap`**, which is not in `LEAKAGE_VALUES`
  (`call_status.py:37`). The call-file loader would raise "measure.notes needs 'leakage: <none|...>'"
  (`call_status.py:96-98`); a module-file loader built on the same code raises or caps at `smoke`.
- D13's `leakage_exception` key fixes the *policy* but not the *parse*: the string still mints a bogus
  leakage value, and any future consumer of `leakage_of()` on R0 gets `"overlap"`.

**Fix.** In the same PR, reword `phasec.LEAKAGE_NOTE` to carry a legal token plus prose, e.g.
`"leakage: unknown (overlap with SignalP 6 training data not measured)"`, and let D13's
`leakage_exception` key be what lifts the cap. Add a test that `leakage_of()` on every R0 entry returns
a member of `LEAKAGE_VALUES`.

### A7. DEFECT (D16 + D7 + D6): one `negative_basis` value cannot describe the pooled Phase C negatives, and the pooled value decides whether R0 stays `estimated`

*Verified.*

- Phase C specificity is `1 - FPR` **over all negative classes of the set** (`phasec.py` module
  docstring, lines 5-6), and `n_neg` comes from `n.all.neg` (`phasec.py:229`). The strata N-int, N-sec
  and **PM-TM** (`phasec.py:39`) are extra per-stratum rates (`phasec.py:259-268`), not separate
  entries.
- D16 says GO non-IEA negatives (N-int, N-sec) count as `curated` and the PM-TM stratum (predicted TM
  topology) is `assumed`. Both are inside the same pooled `n_neg`. A single `negative_basis` on that
  measure is therefore false either way: `curated` over-claims for the PM-TM part; `assumed` caps R0 at
  `smoke` (D7) and breaks D6, `docs/paper/03-status-and-validation-rules.md:91,137` and `CLASSES.md:11`
  ("estimated in A. nidulans").
- PM-TM is also the one stratum whose "negative" status is a prediction, i.e. exactly the case
  COMMON-RULES item 18-19 excludes ("a prediction from a tool is not independent of that tool").

**Fix.** Make `negative_basis` a rule over strata, not a free choice: `curated` only when every
negative stratum in the measure is curated; otherwise `assumed`, unless the owner signs a documented
exception (the D13 key generalises to this). Alternatively exclude PM-TM from the pooled specificity
and report it as a stratum only — but that changes the published R0 specificity number, so it is an
owner decision, not an implementation detail. Either way the owner must state it before Task 3.1.

### A8. CONCERN (D8 + D9): the Kiss truth set cannot support "recall per clade" or leave-clade-out as described

*Verified by counting `docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` (181 data rows).*

| slice | rows | clusters |
|---|---|---|
| adhesion | 51 | 26 |
| cell fusion (both category spellings) | 8 | 6 |
| cell wall (both spellings) | 122 | 44 |

- **No clade column exists.** `node_of_emergence` is the inferred ancestral emergence node of the
  *cluster* ("lca of Dikarya", "lca of Eukaryotes", ...), not the taxon of the protein, and it is blank
  for 37 of 181 rows. The `species` field is empty for all 181 rows. The clade of a row can only be
  derived from the JGI organism code inside `protein_id_jgi`.
- **37 of 181 rows have no `protein_id_jgi` at all** — no sequence is obtainable for them, before any
  mapping loss. Q10 asks about mapping loss; the floor is already 20%.
- **17 distinct organism codes, dominated by four:** Canalb1 36, Aspfu1 26, Neucr2 26, Sacce1 24
  (112 of the 144 ID-bearing rows = 78%). All of Basidiomycota is 13 rows (Ustma22 5, Schco3 3,
  CryneH991 3, Copci 2). Leave-clade-out over Saccharomycotina / Eurotiomycetes / Sordariomycetes is
  possible; holding out Basidiomycota leaves a test fold of ~13 proteins; every other clade is 1-4 rows.
- **46 of 76 clusters contain one organism**; for adhesion, 19 of 26 clusters are single-organism. So
  "recall per cluster" and "per clade" are largely the same number for the adhesion slice, and a
  held-out clade removes whole clusters rather than testing transfer within them.
- **Circularity (asked).** Kiss's gene sets were assembled from orthology plus functional annotation;
  the repo's own reconciliation notes the 23-family set is "42% yeast plus Candida genes"
  (`docs/reviews/2026-10-17-wall-fusion-products-review.md:76`, citing 06 row 17). Measuring *Pfam*
  recall on it is not circular in the training sense (Pfam models were not built from this list), but
  the **seed-like / seed-independent split D8 borrows from `analysis/lectin_validation` is the only
  thing that separates "the model knows this protein family" from "the model generalises"** — and with
  19 of 26 single-organism adhesion clusters, the seed-independent fold will be small. Report it with
  its n; do not pool.
- **Count mismatch to check (CANNOT TELL).** D8 states "23 clusters, 45 curated genes in 06"; this TSV
  gives 51 adhesion rows in 26 clusters. Report 06 may count only the Data 3 'gene list' sheet. The
  numbers in the spec and in the TSV must be reconciled before the set is built, or `counts.md` will be
  wrong on arrival.
- **"Adhesion" validity (asked).** The 51 adhesion rows are Kiss's *functional category* label, mostly
  hydrophobins and yeast flocculins/agglutinins. The 122 "cell wall" rows are not adhesion and must
  never be scored against `cell_wall_adhesion_candidate`; D8 does not say which slice is the truth for
  which call. That is a missing decision (B1).

**Fix.** Before PR 5: (1) add a derived `organism`/`clade` column from `protein_id_jgi` and state the
clade mapping source; (2) decide what the 37 ID-less rows do (drop, with the loss in `counts.md`);
(3) state the minimum fold size for a reported per-clade recall (and that Basidiomycota will be below
it); (4) reconcile 45/23 against 51/26.

### A9. CONCERN (D4 vs the code): the narrow flag needs a per-module column set, and it widens the "config does not cover the family table" hole

*Verified.*

- `pfam.COLUMNS = ["hit", "families"]` (`pfam.py:24`) is used for **all four** Pfam modules in one
  loop (`modules/cli.py:301, 309`). Writing `hit_narrow` for `pfam_adhesion` only means `COLUMNS` must
  become a per-module function; otherwise the other three modules gain a column that is always empty,
  and `_bad_value` / `required_columns` behaviour must be re-checked for them.
- The engine side is fine: `required_columns` (`cli.py:154-180`) walks the config and keys on the
  literal module name, so `{flag: pfam_adhesion.hit_narrow}` makes the column required for
  `pfam_adhesion` alone; an old table stops the run with `missing column 'hit_narrow'`
  (`cli.py:213-216`). Spec 3.4 is right about this.
- `_family_digest` hashes `(pfam_acc, module, second_condition, active)` (`modules/cli.py:225-232`);
  D4/spec 3.4 correctly adds tier and `narrow_ok`. Note what is still uncovered: **`config_sha256` does
  not cover `family_table.tsv`**, so the meaning of the narrow call now depends on a file outside the
  call hash. Today no call file reads `pfam_adhesion`, so nothing goes stale wrongly — but that is a
  coincidence of the current call set, not a protection.
- **Does D12 make a committed number wrong?** Spec section 9 item 3 says D4 "changes the counts in
  `2026-10-07-classifier-runs-and-scale.md:31` (6 to 15 candidates per proteome)". That is not the live
  risk: that report's wall-family row is labelled "wall family domain (PA14 only)" and predates CFEM
  (active 2026-10-08) and the eight families of 2026-10-09. **The published counts are already stale on
  `main`**, by more than D12 will change them. The honest fix is a re-run plus a dated note, not an
  edit framed as a D4 consequence. (`outputs.py:294` still prints the label "wall family domain (PA14,
  CFEM)", also stale — PR #88 touches `outputs.py`; confirm it covers this line.)

### A10. CONCERN (D15 + D10): the blocking-CI claim is not true of the paths the plan actually creates

*Verified.*

- The blocking job runs `pytest tests/surface_glyco tests/cellsurface_sorting_hat tests/class_table`
  (`.github/workflows/build_and_test.yml:46`). The `analysis-tests` matrix is
  `continue-on-error: true` (`:51`) with seven suites.
- D15 says the new guards "live in the already-blocking core suites". Two of them do not, as planned:
  the **L1 control-set validator test** is planned at `tests/controls/test_validate_controls.py`
  (plan Task 4.1) — a new directory in neither job — and the **tier check** imports a script under
  `analysis/family_tests/` from a test in `tests/cellsurface_sorting_hat/` (plan Task 3.5), crossing
  the conftest/sys.path boundary that the workflow comment at `:55-56` exists to avoid.
- Sequencing: D15 says these checks block merges, but D10 puts them in **PR 4**, after the engine
  changes land in PR 3. For the duration of PR 3 the new rules are unguarded, and the golden files are
  regenerated (plan Tasks 3.3 step 7, 3.4 step 5) by the same agent that changes the engine — the
  classic "regenerate the oracle" risk. The plan's own mitigation ("read the diff line by line") is a
  human step, not a check.

**Fix.** Either add `tests/controls` to the blocking job in PR 2 (empty or with the fixture tests
only) and keep the derive script importable as a package module, or state plainly in D15 that L1
blocks only from PR 4 onward. For the golden files: write the expected rows for the six gap proteins
**by hand in PR 3** (plan Task 4.2 already says "written by hand ... before the engine is read again"
— move that task into PR 3).

### A11. CONCERN (D13 + D10): no task re-runs Phase C, so R0 drops to `smoke` between PR 3 and an HPCC run

*Verified / partly inferred.*

- D13 says R0 files are "rewritten once to carry" `leakage_exception`. The R0 status files are not in
  git — `CLASSES.md:11` says "The status files are in git-ignored work directories; the report is
  cited." Only `calibrate phasec` can write them, from the Phase C `metrics.json` and cluster files on
  HPCC.
- Plan Task 3.4 changes `phasec` to emit the key; plan Task 3.9 re-measures only `tandem_repeat_protein`
  and is explicitly "not done by this plan's PR". **Nothing re-runs `phasec`.** Between the PR 3 merge
  and that run, D7's loader sees an R0 entry with no exception key (and, per A6, a bogus leakage token)
  and demotes R0 to `smoke` or refuses the file. `CLASSES.md:11`, the paper rules file (`:91,:137`) and
  every report that says "estimated in A. nidulans" become wrong in the interval.
- The same applies to `ENGINE_SEMANTICS` (D5, `engine.py:30`) and the `categories.yaml` edit (D4): all
  20 committed call-status files lose their `call_hash`/`config_sha256` match
  (`call_status.py:179,194`) and become stale, so `CLASSES.md:12` ("smoke where the call file is
  current") goes false until the HPCC re-measure. The plan acknowledges staleness but assigns no owner
  or date to the re-run.

**Fix.** Make the HPCC re-runs (phasec + the two repeat call files) a numbered task with a stop point
inside PR 3's acceptance, or add a documented interim sentence to `CLASSES.md` that the PR writes in
the same commit.

### A12. CONCERN (D2 + D3): the tier column and the call status use the same three words for different objects

*Verified.* `specificity_tier` is per **family row** of `family_table.tsv`; `status` is per **call and
species** (`status.py:14`). After D11, every `pfam_hydrophobin` family row would be `untested` (all
hydrophobin families score 0/54 E1, `family_tests.tsv`), while `CLASSES.md:14` reports the
`hydrophobin_domain` call as `smoke in 8 proteomes` from eight real measures (verified: the eight
committed files, n_pos 3-7, leakage `partial`). A reader seeing `untested` families behind a `smoke`
call will not know which to believe.

**Fix.** Either use a distinct vocabulary for the family column (D3 considered and rejected `E0-E3`),
or state in one sentence in `CLASSES.md` and in the table header that the tier describes the *family
model's own evidence* and the status describes the *call's measurement*, and that neither implies the
other.

### A13. CONCERN (D7 + D17 + D18 vs `engine.py:225-253`): status bookkeeping for the fourth value is undefined

*Verified (code), consequence inferred.*

- `_leaf` attaches a contributor for every value except `NOT_ASSESSABLE` (`engine.py:174-177`), and a
  record with no contributors gets `UNVALIDATED` with an empty basis (`engine.py:348-349`).
- No decision says whether a `not_applicable` leaf contributes. If it does, the record's status is
  limited by the status of a module that did not apply (for `antigen_lookup`, `unvalidated` — "the
  ranking prints NOT CALIBRATED", `CLASSES.md`). If it does not, then under A1's drop-rule an `and`
  that returns `not_called` purely because of an NA operand has **no** deciding module, i.e. a
  `not_called` record with status `unvalidated` and an empty basis — indistinguishable in
  `calls.long` from a bug.
- D18's "a call whose operands are all not_applicable is itself not_applicable" also needs the
  contributor rule: such a record should carry no status, like `not_assessable` records
  (`engine.py:231-232`).

**Fix.** State it in D18: an NA leaf never contributes to the status basis, and a record whose value is
`not_applicable` takes `unvalidated` with an empty basis, exactly as `not_assessable` does today. Add
the case to the mutation list of Task 3.2.

### A14. CONCERN (plan vs D17): the plan still carries the old Q8 default

*Verified.* D17 decides that an entry without a `measure` can only be `unvalidated`. The plan says
"Entries without `measure`: per Q8 (default: capped at `smoke`)" (`plan:69`) and the spec's Q8 row
(`...design.md:270`) says the same. Minor, but it is the kind of drift that gets implemented from the
plan rather than from the decision.

I verified D17's factual premise: **all 20 files under `docs/reports/data/sorting_hat/call_status/`
carry a `measure` on every entry**, and every entry's status is `smoke`. So D17 breaks no committed
file. `measure.py:79-81` currently allows a measure-less entry at `smoke` — D17 tightens it; that is a
one-line change in `_check_entries` plus the loader.

### A15. CONCERN (D4 + D12): `narrow_ok` is signed per row, but the evidence behind the signature is one number with no denominator

*Verified.* D12's yes/no split tracks `pos_E1_hit` almost exactly: yes for PA14 9, Candida_ALS_N 10,
Candida_ALS 9, Flocculin_t3 8, Flo11 6, GLEYA 12, Flocculin 3; no for CFEM 1, PIR1-like_C 0, Hyr1 0.
Given A4, that number is a mixed-panel coverage count, not a specificity result. Two consequences:
- GLEYA is signed `yes` although it hits an N2 negative (`family_tests.tsv`, `neg_N2_hit = 1`), while
  CFEM is signed `no` on 12 N1 hits — a defensible ordering, but the rule behind it is not written
  down anywhere, so the next family cannot be signed consistently.
- Flocculin (3/54, 26 scan proteins in 2 proteomes) and Flo11 (6/54, 9 proteins in 5 proteomes) are
  signed on very few observations. `active_date` for these rows is 2026-10-09 and the note still says
  "not reviewed" (`family_table.tsv:21-22`).

**Fix.** Write the sign-off criterion into the table header or `pfam.py` (what evidence a `yes`
requires), so `narrow_ok` is auditable rather than a remembered judgement.

---

## B. Decisions the owner has NOT taken that the implementation will force

1. **Which Kiss slice is truth for which call.** D8 names Data 3 and 4 as "the truth set" but the file
   holds three disjoint categories (adhesion 51, cell wall 122, fusion 8). Recall of
   `cell_wall_adhesion_candidate` against the 122 wall rows would be meaningless; recall of a future
   wall-remodeling report against the 51 adhesion rows likewise. Needed before PR 5 Task 5.1.
2. **The `negative_basis` rule for mixed strata** (A7) — a value, or a rule, must be chosen.
3. **Who re-runs `phasec` and the two call-status measures, and when** (A11), plus what the docs say in
   the interval.
4. **The full 4-value truth table, including `not`** (A1, A13), and whether an NA leaf contributes to
   the status.
5. **The `primary_label` residual and whether the antigen rank is a mechanism for labelling** (A5).
6. **Whether a family measurement needs its own leakage field.** Pfam models are built from seed
   alignments that may contain the E1 positives; spec Q11 raises it and defaults to "yes", but D3/D11
   do not include a leakage column in the measurement file. Without it the tier rule has no cap that
   matches the one all other statuses carry (`03-status-and-validation-rules.md:72-75`).
7. **A denominator convention for per-family sensitivity** (A4) — otherwise every future family test
   repeats the "of 54" problem.
8. **What happens to the `active` flag's remaining users.** `NoActiveFamilyError` makes an entire Pfam
   module `unavailable` when no family is active (`pfam.py:136-137`, `modules/cli.py:296-309`). Nothing
   says what happens if every `narrow_ok` is `no`: `wall_family_narrow_domain` would be a call that is
   `not_called` for every protein in every proteome, with no warning in the report.

## C. The three weakest decisions

1. **D18** — it is stated in six words ("neutral in and/or ... dropped from operands") that, read
   literally, invert a call's meaning for thousands of proteins (A1), contradict D5's own recorded
   rule and Q3 answer, leave `not` undefined, and mis-locate the code change (A2).
2. **D11** — it converts a test whose own report says it is descriptive, unclustered and un-reviewed
   into the committed *measurement* file that the tier (and, through D4, the narrow call) is derived
   from (A3, A4). Everything downstream inherits that denominator.
3. **D16** — one `negative_basis` value is applied to a pooled negative set whose strata have
   different provenance, and the value chosen decides whether the repo's single `estimated` claim
   survives (A7).

## D. What is sound

- **D1's core choice** (independent flags + a derived column) keeps `categories.yaml` and every call
  hash unchanged for the labelling change. Verified: precedence in code does not touch
  `config_sha256` (`engine.py:82-95`), so no call file goes stale from D1 alone.
- **D2** is honest and matches the code: `active` really is only an owner switch
  (`pfam.py:66-69` checks `active_by`/`active_date` and nothing else).
- **D4's mechanism** is implementable exactly as specified on the engine side: `required_columns`
  keys on the literal module name, so only `pfam_adhesion` tables must be rebuilt (A9), and
  `_family_digest` is the right place for the new fields.
- **D7's direction** (re-validate on load, allow lower, reject higher) matches what the call-file
  loader already does (`call_status.py:91-100`) and closes a real asymmetry (`status.py:136-164`
  checks only shape). Verified it breaks none of the 20 committed files: every entry is `smoke` with a
  measure.
- **D17's premise is true** — all 20 committed call-status entries carry a measure (verified).
- **D9's staging** (Pfam baseline before custom HMMs, leave-clade-out before any tier above `smoke`,
  no ESM) is the right order and has repo precedent (`analysis/hydrophobin_truth/lco.py`). The
  constraint is data, not method (A8).
- **D10's PR order** is sound apart from the tests-after-engine inversion called out in A10.
- **D6's snapshot check** addresses a real gap: `pfam` records only the R0 module's `params_hash` and
  `artefact_hash` (`modules/cli.py:456-459`) and the SignalP identity is tool+version+mode only
  (`modules/cli.py:240-247`), so two different R0 tables can produce identical identities.
</content>
</invoke>
