> What this is: an independent review of the logic, workflow and data flow of `cellsurface_sorting_hat` (main at c7eadf0). Verbatim copy; the text below the rule is unchanged.
> Date: 2026-10-17. Source: scratchpad copy of sorting_hat_logic_review.md. Findings are cited in the spec as SH-D1..SH-D3, SH-C1..SH-C9, SH-L1..SH-L2.
> Independent review by an Opus agent; not verified by the owner.

---

# cellsurface_sorting_hat: logic review of workflow and data flow

2026-10-17. Independent read-only review of `main` at `c7eadf0`. Reviewer: Claude Code (claude-opus-5-5).

## What I did

**Read (verified):** `AGENTS.md`; `docs/CLASSES.md`; specs `2026-10-04-orchestrator-design.md` (§3.3–3.4), `2026-10-07-per-call-status-design.md`, `2026-10-09-v1-class-list-design.md` (M2), `2026-10-09-sowgp-like-calls-design.md` (M5, header and integration lines), `2026-10-10-lectin-class-design.md` (owner block and §3.4); `docs/paper/03-status-and-validation-rules.md` §5–6. Source files in full: `engine.py`, `logic.py`, `status.py`, `call_status.py`, `cli.py`, `outputs.py`, `categories.yaml`, `modules/{base,pfam,signalp,repeats,cli}.py`, `calibration/{cli,measure,call_files,phasec}.py`. Partly read: `modules/lookups.py` (`antigen_rows`). Also `data/sorting_hat/family_table.tsv`, committed status files under `docs/reports/data/sorting_hat/call_status/`, and the reports `2026-10-07-classifier-runs-and-scale.md` §3–5, `2026-10-07-sorting-hat-run-and-calibration.md` (grep), `2026-10-08-hydrophobin-validation.md` (grep), `2026-10-09-family-tests.md` (grep), `cell-wall-fusion/04` and `07` in full, and the header of `analysis/lectin_class_report/lectin_report.py`.

**Ran:**
- `uv run --no-sync --with pytest pytest tests/cellsurface_sorting_hat -q`: 818 passed, 9 skipped, 2 failed. Both failures are macOS-only (`realpath -m` is GNU) in `test_data_and_scripts.py` (job-script tests). They are not logic failures.
- `shasum -a 256 categories.yaml` = `e4f3d63d…`. This equals `config_sha256` in the committed `…config-with-family-evidence.json` call files, so the files cited in CLASSES.md are current for the config.
- An in-memory probe of `engine.evaluate` with the packaged config and 3 toy proteins, using `ModuleTable` objects only (no files written; the throwaway script was deleted afterwards). The proteins:
  - P1: R0 not_called, repeat02 called;
  - P2: R0 called, every mechanism not_called, `pfam_hsba` row `error`;
  - P3: R0 called, every mechanism not_called, `antigen_lookup` `not_applicable`.

  Results are quoted in D1 and C1.

**Not read:** the reports `01`–`03`, `05`, `06` of cell-wall-fusion; `calibration/{intervals,panel,stability}.py`; `taxonomy.py`; `fasta.py`; module wrappers `allergen`, `cys8`, `hydrophobin_relaxed`; the test files. I make no claim about them.

Verdict scale: SOUND / CONCERN / DEFECT / CANNOT TELL.

---

## Findings, ranked by severity

### D1. DEFECT (logic and wording): the residual `other_*` calls do not close the class list, and two documents call it a partition

*Verified by reading and by the probe.*
- `other_not_surface` needs SP = not_called **and** every assessable mechanism = not_called (`categories.yaml:95-100`; `engine.py:240-253`). Every composite needs SP = called (`categories.yaml:59-75, 79-82`).
- So a protein with **no signal peptide and a mechanism hit** (repeat, wall family, hydrophobin, HsbA, top-15 rank) gets **no** composite and **no** `other_*`. Probe protein P1 (SP not_called, repeat02 called): `cell_wall_adhesion_candidate`, `surface_attachment_candidate`, `other_not_surface` and `other_surface_no_mechanism` are all `not_called`.
- This group is not small. Examples:
  - intracellular repeat proteins (known limit 4, `outputs.py:29-31`);
  - PA14 hits without an SP are removed before the engine sees them (`pfam.py:144-150`), so they count as "no mechanism". CFEM and the 8 families with no condition are not removed;
  - FIG2-like proteins with no R0 call (07 table row FIG2, `07-model-hit-check.md:29`).
- **Wording.** `2026-10-07-classifier-runs-and-scale.md:30` says SP, candidate, no-mechanism and not-surface "partition the assessable proteins". `:70` says the calls are "a partition of the proteome by signal peptide". Both were already untrue when written, because of the no-SP plus repeat case. Today they are also untrue for SP proteins: SP plus hydrophobin only, or SP plus HsbA only, gives `surface_attachment_candidate` called, `cell_wall_adhesion_candidate` not_called and `other_surface_no_mechanism` not_called.
- `CLASSES.md:22` ("No signal peptide, with the mechanism calls not called") is accurate. It does not tell the reader that the complement has no label.
- **M5 widens the gap.** `sowgp_ortholog` would join the mechanism lists of both `other_*` calls but feed no composite (`2026-10-09-sowgp-like-calls-design.md:43`). An SP plus ortholog-only fragment would then have no class either.
- **Fix options:** add an `other_no_surface_with_mechanism` residual, or state in CLASSES.md and the report that the calls are labels, not a partition, and give the unlabeled count.

### D2. DEFECT (doc/code drift on what feeds the narrow call): 8 newly active wall families are absent from every description of `wall_family_domain`

*Verified.*
- `family_table.tsv` has `active=yes` for these rows, activated 2026-10-09 (commits `cbfdc45`, `d94a8ab`):
  - Candida_ALS_N, line 8;
  - Candida_ALS, line 9;
  - Flo11, line 12;
  - GLEYA, line 13;
  - Flocculin, line 21;
  - Flocculin_t3, line 22;
  - Hyr1, line 23;
  - PIR1-like_C, line 24.
- All of them are in module `pfam_adhesion`, so all of them set `wall_family_domain` and hence `cell_wall_adhesion_candidate`.
- Every description still names only CFEM and PA14:
  - `docs/CLASSES.md:13` ("CFEM: no mature TM helix; PA14: R0 signal peptide");
  - the report label `"wall family domain (PA14, CFEM)"` (`outputs.py:294, 317`);
  - orchestrator spec `:191` ("Bys1, ALS, Flo11, GLEYA, Hyr, PIR … are inactive");
  - M2 decision C5, "10 of 24 rows" active (`v1-class-list-design.md:26`). The table now has 18;
  - known limit 4 cites "(CFEM, Bys1, Als)" (`outputs.py:31`). Bys1 is inactive (`family_table.tsv:3`).
- CLASSES.md is "generated … a test checks that this file is current" (`CLASSES.md:3`). The test checks generation from `class_descriptions.yaml`, not agreement with the family table. So the generated table can be current and still wrong about what the call reads.

### D3. DEFECT (premise of the c-i argument and of a docstring): "`active` means passed a specificity test" is false for 8 of 18 active rows

*Verified.*
- `pfam.py:31-32`: `NoActiveFamilyError` = "No family of the module has passed its specificity test and sign-off".
- Report 04 relies on this meaning: `04-repo-coverage-and-gaps.md:20` ("`active` means the family passed a specificity test and an owner sign-off") and `:225` ("`active` already promises a passed specificity test"). So does the lectin spec: `2026-10-10-lectin-class-design.md:266` item (4) and the `:270` (c-ii) cell.
- The rows themselves say otherwise:
  - "specificity test not run": `family_table.tsv:8, 9, 12, 13`;
  - "not reviewed": `:21-24`.
  - PIR1-like_C (`:24`) is active even though the family test found 0 of 54 E1 positives and hits on 3 N1 negatives plus 1 N2 negative (`2026-10-09-family-tests.md:32, 38`). That report recommends "keep in the broader call only, or review" (`:51`).
- **Effect on c-i.** The recommendation for c-i still stands. It rests on two other grounds that hold: no config-hash churn (`call_status.py:179-180`; `categories.yaml` bytes = `config_sha256`, `engine.py:84-92`), and no new call (lectin owner block `:13-14`). But the "keeps `active` meaningful" argument is already lost. The owner should decide what `active` means and record it in the table header or in `pfam.py`.

### C1. CONCERN (call logic): unknown and not-applicable mechanisms count as absent in `other_*`, so the residual can be `called` while the composite is unknown

*Verified (probe and code).*
- By design (orchestrator spec `:208-212`), `other_*` is computed over assessable mechanisms only (`engine.py:243-253`).
- Probe P2 (SP called; `pfam_hsba` row in state `error`): `surface_attachment_candidate` = `not_assessable`, but `other_surface_no_mechanism` = `called` with `other_basis = hsba_domain,cocci_specificity_rank_top15`.
- So a report count of "other, no mechanism" includes proteins whose mechanism status is unknown.
- Probe P3: `antigen_lookup` writes `not_applicable` for every non-Coccidioides protein (`lookups.py:156-157`). The engine treats every non-`ok` state as unknown (`engine.py:185`). So `cocci_specificity_rank_top15` is listed as "left out because not assessable" for **every** protein outside C. immitis RS (`outputs.py:286-289`). "Not applicable to this taxon" and "unknown" share one label.
- **Suggestion:** split `not_applicable` from unknown in `other_basis`, and report how many `other_*` = called records have a non-empty basis.

### C2. CONCERN (status policy): the leakage rule has a documented exception that the per-call rule would not allow

*Verified.*
- Call files cap any leakage other than `none` (including `unknown`) at `smoke` (`call_status.py:96-100`). The `truth` command requires `--leakage` and caps the status (`calibration/cli.py:92-97, 519-521`).
- The Phase C path writes R0 module entries with **no cap** (`phasec.py:269`, no `cap=`). Its notes say "overlap … with the SignalP 6 training data was not measured" (`phasec.py:41, 251`). That is `unknown` in all but name.
- Result: R0 is `estimated` in A. nidulans (`2026-10-07-sorting-hat-run-and-calibration.md:63`; `CLASSES.md:11`).
- This is an explicit owner decision (`docs/paper/03-status-and-validation-rules.md:91`), and the report states the limit (`outputs.py:47-48`). It is not a bug. It is the one place where an `estimated` claim rests on truth (GO-curated secreted proteins) that plausibly overlaps the model's training data.
- **Suggestion:** state the exception beside the status in CLASSES.md, or cap it to match the rule elsewhere.

### C3. CONCERN (status semantics): composite negatives inherit `estimated`, and the planned "derived" marker was never added

*Verified (code); the example is inferred from the code.*
- For a false AND, only the false inputs decide (`engine.py:226-227`). So in A. nidulans every protein without an SP gets `cell_wall_adhesion_candidate[R0] = not_called` with status `estimated`, because R0's status decides it. This holds even when the protein has a wall-family hit.
- No measurement of the composite exists. Owner Q7 says composite statuses are derived (`per-call-status-design.md:151`).
- The marker in `status_basis` was deferred (M2 O5 default "not added", `v1-class-list-design.md:114`). The basis is `step1_rule@R0:taxon:162425`, the same string as for the leaf.
- Only a general sentence in the report says composite statuses are derived (`outputs.py:206`). `calls.long`/`calls.wide` readers see `estimated` on a composite.

### C4. CONCERN (status rule): the measure has no notion of negative quality

*Verified (code and committed files).*
- `status_from_measure` counts every label-0 row as a negative (`status.py:105-133`).
- Committed truth sets use **assumed** negatives:
  - repeat call S288C, 297 negatives "assumed from absent annotation" (`…Scer_S288C/tandem_repeat_protein…family-evidence.json`, notes);
  - hydrophobin Af293, 9,640 negatives "assumed" (`…pfam_hydrophobin.Afum_Af293_UniProt.json`, notes).
- Today small positive counts and the leakage cap keep both at `smoke`, so no wrong claim results.
- But a future truth set with 20 or more positive clusters and assumed negatives, given `--leakage none`, would reach `estimated` with a near-1 specificity interval that measures annotation absence. Leakage is self-declared by the operator. Nothing checks truth IDs against tuning sets or Pfam member lists.
- **Suggestion:** a measure key such as `negatives: assumed|curated` that caps at `smoke` when the negatives are assumed.

### C5. CONCERN (evidence vs label): the narrow `cell_wall_adhesion_candidate` is fed by families whose own data show non-adhesive hits

*Verified (family table and family-tests report).*
- CFEM: 12 N1 negatives in 9 of 76 clusters and 1 of 54 E1 positives (`2026-10-09-family-tests.md:44`). Kept by owner decision O1 (a).
- PIR1-like_C: structural wall family (its own note, `family_table.tsv:24`), 0 E1 positives, N1 hits (above).
- Hyr1: no E1 positive possible (`family-tests.md:39`).
- The engine cannot restrict a family to the broader call. Both composites read the single flag `pfam_adhesion.hit` (`categories.yaml:53-54, 63, 74`).
- The "candidate" suffix and known limit 4 (`outputs.py:29-34`) disclose that domains are not shown to mediate adhesion. But the call name says "adhesion". This matters for the wall/fusion question (Q5). Structural wall proteins (Pir-like, CFEM, Hyr1) **do** reach the narrow call. Report 07 tested enzymes, not these proteins.
- **Suggestion:** the M2 O1 (b)/(c) split, an adhesion-domain flag apart from a wall-family flag, now has a second reason.

### C6. CONCERN (data flow): the R0 table that sets the PA14 condition inside `pfam_adhesion` is never checked against the R0 table the engine gates on

*Verified (code). Whether the two diverge in practice: CANNOT TELL.*
- `pfam` reads the R0 table given as `--sp-module` at module build time (`modules/cli.py:265-269`). It records only that module's `params_hash` and `artefact_hash` (`modules/cli.py:458-459`).
- The SignalP identity is tool+version+mode only, not the input (`modules/cli.py:240-247, 422-426`).
- If R0 is re-run on an edited FASTA, or a different table is passed, the engine's `step1_rule@R0` and the R0 snapshot inside `pfam_adhesion` can disagree, and no check fires. `_check_ids` only requires a subset match (`modules/cli.py:187-200`).
- The same pattern applies to `hydrophobin_relaxed` and `cys8` (`modules/cli.py:352, 369`); both are experimental and read by no call.
- **Minor, safe direction:** `conditions` enters `params` for all four Pfam modules (`modules/cli.py:278-283`). So a SignalP mode change also makes `pfam_hydrophobin`'s status stale, although hydrophobin has no condition.

### C7. CONCERN (enforcement location): the module status loader trusts files that the writer validated

*Verified.*
- The writer path checks the following: status no stronger than the measure (`measure.py:84-88`); species rank and taxon ≥ 2 (`calibration/cli.py:135-180`); and the leakage cap (`cli.py:520`).
- `load_status_source` checks only shape, duplicates and identity (`status.py:136-164`). It does **not** re-check status against measure, the leakage cap or taxon ≥ 2. `_taxa` in the writer even allows 0 (`measure.py:50`).
- The call-file loader does all three (`call_status.py:79-100`).
- So a hand-edited or old module file could claim `estimated` for taxon 1 and apply it to every protein (`status.py:167-177`).
- The rules "one species per entry" and "leakage cap on truth" are enforced only at write time for module files. They close the loophole only for files made by the tool.

### C8. CONCERN (small-n wording in report 07 and CLASSES)

*Verified.*
- Report 07 is careful ("Measurement only", limits at `:112-118`). It derives "Hat labels" from an older committed run with only PA14 active (`:14`, `:115`), so its labels for SAG1/FLO11/FIG2 are not engine output under the current config.
- Its line `:108` is correct as stated: the ESM classifier calls 6 of 13 wall enzymes above 0.5, and no active family-table row hits any of them. But the 13 enzymes are one species, one protein per gene, with no negatives beyond the set (`:114`).
- CLASSES `:14` "smoke in 8 proteomes" is supported by the report (`2026-10-08-hydrophobin-validation.md:24-33`, 3 to 7 positives each). Only 3 of the 8 status files are committed (`docs/reports/data/sorting_hat/call_status/hydrophobin/`), and CLASSES cites one of them plus the report (`CLASSES.md:30-31`). That is acceptable but thin for "status sources are committed files" (M2 §3, `v1-class-list-design.md:35`).
- The repeat interval for S288C, sensitivity 0.435 [0.067, 0.667], is printed as `smoke`. That is correct. Its lower bound shows the number cannot support any claim, and CLASSES reports it as smoke (`CLASSES.md:12`). SOUND as worded.

### C9. CONCERN (specificity tests vs biology)

*Verified (report 07).*
- The committed PF11766 specificity file for S288C counts SAG1 (α-agglutinin) as a false positive because the member list is Candida-derived (`07-model-hit-check.md:79`).
- The specificity tests that gate `active` therefore depend on member lists that are taxon-biased. A non-member hit is not necessarily a non-adhesin.
- This and D3 together mean `active` is a weak signal in both directions.

### L1. LOW: user `--notes` is not escaped before `call=`/`leakage:` parsing

*Verified (code).*
- `_plain` exists to prevent this (`calibration/cli.py:252-254`), but `args.notes` is placed raw into the notes (`:499-502`).
- It is safe today only because the tool appends its own `call=` and `leakage:` later, and both parsers take the last match (`cli.py:308-309`; `call_status.py:61-63`).
- A future reorder of the notes string would silently break this. Apply `_plain(args.notes)`.

### L2. LOW: complete-case sensitivity

*Verified.*
- `not_assessable` truth rows are dropped from y (`calibration/cli.py:480-482`). The pessimistic bound is in the notes only (`:496-502`), not in `measure.sensitivity`. The status is set from the optimistic number.
- Today's committed files have 0 not-assessable rows, so this has no effect yet.

---

## Answers to the five questions

**Q1. Data flow (module tables → engine → labels → report): CONCERN.**
- Schemas are consistent and checked at every hop:
  - module files have `id`/`state` + required columns (`cli.py:131-151, 154-180, 214-216`), and duplicate IDs are refused (`cli.py:146-147`; `base.py:62-65`);
  - FASTA membership is checked in every wrapper (`modules/cli.py:187-200`; `pfam.py:116-122`; `repeats.py:63-67`), plus a no-match guard in the engine (`cli.py:217-223`);
  - bad values become `bad_value` (`cli.py:234-235, 248-253`);
  - missing rows become `error` and the module becomes `partial` (`cli.py:230-242`).
- Species keys come from the user's `--taxon`/`--taxon-map`. They are recorded as given and not checked against sequences (`outputs.py:46`).
- Silent-meaning risks: C6 (two R0 snapshots); C1 (`not_applicable` read as unknown); D2 (the meaning of `wall_family_domain` changed through the family table without any doc change). The config hash does not cover the family table (by design, `lectin spec:266` item 5), so a call-file status measured before 2026-10-09 stays "current" even though the meaning of `wall_family_domain` widened. This affects only call files that read `pfam_adhesion`, and there are none today (`call_eligible`, `engine.py:452-475`). `pfam_adhesion`'s own identity does change, through `_family_digest` (`modules/cli.py:225-232, 290`).

**Q2. Call logic: DEFECT (D1), plus D2 and C1/C5.**
- The Kleene logic and the "deciding inputs" rule are correct and match the spec (`logic.py:28-45`; `engine.py:223-233`; orchestrator `:172-180`).
- `cell_wall_adhesion_candidate` ⊂ `surface_attachment_candidate` by construction (`categories.yaml:59-75`). This nesting is intended, not an overlap bug.
- `serodiagnostic_marker_candidate` and `other_surface_no_mechanism` are mutually exclusive because top-15 is in the mechanism list (`categories.yaml:106`). Allergen calls are outside the mechanism list, so an SP allergen homolog is also `other_surface_no_mechanism`. That is consistent with the spec (`orchestrator :206`).
- Non-exhaustive: D1. Definition drift: D2.

**Q3. Status/calibration logic: CONCERN.**
- The per-call restriction works: `measured_call_of` plus the leaf-call contributors (`engine.py:286-294`; `cli.py:304-309`). `ref` keeps the callee as the leaf, so R0 measured on `signal_peptide_protein` counts inside composites. That is correct.
- Call files close their loopholes: identity, `call_hash`, config, states, leakage, status ≤ measure, taxon ≥ 2 (`call_status.py:66-102, 177-196`).
- One species per entry is enforced at write (`calibration/cli.py:448-454`; `phasec`'s `_set_taxa` `:189-214`). The status for a species applies to all its strains (`status.py:167-177`). That is intended.
- No circularity is hidden: repeat `tuned_on_truth` → smoke; hydrophobin `partial` → smoke.
- Open gaps: C2 (R0 exception), C4 (assumed negatives), C7 (loader trusts files).
- Small-n results are labelled `smoke` everywhere I checked. No `estimated` claim rests on fewer than 20 clusters (`status.py:126-130`).

**Q4. Evidence and wording: CONCERN, with one DEFECT (D1 "partition").**
- "~12% adhesin precision" is the shipped **ESM classifier**, not the hat (`04-repo-coverage-and-gaps.md:88`). The hat never reads ESM scores. `step1_ml@card` is a declared variant (`categories.yaml:15`) with no wrapper in `modules/cli.py`.
- "27 of 28 hard negatives" is R0 alone (`classifier-runs-and-scale.md:51, 56`). `CLASSES.md:11` says R0 claims "nothing more" than a signal peptide. That is consistent.
- Overstatements: D1 (partition), D2/D3 (active and the family list), C5 (the "adhesion" name on CFEM/PIR1-like hits), C3 (`estimated` shown on composites).

**Q5. Wall/fusion work vs the tool: SOUND for enzymes, CONCERN for structural wall proteins, c-i fits.**
- GPI-anchored wall enzymes (GAS1, CRH1, UTR2, CTS1, SUN4, KRE9, ECM33, EXG1) reach `other_surface_no_mechanism`: R0 called, no family hit, repeat not called (`07:39-50, :82`). That matches known limit 3 (`outputs.py:27-28`) and decision C3 (`v1-class-list-design.md:24`). This is logical for v1, and the hat does not repeat the ESM classifier's 6-of-13 error.
- But "Als models are active" sits beside PIR1-like_C, Hyr1 and CFEM in the same flag (C5, D2). So structural wall proteins do reach the narrow call. Report 07 did not test them (`07:33` lists no hits for those models in its 23-protein set).
- c-i fits the engine. The engine loads every table in `modules/` and lists it in run states and calibration rows (`cli.py:189-192, 412`), so a c-ii report-only module would show up there. Any `evidence:` edit changes `config_sha256` and makes call files stale (`call_status.py:179-180`). c-i avoids both.
- Two caveats:
  - (i) D3: the "active stays meaningful" reason no longer holds.
  - (ii) The lectin script takes its own `--signalp` file, matched by ID (`lectin_report.py:17-30`). The 04 proposal adds `tm.n_tm_mature` and R0 from tool outputs (`04:198`). A c-i script that reads hat module tables should check the run identity (`run.json` `module_identities`) or record which R0 it used. Otherwise it repeats C6 outside the tool.

---

## Things that are sound

- **Three-valued logic.** It is small, tested and exact (`logic.py`). Unknown inputs never decide a result; a `not_assessable` record gets `unvalidated` with an empty basis (`engine.py:231-232, 348-349`).
- **Config validation.** It refuses forward refs, ungated-uses-gated, unknown thresholds, `{step1}` misuse, basis calls that are not ungated, and `other` calls without an earlier per-variant surface (`engine.py:98-163`).
- **Per-call status design.** The eligibility, `call_hash`, staleness reasons and grouping by leaf call follow the spec closely (`engine.py:356-373, 452-511`; `call_status.py:177-196`; spec §3.1–3.5). A stale file never raises a status.
- **Status ≤ measure, at write.** This holds for both file kinds (`measure.py:69-93`; `call_status.py:91-95`). `estimated` needs ≥20/20/20 clusters and half-width ≤ 0.10 (`status.py:97-133`).
- **Phase C import.** It is defensive: counts are re-joined to cluster files and refused on mismatch, there is one species per set from names.dmp, SignalP module and mode must match the R0 record, and intervals are widened at 0/1 and by cluster Wilson (`phasec.py:67-95, 98-141, 232-236, 276-299`).
- **Module wrappers fail closed.**
  - Pfam refuses an unfinished domtbl, one run without `--cut_ga`, or targets not in the FASTA (`pfam.py:90-95, 116-122`).
  - A missing condition value gives `error`, not 0 (`pfam.py:144-166`).
  - Repeats refuse short or non-finite rows (`repeats.py:17-27`).
  - An all-error table is refused (`modules/cli.py:203-213`).
  - An empty BLAST file is refused (`modules/cli.py:324-325`).
- **Module identity.** It covers what decides a call: the family digest over accession, module, condition and active, not free text (`modules/cli.py:225-232`); the repeat script sha (`:316`).
- **Report honesty.** The report prints known limits, status definitions, warnings for partial or absent modules, and stale call files (`outputs.py:17-49, 155-180, 199-275`). "candidate" and "evidence" are used as M2 rule 2 says. Hydrophobin and HsbA are explicitly "not adhesion" (`outputs.py:310-311`).
- **Tests.** 818 pass. The 2 failures are platform-only (GNU `realpath -m` on macOS).
