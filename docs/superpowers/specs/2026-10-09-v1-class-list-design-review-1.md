# Review 1 of the M2 spec (v1 class list), revision 1

*2026-10-09. Independent review of `docs/superpowers/specs/2026-10-09-v1-class-list-design.md`. Reviewer: Claude Code (claude-opus-5-5). Branch `fungi-scan-and-specs` at `cdfea39`. No existing file was edited. Nothing was committed. Every finding below was checked with a command or by reading the named file. Statements I did not check are marked "not checked".*

## 0. What I ran

- `sha256sum src/cellsurface_sorting_hat/categories.yaml` gives `8c32064a…`.
- I read all status files under `_workdir/sorting_hat/*/status/` and `docs/reports/data/sorting_hat/call_status/` (Python, `json`).
- I read the latest stored runs: `_workdir/sorting_hat/Scer_S288C/out_attach/` (config `8c32064a…`) and `_workdir/sorting_hat/Cimm_RS/out_hyd/` (config `fd397b94…`).
- I ran a synthetic engine check (scratchpad `rule3.py`, `/usr/bin/python3.12`). It calls `engine.call_eligible`, `engine.evaluate` with status hooks, and `engine.collect_evidence` with `pfam_adhesion.families` added.
- `PYTHONPATH=src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q`: 818 passed, 7 skipped. The run left `git status` unchanged.
- `git rev-list --left-right --count origin/main...HEAD` gives `0 51`. `git branch -vv` shows no upstream for `fungi-scan-and-specs` or `hydrophobin-validation`.

## 1. Findings

### Section 3: the call table and its status column

**1. MAJOR. The table describes an unmerged local branch. No status has a source.**
- Evidence: the branch is 51 commits ahead of `origin/main`, with no upstream. On `origin/main`, `categories.yaml` has 10 calls. It has no `hydrophobin_domain`, `hsba_domain` or `surface_attachment_candidate` (`git show origin/main:src/cellsurface_sorting_hat/categories.yaml`). On `main` only PA14 is active (`git show origin/main:data/sorting_hat/family_table.tsv`).
- No cell of the "Status today" column cites a file.
- Fix: add a line that names the branch and commit the table describes. Add a source column. The sources are:
  - R0: `_workdir/sorting_hat/*/status/step1_rule@R0.json`. *A. nidulans* (taxon 162425) is `estimated`. Taxa 4932, 5476, 746128, 5207 and 5270 are `smoke`.
  - Repeat: `docs/reports/data/sorting_hat/call_status/{Scer_S288C,Calb_SC5314}/tandem_repeat_protein.owner_definition.config-with-surface-attachment.json`.
  - Hydrophobin: `docs/reports/data/sorting_hat/call_status/hydrophobin/pfam_hydrophobin.*.json` (8 files).
- The values in the table match these files. I checked 0.435 / 0.987 and 0.462 / 0.976, and the 8 hydrophobin proteomes.

**2. MINOR. A status holds for each work directory and taxon, not for each call. The latest stored S288C run shows the repeat call as `unvalidated`.**
- Evidence: `out_attach/report.md` ("Call calibration") lists `tandem_repeat_protein … stale: config differs`. All 25 repeat calls in that run are `unvalidated`.
- The status file was rewritten at 23:11:20, after the run. It now carries config `8c32064a…`, so a new run would use it.
- R0 status files exist only in 4 work directories: Af293 UniProt, Calb, S288C and S288C without dubious ORFs. A1163 and W72310 have no R0 status file, so R0 is `unvalidated` there.
- Fix: in the table, say "`smoke` where the status file is installed and current". Do not say "`smoke`" alone.

**3. MINOR. Every call is present, but some rows are wrong.**
- I compared the 13 calls in `categories.yaml` with the table. No call is missing and no call is extra. The kinds are right. These rows are wrong:
  - `iuis_allergen_similarity` reads only `allergen_homology`. The merged row lists `pfam_allergen` for both calls.
  - The meaning of `iuis_allergen_homolog` leaves out its Pfam branch (`categories.yaml` lines 85-91).
  - `signal_peptide_protein`, `other_not_surface` and `other_surface_no_mechanism` have `per_variant: true`. The table does not mark them "gated". The spec's own definition of "gated" applies to them.
  - "NOT CALIBRATED" is not a status value. The module status of `antigen_lookup` is `unvalidated`. "NOT CALIBRATED" is the label that the ranking prints.
- Fix: split the allergen row. Mark every `per_variant` call as gated. Write the status as "`unvalidated` (ranking prints NOT CALIBRATED)".

**4. MAJOR. In practice `iuis_allergen_homolog` never returns `not_called`.**
- Evidence: both allergen Pfam families (PF16541, PF25312) are inactive. So the module `pfam_allergen` is `unavailable` (`modules/cli.py` lines 286-300).
- The call is `or(test, flag pfam_allergen.hit)`. When the test is false, the result is `or(not_called, not_assessable)`, which is `not_assessable`.
- S288C `out_attach`: 15 `called`, 0 `not_called`, 6707 `not_assessable`. Ag2/PRA in RS is also `not_assessable`.
- This is noted only in passing in `docs/reports/2026-10-07-classifier-runs-and-scale.md` line 13.
- Fix: state this in the one-sentence meaning. Add it as an owner question: drop the Pfam branch until a family is active, or keep it and accept `not_assessable`.

**5. MAJOR. The meaning of `wall_family_domain` leaves out the second conditions. One of them makes this ungated call depend on R0.**
- Evidence: `family_table.tsv` gives CFEM `second_condition = no_tm` and PA14 `second_condition = signal_peptide`.
- The S288C module record `modules/pfam_adhesion.json` has `conditions.sp_module.module = step1_rule@R0`. So for a variant other than R0, the PA14 route is still gated by R0.
- Fix: the meaning should say "a hit to an active family that passes that family's condition (CFEM: no mature TM helix; PA14: R0 signal peptide)".

### Section 2: decisions C1 to C9

**6. MINOR. Most decisions have a document source. Two citations should be sharpened. One data file conflicts with C6.**

| # | Supported? | Where |
|---|---|---|
| C1 | yes | orchestrator spec line 449 (D2) and §1 lines 18-22 |
| C2 | yes | `docs/HANDOFF-2026-10-07.md` line 17 |
| C3 | yes | orchestrator spec line 458 (D11). Known limit 3 |
| C4 | yes | `docs/HANDOFF-2026-10-08.md` line 42 |
| C5 | yes against the data | `family_table.tsv` has 24 rows and 10 active: PF05730, PF07691, PF01185, PF06766, PF28987, PF22354, PF29785, PF29802, PF29465 (7 hydrophobin) and PF12296. Cerato-platanin is absent. Its only source is the roadmap line 168 ("PR #75 body"). I found no owner document |
| C6, C7 | yes | `docs/paper/02-training-and-testing-ledger.md` line 143 (H9). H9 says the narrow call is "unchanged". It does not say "not renamed" in words |
| C8 | yes | ledger line 142 (H8); `analysis/hydrophobin_truth/ship_decision.json` |
| C9 | yes | `docs/HANDOFF-2026-10-08-hydrophobin.md` line 27; M5 spec line 3 |

- Conflict with C6: the PF12296 row in `family_table.tsv` has class "2c hydrophobin (HsbA)" and the note "owner 2026-10-08: include as hydrophobins". Known limit 4 says HsbA "may not be a hydrophobin". C6 calls HsbA a separate surface-active call.
- The `source_pmid` cell of PF00624, PF13928, PF15789, PF22799 and PF30910 reads "owner decision 2026-10-08 (hydrophobin category by Pfam models)". These are flocculin, Hyr1, PIR and ALS rows.
- Fix: cite the ledger rows for C6 to C8. Cite the roadmap line for cerato-platanin, or ask the owner. In a data commit, correct the HsbA class text and the five `source_pmid` cells.

**7. MINOR. The O2 family list is wrong.**
- O2 names "Bys1, ALS, Flo11, GLEYA, Hyr1, PIR and the five candidates". Hyr1 (PF15789) is one of the five candidates, so the list counts it twice.
- The list leaves out Hyphal_reg_CWP (PF11765). It also leaves out the two inactive allergen families, PF16541 and PF25312, which finding 4 depends on.
- There are 14 inactive rows: 12 in `pfam_adhesion` and 2 in `pfam_allergen`.
- Fix: list the 14 accessions.

### Section 5: the changes

**8. MAJOR. The 5.1 acceptance check cannot be run as written.**
- `evidence.tsv.gz` is a long table with the columns `protein, module, field, value` (`outputs.py` line 123). It has no column for each field.
- The golden test pins only `calls.long` and `report.md` (`test_cli.py` lines 914-944). The toy `pfam_adhesion` rows have no `families` field (`test_cli.py` line 42).
- So "the golden file test covers it" is false. The only golden change will be the `categories.yaml sha256` line in `report.expected.md`.
- The mechanism itself works. My synthetic check gave `('cfem', 'pfam_adhesion', 'families', 'PF05730')`.
- Rows that are not `ok` are skipped. Rows with an empty field are skipped, which includes every `hit = 0` row (`engine.py` lines 514-526).
- A field name that the module does not write is skipped silently. `_validate` checks only the `MODULE.FIELD` form.
- The value is the accession `PF05730`, not the name "CFEM".
- Fix: change the check to "`evidence.tsv.gz` has rows with `field = families` for the three modules". Add `families` to the toy fixture. Add an assertion in `test_evidence_and_protein_tables`. Change the text "a reader can see CFEM" to "PF05730 (CFEM)", or add a name lookup.

**9. MINOR. Any edit to `categories.yaml`, a comment included, makes the repeat call files stale.**
- Evidence: the config hash is taken over the raw bytes (`engine.py` lines 84-92). A call file is stale on "config differs" (`call_status.py` line 179). `docs/paper/03` §6a says so.
- M2d is in the plan, but the spec does not explain why. It also does not say where the new files are stored.
- Fix:
  - State the reason in M2d.
  - Run M2d right after M2a.
  - Copy the new files to `docs/reports/data/sorting_hat/call_status/*/` with a new suffix.
  - Forbid edits to `categories.yaml` in M2b and M2c. The header comment on lines 1-2 points to the stale §3.4 and is tempting to edit.
  - Note that O1 (b) or (c), and M5, will require the re-measurement again.

**10. MINOR. 5.2: a separate docs yaml is the right place, but "measured where" will go stale and cannot be tested.**
- Keeping descriptions out of `categories.yaml` is correct, by finding 9. The spec should say this is the reason.
- The existing pattern compares spec text with the code (`test_outputs.py` lines 127-151, known limits).
- A docs yaml checked against `load_config()` names in both directions follows that pattern.
- But "measured where" in a static yaml repeats the status files. Those files sit in git-ignored work directories, so CI cannot check them.
- After M2b there are three texts that state meanings: the orchestrator §3.4 table, `_known_limits`, and `docs/CLASSES.md`.
- Fix:
  - "Measured where" should cite committed files, for example `docs/reports/data/sorting_hat/call_status/...`. The test should check that each cited path exists.
  - Name where the generator lives and confirm that CI runs the test.
  - Make the §3.4 table point to `docs/CLASSES.md`. Keep the "Known limits, printed in the report header" block and its 9 items, because `test_outputs.py` parses that block up to the next `### `.

**11. MAJOR. 5.3 names real stale statements but misses several. Its proposed source is itself stale.**
- Confirmed as stale:
  - `STATUS.md` line 26: orchestrator "proposal only, no spec".
  - `STATUS.md` line 22: "15 families, all inactive". There are 24 rows and 10 are active.
  - Orchestrator spec line 5: "No code, data or job exists".
  - `TOOL-ARCHITECTURE.md` lines 72 and 246: PR-AUC 0.94-0.98, a number of the old classifier.
  - `TOOL-ARCHITECTURE.md` line 77 ("solved by HMMs") and line 248 ("no work needed"). The PF01185 note in `family_table.tsv` records that RodD is missed. The hydrophobin reports record the T2 misses.
- Not named in 5.3:
  - Orchestrator §3.4 table: `wall_family_domain` lists PF01185, PF06766, PF28987, PF22354 and the ALS families. `tandem_repeat_protein` is "`unvalidated` in every clade". There are no rows for `hydrophobin_domain`, `hsba_domain` or `surface_attachment_candidate`. The mechanism list leaves out the two surface-active calls.
  - `categories.yaml` lines 1-2 point readers to that table.
  - `docs/paper/03` line 139: "No other module has a status above `unvalidated`".
  - `TOOL-ARCHITECTURE.md` line 238 needs a check in the same pass.
  - The roadmap, which 5.3.1 names as the source for `STATUS.md`, is stale. Lines 15-18 and 34 say PR #76 is open, but `origin/main` is `2d023d6`, "Merge pull request #76".
- Fix: add the §3.4 table and `docs/paper/03` line 139 to 5.3. Use git and the status files as the source, not the roadmap.

### Rule 3: composite status is derived

**12. MAJOR. The machinery derives the status. Rule 3 states it wrongly, and its exception cannot happen.**
- Check: `call_eligible` returns `(False, 'contains ref')` for `cell_wall_adhesion_candidate`, `surface_attachment_candidate` and `serodiagnostic_marker_candidate`.
- In my synthetic run with a `smoke` repeat call file and a `smoke` R0 status:
  - A repeat protein got `surface_attachment_candidate = smoke`, basis `callfile;step1_rule@R0:R0src`.
  - A hydrophobin protein got `smoke` from `pfam_hydrophobin` and R0. Its narrow call was `not_called` with status `unvalidated`.
  - A CFEM protein got `unvalidated` from `pfam_adhesion`.
- So the status is computed per record. It comes from the leaves that decided the value (`engine.py` lines 223-233). It does not come from all the leaves the call reads.
- Rule 3 says "the measured leaf calls it reads". `docs/paper/03` §6 decision 3 says "that decided it".
- Rule 3's exception ("unless it has its own status file") cannot occur. Call files are refused for any call with a `ref` (`docs/paper/03` §6a; `engine.py` lines 452-475).
- The report explains derived statuses only in the "Call calibration" section, and that section appears only when call files exist (`outputs.py` lines 231-239). Most runs, such as the Fungi_5k scan and RS, have no call file, so their reports never say this.
- Owner decision 4 of 2026-10-08 is open (`docs/HANDOFF-2026-10-08.md` line 44; roadmap D3 #16). It asks whether to put a "derived" marker in `status_basis` for every run. Section 6 of the spec does not list it.
- Fix: write rule 3 as "weakest status of the leaf calls that decided this record". Delete the exception, or say that it needs an engine change. Add owner decision 4 as O5.

**13. MAJOR. Rule 6 (status meanings in the report header) has no task. Section 5 says no other code changes.**
- Evidence: `render_report` prints no definition of `unvalidated`, `smoke` or `estimated` (`outputs.py` lines 175-303).
- Rule 6 needs code changes: `outputs.py`, the golden report, and probably the known-limits test, if it becomes a known limit.
- Fix: add a task (M2e) with a golden diff, or move rule 6 to `docs/CLASSES.md` only.

### O1, CFEM, and hydrophobin and HsbA wording

**14. MINOR. O1 is framed correctly. It needs numbers, and the default does not solve problem 2.**
- Code today: a CFEM-only protein with no mature TM helix and an R0 signal peptide is `wall_family_domain = called`. It is then `cell_wall_adhesion_candidate[R0] = called`.
- RS `out_hyd`: `XP_001240075.1` (Ag2/PRA) is called by both, by `cocci_specificity_rank_top15` and by `serodiagnostic_marker_candidate`.
- `wall_family_domain` hits in the module tables:

  | Proteome | CFEM | PA14 |
  |---|---|---|
  | S288C | 1 | 4 |
  | Calb | 6 | 0 |
  | Af293 UniProt | 4 | 0 |
  | RS | 7 | 0 |

- So 18 of the 22 domain hits are CFEM. With CFEM removed, by option (b) or (c), the domain route of the narrow call fires only in S288C.
- Option (b) names ALS and Flo. Both are inactive, so (b) means PA14 only today.
- Default (a) leaves §1 problem 2 as it is. Known limit 4 already says that CFEM is not shown to mediate adhesion.
- Fix: add these counts and the fact about (b). Say plainly that (a) only shows the family and does not change the call.

**15. MAJOR. The report text labels repeat evidence as "adhesin repeat". This breaks rule 2.**
- Evidence: `outputs.py` line 275 maps `tandem_repeat_protein` to "adhesin repeat". Line 276 maps `wall_family_domain` to "adhesion or wall family domain".
- Known limit 4 says repeat proteins include intracellular ones. The roadmap (lines 292-294, from the 2026-10-07 report §3) says the repeat call fires for 5 of 27 S288C hard negatives.
- Fix: rename the labels, for example "tandem repeat" and "wall family domain (PA14, CFEM)". Regenerate the golden report.

**16. MINOR. Other inconsistencies for hydrophobin and HsbA, and for the basis column.**
- The report table puts hydrophobin and HsbA in one row (`outputs.py` lines 277-278). A reader cannot see the HsbA count.
- The table counts records, not proteins (`outputs.py` line 272). If more than one step 1 variant is available, "proteins called" counts each protein once per variant. Today only R0 is available.
- The basis of `surface_attachment_candidate` is only in `calls.long`. `calls.wide` writes `_basis` only for `other_*` (`outputs.py` lines 111-117; the S288C `out_attach` wide header has no `surface_attachment_candidate[R0]_basis`). Spec §1 item 3 does not say this.
- The HsbA class text in `family_table.tsv`: see finding 6.
- The golden files contain the same text as `outputs.py`, so they inherit these labels.
- Fix: split the HsbA row. Count distinct proteins for each variant. Either write the basis to the wide table or say that it is only in `calls.long`.

### Missing items and plan order

**17. MINOR. Gaps in the class list.**
- `hydrophobin_relaxed` and `cys8_pattern` run and appear in the report's module tables (S288C `out_attach`), but no call reads them.
- `tm`, `cys_rich` and `expression` are evidence fields.
- The class list should list these as "evidence or experimental, not a call". Otherwise `docs/CLASSES.md`, which is generated from calls, does not explain them.
- The roadmap M2 row also names Bys1, ALS and SOWgp-like. These are deferred to O2/M4 and M5, which is acceptable.

**18. MINOR. The M5 call name contradicts the M5 spec.**
- Section 3 says M5 adds `secreted_pro_cys_array`. The M5 spec names it `pro_cys_array_protein` (M5 spec lines 41, 50 and 69).
- O4 of M2 overlaps O1 of M5.
- Fix: use the M5 name. Point O4 to M5 O1.

**19. MINOR. Rule 4 cannot be tested.**
- "Or name the measure" lets any name pass.
- Fix: list the allowed suffixes, plus the three named exceptions (`cocci_specificity_rank_top15` and the two `iuis_allergen_*` calls). Then the M2b test can check the names.

**20. MINOR. The M2c check is not defined.**
- "A script lists the replaced statements" names no script, input or output.
- Fix: make the check a table in the commit message or in the PR body. It has one row per replaced line: file, line, old text, new text and supporting file. A reviewer reads the table. No script is needed.

## 2. Verdict

The spec is not ready for code as written. There is no blocker. The design is sound:
- The evidence column mechanism works.
- A separate docs table is the right choice.
- Composite statuses are already derived.

The MAJOR findings need text changes before M2a starts:
- 1: source and branch for the status column.
- 4 and 5: correct meanings for the allergen homolog and wall-domain calls.
- 8: a 5.1 check that can be run.
- 11: the full list of stale statements.
- 12: rule 3 wording, and owner decision 4 added.
- 13: a task for rule 6.
- 15: the "adhesin repeat" label.

M2a can start once finding 8 is fixed and M2d is placed right after it.
