# M2: the v1 class list of `cellsurface_sorting_hat` and how its calls combine

*2026-10-09. Revision 2 (after independent review 1: 0 blockers, 8 major, 12 minor findings, all addressed; the review is `2026-10-09-v1-class-list-design-review-1.md`).
Author: Claude Code (claude-sonnet-5-5). Milestone M2 of `docs/ROADMAP-2026-10-08-DRAFT.md`. Branch `fungi-scan-and-specs`. This spec fixes what the tool reports as a class, what each class
means, and how the calls combine. It also lists the small code and documentation changes that follow. It makes no new scientific claim.*

**What this describes.** The table in section 3 describes branch `fungi-scan-and-specs` (built on `hydrophobin-validation`, PR #77, which contains `main` at `2d023d6`). On `origin/main` itself the calls `hydrophobin_domain`,
`hsba_domain` and `surface_attachment_candidate` do not exist and only PA14 is an active family. The text is written for the state after PR #77 merges.

## 1. Why

1. No single current document says what the tool reports and what each status means. `docs/model-review/STATUS.md` (2026-10-02), the orchestrator spec header ("no code, data or job exists") and `docs/TOOL-ARCHITECTURE.md` are stale (section 5.5).
2. `cell_wall_adhesion_candidate` is called for CFEM proteins (Ag2/PRA `XP_001240075.1` in RS) through `wall_family_domain`. A CFEM domain is a wall family and is not shown to mediate adhesion (known limit 4 says so).
   Domain hits that set `wall_family_domain` today: S288C CFEM 1 and PA14 4; *C. albicans* CFEM 6; Af293 CFEM 4; RS CFEM 7. So 18 of 22 are CFEM, and the domain route of the narrow call would fire only in S288C without CFEM.
3. A reader cannot see which evidence made a composite call. The `basis` of `surface_attachment_candidate` exists, but only in `calls.long`; `calls.wide` writes a basis column only for `other_*` calls.
4. The report's own labels overstate: `outputs.py` calls repeat evidence "adhesin repeat", although the repeat call also fires for intracellular proteins (known limit 4; 5 of 27 S288C hard negatives).

## 2. Decisions already made (owner, 2026-10-04 to 2026-10-08)

| # | Decision | Source |
|---|---|---|
| C1 | Five categories in v1: surface glycoprotein, cell wall and adhesion candidate, antigen, allergen, other. Enzyme classes later. | orchestrator spec D2 (line 449) and section 1 |
| C2 | Step 1 gate is rule R0 (SignalP 6 signal peptide). ML step 1 is deferred. | `docs/HANDOFF-2026-10-07.md` |
| C3 | No `cell_wall_protein` call in v1. GPI-anchored wall proteins are not a class in v1. | orchestrator spec D11, known limit 3 |
| C4 | `tandem_repeat_protein` means a repeating motif or array as in FLO11. A repeating domain counts only if it is a known adhesion-associated domain. | `docs/HANDOFF-2026-10-08.md` |
| C5 | Active Pfam families (10 of 24 rows in `data/sorting_hat/family_table.tsv`): PA14, CFEM, the seven hydrophobin-class models, HsbA. All others inactive until signed off. Cerato-platanin is out of scope (it is absent from the table; the owner's statement is in the conversation of 2026-10-08 and the PR #77 body; no other document). | table; `docs/paper/02` rows H1 to H8 |
| C6 | Hydrophobins and HsbA are surface-active proteins: separate evidence calls, in the broader call, not adhesins. HsbA is kept together with hydrophobin "for now" (it may not be a hydrophobin). | `docs/paper/02` row H9 |
| C7 | The broader call `surface_attachment_candidate` sits beside the narrow `cell_wall_adhesion_candidate`. The narrow call is not renamed. | row H9 |
| C8 | The hydrophobin call is the strict Pfam call only. The relaxed level is not brought in. | row H8; `analysis/hydrophobin_truth/ship_decision.json` |
| C9 | SOWgp-like proteins get both an ortholog call and a repeat-architecture call (M5). | `docs/HANDOFF-2026-10-08-hydrophobin.md`; M5 spec |

## 3. The v1 class list

"Gated" means `per_variant: true` (one value per step 1 variant). Statuses are per work directory and taxon: a status counts only where its file is installed and current.
Status sources are committed files where they exist.

| Call | Kind | Gated | Category | Meaning in one sentence | Reads | Status and source |
|---|---|---|---|---|---|---|
| `signal_peptide_protein` | evidence | yes | surface glycoprotein | SignalP 6 predicts a signal peptide (rule R0); nothing more is claimed. | `step1_rule@R0` | `estimated` in *A. nidulans* (taxon 162425); `smoke` in taxa 4932, 5476, 746128, 5207, 5270; `unvalidated` elsewhere. Files: `_workdir/sorting_hat/*/status/step1_rule@R0.json` (installed in 4 work directories; not committed). Report: `docs/reports/2026-10-07-sorting-hat-run-and-calibration.md` |
| `tandem_repeat_protein` | evidence | no | cell wall and adhesion | One of two period detectors finds a tandem array (period above 0, coverage 0.25 or more, 2.5 copies or more). It also fires for intracellular repeat proteins. | `repeat02`, `repeat14` | `smoke` where the call file is current: S288C 0.435 [0.067, 0.667] / 0.987; *C. albicans* 0.462 [0.000, 0.788] / 0.976. Files: `docs/reports/data/sorting_hat/call_status/{Scer_S288C,Calb_SC5314}/tandem_repeat_protein.owner_definition.config-with-surface-attachment.json`. Stale after any `categories.yaml` edit |
| `wall_family_domain` | evidence | no | cell wall and adhesion | A hit to an active wall-family Pfam model that passes that family's condition (CFEM: no mature TM helix; PA14: R0 signal peptide, which makes this ungated call depend on R0). | `pfam_adhesion` | `unvalidated` |
| `hydrophobin_domain` | evidence | no | cell wall and adhesion (surface-active) | A hit to one of seven hydrophobin-class Pfam models at the gathering cutoff. | `pfam_hydrophobin` | `smoke` in 8 proteomes: `docs/reports/data/sorting_hat/call_status/hydrophobin/pfam_hydrophobin.*.json` |
| `hsba_domain` | evidence | no | cell wall and adhesion (surface-active) | A hit to HsbA (PF12296). | `pfam_hsba` | `unvalidated` (no truth) |
| `cell_wall_adhesion_candidate` | composite | yes | cell wall and adhesion | Signal peptide and (repeat or wall-family domain). | the three calls above through `ref` | derived (rule 3) |
| `surface_attachment_candidate` | composite | yes | cell wall and adhesion | Signal peptide and any of repeat, wall-family domain, hydrophobin, HsbA; the basis column names what held. | as above plus the two surface-active calls | derived (rule 3) |
| `cocci_specificity_rank_top15` | evidence | no | antigen | Top 15% of the *Coccidioides* genus-specificity ranking (RS only). | `antigen_lookup` | `unvalidated` (the ranking prints NOT CALIBRATED) |
| `serodiagnostic_marker_candidate` | composite | yes | antigen | Top 15% ranking and signal peptide. | as above | derived |
| `iuis_allergen_similarity` | evidence | no | allergen | Similarity to a WHO/IUIS fungal allergen of at least 35% over 80 aa. | `allergen_homology` | `unvalidated` |
| `iuis_allergen_homolog` | evidence | no | allergen | At least 70% identity over 80% of an allergen, **or** a hit to an allergen Pfam family. Both allergen families (PF16541, PF25312) are inactive, so the Pfam branch is `not_assessable` and the call is never `not_called` in practice (S288C: 15 called, 0 not_called, 6707 not_assessable). | `allergen_homology`, `pfam_allergen` | `unvalidated` |
| `other_not_surface` | other | yes | other | No signal peptide, with the mechanism calls not called. | R0 and the mechanism calls | derived |
| `other_surface_no_mechanism` | other | yes | other | Signal peptide and no mechanism evidence (repeat, wall family, hydrophobin, HsbA, antigen rank). | as above | derived |

Evidence and experimental modules that no call reads: `tm` (a condition for CFEM), `cys_rich` and `expression` (evidence fields), `cys8_pattern` and `hydrophobin_relaxed` (experimental, not part of any call; C8). They are listed in `docs/CLASSES.md` as "evidence or experimental, not a call".

Added by M5 (not part of this spec): `sowgp_ortholog` and `pro_cys_array_protein`.

## 4. Rules for the class list

1. **One definition per call.** Every call has a one-sentence meaning, the modules it reads, the category and where its status is measured. They live in `docs/class_descriptions.yaml`, outside `categories.yaml` (an edit of `categories.yaml`, a comment included, changes its hash and makes the repeat call files stale). `docs/CLASSES.md` is generated from the two files.
   A test compares the call names of `load_config()` with the names in the yaml in both directions (the pattern of the known-limits test in `tests/cellsurface_sorting_hat/test_outputs.py`).
2. **Evidence is not a claim.** Evidence calls say what a tool found. A composite says what a rule combined. The report and `docs/CLASSES.md` use "candidate" for composites. The report labels for the sub-classes use "tandem repeat", "wall family domain (PA14, CFEM)", "hydrophobin" and "HsbA", not "adhesin".
3. **Composite status is derived per record.** The status of a composite record is the weakest status of the leaf calls that decided that record's value (the existing machinery; `call_eligible` refuses a call file for any call with a `ref`, so a composite cannot have its own status file; giving it one would need an engine change, which this spec does not make).
4. **Names.** Evidence calls end in `_domain`, `_protein`, `_ortholog` or `_similarity`/`_homolog` (allergen). Composites end in `_candidate`. Residual categories start with `other_`. Named exceptions: `cocci_specificity_rank_top15`. The M2c test checks these suffixes.
5. **Adding a call** needs an entry in `docs/class_descriptions.yaml`, a module that writes its input, tests, a golden diff read line by line, and a statement of how it can be measured (or that it cannot).
6. **What a status means** is printed once in the report header (task M2e): `unvalidated` (no measurement), `smoke` (measured on too little truth for a number), `estimated` (enough truth); a leakage cap applies; a status counts only for the species and taxon of its file.

## 5. Changes this spec makes

### 5.1 Evidence rows for the family

`evidence.tsv.gz` is a long table (`protein, module, field, value`). Add `pfam_adhesion.families`, `pfam_hydrophobin.families` and `pfam_hsba.families` to the `evidence:` list of `categories.yaml`. A protein with a hit then has a row `(module, families, PF05730)` etc.
Rows with an empty field (every `hit = 0` row) and modules whose state is not `ok` are skipped by the existing mechanism, and a field name that a module does not write is skipped silently. The value is the accession; `docs/CLASSES.md` and the report give the name lookup (PF05730 CFEM, PF07691 PA14). This shows the family. It does not change any call (option O1 (a)).
**Check:** the toy `pfam_adhesion` rows in `tests/cellsurface_sorting_hat/test_cli.py` get a `families` field, `test_evidence_and_protein_tables` asserts the new rows, and a run on the toy proteome has `field = families` rows for the three modules. The golden files pin only `calls.long` and `report.md`, so the golden change is the `categories.yaml sha256` line.

### 5.2 Report changes

1. Rename the sub-class labels in the surface-attachment table (`outputs.py`): "tandem repeat", "wall family domain (PA14, CFEM)", "hydrophobin", "HsbA" (two rows, not one).
2. Count distinct proteins per variant (today it counts records).
3. Write the basis of `surface_attachment_candidate` to `calls.wide` too (the existing code writes a `_basis` column only for `other_*`).
4. Print the status definitions of rule 6 and one sentence that composite statuses are derived (task M2e). Regenerate the golden report, reading every changed line.

### 5.3 Data fixes in `family_table.tsv` (notes only; the digest covers accession, module, condition and active, so no module identity changes)

1. The PF12296 row says "2c hydrophobin (HsbA)" and "include as hydrophobins". Change the class text to "2c surface-active (HsbA, grouped with hydrophobin)" and keep the owner note.
2. The `source_pmid` cell of PF00624, PF13928, PF15789, PF22799 and PF30910 carries the text "owner decision 2026-10-08 (hydrophobin category by Pfam models)". Replace it with the correct source for each (the owner's classification of 2026-10-08 as flocculin, Hyr1, PIR1-like and ALS_M; no PMID).

### 5.4 Generated class table

`docs/class_descriptions.yaml` (call, meaning, category, reads, status files as committed paths) and `analysis/class_table/make_classes_md.py` produce `docs/CLASSES.md`. "Measured where" cites committed files only (for example `docs/reports/data/sorting_hat/call_status/...`); the test checks that every cited path exists. The orchestrator §3.4 table gets a pointer to `docs/CLASSES.md`; its "Known limits, printed in the report header" block (9 items) stays, because `test_outputs.py` parses it up to the next `### `. CI runs the test (`.github/workflows`; checked in M2c).

### 5.5 Documentation refresh (sources: git history and status files, not the roadmap, which still says PR #76 is open)

| File | Statement to replace |
|---|---|
| `docs/model-review/STATUS.md` line 26 | "Orchestrator: proposal only, no spec" |
| `docs/model-review/STATUS.md` line 22 | "15 families, all inactive" (24 rows, 10 active) |
| orchestrator spec line 5 | "No code, data or job exists for this spec" |
| orchestrator spec §3.4 table | hydrophobins listed in `wall_family_domain`; `tandem_repeat_protein` "unvalidated in every clade"; no rows for `hydrophobin_domain`, `hsba_domain`, `surface_attachment_candidate`; mechanism list lacks the two surface-active calls |
| `docs/TOOL-ARCHITECTURE.md` lines 72 and 246 | PR-AUC 0.94 to 0.98 (a number of the old classifier) |
| `docs/TOOL-ARCHITECTURE.md` lines 77, 248 and 238 | hydrophobin "solved by HMMs" / "no work needed" (T2 misses recorded in the hydrophobin reports; line 238 to be read in the same pass) |
| `docs/paper/03` line 139 | "No other module has a status above `unvalidated`" |
| `categories.yaml` lines 1 and 2 | comment pointing at the stale table (**not edited in M2b to M2e**; see M2d) |

The refresh is recorded as a table in the commit message: file, line, old text, new text, supporting file. A reviewer reads the table.

## 6. Decisions for the owner (stop points)

| # | Question | Default |
|---|---|---|
| O1 | CFEM and the narrow call. A CFEM hit sets `cell_wall_adhesion_candidate`. (a) keep, with the family shown (5.1): this shows the family but does not change the call; (b) split `wall_family_domain` into an adhesion-domain call (PA14 today; ALS and Flo are inactive) and a wall-family call (CFEM), so the narrow call reads only the first; (c) remove CFEM from the narrow call and keep it in the broader call only. Under (b) or (c) the domain route of the narrow call fires only for PA14 (S288C, 4 proteins in the four proteomes counted). | (a) |
| O2 | Which of the 14 inactive families are signed off for v1 (milestone M4): Bys1 PF04681, Candida_ALS_N PF11766, Candida_ALS PF05792, Flo11 PF10182, GLEYA PF10528, Hyphal_reg_CWP PF11765, PIR PF00399, Flocculin PF00624, Flocculin_t3 PF13928, Hyr1 PF15789, PIR1-like_C PF22799, ALS_M PF30910, and the allergen families AltA1 PF16541 and Allergen_Asp_f_4 PF25312. | stay inactive |
| O3 | `iuis_allergen_homolog` is never `not_called` while the allergen families are inactive. Drop its Pfam branch until a family is active, or keep it and accept `not_assessable`. Also whether the antigen and allergen categories stay in the v1 report with status `unvalidated`. | keep, status printed |
| O4 | Does `pro_cys_array_protein` (M5) appear in `surface_attachment_candidate`? Decided in the M5 spec (its O1). | see M5 |
| O5 | A marker in `status_basis` that says the status is derived, for every run (owner decision 4 of 2026-10-08, open). | not added |

## 7. Not in this spec

New calls (M5), new families (M4), measured statuses (other milestones), the repeat detector (M6), Onygenales truth (M7), enzyme classes, GPI as a class, an engine change that lets a composite have its own status file.

## 8. Work plan (order matters: any edit of `categories.yaml` makes the repeat call files stale)

| Task | Content | Check |
|---|---|---|
| M2a | Evidence rows (5.1): `categories.yaml` edit, toy fixture and test, golden `sha256` line read | `pytest tests/cellsurface_sorting_hat`; a toy run has `field = families` rows |
| M2d | **Right after M2a:** re-measure the repeat call files for S288C and *C. albicans* under the new config; copy the files to `docs/reports/data/sorting_hat/call_status/*/` with a new suffix | numbers unchanged (0.435 / 0.987; 0.462 / 0.976) |
| M2b | Report changes (5.2 items 1 to 3) and the data fixes (5.3); golden report regenerated line by line. No edit of `categories.yaml`. | tests; golden diff read |
| M2e | Status definitions in the report header (5.2 item 4). No edit of `categories.yaml`. | tests; golden diff read |
| M2c | `docs/class_descriptions.yaml`, generator, test, `docs/CLASSES.md` (5.4) | the test passes in both directions; cited paths exist |
| M2f | Documentation refresh (5.5) with the table in the commit message | table complete |

Each task is its own commit. M5 and O1 (b) or (c) require the repeat re-measurement again. This revision has not been reviewed again; the review found no blocker, and the fixes are text.
