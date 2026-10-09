# M2: the v1 class list of `cellsurface_sorting_hat` and how its calls combine

*2026-10-09. Revision 1. Author: Claude Code (claude-sonnet-5-5). Milestone M2 of `docs/ROADMAP-2026-10-08-DRAFT.md`. Branch `fungi-scan-and-specs`. Status: draft for independent review. This spec
fixes what the tool reports as a class, what each class means, and how the calls combine. It also lists the small code and documentation changes that follow. It makes no new scientific claim.*

## 1. Why

The tool has grown calls one at a time (R0, repeats, Pfam families, hydrophobin, HsbA, antigen, allergen, composites). Three problems follow (roadmap section E):

1. No single current document says what the tool reports and what each status means. `docs/model-review/STATUS.md` (2026-10-02), the orchestrator spec header ("no code exists") and `docs/TOOL-ARCHITECTURE.md` are stale.
2. The composite `cell_wall_adhesion_candidate` is called for CFEM proteins (Ag2/PRA, PRA2) through `wall_family_domain`. A CFEM domain is a wall family and is not shown to mediate adhesion.
3. A reader cannot see which evidence made a composite call, except for the new `surface_attachment_candidate` basis column.

## 2. Decisions already made (owner, 2026-10-04 to 2026-10-08)

| # | Decision | Source |
|---|---|---|
| C1 | Five categories in v1: surface glycoprotein, cell wall and adhesion candidate, antigen, allergen, other. Enzyme classes later. | orchestrator spec D2 |
| C2 | Step 1 gate is rule R0 (SignalP 6 signal peptide). ML step 1 is deferred. | handoff 2026-10-07 |
| C3 | No `cell_wall_protein` call in v1. GPI-anchored wall proteins are not a class in v1. | orchestrator spec D11 |
| C4 | `tandem_repeat_protein` means a repeating motif or array as in FLO11. A repeating domain counts only if it is a known adhesion-associated domain. | owner, 2026-10-08 |
| C5 | Active Pfam families: PA14, CFEM, the seven hydrophobin-class models, HsbA. Cerato-platanin is out of scope. All other families are inactive until signed off. | owner, 2026-10-07 and 2026-10-08 |
| C6 | Hydrophobins and HsbA are surface-active proteins: separate evidence calls, in the broader call, not adhesins. | owner, 2026-10-08 |
| C7 | The broader call `surface_attachment_candidate` sits beside the narrow `cell_wall_adhesion_candidate`. The narrow call is not renamed. | owner, 2026-10-08 |
| C8 | The hydrophobin call is the strict Pfam call only. The relaxed level is not brought in. | owner, 2026-10-08 |
| C9 | SOWgp-like proteins get both an ortholog call and a repeat-architecture call (specified in M5). | owner, 2026-10-08 |

## 3. The v1 class list

Each row is a call in `src/cellsurface_sorting_hat/categories.yaml`. "Kind" is evidence (one tool), composite (a rule over calls) or other (the residual categories). "Gated" means it has one value per step 1 variant.

| Call | Kind | Category | Meaning in one sentence | Reads | Status today |
|---|---|---|---|---|---|
| `signal_peptide_protein` | evidence | surface glycoprotein | SignalP 6 predicts a signal peptide (rule R0). Nothing more is claimed. | `step1_rule@R0` | measured: `estimated` in *A. nidulans*, `smoke` in five species |
| `tandem_repeat_protein` | evidence | cell wall and adhesion | One of two period detectors finds a tandem array. | `repeat02`, `repeat14` | `smoke` in S288C and *C. albicans* (sensitivity 0.44 to 0.46) |
| `wall_family_domain` | evidence | cell wall and adhesion | A hit to an active wall-family Pfam model (PA14, CFEM). | `pfam_adhesion` | `unvalidated` |
| `hydrophobin_domain` | evidence | cell wall and adhesion (surface-active) | A hit to a hydrophobin-class Pfam model. | `pfam_hydrophobin` | `smoke` in 8 proteomes |
| `hsba_domain` | evidence | cell wall and adhesion (surface-active) | A hit to HsbA (PF12296). | `pfam_hsba` | `unvalidated` |
| `cell_wall_adhesion_candidate` | composite, gated | cell wall and adhesion | Signal peptide and (repeat or wall-family domain). | `signal_peptide_protein`, `tandem_repeat_protein`, `wall_family_domain` | derived (weakest leaf) |
| `surface_attachment_candidate` | composite, gated | cell wall and adhesion | Signal peptide and any of repeat, wall-family domain, hydrophobin, HsbA. The basis column names what held. | as above plus the two surface-active calls | derived |
| `cocci_specificity_rank_top15` | evidence | antigen | Top 15% of the *Coccidioides* genus-specificity ranking (RS only). | `antigen_lookup` | NOT CALIBRATED |
| `serodiagnostic_marker_candidate` | composite, gated | antigen | Top 15% ranking and signal peptide. | as above | derived |
| `iuis_allergen_similarity`, `iuis_allergen_homolog` | evidence | allergen | Similarity to a WHO/IUIS fungal allergen (35% over 80 aa; 70% over 80% of the allergen). | `allergen_homology`, `pfam_allergen` | `unvalidated` |
| `other_not_surface`, `other_surface_no_mechanism` | other | other | No signal peptide; or signal peptide and no mechanism evidence. | R0 and the mechanism calls | derived |

Added by M5 (not part of this spec): `sowgp_ortholog` and `secreted_pro_cys_array`.

## 4. Rules for the class list

1. **One definition per call.** Every call has a one-sentence meaning, the modules it reads, the category it belongs to and where its status is measured. These live in one generated file, `docs/CLASSES.md`, built from `categories.yaml` plus a small docs table (section 5.2). A test fails if a call has no entry or an entry has no call.
2. **Evidence is not a claim.** Evidence calls say what a tool found. A composite says what a rule combined. The report and `docs/CLASSES.md` use "candidate" for composites and "domain" or "repeat" for evidence.
3. **Composite status is derived.** A composite's status is the weakest status of the measured leaf calls it reads (per-call status machinery). No composite gets a status of its own unless it has its own status file and truth.
4. **Names.** Evidence calls end in `_domain`, `_protein`, or name the measure. Composites end in `_candidate`. The residual categories start with `other_`.
5. **Adding a call** needs: an entry in `docs/CLASSES.md`, a module that writes its input, tests, a golden-file diff read line by line, and a statement of how it can be measured (or that it cannot).
6. **What a status means** is stated once in the report header: `unvalidated` (no measurement), `smoke` (measured on too little truth for a number), `estimated` (measured with enough truth). The leakage cap applies.

## 5. Changes this spec makes (no owner decision needed)

### 5.1 Evidence columns: which family hit

`calls.long` cannot say whether `wall_family_domain` came from PA14 or CFEM. `evidence.tsv.gz` can copy a module field. Add `pfam_adhesion.families`, `pfam_hydrophobin.families` and `pfam_hsba.families` to the `evidence:` list of `categories.yaml`. A reader can then see "CFEM" next to an Ag2/PRA row without a new call. Check: `evidence.tsv.gz` has the columns; the golden file test covers it.

### 5.2 A generated class table

`docs/CLASSES.md` is generated by a small script from `categories.yaml` and `docs/class_descriptions.yaml` (call, meaning, category, measured where). A test checks that the file is current and complete. The script and the test are the only code. No call changes.

### 5.3 Documentation refresh

1. `docs/model-review/STATUS.md`: rewrite the head (what exists, what is measured, what is not) from the roadmap table, keeping the 2026-10-02 text as history.
2. `docs/superpowers/specs/2026-10-04-orchestrator-design.md`: header says the code exists, names the tests and the version of the call table.
3. `docs/TOOL-ARCHITECTURE.md`: replace the statements that code and data show to be older (hydrophobin "solved"; the old classifier's repeat PR-AUC 0.94 to 0.98; the family list). Each replaced statement is listed in the commit message with the file that supports the new text.
4. `docs/paper/03` section 7 and the paper notes: update only what the code now shows.

## 6. Decisions for the owner (stop points)

| # | Question | Default if the owner does not answer |
|---|---|---|
| O1 | CFEM and the narrow call. A CFEM hit sets `cell_wall_adhesion_candidate` although adhesion is not shown. Options: (a) keep, with the family shown in evidence (5.1); (b) split `wall_family_domain` into an adhesion-domain call (PA14, ALS, Flo) and a wall-family call (CFEM), so the narrow call reads only the first; (c) remove CFEM from the narrow call and keep it in the broader call only. | (a) |
| O2 | Whether the inactive families (Bys1, ALS, Flo11, GLEYA, Hyr1, PIR and the five candidates) are signed off for v1. This is milestone M4. | stay inactive |
| O3 | Whether allergen and antigen categories stay in the v1 report with their current status (`unvalidated`, NOT CALIBRATED). | stay, with the status printed |
| O4 | Whether `surface_attachment_candidate` should also read the SOWgp-like calls of M5 (SOWgp is a spherule wall protein; its adhesion role is not part of this spec). | decided in the M5 spec |

## 7. Not in this spec

New calls (M5), new families (M4), measured statuses (other milestones), the repeat detector (M6), Onygenales truth (M7), enzyme classes, GPI as a class.

## 8. Work plan

| Task | Content | Check |
|---|---|---|
| M2a | Evidence columns (5.1): config edit, golden regeneration read line by line, tests | `pytest tests/cellsurface_sorting_hat`; `evidence.tsv.gz` has the three columns |
| M2b | `docs/class_descriptions.yaml`, the generator and its test (5.2); `docs/CLASSES.md` | the test passes; the file lists every call |
| M2c | Documentation refresh (5.3) | a script lists the replaced statements and the supporting file for each |
| M2d | Re-measure the repeat call files after the last `categories.yaml` edit | numbers unchanged (0.435 / 0.987; 0.462 / 0.976) |

Each task is its own commit. An independent review precedes M2a.
