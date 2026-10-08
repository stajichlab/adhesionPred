# Review 1: hydrophobin call validation and 8-cysteine rescue rule (design of 2026-10-08)

*2026-10-08. Independent review by Claude Code (claude-opus-5-5). Reviewed:
`docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md` (revision 1, commit d465e4c, branch
`hydrophobin-validation`). The spec was not edited. Every statement below was checked against a file or a
command output. Where I did not check something, I say so.*

## How I checked

- Code on this branch: `categories.yaml`, `modules/pfam.py`, `modules/signalp.py`, `modules/cli.py`,
  `engine.py`, `outputs.py`, tests.
- Per-call status code is **not on this branch**. I read it from local branch `per-call-status` (PR #76, open):
  `engine.py` (`call_eligible`, `call_hash`), `call_status.py`, `docs/paper/03` section 6a,
  `docs/reports/2026-10-08-repeat-call-calibration.md`.
- Re-derived F2 and F3 from `analysis/hydrophobin_truth/` with Python.
- Re-ran `analysis/sorting_hat_run/hydrophobin_discovery/analyse.py` (read-only) on the domain tables in
  `_workdir/sorting_hat/calibration/hydrophobin_discovery/`, and counted hits per proteome with `awk`.
- Fetched PubMed metadata for the seven PMIDs of F7 and `docs/paper/05` section 5.

## Findings

### 1. BLOCKER: the `variants` syntax does not exist, and `pfam_or_cys8` as written cannot get a call status file

Evidence:
- `categories.yaml` has no `variants` key. A call has `name`, `expr`, optional `per_variant`, or `kind: other`.
  `engine._validate` reads `call["expr"]` for every non-`other` call (`engine.py` on `per-call-status`, line
  121). A call with only `variants` fails with a `KeyError`, not a `ConfigError`.
- `per_variant` means the step 1 variant (R0, R1, R2, card). It is not a rule-variant mechanism.
- `pfam_or_cys8` uses `{ref: hydrophobin_domain}` and `{ref: signal_peptide_protein}`. `call_eligible` returns
  `(False, "contains ref")` for any call with a `ref` node (`engine.py` on `per-call-status`, lines 443-466;
  `docs/paper/03` section 6a: "has no `ref`").
- `signal_peptide_protein` is `per_variant: true`. A `ref` to it from an ungated call raises
  `ConfigError: an ungated call cannot use per_variant call` (`engine.py` line 144-145).
- So the only call that the design gives a call status file cannot have one.

Fix:
- Write two named calls. Do not use `variants`.
- Option A: `hydrophobin_protein` with `per_variant: true` and an inline expression
  `or: [{flag: pfam_hydrophobin.hit}, {and: [{flag: cys8_pattern.hit}, {call: "{step1}"}]}]`. It reads three
  modules for each variant and is eligible. It cannot be a `mechanism` of the `other_*` calls (mechanisms must
  be ungated, `engine.py` line 117-119). The spec must say which step 1 variants are measured (one file per
  variant: `status/calls/<call>.<variant>.json`).
- Option B: put the signal-peptide condition inside `cys8_pattern` (as `pfam_adhesion` does with
  `second_condition`). Then `or: [{flag: pfam_hydrophobin.hit}, {flag: cys8_pattern.hit}]` is ungated, reads two
  modules, is eligible, and can be a mechanism. Its status then depends on the R0 module identity.
- Choose one and state it in section 3.3.

### 2. BLOCKER: `cys8_pattern` cannot cut the signal peptide from the step 1 module table

Evidence:
- `modules/signalp.py` keeps only `prediction` and `sp_prob` (line 46). `COLUMNS = ["call", "sp_prob",
  "prediction"]` (line 75). The cleavage site is read from the file but dropped.
- `_condition_table` reads one named column of a module table (`modules/cli.py` lines 360-390). No column holds
  the cleavage position.
- Only `step1_rule@R0` comes from SignalP. The design does not say what `cys8_pattern` does for R1, R2 or
  `step1_ml@card`.

Fix:
- Add a `cs_pos` column to the SignalP module, or read `prediction_results.txt` in the `cys8` wrapper and record
  its hash as an artefact.
- If the SignalP module changes, state the effect on the R0 status (A. nidulans `estimated`, decision 4 in
  `docs/paper/03` section 6). A version bump drops existing R0 statuses.
- Or drop the cut. Run the pattern on the full sequence and keep the SP requirement as a separate leg. The
  regex does not need the mature region. Only `n_cys` and `mature_length` need it.

### 3. BLOCKER: the negative-set rule makes Pfam-only specificity 1.0 by construction

Evidence:
- Section 4.3 item 3: "a hit to any of the eight Pfam models at the gathering cutoff is not allowed in the
  negative set of the Pfam-only variant."
- Every Pfam-only call is a Pfam hit. If Pfam hits are removed from the negatives, the Pfam-only variant has no
  false positives by definition.
- In Af293 this removes 9 HsbA proteins and any non-hydrophobin Pfam hit (finding 7).
- The `pfam_or_cys8` variant uses a different negative set. The "added FP" in section 6 then compares two
  denominators.

Fix:
- Labels must not depend on any prediction. One negative set for both variants.
- A Pfam hit with no label of any tier is a negative for both variants. Report it as a candidate false positive
  and review it by hand (as section 6 already does for rescue-only calls).

### 4. BLOCKER: in the proteomes in scope, the rescue gain cannot be measured

Evidence:
- Section 11 limits the work to proteomes in `_workdir/sorting_hat/`: Afum Af293, A1163, W72310, Bder ER3,
  Calb SC5314, Cimm RS, Cneo H99, Scer S288C. That is five species.
- Of the 174 Swiss-Prot hydrophobins, the only entries from these species are RodA to RodG of A. fumigatus
  (P41746, E9QT94, Q4WBR8, Q4WE22, Q4WC31, Q4WEK0, Q4X055). All 7 are in `Afum_Af293_UniProt.faa`. All 7 carry
  a hydrophobin-class Pfam model in the UniProt cross-reference (RodD by PF29785).
- None of the 9 proteins without a Pfam model is from these species. They are from P. expansum, T. asperellum,
  H. virens, G. zeae (2), F. fulva (3) and F. velutipes.
- The discovery search found 3 secreted, motif-positive, Pfam-negative proteins in Af293. All 3 are annotated
  "GPI anchored protein" (Q4WUX0, Q4WBD4, Q4WBC1). Output of `analyse.py`.
- So with Swiss-Prot truth in Af293: added TP = 0, added FP >= 0. The keep criterion of section 6 can only
  reject the rescue there. This is decided by the scope, not by the measurement.
- Section 4.4 lists A. nidulans, a Trichoderma species, P. ostreatus and B. bassiana as candidates. None is in
  `_workdir/sorting_hat/`. Section 11 excludes them. The two sections contradict.

Fix:
- Either add proteomes of species where Swiss-Prot has a hydrophobin without a Pfam model (G. zeae PH-1, H.
  virens Gv29-8, P. expansum, F. fulva), or state that the rescue gain is measured only at sequence level on the
  Swiss-Prot set (a report table, not a status).
- Remove the contradiction between 4.4 and 11.

### 5. MAJOR: the spec depends on unmerged code, and its branch does not hold that code

Evidence:
- `git merge-base --is-ancestor per-call-status HEAD` fails. `src/cellsurface_sorting_hat/call_status.py` does
  not exist on `hydrophobin-validation`. `docs/paper/03` on this branch has no section 6a.
- PR #76 (`per-call-status`) and PR #75 (`signoff-cfem-hydrophobin`) are open (`gh pr list`).
- F5 says activation "on 2026-10-08 (PR #75)". PR #75 is not merged.

Fix: state in section 10 that H1 starts after PR #75 and PR #76 merge, or that the branch is rebased on
`per-call-status`.

### 6. MAJOR: F3 E-value range is wrong, and F3's regex claim has no stored source

Evidence (`nopfam9.domtbl`, best full-sequence E-value per protein):

| protein | model | E-value | score |
|---|---|---|---|
| HCF2_FULFL | Eas | 2.1e-08 | 23.6 |
| PSH_FLAVE | Hydrophobin | 2.4e-08 | 23.9 |
| HYD1_GIBZE | Eas | 7.6e-07 | 18.6 |
| HFBE_PENEN | Hydrophobin | 7.4e-06 | 15.9 |
| HYD2_GIBZE | Eas | 9.9e-04 | 8.6 |

- The range is 2.1e-8 to 9.9e-4, not "7e-6 to 1e-3". Three of five are below 7e-6.
- PSH_FLAVE has a full-sequence score of 23.9, above the PF01185 sequence GA of 23.00. Its domain score is 22.1,
  below the domain GA of 23.00. It is 0.9 bits from a call. (GA lines from `models.hmm`.)
- Count of 5 with a hit and 4 without is correct.
- "All 9 contain the 8-cysteine pattern" is correct. I re-ran the regex: all 9 match, with 8 to 11 Cys. But the
  cited source (`nopfam9.domtbl`) does not hold this. No file in the repository holds it.

Fix: correct the range. Store the regex result per protein as a table and cite it.

### 7. MAJOR: F6 counts are wrong, and HsbA dominates the new module in the only measurable species

Evidence (proteins hit at `--cut_ga`, discovery domain tables):

| proteome | 7 hydrophobin-class models | with HsbA |
|---|---|---|
| Afum_Af293_UniProt | 7 | 16 |
| Afum_A1163 | 6 | 13 |
| Afum_W72310 | 8 | 14 |
| Cimm_RS | 2 | 2 |
| Bder_ER3 | 2 | 4 |
| Scer_S288C, Calb_SC5314, Cneo_H99 | 0 | 0 |

- F6 says "1 to 6 in the eight proteomes searched". Three have 0. A. fumigatus has 6 to 8. With HsbA, up to 16.
- The eight proteomes are five species (three are A. fumigatus strains). For a per-species status, A. fumigatus
  counts once.
- "More in some Basidiomycota" is not backed by any proteome searched here. It rests on F7 (Xu 2021: 40 genes in
  P. ostreatus).
- In Af293, 9 of 16 `pfam_hydrophobin` hits would be HsbA. HsbA has no 8-Cys hydrophobin label in the truth
  design. `hydrophobin_domain` in A. fumigatus would be mostly non-hydrophobins by the spec's own truth.
  Risk 4 names the issue but does not give this number.

Fix: correct F6. Decide before H1 whether HsbA (and Hydrophobin_like) go into `pfam_hydrophobin`, a separate
module, or a separate call. Report the per-model table as planned, with these counts.

### 8. MAJOR: leakage `none` is not justified

Evidence:
- The author has read the 9 no-Pfam proteins and run a regex on them (F3).
- The author has run the regex on the 8 proteomes that will be measured, and has seen the false-positive list
  (section 6 cites it as the reason for the keep criterion).
- D4 default (no length cap) keeps PSH_FLAVE (515 aa, motif match at residue 445). This is a truth protein.
  The spec cites CFTH1 for the choice, but the choice also fits a protein already seen.
- H4 uses the 9 no-Pfam proteins as unit tests. A rule written to pass these tests is tuned on them. They are
  also truth positives.
- Task order: H2 (build the truth) comes before H3 (fix the spacing). Section 5.1 says the spacing must be fixed
  "before any truth protein is looked at". That is already not true, and the order makes it less true.
- `docs/paper/03` section 5: any leakage except `none` caps at `smoke`. Status is `smoke` anyway (finding 10),
  so the honest label costs nothing.

Fix:
- Record leakage `partial` with a note naming the 9 proteins and the discovery search. Or exclude the 9 proteins
  and their clusters from the measured set, and exclude the 8 discovery proteomes from any `none` claim.
- Do not use the 9 as unit tests. Use synthetic sequences or proteins outside the truth set.
- Move H3 before H2, or freeze the pattern in a commit before H2 starts.

### 9. MAJOR: the dev/test split is conditional and too small to be useful

Evidence:
- Section 5.3 makes the split only "if a change to the pattern is made after seeing truth results". A split
  decided after seeing results is not made "before any result is seen".
- In-scope per-species positives: at most 7 (Af293). A 30% split gives about 2 dev positives. Cluster counts
  are not measured yet (section 4.2 says so).

Fix: make the split unconditional in H2, with the seed, before any measurement. If the test part has fewer
than about 10 positive clusters, say that the split gives no usable test and report the result as `smoke` with
leakage `partial`.

### 10. MAJOR: `estimated` is unreachable in every in-scope species; the goal in section 1 overstates

Evidence:
- `docs/paper/03` section 2: at least 20 positives, 20 negatives, 20 clusters per class, and half-width <= 0.10,
  per species.
- In-scope species with Swiss-Prot hydrophobins: only A. fumigatus, with 7. Cimm and Bder have none in
  Swiss-Prot. Scer, Calb and Cneo have no Pfam hydrophobin hits.
- In the whole Swiss-Prot set, only P. ostreatus has more than 20 named entries (28; next F. velutipes 12,
  T. asperellum 10). Entries with a PubMed citation in the function comment: P. ostreatus 14, B. bassiana 9,
  A. fumigatus 7. This is a rough proxy, not the T1/T2 tier.
- Section 4.4 says "Most will be `smoke`". With these numbers, all will be.
- Section 1 promises "sensitivity X and specificity Y in species Z". Only `smoke` numbers with wide intervals
  are possible.

Fix: say "all status files will be `smoke`" in sections 1, 4.4 and 7.1. State what a `smoke` status adds over
`unvalidated` for the reader.

### 11. MAJOR: bulk assumed negatives conflict with `docs/paper/03`, and differ from the repeat precedent

Evidence:
- `docs/paper/03` section 3: "Specificity is never inferred from absence in a database. A negative needs an
  independent reason."
- Section 4.3 item 2 uses "all other proteins of the measured proteomes" as negatives.
- The repeat work limited the population to curated proteins plus reviewed, secreted UniProt proteins mapped by
  exact sequence (`docs/reports/2026-10-08-repeat-call-calibration.md` section 2). It still called its
  negatives "assumed".
- Hard negatives (section 4.3 item 1) come from Swiss-Prot entries of many species. Most are not in an in-scope
  proteome, so they cannot enter a per-species entry.

Fix: follow the repeat precedent (reviewed secreted proteins of the species), or ask the owner to amend
`docs/paper/03`. Use cross-species hard negatives only in a sequence-level report table.

### 12. MAJOR: the keep criterion and the metric list need work

Evidence and problems:
- "Keep the rescue only if added FP < added TP, per species and pooled." With 0 TP and 0 FP, `0 < 0` is false.
  No tie rule is given.
- The FPs come from assumed negatives. Risk 7.2 says some are probably true hydrophobins. The criterion counts
  them before the manual review.
- No interval. With counts of 0 to 5, the comparison is noise.
- "Per species and pooled" both: in Af293 the criterion fails by design (finding 4).
- Missing metrics: precision of the call and of the rescue-only calls; the number of `not_assessable` records
  (missing SignalP or module rows); FP count per proteome; sensitivity of the pattern on the Pfam-positive
  Swiss-Prot entries (165), which shows how often the fixed spacing misses known members.

Fix: define the criterion on rescue-only calls **after** the manual review, as a precision with a Wilson
interval and a fixed threshold (for example lower bound >= 0.5), and define the result for 0/0. Add the
missing metrics.

### 13. MAJOR: rescue sensitivity on named hydrophobins is partly circular

Evidence:
- F7: the family is defined by the eight cysteines (Kubicek 2008, Seidl-Seiboth 2011 abstracts). Xu 2021 found
  the P. ostreatus genes by this pattern.
- A truth set of named hydrophobins therefore holds proteins that were named because they have the pattern.
  T1/T2 removes Pfam-derived names, but not pattern-derived names.
- T1 lists "a deletion phenotype" as evidence. A generic phenotype (for example slower growth after RNAi of
  FBH1 in Xu 2021) does not show hydrophobin function.

Fix: state that pattern sensitivity measures conformity to the chosen spacing, not discovery. Restrict T1
phenotypes to hydrophobin-specific ones (rodlets, wettability, surface activity).

### 14. MAJOR: H1 makes the existing repeat call files stale; this task is missing

Evidence:
- `call_status.stale_reason` returns `"config differs"` when `config_sha256` differs
  (`call_status.py` on `per-call-status`, lines 177-180). `config_sha256` is the hash of the whole
  `categories.yaml` (`engine.load_config`).
- H1 edits `categories.yaml`. The call files `_workdir/sorting_hat/Scer_S288C/status/calls/tandem_repeat_protein.json`
  and `.../Calb_SC5314/status/calls/tandem_repeat_protein.json` then go stale. The run falls back to module
  statuses.
- `artefact_digest` of every Pfam module is `pfam_sha256 + ":" + _family_digest(families)` over the whole table
  (`modules/cli.py` lines 185-191, 248). Any edit of a hydrophobin row changes the identity of `pfam_adhesion`
  too, and the reverse.

Fix: add a task after H1 to re-run `calibrate truth --call-status` for `tandem_repeat_protein` in both species.
Consider a per-module family digest so that the two Pfam modules do not invalidate each other.

### 15. MAJOR: H1 misses several code and test changes

Evidence:
- `tests/cellsurface_sorting_hat/test_data_and_scripts.py` line 43 asserts the module set is exactly
  `{"pfam_adhesion", "pfam_allergen"}`.
- `outputs.py` lines 29-32: report text lists hydrophobin among the families of `wall_family_domain`.
- `test_cli.py::test_golden_calls` asserts `other_basis == "cocci_specificity_rank_top15"` for ENZ1 and STAR1.
  If `hydrophobin_domain` becomes a mechanism and the fixture has no `pfam_hydrophobin` table, the call is
  `not_assessable` and `other_basis` changes. On `per-call-status` there is also `test_run_output_is_pinned`
  with files in `tests/cellsurface_sorting_hat/golden/`.
- D2 covers only `other_surface_no_mechanism`. Section 3.1 changes `other_not_surface` too. The spec does not
  decide that list.
- Existing module tables have only PA14 active (`modules/pfam_adhesion.json`, `params.families =
  ["PF07691"]`) and the run domain tables predate the new families. H5 must re-run hmmsearch, not only the
  wrapper.

Fix: list these items in H1 and H5. Add `other_not_surface` to D2.

### 16. MAJOR: F2 provenance is incomplete

Evidence:
- F2 counts are correct. I re-derived: 189 rows; 174 match "hydrophobin|rodlet" in the protein name; 165 with
  one of the seven models (PF01185 99, PF06766 47, PF22354 12, PF29802 2, PF28987 2, PF29785 2, PF29465 1; no
  entry has two); 9 with none; 15 rows do not match the name. No named entry has HsbA.
- The `Pfam` column is UniProt's cross-reference. It is not an hmmsearch at Pfam 38.2 GA. The two can differ.
- Section 4.2 says the query is "in `analysis/hydrophobin_truth/`". The directory holds only
  `sp_hydrophobin_query.tsv`, `nopfam9.faa` and `nopfam9.domtbl`. No query string, no UniProt release, no date.
- `analysis/hydrophobin_truth/` is untracked (`git status`: `??`).

Fix: store the query URL, the UniProt release and the fetch date. Say that F2 uses the UniProt cross-reference.
Commit the directory.

### 17. MINOR: section 3.3 contradicts itself

"Both variants have at least two modules" is followed by "`pfam_only` reads one module". Delete the first
sentence.

### 18. MINOR: rule 5.2 item 3 is redundant

The pattern C-CC-C-C-CC-C has eight cysteines. Any match has at least 8 in the span. If the intent is a count in
the mature region or an upper bound, say so.

### 19. MINOR: D4 is not operational

"Count domain-like spans" does not change a per-protein `hit`. State what field records the count and whether
it affects the call.

### 20. MINOR: `cys8_pattern` identity

The pattern, the length window and the spacing source must be in the module `params` so that `params_hash`
changes when they change. Otherwise a call file stays valid after a rule change. The spec does not say this.

### 21. MINOR: family table note for PF01185 disagrees with the discovery tables

`data/sorting_hat/family_table.tsv` line 4: "18 hits in 4 proteomes (A. fumigatus Af293, A1163, W72310 and
C. immitis RS)". The discovery tables give 6 + 5 + 6 + 1 = 18 in those four, plus 2 in Bder_ER3. Twenty hits
in five proteomes.

### 22. MINOR: F7 and the novelty section

- All seven PMIDs resolve to the papers named, and the abstract statements match (PubMed metadata, fetched
  2026-10-08). Yang 2006: doi 10.1186/1471-2105-7-S4-S16. Kubicek 2008: 10.1186/1471-2148-8-4. Seidl-Seiboth
  2011: 10.1007/s00239-011-9438-3. Xu 2021: 10.1016/j.micres.2021.126723. Pitocchi 2026:
  10.1016/j.ijbiomac.2026.153896. De Vries 1999: 10.1046/j.1432-1327.1999.00387.x. Peñas 1998:
  10.1128/AEM.64.10.4028-4034.1998.
- Yang 2006 (PMC1780129), Kubicek 2008 (PMC2253510) and Peñas 1998 (PMC106595) have open full text. "Abstracts
  only" can be lifted now for these three.
- The Yang 2006 abstract says candidates were kept after "filtering by pattern, domain and length", and all 9
  "possess the common pattern and hydrophobin domain". That is a pattern-plus-domain combination. Novelty item 1
  ("sensitivity gain over Pfam alone") must be compared with this before any claim.
- Item 1 says "F2 and F3 show 5% have no Pfam model". 9/174 = 5.2% is correct, but it is a count on
  Swiss-Prot names and UniProt cross-references. It is an upper bound on the gain in that set, not a measured
  gain. Section 8 should say so.
- Item 3 is a property of the tool, not of the rule.
- The section does not claim the pattern. The limits (PubMed only, abstracts only) are stated. Apart from the
  points above, it is honest.

### 23. MINOR: I did not test the class I and II spacing from memory

Section 5.1 gives two spacings "from memory". I did not run them on any truth protein. Running them would add
to the leakage of finding 8. H3 should fix them from the papers first.

## Verdict

**Not ready for a plan.** Four blockers:
1. The call syntax (`variants`, `ref`) cannot produce a call status file with the existing engine.
2. The module cannot cut the signal peptide from the step 1 table.
3. The negative-set rule makes Pfam-only specificity 1.0 by construction.
4. In the proteomes in scope, the rescue cannot gain a true positive, so the keep criterion decides by scope.

Also correct F3 and F6, state leakage `partial`, and state that every status will be `smoke`. Then decide
on HsbA before H1. The Pfam-only part (H1: move rows to `pfam_hydrophobin`, add `hydrophobin_domain`, fix the
F5 side effect) is sound in idea and can go ahead separately once findings 5, 7, 14 and 15 are handled.
