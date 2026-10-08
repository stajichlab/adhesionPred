# Review 2: hydrophobin call validation and 8-cysteine rescue rule (revision 2)

*2026-10-08. Independent review by Claude Code (claude-opus-5-5). Reviewed:
`docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md` (revision 2, commit 78b0739, branch
`hydrophobin-validation`) against review 1 (`...-design-review-1.md`). The spec was not edited. Every
statement below was checked against a file, a command output or code. Where I did not check something, I say
so.*

## How I checked

- Code on branch `per-call-status` (PR #76, open), read with `git show per-call-status:<path>`:
  `engine.py` (`_validate`, `_eval`, `_eval_other`, `call_eligible`, `call_hash`), `call_status.py`
  (`stale_reason`), `modules/cli.py` (`_condition_table`, `_family_digest`, the `pfam` branch of `run`),
  `modules/pfam.py`, `modules/signalp.py`, `modules/base.py` (`write_module`), `calibration/cli.py`
  (`truth`, `_check_call_run`), `categories.yaml`, `docs/reports/2026-10-08-repeat-call-calibration.md`.
- Git: `signoff-cfem-hydrophobin` (PR #75) is an ancestor of `hydrophobin-validation`. `per-call-status` is
  not. PR #75 is not in `per-call-status`. `gh pr list`: #75 and #76 are open.
- `data/sorting_hat/family_table.tsv` on this branch.
- `analysis/hydrophobin_truth/` (tracked since 78b0739): re-counted species and strains of the 174 named
  entries with Python.
- Re-ran `analysis/sorting_hat_run/hydrophobin_discovery/analyse.py` (read-only) on the discovery tables.
- R0 run records (`_workdir/sorting_hat/*/modules/step1_rule@R0.json`): all eight are `ok`.
- Wilson lower bounds computed with z = 1.96.
- PubMed ID converter for PMIDs 21424760, 17217508, 18186925.

## Part 1: status of the 23 findings of review 1

| # | Review 1 finding | Status | Evidence |
|---|---|---|---|
| 1 | `variants` syntax; `pfam_or_cys8` not eligible | RESOLVED | 3.4 uses `or` of two flags, ungated, no `ref`. On `per-call-status`, `call_eligible` returns True for this shape: no `ref` (`_has_ref`), no literal step 1 module, `reads_of_call` gives 2 modules. It passes `_validate_node`. It can be a mechanism if it is placed before the `other_*` calls (new finding 6). |
| 2 | Cleavage site not in the SignalP table | RESOLVED | 3.3: the pattern runs on the full sequence; the SignalP module is not changed. |
| 3 | Pfam-only specificity 1.0 by construction | RESOLVED | 4.3 item 3: one negative set; a Pfam hit with no label is a negative for both variants. |
| 4 | Rescue gain cannot be measured in scope | PARTIAL | 4.5 states it and adds E1 (D7). But the E1 list misses one of the nine species, strain mapping is not checked, and rule 6.3 pools in-scope proteomes where no rescue-only call can be labelled true (new findings 1 and 4). |
| 5 | Depends on unmerged code | RESOLVED | 3.1 states both dependencies. Small wording mismatch with section 10 (new finding 9). |
| 6 | F3 range wrong, no stored regex result | RESOLVED | F3 values equal the review 1 table. `nopfam9_regex.tsv` exists and is tracked: 9 rows, all `regex_match` 1, 8 to 11 Cys. |
| 7 | F6 counts wrong; HsbA dominates | RESOLVED | F6 now matches. My re-run of `analyse.py`: PF01185 plus other-model-only hits are Af293 6+1, A1163 5+1, W72310 6+2, RS 1+1, ER3 2+0. D6 separates HsbA. |
| 8 | Leakage `none` not justified | RESOLVED | 5.3: `partial`, the 9 are not unit tests, H3 before H2, the pattern is frozen in a commit. See new finding 11 for the "seen" list. |
| 9 | Conditional, small split | RESOLVED | 4.4 item 2: unconditional, by cluster, seed, made in H2, with a rule for fewer than 10 test clusters. |
| 10 | `estimated` unreachable | RESOLVED | Sections 1, 4.5 and 7.1 say every status is `smoke`. 4.5 omits the half-width criterion (new finding 12). |
| 11 | Bulk assumed negatives | RESOLVED | 4.3 item 1 follows the repeat precedent, labels them "assumed", and adds D8. Hard negatives are report-only. |
| 12 | Keep criterion and metrics | PARTIAL | 6.3 has a Wilson rule and a 0-case. 6.1 adds the missing metrics. But 6.3 conflicts with 6.2 and pools strains (new finding 1). |
| 13 | Pattern sensitivity partly circular | RESOLVED | 4.1 "What the labels do not remove"; T1 restricted to hydrophobin-specific phenotypes. |
| 14 | H1 makes repeat call files stale | RESOLVED | H1b re-runs both files; H1c is the per-module digest. 0.435 and 0.462 match the repeat report section 6 table. `tandem_repeat_protein` reads only `repeat02` and `repeat14`, so its numbers should not change. But see new finding 3 on timing. |
| 15 | H1 misses code and test changes | RESOLVED | H1 lists `MODULES`, the module-set assertion, `outputs.py`, golden files, `other_basis` strings and both D2 lists. H5 re-runs `hmmsearch`. |
| 16 | F2 provenance | RESOLVED | `README.md` holds the query, the fields, the URL, release 2026_03 and the date, and says the Pfam column is the UniProt cross-reference. The directory is committed. Small defects in new finding 10. |
| 17 | 3.3 contradicts itself | RESOLVED | The sentence is gone. |
| 18 | Rule 5.2 item 3 redundant | RESOLVED | 5.2 item 3 says no separate count is added. |
| 19 | D4 not operational | RESOLVED | `n_spans` is a separate field; `hit` does not change. |
| 20 | `cys8_pattern` identity | PARTIAL | 3.3 puts the pattern, the spacing source and the length window in `params`. It does not put the R0 condition identity there (new finding 2). |
| 21 | PF01185 note disagrees | PARTIAL | H0 schedules the fix. The note in `family_table.tsv` line 4 still says "18 hits in 4 proteomes". |
| 22 | F7 and novelty | RESOLVED | Section 8 item 1 says "upper bound"; H3b reads the open full texts. One new error on open access (new finding 8). |
| 23 | Spacing from memory untested | RESOLVED | 5.1: unverified, in no file, not run. |

Count: 18 RESOLVED, 5 PARTIAL, 0 NOT RESOLVED.

## Part 2: new findings

### 1. MAJOR: rule 6.3 conflicts with rule 6.2, and with in-scope numbers it decides before E1 is measured

Evidence:
- 6.2: "A label changes only with a recorded reason and a PMID."
- 6.3: the rescue is kept only if rescue-only calls are >= 5 and the Wilson lower bound of precision
  "after manual review" is >= 0.5. 6.3 says this "does not count unlabelled true hydrophobins as false".
- These two rules disagree. A rescue-only call is, by definition, a protein that Pfam misses. Most such
  proteins in a proteome have no paper. Under 6.2 a reviewer cannot make one true without a PMID. So under
  6.3 it counts as false. The spec does not say what other review outcome exists.
- Wilson lower bound (z = 1.96) at precision 1.0: n = 5 gives 0.566, n = 6 0.610, n = 7 0.646. For
  n = 5 to 7, every rescue-only call must be true. n = 10 needs 9 true; n = 15 needs 12.
- In-scope proteomes, screening regex of `analyse.py` (secreted by R0, 60 to 400 aa, no hydrophobin-class
  hit): Af293 3, A1163 5, W72310 3, Cimm RS 3, Bder ER3 5, Cneo H99 1, S288C 0, C. albicans 0. Total 20. None
  has a hydrophobin label (F7). The frozen pattern will differ, and the spec has no length cap, so the real
  count can be higher or lower. I did not run any frozen pattern.
- The three *A. fumigatus* strains appear to give the same proteins three times. Example: Q4WBD4 (Af293),
  B0YAH7 (A1163) and KAK9640788.1 (W72310) are each 278 aa with 16 Cys; Q4WUX0 and B0Y4E8 are each 301 aa
  with 18 Cys. I did not align them. Pooling "over all measured proteomes" can count one protein up to three
  times. Wilson assumes independent trials.
- E1 can add at most 8 labelled rescue-only true positives (finding 4 below). With about 20 in-scope calls
  that cannot become true, 8/28 has a Wilson lower bound far below 0.5. So the rule would reject the rescue
  because of the in-scope proteomes, not because of the E1 measurement.

Fix:
- Define the review outcomes in 6.2, for example: `true (PMID)`, `hydrophobin by review` (with the evidence
  that counts: BLAST top hit to a named hydrophobin, spacing class, length), `not hydrophobin`, `unresolved`.
  State that `hydrophobin by review` is a prediction-based label and is reported apart from the T1/T2 result.
- State in 6.3 which outcomes count as true and how `unresolved` is counted (as false, or reported as a
  range with both choices).
- Count orthologs once: for *A. fumigatus*, count rescue-only calls per orthogroup or per species, not per
  strain. Or apply 6.3 per species and require it in at least N species.
- State the expected result of 6.3 if D7 is "no": rejection, decided by scope.

### 2. MAJOR: the R0 dependency of `cys8_pattern` is not in its identity as specified

Evidence:
- 3.4: "Its status depends on the identity of the R0 module through `cys8_pattern`."
- A call status file records `reads` = the identities of the modules the call reads (`call_status.py`
  `build_source`; `stale_reason` compares them). For `hydrophobin_protein`, `reads_of_call` gives
  `pfam_hydrophobin` and `cys8_pattern` only. R0 is not a read module.
- A module identity is (`name`, `version`, `params_hash`, `artefact_hash`) (`base.write_module`). For
  `pfam`, the R0 identity enters only because `run` puts `conditions = {"sp_module": {module, params_hash,
  artefact_hash}}` into `params` (`modules/cli.py`, `pfam` branch and `_condition_table`).
- 3.3 lists the `cys8_pattern` params as the pattern string, the spacing source and the length window. It
  does not list the condition identity. Built as written, a new R0 run would not make a `hydrophobin_protein`
  call file stale.
- The pfam precedent copies `params_hash` and `artefact_hash` of R0, not its `version`. R0's
  `artefact_hash` is `_tool_digest("signalp", version)`. It does not hash the SignalP results.

Fix: state in 3.3 that `cys8_pattern` puts `conditions = {"sp_module": <_condition_table identity>}` in its
params, as `pfam` does. Consider adding the R0 `version` to that identity (and to the pfam one). Add a test
that a changed R0 `params_hash` changes the `cys8_pattern` `params_hash`.

### 3. MAJOR: the order of `categories.yaml` edits makes call files stale again, and the end state after a "drop" is not defined

Evidence:
- `stale_reason` returns "config differs" when `config_sha256` differs. `config_sha256` is the hash of the
  whole `categories.yaml` (`engine.load_config`).
- H1 adds `hydrophobin_protein` and the D2 mechanism changes. The `cys8_pattern` wrapper is H4. Between H1
  and H4 no `cys8_pattern` table exists. Then `hydrophobin_protein` is `not_assessable` for every protein
  without a Pfam hit (`k_or` of `not_called` and `not_assessable`), and `_eval_other` lists it in
  `other_basis` as left out for every surface protein.
- 6.3 may remove the rescue from the default call after H6. That edits `categories.yaml` again. Then the
  `tandem_repeat_protein` files from H1b and the `hydrophobin_protein` file from H6 are stale ("config
  differs").
- If the rescue is dropped, `hydrophobin_protein` becomes `{flag: pfam_hydrophobin.hit}`. That is the same
  as `hydrophobin_domain`, reads one module, and fails `call_eligible` ("reads fewer than two modules"). The
  H6 call status file then names an ineligible call, and `build_source` raises; the run stops
  (`CallStatusResolver` turns it into `CallStatusError`).
- "`cys8_pattern` registration" in H1 has no meaning in the code: modules are not registered in the engine.
  The module is a `cellsurface_sorting_hat_module` subcommand, written in H4.

Fix:
- Fix the end state now. For example: keep `hydrophobin_domain` as the D2 mechanism; add
  `hydrophobin_protein` (and decide its mechanism role) only after H6, and only if 6.3 keeps the rescue. If
  6.3 drops it, `hydrophobin_protein` is not added, and the rescue result is a report table.
- Run H1b after the last `categories.yaml` edit of this work, not after H1. Or re-run it at the end too.
- Remove "`cys8_pattern` registration" from H1, or say what it means.

### 4. MAJOR: E1 does not cover all nine no-model entries, and the strain mapping is not checked

Evidence (`sp_hydrophobin_query.tsv`, the 174 named entries):
- The nine no-model entries are from six species: *P. expansum* 1, *T. asperellum* 1, *H. virens* 1,
  *G. zeae* 2, *F. fulva* 3, *F. velutipes* 1 (PSH_FLAVE).
- 4.5 says "the 9 no-model entries come from five of these". *F. velutipes* is not in E1. Eight of nine come
  from five E1 species. PSH_FLAVE, the protein 0.9 bits below GA (F3, risk 7), cannot enter a proteome
  measurement.
- 4.5 says "0 to 28 named entries each". Counts of named entries in the E1 species: *G. zeae* 5,
  *H. virens* 2, *P. expansum* 7, *F. fulva* 6, *P. ostreatus* 28, *B. bassiana* 9, *T. asperellum* 10. The
  range is 2 to 28.
- 4.3 maps UniProt entries to a proteome by exact sequence. Strain in the UniProt organism field: HYD1_GIBZE
  and HYD2_GIBZE are PH-1; HFB3_HYPVG is Gv29-8. HFBE_PENEN, HYD1_TRIAP, HCF1/2/4_FULFL and PSH_FLAVE have no
  strain. Of the *T. asperellum* named entries, 2 have no strain and 8 are ATCC 204424. I did not check
  whether any of the 8 sequences is in an E1 reference proteome.
- So the labelled rescue gain in E1 is between 0 and 8 proteins, and only 3 are strain-matched by name.

Fix: correct the 4.5 sentence and the range. Add *F. velutipes* to E1 or say PSH_FLAVE is sequence-level
only. Before D7 asks for the download, check by exact sequence (or by a stated identity cut-off, which then
must be named in the label rules) which of the 8 are in the chosen E1 proteomes. Name the strain of each E1
proteome.

### 5. MINOR: `cys8_pattern` row states are not specified, and the pfam precedent would block the call status

Evidence:
- 3.3: `hit` is "Empty (not assessable) if R0 is missing". It does not say the row state.
- The pfam precedent gives a protein with no usable condition value the state `error` (`pfam.pfam_rows`).
  `_write` then sets the run state to `partial`. `calibrate truth --call-status` refuses a read module that is
  not `ok` (`_check_call_run`), and `stale_reason` returns "module state not ok".
- `_condition_table` raises if the R0 table or run record is missing, or its run state is `unavailable`,
  `not_run` or `error`. 3.3 does not say this.
- All eight R0 run records in `_workdir/sorting_hat/` are `ok` now.

Fix: state in 3.3: a protein with no `ok` R0 row gets state `ok` and an empty `hit` (the engine reads an
empty flag as `not_assessable`); an invalid sequence gets `invalid_row`; a missing or unusable R0 table stops
the wrapper. Add a test that a missing R0 row leaves the run state `ok`.

### 6. MINOR: mechanism placement and the embedded R0 gate

Evidence:
- `_validate`: a mechanism "must be an earlier ungated call". `hydrophobin_protein` must stand above
  `other_not_surface` in `categories.yaml`.
- `hydrophobin_protein` is ungated but `cys8_pattern.hit` needs R0 `called`. In the `other_*` calls of
  variants R1, R2 and `card`, a protein secreted by that variant but not by R0 has `cys8_pattern.hit` 0. It
  can be `other_surface_no_mechanism` while it matches the pattern. `wall_family_domain` with PA14
  (`second_condition = signal_peptide`) already behaves this way, so there is a precedent.

Fix: state the placement in H1 and the R0-only gate in 3.4.

### 7. MINOR: `hsba_domain` has no mechanism role and the 2c display step is not specified

Evidence: on this branch HsbA (PF12296) is active in `pfam_adhesion`, so a secreted HsbA protein is
`wall_family_domain` now. D6 moves it to `pfam_hsba` with call `hsba_domain`. D2 adds only
`hydrophobin_protein` to the mechanism lists. A secreted HsbA protein then becomes
`other_surface_no_mechanism`. 3.2 says the report combines the calls "as a display step"; no file or function
is named. `pfam.MODULES` and `SECOND_CONDITIONS` checks also need `pfam_hsba`.

Fix: decide in D2 whether `hsba_domain` is a mechanism. Name where the 2c display combination is made.

### 8. MINOR: open-access claim in 5.1 is wrong for Seidl-Seiboth 2011

Evidence: 5.1 says "full text open for the last two [Seidl-Seiboth 2011, Yang 2006] and for Kubicek". The
PubMed ID converter gives PMC1780129 for Yang 2006 (17217508) and PMC2253510 for Kubicek 2008 (18186925), and
no PMCID for Seidl-Seiboth 2011 (21424760). F8 lists Yang, Kubicek and Peñas as open. Wessels 1994, Linder
2005 and Sunde 2008 have no PMID in the spec; I did not check their access.

Fix: correct 5.1 to match F8. Give PMIDs for the H3 papers and say how H3 gets any closed full text. If a
spacing cannot be read from a source, H3 must record that.

### 9. MINOR: H1 dependency wording differs between 3.1 and section 10

Evidence: 3.1: "after both merge, or after the work branch is rebased on `per-call-status`". Section 10 H1:
"After PR #75 and #76 merge". This branch already contains PR #75 (`signoff-cfem-hydrophobin` is an
ancestor), so a rebase on `per-call-status` gives both. The rest of the order is consistent: H3 and H2 need
no code; H4 needs only the frozen pattern from H3; H5 needs H1 and H4. H3b repeats reading done in H3 (Yang,
Kubicek).

Fix: use the 3.1 wording in H1. Say that H3, H3b and H2 can start before #75 and #76 merge.

### 10. MINOR: H0 is partly done; README has stray lines

Evidence: `analysis/hydrophobin_truth/` is tracked since 78b0739, so "Commit `analysis/hydrophobin_truth/`"
in H0 is done. `README.md` lines 5 to 7 are pasted HTTP response headers (`access-control-expose-headers`,
`x-uniprot-release`, `x-uniprot-release-date`).

Fix: remove the commit item from H0. Replace the header lines with one sentence on the release and date.

### 11. MINOR: the "seen" list must include proteins named in reviews

Evidence: 5.3 names the 9 entries and the discovery search. Review 1 and this review list in-scope
regex-positive, Pfam-negative proteins by ID (for example Bder ER3 `F00FD2C2_006465-T1`, 101 aa, 8 Cys, and
`F00FD2C2_006466-T1`, 127 aa, 9 Cys, both secreted by R0 with no hydrophobin-class hit). I did not check what
these proteins are. They will be candidate rescue-only calls in H6.

Fix: in 5.3 and in the "seen" column, include every protein named in the discovery output and in the review
files.

### 12. MINOR: section 4.5 states the `estimated` rule incompletely

Evidence: 4.5 says "at least 20 positives, 20 negatives and 20 clusters per species". Review 1 cites
`docs/paper/03` section 2 as 20 per class and a half-width <= 0.10. The conclusion (`estimated` not expected)
does not change.

Fix: quote the rule in full.

### 13. MINOR: the keep decision uses the test part

Evidence: 4.4 says the pattern and any threshold are tuned on the development part only. 6.3 decides whether
the rescue is in the default call, using all measured proteomes. Choosing the call is a tuning step. Leakage
is `partial` and status is `smoke` anyway, so the label does not change.

Fix: say in 6.3 that the decision uses all data and that this is part of the `partial` note. Or apply it to
the development part only.

## Verdict

**Not ready for a plan yet. No blockers.** Revision 2 fixes all four blockers of review 1, and the call shape
works with the engine on `per-call-status`. Four MAJOR findings need a spec edit first:

1. Make rules 6.2 and 6.3 agree, and count *A. fumigatus* strains once (finding 1).
2. Put the R0 condition identity in the `cys8_pattern` params (finding 2).
3. Fix the end state of the calls and the order of `categories.yaml` edits, so that call files are not made
   stale or ineligible by a later edit (finding 3).
4. Correct E1 coverage and check strain mapping before D7 (finding 4).

These are text changes. The Pfam-only part (`pfam_hydrophobin`, `hydrophobin_domain`, its module status, H1c)
does not depend on them and is ready for a plan once findings 3 and 7 are decided.
