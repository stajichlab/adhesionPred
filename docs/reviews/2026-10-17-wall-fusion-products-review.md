> What this is: an independent review of the cell wall / cell fusion products (reports 01-07, errata f3b1813, wall_fusion_report 61ecc57). Verbatim copy; the text below the rule is unchanged.
> Date: 2026-10-17. Source: scratchpad copy of wall_fusion_products_review.md. Findings are cited in the spec as WF-F1..WF-F11.
> Independent review by an Opus agent; not verified by the owner.

---

# Review: cell wall / cell fusion products (reports 01-07, errata f3b1813, wall_fusion_report 61ecc57)

Date 2026-10-17. Read-only review. Line refs: `main:` = /Users/jstajich/projects/adhesionPred at c7eadf0;
`err:` = the errata worktree (/Users/jstajich/projects/adhesionPred-errata, f3b1813; line numbers after the errata);
`ws:` = /Users/jstajich/projects/adhesionPred-wallscript (61ecc57).

## What I did (verified) vs inferred

Verified by reading or running:
- Read 05, 06, 07 in full; main 04 in full; the full `git diff main..wall-fusion-errata`; the conclusion sections of 01, 02, 03.
- Spot-checked `analysis/cell_wall_check/summary.tsv`, `full_pfam_ga_hits.tsv`, `fullpfam.domtbl`, `scan.versions.tsv`,
  `docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` (awk counts), and `data/sorting_hat/family_table.tsv` (active column).
- `git log` for `family_table.tsv` (commits cbfdc45 and d94a8ab).
- Read `wall_fusion_report.py`, `wall_fusion_models.tsv`, README, sbatch, `results/*`, the tests; compared to main
  `analysis/lectin_class_report/lectin_report.py`.
- Ran `uv run --no-project --with pytest pytest tests/wall_fusion_report -q` in ws: **15 passed**. No `.venv` or `uv.lock` appeared.
  `__pycache__` dirs under `ws:tests/wall_fusion_report/` and `ws:analysis/wall_fusion_report/` exist (probably from my run, git-ignored);
  `.pytest_cache` and `.ruff_cache` have mtime 10:28, a few minutes before my run (probably the author's). `git status` is clean.
- Checked all 54 errata line numbers in the four logs against the diff hunks.

Inferred (not verified): the biology argument on Gel4 (finding 4), and how HMMER treats domain thresholds under `--cut_ga` (finding 10).

---

## Findings ranked by severity

### F1. DEFECT: the errata misses wrong "main" statuses in the same table it fixed (04)
- `err:04:34` PIR1-like_C PF22799 "inactive / active", `err:04:35` Hyr1 PF15789 "inactive / ...", `err:04:142` "PF00399, PF22799 ... inactive on main".
  On main, `family_table.tsv` rows 23 (PF15789) and 24 (PF22799) are `active=yes`. Commit cbfdc45 ("activate PF00624, PF13928,
  PF15789, PF22799, PF30910 (owner 2026-10-09)") turned them on.
- The errata fixed the next two rows (`err:04:36-37`) but left these. A reader now sees corrected and uncorrected "main" states side by side.
- Also wrong in the errata: `err:04:36` says Flocculin/Flocculin_t3 are active since "commit d94a8ab". d94a8ab turned on only
  Candida_ALS_N, Candida_ALS, Flo11 and GLEYA. **Flocculin PF00624 and Flocculin_t3 PF13928 were turned on in cbfdc45.** The commit
  attribution is not in 07 (07 row 6 names d94a8ab for the Als rows only), so the errata added an unsourced fact, and it is partly wrong.
- Also stale on main and not marked: `err:04:10-12` ("None of these exist on main", "`docs/CLASSES.md` is absent"). `docs/CLASSES.md` and
  `analysis/lectin_class_report/` are on main now (merged with #84, before #85). The same applies to `err:04:171` ("found on branches, not on main").
  05-07 do not state this, so it is outside the errata's stated scope, but the header "Errata applied" suggests the file is current.

### F2. DEFECT: claims contradicted by 07's measurements are left unmarked in 01, 02, 03
07 measured these facts (07:39-51, 07 s4): the shipped ESM model does **not** call GAS1 (0.030), DFG5 (0.020), DCW1 (0.020), ECM33 (0.033),
BGL2 (0.017) or EXG1 (0.011). It calls CRH1, UTR2, SUN4, SCW11, KRE9 and CTS1. In the sorting hat every wall enzyme has the derived label
`other_surface_no_mechanism`, not `cell_wall_adhesion_candidate`.
- `err:03:294`: "Likely to score as adhesion-like: Gas/Gel/Phr family members, Dfg5, Ecm33, ...". **The three named S288C members are
  exactly the ones the ESM model does not call.**
- `err:03:296`: "Gas/Gel and BGT would count as false positives for `cell_wall_adhesion_candidate`". 07 s4 bullet 3 and rows 11-12 say
  they would be `other_surface_no_mechanism`, because no table row hits them and they have no repeat call.
- `err:01:77`: "Many carry signal peptide plus GPI, which a surface-glycoprotein classifier will score high." Of the 6 S288C wall enzymes
  with a UniProt GPI site (GAS1, CRH1, UTR2, ECM33, DFG5, DCW1), 2 are called. Of the 7 without one, 4 are called (07:39-51). n is small,
  but the stated mechanism does not fit the data.
- `err:02:146`: "a Pfam or surface classifier would catch the adhesion step (agglutinins)". 07 s4 bullet 1: on the Pfam route only SAG1 is
  caught; AGA1 and AGA2 hit no row. The ESM model misses AGA2.
- `err:03:297`: "Gel/Gas is glycoside hydrolase family 72 ... I did not verify". 05 Table A (PF03198 description: GPI-anchored
  glucanosyltransferases) and 07 row 17 (Gas1 504.4) now support it. It was left as "unverified", while `err:03:164` marks the same
  GH72 placement as an "unresolved conflict". So report 03 now contradicts itself.
- `err:02:178` (open Q2, "what are their Pfam domains?"): 07 row 23 answers the Pfam half (SCW4: PF00332 45.6 and PF29664 283.2;
  SCW11: PF29664 only). The errata log (`err:02:218`) lists only the GPI half as still open.
- 04, unmarked: `err:04:16` and `err:04:128` (the chitin synthase name discrepancy; 05 settles it: "repo is right", 05:18);
  `err:04:136` (Scw4/Scw10 "none | not in repo"; 07 row 23 "Superseded"); `err:04:111` ("No measurement shows how a proposed enzyme or
  fusion family would behave"; 07 is that measurement, though only on 23 proteins).

### F3. CONCERN: the conclusion "wall enzymes are the main false-positive risk" has to be split by tool
- For the ESM classifier the corrected numbers support it: 6/13 wall enzymes called (07 s4), plus the earlier model review (main `04:88-92`).
- For the sorting hat's family route, 07 refutes it: 0/13 hit an active row (07 row 10, s4 bullet 2), and R0 is an all-secreted filter.
  01 s3b, 03 s3.5 and the script's FALSE_POSITIVE_NOTE (`ws:wall_fusion_report.py:54-60`) say "surface-glycoprotein detection" in general.
  The script note is scoped well enough, because it names the ESM classifier. 03 s3.5 is not (see F2).
- "Fusion proteins mostly missed": supported, but only for 4 S. cerevisiae proteins (07:59-62; shipped p 0.002 to 0.048). Report 02's
  wider list (LFD-1, HAM-6, Nc FIG-1 and others) is untested. The statement is inference and labelled as such (`err:02:146`). SOUND as worded.
- "Multicellularity not driven by adhesin expansion": the claim is attributed to Kiss 2019 ([E] at `err:03:269`) and 06 s2 / s7 confirm the
  text. The "[I] little support" at `err:03:272` is a fair reading. One caveat should travel with it: the 23-family set is 42% yeast plus
  Candida genes (06 row 17), which limits detection in early lineages (03 caveats). SOUND.

### F4. CONCERN: the "Gel4 = CAZy family 17" conflict was escalated, not resolved, and a correct statement in 01 now carries a doubt marker
- `err:01:71` inserts "[unresolved conflict, see 05] (Gel4 is called CAZy family 17 in 03:162; not resolved)" next to "GEL/GH-72".
  The same claim at `err:01:62`, `err:01:77` and `err:01:82` is unmarked, so the marking is inconsistent. "03:162" is a pre-errata line
  number; it is now 03:164.
- The conflict can largely be resolved from data already in the repo. In Data 3, gel4 is in cluster 435 with gas1, phr1 and gel7
  (06 s5; TSV row `gel4 435 ... "GPI protein, Gas/Phr family, CAZy family 17"`). The annotation itself says "Gas/Phr family", and
  S. cerevisiae Gas1 hits Glyco_hydro_72 at 504.4 (07 row 17). Inference: "CAZy family 17" is most likely an annotation slip in the
  supplement. The finding should read "file annotation inconsistent with its own cluster", not "GH72 placement unresolved".
  The errata faithfully copies 05 Summary item 4. The weakness is in 05, and the errata spreads it into 01.

### F5. CONCERN: the errata log is inaccurate for some items, and some markers overstate
- The log says "New: same text with an added marker" for every item. `err:04:123` and `err:04:124` (log items 8 and 9) also **replace**
  the Note cell "unverified: check InterPro/Pfam before use" with new text ("verified ID and name (05); ... hits measured (07)").
  Old text was deleted, which the log does not say.
- `err:04:153` says "Prm1 = PF27758 ... and Fig1 = PF12351 ... verified by hmmsearch". 07 found these with `hmmscan` against all of
  Pfam-A (07:13; `full_pfam_ga_hits.tsv`), not with hmmsearch. "Verified" overstates: 07 s4 says "measured, not yet used"; the names come
  from the HMM NAME line and were never checked in InterPro; n = 1 protein each.
- `err:02:20`: Kiss "defines 23 families" is marked "[corrected ...]", but 06 s7 grades that phrase **CONFIRMED** ("slightly generous").
  A clarification is labelled as a correction. `err:02:159`: "23 families = 45 curated genes" equates two different levels, which is the
  confusion 06 s2 warns about. It should read "23 families containing the 45 curated genes".
- `err:03:162` (cluster 209): the original "Opisthokonta + amoebozoa" was incomplete, not wrong. E. histolytica is an amoebozoan; only
  T. trahens is not. The new text is accurate, but the old text did not need a "corrected" tag. Minor.
- Line numbers: all 54 cited "Line N" values match the pre-errata hunks (checked). Sources cited per item match 05, 06 and 07 sections.

### F6. CONCERN: the script lost the duplicate-ID warning that the lectin script has
- `ws:wall_fusion_report.py:523-537` keeps the best row per (protein, model) and drops the others **silently**. `count_records`
  (`:256-266`) counts records, not unique IDs, so `n_proteins` and `n_hit_proteins` can disagree with no notice.
- main `lectin_report.py:323-339` (`scan_fasta`) returns the number of repeated IDs. The lectin `main` warns per proteome ("records repeat
  an ID inside the file ... counted once per ID") and across proteomes (`:741-748`, `:771-778`). These are the review-era fixes listed in
  04 s4.2 (`test_duplicate_id_inside_one_proteome`).
- `ws:tests/.../test_wall_fusion_report.py:266-273` (`test_best_score_kept_for_repeated_id`) tests the merge but not a warning, so it
  locks in the silent behaviour. There is also no `.zst` test and no test for a model without a GA line, though both paths are in the
  code (`:140-150`, `:188-190`).

### F7. CONCERN: no runtime check of model names; the test's "verified" list is circular
- `extract_models` reads NAME from Pfam (`:176-177`) and the output prints it (`model_name`), but nothing compares it to the TSV `name`.
  A renumbered or renamed accession in a later Pfam release would be scanned silently. Plan 04 s4.4 asked for exactly this drift test
  (main `04:209`).
- `VERIFIED_05_WALL` / `VERIFIED_07_FUSION` (`ws:test...py:27-46`) are hand-copied dicts. `test_model_tsv_only_verified_ids` checks
  the TSV against them. `test_model_sets_and_status_text` checks only that the status text, which the same author wrote, starts with
  "05 MATCH". Neither parses 05 or a Pfam file. The fake Pfam in the tests is built from the TSV itself (`:87-96`), so names always agree.
- I checked by hand: all 16 wall names equal the 05 Table A short names (05:28-42 plus PF02364 at 05:28), and PRM1/Fig1 equal
  `full_pfam_ga_hits.tsv:27-28`. So the list is correct today, but the test guards against edits only, not against error.

### F8. CONCERN: model-list decisions use inconsistent criteria
- Fusion models PF27758 and PF12351 are included on the strength of the 07 Pfam scan alone ("NOT in 05", models TSV rows 18-19).
  PF29664 (BtgE-BtgC: SCW4 283.2, SCW11 266.9, BGL2 224.2) and PF05390 (Kre9_KNH1_C: KRE9 166.2, against weak Kre9_KNH 46.7) are
  excluded because they are "not in the 05 verification" (`ws:README.md:39`). That is the same evidence level with the opposite
  decision. PF29664 is the only family that catches SCW11, a soluble wall glucanase in report 02 (`err:02:84`).
- PF01522 is excluded because the "CE4 mapping is UNVERIFIABLE" (`ws:README.md:30`). But 05 verifies the ID and name (05:35); only the
  CAZy mapping is unverified. GH76 is included even though 05 Table B also marks the Dfg5/Dcw1 link UNVERIFIABLE (05:90); 07 measured
  hits for it. The md note "any ID that 05 marks UNVERIFIABLE" (`ws:wall_fusion_report.py:449-454`) misstates 05: 05 marks claims
  UNVERIFIABLE, not IDs (05:104, all 53 resolve). A defensible reason to drop PF01522 would be "large family, no test protein".
- PF12454 Ecm33 is kept although its namesake S. cerevisiae protein has no hit, even at E 10 (07 row 22). The TSV and README disclose
  this (`ws:wall_fusion_models.tsv:14`, `README.md:42`). Keeping it is acceptable for a report, but the `family` label "Ecm33" reads as
  if a hit means Ecm33; known-member sensitivity is 0 of 1. Suggest the label "Ecm33-like (PF12454; Sc ECM33 not detected)", or dropping it.
- PF07983 X8 is not an enzyme family. Kre9_KNH, SUN, Tos1 and Ecm33 are not shown to be enzymes either, and the chitin and glucan
  synthases are synthesis, not remodeling. All of them sit in set `wall_remodeling`, so the per-set "distinct proteins" count
  (`:420-438`) can count an X8-only protein (X8 is also in plant glucanases) as a "wall_remodeling" protein. Suggest naming the set
  `wall_synthesis_remodeling_context`, or reporting X8 only together with GH72.
- Exclusions that are justified: SH3_1 (07 row 26: weak 23.5, generic domain); WSC (sensor, also in secreted enzymes; main `04:136`);
  GPI-anchor synthesis (ER membrane enzymes; main `04:138`); PF00026 (generic; main `04:135` says it needs curated examples);
  PF00187, PF01607 and LysM (already in the lectin report); PF10333 (the record is Pga1, an ER GPI enzyme, so it falls under the GPI
  exclusion anyway).

### F9. CONCERN: label wording names a function where the lectin precedent names the model
- `model_label` gives `"<set>_family_hit: <family>"` (`ws:wall_fusion_report.py:283-285`), for example `cell_fusion_family_hit: Prm1`
  and `wall_remodeling_family_hit: X8`. The lectin label is `"<Name> (<PFxxxxx>) domain hit"` (main `lectin_report.py:373-376`), which
  names the model only.
- Plan 04 s4.3 item 1 asked for fixed text per row ("It does not show that the protein remodels the wall ..."). The script puts the
  disclaimer only in two `#` comment lines and in the md. A `comment="#"` reader, the documented way to read the table, drops it, and
  then sees "cell_fusion" in the `set` and `label` columns of every fusion row.
- CLASSES.md terms (main `docs/CLASSES.md:5,13`): `_domain` is an evidence call and `_candidate` a composite. `_family_hit` collides with
  neither: good. But `cell_fusion_...` reads like a call name in that scheme. Suggest `label = "PRM1 (PF27758) family hit"`, as in the lectin script.
- `ws:README.md:3` says "It has no threshold", but it applies the Pfam GA cutoff (`:54`, `:393-394`). The intended meaning is "no call
  threshold of its own". Reword.

### F10. MINOR: PRM1 475.9 vs 475.5 is explained: a column difference, not a different run
- 07's "Other Pfam hits" column comes from `full_pfam_ga_hits.tsv`, whose score column is **`best_domain_score`**
  (`analysis/cell_wall_check/make_summary.py:46`). The underlying hmmscan line (`fullpfam.domtbl:50`) gives sequence score **475.9**
  and domain score **475.5**.
- The script reports seq_score 475.9 and domain_score 475.5 (`ws:results/wall_fusion_hits.tsv:17`). Same HMMER 3.4, same Pfam 38.2, same
  values. Only the E-values differ (hmmscan 2.8e-142 vs hmmsearch 2.1e-145), as expected from the different search space.
- Consequences: `ws:README.md:79` ("cause was not checked") can be closed. 07:20 says "Scores are hmmsearch sequence scores", but the
  "Other Pfam hits (full scan)" column holds hmmscan best-domain scores. This mislabel in 07 affects PRM1, the PF29664 values, AGA2 and
  others; it is harmless for single-domain hits. The models TSV `supporting_report` quotes "Prm1 475.5" (the domain score) next to a
  table that shows 475.9.
- Inference on `--cut_ga`: the script counts a protein when its sequence score is at or above GA (sequence). Pfam's own pipeline also
  needs one domain at or above GA (domain). `domain_ge_ga` exposes this per row, but `n_hit_proteins` counts both kinds. This is
  disclosed (`:393-396`) and acceptable.

### F11. MINOR: 07 internal points
- 07 row 10: "PA14 ... does not hit BGL2, EXG1, SCW4 or SUN4". PA14 occurs in GH3 beta-glucosidases, and none of these four is GH3,
  so the negative tests nothing.
- 07's derived sorting-hat labels hold up (Q2). R0 and repeat calls do not depend on the family table, so the older "only PA14 active"
  run is a valid source for those two leaf calls. The derived labels add the fresh Pfam hits and use the current `categories.yaml` rules
  (`categories.yaml:59-63,101-106`). The run's config hash was not compared with the current config, so for an exact match: CANNOT TELL.
  9 of 23 proteins are not in that run, and 07 marks them "not run" / "n/m".

---

## Answers to the five questions

1. Errata fidelity: **CONCERN**, with DEFECT parts. Most edits are supported by the cited 05, 06 or 07 sections, and the spot checks match:
   TSV adhesion 45 rows / 23 clusters; Canalb1 11 + Sacce1 8; cluster 209 node; gel4 annotation; summary.tsv scores; active flags for the
   Als rows. Problems: the wrong commit for the Flocculin rows and the missed PIR1-like_C and Hyr1 statuses (F1); several 07-contradicted
   claims left unmarked (F2); the Gel4 doubt spread into 01 (F4); "verified by hmmsearch" for the fusion IDs, a CONFIRMED claim tagged
   "corrected", and two replacements logged as "added marker" (F5). No open question was silently resolved; the edits lean the other way.
2. Cross-report logic: **CONCERN** (F2, F3). The ESM-classifier risk holds. The sorting-hat false-positive claims in 03 s3.5 do not.
   Fusion and multicellularity conclusions: SOUND as worded. 07's derived labels: SOUND, with the config hash CANNOT TELL (F11).
3. Scripts: **SOUND** on no call, no status and no extra files (only three `write_atomic` calls; test `:230-252`; columns contain no
   call/status). `--cut_ga` on a combined model file is correct, because each model uses its own GA. Empty FASTA, missing model, missing GA,
   duplicate proteome names and exit code 2 are handled. **CONCERN**: no duplicate-ID warning (F6) and no name drift check (F7).
   PRM1 difference: explained (F10).
4. Model list: **CONCERN** (F8). The criteria are inconsistent (fusion IDs vs PF29664/PF05390; the PF01522 reason), Ecm33's label, and
   X8 and the non-enzymes inside `wall_remodeling`. The test list is circular but correct today (F7).
5. Wording: **CONCERN** (F9). The header, README and md disclaimers are clear and say "NOT evidence of adhesion, of cell fusion ...".
   The per-row `label` and `set` values carry function words, unlike the lectin c-i label.

## Things that are sound
- 05: a careful method. All 53 IDs resolve. The chitin synthase swap is placed on the brief, not the repo (05:18). Claims and IDs are kept apart.
- 06: the 23/45 reconciliation is correct, and the spot checks reproduce its counts. The "four distinct sets" table (06:28-35) is useful.
- 07: clear limits section (07:112-118). It states the data source for every column. Measured and derived results are kept apart. Its
  scores reproduce in the script output for every shared model (GAS1 504.4, X8 102.4, DFG5 530.8, DCW1 547.4, BGL2 353.2, SCW4 45.6,
  CRH1 170.8, UTR2 159.5, CTS1 81.6, EXG1 43.7, SUN4 330.1, KRE9 46.7, FIG1 192.3).
- Errata markers are dated, point to a section, keep the old text visible (except F5), and each file ends with a "still unverified" list.
- Script: option c-i is respected. Nothing under `src/`, `family_table.tsv` or `categories.yaml` is touched (diff stat: only
  `analysis/wall_fusion_report/` and `tests/wall_fusion_report/`). Atomic writes, a temp dir under `--tmpdir` or `$SCRATCH`, the sbatch
  takes every path from the environment, and the no-call notice is in both tables and the md. 15 tests pass. They need no network,
  HMMER or Pfam.
