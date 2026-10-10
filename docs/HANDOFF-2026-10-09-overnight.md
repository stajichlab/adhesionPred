# Hand-off: overnight work, 2026-10-08 to 2026-10-09

*Branch `fungi-scan-and-specs` (from `hydrophobin-validation`; local until the owner says otherwise). Owner instructions: scan proteomes from Fungi_5k, 40 to 60 jobs on `short` and `short_gpu`; stop at any owner-decision point and at any failed review gate; commit as tasks complete; write this file. No merges, no deletions, no new PRs.*

## Status log (newest last)

- 2026-10-08 late: `surface_attachment_candidate` call added and pushed to `hydrophobin-validation` (`03dcc23`).
- Descriptive scan set up: 56 proteomes from Fungi_5k (`analysis/fungi_scan/selection.tsv`; 24 Onygenales, 12 Aspergillus, 3 Candida, others). 15 grouped SLURM jobs submitted (`_workdir/fungi_scan/jobs.txt`).
- M2 class list: spec rev 2 (review: 0 blockers, 8 major, all addressed); M2a to M2f implemented and committed (`33aa401`, `4df020d`, `23c923b`, `da8032c`, `bcc8e51`). Repeat call re-measured after the config change: numbers unchanged (S288C 0.435 [0.067, 0.667] / 0.987; *C. albicans* 0.462 [0.000, 0.788] / 0.976). Last full run: 944 passed, 9 skipped, plus 7 class-table tests.
- M5 SOWgp-like spec: rev 1 failed review 1 (2 blockers, 9 major); rev 2 and rev 3 written. Review 2: PASS (0 blockers, 3 major, 8 minor). Rev 3 addresses the 3 major. The 8 minor findings are not applied (`...-sowgp-like-calls-design-review-2.md`). Commits `f444f97`, `df0f767`. No M5 code written.
- Task 07 (shared hard-negative set): an agent built 161 rows in 91 clusters (`data/controls/hard-negatives-shared/`, commit `9962ea7`). Not owner-reviewed. Size target mostly not met: only `laccase_etc` has 20 or more clusters per stratum; Basidiomycota 14, Onygenales 15, Taphrinomycotina 5 clusters. Six rows share a cluster with a positive (flagged, not removed).

- Descriptive scan finished: all 15 jobs completed; converters ran on 56 proteomes without error; report `docs/reports/2026-10-09-fungi-scan-descriptive.md` with tables in `docs/reports/data/sorting_hat/fungi_scan/` (1,053 surface attachment candidates). The SOWgp unit HMM hits only the two *Coccidioides* proteomes.

## Owner decision points reached (work stopped here)

**M5 (needed before M5b):**
- O1: `pro_cys_array_protein` as sub-label of `surface_attachment_candidate`; `sowgp_ortholog` in the `other_*` mechanism lists. Default: yes / yes.
- O2: keep Pro+Cys over 15% and Cys at least 4% (Pro can be low), or add a proline floor. Default: keep. The source of these thresholds is not recorded.
- O2b: detector `repeat02 or repeat14` or `repeat14` only. Default: either; M5b reports both.
- O3: drop BAD1 (11.8% Pro+Cys fails the rule) or change the rule. Default: drop.
- O4: the 20 class 2a candidates and the 74 look-alikes as T4 outcomes only. Default: reported only.

**M2 (O1 to O5):** see `docs/superpowers/specs/2026-10-09-v1-class-list-design.md` section on decisions.

**Task 07 questions from the agent:**
1. Should Msb2 return as a flagged hard negative? (The *A. nidulans* ortholog MsbA affects adhesion, PMID 25294314; the agent rejected the whole family.)
2. Is `serves=step1` read correctly (signal-peptide or membrane protein not at the surface; 9 rows), with `adhesion_level` for the rest?
3. Keep, rename or drop the added strata `secretory_non_surface` and `wall_structural_other`?
4. Move the two PA14 beta-glucosidase positives (A0A1D8PJB5, A0AAW0V6D7) to E3? They are labelled adhesin by the PA14 domain only.
5. Do N2 rows count toward the 20-cluster floor?
6. Is pooling strata into one negative set for one call in one species accepted?

## Failed review gates

None. M5 failed review 1 and passed review 2 after revision.

## Next steps

1. Read the scan report. Say whether to keep the tables in the repository (about 1,500 lines).
2. Owner answers O1 to O4 (M5) and the task 07 questions.
3. Then M5a (plan first, review before code).

## Owner decisions of 2026-10-09 (answered in the session)

| Item | Decision |
|---|---|
| M5 O1 | `pro_cys_array_protein` is a sub-label of `surface_attachment_candidate`; `sowgp_ortholog` joins the `other_*` mechanism lists. |
| M5 O2 | Keep Pro+Cys over 15% and Cys at least 4% (Pro can be low). |
| M5 O2b | Detector: `repeat02 or repeat14` (both). |
| M5 O3 | BAD1 dropped; BAD1-type arrays are outside the call. |
| M5 O4 | The 20 class 2a candidates and the 74 look-alikes are reported only. |
| M2 O1 | (a) CFEM stays in `cell_wall_adhesion_candidate`; the family is shown. |
| M2 O2 | Five families went active: PF00624, PF13928, PF15789, PF22799, PF30910. PF30910 was turned off again the same day after the family test; four stay on. The other nine stay inactive. **Next: test all 14 and decide how to incorporate them (owner, "move next to testing all of these").** |
| M2 O3 | (a) `iuis_allergen_homolog` keeps its Pfam branch; may be `not_assessable`; antigen and allergen categories stay in the report as `unvalidated`. |
| M2 O5 | No "derived" marker in `status_basis` for now. |
| Task 07 Q1 | Msb2 family stays rejected. |
| Task 07 Q2 | `serves=step1` reading accepted (9 rows); other rows use `adhesion_level`. |
| Task 07 Q3 | Keep `secretory_non_surface` and `wall_structural_other`. |
| Task 07 Q4 | A0A1D8PJB5 and A0AAW0V6D7 move to E3 (done in `data/curated/adhesins/adhesins.tsv`; the two rows only). |
| Task 07 Q5 | Only N1 rows count toward the 20-cluster floor; N2 are reported separately. |
| Task 07 Q6 | Pooling strata into one negative set per call per species is accepted. |

**Consequence of Q5 (recount, `controls.tsv` and `clusters.tsv`):** N1 only gives 137 rows in 76 clusters (not 91). Per stratum: laccase_etc 23, wall_hydrolase 13, domain_non_member 9, gpi_wall_enzyme 8, st_linker_enzyme 7, repeat_non_adhesin 6, wall_structural_other 4, mucin_sensor 3, secretory_non_surface 3. Per species: *S. cerevisiae* 16, *A. fumigatus* 15, *A. nidulans* 10, *V. dahliae* 9, *T. rubrum* 8. **No species reaches 20 N1 clusters, so under the accepted rules no pooled species entry can reach `estimated` yet.** More N1 negatives are needed (or the owner reconsiders Q5).

## Next steps (updated)

1. Test all 14 inactive families (descriptive counts on the scan proteomes first, then truth-based tests where truth exists) and decide how to incorporate each (M4).
2. M5a plan (review before code), then M5a.
3. Re-measure the repeat call when `categories.yaml` next changes (M5b2).
4. Hard-negative set: add N1 rows for species near 20 clusters (*S. cerevisiae* 16, *A. fumigatus* 15).

## Family tests (M4) run 2026-10-09 after the owner request

Report: `docs/reports/2026-10-09-family-tests.md`; code `analysis/family_tests/`; table `docs/reports/data/sorting_hat/family_tests/family_tests.tsv`. Key facts for the owner:
- **PF30910 ALS_M (turned on today) does not hit the Candida Als positives** (0 of 14 with a Candida_ALS_N hit, even at E 1e-3). **Decision (owner, 2026-10-09): turned off again** (`family_table.tsv`; the other four of the five stay on).
- **PF22799 PIR1-like_C (turned on today) hits 3 N1 and 1 N2 hard negatives** and no E1 positive.
- **CFEM hits 12 N1 negatives in 9 of 76 clusters** (and 1 of 54 E1 positives): new evidence for M2 O1.
- Hyr1 has no E1 positive to test it; the allergen families need an IUIS test; Candida_ALS_N, Candida_ALS, Flo11 and GLEYA are the next candidates to turn on.

**Owner decision 2026-10-09 (after the family test):** PF30910 ALS_M off; PF11766 Candida_ALS_N, PF05792 Candida_ALS, PF10182 Flo11 and PF10528 GLEYA on. Active families now: 10 + PF00624, PF13928, PF15789, PF22799 + these four = 18 of 24. Inactive: Bys1 PF04681, Hyphal_reg_CWP PF11765, PIR PF00399, ALS_M PF30910, AltA1 PF16541, Allergen_Asp_f_4 PF25312. Known at the time of the decision: GLEYA hits 1 N2 negative; the four families overlap the repeat call in 3 to 15 scan proteins each. The repeat call's status files are not stale (only `categories.yaml` changes the hash), but the pfam family digest changed; the Pfam-reading calls need a new measure when next run in a real work directory (not done).
