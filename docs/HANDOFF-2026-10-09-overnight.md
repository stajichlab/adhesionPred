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
