# Hand-off: hydrophobin validation (2026-10-08)

Branch `hydrophobin-validation` (local, not pushed; it contains `per-call-status` merged in `f005e20`). Read in this order:
`docs/reports/2026-10-08-hydrophobin-validation.md`, the spec and plan in `docs/superpowers/`, and the three review files.

## Done
- Spec rev 3, plan rev 2, three independent reviews (spec x2, plan x1), all findings applied.
- H0 (PF01185 note, commit `1466f90` on `signoff-cfem-hydrophobin`, local), H1c (per-module Pfam digest), H3/H3b (frozen spacing, prior art), H2/H2b/H2c (truth set, mapping, calibration inputs, stored split), H4 (`cys8_pattern`), H1 (module split, `hydrophobin_domain`, `hsba_domain`), H5/H5b (runs, 12 proteomes), H6 (measurement, 8 status files), H1b (repeat call re-measured, same numbers).
- Tests: `tests/cellsurface_sorting_hat`, `tests/hydrophobin_truth`, `tests/calibration_truth` pass.

## Results in one paragraph
All hydrophobin statuses are `smoke`. The 8-cysteine rescue with the published spacing does not help (rule 6.3: no evidence). The published spacing matches 60% of Swiss-Prot hydrophobins. `hydrophobin_protein` was not added.

## Owner decisions waiting
1. Rescue: leave out (default) or build a data-derived spacing from the 122 Pfam-positive T2 entries, freeze it, and test it on the 9 no-model entries (report section 8).
2. PRs: #75 (CFEM and hydrophobin activation) puts hydrophobin hits in `wall_family_domain`. The module split is on this branch. Merge order: #76, #75, then this branch (it already contains #76).
3. Review the 22 unlabelled Pfam hydrophobin calls (`docs/reports/data/sorting_hat/hydrophobin/review_lists.tsv`, list `pfam_call_not_positive`) to add labels.
4. HsbA (`hsba_domain`) has no truth set. Decide whether it is worth one.

## Not done
- Cross-species hard-negative table, wider prior-art search beyond PubMed and the sources named in `docs/paper/05` section 6, Wessels 1994 / Linder 2005 / Sunde 2008 not read.
- No push, PR or merge since the owner authorised #75 and #76.
