# PR 6 ecm33 independent-set evaluation: the custom Ecm33 model against PANTHER

*2026-10-24. Author: Claude Code (claude-sonnet-5-5). Run record: `docs/reports/data/pr6_ecm33_independent_eval_2026-10-24/` (README there lists jobs, pins and files). Freeze: `analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json` (signed 2026-10-23, sha256 `22e6727879f9847b166b3ccdc75ceb702a87589b504983c5208589ee1dffc270`, merged #130). Design: `analysis/wall_fusion_truth/README_pr6_ecm33_eval.md`. Earlier build record: `2026-10-23-pr6-ecm33-build.md`.*

Status: the verdict is `not testable (u < 3): u = 1`; every other cell is `headline` or `descriptive`, none is `estimated`. The custom model stays report-only. No call, tier, panel, `narrow_ok` or status changes with this record.

## 1. What was run

The test set is the UniRef50 clusters of the PANTHER Ecm33 families (187 clusters, 725 proteins), searched with the frozen custom Ecm33 `all` model (`data/custom_models/ecm33_all_2026-10-23/`, cutoff 29.6 bits) and five comparators. PANTHER is itself a profile HMM, so recall here is agreement between two profile models, not truth and not adhesion.

Jobs (code `a8d0876`): extract 30336905 (COMPLETED, 15 s), search 30343079 (COMPLETED exit 0, 1 h 40 min), report run locally. Settings: cutoff 29.6 (frozen), relaxed Pfam 15.0, `min_u` 3, hit rule E <= 1e-5 for the iterative and phmmer comparators, `-Z 70198728`, MMseqs2 `b804fbe3` at 30 percent identity / 50 percent coverage.

**Why the first extract failed.** Job 30324005 (code `e81f3bd`) stopped with FATAL: it found 51 of 187 representatives in the HPCC UniRef50 file (release 2025_03). The 187 ids came from the live UniProt view of 2026-10-23, so 136 do not exist in that older file. No result was produced. The fix (#133): the representatives are fetched live by id from the UniProt REST service (release 2026_03, 187 requested, 187 received), staged as `reps.faa` with sha256 `493d386e31a87e452deeb5327fd20d4850d3abbce977ff6595c7248049605586` (sidecar and provenance committed; the FASTA is not). Caveat: a live representative can differ from the 2025_03 one, so the test set is the live PANTHER / UniRef50 view of 2026-10-23.

## 2. Separation from training

Of the 187 test clusters, 66 are separated from the training members at 30 percent identity / 50 percent coverage (no training member in the same MMseqs2 cluster). At the strict PANTHER cut (E <= 1e-10) there are 109 clusters (563 proteins); only 15 of them (22 proteins at the strict cut; 33 proteins in those clusters in all) are separated. These 15 are the HEADLINE set. The other 94 of the 109 share a 30 percent cluster with training members, because Ecm33 is conserved.

The 15 headline clusters: clades Chytridiomycota 4, Entomophthoromycota 2, and one each in nine others (Dothideomycetes, Eurotiomycetes, Glomeromycota, Leotiomycetes, Mortierellomycota, Mucoromycota, Orbiliomycetes, Sordariomycetes, Zoopagomycota); 13 tier T3_panther_only and 2 T2_named; all subfamily SF3.

## 3. Result

Verdict line (whole headline set, `min_u` 3): `not testable (u < 3): u = 1`. u is the number of headline clusters missed by every comparator. Only one is, so the rule (custom recovers at least ceil(u/2) of u) cannot be applied.

Recall with exact two-sided 95 percent Clopper-Pearson intervals; the headline row and the two descriptive columns:

| Method | Headline 15 | 95% CI | Strict 109 (descriptive) | All 187 (descriptive) |
|---|---|---|---|---|
| Custom `all` model | 14/15 | 0.681-0.998 | 108/109 | 185/187 |
| Pfam GA (PF12454) | 1/15 | 0.002-0.319 | 27/109 | 28/187 |
| Relaxed Pfam (15.0 bits) | 1/15 | 0.002-0.319 | 40/109 | 48/187 |
| Alignment-seeded iteration | 14/15 | 0.681-0.998 | 108/109 | 185/187 |
| Seed-only phmmer | 5/15 | 0.118-0.616 | 99/109 | 124/187 |
| Seed-only jackhmmer | 14/15 | 0.681-0.998 | 108/109 | 185/187 |

Intervals of the descriptive columns: custom and jackhmmer 108/109 0.950-1.000, 185/187 0.962-0.999; Pfam GA 27/109 0.170-0.340, 28/187 0.102-0.209; phmmer seed 99/109 0.838-0.955, 124/187 0.591-0.730.

Headline by tier: T3_panther_only 13 clusters, custom 12/13 (0.640-0.998), Pfam GA 1/13, jackhmmer 12/13, phmmer 3/13; T2_named 2 clusters, custom 2/2, Pfam GA 0/2, phmmer 2/2. Headline by clade (custom, Pfam GA, phmmer seed): Chytridiomycota 4/4, 0/4, 0/4; Entomophthoromycota 2/2, 0/2, 2/2; every one-cluster clade is 1/1 for the custom model except Mortierellomycota 0/1; Pfam GA hits only the Sordariomycetes cluster. Every clade, tier and subfamily row with its interval (headline, strict 109, all 187) is in `report/summary.txt` and `report/groups.tsv`. In the 187-cluster column jackhmmer seed equals the custom model in every clade (Entomophthoromycota 21/21 for both, phmmer seed 7/21, Pfam GA 0/21).

The u clusters: the u = 1 headline cluster is `UniRef50_A0A9P5VCA1` (Mortierellomycota, T3_panther_only, SF3, PANTHER E 5.7e-17), missed by every method. The second u cluster of the 187 is `UniRef50_A0A9W6T7D3` (Saccharomycotina, T3_panther_only, SF3, PANTHER E 0.0011, outside the strict set). Both are one-protein clusters.

## 4. What the numbers say

1. The alignment-independent seed-only jackhmmer (the cheap method the owner asked to test) recovers as much as the custom model: 14/15, 108/109 and 185/187. So the custom model shows no recall advantage over a jackhmmer seeded by the single Kiss protein. The alignment-seeded iteration gives the same counts.
2. One-pass seed-only phmmer is clearly weaker: 5/15, 99/109, 124/187.
3. Pfam PF12454 reaches only 28/187 (GA) and 48/187 (relaxed). So the custom model and jackhmmer both fill the Pfam gap.
4. The custom model is a fixed, versioned, fast artifact (hash-pinned, one `hmmsearch -T 29.6`), whereas the jackhmmer result depends on the target database and the iterations.
5. PANTHER is a profile HMM: agreement with it measures agreement between two profile models, not whether a protein is Ecm33 and not whether it is an adhesin.
6. The headline is only 15 clusters, so intervals are wide (14/15 is 0.681-0.998). One cluster (A0A9P5VCA1) is missed by every method, and u = 1 leaves the freeze's rule untestable.

## 5. Caveats

- A live UniRef50 representative may not carry the PANTHER call of its cluster and may differ from the 2025_03 one.
- The strict and all-187 columns are mostly close relatives of training (94 of 109 strict clusters share a 30 percent cluster), so they test coverage, not extrapolation.
- `purity.json` (descriptive): 164 of 481 training members are in no PANTHER-called cluster, a disagreement between two profile models.
- PANTHER and UniProt counts are of 2026-10-23 and change.

## 6. Questions for the owner (not decisions)

1. Is the custom Ecm33 model worth keeping when seed-only jackhmmer ties it? The value left is a fixed, fast, versioned artifact.
2. Should the 94 non-separated strict clusters be examined? They are close relatives of training, so they test coverage, not extrapolation.
3. Should the same test design be applied to other PANTHER families (Chs7, Chs6, Chs5, Big1, Kre6 from `2026-10-23-panther-family-scoping.md`)?

No tier, panel, call, `narrow_ok` or status change is made by this record.
