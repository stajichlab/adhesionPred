# Baseline recall of the Pfam models on the Kiss 2019 control set (Task 5.4)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`, Task 5.4. Run record: `docs/reports/data/kiss_baseline_2026-10-22/` (job 30233470, code `9472e32`, HMMER 3.4, Pfam 38.2). Every number below was recomputed from the committed tables.*

## 1. What this measures, and what it does not

- It is **family coverage by Pfam models**: is a control protein hit by `hmmsearch --cut_ga` of the models of its slice (adhesion rows: the 18 active family rows; wall rows: the 21 wall models; fusion rows: the fusion models among the 21)? It is **not the sensitivity of the sorting-hat call**. It changes **no tier, no status and no `narrow_ok`** (plan Task 5.4; D4 is an owner decision).
- The unit is the MMseqs2 cluster of Task 5.3 (30 percent identity, 50 percent coverage). A cluster counts as recalled when at least one control member is hit. Rows are given beside clusters.
- **`seed-like` is nearly circular.** A cluster is `seed-like` when a member reaches 30 percent identity to a Pfam seed, so a Pfam model finding it is close to guaranteed. The informative figure is **`seed-independent`**: clusters with no seed relative.
- The 95 percent Wilson interval treats clusters as independent. It is a proxy.
- There is **no pooled number** across slices. Wall rows are a descriptive panel and are never negatives and never truth for the adhesion call (D20). The 35 family representatives are descriptive and never counted.
- 39 models fetched, 39 with a GA line, none without. The Pfam file has the same sha256 as in the earlier jobs.

## 2. Recall per slice

Clusters hit, per slice and seed class (own-slice label: `slice`; label against all 39 models: `all39`). Status: `estimated` needs at least 5 clusters in the cell (decision 4), otherwise `descriptive`; `none` has 0 clusters.

| Slice | Seed class | Label | Clusters hit / n | Rows hit / n | 95% Wilson | Status |
|---|---|---|---|---|---|---|
| adhesion | seed-like | slice | 16 / 16 | 24 / 24 | 0.806 to 1.000 | estimated |
| adhesion | seed-independent | slice | **0 / 12** | 0 / 12 | 0.000 to 0.242 | estimated |
| adhesion | seed-like | all39 | 16 / 16 | 24 / 24 | 0.806 to 1.000 | estimated |
| adhesion | seed-independent | all39 | **0 / 12** | 0 / 12 | 0.000 to 0.242 | estimated |
| wall | seed-like | slice | 14 / 14 | 39 / 39 | 0.785 to 1.000 | estimated |
| wall | seed-independent | slice | **2 / 30** | 3 / 52 | 0.018 to 0.213 | estimated |
| wall | seed-like | all39 | 14 / 17 | 39 / 42 | 0.590 to 0.938 | estimated |
| wall | seed-independent | all39 | 2 / 27 | 3 / 49 | 0.021 to 0.234 | estimated |
| fusion | seed-like | slice and all39 | none | none | none | none |
| fusion | seed-independent | slice and all39 | 0 / 4 | 0 / 5 | 0.000 to 0.490 | descriptive |

- Of the adhesion clusters, 16 have a seed relative and are all hit; the 12 with no seed relative are all missed.
- Of the wall clusters, 14 of 14 `seed-like` are hit and 2 of 30 `seed-independent`. Against all 39 models the `seed-like` count is 17: three wall clusters are close only to a seed of an adhesion-family model, and the own-slice models do not hit them (14 of 17).
- Fusion has 4 clusters (5 rows), all `seed-independent`, none hit. This is descriptive only (below the floor of 5, and decision 5).

## 3. Recall per clade

Fine level (primary), stratum `any`, own-slice label. A cluster that spans several clades counts in each of them, so the clade cells do not add up to the slice totals.

| Slice | Clade | seed-like hit / n | seed-independent hit / n |
|---|---|---|---|
| adhesion | Saccharomycotina | 10 / 10 (estimated) | 0 / 4 (descriptive) |
| adhesion | Sordariomycetes | 3 / 3 (descriptive) | 0 / 4 (descriptive) |
| adhesion | Eurotiomycetes | 1 / 1 (descriptive) | 0 / 2 (descriptive) |
| adhesion | Ustilaginomycetes | 2 / 2 (descriptive) | 0 / 1 (descriptive) |
| adhesion | Tremellomycetes | none | 0 / 1 (descriptive) |
| wall | Saccharomycotina | 10 / 10 (estimated) | 1 / 22 (estimated) |
| wall | Eurotiomycetes | 8 / 8 (estimated) | 1 / 9 (estimated) |
| wall | Sordariomycetes | 7 / 7 (estimated) | 1 / 8 (estimated) |
| wall | Taphrinomycotina | 1 / 1 (descriptive) | 0 / 1 (descriptive) |
| wall | Chytridiomycota, Tremellomycetes, Ustilaginomycetes | none | 0 / 2, 0 / 1, 0 / 1 (descriptive) |
| fusion | Eurotiomycetes, Saccharomycotina, Sordariomycetes, Ustilaginomycetes | none | 0 / 1, 0 / 2, 0 / 1, 0 / 1 (descriptive) |

- At the fine level no adhesion `seed-independent` clade cell reaches 5 clusters: the 12 clusters are spread over five clades (4, 4, 2, 1, 1). The clade-level statement is therefore descriptive for adhesion.
- At the coarse level (second test), adhesion Pezizomycotina `seed-independent` is 0 of 6 (estimated); Saccharomycotina has 0 of 4 (descriptive). Wall coarse: Pezizomycotina 1 of 14, Saccharomycotina 1 of 22 (both estimated).
- `recall_by_clade.tsv` has 190 cells, 60 estimated and 130 descriptive; it also has the `clade_only` and `shared` strata (decision 3).

## 4. The clusters Pfam misses (candidates for PR 6)

`clusters_missed.tsv` has 44 rows: **adhesion 12, wall 28, fusion 4**. All are `seed-independent` by the own-slice label.

- **Adhesion (12, all clade-only):** Saccharomycotina 4 (aga1 in *S. cerevisiae*; ece1, iff4, dfi1 in *C. albicans*); Sordariomycetes 4 (ncl1 to ncl4, *N. crassa*); Eurotiomycetes 2 (hyd1_like, *H. capsulatum*; glfB, *A. fumigatus*); Ustilaginomycetes 1 (rep1, *M. maydis*); Tremellomycetes 1 (cfl1, *C. neoformans*). The evidence flags: row-verified 5, portal-validated 4, locus-confirmed 2, one row without a flag.
- **Wall (28):** 20 clade-only and 8 shared. By clade: Saccharomycotina alone 14, Eurotiomycetes alone 3, Sordariomycetes alone 3, Taphrinomycotina 1 (*S. pombe*), and 7 clusters that span two to five clades. By species, the clusters that are in one species are *C. albicans* 8, *S. cerevisiae* 4, *N. crassa* 3, *A. fumigatus* 3, *E. gossypii* 1, *S. pombe* 1.
- **Fusion (4):** ste2 (*S. cerevisiae*), ste3 with preA (*S. cerevisiae* and *A. fumigatus*), pra1 (*M. maydis*), so_soft (*N. crassa*).

This list names the candidates. It does not say a custom model will succeed; PR 6 needs its own spec and leave-clade-out.

## 5. Wall rows against the adhesion families (descriptive false-hit panel)

The 91 wall control rows (44 clusters) against the 18 active family models: **three rows hit one family each**. They are never negatives and never truth for the adhesion call.

| Family | Pfam | Wall row hit | `narrow_ok` | tier |
|---|---|---|---|---|
| CFEM | PF05730 | Q4WLB9 | no | smoke |
| Flocculin_t3 | PF13928 | P40953 | **yes** | smoke |
| PIR1-like_C | PF22799 | Q59SF7 | no | untested |

The 15 other family models hit no wall row. These three rows are the three wall controls that the other-slice models hit; no adhesion or fusion control is hit by a model of another slice. Q59SF7 (pir1) and Q4WLB9 (cfmA) are hit only by the family model; P40953 (cht2) is also hit by the wall model PF00704. One of the three families, Flocculin_t3, is signed `narrow_ok=yes`. This note reports it and changes nothing; a `narrow_ok` change is the owner's decision.

## 6. The shared negatives (161 sequences, six of them `#neg` copies)

`negatives_false_hits.tsv` gives, for each of the 39 models, the hits on the 128 sequences in negative-only clusters and the 33 in mixed clusters. 17 of the 39 models hit at least one negative; 72 distinct negative sequences are hit by some model.

- For the **18 family models** the hits are false hits for the adhesion families: five models hit a negative. PF05730 (CFEM) hits 16 sequences in 12 negative-only clusters and 1 in a mixed cluster; PF07691 (PA14) 4 in one negative-only cluster; PF22799 (PIR1-like_C) 1 in a negative-only cluster and 3 in one mixed cluster; PF10528 (GLEYA) and PF12296 (HsbA) 1 each in negative-only clusters.
- For the **21 wall and fusion models** a hit on this set is not an error of the model: the shared negatives are hard negatives for adhesion, and the set holds cell wall enzymes. The largest counts are PF00704 (11 and 4), PF00722 (0 and 11), PF00150 (6 and 4) and PF03198 (0 and 4), as negative-only and mixed hits.
- The `#neg` copies are in mixed clusters, as decided. Filter on `cluster_composition` to separate them.

## 7. Family representatives (descriptive, never counted)

| Slice | Rows | Sequences | Hit by own-slice models |
|---|---|---|---|
| adhesion | 4 | 4 | 0 |
| wall | 28 | 23 | 16 |
| fusion | 3 | 3 | 1 |

A representative is the *A. fumigatus* ortholog, not the Kiss protein. It counts toward no floor and is not a test control.

## 8. Limits

- Family coverage is not call sensitivity; the call also uses repeat and other modules.
- Wilson intervals assume independent clusters; a bootstrap over clusters is the upgrade if an interval will be quoted.
- The clade cells overlap (a multi-clade cluster is in each clade).
- A hit is any row of the `--cut_ga` domain table; no other threshold was used.
- The model names were matched in the pressed database by the job (39 of 39 found).
