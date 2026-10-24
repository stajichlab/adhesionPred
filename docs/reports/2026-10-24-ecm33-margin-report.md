# Ecm33 custom model: how far scores sit from the cutoff

*2026-10-24. Author: Claude Code (claude-sonnet-5-5). Built only from tables already committed (no search, no HPCC, no SLURM). Numbers come from `analysis/wall_fusion_truth/ecm33_margin_tables.py` (writes the small tables in `docs/reports/data/ecm33_margin_report_2026-10-24/`) and are recomputed by `tests/wall_fusion_truth/test_ecm33_model_card.py`. Model card: `data/custom_models/ecm33_all_2026-10-23/MODEL_CARD.md`.*

**This report does NOT change the cutoff (29.6 bits), the model, any freeze, call, tier, panel or `narrow_ok`.** A real recalibration needs a new held-out set and a new signed freeze (section 7).

## 1. Words used (each defined once)

| Word | Plain meaning |
|---|---|
| Score (bits) | How well a protein fits the model. Higher is a better fit. A protein is "called" when its score is at or above the cutoff. |
| Cutoff | 29.6 bits for the `all` model. |
| Margin | Score minus cutoff. A large positive margin for a true family member and a negative margin for a false hit mean the cutoff is not delicate. |
| Held-out member | A family member kept out of a model's training, then searched for. |
| Unit | One held-out clade (a fungal group). Each unit has its own model built without that clade, with its own cutoff set by the same rule. The `all` model has no held-out unit. |
| False hit (extra call) | A hit in one of four proteomes (A. fumigatus, C. albicans, C. immitis, S. cerevisiae) that a search from Q5AGC4 does not label as a family member. The freeze allows at most 5 per 10,000 proteins. |

## 2. Summary

| Question | Answer |
|---|---|
| Lowest held-out member score of any unit model, `cluster30` (headline) | 86.9 (Eurotiomycetes), cutoff of that unit 31.1: margin 55.8 |
| Smallest margin of any unit and mode | 41.4 (Kickxellomycota, `cluster30`: lowest 87.3, cutoff 45.9) |
| True Ecm33s among the shared negatives, `all` model | P38248 387.1 (13.1 times the cutoff); Q4WNS8 444.2 (15.0 times) |
| False hits at or above 29.6 in the four proteomes | 4, all in A. fumigatus: 42.7, 34.1, 33.8 and 29.6 (margins +13.1, +4.5, +4.2, 0.0) |
| Highest false-hit score below the cutoff | 29.5, A. fumigatus (0.1 bit below) |
| Shared hard negatives called | 0 of every group (rate 0.0000, limit 0.05); best score 21.1 (8.5 below) |
| Independent evaluation, 187 test clusters | 185 at or above 29.6; lowest 41.7 (+12.1); 0 within 10 bits of the cutoff, 1 within 20; 2 misses |

## 3. Held-out members per unit (build record, unit models)

Columns per mode: the unit's cutoff, the lowest and the median score of its held-out members, and the margin of the lowest. `cluster30` is the headline mode (training members within 30 percent identity of a held-out member removed); `none` removes only the held-out members. Source: `ecm33_build/<unit>/<mode>/hits_custom.tsv` and `cutoff.json`. All held-out members are hits, so the lowest score is the one nearest the cutoff.

| Unit | Held out | cluster30 cutoff | lowest | median | margin | none cutoff | lowest | median | margin |
|---|---|---|---|---|---|---|---|---|---|
| Dothideomycetes | 106 | 32.9 | 93.9 | 270.2 | 61.0 | 19.8 | 98.5 | 274.9 | 78.7 |
| Eurotiomycetes | 55 | 31.1 | 86.9 | 230.2 | 55.8 | 25.3 | 103.2 | 227.3 | 77.9 |
| Glomeromycota | 8 | 37.0 | 118.7 | 273.5 | 81.7 | 36.5 | 125.6 | 282.5 | 89.1 |
| Kickxellomycota | 37 | 45.9 | 87.3 | 241.7 | 41.4 | 35.6 | 87.8 | 252.5 | 52.2 |
| Leotiomycetes | 38 | 36.5 | 92.8 | 257.0 | 56.3 | 28.0 | 116.5 | 255.2 | 88.5 |
| Mucoromycota | 33 | 39.3 | 151.0 | 319.6 | 111.7 | 37.4 | 149.8 | 331.3 | 112.4 |
| Orbiliomycetes | 12 | 38.2 | 156.8 | 242.9 | 118.6 | 33.9 | 165.1 | 247.1 | 131.2 |
| Pezizomycetes | 5 | 40.7 | 232.2 | 383.4 | 191.5 | 31.0 | 248.1 | 404.7 | 217.1 |
| Saccharomycotina | 107 | 40.4 | 103.5 | 285.4 | 63.1 | 39.8 | 116.0 | 308.5 | 76.2 |
| Sordariomycetes | 57 | 41.8 | 87.1 | 260.4 | 45.3 | 34.1 | 104.5 | 265.5 | 70.4 |
| Taphrinomycotina | 11 | 44.9 | 246.3 | 349.0 | 201.4 | 30.7 | 249.1 | 338.2 | 218.4 |

Reading: in every unit and both modes the lowest member sits at least 41 bits above the unit cutoff, and the median is 227 to 405 bits. This is the circular test of the build record (the held-out members come from the same search that supplied the training), so a large margin is expected and says the family is coherent, not that the model finds Ecm33 beyond a Q5AGC4 search. The unit cutoffs (19.8 to 45.9) differ from the `all` model's 29.6 because each is set from the false-hit side of its own model.

## 4. The two true Ecm33s among the shared negatives

P38248 (S. cerevisiae Ecm33) and Q4WNS8 (A. fumigatus ecm33) are in the shared negative set but are true family members, so the freeze leaves them out of the false-hit rate. `all` model (`ecm33_build/all/none/neg_hits.tsv`): P38248 **387.1**, Q4WNS8 **444.2**, against 29.6. In the 22 unit/mode rows (`report/ignored_negatives_line.tsv`) P38248 scores 314.9 to 407.3 and Q4WNS8 389.1 to 454.5. Both are far above any cutoff.

## 5. False hits in the four kept proteomes (`all` model)

The model searched the four proteomes (32,491 proteins; `proteome_hits.tsv` lists every hit down to 0.7 bits). The committed `cutoff.json` says 4 extra calls (A. fumigatus 4; C. albicans, C. immitis and S. cerevisiae 0) and limits of 4.82, 3.11, 4.96 and 3.36 per proteome (16.2455 in all). The labelled list is not committed, so the four extras are inferred: A. fumigatus has 6 hits at or above 29.6, two labelled (444.2 Q4WNS8, 199.4 Q4WMD5) and the 4 lowest are the extras.

| Proteome | Hits at or above 29.6 | Extra calls (cutoff.json) | The extras (score, margin) | Best hit below the cutoff | Hits 0 to 10 bits below |
|---|---|---|---|---|---|
| A. fumigatus Af293 | 6 | 4 (limit 4.82) | Q4WHJ6 42.7 (+13.1), Q4WIT2 34.1 (+4.5), Q4WUZ7 33.8 (+4.2), Q4WP42 29.6 (0.0) | 29.5 (Q4WPI5, -0.1), then 29.1 | 10 |
| C. albicans SC5314 | 3 | 0 (limit 3.11) | none | 20.8 (-8.8) | 1 |
| C. immitis RS | 1 | 0 (limit 4.96) | none | 27.5 (-2.1), 27.0, 26.6 | 8 |
| S. cerevisiae S288C | 4 | 0 (limit 3.36) | none | 14.1 (-15.5) | 0 |

The margin is thin on the false-hit side. The cutoff 29.6 is the score of the fourth extra call; the next hit, 29.5, would be a fifth in A. fumigatus, over that proteome's limit of 4.82, so it was left just below. Nineteen proteome hits score in the 10 bits below the cutoff (A. fumigatus 10, C. immitis 8, C. albicans 1). A different proteome, or a protein set with more diverged cell-wall proteins, could cross 29.6. The extra calls are not family members by the freeze's definition; none was checked by hand here (the earlier PANTHER and Pfam data do not say whether they are Ecm33-like).

## 6. Shared hard negatives and the independent evaluation

**Hard negatives.** Every one of the 9 hard-negative groups has call rate 0.0 for the `all` model (`cutoff.json`, limit 0.05; also 0.0000 in all 198 group rows of `report/hard_negative_calls.tsv`). The shared negatives that were scored at all are G2X962 (21.1, 8.5 below the cutoff) and P43497 (12.6, 17.0 below), besides the two true Ecm33s of section 4.

**Independent evaluation** (`pr6_ecm33_independent_eval_2026-10-24/search/hits_custom.tsv`; 187 test clusters; representatives scored with `hmmsearch -T 29.6`, so only hits are listed). Hits at or above 29.6: 185. Median 262.3, highest 477.1. Lowest scores: 41.7, 51.7, 53.3, 53.4, 55.6, 57.9, 62.7, 68.6.

| Distance above the cutoff | Clusters |
|---|---|
| less than 10 bits (29.6 to 39.6) | 0 |
| less than 20 bits (29.6 to 49.6) | 1 (41.7) |
| less than 45 bits (29.6 to 74.6) | 9 |
| misses (score below 29.6, not printed by `-T`) | 2: UniRef50_A0A9P5VCA1 (Mortierellomycota) and UniRef50_A0A9W6T7D3 (Saccharomycotina); every other method missed both too |

Histogram of the 185 scores (25-bit bins; each # is one cluster):

```
 25- 49 bits | # 1
 50- 74 bits | ######## 8
 75- 99 bits | ###### 6
100-124 bits | ################# 17
125-149 bits | ############ 12
150-174 bits | ################## 18
175-199 bits | ######### 9
200-224 bits | ########### 11
225-249 bits | ###### 6
250-274 bits | ####### 7
275-299 bits | ######### 9
300-324 bits | #### 4
325-349 bits | ####### 7
350-374 bits | ################ 16
375-399 bits | ######################### 25
400-424 bits | ########### 11
425-449 bits | ############# 13
450-474 bits | #### 4
475-499 bits | # 1
```

Two clusters are missed and the distribution has no pile-up just above 29.6: the hits are spread from 41.7 upwards. The two misses sit below the cutoff by an unknown amount. So a lower cutoff would not recover them without evidence, and the false-hit side (section 5) stops at 29.5 within 0.1 bit.

## 7. What this does and does not say

- The model's scores for family members are far from the cutoff (lowest 41.7 in the independent set, 86.9 in the build). The cutoff is set by the false-hit side, and that side is thin: the next A. fumigatus hit is 0.1 bit below.
- These are scores of known or PANTHER-called family members against false hits from four proteomes. They say nothing about whether a protein is an adhesin, and nothing about proteins outside fungi.
- **The cutoff is NOT changed by this report.** It was set from false hits only under the signed freeze `analysis/wall_fusion_truth/pr6_ecm33_freeze.json`. Moving it, even by 0.1 bit, would be a new rule applied to data already seen. A real recalibration needs a new held-out set (see `2026-10-24-ecm33-curated-ortholog-count.md`: no curated independent set exists yet) and a new signed freeze.
