# Test of the 14 inactive Pfam families (and the 10 active ones for reference)

*2026-10-09. Claude Code (claude-sonnet-5-5). Method, inputs and limits: `analysis/family_tests/README.md`. Table: `docs/reports/data/sorting_hat/family_tests/family_tests.tsv`. This is a descriptive test. It gives counts, not estimated sensitivity or specificity. E2 positives are family-defined and are not counted as evidence of sensitivity.*

## What the numbers are

- `pos_E1_hit`: E1 adhesins (direct experimental evidence; 54 with sequence) that hit the family at `--cut_ga`.
- `neg_N1_hit`, `neg_N1_clusters_hit`: N1 hard negatives (137 rows, 76 clusters) that hit.
- `scan_*`: distinct proteins with a raw `--cut_ga` hit in the 56 scan proteomes. This is higher than the counts in the scan report for active families, because that report counts proteins after the module rules (for example, PA14 needs a signal peptide).

## Results for the 14 inactive families

| Family | E1 hit (of 54) | N1 hit (clusters) | N2 hit | Scan proteins (proteomes) | With signal peptide | With repeat call |
|---|---|---|---|---|---|---|
| Bys1 PF04681 | 1 | 0 | 0 | 138 (43) | 128 | 1 |
| Candida_ALS_N PF11766 | 10 | 0 | 0 | 30 (6) | 25 | 13 |
| Candida_ALS PF05792 | 9 | 0 | 0 | 31 (5) | 24 | 15 |
| AltA1 PF16541 | 0 | 0 | 0 | 9 (4) | 9 | 0 |
| Allergen_Asp_f_4 PF25312 | 0 | 0 | 0 | 68 (30) | 57 | 0 |
| Flo11 PF10182 | 6 | 0 | 0 | 9 (5) | 9 | 3 |
| GLEYA PF10528 | 12 | 0 | 1 | 42 (9) | 31 | 6 |
| Hyphal_reg_CWP PF11765 | 3 | 0 | 0 | 63 (7) | 55 | 6 |
| PIR PF00399 | 0 | 4 (2) | 2 | 17 (5) | 17 | 4 |

The five families the owner turned on on 2026-10-09:

| Family | E1 hit (of 54) | N1 hit (clusters) | N2 hit | Scan proteins (proteomes) | With signal peptide | With repeat call |
|---|---|---|---|---|---|---|
| Flocculin PF00624 | 3 | 0 | 0 | 26 (2) | 14 | 8 |
| Flocculin_t3 PF13928 | 8 | 0 | 0 | 54 (6) | 44 | 7 |
| Hyr1 PF15789 | 0 | 0 | 0 | 15 (3) | 15 | 5 |
| PIR1-like_C PF22799 | 0 | 3 (1) | 1 | 57 (25) | 57 | 7 |
| ALS_M PF30910 | 0 | 0 | 0 | 38 (36) | 24 | 0 |

## Findings (facts first, then what each suggests)

1. **ALS_M (PF30910) does not hit Candida Als proteins.** The 14 E1 adhesins with a Candida_ALS_N hit have no ALS_M hit at `--cut_ga` or at E 1e-3. In the scan it hits 38 proteins, one per proteome in 36 proteomes (24 with a signal peptide, none with a repeat call), and only 5 proteomes are Candida-type. It is active as of today. It does not do what the owner expected ("detect Als"). The owner may want to turn it off again until the hits are understood. The Candida_ALS_N (PF11766) and Candida_ALS (PF05792) models do hit the Als positives (10 and 9 of 54 E1; the other E1 positives hit neither).
2. **PIR1-like_C (PF22799) hits 3 N1 negatives (1 cluster) and 1 N2 negative; no E1 positive.** It is active. A hard negative hit means the family fires on a protein with characterised non-adhesive function. PIR (PF00399, inactive) hits 4 N1 negatives (2 clusters) and 2 N2. Pir proteins are structural wall proteins; the data agree with that, not with an adhesion role.
3. **Hyr1 (PF15789) has no E1 hit.** Its 9 E2 hits come from labels that were made from the family. It hits 15 scan proteins in 3 proteomes. No test of its specificity is possible with the current positives.
4. **Flocculin (PF00624) and Flocculin_t3 (PF13928):** 3 and 8 E1 hits. No negative hit. Scan hits are concentrated: PF00624 in 2 proteomes (26 proteins), PF13928 in 6 (54).
5. **Candida_ALS_N, Candida_ALS, Flo11, GLEYA:** E1 hits 10, 9, 6, 12 of 54, no N1 hit (GLEYA hits 1 N2). They are inactive. Candida_ALS_N and Candida_ALS together hit 13 and 15 proteins with a repeat call in 5 to 6 proteomes.
6. **Bys1 (PF04681):** 1 E1 hit (CalA) and 138 scan proteins in 43 proteomes, 128 with a signal peptide. It is common. It also covers CalA's paralogs (README of the curated table). Its adhesion role in the other proteins is not shown.
7. **Allergen families (AltA1, Allergen_Asp_f_4):** no adhesion truth applies. AltA1 hits 9 proteins in 4 proteomes; Allergen_Asp_f_4 hits 68 in 30. Their test would need the IUIS allergen list (`analysis/allergen_scoping/`); I did not run it.
8. **CFEM (PF05730) and PA14 (PF07691), active since 2026-10-07/08, for comparison:** CFEM hits 12 N1 negatives in 9 clusters (of 76) and 5 N2, with 1 of 54 E1 positives. PA14 hits 4 N2 negatives and 9 of 54 E1 positives. CFEM feeds `cell_wall_adhesion_candidate` (owner decision M2 O1 (a), 2026-10-09). On these data CFEM fires on many characterised non-adhesive proteins. The owner decided to keep it with the family shown. These counts are new evidence for that decision.

## Suggested actions (owner decides)

| Family | Suggestion | Based on |
|---|---|---|
| PF30910 ALS_M | Turn off until its hits are examined | Finding 1 |
| PF22799 PIR1-like_C | Keep in the broader call only, or review | Finding 2 |
| PF15789 Hyr1 | Keep on; needs E1 positives (Hyr1/Iff proteins from the literature) to test | Finding 3 |
| PF00624, PF13928 | Keep on | Findings 4 |
| PF11766, PF05792, PF10182, PF10528 | Candidates to turn on next; check the 4 repeat overlaps and the single N2 hit of GLEYA | Finding 5 |
| PF04681 Bys1 | Leave off pending a curated Bys1 positive set | Finding 6 |
| Allergen families | Test against IUIS before any decision | Finding 7 |
| PF05730 CFEM | Revisit M2 O1 with the N1 hit rate in hand | Finding 8 |

## Limits

1. The hard-negative set and the positives were curated by an agent and have not been owner-reviewed.
2. Positives are not clustered. E1 hit counts are protein counts and may count one family several times (for example the FLO proteins).
3. N1 clusters hit are cluster counts within the hard-negative set (76 N1 clusters); they are not rates over all proteins.
4. The scan proteomes have no truth, so the scan columns describe how often a family fires, not whether it is right.
