# Hydrophobin call validation: first measurement

*2026-10-08. Branch `hydrophobin-validation` (local, not pushed). Spec: `docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md`.
Plan: `docs/superpowers/plans/2026-10-08-hydrophobin-validation.md`. Data: `docs/reports/data/sorting_hat/hydrophobin/` and
`docs/reports/data/sorting_hat/call_status/hydrophobin/`. Every number below comes from those files or from the commands named.*

## 1. Summary

1. `pfam_hydrophobin` (the call `hydrophobin_domain`) has a measured status in 8 proteomes. All 8 are `smoke`. None can be `estimated`: each proteome has 3 to 7 labelled positives in 1 to 5 clusters.
2. Sensitivity is 1.00 in four proteomes (*A. fumigatus* x3, *B. bassiana*) and *P. ostreatus* PC9, 0.80 in *P. expansum*, 0.60 in *F. graminearum*, and 0.25 in *F. fulva*. The intervals are wide (lower bounds 0.00 to 0.57).
3. Specificity is 0.9994 to 1.0000. This number is not informative. The negatives are assumed (no annotation), and 19 of the 22 called proteins without a label have 8 or more cysteines (section 5). Many are probably unlabelled hydrophobins.
4. The 8-cysteine rescue rule, with the frozen published spacing, finds one extra labelled hydrophobin in these proteomes (HYD2_GIBZE, *F. graminearum*). It makes 3 more calls with no label. Rule 6.3 gives "no evidence that the rescue helps". `hydrophobin_protein` was not added to `categories.yaml`.
5. The frozen spacing matches only 78 of 131 T2 Swiss-Prot hydrophobins (60%): 75 of 122 with a Pfam cross-reference and 3 of 9 without. As a recall rule, the published spacing is too narrow.

## 2. What was run

- Truth: UniProt release 2026_03, 131 T2 and 43 T3 entries, mapped to the proteomes by BLASTP (at least 95% identity and 90% coverage). Positives are T1/T2 only. T3 mapped proteins are excluded. See `analysis/hydrophobin_truth/README.md`.
- Negatives: every other protein of the proteome. They are **assumed** negatives.
- Calls: `hydrophobin_domain` (module `pfam_hydrophobin`, Pfam 38.2 GA, seven hydrophobin-class models) from the core run `out_hyd`; `cys8_pattern` from `cellsurface_sorting_hat_module cys8` with the spacing file `cys8_spacing.yaml` (14 sets, any one set matches, R0 signal peptide required).
- Status files: `cellsurface_sorting_hat_calibrate truth --module pfam_hydrophobin --leakage partial`, whole truth table (no usable test split, section 6).

## 3. Status of `pfam_hydrophobin` (module status, call `hydrophobin_domain`)

| Proteome | Status | Positives / clusters | Sensitivity [95% CI] | Specificity [95% CI] | Assumed negatives |
|---|---|---|---|---|---|
| *A. fumigatus* Af293 | smoke | 7 / 5 | 1.00 [0.57, 1.00] | 1.0000 [0.9995, 1.0000] | 9640 |
| *A. fumigatus* A1163 | smoke | 6 / 4 | 1.00 [0.51, 1.00] | 1.0000 [0.9995, 1.0000] | 9936 |
| *A. fumigatus* W72310 | smoke | 5 / 3 | 1.00 [0.44, 1.00] | 0.9997 [0.9990, 1.0000] | 10551 |
| *B. bassiana* ARSEF 2860 | smoke | 7 / 4 | 1.00 [0.51, 1.00] | 0.9999 [0.9993, 1.0000] | 9471 |
| *F. fulva* Race5 | smoke | 4 / 3 | 0.25 [0.00, 1.00] | 0.9996 [0.9991, 0.9999] | 13556 |
| *F. graminearum* PH-1 | smoke | 5 / 4 | 0.60 [0.00, 1.00] | 1.0000 [0.9996, 1.0000] | 11188 |
| *P. expansum* MD-8 | smoke | 5 / 4 | 0.80 [0.00, 1.00] | 0.9998 [0.9992, 1.0000] | 10619 |
| *P. ostreatus* PC9 | smoke | 3 / 1 | 1.00 [0.21, 1.00] | 0.9994 [0.9986, 0.9999] | 10479 |

Intervals are the tool's own (cluster bootstrap for sensitivity). With 1 to 5 clusters they are not reliable. Leakage `partial`: the author had seen the 9 no-Pfam entries and run a screening regex on the proteomes.
Proteomes with no truth positives (S288C, *C. albicans*, *C. immitis* RS, *B. dermatitidis* ER3) have no status. Hydrophobin_domain called 0, 0, 2 and 2 proteins there.

## 4. Variants and the rescue (from `summary.tsv`, `review_lists.tsv`)

| Proteome | Pfam only: TP/FN | Pfam or cys8: TP/FN | cys8 only: TP/FN | Rescue-only calls |
|---|---|---|---|---|
| Afum Af293 / A1163 / W72310 | 7/0, 6/0, 5/0 | same | 4/3, 4/2, 4/1 | 0 |
| *B. bassiana* | 7/0 | 7/0 | 3/4 | 1 (unlabelled, 143 aa, 24 Cys, no Swiss-Prot hit at 1e-3) |
| *F. fulva* | 1/3 | 1/3 | 0/4 | 0 |
| *F. graminearum* | 3/2 | 4/1 | 1/4 | 1 (labelled: F0349401_001781 = HYD2_GIBZE, T2) |
| *P. expansum* | 4/1 | 4/1 | 2/3 | 0 |
| *P. ostreatus* PC9 | 3/0 | 3/0 | 3/0 | 0 |
| *C. immitis* RS, *B. dermatitidis* | no truth | | | 1 + 1 (unlabelled) |

- The F. fulva proteins HCF1, HCF2 and HCF4 (T2, no Pfam model) are not matched by the frozen spacing. The screening regex of the discovery search had matched all nine no-model entries. The published sets do not.
- Rescue-only review (Swiss-Prot 2023_03 and the 174 known hydrophobins, BLASTP): F0349401_001781 is identical to HYD2_GIBZE (hydrophobin). *B. dermatitidis* F00FD2C2_006465 (101 aa, 8 Cys) has weak similarity to a hydrophobin (38.9% over 54 aa, E 5e-10) and 49% identity to the uncharacterised secreted protein ARB_06108. *C. immitis* XP_001245267.2 (299 aa, 16 Cys) and *B. bassiana* F1BB8A46_008916 have no hit to a known hydrophobin. Outcomes: 1 hydrophobin, 3 unresolved.
- **Rule 6.3**: 4 rescue-only calls in 4 clusters, 1 resolved as hydrophobin. The rule needs at least 5 resolved clusters. Result: "no evidence that the rescue helps". `cys8_pattern` stays in the tool and is reported. It is not part of a default call.

## 5. Calls without a label

22 proteins have a Pfam hydrophobin call and no T1/T2 label (T3 mapped proteins are excluded from this count). 19 of them have 8 or more cysteines and most are 85 to 170 aa. These are probably unlabelled hydrophobins, not false positives: Swiss-Prot has no entry for them. Examples:
*A. fumigatus* W72310 KAK9636619.1 and KAK9656839.1 (the Hydrophobin_D and DewD hits found earlier), *P. ostreatus* PC9 7 proteins (the species has 40 published hydrophobin genes, PMID 33636611, of which Swiss-Prot has 28 entries).
They are **not** reviewed to label level here. The table in `review_lists.tsv` (`pfam_call_not_positive`) lists each one. Specificity is therefore a lower bound of unknown size, and the text does not claim a specificity.

## 6. Sensitivity of the frozen spacing on Swiss-Prot hydrophobins (`pattern_on_swissprot.tsv`)

| Group | n | Pattern match |
|---|---|---|
| T2 with a Pfam hydrophobin cross-reference | 122 | 75 (61%) |
| T2 without a Pfam cross-reference | 9 | 3 (33%) |
| T3 with a cross-reference | 43 | 32 (74%) |

Class assigned by the matching set: I 82, II 24, other 4 (a protein can match several sets). Not matched, for example: RODL_NEUCR (EAS), MHP1_PYRO7, HYD1C to HYD1F_BEAB2, HYPD_AGABI, HFBF/HFBG_PENEN, SSP1_OIDMZ. This is conformity of known hydrophobins to the published spacing, not discovery (labels were partly named from the pattern).

## 7. Limits

- All statuses are `smoke`, leakage `partial`. The development/test split (30%/70% of clusters, `split.tsv.gz`) is stored, but the test parts hold 2 to 7 positives, so no species has the 10 positive clusters that the spec asks for. The measurement uses the whole table.
- Specificity rests on assumed negatives (section 5).
- The cross-species hard-negative table (spec 4.3 item 2, `hard_negatives_swissprot.tsv`) was not built. Pattern false positives are only seen in whole proteomes.
- No leave-one-species-out table: the rule has no fitted parameters, so a leave-out would give the same numbers as section 6.
- Wessels 1994, Linder 2005 and Sunde 2008 were not read. The spacing sets come from other papers (`cys8_spacing.yaml`).
- `hsba_domain` has no truth set and no status.
- The repeat-call status files (`tandem_repeat_protein`) are stale because `categories.yaml` changed. H1b is not done yet.

## 8. Decision for the owner

The published spacing does not give a useful rescue: it recovers 1 of the 6 labelled proteins that Pfam misses in these proteomes, and 3 of 9 no-model Swiss-Prot entries. Options:
1. Leave the rescue out of the default call (the result of rule 6.3). Default.
2. Build a data-derived spacing from the 122 Pfam-positive T2 entries (not the 9), freeze it, and test it on the 9 no-model entries and on the proteomes. This is a new rule with a new freeze. The 9 entries are not independent of the author, so the result would be `smoke` with leakage `partial` at best.
