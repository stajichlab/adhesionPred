# AGA1 broad flag candidate: signal peptide + GPI + AIUPred disordered run

Date: 2026-10-24. Owner decisions: AIUPred alone is the disorder tool; minimum disordered run 200 residues, kept as a parameter
(`analysis/aga1_class/broad_flag.py --min-disordered-run`, default 200), 100 as the pre-declared sensitivity value.
Tables: `docs/reports/data/aga1_broad_flag_2026-10-24/`. The Ser/Thr stalk rule is not part of this candidate and is still open.

## Counts (flag = yes)

| Group | n | run >= 200 | run >= 100 |
|---|---|---|---|
| tier A | 35 | 23 | 26 |
| tier B | 32 | 18 | 23 |
| hard negatives | 172 | 8 | 23 |
| background assumed negatives | 389 | 7 | 31 |
| other-mechanism adhesins | 49 | 6 | 12 |
| Kiss rows (other clusters) | 70 | 1 | 1 |
| excluded indirect regulators | 60 | 0 | 0 |
| excluded surface/other phenotype | 36 | 0 | 5 |
| uncurated repeat-bearing | 7 | 1 | 3 |

No protein is NA at either cut.

## Reading the misses and the calls (run >= 200)

* **Reference proteins:** AGA1 yes (run 207), DAN4 yes (444), SED1 **no** (run 151, GPI yes), MSB2 no (no GPI call, run 406).
  SED1 is a broad positive by the nested labels, so 200 misses one of the proteins the flag is meant to hold; at 100 it is yes.
* **Tier A misses (12 of 35):** three have GPI yes and a run just under the cut (ALS4498 178, FLO10 175, MAD2 141); eight have no
  PredGPI call (HWP1 and seven fission-yeast map4/Pfl proteins; pfl5 has a run of 518); SAG1 has a run of 11.
  So the three-part rule loses mainly proteins PredGPI does not call, not proteins without disorder.
* **Flagged "hard negatives" (8):** DAN4, TIR4, FIT1, HKR1 (S. cerevisiae) and cfmA, chiA, chiA1, xgeA (Aspergillus). DAN4 and TIR4
  are Ser/Thr-rich wall mannoproteins, which the nested labels count as broad positives, so they are correct calls for the broad
  flag; their "hard negative" label belongs to the narrow flag. The Aspergillus chitinase and glucanase entries are enzyme wall
  proteins with a disordered linker, which the stalk rule or a domain check would have to remove.

## What this is not

An evaluated flag. The cut was chosen from the tier A table before any held-out clade was looked at, but the same proteins
define tier A and the counts above, so they are optimistic. Nothing here tests the flag leaving a clade out.

## Next decisions (owner)

1. Ser/Thr stalk rule (add a composition or window rule, or leave the broad flag as these three parts).
2. Clade unit for narrow classes, glycosylation predictor, ESM-2 tiling for the 74 proteins over 1,022 aa, extra UniProt positives.
