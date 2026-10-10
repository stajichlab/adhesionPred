# Discovery pass on the secretome: Onygenales and Aspergillaceae against a Pezizomycotina background

*2026-10-10. Claude Code (claude-sonnet-5-5). Spec: `docs/superpowers/specs/2026-10-09-secretome-discovery-pass-design.md` (rev 3, passed review 2). Tables: `docs/reports/data/sorting_hat/secretome_discovery/`. This report describes associations between clades and Pfam families. It makes no call, measures no sensitivity or specificity, and does not show function.*

## What was done

- **Genera:** one proteome per named genus (`analysis/fungi_scan/selection_discovery.tsv`): 27 Onygenales, 7 Aspergillaceae, 619 background (Pezizomycotina outside the two targets, including 7 other Eurotiales genera). All genera passed the selection rule after one fix (16 files list 1 protein in `protein_counts.tsv` but hold thousands).
- **Pools:** per proteome, the sorting hat was re-run with the current 18 active families (all 653 proteomes list the same set). Pool S: signal peptide (R0), no mature TM helix, no C-terminal HDEL or KDEL. Pool U: pool S and no mechanism evidence (`other_surface_no_mechanism[R0]`). SignalP 6 and TMHMM came from Fungi_5k (identical to the scan run on 56 proteomes, `2026-10-09-discovery-b2a-fungi5k-agreement.md`). Pfam-A (release 37.x) and PredGPI came from Fungi_5k.
- **Tests:** per Pfam family, genus presence in pool U (and S): two-sided Fisher test, a logistic regression with log pool size, Benjamini-Hochberg over the families tested, minimum 3 target genera. Gains and losses on the Fungi_5k species tree, rooted on Basidiomycota and pruned to the chosen genera. Stability: 5 random subsamples of 150 background genera, stratified by class.

## Numbers

| Target | Pool | Families tested | Enriched, q < 0.05 | In all 5 subsamples | Depleted, q < 0.05 |
|---|---|---|---|---|---|
| Onygenales (27 genera) | U | 3,549 | 49 | 39 | 104 |
| Onygenales | S | 3,577 | 50 | 40 | 106 |
| Aspergillaceae (7 genera) | U | 3,498 | 2 | 0 | 1 |
| Aspergillaceae | S | 3,525 | 3 | 1 | 2 |

"Enriched" uses the full background. "In all 5 subsamples" counts, among the enriched families, those with q < 0.05 in each of 5 random subsamples of 150 background genera (the `seeds_q_lt_0.05_of_5` column). The logistic regression with pool size gives p < 0.05 for all 49 Onygenales pool U families, so pool size alone does not explain them.

**Aspergillaceae: no pattern can be claimed.** With 7 genera, the test finds 2 families (DUF2990, EMC7 beta-sandwich) and neither is stable across subsamples. The review of the spec predicted this: with 7 genera a family must be in nearly all of them and absent from almost all of the background. This report does not say Aspergillaceae lacks patterns. It says this design cannot see them. The within-*Aspergillus* supplement was not run.

## Onygenales families (pool U, q < 0.05, 49 families)

The full table is `families_U_onygenales.tsv` (with Pfam names, genus counts, odds ratio, q, median copies per 1,000 pool proteins, gains and losses on the tree, seeds). The role column holds "no stated role" for every family (see Limits). The descriptions below are the Pfam descriptions. They are not evidence of the protein function in these fungi.

Strongest by q (all stable in 5 of 5 subsamples):

| Family | Pfam description | Onygenales genera / 27 | Background genera / 619 | Gains in target, elsewhere; losses |
|---|---|---|---|---|
| PF01216 | Calsequestrin | 23 | 33 | 1, 23; 8 |
| PF13848 | Thioredoxin-like domain | 26 | 93 | 0, 65; 14 |
| PF01244 | Membrane dipeptidase (M19) | 17 | 21 | 2, 21; 1 |
| PF20138 | DUF6528 | 10 | 2 | 4, 2; 0 |
| PF14862 | Big defensin | 7 | 0 | 4, 0; 0 |
| PF10341 | TPP1/ACD (shelterin) | 3 | 0 | 1, 0; 2 |

Families that are in the three areas the owner named, by description only (a T4 reading, not a role assignment):

- **Cell wall and surface:** PF10287 and PF10290 (Tos1, cell wall protein, C-terminal and N-terminal), PF12454 (Ecm33, GPI-anchored cell wall organization protein), PF04681 (Bys1, *Blastomyces* yeast-phase-specific protein; present in 27 of 27 Onygenales and 448 of 619 background genera, median 6.1 against 2.0 copies per 1,000 pool proteins). Bys1 is an inactive family in the tool.
- **Glycan binding:** PF07938 (fungal fucose-specific lectin), PF18647 (putative glycan binding domain), PF21691 (alpha-galactosyl binding lectin). Present in 6, 15 and 6 Onygenales genera.
- **Secreted degradative enzymes:** M35 deuterolysin and aspzincin (PF02102, PF14521), trypsin (PF00089), dipeptidyl peptidase IV (PF00930), peptidase S41 (PF03572), secretory lipase (PF03583), glycosyl hydrolase family 65 (PF03632, PF03636, PF03633), phospholipase A2 and lysozyme-like domains (PF09056, PF00959). For 29 of the 49 reported families the Onygenales median copy number per 1,000 pool proteins is higher than the background median (the others have a median of 0 in both groups).
- **Unknown function:** DUF6528 (PF20138), DUF4360 (PF14273), DUF2990 (PF11693). DUF6528 is in 10 Onygenales genera and 2 background genera.
- **Not surface-related by description (reported because they pass the test):** calsequestrin, Apc15p (PF05841), TPP1/ACD, EMC7, thioredoxin-like, GATase. Pool U contains proteins that SignalP calls but that may be mis-assigned or intracellular.

## What the tree columns say, and do not say

Gains are counted by a Fitch reconstruction in which ties are resolved to the parent state and a root tie resolves to absent. This favours extra gains for families that are widespread (PF13848 shows 65 gains elsewhere and 14 losses). Read the column as a rough comparison, not as a count of evolutionary events. The one clear contrast is families that are absent from the background and found in several Onygenales genera with gains inside the target clade: Big defensin (PF14862, 7 genera, 4 gains in the target clade) and DUF6528 (4 gains in the target clade). These are candidates for follow-up. This pass does not show what they do.

## Limits

1. **Phylogeny.** Counting genera reduces over-sampling of genera. It does not remove shared ancestry within the Onygenales clade. Fisher p and q are a ranking statistic under an independence assumption that does not hold.
2. **Role table not blind.** The spec required the role table to be written before the statistics were attached. I wrote no roles, so the column says "no stated role" for every family. The descriptions above were read after the results and are labelled as such.
3. **Not run:** the Mann-Whitney test is computed per family but not corrected (BH) in this pass; the Pfam-free cluster analysis; the pathogen contrast inside Onygenales; the within-*Aspergillus* supplement; the second proteome per genus check. The hybrid exclusion list has no literature sources yet.
4. **Pool U depends on the 18 active families.** A family active in one run and inactive in another changes the pool. All 653 runs used the same set.
5. **Pfam-A release.** Full Pfam-A here is release 37.x; the 24 family models are 38.2.
6. **SignalP 6 does not separate wall proteins from secreted enzymes.** `other_surface_no_mechanism` is a residual, not a negative set.
7. **Depletion.** 104 families are depleted in Onygenales pool U at q < 0.05. The Onygenales pool is about half the size of the background pool, but the logistic regression includes pool size and I did not check the depleted list for pool-size effects.
8. **The background contains Eurotiales genera** (7). Results for Onygenales do not change the sister comparison; the sister comparison (descriptive) is not in this report.

## Suggested next steps (owner decides)

1. Add the Mann-Whitney BH correction, the Pfam-free cluster test and the pathogen contrast; run the *Aspergillus* supplement.
2. Look at the candidate proteins of the families in the cell wall, glycan-binding and Big defensin groups (`candidates_U_onygenales.tsv.gz`, 1,053 protein rows) and, for those with literature, add evidence rows.
3. Decide whether to add Bys1 (already inactive in the family table), Tos1 and Ecm33 to the broad-tier families after the family tests of workstream C.
