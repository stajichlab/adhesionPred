# Discovery pass on the secretome: Onygenales and Aspergillaceae against a Pezizomycotina background

*2026-10-10 (revised later the same day with the Mann-Whitney correction, the Pfam-free cluster test, the pathogen contrast, the second-proteome check and the *Aspergillus* supplement). Claude Code (claude-sonnet-5-5). Spec: `docs/superpowers/specs/2026-10-09-secretome-discovery-pass-design.md` (rev 3, passed review 2). Tables: `docs/reports/data/sorting_hat/secretome_discovery/`. This report describes associations between clades and Pfam families. It makes no call, measures no sensitivity or specificity, and does not show function.*

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

Gains are counted by a Fitch reconstruction in which ties are resolved to the parent state and a root tie resolves to absent. This favours extra gains for families that are widespread (PF13848 shows 65 gains elsewhere and 14 losses). Read the column as a rough comparison, not as a count of evolutionary events. Families that are absent from the background and found in several Onygenales genera with gains inside the target clade are Big defensin (PF14862, 7 genera, 4 gains in the target clade) and DUF6528 (4 gains in the target clade). See the pathogen contrast below before treating them as host-interaction candidates.

## Additional checks (run after the first version of this report)

### Copy number (Mann-Whitney with its own BH correction)

Run for families present in at least 3 target genera (Onygenales pool U: 377 families; Aspergillaceae pool U: 239). Table: `U_onygenales_all.tsv`, `U_aspergillaceae_all.tsv` (column `mwu_q`).

| Target (pool U) | Families with mwu_q < 0.05 | With a higher target median copy number | Also enriched by presence (Fisher q < 0.05) |
|---|---|---|---|
| Onygenales | 204 | 85 | 29 |
| Aspergillaceae | 47 | 39 | 2 |

The copy-number test finds families where presence is the same but copies differ. For Aspergillaceae it is the only test with power: 39 families have a higher median copy number, while presence finds 2. By description they include wall-related families (SUN PF03856, Kre9/KNH PF10342, X8 PF07983, LysM PF01476, Ecm33 PF12454, Bys1 PF04681) and many secreted enzymes (glycosyl hydrolases, lipases, pectin lyase, amylase). Most are present in nearly all 40 *Aspergillus* species of the supplement, so the signal is copy number, not presence. These are associations with the clade. Higher copy number can reflect gene family expansion or a different secretome size; the median is per 1,000 pool proteins.

### Pfam-free clusters (MMseqs2 clusters of Fungi_5k, no Pfam hit in the chosen proteomes)

53,022 clusters were tested for Onygenales and 52,327 for Aspergillaceae (pool U; `clusters_*_enriched.tsv`). At q < 0.05: 168 enriched clusters for Onygenales and 24 for Aspergillaceae. The strongest Onygenales cluster (representative `F4EA5946_003492-T1`, *Uncinocarpus reesii*) is in 24 of 27 Onygenales genera and 14 of 619 background genera (q = 5e-26). The representative protein has 173 aa, 0% Cys, 8.1% Pro, 16.8% Ser+Thr and a PredGPI call. It is a Pfam-free, GPI-predicted, Ser/Thr-rich small protein, so a wall-protein candidate by composition. This was not checked against literature. SOWgp-like proteins are expected to be Pfam-free but are present in *Coccidioides* only, so they fall below the 3-genus minimum.

### Pathogen contrast inside Onygenales (descriptive, no test)

Host groups and sources: `analysis/secretome_discovery/onygenales_host_groups.tsv` (PubMed-verified sources). Strict mammal-pathogen group (7 genera): *Coccidioides*, *Histoplasma*, *Blastomyces*, *Emergomyces*, *Paracoccidioides*, *Trichophyton*, *Microsporum*. Broad group (10) adds *Emmonsia*, *Nannizzia*, *Arthroderma*. The other 13 genera are reptile pathogens (3, one source), *Uncinocarpus* (nonpathogenic relative, Sharpton 2009) or not assessed. Table: `pathogen_contrast_U_onygenales.tsv`.

- **Most of the 49 enriched families are clade-wide, not pathogen-specific.** M35 deuterolysin: 7 of 7 pathogen genera and 19 of 20 others. Bys1, Ecm33, Tos1, glycosyl hydrolase family 65: all 7 of 7 against 18 to 20 of 20. This agrees with Sharpton 2009 (PMID 19717792), who found M35 expanded in *Coccidioides* and *Uncinocarpus* (the nonpathogenic relative).
- **Big defensin (PF14862) is absent from all 7 mammal-pathogen genera and present in 7 of 20 other Onygenales genera.** It is clade-specific but the contrast does not point to host interaction.
- **Pathogen-leaning families** (difference of at least 0.4 in the fraction of genera): Defensin_2 (arthropod defensin, PF01097) in 5 of 7 pathogen genera against 1 of 20 others; SARAF (PF06682) in 4 of 7 against 3 of 20. n is 7 genera, so no claim is made.
- **Lectin families** are not pathogen-specific: PF07938 in 2 of 7 against 4 of 20; PF21691 in 2 of 7 against 4 of 20; PF18647 in 2 of 7 against 13 of 20.

### Second proteome per genus (18 target genera have one; 15 Onygenales and 3 Aspergillaceae)

- Presence concordance between the first and second proteome of the same genus: 13,454 of 15,084 family-genus calls agree (0.892; pool U, 838 families).
- Re-running the test with the second proteome substituted where one exists (mixed set): Onygenales pool U gives 52 enriched families against 49; 42 are in both sets, 7 are lost (PF00878, PF01097, PF03583, PF03633, PF06682, PF07938, PF12192) and 10 are new. Aspergillaceae gives 3 against 2; the two original families remain (`second_U_*.txt`).
- Families lost with the second proteome include the fungal fucose-specific lectin PF07938 and Defensin_2. Treat results for families near the 3-genus minimum as fragile.

### *Aspergillus* supplement (40 species in 8 clades, spaced across the species tree)

Clades are average-linkage clusters of patristic distances cut at 8 groups. Table: `supplement_U_clades.tsv` for 39 distinct families of Aspergillaceae pool U (the families with mwu_q < 0.05 and higher target copy number; this set contains the 2 presence-enriched families and the 2 nominal-p families, `asp_families.txt`): species with the family and median copies per 1,000 pool proteins per clade. 21 of the 39 families are present in at least 38 of the 40 species; the 2 presence families are less uniform (EMC7 PF09430 in 26 of 40 species, DUF2990 PF11693 in 27 of 40). GATase_6 and GATase_7 (PF13522, PF13537) are in 3 of 40 species. No test is run.

## Limits

1. **Phylogeny.** Counting genera reduces over-sampling of genera. It does not remove shared ancestry within the Onygenales clade. Fisher p and q are a ranking statistic under an independence assumption that does not hold.
2. **Role table not blind.** The spec required the role table to be written before the statistics were attached. I wrote no roles, so the column says "no stated role" for every family. The descriptions above were read after the results and are labelled as such.
3. **Done after the first version:** the Mann-Whitney correction, the Pfam-free cluster test, the pathogen contrast, the second-proteome check and the *Aspergillus* supplement (see above). **Still open:** the hybrid exclusion list has no literature sources; the pathogen group for 13 genera is "not assessed"; the Pfam-free clusters have no functional annotation beyond composition.
4. **Pool U depends on the 18 active families.** A family active in one run and inactive in another changes the pool. All 653 runs used the same set.
5. **Pfam-A release.** Full Pfam-A here is release 37.x; the 24 family models are 38.2.
6. **SignalP 6 does not separate wall proteins from secreted enzymes.** `other_surface_no_mechanism` is a residual, not a negative set.
7. **Depletion.** 104 families are depleted in Onygenales pool U at q < 0.05. The Onygenales pool is about half the size of the background pool, but the logistic regression includes pool size and I did not check the depleted list for pool-size effects.
8. **The background contains Eurotiales genera** (7). Results for Onygenales do not change the sister comparison; the sister comparison (descriptive) is not in this report.

## Suggested next steps (owner decides)

1. Read the literature for the strongest candidates and add evidence rows (Bys1, Tos1, Ecm33, the top Pfam-free GPI cluster, deuterolysin M35 in the family tables).
2. Look at the candidate proteins of the families in the cell wall and glycan-binding groups (`candidates_U_onygenales.tsv.gz`, 1,053 protein rows) and, for those with literature, add evidence rows.
3. Decide whether to add Bys1 (already inactive in the family table), Tos1 and Ecm33 to the broad-tier families after the family tests of workstream C.
