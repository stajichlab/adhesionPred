# Allergen predictor: status report against the plan (2026-10-10)

Research use only. These are AUROC values from a logistic regression on ESM-2 650M mean embeddings (M1). No result here is a
specificity claim. Every number has its cluster count. Intervals are 95% cluster bootstrap (1,000 resamples of clusters).

## 1. Data built

| Set | Positives | Negatives | Source and notes |
|---|---|---|---|
| Fungal | 200 | 200 length-matched weak negatives, 516 family-matched negatives | Positives: WHO/IUIS, KW-0020, AllergenOnline (fungal and mapped), COMPARE (4 mapped). Negatives: reviewed Swiss-Prot fungi, no allergen annotation. 26 positives from NCBI protein IDs, not UniProt. |
| All-kingdom | 1,679 | 1,669 length-matched (plant 904, animal 765), 1,547 family-matched (plant 654, animal 893) | Positives: AllergenOnline non-fungal rows mapped to UniProt. Negatives: reviewed Swiss-Prot non-allergens from Viridiplantae and Metazoa (human excluded). |

Filters applied to all negatives: no KW-0020, no allergen comment, no "antigen" word, no "Cancer", length 10 to 850,
protein-level evidence, annotation score 4 or 5. No negative shares an accession with a positive or with the other set.

## 2. Pipeline status (plan phases)

| Phase | Status |
|---|---|
| T0.2 positives | Done (fungal 200 proteins in positives.tsv; all-kingdom 1,679 separate) |
| T0.3 candidates | Done (6,615 fungal candidates, 0 overlap with positives) |
| T0.4 pilot | Done (job 29797481) |
| T0.5 pairwise search | Done for the fungal set; pairing dropped 2026-10-10 |
| B6 positive-to-positive search | Done (job 29799979) |
| Length-matched split | Done (fungal, 82 positive clusters; all-kingdom, 1,070 clusters in the shared set) |
| B6 shortcut check | Done, four seeds (see section 3) |
| Embeddings | Done for fungal (916 sequences, job 29807328) and all-kingdom (4,895, job 29821859) |
| Learning curves | Done (sections 3 and 4) |
| Within-kingdom AUROC | Done (section 5) |
| M1b random forest, M1L length-matched run, permutation and seed controls | **Not run** |
| Threshold rule, module integration (T5) | **Not run** |
| COMPARE full download | **Waiting on owner** (form needs an email address) |

## 3. Fungal set

**B6 shortcut check** (nearest-neighbour identity to training positives; mean over four split seeds, stop rule O1 as amended):
- Validation AUROC 0.577 (sd 0.044), 4 seeds. Test AUROC 0.530 (sd 0.042), 4 seeds.
- Below the 0.60 stop threshold. No halt.

**Learning curve, length-matched negatives only** (condition A; 175 training clusters, validation 41 clusters and 62 proteins, test 48 clusters and 62 proteins):

| Train fraction | Train clusters | Validation AUROC (sd) | Test AUROC (sd) |
|---|---|---|---|
| 25% | 44 | 0.812 (0.013) | 0.842 (0.046) |
| 50% | 88 | 0.807 (0.005) | 0.847 (0.011) |
| 75% | 131 | 0.827 (0.014) | 0.856 (0.015) |
| 100% | 175 | 0.814 | 0.860 |

Full training set, 95% CI: test 0.860 (0.735 to 0.954), 48 clusters, 62 proteins. Validation 0.814 (0.679 to 0.932), 41 clusters, 62 proteins.

**Learning curve, family-matched negatives added** (condition B, shared clustering; 272 training clusters, validation 48 clusters and 158 proteins, test 47 clusters and 82 proteins):

| Train fraction | Train clusters | Validation AUROC (sd) | Test AUROC (sd) |
|---|---|---|---|
| 25% | 68 | 0.711 (0.117) | 0.695 (0.026) |
| 50% | 136 | 0.692 (0.058) | 0.655 (0.019) |
| 75% | 204 | 0.678 (0.027) | 0.685 (0.022) |
| 100% | 272 | 0.675 | 0.675 |

Full training set, 95% CI: test 0.675 (0.511 to 0.833), 47 clusters, 82 proteins (14 positives). Validation 0.675 (0.583 to 0.884), 48 clusters, 158 proteins.

Reading: adding family-matched negatives lowers test AUROC from about 0.86 to about 0.68. The earlier figure was partly a family and function effect. The test set has only 14 positives, so the interval is very wide.

## 4. All-kingdom set

**Learning curve** (family and kingdom negatives in training and evaluation; 1,070 training clusters, validation 254 clusters and 837 proteins, test 220 clusters and 694 proteins):

| Train fraction | Train clusters | Validation AUROC (sd) | Test AUROC (sd) |
|---|---|---|---|
| 25% | 268 | 0.799 (0.016) | 0.851 (0.045) |
| 50% | 535 | 0.791 (0.017) | 0.864 (0.021) |
| 75% | 802 | 0.803 (0.005) | 0.871 (0.033) |
| 100% | 1,070 | 0.805 | 0.886 |

Full training set, 95% CI: test 0.886 (0.788 to 0.947), 220 clusters, 694 proteins. Validation 0.805 (0.743 to 0.863), 254 clusters, 837 proteins.

## 5. Within-kingdom AUROC (all-kingdom model, full training set)

Each kingdom is scored against its own negatives only, so a kingdom difference cannot inflate the result.

| Split | Group | Positives | Negatives | Clusters | AUROC | 95% CI |
|---|---|---|---|---|---|---|
| Validation | Plant | 94 | 323 | 133 | 0.805 | 0.733 to 0.860 |
| Validation | Animal | 111 | 309 | 140 | 0.810 | 0.725 to 0.889 |
| Test | Plant | 95 | 226 | 113 | 0.859 | 0.754 to 0.940 |
| Test | Animal | 151 | 222 | 117 | 0.906 | 0.769 to 0.964 |

Reading: the within-kingdom AUROC is close to the pooled value, so the pooled result is not only a plant-versus-animal effect.

## 6. What the results can and cannot support

- The all-kingdom set supports an AUROC claim: about 0.80 on validation and 0.89 on test, with intervals that exclude 0.5.
- The fungal set is too small for a claim. The test set has 14 positives.
- No set has a specificity claim. Every split is below the 500-negative planning floor, and no evidence-grade negative exists (IEDB had none usable).
- Negatives are weak: they come from Swiss-Prot without allergen annotation, and some may be unannotated allergens.

## 7. Next steps (plan items still open)

1. M1b random forest and M1L length-matched run on the all-kingdom set.
2. Five seeds and one label-permutation run for every reported model.
3. Threshold rule (plan O2: validation sensitivity 0.80), then the evidence field for the sorting hat (T5).
4. COMPARE download, when the owner provides an email address for the form.
5. Report to the owner for review. No merges or releases without approval.
