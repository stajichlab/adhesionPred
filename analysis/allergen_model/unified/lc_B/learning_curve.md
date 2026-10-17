# Learning curve, M1 (ESM-2 650M mean embedding + logistic regression)

- Training negatives: ['family_negative', 'weak_negative']; evaluation negatives: ['family_negative', 'weak_negative']
- Run date: 2026-10-10
- Training clusters available: 272 (676 proteins: 140 positives, 536 negatives)
- Validation: 158 proteins in 48 clusters; test: 82 proteins in 47 clusters
- Draws per fraction: 5 (100% has one, the full training set)

| Train fraction | Draws | Train clusters (mean) | Train positives (mean) | Validation AUROC (mean, sd) | Test AUROC (mean, sd) |
|---|---|---|---|---|---|
| 25% | 5 | 68.0 | 33.8 | 0.711 (0.117) | 0.695 (0.026) |
| 50% | 5 | 136.0 | 61.0 | 0.692 (0.058) | 0.655 (0.019) |
| 75% | 5 | 204.0 | 101.4 | 0.678 (0.027) | 0.685 (0.022) |
| 100% | 1 | 272.0 | 140.0 | 0.675 (0.000) | 0.675 (0.000) |

## Full training set, with 95% cluster-bootstrap intervals

- Test AUROC 0.675, 95% CI 0.511 to 0.833, from 47 test clusters (82 proteins)
- Validation AUROC 0.675, 95% CI 0.583 to 0.884, from 48 validation clusters (158 proteins)

No result here is a specificity claim. The split is below the 500-negative planning floor (spec section 3.2).
