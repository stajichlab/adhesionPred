# Learning curve, M1 (ESM-2 650M mean embedding + logistic regression)

- Training negatives: ['weak_negative']; evaluation negatives: ['weak_negative']
- Run date: 2026-10-10
- Training clusters available: 175 (276 proteins: 144 positives, 132 negatives)
- Validation: 62 proteins in 41 clusters; test: 62 proteins in 48 clusters
- Draws per fraction: 5 (100% has one, the full training set)

| Train fraction | Draws | Train clusters (mean) | Train positives (mean) | Validation AUROC (mean, sd) | Test AUROC (mean, sd) |
|---|---|---|---|---|---|
| 25% | 5 | 44.0 | 32.4 | 0.812 (0.013) | 0.842 (0.046) |
| 50% | 5 | 88.0 | 81.8 | 0.807 (0.005) | 0.847 (0.011) |
| 75% | 5 | 131.0 | 103.2 | 0.827 (0.014) | 0.856 (0.015) |
| 100% | 1 | 175.0 | 144.0 | 0.814 (0.000) | 0.860 (0.000) |

## Full training set, with 95% cluster-bootstrap intervals

- Test AUROC 0.860, 95% CI 0.735 to 0.954, from 48 test clusters (62 proteins)
- Validation AUROC 0.814, 95% CI 0.679 to 0.932, from 41 validation clusters (62 proteins)

No result here is a specificity claim. The split is below the 500-negative planning floor (spec section 3.2).
