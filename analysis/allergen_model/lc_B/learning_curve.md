# Learning curve, M1 (ESM-2 650M mean embedding + logistic regression)

- Training negatives: ['family_negative', 'weak_negative']; evaluation negatives: ['family_negative', 'weak_negative']
- Run date: 2026-10-10
- Training clusters available: 277 (604 proteins: 144 positives, 460 negatives)
- Validation: 169 proteins in 69 clusters; test: 143 proteins in 69 clusters
- Draws per fraction: 5 (100% has one, the full training set)

| Train fraction | Draws | Train clusters (mean) | Train positives (mean) | Validation AUROC (mean, sd) | Test AUROC (mean, sd) |
|---|---|---|---|---|---|
| 25% | 5 | 69.0 | 45.6 | 0.654 (0.011) | 0.678 (0.046) |
| 50% | 5 | 138.0 | 66.4 | 0.663 (0.036) | 0.722 (0.039) |
| 75% | 5 | 208.0 | 111.6 | 0.694 (0.020) | 0.721 (0.034) |
| 100% | 1 | 277.0 | 144.0 | 0.691 (0.000) | 0.706 (0.000) |

## Full training set, with 95% cluster-bootstrap intervals

- Test AUROC 0.706, 95% CI 0.547 to 0.837, from 69 test clusters (143 proteins)
- Validation AUROC 0.691, 95% CI 0.553 to 0.819, from 69 validation clusters (169 proteins)

No result here is a specificity claim. The split is below the 500-negative planning floor (spec section 3.2).
