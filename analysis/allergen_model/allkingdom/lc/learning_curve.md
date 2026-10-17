# Learning curve, M1 (ESM-2 650M mean embedding + logistic regression)

- Training negatives: ['family_negative', 'kingdom_negative']; evaluation negatives: ['family_negative', 'kingdom_negative']
- Run date: 2026-10-10
- Training clusters available: 1070 (3364 proteins: 1228 positives, 2136 negatives)
- Validation: 837 proteins in 254 clusters; test: 694 proteins in 220 clusters
- Draws per fraction: 5 (100% has one, the full training set)

| Train fraction | Draws | Train clusters (mean) | Train positives (mean) | Validation AUROC (mean, sd) | Test AUROC (mean, sd) |
|---|---|---|---|---|---|
| 25% | 5 | 268.0 | 278.6 | 0.799 (0.016) | 0.851 (0.045) |
| 50% | 5 | 535.0 | 561.6 | 0.791 (0.017) | 0.864 (0.021) |
| 75% | 5 | 802.0 | 938.2 | 0.803 (0.005) | 0.871 (0.033) |
| 100% | 1 | 1070.0 | 1228.0 | 0.805 (0.000) | 0.886 (0.000) |

## Full training set, with 95% cluster-bootstrap intervals

- Test AUROC 0.886, 95% CI 0.788 to 0.947, from 220 test clusters (694 proteins)
- Validation AUROC 0.805, 95% CI 0.743 to 0.863, from 254 validation clusters (837 proteins)

No result here is a specificity claim. The split is below the 500-negative planning floor (spec section 3.2).
