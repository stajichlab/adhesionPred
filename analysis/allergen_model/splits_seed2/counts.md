# Cluster split with length-matched weak negatives

- Run date: 2026-10-10; seed 2; length bin 50 aa; 1:1 per bin
- Positive clusters: 82; sampled negatives: 200
- Bins where candidates were fewer than positives: none

| Split | Positives | Weak negatives | Negatives meet the 500 planning floor |
|---|---|---|---|
| train | 140 | 138 | NO |
| validation | 31 | 29 | NO |
| test | 29 | 33 | NO |

The 500 floor is a planning value (spec section 3.2). Below it, the split is not feasible for a specificity claim; AUROC is still reported and labelled.
