# Counts: lectin validation truth set

Source: `controls.tsv` (729 rows, 729 sequences, 539 clusters in all). Clusters: MMseqs2 17-b804f, 30% identity, 50% coverage, all rows together. Numbers from `../results/counts.tsv`.

| Class | Label | Rows | Sequences | Clusters | Species | Ascomycota clusters | Basidiomycota clusters | `in_reference` | `none` |
|---|---|---|---|---|---|---|---|---|---|
| structure_backed | positive | 39 | 39 | 23 | 27 | 6 | 17 | 20 | 19 |
| host_binding_h1a | positive | 4 | 4 | 3 | 3 | 3 | 0 | 3 | 1 |
| chitin_masking_h1b | positive | 10 | 10 | 5 | 6 | 5 | 0 | 10 | 0 |
| literature_binding | positive | 2 | 2 | 2 | 2 | 0 | 2 | 1 | 1 |
| known_gap | positive | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 1 |
| **all primary positives** | positive | 56 | 56 | 30 | 37 | 12 | 18 | 34 | 22 |
| pa14_owned_elsewhere (apart) | positive | 10 | 10 | 6 | 4 | 6 | 0 | 0 | 10 |
| er_housekeeping | negative | 18 | 18 | 8 | 13 | 8 | 1 | 18 | 0 |
| enzyme_lectin_like_domain | negative | 12 | 12 | 2 | 7 | 2 | 0 | 8 | 4 |
| n1_shared | negative | 136 | 136 | 74 | 44 | 65 | 13 | 3 | 133 |
| **all reasoned negatives** | negative | 166 | 166 | 84 | 54 | 75 | 14 | 29 | 137 |
| assumed (no reason) | negative | 497 | 497 | 430 | 148 | 381 | 43 | 0 | 497 |

Cluster counts per class are the number of distinct clusters with a member in that class. One cluster (representative Q12713) holds paracoccin (positive) and two N1 endochitinases (negative).

## Size targets (COMMON-RULES rule 7)

| Target | Met? |
|---|---|
| 20 positive clusters | Pooled: 30. Not in any one species (most: 3, *Coprinopsis cinerea*) |
| 20 negative clusters | Pooled reasoned: 84. Not in any one species (most: 19) |
| Interval half-width 0.10 or less | No. Union recall by cluster 19/30, Wilson 0.46 to 0.78 |

I did not pad any class. No species entry reaches `estimated`.
