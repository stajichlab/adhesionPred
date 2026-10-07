# Adjudication of the repeat-mechanism labels

Rows (adhesin, E1 or E2): 96. Rule: highest evidence tier wins; a same-tier conflict goes to an expert.

## How rows were settled

| Rule | Rows |
|---|---|
| agreement | 74 |
| higher_tier_wins | 15 |
| no_claim | 5 |
| conflict_at_same_tier | 2 |

## Final mechanism label

| Label | Rows |
|---|---|
| 2a | 61 |
| 2b-iii | 1 |
| 2c | 3 |
| 2d | 2 |
| other | 22 |
| unknown | 5 |
| unresolved | 2 |

## Independent clusters with a repeat region (30% identity, 50% coverage), cumulative by evidence tier

| Evidence tier included | Clusters | Needed |
|---|---|---|
| paper_stated | 8 | 20 |
| + database_annotation (UniProt repeat features) | 19 | 20 |
| + family_inference | 35 | 20 |

Votes of `other`, `2c` or `2d` are not counted against repeats. They mean that the cited paper gives another mechanism or none.

A cluster count does not make a status. `estimated` also needs a 95% half-width of at most 0.10, which needs about 61 clusters at sensitivity 0.8 (binomial approximation).

## Expert review packet

2 rows. File: `expert_review_packet.tsv`. Each row has both agents' labels, sources and quotes, and blank `expert_*` columns.
