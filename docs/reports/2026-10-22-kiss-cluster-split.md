# Kiss 2019 control set: cluster split (Task 5.3)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Record: `docs/reports/data/kiss_cluster_split_2026-10-22/` (job 30229366). Every number below was recomputed from `clusters_split.tsv`, `members_split.tsv` and `clade_floors.tsv` of that record. No recall, sensitivity or other performance number is made here (that is Task 5.4).*

## 1. What was clustered

293 sequences: 132 positive controls (adhesion 36, wall 91, fusion 5) and the 161 shared hard negatives. MMseqs2 at 30% identity and 50% coverage gave **152 clusters**; 35 of the 76 positive clusters have more than one member (largest 9). The 35 *A. fumigatus* family representatives were not clustered (owner decision 3); they keep their Kiss `cluster_number` (adhesion 4, fusion 3, wall 20 clusters).

| class and slice | clusters |
|---|---|
| positive, adhesion | 28 |
| positive, wall | 44 |
| positive, fusion | 4 |
| negative only | 76 |
| mixed (a control and a negative in one cluster) | 14 |
| all | 152 |

The 14 mixed clusters are 13 wall and 1 adhesion. A mixed cluster is judged by its control members for `seed_class_slice`.

## 2. Seed classes

A cluster is `seed-like` if any member reaches 30% identity (target coverage 0.5) to a Pfam seed sequence of the models used. Two labels (decision 4): `slice` (own-slice models: wall and fusion against the 21 wall and fusion models, adhesion against the 18 active family rows) and `all39`.

| positive clusters | n | seed-like (slice) | seed-independent (slice) | seed-like (all39) | seed-independent (all39) |
|---|---|---|---|---|---|
| adhesion | 28 | 16 | 12 | 16 | 12 |
| wall | 44 | 14 | 30 | 17 | 27 |
| fusion | 4 | 0 | 4 | 0 | 4 |

Three wall clusters are `seed-independent` for the wall models and `seed-like` for all 39: they resemble a seed of an adhesion-family model. Of the 90 clusters that hold a negative, 32 are `seed-like` (all39) and 58 are `seed-independent`.

Five models have very few seed sequences (Hydrophobin_D 3, Glyco_hydro_17 6, Chitin_synth_2 7, Hyd1F 8); a `seed-independent` label that depends on them rests on little. The seeds are the alignments InterPro served on 2026-10-22.

## 3. Which slices reach the floor of 5 clusters

Rule (owner decisions 3 and 4): the floor counts all clusters present in a held-out clade; a stratum is `estimated` only with at least 5 clusters of its own, else `descriptive`; fusion is always `descriptive`. Here `clade_only` means that every positive member of the cluster is in that clade and `shared` that the cluster also has positive members in another clade, taken from the MMseqs2 clusters (a different definition from the Kiss `cluster_number` used in the feasibility report, so the counts differ from it). Full table: `clade_floors.tsv`.

**Fine level (primary), clusters present / clade-only / shared:**

| slice | clade | any | clade-only | shared | status (any; clade-only; shared) |
|---|---|---|---|---|---|
| adhesion | Saccharomycotina | 14 | 14 | 0 | estimated; estimated; none |
| adhesion | Sordariomycetes | 7 | 7 | 0 | estimated; estimated; none |
| adhesion | Eurotiomycetes | 3 | 3 | 0 | descriptive; descriptive; none |
| adhesion | Ustilaginomycetes | 3 | 3 | 0 | descriptive; descriptive; none |
| adhesion | Tremellomycetes | 1 | 1 | 0 | descriptive; descriptive; none |
| wall | Saccharomycotina | 32 | 19 | 13 | estimated; estimated; estimated |
| wall | Eurotiomycetes | 17 | 5 | 12 | estimated; estimated; estimated |
| wall | Sordariomycetes | 15 | 4 | 11 | estimated; descriptive; estimated |
| wall | Chytridiomycota, Taphrinomycotina, Tremellomycetes, Ustilaginomycetes | 1 to 2 each | 0 to 1 | 1 to 2 | descriptive |
| fusion | Saccharomycotina | 2 | 1 | 1 | descriptive |
| fusion | Eurotiomycetes, Sordariomycetes, Ustilaginomycetes | 1 each | 0 or 1 | 0 or 1 | descriptive |

**Coarse level (second, harder test):** adhesion Saccharomycotina 14 and Pezizomycotina 10 (both clade-only, `estimated`); wall Saccharomycotina 32 (19 clade-only, 13 shared) and Pezizomycotina 24 (11 clade-only, 13 shared), all `estimated`; every Basidiomycota, Chytridiomycota and Taphrinomycotina entry and all of fusion stay `descriptive`.

What reaches the floor: adhesion at Saccharomycotina and Sordariomycetes (fine), and Saccharomycotina and Pezizomycotina (coarse); wall at Saccharomycotina, Eurotiomycetes and Sordariomycetes (fine; Sordariomycetes clade-only has 4 and stays `descriptive` there). Everything else, and all of fusion, is descriptive. Compared with the feasibility report (Kiss `cluster_number`), the MMseqs2 clusters split some Kiss clusters, so the adhesion counts are higher (Saccharomycotina 14 instead of 10, Sordariomycetes 7 instead of 6) and the Eurotiomycetes wall clade-only count reaches 5 instead of 4. The eligible clades per rule do not change.

## 4. What this does not say

No recall or specificity. The clusters and seed classes only decide how Task 5.4 reports recall. A control set where about half of the adhesion and a third of the wall clusters are `seed-like` means that a Pfam model recalling them is partly recalling its own seed family; Task 5.4 reports `seed-like` and `seed-independent` apart for that reason.

## 5. Control set

`data/controls/wall-fusion-kiss2019/clusters.tsv` now holds these MMseqs2 cluster ids (it held the Kiss proxy before). `controls.tsv` and `controls_local_sequence.tsv` carry `cluster_id`, `seed_class_slice`, `seed_class_all39` and `leakage` (`none` or `in_reference`). The 112 + 20 rows fall in 76 clusters (67 in `controls.tsv`). The L1 validator passes.
