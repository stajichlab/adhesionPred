# Round 2 of the PR 6 gap search: verdicts re-derived under the fixed clade rule

*2026-10-22. Source: the committed run record `../pr6_gap_uniref50_2026-10-22/` (job 30262455); that record is NOT edited.*

**Why.** The round 2 summary counted the `(no clade)` bucket (fungal hits with no assigned fine clade) as a clade, so the
one-pass set of hyd1_like (3 hits in Sordariomycetes plus 2 hits with no assigned clade) was called `buildable`. The rule
now counts named clades only (`pr6_gap_summary.training_rule`): such hits count in the total of sequences, but they are
not a clade for the 2-clade condition or the strict column. Owner choice c of 2026-10-22.

**How.** `analysis/wall_fusion_truth/pr6_gap_rederive.py` reads `gap_per_query.tsv` and `gap_homologs.tsv.gz` of the
record, applies the fixed rule and writes `rederived_verdicts.tsv` (old and new verdicts side by side). Nothing is
searched again.

| query | gene | one-pass (phmmer): sequences, named clades | iterated (jackhmmer): sequences, named clades | status in the record | status under the fixed rule |
|---|---|---|---|---|---|
| A0ACE4KFX7 | hyd1_like | 5, 1 (2 unassigned) | 43, 3 (3 unassigned) | `buildable` | **`buildable_iterated_only`** |
| P32323 | aga1 | 4, 4 | 64, 19 | `buildable_iterated_only` | unchanged |
| Q07730 | ece1 | 0, 0 | 126, 23 (1 unassigned) | `buildable_iterated_only` | unchanged |
| Q5AGC4 | ecm33 | 379, 18 (5 unassigned) | 953, 29 (20 unassigned) | `buildable` | unchanged |

**What differs from what was expected.** hyd1_like does not become `not_buildable`. The one-pass set no longer meets the
rule (5 sequences, 1 named clade), but the iterated set does (43 sequences in 3 named clades: Sordariomycetes 27,
Dothideomycetes 10, Leotiomycetes 3), so the status is `buildable_iterated_only`: the same drift flag as aga1 and ece1.
hyd1_like is therefore not buildable on the conservative one-pass set; whether the iterated set is usable depends on a
drift check that has not been done. ecm33 stays `buildable` on the one-pass set. The strict reading (5 sequences in each
of 2 clades) is unchanged in every row (`strict_committed` equals `strict_fixed` in `rederived_verdicts.tsv`).
