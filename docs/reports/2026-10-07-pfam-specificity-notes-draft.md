# Pfam specificity notes: draft for sign-off

*2026-10-07. Draft. `data/sorting_hat/family_table.tsv` is **not** changed and no family is `active`. The notes below are in the patch `docs/reports/data/sorting_hat/specificity/family_table_notes.patch` (3 lines changed; `git apply --check` passes; the package loader reads the patched table: 15 families, none active). Counts come from `pfam_proposals.tsv` and the owner-reviewed sheet `review_sheet_first_signoff.tsv`. They describe 8 proteomes and unreviewed draft member lists. They are not measures of sensitivity or specificity.*

## What the classes rest on

- 16 rows were classified by the owner on 2026-10-07 (CFEM 7, PA14 8, Hydrophobin 1).
- All other hits carry classes proposed by `pfam_proposals.py` from table facts. They are not reviewed.

## Notes

| family | hits | classes at the time of writing | `second_condition` |
|---|---|---|---|
| CFEM (PF05730) | 33 in 8 proteomes | 26 `uncharacterised_true_member` (7 owner, 19 proposed), 3 `known_member`, 4 `receptor_like` (proposed) | `no_tm` (unchanged) |
| PA14 (PF07691) | 29 in 6 proteomes | 25 `false_domain_hit` (8 owner, 17 proposed), 4 `known_member` | `signal_peptide` (unchanged) |
| Hydrophobin (PF01185) | 18 in 4 proteomes | 12 `uncharacterised_true_member` (1 owner, 11 proposed), 6 `known_member` | none |

The new `specificity_note` of each family is the text in the patch. Facts checked while writing:

- **CFEM.** `no_tm` uses `n_tm_mature` (helices that start inside the signal peptide are ignored). It removes exactly four hits: W72310 KAK9636691.1 (6 mature helices), W72310 KAK9642401.1 (1), RS XP_001239408.1 (1) and H99 J9VVT8 (1). These are the four proposed `receptor_like` rows. CSA2 (*C. albicans* C4_06920C_A) has one raw helix inside the signal peptide, so it is **called**, although the owner considers it not directly cell wall related.
- **PA14.** The 25 false domain hits have no signal peptide. The only four hits with a signal peptide are *S. cerevisiae* FLO1, FLO5, FLO9 and FLO10. FLO11 is not a PA14 hit (it is a Flo11 family hit).
- **Hydrophobin.** All 18 hits have 8 Cys and a signal peptide. RodD has no hit.

## Consequences of activating each family

An active family makes `wall_family_domain` `called` for a protein that has a hit and passes the second condition. The status stays `unvalidated` (decision 4). The call is domain evidence, not an adhesin call.

- CFEM: 29 of 33 hits would be called. The called set includes RS Ag2/PRA and PRA2, the three Af293 CfmA to C proteins, *C. albicans* PGA7, RBT5, PGA10, CSA1, SSR1 and CSA2.
- PA14: 4 of 29 would be called.
- Hydrophobin: all 18 would be called.
- `cell_wall_adhesion_candidate` (domain or repeat, and a signal peptide) would then call these proteins too.

## Sign-off (owner)

Setting `active = yes`, `active_by` and `active_date` for a family is the owner's decision. Per the plan, that change is committed alone, with the review tables as evidence. This draft does not make it.
