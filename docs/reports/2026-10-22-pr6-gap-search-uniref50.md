# PR 6 gap search, round 2 (UniRef50, job 30262455)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Run record: `docs/reports/data/pr6_gap_uniref50_2026-10-22/`. Round 1 (Swiss-Prot): `2026-10-22-pr6-gap-search-swissprot.md`. Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md` (decisions 3, 5, 7, S3). Scripts: `analysis/wall_fusion_truth/pr6_gap_search.sbatch`, `pr6_gap_summary.py`.*

## 1. What was run

The same four queries as round 1 (the clusters with no Pfam hit) were searched with `phmmer` (one pass) and `jackhmmer -N 3` at `--incE 1e-5` against UniRef50 2025_03 (70,198,728 sequences, sha256 `bc52bdd6…`). The owner's decision was UniRef50 for all four queries. The job ran with `-p batch --time=8:00:00 -c 8 --mem=16G` and finished in 3 h 19 min (exit 0, peak memory about 1.6 GB). Code `e11808b`, HMMER 3.4, Pfam 38.2 with the same sha256 as the earlier jobs. Settings are unchanged from round 1.

Caveats that apply to every number below:
- UniRef50 is about 122 times larger than Swiss-Prot (70,198,728 against 573,661 sequences), so the same E-value cutoff is stricter and round 1 and round 2 counts are not strictly comparable.
- Each hit is the representative of a 50 percent identity cluster, so a count of outside sequences is a count of clusters. It is conservative. The `n=` member counts are not used.
- The 5-in-2-clades rule is the provisional rule of decision 7. The script counts a fungal hit whose fine clade could not be assigned as a clade of its own, called `(no clade)`. See section 3, hyd1_like.

## 2. Result per query

Counted from `gap_homologs.tsv.gz` (hits per query, method and class) and `gap_per_query.tsv`. "Inside" is the held-out fine clade; "self" is the query's own cluster; non-fungal hits are listed apart. "Species" is the column of `gap_per_query.tsv` (species level; `gap_homologs.tsv.gz` holds taxids, so a count of distinct taxids can differ by a few).

| Cluster | Gene | Held-out | Method | Outside seq | Outside clades | Outside species | Inside | Non-fungal | Status |
|---|---|---|---|---|---|---|---|---|---|
| KS0017 | hyd1_like | Eurotiomycetes | phmmer | 5 | 2 (1 named + `(no clade)`) | 4 | 8 | 6 | |
| | | | jackhmmer | 43 | 4 (3 named + `(no clade)`) | 25 | 16 | 0 | `buildable` |
| KS0080 | aga1 | Saccharomycotina | phmmer | 4 | 4 | 4 | 1 | 134 | |
| | | | jackhmmer | 64 | 19 | 52 | 9 | 2,688 | `buildable_iterated_only` |
| KS0107 | ece1 | Saccharomycotina | phmmer | 0 | 0 | 0 | 2 | 6 | |
| | | | jackhmmer | 126 | 24 (23 named + `(no clade)`) | 103 | 15 | 6,896 | `buildable_iterated_only` |
| KS0086 | ecm33 | Saccharomycotina | phmmer | 379 | 19 (18 named + `(no clade)`) | 244 | 106 | 195 | |
| | | | jackhmmer | 953 | 30 (29 named + `(no clade)`) | 427 | 120 | 35,290 | `buildable` |

The status is decided by the one-pass set: `buildable` when `phmmer` meets the rule, `buildable_iterated_only` when only `jackhmmer` does (a flag of possible drift).

Round 1 against round 2 (outside sequences, `phmmer` / `jackhmmer`):

| Gene | Round 1, Swiss-Prot | Round 2, UniRef50 |
|---|---|---|
| hyd1_like | 0 / 1, `not_buildable` | 5 / 43, `buildable` |
| aga1 | 0 / 0, `not_buildable` | 4 / 64, `buildable_iterated_only` |
| ece1 | 0 / 0, `not_buildable` | 0 / 126, `buildable_iterated_only` |
| ecm33 | 4 / 4 (2 clades), `not_buildable` | 379 / 953, `buildable` |

## 3. What the numbers say

- **ecm33 (KS0086, Q5AGC4).** The one-pass set has 379 outside fungal sequences in 19 clade labels (18 named clades) and 244 species. Ten named clades have at least 5 sequences each (largest: Dothideomycetes 106, Sordariomycetes 57, Eurotiomycetes 55, Leotiomycetes 38, Kickxellomycota 37, Mucoromycota 33). The strict reading (at least 5 sequences in each of 2 clades) is also met. This is a widespread fungal family. The shipped model `PF12454` does not reach its GA threshold on either Ecm33 protein (Q5AGC4 score 5.5, P38248 score 0.6, GA 23.2, both `below_GA`; the check is identical to round 1). So Ecm33 is a clear candidate for a custom model and for leave-clade-out in the wall slice. Note that `jackhmmer` adds 35,290 non-fungal hits for this query (the one-pass set has 195), so any alignment should be built from the one-pass set or checked first.
- **hyd1_like (KS0017, A0ACE4KFX7).** The one-pass set has 5 outside fungal sequences, which meets the provisional minimum of 5 exactly. They fall in 2 labels: Sordariomycetes (3 sequences) and `(no clade)` (2 sequences, fungal hits whose fine clade could not be assigned by the script). Only one named clade is in the one-pass set, so the rule is met only because the script counts `(no clade)` as a clade. The strict reading (5 in each of 2 clades) is not met. The iterated set has 43 sequences (Sordariomycetes 27, Dothideomycetes 10, Leotiomycetes 3, `(no clade)` 3), with 3 named clades, and no non-fungal hit. The status `buildable` is therefore marginal. This behaviour of the script was found while writing this report and was not changed here.
- **aga1 (KS0080, P32323) and ece1 (KS0107, Q07730).** Both are `buildable_iterated_only`. The one-pass set has 4 outside sequences for aga1 (4 clades, below the minimum of 5 sequences) and 0 for ece1. `jackhmmer` reaches 64 sequences in 19 clades for aga1 and 126 in 24 clades for ece1, but it also adds 2,688 and 6,896 non-fungal hits, against 64 and 126 outside fungal hits. These are small proteins that are rich in repeats or low complexity or are peptide precursors, and iteration can pull in compositionally similar but unrelated proteins. The flag is a drift flag. **These two are not buildable until a drift check is done.** This report does not claim that a check was made.

## 4. What a drift check would compute (not done)

From the per-iteration `--tblout` tables in `raw/` on HPCC (not committed): the iteration in which each hit and each clade entered the profile; the fraction of the query covered by the alignment of each hit (a hit that aligns only a repeat or a low-complexity stretch is suspect); the same search with a composition-bias filter on, and whether the outside hits survive; whether the outside hits share the query's motif or domain (a domain scan of a sample of them); the order in which clades were added (a clade that enters only at iteration 3 and brings thousands of non-fungal hits with it is suspect). These are descriptions of the computation, not results.

## 5. What this changes and what it does not

- No model was built. No call, status, tier, panel row or `narrow_ok` is changed.
- Decision 5 asked how many clade-only clusters have any outside homolog. For the four gap clusters the answer is: in Swiss-Prot (round 1), 2 of 4 have any outside hit (hyd1_like 1, ecm33 4) and none meets the rule; in UniRef50, all 4 have an outside hit in the one-pass or the iterated set (3 of 4 in the one-pass set; ece1 only after iteration), and 2 of 4 (ecm33 and hyd1_like) have a one-pass set that meets the rule, with hyd1_like marginal.

## 6. Open choices for the owner (not decided here)

1. Which targets go to the build step 6.1: ecm33 clearly; hyd1_like only with the caveat of section 3; aga1 and ece1 only after a drift check.
2. Whether to run the drift check for aga1 and ece1 (one more HPCC job, reading `raw/` and the target database), or to leave them as known gaps.
3. Whether the script should stop counting `(no clade)` as a clade (a small code change and a new test; the verdict for hyd1_like would then be `not_buildable` on the one-pass set).
