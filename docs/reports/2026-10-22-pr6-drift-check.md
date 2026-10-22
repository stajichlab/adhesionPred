# PR 6 drift check of the iterated-only queries (job 30282220)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Run record: `docs/reports/data/pr6_drift_check_2026-10-22/`. Script and freeze: `analysis/wall_fusion_truth/pr6_drift_check.py`. Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md` (decisions 5 and 7). Owner choice c of 2026-10-22.*

Descriptive only: no model is built, and no call, status, tier, `narrow_ok` or panel row changes. Every threshold below is a PROPOSAL of the freeze in `pr6_drift_check.py`, not a decision. Nothing in the script or the freeze is changed by this report.

## 1. What was asked

aga1, ece1 and hyd1_like are `buildable_iterated_only`: the one-pass `phmmer` set does not meet the training rule, the `jackhmmer` set does (`docs/reports/2026-10-22-pr6-gap-search-uniref50.md`, with the clade-bucket re-derivation in `docs/reports/data/pr6_gap_uniref50_2026-10-22_rederived/`). Iteration can pull in unrelated proteins that share only a composition. The check asks of each of the 233 `jackhmmer` outside fungal hits whether a direct pairwise `phmmer` alignment to the query supports it (E <= 1e-5, aligned part >= 0.5 of the query or of the hit, and the same search against the hit with its residues shuffled scores less than half), then applies the training rule (5 sequences in at least 2 named clades) to the supported hits only.

## 2. Result (recomputed from `drift_hits.tsv`)

| Query | Hits | Supported | bias_driven | Unsupported | Supported named clades | Verdict |
|---|---|---|---|---|---|---|
| hyd1_like (A0ACE4KFX7) | 43 | 32 | 0 | 11 | 3 | `supported_rule_met` |
| aga1 (P32323) | 64 | 11 | 0 | 53 | 8 | `supported_rule_met` |
| ece1 (Q07730) | 126 | 0 | 0 | 126 | 0 | `not_supported` |

- **ece1: the drift is confirmed.** Not one of the 126 hits has a pairwise hit to the query that passes the E-value and the coverage test. 123 have no pairwise hit at all, and 3 have a best domain E-value above 1e-5. ece1 is not buildable. The 126 outside hits of the `jackhmmer` run (24 clades) and the 6,896 non-fungal hits in the round 2 record are artefacts of iteration.
- **hyd1_like:** 32 of 43 hits are supported, in 3 named clades (Dothideomycetes 7, Leotiomycetes 3, Sordariomycetes 19) plus 3 hits with no assigned clade, which are not a clade. Their query coverage is 0.531 to 0.99 (median 0.724). 11 hits are unsupported: 2 have no pairwise hit and 9 have a best E-value above 1e-5.
- **aga1:** 11 of 64 hits are supported, in 8 named clades (Agaricomycetes 2, Glomeromycota 2, Kickxellomycota 2, and one each in Malasseziomycetes, Microsporidia, Pezizomycetes, Pucciniomycetes and Sordariomycetes). Of the other 53, 31 have no pairwise hit and 22 do not pass, among them 9 that pass the E-value and fail the coverage.

## 3. A finding for the owner: the aga1 verdict rests on short local matches

The coverage rule of the freeze is "the LARGER of the query coverage and the hit coverage is at least 0.5". For aga1 this is what lets all 11 hits through:

- **Query coverage of the 11 supported hits is 0.148 to 0.31 (median 0.217)** of the 725-aa Aga1 (the length in `analysis/wall_fusion_truth/results/row_mapping.tsv`). Not one reaches 0.5.
- Their hit coverage is 0.549 to 0.946 (median 0.736): each hit is mostly covered by a short stretch of the query.
- The Aga1 query is 64 percent S/T/N/Q/G/P with a low-complexity score of 0.67 (hyd1_like: 50 percent and 0.16; ece1: 31 percent and 0.11). Short Ser/Thr-rich matches are what a composition effect would produce.
- **The shuffled null did not help here.** `bias_driven` is 0 for every query, but this is not evidence of no bias: every residue-shuffled hit scored 0 (no hit even at E 10), so the null ratio is 0.000 for every hit that has a real pairwise hit (32 for hyd1_like, 20 for aga1, 0 for ece1). With one shuffle per hit, the null cannot separate a composition match from a homolog in this data.

So the aga1 verdict `supported_rule_met` rests on short local matches. It is not a family-level homology claim for Aga1.

### Sensitivity: the same table read with QUERY coverage

This is a descriptive re-reading of the committed `drift_hits.tsv`, **not a re-classification and not a new freeze**. It keeps the same E-value (<= 1e-5), the same null threshold (< 0.5) and the same rule on named clades, and only replaces "larger of query and hit coverage" by "query coverage >= 0.5".

| Query | Supported (committed freeze) | Named clades | Supported (query coverage >= 0.5) | Named clades | Verdict under query coverage |
|---|---|---|---|---|---|
| hyd1_like | 32 | 3 | 32 | 3 | `supported_rule_met` (unchanged) |
| aga1 | 11 | 8 | 0 | 0 | `not_supported` (changes) |
| ece1 | 0 | 0 | 0 | 0 | `not_supported` (unchanged) |

**The choice of the coverage definition is an owner decision (made, section 4b).** Recommendation to the owner: use the query coverage. A custom model has to represent the Kiss protein, so a hit that matches only a fifth of it is not the same family member; hyd1_like and ece1 do not change, and aga1 does. Nothing is changed in the script or in the freeze here.

## 4. What this means for the PR 6 targets

| Target | State after the drift check |
|---|---|
| ecm33 (KS0086) | buildable on the one-pass `phmmer` set (round 2); drafts and freeze proposal in `analysis/wall_fusion_truth/README_pr6_ecm33.md` |
| hyd1_like (KS0017) | buildable on its 32 supported hits under either coverage definition; adhesion slice, so leave-cluster-out reported `partial` (decision 4) |
| aga1 (KS0080) | left out of the PR 6 HMM work (owner decision d, below): `not_supported` under the query coverage; needs a different modeling approach |
| ece1 (KS0107) | not buildable: drift confirmed (0 supported hits of 126) |

## 4b. Owner decision (2026-10-22)

The owner chose **option d** for the coverage definition: **the rule is query coverage >= 0.5** (a supported hit must align across at least half of the Kiss protein); the sensitivity table above stays in the record as a descriptive column showing what the looser "larger of query and hit coverage" rule would have given. Consequences:

- **aga1 (KS0080)** is **left out of the PR 6 HMM work** (`not_supported`: 0 supported hits). The owner's reading: aga1 **needs a more involved modeling approach**, not a homology model. Its evidence is short local matches (15 to 31 percent of the protein) in a protein that is 64 percent S/T/N/Q/G/P, which points to sequence-feature or composition-aware modeling (and possibly a domain-level model of one segment, option c of the choice, which would need the alignment coordinates this run did not save). That is outside this PR; it stays a known gap.
- **ece1 (KS0107)** is **left out**: `not_buildable` (drift confirmed, 0 supported of 126).
- **hyd1_like (KS0017)** stays a PR 6 target (32 supported hits, 3 named clades, query coverage 0.531 to 0.99); adhesion slice, so leave-cluster-out reported `partial` (decision 4).
- The PR 6 target set is therefore **ecm33 and hyd1_like**. The committed script (`pr6_drift_check.py`) still classifies with the larger-coverage rule; the owner decision is applied here by reading the committed table (a change to the script's rule is a separate change).

## 5. What the check does not show

- A `supported_homolog` is a pairwise alignment passing three proposed thresholds. It is not a curated family member.
- The shuffled null is one shuffle per hit (seed 1) and scored 0 for every hit, so it carried no information in this run; more shuffles would be a possible upgrade.
- The thresholds (1e-5, 0.5, 0.5, 5 sequences, 2 clades) are proposals, and the sensitivity table shows the one that matters for aga1.
- UniRef50 hits are cluster representatives, so the counts are conservative.
