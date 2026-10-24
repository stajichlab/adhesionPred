# PR 6 family-specificity results: do the candidate Pfam families hit things they should not?

*2026-10-24. Author: Claude Code (claude-sonnet-5-5). Run record: `docs/reports/data/pr6_family_specificity_results_2026-10-24/` (README there lists jobs, commits, hashes). Design: `analysis/wall_fusion_truth/README_pr6_family_specificity.md`.*

Status: a screening step only. Nothing enters the panel, no tier, `narrow_ok` or status changes, and `family_measure.tsv` and `family_table.tsv` are untouched. `narrow_ok` stays owner-signed.

## 1. What was run

Two HPCC jobs searched Pfam 38.2 family models (`hmmsearch` at each family's gathering cutoff, GA) against three sets:

- the family's **own members**: the evidenced Kiss proteins that made the family a candidate;
- **hard negatives**: 161 proteins from the shared negative set (137 N1 plus 24 N2) that should not be adhesion or wall members;
- **four proteomes** (one per species: A. fumigatus, C. albicans, C. immitis, S. cerevisiae; 32,491 proteins). These are not labelled, so a hit there is a count, not an error.

| run | job | code commit | families |
|---|---|---|---|
| run1 | 30423745 | `e36c6e4` | 30 candidate families (Pfam hit covering at least half the protein, or E <= 1e-30) |
| run2 | 30426504 | `f562773` | 3 added by the owner: Big1 (BIG1_C PF27678), KNR4/SMI1 (SMI1_KNR4 PF09346), DFG10 (Steroid_dh PF02544) |

In both runs: the code pin matched, every model was fetched from Pfam 38.2 with a GA line, the inputs were 132 positives and 161 negatives, and no panel or tier changed. Leakage is `tuned_on_truth` (the families were chosen because they hit their own members) and the tier is capped at `smoke`.

**Own-member recall is circular.** It is 100 percent for all 33 families (30 plus 3), and that is by construction, not a result.

Terms: "hits per 10,000" is `10000 * hits / 32,491` over the four proteomes. "Non-member" below means not one of the family's own Kiss members. The reference point is the hydrophobin family from earlier work, which hit 5 per 10,000 (about 16.2 proteins over 32,491).

## 2. Run1: 30 candidate families

**Hard negatives.** Only 3 of the 161 hit: PF00082 Peptidase_S8 hits J9VDC2; PF02055 Glyco_hydro_30 hits Q4WBR2 and Q8J0I9 (the same negative cluster). The other 28 families hit no negative.

**Proteomes.** Families above 5 per 10,000 (hits as reported by the job, `kept_hits_total`):

| family | hits | per 10,000 | note |
|---|---|---|---|
| PF00082 Peptidase_S8 | 37 | 11.39 | also hits one hard negative |
| PF00534 Glycos_transf_1 | 23 | 7.08 | |
| PF13641 Glyco_tranf_2_3 | 20 | 6.16 | also hits 6 other wall-group Kiss positives (3 of them chitin synthases inside the A. fumigatus proteome) |
| PF02366 PMT | 18 | 5.54 | also hits 1 other wall positive; 2 of the 18 are its own members |
| PF00982 Glyco_transf_20 | 17 | 5.23 | 1 of the 17 is its own member |
| PF16192 PMT_4TMC | 17 | 5.23 | 2 of the 17 are its own members |

The job total includes own members that sit in the proteomes. Leaving those out (`summary.tsv`, `kept_per10k_excl_own`), Peptidase_S8 is 36 (11.08), Glycos_transf_1 23 (7.08), Glyco_tranf_2_3 20 (6.16), and PMT 16 (4.92), Glyco_transf_20 16 (4.92), PMT_4TMC 15 (4.62) fall just under 5. So three of the six are above 5 either way and three sit on the line depending on the count used.

Other points:
- Hyphal_reg_CWP (PF11765, already in the family table): 14 hits (4.31 per 10,000), 12 of them in C. albicans and 2 in S. cerevisiae. These are likely paralogs, but that was not checked.
- Rep1 (PF29797) hits nothing in the proteomes and no other protein; only its own member (circular).
- The remaining families give 0.00 to 3.69 per 10,000 (`summary.tsv`).

## 3. Run2: the three extra families

| family | own members (circular recall) | negatives hit (of 161) | other positives hit | proteome hits | per 10,000 |
|---|---|---|---|---|---|
| BIG1_C PF27678 | 1/1 | 0 | 0 | 1 | 0.31 |
| SMI1_KNR4 PF09346 | 2/2 | 0 | 0 | 5 | 1.54 |
| Steroid_dh PF02544 | 1/1 | 0 | 0 | 10 | 3.08 |

All three are under 5 per 10,000, and none of their proteome hits is a Kiss positive.

## 4. What this means

These counts are a screen. A family with few negative and proteome hits is not shown to be specific, and a family with many is not shown to be broad. They say which families deserve a closer look before the owner decides on a panel addition.

## 5. What was NOT established

- Hits were not checked to be true non-members. A proteome hit may be a real family member (a paralog or an orthologue of the adhesion or wall protein), which would be a correct call, not a false one. The 12 C. albicans hits of Hyphal_reg_CWP are an example: likely paralogs, unverified.
- There is no recall test. Families defined by one or two evidenced members have no independent set to measure sensitivity on; own-member recall is circular.
- The four proteomes are one per species and unlabelled; the 161 negatives are a shared set, not chosen per family.
- No interval, sensitivity or specificity was computed, and no family should be called "specific" beyond the counts above.

## 6. Questions for the owner (not decisions)

1. Which of the three families above 5 per 10,000 (Peptidase_S8, Glycos_transf_1, Glyco_tranf_2_3) are worth checking hit by hit for paralogs, before any panel step?
2. Is a general-purpose enzyme domain such as Peptidase_S8 or a glycosyltransferase family wanted as a panel candidate at all, given the hit rates?
