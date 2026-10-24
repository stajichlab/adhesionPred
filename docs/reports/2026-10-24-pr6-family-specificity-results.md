# PR 6 family-specificity results: do the candidate Pfam families hit things they should not?

*2026-10-24. Author: Claude Code (claude-sonnet-5-5). Run record: `docs/reports/data/pr6_family_specificity_results_2026-10-24/` (README there lists jobs, commits, hashes). Design: `analysis/wall_fusion_truth/README_pr6_family_specificity.md`.*

## Correction (2026-10-24, after #141)

What was wrong in #141: `summary.tsv` decided that a proteome hit is an own member by matching the Kiss UniProt accession to the proteome protein id. Only Afum_Af293_UniProt headers carry `sp|ACC|`; Calb_SC5314 uses CGD Assembly 22 ids (`C1_13450W_A`), Scer_S288C uses systematic names (`YNL238W`) and Cimm_RS uses RefSeq ids. So own members in C. albicans and S. cerevisiae were counted as non-member hits, which raised the non-member counts and rates. Nothing else was wrong: the HPCC tables (`run1/`, `run2/`, `hpcc/`) are unchanged byte for byte, `kept_hits_total` and `kept_per10k` come from the job and do not exclude anything, and the `other_positive_hits_*` columns come from the positives FASTA by Kiss accession and were not affected.

What changed: a committed map from every Kiss positive (all 132) to its proteome id (`kiss_positive_proteome_map.tsv`, with provenance: UniProt release 2026_03, queried 2026-10-24), a shared helper (`analysis/wall_fusion_truth/pr6_kiss_proteome_map.py`) and a summary that now gives three counts per family: own-member hits, hits that are Kiss positives of another cluster, and hits that are not Kiss positives. Numbers below are the corrected ones. The three corrections to the text of #141 are in sections 2 and 3 (the "three sit on the line" sentence, the Hyphal_reg_CWP paralog sentence, and the statement that none of the three extra families hits a Kiss positive). Also stated plainly: "tier capped at `smoke`" is recorded only in run2's `provenance.tsv` (`tier_cap`); for run1 it follows from the script design (the script derives `smoke` at most and writes only report tables), not from run1's output.

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

In both runs: the code pin matched, every model was fetched from Pfam 38.2 with a GA line, the inputs were 132 positives and 161 negatives, and no panel or tier changed. Leakage is `tuned_on_truth` (the families were chosen because they hit their own members) and the tier is capped at `smoke` (recorded only in run2's provenance; run1's follows from the script design, see the correction note).

**Own-member recall is circular.** It is 100 percent for all 33 families (30 plus 3), and that is by construction, not a result.

Terms: "hits per 10,000" is `10000 * hits / 32,491` over the four proteomes. "Own" is a hit that is one of the family's own Kiss members (matched by the id its proteome uses, through the committed map). "Other Kiss positive" is a hit that is a Kiss positive of another cluster. "Non-Kiss" is any other hit. "Non-member" means not an own member, so it is other Kiss positives plus non-Kiss hits.

The limit is the hydrophobin precedent, written as "at most 5 unlabelled extra calls per 10,000 proteins" (about 16.2 proteins, so 16 hits pass and 17 do not, over 32,491). Here it is read on **non-member hits** (all hits minus own members; `kept_per10k_excl_own`). The precedent says "unlabelled", which could also mean that Kiss positives of other clusters are not counted; that is ambiguous, so both rates are shown (`kept_per10k_excl_own` and `non_kiss_per10k`). Only Glyco_tranf_2_3 changes side between the two readings (20 vs 16 hits, 6.16 vs 4.92).

## 2. Run1: 30 candidate families

**Hard negatives.** Only 3 of the 161 hit: PF00082 Peptidase_S8 hits J9VDC2; PF02055 Glyco_hydro_30 hits Q4WBR2 and Q8J0I9 (the same negative cluster). The other 28 families hit no negative.

**Proteomes.** Families above 5 per 10,000, corrected. `hits` is the job's `kept_hits_total` (raw); own, other-Kiss-positive (other KP) and non-Kiss are the three groups of a hit; the last two rate columns are non-member hits (own removed, other KP kept) and non-Kiss hits only.

| family | hits | per 10,000 raw | own | other KP | non-Kiss | per 10,000 non-member | per 10,000 non-Kiss | above 5 (non-member) |
|---|---|---|---|---|---|---|---|---|
| PF00082 Peptidase_S8 | 37 | 11.39 | 3 | 0 | 34 | 10.46 | 10.46 | yes |
| PF00534 Glycos_transf_1 | 23 | 7.08 | 0 | 0 | 23 | 7.08 | 7.08 | yes |
| PF13641 Glyco_tranf_2_3 | 20 | 6.16 | 0 | 4 | 16 | 6.16 | 4.92 | yes (no, if other KP are not counted) |
| PF02366 PMT | 18 | 5.54 | 5 | 0 | 13 | 4.00 | 4.00 | no |
| PF00982 Glyco_transf_20 | 17 | 5.23 | 2 | 0 | 15 | 4.62 | 4.62 | no |
| PF16192 PMT_4TMC | 17 | 5.23 | 5 | 0 | 12 | 3.69 | 3.69 | no |

Notes: Peptidase_S8 also hits one hard negative. Glyco_tranf_2_3 hits 6 other wall-group Kiss positives in the positives FASTA (3 chitin synthases inside the A. fumigatus proteome); of those, 4 are in the kept proteomes and are among its 20 proteome hits: P78746, Q4X0H6 and Q4X0H7 in A. fumigatus and CHS3 (P29465, YBR023C) in S. cerevisiae. PMT's one other wall positive hit is not in a kept proteome.

Under the corrected counts, **three families are above 5 per 10,000** on the non-member reading: Peptidase_S8, Glycos_transf_1 and Glyco_tranf_2_3 (Glyco_tranf_2_3 only because its 4 other-Kiss-positive hits count; without them it is 16 hits, 4.92, just under). PMT, Glyco_transf_20 and PMT_4TMC are under 5: with their own members correctly removed they sit at 4.00, 4.62 and 3.69, not "on the line". The #141 sentence that three of the six "sit on the line depending on the count" rested on the broken count (it removed only the own members that happen to carry a UniProt id) and is withdrawn. Whether a hit at 16 (4.92) passes depends on reading "at most 5" as at most 16.2 proteins; 17 hits do not.

Other points:
- Hyphal_reg_CWP (PF11765, already in the family table): 14 hits (4.31 per 10,000 raw), 12 of them in C. albicans and 2 in S. cerevisiae. One of the 14 is its own member: IFF4 (Q5AAL9, `CR_00610W_A`), so 13 non-member hits (4.00). Looking the 13 ids up in UniProt (2026-10-24, by ordered-locus name) gives, in C. albicans: HYR1 (`C1_13450W_A`), IFF6 (`C2_09130C_A`), FLO9 (`C3_00580W_A`, unreviewed entry A0A1D8PIY8), IFF11 (`C3_00600W_A`), IFF5 (`C4_06550C_A`), IFF8 (`C5_00710W_A`; UniProt entry Q59XB0 names it IFF7 in the protein name and IFF8 as the gene), HYR3 (`C5_00730W_A`), RBR3 (`C7_03290C_A`), HYR4 (`CR_00760C_A`), IFF3 (`CR_03630W_A`) and IFF9 (`CR_03880W_A`); in S. cerevisiae: HPF1 (`YOL155C`) and CSS1 (`YIL169C`). So the 11 other C. albicans hits are named members of the Hyr/Iff family (plus RBR3 and FLO9), i.e. likely correct family calls, not false positives, and the two S. cerevisiae hits are the haze protective factor proteins. What was verified: the gene names that UniProt gives for those ids. What was not verified: the domain boundaries or an alignment, whether the proteome sequences equal the UniProt sequences (proteome FASTA files are on HPCC only), or whether these proteins belong in the panel.
- Rep1 (PF29797) hits nothing in the proteomes and no other protein; only its own member (circular).
- The remaining families give 0.00 to 4.00 per 10,000 non-member hits (`summary.tsv`). Corrections to own-member counts: Peptidase_S8 3 (was 1), PMT 5 (2), PMT_4TMC 5 (2), Glyco_transf_20 2 (1), Hyphal_reg_CWP 1 (0), Trehalose_PPase 2 (1), SIS 2 (1), STE3 2 (1), fn3_2 2 (0), ChAPs, Chs7, SKN1_KRE6_Sbg1, STE2, KRE1 and Glyco_transf_24 1 each (0), and the three run2 families 1 each (0). Glycos_transf_1 is unchanged (23, 7.08).

## 3. Run2: the three extra families

| family | own members (circular recall) | negatives hit (of 161) | other positives hit | proteome hits | per 10,000 |
|---|---|---|---|---|---|
| BIG1_C PF27678 | 1/1 | 0 | 0 | 1 | 0.31 |
| SMI1_KNR4 PF09346 | 2/2 | 0 | 0 | 5 | 1.54 |
| Steroid_dh PF02544 | 1/1 | 0 | 0 | 10 | 3.08 |

All three are under 5 per 10,000. Each hits its own member inside a proteome, so the non-member hits are 9 (2.77) for Steroid_dh (own member P40526 DFG10, `YIL049W`), 4 (1.23) for SMI1_KNR4 (own member A0A1D8PE53 SMI1, `C1_07870C_A`) and 0 (0.00) for BIG1_C, whose single proteome hit is its own member BIG1 (P38813, `YHR101C`). None of the three hits another cluster's Kiss positive. The #141 statement that none of their proteome hits is a Kiss positive was wrong: one hit each is.

## 4. What this means

These counts are a screen. A family with few negative and proteome hits is not shown to be specific, and a family with many is not shown to be broad. They say which families deserve a closer look before the owner decides on a panel addition.

## 5. What was NOT established

- Hits were not checked to be true non-members. A proteome hit may be a real family member (a paralog or an orthologue of the adhesion or wall protein), which would be a correct call, not a false one. The 11 other C. albicans hits of Hyphal_reg_CWP are an example: UniProt names them as Hyr/Iff family members, so they are family members by name; no alignment was done.
- The Kiss-positive map rests on UniProt ids and names (release 2026_03). Whether a proteome entry has the same sequence as the UniProt entry, and whether each mapped id exists in its proteome, cannot be checked here (proteome FASTA files are on HPCC only). 34 of the 85 mapped ids appear in the proteome hit tables, which is the only independent evidence for the id format. Cimm_RS has no Kiss positive (none of the 112 UniProt positives is a Coccidioides entry), and no positive is unmapped.
- There is no recall test. Families defined by one or two evidenced members have no independent set to measure sensitivity on; own-member recall is circular.
- The four proteomes are one per species and unlabelled; the 161 negatives are a shared set, not chosen per family.
- No interval, sensitivity or specificity was computed, and no family should be called "specific" beyond the counts above.

## 6. Questions for the owner (not decisions)

1. Which of the three families above 5 per 10,000 on the non-member reading (Peptidase_S8, Glycos_transf_1, Glyco_tranf_2_3) are worth checking hit by hit for paralogs, before any panel step? And should hits that are Kiss positives of other clusters count as extra calls (Glyco_tranf_2_3: 6.16 with them, 4.92 without)?
2. Is a general-purpose enzyme domain such as Peptidase_S8 or a glycosyltransferase family wanted as a panel candidate at all, given the hit rates?
