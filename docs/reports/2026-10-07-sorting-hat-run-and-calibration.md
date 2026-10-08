# Sorting hat: first end-to-end runs and calibration

Research use only. This is not a regulatory allergenicity assessment and not a diagnostic result. No row is supported by an IgE, antibody or T-cell measurement. `signal_peptide_protein` means that SignalP calls a signal peptide and nothing more. Peptide-level calls only; glycan epitopes are not assessed. WHO/IUIS lists no *Coccidioides* allergen; delayed-type hypersensitivity skin-test reactivity is T-cell mediated and is not what the allergen columns address. The antigen ranking covers the *C. immitis* RS reference only (not *C. posadasii*). Non-protein adhesins such as galactosaminogalactan are not detected.

Allergen source: WHO/IUIS Allergen Nomenclature Sub-Committee (allergen.org), extract of 2026-10-04, 111 sequences from 27 species used after cleaning.

*2026-10-07. Plan 2, Tasks 11 to 15. Branch `sorting-hat-run-report`. Everything here was run on the UCR HPCC on 2026-10-06 and 2026-10-07. All data files are under `docs/reports/data/sorting_hat/`. The person who wrote this did not review any Pfam hit or any call as a biologist. Where a section says "not reviewed", a reviewer must do it.*

## 0. Summary

- Seven proteome runs finished end to end with exit code 0: *A. fumigatus* Af293 (UniProt and Fungi_5k annotations), A1163, W72310, *C. immitis* RS, and *S. cerevisiae* S288C (with and without dubious ORFs).
- Every module job finished in under 11 minutes. The longest was TMHMM (5 to 10 minutes). SignalP on a GPU took 46 to 252 seconds.
- The only module with a measured status in these runs is `step1_rule@R0`, and only in the Af293 UniProt run. Its status is `smoke` for Af293. All other module and taxon pairs are `unvalidated`.
- The Pfam module reports `unavailable`. No family is signed off (decision C4). The review tables are ready for an expert.
- The allergen module finds 65 of 111 known allergens by the 35% rule and 41 of 111 by the 70% rule when the only references are allergens of other species.
- The R0 call is stable across the three *A. fumigatus* annotations for sequences that are identical, and mostly stable for non-identical pairs (kappa 0.65 to 0.94).
- FLO11 in S288C is not called by the repeat call, although the panel expects it. Section 8.

## 1. Proteomes

| name | source | retrieved | proteins | unique sequences | invalid | sha256 (first 12) |
|---|---|---|---|---|---|---|
| Afum_Af293_UniProt | UniProt UP000002530 (the Phase C gene models) | 2026-10-06 | 9,647 | 9,647 | 0 | `7cfae982f8d2` |
| Afum_Af293_Fungi5k | Fungi_5k input (a different annotation) | 2026-10-06 | 9,161 | 9,157 | 0 | `252b09e7788f` |
| Afum_A1163 | UniProt UP000001699 | 2026-10-06 | 9,942 | 9,942 | 0 | `b9bd8dc4101b` |
| Afum_W72310 | NCBI GCA_040167795.1 (UCR_Afum_W72310_1.0) | 2026-10-06 | 10,556 | 10,521 | 0 | `dda9946615d9` |
| Cimm_RS | NCBI RefSeq GCF_000149335.2 | 2026-10-06 | 9,910 | 9,875 | 0 | `927d60dfccc6` |
| Scer_S288C | SGD `orf_trans_all.fasta.gz`, with dubious ORFs | 2026-10-06 | 6,722 | 6,638 | 8 | `17e8b47e1ae2` |
| Scer_S288C_nodubious | the same file, 683 `Dubious ORF` entries removed | 2026-10-06 | 6,039 | 5,975 | 8 | `55398250be75` |

Full records: `data/sorting_hat/provenance/*.json`. The eight invalid S288C sequences are excluded from all modules (`na_invalid`). *C. neoformans* H99 (UniProt UP000010091) was searched with Pfam and TMHMM only, for the family review.

Taxa checked against the NCBI taxdump of 2026-10-06: Af293 330879, A1163 451804 and W72310 746128 are all *A. fumigatus* (746128). W72310 has no strain node. The UniProt proteome records give the same IDs for Af293 and A1163.

## 2. Jobs

Table of every job: `data/sorting_hat/calibration/sacct_all.txt`. Summary by job type (elapsed seconds over all proteomes of 6,039 to 10,556 proteins):

| job | n | elapsed (s) | partition |
|---|---|---|---|
| Pfam (15 models, `--cut_ga`) | 8 | 18 to 47 | short |
| repeats (two detectors) | 7 | 97 to 161 | short |
| BLAST against 111 allergens | 7 | 41 to 60 | short |
| TMHMM | 8 | 326 to 612 | short |
| SignalP 6 GPU, fast mode | 7 | 46 to 252 | short_gpu |
| reciprocal best hit BLAST (6 searches) | 1 | 966 | short |

Notes:
- The `exfab` partition has one GPU node. Two SignalP jobs waited there with a start estimated a week ahead. Jobs on `short_gpu` started within hours. The lab account cap on `short_gpu` is four GPUs and other lab jobs also use it. A multi-partition request is refused (the accounts are bound to one partition).
- `hmmsearch` with several models cannot read a gzip FASTA. The Cys-rich Pfam script needs a plain FASTA.
- The module `ncbi-blast/2.14.0+` is not listed on every node. The jobs used 2.14.0+ on compute nodes.
- The plan sizes jobs to about one hour. These jobs run for seconds to minutes. A scatter over many proteomes should bundle several proteomes in one job.

## 3. Calibration status

Only rule R0 has a measured status. The status file comes from the Phase C metrics (run 2026-10-02, widened C grid). Its intervals are the widest of the file's interval and the Wilson interval on the cluster count (decision 2). Leakage: overlap of the Phase C positives with the SignalP 6 training data was not measured (decision 1).

| taxon | Phase C label | calibration set | positives / negatives | clusters (pos / neg) | sensitivity [95% interval] | specificity [95% interval] |
|---|---|---|---|---|---|---|
| *S. cerevisiae* (4932) | smoke | S1:Scer_SGD | 79 / 3,785 | 58 / 3,156 | 0.848 [0.734, 0.943] | 0.966 [0.958, 0.972] |
| *C. albicans* (5476) | smoke | S1:Calb_CGD | 153 / 459 | 113 / 410 | 0.477 [0.356, 0.591] | 0.943 [0.913, 0.968] |
| *A. fumigatus* (746128) | smoke | S3:Afum_ASPFU | 19 / 45 | 17 / 38 | 0.947 [0.738, 1.000] | 1.000 [0.908, 1.000] |
| *A. nidulans* (162425) | **estimated** | S3:Anid_EMENI | 109 / 164 | 100 / 151 | 0.688 [0.592, 0.779] | 0.988 [0.955, 1.000] |
| *C. neoformans* (5207) | smoke | S3:Cneo_H99_GOA | 7 / 32 | 6 / 31 | 0.857 [0.458, 1.000] | 0.938 [0.796, 1.000] |
| *U. maydis* (5270) | smoke | S3:Umay_MYCMD | 9 / 28 | 9 / 24 | 1.000 [0.701, 1.000] | 0.893 [0.712, 1.000]  |

- The Af293 lower bound (0.738) is lower than the 0.833 that the plan expected. The cluster-count interval is wider than the file's own interval.
- The status file was first written to the Af293 UniProt work directory only. On 2026-10-07 the owner decided that it also goes to the two S288C runs, because Phase C measured R0 on the same SGD file, and to the *C. albicans* SC5314 run, because Phase C measured R0 on the same CGD file (after removing 38 duplicate IDs with identical sequences). The S288C runs report R0 as `smoke` (6,714 and 6,031 proteins; the 8 invalid sequences stay `unvalidated`). The *C. albicans* run reports `smoke` with sensitivity 0.477 [0.356, 0.591]. Fungi_5k, A1163, W72310 and RS keep no status source, so R0 is `unvalidated` there. Fungi_5k, A1163 and W72310 are other annotations than the one measured. RS has no truth set.
- Expected false calls and predictive values: not computed. The prevalence of surface proteins is assumed, not measured, in all genomes here.
- Every other module has no calibration: `tm`, `repeat02`, `repeat14`, `allergen_homology`, `antigen_lookup`, `cys_rich`, `expression`, `pfam_adhesion`, `pfam_allergen`. The reason per module is in section 9.

## 4. Allergen homology

**Leave-species-out recall** (`allergen-lso`, BLAST 2.14.0+, 111 sequences from 27 species; the reference for each allergen is every allergen of another species):

| rule | allergens found |
|---|---|
| `iuis_allergen_similarity` (identity at least 35% over at least 80 aa, a single local alignment) | 65 / 111 (59%) |
| `iuis_allergen_homolog` (identity at least 70% and coverage at least 80%) | 41 / 111 (37%) |

This is sensitivity against allergens of other species. It is not a specificity and writes no status entry.

**Counts on the *A. fumigatus* annotations:**

| run | proteins | similarity called | homolog called |
|---|---|---|---|
| Af293 UniProt | 9,647 | 104 | 41 |
| Af293 Fungi_5k | 9,161 | 105 | 41 |
| A1163 | 9,942 | 104 | 39 |
| W72310 | 10,556 | 115 | 43 |
| *C. immitis* RS | 9,910 | 84 | 16 |
| S288C | 6,722 | 84 | 15 |

The Fungi_5k numbers equal the numbers measured on 2026-10-05 (105 and 41).

**Same-species and other-species hits in Af293 UniProt.** Of 104 similarity calls, 51 hit an allergen of *A. fumigatus* and 53 hit another species. Of 41 homolog calls, 30 are same-species and 11 are other-species (all 11 have exposure route "Airway"). 28 proteins reach 95% identity and 80% coverage, all same-species. The Af293 counts are self-recognition: the reference holds the *A. fumigatus* allergens. They are not module performance. The plan quoted "30 reached 95%" for Fungi_5k without a definition. The count here uses identity at least 95% and coverage at least 80%.

**Limits.**
- Five IUIS entries are not searchable and were skipped (111 were written): Epi p 1.0101 (free text) and four peptide fragments of 9 to 29 residues (Asp fl 13, Cur l 1, Tri t 1, Tri t 4). Two of these are Onygenales.
- A hit means sequence similarity. It does not show IgE binding.
- S288C has 15 homolog calls. They match Asp f 19, Asp f 12 and Asp f 23 (*A. fumigatus*), Cand a 1, Mala s 6 and Rho m 1 and 2. Why these proteins are similar was not examined.

## 5. Pfam family review (not signed off)

Pfam 38.2 (sha256 `4b0da6399b97`), HMMER 3.4, `--cut_ga`, 15 models. No family is active. The module writes `unavailable` for `pfam_adhesion` and `pfam_allergen`, and `wall_family_domain` is `not_assessable` for every protein in every run.

Facts per family and proteome (`data/sorting_hat/specificity/<proteome>/review_table.tsv`). "Cys ≥ 8" counts hits that have at least 8 Cys in the envelope of the best domain. "Draft members" counts hits that are in the draft member list.

| family | proteome | hits | with signal peptide (R0) | with TMHMM helix | Cys ≥ 8 | draft members |
|---|---|---|---|---|---|---|
| PF00399 PIR | Scer_S288C | 10 | 10 | 0 | 0 | 5 |
| PF01185 Hydrophobin | Afum_Af293_UniProt | 6 | 6 | 1 | 6 | 6 |
| PF01185 Hydrophobin | Afum_A1163 | 5 | 5 | 1 | 5 | n/a |
| PF01185 Hydrophobin | Afum_W72310 | 6 | 6 | 0 | 6 | n/a |
| PF01185 Hydrophobin | Cimm_RS | 1 | 1 | 0 | 1 | n/a |
| PF04681 Bys1 | Afum_Af293_UniProt | 5 | 5 | 0 | 0 | 1 |
| PF04681 Bys1 | Afum_A1163 / Afum_W72310 | 5 / 5 | 5 / 5 | 0 / 0 | 0 / 0 | n/a |
| PF04681 Bys1 | Cimm_RS | 2 | 2 | 1 | 0 | n/a |
| PF05730 CFEM | Afum_Af293_UniProt | 4 | 4 | 0 | 3 | 3 |
| PF05730 CFEM | Afum_A1163 | 4 | 4 | 0 | 3 | n/a |
| PF05730 CFEM | Afum_W72310 | 7 | 6 | 2 | 6 | n/a |
| PF05730 CFEM | Scer_S288C | 1 | 1 | 0 | 1 | 0 |
| PF05730 CFEM | Cimm_RS | 8 | 7 | 1 | 6 | n/a |
| PF05730 CFEM | Cneo_H99 | 3 | not run | 1 | 3 | n/a |
| PF07691 PA14 | Afum_Af293_UniProt / A1163 / W72310 | 6 / 6 / 5 | 0 / 0 / 0 | 0 / 0 / 0 | 0 | 0 |
| PF07691 PA14 | Scer_S288C | 6 | 4 | 4 | 0 | 4 |
| PF07691 PA14 | Cneo_H99 | 5 | not run | 0 | 0 | n/a |
| PF10182 Flo11 | Scer_S288C | 2 | 1 | 1 | 0 | 1 |
| PF10528 GLEYA | Scer_S288C | 6 | 4 | 4 | 0 | 4 |
| PF11765 Hyphal_reg_CWP | Scer_S288C | 2 | 2 | 0 | 0 | n/a |
| PF11766 Candida_ALS_N | Scer_S288C | 1 | 1 | 1 | 0 | n/a |
| PF25312 Allergen_Asp_f_4 | Afum_Af293_UniProt | 3 | 2 | 0 | 0 | 1 |
| PF25312 Allergen_Asp_f_4 | Afum_A1163 / Afum_W72310 | 3 / 3 | 2 / 2 | 0 / 0 | 0 | n/a |
| PF28987 DewD | Afum_W72310 / Cimm_RS | 1 / 1 | 1 / 1 | 0 / 1 | 1 / 1 | n/a |

What the facts show (not a reviewer's classification):
- **PA14.** All 17 PA14 hits in the three *A. fumigatus* proteomes lack a signal peptide and have 0 or 1 Cys. Four of the six Af293 hits are named beta-glucosidases (BGLH, BGLI, BGLJ, BGLK). Two are unnamed. The 5 H99 hits have 0 to 4 Cys. SignalP was not run on H99, so their signal peptide status is unknown. The `signal_peptide` second condition in the family table removes all of them. In S288C, 4 of 6 have a signal peptide and are the flocculins.
- **CFEM.** The three known Af293 CFEM proteins (CfmA-C) carry 8 or 9 Cys in the domain. A fourth hit (Q4WYI0, 764 aa) has 7. The three *C. neoformans* H99 hits have 8 or more Cys, and one has a TMHMM helix. SignalP was not run on H99. W72310 has 7 CFEM hits against 4 in Af293. Two have a helix.
- **Hydrophobin.** Six Af293 proteins (RodA, B, C, E, F, G) have the Pfam hydrophobin domain with 8 Cys or more. RodD (hydrophobin-like) has no hit. A1163 has 5 hits.
- **Bys1.** CalA (Q4WXJ1) is a hit. Four other proteins also hit. None has 8 Cys. Their function is not known to this report.
- **PIR (S288C).** 10 proteins hit. Five are draft members (PIR1, PIR3, PIR5, HSP150, CIS3). The other five (TIR1, TIR2, CWP1, CWP2, ANS1) are not in the draft list. SGD describes TIR1 and CWP1 as cell wall mannoproteins and CWP2 as a covalently linked cell wall protein. ANS1 is an uncharacterized ORF. All ten hits have a signal peptide.
- A1163 and W72310 have no draft members, because their protein IDs differ. Their tables list hits only.

The draft member lists are in `data/sorting_hat/specificity/*/members.tsv` with the source of each row in `members.rationale.tsv`. They come from UniProt and SGD names and from the papers in `docs/paper/05-literature-verification.md`. They were not reviewed by an expert. The sensitivity and specificity lines in the `*.specificity.tsv` files treat all non-members as negatives. **Do not report them as specificity.** The column `review_class` is empty. A reviewer must classify each non-member hit as false domain hit, uncharacterised true member or receptor-like. Decision C4 (sign-off) is open. No change was made to `data/sorting_hat/family_table.tsv`.

## 6. Antigen lookup, Cys-rich proteins and expression in *C. immitis* RS

- The RefSeq proteins were mapped to the antigen ranking through `protein_map.tsv`. 9,139 proteins are `ok` and 771 are `not_in_reference`. The `idmap_method` is gene-to-best-transcript. 126 RefSeq isoforms inherit values of another transcript. An exact-sequence match is not implemented.
- `cocci_specificity_rank_top15` is called for 1,376 proteins (15% of 9,139). 143 of them have a signal peptide. This is the same 143 as in the call `serodiagnostic_marker_candidate`. The plan expected 143.
- All 14 Tier 1 and all 45 Tier 2 candidates (`analysis/cocci_antigens/`) are called by the top 15% cut. This is not an independent test. The cut was set after the four anchors were seen, and the tiers come from the same ranking. Ag2/PRA, PRA2 and PRA3 are one family. 15% of the proteome is called.
- Cys-rich finder on the RefSeq FASTA (Pfam search of its own, six models): 9,910 proteins, 460 with a signal peptide, 21 `cys_rich_sp_unassigned`, 5 `cys_rich_sp_known_family`, 31 `cys_rich_no_sp`. The same counts as the earlier run on the FungiDB file.
- R0 calls a signal peptide for 460 of 9,910 RS proteins.
- Sequence homology does not transfer epitopes. The ranking is a genus-specificity ranking made from sequence comparison. It is not epitope prediction and was not tested against serology.

## 7. Stability across *A. fumigatus* annotations

Reference: Af293 UniProt. Pairs by best reciprocal BLAST hit (BLAST 2.14.0+, e-value 1e-5). Script: `analysis/sorting_hat_run/rbh_stability.py`. Summaries: `data/sorting_hat/rbh/*.summary.txt`.

| comparison | reciprocal pairs | identical sequences | non-identical pairs |
|---|---|---|---|
| Af293 Fungi_5k | 8,723 | 4,743 | 3,980 |
| A1163 | 9,385 | 5,190 | 4,195 |
| W72310 | 8,923 | 1,409 | 7,514 |

Identical sequences give identical calls for every call and every pair. That tests determinism of the software. It does not test biology.

Non-identical pairs, `signal_peptide_protein[R0]`:

| comparison | agree | called/called | called/not called | not called/called | not called/not called | kappa |
|---|---|---|---|---|---|---|
| Af293 vs Fungi_5k | 3,820 / 3,980 | 164 | 98 | 62 | 3,656 | 0.651 |
| Af293 vs A1163 | 4,157 / 4,195 | 333 | 19 | 19 | 3,824 | 0.941 |
| Af293 vs W72310 | 7,359 / 7,514 | 527 | 78 | 77 | 6,832 | 0.861 |

Other calls on non-identical pairs: `tandem_repeat_protein` kappa 0.889 (Fungi_5k), 0.719 (A1163), 0.785 (W72310); `iuis_allergen_similarity` kappa 0.950, 0.989, 0.977. The differences come from different gene models, mostly different start sites, and were not classified one by one. The bioinformatics review of 2026-10-05 gave 5 of 80 changed R0 calls (6.3%) for pairs that share the C-terminus but differ at the start. The reciprocal-hit comparison gives 160 of 3,980 (4.0%) for Fungi_5k. The two numbers use different pair sets.

Named cases (Af293 accession): CspA (Q4WXC4), RodA (P41746), CalA (Q4WXJ1), Asp f 2 (P79017), Asp f 4 (O60024), Asp f 7 (O42799), Asp f 15 (O60022). In all three comparisons every call of these seven proteins is equal. CspA has 93.7% identity to its A1163 hit and 92.5% identity (coverage 0.82) to its W72310 hit. The repeat region of CspA varies between isolates according to the curated seed. The tandem repeat call did not change.

Limits:
- The Fungi_5k and W72310 gene models come from other pipelines than UniProt. A difference in calls can come from annotation.
- Only Af293 has Phase C truth. A1163 and W72310 are for stability only.
- W72310 has only 1,409 identical sequences to Af293. The reciprocal-hit pairs include proteins with 90% to 95% identity.

## 8. Panel check (report only)

Panel: `data/sorting_hat/panel.tsv`, 11 rows. Results in `data/sorting_hat/calibration/panel_*.txt`. Nothing is gated on it.

| run | agree | disagree | excluded (leakage) | known miss | not in run |
|---|---|---|---|---|---|
| Af293 UniProt | 1 | 0 | 8 | 1 | 1 |
| *C. immitis* RS | 1 | 0 | 8 | 1 | 1 |
| S288C | 0 | 0 | 8 | 1 | 2 |

Most rows are excluded because the panel proteins were used to tune the modules (spec section 4). The excluded rows still show what was observed:
- SOWgp (`XP_001245172.2`): `tandem_repeat_protein` and `cocci_specificity_rank_top15` are both called, as expected.
- Ag2/PRA (`XP_001240075.1`, exact match to UniProt A0A0E1RVD3): `wall_family_domain` is `not_assessable`. No Pfam family covers it, and no family is active.
- RodA and CalA in Af293: `wall_family_domain` is `not_assessable`, for the same reason.
- **FLO11 (`YIR019C`) in S288C: `tandem_repeat_protein` is `not_called`, but the panel expects `called`.** This is a miss of the detectors, checked on 2026-10-07. The call needs coverage at least 0.25 and at least 2.5 copies. `repeat02` finds a weak period-15 signal in residues 753 to 815 only (4.10 copies, coverage 0.045; the protein has 1,368 aa and 50% Ser/Thr). `repeat14` finds no period (z_seq 3.31, no region). The other S288C flocculins: FLO1, FLO5 and FLO9 are called by both detectors (coverage 0.29 to 0.53); FLO10 is called by `repeat14` only (period 63, coverage 0.269). The statement in `analysis/cocci_repeats/REPORT_2026-09-30_pfam_confirmation.md` that FLO11 is a tandem array with coverage near 1.0 is not supported by this measurement. It describes class 2a, not a FLO11 result. FLO11 is also a draft member of the Flo11 Pfam family, so the Flo11 domain would call it once that family is active. The repeat call is therefore a limited evidence source for Ser/Thr-rich, diverged repeats.
- Gel1: R0 calls a signal peptide in Af293 and in RS (agree).
- Asp f 1 (mitogillin) and Asp f 2: `iuis_allergen_homolog` is called in Af293. The expected result is guaranteed, because both are in the reference set.
- The RS Hsp60-like protein (by annotation only): R0 does not call a signal peptide. This matches the expected known miss for moonlighting proteins.

Not done: Asp f 34 (accession not found in the proteome headers), the Als1 and Rbt5 rows, a cross-species allergen row, and Hsp60 of *Histoplasma*.

## 9. What is not measured, and what blocks it

| module or call | state | blocked by |
|---|---|---|
| `step1_rule@R0` | `smoke` for Af293, S288C, *C. albicans*, *C. neoformans*, *U. maydis*. `estimated` for *A. nidulans* only | Basidiomycota and GPI truth (curation). Overlap with the SignalP 6 training set is not measured |
| `step1_rule@R1`, `R2`, `step1_ml@card` | `unavailable` | owner decision on rule, ML or hybrid |
| `tm`, `allergen_homology` | `unvalidated` | no truth table (C2 allergen negatives from IEDB assays) |
| `repeat02`, `repeat14`, `tandem_repeat_protein` | `unvalidated` | per-call status entries (engine change); 8 / 19 / 35 clusters by evidence tier against a floor of 20; decision on family inference |
| `pfam_adhesion`, `pfam_allergen`, `wall_family_domain` | `unavailable` | family sign-off (C4); expert classification of non-member hits; reviewed member lists (agent task 04) |
| `antigen_lookup`, `cocci_specificity_rank_top15`, `serodiagnostic_marker_candidate` | `unvalidated` | leave-species-out antigen test (C3); only four anchors exist |
| `cys_rich`, `expression` | `unvalidated` | no truth table |
| composite calls (for example `cell_wall_adhesion_candidate`) | status from the weakest deciding module | decision 3 (strict reading gives `unvalidated`) |

Other gaps:
- Per-sha256 cache, GPU out-of-memory retry, a sliding-window allergen rule, exact sequence-hash ID mapping and automatic reciprocal-best-hit code in the core tool are not built. The reciprocal-hit script here is an analysis script, not part of the package.
- The CI job for the package and the job scripts under CI were not run. The job scripts were run under SLURM and worked. Two scripts needed a workaround (section 2).
- The draft member lists, the panel proteins and the interpretation lines in sections 5 to 8 need an expert review before any claim goes into a paper.
