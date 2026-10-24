# Scoping a model of the AGA1 class: Ser/Thr-rich, repeat-bearing, GPI-anchored adhesin-like proteins

*2026-10-24. Author: Claude Code (claude-sonnet-5-5). For the project owner. Local study only: no HPCC or SLURM step, no GPU, no change to any call, tier, panel, `narrow_ok` or freeze. Every number was recomputed from committed tables (tables, scripts and exact commands in `docs/reports/data/aga1_class_scoping/`, with its README). The one outside source is read-only UniProt REST, used only to fetch sequences, lineages and feature lines of accessions already named in committed tables (cache outside the repo, no sequence committed). Treat all numbers as descriptive: n is small.*

## Words used (one line each, at first use)

| Word | Plain meaning |
|---|---|
| Class (the "AGA1 class") | Cell-wall adhesin-like proteins that are rich in Ser/Thr (and Asn, Gln, Gly, Pro), carry repeats, and are held in the wall by a GPI anchor. Issue #121: Aga1 itself is one Kiss protein with no reliable relatives, so the target is the class. |
| GPI anchor | A lipid anchor added at the C-terminus that holds a protein in the outer cell wall. Here only the UniProt annotation is available. |
| Hard negative | A protein curated as NOT an adhesin that still looks like one (cell-wall enzymes, mucin-type sensors, other Ser/Thr-rich proteins). |
| ST-rich hard negative | A hard negative with some 30-residue window at least 50 percent Ser/Thr (`hn_st30`, threshold fixed before any model was run). |
| Clade | A fungal group. Fine clade = the rule of `map_kiss.py` (Saccharomycotina as one unit; the class for Pezizomycotina and Basidiomycota). |
| Cluster | A group of proteins at least 30 percent identical over 50 percent of their length (MMseqs2). Two proteins in one cluster are not independent. |
| Leave-clade-out | Test on one clade and train on all the others. The honest test for "does it work in a new lineage". |
| Leave-cluster-out | Hold out whole clusters (5 folds, 20 random repeats), so no close relative is in training. |
| AUC | Probability that a random positive scores above a random negative. 0.5 = chance, 1 = perfect. |
| Cluster bootstrap | An interval made by resampling whole clusters of the tested proteins. It covers sampling of the tested proteins, not training variability. |
| Composition features | Fractions of S, T, N, Q, G, P and Cys. No tool needed. |
| Detector 02 / 14 | The two existing tandem-repeat detectors behind the `tandem_repeat_protein` call. z_seq is detector 14's whole-sequence periodicity score. |
| Frozen embeddings | Protein-language-model vectors used as fixed input; the model is not fine-tuned. |

## 1. Summary

1. **Data exist for a general class finder, but nearly all of it is yeast.** Class positives (curated adhesin and repeat-bearing): **67 proteins in 31 MMseqs clusters; 35 proteins in 19 clusters with independent repeat evidence (tier A).** Tier A by clade: Saccharomycotina 21 (10 clusters), Taphrinomycotina 8 (5), Eurotiomycetes 4 (2), Sordariomycetes 2 (2); no Basidiomycota. Hard negatives: 172 proteins (100 clusters), of which 68 are ST-rich (43 clusters) and 19 carry at least 2 UniProt repeat features. 389 more proteome proteins are assumed negatives. The issue #15 fine-tuning gate (300 clusters) is far away: 19 tier-A clusters.
2. **Composition alone separates the class from hard negatives, less well from ST-rich ones.** Leave-cluster-out, tier A, composition only: AUC 0.923 [0.844, 0.974] against all hard negatives, **0.849 [0.707, 0.944] against ST-rich hard negatives**, 0.805 [0.556, 0.941] against hard negatives with UniProt repeats. The one clade fold with enough data (Saccharomycotina held out, 21 positives, 10 clusters; trained only on 14 non-yeast-clade positives): 0.920 [0.811, 0.987] (all hard negatives, n=28) and 0.907 [0.780, 0.990] (ST-rich, n=23). The other three folds have 8, 4 and 2 positives and too few negatives for an interval; they are listed but not interpreted.
3. **Against the existing repeat call, at the same number of false hits, the feature model recovers more positives, but the intervals overlap.** At 6 of 172 hard negatives called (3.5 percent, what the repeat call does): repeat call 18/35 tier-A positives = 0.514 [0.356, 0.670]; feature model 24/35 = 0.686 [0.520, 0.814]. The model's threshold was picked on the same out-of-fold scores, so its number is optimistic. No claim of improvement is made.
4. **AGA1 looks like the S/T-rich hard negatives SED1 and DAN4, not like FLO1.** AGA1: 725 aa, 55 percent Ser/Thr, 64 percent S/T/N/Q/G/P, a 132-residue S/T-rich run, no period found by either detector, z_seq 3.74. It scores above all but 2 of 28 Saccharomycotina hard negatives (SED1 and DAN4 score higher) in the leave-clade-out fold. Its profile is closer to SED1/DAN4/HPF1/FLO11 than to FLO1/5/9 (z_seq 10 to 18).
5. **Narrow classes: none meets the provisional rule.** Every candidate narrow class with at least 5 positives (Als-like 13, Hyr/Iff-like 12, Hyphal_reg_CWP-type 13, Flo1-type 5) lies in a single fine clade (Saccharomycotina). Pfl-like (8, Taphrinomycotina) and SOWgp (3, Eurotiomycetes) are single-clade too. Strict rule (at least 5 in each of 2 clades): no for all. The data say: one general finder; narrow slices as reports only.
6. **Next stage**: an HPCC step is not needed to decide the design; what is needed first is more non-yeast positives (curation) and owner decisions on the positive and negative definitions (section 7).

## 2. What exists: the inventory

Sources read (all committed): `analysis/calibration_truth/repeat_call_truth/truth_v3.{Scer_S288C,Calb_SC5314}.annotated.tsv`, `data/controls/repeat-mechanism/` (186-row curation table and the adjudicated consensus of 96 adhesins), `data/controls/repeat-adhesins-other-species/`, `data/controls/hard-negatives-shared/`, `data/controls/wall-fusion-kiss2019/`, `analysis/lectin_validation/controls/`, `data/curated/adhesins/adhesins.tsv` (Pfam list), `data/sorting_hat/family_table.tsv`, `docs/CLASSES.md`, `docs/class_descriptions.yaml`. 850 distinct proteins (accession-keyed; `proteins.tsv`).

**Roles** (`source_counts.tsv`; a protein can come from several tables):

| Role | n | Notes |
|---|---|---|
| class_positive (tier A / A+B) | 35 / 67 | adhesin curated E1 to E3 and repeat-bearing (rule in the data README) |
| hard_negative | 172 | N1 or N2 |
| background_assumed_negative | 389 | S288C 272 + C. albicans 117 proteome members with no curated class; negatives assumed, not shown |
| adhesin_other_mechanism | 49 | adhesins whose mechanism is CFEM, hydrophobin, moonlighting, etc.: not this class |
| excluded (indirect regulators 60, adhesion-phenotype only 36) | 96 | not used |
| uncurated_repeat_bearing | 7 | UniProt repeat features, no curated class |
| Kiss 2019 rows with no other role | 70 | wall/fusion/adhesion family rows; no negatives in that table |

**truth_v3 (what the repeat call is scored on).** S288C: 320 rows, 23 repeat-label-1 (UniProt at least 2 array-type repeat features, or a paper statement): adhesin 8, hard negative 10, uncurated 4, adhesion-phenotype 1. C. albicans: 224 rows, 13 label-1: adhesin 10, uncurated 3. Note: truth label 1 is "repeat-bearing", not "adhesin": MSB2, HKR1, EGT2, SED1, DAN4 are label 1 and curated hard negatives, so the repeat call cannot tell them from the class by design. C. albicans has no curated hard negatives in truth_v3.

**By clade (fine) and cluster** (`class_clade_counts.tsv`):

| Group | Saccharomycotina | Taphrinomycotina | Eurotiomycetes | Sordariomycetes | Other clades | All |
|---|---|---|---|---|---|---|
| Positives tier A | 21 (10 cl.) | 8 (5) | 4 (2) | 2 (2) | 0 | 35 (19) |
| Positives A+B | 51 (21) | 10 (6) | 4 (2) | 2 (2) | 0 | 67 (31) |
| Positives with UniProt GPI feature (A+B) | 31 (11) | 2 (2) | 0 | 2 (2) | 0 | 35 (15) |
| Hard negatives (all) | 28 (24) | 5 (5) | 58 (43) | 42 (29) | 39 | 172 (100) |
| ST-rich hard negatives | 23 (21) | 2 (2) | 19 (15) | 16 (10) | 8 | 68 (43) |
| Hard negatives with UniProt repeats | 10 (8) | 0 | 2 (1) | 2 (2) | 5 | 19 (12) |

Hard-negative strata (172): laccase_etc 38, wall_hydrolase 29, domain_non_member 23, gpi_wall_enzyme 22, st_linker_enzyme 21, repeat_non_adhesin 13, wall_structural_other 9, secretory_non_surface 5, mucin_sensor 5, hydrophobin_non_adhesin 4, moonlighting 1, unlabelled 2. Of the 68 ST-rich ones: gpi_wall_enzyme 17, st_linker_enzyme 16, domain_non_member 9, wall_structural_other 8, wall_hydrolase 6, mucin_sensor 5, repeat_non_adhesin 3, laccase_etc 2.

**Other inventory items.** Kiss adhesion slice: 31 rows (Saccharomycotina 18, Eurotiomycetes 5, Sordariomycetes 5, Ustilaginomycetes 2, Tremellomycetes 1); 10 of them (aga1, flo1, flo5, flo9, flo10, flo11, eap1, hwp1, iff4, als3/als1) are class positives; the other 21 are hydrophobins, small cysteine-rich proteins and a few S/T-rich cell-wall proteins (`kiss_adhesion_rows.tsv`; fig2, wap1/csa1, rbt1, dfi1 and bsc1 have an ST-rich window of at least 0.7). The 20 JGI-only Kiss rows have no committed sequence. Lectin controls (733 rows): 153 rows overlap the study, mostly the shared hard negatives (136 of the lectin negatives); 4 lectin positives overlap (stratum `pa14_owned_elsewhere`): FLO1 and FLO5 are class positives, EPA1 and EPA6 are adhesins of another mechanism; one lectin "assumed negative" (IFF4109) is a class positive (`lectin_overlap.tsv`).

**Calls that exist for the class** (`docs/CLASSES.md`): `tandem_repeat_protein` (status smoke: S288C 0.435 sensitivity / 0.987 specificity; C. albicans 0.462 / 0.976), `wall_family_narrow_domain` (7 Pfam families, all `smoke` tier and `narrow_ok=yes` in `family_table.tsv`), composite `cell_wall_adhesion_candidate` = signal peptide and (repeat or narrow domain). I recomputed the repeat call on truth_v3 from the locally re-run detectors: S288C 10/23 and 4/297 false hits; C. albicans 6/13 and 5/211, identical to the committed statuses, and it agrees with all 171 per-protein calls of `curated_rows_mapped.tsv` that could be matched (5 rows have no local feature row).

## 3. Sequences and features

**Sequences** (`sequence_availability.tsv`). Local committed `sequences.faa`: 283 proteins (170 hard negatives, 70 Kiss-only rows, 43 others); the other 567 come from UniProt (cache). All 850 have a sequence. Two ids (ACE80261, XP_045276491.1) are not UniProt accessions; both have local sequences.

**Features computed locally, stdlib only** (`features.tsv`): length; fractions of S, T, N, Q, G, P; S+T, S+T+P, S+T+N+Q+G+P; Cys count; `max_ST30` (highest S+T fraction in any 30-residue window); `lc_frac` (fraction of residues in a 30-window of entropy below 2.5 bits, a crude low-complexity measure); `mean_entropy30`; N-terminal residues before the first S/T-rich window; longest S/T-rich run (maximal segment at least 50 percent S/T); N-X-S/T sequon density (a crude proxy, not a predictor); repeat detectors 02 and 14 period, coverage, copies, z_seq (both ran locally in about 3 s); the existing repeat call and the exploratory z_seq 3.5 rule. UniProt annotation lines (not predictions): signal peptide, GPI-anchor, number of repeat features.

**Needs HPCC or an external tool later:** SignalP 6, a GPI predictor (PredGPI or similar), IUPred/AIUPred, glycosylation predictors (NetOGlyc and similar), ESM-2 embeddings, a repeat-unit profile, OrthoFinder.

**Two confounders to keep in mind.** Length: tier-A positives average 1,153 aa (17 of 35 are over 1,022 aa, the ESM-2 context limit), hard negatives 473 aa (4 of 172 over 1,022 aa). UniProt GPI annotation: 24 of 35 tier-A positives have it, but so do 35 of 170 hard negatives and 28 of the 68 ST-rich ones (GPI wall enzymes); the annotation alone is not a separator.

## 4. Baseline separation (descriptive)

Method: logistic regression (L2, balanced weights, standardised), trained on positives against hard negatives only; background proteins never trained on. Features: `comp` (S, T, N, Q, G, P, Cys fractions), `comp_len_struct` (adds log length, window and run features, N-terminal fraction, sequons), and a third set with repeat-detector numbers (`baseline_*.tsv`). Pre-declared: negative subsets and the 0.5 ST-rich threshold. Intervals only where at least 5 positives, 3 positive clusters, 10 negatives and 3 negative clusters exist; otherwise "n too small".

**Leave-cluster-out, 5 folds, pooled out-of-fold scores (seed 0), cluster bootstrap 95 percent interval**

| Positives | Features | vs all hard negatives (172) | vs ST-rich (68) | vs hard negatives with repeats (19) | vs assumed background (389) |
|---|---|---|---|---|---|
| tier A (35, 19 cl.) | comp | 0.923 [0.844, 0.974] | 0.849 [0.707, 0.944] | 0.805 [0.556, 0.941] | 0.927 [0.861, 0.967] |
| tier A | comp_len_struct | 0.943 [0.884, 0.980] | 0.877 [0.766, 0.953] | 0.809 [0.562, 0.948] | 0.952 [0.905, 0.979] |
| tier A | + repeat detectors | 0.944 [0.880, 0.981] | 0.883 [0.774, 0.956] | 0.786 [0.538, 0.922] | 0.951 [0.893, 0.983] |
| tier A+B (67, 31 cl.) | comp_len_struct | 0.905 [0.821, 0.962] | 0.840 [0.727, 0.924] | 0.771 [0.540, 0.906] | 0.934 [0.869, 0.976] |

Over 20 random fold assignments the mean (min-max) for tier A comp_len_struct is 0.955 (0.920 to 0.968) vs all hard negatives and 0.898 (0.849 to 0.922) vs ST-rich; that range shows split variability only.

**Leave-clade-out (comp_len_struct, tier A; composition-only in brackets)**

| Held-out clade | Positives (clusters) | vs all hard negatives | vs ST-rich | vs hard negatives with repeats |
|---|---|---|---|---|
| Saccharomycotina | 21 (10) | 0.912 [0.798, 0.986] (n=28) ([0.920, 0.811-0.987]) | 0.892 [0.756, 0.985] (n=23) ([0.907, 0.780-0.990]) | 0.757 [0.450, 0.959] (n=10) |
| Taphrinomycotina | 8 (5) | 0.975 (n=5; no interval, n too small) | 0.938 (n=2) | none |
| Eurotiomycetes | 4 (2) | 0.966 (n=58; n too small) | 0.895 (n=19) | 1.0 (n=2) |
| Sordariomycetes | 2 (2) | 0.976 (n=42) | 0.938 (n=16) | 1.0 (n=2) |

Only the Saccharomycotina fold supports a number; there the model was trained on 14 non-yeast-clade positives (8 of them S. pombe proteins; 9 clusters in all). The tier A+B version of the same fold is clearly worse (0.735 [0.578, 0.867] vs all hard negatives; 0.695 vs ST-rich), because tier B is mostly Hyr/Iff/Hyphal_reg_CWP family-inference proteins with short Ser/Thr stretches, which is evidence that tier B is not a clean part of the class.

**Single features, tier A vs ST-rich hard negatives (no training; raw AUC, higher value = more positive-like)** (`single_feature_auc.tsv`): detector-14 z_seq 0.846 [0.704, 0.942]; Thr fraction 0.837; S+T+P fraction 0.830 [0.687, 0.930]; length 0.804; coverage of detector 14 0.796; S/T+N+Q+G+P 0.787; Cys count 0.665; `max_ST30` 0.438 and longest S/T run 0.499 (these two do not separate the class from ST-rich negatives at all, by construction: the negatives are chosen for having S/T windows). Asn, Gln and Gly fractions are lower in positives (AUC 0.27 to 0.38).

**Existing repeat call on the same proteins** (`repeat_call_vs_features.tsv`; Wilson intervals)

| Rule | Positives called, tier A (35) | Hard negatives called (172) | ST-rich called (68) | With UniProt repeats called (19) | Background called (389) |
|---|---|---|---|---|---|
| Repeat call (as shipped) | 18 = 0.514 [0.356, 0.670] | 6 = 0.035 | 3 = 0.044 | 6 = 0.316 | 8 = 0.021 |
| Repeat call or z_seq 3.5 (exploratory) | 27 = 0.771 [0.610, 0.879] | 25 = 0.145 | 13 = 0.191 | 16 = 0.842 | 39 = 0.100 |

For tier A+B the shipped call finds 25/67 = 0.373 [0.267, 0.493]. At the same 6 hard-negative false hits as the shipped call, the feature model (composition + structure, leave-cluster-out seed 0) calls 24/35 tier-A positives = 0.686 [0.520, 0.814] and 44/67 of A+B = 0.657 [0.537, 0.759] (`matched_operating_point.tsv`). Caveats: threshold chosen on the same out-of-fold scores; scores pooled over 5 models; intervals overlap. It is a hint that composition carries more than the period detectors for this class, not a measurement of improvement.

**Where the repeat call stands on the truth tables.** `zseq_rule_scan.tsv` is reproduced (`repeat_rule_reproduction.tsv`, truth_v3): at z_seq 3.5, S288C sensitivity 0.739 (17/23) but 23 of 297 false hits (19 added: FPR 0.077 vs 0.013), and AGA1 is gained only at that threshold together with MSB2, HKR1, EGT2, SED1, HPF1 and DAN4. C. albicans at z_seq 3.5: 0.769 and 33 false hits (28 added).

## 5. Where AGA1 sits (`aga1_vs_named.tsv`)

| Protein | Class | Length | S+T | S/T/N/Q/G/P | Longest S/T run | N-term before S/T window | z_seq | Repeat call | Leave-clade-out score; yeast hard negatives scoring higher (of 27 or 28) |
|---|---|---|---|---|---|---|---|---|---|
| **AGA1** | adhesin E1 (tier A) | 725 | 0.552 | 0.644 | 132 | 62 | 3.74 | no | 4.37; 2 |
| SED1 | hard neg N1 (repeat_non_adhesin) | 338 | 0.417 | 0.610 | 13 | 3 | 3.69 | no | 10.61; 0 of 27 |
| DAN4 | hard neg N2 | 1161 | 0.494 | 0.606 | 154 | 104 | 4.35 | no | 7.62; 1 |
| MSB2 | hard neg N1 (mucin_sensor) | 1306 | 0.417 | 0.590 | 39 | 68 | 4.78 | no | -0.40; 5 |
| HKR1 | hard neg N1 | 1802 | 0.361 | 0.522 | 11 | 22 | 6.07 | no | 1.49; 2 |
| EGT2 | hard neg N1 | 1041 | 0.379 | 0.524 | 16 | 292 | 5.39 | no | 0.76; 4 |
| HPF1 | adhesin E3 (tier A) | 967 | 0.424 | 0.597 | 77 | 15 | 4.16 | no | 1.87; 2 |
| FLO11 | adhesin E1 | 1367 | 0.501 | 0.680 | 205 | 194 | 3.31 | no | 6.64; 2 |
| FLO1 | adhesin E1 | 1537 | 0.411 | 0.586 | 28 | 260 | 17.80 | yes | 6.56; 2 |
| FLO5 / FLO9 / FLO10 | adhesin E1 | 1075 / 1322 / 1169 | 0.40 to 0.41 | 0.54 to 0.58 | 26 to 39 | 258 to 291 | 10.5 / 15.8 / 4.9 | yes | 3.1 to 6.1; 2 |
| SAG1 | adhesin E1 | 650 | 0.299 | 0.471 | 6 | 324 | 2.22 | no | -1.93; 12 |

Reading: AGA1 has the highest S+T fraction of the group (0.552) and one of the longest uninterrupted S/T runs; a feature model trained without any Saccharomycotina protein ranks it above 26 of 28 yeast hard negatives, but SED1 and DAN4 are above it, and the repeat detectors see no period in AGA1, SED1 or MSB2 (the z_seq of 3.7 to 4.8 is the same range as the hard negatives). Its Cys count (11) is within the hard negatives' range. This supports the owner's reading of #121: AGA1 is a member of a composition-defined class, not a repeat-defined one, and the hard negatives SED1 and DAN4 are its nearest lookalikes. Two cautions: SED1 and DAN4 are curated N1 and N2, not proven non-adhesins in every context; and AGA1 is a single protein.

## 6. Narrow-class slices (`narrow_slices.tsv`; reported, not decided)

Rule (provisional): strict = at least 5 positives in each of at least 2 clades; weak = at least 5 in total over at least 2 clades. Clade = fine clade (Saccharomycotina as one unit). Confusing negatives = hard negatives of the same clades; the last column counts those scoring at or above the slice's median positive in the leave-cluster-out run.

| Narrow class | Positives (tier A) | Clades; families | Clusters | Existing Pfam narrow call (family table) | Repeat call | Hard negatives in same clades (ST-rich; scoring at or above slice median) | Strict | Weak |
|---|---|---|---|---|---|---|---|---|
| Als-like (ALS, SAG1) | 13 (11) | Saccharomycotina only; Debaryomycetaceae 9, Metschnikowiaceae 3, Saccharomycetaceae 1 | 4 | Candida_ALS_N, Candida_ALS: active, smoke, narrow_ok yes | 7/13 | 28 (23; 2) | no | no |
| Hyphal_reg_CWP-type (PF11765 members not named Hyr/Iff) | 13 (1) | Saccharomycotina; Metschnikowiaceae 7, Saccharomycetaceae 6 | 7 | Hyphal_reg_CWP: not active, smoke, narrow_ok no | 3/13 | 28 (23; 3) | no | no |
| Hyr/Iff/Rbr-like | 12 (1) | Saccharomycotina; Debaryomycetaceae 11 | 3 | Hyphal_reg_CWP (not active); Hyr1 (active, untested, narrow_ok no) | 1/12 | 28 (23; 3) | no | no |
| Flo1-type (FLO1/5/9/10) | 5 (4) | Saccharomycotina; Saccharomycetaceae 4, Debaryomycetaceae 1 | 2 | Flocculin, Flocculin_t3, PA14: active, smoke, narrow_ok yes | 4/5 | 28 (23; 3) | no | no |
| Pfl-like (S. pombe) | 8 (6) | Taphrinomycotina only | 4 | none at 50 percent of members | 4/8 | 5 (2; 0) | no | no |
| SOWgp (Coccidioides) | 3 (3) | Eurotiomycetes only | 1 | none | 3/3 | 58 (19; 39) | no | no |
| Epa-like | 2 (0) | Saccharomycotina | 1 | GLEYA: active, smoke, narrow_ok yes | 1/2 | 28 | no | no |
| Hwp1/Eap1-like | 2 (2) | Saccharomycotina | 2 | Flocculin_t3 | 0/2 | 28 | no | no |
| Flocculin-repeat type (other) | 2 (0) | Saccharomycotina | 1 | Flocculin, Hyphal_reg_CWP | 1/2 | 28 | no | no |
| Flo11-type | 1 (1) | Saccharomycotina | 1 | Flo11: active, smoke, narrow_ok yes | 0/1 | 28 | no | no |
| AGA1 | 1 (1) | Saccharomycotina | 1 | none (Pfam scan hit nothing) | 0/1 | 28 (23; 2) | no | no |
| No name or domain rule (Mad1/Mad2 of Metarhizium, two Taphrinomycotina, one Eurotiomycetes) | 5 (5) | 3 clades | 5 | none | 1/5 | 105 (37; 2) | no | yes, but a mixed bag, not a class |

Answer to "could the data support a class-specific model?": **no for every candidate** under the strict rule; **no under the weak rule except the residual grab-bag**, which is not a class. The Pfam families ALS, Flo11, GLEYA and Flocculin already give narrow domain calls. The classes with no domain (AGA1, SOWgp, Pfl, Mad) have 1 to 8 positives in one clade each. One observation outside the rule: at family level (Debaryomycetaceae, Metschnikowiaceae, Saccharomycetaceae), only the Hyphal_reg_CWP-type group has at least 5 positives in each of 2 families (7 and 6, mostly tier B and family inference); family-level hold-out inside Saccharomycotina is an owner decision, not a result.

## 7. Proposed next stage (a plan; nothing started)

Order, cheapest and most informative first.

| Step | Where | What | Cost guess (a guess, not measured) |
|---|---|---|---|
| 1. Decide definitions | owner | See decisions below | none |
| 2. Add non-yeast positives and ST-rich negatives by curation | local, literature | The bottleneck: 8 of 35 tier-A positives are S. pombe (6 of them Pfl), 6 more are non-yeast (Coccidioides SOWgp, Aspergillus, Metarhizium), none Basidiomycota; Mad1/Mad2, Coccidioides SOWgp, Aspergillus candidates are the leads | days of curation |
| 3. Tool run on the 850 proteins (then on proteomes) | HPCC, CPU | SignalP 6 (module exists), a GPI predictor (owner picks: PredGPI or similar), AIUPred/IUPred, a glycosylation predictor if licence allows | about 1 to 2 CPU hours for 850 proteins, guess |
| 4. Repeat-unit profile (cheap add-on 1) | local first, HPCC for proteomes | Build unit models from tier-A arrays (existing `analysis/cocci_repeats` unit scripts); test whether the unit finds members in another clade | hours |
| 5. Orthology within a clade (cheap add-on 2) | HPCC | OrthoFinder or reciprocal best hits inside Saccharomycotina and Pezizomycotina proteome sets, to see which FLO/ALS/Hyr groups have orthologs | 4 to 12 hours on 32 cores for 30 to 50 proteomes, guess |
| 6. Frozen ESM-2 embeddings (comparator only) | HPCC | 650M model, mask-aware mean pool (the shipped pooling includes padding, issue from the model review); proteins over 1,022 aa (17 of 35 tier-A positives) need tiling or truncation, which must be decided before the run | 850 proteins on 1 GPU roughly tens of minutes; on CPU roughly 1 to 3 hours (guess). ESM-2 35M or 150M on CPU is the cheap start. ProtT5 only if ESM-2 adds something over composition features |
| 7. Re-run the baseline with step 3 to 6 features | local | Same leave-clade-out and leave-cluster-out; report against the composition baseline and the repeat call | minutes |

Gates: no fine-tuning (issue #15 requires 300 clusters; tier A has 19). Embeddings and tools are comparators against the composition baseline (leave-cluster-out tier A AUC 0.849 vs ST-rich hard negatives; leave-clade-out Saccharomycotina 0.892 to 0.907 with wide intervals). A comparator is worth reporting only if its interval clears that.

**Owner decisions the next stage needs**
1. Positive set: tier A only, or A+B? (Tier B behaves differently: lower AUC, no repeats in most.) Require a GPI feature, or not?
2. Whether Pfl (S. pombe), SOWgp and Mad proteins belong to the class or are separate lines.
3. The ST-rich negative definition (0.5 in a 30-window is mine) and whether SED1 and DAN4 stay hard negatives.
4. The clade unit for narrow classes: keep fine clades (no class qualifies) or allow family level inside Saccharomycotina.
5. Which GPI predictor and which glycosylation predictor to run; ESM-2 size; whether to tile long proteins.
6. Whether to fetch more positives from UniProt/Swiss-Prot by feature search (new accessions, not in a committed table: outside this task's remit).

## 8. What I could not establish

- **Non-yeast generality.** 21 of 35 tier-A positives and 51 of 67 A+B are Saccharomycotina. Eurotiomycetes 4, Sordariomycetes 2, Taphrinomycotina 8 (S. pombe only), Basidiomycota 0. Any claim beyond yeast is untestable from these tables.
- **A leave-clade-out result with an interval except for Saccharomycotina**, and even there training is 14 proteins from 9 clusters, 8 of which are S. pombe.
- **GPI and signal peptide as predictions.** Only UniProt annotation lines exist; these cover 35 of 67 positives and are mostly automatic evidence for the unreviewed ones.
- **Whether the hard negatives are truly non-adhesive.** "Hard negative" means curated N1/N2; SED1, DAN4, MSB2, HKR1, EGT2 sit exactly where AGA1 sits, so the label, not the model, may be the limiting factor.
- **ESM-2 or ProtT5 behaviour.** Not run (no local GPU run asked for); the plan above is a plan.
- **Orthology and repeat-unit routes.** Not tested.
- **Independence of positives.** 30 percent identity clusters may miss remote family relationships (ALS/Hyr/Iff share domains below that cutoff). Cluster bootstrap intervals do not include training variability.
- **Negative background.** 389 background proteins are assumed negatives (some may be unlabelled adhesin-like proteins); they are used only as an extra test set.
- **My definitions were fixed after I looked at the gene list and the median features, before any AUC was computed.** Tier A/B and the 0.5 ST-rich cut are therefore not blind.
- **UniProt version drift.** Sequences and feature lines were fetched on 2026-10-24; later releases may change them.

## 9. Surprises

- The repeat call reproduces the committed statuses exactly from locally re-run detectors (S288C 10/23 and 4/297; C. albicans 6/13 and 5/211; 171 of 171 per-protein calls).
- The truth label that the repeat call is scored on (UniProt at least 2 repeat features) is not the same as "adhesin": 10 of the 23 S288C label-1 proteins are curated hard negatives, so a perfect repeat detector would still call them.
- Over half of the tier-A positives (17 of 35) exceed the ESM-2 window, against 4 of 172 hard negatives: length is the single strongest nuisance variable and is also what ESM tiling decisions will touch.
- Tier-B positives (Hyr/Iff and relatives) are family-inferred "repeat" members with no UniProt repeats in most; including them lowers every AUC.
- `max_ST30` and longest S/T run do not separate the class from ST-rich negatives (AUC 0.44 and 0.50), whereas z_seq, Thr and S+T+P fractions do (0.85, 0.84, 0.83).
