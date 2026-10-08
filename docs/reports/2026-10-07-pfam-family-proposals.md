# Pfam families: proposals for sign-off (decision C4)

*2026-10-07. Proposals only. Nothing is set to `active`, and `data/sorting_hat/family_table.tsv` is not changed. The owner decides each family. Per-hit data: `docs/reports/data/sorting_hat/specificity/pfam_proposals.tsv` (159 family-protein hits in 8 proteomes, including *C. neoformans* H99 with SignalP and *C. albicans* SC5314; script `analysis/sorting_hat_run/pfam_proposals.py`; the rules are at the top of the script). Research use only. A Pfam hit is domain evidence. It does not show adhesion.*

## What "active" means (owner decision 4)

An active family reports its hit as evidence. The family has a `specificity_note` and, where the review shows a need, a `second_condition`. The status stays `unvalidated` until a truth table exists. No family becomes an adhesin call.

## Class counts proposed per family

Classes are `proposed`, not reviewed. `needs_expert` means the table facts do not decide.

| family | hits | known_member (draft list) | uncharacterised_true_member | receptor_like | false_domain_hit | needs_expert |
|---|---|---|---|---|---|---|
| PF00399 PIR | 10 | 5 | 0 | 0 | 0 | 5 |
| PF01185 Hydrophobin | 18 | 6 | 11 | 0 | 0 | 1 |
| PF04681 Bys1 | 17 | 1 | 0 | 0 | 0 | 16 |
| PF05730 CFEM | 33 | 3 | 19 | 4 | 0 | 7 |
| PF05792 Candida_ALS | 9 | 8 | 0 | 0 | 0 | 1 |
| PF07691 PA14 | 29 | 4 | 0 | 0 | 17 | 8 |
| PF10182 Flo11 | 3 | 1 | 0 | 0 | 0 | 2 |
| PF10528 GLEYA | 6 | 4 | 0 | 0 | 0 | 2 |
| PF11765 Hyphal_reg_CWP | 14 | 10 | 0 | 0 | 0 | 4 |
| PF11766 Candida_ALS_N | 9 | 8 | 0 | 0 | 0 | 1 |
| PF25312 Allergen_Asp_f_4 | 9 | 1 | 5 | 0 | 0 | 3 |
| PF28987 DewD | 2 | 0 | 1 | 0 | 0 | 1 |
| PF06766, PF22354, PF16541 | 0 | - | - | - | - | - |

Three families have no hit in any of the eight proteomes (Hydrophobin_2, Eas, AltA1). Nothing can be reviewed for them yet.

## Proposal per family

| family | proposed second condition | evidence from the tables | readiness |
|---|---|---|---|
| PA14 (PF07691) | keep `signal_peptide` | None of the 23 non-yeast PA14 hits has an R0 signal peptide: *A. fumigatus* Af293 6, A1163 6, W72310 5, *C. neoformans* H99 5 (SignalP now run on H99), *C. albicans* 1 (C3_01820W_A, unknown function). The 17 hits in Af293, A1163 and H99 are annotated beta-glucosidases and are proposed `false_domain_hit`. The 5 W72310 hits are "hypothetical protein" and stay `needs_expert`. In S288C, 4 of 6 hits have a signal peptide and are the draft flocculin members. | Ready after an expert check of the 8 `needs_expert` rows. |
| CFEM (PF05730) | keep `no_tm` | 5 of 33 hits have a TMHMM helix (W72310 2, RS 1, H99 1, *C. albicans* 1). Four are proposed `receptor_like` (W72310 2, RS 1, H99 1: two or more helices, or one helix without a signal peptide). The *C. albicans* hit with one helix and a signal peptide (CSA2) stays `needs_expert`. The 3 known Af293 CfmA-C proteins have 8 or 9 Cys and no helix. The 6 *C. albicans* hits (PGA7, RBT5, PGA10, CSA2, CSA1, SSR1 by CGD name) all have 8 Cys and a signal peptide. In RS, the CFEM hits include Ag2/PRA, PRA2 and "proline-rich antigen 7". Two hits with 7 Cys and 680 to 760 aa (Af293 Q4WYI0, RS XP_001244031.1) need an expert. One RS hit is annotated "chitinase 3". | Ready after an expert check of the 7 `needs_expert` rows. Note for the panel: when CFEM is active, Ag2/PRA becomes a `wall_family_domain` call. |
| Hydrophobin (PF01185) | none | All 18 hits (Af293 6, A1163 5, W72310 6, RS 1) have 8 Cys and a signal peptide. RodD (hydrophobin-like, in the draft list) has no Pfam hit, so the model misses at least one known hydrophobin-like protein. One A1163 hit has 1 helix with a signal peptide. | Ready. The RodD miss belongs in the specificity note. |
| DewD (PF28987) | none | 2 hits (W72310, RS), both with 8 Cys and a signal peptide. The RS hit also has 1 helix. | Too few hits to say. Review with the hydrophobin family. |
| Bys1 (PF04681) | none proposed | 17 hits, but only CalA is a known member. 16 are `needs_expert`. Annotations such as "BYS1 domain protein" come from the domain, so they are not independent. The family table names calB and calC as the control, and they are not identified in the tables. | Not ready. Needs the paralog control and an expert. |
| Flo11 (PF10182), GLEYA (PF10528), PIR (PF00399) | none proposed | Reviewed in S288C, and Flo11 also in *C. albicans* (one hit, RBT1). The draft lists cover FLO11, FLO1, FLO5, FLO9, FLO10 and PIR1, PIR3, PIR5, HSP150, CIS3. Non-members: PIR TIR1, TIR2, CWP1, CWP2, ANS1 (all have a signal peptide); GLEYA and PA14 TDA8 and YHR213W (no signal peptide). | Not ready. The PIR non-members need an expert. |
| Candida_ALS_N (PF11766), Candida_ALS (PF05792) | none proposed | *C. albicans* SC5314 was searched on 2026-10-07. Candida_ALS_N has 8 hits in *C. albicans*: ALS1, ALS2, ALS3, ALS4, ALS5, ALS6, ALS7 and ALS9 by CGD gene name (the draft member list), each with a signal peptide and 8 Cys. A CGD lookup of ALS8 returns the ALS3 locus (CR_07070C_A), so ALS8 has no separate hit. Candida_ALS (the repeat region) has the same 8 plus C3_00580W_A (CGD gene FLO9, "adhesin-like cell wall mannoprotein"). S288C adds SAG1 to Candida_ALS_N. | Reviewable. Sensitivity is high on this list, but the list is from gene names. |
| Hyphal_reg_CWP (PF11765) | none proposed | 12 hits in *C. albicans* (plus CSS1 and HPF1 in S288C). Ten are HYR or IFF family proteins by CGD gene name (HYR1, HYR3, HYR4, IFF3, IFF4, IFF5, IFF6, IFF8, IFF9, IFF11) and are the draft members. The other two are C3_00580W_A (FLO9) and C7_03290C_A (RBR3). **HWP1 (C4_03570W_A), HWP2 (C4_03510C_A) and EAP1 (C2_09530W_A) are in the proteome and have no hit from any of the 15 models.** The family table calls PF11765 an "Hwp1-like wall protein". In *C. albicans* its hits are the HYR/IFF group, not HWP1. RBT1 (C4_03520C_A) hits Flo11 (PF10182) and not PF11765. | The class text in the family table needs an expert check before sign-off. |
| Allergen_Asp_f_4 (PF25312), AltA1 (PF16541) | none proposed | Asp f 4: 9 hits in 3 *A. fumigatus* proteomes. 5 are annotated "allergen-like". AltA1: no hit in any proteome (no *Alternaria* proteome searched). | Asp f 4 reviewable. Annotation by similarity is not IgE evidence. |

## What is needed before sign-off

1. ~~SignalP on *C. neoformans* H99.~~ Done 2026-10-07 (job 29605895).
2. ~~A *C. albicans* SC5314 run.~~ Done 2026-10-07 (the Phase C FASTA has 38 duplicate IDs with identical sequences; a deduplicated copy of 6,212 proteins was used; provenance `data/sorting_hat/provenance/Calb_SC5314.json`).
3. An *Alternaria* proteome for AltA1, or leave that family out of the first sign-off.
4. An expert pass on the `needs_expert` rows (51 of 159), starting with CFEM and PA14.
5. Reviewed member lists (agent task 04, issue #69) to replace the draft lists. The draft lists are not reviewed.
6. An expert check of the PF11765 class text (the HWP1 miss above).

Proposed first sign-off, if the owner agrees after step 4: PA14 with `signal_peptide`, CFEM with `no_tm`, and Hydrophobin, each with a `specificity_note` that states the counts above.
