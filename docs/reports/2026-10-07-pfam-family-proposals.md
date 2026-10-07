# Pfam families: proposals for sign-off (decision C4)

*2026-10-07. Proposals only. Nothing is set to `active`, and `data/sorting_hat/family_table.tsv` is not changed. The owner decides each family. Per-hit data: `docs/reports/data/sorting_hat/specificity/pfam_proposals.tsv` (122 family-protein hits in 6 proteomes; script `analysis/sorting_hat_run/pfam_proposals.py`; the rules are at the top of the script). Research use only. A Pfam hit is domain evidence. It does not show adhesion.*

## What "active" means (owner decision 4)

An active family reports its hit as evidence. The family has a `specificity_note` and, where the review shows a need, a `second_condition`. The status stays `unvalidated` until a truth table exists. No family becomes an adhesin call.

## Class counts proposed per family

Classes are `proposed`, not reviewed. `needs_expert` means the table facts do not decide.

| family | hits | known_member (draft list) | uncharacterised_true_member | receptor_like | false_domain_hit | needs_expert |
|---|---|---|---|---|---|---|
| PF00399 PIR | 10 | 5 | 0 | 0 | 0 | 5 |
| PF01185 Hydrophobin | 18 | 6 | 11 | 0 | 0 | 1 |
| PF04681 Bys1 | 17 | 1 | 0 | 0 | 0 | 16 |
| PF05730 CFEM | 27 | 3 | 12 | 3 | 0 | 9 |
| PF07691 PA14 | 28 | 4 | 0 | 0 | 17 | 7 |
| PF10182 Flo11 | 2 | 1 | 0 | 0 | 0 | 1 |
| PF10528 GLEYA | 6 | 4 | 0 | 0 | 0 | 2 |
| PF11765 Hyphal_reg_CWP | 2 | 0 | 0 | 0 | 0 | 2 |
| PF11766 Candida_ALS_N | 1 | 0 | 0 | 0 | 0 | 1 |
| PF25312 Allergen_Asp_f_4 | 9 | 1 | 5 | 0 | 0 | 3 |
| PF28987 DewD | 2 | 0 | 1 | 0 | 0 | 1 |
| PF06766, PF22354, PF05792, PF16541 | 0 | - | - | - | - | - |

Four families have no hit in any of the six proteomes (Hydrophobin_2, Eas, Candida_ALS, AltA1). Nothing can be reviewed for them yet.

## Proposal per family

| family | proposed second condition | evidence from the tables | readiness |
|---|---|---|---|
| PA14 (PF07691) | keep `signal_peptide` | All 17 *A. fumigatus* hits (Af293 6, A1163 6, W72310 5) have no R0 signal peptide call. The 12 in Af293 and A1163 are annotated beta-glucosidases and are the proposed `false_domain_hit` rows. The 5 W72310 hits are "hypothetical protein" and stay `needs_expert`. The 5 *C. neoformans* H99 hits are annotated beta-glucosidases (SignalP not run on H99), so they are also proposed `false_domain_hit`. In S288C, 4 of 6 hits have a signal peptide and are the draft flocculin members. | Ready after SignalP on H99 and expert check of the 7 `needs_expert` rows (S288C TDA8, YHR213W and others). |
| CFEM (PF05730) | keep `no_tm` | 4 of 27 hits have a TMHMM helix (W72310 2, RS 1, H99 1). Three are proposed `receptor_like`. The 3 known Af293 CfmA-C proteins have 8 or 9 Cys and no helix. In RS, the CFEM hits include Ag2/PRA, PRA2 and "proline-rich antigen 7". Two hits with 7 Cys and 680 to 760 aa ("extracellular serine-threonine rich protein", Af293 Q4WYI0 and RS XP_001244031.1) need an expert. One RS hit is annotated "chitinase 3". | Ready after SignalP on H99 and expert check of the 9 `needs_expert` rows. Note for the panel: when CFEM is active, Ag2/PRA becomes a `wall_family_domain` call. |
| Hydrophobin (PF01185) | none | All 18 hits (Af293 6, A1163 5, W72310 6, RS 1) have 8 Cys and a signal peptide. RodD (hydrophobin-like, in the draft list) has no Pfam hit, so the model misses at least one known hydrophobin-like protein. One A1163 hit has 1 helix with a signal peptide. | Ready. The RodD miss belongs in the specificity note. |
| DewD (PF28987) | none | 2 hits (W72310, RS), both with 8 Cys and a signal peptide. The RS hit also has 1 helix. | Too few hits to say. Review with the hydrophobin family. |
| Bys1 (PF04681) | none proposed | 17 hits, but only CalA is a known member. 16 are `needs_expert`. Annotations such as "BYS1 domain protein" come from the domain, so they are not independent. The family table names calB and calC as the control, and they are not identified in the tables. | Not ready. Needs the paralog control and an expert. |
| Flo11 (PF10182), GLEYA (PF10528), PIR (PF00399) | none proposed | Reviewed in S288C only. The draft lists cover FLO11, FLO1, FLO5, FLO9, FLO10 and PIR1, PIR3, PIR5, HSP150, CIS3. Non-members: PIR TIR1, TIR2, CWP1, CWP2, ANS1 (all have a signal peptide); GLEYA and PA14 TDA8 and YHR213W (no signal peptide). | Not ready. S288C alone is thin, and the PIR non-members need an expert. |
| Hyphal_reg_CWP (PF11765), Candida_ALS_N (PF11766), Candida_ALS (PF05792) | none proposed | The family table calls these Candida-specific, but no Candida proteome was searched. S288C has 2 and 1 hits (CSS1, HPF1, SAG1) and no ALS repeat hit. | Not ready. Search *C. albicans* SC5314 first. |
| Allergen_Asp_f_4 (PF25312), AltA1 (PF16541) | none proposed | Asp f 4: 9 hits in 3 *A. fumigatus* proteomes. 5 are annotated "allergen-like". AltA1: no hit in any proteome (no *Alternaria* proteome searched). | Asp f 4 reviewable. Annotation by similarity is not IgE evidence. |

## What is needed before sign-off

1. SignalP on *C. neoformans* H99 (the CFEM and PA14 rows there have no signal peptide call). One short GPU job.
2. A *C. albicans* SC5314 run (Pfam, TMHMM, SignalP) for the Candida-specific families. The FASTA is already in `_workdir/step1_compare/downloads/`.
3. An *Alternaria* proteome for AltA1, or leave that family out of the first sign-off.
4. An expert pass on the `needs_expert` rows (48 of 122), starting with CFEM and PA14.
5. Reviewed member lists (agent task 04, issue #69) to replace the draft lists. The draft lists are not reviewed.

Proposed first sign-off, if the owner agrees after steps 1 and 4: PA14 with `signal_peptide`, CFEM with `no_tm`, and Hydrophobin, each with a `specificity_note` that states the counts above.
