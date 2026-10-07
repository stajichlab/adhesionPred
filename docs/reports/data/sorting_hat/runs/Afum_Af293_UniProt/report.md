# cellsurface_sorting_hat report

- version: 0.2.0.dev421
- proteins: 9647 (0 invalid, excluded from all modules; 0 had a trailing `*`, removed)
- taxa (taxon ID: proteins): 330879: 9647
- taxonomy file sha256: 78392a6926b01fadda65ab819b75e0685e62593fffc0efe09c70b04d26cb7f6b
- categories.yaml sha256: f97cb08508c13100cded294a41d86940d3907419f7d59d76f49bbcc278a08ade
- default gate: step1_rule@R0
- thresholds: allergen_coverage_min = 80, allergen_hit_identity_min = 35, allergen_hit_length_min = 80, allergen_identity_min = 70, antigen_percentile_max = 15

## Warnings

- **WARNING: Modules the rules need and that were not found: antigen_lookup**

## Module run states

| module | run state |
|---|---|
| allergen_homology | ok |
| pfam_adhesion | unavailable |
| pfam_allergen | unavailable |
| repeat02 | ok |
| repeat14 | ok |
| step1_rule@R0 | ok |
| tm | ok |

Unavailable step 1 variants: step1_rule@R1, step1_rule@R2, step1_ml@card

## Module calibration

| module | taxon | status | measured on call | calibration set | positives | negatives | sensitivity [95% CI] | specificity [95% CI] |
|---|---|---|---|---|---|---|---|---|
| allergen_homology | 330879 | unvalidated | - | - | - | - | not measured | not measured |
| repeat02 | 330879 | unvalidated | - | - | - | - | not measured | not measured |
| repeat14 | 330879 | unvalidated | - | - | - | - | not measured | not measured |
| step1_rule@R0 | 330879 | smoke | signal_peptide_protein | S3-Eurotiomycetes:Afum_ASPFU | 19 | 45 | 0.947 [0.738, 1.000] | 1.000 [0.908, 1.000] |
| tm | 330879 | unvalidated | - | - | - | - | not measured | not measured |

## Calls

| call | variant | called | not_called | not_assessable |
|---|---|---|---|---|
| cell_wall_adhesion_candidate | R0 | 8 | 8894 | 745 |
| cocci_specificity_rank_top15 | - | 0 | 0 | 9647 |
| iuis_allergen_homolog | - | 41 | 0 | 9606 |
| iuis_allergen_similarity | - | 104 | 9543 | 0 |
| other_not_surface | R0 | 8879 | 768 | 0 |
| other_surface_no_mechanism | R0 | 745 | 8902 | 0 |
| serodiagnostic_marker_candidate | R0 | 0 | 8894 | 753 |
| signal_peptide_protein | R0 | 753 | 8894 | 0 |
| tandem_repeat_protein | - | 23 | 9624 | 0 |
| wall_family_domain | - | 0 | 0 | 9647 |

## `other_basis` (categories left out because they were not assessable)

- wall_family_domain,cocci_specificity_rank_top15: 9624

## Known limits

1. Research use only. This is not a regulatory allergenicity assessment and not a diagnostic result. No row is supported by an IgE, antibody or T-cell measurement.
2. `signal_peptide_protein` means that SignalP calls a signal peptide (rule R0) and nothing more. Such proteins are secreted, wall-bound or GPI-anchored. They are not shown to be exposed at the cell surface, and glycosylation is not assessed. Plasma membrane mucins can be missed.
3. GPI-anchored and secreted enzymes, and non-adhesive structural wall proteins, get no finer label in version 1. There is no `cell_wall_protein` call.
4. `tandem_repeat_protein` and `wall_family_domain` are evidence. Repeat proteins include intracellular ones (ubiquitin, calmodulin, ankyrin proteins). A domain of a family that is linked to adhesion or wall function in at least one species (CFEM, Bys1, hydrophobin, Als) is not shown to mediate adhesion here. The adhesion call needs a signal peptide.
5. `cocci_specificity_rank_top15` is the top 15% of a fixed Coccidioides immitis ranking (similarity to IEDB antigens, prevalence, absence of orthologs in confounder fungi). It is not epitope prediction. The cut was set after the four anchors were seen. The ranking prints NOT CALIBRATED. `serodiagnostic_marker_candidate` adds a signal peptide. Peptide level only; glycan epitopes are not assessed.
6. `iuis_allergen_similarity` and `iuis_allergen_homolog` are sequence similarity to allergens in the WHO/IUIS fungal set (IgE binding in patients). They suggest possible IgE cross-reactivity at most. A protein with no hit is not thereby a non-allergen. WHO/IUIS lists no Coccidioides allergen.
7. Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, non-protein adhesins, moonlighting proteins and biofilm are not categories.
8. The taxon you give is recorded as given. It is not checked against the sequences.
9. Leakage: overlap between the Phase C positives and the SignalP 6 training data was not measured. The status of rule R0 does not account for it.
