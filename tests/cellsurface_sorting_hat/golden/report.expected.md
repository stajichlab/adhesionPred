# cellsurface_sorting_hat report

- proteins: 5 (1 invalid, excluded from all modules; 1 had a trailing `*`, removed)
- taxa (taxon ID: proteins): 40: 3, 41: 2
- taxonomy file sha256: 492f746e8af7ed537fe0bcf4920a583e583410d1ebf5b12cd6c1284e44d474df
- categories.yaml sha256: fd397b941d9d7eb7bb94b2b48e8afaaabda55d3342dae0e4933c3f242a1fb6d9
- default gate: step1_rule@R0
- thresholds: allergen_coverage_min = 80, allergen_hit_identity_min = 35, allergen_hit_length_min = 80, allergen_identity_min = 70, antigen_percentile_max = 15

## Invalid proteins

- BAD1: internal stop codon

## Module run states

| module | run state |
|---|---|
| allergen_homology | ok |
| antigen_lookup | ok |
| pfam_adhesion | ok |
| pfam_allergen | ok |
| pfam_hsba | ok |
| pfam_hydrophobin | ok |
| repeat02 | ok |
| repeat14 | ok |
| step1_rule@R0 | ok |

Unavailable step 1 variants: step1_rule@R1, step1_rule@R2, step1_ml@card

## Proteins that are not `ok` in a module

| module | state | proteins |
|---|---|---|
| allergen_homology | na_invalid | 1 |
| antigen_lookup | na_invalid | 1 |
| antigen_lookup | not_applicable | 2 |
| pfam_adhesion | na_invalid | 1 |
| pfam_allergen | na_invalid | 1 |
| pfam_hsba | na_invalid | 1 |
| pfam_hydrophobin | na_invalid | 1 |
| repeat02 | na_invalid | 1 |
| repeat14 | na_invalid | 1 |
| step1_rule@R0 | na_invalid | 1 |

## Module calibration

| module | taxon | status | measured on call | calibration set | positives | negatives | sensitivity [95% CI] | specificity [95% CI] |
|---|---|---|---|---|---|---|---|---|
| allergen_homology | 40 | unvalidated | - | - | - | - | not measured | not measured |
| allergen_homology | 41 | unvalidated | - | - | - | - | not measured | not measured |
| antigen_lookup | 40 | unvalidated | - | - | - | - | not measured | not measured |
| antigen_lookup | 41 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_adhesion | 40 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_adhesion | 41 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_allergen | 40 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_allergen | 41 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_hsba | 40 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_hsba | 41 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_hydrophobin | 40 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_hydrophobin | 41 | unvalidated | - | - | - | - | not measured | not measured |
| repeat02 | 40 | unvalidated | - | - | - | - | not measured | not measured |
| repeat02 | 41 | unvalidated | - | - | - | - | not measured | not measured |
| repeat14 | 40 | unvalidated | - | - | - | - | not measured | not measured |
| repeat14 | 41 | unvalidated | - | - | - | - | not measured | not measured |
| step1_rule@R0 | 40 | estimated | - | - | - | - | not measured | not measured |
| step1_rule@R0 | 41 | unvalidated | - | - | - | - | not measured | not measured |

## Calls

| call | variant | called | not_called | not_assessable |
|---|---|---|---|---|
| cell_wall_adhesion_candidate | R0 | 2 | 2 | 1 |
| cocci_specificity_rank_top15 | - | 2 | 0 | 3 |
| hsba_domain | - | 0 | 4 | 1 |
| hydrophobin_domain | - | 0 | 4 | 1 |
| iuis_allergen_homolog | - | 1 | 3 | 1 |
| iuis_allergen_similarity | - | 2 | 2 | 1 |
| other_not_surface | R0 | 1 | 3 | 1 |
| other_surface_no_mechanism | R0 | 1 | 3 | 1 |
| serodiagnostic_marker_candidate | R0 | 2 | 1 | 2 |
| signal_peptide_protein | R0 | 3 | 1 | 1 |
| tandem_repeat_protein | - | 2 | 2 | 1 |
| wall_family_domain | - | 0 | 4 | 1 |

## `other_basis` (categories left out because they were not assessable)

- cocci_specificity_rank_top15: 2

## Known limits

1. Research use only. This is not a regulatory allergenicity assessment and not a diagnostic result. No row is supported by an IgE, antibody or T-cell measurement.
2. `signal_peptide_protein` means that SignalP calls a signal peptide (rule R0) and nothing more. Such proteins are secreted, wall-bound or GPI-anchored. They are not shown to be exposed at the cell surface, and glycosylation is not assessed. Plasma membrane mucins can be missed.
3. GPI-anchored and secreted enzymes, and non-adhesive structural wall proteins, get no finer label in version 1. There is no `cell_wall_protein` call.
4. `tandem_repeat_protein` and `wall_family_domain` are evidence. Repeat proteins include intracellular ones (ubiquitin, calmodulin, ankyrin proteins). A domain of a family that is linked to adhesion or wall function in at least one species (CFEM, Bys1, Als) is not shown to mediate adhesion here. The adhesion call needs a signal peptide. `hydrophobin_domain` (hydrophobin Pfam models, category 2c) and `hsba_domain` (HsbA, PF12296, which may not be a hydrophobin) are separate evidence calls. They do not set `wall_family_domain`.
5. `cocci_specificity_rank_top15` is the top 15% of a fixed Coccidioides immitis ranking (similarity to IEDB antigens, prevalence, absence of orthologs in confounder fungi). It is not epitope prediction. The cut was set after the four anchors were seen. The ranking prints NOT CALIBRATED. `serodiagnostic_marker_candidate` adds a signal peptide. Peptide level only; glycan epitopes are not assessed.
6. `iuis_allergen_similarity` and `iuis_allergen_homolog` are sequence similarity to allergens in the WHO/IUIS fungal set (IgE binding in patients). They suggest possible IgE cross-reactivity at most. A protein with no hit is not thereby a non-allergen. WHO/IUIS lists no Coccidioides allergen.
7. Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, non-protein adhesins, moonlighting proteins and biofilm are not categories.
8. The taxon you give is recorded as given. It is not checked against the sequences.
9. Leakage: overlap between the Phase C positives and the SignalP 6 training data was not measured. The status of rule R0 does not account for it.
