# cellsurface_sorting_hat report

- version: 0+unknown
- proteins: 6039 (8 invalid, excluded from all modules; 6037 had a trailing `*`, removed)
- taxa (taxon ID: proteins): 559292: 6039
- taxonomy file sha256: 78392a6926b01fadda65ab819b75e0685e62593fffc0efe09c70b04d26cb7f6b
- categories.yaml sha256: 68a35bad9dedcb8c34e0368153817edeb1ef489da0e12ed9d9243d4e52cc3721
- default gate: step1_rule@R0
- thresholds: allergen_coverage_min = 80, allergen_hit_identity_min = 35, allergen_hit_length_min = 80, allergen_identity_min = 70, antigen_percentile_max = 15

## What the statuses mean

- `unvalidated`: no measurement exists for this module, call or taxon.
- `smoke`: measured on too little truth for a number (the leakage cap also gives `smoke`).
- `estimated`: measured with enough truth for a number (at least 20 positives, 20 negatives and 20 clusters, and a narrow interval).
- A status holds for one species and taxon: it counts only where its file is installed and current.
- The status of a composite call (a call built from other calls) is the weakest status of the leaf calls that decided that record's value.

## Invalid proteins

- YAR061W: internal stop codon
- YDR134C: internal stop codon
- YER109C: internal stop codon
- YFL056C: internal stop codon
- YIL167W: internal stop codon
- YIR043C: internal stop codon
- YOL153C: internal stop codon
- YOR031W: internal stop codon

## Module run states

| module | run state |
|---|---|
| allergen_homology | ok |
| antigen_lookup | ok |
| pfam_adhesion | ok |
| pfam_allergen | unavailable |
| pfam_hsba | ok |
| pfam_hydrophobin | ok |
| repeat02 | ok |
| repeat14 | ok |
| step1_rule@R0 | ok |
| tm | ok |

Unavailable step 1 variants: step1_rule@R1, step1_rule@R2, step1_ml@card

## Proteins that are not `ok` in a module

| module | state | proteins |
|---|---|---|
| allergen_homology | na_invalid | 8 |
| antigen_lookup | na_invalid | 8 |
| antigen_lookup | not_applicable | 6031 |
| pfam_adhesion | na_invalid | 8 |
| pfam_hsba | na_invalid | 8 |
| pfam_hydrophobin | na_invalid | 8 |
| repeat02 | na_invalid | 8 |
| repeat14 | na_invalid | 8 |
| step1_rule@R0 | na_invalid | 8 |
| tm | na_invalid | 8 |

## Module calibration

| module | taxon | status | measured on call | calibration set | positives | negatives | sensitivity [95% CI] | specificity [95% CI] |
|---|---|---|---|---|---|---|---|---|
| allergen_homology | 559292 | unvalidated | - | - | - | - | not measured | not measured |
| antigen_lookup | 559292 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_adhesion | 559292 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_hsba | 559292 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_hydrophobin | 559292 | unvalidated | - | - | - | - | not measured | not measured |
| repeat02 | 559292 | unvalidated | - | - | - | - | not measured | not measured |
| repeat14 | 559292 | unvalidated | - | - | - | - | not measured | not measured |
| step1_rule@R0 | 559292 | smoke | signal_peptide_protein | S1:Scer_SGD | 79 | 3783 | 0.848 [0.734, 0.943] | 0.966 [0.959, 0.972] |
| tm | 559292 | unvalidated | - | - | - | - | not measured | not measured |

## Calls

| call | variant | called | not_called | not_assessable | not_applicable |
|---|---|---|---|---|---|
| cell_wall_adhesion_candidate | R0 | 15 | 6016 | 8 | 0 |
| cocci_specificity_rank_top15 | - | 0 | 0 | 8 | 6031 |
| hsba_domain | - | 0 | 6031 | 8 | 0 |
| hydrophobin_domain | - | 0 | 6031 | 8 | 0 |
| iuis_allergen_homolog | - | 15 | 0 | 6024 | 0 |
| iuis_allergen_similarity | - | 84 | 5947 | 8 | 0 |
| other_not_surface | R0 | 5727 | 304 | 8 | 0 |
| other_surface_no_mechanism | R0 | 280 | 5751 | 8 | 0 |
| serodiagnostic_marker_candidate | R0 | 0 | 6031 | 8 | 0 |
| signal_peptide_protein | R0 | 296 | 5735 | 8 | 0 |
| surface_attachment_candidate | R0 | 16 | 6015 | 8 | 0 |
| tandem_repeat_protein | - | 19 | 6012 | 8 | 0 |
| wall_family_domain | - | 10 | 6021 | 8 | 0 |
| wall_family_narrow_domain | - | 9 | 6022 | 8 | 0 |

## Primary label

One label per protein for the default gate, from the calls above. The first match wins. The calls stay independent flags; the label is a summary for a reader.

| primary label | proteins |
|---|---|
| cell_wall_adhesion_candidate | 15 |
| surface_attachment_candidate | 1 |
| other_surface_no_mechanism | 280 |
| other_not_surface | 5727 |
| mechanism_without_class | 8 |
| not_assessable | 8 |

## Surface attachment: which evidence held

A protein can hold more than one. Hydrophobin and HsbA are surface-active proteins: they adsorb to surfaces. None of these labels shows adhesion.

| evidence | proteins |
|---|---|
| proteins called | 16 |
| tandem repeat | 13 |
| wall family domain | 7 |
| hydrophobin | 0 |
| HsbA | 0 |

## Known limits

1. Research use only. This is not a regulatory allergenicity assessment and not a diagnostic result. No row is supported by an IgE, antibody or T-cell measurement.
2. `signal_peptide_protein` means that SignalP calls a signal peptide (rule R0) and nothing more. Such proteins are secreted, wall-bound or GPI-anchored. They are not shown to be exposed at the cell surface, and glycosylation is not assessed. Plasma membrane mucins can be missed.
3. GPI-anchored and secreted enzymes, and non-adhesive structural wall proteins, get no finer label in version 1. There is no `cell_wall_protein` call.
4. `tandem_repeat_protein`, `wall_family_domain` and `wall_family_narrow_domain` are evidence. Repeat proteins include intracellular ones (ubiquitin, calmodulin, ankyrin proteins). A domain of a family that is linked to adhesion or wall function in at least one species (CFEM, ALS, Flo11) is not shown to mediate adhesion here. The adhesion call needs a signal peptide. `wall_family_domain` reads every family that the owner enabled (`active`; not a test result). `wall_family_narrow_domain` reads only the families with a specificity tier of smoke or better that the owner signed with narrow_ok=yes. `cell_wall_adhesion_candidate` reads the narrow call; `surface_attachment_candidate` reads the broad one. `hydrophobin_domain` (hydrophobin Pfam models, category 2c) and `hsba_domain` (HsbA, PF12296, which may not be a hydrophobin) are separate evidence calls. They do not set `wall_family_domain`.
5. `cocci_specificity_rank_top15` is the top 15% of a fixed Coccidioides immitis ranking (similarity to IEDB antigens, prevalence, absence of orthologs in confounder fungi). It is not epitope prediction. The cut was set after the four anchors were seen. The ranking prints NOT CALIBRATED. `serodiagnostic_marker_candidate` adds a signal peptide. Peptide level only; glycan epitopes are not assessed.
6. `iuis_allergen_similarity` and `iuis_allergen_homolog` are sequence similarity to allergens in the WHO/IUIS fungal set (IgE binding in patients). They suggest possible IgE cross-reactivity at most. A protein with no hit is not thereby a non-allergen. WHO/IUIS lists no Coccidioides allergen.
7. Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, non-protein adhesins, moonlighting proteins and biofilm are not categories.
8. The taxon you give is recorded as given. It is not checked against the sequences.
9. Leakage: overlap between the Phase C positives and the SignalP 6 training data was not measured. The status of rule R0 does not account for it. This is the owner's recorded exception: the R0 status entry carries `leakage_exception`, and its text is in the basis.
