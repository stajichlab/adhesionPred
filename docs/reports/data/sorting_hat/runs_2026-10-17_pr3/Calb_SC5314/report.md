# cellsurface_sorting_hat report

- version: 0+unknown
- proteins: 6212 (0 invalid, excluded from all modules; 0 had a trailing `*`, removed)
- taxa (taxon ID: proteins): 237561: 6212
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

## Notes

- Module cys8_pattern records no digest of the step1_rule@R0 table it was built against (condition sp_module; an older record). The tables cannot be compared.
- Module hydrophobin_relaxed records no digest of the step1_rule@R0 table it was built against (condition sp_module; an older record). The tables cannot be compared.

## Module run states

| module | run state |
|---|---|
| allergen_homology | ok |
| antigen_lookup | ok |
| cys8_pattern | ok |
| hydrophobin_relaxed | ok |
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
| antigen_lookup | not_applicable | 6212 |

## Module calibration

| module | taxon | status | measured on call | calibration set | positives | negatives | sensitivity [95% CI] | specificity [95% CI] |
|---|---|---|---|---|---|---|---|---|
| allergen_homology | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| antigen_lookup | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| cys8_pattern | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| hydrophobin_relaxed | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_adhesion | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_hsba | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| pfam_hydrophobin | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| repeat02 | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| repeat14 | 237561 | unvalidated | - | - | - | - | not measured | not measured |
| step1_rule@R0 | 237561 | smoke | signal_peptide_protein | S1:Calb_CGD | 153 | 438 | 0.477 [0.356, 0.591] | 0.952 [0.926, 0.969] |
| tm | 237561 | unvalidated | - | - | - | - | not measured | not measured |

## Call calibration

Calls that read several modules, measured as a whole (`status/calls/`). A valid file sets the status of its call. Statuses of composite calls (for example `cell_wall_adhesion_candidate`) are derived: the weakest status of the measured leaf calls that decided them. A stale file is listed with its reason and is not used.

| call | variant | taxon | status | calibration set | positives | negatives | sensitivity [95% CI] | specificity [95% CI] | validity |
|---|---|---|---|---|---|---|---|---|---|
| tandem_repeat_protein | - | 5476 | smoke | repeat_motif_definition_Calb | 13 | 211 | 0.462 [0.000, 0.788] | 0.976 [0.938, 0.996] | valid |

## Calls

| call | variant | called | not_called | not_assessable | not_applicable |
|---|---|---|---|---|---|
| cell_wall_adhesion_candidate | R0 | 29 | 6183 | 0 | 0 |
| cocci_specificity_rank_top15 | - | 0 | 0 | 0 | 6212 |
| hsba_domain | - | 0 | 6212 | 0 | 0 |
| hydrophobin_domain | - | 0 | 6212 | 0 | 0 |
| iuis_allergen_homolog | - | 13 | 0 | 6199 | 0 |
| iuis_allergen_similarity | - | 64 | 6148 | 0 | 0 |
| other_not_surface | R0 | 5834 | 378 | 0 | 0 |
| other_surface_no_mechanism | R0 | 317 | 5895 | 0 | 0 |
| serodiagnostic_marker_candidate | R0 | 0 | 6212 | 0 | 0 |
| signal_peptide_protein | R0 | 367 | 5845 | 0 | 0 |
| surface_attachment_candidate | R0 | 50 | 6162 | 0 | 0 |
| tandem_repeat_protein | - | 22 | 6190 | 0 | 0 |
| wall_family_domain | - | 45 | 6167 | 0 | 0 |
| wall_family_narrow_domain | - | 22 | 6190 | 0 | 0 |

## Primary label

One label per protein for the default gate, from the calls above. The first match wins. The calls stay independent flags; the label is a summary for a reader.

| primary label | proteins |
|---|---|
| cell_wall_adhesion_candidate | 29 |
| surface_attachment_candidate | 21 |
| other_surface_no_mechanism | 317 |
| other_not_surface | 5834 |
| mechanism_without_class | 11 |
| not_assessable | 0 |

## Surface attachment: which evidence held

A protein can hold more than one. Hydrophobin and HsbA are surface-active proteins: they adsorb to surfaces. None of these labels shows adhesion.

| evidence | proteins |
|---|---|
| proteins called | 50 |
| tandem repeat | 15 |
| wall family domain | 41 |
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
