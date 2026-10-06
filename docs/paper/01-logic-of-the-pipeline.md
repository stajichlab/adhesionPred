# 01. Logic of the pipeline

*2026-10-06. Sources: `docs/TOOL-ARCHITECTURE.md`, `docs/model-review/2026-09-27-review-and-framework-plan.md`,
`docs/model-review/STATUS.md`, `docs/step1-plain-language-summary.md`,
`docs/superpowers/specs/2026-10-04-orchestrator-design.md`.*

## 1. The problem

"Adhesin", "antigen", "allergen" and "surface protein" are outcomes. Different proteins reach an
outcome by different molecular routes. These routes have different sequence signatures. One
predictor for all of them does not work.

The first classifier showed this. It was trained on adhesins from Saccharomycotina (FLO and ALS
families) against random proteins. Four measurements from the review of 2026-09-27 explain why it
was stopped:

1. **Its benchmark was saturated.** Amino-acid composition plus length, with no language model,
   reached ROC-AUC 0.995 under leave-family-out cross-validation. The cross-validation could not
   tell a good model from a poor one (`docs/model-review/2026-09-27-review-and-framework-plan.md`,
   section 4.1).
2. **It over-called.** On *S. cerevisiae* S288C the shipped 8M model made 77 calls. It found 8 of 9
   known adhesins. Its precision was about 12%. The other calls were cell wall mannoproteins,
   mucin-like sensors and dubious ORFs (same file, section 4.3).
3. **It detected tandem repeats.** Across 95 curated adhesins, the median repeat coverage was 1.000
   for adhesins that it found and 0.000 for adhesins that it missed (Mann-Whitney p = 0.0037;
   section 4.7). Adhesins that bind through one folded interface have no repeat signal.
4. **Adhesin families are clade-specific.** A model trained on one clade and tested on another
   fell to ROC-AUC 0.60. A size-matched model trained inside the clade reached 0.94 (section 4.5).

The conclusion: split the question. Do not predict "adhesin" in one step.

## 2. The steps

| Step | Question | Implementation | State |
|---|---|---|---|
| 1 | Is the protein outside the plasma membrane (wall, outer face, or secreted)? | Rule R0 (SignalP 6 calls a signal peptide), stricter rules R1 and R2, and ESM-2 plus logistic regression models. | Measured (Phase C). No model ships. |
| 2 | Which adhesion mechanism evidence does the protein have? | Tandem-repeat detectors; HMM family scans (for example CFEM, Bys1, hydrophobin). | Repeat detectors exist. Family scans are wrapped but no family is active. |
| 3 | Does the protein serve a stated purpose (serodiagnostic marker, allergen similarity)? | A ranking of *Coccidioides* proteins by genus specificity; similarity to WHO/IUIS fungal allergens. | Ranking exists for *Coccidioides* only. Allergen module has no calibration. |

Step 1 predicts **location**. It does not predict adhesion, glycosylation or immune response.

## 3. Design rules that follow from the evidence

These rules are the "logic" that the paper should state. Each rule answers a failure found in the
review.

| Rule | Failure that it answers |
|---|---|
| Labels come from GO cell component terms with direct (non-homology) evidence, not from sequence similarity. | Positives that were harvested by similarity were partly not adhesins. |
| An unlabelled protein is never a negative. Negatives need an independent reason. | Random proteins as negatives made the boundary trivial. |
| Proteins are grouped by sequence similarity (MMseqs2, 30% identity, 50% coverage). A group never spans training and test. | Half of the first positive set was near-identical FLO11 alleles (581 positives in 53 clusters). |
| Test on other species and other clades, not only on held-out clusters. | A clade-specific model collapses outside its clade. |
| Intervals are 95% cluster bootstrap, not protein bootstrap. | Independent units are clusters, not proteins. |
| A set is an *estimate* only if it has at least 20 direct positives and a recall half-width of 0.10 or less. Otherwise it is a *smoke test*. | Small test sets give numbers that look better than they are. |
| Accuracy gates are set after measurement, not before. | A gate set before measurement is a guess. |
| Negatives are split into strata (N-int inside the cell; N-sec secretory but not surface; PM-TM plasma membrane with a transmembrane helix) and rates are reported per stratum. A curated hard-negative list (cell wall remodelling enzymes, mucin-like sensors, S/T-rich secreted enzymes) is planned; only a 31-row seed list exists (issue #14 is open). | Random negatives under-estimate the false-positive rate on look-alikes. |

## 4. The sorting tool (`cellsurface_sorting_hat`)

The tool takes a proteome and writes one row per protein. A row has, for each evidence module, the
value, the status of that module for the protein's taxon, and the reason.

Layers:

1. **Evidence modules.** Each module has a name, a version, a `params_hash` (all parameters), an
   `artefact_hash` (database or model file), an `applicable()` test and a `status`.
2. **Category logic.** A rules file (`categories.yaml`) combines module values with three-valued
   (Kleene) logic: `called`, `not_called`, `not_assessable`. A call keeps the weakest status among the
   modules that decided it.
3. **Outputs.** Long and wide call tables, evidence, per-protein table, a report and a run record.

Why three values: a module may not apply to a taxon (for example the antigen lookup outside
*Coccidioides*). "Cannot say" must not become "no".

Why status is separate from the call: a module can give a value before anyone has measured it. The
status says how much to trust the value.

Categories in version 1 (all are evidence for a hypothesis, not a finding):

| Name | Meaning |
|---|---|
| `signal_peptide_protein[v]` | Step 1 variant `v` calls the protein outside the membrane. For R0 this means SignalP calls a signal peptide. |
| `tandem_repeat_protein` | One of two repeat detectors calls a repeat region. Evidence only. |
| `wall_family_domain` | A Pfam family linked to wall function in at least one species is present. Evidence only. |
| `cell_wall_adhesion_candidate[v]` | (repeat OR family domain) AND signal peptide. |
| `cocci_specificity_rank_top15` | *Coccidioides* protein in the top 15% of a genus-specificity ranking. It is not epitope prediction. |
| `serodiagnostic_marker_candidate[v]` | Top-15% rank AND signal peptide. |
| `iuis_allergen_similarity`, `iuis_allergen_homolog` | Sequence similarity to WHO/IUIS fungal allergens at two thresholds. Evidence only. IgE cross-reactivity is a hypothesis at most. |
| `other_*` | Signal peptide protein with none of the three mechanism calls, or not a signal peptide protein. |

Version 1 has five categories. Enzyme classes (synthases, glucanases, degradation enzymes) are later
work.

## 5. What the logic does not claim

- It does not claim that a call is an adhesin, an antigen or an allergen.
- It does not claim that a module is accurate outside the taxa in its status entries.
- It does not claim that agreement between two modules of the same kind has been measured. The spec
  says it will be reported as Cohen's kappa with the 2x2 counts, because kappa is unstable at low
  prevalence. That measurement is planned (Plan 2, Task 14) and has not been run.
