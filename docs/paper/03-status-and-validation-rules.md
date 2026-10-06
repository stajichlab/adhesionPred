# 03. Status and validation rules

*2026-10-06. Sources: `docs/superpowers/specs/2026-10-04-orchestrator-design.md`,
`docs/superpowers/plans/2026-10-05-cellsurface-sorting-hat-modules-and-calibration.md` (Global
Constraints), `src/cellsurface_sorting_hat/calibration/`. Section 6 records decisions from the
conversation of 2026-10-06. Those decisions are implemented in PR #64 (branch `sorting-hat-modules`)
and are in effect once that PR is merged.*

## 1. Why a status exists

A module can give a value for a protein before anyone has measured the module. The value then
looks as firm as a measured one. The status separates the two. It is a property of a module, a
module version, and a taxon.

| Status | Meaning |
|---|---|
| `estimated` | Measured with enough truth to give a number. |
| `smoke` | Measured on too little truth for a number. It shows that something works or fails. |
| `unvalidated` | Not measured, or measured on another version of the module. |

A module cannot set its own status. Without a measurement file it is `unvalidated`.

## 2. The rule for `estimated`

All of these must hold for a module in one species:

1. At least 20 positives and at least 20 negatives.
2. At least 20 independent clusters of each class, when cluster counts are recorded.
3. A 95% interval half-width of at most 0.10 for **both** sensitivity and specificity.
4. The module was tested on negatives. Without negatives the status is at most `smoke`.
5. The truth set did not help to set the rule or its cutoffs, and does not overlap the module's
   reference set (leakage rule, section 5).
6. The Phase C label is `estimate`. The status of a Phase C measurement is the weaker of the Phase C
   label and this rule.

Why 20 and 0.10: these are the Phase C rule (ruling C-8, owner answer of 2026-10-01). The spec says
no Phase C code existed at that time. I did not find a record of how the values were chosen.

## 3. Intervals

- Cluster bootstrap by homology group, 2,000 resamples, fixed seed, percentile 2.5 and 97.5.
  Clusters are the independent units, so resampling proteins would give intervals that are too
  narrow. A test in the code pins this: with 100 positive clusters, five of them holding 60 proteins
  each, resampling proteins gives a half-width of about 0.07 and resampling clusters gives about 0.14.
- The interval used is the widest of the cluster bootstrap and the Wilson interval on the number of
  independent clusters of that class. Wilson uses the unrounded proportion with the number of
  clusters as n. A rounded count made the interval up to 0.018 narrower, which is not conservative.
- A tolerance of 1e-9 is allowed on the half-width check, because decimal sums such as 0.8 - 0.6
  give 0.20000000000000007.
- Specificity of a Phase C measurement comes with per-stratum rates (N-int, N-sec, PM-TM). The pooled
  value depends on the mix of negatives.
- Specificity is never inferred from absence in a database. A negative needs an independent reason.
  If a truth set has no negatives, the module gets no specificity and the report says so.

## 4. One species per status entry

A status entry is one species. The species is resolved from the NCBI taxonomy (`names.dmp`, with
`nodes.dmp` for the rank). A name with zero or several IDs stops the run. A protein's taxon may be a
strain below the species. A status applies to a protein only if its taxon is the tested species or a
descendant. It does not apply because two taxa share a broad label such as "Eurotiomycetes".

Why: a pooled estimate describes neither species. A probe in the review showed the failure. Twenty
*S. cerevisiae* and twenty *C. albicans* positives, pooled, made both species `estimated`, although
neither species alone had 20 positives. The code now refuses more than one taxon per call and
refuses a rank above species.

A status written for an older module identity (module name, version, `params_hash`,
`artefact_hash`) is dropped when a new entry is merged. A status source is refused if its identity
differs from the running module.

## 5. Leakage

`--leakage` is required when a truth set is added. The values are `none`, `partial`,
`tuned_on_truth`, `in_reference` and `unknown`. Any value except `none` caps the status at `smoke`.
A truth set that overlaps a module's reference set, or that helped to tune it, must not
support an `estimated` claim.

Panel check (report only, never changes a status): each expected call is marked `agree`,
`disagree`, `not_in_run`, `known_miss`, `excluded` or `excluded_leakage`. Proteins that tuned a module
(for example PRA3, Ag2/PRA, PRA2, SOWgp for the ranking) and their homologs are excluded from the
metrics.

## 6. Decisions of 2026-10-06 (implemented in PR #64)

The owner agreed with these recommendations in conversation. PR #64 implements them. Sections 2 to 5
describe the behaviour before that PR where they differ: the cluster counts are now recorded
(decision 2) and a status now counts only for the call that was measured (decision 3).

| # | Question | Decision | Consequence |
|---|---|---|---|
| 1 | What counts as leakage for rule R0? | Only tuning of R0's own cutoffs and overlap with R0's own reference set cap a status. Overlap between the Phase C positives and the training data of SignalP 6 is **not measured**. It is recorded as a stated limit in the notes and the report. It does not cap the status. | *A. nidulans* stays `estimated`. `phasec` gets no `--leakage` argument. The plan's constraint line is narrowed to say this. |
| 2 | Should `phasec` record cluster counts and apply the full interval rule? | Yes. `phasec` counts clusters per class from the Phase C cluster file and eval table, records them, applies the widest-of rule, and refuses when protein counts differ from `metrics.json`. | The sentence "cluster floor not checked" is removed. Measured on the real files: *A. nidulans* has half-widths of 0.0935 (sensitivity) and 0.0227 (specificity) and stays `estimated`. The other five species stay `smoke`. |
| 3 | Where does a call get its status? | From a measurement of that call. A call with no measurement of its own reports `unvalidated`, with a note naming the call on which the module was measured. | Until per-call status entries exist (a change to the status file format), `truth` refuses a call that reads more than one module. This covers `tandem_repeat_protein` (repeat02 OR repeat14) and `iuis_allergen_homolog` (allergen homology OR Pfam allergen). Each module is calibrated on a call that reads only that module. |
| 4 | How strictly is a status tied to the data it was measured on? | Enforce it. The engine writes module identities into `run.json`. `truth` refuses when they differ from the work directory's module records. `phasec` takes the SignalP module and mode that Phase C used as required arguments, refuses unless they match the work directory's R0 record, and records both in the notes. | A later re-run of SignalP with another build or mode invalidates the R0 status. This needs a change to the core engine's run record. The repository's Phase C job script uses `signalp/6-gpu` with `--mode fast --organism eukarya`, the same as the new job script. |

Naming note: "Phase C" is the held-out evaluation of step 1. A more descriptive name for future
documents is "step 1 held-out benchmark". The path `_workdir/step1_compare/phasec/` and older
documents keep the old name. No rename is planned now.

## 7. What a reader can rely on today

| Claim | Basis |
|---|---|
| R0 is `estimated` in *A. nidulans*. | Re-derived: 109 positives, 164 negatives, half-widths 0.092 and 0.016 from the file's intervals. |
| R0 is `smoke` in the other five Phase C species. | Re-derived: label or fewer than 20 positives. |
| No other module has a status above `unvalidated`. | No calibration has been run (Plan 2 Tasks 13 and 14). |
