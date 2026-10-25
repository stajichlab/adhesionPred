# AGA1-class modules for the sorting hat, and a cross-proteome aggregator: design

Date: 2026-10-26. Status: DRAFT for review. Decisions: the owner's answers of 2026-10-26 (memory note `aga1-sorting-hat-integration-decisions-2026-10-26`),
recorded here with the reasons. Related: `docs/superpowers/specs/2026-10-04-orchestrator-design.md` (the engine), `docs/CLASSES.md` (current calls),
`analysis/aga1_class/` (the tools this spec turns into modules), `docs/reports/2026-10-24-aga1-broad-flag-candidate.md`,
`docs/reports/2026-10-25-aga1-feature-classifier.md`, `docs/reports/2026-10-26-aga1-nonyeast-and-cpar-test.md`.

## 1. Goal

Reports generated from proteomes and aggregated across groups of taxa: counts of gene categories, with the quality of each count visible. Two parts:

1. Make the AGA1-class predictors (PredGPI, AIUPred, the broad and narrow rule flags, the feature classifier) modules and calls of the sorting hat, so one run
   of `cellsurface_sorting_hat` yields them next to the existing calls, with statuses.
2. Add an aggregator that reads the outputs of many runs and a grouping of taxa, and writes tables, a Markdown report and static figures.

Out of scope here (later specs): orthogroups, diversity and paralog counts; genome context for tandem duplication; phylogeny-aware statistics; a scale-out runner;
interactive HTML; tiers as engine calls; expected false calls per count.

## 2. Part 1: modules, calls, thresholds, statuses

### 2.1 Module wrappers (`cellsurface_sorting_hat_module`)

The engine never runs tools. A wrapper turns one tool's output into `<workdir>/modules/<module>.tsv.gz` (+ `.json` run record) through `write_module`, as the existing
wrappers do (`modules/signalp.py` is the template: parse, validate, one row per protein, `invalid_row` for invalid proteins, `error` for a protein the tool did not return).

| Wrapper | Reads | Module table fields | Rules |
|---|---|---|---|
| `gpi` | PredGPI gff3 (`predgpi_raw.gff3`, one row per protein, `Chain` = no GPI, `GPI-anchor` = omega site, score = confidence class 1.0 / 0.70 / 0.55) | `call` (called = GPI-anchor, not_called = Chain), `omega_site`, `confidence_class` (high / probable / weak), `fpr_upper_bound` | A protein of 40 residues or fewer is never scored by PredGPI: its row is `state=not_assessable` (never `not_called`). Unknown row type, `GPI-anchor` with start != end, a `Chain` that does not span the protein, a duplicate or missing id stop the wrapper. |
| `disorder` | AIUPred 2.x raw (`#>ID` record lines, `pos aa score` rows) | `longest_run`, `frac`, `mean`, `n_res` (threshold 0.5 for a residue) | No `call` column: the calls on it read `longest_run` through `{test: ...}`. Scores up to 0.05 outside [0, 1] are accepted (AIUPred reached 1.023), further out stops. The residue count must equal the FASTA length. The per-residue file carries residues and is read, never copied into the work directory. |
| `aga1_model` | the score table of `feature_classifier.py apply` for the three models (`score_logreg`, `score_rf`, `score_gboost`), and the thresholds table (2.3) | `score_logreg`, `score_rf`, `score_gboost`, `thr_logreg`, `thr_rf`, `thr_gboost`, `n_models_call`, `threshold_clade`, `calibrated` (yes/no), `call` (called = at least 2 of 3 models at or above their threshold) | The thresholds come from the data table by lineage (2.3), so the wrapper needs `--taxon` and `--taxdump`. Missing score for a protein gives `error`. The wrapper never imports scikit-learn. |

The existing `signalp` (`step1_rule@R0`), `repeat14` and `repeat02` modules are reused unchanged. The Pfam and TM modules are those already in the engine (`pfam`, `tm`);
the cluster jobs that produce their input are `TOOL=pfam_family` and `TOOL=tmhmm` of `analysis/aga1_class/proteome_array.sbatch`.

### 2.2 Calls (`categories.yaml`, category "cell wall and adhesion")

| Call | Kind | Definition |
|---|---|---|
| `gpi_anchored_protein` | evidence | `{call: gpi}` |
| `long_disordered_protein` | evidence | `{test: {module: disorder, field: longest_run, op: ">=", value: $disorder_run_min}}` with `thresholds.disorder_run_min: 200` |
| `ser_thr_gpi_wall_candidate` (the broad flag) | composite | `and: [signal_peptide_protein, gpi_anchored_protein, long_disordered_protein]` (gated like the other composites) |
| `flocculin_like_candidate` (the narrow flag) | composite | `and: [ser_thr_gpi_wall_candidate, {call: repeat14}]`. Detector 14 only, as measured; `tandem_repeat_protein` (02 or 14) is NOT used here. |
| `model_surface_candidate` | evidence | `and: [{call: aga1_model}, signal_peptide_protein]` |

Three-valued logic applies as everywhere: a protein with no GPI call (40 residues or fewer) is `not_assessable` for the broad flag unless another part is `not_called`
(a not-called part makes the composite `not_called`). The tiers of the AGA1 candidate list are not engine calls. `docs/class_descriptions.yaml` gets a description and a
"known limits" sentence for each call (the flags do not see short CFEM proteins; enzymes with a Ser/Thr-rich linker are flagged; the model call needs a signal peptide and an
uncalibrated threshold outside yeast), and `docs/CLASSES.md` is regenerated by `analysis/class_table/make_classes_md.py` (a test checks it is current).

### 2.3 Thresholds by clade: `data/sorting_hat/aga1_model_thresholds.tsv`

A signed data table in the style of `family_table.tsv`. Columns: `clade_taxon` (NCBI taxon id, lineage match), `clade_name`, `model` (logreg / rf / gboost), `threshold`,
`basis` (how it was set), `calibrated` (yes / no), `signed_by`, `signed_date`, `evidence` (path of the report or record).

Initial rows: Saccharomycotina (taxon 4891) for each model, value = the score at which the 850 labelled negatives have a 3.5% false-positive rate, out of fold
(logreg 0.56834, rf 0.25448, gboost 0.11414; `docs/reports/data/aga1_proteome_test_2026-10-25/thresholds.tsv`), `calibrated=yes`, basis "yeast 850, leave-cluster-out"; and a fallback row
for taxon 1 (all other taxa) with the same values and `calibrated=no`. The wrapper picks the row of the deepest clade that is an ancestor of the protein's taxon (or the taxon itself);
the engine's `taxonomy.py` supplies the lineage. If no row matches, the wrapper stops. The chosen threshold, clade and `calibrated` flag are written into every module row, so a
report states which calibration produced each call. Adding a calibrated clade after the hand-check is adding signed rows; no code changes.

### 2.4 Statuses

* The rule calls (`gpi_anchored_protein`, `long_disordered_protein`, `ser_thr_gpi_wall_candidate`, `flocculin_like_candidate`) get measured status files for *S. cerevisiae* (taxon 4932 and the strain taxon) and
  *C. albicans* (5476) with `cellsurface_sorting_hat_calibrate truth`: truth = the 850 labelled proteins (tier A positives; hard negatives and background negatives; the excluded roles are not truth),
  status `smoke` (8 and 10 tier A proteins), with a note that the 200-residue cut and the repeat call were chosen on the same proteins (leakage `design` or the closest allowed value of the status schema).
* `model_surface_candidate` has no status file (`unvalidated`) until held-out verdicts exist: the classifier was trained on those proteins, so a measurement on them is in-sample. The leave-cluster-out numbers stay
  in the classifier report.
* A status holds for one species and taxon, as for every call; a species without a file is `unvalidated`.

## 3. Part 2: the aggregator (`cellsurface_sorting_hat_aggregate`)

### 3.1 Inputs

* A list of engine output directories (each with `calls.long.tsv.gz`, `proteins.tsv.gz`, `run.json`), or a manifest TSV (`run_dir`, optional `label`).
* `--taxdump nodes.dmp` and optional `--names names.dmp`.
* `--rank RANK` (default `class`) and optional `--groups TABLE` (TSV `taxon_id`, `group`, optional `note`): a table row overrides the lineage group for that taxon and for taxa below it that have no more specific row.

Stops (FATAL, exit 1) on a run directory without the three files, a duplicate proteome label, a taxon not found in `nodes.dmp`, a malformed group table, or runs made with different `categories.yaml` hashes
(the config hash is in `run.json`) unless `--allow-mixed-config` is given, in which case the report lists the mixture.

### 3.2 Tables

* `counts_by_proteome.tsv`: per proteome and call: `proteins`, `called`, `not_called`, `not_assessable`, `not_applicable`, `assessable` (= called + not_called), `called_per_1000_assessable`, `status` (the call's status for that taxon: unvalidated / smoke / estimated), `status_species`, `measured` (yes/no), `unvalidated_marker`.
* `counts_by_group.tsv`: per group and call: `n_proteomes`, `n_measured`, mean / median / min / max of `called_per_1000_assessable` over all proteomes, the same over measured proteomes only, and the total called.
* `proteomes.tsv`: proteome label, taxon, group, number of proteins, number of invalid proteins, config hash, engine version.
* Both count tables always carry the status columns; nothing is dropped for being unvalidated. `--measured-only` restricts the group summaries to measured proteomes and says so in the report.
* Denominator rule (owner decision): per 1,000 **assessable** proteins for that call and proteome, never per all proteins; the number of proteins is listed next to it.

### 3.3 Report and figures

`report.md` (same style as the engine report): what was aggregated (proteomes, groups, config hash), a table of group summaries for each category of calls, the status coverage (which cells are measured), warnings
(mixed configs, proteomes with many invalid proteins, groups with one proteome), and known limits (a count of an unvalidated call is not a prevalence; the thresholds of the model call are calibrated for Saccharomycotina only).
Figures (PNG, matplotlib, deterministic, saved next to the report and referenced from it): (1) heatmap of called per 1,000 assessable by proteome (rows) and call (columns), unvalidated cells hatched; (2) bar chart per group and call with the number of
proteomes in each bar and error bars only where `n_proteomes >= 3`; (3) status-coverage panel (proteome x call, colour = status). The aggregator never reads sequences.

### 3.4 What is not claimed

Counts are counts of calls. They are not gene numbers per family and not false-positive-corrected; expected false calls from a measured specificity come later and need no change to the table layout (the status columns are already present).

## 4. Tests (strict, tests first)

* Wrapper tests with SYNTHETIC tool output: gpi (each row type, <= 40 residues, malformed rows, duplicates), disorder (layout `#>ID`, slack, length mismatch), aga1_model (threshold choice by lineage, fallback row, stop when no row, 2-of-3 logic, missing score).
* Engine tests: the five calls on a synthetic work directory with every combination of the three-valued inputs; `CLASSES.md` current; config loads; `thresholds.disorder_run_min` read.
* Thresholds table: schema checks, one signed row per clade and model, `calibrated` values, the Saccharomycotina values equal the committed record.
* Aggregator tests: synthetic runs, rank grouping, group-table override, denominators, status columns, `--measured-only`, mixed-config stop, figure files exist and are non-empty.
* End-to-end test on committed small fixtures (a few proteins), plus a real local run on the 14 proteomes as the acceptance check (module tables from the existing SignalP, PredGPI, AIUPred, repeat and classifier outputs; modules without a table are `not_assessable`).

## 5. Work packages (each a PR)

1. This spec and a review; the plan.
2. Module wrappers `gpi`, `disorder`, `aga1_model` + thresholds table + tests.
3. Calls in `categories.yaml`, `class_descriptions.yaml`, regenerated `CLASSES.md`, engine tests.
4. Local run on the 14 proteomes (acceptance), measured status files for the rule calls in *S. cerevisiae* and *C. albicans*.
5. The aggregator and its tests.
6. First aggregate report on the 14 proteomes; README and `docs/TECHNICAL-NOTES.md` updated.

Cluster work (the Pfam and TMHMM arrays) runs in parallel and each submission needs the owner's go-ahead.

## 6. Checked against the code (2026-10-26)

* A numeric threshold on a module field is expressed with `{test: {module, field, op, value}}`, and `"$name"` reads `thresholds` (header of `categories.yaml`); `disorder` therefore needs no `call` column.
* `taxonomy.Lineage.ancestors(taxon)` returns the taxon and its parents from `nodes.dmp` alone (no `names.dmp`); a status entry for a taxon applies to its descendants (`best_entry`), so the species taxa 4932 (*S. cerevisiae*) and 5476
  (*C. albicans*) cover their strains. `MIN_TAXON = 2` applies to status entries and `species_limited`; the thresholds table is its own file with its own validation and uses taxon 1 explicitly for the "all other taxa" fallback.
* The config hash is recorded by the CLI (`config_sha256`) and is therefore available to the aggregator through `run.json`; the aggregator test must confirm the exact key.

## 7. Questions for the review

* Is a `smoke` status justified for the rule calls given 8 and 10 tier A proteins and the in-sample cut (2.4), or should the status schema's leakage field force `unvalidated`?
* Should `flocculin_like_candidate` read `repeat14` only, or both repeat modules? (The measured narrow flag used detector 14 only; `tandem_repeat_protein` uses 02 or 14.)
* Should the aggregator emit proteome-level warnings for proteomes with many `not_assessable` GPI rows (proteins of 40 residues or fewer), since they change the denominators?
