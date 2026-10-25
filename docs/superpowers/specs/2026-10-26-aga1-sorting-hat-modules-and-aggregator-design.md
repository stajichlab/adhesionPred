# AGA1-class modules for the sorting hat, and a cross-proteome aggregator: design

Date: 2026-10-26. Status: revision 2 after an adversarial review (2026-10-26; findings verified against the engine code); the owner accepted the status plan of 2.4 on 2026-10-26. The plan is `docs/superpowers/plans/2026-10-26-aga1-sorting-hat-modules-and-aggregator.md`. Decisions: the owner's answers of 2026-10-26 (memory note `aga1-sorting-hat-integration-decisions-2026-10-26`),
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
| `disorder` | AIUPred 2.x raw (`#>ID` record lines, `pos aa score` rows) | `longest_run`, `frac`, `mean`, `n_res` (threshold 0.5 for a residue) | No `call` column: the calls on it read `longest_run` through `{test: ...}`. Scores up to 0.05 outside [0, 1] are accepted (AIUPred reached 1.023), further out stops. The residue count must equal the FASTA length. The per-residue file carries residues; it is **streamed** (13.76 million rows for 27,811 proteins), never loaded whole and never copied into the work directory. The length check uses the FASTA sequence after the trailing `*` strip (`fasta.py`). |
| `aga1_model` | the score table of `feature_classifier.py apply` for the three models (`score_logreg`, `score_rf`, `score_gboost`), and the thresholds table (2.3) | `score_logreg`, `score_rf`, `score_gboost`, `thr_logreg`, `thr_rf`, `thr_gboost`, `n_models_call`, `threshold_clade`, `calibrated` (yes/no), `call` (called = at least 2 of 3 models at or above their threshold) | The thresholds come from the data table by lineage (2.3), so the wrapper needs `--taxon` (one species-level taxon, or `--taxon-map`) and `--taxdump`. A missing score for a protein gives `error`. The wrapper never imports scikit-learn. **Identity:** the `ModuleSpec` carries `artefact_digest` = hash of the thresholds table rows used + sha256 of the score table, and `params` holds the model names and the 2-of-3 rule, so a changed threshold or score table changes the module identity and makes status files stale (`modules/base.py`). |

Protein ids are the first token of the FASTA header (`fasta.py`); the staged proteome chunks have bare accessions, a raw UniProt FASTA gives `sp|ACC|NAME`. Each wrapper stops on an id that does not match its input (`_check_ids`). The existing `signalp` (`step1_rule@R0`), `repeat14` and `repeat02` modules are reused unchanged. The Pfam and TM modules are those already in the engine (`pfam`, `tm`);
the cluster jobs that produce their input are `TOOL=pfam_family` and `TOOL=tmhmm` of `analysis/aga1_class/proteome_array.sbatch`.

### 2.2 Calls (`categories.yaml`, category "cell wall and adhesion")

The engine's rules shape the call set (verified in `engine.py`): a call with a `ref` node cannot get a status file (`call_eligible`), an ungated call cannot `ref` a per_variant call, a call may not name a step 1 module literally,
and `calibrate truth --module` accepts only a call that reads one module (a call that reads several is measured with `--call-status`). The signal-peptide call is per_variant, so every call that uses it is a gated, derived composite.

| Call | Kind | Gate | Definition |
|---|---|---|---|
| `gpi_anchored_protein` | evidence | no | `{call: gpi}` |
| `long_disordered_protein` | evidence | no | `{test: {module: disorder, field: longest_run, op: ">=", value: "$disorder_run_min"}}`, `thresholds.disorder_run_min: 200` |
| `gpi_long_disorder_protein` | evidence | no | `and: [{call: gpi}, {test: disorder.longest_run >= $disorder_run_min}]` (reads two modules, no `ref`: measurable) |
| `ser_thr_gpi_wall_candidate` (the broad flag) | composite | yes, per_variant | `and: [{ref: signal_peptide_protein}, {ref: gpi_long_disorder_protein}]`; derived status |
| `flocculin_like_candidate` (the narrow flag) | composite | yes, per_variant | `and: [{ref: ser_thr_gpi_wall_candidate}, {call: repeat14}]`. Detector 14 only, as measured; `tandem_repeat_protein` (02 or 14) is NOT used. |
| `model_surface_candidate` | composite | yes, per_variant | `and: [{call: aga1_model}, {ref: signal_peptide_protein}]`; its status is the weakest of its leaves (`aga1_model` unvalidated) |

Three-valued logic applies as everywhere: a protein with no GPI call (40 residues or fewer) is `not_assessable` for the broad flag unless another part is `not_called`. **Primary label:** `labels.py` files a protein that only the new composites
call under `other_surface_no_mechanism`; the new composites join the mechanism lists of the `other_*` calls and `MECHANISM_CALLS` so such a protein gets a mechanism label (decided in the plan; the label tests are updated). The tiers of the AGA1 candidate list are
not engine calls. `docs/class_descriptions.yaml` gets a description and a "known limits" sentence per call (the flags do not see short CFEM proteins; enzymes with a Ser/Thr-rich linker are flagged; the model call has an uncalibrated threshold outside yeast), and
`docs/CLASSES.md` is regenerated by `analysis/class_table/make_classes_md.py`. The thresholds table and the three module names are added to the `evidence:` block (`aga1_model.threshold_clade`, `aga1_model.calibrated`, `gpi.confidence_class`) so the aggregator reads them
from `evidence.tsv.gz`.

### 2.3 Thresholds by clade: `data/sorting_hat/aga1_model_thresholds.tsv`

A signed data table in the style of `family_table.tsv`. Columns: `clade_taxon` (NCBI taxon id, lineage match), `clade_name`, `model` (logreg / rf / gboost), `threshold`,
`basis` (how it was set), `calibrated` (yes / no), `signed_by`, `signed_date`, `evidence` (path of the report or record).

Initial rows: Saccharomycotina (taxon 4891) for each model, value = the score at which the 850 labelled negatives have a 3.5% false-positive rate, out of fold
(logreg 0.56834, rf 0.25448, gboost 0.11414; `docs/reports/data/aga1_proteome_test_2026-10-25/thresholds.tsv`), `calibrated=yes` meaning "set from yeast out-of-fold scores of the labelled 850 (a mixed set: 21 of 35 tier A are Saccharomycotina)", not a validated false-positive rate for the clade; and a fallback row
for taxon 1 (all other taxa) with the same values and `calibrated=no`. The wrapper picks the row of the deepest clade that is an ancestor of the protein's taxon (or the taxon itself);
the engine's `taxonomy.py` supplies the lineage. If no row matches, the wrapper stops. The chosen threshold, clade and `calibrated` flag are written into every module row, so a
report states which calibration produced each call. Adding a calibrated clade after the hand-check is adding signed rows; no code changes.

### 2.4 Statuses (revised after the review)

What the engine can and cannot measure decides this section:

* **Measurable by `calibrate truth`:** calls without `ref` that read exactly one module with `--module` (`gpi_anchored_protein`, `long_disordered_protein`), and calls that read two or more modules with `--call-status` (`gpi_long_disorder_protein`).
  The composites (`ser_thr_gpi_wall_candidate`, `flocculin_like_candidate`, `model_surface_candidate`) contain `ref` and are "derived": the engine reports the weakest status of their leaf calls.
* **Truth must match the property.** The class truth (tier A vs hard and background negatives) measures class membership, not "is GPI-anchored" or "has a 200-residue run": many hard negatives are GPI-anchored. Statuses are therefore assigned as follows:
  * `gpi_long_disorder_protein` (the property combination of the broad flag, without the signal peptide) is measured against the class truth for *S. cerevisiae* (4932) and *C. albicans* (5476), `smoke` (few tier A proteins per species; the actual numbers
    come from `proteins.tsv` of the scoping record, not from the proteome-test report), leakage `tuned_on_truth` (the 200-residue cut was chosen on these proteins). The status file says in its notes that it measures agreement with the AGA1 class labels.
  * `gpi_anchored_protein`: measured against UniProt GPI-anchor features where the 850 table has them (`uniprot_gpi_feature`, about half of tier A), or left `unvalidated`; decided in the plan after counting.
  * `long_disordered_protein`: `unvalidated` (no truth for the property).
  * `model_surface_candidate`: `unvalidated` (the classifier was trained on these proteins; in-sample).
* `calibrate truth --taxa` takes one species-level taxon; a status for 4932 covers its strains by lineage.
* A species without a status file is `unvalidated`.

## 3. Part 2: the aggregator (`cellsurface_sorting_hat_aggregate`)

### 3.1 Inputs

* A list of engine output directories (each with `calls.long.tsv.gz`, `proteins.tsv.gz`, `run.json`), or a manifest TSV (`run_dir`, optional `label`).
* `--taxdump nodes.dmp` and optional `--names names.dmp`.
* `--rank RANK` (default `class`) and optional `--groups TABLE` (TSV `taxon_id`, `group`, optional `note`): a table row overrides the lineage group for that taxon and for taxa below it that have no more specific row.

Stops (FATAL, exit 1) on a run directory without the three files, a duplicate proteome label, a taxon not found in `nodes.dmp`, a malformed group table, or runs made with different `categories.yaml` hashes
(the config hash is in `run.json`) unless `--allow-mixed-config` is given, in which case the report lists the mixture.

Calls are read from `calls.long.tsv.gz` for the **default gate only** (`step1_rule@R0`): the file has one row per step 1 variant for per_variant calls, and counting all variants would multiply the counts. Optional `--variant` selects another.
`threshold_clade` and `calibrated` come from `evidence.tsv.gz` (they are in the `evidence:` block), not from the module table. `nodes.dmp` ranks are read with `calibration/cli.read_nodes` (`taxonomy.Lineage` stores parents only); an unknown taxon stops the aggregator.

### 3.2 Tables

* `counts_by_proteome.tsv`: per proteome and call: `proteins`, `called`, `not_called`, `not_assessable`, `not_applicable`, `assessable` (= called + not_called), `called_per_1000_assessable`, `status` (the call's status for that taxon: unvalidated / smoke / estimated), `status_species`, `measured` (yes/no), `unvalidated_marker`.
* `counts_by_group.tsv`: per group and call: `n_proteomes`, `n_measured`, mean / median / min / max of `called_per_1000_assessable` over all proteomes, the same over measured proteomes only, and the total called.
* `proteomes.tsv`: proteome label, taxon, group, number of proteins, number of invalid proteins, config hash, engine version.
* Both count tables always carry the status columns; nothing is dropped for being unvalidated. `--measured-only` restricts the group summaries to measured proteomes and says so in the report.
* Denominator rule (owner decision): per 1,000 **assessable** proteins for that call and proteome, never per all proteins; the number of proteins is listed next to it. If `assessable == 0` (for example a proteome whose `disorder` table is missing makes every rule call `not_assessable`) the rate is `NA`, never 0, and the report warns, including for calls whose module table is absent.
* Group means of per-1,000 rates over unequal proteome sizes are descriptive, not estimates; `n_proteomes` and the raw called and assessable counts are always printed beside them, and a group of one or two proteomes is flagged.

### 3.3 Report and figures

`report.md` (same style as the engine report): what was aggregated (proteomes, groups, config hash), a table of group summaries for each category of calls, the status coverage (which cells are measured), warnings
(mixed configs, proteomes with many invalid proteins, groups with one proteome), and known limits (a count of an unvalidated call is not a prevalence; the thresholds of the model call are calibrated for Saccharomycotina only).
Figures (PNG, matplotlib with `metadata={"Software": None}` and fixed fonts so the bytes are reproducible, saved next to the report and referenced from it): (1) heatmap of called per 1,000 assessable by proteome (rows) and call (columns), unvalidated cells hatched; (2) bar chart per group and call with the number of
proteomes in each bar and error bars only where `n_proteomes >= 3`; (3) status-coverage panel (proteome x call, colour = status). The aggregator never reads sequences.

### 3.4 What is not claimed

Counts are counts of calls. They are not gene numbers per family and not false-positive-corrected; expected false calls from a measured specificity come later and need no change to the table layout (the status columns are already present).

## 4. Tests (strict, tests first)

* Wrapper tests with SYNTHETIC tool output: gpi (each row type, <= 40 residues, malformed rows, duplicates), disorder (layout `#>ID`, slack, length mismatch), aga1_model (threshold choice by lineage, fallback row, stop when no row, 2-of-3 logic, missing score).
* Engine tests: the five calls on a synthetic work directory with every combination of the three-valued inputs; `CLASSES.md` current; config loads; `thresholds.disorder_run_min` read.
* Identity and status: a changed threshold row or score table changes the `aga1_model` identity; a call-status file for a composite is refused; `gpi_long_disorder_protein` accepts `--call-status`; a stale status file is not used.
* Thresholds table: schema checks, one signed row per clade and model, `calibrated` values, the Saccharomycotina values equal the committed record.
* Aggregator tests: synthetic runs, rank grouping, group-table override, denominators, status columns, `--measured-only`, mixed-config stop, figure files exist and are non-empty.
* End-to-end test on committed small fixtures (a few proteins), plus a real local run on the 10 proteomes as the acceptance check (module tables from the existing SignalP, PredGPI, AIUPred, repeat and classifier outputs; modules without a table are `not_assessable`).

## 5. Work packages (each a PR)

1. This spec and a review; the plan.
2. Module wrappers `gpi`, `disorder`, `aga1_model` + thresholds table + tests.
3. Calls in `categories.yaml`, `class_descriptions.yaml`, regenerated `CLASSES.md`, engine and label tests. This changes the config hash, the golden run and `report.expected.md` (new absent modules appear in the warnings), so the golden files are regenerated in their own commit.
4. Local run on the 10 proteomes (acceptance), status files as in 2.4.
5. The aggregator and its tests.
6. First aggregate report on the 10 proteomes; README and `docs/TECHNICAL-NOTES.md` updated.

Cluster work (the Pfam and TMHMM arrays) runs in parallel and each submission needs the owner's go-ahead.

## 6. Checked against the code (2026-10-26)

* A numeric threshold on a module field is `{test: {module, field, op, value}}` with `"$name"` reading `thresholds` (header of `categories.yaml`); `disorder` needs no `call` column.
* `taxonomy.Lineage.ancestors(taxon)` returns the taxon and its parents from `nodes.dmp` alone; a status entry for a taxon applies to its descendants (`best_entry`), so 4932 and 5476 cover their strains. `MIN_TAXON = 2` applies to status entries and
  `species_limited`; the thresholds table is its own file and uses taxon 1 explicitly for the fallback.
* Engine rules that shaped 2.2 and 2.4: an ungated call cannot `ref` a per_variant call (`_validate_node`); `call_eligible` refuses `ref`, `kind: other` and literal step 1 modules; `calibrate truth --module` refuses a call that reads several modules; `LEAKAGE_VALUES` are none, partial, tuned_on_truth, in_reference, unknown; an absent module table makes a leaf `not_assessable` (`_row`).
* The config hash is recorded by the CLI (`config_sha256`) and reaches `run.json`; the aggregator test confirms the exact key.

## 7. Questions answered by the review

* `smoke` is what the schema yields, but only for calls that can be measured against truth that matches the property (2.4).
* The narrow flag reads `repeat14` only, as measured.
* The aggregator warns on proteomes with many `not_assessable` GPI rows (proteins of 40 residues or fewer) and, more importantly, on absent module tables, which change the denominators most.
