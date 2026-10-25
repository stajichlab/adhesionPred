# adhesionPred: tools for fungal cell-surface and cell-wall proteins

This repository holds a set of tools and analyses for finding and sorting fungal cell-surface proteins: signal-peptide proteins,
GPI-anchored cell-wall proteins, tandem-repeat adhesins (flocculin/ALS/EPA-like), wall-family domain proteins and a few others.

There is **no single "adhesin predictor"**. "Adhesin", "antigen" and "biofilm factor" are outcomes, not sequence features, and the
evidence in this repository showed that one model for all of them is the wrong target (`docs/TOOL-ARCHITECTURE.md`). The
original ESM-2 classifier turned out to be a surface-glycoprotein detector rather than an adhesin predictor, and no
trained model ships with the package. What you can use today are the components below. Each states what it measures and what it
does not.

## What is usable now

| Component | What it does | Status |
|---|---|---|
| [`cellsurface_sorting_hat`](#1-cellsurface_sorting_hat-sort-a-proteome-into-surface-categories) | Reads tool results for a proteome (SignalP, Pfam, repeat detectors, TMHMM, allergen BLAST) and writes per-protein calls with explicit evidence and a three-valued logic (called / not called / not assessable) | Engine, module wrappers, calibration commands and HPCC job scripts exist. Each call has a measured status per species where truth data allowed it (`docs/CLASSES.md`) |
| [AGA1-class tools](#2-aga1-class-tools-flocculin-like-sertthr-rich-gpi-wall-proteins) (`analysis/aga1_class/`) | Rule-based flags and a features-only classifier for Ser/Thr-rich, repeat-bearing, GPI-anchored wall proteins; scripts to run SignalP, PredGPI and AIUPred on whole proteomes and summarise them | Research tools with measured performance on labelled sets (see `docs/reports/`); not a validated predictor |
| Repeat detectors (`analysis/cocci_repeats/02_repeat_profile.py`, `14_repeat_detect_general.py`) | Tandem-repeat arrays per protein (period, copies, coverage); script 14 is the project's general detector | Used by the sorting hat and the AGA1 tools |
| Pfam family table (`data/sorting_hat/family_table.tsv`) | Which Pfam families count as wall, hydrophobin, allergen or HsbA evidence, with the specificity tier and owner sign-off for each | Curated; every row says what was measured |
| [`surface_glyco_*`](#3-surface_glyco-the-older-esm-2-surface-glycoprotein-scorer) (older) | ESM-2 embeddings plus logistic regression that scores a protein as a secreted cell-surface glycoprotein | Still works, but no trained model ships; train your own. Not an adhesin predictor and not used by the newer tools |

Reports with the measurements behind these are in `docs/reports/` and the model review is in `docs/model-review/`.

## Install

```bash
git clone https://github.com/stajichlab/adhesionPred.git
cd adhesionPred
pip install -e ".[dev]"      # Python 3.11 or newer; installs torch, scikit-learn, biopython, pandas ...
```

The external tools are **not** installed by the package, and some have academic licences. Install them yourself where you need them:
SignalP 6, PredGPI, AIUPred (a GPU is much faster), HMMER 3 with Pfam-A, TMHMM, BLAST+, and NCBI `nodes.dmp` (from the taxonomy
`taxdump`) for the sorting hat. The SLURM scripts use cluster-specific module names and partitions that you will need to change
(`docs/TECHNICAL-NOTES.md`).

## 1. `cellsurface_sorting_hat`: sort a proteome into surface categories

Three commands:

| Command | Purpose |
|---|---|
| `cellsurface_sorting_hat_module <tool> ...` | Convert one tool's output into a module table in a work directory. Tools: `signalp`, `pfam`, `repeat02`, `repeat14`, `allergen`, `antigen`, `cys`, `expression`, `hydrophobin_relaxed`, `cys8`, `tm` |
| `cellsurface_sorting_hat --fasta ... --workdir ... --out ...` | The engine. It does **not** run the tools; it reads `<workdir>/modules/<module>.tsv.gz` (plus optional run records and status files), evaluates the rules in `src/cellsurface_sorting_hat/categories.yaml` and writes the calls |
| `cellsurface_sorting_hat_calibrate <task> ...` | Measure a call against a truth table (`truth`), check a Pfam family's specificity (`pfam-specificity`), leave-species-out for the allergen set (`allergen-lso`), validate control sets (`controls`), check against a panel (`panel`), and convert the Phase C metrics (`phasec`) |

Typical run for one proteome:

```bash
# 1. run the tools (on a cluster, scripts/sorting_hat/submit_modules.sh submits SignalP, Pfam, repeats, allergen BLAST, TMHMM)
PROJ_ROOT=$PWD FASTA=proteome.faa WORKDIR=work FAMILY_TABLE=data/sorting_hat/family_table.tsv \
  PFAM_RELEASE=38.2 ALLERGEN_FASTA=iuis_allergens.faa bash scripts/sorting_hat/submit_modules.sh

# 2. convert each tool output to a module table (one command per tool; the arguments differ, see --help)
cellsurface_sorting_hat_module signalp --fasta proteome.faa --workdir work --results signalp_out \
  --signalp-version 6.0h-gpu
cellsurface_sorting_hat_module repeat14 --fasta proteome.faa --workdir work --table repeat14.tsv \
  --script analysis/cocci_repeats/14_repeat_detect_general.py
# ... pfam, tm, allergen, as needed

# 3. sort (TAXON_ID is the NCBI taxon of the proteome)
cellsurface_sorting_hat --fasta proteome.faa --taxon TAXON_ID --taxdump nodes.dmp --workdir work --out out
```

Outputs in `out/`: `calls.long.tsv.gz`, `calls.wide.tsv.gz`, `evidence.tsv.gz`, `proteins.tsv.gz`, `report.md`, `run.json`.

How to read the calls (full list with definitions: `docs/CLASSES.md`, generated from the rules):

* **Evidence calls** say what a tool found (`signal_peptide_protein`, `tandem_repeat_protein`, `wall_family_domain`,
  `hydrophobin_domain`, `iuis_allergen_similarity` ...). **Composites** (`cell_wall_adhesion_candidate`, `surface_attachment_candidate`)
  say what a rule combined. `other_*` calls are the residual categories.
* A call can be **called**, **not called** or **not assessable** (an input is missing or does not apply). Nothing is silently counted as 0.
* Each call carries a **status for one species and taxon**: `unvalidated` (no measurement), `smoke` (too little truth for a number) or
  `estimated` (enough truth). A status only counts where its status file is installed and current. The status of a Pfam family
  (`data/sorting_hat/family_table.tsv`) is separate from the status of a call.
* `narrow_ok` in the family table is signed by the owner per family. A wall family reaches the narrow call only with a specificity tier of
  `smoke` or better and that sign-off.

Design notes: `docs/superpowers/specs/2026-10-04-orchestrator-design.md`; the plans that built it are in `docs/superpowers/plans/`.

## 2. AGA1-class tools: flocculin-like, Ser/Thr-rich, GPI wall proteins

Target class: cell-wall proteins that are Ser/Thr-rich, carry repeats and are held by a GPI anchor (AGA1, FLO1/5/9/11, ALS, HWP1,
EPA ...). The class and its labels are in `docs/reports/2026-10-24-aga1-class-scoping.md`; `analysis/aga1_class/README.md` has the
format notes for every script.

**The two flags** (rules, no training; `broad_flag.py`, and `rule_flags` in `feature_classifier.py`):

* **Broad**: SignalP signal peptide, **and** PredGPI GPI call, **and** an AIUPred disordered run of at least 200 residues
  (the length cut is a parameter). This finds Ser/Thr-rich GPI wall mannoproteins in general
  (SED1/DAN4-like proteins and the ALS/FLO/EPA adhesins).
* **Narrow**: the broad flag **and** the project repeat call (detector 14: coverage at least 0.25 and at least 2.5 copies). Stricter;
  it misses adhesins whose repeats are not periodic.
* An optional Ser/Thr "stalk" rule is available (`stalk_features.py`).

**The classifier** (`feature_classifier.py`): logistic regression, random forest and gradient boosting on per-protein features
(composition, disorder, repeat score, signal peptide and GPI calls, tail hydrophobicity, length). It is evaluated with whole sequence
clusters held out, and `apply` scores a new feature table with a model fit on the labelled set. Thresholds do not carry across
clades, so recalibrate one for each new clade (performance: `docs/TECHNICAL-NOTES.md`).

**Run it on a proteome** (all proteins are keyed by UniProt accession):

```bash
# 1. stage chunks of at most 1,000 proteins (writes outside the repo; refuses an existing --out)
python3 -I analysis/aga1_class/stage_proteomes.py --out /path/NEW_DIR Scer=Scer.faa Calb=Calb.faa

# 2. run the tools, one SLURM array task per chunk (TOOL = signalp | predgpi | aiupred; the header of the script has the
#    sbatch options, including those for AIUPred on a GPU)
sbatch --array=1-<N>%4 --export=ALL,PROJ_ROOT=...,CODE_COMMIT=...,IN=...,OUT=...,TOOL=signalp \
  analysis/aga1_class/proteome_array.sbatch
python analysis/cocci_repeats/14_repeat_detect_general.py NEW_DIR/chunks/*.faa --out rep14.tsv   # local, minutes

# 3. merge the chunks and convert to one table per protein
python3 -I analysis/aga1_class/merge_proteome_chunks.py --chunks NEW_DIR/chunks.tsv --results RESULTS --out MERGED
python3 -I analysis/aga1_class/convert_predictions.py --manifest NEW_DIR/manifest.txt \
  --signalp MERGED/signalp/prediction_results.txt --predgpi MERGED/predgpi/predgpi_raw.gff3 \
  --lengths lengths.tsv --disorder MERGED/aiupred/disorder_raw.txt --disorder-tool aiupred --out predictions.tsv

# 4. features, frozen classifier scores, per-species summary and candidate list
python analysis/aga1_class/feature_classifier.py features --proteins proteins.tsv --predictions predictions.tsv \
  --stalk stalk.tsv --rep14 rep14.tsv --raw MERGED/aiupred/disorder_raw.txt --out features.tsv
python analysis/aga1_class/feature_classifier.py apply --train-features \
  docs/reports/data/aga1_feature_classifier_2026-10-25/features.tsv --new-features features.tsv --model rf --out scores_rf.tsv
python3 -I analysis/aga1_class/proteome_summary.py --features features.tsv --scores rf=scores_rf.tsv ... \
  --oof docs/reports/data/aga1_feature_classifier_2026-10-25/oof_scores.tsv ... --out-dir SUMMARY
```

(`stalk.tsv` comes from `stalk_features.py --raw ... --manifest ... --out stalk.tsv`; `lengths.tsv` and `proteins.tsv` come from the
staging step. The `...` marks arguments you will find in each script's `--help`.)

The scripts refuse to overwrite existing outputs and stop with a clear error instead of writing silent zeros. FASTA files and per-residue tool
output are kept out of the repository (conventions: `docs/TECHNICAL-NOTES.md`).

The flags target flocculin-like proteins; other GPI wall families, such as the CFEM proteins, are not seen by them and are covered by
Pfam families in the family table. Measured performance and its limits: `docs/TECHNICAL-NOTES.md` and the reports in `docs/reports/`.

## 3. `surface_glyco`: the older ESM-2 surface-glycoprotein scorer

This was the first tool in the project and it still works: the commands below train and run it. No trained model ships with the package, so
train one first. It scores "secreted cell-surface glycoprotein", not adhesion, and none of the sorting-hat or AGA1-class tools depend on it
(background: `docs/TECHNICAL-NOTES.md`). The packaged training data (`data/positive`, `data/negative`) are FLO/ALS-related proteins from
*Saccharomyces* and *Candida* against random proteins.

```bash
surface_glyco_train --positive data/positive --negative data/negative      # writes ./models/*.pkl and a JSON model card
export SURFACE_GLYCO_MODELS_DIR=$PWD/models                                 # predict and evaluate look here, not in ./models
surface_glyco_predict --model models/surface_glyco_model_esm2_t6_8M_UR50D.pkl --input proteome.fasta --output scores.csv
surface_glyco_predict --model models/surface_glyco_model_esm2_t6_8M_UR50D.pkl --input query_dir --output combined.csv
```

`--input` may be one FASTA or a folder (`.fasta`, `.fa`, `.pep`, `.aa`, with or without `.gz`). The output column is
`surface_glycoprotein_score`; the 0.5 threshold is not calibrated. It runs on CPU and much faster on a GPU with CUDA. FungiDB raw
downloads need a login; get proteomes from NCBI Datasets or UniProt.

## Repository layout

| Path | Content |
|---|---|
| `src/cellsurface_sorting_hat/`, `src/surface_glyco/` | The two Python packages |
| `analysis/aga1_class/` | AGA1-class scripts (staging, tool arrays, converters, flags, classifier, summaries) |
| `analysis/cocci_repeats/` | Repeat detectors and the *Coccidioides* work they came from |
| `analysis/` (other folders) | Truth sets, calibration, scans and scoping studies, each with a README or report |
| `data/sorting_hat/` | Family table, Cys spacing sets, hydrophobin model, Phase C species |
| `data/controls/`, `data/curated/`, `data/positive`, `data/negative` | Control and truth sets, curated adhesin tables, training sequences |
| `scripts/sorting_hat/` | SLURM job scripts for the module tools |
| `docs/` | `TOOL-ARCHITECTURE.md`, `CLASSES.md`, `reports/` (dated results with data folders and hashes), `model-review/`, `superpowers/` (specs and plans) |
| `tests/` | One folder per component |

## Development and technical notes

Setup, CI, conventions, cluster notes and measured performance are in [`docs/TECHNICAL-NOTES.md`](docs/TECHNICAL-NOTES.md).
`CHANGELOG.md` lists breaking changes and `AGENTS.md` records AI-assisted contributions.

## Author

Jason Stajich, jason.stajich at ucr.edu. Code was developed with support from opencode.ai, Copilot and Claude models.
