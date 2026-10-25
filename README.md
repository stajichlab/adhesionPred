# adhesionPred: tools for fungal cell-surface and cell-wall proteins

This repository holds a set of tools and analyses for finding and sorting fungal cell-surface proteins: signal-peptide proteins,
GPI-anchored cell-wall proteins, tandem-repeat adhesins (flocculin/ALS/EPA-like), wall-family domain proteins and a few others.

There is **no single "adhesin predictor"**. "Adhesin", "antigen" and "biofilm factor" are outcomes, not sequence features, and the
evidence in this repository showed that one model for all of them is the wrong target (`docs/TOOL-ARCHITECTURE.md`). The
original ESM-2 classifier turned out to be a surface-glycoprotein detector with about 12% adhesin precision on *S. cerevisiae*, and no
trained model ships with the package. What you can use today are the components below. Each states what it measures and what it
does not.

## What is usable now

| Component | What it does | Status |
|---|---|---|
| [`cellsurface_sorting_hat`](#1-cellsurface_sorting_hat-sort-a-proteome-into-surface-categories) | Reads tool results for a proteome (SignalP, Pfam, repeat detectors, TMHMM, allergen BLAST) and writes per-protein calls with explicit evidence and a three-valued logic (called / not called / not assessable) | Engine, module wrappers, calibration commands and HPCC job scripts exist. Each call has a measured status per species where truth data allowed it (`docs/CLASSES.md`) |
| [AGA1-class tools](#2-aga1-class-tools-flocculin-like-sertthr-rich-gpi-wall-proteins) (`analysis/aga1_class/`) | Rule-based flags and a features-only classifier for Ser/Thr-rich, repeat-bearing, GPI-anchored wall proteins; scripts to run SignalP, PredGPI and AIUPred on whole proteomes and summarise them | Tested on 850 labelled proteins and on five yeast proteomes (reports below). Research tools, not a validated predictor |
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
`taxdump`) for the sorting hat. The SLURM scripts in `scripts/sorting_hat/` and `analysis/aga1_class/*.sbatch` use module names and
partitions of the UCR HPCC (`signalp/6`, `predgpi/202001`, `aiupred/2.x`); change them for your cluster.

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
  (`--min-disordered-run` sets the cut; 100 is the sensitivity value). This finds Ser/Thr-rich GPI wall mannoproteins in general
  (SED1/DAN4-like proteins and the ALS/FLO/EPA adhesins).
* **Narrow**: the broad flag **and** the project repeat call (detector 14: coverage at least 0.25 and at least 2.5 copies). Stricter;
  it misses adhesins whose repeats are not periodic (AGA1 and FLO11 have no repeat call).
* The Ser/Thr "stalk" rule is available (`stalk_features.py`) but is **not** used as a gate: it removed real positives and no negatives.

**The classifier** (`feature_classifier.py`): logistic regression, random forest and gradient boosting on 30 features (composition,
disorder, repeat score, signal peptide and GPI calls, tail hydrophobicity, length). With whole clusters held out, AUC is about 0.97 on
the labelled tier A positives against hard negatives and background, and at a 3.5% false-positive rate it recovers more positives
than the broad rule (about 0.77 to 0.85 against 0.66), but the confidence intervals overlap, so no improvement is claimed. Composition
carries nearly all the signal. Fixed thresholds do **not** carry across clades (the realised false-positive rate was 12 to 36% with
Saccharomycotina held out), so recalibrate a threshold for each new clade. Report: `docs/reports/2026-10-25-aga1-feature-classifier.md`.

**Run it on a proteome** (the pipeline used for the five-yeast test; all inputs are UniProt accessions):

```bash
# 1. stage chunks of at most 1,000 proteins (writes outside the repo; refuses an existing --out)
python3 -I analysis/aga1_class/stage_proteomes.py --out /path/NEW_DIR Scer=Scer.faa Calb=Calb.faa

# 2. run the tools, one SLURM array task per chunk (TOOL = signalp | predgpi | aiupred; the header of the script has the
#    exact sbatch lines, the pinned code commit and the GPU options for AIUPred)
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

Every script refuses to overwrite an existing output and stops with a `FATAL:` message (exit 1) on a missing accession, a count
mismatch or a malformed row; nothing is written as a silent 0. Per-residue AIUPred files carry residue letters, so they and FASTA
files are kept out of the repository (hashes are recorded instead).

**What has been measured** (`docs/reports/`):

* `2026-10-24-aga1-signalp-predgpi-results.md`, `2026-10-24-aga1-disorder-results.md`, `2026-10-24-aga1-broad-flag-candidate.md`: the predictor
  results on 850 labelled proteins and the first flag.
* `2026-10-25-aga1-feature-classifier.md`: the features-only study (above).
* `2026-10-25-aga1-proteome-test.md`: five yeast proteomes (*S. cerevisiae*, *C. albicans*, *C. glabrata*, *K. marxianus*, *C. auris*).
  In *C. glabrata* the broad rule, or two of three models with a signal peptide, finds 32 of 37 proteins with an adhesin-like UniProt
  name (most of the EPA family, AWP, PWP); in *K. marxianus* 11 of 14. The name reference is only a rough guide. All five are
  Saccharomycotina yeasts, so that test says little about other clades. A test on *A. fumigatus*, *Coccidioides*, *Cryptococcus*
  and *Rhodotorula* is under way and will have its own report.
* The flags do not see **CFEM** proteins (RBT5, PGA7, PGA10, CSA1, CSA2 and relatives): they are short GPI proteins without a long
  disordered run. They are a different class (heme acquisition) and are covered by the Pfam CFEM family (PF05730) in the family table.

## 3. `surface_glyco`: the older ESM-2 surface-glycoprotein scorer

This was the first tool in the project and it still works: the commands below train and run it. It is kept for comparison and for
surface-glycoprotein scoring, and none of the sorting-hat or AGA1-class tools depend on it. An ESM-2 embedding plus logistic regression that scores secreted cell-surface glycoproteins (FLO/ALS-like). It is **not** an adhesin
predictor and no trained model ships with the package (the models of version 0.1.0 were removed after review, issue #25). The
packaged training data (`data/positive`, `data/negative`) are FLO/ALS-related proteins from *Saccharomyces* and *Candida* against
random proteins, so the model learns "surface glycoprotein", mostly through repeats.

```bash
surface_glyco_train --positive data/positive --negative data/negative      # writes ./models/*.pkl and a JSON model card
export SURFACE_GLYCO_MODELS_DIR=$PWD/models                                 # predict and evaluate look here, not in ./models
surface_glyco_predict --model models/surface_glyco_model_esm2_t6_8M_UR50D.pkl --input proteome.fasta --output scores.csv
surface_glyco_predict --model models/surface_glyco_model_esm2_t6_8M_UR50D.pkl --input query_dir --output combined.csv
```

`--input` may be one FASTA or a folder (`.fasta`, `.fa`, `.pep`, `.aa`, with or without `.gz`). The output column is
`surface_glycoprotein_score` and the labels are `surface_glycoprotein` and `other`; the default 0.5 threshold is not calibrated.
The model card records the ESM model, layer, pooling and training data; `predict` and `evaluate` read them from the card and stop
with an error naming the field if it is missing or `--model-name` disagrees. Sequences are truncated at 1,022 residues and `predict`
prints how many. It runs on CPU and much faster on a GPU with CUDA. FungiDB raw downloads need a login; get proteomes from NCBI
Datasets or UniProt.

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

## Development

```bash
pre-commit install
ruff format . && ruff check .      # CI pins ruff 0.3.5; use that version, newer ruff formats differently
pytest tests/cellsurface_sorting_hat tests/aga1_class     # CI runs one pytest process per folder
```

Every analysis result is committed with the command that made it, the code commit, and sha256 hashes of inputs that stay outside
the repository. Proteins are keyed by NCBI or UniProt accession, never by JGI or internal ids. `CHANGELOG.md` lists breaking changes;
`AGENTS.md` records AI-assisted contributions.

## Author

Jason Stajich, jason.stajich at ucr.edu. Code was developed with support from opencode.ai, Copilot and Claude models.
