# Adhesion Protein Predictor

A bioinformatics tool that scores proteins from FASTA files as **secreted cell-surface glycoproteins**
(FLO/ALS-like). Despite the project name, this is step 1 of a multi-step pipeline and not an adhesin
predictor. The removed 0.1.0 model made calls on *S. cerevisiae* S288C of which about 12% were known adhesins. Mechanism-specific
adhesin classes are separate tools; see `docs/TOOL-ARCHITECTURE.md` section 2.0 for the steps and names.

# Requirements

This will run with pytorch on CPUs but will be much faster if on a GPU system with torchvision installed having CUDA bindings.

# Usage

### Training

No trained model ships with this package. Train one before you predict. The training data in `data/` are FLO, ALS1 related proteins from Saccharomyces and Candida. The model that this data produces scores surface glycoproteins, not adhesins. A model trained on a surface-glycoprotein label is planned (see `docs/PLAN-2026-09-30-pipeline-and-decisions.md`).

```
surface_glyco_train --positive data/positive --negative data/negative
```

Training writes the model to `./models` by default. It also writes a JSON model card next to the model. `surface_glyco_predict` and `surface_glyco_evaluate` do not look in `./models`. They look in the packaged model directory, or in the directory in `SURFACE_GLYCO_MODELS_DIR`. After training, pass the model with `--model`, as in the examples below, or run `export SURFACE_GLYCO_MODELS_DIR=$PWD/models`.

### Application 

You can run on a single file at a time and produce a report for each query file.
Note: FungiDB raw-file downloads now require a login (HTTP 401 as of 2026-09); fetch
proteomes from NCBI Datasets or UniProt, or download from FungiDB while logged in.
```
mkdir -p query
pushd query
curl -O https://fungidb.org/a/service/raw-files/release-68/CalbicansSC5314/fasta/data/FungiDB-68_CalbicansSC5314_AnnotatedProteins.fasta
curl -O https://fungidb.org/a/service/raw-files/release-68/CneoformansJEC21/fasta/data/FungiDB-68_CneoformansJEC21_AnnotatedProteins.fasta
curl -O https://fungidb.org/a/service/raw-files/release-68/Spombe972h/fasta/data/FungiDB-68_Spombe972h_AnnotatedProteins.fasta
popd
for qorg in $(ls query/*.fasta)
do
   surface_glyco_predict --model models/surface_glyco_model_esm2_t6_8M_UR50D.pkl --input $qorg --output $(basename $qorg .fasta).surface_glyco.csv
done
```

You can run on a single folder and all results will be combined in a single file. It will look for all .fasta, .fa, .pep, .aa with or without .gz extensions.

```
surface_glyco_predict --model models/surface_glyco_model_esm2_t6_8M_UR50D.pkl --input query --output Combinedquery.surface_glyco.csv
```

## Development Setup

1. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Setup Pre-commit Hooks:**
   ```bash
   pre-commit install
   ```

Training writes a JSON model card next to the model (`surface_glyco_model_<esm>.json`) that records the
ESM model, layer, pooling and training data. `surface_glyco_predict` and `surface_glyco_evaluate` read
the ESM model, layer and pooling from the card. They stop with an error that names the field if the card is
missing, the pooling is not supported, or `--model-name` differs from the card. `surface_glyco_predict`
prints how many sequences it truncates at 1022 residues. The output column is `surface_glycoprotein_score`
and the labels are `surface_glycoprotein` and `other`. The default 0.5 threshold is not calibrated. See `docs/model-review/` for the current model review and validation plan.

## cellsurface_sorting_hat (in development)

`cellsurface_sorting_hat` sorts the proteins of a proteome into cell surface categories. The
current calls are `signal_peptide_protein`, `tandem_repeat_protein`, `wall_family_domain`,
`cocci_specificity_rank_top15`, `serodiagnostic_marker_candidate`, `iuis_allergen_similarity`,
`iuis_allergen_homolog` and `cell_wall_adhesion_candidate`. Design: `docs/superpowers/specs/2026-10-04-orchestrator-design.md`.

The core engine reads module result tables from a work directory:

    cellsurface_sorting_hat --fasta proteome.faa --taxon TAXON_ID --taxdump nodes.dmp \
        --workdir work --out out

`TAXON_ID` is the NCBI taxon ID of the proteome. Look it up in the NCBI taxonomy.

It writes `calls.long.tsv.gz`, `calls.wide.tsv.gz`, `evidence.tsv.gz`, `proteins.tsv.gz`,
`report.md` and `run.json`.

The module wrappers and the calibration commands now exist:

- `cellsurface_sorting_hat_module` turns the output of one tool (SignalP, Pfam `hmmsearch`, the
  repeat detectors, BLASTP against the IUIS allergens, TMHMM) or one lookup table (antigen ranking,
  Cys-rich tiers, spherule expression) into a module table in the work directory.
- `cellsurface_sorting_hat_calibrate` writes status sources from a truth set (`truth`), from the
  Phase C metrics (`phasec`) and runs the Pfam specificity, allergen leave-species-out and panel
  checks.

No calibration has been run on HPCC data yet (Tasks 11 to 15 of the plan are not run). Every module
status is therefore `unvalidated`. The Phase C numbers for rule R0 exist, and `calibrate phasec` can turn them into status entries, but no status file has been written from them in a real work directory.

# Author

Jason Stajich, jason.stajich<at>ucr.edu 

Code was developed with support from opencode.ai, co-pilot and various code models for setting up classifier framework

