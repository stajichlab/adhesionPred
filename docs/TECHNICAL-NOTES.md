# Technical notes

Details that do not belong in the README: development setup and CI, conventions the scripts follow, cluster-specific notes, details
of the older ESM-2 scorer, and a summary of measured performance with pointers to the reports that hold the data. The README describes
what the tools are and how to use them; this file records how they are built and what has been measured.

## 1. Development setup and CI

```bash
pip install -e ".[dev]"
pre-commit install
ruff format . && ruff check .
pytest tests/cellsurface_sorting_hat tests/aga1_class
```

* **Ruff version.** CI pins `ruff==0.3.5`. A newer ruff formats some files differently, so a clean `ruff format --check` with a newer
  version does not mean CI will pass. Lint with the pinned version, for example
  `uv run --no-project --with ruff==0.3.5 ruff format --check .`.
* **Tests.** CI runs one pytest process per folder under `tests/` (each folder's `conftest.py` puts its own `analysis/` directory on the
  path). Tests use synthetic data; a few record tests check the hashes of committed result tables against their README. Running a
  folder without the project installed works with `PYTHONPATH=src` and the packages the folder imports, for example
  `PYTHONPATH=src uv run --no-project --with pytest --with pandas --with numpy --with scikit-learn --with matplotlib python -m pytest tests/aga1_class`.
  `tests/surface_glyco/test_config_paths.py` needs the `esm` package.
* **Generated files.** `docs/CLASSES.md` is generated from `docs/class_descriptions.yaml` and
  `src/cellsurface_sorting_hat/categories.yaml` by `analysis/class_table/make_classes_md.py`; a test checks that it is current.
  `AGENTS.md` records AI-assisted contributions.

## 2. Conventions the scripts follow

* **Accessions.** Proteins are keyed by NCBI or UniProt accession, never by JGI or other internal ids.
* **No silent zeros.** A value that does not apply is `NA` or `not_applicable`, never 0. Examples: PredGPI does not score proteins of 40
  residues or fewer, and the repeat detectors skip proteins under 80 residues; both give a missing value.
* **Fail loudly.** Scripts that write results refuse to overwrite an existing output and stop with a `FATAL:` message and exit code 1 on a
  missing or extra accession, a count mismatch or a malformed row.
* **Sequences stay out of the repository.** FASTA files and per-residue tool output (which carries residue letters) are kept outside
  the repository; the record README holds their sha256 instead. Downloaded data is read with `python3 -I`.
* **Records.** Every committed result has the command that made it, the code commit, the tool versions (from a provenance file written by
  the job) and a README with sha256 hashes that a test checks.
* **Cluster jobs.** The `*.sbatch` files check that the analysis directory of the cluster clone matches the pinned `CODE_COMMIT`, write to a
  new output directory, use `$SCRATCH` for temporary files, and verify input hashes and per-chunk row counts before copying results.
  The scripts are committed unedited after they have run.
* **Frozen decisions.** Family sign-offs (`narrow_ok`), thresholds and label sets are recorded with the date and owner in the data tables or
  report that holds them, and are not changed by a later run.

## 3. Cluster notes

The SLURM scripts in `scripts/sorting_hat/` and `analysis/aga1_class/*.sbatch` were written for the UCR HPCC and use its module names
(`signalp/6`, `predgpi/202001`, `aiupred/2.x`) and partitions (a 2-hour `short` CPU partition; a GPU partition named on the command
line for AIUPred). Change them for another cluster.

* PredGPI needs `PREDGPI_HOME` set explicitly in a non-interactive shell.
* AIUPred needs torch; on a GPU it took about 75 seconds for 850 proteins. IUPred3 on CPU took 21 minutes for 850 proteins.
* Measured times for 850 proteins: SignalP 6 (8 CPUs, fast mode) 12 min 50 s; PredGPI (2 CPUs) 2 min 55 s. In the proteome arrays (chunks of
  up to 1,000 proteins) a task took up to 29 minutes for SignalP, up to 1 min 37 s for PredGPI and 10 to 40 seconds for AIUPred on a GPU.
* Proteome runs use one SLURM array task per chunk of at most 1,000 proteins (`stage_proteomes.py`); the array throttle (`%N`) keeps the run
  within a group CPU limit, and a group at its CPU limit leaves jobs pending with reason `AssocGrpCpuLimit` regardless of the partition.
* AIUPred and IUPred3 scores overshoot 0..1 slightly (up to 1.023 and down to -0.012 were seen); the converter accepts 0.05 outside
  the range and uses the values as written.

## 4. The older ESM-2 scorer (`surface_glyco`)

* **Model card.** Training writes a JSON card next to the model with the ESM model, layer, pooling and training data. `predict` and
  `evaluate` read those values from the card and stop with an error naming the field if the card is missing, the pooling is not supported or
  `--model-name` differs from it. Cards of new models also record library versions.
* **Model location.** `surface_glyco_train` writes to `./models`. `predict` and `evaluate` look in the packaged model directory or in
  `SURFACE_GLYCO_MODELS_DIR`, not in `./models`.
* **Truncation.** Sequences longer than 1,022 residues are truncated; `predict` prints how many. Embeddings average over residue tokens
  only, so a score does not depend on the other sequences in its batch.
* **Why no model ships.** The two models of version 0.1.0 were removed (issue #25). They were trained on FLO/ALS homologs against random
  proteins with padding-inclusive pooling. The model review (`docs/model-review/`) found the classifier behaves as a surface-glycoprotein
  detector (about 12% adhesin precision on *S. cerevisiae* S288C) and that amino-acid composition alone separates the labelled sets
  (leave-family-out ROC-AUC 0.995), so cross-validation on that data is saturated.
* **Output.** `surface_glycoprotein_score` and the labels `surface_glycoprotein` and `other`; the default 0.5 threshold is not calibrated.
  `docs/TOOL-ARCHITECTURE.md` explains why the output must not be reported as an adhesion probability.

## 5. Measured performance (summary; the reports hold the data)

Everything here is research-grade and comes with the limits stated in the cited reports. Counts are proteins.

### Sorting hat
Per-call and per-family statuses (`unvalidated`, `smoke`, `estimated`) and the species they hold for are in `docs/CLASSES.md` and
`data/sorting_hat/family_table.tsv`. Family specificity screens: `docs/reports/2026-10-24-pr6-family-specificity-results.md`.

### AGA1-class tools (labelled set of 850 proteins)
Reports: `docs/reports/2026-10-24-aga1-class-scoping.md`, `...-aga1-signalp-predgpi-results.md`, `...-aga1-disorder-results.md`,
`...-aga1-broad-flag-candidate.md`, `2026-10-25-aga1-feature-classifier.md`, `2026-10-25-aga1-proteome-test.md`.

| Group | n | Broad flag | Narrow flag |
|---|---|---|---|
| tier A positives | 35 | 23 | 14 |
| tier B | 32 | 18 | 4 |
| hard negatives | 172 | 8 | 2 |
| background assumed negatives | 389 | 7 | 0 |

(Broad: signal peptide and GPI and an AIUPred run of at least 200 residues; narrow adds the detector-14 repeat call. A run cut of 100
gives 26, 23, 23 and 31 for the broad flag.) The Ser/Thr stalk rule (a 30-residue window at least 50% Ser/Thr inside the run) changed
tier A from 23 to 18 and removed no negative, so it is an option and not a gate. Median longest disordered run: tier A 302, hard
negatives 59, background 16 residues.

**Features-only classifier** (leave-cluster-out, 5 folds x 10 seeds, cluster-bootstrap 95% intervals):

| Model | ROC-AUC | Recall at 1% FPR | at 3.5% FPR | at 5% FPR |
|---|---|---|---|---|
| Logistic regression | 0.973 [0.936, 0.992] | 0.62 | 0.85 | 0.87 |
| Random forest | 0.975 [0.945, 0.990] | 0.55 | 0.77 | 0.85 |
| Gradient boosting | 0.965 [0.913, 0.988] | 0.59 | 0.79 | 0.87 |
| Broad rule | | 0.66 at FPR 2.7% | | |
| Narrow rule | | 0.40 at FPR 0.4% | | |

The intervals overlap those of the rules, so no improvement is claimed. Composition alone gives AUC about 0.97; predictor calls alone
0.61 to 0.70; N-glycosylation site density carries no signal. Tier B recall at 3.5% FPR is 0.66 to 0.71 for the models against 0.56 for the
broad rule. With Saccharomycotina held out, ranking carries over (AUC 0.87 to 0.97) but a fixed 3.5% threshold gives a realised false-positive
rate of 12 to 36%, so a threshold must be recalibrated per clade.

**Five yeast proteomes** (27,811 proteins; UniProt reference proteomes of *S. cerevisiae*, *C. albicans*, *C. glabrata*, *K. marxianus* FIM1
and *C. auris*): the broad rule flags 11, 25, 25, 10 and 14 proteins and the narrow rule 5, 6, 4, 2 and 5. Against UniProt adhesin-like
names (a rough, partly homology-derived reference) the broad rule or two of three models with a signal peptide finds 32 of 37 in
*C. glabrata* and 11 of 14 in *K. marxianus*. A model score without a signal-peptide condition is not specific: high-scoring proteins
without one are cytoplasmic disordered proteins. All five are Saccharomycotina yeasts.

**Four non-yeast proteomes and *C. parapsilosis*** (frozen models and thresholds, nothing retuned; report
`docs/reports/2026-10-26-aga1-nonyeast-and-cpar-test.md`):

| Species | Proteins | Broad rule | Narrow rule | logreg / rf / gboost calls | 2 of 3 models |
|---|---|---|---|---|---|
| *A. fumigatus* | 9,647 | 7 | 3 | 19 / 15 / 85 | 20 |
| *C. immitis* | 9,712 | 8 | 2 | 28 / 23 / 87 | 27 |
| *C. neoformans* | 7,427 | 2 | 0 | 2 / 5 / 70 | 3 |
| *R. toruloides* | 8,160 | 2 | 1 | 6 / 19 / 117 | 11 |
| *C. parapsilosis* | 5,777 | 28 | 9 | 26 / 29 / 35 | 25 |

The *C. immitis* spherule outer wall glycoprotein (CIMG_04613) passes both rules and all three models; related SOWgp proteins are training positives.
In *C. parapsilosis* the broad rule finds 16 of 37 adhesin-named proteins and the broad rule or 2 of 3 models with a signal peptide finds 18 (ALS1, 3, 6,
7, 11, HYR31, Hyr-like and Flo11-domain proteins). Thresholds are not equally calibrated outside yeast: gradient boosting calls 70 to 117 proteins in the
four non-yeast proteomes against 2 to 28 for logistic regression and 5 to 23 for the random forest, so a single-model call needs recalibration on a new clade;
the "2 of 3 models with a signal peptide" call absorbs most of it. Enzymes with a Ser/Thr-rich linker (chitinases, an endoglucanase) are still flagged by the
broad rule, and secreted Ser/Thr/Pro-rich proteins without a GPI call are scored high by the models but excluded by the broad rule.

**CFEM proteins.** Of the 15 UniProt-annotated CFEM (PF05730) proteins in the five yeast proteomes, all have a signal peptide, 14 have a GPI call and none
is flagged by the rules or the models, as expected for short GPI proteins without a long disordered run. In the non-yeast proteomes 4 of 20 are flagged
(*A. fumigatus* cfmA and an extracellular Ser/Thr-rich protein; *C. immitis* CIMG_00795 and CIMG_03472), and in *C. parapsilosis* 1 of 8 (CFEM5, 918 aa,
a disordered run of 230). The flags miss the short CFEM proteins and catch some long CFEM-domain proteins that carry a long disordered region; the Pfam CFEM
family covers the family.

### Cross-proteome aggregator
`cellsurface_sorting_hat_aggregate` (`src/cellsurface_sorting_hat/aggregate/`) reads `calls.long.tsv.gz`, `proteins.tsv.gz` and `run.json` of
each run at the default gate variant (`R0`), one taxon per run, and refuses mixed `categories.yaml` hashes unless `--allow-mixed-config`.
Status of a (proteome, call) cell is the weakest status among the called records, else `not_called`. A cell with fewer than half of its
proteins assessable is flagged `low_assessable` and left out of the group rate summaries (an OR with an unassessable input makes the rate
1,000 per 1,000 by construction). Figures are deterministic PNGs. First report on 10 proteomes:
`docs/reports/data/aga1_aggregate_2026-10-26/` (limits in its README). Not yet done: orthogroups, diversity and duplication, genome
context, phylogenetic statistics, a scale-out runner, an interactive page.

Second report with the Pfam, TMHMM and repeat02 modules added (`analysis/sorting_hat_run/aga1_domain_modules.py`):
`docs/reports/data/aga1_aggregate_domains_2026-10-26/`. Remaining NA calls: those that need the allergen BLAST and antigen lookup.

