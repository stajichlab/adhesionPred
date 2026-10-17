# Allergen predictor (evidence only): implementation plan

Date: 2026-10-10. Source: `docs/superpowers/specs/2026-10-09-allergen-predictor-spec.md` (draft 5 with
decisions A-G and the Q4 result). Status: draft plan, not reviewed, no code written.

Branch: a new feature branch, `allergen-predictor`, from `main`. No merges without the owner's go-ahead.

## Standing constraints

- Evidence only. Status `unvalidated`. Research use only. No threshold set before measurement.
- Repo Python with `/usr/bin/python3.12`. HPCC jobs sized 1-1.5 h. `${SCRATCH:?}`. No `/tmp`. No `BASH_SOURCE`.
- Publisher PDFs stay out of git (`_workdir/curation_texts_archive/`, hashes only).
- No specificity claim without the floors in spec section 6: 20 evidence-grade negative clusters for the 0.90 row,
  40 for the 0.95 row. **The IEDB count (Q4) shows zero usable protein-level negatives, so no sensitivity-at-specificity
  row is reportable on evidence-grade negatives.** Weak-negative numbers are labelled optimistic.

## Phase 0: gates and counts (no training)

T0.1 Count the non-fungal positive sources (AllerBase, SDAP 2.0, AlgPred2), with a scoping table in the same
     format as `analysis/allergen_scoping/`. Record licence terms. Output: `analysis/allergen_model/positives/`.
T0.2 Build the positive union by accession: IUIS (104) and KW-0020 fungal (111), 63 shared, one row per protein,
     source flags, `evidence_level`. Test: the union count is 152 for fungal alone; the script fails on a duplicate.
T0.3 Weak-negative candidates: reviewed Swiss-Prot with the spec section 3.2 filters (no KW-0020, no allergen
     comment, no "antigen", length 10-850, taxa 9606 and 10239 excluded, annotation score > 3, no "Cancer").
     Count only. Output: candidate table and counts per kingdom and genus.
T0.4 MMseqs2 cluster the positives and the candidate negatives on HPCC (`module load MMseqs2/17-b804f`, 30% identity,
     coverage 0.5). Record cluster counts for fungal and all positives. Job sized to 1 h.
T0.5 **Build gate (decision A).** Count pairable weak negatives at T_C 0.6 and 0.7. An arm with fewer than 500
     paired negatives is dropped. Also check that the length filter, pairing and length-matched run leave a
     feasible set together (M1L). Output: `analysis/allergen_model/gate/gate_counts.md`. **Stop point:** the owner
     reads the counts before any training.

## Phase 1: data build (after the gate)

T1.1 Split the positive clusters 70/15/15 (fixed seed). Test: no cluster in two splits.
T1.2 Pair each weak negative with a positive at identity >= T_C. A negative takes its paired positive's split.
     Unpaired candidates are dropped and logged in `rejected.tsv`. Pair-in-cluster is allowed (decision F).
T1.3 Remove T_S violations (T_S = 0.4, MMseqs2 all-vs-all, decision G(a)) across splits in both classes.
     Record the measured maximum inter-split identity. The split is valid only if that maximum is <= 0.4.
T1.4 Write the control outputs under `data/controls/allergen_model/`: `sequences.faa`, `controls.tsv` (evidence rows,
     14 columns, source and quote on each row), `training_negatives.tsv` (weak rows), `clusters.tsv`, `rejected.tsv`,
     `counts.md`, `unverified_notes.md`, `README.md` (records the schema deviation). Taxon IDs checked against
     `/srv/projects/db/taxonomy/names.dmp`. Large tables compressed (`.tsv.gz` shared, `.zst` intermediates).
T1.5 Tests: a fixture of 20 proteins with known clusters, pairs and violations; the build must reproduce the expected
     splits and reject the planted violations.

## Phase 2: baselines (no GPU)

T2.1 B1 composition (AAC, dipeptides, length), logistic regression. B4 length only. B3 6-mer linear SVM.
T2.2 B2 BLAST to IUIS, leave-species-out (drop the positive's own species from the reference). B6 nearest neighbour to
     training positives (training split only). Both are sanity floors. Expected: near chance on paired negatives.
     The expectation is written down before the run.
T2.3 X1 external predictor: first check that its server allows batch submission and that its terms permit this use.
     If not, X1 is dropped and the report says why. Disclosure of training overlap stays in the report.
T2.4 Each baseline scored on the primary test set, with cluster-bootstrap intervals (1,000 replicates; degenerate
     single-class replicates skipped and counted; rows with more than 5% skipped are not reportable).
T2.5 Tests: each baseline on a toy set with a known answer. A random score gives AUROC 0.5 within tolerance.

## Phase 3: embeddings and models (GPU, after Phase 1 and 2 pass)

T3.1 Embedding job: ESM-2 650M (`esm2_t33_650M_UR50D`), one job per protein set, sized 1-1.5 h, `exfab` queue,
     embeddings compressed to `.zst`. Mean pooling excludes padding; test: batch order does not change a score.
T3.2 M1 logistic regression on embeddings. Hyperparameters chosen on validation only.
T3.3 M1b random forest, same embeddings.
T3.4 M1L: M1 with negatives filtered to the positives' length range (only if T0.5 says it is feasible).
T3.5 Seeds: five per model. Label-permutation run per model, shuffling labels within each fixed split.
     A permutation run that scores well stops the work for review.
T3.6 Leave-family-out (PF16541 and other families): descriptive. Per-family retraining is not budgeted; the plan
     records the job count when it is known.

## Phase 4: evaluation and report

T4.1 Primary table: AUROC and sensitivity at specificity 0.90 (20 evidence-grade negative clusters) and 0.95
     (40 clusters). Under the current evidence the sensitivity rows are **not reportable**; AUROC, MCC and calls are.
     Weak-negative numbers carry the word "optimistic".
T4.2 AUPRC only beside its positive rate.
T4.3 Prevalence scan: score the Fungi_5k proteomes (56 scanned earlier). Report calls per proteome across a range of
     thresholds, not one point. No false positives are claimed.
T4.4 Label-noise check: count paired negatives that match a WHO/IUIS or AllergenOnline entry. Reported as a limit.
T4.5 Report under `docs/reports/`, with the banner. Every number traced to a file in `analysis/allergen_model/`.

## Phase 5: integration (only after the owner accepts Phase 4)

T5.1 New evidence module in `src/cellsurface_sorting_hat/modules/` and CLI wiring in `cli.py` (module names are
     hard-coded there, and `pfam.py` shows the pattern). Evidence field only. No change to existing calls.
T5.2 Tests for the module and the CLI. The sorting-hat test suite passes.

## Open items the plan depends on

- Non-fungal positive sources and their licences (T0.1).
- X1 batch access and terms (T2.3).
- Whether AllergenAI's accuracy (0.94 or 0.81) is right (needs the full text).
- DeepPlantAllergy table values (needs the PDF table).
- Allerdictor has no PDF; B3 stays unverified.
