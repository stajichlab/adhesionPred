# Allergen predictor (evidence only): implementation plan, v3 (job scripts corrected after the fable review)

**Update 2026-10-10: pairing is dropped.** Negatives are length-matched to the positives (spec section 3.2). T1.2 pairing, the T_C arms and the pair-in-cluster exception no longer apply. B6 (nearest-neighbour) is a required shortcut check with the O1 stop margin (spec section 5). Pairing may be revisited if more positives raise the pairable count.

Date: 2026-10-10 (v2 same day, after the fable review of v1). Source: spec draft 5
(`docs/superpowers/specs/2026-10-09-allergen-predictor-spec.md`), decisions A-G, Q4 result
(`analysis/allergen_model/q4_iedb/counts_q4.md`). Status: v2, not yet reviewed again. No code written.

Branch: `allergen-predictor`, from `main`. No merges without the owner's go-ahead.

## Standing constraints

- Evidence only. Status `unvalidated`. Research use only. No threshold set before measurement.
- MMseqs2: use the non-AVX2 build `/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs`. The default module
  loads the AVX2 build, which fails with SIGILL on node c03 (job 29797130 and 29797276, 2026-10-10).
- Python: `/usr/bin/python3.12` for stdlib-only scripts. GPU and ESM jobs use the conda env `adhesionPred`
  (`/rhome/jstajich/.conda/envs/adhesionPred/bin/python`), which has torch and fair-esm. `/usr/bin/python3.12` has no torch.
  No `/tmp`. No `BASH_SOURCE` for paths.
- Every job script: `${SCRATCH:?}` for intermediates, trap cleanup on exit, copy to `/bigdata` under a
  temporary name, then `mv`. Jobs sized 1-1.5 h of real runtime (CLAUDE.md).
- Publisher PDFs stay out of git. T0.6 moves them to `_workdir/curation_texts_archive/` and commits hashes only.
- Q4 result: no usable protein-level negatives. Sensitivity-at-specificity rows are **not reportable**
  (floors 20 clusters for 0.90, 40 for 0.95 are not met by evidence-grade negatives). Weak-negative numbers are
  labelled optimistic.
- Each T_C arm (0.6 and 0.7) runs Phases 1-4 separately. Output paths are keyed by arm (`arm_0.6/`, `arm_0.7/`).
  Every table carries a `T_C` column. An arm below the 500-negative floor is dropped and reported.

## Phase 0: gates and counts (no training)

T0.1 Non-fungal positive sources (AllerBase, SDAP 2.0, AlgPred2): obtainability and licence. Record in
     `analysis/allergen_model/positives/sources.md`. No download until the licence is read.
T0.2 Positive union by accession: **done 2026-10-10** by `build/02_positives.py`. Union of 152 accessions (104 IUIS, 111 KW-0020,
     63 shared). Sequences from UniProt for every member. Result: **147 positive rows** in `analysis/allergen_model/positives/`
     (`positives.faa` with P_ IDs, `positives.tsv`, `counts.md`, `rejected.tsv`). Breakdown: 99 IUIS, 33 KW-0020 curated, 15
     similarity-only. Four accessions are inactive UniProt entries with no sequence, and are rejected with the reason. Three
     accessions were merged into current primaries; the mapping is recorded in `merged_from`. Reconciliation: 152 - 4 rejected
     - 3 merged + 2 new primaries = 147 (one new primary, Q8J0P4, was already in the union).
T0.3 Weak-negative candidates: **done 2026-10-10** by `build/01_candidates.py` (UniProt release 2026_03). Reviewed Swiss-Prot
     fungi with the spec 3.2 filters, and annotation score 4 or 5 (the UniProt range syntax is rejected, so the query uses
     `annotation_score:4 OR annotation_score:5`). Result: **6,627 candidates from 192 genera**, in
     `analysis/allergen_model/candidates/` (`candidates.faa` with N_ IDs, `candidates.tsv`, `counts.md`). The earlier figure of
     10,747 lacked the annotation filter and is no longer used.
T0.4 Positive and candidate **clustering**: `jobs/cluster_mmseqs.sbatch` (MMseqs2 easy-cluster, 30% identity, coverage 0.5,
     `-s 7.5`), partition `stajichlab`, 16 cpus, 32 GB, 6 h. Before the full run, `jobs/pilot_mmseqs.sbatch` times clustering and the
     search on a random sample of 1,000 candidates (fixed seed) and writes `pilot_projection.txt`. The projection is linear, so
     its time limits are a floor and an owner-checked estimate, not a guarantee. Output: `clusters.tsv` in the COMMON-RULES schema.
     Status: the pilot script dry-ran on the fixture (`fixtures/dryrun_out/out/pilot`). It has not run on real candidates.
     Runtime is not yet measured; the job is sized from a timed pilot on a subset before the full run.
T0.5 **Pairwise identity table (new in v2).** `mmseqs search` (or `easy-search`) of every weak-negative candidate
     against the positive set. Parameters recorded: sensitivity, coverage mode, `--min-seq-id 0`, output columns
     including `pident`, `alnlen`, `qcov`, `tcov`. Output: `pairs_candidates_vs_positives.tsv.zst`. This table is the
     only input to T0.6 and T1.2. Job sized from the pilot in T0.4.
T0.6 **Build gate (superseded 2026-10-10: pairing dropped; the gate counts length-matched negatives per split, planning floor 500).** [Previous text:] Count pairable candidates at T_C 0.6 and 0.7 from the T0.5 table. Arms with fewer
     than 500 paired negatives are dropped. Check that the length filter, pairing and the length-matched run (M1L)
     leave a feasible set together. Output: `analysis/allergen_model/gate/gate_counts.md`. **Stop point:** the owner
     reads the counts before Phase 1. Publisher PDFs are moved to the archive here (no commit yet).

## Phase 1: data build (per surviving arm, after the gate)

T1.1 Split positive clusters 70/15/15, fixed seed, at cluster level. Test: no cluster in two splits.
T1.2a Positive-only pass first (spec section 4, step 3). Positive pairs that cross splits at fident >= 0.4 (T_S) are resolved
     by the T1.3 removal rule, applied to positives only, before any negative is paired.
T1.2 (superseded 2026-10-10) Length-matched sampling replaces pairing: negatives sampled to the positives' length distribution, assigned to splits by cluster. Previous pairing text:
     fident >= T_C, qcov >= 0.5 and tcov >= 0.5. The negative takes the split of its paired positive. Unpaired candidates go
     to `rejected.tsv` with counts per condition. Pair-in-cluster is allowed (decision F; task 05 carries the dated exception).
T1.3 Inter-split identity check. Run `cross_split_check.sbatch` (all-vs-all, `--max-seqs 2000`, fident, T_S = 0.4, decision
     G(a)). The job fails if fewer than two split prefixes are present, so it cannot pass vacuously. Record the measured maximum.
     Removal rule, adopted here and marked for owner confirmation: for each crossing pair at identity >= 0.4, remove the
     protein from the smaller split (ties: remove the weak negative before a positive; then the higher accession).
     Re-measure after removal. Repeat until no crossing pair is >= 0.4. The rule is deterministic and logged in `rejected.tsv`.
T1.4 Recount after removal. Per split: positive clusters, negative clusters (evidence-grade and weak counted separately),
     `in_reference` split. A stratum that fails the fungal floor (20 independent test clusters) is marked
     **not reportable** in the tables that use it.
T1.5 Control outputs under `data/controls/allergen_model/arm_<T_C>/`: `sequences.faa`, `controls.tsv` (evidence rows, 14
     columns, source and quote on each row), `training_negatives.tsv` (weak rows, not evidence-grade), `clusters.tsv`,
     `rejected.tsv`, `counts.md`, `unverified_notes.md`, `README.md` (records the schema deviation and the banner). Each
     row carries `in_reference`. Taxon IDs checked against `/srv/projects/db/taxonomy/names.dmp` (T0.7 checks that this file
     exists and its date). Large tables compressed (`.tsv.gz` shared, `.zst` intermediates).
T1.6 Tests: a 20-protein fixture with planted clusters, pairs and crossing violations. The build must reproduce the expected
     splits, remove the planted violations, and log them. A second fixture with a planted leak (a training positive in test)
     must fail the build.
T1.7 Provenance record (spec section 8): data hashes, MMseqs2 version and parameters, aligner parameters for T1.3, ESM-2
     variant, seeds, T_C, commit. Written by each job, not by hand.

## Phase 2: baselines (no GPU, per arm)

T2.1 B1 composition (AAC, dipeptides, length), logistic regression. B4 length only, logistic regression. B3 6-mer linear
     SVM (Allerdictor design, unverified). B5 one-hot + linear head (the encoding floor of DeepPlantAllergy and Applm).
T2.2 B2 BLAST to IUIS, leave-species-out (drop the positive's own species from the reference). B6 nearest neighbour to
     training positives, training split only.
T2.3 X1 external predictor: **dropped** (owner decision O5, 2026-10-10). Reason, for the report: AllerCatPro 2.0 caps batch
     submissions at 50 sequences and its page states no terms of use; the AlgPred 2.0 host did not resolve. No external
     reference is run. The report states this.
T2.4 Each baseline on the primary test set. Cluster-bootstrap intervals, 1,000 replicates. Degenerate single-class
     replicates are skipped and counted; a row with more than 5% skipped is not reportable.
T2.5 **Pass criterion and stop point (new in v2).** Written down before the run: B2 and B6 are near chance on paired negatives.
     Stop if B2 or B6 AUROC is above 0.60 (owner decision O1), or if B4 AUROC is within 0.05 of M1 AUROC once M1 is run in
     Phase 3. A stop halts the arm for owner review.
T2.6 Controls on every baseline: five seeds where the learner is stochastic, and one label-permutation run, shuffling labels
     within each fixed split. A permutation run that scores well stops the arm.
T2.7 Toy tests: a random score gives AUROC 0.5 within tolerance; a planted separable feature gives AUROC near 1.

## Phase 3: embeddings and models (GPU, after Phase 1 and 2 pass, per arm)

T3.1 GPU queue: `exfab`, if the account has access (owner decision O4). A pre-check with `sacctmgr` runs before the first
     GPU job and stops the arm if the account is not in the group. Embedding job:
     ESM-2 650M (`esm2_t33_650M_UR50D`), one job per protein set, sized to 1-1.5 h. Embeddings compressed to `.zst`.
     Test: mean pooling excludes padding, and batch order does not change any score.
T3.2 M1 logistic regression on embeddings. Hyperparameters on validation only.
T3.3 M1b random forest, same embeddings.
T3.4 M1L: M1 with negatives filtered to the positives' length range (only if T0.6 says feasible).
T3.5 Five seeds per model; one label-permutation run per model within fixed splits.
T3.6 Score table `allergen_model.tsv`: protein ID, score, model ID, seed, T_C, split, training-data hash.
T3.7 Leave-family-out (PF16541 and others): descriptive. Per-family retraining is not budgeted; the job count is recorded
     once the families are known.

## Phase 4: threshold and evaluation

T4.1 Threshold (owner decision O2, 2026-10-10): the lowest validation score that gives sensitivity >= 0.80 on validation
     positives. Chosen on validation only, before any test result is read. Recorded as a config value with its rule and date.
T4.2 Primary table: AUROC, MCC, and calls. Sensitivity at specificity rows are **not reportable** under the Q4 result; the table
     says so. Every table is stratified by all positives, fungal positives, leave-family-out, and `in_reference` vs
     `out_of_reference`. Weak-negative numbers carry "optimistic".
T4.3 AUPRC only beside its positive rate.
T4.4 Prevalence scan: score the Fungi_5k proteomes (56 scanned earlier, per the handoff). Calls per proteome across a range
     of thresholds, not one point. No false-positive claim.
T4.5 Label-noise check: count paired negatives matching a WHO/IUIS entry. A matching AllergenOnline check runs only if that
     file is obtained (T0.1); otherwise the report says the check is WHO/IUIS only.
T4.6 Report under `docs/reports/`, with the banner. Every number traces to a file under `analysis/allergen_model/`.

## Phase 5: integration (only after the owner accepts Phase 4)

T5.1 New module file `src/cellsurface_sorting_hat/modules/allergen_model.py` with module ID `allergen_model`. The existing
     `allergen.py` is the BLAST homology module and is not changed. Module names are hard-coded in
     `src/cellsurface_sorting_hat/modules/cli.py` (the `ModuleSpec` entries), so CLI wiring is part of this task.
T5.2 Evidence field only. No change to `iuis_allergen_similarity`, `iuis_allergen_homolog` or any category.
T5.3 Tests for the module and the CLI. The sorting-hat test suite passes.

## Compute annex (v3, 2026-10-10)

Partitions checked with `sinfo` on 2026-10-10. `stajichlab`: 30-day limit, 64+ cpus per node. `exfab`: 30-day limit,
one node with 2x ada6000. `short`/`short_gpu`: 2 h cap. The account has the `exfab` association (`sacctmgr`).
Job scripts are in `analysis/allergen_model/jobs/`. Each follows the repo's sbatch convention: required variables
checked, `${SCRATCH:?}`, trap cleanup, outputs written to a temporary name then `mv`, and a `.done` marker written last.

| Task | Script | Partition | CPUs | Memory | Time | Input size | Notes |
|---|---|---|---|---|---|---|---|
| T0.4 cluster | `cluster_mmseqs.sbatch` | stajichlab | 16 | 32 G | 6 h | 10,747 candidates (upper bound) + 152 positives | Pilot on 1,000 random candidates first; time the real run from it |
| T0.5 pairs | `pairs_search.sbatch` | stajichlab | 16 | 32 G | 6 h | same | All candidate-vs-positive pairs; unfiltered (-e 1000) |
| T1.3 cross-split | `cross_split_check.sbatch` | stajichlab | 16 | 32 G | 4 h | per arm, up to ~11,000 proteins | All-vs-all; removal is done by the build, not here |
| T3.1 embeddings | `embed_esm2.sbatch` | exfab | 8 | 64 G | 4 h | per protein set | 1 GPU; ESM-2 650M; pre-check on the `exfab` association |
| T2.x baselines | (per-task Python, `stajichlab`) | stajichlab | 4 | 16 G | 2 h | per arm | CPU only; no GPU |
| T3.2-3.5 models | (per-task Python, `stajichlab`) | stajichlab | 8 | 32 G | 2 h | per arm | Five seeds and one permutation per model, bundled into one job per model |
| T4.4 scan | (per-task Python, `stajichlab`) | stajichlab | 8 | 16 G | 2 h | 56 Fungi_5k proteomes | Bundle proteomes into jobs of about 1 h |

**Job count (recounted, v3 review).** Shared across both arms, run once: the pilot, the clustering job (T0.4), the pairs
search (T0.5), the embedding jobs for the Fungi_5k proteomes, and the proteome scan. Per arm: one cross-split job after the
build (T1.3), one more if the removal rule changes the sets, the baseline job (one job for B1 to B5, B6 included), three
embedding jobs (train, validation, test) that are shared with the other arm only where the protein set is identical,
and three model jobs (M1, M1b, M1L), each running five seeds and one permutation. Per arm that is about 7 to 8 jobs.
The job count is not checked against any submission-rate limit. CLAUDE.md mentions such limits but gives no number.
Cluster and pairs jobs run for the measured time from the pilot, not the 6 h allocation, which is an upper bound.

**Acceptance criteria for every job script (moved from the standing constraints):** required variables checked and the
job stops if they are unset; `${SCRATCH:?}` for intermediates; `trap` cleanup on exit; outputs by temporary name then
`mv`; `.done` marker written last and removed at start; `/usr/bin/python3.12` for any Python; no `/tmp`; no
`BASH_SOURCE`; a FATAL message and nonzero exit on each failure path.

**Provenance (T1.7) fields** that every job writes: MMseqs2 version (from `mmseqs version`), parameters, ESM-2 model
name, seed, arm (T_C), input file hash (`sha256sum`), and the git commit of the repo.

## Open decisions for the owner

O1. (Decided 2026-10-10: B2/B6 halts when the mean validation AUROC over four split seeds (20261010, 1, 2, 3) is above 0.60; B4 within 0.05 AUROC of M1 halts. Amended 2026-10-10 from a single-seed rule after seed re-runs: validation means 0.577 (sd 0.044); test means 0.530 (sd 0.042). Recorded in T2.5.)
O2. (Decided 2026-10-10: validation sensitivity 0.80. Recorded in T4.1.)
O3. (Decided 2026-10-10: ships with the T4.1 threshold, labelled `unvalidated`, with the banner and the optimistic
    weak-negative limit stated in the output. Recorded in T5.2.)
O4. GPU access: confirm `exfab` access for this account, or name the partition to use (`short_gpu` caps at 2 h; the shared
    `gpu` partition had ~10 h fairshare waits per the handoff).
O5. X1 external tool: whether its batch use is permitted (T2.3). Not checked.
O6. (Decided 2026-10-10: the proposed rule is confirmed. Recorded in T1.3.)

## Open dependencies

- Non-fungal positive sources and licences (T0.1). AllerBase, SDAP 2.0, AlgPred2 not yet checked.
- Full texts for AllergenAI (0.94 vs 0.81 accuracy) and DeepPlantAllergy table values.
- Allerdictor has no PDF; B3 stays unverified.
- MMseqs2 runtime and memory are unmeasured; T0.4 sizing depends on the pilot.
