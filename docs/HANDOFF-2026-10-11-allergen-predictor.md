# Handoff: allergen predictor (evidence only), 2026-10-11

Read first: `docs/reports/2026-10-10-allergen-predictor-status.md` (every result with cluster counts and intervals),
then `docs/superpowers/specs/2026-10-09-allergen-predictor-spec.md` (draft 5 plus amendments) and
`docs/superpowers/plans/2026-10-10-allergen-predictor-plan.md` (v3 plus updates).

Status: research only. Nothing committed. Nothing merged. No specificity claim is supported anywhere.
Branch at handoff: `fungi-scan-and-specs` (new files are untracked). Plan says to use `allergen-predictor` for code;
no code branch has been created yet.

## 1. What the model is

Evidence-only allergen score for `cellsurface_sorting_hat`. Model M1: ESM-2 650M mean-pooled embedding (padding excluded,
layer 33), logistic regression, balanced classes. Status `unvalidated`. Not a call. Not a mechanism category.

## 2. Decisions in force (in order made)

- Owner lifted the "no allergenicity model" rule for an evidence-only model (2026-10-09). The rule came from an agent
  directive on 2026-10-06, not from the owner.
- Positives: all kingdoms, fungal results reported separately (Q1).
- Weak negatives come from reviewed Swiss-Prot with no allergen annotation (Q2). Evaluation negatives need an independent
  reason. None is usable yet (IEDB IgE-negative records are all allergens or peptide-level, see
  `analysis/allergen_model/q4_iedb/counts_q4.md`).
- Pairing (Applm-style, T_C 0.6 and 0.7) is DROPPED (2026-10-10). Only 74 candidates pair at 0.6 and 39 at 0.7, against a
  500 floor. Replaced by length-matched sampling. Can be reopened if more positives raise the pairable count.
- Candidates at fident 0.9 or higher to a positive were removed (6 proteins). O7 coverage rule kept only for the
  label-noise check.
- Split: cluster level, 70/15/15. No cluster is split across sets. Fungal and all-kingdom sets each use one shared clustering.
- Specificity floor: 20 evidence-grade negative clusters for a 0.90 row, 40 for 0.95 (D, G(b)). Planning floor for negatives
  per split: 500 (A). Both are unmet, so no specificity claim.
- Stop rule (O1, amended 2026-10-10): halt if the MEAN validation B6 AUROC over four split seeds exceeds 0.60 (seeds
  20261010, 1, 2, 3). Current mean 0.577, so no halt.
- Threshold rule (O2): lowest validation score with sensitivity >= 0.80. Ship with threshold, labelled unvalidated (O3).
  Not yet computed.
- GPU queue `exfab` (O4). Embedding script checks the account association first.
- External reference tool X1 dropped: AllerCatPro caps batches at 50 and states no terms of use.
- Owner chose to add AllergenOnline (fungal and non-fungal), NCBI protein sequences (source recorded), and the COMPARE
  UniProt-mapped set (4 proteins, evidence label `compare_unverified`).
- Non-fungal positives (1,679) are a SEPARATE set from the fungal set. Their negatives are kingdom-matched.
- Mixed rule for plant or animal (26 positives with two AllergenOnline types): plant if any type is plant.

## 3. Results to date (AUROC only; see the status report for full tables)

Fungal set, 200 positives:
- B6 shortcut (mean of 4 seeds): validation 0.577 (sd 0.044), test 0.530 (sd 0.042). Below the stop threshold.
- Learning curve, length-matched negatives only: full-set test 0.860 (95% CI 0.735 to 0.954; 48 clusters, 62 proteins).
  Validation 0.814 (0.679 to 0.932; 41 clusters, 62 proteins).
- Learning curve with family-matched negatives (shared clustering): full-set test 0.675 (0.511 to 0.833; 47 clusters,
  82 proteins, 14 positives). Validation 0.675 (0.583 to 0.884; 48 clusters, 158 proteins). Honest figure for the fungal set.

All-kingdom set, 1,679 positives, 3,216 negatives (plant 904 length-matched plus 654 family; animal 765 length-matched
plus 893 family). 1,544 clusters in the shared clustering (4,895 sequences).
- Learning curve, full training set: test 0.886 (95% CI 0.788 to 0.947; 220 clusters, 694 proteins). Validation 0.805
  (0.743 to 0.863; 254 clusters, 837 proteins).
- Within-kingdom AUROC (each kingdom scored against its own negatives): validation plant 0.805 (0.733 to 0.860; 133
  clusters, 417 proteins), validation animal 0.810 (0.725 to 0.889; 140 clusters, 420 proteins), test plant 0.859 (0.754
  to 0.940; 113 clusters, 321 proteins), test animal 0.906 (0.769 to 0.964; 117 clusters, 373 proteins).
  Table: `analysis/allergen_model/allkingdom/kingdom_auc.tsv`.

Interpretation: the all-kingdom result is the best-supported claim, and it holds within each kingdom. The fungal set is too
small to support a claim. Negatives are weak (may contain unannotated allergens). This is discrimination, not specificity.

## 4. Where the files are

Fungal pipeline (`analysis/allergen_model/`):
- `positives/positives.tsv`, `positives.faa` (200 rows). Built by `build/02_positives.py`, then `build/04_merge_allergenonline.py`.
  Backups: `*.before_allergenonline`.
- `candidates/candidates.tsv`, `candidates.faa` (6,615 candidates). `build/01_candidates.py` with `--exclude` and `--remove`
  (`candidates/remove_near_identical.txt`).
- `full_cluster/` (shared clustering run of the 200-positive set; ok).
- `splits/splits.tsv` (length-matched split, seed 20261010), `splits_seed1..3/` (other seeds), `splits/family_splits.tsv`
  (family negatives assigned by cluster, first run; superseded by the unified split below for training).
- `unified/` (fungal shared clustering: `clu_cluster.tsv`, `unified_splits.tsv`, `lc_B/`). Build: `build/10_unified_split.py`.
- `embed_split/` (`embeddings.npz`, 916 fungal sequences, ESM-2 650M, layer 33; `embed.done`).
- `b6/`, `b6_seed1..3/` (B6 scores and summaries). `pos_vs_pos/` (positive-to-positive identities).
- `lc_A/` (learning curve, length-matched only).
- `family_matched/` (fungal family-matched negatives: `family_negatives.faa/tsv`, `counts.md`, `README.md`, `cluster/`).
- `tests/species_check.py` and `tests/species_out/` (species tables; fungal set).
- `fixtures/` (T1.6 fixture, dry run, `check_fixture.py`: 18 checks pass).
- `pilot/` (pilot projection), `q4_iedb/` (Q4 IEDB counts), `compare/` (COMPARE mapping, `compare_fungi.tsv`,
  `compare_new.faa`, `compare_new.tsv`, `unmapped.tsv`), `allergenonline/` (fungal mapping, NCBI-resolved 26 in
  `ncbi_unmapped.faa`, `blast_best_names.tsv`).
- `allergenonline_all/` (all-kingdom mapping, 1,679 new proteins; `allergenonline_new.faa/tsv`, `allergenonline_fungi.tsv`).

All-kingdom set (`analysis/allergen_model/allkingdom/`):
- `kingdom_negatives.faa/tsv`, `counts.md` (1,669 length-matched; animal 10 short of the 1:1 target in six bins). Build:
  `build/11_allkingdom_negatives.py`.
- `family/` (family-matched non-allergen negatives: `family_negatives.faa/tsv`, `counts.md`, `README.md`, `cluster/`,
  `raw/`, `scripts/`). Built by an agent; see its README.
- `unified/` (`all_input.faa` 4,895 sequences; `clu_cluster.tsv` (1,544 clusters); `allkingdom_splits.tsv`).
- `embed/` (`embeddings.npz`, `embed.done`; 4,895 sequences; 10 truncated at 1,024 residues).
- `lc/` (learning curve, all-kingdom).
- `kingdom_auc.tsv` (within-kingdom AUROC; built by `build/12_kingdom_auc.py`).

Reports and specs:
- `docs/reports/2026-10-10-allergen-pilot-and-species-check.md`
- `docs/reports/2026-10-10-allergen-predictor-status.md` (status report; read this first)
- `docs/reports/2026-10-10-allergen-predictor-literature.md`
- `docs/superpowers/specs/2026-10-09-allergen-predictor-spec.md` (amended; pairing sections are superseded)
- `docs/superpowers/plans/2026-10-10-allergen-predictor-plan.md` (v3; pairing tasks marked superseded)
- `docs/agent-tasks/05-allergen-controls.md` (dated owner exception for pair-in-cluster; now moot)

Input sources (do NOT commit):
- `to_import/allergen/` holds the PDFs (DeepPlantAllergy, AllerTrans, Applm, fungal review), `AllergenOnline.csv`, and `compare.tsv`.
  Per COMMON-RULES, publisher texts must move to `_workdir/curation_texts_archive/` (hashes only) before any commit.

## 5. Job log (HPCC, 2026-10-10)

- Pilot (T0.4): 29797481, completed. Earlier failures 29797130, 29797276 (AVX2 SIGILL on node c03).
- Full fungal cluster and pairs: 29798025 (cluster), 29798026 (pairs); later 29798806, 29798807; 200-positive run 29799590,
  29799591 (all completed).
- Positive-to-positive (B6 input): 29799979 (completed).
- Fungal embedding: 29807328 (completed). Earlier pending attempt 29806863 cancelled.
- All-kingdom embedding: 29821859 (completed; it waited on the group CPU limit before starting).

## 6. Pitfalls already hit (avoid repeating)

- Default MMseqs2 module is the AVX2 build. It crashes with SIGILL on node c03. Use
  `/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs` (non-AVX2). All job scripts use it.
- `sort | head` under `pipefail` gives exit 141 (SIGPIPE). Write to a file first.
- `/scratch/jstajich` is not writable from the login node. Use a directory under the repo for local debug runs.
- UniProt query syntax rejects range form for annotation score (`[4 TO 5]`). Use `annotation_score:4 OR annotation_score:5`.
- UniProt ID mapping: merged accessions redirect; inactive entries have no sequence. Handle both. NCBI efetch resolves what
  UniProt cannot (for example 2026 *Asp f* entries and 26 AllergenOnline GenBank IDs). Length check against source confirmed match.
- A candidate build that does not exclude positives leaks (P12547 case). `--exclude` is required.
- `conda` env `adhesionPred` has torch 2.10 with CUDA and fair-esm. `/usr/bin/python3.12` has no torch.
- `exfab` jobs can wait on `AssocGrpCpuLimit` for a long time.
- `get_esm_embeddings` default `repr_layer` is 6. Pass 33 explicitly (the embed script does).
- The COMPARE file has no header and 25 placeholder IDs with no sequence. Column order: species, type, name, allergen, accession, length, year, blank.
- The AllergenOnline `unmapped.tsv` rows were rescued through NCBI; the length match is the check, not the accession.

## 7. Open decisions and actions (owner)

1. COMPARE full download: the form emails files to an address you provide. Put them in `to_import/allergen/` and tell me the
   filename. Four mapped COMPARE proteins already in positives (`compare_unverified`). Thirteen COMPARE records for 2026 *Asp f*
   entries came from NCBI and are NOT yet added. Decide whether to add them.
2. Whether to add the 13 NCBI-sourced COMPARE records (2026 *Asp f* and older entries).
3. Whether to reopen pairing if more positives change the pairable count.
4. Whether the 500-negative floor stays a planning value or becomes a rule.
5. Whether to commit the work and create the `allergen-predictor` branch. Publisher PDFs must be archived first.

## 8. Next steps (plan items still open)

1. M1b random forest and M1L length-matched run (all-kingdom set first).
2. Five seeds and one label-permutation run for each reported model.
3. Threshold rule (O2, validation sensitivity 0.80) and evidence field (T5). Module names are hard-coded in
   `src/cellsurface_sorting_hat/modules/cli.py`; a new module needs code, not config.
4. Per-group report for the all-kingdom model that matches `kingdom_auc.tsv`.
5. Owner review before any merge or release. Nothing is merged.

## 9. Known limits (state them in any report)

- Negatives are weak. Some may be unannotated allergens. Evidence-grade negatives do not exist yet.
- Fungal test set has 14 positives. Intervals are very wide.
- Cluster counts for the family-matched fungal learning curve come from a mixed cluster space (family run clusters plus
  shared clusters). The shared-clustering run (`unified/`) is the consistent one.
- The all-kingdom 1:1 negatives are 10 short in six animal length bins.
- B6 is near chance across seeds; that is the expected result for length-matched negatives, not evidence of a clean task.
- Some near-identical positives exist (NCBI and Swiss-Prot duplicates). Clustering keeps them in one split.
