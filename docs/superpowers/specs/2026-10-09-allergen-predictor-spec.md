# Allergen predictor (evidence only): spec, draft 5

Date: 2026-10-09 (draft 5 edited 2026-10-10). Status: draft 5. Not yet approved. No code written.
Supersedes drafts 1-4. Draft 5 applies the fixes from the final review of draft 4. Owner decisions A-E
(section 10) are recorded. Two points are still open owner decisions (section 10, F and G).

Owner decision (2026-10-09): lift the "no allergenicity model" rule for an evidence-only model.
Weak Swiss-Prot negatives may train it. Evaluation negatives still need an independent reason.
The rule came from the agent directive of 2026-10-06 (`docs/agent-tasks/05-allergen-controls.md`),
not from an owner statement.

## 1. Purpose and scope

Give `cellsurface_sorting_hat` one more allergen evidence line. The existing line is BLAST homology to
the WHO/IUIS fungal set (`allergen_homology`). The new model scores a protein for allergen-like sequence
features without an alignment, so it can flag allergens with no close fungal homolog.

- Evidence only. No mechanism category. No change to any existing call.
- Status `unvalidated` until the owner accepts a measured result.
- Research use only. No clinical, regulatory or product claim.
- No threshold is set before measurement. A threshold is chosen on the validation split only, never on
  the test split, and is a config item.

## 2. Literature

Full text read unless stated. Summaries in my own words. Figures as reported in each paper.

| Paper | Data | Split and similarity control | Reported result | Use in this spec |
|---|---|---|---|---|
| DeepPlantAllergy, Briefings Bioinf. 2025 (bbaf605) | 955 plant allergens (AllerBase), MMseqs2 at 85%, coverage 0.8, to 614 clusters. Negatives from UniProt reviewed entries, filtered; 1,817 candidates, 798 after clustering at 20%, 769 after removing any within 20% identity of a positive. About 20% of final negatives are non-plant. | Positives and negatives clustered at 20%; each cluster in one split; 70/20/10 (431/122/61 proteins); five seeds on a fixed split | Test F1 93.9% (ESM-1B) and 93.6% (ProtBert); MCC 0.879 and 0.873; ROC-AUC 97.7-97.8%. One-hot (narrative): F1 80.8%, MCC 0.617, AUC 86.2%. Table values not extracted; check the table before quoting a second decimal. | Design template for negatives and splits. Code: github.com/Lilly-dh/DeepPlantAllergy (licence not checked). |
| AllergenAI, BMC Bioinf. 2025 (s12859-025-06302-1; PMC12625376) | 9,646 allergens (SDAP 2.0, AlgPred2). Negatives: cupin family (PF00190), 27,182 sequences, clustered at 95%, then removed any with an SDAP hit (E <= 0.1); 2,775 kept. | Clusters assigned to train or test; five folds by cluster within training; similarity threshold raised from 40% to 90% | Results section: test accuracy 0.94, AUROC 0.97. A reviewer reports the Methods give test accuracy 0.81. Not resolved: the full text is not local. Training-pool counts do not reconcile between Methods and Results. | Shows the risk of a negative set from one family. |
| AllerTrans, Biol. Methods Protoc. 2025 (bpaf040) | AlgPred 2.0 training and validation sets | 5-fold CV. The methods text describes no identity-based split. The "independent" AlgPred 2.0 result is on the validation set, from the same dataset family | 5-fold: sensitivity 97.91%, specificity 97.69%, MCC 0.956. AlgPred 2.0 validation: sensitivity 93.5%, specificity 95.38%, MCC 0.900. Recombinant check (42 allergen + 42 non-allergen UniProt entries, Table A1): sensitivity 95.24%, specificity 80.95%. | Random CV is optimistic. The recombinant specificity drop is the only out-of-distribution result in the paper. Code: github.com/faezesarlakifar/AllerTrans (licence not checked). |
| Applm, arXiv 2508.10541v1 (2025) | Curated allergen and non-allergen sets | Pipeline with an inter-split threshold T_S and an inter-class threshold T_C. Positives are clustered and split first; each negative is paired with a positive at identity >= T_C and assigned the paired positive's split; then violations of T_S are removed. The paper states that CD-HIT and MMseqs2 cannot by themselves guarantee the similarity separation needed. | Encoding comparison: average AUROC 0.844-0.864 with frozen pLM embeddings, 0.742 one-hot, 0.684 BLOSUM62. External-set average AUROC 0.872; average AUPRC 0.700 against 0.202 for random. AUROC across T_S/T_C scenarios ranges from 0.722 to 0.951 (0.942 is the average in the easiest scenario). Classifier: random forest on frozen embeddings. | Best evaluation design in the set. Length Control (matching length distributions between classes in training) is reported as helping when the length distributions of the two classes are similar in the test set (KS statistic < 0.25). Code: github.com/brianwongsh/Applm (licence not stated in the text). |
| Fukutomi and Taniguchi, invited review, 2015 (`64_5.pdf`) | Not applicable | Not applicable | Not applicable | Background on fungal sensitisation and cross-reactivity. No machine learning. |
| Allerdictor, Bioinformatics 2014 (30(8):1120) | Abstract-level summary only; **no PDF in `to_import/allergen/`** | Three benchmarks A, B, C of varying similarity | AUPRC 0.97 / 0.91 / 0.85 as reported; not verified | 6-mer sparse counts with a linear SVM. Reused as baseline B3 (unverified design). Web server did not resolve; no licence found. |

Lessons carried into this spec:
1. Split by sequence cluster, then check the split by measurement. Random CV is optimistic (AllerTrans).
2. A negative set from one family (AllergenAI) or one kingdom (DeepPlantAllergy) can teach the model that
   family or kingdom. Report the taxonomy of the negatives.
3. Negatives must be matched to positives in similarity, or the model can learn distance to positives
   instead of allergenicity. Applm pairs negatives with positives at a similarity threshold (owner decision).
4. Simple encodings lag language-model embeddings (one-hot F1 80.8% vs about 93.9%, DeepPlantAllergy
   narrative). A composition baseline is still required. The 2026-09-27 model review in `AGENTS.md` found
   composition alone reached ROC-AUC 0.995 under leave-family-out CV on the shipped adhesin model, not on
   allergens. It is a reason to test composition here, not evidence about allergens.
5. Length can be a shortcut. Applm reports that allergens in databases tend to be shorter. A length-only
   baseline is required, plus a length-matched run (M1L).
6. Applm's classifier is a random forest; its alternative is a feed-forward network. Our M1 is logistic
   regression and M1b is a random forest.
7. Random-CV and "independent" validation numbers in the literature are not comparable with each other or
   with our data.

## 3. Data

### 3.1 Positives (evidence-grade)

- WHO/IUIS fungal allergen set: 104 unique UniProt accessions in
  `analysis/allergen_scoping/iuis_fungal_isoallergens.tsv` (column `AccUniProt`), from 121 isoallergen rows.
  The earlier figure of 111 "usable" entries counts isoallergen rows, not accessions. The build step
  re-derives the count.
- UniProt KW-0020 (Allergen), Fungi, reviewed: 111 entries in `analysis/allergen_scoping/uniprot_kw0020_fungi.tsv`.
  15 names contain "homolog". Those are flagged `evidence_level=similarity_only`.
- The two fungal sets share 63 accessions. **Union by accession, one row per protein** (152 proteins
  maximum), with source flags (`iuis`, `kw0020`) and the evidence level. A protein is counted once in any split.
- Non-fungal positives (AllerBase, SDAP 2.0, AlgPred2): included in training (owner decision Q1). **Not yet
  counted**; no scoping table exists for them. The build step writes one with the same columns as the fungal
  tables. Whether the sources can be obtained is not yet checked.
- Every result is reported for all positives and separately for fungal positives.
- Publisher texts (the PDFs in `to_import/allergen/`) must not be committed. COMMON-RULES requires them in
  `_workdir/curation_texts_archive/`, with only hashes in the repo. Move them before any commit.

### 3.2 Negatives for training (weak)

- Source: reviewed Swiss-Prot entries with protein-level evidence, no allergen keyword (KW-0020), no allergen
  comment, and no "antigen" word, length 10-850 aa (DeepPlantAllergy's range). Also exclude human (taxon 9606)
  and viruses (10239), require annotation score above 3, and exclude "Cancer" (DeepPlantAllergy's other filters).
- **Negative selection (owner decision, 2026-10-10: pairing dropped).** Weak negatives are sampled from the
  filtered candidate pool so that their length distribution matches the positives' length distribution
  (length-matched sampling, the Applm Length Control idea). No similarity condition is applied to the choice of
  negatives. Sampling is done with a fixed seed, in length bins, with the bin counts written to `counts.md`.
- **Similarity is measured, not controlled.** Pairing is no longer part of the design. Similarity between
  negatives and positives is measured by B6 (nearest neighbour), which is a required shortcut check
  (section 5). Its result decides whether the model can be read at all (stop margin, section 5).
- **Build gate (planning value 500):** the build counts the length-matched negatives available per split. A
  split with fewer than **500 negatives** is reported as not feasible for that split. 500 is a planning value,
  revised once the count is measured.
- **Label-noise check:** count the candidate negatives that match a WHO/IUIS or AllergenOnline positive at
  identity >= 0.6 (coverage >= 0.5). Such hits are possible unannotated allergens. They are reported as a
  limit, not removed silently. Candidates at identity >= 0.9 are removed from the pool (owner decision, 2026-10-10).
- Kingdom mix: not controlled. The resulting kingdom mix of the sampled negatives is reported as a table.
- Label: `weak_negative`, with the query and date recorded. Not evidence-grade.
- Pair-in-cluster exception (owner decision F) and the T_C sweep (0.6, 0.7) no longer apply. Pairing is
  dropped. The O7 coverage rule is kept only for the label-noise check.
- The count of these negatives is **not yet measured**. No number is assumed here.

### 3.3 Negatives for evaluation (evidence-grade)

- Must have an independent reason under `COMMON-RULES.md`: IgE or serum test negatives, or published binding
  assays with a negative result.
- First step: count IEDB IgE-negative records for fungal proteins, by species and cluster, and report the count
  before building anything else (owner decision Q4). **Done 2026-10-10** (`analysis/allergen_model/q4_iedb/counts_q4.md`):
  328 IgE-negative fungal records, 16 proteins, 9 species. Six proteins have negative records only; five of
  them are named allergens (Sta c 3, Mala f 3, Tri r 2, Asp f 5, Mala f 4). Ten proteins have both negative and
  positive records, all allergens. **No usable protein-level negatives.** The 20-cluster floor cannot be met.
  The model therefore has no specificity claim from IEDB.
- If there are fewer than 20 independent **evidence-grade** negative clusters, the model has **no specificity
  claim**. The report then gives sensitivity on positives and calls on weak negatives, labelled optimistic.
  Weak negatives never count toward the 20-cluster floor.

## 4. Splits, leakage and evaluation sets

Order of operations (follows Applm's pipeline; corrects draft 4, which was circular):

1. Cluster the **positives** with MMseqs2 at 30% identity, coverage 0.5 (COMMON-RULES rule 4). Installed on
   HPCC as `module load MMseqs2/17-b804f` on the login node (`docs/HANDOFF-2026-10-06.md`). Not available on
   this machine.
2. Split the positive clusters: 70% train, 15% validation, 15% test, fixed seed. Fractions are fixed in this
   draft. The fungal floor (section 6) decides whether a fungal table is reported.
3. Remove positive violations of T_S: any positive pair that crosses a split at identity >= T_S is resolved
   by removing one protein. **T_S = 0.4** (owner decision G(a), 2026-10-10), measured with MMseqs2 all-vs-all.
4. Length-matched sampling of weak negatives (section 3.2). Each sampled negative is assigned to the split of
   its **sequence cluster** (MMseqs2 clusters of positives and candidates together, section 4 step 1). A cluster
   is never split across train and test. Candidates in a cluster with no positive are assigned by cluster alone.
5. Remove negative violations of T_S across splits, in both classes (Applm's order).
6. Measure, do not assume, the inter-split identity: all-vs-all alignment across splits (aligner and parameters
   recorded), and report the measured maximum. A split is valid only if the measured maximum is at or below T_S.

Evaluation sets:
   a. Primary test set: the positives and weak negatives **of the test split only**. Training positives are
      never scored in a test table. Any table that uses weak negatives carries the word "optimistic".
   b. Leave-family-out: hold out one Pfam allergen family at a time. Families are small (PF16541: about 3 fungal
      proteins), so these are descriptive only. Per-family retraining needs compute not budgeted in section 8.
   c. Close-homolog set: dropped (owner decision E). The report states that near-homolog performance is not
      measured here.
   d. In-reference and out-of-reference: proteins in the WHO/IUIS reference (`in_reference`) are reported with
      and without.
   e. The hard-negative stratum is deleted (owner decision C). No unpaired candidate enters any table.

- **Leave-species-out for B2 and B6 only** (owner decision): when a positive is scored, the reference set
  excludes its own species. The model is not retrained per species. Its species generalisation is not directly
  tested, and the report states this.
- **Prevalence (realistic scan):** score every protein in one or more Fungi_5k proteomes. Report the score
  distribution and the number of **calls** per proteome across a range of thresholds, not only the validation
  operating point. A threshold chosen on paired data does not transfer to an unpaired proteome. Calls are not
  false positives: without labels only calls are known. Precision at a stated prevalence appears only as a
  conditional sensitivity analysis over an assumed specificity range, labelled as assumed. If section 3.3 gives
  fewer than 20 evidence-grade clusters, this is the only precision-type output, and the report says so.

## 5. Models and baselines

Each baseline is required. No baseline is dropped without a written reason.

| ID | Model | Input | Role |
|---|---|---|---|
| B1 | Composition | amino-acid and dipeptide composition, length | floor |
| B2 | BLAST to IUIS | existing `allergen_homology` output, leave-species-out | sanity floor (see note) |
| B3 | 6-mer linear SVM | 6-mer counts (Allerdictor's design, **unverified**: no PDF read) | cheap classical baseline |
| B4 | Length only | protein length, logistic regression | checks the length shortcut |
| B5 | One-hot + linear head | one-hot per residue, pooled | encoding floor (DeepPlantAllergy, Applm) |
| B6 | Nearest neighbour to training positives | best identity to any training positive (MMseqs2, training split only) | sanity floor (see note) |
| X1 | One published external predictor (owner decision) | **Dropped (owner decision O5, 2026-10-10).** Reason: AllerCatPro 2.0 caps batch submissions at 50 sequences and states no terms of use; the AlgPred 2.0 host did not resolve; Allerdictor's server is not reachable. The report states this reason. | not run |
| M1 | ESM-2 650M + logistic regression | mean-pooled embedding, padding excluded | main candidate |
| M1b | ESM-2 650M + random forest | same embedding | classifier comparison (Applm) |
| M1L | M1 with length-matched negatives | as M1, negatives filtered to the positives' length range | length sensitivity run |

Notes:
- **B6 is a required shortcut check (updated 2026-10-10, pairing dropped).** Negatives are no longer chosen by
  similarity, so similarity to positives may separate the classes. B6 (nearest-neighbour identity to training
  positives) therefore measures whether similarity alone explains the result. The stop margin of owner decision
  O1 applies, as amended 2026-10-10: the stop uses the mean validation B6 AUROC over four split seeds (20261010, 1, 2, 3),
  with the seed standard deviation reported. A mean above 0.60 halts the arm for owner review before any embedding model
  is trained. Single-seed results are reported but do not halt. The
  report records the B6 value whether or not it halts. B2 is kept as a reference to the existing module.
- Use ESM-2 650M (`esm2_t33_650M_UR50D`, the variant AllerTrans names). Record the version.
- Mean pooling excludes padding. Padding in the mean makes scores batch-dependent (`AGENTS.md`).
- Hyperparameters are chosen on the validation split only.
- **Disclosure for X1:** AllerCatPro 2.0 and AlgPred 2.0 were built on allergen databases that contain the
  WHO/IUIS positives. Scoring X1 on our positives is therefore scoring it partly on its own training data.
  X1 is a reference point, not a fair comparison, and the report says so.
- **Controls on every reported model:** five seeds (report mean and spread), and one label-permutation run.
  The permutation shuffles labels **within each fixed split**, so each split keeps its composition and its
  pairs stay intact. Metrics from a permutation run should fall to chance; a permutation run that scores well
  means the pipeline leaks.
- Motif output (Integrated Gradients, optional) is a hypothesis only. The 2026 residue-attribution preprint
  (arXiv 2606.22181) was title-checked only; its finding is not used.

## 6. Metrics and gates

- Comparable pair: AUROC and sensitivity at fixed specificity, with 95% cluster-bootstrap intervals (1,000
  replicates, resampling clusters). MCC is reported as well.
- AUPRC depends on class ratio. It is reported only beside the positive rate of its table.
- **Specificity claims use evidence-grade negatives only.** A sensitivity-at-specificity row is a specificity
  claim. It is reported only when the test split has at least 20 independent evidence-grade negative clusters
  (owner decision D). Weak negatives do not count toward that floor. Below it, the row is marked "not
  reportable", and AUROC and calls are reported instead.
- **Resolution (owner decision G(b), 2026-10-10):** the sensitivity-at-specificity 0.90 row needs 20
  evidence-grade negative clusters (decision D). The 0.95 row needs **40**, because at 20 clusters one step is
  0.05 and 0.95 cannot be resolved finely. Each row shows its cluster count.
- **Bootstrap degeneracy:** a replicate with a single class has no AUROC. Such replicates are skipped and
  their number is reported. A row with more than 5% skipped replicates is marked "not reportable".
- **Interval rule:** test positives are few (the reviewer estimates 10-20 fungal positives in a 15% test split).
  Intervals are reported at their true width. The 0.10 half-width of COMMON-RULES rule 7 is not met and is stated
  as a limit. Adding positives is an owner action (AllergenOnline, COMPARE).
- **Fungal floor (owner decision):** a fungal result is reportable only with at least 20 independent test
  clusters. Below that, the table is marked "not reportable". The cluster count is measured on HPCC before the
  report is written.
- **Label-noise limit:** AUROC and MCC here measure discrimination between allergens and negatives paired by
  similarity, under unquantified systematic noise from unannotated allergens among the negatives. The label-noise
  check in section 3.2 bounds that noise only partly. The report says so.
- Report every table for: all positives, fungal positives, leave-family-out set, in-reference and
  out-of-reference.
- Every baseline appears in the same table.
- Gate values are set after measurement, by the owner. This spec sets none.
- Published numbers are not targets. Their splits differ (section 2).

## 7. Output and integration

- Score table `allergen_model.tsv`: protein ID, score, model ID, seed, T_C, split, training-data hash.
- An evidence field `allergen_model.score` is evidence, not a call, and stays `unvalidated`. **A new evidence
  field needs code.** Module names are hard-coded in the CLI (`src/cellsurface_sorting_hat/modules/cli.py`,
  `ModuleSpec("allergen_homology", ...)`) and the Pfam module lists its families in `pfam.py`. The change is a
  new module and CLI wiring, planned as a code task, not a config edit.
- No change to `iuis_allergen_similarity`, `iuis_allergen_homolog` or any category.
- Control-set outputs under `data/controls/allergen_model/`, following COMMON-RULES:
  - `sequences.faa`: all sequences used.
  - `controls.tsv`: evidence-grade rows only, with the 14 required columns and a source and quote on every row
    (rule 1). Positives and evaluation negatives with an independent reason.
  - `training_negatives.tsv`: weak Swiss-Prot negatives with query and date. Not evidence-grade. Checks run on
    `controls.tsv` only. This is the owner-approved deviation from rule 1; the README records it.
  - `clusters.tsv`: MMseqs2 cluster assignments and parameters.
  - `rejected.tsv`: candidates removed, with reason and count per T_C (pairing failure, split violation, the
    pair-in-cluster conflict if the exception is not granted).
  - `unverified_notes.md`: any claim from memory, marked "from memory, not checked".
  - `README.md`: method, versions, the schema deviation, and the banner below.
  - `counts.md`: the counts reported in the text.
  - Taxon IDs are checked against `names.dmp` (rule 5), at `/srv/projects/db/taxonomy`.
  - Large tables are compressed (CLAUDE.md): `.tsv.gz` for files shared outside the cluster, `.zst` for pipeline
    intermediates.
- Report with the banner: "Research use only. Not a regulatory allergenicity assessment and not a diagnostic
  result."

## 8. Compute and provenance

- Build, embedding and training runs are SLURM jobs on UCR HPCC, each sized to about 1 to 1.5 h (`CLAUDE.md`).
  Use `${SCRATCH:?}` for intermediates (fails loudly if unset), trap cleanup on exit. Copy results to `/bigdata`
  under a temporary name, then `mv` into place. Never write to `/tmp`. Do not use `BASH_SOURCE` for paths.
- Repo Python scripts run with `/usr/bin/python3.12` (CLAUDE.md), not the default `python3`.
- The GPU queue `exfab` was approved for embedding jobs in the Step 1 work. Confirm before use here.
- Embeddings are computed once per protein set and cached.
- Record: data hashes, MMseqs2 version and parameters, aligner and parameters for the T_S check, ESM-2 variant,
  seeds, T_C, code commit.
- Code is written from the method descriptions. DeepPlantAllergy, AllerTrans and Applm code is **not copied**
  until its licence is checked.
- Leave-family-out retraining (section 4b) and five-seed runs multiply the training jobs. The job count is set in
  the plan, once the cluster and pairing counts are known.
- The ESM-2 650M embedding jobs were gated on the IEDB negative count (section 3.3). That count is now known
  (2026-10-10). The gate is cleared, but the result means no evidence-grade specificity claim is possible.
- Git: push over HTTPS and open a pull request, per the repo workflow.

## 9. Decisions recorded (owner, 2026-10-09 and 2026-10-10)

- Q1 Positives: all kingdoms, fungal results reported separately.
- Q2 Weak negatives: Swiss-Prot, protein-level evidence, no allergen keyword or comment, no "antigen" word.
- Pairing (2026-10-10): dropped. Pairing gave too few candidates (74 at T_C 0.6, 39 at 0.7, against a 500
  floor), so length-matched sampling replaces it. Revisit only if more positives change the pairable count.
  The T_C sweep (0.6, 0.7) no longer applies.
- A (2026-10-10): build gate; arm dropped below 500 paired negatives (planning value).
- B (2026-10-10): pairing wins over kingdom mix; the resulting mix is reported.
- C (2026-10-10): hard-negative stratum deleted; unpaired candidates dropped and logged.
- D (2026-10-10): sensitivity-at-specificity rows need at least 20 evidence-grade negative clusters.
- E (2026-10-10): close-homolog set dropped; report states near-homolog performance is not measured.
- Q3 Split clustering: MMseqs2 30% identity, coverage 0.5 (for positives; see section 4 order).
- Q4 Evaluation negatives: count IEDB IgE-negative records first.
- Test power: report wide intervals honestly; add positives as an owner action.
- Prevalence: calls per proteome across thresholds; precision only as a labelled conditional analysis.
- B2 and B6: leave-species-out, and treated as sanity floors (draft 5 note; the reviewer's point 4).
- Length: length filter on negatives, plus the length-matched run M1L.
- Controls: five seeds and a label-permutation run within fixed splits.
- Positives: union by accession, one row per protein, with source flags.
- Classifiers: M1 logistic regression, M1b random forest.
- Data format: `controls.tsv` for evidence rows; `training_negatives.tsv` for weak negatives.
- External reference: one published predictor (X1), with the contamination disclosure.
- Leave-species-out is limited to B2 and B6. The model is not retrained per species.
- Fungal floor: 20 independent test clusters before a fungal result is reportable.

## 10. Owner decisions

Decided:
- A. Pairing feasibility: count first; drop an arm below 500 paired negatives (planning value).
- B. Kingdom mix versus pairing: pairing wins; the resulting mix is reported.
- C. Hard-negative stratum: deleted. Unpaired candidates dropped and logged.
- D. Specificity floor: 20 evidence-grade negative clusters.
- E. Close-homolog set: dropped.

Open (not yet decided):
F. **Pair-in-cluster exception: decided (owner, 2026-10-10).** Granted for this model, recorded in task 05.
- **G(a). T_S: decided (owner, 2026-10-10).** T_S = 0.4, measured with MMseqs2 all-vs-all.
- **G(b). Minimum for a 0.95 specificity row: decided (owner, 2026-10-10).** 40 evidence-grade negative
  clusters. The 0.90 row keeps the 20-cluster floor.

## 11. Still open (not design decisions)

- Whether to add a tree-based model beyond M1b, and an M2 MLP (only if M1 is weak on validation).
- Licence checks for DeepPlantAllergy, AllerTrans and Applm code.
- Owner action to obtain more evidence-grade positives (AllergenOnline; COMPARE needs an email form).
- Allerdictor: no PDF in the repo; its row stays unverified until one is provided.
- AllergenAI: test accuracy (0.94 in Results, 0.81 in Methods, per the review) to be checked against full text.
- DeepPlantAllergy: table values to be read from the PDF table, not the narrative.
- Whether non-fungal positive sources (AllerBase, SDAP 2.0, AlgPred2) can be obtained and under what licence.
- Publisher PDFs to move to `_workdir/curation_texts_archive/` before any commit.
- The length filter, pairing and the length-matched run may not share a feasible negative set. The build gate (A)
  checks this.

## 12. Not in this spec

- Any allergenicity call for a clinical or product decision.
- IgE cross-reactivity claims.
- Vendoring third-party code or model weights before a licence check.
- Changes to the thresholds of existing modules.

## 13. Next steps

1. Owner decides F and G (section 10).
2. Write the plan, with the build gate (A) and the order of operations (section 4) as the first tasks.
3. Independent review of the plan.
4. Build on a feature branch. No merges without the owner's go-ahead.
