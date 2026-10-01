# Design spec: step 1 surface glycoprotein model (ML candidate) and its comparison with the rule

*2026-09-30. Status: draft for owner review. No code exists for this spec. Base: `main` 610918a.*

Inputs: `docs/PLAN-2026-09-30-pipeline-and-decisions.md` (branch `readme-file-index`; it has no
"Changes made later" section, so the later decisions come from the owner's brief: old pickles deleted,
no legacy mode, two step 1 candidates); the review (cited "review n"); `docs/TOOL-ARCHITECTURE.md`;
`docs/plans/2026-09-30-design-review-fable.md` (branch `issue-9-plan`, cited "Fable #n").

Counts marked "GAF" come from a throwaway script run on 2026-09-30 on `sgd.gaf.gz` and `cgd.gaf.gz`
(GO release files dated 2026-05-28) and `go-basic.obo` (2026-08-08), with `is_a` and `part_of`
ancestors and `NOT` rows skipped. Deliverable D1 must reproduce them.

## 1. Purpose and non-goals

**Purpose.** Step 1 asks: does this protein go through the secretory pathway and end up outside the
plasma membrane (cell wall, outer face of the plasma membrane, or secreted)? This spec defines the ML
candidate and one test bench that compares it with the rule (SignalP + GPI anchor + Ser/Thr content).
The owner picks rule, ML or hybrid from that comparison.

**Non-goals.** The model must not be described as:

- an adhesin predictor (old model: about 12% adhesin precision on S288C, review 4.3);
- a glycosylation predictor (GO does not record glycosylation; see section 2);
- a moonlighting-protein detector (Hsp60 and enolase reach the surface without a signal peptide; plan 2d);
- an antigen or epitope predictor;
- validated outside the clades in the truth set (section 3).

The score name is `surface_glycoprotein_score`. "Glycoprotein" stays in the name by owner decision.
The label in this spec is a localisation label. This mismatch is open question Q1.

## 2. Label definition

### 2.1 Facts checked

| Fact | Value | Source |
|---|---|---|
| `surface.tsv` rows | 3,573 (3,540 `uniprot_surface_kw`, 33 `curated_literature`) | `data/curated/surface/surface.tsv` |
| Negatives in `surface.tsv` | 0. Every row has a surface keyword | same; `data/curated/surface/README.md` |
| `signal_peptide` yes / no | 3,338 / 235 | same |
| `gpi_anchor` yes / no | 325 / 3,248 | same |
| Rows longer than 1,022 aa | 164 | same, `length` column |
| `reviewed` (Swiss-Prot) | 982 | same |
| Eurotiomycetes literature rows | 21: 10 Onygenales, 11 Eurotiales. YPS3 has no accession. HSP60 is `moonlighting` | `data/curated/adhesins/eurotiomycetes_seeds.tsv` |
| Old training data | 581 positive (FLO/ALS), 1,500 random negative | `data/positive`, `data/negative` (counted with `zcat | grep -c '>'`) |
| GPI anchoring in GO | no current cellular-component term. GO:0031225, GO:0046658 and GO:0031362 are obsolete | `go-basic.obo` |

The old FLO/ALS-vs-random data are not reused. `surface.tsv` cannot supply negatives. Its
labels come from UniProt keywords, and those keywords share inputs (signal-peptide and GPI predictors)
with the rule. So `surface.tsv` must not be the test truth for a rule-versus-ML comparison.

### 2.2 Candidate positive definitions

| ID | Positive if | S288C (GAF) | *C. albicans* (GAF) |
|---|---|---|---|
| P-wall | GO cell wall (GO:0005618, includes fungal-type cell wall) | 107 | 124 |
| P-ext | P-wall, or GO extracellular region (GO:0005576) | 128 | 261 |
| P-gpi | P-ext, plus plasma-membrane proteins with curated GPI evidence (UniProt reviewed entry or literature) | not counted | not counted |

All three use non-IEA evidence and exclude any gene that also has a non-IEA cytosol, nucleus or
mitochondrion annotation. In this ontology traversal, cell wall terms are descendants of extracellular
region, so P-wall is a subset of P-ext.

The exclusion is needed. With experimental evidence only, 102 of 305 *C. albicans* "extracellular
region" genes and 39 of 125 S288C genes also carry a cytosol, nucleus or mitochondrion term (GAF).
*C. albicans* ENO1, TDH3, PGK1 and ADH1 have IDA annotations to fungal-type cell wall. Under the
exclusion these fall into an **ambiguous** class (48 genes in S288C, 100 in *C. albicans*, non-IEA).
Ambiguous genes are not used for training and are reported as their own stratum.

Evidence tiers for positives: (T-a) curated GO, non-IEA; (T-b) curated literature rows; (T-c) UniProt
keyword only. T-c is used for training at most (Q4). T-c is never test truth.

### 2.3 Candidate negative definitions

| ID | Negative if | S288C (GAF) | *C. albicans* (GAF) |
|---|---|---|---|
| N-int | non-IEA cytosol, nucleus or mitochondrion term, and no endomembrane, plasma membrane, vacuole, membrane, cell periphery, extracellular or cell wall term at any evidence level | 3,358 | 1,795 |
| N-sec | non-IEA endomembrane system, plasma membrane or vacuole term, and no cell wall or extracellular term at any evidence level | 1,674 | 970 |

N-sec holds the hard negatives: ER, Golgi and vacuole residents and membrane proteins. Without N-sec,
a model can score well by learning "has a signal peptide". N-sec also holds GPI-anchored
plasma-membrane proteins (for example the yapsins, assumption), which face outside. P-gpi moves those
to the positive side, but needs a GPI evidence source that is not a predictor (Q2).

Genes with no qualifying annotation are **unlabelled** (3,517 S288C and 4,198 *C. albicans* genes
with some cellular-component term, plus genes with none). They are never negatives. The proteome is
positive-unlabelled (review 6). The labelled truth set is not, because N-int and N-sec need positive
evidence of a non-surface location.

### 2.4 Recommendation

Use **P-ext** for positives and **N-int + N-sec** for negatives, with P-gpi added if Q2 is answered yes.
Reasons:

1. P-wall excludes secreted enzymes such as SUC2 and PHO5. These are secretory, extracellular and
   glycosylated (assumption from the literature, not checked in GO). The rule's Ser/Thr term cannot
   separate them from wall proteins either.
2. GO cannot say whether a protein is glycosylated. A "glycoprotein-only" label cannot be built from GO.
3. N-sec forces both candidates to separate surface proteins from other secretory-pathway proteins.
   That is the hard part of the task.
4. Keep `wall` versus `extracellular-only` as a column. Recall is then reported for each subset.

## 3. Data sources and feasibility

HEAD requests on 2026-09-30 (`current.geneontology.org/annotations/` unless stated):

| File | Result |
|---|---|
| `sgd.gaf.gz` | 200, 2.6 MB, downloaded to `/tmp` |
| `cgd.gaf.gz` | 200, 4.6 MB, downloaded; taxon 237561 used |
| `pombase.gaf.gz` | 200, 2.0 MB, not downloaded |
| `aspgd.gaf.gz` | 403, not available |
| SGD and CGD direct downloads (`downloads.yeastgenome.org`, `www.candidagenome.org`) | 200, not downloaded |

Species with enough manual annotation: S288C and *C. albicans* SC5314 (counts in 2.2 and 2.3).
*S. pombe* has a GAF but is not in the owner's truth set (Q6).

What is missing:

- **No GO truth for Eurotiomycetes.** AspGD returned 403. The only Onygenales/Eurotiales truth is the
  21 literature rows, and about 19 are usable (YPS3 has no accession; HSP60 is not secretory).
  Confidence intervals on about 19 proteins will be wide. A per-clade number there is a smoke test,
  not an estimate.
- **No negatives for Eurotiomycetes.** The literature rows are surface proteins and hard negatives for
  *adhesion*, not negatives for *surface*. Precision and FPR cannot be measured in that clade.
- **Mapping GAF IDs to sequences** (SGD `S000…`, CGD `CAL…`) to protein FASTA is not built.
- **No Basidiomycota, Mucoromycota or chytrid truth** is planned.

## 4. Dataset construction

1. **Extract** labels from the GAFs (deliverable D1). Record the GAF date and file hash in the output.
2. **Attach sequences** from the SGD and CGD protein FASTA. Record the release.
3. **Dedupe.** Hash each cleaned sequence (SHA-256). Keep one copy of exact duplicates within a class
   (for example S288C ASP3-1 to ASP3-4, assumption that they are identical, to be checked). Drop any
   hash present in both classes (plan default 4; Fable #32).
4. **Cluster** all labelled proteins together with MMseqs2: `easy-cluster --min-seq-id 0.3 -c 0.5
   --cov-mode 0` (as in `analysis/model_review/run.sh`).
5. **Class balance.** Train on all labelled proteins with class weights. The truth-set ratio (P-ext
   about 389, negatives about 7,800) is not the genomic prevalence, which is not measured. Report
   precision on the truth set and re-weighted to a PU prior estimated on whole proteomes.
6. **Splits.**
   - S1: homology-grouped 5-fold (`StratifiedGroupKFold` on cluster ID), both species pooled.
   - S2: leave-species-out (train S288C, test *C. albicans*; then the reverse).
   - S3: leave-clade-out (train both yeasts, test the Eurotiomycetes literature rows).
   - Hyperparameters (C for LR, rule thresholds) are tuned in an inner loop only.

| Leakage risk | Control | Test |
|---|---|---|
| Homologs across folds | cluster-grouped folds | `test_no_cluster_spans_train_and_test` |
| Same sequence in both classes | drop hash | `test_no_hash_in_both_classes` |
| IEA labels derived from SignalP/InterPro | IEA excluded from labels | `test_truth_set_has_no_iea` |
| UniProt keywords share inputs with the rule | T-c never in test truth | `test_test_truth_sources` |
| Thresholds tuned on the test fold | nested loop | `test_threshold_fit_uses_train_only` |
| S3 homologs of training yeasts | report max identity of each S3 protein to the training set | written into the report |

## 5. Model candidates and baselines

| Candidate | Features | Cannot see |
|---|---|---|
| B0 length | length | everything else |
| B1 composition | 20 AA fractions + length, StandardScaler + LR | order, motifs, signal position |
| R rule | SignalP 6 SP call, GPI call (PredGPI module `predgpi/202001`; NetGPI not installed, checked with `module avail`), Ser+Thr fraction | anything not in the three features. SignalP calls 4.6% of *C. immitis* RS proteins (`analysis/cocci_repeats/signalp_summary.tsv`); whether that is an under-call has no truth data |
| M8, M35 | ESM-2 8M layer 6 and 35M layer 6 (as now), residue-mean pooling, StandardScaler + LR | residues past 1,022. The current code keeps the first 1,022 (`embeddings.py:14,163`), so the C-terminal GPI signal of long proteins is lost (164 `surface.tsv` rows are longer) |
| H hybrid | M8 or M35 embedding plus SignalP probability, GPI score and Ser+Thr fraction in one LR | same limits as its parts |

**Why the old CV is useless.** On the old task B1 alone reached ROC-AUC 0.995 under leave-family-out
(review 4.1). A CV that a composition baseline saturates cannot rank models. The new CV must show:
(a) B1 does not saturate it; (b) ML beats B1 and R on N-sec, with the difference outside the bootstrap
CI; (c) this holds under S2. If ML does not beat B1, the report says so.

**Long proteins.** Add a C-terminal window (last 1,022 residues) variant of M8/M35 so the GPI signal is
visible. Report the >1,022 aa subset on its own.

**Rule threshold.** No step 1 Ser+Thr cut-off exists in the repo. The 25% in
`analysis/cocci_repeats/03_repeat_surface_candidates.py` is for class 2a. The cut-off is fitted in the
inner loop.

## 6. Evaluation

- **Metrics per clade** (Saccharomycotina: S288C, *C. albicans*; Eurotiomycetes: literature rows):
  recall, precision, FPR, ROC-AUC, PR-AUC. 95% CI from 2,000 bootstrap resamples, resampled by
  cluster, not by protein.
- **Per stratum:** wall, extracellular-only, N-int, N-sec, ambiguous, >1,022 aa.
- **Precision at fixed recall** (0.8 and 0.9) and **recall at fixed FPR** (0.01).
- **Calibration:** reliability curve and Brier score on S2 test folds.
- **Agreement matrix:** rule call versus ML call on the truth set and on whole S288C, *C. albicans* and
  *C. immitis* RS proteomes. Each discordant cell lists its protein IDs.
- **Calls per proteome** for each candidate (T5).

**Named error-analysis panel.** All IDs were found in repo files on 2026-09-30, except the GO-only IDs.

| Protein | ID | Where found | Why it is in the panel |
|---|---|---|---|
| FLO1 | P32768 | `surface.tsv` | long GPI adhesin, 1,537 aa (truncation case) |
| SAG1 | P20840 | `surface.tsv` | GPI wall protein |
| CWP1 | P28319 | `surface.tsv` | small GPI wall protein |
| CCW12 | Q12127 | `surface.tsv` | 133 aa GPI wall protein |
| GAS1 | P22146 | `surface.tsv` | GPI enzyme |
| PIR1 | Q03178 | `surface.tsv` | wall protein without GPI |
| MSB2 | P32334 | `surface.tsv` | mucin-like plasma-membrane sensor (Q3) |
| HKR1 | P41809 | `surface.tsv` | 1,802 aa mucin-like sensor (Q3) |
| SUC2 | P00724 | `surface.tsv` | secreted enzyme, not wall (Q1) |
| PHO5 | P00635 | `surface.tsv` | secreted enzyme |
| ALS3 | Q59L12 | `surface.tsv` | *C. albicans* GPI adhesin |
| HWP1 | P46593 | `surface.tsv` | *C. albicans* GPI wall protein |
| SAP9 | Q59SU1 | `surface.tsv` | GPI protease |
| ECM33 | A0A1D8PCY4 | `surface.tsv` | `gpi_anchor=no` in UniProt keywords; tests keyword gaps |
| ENO1, TDH3 | CAL0000185645, CAL0000197744 | `cgd.gaf.gz` | IDA cell wall plus cytosol: ambiguous stratum |
| SOWgp58 | Q8NK60 | `surface.tsv` | Onygenales Pro-rich surface protein |
| SOWgp (RS) | CIMG_04613 | local SignalP output | see note below |
| CTS1 | Q1E3R8 | `surface.tsv` | *Coccidioides* secreted chitinase |
| CspA | Q4WXC4 | `surface.tsv` | *A. fumigatus* repeat wall protein |
| cfmA | Q4WLB9 | `surface.tsv` | CFEM GPI protein |
| HSP60 | P50142 | `surface.tsv` | moonlighting; must not be a positive |

Note on SOWgp. `docs/TOOL-ARCHITECTURE.md` says SOWgp "has no call at all". The SignalP 6 output in
`analysis/cocci_repeats/signalp/CimmitisRS_FungiDB/prediction_results.txt` (git-ignored, primary
checkout only) calls CIMG_04613 SP with probability 0.9998. The "no call" statement probably refers
to the earlier annotation-based keywords (assumption, not checked). The SignalP run must be re-done
from a tracked script before the rule is measured.

**Measure-first gating.** (1) Run all candidates on all splits. (2) Write
`results/step1_compare/metrics.json` with point values and CIs. (3) The owner reviews and chooses rule,
ML or hybrid. (4) The owner sets each gate at or just below the measured value. (5) Gates go into
`tests/gates/step1_gates.json`. A regression test compares each new harness run with that file. No
gate value exists before step 3.

## 7. Test framework

| Risk | Test | Type | Runs in |
|---|---|---|---|
| Padding enters the mean (Fable #1) | `test_pool_ignores_padding` | unit | CI |
| Batch order changes scores (T0) | `test_scores_independent_of_batch_order_and_size` (shuffle, batch 1 vs 8, tolerance 1e-5) | unit | CI (8M, CPU) |
| Repeat runs differ (T0) | `test_embedding_deterministic` | unit | CI |
| GPU vs CPU difference (not measured, Fable 5) | `gpu_cpu_diff.py`, report max difference | harness | HPCC |
| Skipped batch misaligns labels (Fable #2) | `test_skip_retry_keeps_alignment` | unit | CI |
| Truncation not counted (Fable #4) | `test_truncation_count_reported`; C-terminal window test | unit | CI |
| Card not enforced (Fable #3, #6) | `test_predict_refuses_card_mismatch` | integration | CI |
| Model loaded from cwd (Fable #31) | `test_model_path_not_from_cwd` | unit | CI |
| Leakage (section 4) | the leakage tests, on fixture tables | unit | CI |
| GO parsing; rule logic | `test_gaf_extract_golden` (50-line fixture); `test_rule_truth_table` | unit | CI |
| Output drift | `test_golden_scores` on a 5-sequence FASTA | integration (`slow`) | CI |
| Accuracy regression | `step1_gates` against `metrics.json` | tier harness T1, T3 to T5 | HPCC |
| Calibration across the kingdom | T6 | tier harness | HPCC, later |

CI uses CPU and Python 3.12 (plan default 2). The 8M weights are cached in CI. Tier T2 (out-of-family)
is replaced by S2 and S3 here, because step 1 has no families. T4 is the N-sec stratum.

## 8. Compute plan

Sizes: S288C 6,722 proteins (plan section 5); *C. immitis* RS 9,910 (`signalp_summary.tsv`);
*C. albicans* not counted. The labelled set is at most about 8,200 genes before dedupe (2.2 plus 2.3).

Review 9.1 has no ESM-2 8M or 35M throughput. It reports ESM-2 150M at 60 proteins/s on one RTX 6000
Ada. Assumption: 8M and 35M are at least that fast, so 30,000 proteins take under 10 minutes per model.
J0 measures the real rate before J1 to J3 are sized.

| Job | Resource | Content | Length |
|---|---|---|---|
| J0 pilot | `exfab`, 1 GPU | 2,000 proteins, 8M and 35M, throughput JSON | under 15 min (assumption) |
| J1 features | `exfab`, 1 GPU | SignalP 6 (GPU build), PredGPI on three proteomes | sized to 1 to 1.5 h |
| J2 embed | `exfab`, 1 GPU | 8M, 35M, C-terminal windows; truth set and three proteomes | sized from J0 to 1 to 1.5 h |
| J3 evaluate | CPU | MMseqs2, CV, bootstrap, metrics JSON | not measured |

Jobs write to `${SCRATCH:?}`, compress large tables with `zstd`, and copy results to
`analysis/step1_compare/` on `/bigdata` before exit. Scripts take `PROJ_ROOT` from the environment, not
from `BASH_SOURCE`. This spec submits no job.

## 9. Risks and unknowns

| Risk | How it is detected |
|---|---|
| GO cell-wall terms include moonlighting enzymes | ambiguous stratum count; ENO1/TDH3 in panel |
| Yeast-only training does not transfer; about 19 Eurotiomycetes proteins give wide CIs | S3 result with CI width; marked smoke test |
| No negatives outside yeasts, so FPR there is unknown | stated in every report; Q6 |
| Truncation hides GPI signals | >1,022 aa stratum |
| SignalP under-calls Onygenales (unproven) | rule recall on literature rows |
| ML learns "signal peptide" only, or composition again explains everything | N-sec FPR of M8/M35 versus R and B1 |
| GAF releases change labels | GAF hash in output; golden extract test |

## 10. Open questions for the owner

1. **Q1. Positive scope.** Recommended: P-ext (wall + secreted, incl. secreted enzymes such as SUC2).
   Alternatives: P-wall only; GPI-CWPs only. Changes: positive count (128 versus 107 in S288C) and
   whether "glycoprotein" in the name is accurate; GPI-CWP-only gives far fewer positives (not counted).
2. **Q2. GPI-anchored plasma-membrane proteins.** Recommended: positive if GPI evidence is curated
   (Swiss-Prot or literature). Alternatives: always negative (in N-sec); exclude. Changes: N-sec
   contents and the source of GPI truth, because GO has no current GPI term.
3. **Q3. Membrane proteins with large extracellular domains** (MSB2, HKR1). Recommended: exclude from
   training, report as a stratum. Alternatives: positive; negative. Changes: whether a TM helix is a
   negative signal for the model.
4. **Q4. UniProt keyword tier (T-c).** Recommended: not used for training in the first run; never test
   truth. Alternative: add T-c positives to training. Changes: about 3,500 more positives, with input
   overlap with the rule.
5. **Q5. GO evidence codes.** Recommended: all non-IEA codes. Alternatives: experimental only (EXP,
   IDA, IMP and the high-throughput codes); include IEA. Changes: *C. albicans* N-int drops from 1,795
   to 168 with experimental codes only (GAF).
6. **Q6. Species list.** Recommended: S288C and *C. albicans* now; add *S. pombe* (PomBase GAF
   available) as a third test species. Alternative: two species only. Changes: S2 gets a species
   outside Saccharomycotina with real negatives.
7. **Q7. Onygenales/Eurotiales rows.** Recommended: test only (S3), never training, reported as a
   smoke test. Alternative: add to training. Changes: S3 becomes impossible if they are in training.
8. **Q8. Clades outside Ascomycota.** Recommended: out of scope for validation; scores outside the
   validated clades carry a "not validated" flag. Alternative: curate Basidiomycota truth first.
9. **Q9. Ambiguous genes.** Recommended: exclude and report. Alternative: treat as positive (GO says
   wall). Changes: moonlighting enzymes enter the positive class.
10. **Q10. Long proteins.** Recommended: add the C-terminal window variant. Alternative: sliding
    windows over the full length. Changes: compute (small at these sizes; not measured).

## 11. Deliverables and definition of done

| ID | Deliverable |
|---|---|
| D1 | `analysis/step1_compare/01_extract_go_truth.py` and `truth_set.tsv.gz` (labels, tier, stratum, cluster, GAF hash) |
| D2 | tracked SignalP 6 and PredGPI job and feature table for the three proteomes |
| D3 | embedding job (8M, 35M, C-terminal windows) with throughput JSON |
| D4 | evaluation script writing `metrics.json` and the agreement matrix with protein IDs |
| D5 | report: per-clade, per-stratum table with CIs and the named-panel table |
| D6 | unit and integration tests in section 7, passing in CI |
| D7 | after owner review: `tests/gates/step1_gates.json` and a model card that reads `metrics.json` |

Done means: D1 to D6 exist and pass; the owner has the comparison; no accuracy statement appears in a
README or card unless it comes from `metrics.json`.

Not in scope: steps 2 and 3; adhesion labels; fine-tuning (review 8.3 gate is not met); ESM C 300M;
Fungi_5k re-scoring; any gate value chosen before measurement.
