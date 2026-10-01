# Design spec: step 1 Phase C, the rule-versus-ML comparison

*2026-10-01. Status: draft for owner review. Parent spec:
`docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md` (called "the parent spec";
sections 4 to 6 and 8 to 11 apply and are not repeated here). No code exists for Phase C. This spec
submits no job.*

## 1. Purpose, scope, non-goals

**Purpose.** Phase A made the truth set. Phase B made the features and embeddings. Phase C answers one
question for the owner: on held-out proteins, how well does the rule (SignalP, GPI call, Ser+Thr) find
cell-wall and secreted proteins, how well do the ESM-2 models find them, and where do they disagree? The
owner then picks rule, ML or hybrid and one training variant (parent spec section 6, "Measure-first
gating").

**Output.** `metrics.json` (values with 95% intervals), `proteome_calls.tsv.gz` (the agreement matrix
and discordant protein IDs), `findings.json` (machine-checked statements), and `report.md` generated
from them.

**In scope.** Training-table construction, clustering and splits (parent section 4), the candidates of
parent section 5 with fixed definitions (section 3 here), evaluation with cluster-bootstrap intervals
(parent section 6), and the tests of parent section 7 that belong to evaluation.

**Not in scope.**

- Gate values. No gate value exists before the owner reads `metrics.json` (parent section 6, last
  paragraph). `tests/gates/step1_gates.json` is deliverable D7, after this phase.
- A shippable model, a model card, a version bump or a release.
- The Basidiomycota truth curation (parent Q8) and the `curated_gpi.tsv` literature curation (R-B).
  Phase C reports what the current truth set supports and says what it cannot.
- Steps 2 and 3, fine-tuning, ESM models larger than 35M.
- Any statement about the whole-proteome prevalence of surface proteins. It is not measured.

## 2. Inputs (all exist; checked 2026-10-01)

| Input | Path under `$STEP1_WORKDIR` | Content |
|---|---|---|
| Labels | `truth_set_triaged.tsv.gz` | label, subset, stratum, `d8_class`, `homology_only`, role, per gene |
| Features | `phaseb/features.tsv.gz` | per member: `ser_thr_frac`, SignalP columns, PredGPI columns, `emb_row`, `emb_cterm_row`, plus label columns for `truth` members |
| Embeddings | `phaseb/emb/<model>.nterm.npy`, `<model>.cterm.npy` | float32; 8M 69,941 x 320 (C-terminal 5,541 x 320); 35M 69,941 x 480 (5,541 x 480) |
| Unique sequences | `phaseb/unique_sequences.tsv.gz` | `seq_sha256`, length, sequence |
| Keyword tier | `keyword_tier.tsv.gz`, `keyword_tier_removed.tsv` | T-c rows (3,092 kept, 448 removed) |
| Provenance | `phaseb/features_run.json`, `phaseb/emb/embedding_run.json`, `d8_run.json`, `keyword_tier_run.json` | hashes, `all_sources` |

**Counts that size the problem** (`truth_set_triaged.tsv.gz`, all non-IEA labels, before sequence
dedupe and direct-evidence filtering; recomputed 2026-10-01):

| Source | P-ext | P-ext, PM-TM | P-ext, pm-unresolved | N-int | N-sec | ambiguous |
|---|---|---|---|---|---|---|
| S288C (train) | 113 | 3 | 9 | 2,645 | 1,548 | 48 |
| *C. albicans* (train) | 199 (+1 P-gpi) | 22 | 37 | 1,772 | 961 | 100 |
| *S. pombe* | 96 | 8 | 12 | 5,928 | 2,426 | 28 |
| *A. fumigatus* | 123 | 1 | 8 | 2,018 | 1,062 | 1 |
| *A. nidulans* | 205 (+1 P-gpi) | 2 | 3 | 2,054 | 1,066 | 38 |
| *U. maydis* | 59 | 1 | 2 | 1,712 | 827 | 1 |
| *C. neoformans* H99 | 7 | 0 | 4 | 23 | 25 | 0 |
| *C. deneoformans* JEC21 | 60 | 2 | 2 | 3,626 | 1,634 | 0 |

The V-go training pool is therefore about 313 positives and about 6,950 negatives from the two
training yeasts (before dedupe and the direct-evidence rule, which applies to test sets only). Direct
evidence counts are smaller (parent section 3.3). H99 has 7 positives: any H99 recall interval is wide.
JEC21 positives are IBA-only (parent section 3.3), so JEC21 is not a headline test set.

## 3. Definitions that Phase C fixes

### 3.1 Class mapping (one function, one test)

`labelmap.class_of(label, d8_class)` returns `pos`, `neg` or `excluded`:

| `label` | `d8_class` | Class | Stratum for reporting |
|---|---|---|---|
| P-ext | empty or P-gpi | pos | `wall` or `extracellular-only` (column `subset`); P-gpi is a list until literature rows exist (R-B) |
| P-ext | PM-TM | neg | PM-TM |
| P-ext | pm-unresolved | excluded | pm-unresolved (list) |
| N-int | any | neg | N-int |
| N-sec | any | neg | N-sec |
| ambiguous | any | excluded | ambiguous (score distribution only) |
| unlabelled | any | excluded | never a negative |

*Ruling C-1: pm-unresolved is excluded.* The parent spec places no class on it. These are P-ext genes
with a plasma-membrane term where D8 found neither curated GPI evidence nor a TM segment, so either
label would be a guess. Cost if wrong: 9 (S288C) and 37 (*C. albicans*) genes missing from training;
the lists stay visible.

### 3.2 Test truth

Headline metrics use **direct-evidence** rows: `homology_only == "no"` and the row's own label (parent
section 3.3). Metrics on all non-IEA labels are reported beside them, in the same table, never instead.
Strata and metrics follow parent section 6. The literature rows (parent Q7) have positives only: report
recall, no precision, no AUC.

### 3.3 Candidates (fixed here so the comparison cannot drift)

| ID | Definition |
|---|---|
| B0 | logistic regression on log(length) |
| B1 | 20 amino-acid fractions + log(length), StandardScaler + logistic regression |
| R0 | rule: SignalP call is `SP` |
| R1 | rule: SP and GPI call at or above class `g` |
| R2 | rule: SP and (GPI call at or above class `g` or `ser_thr_frac` >= `t`) |
| M8, M35 | ESM-2 8M or 35M, N-terminal window embedding (`emb_row`), StandardScaler + logistic regression |
| M8-C, M35-C | as M8 and M35, but proteins longer than 1,022 aa use the C-terminal window (`emb_cterm_row`); shorter proteins use the same row as M8 and M35 |
| H | one logistic regression on [best ESM variant embedding, `sp_prob`, `gpi_prob`, `ser_thr_frac`] |

- **Why R0, R1, R2.** The plan names the three rule inputs but gives no logic. Three nested rules show
  what each input adds. *Ruling C-2: R2 is "the rule"* for the agreement matrix; R0 and R1 are reported
  for context. The owner can change the logic (open question 1).
- **GPI class `g`** takes the values `highly_probable`, `probable`, `weakly` (PredGPI scores 1.0, 0.70,
  0.55; proteins of 40 aa or less are never GPI). **`t`** takes 0.10 to 0.40 in steps of 0.05.
  Both are fitted on training data only (parent section 5, "Rule threshold").
- **M8-C and M35-C** differ from M8 and M35 only on proteins longer than 1,022 aa (5,541 unique
  sequences in the sets of Phase B). The long subset is reported alone (parent section 5).
- **Hyper-parameters.** Logistic regression `C` in {0.01, 0.1, 1, 10}, `class_weight="balanced"`,
  chosen in an inner 3-fold `StratifiedGroupKFold` on training folds only (parent section 4, last
  bullet of "Splits"). Seeds are fixed and recorded.
- **H's ESM variant** is the one of {M8, M35, M8-C, M35-C} with the highest out-of-fold PR-AUC on S1
  with V-go, chosen once and recorded. It is never chosen on S2 or S3 data.
- **Operating point.** Every binary call (rule and ML) uses one criterion: the setting or threshold that
  maximises recall minus FPR (Youden's J) on the out-of-fold training predictions. *Ruling C-3.* The
  criterion does not depend on class prevalence, which is unknown. Cost if wrong: the binary numbers
  move; the threshold-free numbers (ROC-AUC, PR-AUC, precision at fixed recall) do not.
- **Training variants** V-go and V-kw (parent section 2.3). Each trained candidate runs under both.

### 3.4 Datasets, clusters, splits

1. **Table.** One row per unique sequence of the `truth` set and the T-c set, with class, stratum,
   `homology_only`, species, role, length, `seq_sha256`, `emb_row`, `emb_cterm_row`.
2. **Dedupe** by `seq_sha256` as in parent section 4 step 4: keep one copy inside a class, drop a hash
   present in both classes (the existing test `test_dedupe_drops_sequences_present_in_both_classes`
   covers the rule). Dropped hashes are listed.
3. **Clusters.** MMseqs2 `easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0` over all labelled
   proteins of all species and the T-c proteins together (module `MMseqs2/17-b804f`; `mmseqs` is not
   on `PATH` without it). The cluster table is hash-pinned.
4. **S1, S2, S3** as in parent section 4 step 8. S2 and S3 test sets are fixed by species (roles in
   `species.tsv`: train, test_species, test_clade, undecided; alternate files are not test sets).
5. **Max identity to training.** For every S2 and S3 test protein, the highest sequence identity to the
   training set (MMseqs2 `easy-search`, all training proteins as the target). Metrics are reported for
   all test proteins and for those with maximum identity below 0.3. *Ruling C-4.* The parent spec asks
   to report the maximum identity; this turns it into a stratum so a reader sees how much of a result
   comes from near-homologs.
6. **V-kw leakage.** T-c rows are removed when they are a test protein (D10, done in Phase A for the
   fixed test sets), **and** when their cluster contains any protein of the test fold of the current
   split. *Ruling C-5.* The parent spec removes by accession and exact hash only. A homolog of a test
   protein would leak the label. Cost if wrong: fewer V-kw training rows; the removal count per split
   is logged.

## 4. Evaluation

Follows parent section 6. Phase C adds these exact definitions.

- **Bootstrap.** 2,000 resamples. The unit is the cluster: all members of a cluster enter or leave a
  resample together. The same resample indices serve every candidate and both variants, so differences
  are paired. Seed fixed and recorded.
- **Metrics per test set and stratum:** recall, precision, FPR (binary call), ROC-AUC, PR-AUC
  (scored candidates), precision at recall 0.8 and 0.9, recall at FPR 0.01. R0 to R2 have a binary call
  only. For comparison the report gives each ML candidate's precision at the rule's recall and recall at
  the rule's FPR.
- **Estimate or smoke test** (parent section 6, rule fixed 2026-10-01): a test set is an "estimate" when
  the 95% cluster-bootstrap half-width of recall is 0.10 or less **for every candidate**; otherwise it
  is a "smoke test". The classification is computed by the script and stored in `metrics.json`, and
  `report.md` prints it beside every number from that test set.
- **Precision and prevalence.** Training prevalence is not genomic prevalence and the genomic value is
  not measured. Precision is reported at the test set's own prevalence, and in a sensitivity table at
  assumed prevalences of 1%, 3%, 5% and 10%. The table is labelled "assumed, not measured".
- **Calibration** (ML candidates, S2 test sets): reliability curve data and Brier score.
- **Variant effect:** V-kw minus V-go per metric with the paired interval.
- **Agreement.** `proteome_calls.tsv.gz` has one row per protein of every proteome set in
  `sequence_sets.tsv` and one column per candidate call. Proteins that are in a training fold use their
  out-of-fold score (`score_source = oof`); all others use the model trained on the full training pool
  (`score_source = final`). The agreement matrices and discordant lists are views of this table.
- **Named panel** (parent section 6): each listed protein with its label, stratum and every
  candidate's call and score.
- **Findings** (`findings.json`). The script computes the statements the parent spec asks for
  (section 5, "Why the old CV is useless") as booleans with the numbers behind them: (a) B1 does not
  saturate S1 and S2; (b) ML beats B1 and R on N-sec outside the bootstrap interval; (c) (b) also holds
  under S2. If ML does not beat B1, the field says so. `report.md` quotes the field, not free text.

## 5. Code layout

New directory `analysis/step1_compare/eval/`, run with the conda environment Python (numpy, scipy,
scikit-learn). The stdlib-only rule of `test_every_module_imports_only_stdlib_or_local` stays for the
top-level modules; `eval/` gets its own allowed list (numpy, scipy, sklearn, stdlib, local), pinned by a
test that fails on any other import. Scripts continue the numbering:

| Script | Content | Runs on |
|---|---|---|
| `08_build_eval_tables.py` | class mapping, dedupe, training and test tables, input checks | login |
| `09_cluster_and_split.sh` + `09_make_splits.py` | MMseqs2 clustering and max-identity search; S1, S2, S3 fold tables | SLURM CPU |
| `10_fit_and_score.py` | fit R0-R2, B0, B1, M*, H, both variants, all splits; out-of-fold and test scores | SLURM CPU |
| `11_evaluate.py` | metrics, bootstrap, strata, agreement, calibration, findings; writes `metrics.json` | SLURM CPU |
| `12_report.py` | `report.md` from `metrics.json` and `findings.json` | login |

Library modules: `eval/labelmap.py`, `eval/splits.py`, `eval/rules.py`, `eval/models.py`,
`eval/metrics.py`, `eval/bootstrap.py`. Each script follows the Phase B contract: `STOP:` on stderr and
exit 2 with no partial output, temp name then `os.replace`, a run JSON with `all_sources`,
`truth_set_sha256`, input hashes, git commit, library versions, MMseqs2 version and seeds. The first
script checks that `features_run.json` `input_sha256["unique_sequences.tsv.gz"]` equals
`embedding_run.json` `unique_sequences_sha256` (Phase B review item M3) and STOPs otherwise.

**Compute.** No GPU. Sizes: about 45,000 sequences to cluster, about 7,300 training rows per variant,
embedding width 320 or 480. These steps are small. *No run time is measured yet.* The first execution
of each script logs wall time. The scripts form one chain that the plan runs in as few SLURM jobs as
the dependencies allow, on `$SCRATCH`, with results copied to `/bigdata`. They are not split into many
small jobs (global rule on job sizing).

## 6. Test framework (deliverable D6)

Every risk gets a test that can fail. Tests named in the parent spec keep their names.

| Risk | Test | How it fails |
|---|---|---|
| Wrong class for a label | `test_class_of_truth_table` | all `label` x `d8_class` pairs asserted |
| Cluster spans train and test | `test_no_cluster_spans_train_and_test` | plant one spanning cluster in a fixture; the check must reject it |
| IEA label in truth | `test_truth_set_has_no_iea` | fixture with an IEA row |
| T-c is test truth | `test_test_truth_sources` | fixture that routes T-c into a test set |
| T-c holds a test protein or its cluster mate | `test_tc_excludes_test_proteins`, `test_tc_excludes_test_cluster_mates` | plant an accession, a hash and a cluster mate |
| Threshold fitted on test data | `test_threshold_fit_uses_train_only` | change test labels; fitted `t`, `g` and ML threshold must not change |
| Metric formulas wrong | `test_metrics_hand_computed` | recall, precision, FPR, ROC-AUC, PR-AUC on a 10-row case computed by hand |
| Bootstrap unit wrong | `test_bootstrap_resamples_whole_clusters`, `test_bootstrap_seeded` | cluster members move together; same seed, same result |
| Pairing lost | `test_paired_bootstrap_uses_same_indices` | two candidates, same indices |
| A leak inflates CV | `test_shuffled_labels_give_chance_auc` | labels shuffled inside a fixture; out-of-fold AUC must stay near 0.5 |
| Model cannot learn | `test_planted_signal_is_recovered` | synthetic embeddings with a linear signal; AUC must be high |
| Rule logic wrong | `test_rule_truth_table` | R0 to R2 on every SP, GPI, Ser+Thr combination |
| Report states a number that is not in `metrics.json` | `test_report_numbers_come_from_metrics` | alter one number in a copy; the check must fail |
| Result files drift | `test_golden_metrics` | fixed small fixture, frozen expected `metrics.json` |
| Stale inputs | `test_stale_input_stops` | change one input hash; each script must STOP |

## 7. Deliverables and definition of done

| ID | Deliverable |
|---|---|
| E1 | `eval/` library and scripts 08 to 12 with the tests of section 6, passing |
| E2 | `metrics.json`, `proteome_calls.tsv.gz`, `findings.json`, run JSONs for every script |
| E3 | `report.md`: per-species, per-clade, per-stratum tables with intervals; estimate or smoke-test label on every test set; variant effect; calibration; agreement; named panel; what the data cannot show |
| E4 | `analysis/step1_compare/COLUMNS.md` and `README.md` sections for every Phase C output (strict doc tests, as in Phase B) |

Done means: E1 to E4 exist, the full `tests/step1_compare` suite passes, an independent reviewer has
checked the code and one re-run on the real data, and the owner has `report.md`. No accuracy statement
goes into a README or card unless it comes from `metrics.json`.

## 8. Risks

| Risk | Detection |
|---|---|
| Few positives in held-out sets (H99 has 7; direct-evidence sets of 9 to 113 P-ext genes) | estimate or smoke-test label; interval widths |
| The two training yeasts share the same biology; transfer to other clades fails | S2 and S3; metrics split by maximum identity to training |
| Ser+Thr and SignalP outputs explain most of the ML signal | B1, H, R2 in one table; N-sec FPR |
| The rule is defined wrongly (plan gives inputs, not logic) | R0 to R2 nested; owner decides (open question 1) |
| Youden's J picks an operating point with high FPR on N-sec | FPR by stratum is always printed |
| SignalP under-calls *C. immitis* RS (460 of 9,910, 4.6%) and no truth shows whether that is wrong | the Onygenales literature rows (recall only) and the named panel; stated as unverified |
| T-c helps because keywords share inputs with the rule | V-kw minus V-go; T-c never test truth |

## 9. Open questions for the owner

1. **The rule's logic.** R2 = SP and (GPI call at or above `g`, or Ser+Thr at or above `t`). Is this what
   you mean by "SignalP + GPI + Ser/Thr"? Recommended: yes, keep R0 and R1 as context.
2. **Basidiomycota as a test clade.** Phase C scores *U. maydis* (59 P-ext, role undecided) and H99 (7)
   but labels them smoke tests and does not call the result the Q8 validation. Recommended: yes.
3. **Operating point.** Youden's J for all binary calls (ruling C-3). Alternative: fixed recall 0.9.
   Recommended: Youden's J, with precision at recall 0.8 and 0.9 reported beside it.
4. **Rulings C-1 (pm-unresolved excluded), C-4 (identity strata) and C-5 (cluster-mate removal from
   T-c).** Each is cheap to reverse.

## 10. Rulings in this draft

C-1 pm-unresolved excluded; C-2 R2 is "the rule"; C-3 Youden's J operating point; C-4 maximum-identity
strata for S2 and S3; C-5 T-c rows removed when any cluster mate is in the test fold. Each states its
cost if wrong in the section where it appears.
