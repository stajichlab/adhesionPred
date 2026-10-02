# Design spec: step 1 Phase C, the rule-versus-ML comparison

*2026-10-01. Status: revised after the independent review (Fable, 2026-10-01, verdict "fix first"; all
findings addressed below), then after the independent plan review (Fable, 2026-10-01). The owner
answers of 2026-10-01 settled three items: the Ser+Thr grid is 0.20 to 0.40 (section 3.3), the
estimate rule has a floor of 20 direct-evidence positives (section 4, ruling C-8), and the
`hard_negative` literature rows count as positives (section 3.2, ruling C-9). Rulings C-12 to C-15
come from the plan review. Parent spec:
`docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md` ("the parent spec"; its
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
parent section 5 with fixed definitions (section 3), evaluation with cluster-bootstrap intervals
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

| Input | Path | Content |
|---|---|---|
| Labels | `$STEP1_WORKDIR/truth_set_triaged.tsv.gz` | label, subset, stratum, `d8_class`, `homology_only`, `internal_evidence_htp_only`, role, per gene |
| Features | `$STEP1_WORKDIR/phaseb/features.tsv.gz` | per member: `ser_thr_frac`, `sp_prediction`, `sp_prob`, `gpi_call`, `gpi_prob`, `emb_row`, `emb_cterm_row`; label columns for `truth` members |
| Embeddings | `$STEP1_WORKDIR/phaseb/emb/<model>.nterm.npy`, `.cterm.npy` | float32; 8M 69,941 x 320 (C-terminal 5,541 x 320); 35M 69,941 x 480 (5,541 x 480) |
| Unique sequences | `$STEP1_WORKDIR/phaseb/unique_sequences.tsv.gz` | `seq_sha256`, length, sequence |
| Keyword tier | `$STEP1_WORKDIR/keyword_tier.tsv.gz`, `keyword_tier_removed.tsv` | T-c rows (3,092 kept, 448 removed); columns `accession`, `gene`, `genome`, `taxon_id`, `length`, `seq_sha256`, `tier` |
| Literature rows | `data/curated/adhesins/eurotiomycetes_seeds.tsv`, `data/curated/surface/surface.tsv` | Eurotiomycetes seeds (21 rows) and the surface keyword table |
| Species roles | `analysis/step1_compare/species.tsv` | train, test_species, test_clade, undecided, alternate_file |
| Provenance | `phaseb/features_run.json`, `phaseb/emb/embedding_run.json`, `d8_run.json`, `keyword_tier_run.json` | hashes, `all_sources` |

**Counts that size the problem.** Source: `truth_set_triaged.tsv.gz`, grouped by `source_id`, all
non-IEA labels, before sequence dedupe, recomputed 2026-10-01. Sources with role `alternate_file`
(`Spom_SCHPO-mod`, `Cneo_JEC21_GOA`, `Cneo_CRYD1`) are not test sets and are left out.

| Source (role) | P-ext | P-ext, PM-TM | P-ext, pm-unresolved | N-int | N-sec | ambiguous | direct-evidence P-ext (pos class) |
|---|---|---|---|---|---|---|---|
| `Scer_SGD` (train) | 113 | 3 | 9 | 2,645 | 1,548 | 48 | 79 |
| `Calb_CGD` (train) | 199 (+1 P-gpi) | 22 | 37 | 1,772 | 961 | 100 | 154 |
| `Spom_PomBase` (test_species) | 49 | 4 | 6 | 2,962 | 1,214 | 14 | 33 |
| `Afum_ASPFU` (test_clade) | 123 | 1 | 8 | 2,018 | 1,062 | 1 | 19 |
| `Anid_EMENI` (test_clade) | 205 (+1 P-gpi) | 2 | 3 | 2,054 | 1,066 | 38 | 109 |
| `Cneo_H99_GOA` (test_clade) | 7 | 0 | 4 | 23 | 25 | 0 | 7 |
| `Umay_MYCMD` (undecided) | 59 | 1 | 2 | 1,712 | 827 | 1 | 9 |

The V-go training pool is 313 positives and 6,951 negatives (7,264 members, before sequence dedupe).
After the table step (section 3.4) the pool has 308 positive and 6,842 negative unique sequences
(plan prototype, 2026-10-01). The counts differ for these reasons (ruling C-14). Positives:
313 members, minus 2 that share a hash with an ambiguous gene and become `excluded` (1
*C. albicans*, 1 S288C), minus 3 merged members with one hash, give 308. Negatives: 6,951 members,
minus 75 labelled negatives without a row in `features.tsv.gz` (`Calb_CGD` 72, `Scer_SGD` 3; no
sequence, so no table row), minus 34 merged members, give 6,842. Source: `build_run.json`
(`labelled_genes_without_sequence`) and `eval_dedupe_log.tsv` of the plan prototype.
Direct-evidence positives in test sets range from 7 to 109. The JEC21 sources have 30 P-ext genes each
and all of them are IBA-only (parent section 3.3), so they are not headline test sets.

## 2a. Terms used in this spec

| Term | Meaning |
|---|---|
| **PM** | Plasma membrane. The GO term is GO:0005886 (`labels.PLASMA_MEMBRANE`). |
| **PM candidate** | A P-ext gene that also has a non-IEA PM term. Column `pm_candidate == "yes"`. Such a gene is a surface protein in the GO sense (wall or extracellular) and sits in the plasma membrane. It may be a GPI-anchored wall protein (positive) or a transmembrane protein (negative under decision Q3). |
| **TM** | Transmembrane segment, from the UniProt Transmembrane feature (any evidence code). |
| **GPI** | Glycosylphosphatidylinositol anchor. It holds a protein on the outer face of the PM or in the wall. |
| **P-gpi** | PM candidate with curated GPI evidence (reviewed UniProt entry with ECO:0000269, or a literature row). Positive. |
| **PM-TM** | PM candidate with no curated GPI evidence and at least one TM segment (for example MSB2, HKR1). Negative. |
| **pm-unresolved** | PM candidate with neither curated GPI evidence nor a TM segment. Excluded (ruling C-1). |
| **SP** | Signal peptide call by SignalP 6. |
| **P-ext, N-int, N-sec** | Positive (wall or extracellular), internal negative, secretory-pathway negative (parent spec section 2.2). |
| **T-c** | Keyword-only training tier: proteins labelled surface only by UniProt keywords, not by GO (3,092 rows kept). Training only, never test truth (parent spec section 2.3). |
| **V-go** | Training variant 1. Positives are GO-labelled proteins (P-ext and P-gpi). No T-c rows. |
| **V-kw** | Training variant 2. V-go plus the T-c rows that survive the leakage controls of section 3.4. Both variants are scored on the same GO test truth. |
| **`C`** | Regularisation strength of the logistic regression (smaller means stronger regularisation). Values tried: 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.1, 1, 10 (rulings C-12 and C-17). Chosen on the inner folds by the highest inner out-of-fold PR-AUC. |
| **`g`** | GPI class cut for the rule. A protein counts as GPI-anchored when its PredGPI class is at or above `g`. Values: `highly_probable`, `probable`, `weakly`. Fitted on the inner folds. |
| **`t`** | Ser+Thr cut for the rule. A protein counts as Ser/Thr-rich when `ser_thr_frac` is at or above `t`. Values: 0.20 to 0.40 in steps of 0.05 (owner, 2026-10-01). Fitted on the inner folds. |
| **R0, R1, R2** | The three nested rules of section 3.3. R2 (SP and either GPI at or above `g` or Ser+Thr at or above `t`) is "the rule". |
| **B0, B1** | Baselines: length only; amino-acid composition plus length. |
| **M8, M35; M8-C, M35-C** | ESM-2 8M and 35M embeddings with logistic regression. The "-C" versions use the last 1,022 residues for proteins longer than 1,022 aa. |
| **H** | Hybrid: embedding plus SignalP probability, GPI score and Ser+Thr in one logistic regression. |
| **S1, S2, S3** | Split schemes: homology-grouped 5-fold cross-validation; leave-species-out; leave-clade-out (parent spec section 4). |
| **Inner / outer folds** | Outer folds give the test result. Inner folds sit inside each outer training set and pick `C`, `g`, `t` and thresholds, so the test labels never influence a choice. |
| **Out-of-fold (oof) score** | A score for a protein from a model that did not train on that protein or its cluster. |
| **Youden's J** | Recall minus false-positive rate. The operating point is the setting that maximises it. |
| **Direct evidence** | A label whose supporting GO codes are not homology-transfer codes (`homology_only == "no"`). Used for headline test numbers. |
| **Estimate / smoke test** | A test set is an "estimate" when the 95% interval of recall has half-width 0.10 or less for every listed candidate and the set has at least 20 direct-evidence positives; otherwise it is a "smoke test" (section 4, ruling C-8). |

## 3. Definitions that Phase C fixes

### 3.1 Class mapping (one function, one test)

`labelmap.class_of(label, d8_class)` returns `pos`, `neg` or `excluded`:

| `label` | `d8_class` | Class | Stratum for reporting |
|---|---|---|---|
| P-ext | empty or P-gpi | pos | column `subset`: `wall` or `extracellular-only`; P-gpi is a list until literature rows exist (R-B) |
| P-ext | PM-TM | neg | PM-TM |
| P-ext | pm-unresolved | excluded | pm-unresolved (list) |
| N-int | any | neg | N-int |
| N-sec | any | neg | N-sec |
| ambiguous | any | excluded | ambiguous (score distribution only); sub-stratum `internal_evidence_htp_only` |
| unlabelled | any | excluded | never a negative |

Read `subset`, not `stratum`, for positives: `03_triage_pm.py` writes `d8_class` over `stratum` for PM
candidates (COLUMNS.md, `truth_set_triaged.tsv.gz`).

*Ruling C-1: pm-unresolved is excluded.* The parent spec places no class on it, and the `d8_triage.py`
docstring says it "stays P-ext"; this ruling overrides that for Phase C. These are P-ext genes with a
PM term where D8 found neither curated GPI evidence nor a TM segment, so either label would be a guess.
Cost if wrong: 9 (S288C) and 37 (*C. albicans*) genes missing from training; the lists stay visible.

### 3.2 Test truth

Headline metrics use **direct-evidence** rows: `homology_only == "no"` and the row's own label (parent
section 3.3). The filter applies to test rows only, for positives and negatives. Metrics on all non-IEA
labels are reported beside them, in the same table, never instead.

**Literature rows (parent Q7).** The file `eurotiomycetes_seeds.tsv` has a `class` column that describes
adhesion (`adhesin` 12, `hard_negative` 9), not surface location. Step 1 therefore defines literature
positives as the rows that have an accession, a sequence in the keyword set and `moonlighting` not
`YES`. The `adhesin` and the `hard_negative` rows both count (owner, 2026-10-01): 19 literature
positives in the Phase A data. The script logs the count. HSP60 (P50142) is scored but is never a
positive. *Ruling C-9 (amended 2026-10-01): `hard_negative` rows are positives, not negatives* (the
label says "not an adhesin", which does not mean "not on the surface"). The named panel also lists
every `hard_negative` row with its scores. Literature rows have positives only: report recall only.
`metrics.json` stores precision, PR-AUC, precision at a recall level and precision at the rule's
recall as `null` for this set, and `report.md` does not print them.

### 3.3 Candidates (fixed here so the comparison cannot drift)

| ID | Definition |
|---|---|
| B0 | logistic regression on log(length) |
| B1 | 20 amino-acid fractions + log(length), StandardScaler + logistic regression |
| R0 | rule: SignalP call is `SP` |
| R1 | rule: SP and GPI call at or above class `g` |
| R2 | rule: SP and (GPI call at or above class `g` or `ser_thr_frac` >= `t`) |
| M8, M35 | ESM-2 8M or 35M, N-terminal window embedding (`emb_row`), StandardScaler + logistic regression |
| M8-C, M35-C | as M8 and M35, but proteins longer than 1,022 aa use the C-terminal window (`emb_cterm_row`); shorter proteins use `emb_row`, as in M8 and M35 |
| H | one logistic regression on [embedding of the chosen ESM variant, `sp_prob`, `gpi_prob`, `ser_thr_frac`] |

- **Rule (owner, 2026-10-01).** A protein needs a signal peptide. It is then a surface candidate if it
  has a GPI call at or above class `g`, or a Ser+Thr fraction at or above `t`. This is R2 and it is
  "the rule" in the agreement matrix (ruling C-2). R0 and R1 are context.
- **GPI class `g`** takes `highly_probable`, `probable`, `weakly` (PredGPI scores 1.0, 0.70, 0.55;
  proteins of 40 aa or less are never GPI). **`t`** takes 0.20 to 0.40 in steps of 0.05 (owner,
  2026-10-01). The plan review found that 96.7% of the SP proteins have `ser_thr_frac` at or above
  0.10, so at `t` = 0.10 R2 equals R0 on the real data; 0.10 and 0.15 are therefore not in the grid.
- **M8-C and M35-C** differ from M8 and M35 only on proteins longer than 1,022 aa (5,541 unique
  sequences in the Phase B sets). The long subset is reported alone (parent section 5).
- **Hyper-parameters.** Logistic regression `C` in {0.0001, 0.0003, 0.001, 0.003, 0.01, 0.1, 1, 10},
  `class_weight="balanced"`. `C` is the value with the highest inner out-of-fold PR-AUC (ties: the
  smaller `C`). *Ruling C-12:* the grid goes down to 0.001, because on the real S1 fold 0 every
  logistic-regression candidate chose 0.01, the lowest value of the old grid. The choice stays on the
  inner folds, so the extension adds no leakage. It costs two more `C` values per candidate and unit
  (six more inner fits; for H, eight more (variant, `C`) settings).
  *Ruling C-17 (owner, 2026-10-02):* the grid goes down to 0.0001, because in the first real run (grid
  0.001 to 10) every logistic-regression candidate chose 0.001, the new lower edge. The choice stays on
  the inner folds.
- **Nested protocol (one rule for every fitted quantity).** The fitted quantities are: scaler, `C`,
  `g`, `t`, the ML decision threshold and H's ESM variant. For each outer training set they are fitted
  on the **inner** 3-fold `StratifiedGroupKFold` (groups = clusters) of that training set only. The
  decision threshold and `g`, `t` use the inner out-of-fold predictions. For S1 the outer training set
  is the four training folds of each outer fold. For S2 and S3 it is the whole training pool.
  Outer test-fold or test-species labels are never read before the final scoring.
- **H's ESM variant** is one of {M8, M35, M8-C, M35-C}, chosen by inner out-of-fold PR-AUC inside each
  outer training set, and recorded per outer fold. The S1 numbers of H are then free of selection
  bias from the outer folds.
- **Operating point.** Every binary call (rule and ML) uses one criterion: the setting or threshold that
  maximises recall minus FPR (Youden's J) on the inner out-of-fold predictions (ruling C-3). The
  criterion does not depend on class prevalence, which is unknown. Cost if wrong: the binary numbers
  move; the threshold-free numbers do not.
- **Training variants** V-go and V-kw (parent section 2.3). Each trained candidate runs under both.

### 3.4 Datasets, clusters, splits

1. **Table.** Start from the `truth` members whose source role is not `alternate_file` and whose label
   is not `unlabelled`, and the T-c rows. One row per unique sequence, with class, `subset`,
   `internal_evidence_htp_only`, `homology_only`, source, species, role, length, `seq_sha256`,
   `emb_row`, `emb_cterm_row`. Match by `seq_sha256`: proteome IDs differ from truth IDs.
2. **Dedupe and precedence.** Inside the GO truth: keep one copy of a hash within a class and drop a
   hash present in both classes (parent section 4 step 4; no such hash exists in the Phase A data, so
   the rule is a guard). *Ruling C-14:* a hash with a positive or negative member and an excluded
   member (for example a direct-evidence positive that shares its sequence with an ambiguous gene)
   becomes one `excluded` row; in the Phase A data this applies to 2 hashes (1 *C. albicans*, 1
   S288C positive). **A T-c row is dropped when its hash has any GO row (positive, negative or
   excluded) in the training table; the GO label wins** (ruling C-6). In the Phase A data 224 GO
   negatives (207 N-sec, 8 N-int, 9 PM-TM) and 65 excluded genes (27 pm-unresolved, 38 ambiguous) share
   a hash with a kept T-c row (reviewer count, to be re-checked by the script). Without this rule V-kw
   would lose those negatives and would admit excluded genes as positives. The script logs the counts
   per class. The dedupe function is implemented in `phasec/` with its own test; it does not import the
   `surface_glyco` package function.
3. **Clusters.** MMseqs2 `easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0` over all table sequences
   (module `MMseqs2/17-b804f`; `mmseqs` is not on `PATH` without it). The cluster table is hash-pinned.
4. **S1, S2, S3** as in parent section 4 step 8. S2 and S3 test sets are fixed by species role in
   `species.tsv`: training species have role `train`; test sets are `test_species`, `test_clade` and
   `undecided`. Alternate files are never test sets.
5. **T-c leakage controls for V-kw** (in addition to D10 of Phase A):
   a. T-c rows that are test proteins (by accession or hash) are removed.
   b. *Ruling C-5:* T-c rows whose cluster contains any protein of the test fold (S1) or any test
      protein (S2, S3) are removed.
   c. *Ruling C-7:* for S2 and S3 the T-c rows whose `taxon_id` is the test species (S2) or lies in the
      test clade (S3) are removed. In the Phase A data the kept T-c rows hold, for example, 666
      *A. fumigatus*, 413 *C. albicans*, 350 S288C, 555 *C. immitis* and 511 *C. posadasii* proteins
      (reviewer count). Without this rule "train both yeasts, test Eurotiomycetes" under V-kw would
      train on 666 *A. fumigatus* keyword positives. The clade of each T-c taxon comes from a table
      `phasec/tc_taxon_clades.tsv` (one row per `taxon_id` in `keyword_tier.tsv.gz`, with a clade or
      `other`). A test fails when a `taxon_id` has no row.
   The number of rows removed by each rule is logged per split.
6. **Maximum identity to training** (S2 and S3 test proteins). MMseqs2 `easy-search` with `-s 7.5`,
   query = test protein, target = the GO training proteins of the split (the same targets under V-go
   and V-kw, so the stratum does not depend on the variant). Identity is `fident`, coverage as in the
   clustering (`-c 0.5 --cov-mode 0`). No hit means below 0.3. A protein is never compared with
   itself. Metrics are reported for all test proteins and for those with identity below 0.3
   (ruling C-4).

## 4. Evaluation

Follows parent section 6. Phase C adds these definitions.

- **Bootstrap.** 2,000 resamples. The unit is the cluster: all members of a cluster enter or leave a
  resample together. The same resample indices serve every candidate and both variants, so differences
  are paired. Seed fixed and recorded. The intervals describe sampling of the test set only. The
  models are fitted once per split and are not refitted inside the bootstrap. This is stated beside the
  variant-effect interval in the report.
- **Metrics per test set and stratum:** recall, precision, FPR (binary call), ROC-AUC, PR-AUC (scored
  candidates), precision at recall 0.8 and 0.9, recall at FPR 0.01. Precision at a recall level reads
  the step-function precision-recall curve (the highest precision at any threshold with recall at or
  above the level). R0 to R2 have a binary call only. For comparison the report gives each ML
  candidate's precision at the rule's recall and recall at the rule's FPR.
- **Estimate or smoke test** (parent section 6, written before any run). Truth = direct-evidence rows.
  Variant = V-go. Candidates = R2 and every ML candidate (M8, M35, M8-C, M35-C, H); B0, B1, R0 and R1
  are baselines and are left out. The half-width is (97.5th percentile minus 2.5th percentile) / 2 of
  recall. A test set is an "estimate" when the half-width is 0.10 or less for every listed candidate
  **and** the set has at least 20 direct-evidence positives; otherwise it is a "smoke test" (ruling
  C-8, count floor added by the owner on 2026-10-01). The floor exists because a percentile interval
  has zero width when every resample gives the same recall (for example 7 of 7 positives found), and
  such a set would otherwise pass. `metrics.json` records, per test set, `n_direct_positives`,
  `floor_met` and `zero_width_recall_interval` (the label candidates whose recall interval has lo
  equal to hi) beside the label. The script also prints the half-width of FPR beside the label; it
  does not enter the label, because the parent rule is about recall. The classification is stored in
  `metrics.json` and `report.md` prints it beside every number from that test set.
- **Precision and prevalence.** Training prevalence is not genomic prevalence and the genomic value is
  not measured. Precision is reported at the test set's own prevalence, and in a sensitivity table at
  assumed prevalences of 1%, 3%, 5% and 10%, labelled "assumed, not measured" (ruling C-11; it replaces
  the parent's "prior estimated on whole proteomes", which no measurement supports).
- **Calibration** (ML candidates, S2 test sets). `class_weight="balanced"` moves probabilities toward
  0.5 prevalence by construction. Probabilities are therefore calibrated by Platt scaling fitted on the
  inner out-of-fold predictions of the training pool only (ruling C-10). The report gives reliability
  data and the Brier score after this step, and says so.
- **Variant effect:** V-kw minus V-go per metric with the paired interval.
- **Agreement.** `proteome_calls.tsv.gz` has one row per protein of every proteome set in
  `sequence_sets.tsv`, joined to the labelled tables by `seq_sha256`. It has one call column and one
  score column per candidate and variant. `score_source` is `oof` when the protein's hash is in an S1
  fold (T-c rows are assigned to S1 folds by cluster, so they also get out-of-fold scores), `final` when
  it is not in any training table, and `in_sample` for any remaining training protein. The agreement
  matrices and discordant lists are views of this table.
- **Named panel** (parent section 6): each listed protein, found by accession in the keyword set and
  joined to truth rows by hash, with its label, stratum and every candidate's call and score. The
  `hard_negative` literature rows are listed here (ruling C-9).
- **Findings** (`findings.json`). The script computes booleans and the numbers behind them:
  (a) B1 does not saturate S1 or S2: its ROC-AUC is below 0.99 (the old saturation value was 0.995,
  parent section 5) on direct-evidence truth; the value is printed.
  (b) ML beats B1 and R2 on N-sec: for V-go on direct-evidence truth, the paired cluster-bootstrap 95%
  interval of FPR_N-sec(comparator) minus FPR_N-sec(ML), both read at the recall of R2, excludes 0.
  (c) the same statement holds under S2.
  *Ruling C-13:* no headline ML candidate is named. (b) and (c) are evaluated for each ML candidate
  (M8, M35, M8-C, M35-C, H), and `holds_for` lists the candidates for which the statement holds; (c)
  holds for a candidate only when it holds on all three S2 test sets. "ML beats" therefore means "any
  listed ML candidate beats", and `report.md` states this and prints the per-candidate table.
  If ML does not beat B1, the field says so. `report.md` quotes the field, not free text.
- **Context note.** The T-c rows for *C. immitis* (555) and *C. posadasii* (511) stay in V-kw training
  for S1 and S2. The report prints this beside the Onygenales literature recall and the *C. immitis* RS
  agreement matrix, because both are read as held-out numbers.

## 5. Code layout

New directory `analysis/step1_compare/phasec/` (not `eval/`, which shadows a built-in name), run with
the conda environment Python (numpy, scipy, scikit-learn). The existing test
`test_every_module_imports_only_stdlib_or_local` scans only `analysis/step1_compare/*.py`. Scripts 08 to
12 therefore live **inside** `phasec/` so that numpy imports do not fail it. `phasec/` gets its own
allowed list (numpy, scipy, sklearn, stdlib, local), pinned by a test that fails on any other import.
The scripts import the top-level modules (`truth_table`, `labels`, `runinfo`) through `PYTHONPATH`, set
as in the README for `jobs/`.

| Script | Content | Runs on |
|---|---|---|
| `phasec/08_build_eval_tables.py` | class mapping, dedupe and precedence, tables, `tc_taxon_clades.tsv` check, input checks | login |
| `phasec/09_cluster_and_split.sh` + `phasec/09_make_splits.py` | MMseqs2 clustering and max-identity search; S1, S2, S3 fold tables; T-c removals | SLURM CPU |
| `phasec/10_fit_and_score.py` | fit R0-R2, B0, B1, M*, H, both variants, all splits; out-of-fold and test scores | SLURM CPU |
| `phasec/11_evaluate.py` | metrics, bootstrap, strata, agreement, calibration, findings; writes `metrics.json` | SLURM CPU |
| `phasec/12_report.py` | `report.md` from `metrics.json` and `findings.json` | login |

Library modules: `phasec/labelmap.py`, `splits.py`, `dedupe.py`, `rules.py`, `models.py`, `metrics.py`,
`bootstrap.py`. Each script follows the Phase B contract: `STOP:` on stderr and exit 2 with no partial
output, temp name then `os.replace`, a run JSON with `all_sources`, `truth_set_sha256`, input hashes,
git commit, library versions, MMseqs2 version and seeds. The first script checks that
`features_run.json` `input_sha256["unique_sequences.tsv.gz"]` equals `embedding_run.json`
`unique_sequences_sha256` (Phase B review item M3) and STOPs otherwise.

**Compute.** No GPU. Sizes: about 23,400 unique sequences to cluster (labelled non-alternate truth plus
T-c; reviewer count 23,413, the script logs the exact value), 7,264 V-go training rows (before dedupe)
and up to 3,092 more T-c rows for V-kw before removals, embedding width 320 or 480. These steps are
small. *No run time is measured yet.* The first execution of each script logs wall time. *Ruling
C-15:* the SLURM job requests partition `epyc` with `--constraint=ryzen`. The MMseqs2 module build
needs AVX2 and stops with exit 132 on the abu_dhabi Opterons of partition `batch`; the epyc nodes
carry the features ryzen, amd and milan. A job that runs an AVX2 tool must request a node feature
that has AVX2; a test scans every shell script of `analysis/step1_compare/` for this. The job
copies its wall-time lines to the log directory on every exit path, also after a failed step. The plan runs
the chain in as few SLURM jobs as the dependencies allow, on `$SCRATCH`, with results copied to
`/bigdata`. They are not split into many small jobs (global rule on job sizing).

## 6. Test framework (deliverable E1)

Every risk gets a test that can fail. Tests named in the parent spec keep their names.

| Risk | Test | How it fails |
|---|---|---|
| Wrong class for a label | `test_class_of_truth_table` | all `label` x `d8_class` pairs asserted |
| Alternate-file sources enter the table | `test_alternate_files_excluded` | fixture with an alternate source |
| T-c row overrides a GO label | `test_tc_row_yields_to_go_label` | T-c hash equal to an N-sec, a pm-unresolved and an ambiguous gene; the GO class must stay and no T-c positive may appear |
| Cluster spans train and test | `test_no_cluster_spans_train_and_test` | plant one spanning cluster; the check must reject it |
| IEA label in truth | `test_truth_set_has_no_iea` | fixture with an IEA row |
| T-c is test truth | `test_test_truth_sources` | fixture that routes T-c into a test set |
| T-c holds a test protein or its cluster mate | `test_tc_excludes_test_proteins`, `test_tc_excludes_test_cluster_mates` | plant an accession, a hash and a cluster mate |
| T-c holds the test species or clade | `test_tc_excludes_test_taxa` | plant rows of the test taxon in S2 and S3 |
| A T-c taxon has no clade | `test_tc_taxon_table_is_complete` | remove one row from the table |
| Threshold fitted on test data | `test_threshold_fit_uses_train_only` | permute the outer test-fold labels; fitted `t`, `g`, scaler, `C`, H variant and ML threshold must not change |
| Metric formulas wrong | `test_metrics_hand_computed` | recall, precision, FPR, ROC-AUC, PR-AUC on a 10-row case computed by hand; precision at fixed recall on a case with ties |
| Bootstrap unit wrong | `test_bootstrap_resamples_whole_clusters`, `test_bootstrap_seeded` | cluster members move together; same seed, same result |
| Pairing lost | `test_paired_bootstrap_uses_same_indices` | two candidates, same indices |
| Test-set contamination inflates CV | `test_shuffled_labels_give_chance_auc` | labels shuffled inside a fixture; out-of-fold AUC must fall within 0.5 +/- 0.05 for a fixed seed. The test detects contamination of test folds. It cannot detect homology leakage, which `test_no_cluster_spans_train_and_test` covers |
| Model cannot learn | `test_planted_signal_is_recovered` | synthetic embeddings with a linear signal; AUC must be high |
| Rule logic wrong | `test_rule_truth_table` | R0 to R2 on every SP, GPI, Ser+Thr combination, with a value at the cut-off |
| C-terminal row wrong | `test_cterm_variant_row_selection` | proteins at 1,022 and 1,023 aa; row choice asserted |
| Direct-evidence filter wrong | `test_direct_filter_applies_to_test_rows_only` | training rows with `homology_only == yes` stay in training; test rows are dropped |
| Estimate label wrong | `test_estimate_label_rule` | planted recall intervals around 0.10; a zero-width interval with 7 positives must be a smoke test; small widths with 20 positives must be an estimate |
| Literature set prints precision | `test_undefined_metrics_are_null`, `test_literature_section_reports_recall_only` | precision, PR-AUC and precision at recall of the literature set must be `null`; its report section must hold no precision column |
| Recall at FPR wrong when the top score is a negative | `test_recall_at_fpr` | the empty call set must give recall 0.0, not NaN |
| Wall times lost when a C1 step fails | `test_c1_keeps_the_wall_times_after_a_failed_step` | step 10 stops; `wall.<job>.txt` must hold the step 09 line |
| AVX2 tool on a node without AVX2 | `test_avx2_tools_have_a_cpu_constraint` | delete the `--constraint` line of a script that runs MMseqs2 |
| Report caveats missing | `test_report_states_the_caveats` | remove a caveat line from the report |
| Prevalence table wrong | `test_prevalence_table_formula` | hand-computed precision at assumed prevalence |
| Max-identity stratum wrong | `test_max_identity_excludes_self_hits` | a test protein that is also in the target set |
| Proteome score source wrong | `test_proteome_calls_use_oof_for_training_hashes` | a training hash must never get `final` |
| Findings wrong | `test_findings_booleans` | fixture where B1 saturates and where it does not |
| Report states a number that is not in `metrics.json` | `test_report_numbers_come_from_metrics` | alter one number in a copy; the check must fail |
| Result files drift | `test_golden_metrics` | fixed small fixture, frozen expected `metrics.json`; compared with a tolerance of 1e-9, library versions pinned in the run JSON |
| Stale inputs | `test_stale_input_stops` | change one input hash; each script must STOP |
| Import rule | `test_phasec_imports_only_allowed` | add a forbidden import to a copy |

## 7. Deliverables and definition of done

| ID | Deliverable |
|---|---|
| E1 | `phasec/` library and scripts 08 to 12 with the tests of section 6, passing |
| E2 | `metrics.json`, `proteome_calls.tsv.gz`, `findings.json`, run JSONs for every script |
| E3 | `report.md`: per-species, per-clade, per-stratum tables with intervals; estimate or smoke-test label on every test set; variant effect; calibration; agreement; named panel; what the data cannot show |
| E4 | `analysis/step1_compare/COLUMNS.md` and `README.md` sections for every Phase C output (strict doc tests, as in Phase B) |

Done means: E1 to E4 exist, the full `tests/step1_compare` suite passes, an independent reviewer has
checked the code and one re-run on the real data, and the owner has `report.md`. No accuracy statement
goes into a README or card unless it comes from `metrics.json`.

## 8. Risks

| Risk | Detection |
|---|---|
| Few positives in held-out sets (direct-evidence P-ext: 7 to 109) | estimate or smoke-test label; interval widths |
| The two training yeasts share the same biology; transfer to other clades fails | S2 and S3; metrics split by maximum identity to training |
| Ser+Thr and SignalP outputs explain most of the ML signal | B1, H, R2 in one table; N-sec FPR |
| The rule is defined wrongly (the plan gave inputs, not logic) | R0 to R2 nested; the owner confirmed R2 |
| Youden's J picks an operating point with high FPR on N-sec | FPR by stratum is always printed |
| SignalP under-calls *C. immitis* RS (460 of 9,910, 4.6%) and no truth shows whether that is wrong | the Onygenales literature rows (recall only) and the named panel; stated as unverified |
| T-c helps because keywords share inputs with the rule | V-kw minus V-go; T-c never test truth; precedence rule C-6 |
| Recall passes the estimate rule while FPR rests on 28 to 44 negatives (H99, *U. maydis*, *A. fumigatus*) | FPR half-width printed beside the label |
| A zero-width recall interval of a small set passes the half-width rule | count floor of 20 direct-evidence positives (ruling C-8); `zero_width_recall_interval` printed |
| Pooled S1 ROC-AUC and PR-AUC mix decision values from five fold models, each with its own score scale (plan review M-3) | fixed caveat line in `report.md`; the per-fold diagnostic is deferred |
| The ML threshold and the Platt scaling are fitted on the inner out-of-fold predictions and applied to a model refitted on all training rows, which can have another score scale (plan review M-4) | fixed caveat line in `report.md`; not measured |
| The report-number check is set membership: a number can come from `metrics.json` and still sit in the wrong cell (plan review M-7) | fixed caveat line in `report.md`; the tables are generated from the same keys |

## 9. Decisions (owner, 2026-10-01)

| Question | Decision | Consequence |
|---|---|---|
| 1. Rule logic | A protein must have a signal peptide. It is then a candidate if it has a GPI call at or above class `g`, or Ser+Thr at or above `t` | R2 is "the rule" (ruling C-2). R0 and R1 stay as context |
| 2. Basidiomycota | *U. maydis* and H99 are scored as smoke tests. They are not the Q8 validation | Section 4 label rule; the report says so |
| 3. Operating point | Recall minus FPR (Youden's J) for all binary calls | Ruling C-3. Precision at recall 0.8 and 0.9 is reported beside it |
| 4. Rulings C-1, C-4, C-5 | Accepted | pm-unresolved excluded; maximum-identity strata; T-c cluster-mate removal |
| 5. Ser+Thr grid (after the plan review) | `t` from 0.20 to 0.40 in steps of 0.05 | Section 3.3; 0.10 and 0.15 are dropped |
| 6. Estimate rule (after the plan review) | Add a floor: at least 20 direct-evidence positives | Section 4; ruling C-8 amended |
| 7. `hard_negative` literature rows (after the plan review) | They count as positives (19 literature positives); HSP60 is never positive | Section 3.2; ruling C-9 amended |

Open items outside Phase C (Basidiomycota truth curation, `curated_gpi.tsv` curation, merge order of the
open PRs) are unchanged.

## 10. Rulings

| ID | Ruling | Origin | Cost if wrong |
|---|---|---|---|
| C-1 | pm-unresolved excluded | draft | 9 (S288C) and 37 (*C. albicans*) genes missing from training |
| C-2 | R2 is "the rule" | owner | none; R0 and R1 are reported |
| C-3 | Youden's J operating point | owner | binary numbers move; threshold-free numbers do not |
| C-4 | maximum-identity strata for S2 and S3 | owner | one extra column and one extra table |
| C-5 | T-c rows removed when a cluster mate is a test protein | owner | fewer V-kw training rows |
| C-6 | A GO label wins over a T-c row with the same hash | review C1 | V-kw would otherwise lose 224 negatives and admit 65 excluded genes |
| C-7 | T-c rows of the test species or clade are removed for S2 and S3 | review C2 | fewer V-kw rows; without it S2 and S3 are not leave-species-out under V-kw. The alternative is to keep them and label S2 and S3 V-kw as "same-species keyword rows in training" |
| C-8 | Estimate rule: candidates R2 and ML, direct truth, V-go, recall half-width at most 0.10, and at least 20 direct-evidence positives (floor added by the owner, 2026-10-01) | review I7; owner | label changes for a few test sets; the parent rule did not name the candidates; the direct-evidence positives of H99 (7), *U. maydis* (9), *A. fumigatus* (19) and the literature set (19) are below the floor, so these sets are smoke tests |
| C-9 | `hard_negative` literature rows are positives, not negatives (amended by the owner, 2026-10-01); HSP60 is never positive | review I3; owner | literature recall is computed on 19 rows, not on the 10 `adhesin` rows; the named panel lists the `hard_negative` rows |
| C-10 | ML probabilities calibrated by Platt scaling on inner out-of-fold predictions | review I5 | calibration plots differ; ranking metrics do not |
| C-11 | Assumed prevalence grid replaces the parent's whole-proteome prior | review M5 | none; no measurement of the prior exists |
| C-12 | Logistic-regression `C` grid {0.001, 0.003, 0.01, 0.1, 1, 10}, chosen by inner out-of-fold PR-AUC | plan review (S1 fold 0 chose the old edge 0.01) | two more `C` values per candidate and unit (six more inner fits); no leakage, the choice stays on the inner folds |
| C-13 | No headline ML candidate is named; findings (b) and (c) report `holds_for` per candidate, and the report says "any" | plan review | the owner reads a per-candidate table |
| C-14 | A positive or negative that shares a sequence with an excluded gene becomes `excluded`; the pool counts 308 / 6,842 differ from 313 / 6,951 by this rule, by merges and by 75 labelled negatives without a feature row | plan review | 2 positives missing from training (1 *C. albicans*, 1 S288C) |
| C-15 | The SLURM job requests `-p epyc --constraint=ryzen`; every job that runs an AVX2 tool requests an AVX2 node feature | plan review; owner suggestion | fewer eligible nodes (19 epyc nodes) |
| C-17 | The logistic-regression `C` grid goes down to 0.0001: 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.1, 1, 10 | owner, 2026-10-02 | in the first real run every ML candidate chose 0.001, the old lower edge; two more `C` values per candidate and unit |

C-6 to C-11 came from the first independent review. The owner asked on 2026-10-01 to go forward with
the plan; C-6 to C-11 were not reviewed one by one.
C-12 to C-15 came from the plan review (2026-10-01). The owner answers of 2026-10-01 amended C-8
and C-9.
