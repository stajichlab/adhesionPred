# Adhesion predictor: code/model review, validation experiments, and a framework plan

*2026-09-27. Review by Claude Code (Opus 5.5) at the request of J. Stajich.
Experiments are reproducible with `analysis/model_review/run.sh`.*

## 1. Bottom line

1. **The current classifier discriminates its training classes almost perfectly, but that benchmark is saturated and uninformative.**
   Homology-aware and leave-family-out cross-validation *do not* reduce performance
   (ROC-AUC ≈ 0.999). A logistic regression on **20 amino-acid frequencies + length,
   with no language model at all**, reaches ROC-AUC 0.995 under leave-family-out CV.
   The task "FLO/ALS vs. random proteins from three proteomes" reduces to "long,
   Ser/Thr-rich secreted protein vs. typical protein". Better CV alone will not tell
   you whether the kingdom-wide calls are right.
2. **Genome-wide, the model over-calls by roughly 3–10×.** In the Fungi_5k screen (the
   `esm2_t12_35M` model) the median species has 113 calls (≈1.3% of the proteome). *S. pombe*
   has 38 and *C. neoformans* has 60, even though 500 proteins from each of those proteomes
   were training negatives. *N. crassa* has 136 and *A. fumigatus* 114. Published adhesin inventories are roughly 10–30 per genome
   in most species (larger in *Nakaseomyces glabratus* and *C. auris*). The top CAZy family
   among the calls is AA1 multicopper oxidases/laccases (12.5k calls), which are not adhesins.
   On *S. cerevisiae* S288C the shipped 8M model makes 77 calls. It recovers 8/9 known adhesins,
   but **precision is ≈12%**: the other calls are GPI cell-wall mannoproteins (SED1, TIR4, DAN1/4, CCW12…),
   mucin-like sensors (MSB2, HKR1, WSC2–4) and dubious ORFs. It also calls 12/20 curated hard negatives.
   **The model detects Ser/Thr-rich fungal cell-surface glycoproteins, of which adhesins are a subset** (§4.3).
3. **The classifier is a tandem-repeat detector, and the limit is the label, not the model (§4.7).**
   Every adhesin it misses has zero tandem repeats (median repeat coverage 1.000 for found vs
   0.000 for missed, p=0.0037). "Adhesin" spans avidity-mediated binding (many weak sites on
   repeats — blatant in sequence, PR-AUC 0.98, transfers across clades) and affinity-mediated
   binding (one folded interface — almost no sequence signal, scores ~0.000). **Do not buy a
   bigger model**: ESM C + Pfam is no better than ESM C alone. Split the label by mechanism,
   predict interpretable properties, and use structure only for the globular class (§8).
4. **Adhesin families are clade-specific, so "one fungal adhesin model" is the wrong target (§4.5).**
   The FLO/ALS/Hyr1 families behind the current model are Saccharomycotina-specific (Flocculin:
   330 proteins there, **0** in every other subphylum). Basidiomycota use CPL1-like (1,640
   proteins, 0 in any Ascomycota); chytrids use VWD and CBM18. A model trained on one clade and
   tested on another collapses to ROC 0.60 while a size-matched within-clade model reaches 0.94.
   Plan for a general stage 1 and **clade-specific stage 2** models.
5. **Stage 2 is learnable where the mechanism is visible (§4.4).** 75 curated adhesins vs 63
   curated non-adhesive surface proteins separate at PR-AUC 0.94 across held-out genomes with
   ESM C 300M; architecture alone (GPI + signal peptide + length) reaches only 0.70, so the
   signal is real and the PLM earns its place. Read together with item 3, this holds *within
   the avidity class*. **Fine-tuning is not warranted**: frozen embeddings plus explicit
   properties are the right next model, and fine-tuning is a decision gate for later (§8.3).
6. **Several code bugs affect reproducibility.** Padding tokens are included in the mean-pool,
   so results depend on batch composition. Failed batches silently shift labels. A mid-network
   layer is hard-coded for the 12-layer model. Sequences are truncated at 1,022 aa, which drops
   the C-terminus of 57% of positives. Only positive calls are written, so thresholds cannot be
   revisited. Details in §3.
7. **Compute:** embed every Fungi_5k proteome **once** with the chosen PLM on GPUs, store the
   embeddings in S3, and then do all classifier training, CV and re-scoring on CPU. NRP Nautilus
   is the right place for the bulk embedding run (sharded GPU Jobs driven by Nextflow, following the
   `nf_funannotate1` k8s pattern). HPCC is the right place for pilots, the `function.duckdb`
   joins, and classifier work (§9).

## 2. What exists (inventory)

| Piece | State |
|---|---|
| `src/adhesion_predict/` | ESM-2 (8M or 35M) mean-pooled embeddings → sklearn `LogisticRegression(C=1)`, no scaling. CLIs `adhesion_train/predict/evaluate`. |
| `models/*.pkl` | Two pickled LR models (320-d t6_8M; 480-d t12_35M). No metadata on layer, pooling or truncation. |
| `data/positive` | 581 sequences: `FLO11_Scer` 333 (**178 unique**), `ALS1_homologs` 128, `FLO_homologs` 118, plus FLO11 (S288C) and one *C. albicans* protein. |
| `data/negative` | 1,500 sequences: 500 random proteins each from FungiDB-68 *C. albicans*, *C. neoformans* JEC21 and *S. pombe*, after removing FASTA36 hits (E<1e-5) to positives. |
| `stage/` | FASTA36/phmmer homolog harvesting at ≥90% id (FLO11 in *S. cerevisiae* strains) and negative sampling. `--seed` default is `None`, so the negative sample is not reproducible. |
| `analysis/kingdom_survey` | 5,802 species; 749,697 positive calls (t12_35M, p>0.5 only). Clade stats. |
| `analysis/adhesion_properties` | Calls vs. size-matched background: 70.7% of calls have a signal peptide (vs. 6.3%); 40% have any Pfam; top background-absent Pfams are Flocculin, Hyr1, Msg2_C, Cadherin_4; AA1 laccases are enriched. |
| `analysis/embedding_clustering` | ESM2-classifier space vs. ESM C 300M space for the 750k calls; ESM C resolves 512 clusters. An ESM C venv (`.venv_esmc`) already exists on HPCC. |

Positive-set composition, measured here (MMseqs2 `easy-cluster --min-seq-id 0.3 -c 0.5`):
**581 positives collapse to 53 homology clusters.** The largest cluster (FLO11 across
*S. cerevisiae* strains) holds 284 sequences, about half the positive set. The effective
number of independent positives is ≈50. No cluster mixes positives and negatives.

Length: positives have median length 715–1,414 aa (by source file), negatives 370–455 aa.
57% of positives (330/581) are longer than the 1,022-residue truncation.

## 3. Code review findings

| # | Severity | Location | Finding | Fix |
|---|---|---|---|---|
| C1 | High (reproducibility) | `embeddings.py:116-117` | `token_representations.mean(dim=1)` averages over BOS, EOS **and padding** tokens, so a protein's embedding depends on the longest sequence in its batch (i.e. on file order and batch size, which differ between CPU=4 and GPU=auto). Measured: 7/2,081 training sequences change call with the 8M model; probability swings up to 0.57. | Mask to residue tokens (1…L) before averaging. Sort by length when batching. |
| C2 | High (silent corruption) | `embeddings.py:123-125`, `train.py:63`, `evaluate.py:44` | A failed batch is logged and skipped. The labels are then recovered with `sequences[:len(embeddings)]`, which misaligns every label after the failed batch. | Return ids; join labels by id; fail loudly or retry on CPU. |
| C3 | Medium | `embeddings.py:113-115` | `repr_layers=[6]` is hard-coded. For `esm2_t12_35M` (the model used for the Fungi_5k screen) that is the middle layer, not the final one. This is valid but undocumented, and it is not what the "ESM-2 embedding" label implies. | Use `model.num_layers`, or make the layer an explicit, saved hyperparameter. |
| C4 | Medium | `embeddings.py:106` | Truncation at 1,022 aa keeps only the N-terminal adhesion domain plus part of the repeats. The C-terminal GPI signal and most of the tandem repeat region are never seen for 57% of positives. | Sliding windows (e.g. 1,022 aa, stride 512) with mean or max across windows. Add explicit C-terminal features (GPI). |
| C5 | Medium | `model.py`, `models/*.pkl` | Bare pickled sklearn object: no record of ESM model, layer, pooling, truncation, sklearn version or training-data hash. No feature scaling. `C=1` is not tuned. | Save a model bundle (joblib + JSON model card). Use a `Pipeline(StandardScaler, LR)` with nested-CV-tuned `C`. |
| C6 | Medium | `scripts/predict.py:88-99` | Only positive calls are written by default, so the Fungi_5k run cannot be re-thresholded (acknowledged in the survey report). IDs containing commas are not CSV-escaped. | Always write all scores (parquet); filter downstream. Use the `csv` module. |
| C7 | Medium | `scripts/train.py` / `model.py` | Reported "test"/CV accuracy uses random splits of a heavily redundant set. `evaluate.py` defaults to the training directories, so it reports training-set performance. | See §6 validation framework. |
| C8 | Low | `features.py:17` | `seq.replace("J","L")` result discarded (no-op); module unused by the pipeline. | Remove or use as a baseline feature block. |
| C9 | Low | `io.py:112-121,148-157` | File handles never closed; `find_fasta_files` also globs upper-case extensions (duplicates on case-insensitive filesystems). | `with` blocks; de-duplicate paths. |
| C10 | Low | `README.md` | `scripts/training.py` does not exist (`adhesion_train`). FungiDB raw-file URLs now return **HTTP 401** (tested 2026-09-27), so the README and `stage/negtraining/download.sh` cannot be run as written. | Update docs; fetch via NCBI Datasets/UniProt or authenticated FungiDB. |
| C11 | Low | `stage/negtraining/sample_ids_from_fasta.py` | Default `--seed None`; the seed used for the committed negatives is unknown. | Record the seed; commit ID lists. |
| C12 | Info | `get_optimal_batch_size` | Batch sizing ignores sequence length (the real memory driver). | Token-budget batching after length sort. |

## 4. Experiments run for this review

§4.1–4.3 are on CPU (M-series Mac) with ESM-2 t6_8M; §4.4–4.7 use ESM C 300M embeddings
computed on HPCC GPUs. Scripts: `analysis/model_review/`.
Classifier for all rows: `StandardScaler → LogisticRegression(C=1)`. Metrics are pooled
out-of-fold predictions. "fpr" and "recall" are at p>0.5.

### 4.1 Cross-validation schemes × feature sets

| features | CV scheme | ROC-AUC | PR-AUC | recall | FPR |
|---|---|---|---|---|---|
| length only | random 5-fold (current practice) | 0.841 | 0.709 | 0.596 | 0.085 |
| length only | homology-grouped 5-fold | 0.796 | 0.599 | 0.294 | 0.082 |
| length only | leave-family-out | 0.820 | 0.659 | 0.428 | 0.075 |
| AA composition + length | random 5-fold | 0.999 | 0.998 | 0.981 | 0.006 |
| AA composition + length | homology-grouped 5-fold | 0.997 | 0.994 | 0.948 | 0.007 |
| AA composition + length | leave-family-out | 0.995 | 0.990 | 0.945 | 0.014 |
| ESM-2 8M, current pooling | random 5-fold | 1.000 | 0.999 | 0.993 | 0.005 |
| ESM-2 8M, current pooling | homology-grouped 5-fold | 0.999 | 0.998 | 0.985 | 0.005 |
| ESM-2 8M, current pooling | leave-family-out | 0.999 | 0.998 | 0.981 | 0.009 |
| ESM-2 8M, masked pooling | random 5-fold | 1.000 | 1.000 | 0.995 | 0.004 |
| ESM-2 8M, masked pooling | homology-grouped 5-fold | 0.999 | 0.996 | 0.983 | 0.005 |
| ESM-2 8M, masked pooling | leave-family-out | 0.999 | 0.998 | 0.990 | 0.013 |

*Homology-grouped* = `StratifiedGroupKFold` on MMseqs2 30%-identity clusters.
*Leave-family-out* = three folds, each holding out one positive family (ALS, FLO, FLO11)
together with one negative species' 500 proteins.

**Interpretation.** Homology leakage is real in the data (53 clusters), but it does not
inflate the metrics, because the classes are separable by composition alone. Even
leave-family-out, the strictest split possible with these data, stays at 0.995–0.999:
ALS and FLO share the same "secreted, Ser/Thr-rich, long" architecture. A 0.5–1% FPR looks
small but matters a lot at proteome scale. With ~6,000 proteins and a true prevalence
of ~0.3%, 0.5% FPR is ~30 false calls per genome, about the same as the number of true adhesins.
And FPR on *random* negatives underestimates FPR on the proteins that actually resemble adhesins.

### 4.2 Existing Fungi_5k calls for reference species (t12_35M model, p>0.5)

| species | proteins | calls | fraction |
|---|---|---|---|
| *Schizosaccharomyces pombe* (training-negative source) | 4,892 | 38 | 0.77% |
| *Cryptococcus neoformans* (training-negative source) | 5,987 | 60 | 1.00% |
| *Candida albicans* (training-negative source) | 5,700 | 84 | 1.47% |
| *Candidozyma auris* | 5,297 | 40 | 0.75% |
| *Nakaseomyces glabratus* | 5,092 | 106 | 2.08% |
| *Batrachochytrium dendrobatidis* | 7,021 | 97 | 1.38% |
| *Aspergillus fumigatus* | 9,161 | 114 | 1.24% |
| *Neurospora crassa* | 8,378 | 136 | 1.62% |

For scale: *C. albicans* has 8 Als and ~12 Hyr/Iff proteins plus Hwp1, Eap1, Rbt1 and a few
others (de Groot 2013; Hoyer & Cota 2016; Smoak 2023). *N. glabratus* has an unusually large
Epa/Pwp/Awp repertoire (~60–70). *S. pombe* has a handful (Linder & Gustafsson 2008).
Calls in *S. pombe* and *C. neoformans* at 38–60 per genome suggest the FP rate on real
proteomes is well above what CV reports.

### 4.3 Genome-wide check on *S. cerevisiae* S288C (shipped 8M model)

SGD `orf_trans_all` (6,722 ORFs, including dubious ones), embedded exactly as `predict.py` does (batch 4, file order, current pooling),
scored with `models/adhesion_model_esm2_t6_8M_UR50D.pkl`. For comparison: an LR retrained with masked pooling ("retrained"), and
the AA-composition + length LR ("aacomp"). Known adhesins (9): FLO1/5/9/10/11, AGA1, AGA2, SAG1, FIG2. Hard negatives (20): §7.2 list.

| model | calls | known adhesins called | hard negatives called | median length of calls |
|---|---|---|---|---|
| shipped (8M, current pooling) | 77 (1.15%) | 8/9 (misses 87-aa AGA2) | **12/20** | 341 aa (proteome 357) |
| retrained, masked pooling | 64 (0.95%) | 9/9 | 8/20 | 368 aa |
| AA composition + length | 124 (1.84%) | 9/9 | 8/20 | 153 aa |

**What the shipped model calls in S288C.** Of the 77 calls, only ≈9 are adhesins or FLO pseudogenes
(AGA1, SAG1, FIG2, FLO1/5/9/10/11, YHR213W), so **precision ≈ 12%** at recall 8/9. The rest are almost entirely
cell-wall and surface glycoproteins:
- GPI/cell-wall mannoproteins: SED1, DAN1, DAN4, TIR4, CCW12/14/22, FIT1–3, SPI1, CRH1, UTR2, PIR1, PIR5, CIS3, HSP150, HPF1, EGT2, SRL1, SUN4, SCW11, DSE2, CSS1/3, SVS1, KRE1/9, NCW2, CTS1
- mucin-like sensors: MSB2, HKR1, WSC2/3/4, MID2, SLG1, MTL1
- PRY1/2/3, and 11 dubious ORFs (including five identical subtelomeric 191-aa ORFs)

The call set is not length-driven (Spearman score-vs-log-length ≈ −0.07). **The model is a good detector of
fungal Ser/Thr-rich, O-mannosylated cell-surface glycoproteins, of which adhesins are a subset.** The kingdom
survey's `adhesion_fraction` is therefore better read as a "surface glycoprotein fraction" until the model is
retrained with hard negatives. Correcting the pooling (C1) changes individual calls but does not fix this.

### 4.4 Stage-2 proof of concept: adhesins vs other surface glycoproteins (2026-09-27)

The central question of this review is whether the *hard* discrimination is learnable at
all. Using the curated labels (§7, `data/curated/surface/surface.tsv`) and ESM C 300M
embeddings of those same proteins, `analysis/model_review/stage2_proof_of_concept.py`
answers it. Groups for CV are MMseqs2 30%-identity clusters, so no protein family spans a
fold; leave-genome-out additionally holds out a whole species.

**Hard version — 75 curated adhesins vs 63 curated non-adhesins** (Sed1, Gas1, Cwp1/2, Msb2,
Hkr1, Pir/Tir family and the rest of §7.2, i.e. exactly the proteins the shipped model
false-positives on):

| features | homology-grouped PR-AUC | leave-genome-out PR-AUC | P@k (leave-genome-out) |
|---|---|---|---|
| **ESM C 300M, full-length mean** | **0.981** | **0.939** | **0.91** |
| ESM C 300M, full + N-terminal | 0.980 | 0.955 | 0.91 |
| ESM C 300M, N-terminal only | 0.972 | 0.943 | 0.85 |
| AA composition + length | 0.939 | 0.907 | 0.87 |
| architecture only (GPI, signal peptide, length) | 0.702 | 0.743 | 0.80 |

Easier version, adding the 1,257 enzyme-annotated `non_adhesin_putative` proteins as
negatives: ESM C full-length reaches PR-AUC 0.976 (homology-grouped) / 0.951
(leave-genome-out) at 92–93% precision among its top-ranked 75.

**What this changes.**
1. **Stage 2 is learnable.** The shipped model's ~12% precision is a *label* problem, not a
   representation problem. With the right negatives, the same class of model separates
   adhesins from their look-alikes across held-out genomes.
2. **The PLM earns its place here, unlike in stage 1.** Architecture alone (GPI + signal
   peptide + length) gets PR-AUC 0.70, so the signal is not simply "GPI-anchored and long".
   ESM C beats amino-acid composition by a real margin (0.94 vs 0.91 leave-genome-out),
   the opposite of the saturated stage-1 benchmark in §4.1 where composition matched the PLM.
3. **75 positives are already enough to get started.** More curation will still help most,
   but this no longer blocks a usable v2 model.

**Caveats.** n=138 in the hard test, so the confidence intervals are wide and small
differences between the top rows are not meaningful. Negatives are mostly *S. cerevisiae*
and *C. albicans*. There are **no adhesin labels for Coccidioides, chytrids or most
Basidiomycota**, so "leave-genome-out" here means across six well-studied yeasts and
*A. fumigatus*, not across the fungal kingdom. And the labels themselves are a draft that
has not had expert review.

### 4.5 Adhesin families are clade-specific, and so is the model (2026-09-27)

Two independent lines of evidence say a single kingdom-wide adhesin model is the wrong
target. Script: `analysis/model_review/stage2_clade_transfer.py`.

**Pfam census.** UniProt protein counts per subphylum for the families that define known
adhesins. AMP-binding is a housekeeping control for how deeply each clade is sequenced;
read every row against it, not in absolute terms.

| family | Saccharomycotina | Pezizomycotina | Taphrinomycotina | Basidiomycota | Mucoromycota | Chytridiomycota |
|---|---|---|---|---|---|---|
| Flocculin (PF00624) | **330** | 0 | 0 | 0 | 0 | 0 |
| Candida_ALS_N (PF11766) | **313** | 0 | 0 | 6 | 0 | 0 |
| Hyr1 (PF11765) | **639** | 10 | 6 | 4 | 0 | 0 |
| **Flo11 (PF10182)** | **213** | **0** | 2 | 5 | 0 | 0 |
| GLEYA (PF10528) | 529 | 1,136 | 9 | 8 | 0 | 0 |
| Flocculin_t3 (PF13928) | **788** | **0** | 0 | 16 | 0 | 0 |
| DIPSY (PF11763) | 0 | 0 | **8** | 0 | 0 | 0 |
| Candida_ALS (PF05792) | 262 | 230 | 0 | 3 | 0 | 0 |
| PA14 (PF07691) | 289 | 3,043 | 2 | 675 | 7 | 7 |
| **CPL1-like (PF21671)** | **0** | **0** | **0** | **1,640** | **0** | **0** |
| CFEM (PF05730) | 369 | 11,592 | 25 | 1,615 | 1 | 5 |
| Hydrophobin_1 (PF01185) | 0 | 1,061 | 0 | **4,340** | 0 | 21 |
| Chitin_bind_1/CBM18 (PF00187) | 33 | 4,745 | 7 | 131 | 9 | **86** |
| VWD (PF00094) | 0 | 0 | 0 | 0 | 0 | **22** |
| *AMP-binding (control)* | *1,441* | *46,441* | *123* | *9,614* | *1,234* | *443* |

> **Correction (applied after first commit):** an earlier version of this table labelled
> PF10528 as "Flo11". PF10528 is **GLEYA**; the Flo11 domain is **PF10182**. The corrected
> rows are above, and the correction *strengthens* the conclusion: real Flo11 is absent from
> Pezizomycotina (0, not 1,136), as is Flocculin_t3. GLEYA is the family that actually spans
> Saccharomycotina and Pezizomycotina.

The FLO/ALS/Hyr1/Flo11 families the current model is trained on are
**Saccharomycotina-specific**: Flocculin, Flocculin_t3, Candida_ALS_N and Flo11 are all
absent from Pezizomycotina, and Hyr1 is ~1,000x rarer there after normalization. DIPSY is
Taphrinomycotina-only (the *S. pombe* pfl adhesins). Basidiomycota instead have CPL1-like (1,640 proteins,
**zero** in any Ascomycota) and hydrophobins; Chytridiomycota have VWD and CBM18. Three families do span clades at comparable normalized frequency — **CFEM** (all three major
clades), **PA14** (Pezizomycotina + Basidiomycota) and **GLEYA** (Saccharomycotina +
Pezizomycotina) — and those are the only plausible foundations for anything clade-transcending.

**Transfer test.** Using the curated labels and ESM C 300M embeddings:

| test | ROC-AUC | PR-AUC |
|---|---|---|
| leave-one-genome-out *within* Saccharomycotina (*C. albicans*, *N. glabratus*, S288C) | 0.948–1.000 | 0.830–1.000 |
| **leave-one-CLADE-out, Saccharomycotina held out** | **0.598** | **0.557** |
| size-matched control: same 13 pos + 5 neg, drawn from *within* Saccharomycotina | **0.940** (5–95%: 0.908–0.972) | 0.906 |

The cross-clade training set is tiny, so the control in row 3 is the one that matters: a
training set of *identical size* drawn from within the clade reaches 0.940, and **0 of 200
resamples scored as low as the cross-clade 0.598**. The collapse is about clade, not about
sample size.

**Caveats.** Pezizomycotina has only 2 positives and 5 negatives here, so its 0.900 is
meaningless. Taphrinomycotina has no negatives at all, so it cannot be scored and its 11
positives enter the cross-clade training set unbalanced — "clade" and "negative composition"
are therefore not fully separable in test B. The control in C is what carries the argument.
Basidiomycota and Chytridiomycota have **no labels at all**, so neither appears here.

**What this changes.** The plan in §8 should be **clade-aware**, not one kingdom-wide stage-2
model:
1. **Saccharomycotina** — labels exist, the model works (PR-AUC 0.94–0.98). Ship this first,
   scoped honestly to the clade it was trained on.
2. **Pezizomycotina** — Flo11, PA14, CFEM and CBM18 dominate; needs its own positives. The
   filamentous adhesins already seeded (Mad1/2, MPG1, CspA) are the starting point.
3. **Basidiomycota** — needs a CPL1-like-centred label set; see
   `analysis/curation/BASIDIOMYCETE_NOTES.md`. Expect near-zero recall from any
   Saccharomycotina-trained model, since the composition is inverted (Cfl1 family: 199–409 aa,
   5–10% Cys, secreted, not GPI; vs Als: 1,155–1,260 aa, ~1% Cys, GPI-anchored).
4. **Chytridiomycota** — no labels; VWD and CBM18 are leads only (`analysis/chytrid_batrach/`).

Stage 1 (surface glycoprotein) may still generalize, since signal peptides, GPI anchors and
Ser/Thr enrichment are universal. That split — **general stage 1, clade-specific stage 2** —
is the design the evidence supports.

### 4.6 Does the language model find anything HMMs cannot? Yes (2026-09-27)

The question that decides whether a PLM earns its place at all: if a Pfam rule matches it,
the cheaper and more interpretable method should win. Script:
`analysis/model_review/esm_vs_hmm.py`. Homology-grouped CV, 75 curated adhesins vs 63
curated non-adhesins.

| features | ROC-AUC | PR-AUC | ROC on domain-blind subset | PR on domain-blind subset |
|---|---|---|---|---|
| HMM rule (carries a known adhesin-family Pfam) | 0.857 | 0.829 | 0.410 | 0.115 |
| all Pfam domains (one-hot, the generous baseline) | 0.928 | 0.946 | 0.763 | 0.246 |
| **ESM C 300M** | **0.977** | **0.983** | **0.879** | **0.716** |
| ESM C + Pfam | 0.979 | 0.986 | 0.877 | 0.723 |

The "domain-blind subset" is the 66 proteins carrying **no** adhesin-family domain (9
adhesins, 57 non-adhesins) — where HMMs have nothing to go on by construction.

Three results:
1. **A pure HMM rule recovers only 73% of curated adhesins** (55/75) at 93% precision. It
   misses 20 outright, including Hwp1, Hwp2, Eap1, Scf1, Aga1, Aga2, and the *S. pombe*
   gsf2/pfl proteins. High precision, and a hard ceiling on recall.
2. **On exactly those proteins, ESM is ~3x better than the best domain baseline**
   (PR-AUC 0.716 vs 0.246). Of the 9 domain-blind adhesins, ESM scores 7 above 0.5, median
   0.693, against a median of 0.001 for curated non-adhesins.
3. **Pfam adds essentially nothing on top of ESM** (0.986 vs 0.983). The embedding already
   encodes what the domain annotation encodes, plus whatever lets it rank the domain-blind
   cases.

**So yes — this is moving in the direction of cataloguing beyond HMMs.** The honest framing
is that the PLM is not replacing domain annotation, it is extending it into the
low-complexity, repeat-rich, poorly-annotated fraction of the surface proteome where fungal
adhesins disproportionately live. That is also precisely the fraction where we have the
least ground truth, so the claim needs independent validation before it is leaned on.

**Caveat:** only 9 adhesins in the domain-blind subset. The direction is clear, the effect
size is not well estimated.

### 4.7 What the model actually detects: tandem repeats (2026-09-27)

The Eurotiomycetes labels (§7.3) exposed a consistent failure pattern: BAD1, SOWgp and CspA
are found, while rodA, CalA and gp43 are missed entirely. Testing whether that tracks repeat
content rather than clade or training composition — out-of-fold scores under homology-grouped
CV, correlated against sequence properties of the 95 curated adhesins. Script:
`analysis/model_review/what_the_model_detects.py`.

| property | Spearman rho | p |
|---|---|---|
| **8-mer repeat coverage** | **+0.302** | **0.003** |
| low-complexity fraction (top-3 aa) | +0.140 | 0.18 |
| length | +0.085 | 0.41 |

| | n | median repeat coverage |
|---|---|---|
| **found** (p>0.5) | 85 | **1.000** |
| **missed** (p<0.5) | 10 | **0.000** |

Mann-Whitney, found > missed: repeat coverage p=0.0037, low-complexity p=0.047, length p=0.044.

**Every adhesin the model misses has zero tandem repeats**: rodA (159 aa hydrophobin), CalA
(177 aa invasin), gp43 (416 aa moonlighting glucanase), PGA1 (132 aa), SAG1, and two PA14
proteins.

So the classifier is, functionally, a **tandem-repeat surface protein detector**. That is a
sharper and less flattering description than "surface glycoprotein detector" (§4.3), and it
supersedes it.

**This locates the limit in the label, not the model.** "Adhesin" spans at least two
physically distinct mechanisms:

| mechanism | how it binds | sequence signal | model performance |
|---|---|---|---|
| **avidity-mediated** | many weak sites on tandem repeats; strength from multivalency | blatant (repeat periodicity) | PR-AUC 0.98, transfers across clades — SOWgp and BAD1 score ~1.0 despite being Onygenales |
| **affinity-mediated** | one folded domain, one high-affinity interface | very weak — "this fold has a binding site" is not a sequence feature | ~0.000 |

Two lines of evidence say model capacity is not the constraint: ESM C + Pfam is no better than
ESM C alone (0.986 vs 0.983, §4.6), and within the repeat-rich class performance is already
0.98. Scaling parameters will not manufacture a signal that is not in the sequence.

## 5. Why the model over-calls: diagnosis

1. **Negatives are random proteins**, so the easiest decision boundary is secreted +
   low-complexity Ser/Thr-rich + long. Every fungal GPI cell-wall protein, mucin-like sensor
   (Msb2, Hkr1), secreted glycosidase with an S/T-rich linker, and heavily O-glycosylated
   extracellular enzyme lands on the positive side. This is consistent with the adhesion_properties
   report (70% signal peptide; AA1 laccases and CBMs enriched).
2. **Positives are one architecture from one clade**: FLO/ALS/FLO11 from Saccharomycotina,
   with half the set being near-identical FLO11 alleles. Adhesins with other architectures
   (Epa/PA14 lectins, Hyr/Iff, Msg, Mad1, CotH, Bad1 tandem repeats, Scf1, hydrophobin-like
   and fibrillar proteins) are absent. The model cannot learn "adhesion", only "looks like FLO/ALS".
3. **The label means different things.** "Adhesin" in the literature ranges from demonstrated
   ligand binding to GPI-CWP family membership. The training positives were harvested by
   sequence similarity (E<1e-10 hits to FLO11 can be driven by Ser/Thr-rich repeats rather than
   domain homology), so some positives are likely not adhesins either. Examples to check:
   *Fusarium*, *Colletotrichum* and *Zymoseptoria* members of `FLO_homologs`.
4. **Threshold 0.5 on a model trained at ~1:2.6 class ratio** implies a prior far above the true
   genomic prevalence (~0.2–0.5%). Probabilities are not calibrated for the proteome setting.

## 6. Proposed validation framework (tiers)

Every candidate model is scored on all tiers. Report per-tier numbers in a model card.

| Tier | Question | Data | Metric |
|---|---|---|---|
| T0 sanity | Is the pipeline deterministic? | fixed set, shuffled batch orders | identical scores |
| T1 in-distribution | Does it separate curated positives from random negatives without leakage? | curated set, MMseqs2 30%/50%-cov clusters (or GraphPart) | homology-grouped CV ROC/PR-AUC; nested CV for hyperparameters |
| T2 out-of-family | Does it recognize adhesin *families it never saw*? | leave-one-family-out over ≥8 families (§7.1) | recall per held-out family at a fixed FPR |
| T3 out-of-clade | Does it transfer across lineages? | leave-one-subphylum/class-out (Saccharomycotina, Pezizomycotina, Basidiomycota, Mucoromycota, Chytridiomycota) | recall/precision per clade |
| T4 hard negatives | Does it reject look-alikes? | curated non-adhesive GPI-CWPs, mucin-like sensors, S/T-rich secreted enzymes, flagged laccases (§7.2) | FPR on hard set; PR-AUC positives vs. hard |
| T5 proteome precision | How many calls per genome, and are they right? | complete proteomes with curated inventories: S288C, *C. albicans* SC5314, *N. glabratus* CBS138, *C. auris* B8441, *S. pombe*, *A. fumigatus* Af293 | precision/recall vs. inventory; calls per genome; **threshold chosen here**, not at 0.5 |
| T6 calibration and consistency | Are scores meaningful across the kingdom? | Fungi_5k | reliability curve on T5; Microsporidia and other expected-low clades as negative controls; call-count vs. proteome size |

Implementation notes:
- Clustering: MMseqs2 `easy-cluster` (fast, installed on HPCC as a module or via conda) or
  GraphPart (Teufel et al. 2023), which guarantees a maximum identity between partitions.
- Use `StratifiedGroupKFold` on cluster IDs everywhere. Tune only in the inner loop.
- Because unlabeled proteomes contain unknown adhesins, treat the proteome-scale data as
  **positive-unlabeled (PU)**. The class prior can be estimated and used to calibrate
  (Bekker & Davis 2020).

## 7. Is more literature curation warranted? Yes. It is the highest-value next step.

### 7.1 Positive families to curate (with anchor references)

| Family / protein | Lineage | Architecture | Anchor refs |
|---|---|---|---|
| Flo1/5/9/10 (Flocculin/PA14), Flo11 | *Saccharomyces*, Saccharomycotina | N-terminal lectin (PA14) or Flo11 domain; S/T repeats; GPI | Verstrepen & Klis 2006 |
| Als1–9 | *Candida* spp. | Als N-terminal Ig-like + T domain; amyloid; 36-aa repeats; GPI | Hoyer & Cota 2016 |
| Hyr/Iff | *Candida*, *C. auris*, *N. glabratus* | Hyr1 N-terminal domain; S/T; GPI; parallel expansions | Smoak et al. 2023; Santana et al. 2023 (Iff4109) |
| Epa (PA14), Pwp, Awp | *Nakaseomyces glabratus* | PA14 lectin; subtelomeric | Cormack et al. 1999; de Groot 2013 |
| Scf1 | *C. auris* | novel; cationic-residue adhesion | Santana et al. 2023 |
| Als-like / Scf / Iff in *C. auris* | *C. auris* | redundancy in aggregation | Wang et al. 2024 |
| Map4 / Mam3 / Pfl family | *S. pombe*, Taphrinomycotina | fission-yeast adhesin family | Linder & Gustafsson 2008 |
| Msg (major surface glycoprotein) | *Pneumocystis* | Msg domains; large subtelomeric family | de Groot 2013 (review) |
| Mad1 / Mad2 | *Metarhizium* (Sordariomycetes) | insect vs. plant adhesion | Wang & St Leger 2007 |
| CotH3 (binds GRP78) | Mucorales (*Rhizopus*) | spore-coat homolog; host-receptor ligand | Gebremariam et al. 2019 |
| Bad1, SOWgp, Pb gp43, Aspergillus CspA and others | dimorphic Onygenales, Aspergillus | mixed: tandem repeats, GPI, secreted | de Groot 2013; Chaudhuri 2011 (FungalRV positive set) |
| Diverse GPI-anchored adhesins (structural classes) | Ascomycota | PA14, Flo11, Als, Hyr, CFEM, etc. | Essen et al. 2020 |

Existing labelled sets to harvest and de-duplicate: **FungalRV** (Chaudhuri et al. 2011; positive set from
10 human pathogens), **FaaPred** (Ramana & Gupta 2010) and Nath 2019. All three used compositional SVMs and
reported ≥0.87 MCC under random CV, which is the same saturated-benchmark problem described here.

Curation deliverable: `data/curated/adhesins.tsv` with UniProt/NCBI accession, species, family,
**evidence level** (E1 = experimental adhesion/binding phenotype; E2 = member of an experimentally
characterized family with a conserved adhesion domain; E3 = computational only), and PMID.
Train on E1+E2, report E3 separately.

### 7.2 Hard negatives to curate

Experimentally characterized, non-adhesive proteins that share the adhesin architecture:
- S. cerevisiae GPI cell-wall proteins without adhesion function: Gas1, Crh1, Cwp1/2, Ccw12, Ccw22,
  Sed1, Tir1/3/4, Dan1, Ecm33, Yps1, Ncw2; Pir family (Pir1, Cis3, Hsp150).
- Mucin-like signaling proteins: Msb2, Hkr1 (and *C. albicans* Msb2).
- Secreted, S/T-rich-linker enzymes: glucanases (Exg1), chitinases, and laccases/AA1
  (currently enriched among the calls).
- Orthologs of the above across the reference proteomes (via OrthoFinder or MMseqs2 reciprocal hits).

### 7.3 Onygenales and Eurotiales (curated 2026-09-27)

`data/curated/adhesins/eurotiomycetes_seeds.tsv` — 21 rows, 10 Onygenales and 11 Eurotiales,
12 adhesins (4 at E1) and 9 hard negatives, all with PMIDs and resolved accessions. Trainable
Eurotiomycetes labels went 7 → 22; Onygenales went from zero to seven.

Adhesins: SOWgp (three size alleles differing in tandem-repeat number; rSOWgp binds
laminin/fibronectin/collagen IV and deletion cuts virulence), BAD1, CalA, CspA, rodA/rodB,
gp43, Ag2/PRA (E3 — surface antigen, adhesion not directly shown).

Hard negatives are the scarcer and more valuable half, and were chosen to be genuinely hard:
Mp1p (abundant *Talaromyces* surface mannoprotein whose demonstrated mechanism is
arachidonic-acid sequestration), cfmA/cfmC (CFEM GPI proteins whose triple deletion affects
only cell-wall stability — and CFEM is the clade-spanning family), CTS1 (immunodominant
*Coccidioides* antigen, but an enzyme), gel1, ecm33, abr2 laccase, Cbp1.

**Result: the added labels did not improve prediction** (ROC 0.732 → 0.670 leave-one-out),
but they made the clades measurable for the first time, and the per-protein pattern is what
led to §4.7. Three pipeline bugs surfaced: `surface.tsv` was silently dropping curated
proteins from non-reference isolates; moonlighting proteins (Hsp60) were entering training as
positives despite having no secretion signal; and E3 domain-only guesses were being trained
on. All three are fixed. Full detail in `docs/model-review/STATUS.md` §7.

## 8. Models: what to use, and what not to buy

**Conclusion from §4.4–4.7: do not buy a bigger model, and do not reach for structure first.
Split the label by adhesion mechanism.** The evidence is that ESM already saturates the
mechanism it can see, and is structurally blind to the others.

### 8.0 Predict mechanisms and properties, not "adhesin"

"Adhesin" is an outcome (the cell sticks to something), not a molecular feature, and the
classifier has quietly reduced it to "has tandem repeats" (§4.7). The design that follows
from the evidence is to predict **interpretable properties**, each mapping to a mechanism,
and treat "adhesin" as an inference over that profile:

| property | mechanism it indicates | predictable from sequence? | status here |
|---|---|---|---|
| tandem repeat content / periodicity | avidity-mediated adhesion | yes, and already discriminative (p=0.003) | **works now** |
| GPI anchor + cell-wall retention | surface display | yes (NetGPI/PredGPI) | available, not yet used directly |
| beta-aggregation / amyloid propensity | Als-type amyloid-mediated adhesion, a documented mechanism | yes (TANGO/Waltz-style) | **unexploited — the clearest gap** |
| hydrophobin Cys pattern | surface-hydrophobicity attachment (rodA, Rep1 repellents) | yes, trivially (PF01185) | a solved problem currently misfiled as a model failure |
| surface charge distribution | cation-mediated adhesion (*C. auris* Scf1) | partly | not attempted |
| folded binding interface | affinity-mediated receptor binding (CalA) | **no** | needs structure (§8.1) |

The gain is honesty as much as accuracy: a per-mechanism output states *why* a protein is
called, instead of one opaque score that silently means "has repeats".

### 8.1 Where structure genuinely helps — and why it is complementary, not a replacement

Structure is the right tool for the **affinity-mediated** class specifically: exposed binding
interfaces, fold recognition, electrostatic patches. It is not a general fix.

The useful asymmetry: **AlphaFold/ESMFold are weakest on the repeat-rich, disordered regions
that define FLO/ALS** — exactly where sequence already works — and strongest on compact
globular domains, exactly where sequence fails. The two approaches are complementary by
construction, so the question is not sequence *or* structure but which class each is applied to.

Not a starting point, though: folding thousands of candidates is expensive, and without the
mechanism split there is no principled way to choose which proteins warrant it. Do the split
first, then fold the globular candidates.

### 8.2 Orthogonal signal both HMMs and PLMs ignore

Adhesins sit disproportionately in subtelomeric, repeat-rich, copy-number-variable regions.
This is the same two-speed genome architecture documented for *Batrachochytrium* virulence
factors (Wacker et al. 2023; see `analysis/chytrid_batrach/`). **Genomic context — subtelomeric
position, TE proximity, copy-number variation between isolates — is information no
sequence-only model uses.** It is cheap to compute from existing assemblies and is worth
testing as a feature.

### 8.3 Frozen embeddings, and the fine-tuning gate (unchanged)

**Two-stage framing (suggested by the S288C result).** The current model is already good at stage 1.
- **Stage 1: cell-surface glycoprotein** (signal peptide + GPI or S/T-rich mucin-like, O-mannosylated).
  This is well defined, easy to label from SignalP/NetGPI/composition, and is roughly what `adhesion_fraction`
  measures today. Report it as its own, honestly named quantity per genome.
- **Stage 2: adhesin vs. other surface protein**, trained *only within stage-1 proteins*, with curated adhesins
  as positives and non-adhesive cell-wall proteins as negatives. This is where N-terminal-domain embeddings,
  adhesion-domain hits and family-held-out validation (T2, T4) matter.
Scores from both stages are kept, so kingdom-wide analyses can use either level.

1. **Embeddings.** Switch to ESM C 300M (open licence; about ESM-2 650M quality at ~½ the cost)
   or ESM C 600M (non-commercial licence; fine for academic use). The embedding_clustering work
   already has an ESM C 300M environment on HPCC. Vieira, Handojo & Wilke 2025 found medium-sized
   models (ESM-2 150M, ESM C 600M) with **mean pooling** perform about as well as the largest
   models for transfer learning on realistic dataset sizes. Pool over residue tokens only.
   Use sliding windows for proteins >1,022 aa (ESM-2) or ESM C's longer context.
   Also keep a separate **N-terminal-domain embedding** (first ~300 aa after the signal peptide),
   because adhesion specificity sits in the N-terminal domain, while the repeats and GPI dominate a full-length mean.
2. **Architecture features.** Signal peptide (SignalP 6.0; Teufel 2022). GPI-anchor
   (NetGPI 1.1; Gíslason et al. 2021, or PredGPI, which you already have in `~/projects/predgpi`).
   TM helices. S/T/P fraction. Tandem-repeat content (e.g. T-REKS/XSTREAM). Intrinsic disorder.
   Presence of known adhesion domains (Pfam/InterPro hits). Most of these are already computed in
   the Fungi_5k `function.duckdb`.
3. **Classifier.** Regularized logistic regression or gradient boosting on the concatenated blocks,
   with class weights. Probability calibration (isotonic or Platt) on T5 genomes. Threshold chosen for a target
   precision per genome. Report the PU-estimated class prior.
4. **Fine-tuning gate.** Consider parameter-efficient fine-tuning (LoRA) of ESM-2 150M/650M or ESM C 300M only if:
   (a) curated E1+E2 positives exceed ~300 homology clusters; (b) frozen-embedding models plateau on T2/T4;
   and (c) the gain is shown on held-out families. Schmirler, Heinzinger & Rost 2024 show fine-tuning
   usually helps, but gains are task-dependent and largest with sufficient data. With ~50
   independent positive clusters today, fine-tuning would mostly memorize FLO/ALS.
5. **Things *not* to do:** do not retrain on the Fungi_5k calls (circular). Do not cluster calls and
   treat clusters as new positives without independent evidence.

## 9. Compute plan (HPCC vs. NRP Nautilus)

### 9.1 Sizing — **measured**, not estimated (pilot run 2026-09-27, HPCC `short_gpu` gpu09, RTX 6000 Ada)

`analysis/model_review/pilot/` embedded all six reference proteomes with three models
(length-sorted token-budget batching, bf16, residue-only mean pooling, sliding windows for
long proteins plus a separate N-terminal embedding). Raw numbers: `pilot/throughput.jsonl`.

| model | proteins/s | residues/s | peak GPU mem | dim | projected GPU-hours for Fungi_5k (~58M proteins) |
|---|---|---|---|---|---|
| **ESM C 300M** | **100** | 48,800 | 5.5 GB | 960 | **~162** |
| ESM-2 150M | 60 | 29,500 | 9.8 GB | 640 | ~267 |
| ESM-2 650M | 35 | 17,000 | 13.6 GB | 1280 | ~462 |

Mean GPU utilization during the run was **91%** (NRP requires >40%), so the batching
strategy is already efficient enough for the cluster's policy.

**ESM C 300M is the clear choice**: fastest, smallest memory footprint, open licence, and
the literature puts it near ESM-2 650M in representation quality (Vieira et al. 2025).
At 162 GPU-hours, the full kingdom re-screen is ~20 hours of wall time across 8 concurrent
NRP GPUs, or a few days on 2-3. This is a one-time cost: embeddings are then reused for
every retrain.

- Storage: 58M x 960-d float16 = **~110 GB** for the full-length pooled embeddings, doubled
  if the N-terminal embeddings are kept too (they are small and worth keeping).
- Everything downstream (training, CV, re-scoring 58M proteins with a linear model) is CPU
  work: minutes to hours on HPCC.
- For comparison, the original screen ran the shipped classifier on CPU across a SLURM array
  (`run_fungi5k*.sh`), which is why it took days.

### 9.2 Where to run what
| Work | Where | Why |
|---|---|---|
| Curation, training, CV, figures | laptop / HPCC CPU | small data |
| Pilot embedding (T5 reference genomes, ~60k proteins) and throughput benchmark | HPCC `exfab` partition (gpu12, 2× ada6000) | already configured; high PriorityTier; the `gpu` partition had 10 h fairshare waits |
| Bulk embedding of Fungi_5k | **NRP Nautilus**, namespace `ucr-stajichlab` | many GPUs; S3 next to compute; Nextflow k8s pattern already working for `nf_funannotate1` |
| Joins with `function.duckdb` (278 GB) | HPCC | the database lives on `/bigdata` |

### 9.3 NRP design (fits the NRP policies)
- **Nextflow on k8s**, reusing `nf_funannotate1`'s run-as-Job head pod, RBAC `nextflow-runner`, PVC and
  `site.config`. One process `EMBED(asmid)` per proteome shard with `accelerator 1` and a GPU image
  (PyTorch + `esm`/`fair-esm`). Inputs are pulled from and outputs pushed to `s3://stajichlab`
  (`https://s3-west.nrp-nautilus.io`). Use `-resume` for restarts.
- **GPU rules** (NRP docs): request `nvidia.com/gpu: 1` per pod (max 2 per pod). Select the card via node
  affinity on `nvidia.com/gpu.product`, and match the CUDA runtime via `nvidia.com/cuda.runtime.major`.
  A100/H100 need quota or `priorityClassName: opportunistic` (preemptible, so shards must be idempotent).
  Keep **GPU utilization >40%**, ideally ~100% (Grafana namespace GPU dashboard). This is why large
  length-sorted batches and pre-staged inputs matter: do not let a GPU idle on S3 I/O.
- **Pod rules:** limits within 20% of requests (equal if >100 pods). No `sleep infinity`. Jobs must exit
  on their own. There is no fair queue, so cap concurrency (`executor.queueSize` ≈ 8–16 GPU tasks) rather than
  submitting 5,800 pods at once. Set `backoffLimit` > 0 and `/dev/shm` as a memory `emptyDir` for DataLoader workers.
- **Shard size:** ~50 proteomes (~500k proteins) per task, so ~120 tasks of ~30–60 min each. That size
  tolerates preemption and keeps per-task overhead (model load ~seconds) negligible.
- Monitoring: https://grafana.nrp-nautilus.io (namespace pods and GPU dashboards).

## 10. Roadmap

Revised 2026-09-27 after §4.7. The change from the original plan: **P1 is no longer "curate
more adhesins" but "split the label by mechanism"**, because the Eurotiomycetes curation
showed more labels of the same architecture do not help.

| Phase | Deliverable | Compute |
|---|---|---|
| P0 fixes | C1–C6 fixed; deterministic embeddings; model card; tests | laptop — **done**, PR #18 |
| P1 **mechanism split** | annotate every curated adhesin with its mechanism class (avidity / affinity / hydrophobin / moonlighting); ship the avidity-class model, scoped and named honestly | laptop |
| P1b property predictors | beta-aggregation/amyloid propensity, GPI, repeat periodicity, hydrophobin Cys pattern, surface charge — as explicit features, not one opaque score | laptop / HPCC CPU |
| P2 benchmark | T0–T6 harness, reported **per mechanism class and per clade**, never as one pooled number | HPCC CPU |
| P3 structure, selectively | fold the globular (affinity-class) candidates only; test whether binding-interface features recover CalA-type invasins | HPCC / NRP GPU |
| P3b genomic context | subtelomeric position, TE proximity, copy-number variation as features (§8.2) | HPCC CPU |
| P4 kingdom re-screen | NRP Nextflow embedding of Fungi_5k to S3; re-score; re-run the surveys with all-score files | NRP GPU + HPCC CPU |
| P5 (gated) | fine-tuning only if a mechanism-specific model plateaus with adequate labels | NRP / HPCC GPU |

Tracking: GitHub issues on `stajichlab/adhesionPred` (§12).

## 11. Sources

### Literature (checked via PubMed or publisher pages on 2026-09-27)
- Chaudhuri R, et al. FungalRV: adhesin prediction and immunoinformatics portal for human fungal pathogens. *BMC Genomics* 2011;12:192. PMID 21496229. https://doi.org/10.1186/1471-2164-12-192
- Ramana J, Gupta D. FaaPred: a SVM-based prediction method for fungal adhesins and adhesin-like proteins. *PLoS One* 2010;5:e9695. PMID 20300572. https://doi.org/10.1371/journal.pone.0009695
- Nath A. Prediction and molecular insights into fungal adhesins and adhesin like proteins. *Comput Biol Chem* 2019;80:333-340. PMID 31078912. https://doi.org/10.1016/j.compbiolchem.2019.05.001
- de Groot PWJ, et al. Adhesins in human fungal pathogens: glue with plenty of stick. *Eukaryot Cell* 2013. PMID 23397570. https://doi.org/10.1128/EC.00364-12
- Verstrepen KJ, Klis FM. Flocculation, adhesion and biofilm formation in yeasts. *Mol Microbiol* 2006. PMID 16556216. https://doi.org/10.1111/j.1365-2958.2006.05072.x
- Hoyer LL, Cota E. *Candida albicans* Agglutinin-Like Sequence (Als) family vignettes. *Front Microbiol* 2016;7:280. PMID 27014205. https://doi.org/10.3389/fmicb.2016.00280
- Linder T, Gustafsson CM. Molecular phylogenetics of ascomycotal adhesins, a novel family of putative cell-surface adhesive proteins in fission yeasts. *Fungal Genet Biol* 2008 (epub 2007). PMID 17870620. https://doi.org/10.1016/j.fgb.2007.08.002
- Smoak RA, et al. Parallel expansion and divergence of an adhesin family in pathogenic yeasts. *Genetics* 2023. PMID 36794645. https://doi.org/10.1093/genetics/iyad024
- Cormack BP, Ghori N, Falkow S. An adhesin of the yeast pathogen *Candida glabrata* mediating adherence to human epithelial cells. *Science* 1999;285:578. PMID 10417386. https://doi.org/10.1126/science.285.5427.578
- Santana DJ, et al. A *Candida auris*-specific adhesin, Scf1, governs surface association, colonization, and virulence. *Science* 2023;381:1461-1467. PMID 37769084. https://doi.org/10.1126/science.adf8972
- Wang et al. Functional redundancy in *Candida auris* cell surface adhesins crucial for cell-cell interaction and aggregation. *Nat Commun* 2024. PMID 39455573. https://doi.org/10.1038/s41467-024-53588-5
- Wang C, St Leger RJ. The MAD1 adhesin of *Metarhizium anisopliae* links adhesion with blastospore production and virulence to insects, and the MAD2 adhesin enables attachment to plants. *Eukaryot Cell* 2007. PMID 17337634. https://doi.org/10.1128/EC.00409-06
- Gebremariam T, et al. Anti-CotH3 antibodies protect mice from mucormycosis by prevention of invasion and augmenting opsonophagocytosis. *Sci Adv* 2019. PMID 31206021. https://doi.org/10.1126/sciadv.aaw1327
- Essen LO, et al. Diversity of GPI-anchored fungal adhesins. *Biol Chem* 2020. PMID 33035180. https://doi.org/10.1515/hsz-2020-0199
- Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model (ESM-2). *Science* 2023. PMID 36927031. https://doi.org/10.1126/science.ade2574
- Schmirler R, Heinzinger M, Rost B. Fine-tuning protein language models boosts predictions across diverse tasks. *Nat Commun* 2024. PMID 39198457. https://doi.org/10.1038/s41467-024-51844-2
- Vieira LC, Handojo ML, Wilke CO. Medium-sized protein language models perform well at transfer learning on realistic datasets. *Sci Rep* 2025;15:21400. https://doi.org/10.1038/s41598-025-05674-x
- Teufel F, et al. GraphPart: homology partitioning for biological sequence analysis. *NAR Genom Bioinform* 2023. PMID 37850036. https://doi.org/10.1093/nargab/lqad088
- Teufel F, et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nat Biotechnol* 2022. PMID 34980915. https://doi.org/10.1038/s41587-021-01156-3
- Gíslason MH, Nielsen H, Almagro Armenteros JJ, Johansen AR. Prediction of GPI-anchored proteins with pointer neural networks (NetGPI). *Curr Res Biotechnol* 2021. https://doi.org/10.1016/j.crbiot.2021.01.001 ; server: https://services.healthtech.dtu.dk/services/NetGPI-1.1/
- Pierleoni A, Martelli PL, Casadio R. PredGPI: a GPI-anchor predictor. *BMC Bioinformatics* 2008;9:392 (local copy `~/projects/predgpi`).
- EvolutionaryScale. ESM Cambrian (ESM C 300M/600M/6B), 2024. https://www.evolutionaryscale.ai/blog/esm-cambrian ; licence: https://github.com/evolutionaryscale/esm/blob/main/LICENSE.md
- *Cited from memory, not re-verified this session:* Bekker J, Davis J. Learning from positive and unlabeled data: a survey. *Mach Learn* 2020;109:719-760.

### Infrastructure documentation
- NRP GPU pods: https://nrp.ai/documentation/userdocs/running/gpu-pods/
- NRP Jobs: https://nrp.ai/documentation/userdocs/running/jobs/
- NRP scheduling (Armada priority classes): https://nrp.ai/documentation/userdocs/running/scheduling/
- NRP monitoring: https://nrp.ai/documentation/userdocs/running/monitoring/
- NRP cluster policies: https://nrp.ai/documentation/userdocs/start/policies/
- Lab NRP configs: `~/projects/nrp-deploy` (namespace `ucr-stajichlab`, project `stajichlab-fungi-bfd`,
  `executor.queueSize = 25`, S3 endpoint `https://s3-west.nrp-nautilus.io`); Nextflow CPU run pattern in
  `~/projects/Bd_pangenome` + `nf_funannotate1` k8s base.
- HPCC: `analysis/embedding_clustering/run_extract_esm2.sbatch` (partition notes: `exfab` / gpu12 ada6000;
  P100 nodes incompatible with current PyTorch; `gpu` partition fairshare waits).
- Data: `/bigdata/stajichlab/shared/projects/Fungi_5k/{samples.csv,input/,functionalDB/function.duckdb}`.
- S288C proteome: http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans_all.fasta.gz
  (FungiDB raw files returned HTTP 401 on 2026-09-27).

### In-repo evidence used
- `analysis/kingdom_survey/REPORT.md`, `tables/species_adhesion_summary.csv`
- `analysis/adhesion_properties/REPORT.md` (domain/topology enrichment, AA1 laccases)
- `analysis/embedding_clustering/REPORT.md`
- `analysis/model_review/` (this review's scripts)

## 12. Tracking issues (filed 2026-09-27)

| # | Title |
|---|---|
| [#7](https://github.com/stajichlab/adhesionPred/issues/7) | Mean-pooling includes BOS/EOS/padding tokens (C1) |
| [#8](https://github.com/stajichlab/adhesionPred/issues/8) | Failed batches skipped silently → labels misaligned (C2) |
| [#9](https://github.com/stajichlab/adhesionPred/issues/9) | Model bundle + card; stop hard-coding layer 6 (C3, C5) |
| [#10](https://github.com/stajichlab/adhesionPred/issues/10) | Proteins >1022 aa: sliding windows + N-terminal embedding (C4) |
| [#11](https://github.com/stajichlab/adhesionPred/issues/11) | predict: write all scores; CSV escaping (C6) |
| [#12](https://github.com/stajichlab/adhesionPred/issues/12) | Validation framework tiers T0–T6 |
| [#13](https://github.com/stajichlab/adhesionPred/issues/13) | Curate expanded adhesin positives with evidence levels |
| [#14](https://github.com/stajichlab/adhesionPred/issues/14) | Curate hard negatives |
| [#15](https://github.com/stajichlab/adhesionPred/issues/15) | Model v2 (frozen PLM + architecture features; fine-tuning gated) |
| [#16](https://github.com/stajichlab/adhesionPred/issues/16) | NRP Nextflow GPU embedding of Fungi_5k → S3 |
| [#17](https://github.com/stajichlab/adhesionPred/issues/17) | Docs/data hygiene |

## Appendix: files from this review
- `analysis/model_review/run.sh`: end-to-end reproduction (downloads S288C from SGD, MMseqs2 clustering, embeddings, CV, proteome check)
- `analysis/model_review/01_embed_esm2_8M.py`: legacy vs. masked pooling embeddings
- `analysis/model_review/02_cv_and_proteome_eval.py`: CV schemes × feature sets; S288C adhesin and hard-negative panel
- `analysis/model_review/stage2_proof_of_concept.py`: stage-2 separability (§4.4)
- `analysis/model_review/stage2_clade_transfer.py`: clade-specificity and the size-matched control (§4.5)
- `analysis/model_review/esm_vs_hmm.py`: PLM vs domain-annotation baselines (§4.6)
- `analysis/model_review/what_the_model_detects.py`: repeat-content analysis (§4.7)
- `analysis/model_review/clade_status_report.py`: per-clade status tables (STATUS.md)
- `analysis/model_review/pilot/`: measured PLM throughput on HPCC GPUs (§9.1)
- `analysis/model_review/cv_results.tsv`, `s288c_scores.tsv`, `results.txt`: outputs quoted above
