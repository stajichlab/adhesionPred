# Design spec: step 1 surface glycoprotein model (ML candidate) and its comparison with the rule

*Drafted 2026-09-30. Revised 2026-10-01 with the owner's decisions (section 10). Revised again on
2026-10-01 after an independent review (fix brief S1 to S13, controller rulings R-A, R-B, R-C).
No code exists for this spec; the Phase A plan is `docs/superpowers/plans/2026-10-01-step1-truth-set.md`.*

Inputs: `docs/PLAN-2026-09-30-pipeline-and-decisions.md` (cited "plan n", sections 4 and 4a); the
review `docs/model-review/2026-09-27-review-and-framework-plan.md` ("review n");
`docs/TOOL-ARCHITECTURE.md`; the Fable design review `docs/plans/2026-09-30-design-review-fable.md`
(branch `issue-9-plan`, "Fable #n"; summarised in plan 4 and 7).

State of `main` that this spec builds on (plan 4a): the package is `surface_glyco`. **No model
ships.** The old pickles are deleted; there is no legacy mode. The model card framework exists in
`src/surface_glyco/card.py` (card version 2, `MAX_RESIDUES = 1022`). Version 0.2.0 is cut only when a
validated model ships. Plan 4a requires this spec and an independent review before any plan or code.

## 1. Purpose and non-goals

**Purpose.** Step 1 asks: does this protein go through the secretory pathway and end up outside the
plasma membrane (cell wall, outer face of the plasma membrane, or secreted)? Step 1 has two candidates
(plan 4a): the rule (SignalP + GPI anchor + Ser/Thr content) and an ESM-based model retrained on a
surface label. This spec defines the ML candidate and one test bench for both. The owner picks rule,
ML or hybrid from the comparison.

**Non-goals.** The model must not be described as:

- an adhesin predictor (old model: about 12% adhesin precision on S288C, review 4.3);
- a glycosylation predictor (section 2.1);
- a moonlighting-protein detector (Hsp60 and enolase reach the surface without a signal peptide);
- an antigen or epitope predictor;
- validated outside the species in the truth set (section 3).

**The name.** The score stays `surface_glycoprotein_score` (decision Q1). The label is a
**localisation** label: GO cell wall or extracellular region. "Glycoprotein" is a loose fit. GO has
no cellular-component term for glycosylation, so the truth set cannot test glycosylation. Most fungal
wall and secreted proteins are glycosylated (assumption, from the literature; not checked here). Every
report and the model card must state that the score predicts location, not glycosylation.

## 2. Label definition

### 2.1 Facts checked

| Fact | Value | Source |
|---|---|---|
| `surface.tsv` rows | 3,573 (3,540 `uniprot_surface_kw`, 33 `curated_literature`) | `data/curated/surface/surface.tsv` |
| Negatives in `surface.tsv` | 0 | same; `data/curated/surface/README.md` |
| Rows longer than 1,022 aa | 164 | same, `length` column |
| Eurotiomycetes literature rows | 21: 10 Onygenales, 11 Eurotiales. YPS3 has no accession. HSP60 is `moonlighting` | `data/curated/adhesins/eurotiomycetes_seeds.tsv` |
| GPI anchoring in GO | no current cellular-component term. GO:0031225, GO:0046658, GO:0031362 are obsolete | `go-basic.obo` releases/2026-07-26 (HTTP Last-Modified 2026-08-08) |
| Glycosylation in GO cellular component | no term | same |

`surface.tsv` labels come from UniProt keywords. Those keywords share inputs (signal-peptide and GPI
predictors) with the rule. So `surface.tsv` is never test truth (decision Q4).

### 2.2 Labels (decided)

All labels use **non-IEA GO evidence** (Q5) and the D1 definition (section 3.1).

| Label | Rule | Decision |
|---|---|---|
| **P-ext** (positive) | non-IEA cell wall (GO:0005618) or extracellular region (GO:0005576), and no non-IEA cytosol, nucleus or mitochondrion term. Column `subset` = `wall` or `extracellular-only` | Q1 |
| **P-gpi** (positive) | P-ext gene with a non-IEA plasma-membrane term and **curated** GPI evidence: a reviewed UniProt entry with experimental evidence (ECO:0000269) for the GPI-anchor lipidation, or a literature row. Predictor output never counts. Reported as a list until literature rows exist (R-B) | Q2 |
| **N-int** (easy negative) | non-IEA cytosol, nucleus or mitochondrion; no endomembrane, plasma membrane, vacuole, membrane, cell periphery, wall or extracellular term at any evidence level | Q5 |
| **N-sec** (hard negative) | non-IEA endomembrane system, plasma membrane or vacuole; no wall or extracellular term at any evidence level | Q3 |
| **PM-TM** (negative stratum) | P-ext gene with a non-IEA plasma-membrane term, a transmembrane segment, and no curated GPI evidence (for example MSB2, HKR1). Labelled negative | Q3 |
| **Ambiguous** | surface term and internal term together (non-IEA) | Q9 |
| Unlabelled | everything else. Never a negative | review 6 |

**Ambiguous genes** (Q9) are excluded from training and from precision. They are a reported stratum.
Example: *C. albicans* ENO1, TDH3, PGK1 and ADH1 have experimental extracellular-region evidence and
also an internal term. Of the extracellular-region genes with an experimental code (305 *C. albicans*,
125 S288C), 102 and 39 carry an internal term at any evidence level (IEA and IBA included). With
experimental internal evidence the numbers are 13 of 305 and 32 of 125 (recomputed on 2026-10-01 from
`cgd.gaf.gz` and `sgd.gaf.gz` with the Phase A filter code).

**High-throughput internal evidence (ruling R-A).** The D1 rule stays unchanged, so it reproduces
3.3. The truth table stores `internal_evidence_htp_only`: true when every non-IEA internal-term
annotation of the gene has code HDA, HMP, HEP, HGI or HTP. Ambiguous genes with this flag are a
reported sub-stratum of "ambiguous": S288C 18 (for example TIP1, CWP2, CCW14, HSP150, ECM33, CTS1),
*S. pombe* 4, other sources 0. They stay ambiguous and outside training until the owner decides.

**Consequence of Q3.** MSB2 and HKR1 have GO wall or extracellular annotations, so D1 places them in
P-ext first (both are in the S288C P-ext list from `/tmp/glyco_spec/d1_count.py`). The PM-TM rule
moves them to the negative side. The model then learns that a transmembrane helix means "not surface".
Both candidates are scored on PM-TM as its own stratum, so this effect is visible. HKR1's
extracellular annotation is IBA only; under direct evidence (3.3) HKR1 is unlabelled, not PM-TM.

**RNA genes.** CGD types RNA genes as `gene_product`, so they pass the object-type filter. 72 of the
1,772 *C. albicans* N-int genes have no protein (tRNA, rRNA and snoRNA genes). The claim that the
filter drops ncRNAs holds for SGD only. SGD has 3 N-int genes without a protein (RAF1, REP1, REP2 of
the 2-micron plasmid); PomBase (`pombase.gaf.gz`) has 1 N-int and 1 N-sec gene without a protein.
Phase A counts these genes as unmatched and keeps them out of sequence-level sets.

**PM-TM and P-gpi need a triage step (D8).** P-ext genes with a non-IEA plasma-membrane term: S288C
12, *C. albicans* 60, *S. pombe* 10, *A. fumigatus* 9. The S288C list mixes GPI proteins from the
literature (GAS3, YPS1; both carry only ECO:0000255, predicted, GPI sites in UniProt) with TM proteins
(MSB2, HKR1). D8 splits them with curated GPI evidence and the UniProt
transmembrane feature. For unreviewed entries that feature probably comes from a predictor
(assumption, not checked). MSB2 and HKR1 TM features are ECO:0000255 (predicted) even in reviewed
entries, so "PM-TM" rests on predicted topology.

**Curated GPI evidence is nearly absent.** UniProt 2026_03 reviewed entries with ECO:0000269 GPI-anchor
evidence: S288C 3 (GAS1, TIP1, YNL190W), *C. albicans* 3 (HWP1, YWP1, ECM331), *S. pombe* 0,
*A. fumigatus* 0, *A. nidulans* 1 (chiA). GAS1 and TIP1 are ambiguous under D1. A prototype D8 run
against UniProt 2026_03 on 2026-10-01 gave P-gpi / PM-TM / unresolved: S288C 0 / 3 / 9,
*C. albicans* 1 / 22 / 37, *S. pombe* 0 / 4 / 6, *A. fumigatus* 0 / 1 / 8, *A. nidulans* 1 / 2 / 3.
The independent reviewer re-ran S288C and *S. pombe* only.

**P-gpi rule (one rule, ruling R-B).** D8 triages only P-ext genes with a non-IEA plasma-membrane term.
Genes outside P-ext with curated GPI evidence (for example GAS1 and TIP1, both ambiguous) are listed in
`d8_gpi_outside_pext.tsv` and are not relabelled. P-gpi is reported as a **list, not a scored
stratum**, until `curated_gpi.tsv` holds literature rows (PMIDs for GAS1, TIP1, YPS1, SAG1, CWP1,
CCW12 and others). That curation is separate work. Owner decision needed: fund it or not. Until then
Q2 yields almost no positives.

**Evidence tiers for positives.** T-a: curated GO, non-IEA. T-b: literature rows. T-c: UniProt
keyword only. T-c is training-only and never test truth (Q4).

### 2.3 Training variants (Q4)

| Variant | Training positives | Test truth |
|---|---|---|
| V-go | T-a + P-gpi | GO truth (same for both) |
| V-kw | V-go + T-c (`surface.tsv` keyword rows) | GO truth (same for both) |

Before V-kw training, D10 removes from T-c every accession that is a test protein: Q7 literature rows,
all *S. pombe* proteins, Basidiomycota truth proteins, and the test fold in each CV split. Removal is
by accession and by exact-sequence hash. The variants apply to every trained candidate (B1, M8, M35,
their C-terminal variants, H, and the rule's fitted threshold). This **doubles** the evaluation and
the training compute (sections 6 and 8).

## 3. Data sources and feasibility

### 3.1 D1 definition (the counting rule)

The counts below come from `/tmp/glyco_spec/d1_count.py` (SHA-256 `25021b6a…`), run with
`/usr/bin/python3.12` on 2026-10-01. D1 must implement exactly this rule and reproduce these numbers on
the same files:

1. Parse `go-basic.obo` releases/2026-07-26 (HTTP Last-Modified 2026-08-08). Ancestors = transitive closure over `is_a` and
   `relationship: part_of`.
2. Keep GAF rows whose column 1 is the file's most frequent database (SGD, CGD, PomBase, UniProtKB)
   and whose object type is `protein` or `gene_product`. This drops ComplexPortal complexes,
   RNAcentral RNAs, SGD ncRNAs and the 188 UniProtKB IBA rows in `cgd.gaf.gz` that duplicate CGD genes.
   For `cgd.gaf.gz`, keep taxon 237561 only.
3. Keep aspect `C`. Skip rows with `NOT` in the qualifier. Skip obsolete terms (none occurred).
4. Gene ID = column 2. Label terms = ancestors of non-IEA rows. "Any evidence" = ancestors of all rows.
5. Apply the rules in 2.2. IEA fraction = IEA share of distinct (gene, term, evidence) triples, aspect C.

**Correction to the earlier draft.** The draft counted complexes, ncRNAs and duplicate CGD UniProtKB
rows. S288C P-ext changes from 128 to 125, N-int from 3,358 to 2,645, N-sec from 1,674 to 1,548.
*C. albicans* P-ext changes from 261 to 259, N-int from 1,795 to 1,772, N-sec from 970 to 961.

### 3.2 Files checked (HEAD and download on 2026-10-01; data in `/tmp/glyco_spec/`, not in the repo)

The GO index `https://current.geneontology.org/annotations/gaf/index.html` lists UniProt-based GAFs.
Fungal files include ASPFU, EMENI (*A. nidulans*), CRYD1 (*Cryptococcus* JEC21), MYCMD (*Ustilago
maydis*), PUCGT, BATDJ, NEUCR, SCHPO, YEAST and CANAL. No H99 (CRYNH) file and no Onygenales file is
listed. `aspgd.gaf.gz` still returns 403. The
earlier claim "no GO truth for Eurotiomycetes" was wrong: ASPFU and EMENI exist.

| File | HTTP | Generated |
|---|---|---|
| `annotations/sgd.gaf.gz`, `cgd.gaf.gz`, `pombase.gaf.gz` | 200 | 2026-05-21 |
| `gaf/SCHPO-mod.gaf.gz` | 200 | 2026-07-28 |
| `gaf/ASPFU-uniprot.gaf.gz`, `EMENI-uniprot.gaf.gz`, `CRYD1-uniprot.gaf.gz`, `MYCMD-uniprot.gaf.gz` | 200 | 2026-07-28 |
| EBI GOA `20846.C_neoformans_JEC21.goa`, `313589.C_neoformans_var_grubii_H99.goa` | 200 | 2026-07-28 |

The root-level `sgd`/`cgd`/`pombase` files are older than the `gaf/` files. D1 pins one file per
species and records its hash. `SCHPO-mod` gives 57 / 2,966 / 1,212 / 14, close to `pombase.gaf.gz`.

### 3.3 Measured counts (D1 rule)

| Species (file) | P-ext (wall / ext-only) | N-int | N-sec | Ambig. | IEA fraction (C) | Direct evidence: P-ext / N-int / N-sec |
|---|---|---|---|---|---|---|
| *S. cerevisiae* S288C (sgd) | 125 (104 / 21) | 2,645 | 1,548 | 48 | 0.262 | 88 / 2,360 / 1,457 |
| *C. albicans* SC5314 (cgd) | 259 (122 / 137) | 1,772 | 961 | 100 | 0.158 | 211 / 179 / 263 |
| *S. pombe* (pombase) | 59 (40 / 19) | 2,962 | 1,214 | 14 | 0.292 | 42 / 2,799 / 886 |
| *C. neoformans* H99 (GOA) | 11 (3 / 8) | 23 | 25 | 0 | 0.986 | 9 / 17 / 15 |
| *C. deneoformans* JEC21 (GOA; CRYD1 identical) | 32 (10 / 22) | 1,813 | 817 | 0 | 0.520 | 0 / 3 / 1 |
| *A. fumigatus* Af293 (ASPFU) | 132 (35 / 97) | 2,018 | 1,062 | 1 | 0.533 | 23 / 18 / 26 |
| *A. nidulans* (EMENI) | 211 (43 / 168) | 2,054 | 1,066 | 38 | 0.517 | 113 / 97 / 66 |
| *U. maydis* (MYCMD) | 62 (14 / 48) | 1,712 | 827 | 1 | 0.519 | 10 / 6 / 21 |

"Homology codes" = IBA, IBD, IKR, IRD, ISS, ISO, ISA, ISM, RCA.

**Direct-evidence truth is an intersection.** It is the subset of genes whose label is the same under
all non-IEA evidence and under evidence without homology codes (`homology_only == no`, that is
`label == X and homology_only == "no"`). The last column gives these counts (recomputed on 2026-10-01
from `truth_set.tsv.gz` built by the Phase A extractor on the pinned files). Recomputing labels on
direct evidence alone is NOT used: it moves genes whose only internal term is IBA or ISS into P-ext
(83 in *C. albicans*, including TDH3, PGK1, ADH1, FBA1, TPI1, SSA2 and ribosomal proteins), which
contradicts Q9. Headline metrics and gates use the direct-evidence subset. All non-IEA truth is
reported beside it.

**IBA dominates outside the model yeasts.** In JEC21, all 33 non-IEA surface evidence rows on P-ext
genes are IBA, and no P-ext gene has experimental evidence. In ASPFU, 100 of 146 such rows are IBA
(EXP 14, IDA 14, ISS 18). IBA labels are transferred along gene trees from annotated homologs
(GO evidence definition). Which species the source annotations come from is not verified. If they come
from yeasts, a leave-clade-out test on IBA truth partly measures agreement with yeast annotation. So
held-out clades report two numbers: all non-IEA truth, and direct-evidence truth (the intersection
above).

**ID-to-sequence mapping.**

| Species | GAF gene ID | Sequence source (same provider) |
|---|---|---|
| S288C | SGD `S000…` | SGD protein FASTA (`downloads.yeastgenome.org`); header field `SGDID:`; 6,052 of 6,056 genes match |
| *C. albicans* | CGD `CAL…`; the FASTA key is the first synonym in GAF column 11 (`C1_00010W_A`) | CGD protein FASTA (`www.candidagenome.org`). Its headers hold only the systematic name, no `CAL` ID. 6,060 of 6,313 genes match; the 253 unmatched are RNA genes and unlabelled genes, none P-ext or ambiguous |
| *S. pombe* | PomBase systematic ID `SP…` | PomBase `peptide.fa.gz` (IDs carry a `.1:pep` transcript suffix); 5,020 of 5,025 `pombase.gaf.gz` genes match |
| H99, JEC21, Af293, *A. nidulans*, *U. maydis* | UniProtKB accession (column 3 holds locus tags such as `CNAG_`, `CNA…`, `AFUA_`) | UniProt reference proteomes UP000010091 (H99, 7,427 proteins), UP000002149 (JEC21, 6,740), UP000002530 (Af293, 9,647); queried at `rest.uniprot.org/proteomes` on 2026-10-01 |

Match counts come from the Phase A prototype run on the FASTA files downloaded on 2026-10-01. Version
skew between GAF and FASTA is checked in D1 (unmatched IDs are counted and reported).

### 3.4 Feasibility verdict

- **S288C, *C. albicans*:** usable (Q6 training species).
- ***S. pombe*:** usable as a leave-species-out test species: 59 positives, 4,176 negatives.
- ***A. fumigatus*:** partly usable. Non-IEA gives 132 / 2,018 / 1,062, but only 23 / 18 / 26 with
  direct evidence. *A. nidulans* has more direct truth (113 / 97 / 66) and is a second S3 test species. The literature rows (about 19
  usable) remain a smoke test (Q7).
- ***Cryptococcus*:** **not usable from GO alone.** H99 has non-IEA cellular-component labels for 77
  genes (11 P-ext). JEC21 has many non-IEA labels, but all surface evidence is IBA and there is no
  experimental P-ext gene. JEC21 is *C. deneoformans* (serotype D), not *C. neoformans* sensu stricto
  (UniProt proteome name).
- **No Onygenales GO file** exists. Onygenales truth is the literature rows only.

**Options for Basidiomycota truth (Q8, for the owner; not chosen here):**

1. Literature curation of *C. neoformans* H99 cell-wall and secreted proteins, plus curated
   negatives. The 11 H99 P-ext genes (CDA1-3, LAC1, LAC2, PLB1, CPL1, QSP1, PQP1, YOR1, cnap1) are a
   start. Size and sources not estimated.
2. Use JEC21 IBA truth (32 / 1,813 / 817) as a separately reported "homology-transfer" stratum.
3. Combine 1 and 2: literature positives from H99, IBA negatives from JEC21.
4. Add *U. maydis* (MYCMD) direct truth (10 / 6 / 21) as a second, small Basidiomycota set.
5. Check FungiDB community GO curation for H99 (not checked).

## 4. Dataset construction

1. **Extract** labels with D1 and D9. Record each file's date and SHA-256.
2. **Attach sequences** (3.3 mapping). Record the release. Count unmatched IDs.
3. **Triage** P-ext ∩ plasma membrane with D8 (P-gpi, PM-TM).
4. **Dedupe.** SHA-256 of each cleaned sequence. Keep one copy within a class. Drop any hash in both
   classes (plan default 4; existing test `test_dedupe_drops_sequences_present_in_both_classes`).
5. **Overlap removal** for V-kw (D10, section 2.3).
6. **Cluster** all labelled proteins from all species together with MMseqs2 `easy-cluster
   --min-seq-id 0.3 -c 0.5 --cov-mode 0` (as in `analysis/model_review/run.sh`).
7. **Class balance.** Class weights. Training ratio (P-ext 384, negatives 6,926 before dedupe and
   triage) is not genomic prevalence, which is not measured. Report precision on the truth set and
   re-weighted to a prior estimated on whole proteomes.
8. **Splits.**
   - **S1:** homology-grouped 5-fold (`StratifiedGroupKFold` on cluster ID), S288C and *C. albicans*
     pooled.
   - **S2:** leave-species-out: train S288C, test *C. albicans*; the reverse; train both, test
     *S. pombe*. *S. pombe* is never in training (Q6).
   - **S3:** leave-clade-out: train both yeasts, test Eurotiomycetes (literature rows; ASPFU, and
     EMENI, decided 2026-10-01) and Basidiomycota (the truth chosen under Q8). Report the maximum
     identity of each test protein to the training set.
   - Hyperparameters and rule thresholds are tuned in an inner loop on training folds only.

| Leakage risk | Control | Test |
|---|---|---|
| Homologs across folds | cluster-grouped folds | `test_no_cluster_spans_train_and_test` |
| Same sequence in both classes | drop hash | existing dedupe test |
| IEA labels from SignalP/InterPro | IEA excluded | `test_truth_set_has_no_iea` |
| UniProt keywords share inputs with the rule | T-c never test truth | `test_test_truth_sources` |
| T-c contains test proteins (Q4) | D10 removal by accession and hash | `test_tc_excludes_test_proteins` |
| IBA truth transferred from training species | headline metrics on direct-evidence truth (intersection, 3.3); all non-IEA truth beside it | `test_direct_evidence_counts_on_real_files` (Phase A) |
| Thresholds tuned on test data | nested loop | `test_threshold_fit_uses_train_only` |

## 5. Model candidates and baselines

| Candidate | Features | Cannot see |
|---|---|---|
| B0 length | length | everything else |
| B1 composition | 20 AA fractions + length, StandardScaler + LR | order, motifs, signal position |
| R rule | SignalP 6 SP call, GPI call (PredGPI `predgpi/202001`; NetGPI not installed, `module avail`), Ser+Thr fraction | anything else |
| M8, M35 | ESM-2 8M and 35M, layer 6, residue mean (padding excluded), StandardScaler + LR | residues past 1,022: `embeddings.py` keeps the first `MAX_RESIDUES` |
| M8-C, M35-C | same, on the last 1,022 residues (Q10) | the N-terminus of long proteins, including the signal peptide |
| H hybrid | best ESM variant + SignalP probability, GPI score, Ser+Thr fraction in one LR | limits of its parts |

Each trained candidate runs as V-go and V-kw (section 2.3). M8-C and M35-C are separate candidates,
not a merged score (Q10). The >1,022 aa subset is reported on its own.

**SignalP coverage in Onygenales.** A direct SignalP 6 run calls 460 of 9,910 (4.6%) *C. immitis* RS
proteins, and calls SOWgp (CIMG_04613) as a signal peptide with probability 0.9998
(`docs/TOOL-ARCHITECTURE.md`, known-gap row; `analysis/cocci_repeats/signalp_summary.tsv`). No truth
data show whether 4.6% is an under-call. The SignalP run must come from a tracked script (D2).

**Why the old CV is useless.** On the old task B1 reached ROC-AUC 0.995 under leave-family-out (review
4.1). The new CV must show: (a) B1 does not saturate it; (b) ML beats B1 and R on N-sec, with the
difference outside the bootstrap CI; (c) this holds under S2. If ML does not beat B1, the report says so.

**Rule threshold.** No step 1 Ser+Thr cut-off exists. The 25% in
`analysis/cocci_repeats/03_repeat_surface_candidates.py` is for class 2a. The cut-off is fitted in the
inner loop.

## 6. Evaluation

- **Runs:** every candidate × {V-go, V-kw} × {S1, S2, S3}. The rule and B0 have no V-kw difference
  except the fitted threshold. All runs are scored on the same GO test truth.
- **Metrics per species and clade:** recall, precision, FPR, ROC-AUC, PR-AUC. 95% CI from 2,000
  bootstrap resamples by cluster. Clades: Saccharomycotina (S288C, *C. albicans*),
  Taphrinomycotina (*S. pombe*), Eurotiomycetes, Basidiomycota.
- **Strata:** wall; extracellular-only; N-int; N-sec; PM-TM; ambiguous (score distribution
  only, no precision) with its sub-stratum `internal_evidence_htp_only` (R-A); >1,022 aa;
  direct-evidence truth (headline). P-gpi is a list, not a scored stratum, until `curated_gpi.tsv`
  has literature rows (R-B).
- **Variant effect (Q4):** V-kw minus V-go for each metric, with paired bootstrap CI.
- **Precision at fixed recall** (0.8, 0.9) and **recall at FPR 0.01**.
- **Calibration:** reliability curve and Brier score on S2 test sets.
- **Agreement matrix:** rule versus ML on the truth set and on whole S288C, *C. albicans*,
  *S. pombe* and *C. immitis* RS proteomes. Discordant cells list protein IDs.

**CI width (not measured).** *S. pombe* has 59 positives (42 with direct evidence). S3 direct-evidence
sets have 9 to 113 P-ext genes. Intervals there will be wide. **Estimate or smoke test (rule fixed on
2026-10-01, before any run):** a test set is an "estimate" when the 95% cluster-bootstrap interval for
recall has half-width <= 0.10; otherwise it is a "smoke test".

**Named error-analysis panel.** IDs found in repo files on 2026-09-30, except GO-only IDs.

| Protein | ID | Why |
|---|---|---|
| FLO1 | P32768 | 1,537 aa GPI adhesin (truncation; C-terminal variant) |
| SAG1, CWP1, CCW12 | P20840, P28319, Q12127 | GPI wall proteins (CCW12 is 133 aa) |
| GAS1 | P22146 | GPI enzyme with ECO:0000269 GPI evidence; ambiguous under D1 (mitochondrion HDA, nuclear periphery IDA), so listed outside P-ext, not P-gpi |
| PIR1 | Q03178 | wall protein without GPI; ambiguous under D1 (internal IDA) |
| MSB2, HKR1 | P32334, P41809 | PM-TM negatives (Q3); HKR1 is 1,802 aa; HKR1 surface evidence is IBA only (outside direct evidence) |
| SUC2, PHO5 | P00724, P00635 | secreted enzymes; SUC2 is ambiguous (cytosol, HDA and RCA evidence); PHO5 is P-ext, subset wall |
| ALS3, HWP1, SAP9 | Q59L12, P46593, Q59SU1 | *C. albicans* GPI proteins; P-ext wall |
| ECM33 | A0A1D8PCY4 | *C. albicans*; `gpi_anchor=no` in keywords; P-ext wall. S288C ECM33 is ambiguous (HDA-only internal evidence) |
| ENO1, TDH3 | CAL0000185645, CAL0000197744 | *C. albicans*; ambiguous stratum; TDH3's internal term is IBA only (`homology_only == yes`). S288C ENO1 is N-sec |
| SOWgp58, SOWgp (RS) | Q8NK60, CIMG_04613 | Onygenales Pro-rich surface protein |
| CTS1, CspA, cfmA | Q1E3R8, Q4WXC4, Q4WLB9 | Eurotiomycetes literature rows; CspA and cfmA are also P-ext wall in ASPFU |
| HSP60 | P50142 | moonlighting; must not be positive |

**Measure-first gating** (plan decision 6). Run all candidates. Write
`results/step1_compare/metrics.json` with values and CIs. The owner chooses rule, ML or hybrid and one
training variant. The owner sets each gate at or just below the measured value. Gates go into
`tests/gates/step1_gates.json`. No gate value exists before that.

## 7. Test framework

Already on `main` (`tests/surface_glyco/`): batch-mate independence
(`test_embedding_does_not_depend_on_batch_mates`), input-order restore, skipped-sequence alignment,
truncation count, card refusal before unpickling, model path not from cwd, dedupe across classes.

| Risk | New test | Type | Runs in |
|---|---|---|---|
| C-terminal window wrong | `test_cterm_window_takes_last_1022` | unit | CI |
| GO parsing drifts from 3.1 | `test_gaf_extract_golden` (50-line fixture, incl. complex, ncRNA, duplicate-DB and NOT rows) | unit | CI |
| Label rules wrong | `test_label_rules_truth_table` (P-ext, PM-TM, P-gpi, ambiguous, N-int, N-sec) | unit | CI |
| Rule logic | `test_rule_truth_table` | unit | CI |
| Leakage (section 4) | leakage tests on fixture tables | unit | CI |
| Output drift | `test_golden_scores` on a 5-sequence FASTA | integration (`slow`) | CI |
| GPU vs CPU difference (not measured) | `gpu_cpu_diff.py`, max difference | harness | HPCC |
| Accuracy regression | `step1_gates` against `metrics.json` | tier harness | HPCC |

CI uses CPU and Python 3.12 (plan default 2). The analysis-tests CI job has `continue-on-error: true`,
so the `tests/step1_compare` tests are **informational in CI**: they cannot fail the build. Moving them
into a required job is an owner decision.

## 8. Compute plan

Sizes: S288C 6,722 proteins (plan 5); *C. immitis* RS 9,910; *S. pombe* 5,130, H99 7,427, Af293 9,647
(UniProt proteome counts); *C. albicans* not counted. Labelled set before dedupe: about 7,300 genes in
the two training yeasts.

Review 9.1 gives ESM-2 150M at 60 proteins/s on one RTX 6000 Ada; it has no 8M or 35M rate.
Assumption: 8M and 35M are at least that fast. J0 measures the real rate before J1 to J3 are sized.
Embeddings do not depend on the training variant, so Q4 doubles only J3 (training and evaluation).
The C-terminal variants (Q10) double J2 for proteins longer than 1,022 aa only.

| Job | Resource | Content | Length |
|---|---|---|---|
| J0 pilot | `exfab`, 1 GPU | 2,000 proteins, 8M and 35M, throughput JSON | under 15 min (assumption) |
| J1 features | `exfab`, 1 GPU | SignalP 6, PredGPI on all truth species and *C. immitis* RS | sized to 1 to 1.5 h |
| J2 embed | `exfab`, 1 GPU | 8M, 35M, C-terminal windows | sized from J0 to 1 to 1.5 h |
| J3 evaluate | CPU | MMseqs2, CV, 2 variants, bootstrap, metrics JSON | not measured |

Jobs write to `${SCRATCH:?}`, compress large tables with `zstd`, and copy results to
`analysis/step1_compare/` before exit. Scripts take `PROJ_ROOT` from the environment, not from
`BASH_SOURCE`. This spec submits no job.

## 9. Risks and unknowns

| Risk | How it is detected |
|---|---|
| GO cell-wall terms include moonlighting enzymes | ambiguous stratum; ENO1, TDH3 in panel |
| *C. albicans* P-ext includes cytosolic-looking proteins without an internal term (for example CDC19, BMH1, CEF3, UGP1 in the P-ext ∩ plasma-membrane list) | D8 triage list; per-species precision |
| IBA truth in held-out clades is homology transfer | direct-evidence stratum (intersection, 3.3) |
| PM-TM rule teaches "TM helix = not surface" (Q3) | PM-TM stratum score distribution |
| TM feature for PM-TM comes from a predictor | D8 records evidence code per protein |
| Yeast-only training does not transfer | S2 (*S. pombe*) and S3 with CI width |
| ML learns "signal peptide" only, or composition explains all | N-sec FPR of M8/M35 versus R and B1 |

## 10. Decisions (owner, 2026-10-01)

| Question | Decision | Consequence for the design | Alternatives rejected |
|---|---|---|---|
| Q1 Positive scope | GO wall or extracellular region (P-ext); name stays `surface_glycoprotein_score` | `subset` column; recall per subset; reports say "localisation label" | P-wall only: drops secreted enzymes (SUC2, PHO5) that the rule cannot separate either |
| Q2 GPI plasma-membrane proteins | Positive with curated GPI evidence | P-gpi; D8; nearly empty until literature rows exist (2.2, R-B) | predictor GPI calls: circular with the rule |
| Q3 TM proteins with large extracellular domains | Negative | PM-TM stratum; model learns TM = not surface | positive; exclude |
| Q4 UniProt-keyword tier | Both variants; T-c never test truth; test proteins removed from T-c | V-go, V-kw; D10; J3 and evaluation doubled | one variant only |
| Q5 GO evidence | All non-IEA codes | D1 rule (3.1); direct-evidence truth = intersection `label == X and homology_only == no` (3.3) for headline metrics; all non-IEA truth reported beside it | experimental only: *C. albicans* N-int falls from 1,772 to 168 |
| Q6 Species | S288C and *C. albicans* train; *S. pombe* test only | S2 adds *S. pombe* (59 / 2,962 / 1,214 / 14) | two species only |
| Q7 Eurotiomycetes literature rows | Test only (smoke test), Q4 overlap rule | S3; removed from T-c | add to training: S3 impossible |
| Q8 Outside Ascomycota | Curate Basidiomycota (*C. neoformans*) truth before the first release | GO alone insufficient (3.4); options; D9 | "not validated" flag only |
| Q9 Ambiguous genes | Excluded from training and precision; own stratum | 2.2, 6 | positive: moonlighting enzymes enter the positive class |
| Q10 Long proteins | C-terminal window (last 1,022 aa), separate candidates | M8-C, M35-C; >1,022 aa reported alone | sliding windows |

## 11. Deliverables and definition of done

| ID | Deliverable |
|---|---|
| D1 | `analysis/step1_compare/01_extract_go_truth.py` implementing 3.1; `truth_set.tsv.gz` (gene, species, label, subset, tier, stratum, evidence codes, homology-code flag `homology_only`, `internal_evidence_htp_only`, file hash). The `cluster` column is added by the later dataset plan, not by D1 |
| D2 | tracked SignalP 6 and PredGPI job and feature table for all truth species and *C. immitis* RS |
| D3 | embedding job (8M, 35M, C-terminal windows) with throughput JSON |
| D4 | evaluation script writing `metrics.json` (both variants) and the agreement matrix with IDs |
| D5 | report: per-species, per-clade, per-stratum table with CIs; variant effect; named panel |
| D6 | the tests in section 7, passing (informational in CI; see section 7) |
| D7 | after owner review: `tests/gates/step1_gates.json` and a card (`card.py` framework) whose `validation` entries come from `metrics.json` |
| D8 | GPI and TM triage of P-ext ∩ plasma membrane: curated-GPI list with source per protein (P-gpi, a list until literature rows exist), PM-TM list with TM evidence code; counts; `d8_gpi_outside_pext.tsv` |
| D9 | *S. pombe* truth extraction (D1 rule) and Basidiomycota truth set per the owner's Q8 option, with sources |
| D10 | T-c training table for V-kw with the test-protein removal log (accession and hash) |

Done means: D1 to D6 and D8 to D10 exist and pass; the owner has the comparison; no accuracy
statement appears in a README or card unless it comes from `metrics.json`. No model ships before the
Basidiomycota truth is in the evaluation (Q8).

Not in scope: steps 2 and 3; adhesion labels; fine-tuning (review 8.3 gate not met); ESM C 300M;
Fungi_5k re-scoring; any gate value chosen before measurement.

## 12. Questions for the owner

Closed on 2026-10-01 (review fix brief):

- Homology-code labels: closed by the direct-evidence intersection (3.3).
- *A. nidulans* (EMENI) is a second Eurotiomycetes test species in S3.
- H99 is the *Cryptococcus* reference. A FungiDB H99 GO check comes first, then literature curation;
  both are separate later work.
- Estimate versus smoke test: the half-width rule in section 6.

Controller rulings (the owner can reverse them): R-A (keep the D1 rule; report
`internal_evidence_htp_only` as a sub-stratum of ambiguous), R-B (P-gpi is a list until
`curated_gpi.tsv` has literature rows), R-C (the shipping gate in section 11 is unchanged).

Open:

1. Which Basidiomycota truth option (3.4, options 1 to 5)?
2. Fund the GPI literature curation for `curated_gpi.tsv` (R-B)?
3. Move the ambiguous genes with high-throughput-only internal evidence (R-A) into training or keep
   them as a reported sub-stratum?
4. Move `tests/step1_compare` into a required CI job?
