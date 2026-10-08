# 02. Training and testing ledger

*2026-10-06. One row per training, test or validation activity. Use this table to build the Methods
and Results. Column "How checked": **re-derived** = recomputed in the 2026-10-05/06 sessions;
**copied** = taken from the cited report without recomputing. "State" is measured, planned, or not
measured. No row describes wet-lab validation, because none exists.*

## A. Review of the first classifier (2026-09-27)

Source for all A rows: `docs/model-review/2026-09-27-review-and-framework-plan.md` (section in brackets). Scripts: `analysis/model_review/`. How checked: copied.

| ID | Activity | Data and split | Result | State |
|---|---|---|---|---|
| A1 | Cross-validation schemes against feature sets (length; amino-acid composition; ESM-2 8M with two poolings) | 581 positives and 1,500 negatives. Random 5-fold, homology-grouped 5-fold (MMseqs2 30%), leave-family-out | Composition plus length: ROC-AUC 0.999, 0.997, 0.995. ESM-2 8M: 0.999 to 1.000. Length only: 0.796 to 0.841. The benchmark is saturated [4.1] | Measured |
| A2 | Homology structure of the positives | MMseqs2 `easy-cluster --min-seq-id 0.3 -c 0.5` | 581 positives form 53 clusters. The largest has 284 sequences. About 50 independent positives [2] | Measured |
| A3 | Genome-wide check on *S. cerevisiae* S288C | Shipped 8M model; known adhesins and 20 curated hard negatives | 77 calls; 8 of 9 known adhesins found; precision about 12%; 12 of 20 hard negatives called [4.3] | Measured |
| A4 | Fungi_5k screen of the 35M model | 5,802 species, 749,697 positive calls | Median 113 calls per species (about 1.3%). *S. pombe* 38, *C. neoformans* 60, *N. crassa* 136, *A. fumigatus* 114 [1, 4.2] | Measured |
| A5 | Dependence on tandem repeats | 95 curated adhesins | Spearman 0.302 (p = 0.003) between score and repeat coverage. Median coverage 1.000 (found) against 0.000 (missed), p = 0.0037 [4.7] | Measured |
| A6 | Clade transfer | Train on one clade, test on another | ROC-AUC 0.60 across clades; 0.94 for a size-matched model inside the clade [4.5] | Measured |
| A7 | Stage-2 proof of concept (not packaged) | 75 curated adhesins against 63 non-adhesive surface proteins; ESM C 300M frozen embeddings; held-out genomes | PR-AUC 0.94. Architecture alone (GPI, signal peptide, length): 0.70. Holds inside the repeat-mediated class [4.4] | Measured, prototype only |
| A8 | Code review for reproducibility | 2,081 training sequences | Padding tokens were in the mean pool. With the 8M model 7 of 2,081 sequences changed call by batch composition; probability moved by up to 0.57 [3] | Measured; fixed on `main` |

**Logic that A1 to A8 give the paper:** the first benchmark could not detect the problem (A1);
the problem showed at proteome scale (A3, A4); its cause was the label and the data (A2, A5, A6).

## B. Step 1: surface protein prediction

Sources: `docs/step1-plain-language-summary.md`, `docs/model-review/STATUS.md`,
`docs/superpowers/specs/2026-10-01-step1-phase-c-evaluation-design.md`, run outputs in
`_workdir/step1_compare/phasec/` (git-ignored). Phase C was run 2026-10-01; the widened C grid run on 2026-10-02 (job 29350110).

| ID | Activity | Data and split | Result | How checked | State |
|---|---|---|---|---|---|
| B1 | **Phase A:** truth set from GO cell component terms | GAF files per species; IEA dropped; labels P-ext, P-gpi, N-int, N-sec, PM-TM; ambiguous genes left out; UniProt keywords never test truth | Direct-evidence surface genes: *S. cerevisiae* 88, *C. albicans* 211, *S. pombe* 42, *A. fumigatus* 23, *A. nidulans* 113, *C. neoformans* 9, *U. maydis* 10 | copied | Measured |
| B2 | **Phase B:** features and embeddings | SignalP 6, GPI call, amino-acid fractions, ESM-2 8M and 35M (layer 6, mean over residues, first 1,022 aa; last 1,022 aa for "-C" variants) | Feature tables | copied | Measured |
| B3 | **Phase C:** held-out comparison of candidates B0, B1, R0, R1, R2, M8, M35, M8-C, M35-C, H | Training on *S. cerevisiae* and *C. albicans* only. S1: cluster-grouped 5-fold. S2: train on one species, test on another. S3: leave one clade out. 95% cluster bootstrap, 2,000 resamples | See table B3 below | copied | Measured |
| B4 | ML against each rule at the rule's own operating point | S1:all, direct truth, variant V-go (job 29351331) | At R0's recall (0.603) ML FPR is 0.013 to 0.018 against R0's 0.037. At R2's recall (0.418) ML and R2 intervals overlap. B1 is not better than R0 at R0's recall | copied | Measured |
| B5 | Literature rows | 19 Onygenales and Eurotiales adhesins; positives only | Recall 1.00 for R0, B1, M8, M8-C, H; 0.947 for M35 and M35-C; 0.684 for R2; 0.632 for R1. R1 and R2 miss CBP1, CTS1, abr2, aspf2 | copied | Smoke test |
| B6 | Sensitivity to the regularisation grid | C grid 0.001 to 10, then 0.0001 to 10 | Recall and FPR moved by a few points. Largest move: Basidiomycota M8 recall 0.875 to 0.750 on 16 positives | copied | Measured |
| B7 | **R0 per species** (basis of the rule's status) | Phase C per-species sets | See table B7 below | **re-derived** (independent recomputation from `metrics.json`; `calibrate phasec` run with the real NCBI names.dmp and nodes.dmp) | Measured |

### Table B3: Phase C headline results (R0, R2, B1, M8)

*Recomputed from `_workdir/step1_compare/phasec/metrics.json` (test_sets, direct truth, stratum all, variant V-go) on 2026-10-07. All 24 recall and FPR values and the positive and negative counts match. The earlier "not read" cell for B1 on *A. nidulans* is now filled.*

| Test set | Label | Positives / negatives | R0 recall / FPR | R2 recall / FPR | B1 recall / FPR | M8 recall / FPR |
|---|---|---|---|---|---|---|
| S1:all (two yeasts, cross-validated) | estimate | 232 / 4,244 | 0.603 / 0.037 | 0.418 / 0.006 | 0.763 / 0.130 | 0.772 / 0.084 |
| Eurotiomycetes (leave one clade out) | estimate | 128 / 208 | 0.727 / 0.010 | 0.227 / 0.005 | 0.859 / 0.168 | 0.898 / 0.034 |
| *A. nidulans* alone | estimate | 109 / 164 | 0.688 / 0.012 | 0.165 / 0.006 | 0.862 / 0.146 | 0.881 / 0.024 |
| Basidiomycota | smoke test | 16 / 60 | 0.938 / 0.083 | 0.125 / 0.017 | 0.625 / 0.283 | 0.750 / 0.133 |

Findings stored in `findings.json` (copied):
- B1 is not saturated: ROC-AUC 0.894 on S1:all and 0.736 to 0.957 on the three S2 sets.
- No ML candidate beats both B1 and R2 on the false-positive rate for non-secreted proteins at
  R2's recall. This holds on S1:all and on each S2 set.

### Table B7: R0 per species (basis of the status entries)

Sensitivity and specificity are with 95% intervals. The converter reads the intervals from
`metrics.json`. It widens a rate of exactly 0 or 1 with the Wilson interval on protein counts.
Cluster counts exist in `clusters.tsv.gz`. They are not yet recorded in the status entries (open
decision, file `03`).

| Species set | Positives / negatives | Clusters (pos / neg) | Sensitivity | Specificity | Phase C label | Status from the rule |
|---|---|---|---|---|---|---|
| *S. cerevisiae* | 79 / 3,785 | 58 / 3,156 | 0.848 [0.738, 0.943] | 0.966 [0.958, 0.972] | smoke test | `smoke` (label) |
| *C. albicans* | 153 / 459 | 113 / 410 | 0.477 [0.356, 0.591] | 0.943 [0.913, 0.968] | smoke test | `smoke` (label) |
| *A. fumigatus* | 19 / 45 | 17 / 38 | 0.947 [0.833, 1.000] | 1.000 [0.921, 1.000] | smoke test | `smoke` (19 positives) |
| *A. nidulans* | 109 / 164 | 100 / 151 | 0.688 [0.595, 0.779] | 0.988 [0.968, 1.000] | estimate | `estimated` |
| *C. neoformans* | 7 / 32 | 6 / 31 | not checked | not checked | smoke test | `smoke` |
| *U. maydis* | 9 / 28 | 9 / 24 | not checked | not checked | smoke test | `smoke` |

Only *A. nidulans* meets the `estimated` rule. Pooled values hide the spread between species
(sensitivity 0.85 in *S. cerevisiae*, 0.48 in *C. albicans*). This is why a status entry is one
species.

## C. Step 2: mechanism evidence

| ID | Activity | Data and result | Source | How checked | State |
|---|---|---|---|---|---|
| C1 | Repeat detector run blind on known proteins | Period-based, composition-agnostic detector. SOWgp58: period 47, 3.8 copies (published 41 to 47 aa, 4 copies). SOWgp82: period 47, 5.8 copies (published 6) | `docs/reports/2026-09-27-cocci-repeat-surface-proteins.md` | copied | Measured |
| C2 | Repeat surface proteins in *Coccidioides* from long reads | 5 UArizona assemblies. About 9 repeat proteins per long-read strain against about 6 per short-read reference | same | copied | Measured, no experiment |
| C3 | Transfer of repeat features | Composition-agnostic features helped only marginally (2 of 5 to 3 of 5 recovered). SOWgp is short and Pro/Cys-rich; training positives are long and Ser/Thr-rich | same, section 1 | copied | Measured |
| C4 | CFEM, Bys1, hydrophobin HMM scans | HMMs exist. No scan wrapper had a specificity test before the sorting tool | `docs/model-review/STATUS.md` | copied | Not measured (all families inactive in the sorting tool) |
| C5 | PRA3 full-length structure | AlphaFold model, pLDDT 62.2; no Foldseek hit at TM 0.5. The earlier "no fold" result is not explained by truncation. Not proof of no fold | `docs/reports/2026-10-02-pra3-fulllength-structure.md` | copied | Measured, computational |
| C6 | PF28404 family | 831 proteomes searched; 408 proteins in 283 proteomes; 217-protein tree; not specific to *Coccidioides*; function unknown | `docs/reports/2026-10-02-pf28404-family.md` | copied | Measured, computational |

## D. Step 3: purpose layers (*Coccidioides*)

| ID | Activity | Data and result | Source | How checked | State |
|---|---|---|---|---|---|
| D1 | Genus-specificity ranking | 488-proteome pangenome (15,857 orthogroups); 8,542 orthogroup representatives; confounders *Histoplasma*, *Blastomyces*, *Paracoccidioides*, *A. fumigatus*, human. Two axes: antigenicity and specificity. Acceptance test on known antigens: 3 of 4 anchors pass the top-decile test (PRA3, Ag2/PRA, SOWgp pass; PRA2 fails), so the script prints `NOT CALIBRATED`. With the 15% cut all four anchors are inside | `docs/reports/2026-09-27-coccidioides-antigen-findings.md`; orchestrator spec section 3.4, D12 | copied | Measured. **No serology** |
| D2 | Cross-reactivity | Ag2/PRA and PRA2 have orthologs at 55 to 70% identity in the dimorphic confounders. SOWgp and PRA3 have no hit | same, section 4.2 | copied | Computational prediction, untested |
| D3 | Spherule expression | Carlin et al. 2021 RNA-seq (2 replicates per state). SOWgp 13 TPM (mycelia) to 15,000 TPM (spherule 48 h). Ag2/PRA, PRA2, PRA3 go down | same; `docs/reports/2026-10-03-cocci-spherule-surface-table.md` | copied | Measured |
| D4 | Spherule surface table | 564 of 9,630 genes up at 48 h; 23 of the 564 have a signal peptide (4.1%) against 4.7% for all genes; 4 up genes pass the "specific" flag | same | copied | Measured |
| D5 | SOWgp presence and absence | Read depth in 559 CRAMs. The 39 strains with no gene model have normal depth (ratio 0.92 to 1.43). The 8% absence is a gene-model gap | `analysis/cocci_repeats/REPORT_2026-10-01_sowgp_depth_vs_repeats.md` | copied | Measured, computational |
| D6 | Allergen scoping | WHO/IUIS fungal set: 111 usable of 116 rows with a sequence. Rule: 35% identity over 80 aa. No *Coccidioides* allergen is listed. No model was trained | `docs/reports/2026-10-04-fungal-allergen-scoping.md` | copied | Planned (module wrapper only) |
| D7 | Map of cell wall gene classes to tools | Three reviews read from abstracts and headings. No class was run through a tool | `docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md` | copied | Not a test |

## E. Software tests of the sorting tool

| ID | Activity | Result | How checked | State |
|---|---|---|---|---|
| E1 | Core engine unit tests (Kleene tables, status propagation, taxonomy scope, atomic writes, FASTA edge cases) | 184 tests pass | re-derived (pytest run 2026-10-05) | Measured |
| E2 | Module wrappers, calibration commands, job-script checks | 626 pass, 7 skipped (the 7 need `shellcheck`, not installed) | re-derived | Measured |
| E3 | Real-data loader checks | SignalP RS FungiDB file: 9,910 rows, 460 signal peptides (4.6%). Antigen lookup on the RS RefSeq FASTA: 9,139 ok and 771 not in reference (total 9,910). Expression lookup: 9,910 ok. BLAST 2.14.0+ keeps ID forms | re-derived (final review probes) | Measured |
| E4 | Toy end-to-end run (module commands, `calibrate phasec` with real names.dmp, `calibrate truth`, engine) | Engine reads all tables and status files; renamed calls appear; *A. nidulans* `estimated`, others `smoke` | re-derived | Measured on toy data |
| E5 | Real-proteome end-to-end run (*A. fumigatus* Af293), calibrations, strain stability | Plan 2 Tasks 11 to 15 | — | **Not run** (needs HPCC jobs) |

Software correctness is not accuracy. E1 to E4 show that the tool does what its rules say. They do
not show that its calls are right.

## F. Control-set counts (2026-10-06)

| ID | Activity | Result | Source | How checked | State |
|---|---|---|---|---|---|
| F1 | Count of the curated truth for the repeat call (decision C1) | `adhesin` E1: 55 rows, 54 sequences, 32 clusters. E1+E2: 98 rows, 96 sequences, 50 clusters. Hard negatives N1+N2: 42 rows, 33 clusters. Rows with "repeat" or "tandem" in family, summary or name: 8 rows, 5 clusters. 47 of 98 E1/E2 rows have an empty family column. The table has no mechanism label | `analysis/calibration_truth/c1_truth_count.py`, output `c1_counts_2026-10-06.txt` | **re-derived** (sequences fetched live from UniProt; MMseqs2 30% identity, 50% coverage) | Measured. A count, not a calibration |

The count shows that the 20-cluster floor is not met for repeat-mediated adhesins. The missing piece is
a mechanism label, not more rows. See `docs/agent-tasks/03-repeat-mechanism-controls.md`.

## G. Per-call status and the first calibration of the repeat call (2026-10-08)

| ID | Activity | Result | Source | How checked | State |
|---|---|---|---|---|---|
| G1 | Per-call status entries (a status for a call that reads several modules) | Spec, plan, two independent reviews, seven tasks built test first. 782 tests pass (683 before). A call file is used only when config, modules read, identities, run states and the call definition match. | `docs/superpowers/specs/2026-10-07-per-call-status-design.md`, `plans/2026-10-08-per-call-status*.md`, `docs/paper/03` section 6a | tests, about 40 mutation checks | Measured (software) |
| G2 | Repeat call (`repeat02` OR `repeat14`), *S. cerevisiae* S288C | 31 positives (22 clusters), 297 negatives (227): sensitivity 0.323 [0.043, 0.535], specificity 0.987 [0.962, 1.000]. Status `smoke` (cluster bootstrap; leakage `tuned_on_truth`). | `docs/reports/2026-10-08-repeat-call-calibration.md`, `docs/reports/data/sorting_hat/call_status/` | recomputed from the status file | Measured, exploratory |
| G3 | Repeat call, *C. albicans* SC5314 | 16 positives (9 clusters), 210 negatives (155): sensitivity 0.375 [0.000, 0.680], specificity 0.976 [0.938, 0.995]. Status `smoke`. | same | same | Measured, exploratory |
| G4 | Repeat call, *A. fumigatus* Af293 | 8 positives (4 clusters), all enzyme repeat domains (PbH, BNR); 159 negatives (117), none called. No status written: sensitivity for arrays is not measurable from UniProt here. | same, section 4a | same | Not measured for arrays |
| G5 | Labels of the truth set | Positive: paper statement or at least 2 UniProt `Repeat` features. Negative: no `Repeat` feature and no text mention (an **assumed** negative). Family inference never counts. | `analysis/calibration_truth/repeat_call_truth/` | scripts, UniProt release 2026_03 | Assumption stated |
| G6 | Why sensitivity is low | 18 of 31 (S288C) and 8 of 16 (*C. albicans*) positives have no period in either detector. About half are arrays (AGA1, HPF1, SED1, EGT2, MSB2, EAP1, PGA18, ALS5 to ALS7) and half repeat domains. Coverage-cutoff and z-score scans do not fix it. | report section 4; `threshold_scan_v2.tsv`, `zseq_rule_scan.tsv` | recomputed | Measured, exploratory |
| G7 | ALS7 adhesion label | GO IMP cites PMID 17510860 (Als1p and Als5p only); Sheppard 2004 shows no adherence for Als7p. Override E1 to E3. | `data/curated/adhesins/manual_overrides.tsv`, errata | PubMed abstract, extracted paper facts | Owner to confirm |
| G8 | Repeat call under the owner's definition (motif or array; adhesion-associated repeat domains only) | S288C 23 positives (16 clusters), 297 negatives: sensitivity 0.435 [0.067, 0.667], specificity 0.987; *C. albicans* 13 (6), 211: 0.462 [0.000, 0.788], 0.976. `smoke`. All 20 remaining misses are arrays (11 at the z gate, 5 at the score pre-filter, 4 at the region test or coverage). | report section 8; `truth_v3.*`, `detector14_stopping_points_v3.tsv` | recomputed from the status files and scripts | Measured, exploratory |
