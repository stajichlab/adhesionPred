# 04. Limits and open questions

*2026-10-06. Sources: `docs/model-review/STATUS.md` section 3, the orchestrator spec sections 4
and 7, the Plan 2 review notes of 2026-10-05, and the conversation of 2026-10-06.*

## 1. What the paper must not claim

| Do not claim | Reason |
|---|---|
| That any candidate is an adhesin, antigen or allergen. | All lists are hypotheses. No wet-lab validation exists for any candidate. |
| That step 1 is an adhesin predictor. | It predicts location. The first classifier detected tandem repeats and secreted Ser/Thr-rich proteins (ledger A3, A5). |
| That a trained model ships or that the rule/ML choice is final. | No model ships. A hybrid gate was chosen in principle on 2026-10-02. No gate values are set. |
| That `serodiagnostic_marker_candidate` or `cocci_specificity_rank_top15` predicts epitopes or serological performance. | It is a genus-specificity ranking made from sequence comparison. The testable serology prediction (sera from histoplasmosis, blastomycosis and paracoccidioidomycosis patients should react with Ag2/PRA and the CF antigen, and not with SOWgp or PRA3) has not been tested. |
| That allergen similarity means IgE cross-reactivity. | It is similarity to a reference set at two thresholds. IgE binding is a hypothesis at most. WHO/IUIS lists no *Coccidioides* allergen. |
| That a module is accurate in a taxon without a status entry. | Status is per species. Unlisted taxa are `unvalidated`. |
| A pooled accuracy across species or clades. | Species differ (R0 sensitivity 0.85 in *S. cerevisiae*, 0.48 in *C. albicans*). The spec allows no pooled number unless each clade passes the `estimated` rule. |
| That the Phase C numbers apply to *Coccidioides*. | Step 1 has no Onygenales truth set. Whether SignalP under-calls in *C. immitis* is untested. |

## 2. What is not measured

From the Phase C report (copied from `docs/model-review/STATUS.md`):
- Basidiomycota performance. *Cryptococcus* (7 positives) and *Ustilago* (9) are smoke tests. No curated truth.
- Performance on GPI-anchored proteins. `curated_gpi.tsv` has no literature rows.
- Whole-proteome prevalence of surface proteins. The prevalence table uses assumed values.
- Variance from refitting the models. The intervals describe test-set sampling only.
- The effect of pooling five fold models in the S1:all ROC-AUC.

From the sorting tool (2026-10-05 reviews):
- Overlap between the Phase C positives and the SignalP 6 training data (decision 1 records it as a stated limit).
- Accuracy of the repeat detectors, family scans, antigen lookup and allergen similarity outside the
  cases listed in the ledger. No module has HPCC calibration yet.
- Strain-to-strain stability of calls for *A. fumigatus* (A1163, W72310). Planned (Plan 2, Task 14).
- Agreement between two modules of the same kind (kappa). Planned.
- The cause of 15 failing repository tests and 91 test errors elsewhere in the repository. The
  lists are identical on the base commit and the branch. The environment lacks several optional
  packages (`esm`, `seaborn`, `statsmodels`, `umap`, `scikit_posthocs`). I did not diagnose them.
- The new CI job, and the job scripts under SLURM. Neither has run.
- Whether TMHMM returns a non-zero exit code on bad input.

Known data limits that affect interpretation:
- Outside the yeasts, direct-evidence positives are few: 9 in *C. neoformans*, 10 in *U. maydis*, 23 in *A. fumigatus*.
- The truth set is built from GO annotations. Most labels outside yeast are homology-transferred. The headline numbers use direct evidence only.
- MSB2, HKR1 and SAP9 are labelled negative and most candidates call them positive. This was not analysed.
- 126 *C. immitis* RS isoforms inherit the antigen and expression values of another transcript
  (`idmap_method` is `gene_best_transcript`). An exact-sequence match is not implemented.
- The taxon IDs 330879, 451804 and 746128 were checked on 2026-10-06 against `names.dmp` and `nodes.dmp`:
  Af293, A1163 and the species *A. fumigatus*. UniProt agrees for Af293 and A1163. See `05-literature-verification.md`.
- Composite calls (for example `cell_wall_adhesion_candidate`) take the weakest status of the measured
  leaf calls that decided them. That status is a derived value. It is not a measurement of the composite
  call. The owner confirmed this reading on 2026-10-07 (`03-status-and-validation-rules.md`, decision 3).
  The report should label such statuses as derived (planned with the per-call status work).
- A module with a bare SignalP version string (no `module=` token) is checked only by major
  version and the `gpu` tag.

## 3. Decisions still open for the owner

| # | Decision | Recommendation | Where |
|---|---|---|---|
| D-B, D-C, D-D | Basidiomycota label unit; whether a curated row may outweigh GO transfer evidence; how to build JEC21 and *U. maydis* transfer strata | The specs give a recommendation for each | `docs/step1-plain-language-summary.md` section 9 |
| Gate values for step 1 | Set after the owner reads `metrics.json` | None before measurement | Phase C spec section 1 |
| Whether to fund GPI and Basidiomycota curation | Literature work; no effort estimate | | `docs/model-review/STATUS.md` section 4 |
| Whether `phasec` should run on A1163, W72310 and the Fungi_5k *A. fumigatus* annotation | Proposed: no. Show the Af293 UniProt numbers as a cross-annotation reference only | Plan 2, Task 13 text marked "proposed, owner to confirm" | |
| Whether an ML model enters the sorting tool | Needs a model card and a freeze | Plan 2 | |
| Per-call status entries (needed to calibrate calls that read several modules) | Change the status file format in the core engine | Decision 3, file `03` | |

Decisions 1 to 4 of 2026-10-06 are in `03-status-and-validation-rules.md`, section 6.

## 4. Items that must be verified before they enter the paper

1. Every "copied" number in `02`. Table B3 was recomputed from `metrics.json` on 2026-10-07 and matches
   (`05-literature-verification.md`, section 4). The other "copied" rows are still unchecked.
2. The literature statements (Vaknin 2014, Liu 2016, Balajee 2007, Fedorova 2008, Gravelat 2013, the A1163
   lineage, gp43) were checked against PubMed, NCBI and UniProt on 2026-10-06. All papers exist. The sentence each
   one supported was not recorded, so compare your sentence with `05-literature-verification.md`. Gravelat 2013
   describes a polysaccharide adhesin. The gp43 glucanase claim is not checked.
3. The three cell wall reviews (Gow 2017, Riquelme 2020, Gow 2023) were read from abstracts,
   editorial text and heading summaries. Read the full text before citing a gene class.
4. The WHO/IUIS counts (120 molecules, 31 species) are from downloads of 2026-10-04.
5. Whether the study design (training on two yeasts only, testing elsewhere) is the design the paper
   should present, or whether the paper should present Phase C as exploratory.
6. The 8-cysteine rescue rule for hydrophobins (`docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md`).
   The pattern is not novel (docs/paper/05, sections 5 and 6). Updated 2026-10-08 after full-text reads and a
   wider search (PubMed, bioRxiv by web search, InterPro and Pfam entries; docs/paper/05 section 6.4). Found:
   Jensen 2010 (PMID 21182770) screens genomes with the eight-cysteine pattern plus a signal sequence and
   reports 5 of 50 hydrophobins with no Pfam hit. Li 2021 (PMID 33440688) takes the union of Pfam hits and
   class I and II cysteine-spacing signatures. Lovett 2022 (bioRxiv, not peer reviewed, doi
   10.1101/2022.08.19.504535) combines a 6-cysteine pattern, Pfam hmmsearch and SignalP, and reports 15
   candidates found only by the pattern, plus 996 other proteins with the pattern. Yang 2006 used the pattern
   as a filter, not a rescue. No source found that measures the sensitivity and false-positive cost of a rule
   "8-cysteine pattern plus signal peptide added to a Pfam call" against a truth set. This covers only the
   sources named. Not searched: Google Scholar, Scopus, and the full text of Wessels 1994, Linder 2005 and
   Sunde 2008 (not readable). Do not claim novelty of the pattern, of the pattern-plus-Pfam union or of the
   signal-peptide condition. Claim only measured sensitivity and false-positive cost. The published spacings
   disagree in places (docs/paper/05 section 6.3). The paper must state which spacing was used and why.
7. On `hydrophobin-validation` the hydrophobin Pfam models have their own module (`pfam_hydrophobin`) and a hit
   no longer sets `wall_family_domain`. PR #75, as pushed, still has them in `pfam_adhesion`, so a hydrophobin hit
   there sets `wall_family_domain` and can set `cell_wall_adhesion_candidate`. Report class counts with this in mind
   until the branch is merged.
8. Hydrophobin measurement (docs/reports/2026-10-08-hydrophobin-validation.md): every status is `smoke`, specificity
   rests on assumed negatives and is a lower bound of unknown size (19 of 22 unlabelled calls have 8 or more
   cysteines), and the published spacing recovers only 60% of Swiss-Prot hydrophobins. Do not report a hydrophobin
   specificity. Do not report the rescue as an improvement.
9. The hydrophobin call in the tool is the strict Pfam call only (`hydrophobin_domain`). The relaxed level (`hydrophobin_relaxed`) recovers
   most hydrophobins that the strict call misses, but in whole proteomes its extra calls could not be shown to be hydrophobins and the cost
   limit failed in 2 of 7 test proteomes. The owner decided on 2026-10-08 not to bring it into the tool. The curation behind this is tier T4
   (assistant, from evidence, not owner-reviewed). Do not describe the relaxed level as part of the tool's hydrophobin detection.
