# Plan: a multi-step protein classification pipeline, and the decisions behind it

*2026-09-30. Summary of the review session and the owner's decisions. Written to be read and shared
without the conversation. Status: decisions made, no code changed yet for them.*

## 1. Where this came from

The shipped `adhesion_predict` tool is not an adhesin predictor. Evidence:

- On *S. cerevisiae* S288C about 12% of its calls are known adhesins
  (`docs/model-review/2026-09-27-review-and-framework-plan.md`).
- Its score tracks tandem-repeat content (Spearman 0.302, p = 0.003 across 95 curated adhesins).
- It was trained on FLO/ALS homologs versus random proteins, so it learned "FLO/ALS-like surface
  glycoprotein".

So the project is split into steps, each with a name that matches what it measures. The overall
goal stays the same: take a proteome and sort its proteins into categories.

## 2. The flow of one protein

```
                      PROTEOME (FASTA)
                            |
             +--------------v---------------+
             | STEP 1  surface / secreted?  |   rule: SignalP + GPI + Ser/Thr
             +------+----------------+------+
                    |                |
                 yes|                |no
                    |                +--------> NOT SURFACE
                    |                           (cytosolic, membrane, organellar)
     +--------------v-----------------------------------------+
     | STEP 2  which adhesion mechanism? (checks run together) |
     |                                                         |
     |  2a  repeat / avidity    tandem-repeat detector         |
     |  2b  small secreted      HMM scans:                     |
     |        CFEM  PF05730   Bys1 PF04681                     |
     |        Cys-knot (PRA3): not a classifiable class        |
     |  2c  hydrophobin         HMM scans PF01185 / PF06766    |
     +--------------+------------------------------------------+
                    |
          matched >= 1 class?  -- yes --> ADHESION CLASS(ES)
                    |
                    no --------------> SURFACE, OTHER
                    |
     +--------------v-----------------------------------------+
     | STEP 3  purpose layers (any surface protein)            |
     |  antigen   orthology vs confounder fungi, prevalence,   |
     |            phase expression (Coccidioides only today)   |
     |  biofilm   not built (needs phenotype data)             |
     +---------------------------------------------------------+

 2d moonlighting (Hsp60, gp43): cytosolic or enzymatic, so they fail step 1.
    Not predictable from sequence; reported as out of scope.
```

Points that are easy to get wrong:

- **Antigen is a purpose layer (step 3), not a sibling of adhesin.** It runs on surface proteins and
  can co-occur with any step 2 class. The current antigen tool is *Coccidioides*-specific and is a
  comparative-genomics ranking. It is not epitope prediction. No B-cell epitope predictor exists here.
- **Categories overlap.** One protein can be surface, class 2a and an antigen candidate. The output is
  a table with one column per tool, not one category per protein.
- **"Other" has two meanings.** "Not surface" and "surface, no mechanism matched" are different and
  must stay separate in the output.
- **Step 2 is not one classifier.** 2a is a detector. CFEM, Bys1 and hydrophobin are HMM scans. The
  Cys-knot (PRA3) is one protein and cannot be a tool. Only 2a has validation today.

## 3. Status of each piece

| Piece | State |
|---|---|
| Step 1: SignalP + GPI + Ser/Thr rule | decided, to build. A direct SignalP 6 run calls 460 of 9,910 *C. immitis* RS proteins (4.6%). It calls SOWgp (CIMG_04613) with probability 0.9998. No truth data show whether 4.6% is an under-call. See `docs/TOOL-ARCHITECTURE.md`, Stage 1 |
| 2a repeat detector | works. Two detectors; `02` reaches 50% recall at ~75% unit identity, `14` at ~35%. Neither supersedes the other. Validated only in Saccharomycotina |
| 2b-i CFEM (PF05730), 2b-iii Bys1 (PF04681), 2c hydrophobin | HMMs exist. A scan wrapper and a specificity test are missing (for Bys1: the *A. fumigatus* paralogs calB and calC) |
| 2b-ii Cys-knot (PRA3) | one protein; not a tool |
| 2d moonlighting | out of scope by construction |
| Step 3 antigen | built for *Coccidioides* only; calibration 3 of 4 |
| Step 3 biofilm | not built; blocked on phenotype data |
| Legacy ESM + logistic regression | frozen; renamed `surface_glyco` (section 4) |

## 4. Decisions (owner, 2026-09-30)

| # | Decision | Reason |
|---|---|---|
| 1 | Step 1 is a **rule**: SignalP + GPI anchor + Ser/Thr content. Not an ML model. | The review concluded stage 1 "does not need ML". The ESM model's labels (FLO/ALS vs random) do not define surface glycoproteins. |
| 2 | The **ESM + LR model is frozen as legacy**. Safety fixes only: card, explicit legacy-pooling flag, refuse on mismatch. No retraining, no tuning. | Step 1 measurement: the scaler barely matters for ESM features (section 5). Keep the embedding code for a possible step 2. |
| 3 | Step 1 truth set: **curated GO annotation** (SGD/CGD cell wall and GPI-anchored vs cytosolic/nuclear) in S288C and *C. albicans*, plus the curated Onygenales/Eurotiales literature rows. Report per clade. | SignalP cannot check itself. `surface.tsv` labels come from UniProt keywords and share inputs with the rule. The literature rows cover the known Onygenales under-call. |
| 4 | **Breaking rename now, no aliases.** Package `surface_glyco`, holding only the legacy code. Entry points `surface_glyco_predict`, `surface_glyco_train`, `surface_glyco_evaluate`. Column `surface_glycoprotein_score`. Labels `surface_glycoprotein` and `other`. Model files renamed to match. Repo name `adhesionPred` stays. | The old names assert something the tool does not do. |
| 5 | New tools (step 1 rule, step 2) get **their own packages or modules**. | Each name then matches what it does. |
| 6 | Accuracy gates: **measure first, then freeze as regression tests.** Report recall, precision and false-positive rate per clade with confidence intervals; the owner reviews; gates are set at or just below the measured values. | No target is guessed in advance. |

Defaults adopted for the remaining questions in the Fable review (change on request):

1. The legacy threshold stays at 0.5 and is documented as uncalibrated.
2. Python 3.12 is the canonical version. CI uses CPU only.
3. Both packaged pickles (8M and 35M) are tracked in one place, `src/surface_glyco/models/`. The duplicate root `models/` copies are dropped.
4. A sequence present in both classes is dropped from training. This only matters if the legacy model is retrained.
5. Curated labels may be smoke-test fixtures, not accuracy gates, until expert review.
6. Old Fungi_5k outputs are frozen and labelled "legacy pooling, 35M".

## 4a. Changes made later on 2026-09-30

| Was | Now | Reason |
|---|---|---|
| ESM + LR model frozen as legacy; pickles kept, legacy pooling emulated behind a flag (decision 2) | The old pickles are **deleted**. No legacy mode, no `--allow-legacy-pooling`. | The tool is not in circulation. The emulation of the original padding-inclusive pooling could only be approximate and could not be verified. |
| Old result files stay readable by downstream code | Old-schema files are **not** read; result discovery uses `*.surface_glyco.csv`. | Not in circulation; the Fungi_5k outputs came from the pre-fix model and are already due for a re-run (#16). |
| Step 1 is a SignalP + GPI + Ser/Thr rule (decision 1); ESM model only a possible step 2 base | Step 1 is **two candidates**: the rule, and an ESM-based model **retrained on a surface-glycoprotein label** (surface vs non-surface). Both are tested on the same GO truth set (decision 3); rule, ML or hybrid is decided from that comparison. | Owner is open to retraining on the correct target. This also measures the rule-vs-ESM agreement (R3). |
| 0.2.0 ships with the rename | **0.2.0 is cut when a validated model ships.** The rename and card framework merge to `main` earlier, unreleased. | The package has no useful default model until the new one exists. |
| Simple retrain on the current labels in 0.2.0 | **Skipped.** | Replaced by the surface-glycoprotein-label model. |

The surface-glycoprotein-label model needs its own spec (positive and negative definitions,
homology-grouped CV, GO truth set, head-to-head with the rule) and independent review before
any plan or code.

## 5. What was measured (step 1 job 29301182, 8M model)

Records: `analysis/model_review/results/step1/`, issue #25.

- The shipped pickles (committed 2026-02-14) were trained on pooling that averaged padding tokens.
  `predict` has used residue-only pooling since 2026-09-27. "Legacy pooling" in the job is an
  emulation (batch 4, file order); the real training batches were not recorded.
- On S288C (6,722 proteins): 77 calls under emulated legacy pooling, 104 under residue-only. All 77
  remain; 27 are new and unchecked. Spearman rho 0.973, median score difference 0.0002, maximum 0.716.
- Scaler on versus off: ESM features change by 0.001 ROC-AUC or less. Composition + length features
  lose more without scaling (leave-family-out recall 0.945 vs 0.801).
- Caveats: 8M model only; the training embeddings include 167 duplicate sequences; two values
  differ slightly from the 2026-09-27 review table (for example 0.818 vs 0.796), cause not checked.

## 6. The high-level interface (proposal, not built)

The owner wants one entry point that takes a proteome and sorts it, while the sub-tools stay separate.

- A thin **orchestrator** takes a proteome and runs the registered tools.
- It writes **one row per protein** with a column per tool, for example `surface_call`, `mech_2a`,
  `mech_cfem`, `mech_hydrophobin`.
- Each tool returns a **call, its evidence** (domain hit, repeat period) and a **validation status**,
  so an unvalidated call cannot look like a validated one.
- Each tool reports its **clade scope**. The 2a detector is validated only in Saccharomycotina; the
  antigen tool only in *Coccidioides*.
- Every tool also runs alone. A new tool is another plugin with the same small interface.

This is new code, so it needs a design spec and an independent review before any implementation.

## 7. Work, in order

1. **Legacy rename and safety fixes.** In progress (branch `surface-glyco-rename`). Touches 46 files: package, tests, `kingdom_survey`, analysis
   scripts, HPCC scripts, docs. The CHANGELOG entry exists. No version bump until a validated model ships (section 4a). Existing Fungi_5k CSVs are
   labelled legacy. Shrinks the #9 plan to the card, the pooling flag, refusal on mismatch, and tests.
2. **GO truth-set extraction** for step 1 (S288C, *C. albicans*), plus the literature rows.
3. **Step 1 rule and its baseline measurement**, per clade with confidence intervals. Then freeze the gates.
4. **Orchestrator design spec** and independent review.
5. **Step 2 scan wrappers** (CFEM, Bys1, hydrophobin) with specificity tests; decide how 2a's two
   detectors are exposed.
6. Later: antigen beyond *Coccidioides* (needs confounder genomes); biofilm (needs phenotype data).

Testing framework recommended by the Fable review: fast unit tests in `tests/` run in CI; slow
integration tests with golden scores; a tier harness under `analysis/` that runs on HPCC and writes
JSON the model card reads. Accuracy results enter a card only through that harness.

## 8. Open items and risks

- The release flow (`version_bump`, then `publish_release`) is untested.
- Not measured: agreement between the ESM+LR score and a SignalP + GPI call; GPU versus CPU
  embedding differences.
- The legacy-pooling emulation cannot be shown equal to the original input.
- Curated labels are a draft (`needs_review=yes` on nearly all rows; one row has no accession).
- `config.ESM2_MODELS` lists models that do not exist in fair-esm; delete it.
- Untracked scratch data is logged in issue #26.

## 9. Where things are

| Item | Location |
|---|---|
| This plan | `docs/PLAN-2026-09-30-pipeline-and-decisions.md` |
| Steps and names (implemented, unreleased) | `docs/TOOL-ARCHITECTURE.md` section 2.0 |
| #9 plan, with Opus review (section 7) | `docs/plans/2026-09-30-model-bundle-issue-9.md` (branch `issue-9-plan`) |
| Independent design and test review (Fable) | `docs/plans/2026-09-30-design-review-fable.md` (branch `issue-9-plan`) |
| Step 1 measurement records | `analysis/model_review/results/step1/` (branch `issue-9-plan`) |
| Pooling issue | #25 |
| Repo hygiene | #26 |
| CI and release workflow fixes | PR #27 |
| `cocci_repeats`, class 2b, naming docs | PR #28 |
