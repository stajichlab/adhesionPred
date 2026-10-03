# Where the tools stand — 2026-10-02

*Rewritten 2026-10-02. The 2026-09-27 text is kept at the end as history. It describes the
tool before the rename and before step 1 was measured.*

Short answer: **the shipped CLI is step 1 (surface glycoprotein prediction), not an adhesin
predictor. No trained model ships. Step 1 has now been measured, but the choice between a rule
and an ML model is not made. Step 2 (adhesion mechanism classes) is mostly unbuilt.**

Plan and decisions: `docs/PLAN-2026-09-30-pipeline-and-decisions.md`.
Spec: `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md`.

## 1. What exists

| component | state | where |
|---|---|---|
| `surface_glyco` package (ESM-2 + logistic regression, legacy code only) | renamed from `adhesion_predict` (PR #31). The two old pickles are deleted. `predict` and `evaluate` read the model, layer and pooling from a model card and refuse on mismatch | `src/surface_glyco/` on `main` |
| Residue-only pooling, no silent batch drops, all scores written, duplicate positives dropped | on `main` (PRs #18, #23, #24) | |
| CI: lint (ruff 0.3.5) and unit tests | on `main` (PR #27). The Lint job passes on PR #36 after `zip(strict=True)` | `.github/workflows/` |
| Step 1 truth set (Phase A), features and embeddings (Phase B), rule-versus-ML evaluation (Phase C) | built, run on HPCC, and on `main` since PR #38 (2026-10-02). The reported numbers come from the run with the widened C grid | `analysis/step1_compare/`, outputs in `_workdir/step1_compare/` (git-ignored) |
| Step 2a repeat detectors | two detectors, validated only in Saccharomycotina | `analysis/cocci_repeats/` |
| Step 2b/2c HMM scans (CFEM, Bys1, hydrophobin) | HMMs exist. No scan wrapper, no specificity test | |
| Step 3 antigen layer | *Coccidioides* only | `analysis/cocci_antigens/` |
| Cysteine-rich secreted candidates (PRA3-like) | on `main` (PR #36). Full-length PRA3 structure re-run found no fold (PR #42, merged) | `analysis/cys_candidates/` |
| PF28404 (ARB_05178) family | search of 831 proteomes and a 217-protein tree (PR #44, merged). Not specific to *Coccidioides*; four paralog groups older than the genus. Function unknown | `analysis/pf28404_family/`, `docs/reports/2026-10-02-pf28404-family.md` |
| Orchestrator (one table, one column per tool) | proposal only, no spec | plan §6 |
| Stage-2 adhesin classifier prototype (ESM C 300M) | prototype, not packaged | `analysis/model_review/stage2_proof_of_concept.py` |

No model ships. Version 0.2.0 is cut only after a validated surface-glycoprotein model exists.

## 2. Step 1 measurement (Phase C, run 2026-10-01)

Report: `_workdir/step1_compare/phasec/report.md`. Intervals are 95% cluster bootstrap
(2,000 resamples). Truth is direct GO evidence. Candidates: B0 and B1 (baselines), R0, R1, R2
(rules), M8, M35, M8-C, M35-C and H. M8 and M35 are logistic regressions on ESM-2 8M and 35M
embeddings of the N-terminal window. The -C variants use the C-terminal window for proteins over
1,022 aa. H is a hybrid: embedding, SignalP probability, GPI score and Ser+Thr in one regression.
B0 uses log length only. B1 uses amino-acid fractions and log length. R0 is "SignalP calls a signal
peptide". R1 is R0 plus a GPI call. R2 is "signal peptide and (GPI call or Ser+Thr fraction at or above t)".
A test set is an *estimate* when the recall half-width is at most 0.10 and it has at least
20 direct positives. Otherwise it is a *smoke test*.

The S1 set is homology-grouped cross-validation. The Eurotiomycetes and Basidiomycota rows come
from models that did not train on that clade (leave one clade out).

| Test set | Label | Positives / negatives | R0 recall / FPR | R2 recall / FPR | B1 recall / FPR | M8 recall / FPR |
|---|---|---|---|---|---|---|
| S1:all (S288C + *C. albicans*, cross-validated) | estimate | 232 / 4,244 | 0.603 / 0.037 | 0.418 / 0.006 | 0.763 / 0.130 | 0.759 / 0.088 |
| Eurotiomycetes clade (*A. fumigatus* + *A. nidulans*) | estimate | 128 / 208 | 0.727 / 0.010 | 0.227 / 0.005 | 0.859 / 0.168 | 0.914 / 0.058 |
| *A. nidulans* alone | estimate | 109 / 164 | 0.688 / 0.012 | 0.165 / 0.006 | not read | 0.899 / 0.037 |
| Basidiomycota (*Cryptococcus*, *Ustilago*) | smoke test | 16 / 60 | 0.938 / 0.083 | 0.125 / 0.017 | 0.625 / 0.283 | 0.875 / 0.167 |

Findings stored in `findings.json`:
- **(a)** B1 is not saturated: ROC-AUC 0.894 on S1:all, 0.736 to 0.957 on the three S2 sets.
  The comparison is not trivially won by a baseline.
- **(b)** No ML candidate beats both B1 and R2 on the false-positive rate for non-secreted
  proteins at the recall of R2 (S1:all). M8 beats B1 only.
- **(c)** The same holds on each of the three S2 sets (train on one species, test on another).

What this means, in plain terms:
- The rule R2 is precise and misses more than half of the surface proteins in S1:all, and about
  three quarters in Eurotiomycetes. R0 recalls more at a higher FPR.
- The ML models recall more than R2, with an FPR of 6 to 9% in the estimate sets. Within the
  confidence intervals, they are not shown to be better than R2 on non-secreted false positives.
- The choice of rule, ML or hybrid is the owner's. No gate values are set.

Literature rows (19 Onygenales and Eurotiales adhesins, positives only, smoke test): recall 1.00
for R0, B1, M8, M8-C and H; 0.947 for M35 and M35-C; 0.684 for R2; 0.632 for R1. SOWgp, CspA
and CBP1 are called by nearly all candidates. R1 and R2 miss CBP1, CTS1, abr2 and aspf2. HSP60
is missed by all candidates, as expected, since it has no signal peptide.

In the named panel, MSB2, HKR1 and SAP9 are labelled negative and most candidates call them
positive. This was not analysed.

## 3. What the data cannot show

From the report:
- Basidiomycota performance. *Cryptococcus* (7 positives) and *Ustilago* (9) are smoke tests.
  Curated truth does not exist.
- Performance on GPI-anchored proteins. `curated_gpi.tsv` has no literature rows.
- Whether SignalP under-calls in *C. immitis*. No truth set tests it.
- Whole-proteome prevalence of surface proteins. The prevalence table uses assumed values.
- Variance from refitting the models. The intervals describe test-set sampling only.
- Pooled S1:all ROC-AUC mixes the score scales of five fold models. The ML threshold and
  Platt scaling are fitted on out-of-fold predictions and applied to a refitted model. The
  effect was not measured.

Also, from the fitted-settings table: the C value of every ML candidate is at the lower edge of
the grid (0.001). Per the owner's decision of 2026-10-02 the grid was widened to 0.0001 and
0.0003. The report on disk was written before that change, so the ML numbers above come from
the narrower grid. They have not been re-run.

## 4. Remaining work

1. ~~Get the Phase A–C code onto `main`.~~ Done, PR #38.
2. **Owner decisions:** rule, ML or hybrid for step 1; whether to change the Youden criterion
   (R2 cannot beat R0 on it); gate values after review; whether to fund curation of GPI and
   Basidiomycota truth.
3. Re-run Phase C with the widened C grid, if the owner wants ML numbers from the wider grid.
4. Train and validate the chosen step 1 model, then ship it with a card (0.2.0).
5. Orchestrator design spec and independent review.
6. Step 2 scan wrappers (CFEM, Bys1, hydrophobin) with specificity tests.
7. Step 3: antigen beyond *Coccidioides*; biofilm (blocked on phenotype data).

Open issues: #9 to #17, #19, #25 and #26. See the issue tracker for the current order.

---

# History: where the tools stood on 2026-09-27

*This section is the earlier text, kept unchanged. It uses the old names (`adhesion_predict`,
"stage 1", "stage 2") and the old shipped models, which no longer exist.*

*Updated after the Onygenales/Eurotiales curation (§7). The label gap is partly closed;
the measurement it enabled is the important part.*

Short answer to "what can this predict today, and how well in Onygenales and *Aspergillus*":
**stage 1 works broadly, stage 2 works only in Saccharomycotina, and in Eurotiomycetes
(*Aspergillus* + Onygenales) the tool is effectively unvalidated.** Numbers below.

## 1. What exists

| component | state | where |
|---|---|---|
| `adhesion_predict` CLI (ESM-2 + logistic regression) | **shipped, in use** — but it is a *surface glycoprotein* detector, ~12% adhesin precision on S288C | `src/adhesion_predict/` |
| Embedding bugs (padding in mean-pool, silent batch drops) | **fixed**, not yet merged to `main` | PR #18 |
| Stage-2 adhesin classifier (ESM C 300M) | **prototype only** — not packaged, no CLI, no model file | `analysis/model_review/stage2_proof_of_concept.py` |
| Curated labels | 75 adhesins / 63 hard negatives, ~87% Saccharomycotina | `data/curated/` |
| Kingdom screen (5,802 species) | done, but reads as *surface glycoprotein* counts, and used the pre-fix embedding code | `analysis/kingdom_survey/` |

So there is no deployable adhesin predictor today. There is a deployable surface-protein
predictor, and a stage-2 prototype that works in one clade.

## 2. Label coverage for the clades in question

From `data/curated/surface/surface.tsv`:

| genome | surface proteins | adhesin | non-adhesin | labeled |
|---|---|---|---|---|
| *C. albicans* SC5314 | 413 | 26 | 21 | 11% |
| *N. glabratus* CBS138 | 270 | 13 | 8 | 8% |
| S288C | 350 | 10 | 29 | 11% |
| **A. fumigatus Af293** | **895** | **2** | **5** | **0.8%** |
| ***C. immitis* RS** | **557** | **0** | **0** | **0%** |
| ***C. posadasii* C735** | **512** | **0** | **0** | **0%** |

*Aspergillus* has 7 labeled proteins out of 895. Onygenales has **none**. Nothing can be
validated in Onygenales at all — not precision, not recall, not calibration.

## 3. What a Saccharomycotina-trained stage-2 model does in these genomes

Trained on 62 positives / 58 negatives from the four Saccharomycotina genomes, applied cold:

| genome | surface proteins | called adhesin | rate |
|---|---|---|---|
| *A. fumigatus* Af293 | 895 | 43 | 4.8% |
| *C. immitis* RS | 557 | 28 | 5.0% |
| *C. posadasii* C735 | 512 | 27 | 5.3% |

The rates look plausible, which is the trap: **in Onygenales there is no ground truth to
check any of those 27–28 calls against.** They are unvalidated output, not predictions with
known error rates.

Where ground truth does exist (*A. fumigatus*, n=7), it recovers 1 of 2 adhesins and
correctly rejects 5 of 5 non-adhesins — too small to mean anything.

## 4. Protein-level reality check

Scoring the characterized surface proteins of these clades:

| protein | organism | length | GPI | score | result |
|---|---|---|---|---|---|
| CspA (repeat-rich GPI cell-wall protein) | *A. fumigatus* | 430 | no | **0.872** | called |
| CalA (invasin, binds integrin α5β1) | *A. fumigatus* | 177 | no | **0.000** | **missed** |
| Ag2/PRA (proline-rich antigen) | *Coccidioides* | 194 | yes | **0.000** | **missed** |
| Gel1 (β-1,3-glucanosyltransferase) | *Coccidioides* | 447 | yes | 0.000 | missed (an enzyme — arguably correct) |

This is the clade-specificity result playing out one protein at a time. **CspA is found
because it looks like a yeast adhesin** — repeat-rich, cell-wall associated. **CalA is
invisible** because it is a 177-aa thaumatin-like invasin with no resemblance to FLO/ALS, and
Ag2/PRA is invisible because it is short where the model expects long. The model has learned
"long, Ser/Thr-rich, GPI-anchored", and Pezizomycotina/Onygenales adhesins frequently are not
that.

## 5. Gaps, ranked by how much they block progress

1. **No Onygenales labels at all.** This is the binding constraint. Nothing about
   *Coccidioides*, *Histoplasma* or *Blastomyces* can be validated. Needed: a curated set from
   BAD1, SOWgp, Ag2/PRA, the Hyr/Iff-like families, plus non-adhesive Onygenales cell-wall
   proteins as negatives. Realistically 20–40 proteins with PMIDs would change the picture.
2. **Only 7 labeled proteins in *Aspergillus*.** CalA, CspA, and the CFEM/hydrophobin
   families are the obvious seeds; several are already in `literature_seeds.tsv` but unresolved
   or E2/E3.
3. **Architectural blind spot.** Short (<250 aa), non-GPI, cysteine-rich surface adhesins
   score ~0 regardless of clade. Both CalA and Ag2/PRA fail this way, as does the whole
   basidiomycete CPL1-like family. Any fix has to include short cysteine-rich proteins in
   training, not just more of the same architecture.
4. **Stage 2 is not packaged.** No CLI, no saved model, no calibration. It exists as an
   evaluation script.
5. **Stage 1 is keyword-based.** The surface population comes from UniProt keywords, which
   are largely automatic for these genomes. Running SignalP 6.0 + NetGPI directly would give
   a better and genome-independent stage-1 set (issue #15).
6. **The kingdom screen predates the bug fixes** and its `adhesion_fraction` is a surface
   glycoprotein fraction. It needs re-running after the model is retrained (#16).

## 6. What would move Eurotiomycetes fastest

In order:
1. Curate ~30 Onygenales + ~20 *Aspergillus* adhesins and an equal number of characterized
   non-adhesive cell-wall proteins. This is a literature task, not a compute task, and it
   unblocks everything else.
2. Retrain stage 2 per clade and report per-clade performance honestly, rather than one
   kingdom-wide number (§4.5 of the review).
3. Add a short-protein / cysteine-rich axis to the feature set and check whether CalA-type
   invasins become findable.
4. Only then re-screen the kingdom.

Until step 1 exists, any *Coccidioides* or *Aspergillus* adhesin call from this tool should be
labelled explicitly as **a hypothesis with unknown error rate**.

*Reproduce the tables here with `analysis/model_review/clade_status_report.py`.*


---

## 7. After curating Onygenales + Eurotiales (added 2026-09-27)

`data/curated/adhesins/eurotiomycetes_seeds.tsv` adds 21 literature-curated rows (10
Onygenales, 11 Eurotiales). Trainable Eurotiomycetes labels went **7 → 22** (8 adhesin,
14 non-adhesin). Onygenales went from zero labels to seven.

**These clades can now be measured at all** — that is the main gain. The measurement:

| training set | ROC | PR |
|---|---|---|
| Saccharomycotina only, applied to the 22 Eurotiomycetes labels | 0.732 | 0.763 |
| + the Eurotiomycetes labels themselves (leave-one-out) | 0.670 | 0.760 |

**Adding the labels did not help.** With 22 labels against 145 from Saccharomycotina, they
are swamped, and leave-one-out on a 22-protein set is noisy. But the per-protein result shows
something more specific than "not enough data":

| found (p > 0.9) | missed (p ≈ 0.000) |
|---|---|
| BAD1 (Blastomyces, tandem repeat) | **rodA** (hydrophobin, ~16 kDa, Cys-rich) |
| SOWgp58/66/82 (Coccidioides, Pro-rich tandem repeat) | **CalA** (177 aa thaumatin-like invasin) |
| CspA (A. fumigatus, repeat-rich CWP) | **gp43** (Paracoccidioides, glucanase that binds laminin) |

5 of 8 adhesins are recovered and 13 of 14 hard negatives are correctly rejected (Cbp1 at
0.475 is the one near-miss). The three failures are **three different mechanism classes**,
none of which resembles FLO/ALS: a hydrophobin that works by surface hydrophobicity, a small
receptor-binding invasin, and a moonlighting enzyme.

**The conclusion is architectural, not quantitative.** The model finds repeat-rich surface
adhesins in *any* clade — SOWgp and BAD1 score ~1.0 despite being Onygenales, because they
look like what it was trained on. It is blind to hydrophobins, small invasins and
moonlighting enzymes *regardless of clade or training data*. More labels of the same
architecture will not fix that.

So the ranked gaps in §5 change: **#3 (the architectural blind spot) is now the binding
constraint, not #1 (labels)**. Concretely:
- Adhesion by surface hydrophobicity (hydrophobins, repellents) is a different physical
  mechanism and probably needs its own model, not more training data.
- Small receptor-binding invasins (CalA) may be unreachable from sequence alone without
  structure.
- Moonlighting adhesins (gp43, Hsp60) violate the stage-1 premise entirely: Hsp60 has no
  signal peptide, and gp43 is an annotated enzyme. Both are excluded by construction from a
  secretion-defined surface population, which is a design decision worth revisiting.

**Followed up in §4.7 of the review**: the split is driven by tandem-repeat content
(repeat coverage 1.000 for found adhesins vs 0.000 for missed, p=0.0037), so the classifier is
functionally a tandem-repeat surface protein detector. The design consequences are in §8 of
the review. A reasonable next step is to split "adhesin" by **mechanism class** — repeat-rich
GPI/cell-wall adhesins, hydrophobin-type, receptor-binding invasins, moonlighting — and model
the first class well rather than pretending one classifier covers all four.
