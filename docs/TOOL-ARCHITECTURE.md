# Tool architecture: different predictors for different intentions

*2026-09-27. Consolidates what this project has actually built, after the evidence showed that
a single "fungal adhesin predictor" is the wrong target. Written in response to the framing:
we are building several tools with different intentions, not one tool.*

## 1. The core argument, with the measurement behind it

"Adhesin", "antigen", "biofilm factor" and "surface protein" are **outcomes**, not molecular
features. Proteins reach each outcome by physically different routes, and those routes have
different sequence signatures — so they need different predictors, different training sets,
different validation, and different scope statements.

The measurement that forced this: the existing classifier is functionally a **tandem-repeat
surface protein detector**. Across 95 curated adhesins, model score tracks 8-mer repeat
coverage (Spearman ρ = +0.302, p = 0.003); median repeat coverage is **1.000 for the adhesins
it finds and 0.000 for the ones it misses** (Mann-Whitney p = 0.0037). It is not detecting
adhesion. It is detecting repeats.

## 2. The architecture that follows

### Stage 1 — general surface/secreted protein (clade-general)

| | |
|---|---|
| question | is this protein displayed on the cell surface or secreted? |
| features | signal peptide, GPI anchor, cell-wall retention, TM topology |
| tools | SignalP 6.0, NetGPI/PredGPI — **not** a language model |
| scope | all fungi; these features are universal |
| status | **works**; currently keyword-derived, should be run directly |
| known gap | annotation coverage. Only 371 of ~9,910 *C. immitis* RS proteins carry a SignalP call (~4%, vs ~10% expected), and **SOWgp — a known surface antigen — has no call at all**. Stage 1 is the foundation of everything downstream and it is under-called. |

This is the one genuinely general tool, and it does not need ML.

### Stage 2 — split by mechanism, not by outcome

| class | exemplars | architecture | detectable from sequence? | status |
|---|---|---|---|---|
| **2a. Repeat/avidity surface proteins** | FLO11, ALS1, AGA1, SOWgp, BAD1, CspA | tandem repeats, repeat coverage ≈ 1.0 | **yes, easily** | **works**: PR-AUC 0.98 homology-grouped, 0.94 leave-genome-out |
| **2b. Small receptor-binding invasins** | CalA (177 aa), Ag2/PRA, PRA3 | compact fold, one binding interface | **no** — scores ~0.000 | not achievable from sequence; needs structure |
| **2c. Hydrophobin/repellent attachment** | RodA, RodB, Ustilago Rep1 | 8-Cys hydrophobin pattern | **yes, trivially** — PF01185/PF06766 | *solved by HMMs*; was misfiled as an ML failure |
| **2d. Moonlighting surface proteins** | Histoplasma Hsp60, *Paracoccidioides* gp43 | cytoplasmic or enzymatic proteins on the surface | **no, by construction** | violates the stage-1 premise; must be held out, not predicted |

### A finding that refines class 2a

Repeat coverage alone does **not** define one class. Composition splits it in two:

| protein | len | repeat cov. | %Ser+Thr | %Pro | %Cys |
|---|---|---|---|---|---|
| FLO11 (*S. cerevisiae*) | 1367 | 1.00 | **50.1** | 10.5 | 1.5 |
| ALS1 (*C. albicans*) | 1260 | 1.00 | **36.8** | 6.7 | 0.8 |
| AGA1 (*S. cerevisiae*) | 725 | 0.75 | **55.2** | 4.3 | 1.5 |
| **SOWgp (*Coccidioides*)** | 328–422 | **1.00** | 11.4 | **16.6** | **6.6** |
| **BAD1 (*Blastomyces*)** | 1146 | **1.00** | 8.5 | 3.8 | **8.0** |

**SOWgp and BAD1 are in the FLO11 repeat class but have inverted composition** — proline- and
cysteine-rich rather than Ser/Thr-rich. So:
- a *repeat* detector generalises across both (and should find SOWgp)
- a classifier trained on FLO/ALS *composition* will not

This is the mechanistic reason behind the clade-transfer collapse (ROC 0.60 cross-clade vs
0.94 size-matched within-clade): the Onygenales repeat antigens are the same structural class
with different chemistry.

And the small *Coccidioides* antigens are a different problem entirely — Ag2/PRA (repeat 0.25),
PRA3 (0.00), PRA2 (0.00), CF antigen (0.00) — alongside AGA2 (0.00), RodA (0.00), CalA (0.00).
**No repeat-based tool will ever find these.**

### Purpose-specific predictors (built on stage 1, not on stage 2)

| tool | question | positives | key discriminator | validation | status |
|---|---|---|---|---|---|
| **Antigen (per clade)** | would a host antibody response see this, and is it species-specific? | known antigens of that clade | **orthology vs confounder fungi** + phase-specific expression | known antigens as controls, with a hard acceptance test | **built for *Coccidioides*** (`analysis/cocci_antigens/`) |
| **Biofilm** | is this a biofilm factor? | curated biofilm genes | unknown — half of curated biofilm genes are regulators, not surface proteins | **genotype–phenotype association** | not built; needs phenotype (Rhodotorula SBF data exists, n=11) |
| **Adhesin (2a only)** | repeat/avidity adhesin? | 96 curated adhesins | repeat content + composition | homology-grouped and leave-genome-out CV | **works within Saccharomycotina**; scope must be stated |

**The antigen predictor is the clearest case that these are different tools.** Its
discriminating signal is not a sequence property at all — it is *absence of orthologs in
Histoplasma/Blastomyces/Paracoccidioides/Aspergillus*, which no single-sequence model can
compute. ML contributed the search space; comparative genomics produced the result.

## 3. What this means for *Coccidioides* antigens specifically

Building a classifier on sequence/structural expectations is the right instinct, and the
property profile of a good *Coccidioides* serodiagnostic target is now explicit:

| property | target value | why | source |
|---|---|---|---|
| signal peptide | present | must be secreted/surface to be seen | SignalP |
| ortholog in confounder fungi | **absent** | cross-reaction is the dominant clinical failure mode | orthology vs 4 genomes |
| human ortholog | absent | self-tolerance and assay background | UniProt |
| prevalence across isolates | ≥95% of 488 | must detect every isolate | pangenome |
| phase expression | spherule-induced *(class-dependent)* | the host meets spherules, not mycelia | RNA-seq |
| repeat content | high *(SOWgp class)* or low *(PRA class)* | **the two classes are different tools** | this document |

Note the last row. SOWgp and the PRA family are **not the same kind of target**: SOWgp is a
Pro/Cys-rich tandem-repeat protein, massively spherule-induced (13 → 15,000 TPM); the PRAs are
small, non-repetitive and *mycelia*-high. One classifier covering both would be the same
mistake as one classifier covering all adhesins.

## 4. Rules this project now follows

1. **Name the tool after what it detects, not what you want.** "Tandem-repeat surface protein
   detector", not "adhesin predictor".
2. **State the clade scope.** Adhesin families are phylogenetically disjoint: Flocculin is
   Saccharomycotina-only (330 proteins there, 0 elsewhere); CPL1-like is Basidiomycota-only
   (1,640 vs 0 in Ascomycota).
3. **Hold every predictor to an acceptance test on known positives**, and let it fail loudly.
   The *Coccidioides* scorer prints `NOT CALIBRATED` and currently does so at 3/4.
4. **Keep axes separate; do not merge into one score.** Merging antigenicity with specificity
   hid that the CF antigen is immunogenic *and* cross-reactive. Merging spherule expression in
   would have demoted PRA3.
5. **Validate against something independent of the labels** where labels are thin: phenotype,
   pangenome variability, or phase expression.
6. **Say when the right tool is not ML.** Hydrophobins are an HMM problem. Stage 1 is a
   SignalP problem. Antigen specificity is a comparative-genomics problem.

## 5. Status

| tool | state |
|---|---|
| Stage 1 surface/secreted | works; needs SignalP/NetGPI run directly (issue #15) |
| 2a repeat/avidity adhesin | works within Saccharomycotina; PR-AUC 0.94–0.98 |
| 2c hydrophobin | solved by existing HMMs; no work needed |
| *Coccidioides* antigen | built, 3/4 calibration, 14 Tier-1 / 45 Tier-2 candidates |
| Biofilm | **not built** — blocked on phenotype linkage, not on modelling |
| 2b invasin / 2d moonlighting | **not achievable from sequence**; documented as out of scope |

## 6. Related documents

- `docs/reports/2026-09-27-coccidioides-antigen-findings.md` — the *Coccidioides* results
- `docs/model-review/2026-09-27-review-and-framework-plan.md` §4.7, §8 — the repeat-detector evidence and the mechanism-split argument
- `docs/model-review/STATUS.md` — per-clade status and what is validatable
- `docs/model-review/SEARCH-FRAMEWORK.md` — the Rhodotorula/Coccidioides search design and its independent review
- `analysis/curation/BASIDIOMYCETE_NOTES.md` — why Basidiomycota need their own tool
