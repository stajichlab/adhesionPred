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
| **2b. Small secreted proteins** — *see below; this was one row and is now three* | CalA (177 aa), Ag2/PRA, PRA3 | **not one architecture** | **no** — scores ~0.000 | **split 2026-09-30.** The single row was a residual bucket, not a mechanism class |
| **2b-i. CFEM-domain surface proteins** | Ag2/PRA, PRA2, ~7 per Onygenales genome | CFEM hemophore fold, structurally confirmed | **yes** — PF05730 | *HMM problem, not ML*. Same standing as 2c |
| **2b-ii. Small Cys-knot secreted proteins** | PRA3 | 38 aa 7-Cys knot on a disordered stalk; no fold assignment | **no** | **open structural question**, one protein. Not a classification task |
| **2b-iii. Bys1-domain invasins** | CalA | Bys1 domain (thaumatin-*like* fold) | **yes** — PF04681 | *HMM problem, not ML* |
| **2c. Hydrophobin/repellent attachment** | RodA, RodB, Ustilago Rep1 | 8-Cys hydrophobin pattern | **yes, trivially** — PF01185/PF06766 | *solved by HMMs*; was misfiled as an ML failure |
| **2d. Moonlighting surface proteins** | Histoplasma Hsp60, *Paracoccidioides* gp43 | cytoplasmic or enzymatic proteins on the surface | **no, by construction** | violates the stage-1 premise; must be held out, not predicted |

### Why class 2b was split (2026-09-30)

Source: `docs/reports/2026-09-29-class2b-structure.md`, `analysis/class2b_structure/`.
AlphaFold DB models, confident cores only (pLDDT ≥ 70), compared by TM-align and read together
with LDDT.

**The positive set is 2 proteins, or 6 drawn as generously as the evidence allows** — of 110
curated adhesins, 12 are ≤ 260 aa, of which 4 are hydrophobins (2c) and 4 are flocculin/PA14
fragments (2a). Those 6 already span five Pfam families. **This cannot support a classifier**:
no homology-grouped CV, no precision estimate, no train/test split. A fold survey was run
instead, and no classifier was built.

What the structures show:

| member | confident core | fold | evidence |
|---|---|---|---|
| Ag2/PRA, PRA2 | 64, 67 aa | **CFEM hemophore** | TM 0.80 / LDDT 0.77 to Csa2 and Rbt5; and TM 0.807 / LDDT 0.77 to the **experimental** Csa2 crystal `4Y7S` |
| PRA3 | **38 aa**, 7 Cys | none assigned | best score 0.477 (LDDT 0.47), below threshold. Rest of the protein is Pro/Thr/Glu at pLDDT 35–50 |
| CalA | 145 aa | **Bys1 / thaumatin-like** | TM 0.83 / LDDT 0.74 to the experimental *Magnaporthe* elicitor MoHrip2 `5FID`, then plant thaumatins |

Three consequences:

1. **The old row's name asserted a mechanism the structures do not support.** The only 2b fold
   assignable with confidence is a *hemophore* — an iron-acquisition architecture, not a
   host-receptor-binding one. Whether Ag2/PRA binds heme was not tested and is not claimed.
2. **Two of the three resulting classes are HMM problems** (PF05730, PF04681), the same
   conclusion already reached for hydrophobins. Only PRA3 is an open structural question, and
   it concerns one protein.
3. **Class 2a and class 2b differ in kind, not degree.** The 2a control SOWgp has **zero**
   residues above pLDDT 70. One class has a fold; the other does not.

Two methodological cautions, both of which fired in practice:

- **TM-score alone is not sufficient evidence for a short query.** Ag2/PRA scores TM 0.518
  against a 392 aa TIM barrel, but LDDT 0.41. Every genuine relationship here sits at LDDT
  0.76–0.77 and every artifact at 0.33–0.48. Report LDDT with TM.
- **A whole-PDB Foldseek search missed the real answer.** It returns no hit to `4Y7S` for
  Ag2/PRA in either mode, while direct pairwise TM-align with the prefilter disabled finds it
  at TM 0.807. The 3Di prefilter drops true neighbours of ~60-residue queries. For a class
  *defined* by being small, "no Foldseek hit" is not evidence of no fold relationship — run
  both, and do not read a null as negative.

**Naming caution.** PF04681 is `Bys1` (*Blastomyces* yeast-phase-specific protein, IPR006771).
It is **not** a thaumatin Pfam — thaumatin proper is PF00314. UniProt names `Q4WXJ1`
"Extracellular thaumatin domain protein, putative" and the structural neighbours are
thaumatin-like, so "thaumatin-like fold, Bys1 family" is the accurate phrasing. CalA's two
*A. fumigatus* paralogs, calB (`Q4WBB5`, AFUA_8G01710) and calC (`Q4WFZ4`, AFUA_3G00510), carry
the same domain and are **not** curated or tested; they are the obvious specificity control for
any PF04681 rule.

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
**No repeat-based tool will ever find these.** The 2b structure survey now gives this a
structural reading: these proteins have a folded domain and SOWgp has none, so they are not
two points on a size axis.

### What the repeat detector actually measures (2026-09-30)

Source: `analysis/cocci_repeats/REPORT_2026-09-29_repeat_detector_divergence.md`, scripts 14-18.

The class-2a detector was measured against a synthetic divergence series for the first time.
`02_repeat_profile.py` scores a period by **exact** residue matches, and that has a floor:

| | old (`02`) | new (`14`) |
|---|---|---|
| 50% recall reached at unit identity of about | **75%** | **35%** |
| copy-number error, mean | **−1.43** (grows with array length: −3.43 at 15 copies) | **+0.38** (roughly flat) |
| calls on the same 71,044 *Coccidioides* proteins | 83 | 189 |

**So the 58 proteins reported in the 2026-09-27 survey are not the *Coccidioides* proteins with
tandem repeats. They are the ones whose repeat units are more than about 75% identical.** That
is a statement about the detector, not about the biology, and every downstream count inherits
it.

Four cautions before anyone adopts `14`:

- **It is not a replacement.** Above ~85% unit identity `02` has the higher recall (96–99% vs a
  ~90% plateau), cleaner period calls, and its fractional count rounds to the anchored SOWgp
  truth more often (87.1% vs 24.3%).
- **It loses 4 Pro/Cys-rich curated class-2a candidates.** Pro/Gly-rich units match at *every*
  period, so the significance test that suppresses low-complexity false positives also rejects
  them. Unresolved, and it hits exactly the SOWgp/BAD1 composition class this section is about.
- **Most of the gain is not the substitution matrix.** It is replacing a flat 0.3 cut on the
  raw match rate with a significance test. Similarity scoring adds ~5 points and costs a
  27-fold rise in low-complexity false positives that must then be controlled.
- **34 of the 123 new calls are ankyrin repeats** — real repeats, but intracellular. A new
  false-positive class for 2a curation at the biology level, which stage 1 should filter.

The divergence axis is synthetic, single substitution model, **no indels**, so ~35% is an upper
bound. There is no real benchmark between 40–60% unit identity.

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
| 2a repeat *detector* | divergence floor measured 2026-09-30: `02` ≈ 75% unit identity, `14` ≈ 35%. Neither supersedes the other |
| 2c hydrophobin | solved by existing HMMs; no work needed |
| *Coccidioides* antigen | built, 3/4 calibration, 14 Tier-1 / 45 Tier-2 candidates |
| Biofilm | **not built** — blocked on phenotype linkage, not on modelling |
| 2b-i CFEM (Ag2/PRA, PRA2) | fold confirmed against experimental `4Y7S`; **HMM problem** (PF05730), not ML |
| 2b-ii Cys-knot (PRA3) | **open structural question**, one protein. No fold assignment |
| 2b-iii Bys1 (CalA) | fold confirmed against experimental `5FID`; **HMM problem** (PF04681), not ML |
| 2d moonlighting | **not achievable from sequence**; documented as out of scope |

## 6. Related documents

- `docs/reports/2026-09-27-coccidioides-antigen-findings.md` — the *Coccidioides* results
- `docs/model-review/2026-09-27-review-and-framework-plan.md` §4.7, §8 — the repeat-detector evidence and the mechanism-split argument
- `docs/model-review/STATUS.md` — per-clade status and what is validatable
- `docs/model-review/SEARCH-FRAMEWORK.md` — the Rhodotorula/Coccidioides search design and its independent review
- `analysis/curation/BASIDIOMYCETE_NOTES.md` — why Basidiomycota need their own tool
