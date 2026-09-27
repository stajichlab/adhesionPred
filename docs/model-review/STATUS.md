# Where the tools stand — 2026-09-27

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

A reasonable next step is to split "adhesin" by **mechanism class** — repeat-rich
GPI/cell-wall adhesins, hydrophobin-type, receptor-binding invasins, moonlighting — and model
the first class well rather than pretending one classifier covers all four.
