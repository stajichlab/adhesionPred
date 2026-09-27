# Is the CBM18 expansion in *Batrachochytrium* what the classifier calls "adhesion"?

*2026-09-27. Prompted by the question of whether CBM18 genes relate to adhesion or to immune
"cloaking". Data: Fungi_5k `function.duckdb` (Pfam, SignalP) joined to the existing
`adhesion_predict` result files. Scripts: `01_cbm18_vs_calls.py`, `02_cbm18_architecture.py`.
Per-protein table: `cbm18_batrachochytrium.tsv`.*

## What the classifier does with CBM18

CBM18 is `Chitin_bind_1` (PF00187) in Pfam, **not** a CAZy-overview family in this database —
querying `cazy_overview` for "CBM18" returns nothing, which is a trap worth recording.

| genome | proteome | adhesion calls | CBM18 proteins | CBM18 called | odds ratio | Fisher *p* |
|---|---|---|---|---|---|---|
| *B. dendrobatidis* JAM81 | 7,021 | 97 (1.38%) | 10 | **5 (50%)** | 75.2 | 1.1e-07 |
| *B. salamandrivorans* AMFP13 | 16,260 | 331 (2.04%) | 13 | **4 (31%)** | 21.6 | 1.0e-04 |
| *Homolaphlyctis polyrhiza* | 6,632 | 33 (0.50%) | 1 | 0 | — | 1 |
| *Spizellomyces punctatus* | 8,319 | 55 (0.66%) | 3 | 0 | — | 1 |

So CBM18 proteins are called at 15–36x the genome background in both *Batrachochytrium*
species, and not at all in the two non-pathogenic chytrid outgroups (which have few CBM18
genes to begin with).

## The architecture split is the interesting part

| species | architecture | n | signal peptide | called |
|---|---|---|---|---|
| Bd | CBM18 only | 4 | 2 | **4/4** |
| Bd | CBM18 + chitin deacetylase (`Polysacc_deac_1`) | 5 | 5 | 1/5 |
| Bd | CBM18 + tyrosinase | 1 | 1 | 0/1 |
| Bsal | CBM18 only | 5 | 3 | **3/5** |
| Bsal | CBM18 + chitin deacetylase | 7 | 5 | 1/7 |
| Bsal | CBM18 + tyrosinase | 1 | 0 | 0/1 |

These are the three categories described in Abramyan & Stajich 2012 (CBM18 alone, CBM18 with
a chitin-active enzyme, CBM18 with tyrosinase). The classifier picks out the **binding-only**
class almost perfectly and largely ignores the ones carrying a catalytic partner domain.

Copy number and length vary as expected for this family: 1–6 CBM18 repeats per protein,
161–1,197 aa.

## How to read this

**This is the model's opinion, not evidence of adhesion.** The same model calls the whole
Ser/Thr-rich surface-glycoprotein class at ~12% adhesin precision (see
`docs/model-review/`). A repeat-rich, signal-peptide-bearing, non-catalytic cell-surface
protein is exactly what it flags, whatever that protein actually binds. The enrichment is
real and specific to the two pathogens, but it is a hypothesis generator.

**Adhesion to host via CBM18 is mechanistically awkward**: amphibian skin is keratinized and
contains no chitin. A hevein-like chitin-binding module is better placed to act on fungal
cell-wall chitin than on host tissue. Two readings fit better than direct adhesion:
- **self-masking / "cloaking"**: shielding cell-wall chitin (or sequestering released chitin
  oligomers) from host chitinases and pattern-recognition receptors, by analogy with LysM
  (CBM50) effectors such as *Cladosporium fulvum* Ecp6 in plant pathogens. Different module,
  same logic.
- **binding host GlcNAc-bearing glycoconjugates** (e.g. skin mucins), which would blur the
  line between cloaking and adhesion.

Distinguishing these needs binding assays and knockouts, not sequence analysis.

## Caveat on the count

Only 10–13 CBM18 proteins are recovered here, which is **an undercount** relative to the
expansion originally described. This analysis uses the Fungi_5k annotation of JAM81 with
Pfam gathering thresholds; diverged repeats and copies below threshold will be missed, and
the original work used a different assembly plus manual HMM searching. Treat the counts as a
floor. A dedicated HMM search over both genomes would give the real family size.

## The published adhesion-candidate claims do not survive contact with this data

Van Rooij et al. 2015 (*Vet Res* 46:137) lists vinculin, fibronectin and fasciclin as
*B. dendrobatidis* adhesion-related genes, citing transcriptomic work (Rosenblum et al. 2008
PNAS; 2012 *Mol Ecol*). In this annotation:

| domain | Bd copies | called adhesion | Bsal copies | called |
|---|---|---|---|---|
| vinculin | 1 | 0 | 1 | 0 |
| fibronectin (fn3) | 1 | 0 | 3 | 0 |
| fasciclin | 3 | 0 | 2 | 0 |
| lectin (any) | 4 | 0 | 6 | 0 |

None are called, and vinculin in particular is a cytoskeletal protein, not a surface adhesin;
fn3 is a generic Ig-like fold. These look like domain-name transfers from animal cell-adhesion
vocabulary rather than identified fungal adhesins. The review itself concedes "the exact
factors mediating adhesion remain uncertain". **There is still no molecularly identified
adhesin in either *Batrachochytrium*.**

## What the calls actually are

Most have no domain annotation at all: 22/97 (Bd) and 64/331 (Bsal) carry any Pfam hit.
Among those that do:
- **Bd**: `Chitin_bind_1` (5), `CAP` (3), Cu-oxidase/AA1 laccase (3)
- **Bsal**: `VWD` (28 of 64 — von Willebrand factor type D, a mucin-associated domain),
  `CAP` (5), `Chitin_bind_1` (4), Cu-oxidase (4)

The Cu-oxidase/AA1 hits are the same laccase over-calling seen kingdom-wide. The Bsal VWD
enrichment is a separate lead worth its own look: VWD domains occur in gel-forming mucins.

## Consequence for the curation

None of these proteins can be added to `data/curated/adhesins/` as positives — there is no
functional evidence. They belong in a hypothesis list. The honest position stays as recorded
in `data/curated/surface/README.md`: **there are no curated adhesin labels for chytrids**, so
the classifier's behaviour in Chytridiomycota (the phylum with the highest predicted
`adhesion_fraction` in the kingdom survey) is untested and unvalidated.
