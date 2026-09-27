# *Batrachochytrium* surface/adhesion factor investigations

> **Curation decision (2026-09-27, JS): CBM18 is OFF the table as an adhesion candidate.**
> The chitin-binding rationale does not support adhesion to keratinized amphibian skin, and
> there is no functional evidence. Section 1 is kept as the record of why — it is a negative
> result about the classifier's behaviour, not a lead. Do not enter CBM18 proteins in any
> adhesin table. The open question for CBM18 is self-masking / cell-wall biology, which is a
> different project.
>
> **VWD (section 2) remains an open lead**, now being followed up in the separate
> `Bsal_VWD` repository.

## 1. Is the CBM18 expansion what the classifier calls "adhesion"?

*2026-09-27. Section 1 prompted by the question of whether CBM18 genes relate to adhesion or to immune
"cloaking". Data: Fungi_5k `function.duckdb` (Pfam, SignalP) joined to the existing
`adhesion_predict` result files. Scripts: `01_cbm18_vs_calls.py`, `02_cbm18_architecture.py`.
Per-protein table: `cbm18_architecture.tsv`.*

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

---

## 2. A *B. salamandrivorans*-specific expansion of secreted VWD proteins

*Scripts: `03_domain_enrichment.py`, `04_vwd_profile.py`. Tables:
`domain_enrichment_chytrids.tsv`, `vwd_chytrids.tsv`.*

Systematic Pfam enrichment among the adhesion calls (Fisher, BH-corrected, per genome)
turned up **VWD (von Willebrand factor type D, PF00094) as the strongest signal anywhere in
these four genomes**: 28 of 41 *Bsal* VWD proteins are called, OR = 113, q = 3.3e-37.

The family itself is what is striking, independent of the classifier:

| genome | proteome | VWD proteins | per 1,000 proteins | called | signal peptide | TM |
|---|---|---|---|---|---|---|
| *B. salamandrivorans* AMFP13 | 16,260 | **41** | 2.52 | 28 | 30 | 4 |
| *B. dendrobatidis* JAM81 | 7,021 | 1 | 0.14 | 1 | 1 | 0 |
| *Homolaphlyctis polyrhiza* | 6,632 | 0 | 0 | — | — | — |
| *Spizellomyces punctatus* | 8,319 | 1 | 0.12 | 1 | 1 | 0 |

41 copies in *Bsal* against 1 in *Bd* and 0–1 in the outgroups. *Bsal*'s proteome is 2.3x
larger than *Bd*'s, so genome expansion alone does not explain a 41-fold difference. The
proteins are **VWD-only** (no partner Pfam domain in any of the 41), **secreted** (30/41 with
a signal peptide, only 4 with a TM helix), and range 170–1,163 aa (median 670).

**Why this is interesting.** In animals the VWD domain is the multimerization module of
gel-forming mucins (MUC2, MUC5AC) and von Willebrand factor; it drives disulphide-linked
oligomerization into gels. A secreted, VWD-only family expanded specifically in the
salamander-skin pathogen is a plausible host-interface family — either forming an
extracellular matrix/gel itself, or interacting with host skin mucus, which salamanders
produce heavily. That is a hypothesis. We found **no literature on VWD domains in
*Batrachochytrium***; searches returned nothing, so as far as we can tell this expansion has
not been described.

**Caveats.** One annotation, one *Bsal* assembly (GCA_002006685.2). Pfam VWD hits should be
confirmed by HMM search and the family checked for repeat/TE association before any claim of
expansion. It is also exactly the kind of low-complexity, cysteine-rich secreted protein that
annotation pipelines over- and under-call. Worth a dedicated look; not yet a finding.

---

## 3. What the classifier ignores: the real virulence factors

*Script: `05_known_virulence_families.py`. Counts are proteins (called, signal peptide).*

| family | *Bd* JAM81 | *Bsal* AMFP13 | *H. polyrhiza* | *S. punctatus* |
|---|---|---|---|---|
| M36 fungalysin | 37 (**0** call, 27 SP) | **342** (**0** call, 216 SP) | 3 (0, 1) | 3 (0, 3) |
| S41 protease | 36 (**0** call, 21 SP) | **161** (**0** call, 86 SP) | 6 (0, 3) | 1 (0, 1) |
| aspartyl protease | 40 (0 call, 2 SP) | 19 (0 call, 2 SP) | 8 (0, 4) | 5 (0, 4) |
| tyrosinase | 13 (0 call, 10 SP) | 26 (0 call, 10 SP) | 13 (0, 9) | 3 (0, 3) |
| chitin deacetylase | 28 (1 call, 23 SP) | 40 (2 call, 24 SP) | 12 (0, 7) | 12 (0, 7) |
| CBM18 | 10 (**5** call, 8 SP) | 13 (**4** call, 8 SP) | 1 (0, 0) | 3 (0, 2) |
| VWD | 1 (**1** call, 1 SP) | 41 (**28** call, 30 SP) | 0 | 1 (1, 1) |
| CAP/PR-1 | 11 (3 call, 5 SP) | 13 (5 call, 8 SP) | 2 (2, 1) | 2 (1, 2) |
| Cu-oxidase/AA1 | 3 (**3** call, 2 SP) | 4 (**4** call, 3 SP) | 4 (3, 1) | 2 (2, 1) |

**The M36 metalloproteases — the best-documented chytrid virulence family — are never
called.** Not one of 342 secreted M36 proteins in *Bsal*, nor 37 in *Bd*. Same for the S41
protease expansion (161 in *Bsal*). This is the expected behaviour and a useful negative
control: the model detects Ser/Thr-rich surface glycoproteins and binding modules, not
secreted enzymes, so it is orthogonal to the protease-centred virulence literature. It is
**not** a virulence-factor detector, and should never be presented as one.

Conversely the AA1 laccases are called in every genome (3/3, 4/4, 3/4, 2/2), which is the same
kingdom-wide laccase false-positive documented in `analysis/adhesion_properties/REPORT.md`.

The M36 count deserves a flag: Wacker et al. 2023 report n=177 M36 genes in *Bsal* from a
nanopore assembly; this annotation of GCA_002006685.2 gives 342. Same direction, different
magnitude — an annotation/assembly difference, not a contradiction, but do not quote either
number as settled.

---

## 4. Summary of what would actually move this forward

1. **HMM-based recount** of CBM18 and VWD in both species against current assemblies. Pfam
   gathering thresholds undercount diverged repeats, and both families are repeat-rich.
2. **Are these in the repeat-rich compartment?** Wacker et al. 2023 show *Batrachochytrium*
   has a two-speed genome with virulence factors in gene-sparse/repeat-rich regions. Testing
   whether the VWD and CBM18 families sit there too is a direct, cheap test of whether they
   behave like the known virulence families.
3. **Expression during infection.** Neither family has been checked against in-host
   expression. Existing *Bsal* infection RNA-seq would answer whether they are induced on
   salamander skin.
4. **Binding assays** are the only way to separate adhesion from cloaking for CBM18, and to
   test whether VWD proteins form gels or bind host mucus.

None of these proteins should enter `data/curated/adhesins/` as positives: there is no
functional evidence for any of them. The chytrid label gap recorded in
`data/curated/surface/README.md` stands.

## Sources

- Abramyan J, Stajich JE. 2012. Species-specific chitin-binding module 18 expansion in the amphibian pathogen *Batrachochytrium dendrobatidis*. *mBio* 3(3):e00150-12. https://doi.org/10.1128/mBio.00150-12
- Wacker T, Helmstetter N, Wilson D, Fisher MC, Studholme DJ, Farrer RA. 2023. Two-speed genome evolution drives pathogenicity in fungal pathogens of animals. *PNAS* 120(2):e2212633120. https://doi.org/10.1073/pnas.2212633120 — *Bsal* is the most repeat-rich Chytridiomycota genome (40.9% repeats), ~3x the length of *Bd*; M36 metalloproteases highly expanded (n=177), 53% flanked by TEs; virulence factors enriched in gene-sparse/repeat-rich compartments.
- Van Rooij P, Martel A, Haesebrouck F, Pasmans F. 2015. Amphibian chytridiomycosis: a review with focus on fungus-host interactions. *Vet Res* 46:137. https://doi.org/10.1186/s13567-015-0266-0 — source of the vinculin/fibronectin/fasciclin adhesion claims assessed in §1.
- Rosenblum EB, Stajich JE, Maddox N, Eisen MB. 2008. Global gene expression profiles for life stages of the deadly amphibian pathogen *Batrachochytrium dendrobatidis*. *PNAS* 105(44):17034-9. https://doi.org/10.1073/pnas.0804173105
- Rosenblum EB, Poorten TJ, Joneson S, Settles M. 2012. Substrate-specific gene expression in *Batrachochytrium dendrobatidis*, the chytrid pathogen of amphibians. *Mol Ecol* 21(13):3110-20. https://doi.org/10.1111/j.1365-294X.2012.05481.x
- No literature found describing VWD-domain proteins in *Batrachochytrium* (searched 2026-09-27).
