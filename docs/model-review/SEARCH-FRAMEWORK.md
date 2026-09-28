# Framework: using ESM as a search tool for Rhodotorula adhesion and Coccidioides antigens

*Drafted 2026-09-27. The stated goal is **to build tools and test whether the classifier/ESM
is actually useful in this search** — not to assume it is. The two biological searches are the
test cases.*

## Independent review (Fable, 2026-09-27) — verdict and verified defects

An independent feasibility review found four defects in this plan. **I verified each against
live APIs and our own files before accepting them**; all four hold.

| finding | verified |
|---|---|
| **The *C. posadasii* reference is a SOWgp DELETION strain.** `surface.tsv` genome label is literally `Coccidioides posadasii C735 delta SOWgp` (taxon 443226) | yes — it was in our own output and went unnoticed |
| **SOWgp and Ag2/PRA are absent from the antigen ranking.** Q8NK60/Q8NK61/Q96V71/Q12295/Q2TVJ9 all missing from `coccidioides_candidates.tsv` | yes — all five absent. The ranking built to extend SOWgp does not contain SOWgp |
| **The ranking is uncalibrated.** Endochitinase-1 (CiX1, `P0CB51`/`Q1E3R8`) — the complement-fixation antigen in actual clinical serology — scores 4/8, **tied with 678 other proteins** (63% of rows) | yes, exactly |
| **Rhodotorula: wrong species is annotated.** The clinical biofilm organism *R. mucilaginosa* has 125 assemblies but only **3 annotations**; the PF21671 proteins cited are mostly *R. toruloides* (26 vs 13) | yes |

Two further points, verified and consequential:
- **Genus *Rhodotorula* has ZERO hydrophobins** (PF01185 = 0, PF06766 = 0; AMP-binding control
  = 110). If its notable cell-surface hydrophobicity were protein-mediated by the known
  mechanisms, those genes would be present. The likely explanation is carotenoid/lipid/capsular
  — **non-protein**. This is the *Malassezia*/GXM trap recurring, and it is a framing problem,
  not a data problem.
- **The ranking has no fungal cross-reactivity penalty.** It penalizes human homology only, but
  the actual clinical failure mode of coccidioidomycosis serology is cross-reaction with
  *Histoplasma*, *Blastomyces*, *Paracoccidioides* and *Aspergillus*. The current top of the
  list — Gel1 glucanosyltransferases, pepsins, subtilisins, chitinases — is pan-fungal
  conserved, i.e. maximally cross-reactive. Adding this term should substantially reorder it.

### Revised verdict per goal

| goal | verdict |
|---|---|
| 1. *Rhodotorula* adhesins | **Not achievable as stated.** The classifier is blind to the entire plausible candidate space (short, Cys-rich, non-GPI, repeat-poor → ~0.000), the clinical species is unannotated, and the phenotype may not be protein-mediated at all. Reframe as a surface-proteome catalog + genotype–phenotype association, labelled a hypothesis pool. Do not publish a ranked "adhesin" list. |
| 2. *Coccidioides* antigens | **Achievable now and the strongest goal** — after fixing the four defects above. |
| 3. Pangenome variability | **Achievable, but not on the annotated assemblies and not for the genes of most interest.** 9 annotations of the same ~28 Mb genome differ by 48% in gene count (7,253–10,724) — the *Bsal* 0/24/41 failure in advance. Only ~6 assemblies have contig N50 > 150 kb. Do it **from reads** (402 Coccidioides WGS runs in SRA): read depth over target loci gives presence/absence *and* copy number, immune to assembly and annotation. |

**The honest headline for goal 3 is likely repeat-copy-number variation, not gene
presence/absence** — SOWgp alleles are *defined* by repeat number (58/66/82 kDa = 3/5/6
repeats), and tandem arrays are exactly what short-read assemblies collapse.

Also flagged and adopted: **spherule-phase expression** is the most discriminating filter
available and costs only CPU — 249 *Coccidioides* RNA-seq runs are in SRA, and SOWgp and
Ag2/PRA matter precisely because they are spherule-phase. And B-cell (antibody) evidence should
be scored separately from T-cell, since a serodiagnostic needs the former; the current score
merges them.

## The three goals

1. Candidate proteins in *Rhodotorula* involved in biofilm formation and stickiness to
   surfaces or human cells.
2. Antigenic candidates in *Coccidioides*, building on SOWgp and the PRA family.
3. Whether those proteins/genes vary in sequence or presence/absence across the *Coccidioides*
   pangenome.

## Why this is a good test of ESM, and the trap to avoid

The honest problem with goals 1 and 2 is that **neither has enough labels to validate a
prediction**. *Rhodotorula* has **zero** named adhesin proteins in the literature
(`analysis/curation/BASIDIOMYCETE_NOTES.md`); *Coccidioides* had zero curated adhesins until
today and has 5 IEDB antigens. Any ranked candidate list will therefore look plausible and be
unfalsifiable — the same trap already documented for the 27 *Coccidioides* calls in
`STATUS.md` §3.

**The way out is orthogonal validation**: score candidates with ESM, then ask whether the
ranking is enriched for something measured independently.

| search | independent signal available | where |
|---|---|---|
| *Rhodotorula* adhesion | **biofilm phenotype per strain** (SBF index, 11 strains + environment metadata) | `shared/projects/Rhodotorula/Biofilm/00_input/Biofilm_Results_Table.csv` |
| *Coccidioides* antigens | **pangenome presence/absence + CNV + variant calls** | `shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/` |

This turns goal 3 from a separate question into **the validation method for goals 1 and 2**.
Adhesins and antigens are characteristically variable — copy number, repeat length,
presence/absence. If ESM-ranked candidates are enriched for pangenome variability relative to
matched surface proteins, that is evidence the ranking carries signal, *without needing
labels*. If they are not, that is an equally useful negative result.

## Data already on HPCC (found 2026-09-27; confirm before relying on it)

**Coccidioides** — `/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/`
- `Pangenome/` with `OrthoFinder_diamond/`, `Functional_ann/`, `results/`, `protein_count.txt`
  → **orthogroups give presence/absence directly**; no need for the 6-frame method used in `Bsal_VWD`
- `Assembly/` (36 entries), `genome/` (18) incl. *C. posadasii* Silveira annotations
- Variant calling against several references (RS, WA211, Silveira 2022), **CNV plots already computed**
- Sample sheets: `Genotyping/C_immitis.samples.csv`, `C_posadasii.samples.csv`

**Rhodotorula** — `/bigdata/stajichlab/shared/projects/Rhodotorula/`
- `Biofilm/00_input/Biofilm_Results_Table.csv` — `S_ID, OD_A2, Env, OD_A3, Growth_OD, Strain_Name, SBF`
  (11 strains; SBF = specific biofilm formation, normalized to growth)
- `Biofilm/Strain_Metadata/`, growth curves, `Genomics_PCA/`
- `FunFinder_Pangenome/`, `ExtremeRhodotorula_DraftGenomes_v2/`
- Also `ctsai085/projects/Rhodotorula_comparative_genomics/` (metadata.tsv, functional annotation workflow)

Fungi_5k covers 21 *Rhodotorula* and 4 *Coccidioides* genomes, so the 278 GB
`function.duckdb` can supply Pfam/SignalP/CAZy for those but **not** for arbitrary new
assemblies — those need their own annotation.

## Method, per goal

### Goal 1 — *Rhodotorula* adhesion candidates
1. Build the stage-1 surface set per genome (SignalP 6.0 + NetGPI directly, **not** UniProt
   keywords — these are draft genomes without curation).
2. Rank with three independent scores, kept separate rather than merged:
   - ESM C 300M embedding score from the stage-2 model *(expected to underperform here — see
     Risks)*
   - explicit property profile: tandem-repeat content, Cys pattern, GPI, Ser/Thr, β-aggregation
     propensity, surface charge (the §8.0 design in the review)
   - family membership: CPL1-like (PF21671), hydrophobins, CFEM, PA14
3. **Validate against phenotype**: test whether candidate presence/absence or sequence
   variation correlates with SBF across the 11 phenotyped strains.

### Goal 2 — *Coccidioides* antigen candidates
Extend `data/curated/antigens/coccidioides_candidates.tsv` (1,069 surface proteins already
ranked, and it independently recovered Ag2/PRA, Gel1, CFEM and subtilisins). Additions needed:
- the full PRA family as anchors: Ag2/PRA `Q12295`, PRA2 `Q6K1L8`, PRA3 `Q2TVJ9`, plus SOWgp
  `Q8NK60/Q8NK61/Q96V71`. Note UniProt already holds many short (69–77 aa) `pra` entries from
  different *C. posadasii* isolates — **allelic variation at this locus is already visible in
  the databases**, which is a hint for goal 3.
- B-cell epitope surface accessibility, not just whole-protein homology.

### Goal 3 — pangenome variability (and the validation engine)
Using the existing OrthoFinder pangenome and VCFs:
- **presence/absence**: orthogroup membership per strain → core / accessory / strain-specific
- **copy number**: orthogroup size per strain, plus the existing CNV calls
- **sequence variability**: dN/dS and π per gene from the VCFs; repeat-length variation needs
  assembly-level inspection, not short-read calls
- **the test**: are ESM-ranked candidates enriched for accessory/variable status relative to
  length- and expression-matched surface proteins? Report effect size, not just a p-value.

## Risks and where this could go wrong

1. **ESM is a tandem-repeat detector** (§4.7). *Rhodotorula* CPL1-like proteins are short,
   Cys-rich and repeat-poor — the exact architecture ESM scores ~0. **Expect ESM to
   underperform the property profile on goal 1.** That is a legitimate result for "is ESM
   useful here", and the reason all three scores are kept separate.
2. **Presence/absence from fragmented assemblies is the classic false-positive generator.** A
   gene absent from a draft assembly is usually a gap, not a deletion. Any absence call must be
   checked against read coverage at the locus; contig-level assemblies cannot support absence
   claims on their own. This project has already been bitten once by annotation heterogeneity
   (same *Bsal* strain: 0 vs 24 vs 41 copies under three annotations).
3. **Repeat-length variation is invisible to short reads.** SOWgp's defining feature is tandem
   repeat number, and that is exactly what Illumina assemblies collapse. Long-read assemblies
   only.
4. **n=11 for the biofilm phenotype** supports correlation at best, never association testing
   with genome-wide correction.
5. **Ag2/PRA is small and non-repetitive** — ESM scored it 0.000 (`STATUS.md` §4). The tool
   that found SOWgp will not find PRA-like antigens.

## What the honest deliverable is

- Goal 1: a ranked candidate list with explicit per-property evidence, labelled a hypothesis
  set, plus a statement of whether ESM added anything over the property profile.
- Goal 2: an extended candidate table anchored on SOWgp + PRA, with the same caveat.
- Goal 3: a variability table per candidate, and the enrichment test that says whether the
  ESM ranking carries independent signal.

**The scientific result of this exercise is the answer to "is ESM useful here", and a negative
answer is publishable and useful.** It should not be buried if it comes out that way.
