# *Coccidioides* surface antigens: specificity, pangenome variability, and what the ML tools can and cannot do

**Report, 2026-09-27.** Stajich lab, UC Riverside. Analysis by Claude Code (Opus 5) with
J. Stajich. All numbers are reproducible from `analysis/cocci_antigens/` and
`analysis/model_review/`; every claim traces to a script, a database record, or a PMID.

**Status: internal working report. Provisional, not peer-reviewed, no wet-lab validation.**

---

## 1. Summary

Working from the 488-proteome *Coccidioides* pangenome, three results stand out.

1. **The PRA family splits on specificity.** SOWgp and PRA3 have **no detectable homolog** in
   *Histoplasma*, *Blastomyces*, *Paracoccidioides* or *Aspergillus fumigatus*. Ag2/PRA and
   PRA2 both have clear orthologs in all three dimorphic confounders (55–70% identity, full
   query coverage). So does the complement-fixation antigen in clinical use. This predicts
   that Ag2/PRA-based and CF-based assays cross-react and that **SOWgp and PRA3 are the
   species-specific markers** in this set.
2. **SOWgp is not universal.** It is present in **92.0%** of 488 proteomes, with copy-number
   variation (mean 1.09, CV 0.26). Every other anchor is at 98–99.8% and strictly single-copy.
   For a diagnostic aiming at complete sensitivity this is a real liability, and it is
   invisible unless you look across the pangenome.
3. **Each anchor has a distinct and diagnostically meaningful profile** (§3). PRA3 is the only
   one that is simultaneously near-universal, single-copy and specific.
4. **Spherule expression validates SOWgp emphatically but does not generalise.** SOWgp goes
   from 13 TPM in mycelia to 15,000 TPM in 48 h spherules (log2FC +10.05, top percentile of the
   genome), while Ag2/PRA, PRA2 and PRA3 are all *down*-regulated in spherules. Intersecting
   specificity with spherule induction yields **40 uncharacterized candidates** with the SOWgp
   profile (§3.7).

A fourth result is methodological and matters for how this work is presented: **the machine
learning tools contributed the search space, not the discrimination.** The specificity signal
came from orthology and cross-genome comparison, not from the language model — which scores
Ag2/PRA at 0.000 and is blind to this entire protein class (§5).

---

## 2. Background and goals

The goal set by the PI: *"find antigenic candidate proteins in Coccidioides to build on SOWgp
and PRA1/PRA3 to see if additional candidates exist"*, and *"test if these proteins/genes are
variable in terms of sequence or presence/absence in the pangenome"*.

An earlier ranking of 1,069 surface proteins (`data/curated/antigens/coccidioides_candidates.tsv`)
was superseded because an independent review found it unusable for this purpose:

| defect | consequence |
|---|---|
| The *C. posadasii* reference was strain **C735 ΔSOWgp — a SOWgp deletion strain** | the protein the whole exercise builds on was deleted from half the reference data, and SOWgp was absent from the ranking entirely |
| No fungal cross-reactivity term (human only) | the top of the list was pan-fungal conserved — Gel1 glucanosyltransferases, pepsins, subtilisins, chitinases — i.e. maximally cross-reactive |
| 63% of rows tied at one integer score | no usable ranking |
| Known antigens not used as controls | the complement-fixation antigen in clinical use ranked 4/8, tied with 678 others |

---

## 3. Findings

### 3.1 Anchor profiles across 488 proteomes

| protein | locus (*C. immitis* RS) | prevalence | copy no. (mean, CV) | max identity to a confounder | interpretation |
|---|---|---|---|---|---|
| **PRA3** | CIMG_02492 | 98.0% | 1.00, 0.00 | **0%** | near-universal, single-copy, specific — **the cleanest profile** |
| **SOWgp** | CIMG_04613 | **92.0%** | 1.09, **0.26** | **0%** | specific, but **absent from 8% of isolates** and copy-variable |
| Ag2/PRA | CIMG_09696 | 98.8% | 1.00, 0.00 | 60% | universal and single-copy, but **cross-reactive** |
| PRA2 | CIMG_09560 | 98.8% | 1.00, 0.00 | 69% | as Ag2/PRA |
| CF antigen (CiX1 / CTS1) | CIMG_02795 | 99.8% | 1.00, 0.04 | 64% | the clinical antigen: universal, **cross-reactive** |

### 3.2 The cross-reactivity result is orthology, not artifact

Ag2/PRA and PRA2 hit the **same three genes** in all three dimorphic confounders —
*Histoplasma* `FBAD1291_004826`, *Blastomyces* `F00FD2C2_005998`, *Paracoccidioides*
`FD6225C6_004143` — at 55–70% identity with **full query coverage (qcov = 1.000)**. Full-length
alignment at that identity across three genera is orthology, not the low-complexity
compositional matching that proline-rich proteins are prone to. PRA2 additionally hits
*A. fumigatus*.

SOWgp and PRA3 return **no hit at all** in any of the four confounder proteomes.

**This is a testable prediction**: sera from histoplasmosis, blastomycosis and
paracoccidioidomycosis patients should react with Ag2/PRA and with the CF antigen, and should
not react with SOWgp or PRA3.

### 3.3 Pangenome structure

8,542 orthogroups with a *C. immitis* RS representative, across 488 proteomes:

| compartment | orthogroups | fraction |
|---|---|---|
| core (≥99% of proteomes) | 5,657 | 66.2% |
| soft-core (95–99%) | 1,305 | 15.3% |
| shell (15–95%) | 841 | 9.8% |
| cloud (<15%) | 739 | 8.7% |

**1,580 orthogroups are accessory** (<95%), and **286 have copy-number CV > 0.3**. So there is
ample variable gene content, and the anchors sit in it: SOWgp is the only anchor outside the
soft-core.

### 3.4 Sequence variability at the SOWgp locus

The three known SOWgp alleles (58, 66 and 82 kDa; UniProt Q8NK60, Q8NK61, Q96V71) all map to
the **single locus** CIMG_04613 at **96.0%, 85.3% and 74.4%** identity with full coverage.
That descending identity is the published tandem-repeat-number difference (4, 5 and 6 copies
of a 41–47 aa Pro/Asp-rich repeat; Hung et al. 2002) showing up directly in the alignment.

**Caveat, and it is a serious one for this specific locus**: the prevalence and copy-number
figures above come from an assembly-derived orthogroup table, and tandem arrays are exactly
what short-read assemblies collapse or fragment. SOWgp's 92% prevalence and CV 0.26 must be
confirmed from **read depth** over the locus before being relied on. An absence call from a
fragmented assembly is usually a gap, not a deletion.

### 3.5 A 595-orthogroup candidate shortlist — and its honest limitation

`analysis/cocci_antigens/shortlist_specific_universal.tsv`: orthogroups present in ≥95% of 488
proteomes with **no detectable homolog** in any of the four confounders or in human. 595
orthogroups, **557 of them single-copy**, 38 additionally copy-number variable.

For context, across all 8,542 orthogroups only **22% have no detectable confounder homolog**
and **28% are ≥70% identical** to one. Specificity is genuinely scarce.

**The limitation: zero of the 595 have an IEDB antigen homolog.** That is not a coincidence,
it is structural. The antigenicity axis measures homology to curated IEDB antigens, which are
dominated by other fungi — so *by construction* it penalises exactly the species-specific
proteins we are looking for. SOWgp and PRA3 themselves sit at the 75th–80th percentile on it.

**The shortlist is therefore specific but immunologically uncharacterized.** Specificity is
doing all the work; antigenicity is not yet evidenced.

### 3.6 Spherule-phase expression — tested, and it does not do what was expected

The lab's own RNA-seq (`jstajich/projects/Coccidioides_UCSD_SpheruleMycelium`, *C. immitis* RS,
kallisto TPM, 2 replicates each of mycelia / spherule 48 h / spherule 8 d) joined to **all
8,542** ranked orthogroups by gene ID.

| protein | mycelia TPM | spherule 48 h | spherule 8 d | log2FC (48 h) | percentile |
|---|---|---|---|---|---|
| **SOWgp** | 13.1 | **15,000.4** | 4,839.2 | **+10.05** | **100th** |
| CF antigen (CiX1) | 2.1 | 8.9 | 9.1 | +1.68 | 79th |
| Ag2/PRA | 2,330.3 | 694.0 | 564.3 | **−1.75** | 4th |
| PRA3 | 24.2 | 5.5 | 3.7 | **−1.96** | 4th |
| PRA2 | 207.6 | 15.0 | 9.9 | **−3.71** | 1st |

**SOWgp is spectacularly confirmed**: essentially off in mycelia and among the most abundant
transcripts in the 48 h spherule — a ~1,000-fold induction, at the very top of the genome. The
48 h > 8 d ordering matches the published "elevated during early spherule development"
(Hung et al. 2002). This is independent validation that the data and the ID join are correct.

**But the PRA family goes the other way.** Ag2/PRA, PRA2 and PRA3 are all *mycelia*-high and
down-regulated in spherules in this dataset. Ag2/PRA is abundant (2,330 TPM) — but in the
wrong phase.

**Consequence, and it is a warning rather than a win:** spherule induction is *not* a general
antigenicity filter for this protein set. Applied naively it would have promoted SOWgp and
**demoted PRA3, the best specificity candidate in the panel**. It identifies one specific
class — parasitic-phase surface antigens of the SOWgp type — and should be used as a labelled
axis, not folded into a single score. This is the same mistake as v1's merged score, and the
report avoids repeating it.

### 3.7 The intersected shortlist: 40 candidates

Requiring **all four** criteria — ≥95% prevalence across 488 proteomes, no confounder homolog,
no human homolog, and spherule-48 h induction (log2FC > 1, TPM > 50) — gives **40 orthogroups**
(`shortlist_spherule_induced.tsv`). Top by spherule abundance:

| protein | mycelia | spherule 48 h | log2FC | prevalence |
|---|---|---|---|---|
| CIMG_06250 | 3.8 | 691.9 | +7.18 | 1.00 |
| CIMG_03452 | 48.3 | 647.6 | +3.72 | 0.99 |
| CIMG_04662 | 0.5 | 480.2 | **+8.32** | 0.99 |
| CIMG_01584 | 8.6 | 454.5 | +5.56 | 1.00 |
| CIMG_05599 | 2.9 | 257.6 | +6.04 | 1.00 |
| CIMG_06249 | 14.1 | 257.5 | +4.10 | 1.00 |

These have the SOWgp *profile* — near-silent in mycelia, strongly induced in early spherules,
universal across isolates, no homolog in the confounding fungi — without being SOWgp.
**CIMG_06249 and CIMG_06250 are adjacent**, suggesting a locus worth looking at directly.

SOWgp itself is excluded only by the prevalence filter (92%, see §3.1).

**None of these 40 has been characterized.** They are a hypothesis set with a defined and
unusually specific profile, not validated antigens.

---

## 4. What would make this decisive

In priority order:

1. ~~Spherule-phase expression~~ — **done** (§3.6–3.7), using the lab's own RS RNA-seq. It
   validated SOWgp emphatically and produced a 40-candidate intersected shortlist, but it did
   *not* generalise to the PRA family. Worth extending to *C. posadasii* and to the 249 SRA
   runs for replication, and to endospore stage — SOWgp is reportedly depleted on endospores
   (Hung et al. 2007), which is an immune-evasion angle this dataset cannot address.
2. **Characterize the 40.** Domain/signal-peptide annotation, structure, and whether any are
   already in the IEDB or proteomics data. This is cheap and immediately informative.
3. **Read-depth confirmation of SOWgp prevalence and copy number** (402 *Coccidioides* WGS runs
   in SRA), which is immune to assembly fragmentation and annotation heterogeneity.
4. **B-cell epitope surface accessibility**, scored separately from T-cell evidence — a
   serodiagnostic needs antibody epitopes.
5. **Serological test of the cross-reactivity prediction** in §3.2, which is a direct wet-lab
   experiment with an unambiguous outcome.

---

## 5. What the ML tools contributed — and did not

This project set out partly to test whether a protein language model is useful in a search
like this. On this evidence:

- **The language model did not find these antigens.** It scores Ag2/PRA at **0.000** and is
  documented (`docs/model-review/`, §4.7) to be functionally a *tandem-repeat surface protein
  detector*: every adhesin it misses has zero tandem repeats, and it is blind to short,
  cysteine-rich, non-GPI proteins — which is what most of this candidate space looks like.
- **Where it does earn its place** is on proteins with no domain annotation: PR-AUC 0.716 vs
  0.246 for the best domain-based baseline on the domain-blind subset (§4.6). That is a real
  capability, but it is orthogonal to the specificity question that drove these findings.
- **The discriminating signal here came from comparative genomics**, not ML: orthology against
  confounder proteomes, and prevalence/copy number across a 488-proteome pangenome.

The honest framing for any write-up: **ML narrowed and organized the search space; cross-genome
comparison produced the result.** Presenting the language model as the discovery engine would
misstate what happened.

---

## 6. Limitations

1. **Four usable positive controls.** The ranking's acceptance test passes 3/4, and the script
   prints `NOT CALIBRATED` whenever it is not 4/4. That warning is left switched on deliberately.
2. **Assembly-derived counts** for prevalence and copy number (§3.4).
3. **One reference for the universe.** Proteins absent from *C. immitis* RS are not scored,
   so *C. posadasii*-specific antigens are systematically missed. A second pass anchored on
   the Silveira proteome would close this.
4. **Cross-reactivity was assessed against four confounder genomes**, one strain each. Broader
   sampling, especially within *Histoplasma*, would firm up the specificity calls.
5. **No experimental validation of anything in this report.**

---

## 7. Sources and provenance

**Data**
- *Coccidioides* pangenome: `/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome` — 493 input proteomes, OrthoFinder result `Cocci_496_OG2_5_5`, 15,857 orthogroups
- Confounder proteomes from `Fungi_5k/input`: *Histoplasma capsulatum* G186AR, *Blastomyces dermatitidis* ER-3, *Paracoccidioides brasiliensis* Pb18, *Aspergillus fumigatus* Af293
- Human proteome UP000005640 (UniProt); IEDB Query API for curated fungal antigens
- Anchors from UniProt: Q8NK60, Q8NK61, Q96V71 (SOWgp); Q12295, A0A0E1RVD3 (Ag2/PRA); Q6K1L8 (PRA2); Q2TVJ9 (PRA3); Q1E3R8, P0CB51 (CF antigen)

**Literature**
- Hung CY, Yu JJ, Seshan KR, Reichard U, Cole GT. 2002. A parasitic phase-specific adhesin of *Coccidioides immitis* contributes to the virulence of this respiratory fungal pathogen. *Infect Immun* 70:3443-56. https://doi.org/10.1128/IAI.70.7.3443-3456.2002 — SOWgp binds laminin > fibronectin > collagen IV; deletion reduces ECM binding and virulence; 4–6 tandem repeats, size varies by isolate
- Hung CY, Xue J, Cole GT. 2007. Virulence mechanisms of *Coccidioides*. *Ann N Y Acad Sci* 1111:225-35. https://doi.org/10.1196/annals.1406.020

**Code** (all in `adhesionPred`)
- `analysis/cocci_antigens/01_build_inputs.sh` — orthology and cross-reactivity searches
- `analysis/cocci_antigens/02_score_antigens.py` — scoring and acceptance test
- `analysis/cocci_antigens/NOTES.md` — method detail and the v1 → v2 defect list
- `docs/model-review/` — the model evaluation underlying §5
