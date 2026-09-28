# *Coccidioides* surface antigens: specificity, spherule expression, and pangenome variability

**Working report — 2026-09-27.** Stajich lab, UC Riverside.
Analysis: Claude Code (Opus 5) with J. Stajich.

> **Status: internal working report, shared for comment. Provisional. No experimental
> validation of anything below. Every number is reproducible from `analysis/cocci_antigens/`;
> every claim traces to a script, a database record, or a PMID. Please read §7 (limitations)
> before using any candidate list.**

---

## 1. Executive summary

Working from the **488-proteome *Coccidioides* pangenome**, plus spherule/mycelium RNA-seq and
cross-genome comparison against the fungi that confound coccidioidomycosis serology:

1. **The PRA/SOWgp panel splits cleanly on species specificity.** SOWgp and PRA3 have **no
   detectable homolog** in *Histoplasma*, *Blastomyces*, *Paracoccidioides* or *Aspergillus
   fumigatus*. Ag2/PRA, PRA2 and the complement-fixation antigen all have clear orthologs in
   the dimorphic confounders (55–70% identity, full-length). This predicts that Ag2/PRA-based
   and CF-based assays cross-react, and identifies **SOWgp and PRA3 as the species-specific
   markers** in this panel.
2. **SOWgp is not universal.** It is present in **92.0%** of 488 proteomes, and copy-variable
   (mean 1.09, CV 0.26). Every other anchor is at 98–99.8% and strictly single-copy. For an
   assay aiming at complete sensitivity this is a real liability, invisible without a pangenome.
3. **Spherule-phase expression validates SOWgp spectacularly but does not generalise.** SOWgp
   goes from 13 TPM in mycelia to **15,000 TPM** at spherule 48 h (log2FC **+10.05**, top
   percentile of the genome). But Ag2/PRA, PRA2 and PRA3 are all *down*-regulated in spherules.
   Spherule induction is therefore a **labelled axis, not a universal antigenicity filter** —
   using it naively would have demoted PRA3, the best specificity candidate in the panel.
4. **14 Tier-1 candidates** satisfy all four criteria simultaneously (secreted, universal,
   *Coccidioides*-specific, spherule-induced), and **45 Tier-2** satisfy the first three (§5).
   None is characterized. All are hypotheses.
5. **Methodological, and important for how this is written up:** the discriminating signal came
   from **comparative genomics, not machine learning**. The protein language model scores
   Ag2/PRA at 0.000 and is blind to this protein class (§6).

---

## 2. Goals and what preceded this

Goals set by the PI: *find antigenic candidates in* Coccidioides *building on SOWgp and
PRA1/PRA3*, and *test whether those genes vary in sequence or presence/absence across the
pangenome*.

An earlier ranking of 1,069 surface proteins was **superseded** after independent review found
four defects. Recorded here because they are instructive:

| defect | consequence |
|---|---|
| The *C. posadasii* reference used was strain **C735 ΔSOWgp — a SOWgp deletion strain** | the protein the exercise was built on was deleted from half the reference data; SOWgp was absent from the ranking entirely |
| No fungal cross-reactivity term (human only) | the top of the list was pan-fungal conserved — Gel1 glucanosyltransferases, pepsins, subtilisins, chitinases — i.e. maximally cross-reactive |
| 63% of rows tied at one integer score | no usable ranking |
| Known antigens not used as controls | the CF antigen in clinical use ranked 4/8, tied with 678 others |

The rebuild scores the **pangenome**, adds a fungal cross-reactivity penalty, is continuous,
deduplicates by orthogroup, and is held to an acceptance test that prints `NOT CALIBRATED`
when known antigens are missed.

---

## 3. Data and methods

| input | source |
|---|---|
| 488-proteome pangenome, 15,857 orthogroups | `shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome` (OrthoFinder `Cocci_496_OG2_5_5`) |
| scoring universe | *C. immitis* RS (FungiDB), 9,139 proteins assigned to orthogroups → **8,542 orthogroup representatives** |
| confounder proteomes | `Fungi_5k/input`: *Histoplasma capsulatum* G186AR, *Blastomyces dermatitidis* ER-3, *Paracoccidioides brasiliensis* Pb18, *Aspergillus fumigatus* Af293 |
| human proteome | UniProt UP000005640 |
| secretion / domains | Fungi_5k `function.duckdb` (SignalP, TMHMM, Pfam), joined via an MMseqs2 id map |
| spherule expression | `jstajich/projects/Coccidioides_UCSD_SpheruleMycelium`, *C. immitis* RS, kallisto TPM, 2 replicates each of mycelia / spherule 48 h / spherule 8 d. **Data from Carlin et al. 2021** (mycelia, young and mature spherules); see §10. Related regulatory work: Duttke et al. 2022 (csRNA-seq of the phase transition). |
| anchors | UniProt Q8NK60, Q8NK61, Q96V71 (SOWgp); Q12295, A0A0E1RVD3 (Ag2/PRA); Q6K1L8 (PRA2); Q2TVJ9 (PRA3); Q1E3R8, P0CB51 (CF antigen) |

**Score.** Two axes, deliberately kept separate:
- *antigenicity* = 2.5·prevalence + 2.0·(homology to curated IEDB fungal antigens)
- *specificity* = −3.0·(max % identity to a confounder) − 0.5·(breadth of confounders) − 1.5·(human homology)

Copy-number variability is **reported, not scored** — it answers the pangenome question but is
a liability for a diagnostic.

---

## 4. Findings

### 4.1 Anchor profiles across 488 proteomes

| protein | locus (*C. immitis* RS) | prevalence | copy no. (mean, CV) | max ident. to confounder | signal peptide | interpretation |
|---|---|---|---|---|---|---|
| **PRA3** | CIMG_02492 | 98.0% | 1.00, 0.00 | **0%** | yes | near-universal, single-copy, specific — **cleanest profile** |
| **SOWgp** | CIMG_04613 | **92.0%** | 1.09, **0.26** | **0%** | *not annotated* | specific, but **absent from 8% of isolates**, copy-variable |
| Ag2/PRA | CIMG_09696 | 98.8% | 1.00, 0.00 | 60% | yes | universal, single-copy, but **cross-reactive** |
| PRA2 | CIMG_09560 | 98.8% | 1.00, 0.00 | 69% | yes | as Ag2/PRA |
| CF antigen (CiX1/CTS1) | CIMG_02795 | 99.8% | 1.00, 0.04 | 64% | yes | the clinical antigen: universal but **cross-reactive** |

### 4.2 Cross-reactivity is orthology, not artifact

Ag2/PRA and PRA2 hit the **same three genes** in all three dimorphic confounders —
*Histoplasma* `FBAD1291_004826`, *Blastomyces* `F00FD2C2_005998`, *Paracoccidioides*
`FD6225C6_004143` — at 55–70% identity with **full query coverage (qcov = 1.000)**.
Full-length alignment at that identity across three genera is orthology, not the
low-complexity compositional matching that proline-rich proteins are prone to. PRA2
additionally hits *A. fumigatus*. SOWgp and PRA3 return **no hit at all** in any confounder.

> **Testable prediction.** Sera from histoplasmosis, blastomycosis and paracoccidioidomycosis
> patients should react with Ag2/PRA and with the CF antigen, and should **not** react with
> SOWgp or PRA3. This is a direct wet-lab experiment with an unambiguous outcome.

Context: across all 8,542 orthogroups, only **22%** have no detectable confounder homolog and
**28%** are ≥70% identical to one. Specificity is genuinely scarce.

### 4.3 Pangenome structure

| compartment | orthogroups | fraction |
|---|---|---|
| core (≥99% of 488 proteomes) | 5,657 | 66.2% |
| soft-core (95–99%) | 1,305 | 15.3% |
| shell (15–95%) | 841 | 9.8% |
| cloud (<15%) | 739 | 8.7% |

**1,580 orthogroups are accessory** (<95%); **286 have copy-number CV > 0.3**. SOWgp is the
only anchor outside the soft-core.

### 4.4 Sequence variability at the SOWgp locus

The three known SOWgp alleles (58, 66, 82 kDa; Q8NK60/Q8NK61/Q96V71) all map to the **single
locus** CIMG_04613 at **96.0%, 85.3% and 74.4%** identity with full coverage. That descending
identity is the published tandem-repeat-number difference (4, 5 and 6 copies of a 41–47 aa
Pro/Asp-rich repeat; Hung et al. 2002) appearing directly in the alignment.

> **Warning specific to this locus.** Prevalence and copy number here come from an
> assembly-derived orthogroup table, and tandem arrays are exactly what short-read assemblies
> collapse or fragment. **SOWgp's 92% prevalence and CV 0.26 must be confirmed from read depth**
> (402 *Coccidioides* WGS runs are in SRA) before being relied on. An absence call from a
> fragmented assembly is usually a gap, not a deletion.

### 4.5 Spherule-phase expression

| protein | mycelia TPM | spherule 48 h | spherule 8 d | log2FC (48 h) | percentile |
|---|---|---|---|---|---|
| **SOWgp** | 13.1 | **15,000.4** | 4,839.2 | **+10.05** | **100th** |
| CF antigen | 2.1 | 8.9 | 9.1 | +1.68 | 79th |
| Ag2/PRA | 2,330.3 | 694.0 | 564.3 | **−1.75** | 4th |
| PRA3 | 24.2 | 5.5 | 3.7 | **−1.96** | 4th |
| PRA2 | 207.6 | 15.0 | 9.9 | **−3.71** | 1st |

SOWgp is essentially off in mycelia and among the most abundant transcripts in the 48 h
spherule — ~1,000-fold induction, top of the genome. The 48 h > 8 d ordering matches the
published "elevated during early spherule development" (Hung et al. 2002), which is independent
evidence the dataset and the ID join are correct.

**Independently corroborated by the source study.** Carlin et al. 2021, reporting this
dataset, note that genes highly upregulated in young spherules include "a spherule surface
protein" and that "genes that are unique to *Coccidioides* spp. are also overrepresented in
this group". Both observations match what falls out of the analysis here: SOWgp is the
top-percentile spherule-induced surface protein, and the Tier-1 candidates (§5.2) are by
construction *Coccidioides*-specific **and** spherule-induced — the intersection that paper
flags as enriched. Our candidate criteria were derived independently, so the agreement is
support for the approach rather than a circular result.

> **Warning.** The PRA family goes the *other* way: Ag2/PRA, PRA2 and PRA3 are all mycelia-high
> and down in spherules here. **Spherule induction is not a general antigenicity filter for this
> protein set.** Applied naively it promotes SOWgp and demotes PRA3 — the best specificity
> candidate in the panel. It identifies one class (parasitic-phase surface antigens of the
> SOWgp type) and is reported as a labelled axis, never folded into a single score.

---

## 5. Candidate lists

### 5.1 A correction worth recording

A first shortlist of 40 (specific + universal + spherule-induced) turned out on annotation to
contain **zero proteins with a predicted signal peptide** and 36 of 40 with no Pfam domain at
all. The filter had never required secretion. **Those 40 are not serodiagnostic candidates** —
they may be interesting spherule biology, but a serodiagnostic antigen must be secreted or
surface-exposed. The file is retained as `shortlist_spherule_induced.tsv` with this caveat.
The lists below add the secretion requirement.

### 5.2 Tier 1 — secreted + universal + *Coccidioides*-specific + spherule-induced (n = 14)

`analysis/cocci_antigens/TIER1_candidates.tsv`

| protein | len | TM | mycelia TPM | spherule 48 h | log2FC | prevalence | Pfam |
|---|---|---|---|---|---|---|---|
| CIMG_00143 | 124 | no | 12.3 | 40.7 | +1.65 | 0.97 | — |
| CIMG_02078 | 124 | no | 10.5 | 22.2 | +1.01 | 1.00 | — |
| CIMG_09680 | 385 | yes | 8.8 | 19.3 | +1.05 | 0.99 | — |
| CIMG_10280 | 628 | yes | 7.2 | 19.3 | +1.30 | 0.98 | — |
| CIMG_09000 | 199 | yes | 6.9 | 17.8 | +1.24 | 1.00 | — |
| CIMG_00211 | 149 | no | 4.6 | 16.7 | +1.66 | 0.99 | — |
| CIMG_09828 | 377 | no | 2.9 | 14.3 | +1.97 | 0.99 | PAN_1; PAN_4 |
| CIMG_06630 | 431 | no | 3.0 | 14.2 | +1.93 | 0.99 | DA_C |
| CIMG_10135 | 444 | no | 3.1 | 13.2 | +1.81 | 0.98 | — |
| CIMG_02158 | 316 | no | 4.0 | 11.3 | +1.29 | 0.96 | — |
| CIMG_02073 | 759 | no | 1.5 | 8.0 | +1.85 | 1.00 | DNase_NucA_NucB |
| CIMG_05465 | 466 | yes | 1.7 | 7.0 | +1.56 | 1.00 | — |
| CIMG_02833 | 522 | no | 1.5 | 4.1 | +1.03 | 1.00 | — |
| CIMG_09495 | 214 | no | 1.1 | 3.4 | +1.07 | 0.99 | — |

> **Warning: these are modestly expressed.** The highest is 40.7 TPM against SOWgp's 15,000.
> None has the abundance profile that makes SOWgp a good serological target. Ten of fourteen
> have no Pfam domain, so function is unknown. `PAN_1/PAN_4` (CIMG_09828) is a
> protein-interaction/adhesion module and is the most interesting on architecture alone.

### 5.3 Tier 2 — secreted + universal + *Coccidioides*-specific (n = 45)

`analysis/cocci_antigens/TIER2_candidates.tsv`. Superset of Tier 1 without the expression
requirement. Use this if spherule-stage expression is not a requirement, e.g. for antigens
detectable in mycelial-phase laboratory exposure.

### 5.4 Acceptance test — the ranking is only partly calibrated

| anchor | combined percentile | antigenicity | specificity | verdict |
|---|---|---|---|---|
| PRA3 | 6.1% | 75.0% | 6.0% | PASS |
| Ag2/PRA | 7.2% | 0.0% | 50.0% | PASS |
| SOWgp | 7.3% | 79.7% | 7.3% | PASS |
| PRA2 | 10.7% | 0.5% | 71.7% | FAIL (just outside) |
| CF antigen | 56.1% | 17.9% | 63.4% | PASS as a *negative* specificity control |

**3/4 *Coccidioides*-specific anchors in the top decile.** The script prints `NOT CALIBRATED`
whenever that is not 4/4, and that warning is left switched on. Four usable controls is thin.

---

## 6. What the ML tools contributed — and did not

This work doubled as a test of whether a protein language model helps in a search like this.

- **The language model did not find these antigens.** It scores Ag2/PRA at **0.000**. It is
  documented (`docs/model-review/` §4.7) to be functionally a *tandem-repeat surface protein
  detector*: every adhesin it misses has zero tandem repeats, and it is blind to short,
  cysteine-rich, non-GPI proteins — which is most of this candidate space.
- **Where it does earn its place** is on proteins lacking domain annotation: PR-AUC 0.716 vs
  0.246 for the best domain-based baseline on the domain-blind subset. Real, but orthogonal to
  the specificity question that drove these findings.
- **The discriminating signal came from comparative genomics**: orthology against confounder
  proteomes, prevalence across 488 proteomes, and phase-specific expression.

> **For any write-up:** ML narrowed and organised the search space; cross-genome comparison
> produced the result. Presenting the language model as the discovery engine would misstate
> what happened.

---

## 7. Limitations — please read before using any list

1. **Annotation coverage is incomplete.** Only **5,817 of 9,139** reference proteins (64%) map
   into the Fungi_5k annotation, so signal-peptide status is *unknown*, not negative, for 36%.
2. **The secreted set is almost certainly under-called.** Only **371 of ~9,910** RS proteins
   (~4%) carry a SignalP annotation, well below the ~10% typical for a fungal proteome.
   Tier 1/2 are therefore **conservative and incomplete**.
3. **SOWgp itself has no Fungi_5k match**, so it would be *excluded* by the secretion filter
   that defines Tier 1/2 — a direct demonstration of limitation 1. The absence of a known
   surface antigen from the filtered lists shows the filter misses real antigens.
4. **Assembly-derived prevalence and copy number** (§4.4), unconfirmed by read depth.
5. **One reference genome defines the universe** (*C. immitis* RS). *C. posadasii*-specific
   antigens are systematically missed. A second pass anchored on the Silveira proteome would
   close this.
6. **Cross-reactivity assessed against four confounder genomes, one strain each.** Broader
   sampling, especially within *Histoplasma*, would firm up the specificity calls.
7. **Expression is one experiment, one strain, *C. immitis* RS, n = 2 per condition.** No
   statistical testing was applied; log2FC is a descriptive ratio of replicate means.
8. **The endospore stage is not covered.** SOWgp is reportedly depleted on endospores (Hung et
   al. 2007) — an immune-evasion angle this dataset cannot address.
9. **No experimental validation of anything in this report.**

---

## 8. Recommended next steps

1. **Serological test of the §4.2 prediction** — cross-reactivity of Ag2/PRA and CF antigen vs
   specificity of SOWgp/PRA3, using heterologous patient sera. Direct, unambiguous, and the
   highest-value experiment here.
2. **Read-depth confirmation of SOWgp prevalence and copy number** across the 402 SRA WGS runs
   — resolves limitation 4 and tests whether the 8% absence is real or assembly artifact.
3. **Repair the secretion annotation**: run SignalP 6.0 and NetGPI directly over the RS and
   Silveira proteomes rather than relying on partial Fungi_5k coverage. This is cheap and
   directly widens Tiers 1–2.
4. **Extend to *C. posadasii*** — second reference universe, and replicate the expression
   analysis.
5. ***C. posadasii* host-outcome phenotypes.** If strain-level variation in infection outcome
   is available, the accessory/copy-variable gene set (§4.3: 1,580 accessory orthogroups, 286
   copy-variable) can be tested for association with outcome. That would convert a descriptive
   pangenome into a hypothesis-generating genotype–phenotype screen, and it is the most
   promising route to *functional* significance for any of these candidates. Worth scoping
   what phenotype data exists and how many strains overlap the 488 proteomes.
6. **Characterize Tier 1** — structure prediction and B-cell epitope accessibility for the ten
   with no Pfam domain.

---

## 9. Files

All under `analysis/cocci_antigens/`:

| file | contents |
|---|---|
| `TIER1_candidates.tsv` | 14 secreted + universal + specific + spherule-induced |
| `TIER2_candidates.tsv` | 45 secreted + universal + specific |
| `cocci_antigen_ranking.tsv` | all 9,139 proteins / 8,542 orthogroup representatives, full evidence columns |
| `shortlist_specific_universal.tsv` | 595 universal + specific (no secretion filter) |
| `shortlist_spherule_induced.tsv` | the superseded 40 — **see §5.1** |
| `candidates_annotated.tsv` | domain/secretion annotation of the 40 |
| `NOTES.md` | method detail and the v1 → v2 defect list |
| `01`–`06_*.py/.sh` | pipeline, numbered in run order |

## 10. Sources

**Literature**
- Hung CY, Yu JJ, Seshan KR, Reichard U, Cole GT. 2002. A parasitic phase-specific adhesin of *Coccidioides immitis* contributes to the virulence of this respiratory fungal pathogen. *Infect Immun* 70:3443-56. https://doi.org/10.1128/IAI.70.7.3443-3456.2002 — SOWgp binds laminin > fibronectin > collagen IV; deletion reduces ECM binding and virulence; 4–6 tandem repeats; parasitic-phase specific, elevated in early spherule development.
- Hung CY, Xue J, Cole GT. 2007. Virulence mechanisms of *Coccidioides*. *Ann N Y Acad Sci* 1111:225-35. https://doi.org/10.1196/annals.1406.020 — SOWgp depletion on endospores as immune evasion.

**Spherule/mycelium RNA-seq**
- Carlin AF, Beyhan S, Peña JF, Stajich JE, Viriyakosol S, Fierer J, Kirkland TN. 2021. Transcriptional analysis of *Coccidioides immitis* mycelia and spherules by RNA sequencing. *J Fungi (Basel)* 7(5):366. PMID 34067070. https://doi.org/10.3390/jof7050366 — **the source of the expression data used in §3.6/§4.5.** Reports young- and mature-spherule upregulation including a spherule surface protein, and overrepresentation of *Coccidioides*-unique genes among spherule-upregulated genes.
- Duttke SH, Beyhan S, Singh R, Neal S, Viriyakosol S, Fierer J, Kirkland TN, Stajich JE, Benner C, Carlin AF. 2022. Decoding transcription regulatory mechanisms associated with *Coccidioides immitis* phase transition using total RNA. *mSystems* 7(1):e0140421. PMID 35076277. https://doi.org/10.1128/msystems.01404-21 — csRNA-seq of the phase transition; identifies alternative promoter usage and a WOPR-family transcription factor (CIMG_02671) as critical for pathogenic growth. Relevant if the Tier-1 candidates are followed up for regulatory control.

**Databases** — UniProtKB, IEDB Query API, Pfam/InterPro, NCBI Datasets, SRA.

**Related internal documents** — `docs/model-review/2026-09-27-review-and-framework-plan.md`
(model evaluation behind §6), `docs/model-review/STATUS.md`, `docs/model-review/SEARCH-FRAMEWORK.md`.
