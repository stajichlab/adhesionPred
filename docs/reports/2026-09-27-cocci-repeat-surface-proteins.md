# Class 2a curation in *Coccidioides*: repeat surface proteins from long-read genomes

**Working report — 2026-09-27.** Stajich lab, UC Riverside.
Companion to `docs/TOOL-ARCHITECTURE.md` (class 2a = repeat/avidity surface proteins) and
`docs/reports/2026-09-27-coccidioides-antigen-findings.md`.

> **Status: provisional, no experimental validation.** Reproducible from
> `analysis/cocci_repeats/`. Read §5 before using the candidate list.

---

## 1. Why this was done

The class-2a model (repeat/avidity surface proteins) works within Saccharomycotina but misses
SOWgp. A transfer test (`analysis/model_review/repeat_structure_transfer.py`) showed
composition-agnostic repeat features help only marginally (2/5 → 3/5 recovered), and suggested
the real gap is that **training positives are long and Ser/Thr-rich (FLO11 1367 aa, ALS1 1260 aa)
while SOWgp is short and Pro/Cys-rich (324–422 aa)**. Length and composition are confounded
with clade in the current label set.

This is the curation step: find more repeat surface proteins in *Coccidioides* directly, using
**long-read assemblies**, where tandem arrays survive.

## 2. Data and method

| input | detail |
|---|---|
| long-read strains | 5 UArizona assemblies — *C. immitis* CiB10637, CiB10992, VFC140; *C. posadasii* Cpos1038, Cpos3700 (`shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc`) |
| references | *C. immitis* RS and *C. posadasii* Silveira 2022 (FungiDB, from the pangenome) |
| repeat detection | periodicity-based and **composition-agnostic**: for each period *p*, score positions where s[i] == s[i+p], take the best *p*, delimit the region sustaining it. Reports period, copy number, extent and the repeat unit. |
| secretion | **SignalP 6** (GPU build, `short_gpu`) — 436–686 secreted proteins per proteome |
| composition | Ser+Thr, Pro, Cys, Gly+Ala, charged, top-3 residue fraction |

**Why long reads matter:** tandem arrays collapse or fragment in short-read assemblies, and
copy number is precisely the feature being measured. Long-read strains yielded ~9 repeat
proteins each against ~6 per short-read reference.

## 3. Validation — the detector recovers published repeat biology

Run blind on known proteins, it reproduces the literature:

| protein | detected period | detected copies | published |
|---|---|---|---|
| SOWgp58 | **47** | 3.8 | 41–47 aa repeat, 4 copies (Hung et al. 2002) |
| SOWgp82 | **47** | 5.8 | 6 copies |
| BAD1 | **24** | 26.9 | its characterized 24-aa tandem repeat |
| ALS1 | **36** | 9.9 | Als 36-aa tandem repeat |
| Ag2/PRA, CF antigen, RodA | — | 0 | correctly, these are not repeat proteins |

## 4. Findings

### 4.1 41 class-2a candidates

71,044 proteins profiled → 58 with a substantial tandem repeat → **41 that are also secreted**
(`analysis/cocci_repeats/class2a_candidates.tsv`).

| composition class | n | exemplar |
|---|---|---|
| **Pro/Cys-rich (SOWgp / BAD1 type)** | **20** | SOWgp itself |
| Ser/Thr-rich (FLO / ALS type) | 10 | the period-17 family below |
| other | 11 | — |

**Pro/Cys-rich repeats outnumber Ser/Thr-rich 2:1 in *Coccidioides*.** This is the concrete
reason a FLO/ALS-trained class-2a model fails here, and these 20 are the training examples the
model was missing.

### 4.2 SOWgp copy-number variation, measured directly

| genome | protein | length | repeat copies |
|---|---|---|---|
| *C. immitis* RS | CIMG_04613 | 324 aa | **3.4** |
| *C. immitis* CiB10637 | CIB10637_003943 | 371 aa | **4.4** |
| *C. immitis* CiB10992 | CIB10992_003451 | 371 aa | **4.4** |
| *C. posadasii* Silveira | QVM09276.1 | 328 aa | **3.8** |

This is the repeat-length variation *measured on assemblies that can represent it*, rather
than inferred from allele sizes. It is consistent with the published 58/66/82 kDa size series.

> **Not detected in VFC140, Cpos1038 or Cpos3700.** Do **not** read this as absence. It may be
> a real deletion, a repeat array below the detection threshold (coverage ≥ 0.25, ≥ 2.5
> copies), a collapsed array, or a gene-model failure. Resolving this requires read-depth
> analysis at the locus (§6). Given SOWgp's 92% pangenome prevalence reported in the antigen
> study, some genuine absence is plausible — but it is not established here.

### 4.3 A Ser/Thr-rich period-17 family, with copy-number variation

10 members across 6 strains, 423–462 aa, ~31% Ser+Thr, 7–9 repeat copies — i.e. genuinely
FLO/ALS-like in composition. Distribution is uneven:

| strain group | members |
|---|---|
| *C. posadasii* Cpos1038 | **5** |
| *C. immitis* CiB10637 / CiB10992 | 2 |
| *C. immitis* VFC140 | 2 |
| *C. posadasii* Silveira | 1 |

Five copies in one strain against one in another is an expanded family with apparent
copy-number variation, and it is the closest thing *Coccidioides* has to a classical
Ser/Thr-rich adhesin family. It has not been characterized here and deserves its own look.

### 4.4 Other recurring signals

- A **period-76, 305 aa protein at coverage 1.00** appears in five of seven genomes.
- A **period-34, Pro-rich (~27%) family** recurs across strains at 355–534 aa.

Neither is identified; both are single-family leads worth a targeted look.

## 5. Limitations

1. **Secretion is SignalP-only.** No GPI-anchor prediction was run, so genuinely cell-wall
   anchored proteins without a classical signal peptide are missed. NetGPI/PredGPI is the
   obvious addition.
2. **Absence is not established** (§4.2) — repeat arrays collapse, gene models fail, and the
   detection threshold is arbitrary.
3. **Thresholds are arbitrary**: coverage ≥ 0.25 and ≥ 2.5 copies were chosen to recover the
   known positives, not derived from anything.
4. **Degenerate repeats are under-detected.** FLO11 itself scores only 0.045 coverage under
   this detector because its repeats are diverged — so this method finds *regular* arrays and
   will miss diverged ones. Real *Coccidioides* adhesins with degenerate repeats are likely
   being missed.
5. **No functional evidence** for any candidate. These are training examples and hypotheses.
6. Gene models differ between the long-read annotations and the FungiDB references, so
   cross-strain comparisons of counts carry annotation noise.

## 6. Next steps

1. **Feed the 20 Pro/Cys-rich candidates back into the class-2a model** as training positives
   and re-run the transfer test. This is the direct test of whether curation closes the
   SOWgp gap identified in `repeat_structure_transfer.py`.
2. **Read-depth at the SOWgp locus** across the strains lacking a call, and across the 402
   SRA WGS runs, to separate real absence from assembly and annotation artifact.
3. **Add NetGPI** to the surface filter (limitation 1).
4. **Identify the period-17 and period-76 families** — BLAST/HMM against characterized
   proteins, check expression and pangenome prevalence.
5. **Extend to Onygenales** — *Histoplasma*, *Blastomyces*, *Paracoccidioides* long-read
   genomes, to test whether the Pro/Cys-rich repeat class is an Onygenales-wide feature
   (BAD1 suggests it is).

## 7. Files

`analysis/cocci_repeats/`:

| file | contents |
|---|---|
| `01_signalp.sh` | SignalP 6 (GPU) over long-read + reference proteomes |
| `02_repeat_profile.py` | periodicity-based repeat detection + composition |
| `03_repeat_surface_candidates.py` | join to class-2a candidates |
| `class2a_candidates.tsv` | **the 41 candidates** |
| `repeat_profile_{longread,reference}.tsv` | full per-protein profiles (untracked, regenerable) |

## 8. Sources

- Hung CY, Yu JJ, Seshan KR, Reichard U, Cole GT. 2002. *Infect Immun* 70:3443-56. https://doi.org/10.1128/IAI.70.7.3443-3456.2002 — SOWgp tandem repeats and size alleles.
- Teufel F, et al. 2022. SignalP 6.0. *Nat Biotechnol*. https://doi.org/10.1038/s41587-021-01156-3
- Long-read assemblies: `shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc`
