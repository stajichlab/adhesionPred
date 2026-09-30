# Are there other proteins with SOWgp-like architecture?

**2026-09-30.** Scripts 35-39 in `analysis/cocci_repeats/` (written by a subagent that was
stopped before it produced results; the pipeline was run directly afterwards).

The question is one of protein architecture: does the *Coccidioides* genome, or any related
genome, contain other proteins built like SOWgp? SOWgp is a secreted protein with a signal
peptide, a Pro/Cys-rich tandem array of 47 aa units, and an unrepeated N-terminus.

> **Status: computational only. No experimental validation.**
>
> **Provisional.** Search job 29300347 was still running when this was written. The numbers
> below come from output files stamped 15:08. The BAD1 control (script 38) and the copy-number
> summary (script 39) have not run. Re-check section 3 against the final files.

---

## 1. Summary

1. **The SOWgp 47 aa unit has no detectable relative outside *Coccidioides*.** A profile HMM
   built from 34 SOWgp units scored 2,055 domain hits across 6,313 proteomes. All hits at
   score >= 26 (1,293 of them) are in *Coccidioides*. The best hit anywhere else scores 19.7.
2. **Inside *Coccidioides*, nothing else has SOWgp's composition.** Among 26 distinct repeat
   families, only the SOWgp orthologs lie close to SOWgp in composition (distance 1.1-3.1). The
   next family is at 8.6.
3. **Two families are partial matches.** A period-43 family has BAD1-like composition. A
   period-34 family is Pro-rich with many copies (9-10). Neither shares the SOWgp unit.
4. **Across Onygenales, repeat proteins with Pro/Cys-rich composition exist, but as many
   separate families.** 74 proteins beyond SOWgp fall into 27 different repeat periods. They
   are architectural look-alikes, not unit-level relatives.
5. **The test is coarse.** It compares composition, not layout. See section 6.

## 2. The HMM search

`sowgp_unit.hmm` is a profile built from 34 aligned SOWgp units (model length 47 aa).
`unit_genomes.tsv` lists 6,313 proteomes: 7 long-read *Coccidioides*, 493 pangenome strains
and 5,813 Fungi_5k proteomes. Of these, 571 are Onygenales (500 *Coccidioides*, 71 other).

### 2.1 Calibration

Two questions were answered before the search was trusted.

**Null.** The HMM was run on proteomes shuffled per protein (length and composition kept,
order destroyed). The null table is empty: zero hits.

**Sensitivity.** Arrays of four units were built from the modal unit at fixed identity and
embedded in random flanks, 300 arrays per cell. Natural and Pro-rich flanks behaved the same.

| unit identity | arrays with >= 1 unit detected | all copies found | median best score |
|---|---|---|---|
| 60% | 100% | 100% | 50.6 |
| 50% | 100% | 95% | 40.0 |
| 40% | 92% | 6% | 26.0 |
| 30% | 1% | 0% | 11.6 |

The floor is sharp. No indels were simulated, so real diverged units would fall off sooner.

### 2.2 Result

| | hits | min score | median | max | score >= 26 | score >= 40 |
|---|---|---|---|---|---|---|
| *Coccidioides* | 1,748 | 18.3 | 68.5 | 99.3 | 1,293 | 1,293 |
| all other genomes | 307 | 11.4 | 13.0 | 19.7 | 0 | 0 |

The two groups do not overlap. No hit scores between 26 and 40 in either group. The 307
non-*Coccidioides* hits are singletons spread across unrelated orders (Pleosporales,
Glomerellales, Eurotiales, Boletales and others), which is the pattern expected from noise.

**What this supports:** no protein in 5,813 Fungi_5k proteomes carries SOWgp-like units at
about 40% identity or better. **What it cannot support:** a relative more diverged than that
would be invisible. The negative has a stated reach, not an unlimited one.

## 3. Architecture without unit homology

The third search ran `14_repeat_detect_general.py` on the Onygenales proteomes, giving 644,399
protein rows. Of these, 1,457 carry a repeat call (coverage >= 0.25, copies >= 2.5).

Composition filter for SOWgp-like: Pro + Cys > 15% **and Cys >= 4%**. The cysteine floor is
deliberate. `03_repeat_surface_candidates.py` uses a sum without a cysteine floor, which lets
cysteine-free proteins through (see `REPORT_2026-09-30_ptgiptewp_family.md`, section 2).

83 proteins pass, in 50 proteomes. Removing SOWgp itself (period 47, *Coccidioides*) leaves 74.

| period | n | genera | median length | median %Cys | median %Pro |
|---|---|---|---|---|---|
| 34 | 11 | *Coccidioides* 9, *Nannizziopsis* 2 | 501 | 4.8 | 22.7 |
| 40 | 10 | *Trichophyton* 5, *Microsporum* 3, *Blastomyces* 1 | 251 | 4.0 | 12.7 |
| 41 | 6 | *Trichophyton* 3, *Blastomyces* 1, *Arthroderma* 1 | 324 | 10.4 | 5.7 |
| 80 | 5 | mixed | 273 | 17.6 | 7.3 |
| 27 | 5 | *Coccidioides* 5 | 199 | 12.1 | 3.5 |
| 55 | 4 | *Blastomyces* 2, *Emergomyces* 1, *Emmonsia* 1 | 382 | 7.0 | 10.9 |
| 5 | 4 | *Coccidioides* 3, *Histoplasma* 1 | 86 | 14.0 | 20.9 |
| 56 | 3 | *Emergomyces* 2, *Histoplasma* 1 | 359 | 7.1 | 12.6 |
| 25 | 3 | *Microsporum* 3 | 339 | 15.3 | 5.9 |

The table lists the nine largest of 27 periods. The rest have one or two members each.

**Reading.** The 74 proteins split across 27 different repeat periods, none at 47 aa. The HMM
found no SOWgp unit in them. They share a composition class with SOWgp but not its repeat. They
are therefore convergent in architecture, or diverged beyond the HMM floor. This analysis
cannot tell those two apart.

## 4. Inside *Coccidioides*: closest families to SOWgp and BAD1

The 58 secreted repeat candidates (`class2a_candidates_general.tsv`) collapse to 26 distinct
repeat families. Each was ranked by Euclidean distance in (%Ser+Thr, %Pro, %Cys) to SOWgp
(11.4, 16.6, 6.6) and BAD1 (8.5, 3.8, 8.0). These reference values come from
`docs/TOOL-ARCHITECTURE.md`; SOWgp is not recomputed here. BAD1 (A4D962) was recomputed from
sequence for the classifier report and agrees.

| family (example protein) | period | len | %S+T | %Pro | %Cys | copies | dist SOWgp | dist BAD1 |
|---|---|---|---|---|---|---|---|---|
| SOWgp orthologs (`CIMG_04613`) | 47 | 324 | 10.8 | 15.7 | 6.5 | 4.0 | **1.1** | 12.2 |
| SOWgp, VFC140 allele (`VFC140_004201`) | 47 | 234 | 11.5 | 13.7 | 5.6 | 3.0 | 3.1 | 10.6 |
| period-43 (`CPOS3700_005551`) | 43 | 274 | 5.8 | 5.5 | 6.6 | 4.0 | 12.4 | **3.5** |
| period-43 (`VFC140_000791`) | 43 | 281 | 6.0 | 5.3 | 5.7 | 5.0 | 12.6 | 3.7 |
| period-43 (`CPOS1038_001811`) | 43 | 253 | 6.3 | 5.5 | 5.1 | 4.0 | 12.3 | 4.0 |
| period-39 (`CIMG_01182`) | 39 | 216 | 14.8 | 5.6 | 3.2 | 4.0 | 12.0 | 8.1 |
| period-34 (`CIMG_00195`) | 34 | 501 | 18.0 | 21.8 | 4.8 | 9.0 | 8.6 | 20.6 |
| period-7 (`CIMG_07912`) | 7 | 105 | 17.1 | 20.0 | 1.0 | 7.0 | 8.7 | 19.6 |

- **Only the SOWgp orthologs are close to SOWgp.** The gap to the next family (8.6) is large.
- **The period-43 family resembles BAD1.** Low Pro, moderate Cys, four to five copies. The
  current classifier files these proteins under `other`. BAD1 itself is filed there too.
- **The period-34 family shares the Pro/Cys character of SOWgp** but is about 500 aa with 9-10
  copies, against SOWgp's 324 aa and 4 copies.
- The remaining families are short Pro/Thr repeats, or low-complexity arrays (`GGHKH`).

## 5. Answer to the question

Among the proteins this project has examined, **SOWgp-like architecture is not repeated inside
the *Coccidioides* genome**. Its closest relatives there are partial: one family that matches
BAD1 better than SOWgp, and one Pro-rich family with a longer array. Across Onygenales,
proteins with similar composition and tandem repeats exist in several genera, but they use
different repeat periods and show no unit homology to SOWgp.

## 6. Limits

1. **The candidate pool is restricted.** Section 4 uses proteins that passed the detector
   thresholds and carry a signal peptide. It is not a whole-proteome search. Proteins with short
   arrays or no predicted signal peptide are not covered.
2. **Composition distance is crude.** It does not test layout: position of the signal peptide,
   length of the non-repeat N-terminus, or the ordering of domains. Two proteins can match in
   composition and differ in organisation.
3. **The 83 and 74 counts are composition filters.** They are not validated family assignments.
4. **The BAD1 control has not run.** Without it, the method has not been shown to separate
   architecture-level relatedness from unit-level homology. It has only been shown to find no
   unit homology.
5. **The HMM floor is optimistic.** The calibration simulated no indels. Real diverged units
   would be lost at higher identity than 40%.
6. **The search job was unfinished** at the time of writing (section status note).
7. **Number reuse.** `34_sowgp_unit_hmm.py` (builds `sowgp_unit.hmm`) and
   `34_newfam_pangenome.sh` share the number 34. The file names differ, so nothing is
   overwritten, but the numbering is ambiguous. An earlier draft of this report wrongly said the
   HMM builder did not exist; it does, and it is in the repository.

## 7. Suggested next step

Compare layout directly. For SOWgp and each candidate family, record signal-peptide length,
N-terminal non-repeat length, repeat start, repeat end and C-terminal length, then test whether
those five numbers cluster. That would show whether any family shares SOWgp's organisation,
rather than only its composition.

## 8. Files

| file | contents |
|---|---|
| `34_sowgp_unit_hmm.py` | builds `sowgp_unit.hmm` from 34 aligned SOWgp units |
| `35_unit_genomes.py` | builds `unit_genomes.tsv`, 6,313 proteomes with taxonomy |
| `36_unit_calibrate.py`, `.sh` | null and sensitivity calibration of the unit HMM |
| `37_unit_search.sh` | HMM, anchor and architecture searches over all proteomes |
| `38_unit_bad1_control.py`, `.sh` | BAD1 control (not yet run) |
| `39_unit_distribution.py`, `.sh` | copy-number summary (not yet run) |
| `sowgp_unit.hmm` | profile HMM, 34 units, model length 47 |
| `unit_hmm_hits.tsv.gz`, `unit_anchor_hits.tsv.gz`, `unit_architecture_onygenales.tsv.gz` | search output |
| `unit_calibration_null.tsv`, `unit_calibration_sensitivity.tsv` | calibration output |
