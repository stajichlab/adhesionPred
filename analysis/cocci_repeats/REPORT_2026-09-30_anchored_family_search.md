# Anchored family search: SOWgp's 92% prevalence is real, and most short models are real alleles

**2026-09-30.** Scripts 30-33 in `analysis/cocci_repeats/`. Generalises
`20_sowgp_anchor_search.py` into a motif-anchored family search with a derived anchor, a
measured specificity, and a classification of the failure modes that hid three SOWgp models.

> **Status: computational only. No experimental validation.**

---

## 1. Summary

1. **The anchor can be derived automatically.** From the SOWgp seed units the tool picks
   `KKYGDC` (conservation 0.79, **0 hits in a per-protein shuffled null**) and recovers exactly
   the same 10 proteins as the hand-picked `PTDCYGDC`.
2. **A shuffle null is not sufficient on its own** — a motif with zero chance hits still lands
   in unrelated real proteins. An alignment filter, not the anchor, is what defines membership.
3. **SOWgp's 92% prevalence is NOT an annotation artifact.** This was the hypothesis that
   motivated the run, and it is **wrong**. Anchored search gives 446/493 = 90.5%, the
   orthogroup approach 449/493 = 91.1%, and even the **union of both methods reaches only
   92.3%**. SOWgp still fails the >=95% antigen-tier filter.
4. **Most short SOWgp models are genuine short alleles, not truncations.** 383 of 474 member
   hits (81%) cover the reference completely with an internal deletion of a whole number of
   repeat units. Only 91 (19%) are model problems.
5. **26 split gene models exist across the pangenome**, GFF-confirmed. But the 80 bp gap seen
   in the two long-read strains **does not generalise** — observed gaps run 79-454 bp.

## 2. Deriving the anchor, and why the null was not enough

`30_anchor_family_search.py derive` takes the family's repeat units, enumerates k-mers
(k = 5-12), and ranks them by how often they occur in a **per-protein residue shuffle** of the
same proteomes. That shuffle keeps each protein's own composition and length, so a motif that is
merely a composition artifact (a Pro/Cys-rich string in a Pro/Cys-rich protein) hits the shuffled
copy about as often as the real one.

For SOWgp, from 42 units across 8 seed sequences (modal period 47 aa):

| motif | k | conservation | proteins hit | hits in shuffled null |
|---|---|---|---|---|
| **KKYGDC** | 6 | 0.786 | **10** | **0** |
| KKYGDCD | 7 | 0.786 | 10 | 0 |
| PPPPKKY | 7 | 0.762 | 10 | 0 |
| KYGDC | 5 | 0.786 | 10 | 2 |
| PPPPKK | 6 | 0.762 | **35** | 2 |

`KKYGDC` was chosen automatically and finds **the same 10 proteins** as the hand-picked
`PTDCYGDC` in `20_sowgp_anchor_search.py`. Note `PPPPKK` in the last row: a nearly as conserved
motif that hits 35 proteins. The specificity measurement is doing real work.

**The negative result that shaped the tool.** A zero-chance-hit motif is still not enough. The
anchor `LAAKIS`, derived for a different class-2a family, has zero null hits and still lands in
eleven 1589 aa proteins that align to their reference at 44.8% identity over 11% of it. The
shuffle destroys homology and convergent sequence, so it understates the false-positive rate for
an ordinary-composition motif. **Membership is therefore decided by the alignment**
(`MIN_IDENTITY = 60%` over ungapped aligned columns, `MIN_ALN_LEN = 30 aa`), with the anchor
used only to select candidates cheaply. All 10 SOWgp hits pass; the `LAAKIS` proteins do not.

## 3. The seven-proteome run reproduces the 2026-09-30 correction

| proteome | protein | len | anchors | ref span | id% | call |
|---|---|---|---|---|---|---|
| CiB10637 | `CIB10637_003943-T1` | 371 | 5 | 1-422 | 96.8 | short_allele, deletion **1.0 units** |
| CiB10992 | `CIB10992_003451-T1` | 371 | 5 | 1-422 | 96.8 | short_allele, **1.0 units** |
| Cpos1038 | `CPOS1038_003233-T1` | 164 | 2 | 1-164 | 95.7 | **split_candidate**, gap 80 bp |
| Cpos1038 | `CPOS1038_003234-T1` | 202 | 3 | 221-422 | 96.5 | **split_candidate**, gap 80 bp |
| Cpos3700 | `CPOS3700_010247-T1` | 115 | 1 | 1-117 | 98.2 | **split_candidate**, gap 80 bp |
| Cpos3700 | `CPOS3700_010248-T1` | 108 | 1 | 315-422 | 99.1 | **split_candidate**, gap 80 bp |
| Cpos3700 | `CPOS3700_003933-T1` | 202 | 3 | 221-422 | 95.5 | truncated (partner is on the haplotig) |
| VFC140 | `VFC140_004201-T1` | 234 | 2 | 1-422 | 99.1 | short_allele, deletion **4.0 units** |
| RS | `CIMG_04613-t26_1-p1` | 324 | 4 | 1-422 | 97.5 | short_allele, 1.7 units |
| Silveira | `QVM09276.1` | 328 | 4 | 1-422 | 99.1 | short_allele, **2.0 units** |

10/10 pass the member filter; 7/7 proteomes have a hit. The tool independently reproduces the
manual 2026-09-30 correction, and adds that **VFC140's internal deletion is exactly 4.0 repeat
units** — a whole number, which is what a unit-count allele should be and not what a random
truncation would give.

Correction to an earlier statement: the intergenic gap is **80 bp**, not 81. The earlier figure
was the difference of the two GFF coordinates; the gap between the genes is 80.

## 4. The pangenome: prevalence is not rescued

`31_anchor_pangenome.sh` over 493 pangenome proteomes. 474 member hits in **446 proteomes**.

| method | strains with SOWgp | prevalence |
|---|---|---|
| anchored search (this work) | 446 | **90.5%** |
| orthogroup, `05_sowgp_pangenome.py` | 449 | 91.1% |
| **union of both** | **455** | **92.3%** |

Agreement is high (440 strains found by both); each method finds a handful the other misses
(6 anchored-only, 9 orthogroup-only). Neither is complete, and **the union still does not reach
95%**.

**This refutes the hypothesis the run was designed to test.** The reasoning was: SOWgp is
excluded from the *Coccidioides* antigen Tier 1/Tier 2 lists
(`docs/reports/2026-09-27-coccidioides-antigen-findings.md`) solely because its prevalence is
0.9201, below the >=0.95 filter; the 2026-09-30 finding that SOWgp models are split, truncated
and missed by orthogroup assignment suggested that 92% was an annotation artifact. **It is not.**
Two independent detection methods and their union all land at 90-92%.

What remains is **assembly**, not annotation. A protein-level search cannot find a locus that was
never assembled, and tandem arrays are what short-read assemblies collapse. The 2026-09-27
report's own recommendation stands unaltered and is now the only route left: confirm SOWgp
presence and copy number **from read depth** across the WGS runs, not from assembled proteomes.

**SOWgp's exclusion from the antigen tiers should be treated as correct until read-depth data
says otherwise.**

## 5. Most short models are real alleles, not fragments

This corrects the framing in `README.md` section 3 and the 2026-09-29 report, which recorded
that "only 205 of the 489 SOWgp models are >= 250 aa. The rest are fragments (72-249 aa).
Short-read assemblies and gene models can truncate the repeat array."

Classifying the 474 member hits by whether they cover the reference completely:

| call | n | % | meaning |
|---|---|---|---|
| `full_length` | 107 | 23% | complete, 6 units |
| `short_allele` | 276 | 58% | **complete reference coverage**, internal deletion of a whole number of units |
| `truncated` | 65 | 14% | partial span, no partner — a real model problem |
| `split_candidate` / `split_model` | 26 | 5% | one gene called as two |

**383 of 474 (81%) are genuine alleles.** Only 91 (19%) are model problems.

The unit spectrum among complete-coverage alleles:

| units | strains |
|---|---|
| 6 | 107 |
| 5 | 97 |
| 4 | 164 |
| 3 | 15 |

And **205 of those 383 are >= 250 aa** — the same 205 the README reports. The number is right;
its interpretation was wrong. The 178 complete alleles below 250 aa are **3- and 4-unit alleles**,
not truncations. The array-length variation is largely real biology, not assembly damage.

This does not overturn the assembly caveat for *copy number* — a collapsed array would look like
a short allele here, and this method cannot tell those apart. It does mean the "most models are
fragments" reading is not supported.

## 6. Split gene models across the pangenome

26 GFF-confirmed split pairs. The 80 bp intergenic gap found in Cpos1038 and Cpos3700
**does not generalise**: observed gaps are 79, 80, 110, 112, 116, 203, 241, 257, 277 and 454 bp.
The earlier observation that two strains shared an identical gap was **n = 2 and coincidental**.
What generalises is the pattern — consecutive locus tags, collinear on one contig, ref spans that
overlap or abut, pair coverage 1.0 — not the gap size.

26 proteins are about 13 genes. That is enough to matter for per-strain copy-number and length
statistics, and not enough to move prevalence.

## 7. Limits

1. **The anchored search has its own blind spot.** `KKYGDC` sits inside the repeat unit, so a
   model that lost the entire array is invisible to it. That is exactly the failure mode most
   likely to affect the 47 strains with no hit, so 90.5% is a floor, not an estimate.
2. **`short_allele` cannot distinguish a real short allele from a collapsed assembly.** Both
   present as complete reference coverage with a whole number of units missing. Only read depth
   separates them.
3. The unit-count spectrum in section 5 is measured against the 6-unit SOWgp82 reference. A
   different reference would shift every count by a constant.
4. Only SOWgp was run at pangenome scale. The class-2a sweep (`32`/`33`) over the other families
   is written but not run here.
5. Species is not separated in section 5; *immitis* and *posadasii* have different modal alleles,
   so the spectrum mixes two distributions.

## 8. Files

| file | contents |
|---|---|
| `30_anchor_family_search.py` | the tool: `derive`, `search`, `family` subcommands |
| `31_anchor_pangenome.sh` | 493-proteome SLURM run -> `sowgp_anchored_pangenome.tsv` |
| `32_class2a_anchor_sweep.py`, `33_anchor_sweep.sh` | class-2a family sweep (written, **not run**) |
| `sowgp_anchor_family.tsv` | the 7-proteome run |
| `sowgp_anchored_pangenome.tsv` | 474 member hits across 446 proteomes |

Run: `./30_anchor_family_search.py family --seed sowgp_seed.fa --proteomes longread --out
sowgp_anchor_family.tsv`, then `sbatch 31_anchor_pangenome.sh`.

---

## 9. The class-2a sweep (scripts 32/33)

Applied to every class-2a candidate family, not only SOWgp. 63 families from the 41 + 58
candidate sets; **40 have an anchor** with <= 5 chance hits; 12 distinct anchors after
collapsing duplicates.

### 9.1 Which families can be anchored at all

| composition class | families with an anchor | families with none |
|---|---|---|
| Pro/Cys-rich (SOWgp / BAD1 type) | 18 | 14 |
| other | 15 | 9 |
| **Ser/Thr-rich (FLO / ALS type)** | **2** | 0 |

Only 2 Ser/Thr-rich families needed an anchor and both got one (`TTECEE` conservation 1.00,
`TTTEQP` 0.90), which is a better outcome than expected — a Ser/Thr array was the case where a
composition-based rule was predicted to fail. 23 of 63 families have no usable anchor at all,
and those are mostly Pro/Cys-rich, where the repeat is too degenerate to yield a specific k-mer.
**An anchored search is not a general substitute for a detector; it works for 2 families in 3.**

### 9.2 The gain is real but much smaller than the raw counts suggest

124 hits over 86 distinct proteins. **28 proteins are rejected by the alignment filter as
`unrelated`**, leaving **58 members**.

| | members only | all hits (misleading) |
|---|---|---|
| anchored proteins the old detector (02) does not call | **24** | 52 |
| anchored proteins the new detector (14) does not call | **14** | 42 |

The first pass of this script printed the right-hand column. It is wrong to quote: a hit that
fails the alignment filter is not a family member, so counting it as "found by anchoring,
missed by the detector" inflates the gain three-fold. The script now prints the member-only
totals; the per-family `n_missed_by_*` columns remain all-hit counts and must be read together
with the `unrelated` column.

**The alignment filter earns its place here.** The `LAAKIS` family has 25 hits and **19 are
rejected as unrelated** — the exact false positive documented in section 2. Without the filter
that one family alone would have contributed 21 spurious "novel" proteins.

### 9.3 What the 14 genuinely new members are

| what | n | detail |
|---|---|---|
| **a 9 aa-period Pro/Cys family** (anchor `IPTEWP`) | 6 | `CIMG_04070`, `CPOS3700_003468`, `VFC140_003739`, `QVM08828.1`, `CIB10637_003487`, `CIB10992_003902`. 142-163 aa, 2-4 anchored units, 98.6-100% identity. **One orthologue per proteome, in all 7.** Neither detector calls any of them. |
| **a 17 aa Ser/Thr family** | 3 | `CIMG_03231`, `CPOS3700_013961`, `CPOS3700_002743`. 330-372 aa, 4-5 units. FLO/ALS-type. |
| **SOWgp split fragments** | 3 | `CPOS1038_003233`, `CPOS3700_010247`, `CPOS3700_010248`. Each carries 1-2 units, below `03`'s `copies >= 2.5`, so no periodicity detector can reach them. |
| borderline "other" 43 aa hits | 2 | `CIB10637_000756`, `CIB10992_000764`, 120 aa, identity 70-71% — just above the 60% member threshold. **Weakest of the set; treat as unresolved.** |

The 9 aa Pro/Cys family is the most interesting result of the sweep. It is a complete orthologue
set, present once in every one of the seven proteomes at 98.6-100% identity, with a regular
short-period Pro/Cys repeat, and **neither periodicity detector calls it**. Short proteins with
a 9 aa period do not reach the coverage and copy-number thresholds `03` applies.

**All six carry a signal peptide.** Joined against the 2026-09-30 SignalP run, every one is
called `SP` with an identical cleavage site (19-20) and an identical probability (0.9869):

| protein | proteome | len | units | SignalP |
|---|---|---|---|---|
| `CIMG_04070-t26_1-p1` | *C. immitis* RS | 151 | 3 | SP, CS 19-20, 0.987 |
| `CPOS3700_003468-T1` | Cpos3700 | 152 | 4 | SP, CS 19-20, 0.987 |
| `VFC140_003739-T1` | VFC140 | 152 | 4 | SP, CS 19-20, 0.987 |
| `QVM08828.1` | Silveira | 163 | 4 | SP, CS 19-20, 0.987 |
| `CIB10637_003487-T1` | CiB10637 | 142 | 2 | SP, CS 19-20, 0.987 |
| `CIB10992_003902-T1` | CiB10992 | 142 | 2 | SP, CS 19-20, 0.987 |

So this is a **secreted, Pro/Cys-rich, tandem-repeat protein family, universal across the seven
genomes, with variable unit count (2-4), that no detector in this project has ever called.** By
the class-2a definition in `docs/TOOL-ARCHITECTURE.md` it is a class-2a candidate, and it is a
new one. Its unit count varies between strains in the same way SOWgp's does, which makes it a
candidate for the same antigenic-variation question.

Nothing is known here about its function, expression or antigenicity. It has not been checked
against the spherule RNA-seq, the pangenome, or any confounder genome.

The 17 aa Ser/Thr family is more mixed: `CIMG_03231` and `CPOS3700_002743` are `SP` but at a
lower probability (0.64), and `CPOS3700_013961` is `OTHER`.

The SOWgp split fragments confirm the tool's design premise directly: **an anchored search
recovers gene-model failures that a periodicity detector cannot**, because it needs one motif
occurrence rather than a detectable period.

### 9.4 Limits of the sweep

1. Families were defined from candidate proteins, not from curated families, so "family" here
   means "this protein and things that align to it".
2. 23 of 63 families have no anchor and are invisible to this method. They are not reported as
   absent, they are untested.
3. The two borderline 43 aa hits sit close to the identity threshold; a stricter cut would drop
   them and the member count would be 12, not 14.
4. Only the 7 long-read/reference proteomes were swept. The pangenome sweep was run for SOWgp
   only.
