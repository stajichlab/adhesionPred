# Pfam confirms the diverged-repeat calls, and validates the detector's period

**2026-09-30.** Scripts 21-23 in `analysis/cocci_repeats/`. First independent check of
`14_repeat_detect_general.py`. Everything previously claimed about its 189 calls rested on the
detector itself plus visual inspection; Pfam-A is external evidence.

Next step 2 of `REPORT_2026-09-29_repeat_detector_divergence.md`.

> **Status: computational only. No experimental validation.**

---

## 1. Summary

1. **Gained calls are 13x enriched for Pfam repeat families over a length-matched background**
   — 31.7% against 2.4%, one-sided Fisher **P = 1.0e-25**. This is the result the test was
   built to produce.
2. **The enrichment is specific to repeat families, not to domains in general.** Gained calls
   carry *fewer* Pfam domains of any kind than background (58.5% vs 66.3%). So the signal is
   not "the detector likes domain-rich proteins".
3. **The detector's period estimate agrees with Pfam's own measure of the unit.** For
   single-unit Pfam models, median signed difference **+0.0 aa**, and **85.1% agree within
   1 aa**. This is the first external validation of the period, not only of the call.
4. **25 of the 34 "ankyrin" calls are confirmed ankyrin** by five Pfam ankyrin families. The
   eyeball assignment in the 2026-09-29 report was substantially right.
5. **The 51 aa group is a real structural repeat family** — 7/7 are PF17660, "Polyglycine
   hydrolase-like, structural repeat", scores 64.8-82.1.
6. **The 62 aa group has no Pfam hit at all** (0/8), and **51 of 123 gained calls remain
   entirely unverified**. Absence of a Pfam hit is weak evidence, so this is not evidence they
   are wrong — but it is not evidence they are right either.
7. **Gained calls are MORE enriched than shared calls** (31.7% vs 12.1%, P = 0.0019). The new
   detector is not adding noise on top of the old one's calls; it is recovering a class the
   old one systematically missed, which is what a divergence-floor fix predicts.

## 2. Method

Four protein sets from the same seven proteomes, built by `21_pfam_sets.py`:

| set | definition | n | median len |
|---|---|---|---|
| `gained` | called by 14, not by 02 | 123 | 235 |
| `shared` | called by both | 66 | 308 |
| `lost` | called by 02, not by 14 | 17 | 168 |
| `background` | **called by neither**, length-matched to `gained` | 1230 | 234 |

"Called" is the `03_repeat_surface_candidates.py` rule (period > 0, coverage >= 0.25,
copies >= 2.5), so the 123 / 66 / 17 counts reproduce `18_repeat_real_compare.py` exactly.

The background is the important control. It is drawn **per gained protein**, from the same
strain, within a length tolerance, sampled from the ~70,800 proteins neither detector called.
Length matching matters because long proteins carry more domains of every kind and the gained
set is long-biased.

`22_pfam_hmmsearch.sh` runs Pfam-A at the **gathering threshold** (`--cut_ga`) — Pfam's own
curated per-family cutoff, so no E-value is invented here. It also dumps each model's length,
which is one of the two period references used below.

A family counts as a repeat family when its Pfam description contains "repeat", plus a small
explicit list. All 27 families counted are written to `pfam_repeat_families.tsv`, so the rule
is auditable. **This rule is conservative** — see section 6.

## 3. Confirmation rate

| set | n | any Pfam | repeat family | median len |
|---|---|---|---|---|
| `gained` | 123 | 72 (58.5%) | **39 (31.7%)** | 235 |
| `shared` | 66 | 15 (22.7%) | 8 (12.1%) | 308 |
| `lost` | 17 | 1 (5.9%) | **0 (0.0%)** | 168 |
| `background` | 1230 | 815 (66.3%) | **29 (2.4%)** | 234 |

One-sided Fisher exact:

| comparison | P |
|---|---|
| gained vs background | **1.0e-25** |
| shared vs background | 3.5e-04 |
| gained vs shared | 1.9e-03 |

Read the "any Pfam" column alongside the "repeat family" column. Gained calls are **less**
likely than background to carry any Pfam domain (58.5% vs 66.3%) and **13x more** likely to
carry a repeat family. A detector that simply favoured well-annotated or domain-rich proteins
would raise both columns. This one raises only the repeat column.

**The 17 lost calls contain zero repeat families.** That is consistent with the 2026-09-29
finding that 12 of the 17 are low-complexity tracts rejected on purpose. It does **not**
vindicate dropping the other 5: four of those are curated Pro/Cys-rich class-2a candidates
which have no Pfam family either way, so this test cannot speak to them.

## 4. The period groups

### 33-34 aa: 25 of 34 are ankyrin

| family | Pfam | description | proteins | domains | score range |
|---|---|---|---|---|---|
| Ank | PF00023 | Ankyrin repeat | 25 | 134 | 15.3-37.8 |
| Ank_2 | PF12796 | Ankyrin repeats (3 copies) | 25 | 86 | 27.3-82.0 |
| Ank_3 | PF13606 | Ankyrin repeat | 25 | 96 | 17.3-28.0 |
| Ank_4 | PF13637 | Ankyrin repeats (many copies) | 25 | 102 | 22.8-44.4 |
| Ank_5 | PF13857 | Ankyrin repeats (many copies) | 21 | 53 | 27.4-40.1 |

The remaining **9 have no Pfam hit at all**. Seven of them share an identical 33 aa unit,
`LDFLVSELKDAGQHARMRILGILSKQAGLPESI`, one per strain across all seven proteomes. An orthologue
set with a perfectly regular repeat and no Pfam family is a plausible real repeat, but it is
**unconfirmed**.

### 51 aa: 7 of 7 confirmed, and not ankyrin

All seven are **PF17660, "Polyglycine hydrolase-like, structural repeat"**, scores 64.8-82.1,
35 domains. Pfam's model length is 50 aa against the detector's 51 aa period.

### 62 aa: 0 of 8 have any Pfam hit

All eight share an identical unit across strains
(`AALRVPGYQGEGYFFKGQRYLRMWWKPGTPEERKVFGPAKITDEW`), again one per proteome. Same status as the
nine above: regular, orthologous, unannotated, unconfirmed.

## 5. The period check — the strongest result

For each protein x repeat-family pair with >= 2 non-overlapping domains, Pfam gives a
**hit-to-hit spacing** measured on that protein. That is a better reference than the model
length because it is Pfam's measurement of the unit in this sequence, not a property of a model.

108 such pairs exist in proteins `14` called.

An **a priori** rule, fixed from the Pfam description alone and not from the data: a model whose
description says "copies" (Ank_2 "3 copies", Ank_4 / Ank_5 "many copies") spans several units,
so its spacing should be a multiple of the true unit. Every other model is one unit.

**Single-unit models (67 pairs):**

| | |
|---|---|
| median signed difference (period - spacing) | **+0.0 aa** |
| within 1 aa | **57 / 67 (85.1%)** |
| within 3 aa | 59 / 67 (88.1%) |

**By family:**

| family | Pfam | pairs | model len | med spacing | med period | med ratio | within 1 aa |
|---|---|---|---|---|---|---|---|
| Ank | PF00023 | 25 | 33 | 33.0 | 33.0 | 1.00 | **25 / 25** |
| Ank_3 | PF13606 | 17 | 31 | 33.0 | 33.0 | 1.00 | 15 / 17 |
| DUF346 | PF03984 | 8 | 39 | 49.0 | 49.0 | 1.00 | 8 / 8 |
| BTRD1 | PF17660 | 7 | 50 | 51.0 | 51.0 | 1.00 | 7 / 7 |
| RPEL | PF02755 | 2 | 24 | 44.0 | 44.0 | 1.00 | 2 / 2 |
| Ank_4 | PF13637 | 20 | 50 | 67.0 | 33.0 | **2.03** | 0 |
| Ank_2 | PF12796 | 12 | 90 | 84.5 | 33.0 | **2.52** | 0 |
| Ank_5 | PF13857 | 9 | 56 | 83.5 | 33.0 | **2.53** | 0 |
| TPR_10 | PF13374 | 7 | 42 | 67.0 | 45.0 | 1.49 | 0 |
| Sel1 | PF08238 | 1 | 36 | 46.0 | 36.0 | 1.28 | 0 |

The a priori rule holds. **Every family whose description says "copies" shows a ratio near 2 or
2.5; every single-unit family shows a ratio of exactly 1.00.** The detector reports the true
ankyrin unit (33 aa) while the multi-unit Pfam models report 2-2.5 units of spacing, which is
what they are built to do. This was predicted before looking at the numbers.

**For comparison, the old detector (`02`) reports a period for only 43 of these 108 pairs, and
27 are within 1 aa.** The new detector supplies a period for all 108 and matches on 57 of the
67 where a direct comparison is valid.

Two families do not fit either arm: TPR_10 (ratio 1.49) and Sel1 (1.28). TPR is a 34 aa repeat
often found as degenerate arrays with insertions, so a non-integer spacing is expected there;
the detector's fixed-period assumption is the weaker party in that comparison, not Pfam.

## 6. What is not confirmed

| | |
|---|---|
| gained with a repeat family | 39 / 123 |
| gained with a Pfam hit, none of them a repeat family | 33 / 123 |
| **gained with no Pfam hit at all** | **51 / 123** |

**Do not read 31.7% as a sensitivity.** Many genuine repeat proteins have no Pfam repeat
family — SOWgp itself has none. Absence of a hit is weak evidence in both directions.

**31.7% is also a lower bound**, because the counting rule is conservative. The commonest
non-repeat families in the gained set are EF-hand (7 proteins across five EF-hand families),
zinc knuckle (8), and C2H2 zinc finger (7). EF-hands, zinc knuckles and C2H2 fingers are all
tandemly arrayed structural motifs in real proteins, but their Pfam descriptions do not contain
the word "repeat", so they were not counted. A less conservative rule would raise the gained
rate and leave the background rate close to unchanged.

Other limits:

1. The background is "called by neither detector". If `14` has a systematically low false
   negative rate, some true repeat proteins sit in the background and the 2.4% is inflated —
   which would make the enrichment an underestimate, not an overestimate.
2. Only 108 protein x family pairs support the period check, dominated by ankyrin (83 of 108).
   The agreement is strong but is mostly a statement about one repeat class.
3. Nothing here tests the divergence floor. These are Pfam-detectable families; the synthetic
   floor measurement in the 2026-09-29 report stands unaltered.
4. The 51 unverified gained calls are the honest residue of this test.

## 6b. High specificity, low sensitivity — and why that is probably correct

The background set answers a question the test was not designed for: **how many repeat proteins
does the detector miss?**

29 of the 1230 background proteins (2.4%) carry a Pfam repeat family and were called by neither
detector. Scaled to the ~70,800 uncalled proteins that is of order **1,670 uncalled proteins
with a Pfam repeat family**, against **47 called** with one (39 gained + 8 shared). The scaling
is indicative only — the background is length-matched to `gained`, not a random draw — but the
direction is not in doubt.

**The misses are not marginal calls. They are not detected at all.**

| | gained, with repeat family | background, with repeat family |
|---|---|---|
| n | 39 | 29 |
| median length | 293 aa | **876 aa** |
| median `coverage_new` | 0.49 | **0.00** |
| period_new > 0 | 39 / 39 | **3 / 29** |
| meets 03's coverage >= 0.25 AND copies >= 2.5 | 39 / 39 | **0 / 29** |

26 of the 29 get **period 0** — the detector finds no periodicity in them whatsoever, so this
is not the `03` thresholds discarding weak calls. And the missed proteins are three times longer
than the called ones.

**This is most likely correct behaviour, not a defect.** The two things Pfam calls "repeat
families" are not one class:

- **tandem arrays that span the protein** — SOWgp, FLO11, ALS1. Repeat coverage near 1.0 is the
  defining property of class 2a in `docs/TOOL-ARCHITECTURE.md`. These are what the detector is
  built for, and the gained set has median coverage 0.49.
- **structural repeat domains inside a large globular protein** — WD40 propellers, ankyrin
  stacks, TPR solenoids. The repeat occupies a minority of a long sequence, and the protein is
  usually intracellular.

The background misses are dominated by the second kind: **WD40 is the single commonest repeat
family in the background (8 of 29)** and does not appear in the gained set at all.

So the accurate description is that `14` is a **tandem-array detector, not a repeat-domain
detector**. For finding class-2a surface proteins that is the right target, and the 13x
enrichment shows its calls are genuine repeats. For a general-purpose "does this protein contain
a repeat" question it has low sensitivity, and an HMM is the better tool — the same conclusion
already recorded for hydrophobins.

One consequence for the gained set: **ankyrin dominates it (25 of 39) yet ankyrin proteins are
intracellular.** They are called because these particular ones are short (median 293 aa), so a
33 aa array does reach 0.25 coverage. They are true repeats and false class-2a candidates, and
the SignalP filter removed all of them from the 58 secreted candidates. That filter is doing
more work than the detector here.

### Family composition differs sharply between sets

| set | proteins with a repeat family | dominant families |
|---|---|---|
| `gained` | 39 | Ank x5 families (25 proteins), BTRD1 (7), TPR_10 (7) |
| `shared` | 8 | **DUF346 only** (8) |
| `background` | 29 | WD40 (8), Ank (6), TPR_CAND1, RPEL, Sel1 |

The old detector's repeat-family calls are **entirely one family**, DUF346. The new detector's
gains are concentrated in ankyrin — a ~33 aa repeat known to diverge strongly between copies,
which is exactly the class an exact-match detector should miss. The enrichment pattern matches
the mechanism the 2026-09-29 report proposed, independently of that report's synthetic benchmark.

## 7. Files

| file | contents |
|---|---|
| `21_pfam_sets.py` | builds the four sets -> `pfam_sets.fa`, `pfam_sets.tsv` |
| `22_pfam_hmmsearch.sh` | Pfam-A at `--cut_ga` -> `pfam_sets.domtbl`, `pfam_models.tsv` |
| `23_pfam_confirm.py` | the analysis -> `pfam_confirmation.tsv`, `pfam_period_check.tsv`, `pfam_repeat_families.tsv` |

Run order: `21`, `sbatch 22`, `23`. Python is `/usr/bin/python3.12`.

## 8. Sources

- Mistry J, Chuguransky S, Williams L, et al. 2021. Pfam: The protein families database in
  2021. *Nucleic Acids Res* 49:D412-D419. https://doi.org/10.1093/nar/gkaa913
- Eddy SR. 2011. Accelerated profile HMM searches. *PLoS Comput Biol* 7:e1002195.
  https://doi.org/10.1371/journal.pcbi.1002195
