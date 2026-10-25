# AGA1-class flags and classifier on five yeast proteomes

Date: 2026-10-25. Record, inputs and hashes: `docs/reports/data/aga1_proteome_test_2026-10-25/`. Scripts: PR 147 (staging and the HPCC
array), PR 148 (classifier), and `merge_proteome_chunks.py`, `proteome_summary.py` here.

## What was done

SignalP 6, PredGPI and AIUPred ran on the UniProt reference proteomes of *S. cerevisiae* S288C, *C. albicans* SC5314,
*C. glabrata* CBS138, *K. marxianus* (strain FIM1; the DMKU3-1042 proteome has no UniProtKB entries) and *C. auris*, 27,811
proteins in 32 chunks (SLURM arrays 30483377, 30483378, 30483379; all 96 tasks completed on the pinned commit). The project's
tandem-repeat detector ran locally. Then per protein: the broad rule (signal peptide and GPI and a disordered run of 200 or more
residues), the narrow rule (broad and the repeat call), and the three frozen classifiers of PR 148 (logistic regression,
random forest, gradient boosting, fit on all tier A positives and negatives of the 850 proteins). A classifier call is a score at or
above the score where the 850 negatives (out of fold) have a 3.5% false-positive rate.

## Counts per species

| Species | Proteins | Signal peptide | GPI | Broad rule | Narrow rule | Called by 2 of 3 models | Of the 850 present |
|---|---|---|---|---|---|---|---|
| *S. cerevisiae* | 6,067 | 297 | 82 | 11 | 5 | 17 | 298 |
| *C. albicans* | 6,041 | 365 | 135 | 25 | 6 | 35 | 239 |
| *C. glabrata* | 5,200 | 245 | 117 | 25 | 4 | 48 | 21 |
| *K. marxianus* | 5,070 | 212 | 80 | 10 | 2 | 28 | 0 |
| *C. auris* | 5,433 | 269 | 108 | 14 | 5 | 19 | 14 |

*S. cerevisiae* and *C. albicans* are mostly the species the 850 came from, so their counts include training proteins: no *S. cerevisiae*
protein passes the broad rule that is not already among the 850 (0 novel), which says the known broad class of S288C is covered.
*C. glabrata* (21 of the 850 present, none tier A) and *K. marxianus* (none) are the real tests.

## Does it find the annotated adhesin-like proteins? (rough, not independent)

Reference: UniProt protein names that mark an adhesin-like wall protein (PA14, adhesin, flocculation, agglutinin-like,
Hyphally-regulated, Flo11, Als, Hwp, Iff, mucin). Names of unreviewed entries are often inferred from homology, so this is a rough
guide, not a truth set. A flag here is the broad rule, or at least 2 of 3 models together with a signal peptide.

| Species | Adhesin-named | Found by the broad rule | Found by broad rule or 2 of 3 models + SP | Other proteins called |
|---|---|---|---|---|
| *C. glabrata* | 37 | 17 | **32** | 15 of 5,163 |
| *K. marxianus* | 14 | 4 | **11** | 12 of 5,056 |
| *C. auris* | 16 | 11 | 13 | 6 of 5,417 |
| *C. albicans* | 24 | 14 | 17 | 20 of 6,017 |
| *S. cerevisiae* | 11 | 6 | 8 | 7 of 6,056 |

In *C. glabrata* the broad-rule-or-model route recovers most of the EPA family (EPA1, EPA3, EPA6, EPA8, EPA9, EPA10, EPA15, EPA23),
AWP1, AWP2, AWP3b, PWP5 and PWP6 and many of the PA14-domain and Hyr/Iff-like proteins, but not EPA2 or EPA20 (EPA2 has no signal
peptide call); in *K. marxianus* it recovers FLO1, FLO5, FLO9 and
FLO11-like (MUC1) proteins (11 of 14 named; one FLO9 and two FLO5 entries are missed). It does this from features alone, with no sequence search, in species whose adhesins the
training set barely contained. The broad rule alone finds about half of the named proteins, because many lack a PredGPI call or a
200-residue run (several EPA proteins have GPI no, and EPA2 has no signal peptide call).

## Candidates beyond the annotated names

The "other proteins called" are the interesting part and the weak part at once. Among them are real wall proteins (SED1, HKR1,
MSB2-like, SAG1, PIR/HSP150-type) and many "uncharacterized" proteins with a signal peptide, a GPI call and a long disordered
run, for example *C. glabrata* CAGL0J0515, CAGL0L0015 and CAGL0G1021 (the last two pass the narrow rule) and *K. marxianus*
FIM1_3170 (narrow). They are listed in `candidates.tsv` with the reason each is listed. Nothing here says which are true cell-wall
proteins; the upper bound on false calls is the "other" column (0.2 to 0.3 % of each proteome), and it includes real hits.

## What to keep in mind

* **All five species are Saccharomycotina yeasts** and the classifier was trained mostly on yeast positives (21 of 35 tier A). This is
  a test across yeast species, not across fungi. A non-yeast test (Pezizomycotina, Basidiomycota) is still missing and is where the
  fixed-threshold problem found in the clade hold-out (realised false-positive rate 12 to 36 % at a 3.5 % target) would show up.
* **A model score alone is not specific.** Without the signal-peptide condition, high-scoring novel proteins in *S. cerevisiae* are
  cytoplasmic disordered proteins (VRP1, PRM7, TDA7, LAS17, NUP1). The calls above require a signal peptide with the models.
* The thresholds come from out-of-fold scores of the 850 and are applied to scores of a model fit on all 850, so their scales differ
  somewhat; the call counts are a rough calibration and not a precision estimate.
* The name-based reference is partly homology-derived and not independent.
* Of PredGPI's limits: several named adhesins have no GPI call and are found only by the model route.

## Next decisions (owner)

1. Treat "at least 2 of 3 models and a signal peptide" as the working candidate call for further species, with the narrow rule as the
   stricter flag; or choose another rule.
2. Add a non-yeast proteome test (for example *Aspergillus fumigatus*, *Coccidioides*, *Cryptococcus*) to check thresholds across clades.
3. Hand-check the candidate list (starting with *C. glabrata* and *K. marxianus* narrow-rule proteins) against the literature and
   Pfam, or add them as a candidate positive set.
4. ESM-2 comparison on the same folds (separate comparator, per the earlier decision).
