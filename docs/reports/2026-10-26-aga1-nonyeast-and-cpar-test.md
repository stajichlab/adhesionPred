# AGA1-class flags and classifiers on four non-yeast proteomes and *C. parapsilosis*

Date: 2026-10-26. Records and hashes: `docs/reports/data/aga1_nonyeast_test_2026-10-26/` and `docs/reports/data/aga1_cpar_test_2026-10-26/`.
Same pipeline as the five-yeast test (`2026-10-25-aga1-proteome-test.md`).

## Declared before the run (nothing was retuned)

The three classifiers are the ones fit on the 850 labelled proteins; their score thresholds are the out-of-fold 3.5% false-positive scores
already recorded for the yeast run (`thresholds.tsv` is identical in every record); the rules are the broad rule (signal peptide, GPI call,
AIUPred run of 200 or more) and the narrow rule (broad plus the detector-14 repeat call); the working candidate call is "at least 2 of 3
models and a signal peptide". The question was whether these carry across clades, which the clade hold-out of the classifier study
predicted they might not do evenly.

## What was run

UniProt reference proteomes of *A. fumigatus* Af293, *C. immitis* RS, *C. neoformans* H99, *R. toruloides* NP11 (34,946 proteins, 37 chunks;
SLURM arrays 30547662, 30547663, 30547664, all 111 tasks completed) and *C. parapsilosis* CDC317 (5,777 proteins, 6 chunks; arrays 30560275,
30560276, 30560277, all 18 tasks completed), on the pinned commit of PR 147.

## Counts

| Species | Proteins | Signal peptide | GPI | Broad rule | Narrow rule | logreg | rf | gboost | 2 of 3 models | Of the 850 present |
|---|---|---|---|---|---|---|---|---|---|---|
| *A. fumigatus* | 9,647 | 753 | 150 | 7 | 3 | 19 | 15 | 85 | 20 | 52 |
| *C. immitis* | 9,712 | 446 | 91 | 8 | 2 | 28 | 23 | 87 | 27 | 2 |
| *C. neoformans* | 7,427 | 331 | 89 | 2 | 0 | 2 | 5 | 70 | 3 | 9 |
| *R. toruloides* | 8,160 | 415 | 120 | 2 | 1 | 6 | 19 | 117 | 11 | 0 |
| *C. parapsilosis* | 5,777 | 337 | 127 | 28 | 9 | 26 | 29 | 35 | 25 | 1 |

(For comparison, the five yeast proteomes had 10 to 25 broad-rule proteins and 17 to 48 proteins called by 2 of 3 models.)

## What carried across

* **The rules found the known adhesin-like proteins they should.** In *C. immitis* the spherule outer wall glycoprotein (SOWgp, CIMG_04613) passes the broad
  rule, the narrow rule and all three models (mean rank 0.9998). Related SOWgp proteins from other *Coccidioides* entries are tier A
  training positives, so this is a check on a close relative, not an independent protein. In *A. fumigatus* the one tier A protein of the 850
  (AFUA_3G08990) is recovered, and "cell wall mannoprotein 1" (AFUA_2G05150), AFUA_4G09600 ("GPI anchored protein, putative", run 929, narrow
  rule), an agglutinin-like middle-domain protein (AFUA_4G00500) and several extracellular Ser/Thr- or proline-rich proteins are called.
* ***C. parapsilosis* behaves like the yeasts.** 37 proteins carry an adhesin-like UniProt name (Als, Hyr, Flo11-domain, hyphal wall); the broad
  rule finds 16 and the broad rule or 2 of 3 models with a signal peptide finds 18, including ALS1, ALS3, ALS6, ALS7, ALS11, HYR31 and
  several Hyr-like and Flo11-domain proteins (the narrow rule flags 9, among them ALS6, ALS7 and RBT1). The name reference is rough and partly homology-derived.
* ***Cryptococcus* and *Rhodotorula* have few calls** (broad 2 and 2, narrow 0 and 1). Flagged there: *C. neoformans* CNAG_05915 and CNAG_00759
  (uncharacterized) and *R. toruloides* proteophosphoglycan-like proteins (ppg4 passes the narrow rule). Whether these are wall proteins is not tested here.

## What did not carry evenly

* **Thresholds are not equally calibrated.** At the same 3.5% false-positive threshold the gradient-boosting model called 70 to 117 proteins in the four
  non-yeast proteomes (0.9 to 1.4% of each), against 2 to 28 for logistic regression and 5 to 23 for the random forest, and 27 to 57 in the yeasts.
  Its calls outside yeast are mostly not supported by the other two models (2 of 3 models gives 3 to 27). This is the clade effect the hold-out
  predicted. The 2-of-3 rule absorbs it here, but a single-model call should not be used on a new clade without recalibration.
* **Models call secreted Ser/Thr/Pro-rich proteins without a GPI call**, e.g. *A. fumigatus* "extracellular proline-rich protein" and
  "extracellular threonine-rich protein": the models score them high, the broad rule excludes them (no GPI call). Whether they are wall-bound is not
  shown by these tools.
* **Enzymes with a Ser/Thr-rich linker are still flagged by the broad rule** (*A. fumigatus* endochitinase A1 and an endo-1,3(4)-beta-glucanase,
  *C. immitis* endochitinase CTS2, *C. parapsilosis* GH18 chitinases CHT21 and CHT22), as in the labelled set. A domain check would be needed to remove them.
* The name-based reference is not usable in the four non-yeast proteomes (2 to 8 matching names each).

## CFEM proteins across all ten proteomes

UniProt-annotated CFEM (PF05730) proteins: 15 in the five yeasts, none flagged; in the non-yeast proteomes 20 (*A. fumigatus* 4, *C. immitis* 8,
*C. neoformans* 3, *R. toruloides* 5) of which 4 are flagged: *A. fumigatus* cfmA (305 aa, a labelled hard negative in the 850) and an
extracellular Ser/Thr-rich protein with a CFEM domain; *C. immitis* CIMG_00795 and an extracellular Ser/Thr-rich protein (CIMG_03472). In
*C. parapsilosis* 1 of 8 is flagged, CFEM5 (918 aa, a disordered run of 230); the other long one, CFEM6/CSA1 (924 aa), and the short hemophore-type
CFEM proteins are not. So the flags miss the short CFEM proteins as intended, and catch some long CFEM-domain proteins that also carry a long
disordered region; the Pfam CFEM call remains the route for the family.

## Limits

The flags and models still rest on the same 850 labelled proteins, mostly yeast; the SOWgp and one *A. fumigatus* positive are the only
non-yeast tier A members. The test shows the rules and the 2-of-3 call keep a low call rate and find known adhesin-like proteins in these
proteomes; it does not give a false-positive rate, because "other proteins called" includes proteins that may be real. Name matching is rough.

## Next decisions (owner)

1. Keep "2 of 3 models and a signal peptide" as the working candidate call and treat single-model calls as clade-specific.
2. Hand-check the candidates (`candidates.tsv` in each record), starting with the narrow-rule proteins.
3. Decide whether to add a domain check (catalytic domain) to remove enzyme linkers from the broad flag.
4. ESM-2 comparison on the same folds (separate comparator).
