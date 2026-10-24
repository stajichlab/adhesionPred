# AGA1-class prediction step: SignalP 6 and PredGPI results (2026-10-24)

Data and provenance: `docs/reports/data/aga1_prediction_step_2026-10-24/` (README with job ids, code commit `424f5e8`, hashes).
Plan: `docs/reports/2026-10-24-aga1-prediction-step-plan.md`. The AIUPred disorder job had not finished; it is not part of this report.

## What ran

Two tools, on the same 850 proteins (class positives tier A and B, Kiss rows, hard negatives, background, other-mechanism and
excluded proteins; see `docs/reports/data/aga1_class_scoping/proteins.tsv`).

* SignalP 6 (job 30451696; the path is the 6.0i install, the tool prints 6.0h), eukarya, fast mode. For each protein it says SP
  (a signal peptide for secretion is predicted, with a cleavage site) or OTHER. 850 rows: 619 SP, 231 OTHER.
* PredGPI 202001 (job 30453935), gff3 output. For each protein it gives either no GPI anchor (a `Chain` row) or a GPI anchor with the
  omega site (the residue where the anchor is attached). 850 rows: 250 GPI-anchor, 597 no, and 3 not scored.

How to read PredGPI's confidence: the gff3 score column is a confidence class of the tool, not a probability. 1.0 (called `high`
in the table) means an estimated false positive rate of 0.0015 or less, 0.70 (`probable`) means 0.005 or less, 0.55 (`weak`) means
0.01 or less. The tool calls a GPI anchor at 0.01 or less. Here 232 are high, 8 probable, 10 weak. The file has no false positive
rate and no raw score, so the table has the class and its upper bound only; nothing was invented. PredGPI never scores a sequence of
40 aa or fewer (it returns a silent "no"). Q3E837 (39 aa), Q8TGS4 and Q8TGS5 (35 aa) are `NA`, not `no`.

## SignalP signal peptide against UniProt annotation

Counts are recomputed from the committed table.

| Group (n) | SignalP predicts SP | UniProt signal annotation |
|---|---|---|
| class_positive tier A (35) | 35 | 35 |
| class_positive tier B (32) | 29 | 31 |
| hard_negative (172) | 160 | 163 (2 not annotated) |
| background_assumed_negative (389) | 303 | 389 |
| kiss_row_positive (70) | 20 | none (no UniProt annotation) |
| adhesin_other_mechanism (49) | 35 | 36 |
| excluded_surface_other_adhesion_phenotype (36) | 30 | 30 |
| excluded_indirect_regulator (60) | 1 | 0 |

Among the 778 proteins with a UniProt signal annotation (yes or no), SignalP agrees for 682 (596 yes/yes, 86 no/no) and differs
for 96: 95 have a UniProt signal but SignalP says OTHER (86 of them in the background group, 5 hard negatives, 2 tier B, 1
other-mechanism, 1 uncurated), and 1 the other way (an indirect regulator). So the agreement is complete for tier A and good
overall, but it is not perfect: SignalP misses about 22% of the background proteins that UniProt gives a signal (86 of 389). The
reason was not examined (UniProt may list signal peptides that SignalP fast mode does not call).

The 74 proteins over 1,022 aa all got a call: 53 SP and 21 OTHER. Of those, 54 have a UniProt signal annotation and 50 of the 54
were called SP (4 missed, about the same miss rate as the shorter proteins); the 21 OTHER include 12 Kiss rows and 5 indirect
regulators. Whether SignalP truncated the long sequences internally was not checked.

## PredGPI against UniProt GPI annotation

| Group (n) | PredGPI predicts GPI | UniProt GPI feature | predicted but not in UniProt |
|---|---|---|---|
| class_positive tier A (35) | 27 | 24 | 4 |
| class_positive tier B (32) | 28 | 11 | 18 |
| hard_negative (172) | 53 | 35 | 18 |
| background_assumed_negative (389) | 99 | 95 | 13 |
| adhesin_other_mechanism (49) | 20 | 9 | 11 |
| excluded_surface_other_adhesion_phenotype (36) | 16 | 7 | 9 |
| excluded_indirect_regulator (60) | 0 | 0 | 0 |
| kiss_row_positive (70) | 3 | none (no UniProt annotation) | not applicable |

Reference proteins: AGA1 (P32323) GPI, omega 700; SED1 (Q01589) omega 307; DAN4 (P47179) omega 1141, all class 1.0 (high);
MSB2 (P32334) no GPI. Over the 775 proteins that have a UniProt GPI feature value and a PredGPI call, 173 agree on yes, 517 agree
on no, 74 are PredGPI-only and 11 are UniProt-only.

## What this says

* The signal peptide agrees with UniProt wherever it is annotated for the class positives (35 of 35 in tier A) and does not
  separate the class from other secreted proteins: 303 of 389 background proteins and 160 of 172 hard negatives also have one.
* PredGPI adds calls UniProt lacks, mostly in tier B (18 of 28), but it is not decisive alone: 27 of 35 tier A (77%) against
  53 of 172 hard negatives (31%) and 99 of 389 background proteins (25%).
* The broad flag is NOT defined here. The owner decides after seeing these counts, the disorder run, and the Ser/Thr stalk rule.

## What was not established

* PredGPI calls are not validated beyond agreement with UniProt; UniProt's own GPI feature is incomplete, so "predicted only" is
  not the same as "wrong".
* The 1% cut is a threshold per sequence, not a false discovery rate over a proteome: 850 proteins at 0.01 can hold many false
  calls, and the background group (25% predicted) shows the order of magnitude of how many proteins are called.
* Omega site accuracy is unknown.
* No check of PredGPI on Candida or other clades was made; the 850 proteins are the scoping set only (clade counts are in
  `docs/reports/data/aga1_class_scoping/`).
* Whether SignalP truncated proteins over 1,022 aa is not known; only that all 74 got a call.
* The SignalP version string (6.0h) and the install path (6.0i) disagree; the tool output is kept as printed.

## Housekeeping

The plan note was said to claim that `sha256sum -c` works on the `.sha256` file. No sentence like that was found in the plan note
or in `analysis/aga1_class/README.md` (searched for `sha256sum`); the `.sha256` file holds the hash only. The README now says so
and gives the hand comparison (the sbatch pin checks already compare the hash this way).
