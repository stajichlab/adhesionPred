# AGA1-class prediction step: disorder results (AIUPred and IUPred3)

Date: 2026-10-24. Third and last tool run of the prediction step planned in `2026-10-24-aga1-prediction-step-plan.md`.
Record and hashes: `docs/reports/data/aga1_disorder_2026-10-24/`. SignalP and PredGPI: `2026-10-24-aga1-signalp-predgpi-results.md`.
Nothing here defines the broad flag; that is the owner's call (decision 4 of the AGA1 grilling).

## What ran

| Tool | Job | Where | Time |
|---|---|---|---|
| AIUPred 2.0 (`aiupred/2.x`, defaults, GPU) | 30453936 | `exfab`, node gpu12 | 1 min 15 s |
| IUPred3 (`iupred/3`, `long` mode, CPU) | 30473264 | `short`, node c10 | 21 min |

Both covered all 850 proteins (445,655 residues, row counts exact). The AIUPred job had waited behind an account CPU cap
(`AssocGrpCpuLimit`, `exfab` account limit 32 CPUs), not behind a GPU shortage; once it ran it took about a minute. IUPred3 on CPU
was the fast answer while it waited. The scoping and SignalP/PredGPI numbers are unchanged.

## Two things the real output changed

1. **AIUPred writes its record headers as `#>ID`.** The converter treated them as comments and would have stopped with
   "score row before any '>' record". It now reads `#>ID` as a record line (test added).
2. **Both tools slightly overshoot [0, 1]** (AIUPred 12 rows, maximum 1.0036; IUPred3 142 rows, from -0.0123 to 1.0096). The
   range check was strict and stopped the first conversion. It now allows 0.02 outside the range and uses the values as written (test added);
   anything further out still stops it.

The per-residue raw files are not committed because each row carries the residue letter (they would reproduce the sequences).

## Results

Disorder is the fraction of residues with a score of at least 0.5 (`frac`) and the longest run of such residues (`longest`).

| Group | n | median length (aa) | AIUPred median frac | AIUPred median longest run | IUPred3 median longest run |
|---|---|---|---|---|---|
| tier A positives | 35 | 973 | 0.61 | 302 | 191 |
| tier B | 32 | 970 | 0.56 | 403 | 278 |
| hard negatives | 172 | 444 | 0.19 | 59 | 26 |
| background assumed negatives | 389 | 318 | 0.09 | 16 | 9 |

Per-group means are in the two tables. The two tools agree well per protein (Pearson r over 850 proteins: mean score 0.83,
disordered fraction 0.89, longest run 0.92). IUPred3 gives lower values throughout; the order of the groups is the same.

**Length is a confound.** Tier A proteins are about twice as long as the hard negatives, and a long protein tends to hold a longer run.
Keeping only proteins of 500 aa or more, tier A still has a median longest run of 310 (27 proteins) against 57 for hard negatives (61)
and 48 for background (105). So disorder is not only a length proxy, but a threshold has not been tested for length.

Combined with the earlier calls (signal peptide yes AND predicted GPI yes), the number of proteins that also hold a disordered
run of at least 100 residues is:

| Group | n | SP + GPI | and AIUPred run >= 100 | and AIUPred run >= 200 | and IUPred3 run >= 100 |
|---|---|---|---|---|---|
| tier A | 35 | 27 | 26 | 23 | 21 |
| tier B | 32 | 27 | 23 | 18 | 19 |
| hard negatives | 172 | 52 | 23 | 8 | 14 |
| background assumed negatives | 389 | 92 | 31 | 7 | 12 |
| Kiss rows (other clusters' positives) | 70 | 3 | 1 | 1 | 1 |
| other-mechanism adhesins | 49 | 20 | 12 | 6 | 9 |

These are counts for a candidate rule, not an evaluated flag: the 100 and 200 cuts were chosen to read the table, nothing was fitted,
and tier A (35) is the group the rule would be trained on. For example, SP + GPI alone keeps 27/35 tier A against 52/172 hard
negatives; adding "AIUPred run >= 200" keeps 23/35 against 8/172 hard negatives and 7/389 background.

## Next decisions (owner)

1. Define the broad flag, now that all three counts exist. Candidate: SP + GPI + a disordered run, with the run cut and the tool
   (AIUPred, IUPred3, or both) to be chosen. The Ser/Thr-rich stalk rule (my 0.5-in-30-residue-window cut) is still undecided.
2. The still-open items from the scoping study: clade unit for narrow classes, glycosylation predictor, ESM-2 tiling for the 74
   proteins over 1,022 aa, extra UniProt positives.
