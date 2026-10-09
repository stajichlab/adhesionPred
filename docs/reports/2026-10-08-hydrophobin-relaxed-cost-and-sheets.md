# Hydrophobin relaxed level: proteome cost (L7a) and evidence sheets (L8)

*2026-10-08. Branch `hydrophobin-validation` (local, not pushed). Module `hydrophobin_relaxed` (L6a, commit `1201bc4`) run on the 12 proteomes. Freeze commit `f1986e2`.
Tables: `docs/reports/data/sorting_hat/hydrophobin_ext/relaxed_cost.tsv`, `analysis/hydrophobin_truth/extra_calls.tsv`, `analysis/hydrophobin_truth/evidence_sheets.tsv`.*

## 1. Cost on the test proteomes (limit 5 unlabelled extra calls per 10,000 proteins, frozen)

Extra call = `hydrophobin_relaxed` hit that `pfam_hydrophobin` and `pfam_hsba` do not call and that has no T1/T2/T3/LP label. Unlabelled true hydrophobins count as cost (spec 6.3).

| Proteome | Part | Proteins | Relaxed hits | Also strict | Labelled extra | Unlabelled extra | Per 10,000 | Limit |
|---|---|---|---|---|---|---|---|---|
| *A. fumigatus* Af293 | test | 9647 | 7 | 6 | 0 | 1 | 1.04 | ok |
| *A. fumigatus* A1163 | test | 9942 | 6 | 5 | 0 | 1 | 1.01 | ok |
| *A. fumigatus* W72310 | test | 10556 | 9 | 7 | 0 | 2 | 1.89 | ok |
| *C. immitis* RS | test | 9910 | 3 | 2 | 0 | 1 | 1.01 | ok |
| *B. bassiana* ARSEF 2860 | test | 9478 | 17 | 7 | 0 | 10 | **10.55** | **over** |
| *F. graminearum* PH-1 | test | 11193 | 17 | 3 | 2 | 12 | **10.72** | **over** |
| *P. expansum* MD-8 | test | 10624 | 8 | 5 | 1 | 2 | 1.88 | ok |
| S288C | tuning | 6722 | 0 | 0 | 0 | 0 | 0.00 | ok |
| *C. albicans* | tuning | 6212 | 1 | 0 | 0 | 1 | 1.61 | ok |
| *B. dermatitidis* ER3 | tuning | 8107 | 6 | 2 | 0 | 4 | 4.93 | ok |
| *F. fulva* Race5 | tuning | 13560 | 11 | 6 | 3 | 2 | 1.47 | ok |
| *P. ostreatus* PC9 | tuning | 10483 | 13 | 9 | 0 | 4 | 3.82 | ok |

No relaxed hit overlaps `pfam_hsba` (E11 column is 0 everywhere). Tuning proteomes are in sample for the cutoff.

**Cost rule result: failed in 2 of 7 test proteomes** (*B. bassiana* and *F. graminearum*, about 10.6 and 10.7 per 10,000, twice the limit). By the ship rule frozen in `freeze.json`, the relaxed level is **not added to `categories.yaml`** while any condition fails. `hydrophobin_relaxed` stays as a module.

## 2. What the over-limit calls are (triage only; no label is written)

The 22 extra calls in the two over-limit proteomes include proteins that look like non-hydrophobins by a Swiss-Prot hit (for example *F. graminearum* F0349401_003565, 62.7% identical to a cutinase; F0349401_005603, 83% identical to an endoglucanase; F0349401_004011, a PR5-like receptor kinase; *B. bassiana* F1BB8A46_006071, an endochitinase fragment match) and small cysteine-rich proteins with no hit at all (for example F1BB8A46_000663, 74 aa, 8 Cys; F0349401_010222, 70 aa, 8 Cys). Whether the second group are hydrophobins is the owner's call (L8).

## 3. Evidence sheets (L8)

`evidence_sheets.tsv` has 61 rows, one per protein with a Pfam or relaxed call and no label: 40 relaxed-only, 18 strict and relaxed, 3 strict-only. Columns: length, cysteine count and gaps, R0, TMHMM, strict families, relaxed score and model, frozen-spacing class, BLAST top hits against the 174 known hydrophobins and Swiss-Prot 2023_03 (E <= 1e-3), and empty `owner_decision` and `owner_reason` columns.
Triage counts (a heuristic, not a label): among the 40 relaxed-only rows, 7 have a strong Swiss-Prot hit (at least 40% identity over at least 70% of the query) to a non-hydrophobin, 4 are longer than 250 aa, 3 have a weak hit and 26 have no Swiss-Prot hit. None of the 61 has a BLAST hit to the 174 known hydrophobins at E <= 1e-3, including the 21 with a Pfam call, which shows how divergent the family is.

## 4. What this means for the ship rule

| Condition (spec 6.2 to 6.5) | Status |
|---|---|
| Recall (at least 3 of 6 Pfam-missed clusters) | passed (6 of 6, report L5) |
| Cost in every test proteome | **failed** (2 of 7) |
| Hard negatives per group | passed (0 calls in every group) |
| Precision (at least 5 resolved extra-call clusters, Wilson lower bound at least 0.5) | pending the owner's decisions on the sheets |

If the owner finds that most extra calls in *B. bassiana* and *F. graminearum* are real hydrophobins, the cost definition (which counts unlabelled true hydrophobins as cost) is mis-specified for species with many unannotated hydrophobins. Changing it is a new pre-registration (v1.1), not a reinterpretation of this result.
