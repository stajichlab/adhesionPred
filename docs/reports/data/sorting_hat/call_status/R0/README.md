# R0 status file (refreshed 2026-10-17)

`step1_rule@R0.json` is the status file of the module `step1_rule@R0` (the signal peptide call `signal_peptide_protein`). It was written again on HPCC with `calibrate phasec` (Task 3.10, script `scripts/sorting_hat/refresh_3_10/01_phasec_r0.sbatch`).

- There is one copy here. The file is byte-identical in the four work directories that had an R0 status file (Afum_Af293_UniProt, Calb_SC5314, Scer_S288C, Scer_S288C_nodubious). Size 12,889 bytes. sha256 `285818ddb8d9466d761e56ecb5509f61d157521e01408e186518f83c8513bf7d`.
- It is the only committed record of the one `estimated` claim for `signal_peptide_protein`.

## Statuses in the file

| taxon | calibration set | status | n_neg | clusters (neg) |
|---|---|---|---|---|
| 4932 | S1:Scer_SGD | smoke | 3783 | 3154 |
| 5476 | S1:Calb_CGD | smoke | 438 | 399 |
| 746128 | S3-Eurotiomycetes:Afum_ASPFU | smoke | 45 | 38 |
| 162425 (A. nidulans) | S3-Eurotiomycetes:Anid_EMENI | estimated | 163 | 150 |
| 5207 | S3-Basidiomycota:Cneo_H99_GOA | smoke | 32 | 31 |
| 5270 | S3-Basidiomycota:Umay_MYCMD | smoke | 27 | 23 |

A. nidulans has 109 positives in 100 clusters.

## The leakage exception

Every entry carries a `leakage_exception`: date 2026-10-06; reason "leakage unmeasured (SignalP training overlap), uncapped by owner decision"; signer "owner (docs/paper/03-status-and-validation-rules.md section 6, item 1)". Without it the loader caps each entry at `smoke`, because the overlap between the positives and the SignalP 6 training data was not measured. The owner chose to leave `estimated` uncapped for this call.

## Negatives

`negative_basis` is `curated` for the strata N-int and N-sec and `assumed` for PM-TM (where present). Only `curated` strata count toward `n_neg`, so the numbers above are curated negatives.
