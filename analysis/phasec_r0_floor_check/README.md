# R0 floor check on the curated negative strata (S1 gate, D16, Q4)

`r0_stratum_floor_check.py` answers one question of the 2026-10-17 review plan (stop point S1): do the
negatives of the curated strata alone (N-int and N-sec; GO location terms, non-IEA) meet the status floors
for the R0 call in each Phase C test species? PM-TM (predicted topology) is `assumed` and is left out
(decision D16).

Floors (`docs/agent-tasks/COMMON-RULES.md` item 7; `src/cellsurface_sorting_hat/status.py`): at least 20
negatives, at least 20 independent clusters of negatives, and a 95% interval half-width of at most 0.10 for
the specificity.

The script is read-only. It uses the standard library only (it imports without HPCC). It reads three files
of the Phase C directory (`metrics.json`, `clusters.tsv.gz`, `eval_table.tsv.gz`) and writes one TSV to the
path of `--out`. How the numbers are made is in the module docstring: the same rules as
`src/cellsurface_sorting_hat/calibration/phasec.py` (`_curated_headline`), written out.

## Result of the HPCC run (2026-10-17)

Curated strata only. `halfw` is the half-width of the 95% interval of the specificity (Wilson on the protein
count, widened with Wilson on the cluster count). "All" is the pooled value that the old R0 status files carry.

| set | negatives (all) | clusters (all) | spec (all) | halfw (all) | curated negatives | curated clusters | spec (curated) | halfw (curated) | floors met |
|---|---|---|---|---|---|---|---|---|---|
| S1:Scer_SGD | 3785 | 3156 | 0.9657 | 0.0072 | 3783 | 3154 | 0.9659 | 0.0064 | yes |
| S1:Calb_CGD | 459 | 410 | 0.9434 | 0.0278 | 438 | 399 | 0.9521 | 0.0213 | yes |
| S3-Eurotiomycetes:Afum_ASPFU | 45 | 38 | 1.0000 | 0.0459 | 45 | 38 | 1.0000 | 0.0459 | yes |
| S3-Eurotiomycetes:Anid_EMENI | 164 | 151 | 0.9878 | 0.0227 | 163 | 150 | 0.9939 | 0.0174 | yes |
| S3-Basidiomycota:Cneo_H99_GOA | 32 | 31 | 0.9375 | 0.1022 | 32 | 31 | 0.9375 | 0.0937 | yes |
| S3-Basidiomycota:Umay_MYCMD | 28 | 24 | 0.8929 | 0.1442 | 27 | 23 | 0.9259 | 0.1163 | no (half-width) |

A. nidulans (the only species with R0 `estimated`) meets all floors on the curated strata. Umay_MYCMD fails the
half-width floor (0.1163) on the curated strata and also failed it on the pooled negatives (0.1442); its status is
`smoke` either way. No species has fewer than 20 curated negatives or fewer than 20 curated clusters.

## How to run it (HPCC, read-only, small job)

Do not run it on the login node. Use a compute node, and write the output to `$SCRATCH`, not to the repository:

```bash
srun -p short -c 1 --mem=4G --time=00:10:00 \
  python3 analysis/phasec_r0_floor_check/r0_stratum_floor_check.py \
    --workdir _workdir/step1_compare/phasec \
    --out "${SCRATCH:?}/r0_stratum_floor_check.tsv"
```

Options: `--variant` (default `V-go`), `--candidate` (default `R0`), `--truth` (default `direct`).
The script prints the table and the line `FLOORS_MET for every set: True|False`.
