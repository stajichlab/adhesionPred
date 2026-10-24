# AGA1-class prediction step: SignalP 6 and PredGPI results (2026-10-24)

Record of the first two tool runs of the AGA1-class prediction step (plan: `docs/reports/2026-10-24-aga1-prediction-step-plan.md`;
report: `docs/reports/2026-10-24-aga1-signalp-predgpi-results.md`). The AIUPred disorder run was not finished when this record
was made and is NOT included; the table has no disorder columns.

| Item | Value |
|---|---|
| Code commit run on HPCC | `424f5e8861dfd4c0cb56fd7533c5e69610120086` (merge of PR 143) |
| SignalP job | SLURM 30451696, `signalp/6` module, `--organism eukarya --mode fast --format none` |
| PredGPI job | SLURM 30453935, `predgpi/202001` module, `-m gff3` (default omega site) |
| SignalP version | the resolved path is the 6.0i install, the tool prints `6.0h` (see `raw/signalp/provenance.txt`) |
| PredGPI | `predgpi.py` sha256 `6b2a92d64ef7f37c55e9594ab08669e207150a63d096360f6439b0d268af4853` (see `raw/predgpi/provenance.txt`) |
| Input | 850 proteins (`analysis/aga1_class/manifest_accessions.txt`); FASTA sha256 `708384fbdbe4a2b24714c8547a12ef23671faff7e1bd75fd5c9cf8475ac63a0e` (the FASTA holds sequences and is not committed) |
| HPCC work directories | `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/aga1_prediction/signalp` and `.../predgpi` |

## Files

| File | Content |
|---|---|
| `aga1_signalp_predgpi.tsv` | one row per protein: accession, gene, organism, clade, `role_group` and UniProt features (from `docs/reports/data/aga1_class_scoping/proteins.tsv`), SignalP columns, PredGPI columns |
| `raw/signalp/prediction_results.txt`, `raw/signalp/provenance.txt` | raw SignalP output and the job's provenance file, unedited |
| `raw/predgpi/predgpi_raw.gff3`, `raw/predgpi/provenance.txt` | raw PredGPI output and the job's provenance file, unedited |

## How the table was made

```
python3 -I analysis/aga1_class/convert_predictions.py --manifest analysis/aga1_class/manifest_accessions.txt \
  --signalp raw/signalp/prediction_results.txt --predgpi raw/predgpi/predgpi_raw.gff3 \
  --lengths <lengths.tsv of stage_inputs.py> --no-disorder --out converted.tsv
python3 -I analysis/aga1_class/build_record_table.py --predictions converted.tsv --out aga1_signalp_predgpi.tsv
```

`lengths.tsv` (accession, length) is made by `stage_inputs.py` from the out-of-repo sequence cache and is not committed. PredGPI
columns: `gpi_pred` yes / no / NA, `gpi_omega_site`, `gpi_confidence_class` (high = estimated false positive rate <= 0.0015,
probable <= 0.005, weak <= 0.01), `gpi_fpr_upper_bound`, `gpi_note`. The gff3 has no false positive rate; none is invented.
Proteins of 40 aa or fewer (Q3E837, Q8TGS4, Q8TGS5) are not scored by PredGPI and are `NA`.

## sha256 (a test checks every line)

```
32f5ac7ef7a0b3745ba958b0c7fa7125feae135e6ad7d81a0f5ae10b59861750  raw/predgpi/predgpi_raw.gff3
1ab8f44182b4de05aed3a119cbab7ee056faac3bb14dfd2a14f993c1e406df2b  raw/predgpi/provenance.txt
34b0fc0db5f7527eeaf5cda8dcceea692af3ed537f4b3058beaab9651bcceee8  raw/signalp/prediction_results.txt
d71be5f2a1ff51054892eebc6224d9178129d46a2e662f3ccdff8fd58b13faaa  raw/signalp/provenance.txt
1465f31b1ff44d2922969c7060ad708253830dce632b855ef09d2dee5bfc9d35  aga1_signalp_predgpi.tsv
```

The raw files were copied back from HPCC by the owner (to a folder outside the repo) and committed here unedited; the hashes
above are those of the committed copies.
