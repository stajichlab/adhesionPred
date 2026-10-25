# AGA1-class prediction step: AIUPred and IUPred3 disorder results (2026-10-24)

Record of the third tool run of the AGA1-class prediction step (plan: `docs/reports/2026-10-24-aga1-prediction-step-plan.md`;
report: `docs/reports/2026-10-24-aga1-disorder-results.md`). SignalP and PredGPI are in
`docs/reports/data/aga1_prediction_step_2026-10-24/`; their columns are repeated here so each table is complete on its own.

| Item | Value |
|---|---|
| Code commit run on HPCC | `424f5e8861dfd4c0cb56fd7533c5e69610120086` (merge of PR 143; `disorder.sbatch` is unchanged since) |
| AIUPred job | SLURM 30453936, partition `exfab`, node gpu12, 1 GPU + 4 CPUs, 1 min 15 s, `aiupred/2.x`, defaults (`-g 0`), torch 2.6.0+cu124 |
| IUPred3 job | SLURM 30473264, partition `short`, node c10, 4 CPUs, 21 min, `iupred/3`, `long` mode |
| Input | 850 proteins, 445,655 residues; FASTA sha256 `708384fbdbe4a2b24714c8547a12ef23671faff7e1bd75fd5c9cf8475ac63a0e` (sequences, not committed) |
| HPCC work directories | `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/aga1_prediction/aiupred` and `.../iupred3` |

## Files

| File | Content |
|---|---|
| `aga1_signalp_predgpi_aiupred.tsv` | one row per protein: scoping columns, SignalP, PredGPI, and AIUPred disorder columns |
| `aga1_signalp_predgpi_iupred3.tsv` | the same with IUPred3 disorder columns |
| `raw_provenance/aiupred_provenance.txt`, `raw_provenance/iupred3_provenance.txt` | the jobs' provenance files, unedited |

The raw per-residue files (`disorder_raw.txt`, 5.7 MB each) are NOT committed: every row carries the residue letter, so they
reproduce the protein sequences. Their hashes are below; the files stay in the HPCC work directories and in the owner's local
copy (`adhesionPred-specificity-results/aga1_aiupred`, `.../aga1_iupred3`).

## How the tables were made

```
python3 -I analysis/aga1_class/convert_predictions.py --manifest analysis/aga1_class/manifest_accessions.txt \
  --signalp ../aga1_prediction_step_2026-10-24/raw/signalp/prediction_results.txt \
  --predgpi ../aga1_prediction_step_2026-10-24/raw/predgpi/predgpi_raw.gff3 --lengths <lengths.tsv of stage_inputs.py> \
  --disorder <disorder_raw.txt> --disorder-tool <aiupred|iupred3> --out converted.tsv
python3 -I analysis/aga1_class/build_record_table.py --predictions converted.tsv --out <table>
```

Disorder columns: `disorder_tool`, `disorder_threshold` (0.5), `disorder_n_res`, `disorder_mean`, `disorder_frac` (fraction of
residues at or above 0.5), `disorder_longest_run` (longest run of residues at or above 0.5). Real layouts found on 2026-10-24:
AIUPred writes a banner of `#` lines, then per record a `#>ID` line, `pos<TAB>aa<TAB>score` rows and a blank line; the converter
now reads `#>ID` as a record line. IUPred3 rows are written by the job under its own `>ID` line. Both tools overshoot [0, 1]
slightly (AIUPred 12 rows up to 1.0036; IUPred3 142 rows from -0.0123 to 1.0096); the converter accepted up to 0.02 outside at that time (widened to 0.05 for the proteome run) and uses
the values as written. The IUPred3 provenance file shows a `torch_version` traceback only because IUPred3 does not use torch.

## sha256 (a test checks every line that names a committed file)

```
847842ea44dea074683a43ba6b335be36b045ede8085e5040e5167fe8c1e0c62  aga1_signalp_predgpi_aiupred.tsv
afb4f2bc948429af4b9415b9d3eb2d917a1f461f1f71479b89038f1a9e253921  aga1_signalp_predgpi_iupred3.tsv
3148e6590fa6d4af768c193d4911f6ee9b1a43b8331c1e3f52618a8043b6554b  raw_provenance/aiupred_provenance.txt
59c605f8e5e42c001933f684a1fb490bd9c4247880b8d7e1152c78908059580b  raw_provenance/iupred3_provenance.txt
```

Raw files not committed (hash of the file in the HPCC work directory, copied back by the owner):

```
5810381d06e5fd89e76e9f04b4b18a7e7283428712dbb15ed90ebda9e97b60b3  aiupred/disorder_raw.txt
93d415a8a94c732ac6bbfe63d7cae306073ddfdf7956092ef49c09d915bbfb24  iupred3/disorder_raw.txt
```
