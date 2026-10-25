# AGA1-class proteome test on Candida parapsilosis CDC317, 2026-10-26

Report: `docs/reports/2026-10-26-aga1-nonyeast-and-cpar-test.md`. Everything is keyed by UniProt accession; no sequence, per-residue score or repeat-unit text is in this folder. Same pipeline, frozen models and thresholds as `docs/reports/data/aga1_proteome_test_2026-10-25/`; nothing was retuned.

## What was run

| Item | Value |
|---|---|
| Proteomes | UniProt reference proteome *C. parapsilosis* CDC317 UP000005221 (5,777 proteins in 6 chunks); downloaded from rest.uniprot.org on 2026-10-25/26 |
| Staging | `analysis/aga1_class/stage_proteomes.py` (`chunks.tsv`) |
| HPCC code commit | `87ee2da4a08b3218013d6631beb0ef7997118afc`, `analysis/aga1_class/proteome_array.sbatch` |
| Jobs | SignalP 6 eukarya fast: 30560275; PredGPI `-m gff3`: 30560276; AIUPred 2.x on 1 GPU (`exfab`): 30560277; all 18 tasks COMPLETED; every provenance file names the commit above |
| Merge | `analysis/aga1_class/merge_proteome_chunks.py` (`merge_provenance.txt`) |
| Repeat detector | `analysis/cocci_repeats/14_repeat_detect_general.py` run locally (5,671 proteins of 80 aa or more; the `rep_unit` column is not kept) |
| Features and scores | `feature_classifier.py features`, then `apply` with logreg, rf, gboost fit on all tier A + negatives of `docs/reports/data/aga1_feature_classifier_2026-10-25/features.tsv` |
| Summary | `analysis/aga1_class/proteome_summary.py`, false-positive target 3.5% (thresholds in `thresholds.tsv`, identical to the five-yeast run) |

## Files

`summary_by_species.tsv`, `name_based_recovery.tsv`, `thresholds.tsv`, `candidates.tsv`, `known_proteins.tsv`, `all_proteins_flags.tsv`, `chunks.tsv`, `merge_provenance.txt` (contents as in the five-yeast record).

## Inputs not in the repository (sha256)

```
63329781f124f9431640a9e9f34d944e9e9bb32a52bcb574935258372ad55749  features_cpar.tsv
ecb48dce92f729ce2f6e0bb9dba09851ed47183fbed965464289ee2da27d435e  predictions_cpar.tsv
651fd72b458a45a847646d38edbe65d0b8f6049d849837b11649cd5940c6e528  stalk_cpar.tsv
e07178390d31ff1efc8c28f3137b434154d0ec358dce6c6848429c1afd1c0b27  rep14_cpar.tsv
351d7c4558d4d87c6e949a1a2dfab6c13e0f23e223fbd4260a553287f46adbc5  scores_logreg.tsv
30d25cf8274b485884a62e670e3c431220524f74515875d2c3a97e573b4c090c  scores_rf.tsv
2394c0703232972797f5bd837363d0e8a61fcf721a548711063d7e1d1a7fea16  scores_gboost.tsv
merged raw outputs (aiupred carries residue letters):
a0a58a79d346e03310b2d5652556036581d8f93ba1e496c69236672c81a16272  signalp/prediction_results.txt
412993e0da82416ed20ea99712f39d0e9fe7f7fdc2706b1949b55ff4699f217c  predgpi/predgpi_raw.gff3
83343b958d3016f231f4c0db530789312000e3c695dc04bdc7e0af6f65e3f195  aiupred/disorder_raw.txt
```

UniProt FASTA files (sha256 of the downloaded files):

```
dd103b0807145dba15a78c05ac16b335b6715e438a10359e51906c9dba9c7c80  Cpar_CDC317.UP000005221.faa
```

## sha256 of the committed files (a test checks every line)

```
0678ce7c63454484c28badc7e5793cafe933fdbca7c601a4447d3d6eee4a1f31  all_proteins_flags.tsv
d62b18e88f8c249e64c3cfba9ec5e3b567b2d4ba36ffb0666670043230312989  candidates.tsv
51d0ee88a0c589ca96d5ed1f0866723f525f4bcd1f5686f6bf4d30b265420ce3  chunks.tsv
243a1b45ae52f04cb3bdfc99f5b4a47bbf4d494cb2a62d1da24d1c205ecad6dc  known_proteins.tsv
6f4065752a36261265550fd236550a27a65c9d3b5b53c537d94309560f06b389  merge_provenance.txt
3c919340876fe44e1873b725c6576725238b26b9edc9cfad57b63ac4f6fa2a4e  name_based_recovery.tsv
32fa342fc22cc487ac290b1b89eb6b6ab4db5614e7a7d1840bb8d5158a3e74b4  summary_by_species.tsv
5bda0e7e8a432dbee140e0004e5db01f6078ce528491365215ef42582bf50d5f  thresholds.tsv
```
