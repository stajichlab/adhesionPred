# AGA1-class proteome test on four non-yeast fungal proteomes, 2026-10-26

Report: `docs/reports/2026-10-26-aga1-nonyeast-and-cpar-test.md`. Everything is keyed by UniProt accession; no sequence, per-residue score or repeat-unit text is in this folder. Same pipeline, frozen models and thresholds as `docs/reports/data/aga1_proteome_test_2026-10-25/`; nothing was retuned.

## What was run

| Item | Value |
|---|---|
| Proteomes | UniProt reference proteomes: *A. fumigatus* Af293 UP000002530 (9,647), *C. immitis* RS UP000001261 (9,712), *C. neoformans* H99 UP000010091 (7,427), *R. toruloides* NP11 UP000016926 (8,160); 34,946 proteins in 37 chunks; downloaded from rest.uniprot.org on 2026-10-25/26 |
| Staging | `analysis/aga1_class/stage_proteomes.py` (`chunks.tsv`) |
| HPCC code commit | `87ee2da4a08b3218013d6631beb0ef7997118afc`, `analysis/aga1_class/proteome_array.sbatch` |
| Jobs | SignalP 6 eukarya fast: 30547662; PredGPI `-m gff3`: 30547663; AIUPred 2.x on 1 GPU (`exfab`): 30547664; all 111 tasks COMPLETED; every provenance file names the commit above |
| Merge | `analysis/aga1_class/merge_proteome_chunks.py` (`merge_provenance.txt`) |
| Repeat detector | `analysis/cocci_repeats/14_repeat_detect_general.py` run locally (34,293 proteins of 80 aa or more; the `rep_unit` column is not kept) |
| Features and scores | `feature_classifier.py features`, then `apply` with logreg, rf, gboost fit on all tier A + negatives of `docs/reports/data/aga1_feature_classifier_2026-10-25/features.tsv` |
| Summary | `analysis/aga1_class/proteome_summary.py`, false-positive target 3.5% (thresholds in `thresholds.tsv`, identical to the five-yeast run) |

## Files

`summary_by_species.tsv`, `name_based_recovery.tsv`, `thresholds.tsv`, `candidates.tsv`, `known_proteins.tsv`, `all_proteins_flags.tsv`, `chunks.tsv`, `merge_provenance.txt` (contents as in the five-yeast record).

## Inputs not in the repository (sha256)

```
bc3417c0c08e36e293822a3f2bcf26046b8e570aa478f894a0bfc70fd4bb71c6  features_nonyeast.tsv
312d51155461e59a79fb4c29ab40307884a2a822fc2c0db8eae58f86db6c569d  predictions_nonyeast.tsv
cbbfbf430c2438e0cb750285fc894ace3afac5e0cca69d76d51a248b966b6519  stalk_nonyeast.tsv
35c3de00ceef154587a0faacfd044a872f56facdf09200ae9ab36df57b341497  rep14_nonyeast.tsv
0cdb2abb77357521cc309ad5b9e23cfda6bdd70d85bba0d3f2693f92c8f9f99b  scores_logreg.tsv
f243a9f8a6d565487ffa3932f0c0739546c0c78e4117c96129e2e42ba2b0d285  scores_rf.tsv
3fcf70c99732a900228cfba85b183fbbfba552292da582353504ad6cb1a39879  scores_gboost.tsv
merged raw outputs (aiupred carries residue letters):
bc7837266e76601e0455d61ec9f452ebdaa891c70a4c4aebe8e0265a81bb5209  signalp/prediction_results.txt
be849a35a65420aba4c76fbc3e4131136331c63b6f83744701a412f2cececdf9  predgpi/predgpi_raw.gff3
f2567b67df1dad5bb29b166a8678c88cebc5c3ae0861c865c350a81afb7c6fad  aiupred/disorder_raw.txt
```

UniProt FASTA files (sha256 of the downloaded files):

```
155a063bbe25057a1ab30107ea80f76a583e57538b2efb76b9f07059b1d255fd  Afum_Af293.UP000002530.faa
bdf83d4cfdaeda9b74e403a12f928155f7c57689746e8cafc9018fd886e669fc  Cimm_RS.UP000001261.faa
12d99a40f6a3bc00f48fca9bbaa16ab21373ed2c52c2c5b3cc5c4aaf1f2e2049  Cneo_H99.UP000010091.faa
2304a94b3bb9d682d651b43b37e15f80a3e62d6c2d322aaac9bf2f8c8b23cecc  Rtor_NP11.UP000016926.faa
```

## sha256 of the committed files (a test checks every line)

```
142ea507f60ad0c6248f4739c3ffc0c0b679a6ae5d704869295262026ffbdbd5  all_proteins_flags.tsv
41059d74599f0af85431fd2c914ed357857de9de7068d25acc4b0b2e2cb6fddb  candidates.tsv
dad685aca76ab60e0ec69283591c4b903f1dcb8a3b39aa374a97177326afe563  chunks.tsv
9b65bb4dfc08e3b5c77f4dba3cddd9146abb095d7286ce30244eb908488285d3  known_proteins.tsv
8265e5886e81aacdbff6142ef7bbfa9788d698076fae40d2548c9e4e2b0399c6  merge_provenance.txt
a2427996fe16c6a830ed433b779a4e3f55a6a02d58b46f0aa54124ec1718c56f  name_based_recovery.tsv
d755fa375110090a82e8ae430b44ddfdd1a34642925c750dbe73e15a72f25e4e  summary_by_species.tsv
5bda0e7e8a432dbee140e0004e5db01f6078ce528491365215ef42582bf50d5f  thresholds.tsv
```
