# AGA1-class proteome test on five yeast proteomes, 2026-10-25

Report: `docs/reports/2026-10-25-aga1-proteome-test.md`. Everything is keyed by UniProt accession; no sequence, per-residue score or repeat-unit text is in this folder.

## What was run

| Item | Value |
|---|---|
| Proteomes | UniProt reference proteomes: S. cerevisiae S288C UP000002311 (6,067), C. albicans SC5314 UP000000559 (6,041), C. glabrata CBS138 UP000002428 (5,200), K. marxianus UP000422736 (5,070; strain FIM1, assembly GCA_001854445.2; the DMKU3-1042 proteome UP000065495 has no UniProtKB entries), C. auris UP000230249 (5,433); 27,811 proteins; downloaded from rest.uniprot.org on 2026-10-24/25 |
| Staging | `analysis/aga1_class/stage_proteomes.py`, 32 chunks of at most 1,000 (`chunks.tsv`) |
| HPCC code commit | `87ee2da4a08b3218013d6631beb0ef7997118afc` (PR 147), `analysis/aga1_class/proteome_array.sbatch`, array tasks 1-32 each |
| Jobs | SignalP 6 eukarya fast: 30483377; PredGPI `-m gff3`: 30483378; AIUPred 2.x on 1 GPU (`exfab`): 30483379; all 96 tasks COMPLETED, 32 chunks per tool, every provenance file names the commit above |
| Merge | `analysis/aga1_class/merge_proteome_chunks.py` (`merge_provenance.txt`) |
| Repeat detector | `analysis/cocci_repeats/14_repeat_detect_general.py` run locally on the chunks (27,224 proteins of 80 aa or more; the `rep_unit` column is not kept) |
| Features and scores | `feature_classifier.py features`, then `apply` with logreg, rf, gboost fit on all tier A + negatives of `docs/reports/data/aga1_feature_classifier_2026-10-25/features.tsv` |
| Summary | `analysis/aga1_class/proteome_summary.py` (false-positive target 3.5%) |

```
python3 -I analysis/aga1_class/proteome_summary.py --features features_proteomes.tsv \
  --scores logreg=scores_logreg.tsv rf=scores_rf.tsv gboost=scores_gboost.tsv \
  --oof docs/reports/data/aga1_feature_classifier_2026-10-25/oof_scores.tsv \
  --train-features docs/reports/data/aga1_feature_classifier_2026-10-25/features.tsv \
  --proteins proteins.tsv --names <the five UniProt FASTA files> --fpr 0.035 --out-dir <new dir>
```

## Files

| File | Content |
|---|---|
| `summary_by_species.tsv` | per species: proteins, signal peptide, GPI, broad and narrow rule, calls per model, known (850) proteins and novel ones |
| `name_based_recovery.tsv` | proteins with an adhesin-like UniProt name and how many each flag finds; the same for all other proteins |
| `thresholds.tsv` | score of each model at which the 850 negatives (out of fold) have a 3.5% false-positive rate |
| `candidates.tsv` | proteins not among the 850 that pass the broad or narrow rule or 2 of 3 models, plus the top 25 per species by mean rank |
| `known_proteins.tsv` | the 850-protein members present in each proteome and their flags |
| `all_proteins_flags.tsv` | every protein: calls, rule flags, scores |
| `chunks.tsv`, `merge_provenance.txt` | the staged chunks (sha256) and the merge record |

## Inputs not in the repository (sha256)

```
7f1eb913d176824eb55eeab8a80dda4200c6019ebbe817f6f3ff1facd4988c68  features_proteomes.tsv
e700f4e07ffa4740efb77809b5fb2e1b1c6f7290ee3800b96b408916e11b3816  predictions_proteomes.tsv
1a3395187f359fa4b994d4ba91f8029e246ca7bd5a95b4c0cf8e3d727f8cd7a6  stalk_proteomes.tsv
3901d5f9ba49f6af6f410bfa08801bd03b1e918eac25013b9da6eadc12ec19a3  rep14_proteomes.tsv
39ae6d937c0c3694566a0d87c7b0c19bad3d4dfe8cc9393cf3f220204f878638  scores_logreg.tsv
5ee861840e73807c324b2a59307dfbbdf925be910dc1ef2870a346e7f6523bdd  scores_rf.tsv
501f13fb456085049c4309e30e90c6e5a9a41599c563f7eb1c43a0048f2e59c4  scores_gboost.tsv
aiupred/disorder_raw.txt (merged raw AIUPred output with residue letters)  75fab7fc4000b05d3ade81d08e98bfcd17ef54ef5faa88580429d177371a5b29
signalp/prediction_results.txt  3c5909d57334d2aeab0094b211c4ebc4c483f325b5b06595bc79e88a5fafc458
predgpi/predgpi_raw.gff3  a9f4479c77199155012e754706d8ffb468fa3e6bcbdf33d5106369152a0e92b4
```

UniProt FASTA files (sha256 of the downloaded files):

```
3b922074c272716908e40e64e87da7b373a5252c6b23e9c8609b590d244a1ea3  Calb_SC5314.UP000000559.faa
472ca2f708a4035ccc46795030762a1610890fe81dbdb9167856ff21e07b3ca7  Caur_UniProt.UP000230249.faa
e8170003d6101d0abd454f788822a8ccef8ab7a01614a8f294e155f5cb335fc0  Cgla_CBS138.UP000002428.faa
214e2f3c0199a9bb4bb9381226cee2eced9c4a60daf464eb9b582deab9352bc9  Kmar_UniProt.UP000422736.faa
1aaa81860c8bda68e81e4290703d871c7c8bd12a928038772fbc91b7e17b342a  Scer_S288C.UP000002311.faa
```

## sha256 of the committed files (a test checks every line)

```
0c5a7f074014ec757aad4982cef34a794da30f152d47ec616178138de2b29ff9  all_proteins_flags.tsv
f19f29afd32b04475e748c39935c762a1defbd5d2e064abb1435f2935df7772d  candidates.tsv
ec10ef7d42026e7833db731e3519117666270bdea8f3997f145e97b530d8e0c8  chunks.tsv
948a3101d96270eb1981bd1dd9b5fc181d055f218fba646e90d8e0ae56e49c81  known_proteins.tsv
b0f30e0296569548ab45962da0816370f9941622f68ae907ec8fcc0f9edff4c9  merge_provenance.txt
b052c3683d9c83abbe9509e4e8f260e160e528f85547a802df2563a40098ec0c  name_based_recovery.tsv
7a54e561c5b9346f8d1e5500203da63a3387a9026538051411f954c2ba374a86  summary_by_species.tsv
5bda0e7e8a432dbee140e0004e5db01f6078ce528491365215ef42582bf50d5f  thresholds.tsv
```
