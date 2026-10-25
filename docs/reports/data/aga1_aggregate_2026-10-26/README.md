# First aggregate report: 10 proteomes (2026-10-26)

Output of `cellsurface_sorting_hat_aggregate` on the engine outputs of 10 fungal proteomes, grouped by NCBI class.
The input runs (engine outputs, FASTA files, NCBI taxdump) are outside the repository; no sequences are stored here.

```
cellsurface_sorting_hat_aggregate --manifest manifest.tsv --taxdump nodes.dmp --names names.dmp --rank class --out <new folder>
```

Manifest labels: Scer_S288C, Calb_SC5314, Cgla_CBS138, Kmar_FIM1, Caur_UniProt, Cpar_CDC317, Afum_Af293, Cimm_RS, Cneo_H99, Rtor_NP11
(each `run_dir` is the engine `--out` folder of that proteome). Code commit: the aggregator commit of the branch `aga1-aggregator`.
`categories.yaml` sha256 of all runs: `d12d160697632ff5f92d673205c691a54714d2e6b13eec816ee10b6a515fe50b`. NCBI `nodes.dmp` sha256 starts `171b54b2`.

Files: `report.md` (read this first), `proteomes.tsv`, `counts_by_proteome.tsv`, `counts_by_group.tsv`, `heatmap.png`, `groups.png`, `status_coverage.png`.

| file | sha256 |
|---|---|
| counts_by_group.tsv | 1282abc420cd99cb3a1522c4c6920b078bb747abdea27fe903bd46cebf0b811a |
| counts_by_proteome.tsv | 64dfd498412cdb70ac94daf362ccaaec6e82bf065686db5f42d96b9f2bf0c733 |
| proteomes.tsv | 97cc6391791ce7dd581f512e978fc4dbca920930d39a2e88082fe1473734cf93 |
| report.md | 95d57d8ff2c7a924cf1c526214a18728af09fc1b0bdeca08aeaf36a511124b3e |

## Known limits of this first report

* Seven modules have no table in these runs (`pfam_*`, `antigen_lookup`, `allergen_homology`, `repeat02`). Calls that read only them are
  NA (not drawn as 0). The Pfam and TMHMM scans have finished on the cluster but are not yet converted to module tables.
* `tandem_repeat_protein` is an OR of two detectors, and `repeat02` has no table, so it is not assessable for most proteins. Its raw rate
  is 1,000 per 1,000 assessable by construction; the aggregator flags every cell with fewer than half of its proteins assessable
  (`low_assessable`) and leaves it out of the group summaries.
* Only `gpi_anchored_protein` and `gpi_long_disorder_protein` have a (smoke) status, and only for *S. cerevisiae* and *C. albicans*. All
  other cells are `unvalidated` and marked `*`. Frozen yeast thresholds are used for the model call in all clades; per-clade
  recalibration waits for the hand-check verdicts.
* Counts are counts of calls, not prevalences. Groups of one proteome carry no spread.
