# Baseline recall of the Pfam models on the Kiss 2019 control set (Task 5.4, run 2026-10-22)

Run record of SLURM job 30233470 (`analysis/wall_fusion_truth/submit_baseline.sbatch`, state COMPLETED, 13 s, exit 0, 35 MB). Result note: `docs/reports/2026-10-22-kiss-baseline-recall.md`. The job writes no call, no status and no tier.

| item | value |
|---|---|
| code pin | `9472e32e319a213211a69051f87208c8ff2a3f5f` (merge commit of PR #104; HPCC clone `adhesionPred-kiss53`, clean tree) |
| HMMER | 3.4 (`hmmer/3.4`) |
| Pfam | release 38.2, `/srv/projects/db/pfam/current/Pfam-A.hmm` resolved to `/bigdata/operations/pkgadmin/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`; sha256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`, the same as in the earlier jobs |
| settings | `hmmsearch --cut_ga` for every model; `-E 1e-5 --domE 1e-5` would apply only to a model without a GA line (there is none) |
| models | 39 fetched, 39 with a GA line, 0 without (21 wall and fusion models of `analysis/wall_fusion_report/wall_fusion_models.tsv`, 18 active family rows of `data/sorting_hat/family_table.tsv`) |
| inputs (sha256) | `kiss_sequences.faa` `66af609eeea24212f7813918fd5957489ded4e4e45792deadd955a0a93bacdf6` (293 sequences: 132 controls, 161 shared negatives); `family_rep_sequences.faa` `1fa26c78a81f7bb0bf6bbd94891cadccb0d05d953e5f77e71b0722a6b5f79a8f` |
| other hashes | see `provenance.tsv` (domain table, members table and cluster table of Task 5.3) |

## Files

| file | what it is |
|---|---|
| `provenance.txt`, `provenance.tsv` | written by the job |
| `models_searched.tsv`, `no_ga_models.tsv` | the 39 models; the models without a GA line (empty) |
| `recall_by_slice.tsv` | per slice and seed class (`slice` and `all39` labels): rows and clusters hit, Wilson interval, status. No pooled row. |
| `recall_by_clade.tsv` | the same per clade (fine and coarse) and stratum (`any`, `clade_only`, `shared`) |
| `clusters_missed.tsv` | the clusters of positive controls that no model of their slice hits (44 rows); the candidate list of PR 6 |
| `members_hits.tsv` | per control sequence: cluster, clade, seed classes, which models hit it |
| `wall_rows_false_hits.tsv` | the 91 wall rows against the 18 family models, descriptive only |
| `negatives_false_hits.tsv` | hits of each of the 39 models on the 161 shared negatives, by cluster composition |
| `family_rep_descriptive.tsv` | the family representatives, descriptive, never counted |

## Left on HPCC, not committed

`domtbl_kiss.txt` (48 KB) and `domtbl_famrep.txt` (5.9 KB) in `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_cluster_split/baseline_out/`, the raw `hmmsearch --cut_ga` domain tables (their sha256 is in `provenance.tsv`). The committed tables hold accessions, ids and counts; they hold no sequence.
