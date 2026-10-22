# Cluster split of the Kiss 2019 control set (Task 5.3, run 2026-10-22)

Run record of SLURM job 30229366 (`analysis/wall_fusion_truth/cluster_split.sh`, state COMPLETED, 15 s, exit 0, 0.9 GB). Result note: `docs/reports/2026-10-22-kiss-cluster-split.md`.

| item | value |
|---|---|
| code pin | `df6367d7e00c6664d6a10058dce7e99a014403ea` (merge commit of PR #103) |
| MMseqs2 | `MMseqs2/17-b804f` (build `b804fbe384e6f6c9fe96322ec0e92d48bccd0a42`) |
| settings | `easy-cluster --min-seq-id 0.3 -c 0.5`; `easy-search --min-seq-id 0.3 --cov-mode 1 -c 0.5 -s 7.5 -e 1e-3` |
| inputs (sha256) | `kiss_sequences.faa` `66af609eeea24212f7813918fd5957489ded4e4e45792deadd955a0a93bacdf6`; `members.tsv` `74371cc4008c37539f704dc0bc44bd21924e22f6e2edcdab23292963bf1a5225`; `seeds.faa.gz` `e0254d76bcb7ddc6e95479904b56bc3a857945c610a14f4c2651761daf456888` |
| inputs, size | 293 sequences (132 positive controls, 161 shared negatives, six of them `<accession>#neg` copies), 3851 Pfam seed sequences of 39 models |

## Files

| file | what it is |
|---|---|
| `provenance.txt`, `mmseqs_version.txt` | written by the job |
| `clusters.tsv` | `accession<TAB>cluster_id` (MMseqs2 clusters) |
| `members_split.tsv` | per sequence: cluster, class, slice, overlap flags, `seed_class_slice`, `seed_class_all39`, best identity and nearest model, `leakage_suggested` |
| `clusters_split.tsv` | per cluster: members, classes, slices, `mixed_pos_neg`, both seed classes, clades and strata (fine and coarse) |
| `summary.tsv` | clusters per class and slice; the 35 family representatives are listed as not clustered |
| `clade_floors.tsv` | made locally from `clusters_split.tsv` by `analysis/wall_fusion_truth/clade_floors.py`: clusters per slice, clade and stratum, with the floor status |

The tables hold accessions, ids and flags only. They hold no sequence, and none of the 20 JGI-only rows appears with a sequence (the local keys are `jgi:<portal>:<id>`). `seed_hits.m8` (485 KB, the MMseqs2 hits to the seed sequences) and `clu_cluster.tsv` stay on HPCC under `adhesionPred/_workdir/kiss_cluster_split/out/` and are not committed. The HPCC clone is `adhesionPred-kiss53`.
