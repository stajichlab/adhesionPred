# PR 6 step 6.0: whole-Pfam scan of the 44 clusters the panel misses (run 2026-10-22)

Run record of SLURM job 30235946 (`analysis/wall_fusion_truth/pr6_pfam_scan.sbatch`, state COMPLETED, 3 min 35 s, exit 0, 0.4 GB). Result note: `docs/reports/2026-10-22-pr6-pfam-scan.md`. The job and the summary write no call, no status and no tier.

| item | value |
|---|---|
| code pin | `ebb29a9a4ceaee9caf0453f74ed411b2508cce82` (merge commit of PR #107; HPCC clone `adhesionPred-kiss53`, clean tree) |
| HMMER | 3.4 (`hmmer/3.4`) |
| Pfam | release 38.2, `/srv/projects/db/pfam/current/Pfam-A.hmm` resolved to `/bigdata/operations/pkgadmin/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`; sha256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`, the same as in the earlier jobs (`provenance.txt`) |
| settings | `hmmscan --cut_ga` (hits); `hmmscan -E 1e-3 --domE 1e-3` without `--cut_ga` (weak hits, descriptive only) |
| queries | 66 member sequences of the 44 missed clusters, extracted from the staged `kiss_sequences.faa` (sha256 `66af609eeea24212f7813918fd5957489ded4e4e45792deadd955a0a93bacdf6`); sha256 of the extract `bd634aa2b46910f141a0a9141b711b9e88a731a54a5155aeb013ea1c16284f95` |
| other hashes | the domain tables, the missed-cluster table and the group table are in `provenance.tsv` |
| provisional setting | `proposal_fraction` 0.5 (cover of the sequence). It is a proposal for the owner. Nothing tuned it. |

## Files

| file | what it is |
|---|---|
| `provenance.txt`, `provenance.tsv` | written by the job and the summary |
| `summary.txt` | the counts and the loud report line |
| `scan_counts.tsv` | per slice and per function group of the PR 6 spec: `hit_panel_model`, `hit_other_pfam`, `no_pfam_hit` |
| `scan_per_cluster.tsv` | the 44 clusters: class, proposal, panel models, other families (cover, score and E-value of the best member), weak families |
| `custom_targets.tsv` | the clusters with no Pfam hit at the GA threshold (4); the input of step 6.1 |
| `panel_proposals.tsv` | per Pfam family: clusters, slices and cover range, for the clusters with a hit covering at least 0.5 of the sequence. Proposals for the owner. |
| `panel_hits.tsv` | the clusters that a panel model hits (2) |

## Left on HPCC, not committed

`domtbl_ga.txt` (34,934 B) and `domtbl_weak.txt` (44,606 B) in `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_cluster_split/pfam_scan_out/`, the raw `hmmscan` domain tables (their sha256 is in `provenance.tsv`). No file here holds a sequence (a test scans for it).
