# PR 6 step 6.0: the whole-Pfam scan of the clusters the panel misses (DRAFT, not submitted)

Plan: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md`, decision 1 and section 6. Nothing here has run on HPCC. Decision 1 (owner, 2026-10-22): the whole-Pfam scan comes first; only a cluster with no Pfam hit becomes a custom-HMM target; a cluster that hits is dropped from the custom list and flagged (the proposal for it is to add that Pfam model to the panel, a separate later decision per model). Fusion clusters are scanned too, descriptively.

## Files

| File | Purpose |
|---|---|
| `pr6_pfam_scan.sbatch` | The HPCC job: extracts the 66 member sequences of the 44 missed clusters from the staged `kiss_sequences.faa`, runs `hmmscan --cut_ga` against all of Pfam-A 38.2 and a separate weak-hit `hmmscan -E 1e-3`, then the summary script. Partition `short`, `${SCRATCH:?}`, no `/tmp`, pin check, refuses an existing `$OUT`. |
| `pr6_pfam_scan_summary.py` | Standard library only. Turns the two domain tables into the per-cluster tables. Writes no call, status or tier. |
| `tests/wall_fusion_truth/test_pr6_pfam_scan.py` | Fixture tests (hand-counted classes, the loud report, refusals, determinism) and static checks of the sbatch. |

`hmmscan` and not `hmmsearch`: 66 query sequences against about 20,000 models is the shape hmmscan is made for (it reads the pressed binary database; the `.h3f/.h3i/.h3m/.h3p` files exist next to `Pfam-A.hmm`, confirmed read-only on HPCC on 2026-10-22). `--cut_ga` works for hmmscan: a hit passes the GA threshold of the model it hits, so no model is skipped.

## Classes (per missed cluster; a cluster is hit when any member is)

| Class | Meaning | What happens to it |
|---|---|---|
| `no_pfam_hit` | no GA hit to any Pfam family | stays a custom-HMM target (`custom_targets.tsv`); its weak hits (E <= 1e-3, failing GA) are listed, descriptive, and change nothing |
| `hit_other_pfam` | a GA hit to a family outside the 39-model panel | dropped from the custom list. Proposal column: `add_to_panel_proposal` if the best family covers at least 0.5 of one member's sequence (union of the alignment intervals / sequence length; **a provisional proposal for the owner, not a decision**), else `partial_domain_for_owner` (the hit may be one generic domain of a larger protein) |
| `hit_panel_model` | a GA hit to one of the 39 panel models | reported LOUDLY (stderr banner, `panel_hits.tsv`, `summary.txt`). Not an error: a cluster is missed by the models of its own slice and a model of another slice can still hit it |

**Expected from the baseline (check, not assume).** `docs/reports/data/kiss_baseline_2026-10-22/members_hits.tsv` shows two members of missed clusters hit by a panel model of another slice: Q59SF7 (pir1, cluster KS0082, PF22799 PIR1-like_C) and Q4WLB9 (cfmA, cluster KS0116, PF05730 CFEM). The hmmscan run must therefore report at least these two as `hit_panel_model`; if it does not, the two tools disagree and that is a finding to report.

## Outputs (a new `$OUT`)

`scan_per_cluster.tsv` (the 44 clusters: class, proposal, panel models, other families with cover, score and E-value of the best member, weak families), `scan_counts.tsv` (per slice and per function group of the spec: `no_pfam_hit`, `hit_other_pfam`, `hit_panel_model`), `custom_targets.tsv`, `panel_proposals.tsv` (per Pfam family: clusters, slices, cover range), `panel_hits.tsv`, `summary.txt`, `provenance.tsv`, and the raw `domtbl_ga.txt`, `domtbl_weak.txt` and `provenance.txt` (job, commit, HMMER version, Pfam path, resolved path, sha256 and whether it equals the one of the earlier jobs, the sequence extract's sha256). The extracted sequences are not an output. No file written by the summary holds a sequence (a test scans for it).

## Run order (exact commands; none has been run)

Each step needs the owner's go-ahead (D21). **No new staging is expected**: the job reads `kiss_sequences.faa` from the directory that was staged on 2026-10-22 and everything else from the clone. This is true only once the clone holds the merge of the PRs with these scripts and these record tables (the baseline record, the cluster record and the spec data), because the job's pin check and file checks read them from the clone.

1. Update the clone to the merge commit that holds this PR:
   ```
   B=/bigdata/stajichlab/jstajich/projects
   git -C $B/adhesionPred-kiss53 pull --ff-only origin main
   git -C $B/adhesionPred-kiss53 rev-parse HEAD     # the merge commit
   git -C $B/adhesionPred-kiss53 status --short     # must print nothing
   ```
2. Submit (one job, about a minute of search plus the hash of the 2.2 GB `Pfam-A.hmm`):
   ```
   W=$B/adhesionPred/_workdir/kiss_cluster_split
   sbatch --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<merge commit>,IN=$W/in,OUT=$W/pfam_scan_out,PFAM_RELEASE=38.2 \
     -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log \
     $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_pfam_scan.sbatch
   ```
3. Read the log: `pin ok`, `66 member sequences`, `inputs ok`, then the `summary.txt` lines. The job stops (and says why) if a member is missing from the FASTA, the Pfam release name does not match, a pressed index file is missing, or a domain table has no `# [ok]` trailer. A different `Pfam-A.hmm` sha256 is recorded and reported, not fatal.

## Checks after the job

- `pfam_sha256_same_as_earlier_jobs: yes` in `provenance.txt` (Pfam 38.2, `4b0da639…`).
- 44 clusters in `scan_per_cluster.tsv` (adhesion 12, wall 28, fusion 4) and 66 members.
- KS0082 and KS0116 are `hit_panel_model` (see above).
- The weak table is descriptive: it must not move any cluster out of `no_pfam_hit`.

## How the result feeds 6.1

Only the `no_pfam_hit` clusters (`custom_targets.tsv`) go to step 6.1 (one model per homologous family, found by a sensitive sequence search; decision 2). The `panel_proposals.tsv` families are a list for the owner, one decision per model; none is added to the panel by this job. If no cluster is left without a Pfam hit, the PR ends at the record PR (6c) with the list of Pfam models proposed for the panel. The next owner go-aheads are: the clone update, this job, then (after the owner reads the tables) the record PR and the 6.1 scripts.

## Result of the run

The run of job 30235946 (2026-10-22, code `ebb29a9`) is recorded in `docs/reports/data/pr6_pfam_scan_2026-10-22/` (tables, provenance and a README; the raw domain tables stay on HPCC) and described in `docs/reports/2026-10-22-pr6-pfam-scan.md`. It found 2 clusters hit by a panel model (KS0082 and KS0116, the two of the baseline), 38 hit by another Pfam family and 4 with no Pfam hit. The use of the 38 is an open owner decision.
