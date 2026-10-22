# Kiss 2019 control set and cluster split: Tasks 5.2 and 5.3 (DRAFT)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`, section PR 5. Decisions: section 10 of `docs/reports/2026-10-18-kiss-feasibility.md` (owner, 2026-10-22). **Nothing here has run on HPCC and nothing is submitted.** `cluster_split.sh` is a draft that needs the owner's go-ahead for each step (D21).*

## What this adds

| File | Task | Where it runs | What it does |
|---|---|---|---|
| `build_control_set.py` | 5.2 | local, offline, stdlib | Writes `data/controls/wall-fusion-kiss2019/` from `results/row_mapping.tsv` and the UniProt answers in `results/query_cache.json`. Never reads a JGI file. |
| `build_kiss_fasta.py` | 5.2 / 5.3 | **local only** | Joins the control sequences (UniProt from the repository, JGI from the owner's downloaded files) into one FASTA outside the repository. Checks every sha256. |
| `fetch_seeds_models.py` | 5.3 | local, needs network | Fetches the Pfam seed alignments of 39 models (21 wall and fusion models plus 18 active family rows) from InterPro into `seeds.faa.gz`. |
| `cluster_split.sh` | 5.3 | HPCC, SLURM `short` | MMseqs2 clusters (30% identity, 50% coverage) and the seed search; then runs `seed_split.py`. |
| `seed_split.py` | 5.3 | HPCC job or local | Labels each cluster `seed-like` (any member at 30% identity or more to a seed) or `seed-independent`. Writes `clusters.tsv`, `members_split.tsv`, `clusters_split.tsv`, `summary.tsv`. |

Tests: `tests/wall_fusion_truth/test_build_control_set.py`, `test_build_kiss_fasta.py`, `test_seed_split.py`, `test_cluster_split_script.py`. The control set is also checked by the L1 validator (`tests/cellsurface_sorting_hat/test_controls_validator.py` finds every committed `controls.tsv`).

## The control set (Task 5.2), as committed

`data/controls/wall-fusion-kiss2019/` (README there). Positives only, keyed by UniProt accession.

| | adhesion | wall | fusion | total |
|---|---|---|---|---|
| `controls.tsv` (public UniProt sequence; L1 validator: PASS) | 29 | 76 | 5 | 110 |
| `controls_local_sequence.tsv` (key `jgi:<portal>:<id>`, no public sequence) | 7 | 15 | 0 | 22 |
| both | 36 | 91 | 5 | 132 |

The 132 rows are the 133 strictly mapped Kiss rows with one merge (cot-2 and gs-1 are the same JGI protein 8184_Neucr2). In `controls.tsv` there are 56 Kiss clusters (per slice); with the local rows 66. `family_reps.tsv` has 35 rows (30 sequences; 12 of those accessions are also controls). `rejected.tsv` has 13 rows: 11 JGI-ID rows that are not mapped (hyr1, bhp2 and bhp3 as `name_conflict`) and 2 family rows without a locus.

**Why two control files.** The L1 validator requires every accession to have a sequence in `sequences.faa`. The 22 rows with no UniProt accession have a sequence only in a JGI file, and JGI sequences are never committed. They are therefore in a second table with the same columns, kept out of `controls.tsv`. They count as controls; `build_kiss_fasta.py` reads them back. This is a trade-off the owner may change (open question 2).

**Clusters are a proxy until Task 5.3.** `clusters.tsv` uses `KISS-<slice>-<cluster_number>`. COMMON-RULES rule 4 asks for MMseqs2. After the job, `build_control_set.py --clusters-tsv ... --split-members ...` replaces them and fills `leakage` and `seed_class`.

## Run order (all commands are NOT run; each HPCC step needs a go-ahead)

1. **Local: the seed FASTA** (network): `python3 -I analysis/wall_fusion_truth/fetch_seeds_models.py --out-dir ~/kiss_5_3_inputs`. It writes `seeds.faa.gz`, `seeds_summary.tsv` (the Pfam version of every alignment: compare with Pfam 38.2) and `seeds.faa.gz.sha256`. A model with no seed alignment on InterPro stops the script; report it, do not skip it.
2. **Local: the sequence FASTA**, outside the repository (the script refuses a path inside it):
   `python3 -I analysis/wall_fusion_truth/build_kiss_fasta.py --jgi-dir to_import/kiss2019_jgi --out-dir ~/kiss_5_3_inputs --include-family-reps --negative-set <DIR> [--negative-set <DIR2>] --on-overlap keep-both`.
   It writes `kiss_sequences.faa`, `members.tsv` and the `.sha256` files. **Do not pass a negative set until the owner confirms the sources (open question 1).** Without `--negative-set` the job clusters positives only.
3. **Stage to HPCC** (owner's go-ahead, then): `rsync -a ~/kiss_5_3_inputs/ hpcc:/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_5_3/in/`. The inputs hold JGI sequences: they stay under `_workdir` and are not committed.
4. **HPCC: a clone at the reviewed commit** (owner's go-ahead): the previous clone was deleted on 2026-10-22. `git clone` the repository to `/bigdata/stajichlab/jstajich/projects/adhesionPred-kiss_5_3` and check out the merge commit that holds these scripts. `CODE_COMMIT` is that commit.
5. **HPCC: the job** (owner's go-ahead):
   ```
   W=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_5_3
   sbatch --export=ALL,PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred-kiss_5_3,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/out \
     -o $W/cluster_split.log -e $W/cluster_split.log \
     /bigdata/stajichlab/jstajich/projects/adhesionPred-kiss_5_3/analysis/wall_fusion_truth/cluster_split.sh
   ```
   `$W` must exist and `$W/out` must not. Estimated size: under 400 sequences and 39 models: minutes. The script asks for 4 CPUs, 8 GB and one hour.
6. **Checks after the job** (read-only): `squeue`/`sacct` state COMPLETED; `out/provenance.txt` shows the pinned commit, the MMseqs2 version and the input sha256; `out/summary.tsv` has no pooled line; every sequence of `members.tsv` is in `members_split.tsv`; the cluster count matches `clusters.tsv`.
7. **Local: bring back the small tables** (they hold accessions, no sequences) and rebuild the set:
   `python3 -I analysis/wall_fusion_truth/build_control_set.py --clusters-tsv <out>/clusters.tsv --split-members <out>/members_split.tsv`. The L1 validator then checks the MMseqs2 clusters. Only then open a PR with `data/controls/wall-fusion-kiss2019/` and `<out>/clusters_split.tsv`, `summary.tsv`.

## What the job does not do

It does not choose the hard negatives, change a status, tier or call, or run `hmmsearch` (that is Task 5.4). It uses `module load MMseqs2/17-b804f` as `analysis/lectin_validation/run_cluster_seed.sh` does; it does not need `hmmer`.

## Open questions for the owner

1. **Which hard negatives?** Spec section 9 item 12: `data/curated/adhesins/hard_negative_seeds.tsv` (30 rows plus a comment line and header) and the shared set `data/controls/hard-negatives-shared/controls.tsv` (161 rows). The confirmation is not given, so no script reads either by default. The seeds file holds UniProt queries (column `uniprot_query`), not accessions or sequences, so using it needs a fetch step that is not written. The shared set has sequences.
2. **Two control files or one?** To keep every accepted row in one `controls.tsv`, the L1 validator would need a rule that a row with a local key may have no sequence in `sequences.faa` (a check of `sequence_sha256` instead). I did not change `src/`.
3. **Kiss wall controls that are also hard negatives.** The shared negative set holds 8 sequences that are also Kiss controls or family representatives: P23776 (exg1), P53301 (crh1), P32623 (crh2/utr2), P22146 (gas1), P29029 (cts1), Q4WLB9 (cfmA), Q4WXE6 (crh3), P0C954 (gel2). In the shared set they are negatives for the adhesion call (wall enzymes); in the Kiss set they are positives of the wall slice. D20 allows both, because the calls differ. A FASTA holds one sequence per accession, so `build_kiss_fasta.py` stops on them. `--on-overlap keep-both` writes the negative copy as `<accession>#neg`; the two copies land in one cluster, which `seed_split.py` marks `mixed_pos_neg`. Options: keep both (default suggestion), remove them from the negative set, or leave them out of the clustering.
4. **Family representatives in the clustering?** With `--include-family-reps` the 18 representatives that are not already controls join the clustering as class `family_rep`. They may bridge a control cluster to a training seed. Without the flag they are not clustered. The decision says they are training seeds and baseline report only; I did not decide whether they should shape the clusters.
5. **Which models are "the seeds"?** The 39 models are the 21 wall and fusion models plus the 18 active rows of `data/sorting_hat/family_table.tsv`. `analysis/lectin_validation` used 20 lectin models for its own set. Add models by editing `fetch_seeds_models.model_list`.
6. **`names.dmp` check.** The taxon names in the control set come from NCBI answers cached in 5.1. They were not checked against `names.dmp` (`python3 -m cellsurface_sorting_hat.controls <dir> --names-dmp names.dmp`).

## Could not confirm (nothing was run on HPCC)

- That `MMseqs2/17-b804f` is still installed under that name (it is the name in COMMON-RULES and in the repository scripts).
- That a compute node has no need for the network (the job reads only staged files; `fetch_seeds_models.py` is local on purpose).
- The `pident` scale of this MMseqs2 build. `finalize_truth.py` treats a file whose values are all 1.0 or less as a fraction, and so does `seed_split.py`.
- That InterPro serves seed alignments for all 39 models.

## How Task 5.4 uses the outputs

`baseline_recall.py` (Task 5.4, not written) reads `clusters_split.tsv` and `members_split.tsv`. Per slice (adhesion, wall, fusion; never pooled) it reports recall of the Pfam models per cluster and per held-out clade, with these columns of `controls.tsv`: `clade_fine` (primary level) and `clade_coarse` (second test); `cluster_stratum_fine` and `cluster_stratum_coarse` (shared or clade-only); `jgi_accept` (filter later); `seed_class` (report `seed-like` and `seed-independent` apart). A stratum is `estimated` only with 5 or more clusters of its own, otherwise `descriptive`; fusion is always descriptive. Family representatives appear only in the recall report and the training seeds, flagged `family_rep`; they never count toward a floor. The wall rows are a descriptive false-hit panel for the narrow call, never negatives, never truth for the adhesion call. The table of clusters that Pfam misses names the PR 6 candidates.
