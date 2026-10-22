# Kiss 2019 baseline recall: Task 5.4 (DRAFT)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`, Task 5.4. Decisions: section 10 of `docs/reports/2026-10-18-kiss-feasibility.md` and `docs/reports/2026-10-22-kiss-cluster-split.md`. **Nothing here has run on HPCC. No job is submitted.** The result of the run is not known; this file says what will be produced and how to check it.*

## What this adds

| File | Where it runs | What it does |
|---|---|---|
| `submit_baseline.sbatch` | HPCC, SLURM `short`, 4 CPUs, 8 GB, 1 h | Pin check; fetches the 39 models from Pfam-A.hmm with `hmmfetch -f` (by name, the method of `scripts/sorting_hat/pfam_hmmsearch.sbatch`); `hmmsearch --cut_ga` against the 293 Kiss sequences and against the 30 sequences of the family representatives; runs `baseline_recall.py`; writes everything to a new `$OUT`. |
| `baseline_recall.py` | HPCC job or local, stdlib | Reads the domain tables and the Task 5.3 tables; writes the tables below. Writes no call, no status, no tier. |

Tests: `tests/wall_fusion_truth/test_baseline_recall.py` (hand-counted fixture), `test_submit_baseline_script.py` (static checks and the GA-split step).

## Pfam and HMMER (confirmed read-only on HPCC, 2026-10-22)

`/srv/projects/db/pfam/current/Pfam-A.hmm`: `current` points to `2026-01-27-Pfam38.2`, which resolves to `/bigdata/operations/pkgadmin/srv/projects/db/pfam/2026-01-27-Pfam38.2/` (Pfam-A.hmm 2,246,909,846 B with the `.h3f`, `.h3i`, `.h3m`, `.h3p` index files). It is Pfam 38.2, the release the earlier jobs used (earlier provenance: sha256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`, HMMER 3.4); module `hmmer/3.4` exists (the default). The job, not the login node, hashes the 2.2 GB file, records the resolved path and the sha256 in `provenance.txt`, and checks only the release name (`*Pfam38.2/`). A different sha256 is recorded and printed as a `REPORT:` line; it does not stop the job. The 39 model names of the two repository tables equal the Pfam names of the InterPro entries (checked against `seeds_summary.tsv` of the input build, 39 of 39).

## Run order (all commands are NOT run; each HPCC step needs a go-ahead)

1. **Clone** at the merge commit that holds these scripts. The clone of the cluster job (`adhesionPred-kiss53`) is at an older commit; use a new clone or `git fetch` and check out the new commit there. `CODE_COMMIT` is that commit.
2. **Input.** `kiss_sequences.faa` and its `.sha256` are already staged at `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_cluster_split/in/`. Everything else the job reads is committed in the clone: the run record `docs/reports/data/kiss_cluster_split_2026-10-22/{members_split,clusters_split}.tsv`, the control set (both files) and `family_rep_sequences.faa`, the two model lists. **No new staging is needed.**
3. **Job** (after the go-ahead):
   ```
   B=/bigdata/stajichlab/jstajich/projects
   W=$B/adhesionPred/_workdir/kiss_cluster_split
   sbatch --export=ALL,PROJ_ROOT=<clone>,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/baseline,PFAM_RELEASE=38.2 \
     -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log <clone>/analysis/wall_fusion_truth/submit_baseline.sbatch
   ```
   `$W/baseline` must not exist. Expected run time: minutes (39 models against 323 sequences; most of the time is hashing Pfam-A.hmm).
4. **Checks after the job** (read-only): `sacct` COMPLETED; the log shows `pin ok`, `inputs ok: 293 sequences`, `models: 39 fetched, 39 with a GA line, 0 without` (a line that says otherwise is reported, see below) and `wrote $W/baseline`; `provenance.txt` shows the pinned commit, hmmer 3.4, Pfam 38.2 and the sha256 line; every `domtbl_*.txt` ends with `# [ok]` and shows `--cut_ga` in its option line (the tables are refused otherwise); `recall_by_slice.tsv` has no line without a slice; `members_hits.tsv` has 132 rows.
5. **Bring the tables back** (they hold accessions and flags, no sequence) and write the report `docs/reports/2026-10-XX-wall-fusion-baseline.md` with the numbers recomputed from them. Only then a results PR.

## Outputs (`$OUT`)

| file | contents |
|---|---|
| `recall_by_slice.tsv` | per slice, per `seed_class_kind` (`slice`, `all39`) and per seed class (`seed-like`, `seed-independent`, always apart): n rows, k rows hit, n clusters, k clusters with at least one member hit, k clusters with every member hit, row and cluster recall, Wilson 95% interval on the clusters, status. There is no line that pools slices and no line that pools the two seed classes. |
| `recall_by_clade.tsv` | the same per held-out clade (fine = primary, coarse = second test) and stratum (`any`, `clade_only`, `shared`). `estimated` only when the cell has at least 5 clusters of its own; fusion is always `descriptive`. |
| `clusters_missed.tsv` | the positive clusters with no hit: slice, members, genes, species, clades, strata, both seed classes, best identity to a seed, `jgi_accept` flags. **This table, not a guess, names the PR 6 candidates.** |
| `wall_rows_false_hits.tsv` | descriptive false-hit panel: how many wall rows each of the 18 family models hits. Wall rows are never negatives and never truth for the adhesion call. |
| `negatives_false_hits.tsv` | false hits of each of the 39 models on the 161 shared negatives, per cluster composition (`negative_only`, `mixed`). The six `<accession>#neg` copies are among them. |
| `family_rep_descriptive.tsv` | the 35 family representatives (30 sequences): hits on their own slice's models, flagged `family_rep`, `counts_toward_floor=no`, `may_be_test_control=no`. They are in no recall table. |
| `no_ga_models.tsv` | models without a GA line (see below); empty when all 39 have one. |
| `members_hits.tsv` | per control: hit on its own slice's models, on the other slice's models, on the 2 fusion models, and the list of models hit. |
| `domtbl_*.txt`, `models_searched.tsv`, `provenance.txt`, `provenance.tsv` | raw tables and provenance. |

## Rules the numbers follow

- **A hit** is a row in the `--cut_ga` domain table of the model. No other threshold decides a recall number.
- **Adhesion** rows: recall of the 18 active family rows (family coverage). It is not the sensitivity of the sorting-hat call and changes no tier. **Wall** rows: recall of the 21 wall and fusion models. **Fusion**: descriptive only (5 rows, at most 2 clusters per clade).
- **A model with no GA line** is reported in `no_ga_models.tsv` and in the log (`REPORT:`), not skipped silently. `hmmsearch --cut_ga` stops on such a model, so the job splits it off and searches it with an E-value cut-off (1e-5) into its own table; its hits are counted in `no_ga_models.tsv` only and enter no recall number. Pfam-A has a GA line for every family, so none is expected.
- **Units.** The independent unit is the MMseqs2 cluster of Task 5.3. A cluster is recalled when at least one control member is hit; the table also gives clusters whose members are all hit. Every number carries its n.
- **`seed-like` clusters** partly measure the seed family the model was built from; they are reported apart from `seed-independent` clusters, under both seed class labels.
- **No performance claim is made by this draft.** No status file, tier or `narrow_ok` changes (D4 is an owner decision).

## Open questions

1. Whether the job's `hmmfetch -f` finds all 39 names in the pressed database. The names equal the InterPro IDs of the entries, but the pressed database itself was not opened here; the job stops with `hmmfetch found N of 39` if one is missing.
2. A sequence name in `kiss_sequences.faa` with a `#` (the six `#neg` copies) and the local keys `jgi:<portal>:<id>`: hmmsearch reads the name up to the first space; both were checked against `members_split.tsv` by the job, not by running hmmsearch on them.
3. The Wilson interval treats clusters as independent. It is a proxy; a bootstrap over clusters (as `calibration`) is the upgrade if the interval will be quoted.

## Run record and result (2026-10-22)

Job 30233470 ran this script at `9472e32` (COMPLETED, 13 s). The small output tables and the provenance are committed in `docs/reports/data/kiss_baseline_2026-10-22/` (README there); the two raw domain tables stay on HPCC. The result note is `docs/reports/2026-10-22-kiss-baseline-recall.md`. It changes no tier, status or `narrow_ok`.
