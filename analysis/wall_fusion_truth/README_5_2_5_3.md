# Kiss 2019 control set and cluster split: Tasks 5.2 and 5.3 (DRAFT)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`, section PR 5. Decisions: section 10 of `docs/reports/2026-10-18-kiss-feasibility.md` and the owner's answers of 2026-10-22 (below). **Nothing here has run on HPCC and nothing is submitted.** `cluster_split.sh` is a draft that needs the owner's go-ahead for each step (D21).*

## Owner decisions that shape these scripts (2026-10-22)

1. **Hard negatives = the shared set of 161 rows** (`data/controls/hard-negatives-shared/controls.tsv`), and only that set. It is the default negative input. `data/curated/adhesins/hard_negative_seeds.tsv` is excluded by this decision (it may be added later); there is no fetch step for it.
2. **A protein that is both a Kiss control and a shared negative keeps both roles.** Eight accessions overlap the shared set: P23776, P53301, P32623, P22146, P29029, Q4WLB9 (Kiss controls) and Q4WXE6, P0C954 (family representatives only). The six controls get a negative copy `<accession>#neg` in the clustering input (`--on-overlap keep-both` is the default); the two copies fall in one cluster, which `clusters_split.tsv` flags `mixed_pos_neg=yes` (and `has_overlap_copy=yes`). Filter on `mixed_pos_neg` to leave them out. The two representatives are not clustered, so they get no copy.
3. **Family representatives are not clustered.** The default of `build_kiss_fasta.py` leaves them out; `--include-family-reps` stays as an option and is not used. They keep their Kiss `cluster_number` (`family_reps.tsv`) and are a flagged group in the report (`summary.tsv` lists them as `family_rep (not clustered)`).
4. **Two seed classes.** `seed_split.py` writes `seed_class_slice` (wall and fusion rows against the 21 wall and fusion models; adhesion rows against the 18 active family rows) and `seed_class_all39` (against all 39). A cluster is `seed-like` if any member has 30% identity or more (target coverage 0.5 or more, E <= 1e-3) to a seed. A cluster of only negatives has `seed_class_slice=not_applicable`; a mixed cluster is judged by its control members.
5. **The set is two files.** `controls.tsv` and `controls_local_sequence.tsv`. Anything that reads the set must read both.

## What this adds

| File | Task | Where it runs | What it does |
|---|---|---|---|
| `build_control_set.py` | 5.2 | local, offline, stdlib | Writes `data/controls/wall-fusion-kiss2019/` from `results/row_mapping.tsv` and the UniProt answers in `results/query_cache.json`. Never reads a JGI file. |
| `build_kiss_fasta.py` | 5.2 / 5.3 | **local only** | Joins the control sequences (UniProt from the repository, JGI from the owner's downloaded files) and the 161 shared negatives into one FASTA outside the repository. Checks every sha256. |
| `fetch_seeds_models.py` | 5.3 | local, needs network | Fetches the Pfam seed alignments of the 39 models (21 wall and fusion models plus 18 active family rows) from InterPro into `seeds.faa.gz`. |
| `cluster_split.sh` | 5.3 | HPCC, SLURM `short` | MMseqs2 clusters (30% identity, 50% coverage) and the seed search; then runs `seed_split.py`. |
| `seed_split.py` | 5.3 | HPCC job or local | Writes `members_split.tsv`, `clusters_split.tsv`, `clusters.tsv`, `summary.tsv` with the two seed classes and the mixed-cluster flag. |

Tests: `tests/wall_fusion_truth/test_build_control_set.py`, `test_build_kiss_fasta.py`, `test_seed_split.py`, `test_cluster_split_script.py`. The control set is also checked by the L1 validator (`tests/cellsurface_sorting_hat/test_controls_validator.py` finds every committed `controls.tsv`).

## The control set (Task 5.2), as committed

`data/controls/wall-fusion-kiss2019/` (README there). Positives only, keyed by UniProt accession.

| | adhesion | wall | fusion | total |
|---|---|---|---|---|
| `controls.tsv` (public UniProt sequence; L1 validator: PASS) | 31 | 76 | 5 | 112 |
| `controls_local_sequence.tsv` (key `jgi:<portal>:<id>`, no public sequence) | 5 | 15 | 0 | 20 |
| both | 36 | 91 | 5 | 132 |

The 132 rows are the 133 strictly mapped Kiss rows with one merge (cot-2 and gs-1 are the same JGI protein 8184_Neucr2). In `controls.tsv` there are 57 Kiss clusters (per slice); with the local rows 66. `family_reps.tsv` has 35 rows (30 sequences; 12 of those accessions are also controls). `rejected.tsv` has 13 rows: 11 JGI-ID rows that are not mapped (hyr1, bhp2 and bhp3 as `name_conflict`) and 2 family rows without a locus.

**Accessions found by sequence.** Of the 22 rows that map without an accession by name, two (pga62_flo1 and als3_als1, accessions Q5AF41 and Q59L12) have a sequence identical (sha256) to exactly one UniProt entry in the answers cached in Task 5.1. They are keyed by that accession and are in `controls.tsv` (the note says how). A row whose sequence matches two cached entries keeps its local key. One more (hum3) has a sequence identical to a hydrophobin sequence in `analysis/hydrophobin_truth/`; no accession was looked up for it (that needs a UniProt query).

**Why two control files.** The L1 validator requires every accession to have a sequence in `sequences.faa`. The 20 rows with no usable UniProt accession have a sequence only in a JGI file, and JGI sequences are never committed. They are in a second table with the same columns, kept out of `controls.tsv`. They are controls like the others; `build_kiss_fasta.py` reads them back.

**Clusters are a proxy until Task 5.3.** `clusters.tsv` uses `KISS-<slice>-<cluster_number>`. COMMON-RULES rule 4 asks for MMseqs2. After the job, `build_control_set.py --clusters-tsv ... --split-members ...` replaces them and fills `leakage`, `seed_class_slice` and `seed_class_all39`.

## Run order (all commands are NOT run; each HPCC step needs a go-ahead)

1. **Local: the seed FASTA** (network): `python3 -I analysis/wall_fusion_truth/fetch_seeds_models.py --out-dir ~/kiss_5_3_inputs`. It writes `seeds.faa.gz`, `seeds_summary.tsv` (the Pfam version of every alignment: compare with Pfam 38.2) and `seeds.faa.gz.sha256`. A model with no seed alignment on InterPro stops the script; report it, do not skip it.
2. **Local: the sequence FASTA**, outside the repository (the script refuses a path inside it):
   `python3 -I analysis/wall_fusion_truth/build_kiss_fasta.py --jgi-dir to_import/kiss2019_jgi --out-dir ~/kiss_5_3_inputs`.
   With the defaults it writes 293 sequences: 132 controls and the 161 shared negatives (six of the negatives as `<accession>#neg`). It also writes `members.tsv` and the `.sha256` files.
3. **Stage to HPCC** (owner's go-ahead, then): `rsync -a ~/kiss_5_3_inputs/ hpcc:/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_5_3/in/`. The inputs hold JGI sequences: they stay under `_workdir` and are not committed.
4. **HPCC: a clone at the reviewed commit** (owner's go-ahead): the previous clone was deleted on 2026-10-22. `git clone` the repository to `/bigdata/stajichlab/jstajich/projects/adhesionPred-kiss_5_3` and check out the merge commit that holds these scripts. `CODE_COMMIT` is that commit.
5. **HPCC: the job** (owner's go-ahead):
   ```
   W=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_5_3
   sbatch --export=ALL,PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/adhesionPred-kiss_5_3,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/out \
     -o $W/cluster_split.log -e $W/cluster_split.log \
     /bigdata/stajichlab/jstajich/projects/adhesionPred-kiss_5_3/analysis/wall_fusion_truth/cluster_split.sh
   ```
   `$W` must exist and `$W/out` must not. About 300 sequences and 39 models: minutes. The script asks for 4 CPUs, 8 GB and one hour. The pin check covers `analysis/wall_fusion_truth`, the control set and the two model lists.
6. **Checks after the job** (read-only): `sacct` state COMPLETED; `out/provenance.txt` shows the pinned commit, the MMseqs2 version and the input sha256; `out/summary.tsv` has no pooled line; every sequence of `members.tsv` is in `members_split.tsv`; `clusters_split.tsv` has the `mixed_pos_neg` column and the six overlap clusters are `yes`.
7. **Local: bring back the small tables** (accessions only, no sequences) and rebuild the set:
   `python3 -I analysis/wall_fusion_truth/build_control_set.py --clusters-tsv <out>/clusters.tsv --split-members <out>/members_split.tsv`. The L1 validator then checks the MMseqs2 clusters. Only then open a PR with `data/controls/wall-fusion-kiss2019/` and `<out>/clusters_split.tsv`, `summary.tsv`.

## What the job does not do

It does not change a status, tier or call, or run `hmmsearch` (that is Task 5.4). It uses `module load MMseqs2/17-b804f` as `analysis/lectin_validation/run_cluster_seed.sh` does; it does not need `hmmer`.

## Open questions

Only what is left (all need HPCC or the network, so none could be settled here):

1. **MMseqs2 module name and `pident` scale.** `MMseqs2/17-b804f` is the name in COMMON-RULES and the repository scripts; not confirmed for today's HPCC. The `pident` scale of this build is unconfirmed: `finalize_truth.py` treats a file whose values are all 1.0 or less as a fraction, and `seed_split.py` does the same.
2. **InterPro seed availability.** That InterPro serves a seed alignment for all 39 models is unconfirmed; the fetch stops on a missing one, and the Pfam release it serves may differ from 38.2 (`seeds_summary.tsv` records it).
3. **Taxon names against `names.dmp`.** The names in the control set come from NCBI answers cached in 5.1 and were not checked against `names.dmp` (`python3 -m cellsurface_sorting_hat.controls <dir> --names-dmp names.dmp`).

## How Task 5.4 uses the outputs

`baseline_recall.py` (Task 5.4, not written) reads `clusters_split.tsv` and `members_split.tsv`. Per slice (adhesion, wall, fusion; never pooled) it reports recall of the Pfam models per cluster and per held-out clade, with these columns of the set (both files): `clade_fine` (primary level) and `clade_coarse` (second test); `cluster_stratum_fine` and `cluster_stratum_coarse` (shared or clade-only); `jgi_accept` (filter later); `seed_class_slice` and `seed_class_all39` (report `seed-like` and `seed-independent` apart). A stratum is `estimated` only with 5 or more clusters of its own, otherwise `descriptive`; fusion is always descriptive. Clusters with `mixed_pos_neg=yes` are flagged and can be filtered. The family representatives appear as a flagged group (`family_reps.tsv`), never count toward a floor, and are never test controls. The wall rows are a descriptive false-hit panel for the narrow call, never negatives, never truth for the adhesion call. The table of clusters that Pfam misses names the PR 6 candidates.
