# Independent-set evaluation of the frozen custom Ecm33 'all' model (DRAFT, nothing run)

Drafted 2026-10-23 by Claude Code (claude-sonnet-5-5). **Nothing was run: no HPCC job, no SLURM submission, no model scored.** The custom model stays **report-only**: nothing here
writes or changes a call, a status, a tier, a panel, `narrow_ok`, `family_measure.tsv` or `family_table.tsv`. Scoping: `docs/reports/2026-10-23-ecm33-independent-set-scoping.md`.
Build record (why this evaluation exists: u = 0 because the alignment-seeded iteration comparator tied the model): `docs/reports/2026-10-23-pr6-ecm33-build.md`.

The design is written in `pr6_ecm33_eval_freeze_PROPOSED.json` (`"frozen": false`). **Decided by the owner on 2026-10-23** (applied here): (1) the test unit is the UniRef50 cluster representative; (2) separation from the 481 training members by the shared-MMseqs2-cluster rule (cluster30, `--min-seq-id 0.3 -c 0.5`; no order-independent easy-search is added); (3) `-Z 70198728` on **all three** iterative comparators, including the alignment-seeded jackhmmer; (4) cluster attributes: clade = modal clade, tier = T2 if any protein of the cluster is T2, subfamily = SF12 if any protein is SF12 (they only decide which row a cluster appears in); (5) the verdict applies to the whole headline set (below). Decided earlier: custom model report-only; `min_u` 3 and the ceil(u/2) rule; the model and the cutoff 29.6 bits of the build freeze; comparators that do not use the training alignment. The remaining items stay **proposals** (listed under `proposals` in the file) until the owner signs.

## What PANTHER is, and what recall here means

**PANTHER PTHR31018 is itself a profile HMM** (built on a phylogenetic tree of reference proteomes, scored by E-value). Recall on this set therefore measures **agreement between two profile
models** (the custom Ecm33 model and PANTHER), not truth, not Ecm33 function and not adhesion. A protein that both call could still be a non-Ecm33 member of the wider family (the scoping lists about 330
"Receptor L-domain domain-containing protein" rows). The independence is from the training search (phmmer of Q5AGC4 against UniRef50), not from sequence similarity as a method.

## Files

| File | Purpose |
|---|---|
| `pr6_ecm33_eval_members.py` | stdlib; `select` (test clusters and id lists from `candidates.tsv`), `separate` (cluster30 table from the MMseqs2 clustering), `purity` (the 481 training members versus PANTHER-called clusters). Ids and metadata only. |
| `pr6_ecm33_eval_fetch_reps.py` | stdlib, LOCAL, read-only: fetches the 187 UniRef50 representatives by id from the public UniProt REST service into a directory OUTSIDE the repo (`reps.faa`, `reps.faa.sha256`, `provenance.txt`); fails loudly on any missing, merged or obsolete id. |
| `pr6_ecm33_eval_extract.sbatch` | DRAFT job: selects, verifies the staged `reps.faa` (sha256 sidecar, 187 ids exactly once), runs `mmseqs easy-cluster` of training plus test sequences, writes `test_clusters.faa` (HPCC only), `separation.tsv`, id lists, `extract.sha256`, `run.json`. |
| `pr6_ecm33_eval_search.sbatch` | DRAFT job: all six searches with fixed thresholds; writes the hit tables, `run.json` (freeze sha256 and sha256 of each result file) and `provenance.txt`. |
| `pr6_ecm33_eval_report.py` | stdlib; refuses without a signed freeze and matching sha256 values; the verdict line (whole headline set), per clade, tier, subfamily descriptive tables, u, comparators side by side, Clopper-Pearson intervals in every row. |
| `pr6_ecm33_eval_freeze_PROPOSED.json` | the proposed freeze (`"frozen": false`). The signed `pr6_ecm33_eval_freeze.json` does **not** exist; only the owner creates it. |
| `tests/wall_fusion_truth/test_pr6_ecm33_eval.py` | tests (including Clopper-Pearson reference values, the verdict row and `-Z` on all three iterative commands): counts from the real `candidates.tsv`, separation, u and the rule at the boundaries, freeze and sha256 refusals, static checks of the jobs. |

## The design

- **Test set**: `candidates.tsv` rows with `overlap_with_training == new_cluster`: **725 proteins in 187 UniRef50 clusters** (recounted by the script and the tests).
- **Unit**: one UniRef50 cluster, tested by the cluster's **representative sequence** (the UniRef50 FASTA record of the cluster id), exactly like the training members. The representative need not be one of the
  PANTHER proteins of the cluster; recall is of the representative. Cluster attributes: clade = most frequent clade of its proteins (2 clusters have two clades, so a clade count can differ by one
  from the scoping table, for example Eurotiomycetes 18 against 19); tier = best tier (35 clusters mix T2 and T3); subfamily = SF12 if any protein is SF12 (4 clusters; 7 proteins).
- **Columns**: HEADLINE = clusters with a protein of PANTHER E <= 1e-10 (**563 proteins, 109 clusters**) that are also **separated** (no training member in the same MMseqs2 30 percent identity /
  50 percent coverage cluster); descriptive 1 = the same E cut without the separation; descriptive 2 = all 187 clusters (725 proteins). Not separated clusters are listed (`not_separated.tsv`).
  A different UniRef50 cluster does not guarantee under 30 percent identity, hence the separation. Its limit: cluster membership is a greedy heuristic, a pair under the threshold can sit in different clusters.
- **Strata** (reported for every column): fine clade, tier (T2 named, T3 PANTHER-only; T1 has no new cluster), PANTHER subfamily (SF12 apart from SF3).
- **Methods** (fixed thresholds, one pass, nothing retuned): `custom` = frozen all-model, hit at full-sequence score >= 29.6 (`hmmsearch -T 29.6`; the job keeps hits down to 0.0 and the report applies the
  cutoff, as the build did, so near-miss scores can be read but retune nothing); `pfam_ga` = PF12454 `--cut_ga`; `pfam_relaxed` = PF12454 at 15.0 bits (the build's relaxed cutoff);
  `iteration_alignment_seeded` = the build job's exact `jackhmmer -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5` seeded by the build alignment of the all-model (with `-Z 70198728`, as the other two); **new, independent of the training alignment**:
  `phmmer_seed_only` (`-E 1e-5 --domE 1e-5 -Z 70198728`) and `jackhmmer_seed_only` (`-N 3 ... -Z 70198728`), both of Q5AGC4 alone.
- **-Z 70198728** (UniRef50 2025_03 sequence count) makes a seed-only E-value the one a UniRef50 search would give, as in the hyd1_like design (decision 10). The iterative and phmmer methods search the test sequences plus a
  background (161 shared negatives and the four kept proteomes, the build's iteration target); the training members are in **no** target. A hit is E <= 1e-5. All three iterative comparators use `-Z 70198728` (owner decision 2026-10-23), so their E-values are on the same database size.
- **u** = clusters missed by **every** comparator (`pfam_ga`, `pfam_relaxed`, `iteration_alignment_seeded`, `phmmer_seed_only`, `jackhmmer_seed_only`).
- **Verdict (owner decision 2026-10-23)**: one row only, the **whole headline set** (PANTHER E <= 1e-10 and separated, about 109 clusters before separation), `min_u` 3. The model **passes** when u >= 3 and it recovers at least
  ceil(u/2) of the u missed clusters (u 3 needs 2; u 4 needs 2; u 5 needs 3). With u < 3 the verdict line says `not testable (u < 3)`. Per-clade, per-tier and per-subfamily rows, and the two descriptive columns, are **descriptive** and never a verdict.
- **Every row** (verdict row and descriptive rows) prints n clusters, u, the number of u recovered, the plain per-comparator recall, and the exact two-sided 95 percent Clopper-Pearson interval for the recovered fraction (recovered of u)
  and for each recall (hits of n), 3 decimals (stdlib only; lower bound 0 for k = 0, upper bound 1 for k = n). **At small u the interval is wide and a wide row is not a result** (2 of 3 gives about 0.094 to 0.992).
  No equal-false-hit-budget view (not asked).
- **False hits**: the frozen cutoff only; the build's numbers (4 extra calls in 32,491 proteome proteins, call rate 0 in every hard-negative group; P38248 and Q4WNS8 are true members among the shared negatives and are ignored by the build freeze) are quoted, not rerun.
- **Purity check** (descriptive, no job): `purity` reads `candidates.tsv` and the committed `members.tsv`; **164 of 481** training members are in no PANTHER-called UniRef50 cluster. It is also run inside the extract job
  and written to `purity.json`, `purity_members.tsv`.
- Every cell is labelled `headline` or `descriptive`; the intervals are Clopper-Pearson and treat clusters as independent, which clusters of one clade are not.

## Why the sequences are fetched live

The first extract job (30324005, code e81f3bd) failed safely: `FATAL: found 51 of 187 representatives in /srv/projects/db/Uniprot/2025_03/uniref50.fasta`. The 187 UniRef50 cluster ids of `candidates.tsv` were taken from the
LIVE UniProt on 2026-10-23, and the HPCC file is the 2025_03 release, so 136 ids do not exist there. The signed freeze defines the test sequence as the cluster representative (UniRef50 FASTA record) and does not pin a sequence file, so the
representative is fetched by cluster id from the live UniRef50 service. Caveat: the representative of a live cluster may differ from the 2025_03 one (and a live cluster may have been merged since its id was listed). The test set is
therefore the live PANTHER / UniRef50 view of 2026-10-23, recorded with the fetch date, the service release header and the sha256 of `reps.faa`. The `-Z 70198728` value is the 2025_03 count and stays as frozen.

## Run order (each step needs the owner's go-ahead, D21; exact commands, NOT run)

0. Merge the PR that holds this work. Owner: `cp pr6_ecm33_eval_freeze_PROPOSED.json pr6_ecm33_eval_freeze.json`, edit what is not accepted, set `"frozen": true` and `"owner_signature"`, commit (`freeze: ecm33 independent-set evaluation (owner)`).
1. Fetch locally (public read-only API; only the ids of the committed table are sent; no owner go-ahead is needed for this read; the output directory must be outside the repository):
   ```
   python3 -I analysis/wall_fusion_truth/pr6_ecm33_eval_fetch_reps.py --out-dir ~/ecm33_eval_reps
   ```
   It stops with a list if any id does not return its own record. Read `provenance.txt` (187 requested and received, service release) first.
2. Stage to HPCC (needs the owner's go-ahead; sequences go to the work directory only, never into the repository):
   ```
   scp ~/ecm33_eval_reps/reps.faa ~/ecm33_eval_reps/reps.faa.sha256 ~/ecm33_eval_reps/provenance.txt <hpcc>:$WORK/ecm33_eval_reps/
   ```
3. Update the HPCC clone to the merge commit (`git -C $PROJ_ROOT pull --ff-only`), then the extract job:
   ```
   sbatch --export=ALL,PROJ_ROOT=$PROJ_ROOT,CODE_COMMIT=<sha>,OUT=$WORK/ecm33_eval_extract,REPS_FASTA=$WORK/ecm33_eval_reps/reps.faa \
     -o $LOGS/%x.%j.log -e $LOGS/%x.%j.log analysis/wall_fusion_truth/pr6_ecm33_eval_extract.sbatch
   ```
   Read `provenance.txt` (counts 187 and 109; how many separated) and `separation.tsv` before the search.
4. ```
   sbatch --export=ALL,PROJ_ROOT=$PROJ_ROOT,CODE_COMMIT=<sha>,EXTRACT=$WORK/ecm33_eval_extract,OUT=$WORK/ecm33_eval_search,PFAM_RELEASE=38.2 \
     -o $LOGS/%x.%j.log -e $LOGS/%x.%j.log analysis/wall_fusion_truth/pr6_ecm33_eval_search.sbatch
   ```
5. Copy the small tables to a local clone (never the FASTA), then:
   ```
   python3 -I analysis/wall_fusion_truth/pr6_ecm33_eval_report.py --results <search_out> --clusters <extract_out>/test_clusters.tsv \
     --separation <extract_out>/separation.tsv --freeze analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json \
     --cutoff-json data/custom_models/ecm33_all_2026-10-23/cutoff.json --out <report_out>
   ```
   Outputs: `groups.tsv` (column x group type x group: n, recall and interval of the custom model and of each comparator, any comparator, u, recovered of u with interval, testable and rule on the verdict row only), `clusters.tsv`, `u_clusters.tsv`, `not_separated.tsv`, `summary.txt` (verdict line first, then every row with its intervals).

## Checks the scripts make

Pinned code (`git diff --quiet $CODE_COMMIT HEAD`), signed freeze, staged `reps.faa` against its sha256 sidecar (recorded in `run.json`), training FASTA sha256 against the freeze, 187 ids found exactly once and equal to the selection, no id shared by training and test, model and cutoff sha256 against the freeze,
relaxed cutoff and `-Z` against the freeze, alignment ids equal to the 481 committed training ids, PF12454 fetched by accession, `extract.sha256` verified by the search job, `run.json` (freeze sha256, sha256 of every hit file) verified by the report,
no result id outside the test set.

## Unconfirmed items

- Resource requests (partition `batch`, time, memory, cores) are guesses; no run time was measured. The extract job no longer streams UniRef50; its resources are lowered guesses. The live fetch (rate, retries) was not run by the author; the URL form was checked once on one id.
- `mmseqs easy-cluster` on 668 sequences and `jackhmmer -Z` are not run by the author; the alignment-seeded `jackhmmer` and `esl-reformat` step is the build job's, which ran.
- The build alignment `all/none/train.aln.fa` is on HPCC only and has no committed hash: its ids are checked against the committed manifest and its sha256 is recorded, not pinned.
- Whether a representative of a PANTHER cluster is a good stand-in for its PANTHER proteins is not known (a representative can lack a PANTHER call).
- PANTHER and UniProt counts are of 2026-10-23 and change.

## Owner choices still open (proposals, with recommendations)

Decisions 1 to 5 (unit, cluster30 separation, `-Z` on all three iterative comparators, attribute rules, verdict on the whole headline set) are made and no longer listed. Open until the owner signs:

1. Headline PANTHER E cut 1e-10 (recommended), with the descriptive columns (E cut only; all 187).
2. The background of the seed-only searches (negatives plus kept proteomes) and `evalue_db_size` 70198728.
3. Hit rule E <= 1e-5 for the iterative and phmmer comparators; relaxed Pfam cutoff 15.0 bits from the build record.
4. The strata (clade, tier, subfamily) as descriptive rows.
