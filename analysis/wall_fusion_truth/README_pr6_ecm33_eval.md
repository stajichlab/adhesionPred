# Independent-set evaluation of the frozen custom Ecm33 'all' model (DRAFT, nothing run)

Drafted 2026-10-23 by Claude Code (claude-sonnet-5-5). **Nothing was run: no HPCC job, no SLURM submission, no model scored.** The custom model stays **report-only**: nothing here
writes or changes a call, a status, a tier, a panel, `narrow_ok`, `family_measure.tsv` or `family_table.tsv`. Scoping: `docs/reports/2026-10-23-ecm33-independent-set-scoping.md`.
Build record (why this evaluation exists: u = 0 because the alignment-seeded iteration comparator tied the model): `docs/reports/2026-10-23-pr6-ecm33-build.md`.

**Every design choice below is a PROPOSAL for the owner** and is written in `pr6_ecm33_eval_freeze_PROPOSED.json` (`"frozen": false`). Decided already (not reopened): custom model report-only; `min_u` 3 and the
ceil(u/2) rule; the model and the cutoff 29.6 bits of the build freeze; the wish for comparators that do not use the training alignment.

## What PANTHER is, and what recall here means

**PANTHER PTHR31018 is itself a profile HMM** (built on a phylogenetic tree of reference proteomes, scored by E-value). Recall on this set therefore measures **agreement between two profile
models** (the custom Ecm33 model and PANTHER), not truth, not Ecm33 function and not adhesion. A protein that both call could still be a non-Ecm33 member of the wider family (the scoping lists about 330
"Receptor L-domain domain-containing protein" rows). The independence is from the training search (phmmer of Q5AGC4 against UniRef50), not from sequence similarity as a method.

## Files

| File | Purpose |
|---|---|
| `pr6_ecm33_eval_members.py` | stdlib; `select` (test clusters and id lists from `candidates.tsv`), `separate` (cluster30 table from the MMseqs2 clustering), `purity` (the 481 training members versus PANTHER-called clusters). Ids and metadata only. |
| `pr6_ecm33_eval_extract.sbatch` | DRAFT job: selects, streams UniRef50 once with `seqkit grep -f`, runs `mmseqs easy-cluster` of training plus test sequences, writes `test_clusters.faa` (HPCC only), `separation.tsv`, id lists, `extract.sha256`. |
| `pr6_ecm33_eval_search.sbatch` | DRAFT job: all six searches with fixed thresholds; writes the hit tables, `run.json` (freeze sha256 and sha256 of each result file) and `provenance.txt`. |
| `pr6_ecm33_eval_report.py` | stdlib; refuses without a signed freeze and matching sha256 values; per clade, tier, subfamily tables, u, the rule, comparators side by side. |
| `pr6_ecm33_eval_freeze_PROPOSED.json` | the proposed freeze (`"frozen": false`). The signed `pr6_ecm33_eval_freeze.json` does **not** exist; only the owner creates it. |
| `tests/wall_fusion_truth/test_pr6_ecm33_eval.py` | 28 tests: counts from the real `candidates.tsv`, separation, u and the rule at the boundaries, freeze and sha256 refusals, static checks of the jobs. |

## The design (PROPOSALS)

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
  `iteration_alignment_seeded` = the build job's exact `jackhmmer -N 3 -E 1e-5 --incE 1e-5 --incdomE 1e-5` seeded by the build alignment of the all-model (no `-Z`); **new, independent of the training alignment**:
  `phmmer_seed_only` (`-E 1e-5 --domE 1e-5 -Z 70198728`) and `jackhmmer_seed_only` (`-N 3 ... -Z 70198728`), both of Q5AGC4 alone.
- **-Z 70198728** (UniRef50 2025_03 sequence count) makes a seed-only E-value the one a UniRef50 search would give, as in the hyd1_like design (decision 10). The iterative and phmmer methods search the test sequences plus a
  background (161 shared negatives and the four kept proteomes, the build's iteration target); the training members are in **no** target. A hit is E <= 1e-5. The alignment-seeded comparator keeps the build's options
  (no `-Z`): its E-values come from a smaller database, so it is the more permissive comparator, which is the stricter test for the custom model (open choice 5).
- **u** = clusters missed by **every** comparator (`pfam_ga`, `pfam_relaxed`, `iteration_alignment_seeded`, `phmmer_seed_only`, `jackhmmer_seed_only`). A group (a clade, or the whole headline set) is **testable**
  when u >= 3 and the custom model recovers >= ceil(u/2) of u. The report prints the rule, the plain per-comparator recall (counts) and the per-clade table. No equal-false-hit-budget view (not asked).
- **False hits**: the frozen cutoff only; the build's numbers (4 extra calls in 32,491 proteome proteins, call rate 0 in every hard-negative group; P38248 and Q4WNS8 are true members among the shared negatives and are ignored by the build freeze) are quoted, not rerun.
- **Purity check** (descriptive, no job): `purity` reads `candidates.tsv` and the committed `members.tsv`; **164 of 481** training members are in no PANTHER-called UniRef50 cluster. It is also run inside the extract job
  and written to `purity.json`, `purity_members.tsv`.
- Every cell is labelled `headline` or `descriptive`; no interval, no `estimated`. Clusters of one clade are related, not independent.

## Run order (each step needs the owner's go-ahead, D21; exact commands, NOT run)

0. Merge the PR that holds this work. Owner: `cp pr6_ecm33_eval_freeze_PROPOSED.json pr6_ecm33_eval_freeze.json`, edit what is not accepted, set `"frozen": true` and `"owner_signature"`, commit (`freeze: ecm33 independent-set evaluation (owner)`).
1. Update the HPCC clone to the merge commit (`git -C $PROJ_ROOT pull --ff-only`), then:
   ```
   sbatch --export=ALL,PROJ_ROOT=$PROJ_ROOT,CODE_COMMIT=<sha>,OUT=$WORK/ecm33_eval_extract \
     -o $LOGS/%x.%j.log -e $LOGS/%x.%j.log analysis/wall_fusion_truth/pr6_ecm33_eval_extract.sbatch
   ```
   Read `provenance.txt` (counts 187 and 109; how many separated) and `separation.tsv` before step 2.
2. ```
   sbatch --export=ALL,PROJ_ROOT=$PROJ_ROOT,CODE_COMMIT=<sha>,EXTRACT=$WORK/ecm33_eval_extract,OUT=$WORK/ecm33_eval_search,PFAM_RELEASE=38.2 \
     -o $LOGS/%x.%j.log -e $LOGS/%x.%j.log analysis/wall_fusion_truth/pr6_ecm33_eval_search.sbatch
   ```
3. Copy the small tables to a local clone (never the FASTA), then:
   ```
   python3 -I analysis/wall_fusion_truth/pr6_ecm33_eval_report.py --results <search_out> --clusters <extract_out>/test_clusters.tsv \
     --separation <extract_out>/separation.tsv --freeze analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json \
     --cutoff-json data/custom_models/ecm33_all_2026-10-23/cutoff.json --out <report_out>
   ```
   Outputs: `groups.tsv` (column x group type x group: n, recall of the custom model and of each comparator, any comparator, u, recovered of u, testable, rule), `clusters.tsv`, `u_clusters.tsv`, `not_separated.tsv`, `summary.txt`.

## Checks the scripts make

Pinned code (`git diff --quiet $CODE_COMMIT HEAD`), signed freeze, UniRef50 file size, training FASTA sha256 against the freeze, 187 ids found, no id shared by training and test, model and cutoff sha256 against the freeze,
relaxed cutoff and `-Z` against the freeze, alignment ids equal to the 481 committed training ids, PF12454 fetched by accession, `extract.sha256` verified by the search job, `run.json` (freeze sha256, sha256 of every hit file) verified by the report,
no result id outside the test set.

## Unconfirmed items

- Resource requests (partition `batch`, time, memory, cores) are guesses; no run time was measured. The extract pass over 28.8 GB follows the build's extract job.
- `mmseqs easy-cluster` on 668 sequences and `jackhmmer -Z` are not run by the author; the alignment-seeded `jackhmmer` and `esl-reformat` step is the build job's, which ran.
- The build alignment `all/none/train.aln.fa` is on HPCC only and has no committed hash: its ids are checked against the committed manifest and its sha256 is recorded, not pinned.
- Whether a representative of a PANTHER cluster is a good stand-in for its PANTHER proteins is not known (a representative can lack a PANTHER call).
- PANTHER and UniProt counts are of 2026-10-23 and change.

## Owner choices (with recommendations)

1. Unit and sequence: the cluster representative, one per cluster (recommended; it matches the training members). Alternative: one PANTHER protein per cluster, needing a fetch by accession.
2. Headline: PANTHER E <= 1e-10 and cluster30-separated (recommended); descriptive: E cut only, and all 187.
3. Separation by shared MMseqs2 cluster (the freeze's cluster30 rule, recommended) versus an order-independent `mmseqs easy-search` of test against training at 30 percent / 50 percent coverage (not drafted; add if you want the stricter rule).
4. `-Z 70198728` and the background (negatives plus kept proteomes) for the seed-only searches (recommended).
5. Alignment-seeded comparator with the build's exact options, no `-Z` (recommended for comparability and because it is the stricter test); alternative: add `-Z` for all three iterative comparators.
6. Cluster attribute rules (modal clade, best tier, SF12 if any) (recommended; 2 and 35 clusters are affected).
7. Which numbers count as the verdict: the rule at the whole-headline level first, then per clade only where u >= 3; with 109 headline clusters in 18 clades expect most clades to be untestable. Recommended: read the whole-headline row and the u list first.
