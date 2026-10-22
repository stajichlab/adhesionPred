# PR 6 drift check of the iterated-only queries (DRAFT, not run on HPCC)

Owner choice c of 2026-10-22: build ecm33 only, treat hyd1_like as a known gap under the strict reading, and check whether
aga1 and ece1 can ever be built. After the clade-bucket fix (`README_pr6_gap.md`, section "Clade rule fix") hyd1_like is
`buildable_iterated_only` too, so the check runs on **three queries**: aga1 (KS0080, P32323), ece1 (KS0107, Q07730) and
hyd1_like (KS0017, A0ACE4KFX7).

## The question

The one-pass `phmmer` set does not meet the training rule for these three, the `jackhmmer` set does (aga1 64 outside
sequences in 19 clades, ece1 126 in 23 named clades, hyd1_like 43 in 3 named clades; `jackhmmer` also adds thousands of
non-fungal hits for aga1 and ece1). Iteration can pull in unrelated proteins that share a composition: Ser/Thr-rich
protein, repeat-rich proteins and short peptide precursors are the usual cases. The raw output holds the final iteration
only (no per-iteration data, no domain coordinates), so the check does **not** read iterations and does **not** run
`jackhmmer` again. It asks of each `jackhmmer` outside fungal hit: does a direct pairwise alignment to the query support it,
or is its score explained by composition?

## Design (`pr6_drift_check.py`; freeze in its docstring, written before any result; all PROPOSALS for the owner)

| step | what | where |
|---|---|---|
| 1 `extract-ids` | the `jackhmmer` outside fungal hits of the three queries, read from the committed round 2 record, with a flag for hits that `phmmer` also found | local or job |
| 2 `seqkit grep -f` | one streaming pass over `uniref50.fasta` (28.8 GB) pulls the sequences of those hits; every id must come back exactly once | job |
| 3 `split` | per query: the query, its hits, and the hits with residues shuffled (seed 1, one stream per query) | job |
| 4 `phmmer -E 10 --domE 10 --domtblout` | query against its hits (real) and against the shuffled hits (null); seconds | job |
| 5 `summarize` | class of every hit; the training rule on **supported** hits only | job |

Classes of a hit: `supported_homolog` (best domain i-Evalue <= 1e-5 AND the aligned part of the query OR of the hit >= 0.5 of
its length AND the shuffled-null full-sequence score < 0.5 of the real score), `bias_driven` (E passes, the null score is
>= 0.5 of the real one: a shuffled hit scores about as well, so composition explains it), `unsupported` (no pairwise hit at
E 1e-5, or coverage below 0.5). Composition columns (fraction of S/T/N/Q/G/P and a low-complexity fraction) are
descriptive only: the classifier does not use them. The verdict of a query is `supported_rule_met` (the supported hits
meet decision 7: at least 5 sequences in at least 2 **named** clades) or `not_supported`. The script never promotes a
query: the owner decides.

## Files

`pr6_drift_check.py` (standard library only), `pr6_drift_check.sbatch` (DRAFT), `tests/wall_fusion_truth/test_pr6_drift_check.py`.

## Run order (exact commands, NOT run; each HPCC step needs the owner's go-ahead)

1. Merge the PR that holds these scripts (and #114, the round 2 record, if not yet merged); update the HPCC clone
   (`git -C $B/adhesionPred-kiss53 pull --ff-only origin main`; the job's pin check needs the clone to hold the record).
2. Submit (no new staging: `kiss_sequences.faa` is already staged; UniRef50 is read in place):

```
B=/bigdata/stajichlab/jstajich/projects
W=$B/adhesionPred/_workdir/kiss_cluster_split
sbatch -p batch --time=4:00:00 -c 4 --mem=8G \
  --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<merge sha>,IN=$W/in,OUT=$W/drift_check_out,UNIPROT_FASTA=/srv/projects/db/Uniprot/2025_03/uniref50.fasta \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log \
  $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_drift_check.sbatch
```

The job refuses to run if the pin check fails, the staged FASTA differs from its hash, the UniRef50 file is not the round 2
file (size check; its sha256 is not recomputed), `seqkit` returns a different number of sequences from the number of ids,
or any search lacks the `# [ok]` trailer. It never overwrites `$OUT`.

## Expected outputs (`$OUT`)

`drift_summary.tsv` (one row per query: counts per class, supported named clades, verdict), `drift_hits.tsv` (one row per
hit: class, domain E-value, query and hit cover, real and null score, composition), `drift_per_clade.tsv`, `summary.txt`,
`hits_manifest.tsv`, `provenance.txt`. No sequence is written. All tables are small; a results record PR follows the pattern
of #112 and #114.

## What the verdict means

aga1, ece1 and hyd1_like stay **not buildable** unless their supported hits meet the rule. `supported_rule_met` would make
the query a candidate for the build step (6.1) with the supported hits as members, not a model. A `bias_driven` majority
confirms the drift suspicion. Neither outcome changes any call, status, tier or `narrow_ok`.

## Open items (unconfirmed without HPCC)

- `seqkit grep -f` on a 70-million-record file: time and memory are not measured; the 4 h, 4 core, 8 GB request is a guess.
- The `phmmer --domtblout` columns are parsed by the documented HMMER 3.4 layout (target name, tlen, query, qlen, full
  score, domain i-Evalue, hmm from/to = query coordinates for phmmer, ali from/to = target coordinates); tested on fixtures
  only, since HMMER is not installed locally. The first real run confirms it.
- The null is one shuffle per hit (seed 1). A hit whose real and null scores are both high is called `bias_driven`; a hit
  with a low real score can slip below the ratio by chance. More shuffles are a possible upgrade.
- The 0.5 coverage, the 1e-5 E-value and the 0.5 null ratio are proposals; the owner may change them before the run (they
  are constants at the top of `pr6_drift_check.py` and a test pins them).

## Result of the run (job 30282220, 2026-10-22)

The run record is in `docs/reports/data/pr6_drift_check_2026-10-22/` and the note in `docs/reports/2026-10-22-pr6-drift-check.md` (a results record PR in the pattern of #112 and #114). hyd1_like: 32 of 43 hits supported in 3 named clades (`supported_rule_met`); aga1: 11 of 64 supported in 8 named clades (`supported_rule_met`) but all 11 pass only through the "larger of query and hit coverage" rule, with query coverage 0.148 to 0.31; ece1: 0 of 126 supported (`not_supported`, drift confirmed). Read with query coverage >= 0.5, aga1 and ece1 are `not_supported` and hyd1_like is unchanged. The choice of the coverage definition is an owner decision; the script and its freeze are not changed by the record.
