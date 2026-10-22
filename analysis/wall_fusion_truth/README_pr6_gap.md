# PR 6 gap search (steps 6.1 and 6.2) and the Ecm33 check: DRAFT, not run

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md` (decisions 3, 5, 7). Nothing here has run on HPCC; no SLURM job was submitted.*

The whole-Pfam scan (step 6.0, record in `docs/reports/data/pr6_pfam_scan_2026-10-22/`) left four of the 44 missed clusters with no Pfam hit. Decision 3 of 2026-10-22: search for outside homologs of these four, build a model only where the training rule is met, and list the rest as `not_buildable` known gaps. Decision 7 sets the rule: at least 5 sequences in at least 2 clades outside the held-out unit (provisional).

| Cluster | Gene | Query | Slice | Held-out clade (fine) |
|---|---|---|---|---|
| KS0017 | hyd1_like | A0ACE4KFX7 | adhesion | Eurotiomycetes |
| KS0080 | aga1 | P32323 | adhesion | Saccharomycotina |
| KS0107 | ece1 | Q07730 | adhesion | Saccharomycotina |
| KS0086 | ecm33 | Q5AGC4 | wall | Saccharomycotina |

## Files

| File | What it does |
|---|---|
| `pr6_gap_search.sbatch` | The HPCC job (partition `short`, `${SCRATCH:?}`, pin check, Pfam release check). Ecm33 check first, then `phmmer` and `jackhmmer` of each query against the UniProt set, then the summary. |
| `pr6_gap_summary.py` | Standard library only. Reads the `--tblout` files, maps each hit's `OX=` taxon to a clade with the NCBI taxdump, applies the training rule, writes the tables. Writes no call, no status, no tier, no sequence. |
| `tests/wall_fusion_truth/test_pr6_gap_search.py` | Fixtures and hand counts; static checks of the sbatch. |

## Order of work (each HPCC step needs the owner's go-ahead, D21)

1. Merge the PR that holds these scripts. Update the HPCC clone `adhesionPred-kiss53` (fast-forward) and read back `HEAD`.
2. **Stage the UniProt set** (open item 1 below). The staged `kiss_sequences.faa` is already in `.../_workdir/kiss_cluster_split/in/`; nothing else is staged. `P38248` (the *S. cerevisiae* Ecm33 for the check) is read from the committed shared negatives in the clone.
3. Submit:
   ```
   B=/bigdata/stajichlab/jstajich/projects
   W=$B/adhesionPred/_workdir/kiss_cluster_split
   sbatch --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<merge commit>,IN=$W/in,OUT=$W/gap_search_out,PFAM_RELEASE=38.2,UNIPROT_FASTA=<staged path> \
     -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_gap_search.sbatch
   ```
   (not run). `OUT` must not exist.
4. Read `summary.txt`, `ecm33_check.tsv`, `gap_per_query.tsv`; commit the small tables as a run record (not the raw `--tblout` files if large; the per-query `.tbl` files are listed in provenance).

## What the job does

1. **Ecm33 check.** `hmmfetch` of the model `Ecm33` (looked up by name; the job stops if its `ACC` is not `PF12454`), `hmmsearch --max -E 10` against Q5AGC4 and P38248. `ecm33_check.tsv` has, per sequence: its length, the ratio to the 406 aa of P38248, the model's GA thresholds (read from the model file), the full-sequence and best-domain scores, the share of the model covered, and a verdict `pass_GA`, `below_GA` or `no_hit_at_E10`. Report 07 row 22 found no hit of P38248 even at E 10; the check shows whether Q5AGC4 behaves the same or sits just below GA.
2. **Searches.** `phmmer` (one pass) and `jackhmmer -N 3`, both with `--incE 1e-5 --incdomE 1e-5` and reporting at `-E 1e-5 --domE 1e-5`, fixed in the script before any result is read. `jackhmmer`'s `--tblout` holds the last iteration.
3. **Classification of every hit** by the taxon in `OX=` (UniProt headers) or `TaxID=` (UniRef headers): `self`, `no_taxid`, `non_fungal`, `inside_clade` (the held-out fine clade; used only for the descriptive leave-species-out question, decision 10) and `outside_clade`.
4. **The rule** on the `outside_clade` fungal hits (distinct accessions): at least 5 sequences in total, over at least 2 clades. Status: `buildable` (the one-pass set meets it), `buildable_iterated_only` (only the iterated set does; iteration can drift, a flag for the owner) or `not_buildable` (with the reason).

Outputs: `gap_per_query.tsv` (counts per method, status, reason, mode `leave-clade-out` for wall and `leave-cluster-out (partial)` for adhesion), `gap_clades.tsv`, `gap_homologs.tsv` (one line per hit: accessions, E-value, score, taxon, clade; no sequences), `ecm33_check.tsv`, `summary.txt`, `provenance.tsv`, plus `provenance.txt`, `ecm33.domtbl` and `raw/{phmmer,jackhmmer}/<query>.tbl`.

## Round 2 on UniRef50 (the owner's decision: Swiss-Prot first, then UniRef50)

`/srv/projects/db/Uniprot/2025_03/uniref50.fasta` (28.8 GB). Its headers are `>UniRef50_<id> <cluster name> n=<members> Tax=<species> TaxID=<taxid> RepID=<representative>`, not UniProt `OX=` headers. `pr6_gap_summary.py` therefore reads `OX=` or `TaxID=` (if a hit has both and they disagree it gets no taxon, is not counted, is listed in `taxid_conflicts.tsv` and the summary prints a loud banner). The self hit of a query is the cluster `UniRef50_<query>` or `RepID=<query>`.

- **Each hit is the representative of a 50 percent identity cluster (one sequence per cluster).** Counts of outside sequences are therefore conservative against the same search of all UniProt. The header's `n=` member counts exist but are NOT used for the rule.
- A query that is only a member of a cluster whose representative is another protein is not "self": that representative is a homolog at 50 percent identity or more and is classified by its own taxon.
- E-values depend on the size of the target set, so the same cut (`--incE 1e-5`, unchanged) is stricter against UniRef50 than against Swiss-Prot; counts of round 1 and round 2 are not strictly comparable.
- The job counts the sequences of the target set once and hashes it once (minutes on 28.8 GB). Partition, time, cores and memory are defaults on the `#SBATCH` lines; override them on the `sbatch` command line for the longer run (`-p batch --time=... -c ... --mem=...`); nothing in the script fixes them. If the summary fails after the searches, the raw tables are kept in `$OUT.raw_after_summary_failure` (the job's scratch directory is deleted on exit).

## Why phmmer and jackhmmer

The spec allows either. `phmmer` is one pass and cannot drift into another family; it is the set the rule is judged on. `jackhmmer` can reach more distant homologs, which matters for small, fast-evolving proteins such as these, but it can also pull in unrelated proteins. Running both makes the drift visible: a target that passes only after iteration is flagged, not silently accepted. This is my reading; the owner can choose otherwise.

## Assumptions (mine, for the owner to confirm or change)

- **The rule reading.** "At least 5 sequences in at least 2 clades" is read as 5 outside sequences in total, spread over 2 clades. The stricter reading (5 in each of 2 clades) is a descriptive column (`*_strict_5_in_2_clades`).
- **Only fungi count** toward the rule (taxon 4751 among the ancestors). Hits from other kingdoms are counted apart (`*_non_fungal_seq`). A Candida-specific peptide with bacterial look-alikes should not become "buildable" by them.
- **The held-out unit is the cluster's own fine clade** (all four are clade-only), as in decision 4. The query itself is excluded. A hit from the held-out clade is not an outside homolog.
- **Decision 5 asked how many clade-only clusters have any outside homolog.** This job answers it for the four gap clusters only. The other 28 clade-only missed clusters hit a Pfam family and are not searched here. Running the same search for all 32 clade-only clusters would be a separate, owner-approved job.
- **Distinct sequences, not species.** The count is of sequences; `*_outside_species` is a column. Sequences are not independent clusters (a proxy, as everywhere in this plan). The one-proteome-per-species rule (decision 9b) is about proteomes and does not apply to UniProt hits.
- `P38248` comes from the committed shared negatives (it is there as `ECM33 YBR078W`); it is not in the staged `kiss_sequences.faa`.

## Open items (could not be confirmed without HPCC)

1. **The UniProt database path.** No UniProt FASTA path is documented in the repository, so `UNIPROT_FASTA` is a required argument and none is guessed. The owner should look (read-only `ls`) for a UniProt FASTA on HPCC. If there is none, the owner downloads one and stages it. The spec names the reviewed fungal proteins as the first choice (Swiss-Prot restricted to taxon 4751); as I know the public REST API, that is the query `taxonomy_id:4751 AND reviewed:true` on the UniProtKB stream endpoint with `format=fasta` (check this against the live API before relying on it). Fungal TrEMBL is much larger and would need a longer job. Whatever is used, the job records its size, sequence count and sha256, and the report records the download date. The summary needs `OX=` in the headers, which standard UniProt FASTA files carry.
2. **The taxdump.** The default is `$DATA_ROOT/_workdir/sorting_hat/taxdump` (the `TAXDUMP` of `scripts/sorting_hat/refresh_3_10/lib.sh`). The job stops with a clear message if `nodes.dmp` or `names.dmp` is missing there.
3. **Partition limits and memory.** `short`, 8 cores, 16 GB, 2 hours are my estimate (the earlier jobs used 1 hour and 8 GB for much smaller work). Whether `short` allows 2 hours is not confirmed. The time for `jackhmmer` against the staged set is not measured.
4. **`hmmfetch -f` by name for `Ecm33`** follows `scripts/sorting_hat/pfam_hmmsearch.sbatch` and the Task 5.4 job, which fetched all 39 models by name; the job stops if the fetched model is not `PF12454`.
5. **The parsers** follow the documented `--tblout` and `--domtblout` layouts of HMMER 3.4. They were tested on fixtures only (HMMER is not installed locally). The first real run confirms them; the `# [ok]` trailer is required.

## What this does not do

No model is built (that is step 6.3, only for a `buildable` target, with its own scripts and its own go-ahead). No call, status, tier, `narrow_ok` or panel row changes. The comparators of decision 9 are step 6.4.

## Result of round 1 (job 30260197, Swiss-Prot 2025_03)

Run on 2026-10-22 with the owner's go-ahead: code `6f8aa83`, `UNIPROT_FASTA=/srv/projects/db/Uniprot/2025_03/uniprot_sprot.fasta` (reviewed entries only, 573,661 sequences; the open item 1 above is settled for this round), the default taxdump `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat/taxdump` (it existed; open item 2 settled), `short`, 8 cores, 16 GB: COMPLETED in 2 min 16 s, so the 2-hour request was not needed (open item 3 settled for these inputs). The parsers read the real HMMER 3.4 output (open item 5 settled). Result: all four targets `not_buildable` (hyd1_like 1 outside sequence, aga1 0, ece1 0, ecm33 4 outside sequences in 2 clades, one short of the minimum of 5). The Ecm33 check: `PF12454` is below its GA for both Q5AGC4 (score 5.5) and P38248 (0.6) against a GA of 23.2. The record is `docs/reports/data/pr6_gap_swissprot_2026-10-22/`, the note `docs/reports/2026-10-22-pr6-gap-search-swissprot.md`. A second round on UniRef50 2025_03 is the owner's decision; its query set is not chosen yet.

## Result of round 2 (job 30262455, UniRef50 2025_03)

Run on 2026-10-22 with the owner's go-ahead: code `e11808b`, `sbatch -p batch --time=8:00:00 -c 8 --mem=16G`, `UNIPROT_FASTA=/srv/projects/db/Uniprot/2025_03/uniref50.fasta` (70,198,728 sequences, sha256 `bc52bdd6…`), COMPLETED in 3 h 19 min. Statuses: ecm33 `buildable` (379 outside fungal sequences with `phmmer`, 19 clade labels), hyd1_like `buildable` but marginal (5 outside sequences, 1 named clade plus the `(no clade)` bucket that the script counts as a clade), aga1 and ece1 `buildable_iterated_only` (drift flag: not buildable until a drift check is done). Record: `docs/reports/data/pr6_gap_uniref50_2026-10-22/`; report: `docs/reports/2026-10-22-pr6-gap-search-uniref50.md`. UniRef50 is about 122 times larger than Swiss-Prot, so round 1 and round 2 counts are not strictly comparable.

## Clade rule fix after round 2 (2026-10-22, owner choice c)

Round 2 counted the `(no clade)` bucket (fungal hits with no assigned fine clade) as a clade, so hyd1_like was called
`buildable` on a one-pass set of 3 hits in Sordariomycetes plus 2 hits with no assigned clade. `pr6_gap_summary.py` now
counts named clades only for the 2-clade condition and for the strict column; sequences without a clade still count in
the total. `pr6_gap_rederive.py` re-derives the round 2 verdicts from the committed record without searching again; the
result is in `docs/reports/data/pr6_gap_uniref50_2026-10-22_rederived/`. Only hyd1_like changes: the one-pass set no
longer meets the rule (5 sequences, 1 named clade) but the iterated set does (43 sequences in 3 named clades), so it
becomes `buildable_iterated_only`, with the same drift flag as aga1 and ece1. ecm33 stays `buildable`. The committed run
record of job 30262455 is not edited. The drift check of aga1, ece1 and hyd1_like is in `README_pr6_drift.md`.
