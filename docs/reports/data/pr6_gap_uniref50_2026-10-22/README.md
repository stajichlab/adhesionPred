# PR 6 gap search, round 2 on UniRef50: run record of job 30262455

Small tables written by `analysis/wall_fusion_truth/pr6_gap_search.sbatch` and `pr6_gap_summary.py` on HPCC (SLURM job 30262455, COMPLETED, 3 h 19 min 17 s, exit 0, peak memory about 1.6 GB). Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md` (decisions 3, 5, 7, S3). Round 1 (Swiss-Prot): `../pr6_gap_swissprot_2026-10-22/`. Report: `docs/reports/2026-10-22-pr6-gap-search-uniref50.md`.

## What the run is tied to

| item | value |
|---|---|
| code and data commit pin | `e11808b2bf297cc0aa1dfb345dfb3bf568a2c937` (HPCC clone `adhesionPred-kiss53`, clean tree; includes the UniRef50 header patch of #113) |
| submitted as | `sbatch -p batch --time=8:00:00 -c 8 --mem=16G` (the `#SBATCH` defaults of the script, partition `short`, are overridden on the command line) |
| HMMER | 3.4 (module `hmmer/3.4`) |
| Pfam | 38.2, `Pfam-A.hmm` sha256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`, same as the earlier jobs |
| target | `/srv/projects/db/Uniprot/2025_03/uniref50.fasta` (UniRef50, one representative per 50 percent identity cluster), 70,198,728 sequences, 28,781,193,344 bytes, sha256 `bc52bdd6c97ba2cc365b7e42ae7abb34865e0a1a82a2f95a74590aae040c81ff` |
| taxdump | `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat/taxdump` (`nodes.dmp` sha256 `78392a69…`, `names.dmp` sha256 `57a92ede…`) |
| input | `kiss_sequences.faa` sha256 `66af609eeea24212f7813918fd5957489ded4e4e45792deadd955a0a93bacdf6` |
| settings | `phmmer -E 1e-5`; `jackhmmer -N 3 --incE 1e-5 --incdomE 1e-5`; Ecm33 check `hmmsearch --max -E 10` (unchanged from round 1) |
| rule (provisional, decision 7) | at least 5 outside fungal sequences in total over at least 2 clades |

## Files

| file | what it is |
|---|---|
| `summary.txt` | the summary the script printed |
| `gap_per_query.tsv` | one line per query: counts per method, status, reason |
| `gap_clades.tsv` | outside hits per query, method and fine clade (`(no clade)` = a fungal hit whose fine clade could not be assigned) |
| `gap_homologs.tsv.gz` | one line per hit (47,072 lines, gzipped, 0.56 MB): accession (`UniRef50_<representative>`), E-value, score, taxid, clade, where (`self`, `inside_clade`, `outside_clade`, `non_fungal`) |
| `ecm33_check.tsv` | PF12454 against Q5AGC4 and P38248, as in round 1 (the same numbers) |
| `provenance.txt`, `provenance.tsv` | job, pin, hashes, settings |

All files hold accessions, numbers and labels only. No sequence is committed (a test checks it). Each hit is the representative of a UniRef50 cluster; the `n=` member counts of the headers are not used by the rule.

## Left on HPCC (not committed)

`.../_workdir/kiss_cluster_split/gap_uniref50_out/raw/` (the `phmmer` and `jackhmmer` outputs of the four queries, including the per-iteration tables that a drift check would read) and `ecm33.domtbl` (1.7 KB).
