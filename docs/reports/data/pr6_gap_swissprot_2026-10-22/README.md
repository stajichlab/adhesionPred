# PR 6 gap search, round 1 on Swiss-Prot: run record of job 30260197

Small tables written by `analysis/wall_fusion_truth/pr6_gap_search.sbatch` and `pr6_gap_summary.py` on HPCC (SLURM job 30260197, COMPLETED, 2 min 16 s, exit 0, peak memory about 1.0 GB). Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md` (decisions 3, 5, 7, S3). Report: `docs/reports/2026-10-22-pr6-gap-search-swissprot.md`.

## What the run is tied to

| item | value |
|---|---|
| code and data commit pin | `6f8aa839d10159549c7240bb0e28a65e4afac2cb` (HPCC clone `adhesionPred-kiss53`, clean tree) |
| HMMER | 3.4 (module `hmmer/3.4`) |
| Pfam | 38.2, `Pfam-A.hmm` sha256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`, same as the earlier jobs |
| UniProt target | `/srv/projects/db/Uniprot/2025_03/uniprot_sprot.fasta` (Swiss-Prot, reviewed entries only), 573,661 sequences, 286,776,098 bytes, sha256 `6a9a6c21d0d3f5074a0d6487b6f07e0861d915e6783167b3ff4a73f50f6c7ed4` |
| taxdump | `/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/sorting_hat/taxdump` (`nodes.dmp` sha256 `78392a69…`, `names.dmp` sha256 `57a92ede…`) |
| input | `kiss_sequences.faa` sha256 `66af609eeea24212f7813918fd5957489ded4e4e45792deadd955a0a93bacdf6` |
| settings | `phmmer -E 1e-5`; `jackhmmer -N 3 --incE 1e-5 --incdomE 1e-5`; Ecm33 check `hmmsearch --max -E 10` |
| rule (provisional, decision 7) | at least 5 outside fungal sequences in total over at least 2 clades |

## Files

| file | what it is |
|---|---|
| `summary.txt` | the summary the script printed |
| `gap_per_query.tsv` | one line per query: counts per method, status, reason |
| `gap_clades.tsv` | outside hits per query, method and fine clade |
| `gap_homologs.tsv` | one line per hit: accession, E-value, score, taxid, clade, where (`self`, `inside_clade`, `outside_clade`, `non_fungal`) |
| `ecm33_check.tsv` | PF12454 against Q5AGC4 and P38248: score against the model's GA thresholds, verdict |
| `provenance.txt`, `provenance.tsv` | job, pin, hashes, settings |

All files hold accessions, numbers and labels only. No sequence is committed (a test checks it).

## Left on HPCC (not committed)

`.../_workdir/kiss_cluster_split/gap_swissprot_out/raw/` (the `phmmer` and `jackhmmer` outputs of the four queries) and `ecm33.domtbl` (1.7 KB). The tables above are what the summary script made from them.
