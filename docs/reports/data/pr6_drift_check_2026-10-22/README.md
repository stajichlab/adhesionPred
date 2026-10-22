# PR 6 drift check: run record of job 30282220 (2026-10-22)

Small tables of the drift check of the three queries whose outside homologs come only from `jackhmmer` iteration:
aga1 (KS0080, P32323), ece1 (KS0107, Q07730) and hyd1_like (KS0017, A0ACE4KFX7). The script, its freeze (thresholds) and its
classes are in `analysis/wall_fusion_truth/pr6_drift_check.py`; the job is `pr6_drift_check.sbatch`; the design is in
`analysis/wall_fusion_truth/README_pr6_drift.md`. The report is `docs/reports/2026-10-22-pr6-drift-check.md`.

| Item | Value |
|---|---|
| Job | 30282220 (`pr6_drift_check`), COMPLETED, exit 0, 3 min 31 s, MaxRSS 75 MB |
| Submit | `sbatch -p batch --time=4:00:00 -c 4 --mem=8G` (partition batch; the request was a guess) |
| Code pin | `e333a66565fdd21c0fc95743958b08d013a5efed` (HPCC clone `adhesionPred-kiss53`, clean tree) |
| Target | `/srv/projects/db/Uniprot/2025_03/uniref50.fasta`, 28,781,193,344 bytes (equals the round 2 file by size; sha256 of round 2 `bc52bdd6c97ba2cc365b7e42ae7abb34865e0a1a82a2f95a74590aae040c81ff`, not recomputed) |
| Input | `kiss_sequences.faa` sha256 `66af609eeea24212f7813918fd5957489ded4e4e45792deadd955a0a93bacdf6` |
| Tools | HMMER 3.4, seqkit 2.4.0, Python 3.12.13 |
| Settings | `phmmer -E 10 --domE 10` (pairwise, real and shuffled hits); shuffle seed 1 |
| Freeze (PROPOSALS) | E <= 1e-5, cover >= 0.5 (the LARGER of query and hit coverage), null ratio < 0.5, rule 5 supported sequences in >= 2 named clades |
| Hits | 233 jackhmmer outside fungal hits (43 + 64 + 126), read from the committed round 2 record |

## Files

| File | Content |
|---|---|
| `drift_summary.tsv` | per query: hits, supported, bias_driven, unsupported, supported named clades, verdict |
| `drift_per_clade.tsv` | per query and fine clade: the three classes |
| `drift_hits.tsv` | one line per hit: class, pairwise E-value, query and hit coverage, real and shuffled score, composition columns |
| `hits_manifest.tsv` | the 233 hits taken from the round 2 record (ids, clades, taxid, E-value, score) |
| `summary.txt`, `provenance.txt` | the job's own summary and provenance |

SHA-256 of the committed files:

| File | sha256 |
|---|---|
| `drift_hits.tsv` | `6d78bbbb412bf0d11d98904c7aff4e48c297a826a72e0069ea32b168fe1eb4c9` |
| `drift_per_clade.tsv` | `b9bb2d8bb7c98a885add2d49c3b0d95da7266ee90feb94e25ed207bd516b6df2` |
| `drift_summary.tsv` | `51d1855967caf8f5b548d468adf5daac7ba30519ff37152a0e97d392065cfa2c` |
| `hits_manifest.tsv` | `f36df45bef919ee12c978cfb575265b724e3451d38c28f2e53b1b8a2e0b4fba2` |
| `provenance.txt` | `443fcebaedde592d4633151e64114f235724a90f9036995d966100c0e4c2c9c7` |
| `summary.txt` | `c36d2e76903b6f31ef3b8f79b4cfd939de2fe1575cb80d42b0dc4896b7e6fc00` |

No sequence is in any file. The hit sequences were extracted on HPCC into the job's scratch directory and deleted with it.
