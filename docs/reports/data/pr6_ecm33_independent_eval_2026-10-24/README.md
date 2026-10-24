# PR 6 ecm33 independent-set evaluation: run record of jobs 30336905 (extract) and 30343079 (search)

Small tables copied back from HPCC after the independent-set evaluation of the frozen custom Ecm33 `all` model (`data/custom_models/ecm33_all_2026-10-23/`, cutoff 29.6 bits) on the independent PANTHER test set, plus the local report tables. Freeze: `analysis/wall_fusion_truth/pr6_ecm33_eval_freeze.json`. Design: `analysis/wall_fusion_truth/README_pr6_ecm33_eval.md`. Report: `docs/reports/2026-10-24-pr6-ecm33-independent-eval.md`.

## What the run is tied to

| item | value |
|---|---|
| code pin | `a8d0876afc2bf9a56f52c87d616f7cc8feb848f9` (extract and search jobs, a second clone on HPCC) |
| freeze sha256 | `22e6727879f9847b166b3ccdc75ceb702a87589b504983c5208589ee1dffc270` (signed 2026-10-23, merged #130) |
| first extract job | 30324005, FAILED safely (code `e81f3bd`): only 51 of 187 UniRef50 ids exist in the HPCC 2025_03 file |
| representatives | 187 sequences fetched live by id from the UniProt REST service (release 2026_03, 2026-10-23); `reps.faa` sha256 `493d386e31a87e452deeb5327fd20d4850d3abbce977ff6595c7248049605586` |
| extract job | 30336905, COMPLETED, 15 s |
| search job | 30343079, COMPLETED exit 0, 1 h 40 min |
| settings | custom cutoff 29.6 bits (frozen), relaxed Pfam 15.0 bits, Pfam 38.2 PF12454, HMMER 3.4, MMseqs2 `b804fbe3` (30 percent identity, 50 percent coverage), `min_u` 3, `-Z 70198728`, hit rule E <= 1e-5 for the iterative and phmmer comparators, headline PANTHER E <= 1e-10 |
| model | `model.hmm` sha256 `9a6e2c812540b0b1892f081611bf4044ec674bc524d68f5c6fc3d16001c50937` |
| training FASTA | `ecm33_members.faa` sha256 `b5e9110d7e315d4035ab16ab1ed7f0f294d410bf603affb429f6d8de3f3eb8ee`; training alignment `train.aln.fa` sha256 `037e252fdbab66eb112e93133d2516445c1b61084b753b49d1359aabc7d2550a` (recorded, not pinned) |
| report | run locally from the copied-back tables (draft, report-only) |

## Files

| path | what it is |
|---|---|
| `reps/provenance.txt`, `reps/reps.faa.sha256` | the live fetch (ids requested 187, received 187, service release header) and the hash of `reps.faa` |
| `extract/counts.json` | 187 test clusters (725 proteins), 109 at the strict PANTHER E cut (563 proteins), strata counts |
| `extract/purity.json` | descriptive: 164 of 481 training members are in no PANTHER-called cluster |
| `extract/separation.tsv` | per test cluster: its 30 percent MMseqs2 cluster, training members in it, `separated` |
| `extract/test_clusters.tsv`, `test_cluster_ids.txt`, `headline_pe_cluster_ids.txt` | the 187 test clusters with clade, tier, subfamily, best PANTHER E; id lists |
| `extract/run.json`, `extract.sha256`, `provenance.txt` | pin, freeze and reps hashes; hashes of the extract outputs; job provenance |
| `search/hits_*.tsv` | the six hit tables (custom, Pfam GA, relaxed Pfam, alignment-seeded iteration, seed-only phmmer, seed-only jackhmmer); byte-identical to the HPCC files, hashes in `run.json` |
| `search/run.json`, `provenance.txt` | pin, sha256 of every hit file, freeze hash, HMMER version, Pfam path; job provenance |
| `report/summary.txt` | verdict, headline row, and every descriptive row with exact 95 percent Clopper-Pearson intervals |
| `report/groups.tsv`, `clusters.tsv`, `u_clusters.tsv`, `not_separated.tsv` | group rows; one row per cluster with all hits; the u clusters; the strict clusters that share a 30 percent cluster with training |

All files hold identifiers, scores and labels only. No sequence is committed (a test checks it). The record is about 2.6 MB; `hits_iter_aln.tsv` is kept uncompressed.

## Left on HPCC (not committed)

The test FASTA `test_clusters.faa` (the 187 representatives, same bytes as `reps.faa`), `purity_members` and `test_proteins` tables, `train.aln.fa`, `model.hmm`, and the job logs. `reps.faa` is also in the local copy outside any repo (`/Users/jstajich/projects/adhesionPred-ecm33-results/eval_reps/`), not committed.

Folder names avoid the ignored patterns (`build/`, `results/*`, `input/`, `output/`, `dist/`; `*.csv`, `*.log`).
