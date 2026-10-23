# PR 6 hyd1_like build: run record of jobs 30321490 (pre-check), 30325384 (extract) and 30335767 (build)

Small tables copied back from HPCC after the pre-check, the extract job and the build array of the hyd1_like custom-HMM experiment (query KS0017, A0ACE4KFX7), plus the local report tables. Freeze: `analysis/wall_fusion_truth/pr6_hyd1_freeze.json` (signed 2026-10-23, merged in #132). Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md`. Report: `docs/reports/2026-10-23-pr6-hyd1-like-build.md`. Earlier record (members come from it): `../pr6_drift_check_2026-10-22/`. Same layout as `../pr6_ecm33_build_2026-10-23/`.

## What the run is tied to

| item | value |
|---|---|
| code pin (extract and build) | `9f4843b1f7d3111c1f2c5db8544076b78d2b8a67` (`9f4843b`, in every `units/*/*/provenance.txt` and `extract/provenance.txt`) |
| code pin (pre-check) | `5fb6b4864ce0e573ecd46fc2baf24b7ac6884d24` |
| freeze sha256 | `65ed9bb6dbcc6035a4093537a5afc0f92a179d3f4356902c1aba83a12e552ff2` (file as read by the report; in every `cutoff.json` and build `provenance.txt`) |
| pre-check job | 30321490 (`pr6_hyd1_negcheck`), 3 min 24 s |
| extract job | 30325384 (`pr6_hyd1_extract`), 3 min 37 s: 30 members (query + 29 UniRef50 representatives), 4 named clades, 3 excluded (no named clade), 6 external proteins |
| build job | 30335767, `sbatch --array=0-9%4`, partition `batch`, 8 cores, 16 GB, 4 h per task. 10 tasks (0-9), all COMPLETED exit 0, 11 to 37 s each |
| settings | MAFFT 7.505 L-INS-i, `hmmbuild --amino`, HMMER 3.4, MMseqs2 `b804fbe3` (cluster30: 30 percent identity, 50 percent coverage), Pfam 38.2; `min_u` 3; extra-call limit 5 per 10,000 (16.2455 over 32,491 proteins); hard-negative rate 0.05; `negatives_ignored` empty; `panel_ignored` 19 ids |
| report | run locally: `python3 -I analysis/wall_fusion_truth/pr6_hyd1_report.py --results <R>/hyd1_build --freeze analysis/wall_fusion_truth/pr6_hyd1_freeze.json --out <R>/hyd1_report` |

## Files

| path | what it is |
|---|---|
| `negcheck/` | pre-check tables: `summary.txt`, `negatives_family_table.tsv` (header only: 0 family members among the 161 shared negatives), `negatives_panel_table.tsv` (B6QUR1, HsbA), `related.tsv` (7 curated truth proteins at E <= 1e-3), `proteome_panel_ga.tsv` (19 rows, 18 distinct proteome proteins), `proteome_labelled_phmmer.tsv` (header only), `weak_E1e-3.tsv`, `proposed_negatives_ignored.json`, `provenance.txt` |
| `extract/` | `members/` (`members.tsv` 30 rows, `excluded.tsv` 3, `units.tsv`, `member_ids.txt`, `summary.txt`), `units_cluster.tsv`, `clu_cluster.tsv`, `ext_cluster.tsv`, `cysteines.tsv` (descriptive), `provenance.txt`, `hyd1_members.faa.sha256` |
| `report/` | `recall_by_unit.tsv` (9 unit/mode rows, all columns including the exact 95 percent intervals and the panel_ignored counted columns), `notes.txt` |
| `units/<unit>/<mode>/` | one folder per unit and mode (`cluster30` headline, `none` descriptive); the `all` model is `units/all/none/`; Sordariomycetes under `cluster30` is not buildable and holds only `manifest.json` |
| `.../cutoff.json` | frozen-rule cutoff of the custom model (with the freeze sha256, `positives_used: false`, the cutoff if the panel_ignored ids were counted) |
| `.../manifest.json`, `train_ids.txt`, `heldout_ids.txt` | held-out ids, training ids with clade and cluster, excluded 30 percent neighbours |
| `.../hits_custom.tsv`, `hits_pfam_ga.tsv`, `hits_relaxed_panel.tsv` | held-out member hits (custom and relaxed hydrophobin panel with scores; Pfam GA ids only) |
| `.../hits_phmmer_seed.tsv`, `hits_jackhmmer_seed.tsv` | seed-only phmmer and jackhmmer (training-alignment seeded) comparator hits |
| `.../neg_hits.tsv`, `proteome_hits.tsv`, `provenance.txt` | shared-negative and kept-proteome hits, per-task provenance (the `hmmbuild.log` files are not committed: the repository ignores `*.log`; they stay on HPCC with the alignments and models) |

All files hold identifiers, scores and labels only. No sequence is committed (a test checks it). The whole record is about 0.5 MB, so no hits table is gzipped or left out.

## Left on HPCC (not committed)

Per unit and mode `train.aln.fa` and `model.hmm` (the `all` model.hmm sha256 is in `units/all/none/provenance.txt`), and the FASTA files of the extract job (`hyd1_members.faa`, sha256 `34e11fc855b4576461ff4b34b24ef061a061e13398875cc8811af47e4ad450d4`, and `external.faa`). Nothing else stays only on HPCC.

Local copy of the fetched results (outside any repo): `/Users/jstajich/projects/adhesionPred-ecm33-results/`.
