# PR 6 family-specificity results: run record of jobs 30423745 (30 candidate families) and 30426504 (3 extra families)

Small tables copied back from HPCC after the two family-specificity screens of the PR 6 candidate Pfam families. Report: `docs/reports/2026-10-24-pr6-family-specificity-results.md`. Design: `analysis/wall_fusion_truth/README_pr6_family_specificity.md` and `README_pr6_family_specificity_extras.md`. This record changes no tier, panel, `narrow_ok` or status, and `data/sorting_hat/family_measure.tsv` and `family_table.tsv` are untouched (a test compares them with `origin/main`).

## What the runs are tied to

| item | run1 | run2 |
|---|---|---|
| job | 30423745 (`pr6_family_specificity`) | 30426504 (`pr6_family_specificity_extras`) |
| code commit (pin ok in the log) | `e36c6e4a02031f30e1c64a74f581ee0f7389e8d7` | `f56277354f174f0a852f86cf7fc61234c84995e0` |
| script | `analysis/wall_fusion_truth/pr6_family_specificity.py` | `analysis/wall_fusion_truth/pr6_family_specificity_extras.py` |
| families | 30 (candidate rule: coverage >= 0.5 or E <= 1e-30) | 3, the owner list of 2026-10-24: BIG1_C PF27678, SMI1_KNR4 PF09346, Steroid_dh PF02544 |
| models | 30 fetched from Pfam 38.2, 30 with a GA line | 3 fetched, 3 with a GA line |
| Pfam | 38.2, `Pfam-A.hmm` sha256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e` | same |
| HMMER | 3.4 (`hmmsearch --cut_ga`) | same |
| inputs | 132 positives, 161 hard negatives (137 N1 + 24 N2), four kept proteomes, 32,491 proteins | same |
| leakage / tier | `tuned_on_truth`; tier capped at `smoke` (the script writes only report tables) | same |

Kept proteomes: Afum_Af293_UniProt 9,647; Calb_SC5314 6,212; Cimm_RS 9,910; Scer_S288C 6,722. Input FASTA hashes are in `run1/provenance.txt` (the same six hashes are in `run2/provenance.txt`).

## Files

| path | what it is |
|---|---|
| `run1/`, `run2/` `family_specificity.tsv` | one row per family: own-member recall (circular), other-positive hits, negative hits, hits per proteome and per 10,000 |
| `positives_hits.tsv`, `negatives_hits.tsv` | which Kiss positives (own or other) and which hard negatives each family hit |
| `proteome_hits.tsv` | every hit in the kept proteomes with its FASTA header text (an identifier and a name, no sequence) |
| `provenance.tsv`, `provenance.txt` | pin, Pfam release and hash, input hashes, proteome sizes, `tier_changes none`, `panel_changes none` |
| `hpcc/*.log.txt` | the two job logs (renamed from `.log`, which is git-ignored); `hpcc/hpcc_sha256.txt` is the sha256 list of the HPCC-side originals |
| `summary.tsv` | both runs in one table, made by `analysis/wall_fusion_truth/pr6_family_specificity_results_summary.py`; adds `own_in_proteomes`, `kept_hits_excl_own`, `kept_per10k_excl_own` |
| `SHA256SUMS.txt` | sha256 of every committed file except this README and itself (`cd` here and run `sha256sum -c SHA256SUMS.txt`) |

Meaning of the last columns of `family_specificity.tsv`: `kept_hits_total` is the number of proteins in the four kept proteomes that the family model hits at its GA cutoff; `kept_n_proteins` is 32,491; `kept_per10k` is `10000 * kept_hits_total / 32491`. The total does NOT exclude Kiss positives that happen to be in those proteomes. `summary.tsv` counts those separately (`own_in_proteomes`, 0 to 2 per family in run1, 0 in run2).

## sha256 of the HPCC originals

The committed run tables are byte-identical to the HPCC files (a test checks the list in `hpcc/hpcc_sha256.txt` against the committed copies). Paths are relative to `_workdir/kiss_cluster_split/`; run1 is `family_specificity/`, run2 is `family_specificity_extras_out/`.

| file | run1 sha256 | run2 sha256 |
|---|---|---|
| `family_specificity.tsv` | `ee1a2b9a8fffad18bab33db79aa79eb489a5df01d03c78b05e997303352f5090` | `c947d137002b4e1683d4d31fbe440b30d9df3e6ad8e006f6d65a11cd3390c58a` |
| `negatives_hits.tsv` | `3af8849c8fdb4be6258f853bcdd52a94a9de4e293fdf5629a3696aa1999555c3` | `15efdf7f02b19d1cc3c9bd51e58eaefecce17f6b7493b64f5de9ba39e6bba660` |
| `positives_hits.tsv` | `9af6361749587f486d7cd3c5e3ffa55e921dac88a4e1509ebf1c84d00fffaa3d` | `f24eaf0b31044c160446928cf4a7b604c1b55531b63f19ee37cd35a61a68faf4` |
| `proteome_hits.tsv` | `a78a689324739bab224a3a67ad6728f26203d731d86642fa1980731628c7af63` | `5e5826609bb51e24996568f8ce5fd37f3f2ff402cd6e5a78fd5e34152e08d383` |
| `provenance.tsv` | `73c81c79799cb96fd4935a1581b3028d13cac53e711d96a88d883b2f2c86640c` | `79badc8abbfd48a4f5111d278e48882bcae158d4f92be8ede120bd52b231a440` |
| job log | `9d166eb479bf57467bd76d5cfc018f6b4edf71fd3a40fecda80571e50f2ddbbb` | `a7c2138a47084b9d45d5c7199e0de0c38f140b21a970f014cd2224c44e00d2d1` |

## Left on HPCC (not committed)

The domain tables `domtbl_*.txt` (the raw `hmmsearch` output) stay on HPCC:
`/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/kiss_cluster_split/family_specificity` (run1) and `.../family_specificity_extras_out` (run2). `provenance.txt` was copied back without a recorded hash (it is not in the HPCC list); the committed hash is in `SHA256SUMS.txt`.

## How to reproduce

1. The tables: re-run the two sbatch scripts at the pinned commits (needs the owner's go-ahead per job): `analysis/wall_fusion_truth/pr6_family_specificity.sbatch` at `e36c6e4` and `pr6_family_specificity_extras.sbatch` at `f562773`.
2. `summary.tsv`: `python3 -I analysis/wall_fusion_truth/pr6_family_specificity_results_summary.py` (deterministic; a test regenerates it and compares).
3. Tests: `pytest tests/wall_fusion_truth/test_pr6_family_specificity_results.py`.

All files hold identifiers, scores, names and labels only. No sequence is committed. Folder and file names avoid the ignored patterns (`build/`, `results/`, `input/`, `output/`, `dist/`; `*.csv`, `*.log`).
