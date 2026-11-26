# Cross-clade domain search, round 1 results (2026-11-26)

Design, cut-offs and the rules for what counts as a result were declared before the run in
`analysis/crossclade_domain_search/README.md`; this note reports against them.

## Run

| Step | Where | Record |
|---|---|---|
| inputs | HPCC login node, `build_inputs.py` at 828c3b4, UniProt release 2026_04 | `inputs_provenance.txt`, `inputs_SHA256SUMS`, `query_domains.tsv`, `truth_pf21671.tsv` |
| search | SLURM job 32482071, `exfab` (gpu12, RTX 6000 Ada), 15 min 54 s; Foldseek 10-941cd33 (GPU), MMseqs2 17-b804f; ProstT5 f16 weights | `search_provenance.txt` |
| summary | `summarize.py` at 828c3b4 with `/usr/bin/python3.12` (the default `python3` on the HPCC is 3.9; the repo needs ≥ 3.11) | `summary.tsv`, `per_query.tsv`, `hits.tsv`, `decoy_hits.tsv` |

The code ran from a detached worktree at 828c3b4 (`_workdir/crossclade_domain_search/code`), because the main HPCC clone holds
another branch's work in progress. Targets: 15,587 proteins (*R. toruloides* NP11 8,160, *C. neoformans* H99 7,427). 14 queries,
714 decoys (51 per query). Raw search tables and the target 3Di FASTA stay on the HPCC in `_workdir/crossclade_domain_search/round1/`.

## Panel `basidio` (known-positive check)

| Method | E ≤ | NP11 PF21671 recovered | Other NP11 targets | Decoy hits (decoys with a hit / 357) |
|---|---|---|---|---|
| Foldseek | 1e-3 | 12/12 | 1 | 0 (0) |
| Foldseek | 1e-2 | 12/12 | 5 | 0 (0) |
| Foldseek | 0.1 | 12/12 | 8 | 2 (1) |
| MMseqs2 | 1e-3 | 12/12 | 1 | 12 (10) |
| MMseqs2 | 1e-2 | 12/12 | 2 | 27 (23) |
| MMseqs2 | 0.1 | 12/12 | 3 | 94 (79) |

- Both methods recover all 12 NP11 CPL1-like proteins at the primary cut-off, so sensitivity does not separate them here. The
  difference is the null: Foldseek has no decoy hit at E ≤ 1e-2 (best decoy E-value per query 0.07 to 0.57), MMseqs2 has decoy
  hits at E ≤ 1e-3 for 5 of 7 queries (best 4.5e-5 to 1.4e-2). The H99 PF21671 domains are short and Ser/Thr-biased, which a
  composition-keeping shuffle reproduces; the 3Di alphabet does not.
- NP11 targets beyond the Pfam reference (all hit by several independent H99 queries):
  - **M7WU28**: E 1.9e-5 (Foldseek) and 9.2e-5 (MMseqs2), hit by all 7 queries with Foldseek; **no signal peptide** in the flags
    table. Likely a CPL1-like protein that Pfam misses; check the gene model's N terminus.
  - Foldseek only, secreted, at E ≤ 1e-2 where its decoys give no hit: **M7WYW5** (2.2e-3; MMseqs2 7e-3 for CPL1), **M7XG46**
    (3.3e-3, 7 queries), **M7WND4** (6.1e-3, 5 queries), **M7WT56** (9.6e-3).
  - Weaker (0.01 < E ≤ 0.1), secreted: M7WHW1, M7XWK5, M7WS82.
  - These are round 2 candidates (Reseek on AFDB models, domain alignment), not results: the primary cut-off is 1e-3.

## Panel `asco` (discovery)

| Method | E ≤ | Targets hit | Decoy hits (decoys with a hit / 357) | Targets with SP |
|---|---|---|---|---|
| Foldseek | 1e-3 | 13 | 42 (4) | 1 |
| Foldseek | 1e-2 | 39 | 230 (33) | 4 |
| MMseqs2 | 1e-3 | 0 | 0 (0) | 0 |
| MMseqs2 | 0.1 | 1 | 20 (19) | 1 |

- The null fails for three queries: HYR1 and *C. auris* IFF4109 (Hyr/Iff β-helix) and FLO11 have decoys hitting at E 5e-6 to
  2.4e-4 (all 42 decoy hits at 1e-3: HYR1 34, FLO11 6, IFF4109 2). Their hits cannot be read; repetitive folds need a stricter null or Reseek.
- Declared rule (secreted target, E ≤ 1e-3, and below the query's best decoy E-value): **one target passes**, NP11 **M7X8Q7**
  for SAG1 (E 4.9e-4 against best decoy 1.05e-3, so a margin of only 2-fold). Hand-check list.
- Not secreted, so outside the rule but worth noting:
  - FLO5 (PA14) and EPA1 (GLEYA, a PA14-like fold) hit the same five H99 proteins (J9VNL1, J9VVK7, J9VP97, J9VZB3, J9VNK7),
    target coverage 0.13 to 0.16, i.e. a PA14-like domain inside larger intracellular proteins; agreement between two queries
    supports the fold, not an adhesin function.
  - SAG1 and ALS3 both hit NP11 M7WQR1 (2.5e-5 and 3.2e-4) and M7XGB3; M7WQR1 has a 494-residue disordered run and no signal
    peptide, so a gene-model check is worth doing.
- MMseqs2 finds nothing below 1e-3 in this panel, as expected for cross-phylum distances.

## Conclusions for round 2

1. For the CPL1-like family, ProstT5/Foldseek is as sensitive as sequence search and has a much cleaner null; it adds 4–5
   secreted NP11 candidates at E ≤ 1e-2 with zero decoy hits. Worth confirming with Reseek on AlphaFold models.
2. Ascomycete adhesin domains give one weak secreted candidate (M7X8Q7) in NP11; β-helix and Flo11 queries need a null that
   handles repeats before their hits mean anything.
3. M7WU28 and M7WQR1 lack a signal peptide in the flags; check their gene models before reading that as "not secreted".
