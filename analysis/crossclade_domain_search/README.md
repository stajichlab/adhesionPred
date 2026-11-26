# Cross-clade domain search, round 1 (strategy A)

Plan and background: `docs/reports/2026-11-25-plm-survey-rhodotorula-adhesion.md` (section 2A). Round 1 ran on 2026-11-26 (results: `docs/reports/data/crossclade_domain_search_round1_2026-11-26/`); every HPCC step
needs the owner's go-ahead.

## Question

Can structure-alphabet search (ProstT5-predicted 3Di + Foldseek) find adhesion-domain proteins in *Rhodotorula toruloides* NP11 that sequence
search and the current sorting-hat rules miss?

## Design (declared before the run)

**Targets:** *R. toruloides* NP11 (8,160 proteins) and *C. neoformans* H99 (7,427), taken from the staged proteomes of the non-yeast test
(`_workdir/nonyeast_test/in`), so the accessions match the committed surface flags
(`docs/reports/data/aga1_nonyeast_test_2026-10-26/all_proteins_flags.tsv`).

**Two query panels** (`queries.tsv`; each query is its Pfam domain, InterPro envelope ± 10 residues):

| Panel | Queries | What it tests | Reference |
|---|---|---|---|
| `asco` | 7 Saccharomycotina adhesin binding domains: Als N-terminal (ALS3, SAG1), PA14 (FLO5), Flo11 domain (FLO11), GLEYA (EPA1), Hyr/Hil β-helix (HYR1, *C. auris* IFF4109) | Discovery: do ascomycete adhesin folds have relatives in basidiomycete yeasts? | None exists (no basidiomycete homolog of these families is known), so this panel reports candidates and decoy rates only |
| `basidio` | the 7 CPL1-like (PF21671) proteins of *C. neoformans* H99, including CPL1 (Cfl1/Cpl1 family, PMID 28039134) | Sensitivity check across Tremellomycetes → Microbotryomycetes, and whether structure search adds CPL1-like proteins beyond Pfam | The 12 NP11 proteins with PF21671 in UniProt (fetched at build time into `truth_pf21671.tsv`) |

Why two panels: the Cfl1 family is not homologous to Als, Flo or Epa. A check "ascomycete queries recover Cfl1" would score 0 by design, so the
known-positive check is done within the basidiomycetes (panel `basidio`), and the cross-phylum question (panel `asco`) is reported as discovery.

**Methods:** (1) Foldseek 3Di+AA with ProstT5 3Di on GPU, exhaustive; (2) MMseqs2 sequence baseline with the prefilter off
(`--prefilter-mode 2`, all-vs-all). Masking is off in both (`--mask 0`) so the comparison is not a masking difference. Reseek (needs 3D models) and TEA (needs
installation) are round 2; the "2 of 3 methods agree" rule of the plan applies only once all three have run.

**Null:** 51 decoys per query (one reversal and 50 shuffles within 10-residue windows, which keep length and local composition, so
S/T-rich and repeat-rich domains keep their bias). All decoy hits count, including a decoy hitting its own forward protein. Reported:
decoy hits and decoys with a hit per cut-off, and the best decoy E-value per query (`per_query.tsv`). With 357 decoys per panel, an
empirical false-hit rate near 1 in 300 is resolvable; 7 reversals alone would not have been.

**Cut-offs:** E ≤ 1e-3 primary; 1e-2 and 0.1 reported beside it. Query hits to any query accession are dropped (self and paralog
queries), and PF21671 recovery counts only non-query members: all 7 H99 members are queries, so H99 recovery is 0/0 by design and only
the 12 NP11 members are scored.

**Surface context:** hits are joined with the committed per-protein flags (signal peptide, GPI, longest disordered run, broad/narrow flags).
For basidiomycetes the relevant context is a **signal peptide, not a GPI anchor**: all 19 CPL1-like proteins of H99 and NP11 have a signal
peptide and none has a GPI call, which is why the sorting hat flags none of them.

**What would count as a result:** panel `basidio`: Foldseek recovers more of the 12 NP11 PF21671 members than MMseqs2 at the same decoy rate, or
finds secreted non-Pfam proteins that a later check (Reseek, domain alignment) confirms. Panel `asco`: a secreted NP11 or H99 target hit at
E ≤ 1e-3 whose E-value is also below the best decoy E-value of that query goes to round 2 and to the hand-check list.

## Files

| File | Where it runs | Role |
|---|---|---|
| `queries.tsv` | committed | the 14 queries, panel, Pfam domain |
| `build_inputs.py` | HPCC login node (about 30 small HTTPS requests) | crops the query domains, writes 51 decoys per query (reversal + window shuffles), extracts the target proteomes, fetches the PF21671 reference; output outside the repo |
| `search.sbatch` | SLURM, `exfab` (1 Ada 6000) | ProstT5 3Di, Foldseek and MMseqs2 searches, provenance |
| `summarize.py` | HPCC (seconds) | `hits.tsv`, `summary.tsv`, `per_query.tsv` (no sequences); copied to `docs/reports/data/` afterwards |
| `tests/crossclade_domain_search/` | anywhere | fixture tests (offline) |

## Commands (one go-ahead per step)

Preconditions: the commit holding this folder is on the HPCC clone, which is clean for `analysis/crossclade_domain_search`.

```
P=/bigdata/stajichlab/jstajich/projects/adhesionPred
W=$P/_workdir/crossclade_domain_search
cd $P && git fetch && git checkout <branch-or-merge-commit>
mkdir -p $W/logs

# 1. inputs (login node, tiny)
python3 -I analysis/crossclade_domain_search/build_inputs.py --queries analysis/crossclade_domain_search/queries.tsv \
  --staged $P/_workdir/nonyeast_test/in --species Rtor_NP11=1130832 --species Cneo_H99=235443 --out $W/in

# 2. ProstT5 weights, once (login node, download only)
module load foldseek/10-941cd33 && foldseek databases ProstT5 $W/prostt5_weights $W/prostt5_tmp && rm -rf $W/prostt5_tmp

# 3. search (exfab GPU)
sbatch --export=ALL,PROJ_ROOT=$P,CODE_COMMIT=$(git -C $P rev-parse HEAD),IN=$W/in,WEIGHTS=$W/prostt5_weights,OUT=$W/round1 \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log analysis/crossclade_domain_search/search.sbatch

# 4. summary (seconds; needs Python >= 3.10, the HPCC default python3 is 3.9)
/usr/bin/python3.12 -I analysis/crossclade_domain_search/summarize.py --inputs $W/in --results $W/round1 \
  --flags docs/reports/data/aga1_nonyeast_test_2026-10-26/all_proteins_flags.tsv --out $W/round1_summary
```

## Later rounds (not built)

- Round 2: Reseek on AFDB models of the hits (NP11 is in AFDB, UP000016926) and TEA; apply the 2-of-3 rule; add the reversed-decoy rate per method.
- More targets: other Microbotryomycetes, *Ustilago maydis* (Rep1/Hum3), *Trichosporon asahii*.
- Strategy C link: the Rhodotorula biofilm assay (`config/site.yaml` `rhodotorula_biofilm`, 23 strains) could test whether candidate-family
  copy number follows biofilm OD across strains with genomes.
