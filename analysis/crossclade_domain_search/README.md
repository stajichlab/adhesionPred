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

## Round 2: real AlphaFold models, three structure methods (built 2026-11-26)

Design, declared before the run (revised after an independent review, before any round 2 run):

- **Structures:** AlphaFold DB models (`round2_inputs.py`), used only if the model sequence equals the round 1 sequence. Queries
  are cut to the round 1 crop (the crop's sequence must equal `queries.faa`). Target models lose residues with pLDDT < 50 (long
  low-confidence disordered regions give spurious structure similarity); fewer than 40 residues left = not used. Targets:
  (a) every round 1 target hit at E <= 0.1 by either method (`candidate`), (b) the PF21671 reference (`pf21671`), and
  (c) background: per species, proteins drawn in seeded random order from the **whole proteome** (minus queries and PF21671
  members; candidates can be drawn) until 500 models are usable, alternating into half A (null) and half B (held out).
- **Methods:** `prostt5` = round 1 Foldseek with predicted 3Di (not re-run); `af_fs` = Foldseek on the AFDB models;
  `reseek` = Reseek v3.01 `-verysensitive -stats superfamily`. TEA is left out (not installed); MMseqs2 is context only.
- **One scale:** every score is put on the round 1 E-value scale (expected chance hits among the 15,587 round 1 targets):
  `af_fs` E x 15,587 / number of target models searched; `reseek` P x 15,587; `prostt5` as is.
- **Null:** for `af_fs` and `reseek`, the best scaled score of the query against half A, leaving out the target itself; for
  `prostt5`, the best of its 51 sequence decoys. `nulls.tsv` gives how many half-A proteins each null rests on.
- **Support:** scaled score <= 1e-3 and better than the query's null. A target is **confirmed** when one query pair has support
  from 2 of the 3 methods. Background proteins are never scored as hits; `background_hits.tsv` lists those at <= 1e-3 for a look.
- **False-pass rate:** half B is scored by the same rule (null from half A); `summary.tsv` `half_B_passing` per method and for
  the 2-of-3 rule is the rate at which randomly drawn proteins pass.
- **Caveats, stated before the run:** `af_fs` and `reseek` read the same model and `af_fs` shares Foldseek scoring with
  `prostt5`, so the votes are not independent (half B measures the rule as a whole); every scored target was selected by a round 1
  hit; Reseek keeps a k-mer prefilter even in `-verysensitive` mode, so some pairs have no score (= no hit). AFDB has no model for
  proteins over 2,700 residues, so *C. auris* IFF4109 (2,946 aa) is expected to be unscorable by `af_fs`/`reseek`;
  unscorable queries and PF21671 members are counted in `summary.tsv`.
- **Gene models:** `gene_model_check.py` lists, for targets without a signal peptide (M7WU28, M7WQR1) and M7X8Q7, the UniProt
  record and the other Rhodotorula members of the same UniRef50 cluster (length, signal peptide, N terminus).

```
P=/bigdata/stajichlab/jstajich/projects/adhesionPred
W=$P/_workdir/crossclade_domain_search
C=$W/code   # detached worktree at the reviewed commit (git -C $P worktree add --detach $C <commit>, or checkout there)

# 5. AFDB models (login node, about 1,200 small downloads)
/usr/bin/python3.12 -I $C/analysis/crossclade_domain_search/round2_inputs.py --round1-inputs $W/in \
  --round1-hits $W/round1_summary/hits.tsv --background-per-species 250 --out $W/in2

# 6. Reseek binary, once (login node)
mkdir -p $W/tools && curl -fsSL -o $W/tools/reseek-v3.01 \
  https://github.com/rcedgar/reseek/releases/download/v3.01/reseek-v3.01-linux-x86 && chmod +x $W/tools/reseek-v3.01

# 7. search (exfab, CPU)
sbatch --export=ALL,PROJ_ROOT=$C,CODE_COMMIT=$(git -C $C rev-parse HEAD),IN2=$W/in2,RESEEK=$W/tools/reseek-v3.01,OUT=$W/round2 \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log $C/analysis/crossclade_domain_search/round2_search.sbatch

# 8. summary and gene-model check (seconds)
/usr/bin/python3.12 -I $C/analysis/crossclade_domain_search/round2_summarize.py --round1-inputs $W/in \
  --round1-summary $W/round1_summary --inputs2 $W/in2 --results2 $W/round2 --out $W/round2_summary
/usr/bin/python3.12 -I $C/analysis/crossclade_domain_search/gene_model_check.py --acc M7WU28 --acc M7WQR1 --acc M7X8Q7 \
  --out $W/gene_model_check.tsv
```

## Later rounds (not built)

- TEA as a fourth method; a repeat-aware null for the β-helix and Flo11 queries (HYR1, IFF4109, FLO11).
- More targets: other Microbotryomycetes, *Ustilago maydis* (Rep1/Hum3), *Trichosporon asahii*.
- Strategy C link: the Rhodotorula biofilm assay (`config/site.yaml` `rhodotorula_biofilm`, 23 strains) could test whether candidate-family
  copy number follows biofilm OD across strains with genomes.
