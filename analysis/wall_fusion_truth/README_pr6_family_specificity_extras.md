# PR 6: the second specificity run, three extra families (DRAFT, not run on HPCC)

Sibling of `README_pr6_family_specificity.md` (the first run, 30 families by the rule "coverage >= 0.5 or E <= 1e-30"). The first run's script,
sbatch, outputs and tests are unchanged; this run imports its code (`pr6_family_specificity.py`) and does not copy it.

**Owner decision of 2026-10-24:** pursue three families that the Pfam rule left out because their Pfam hits were weak partial-domain hits, each
with a PANTHER family of 1,000+ fungal proteins. **This work changes no tier, no panel, no status, no `narrow_ok`, and writes neither
`data/sorting_hat/family_measure.tsv` nor `family_table.tsv`.** `convert` has no `--out`: it is a dry run only.

## The list (`pr6_family_specificity_extras_list.tsv`, committed, in this order)

| Order | Family | Pfam (name for `hmmfetch -f`) | Scan hit (cluster, E-value, coverage of the query) | PANTHER (descriptive) | Own members | Expected |
|---|---|---|---|---|---|---|
| 1 | Big1 (`big1`) | PF27678 `BIG1_C` | KS0087, 8.2e-30, 0.170 | PTHR28285, 1,349 fungal, scoping call "specific" | 1 (P38813) | specific |
| 2 | KNR4/SMI1 (`cot-2/gs-1`, `smi1`) | PF09346 `SMI1_KNR4` | KS0006, 5.2e-19, 0.323 | PTHR47432, 1,432 fungal, "specific" | 2 (A0A1D8PE53, P38678) | specific |
| 3 | DFG10 (`dfg10`) | PF02544 `Steroid_dh` | KS0090, 7.4e-09, 0.419 | PTHR14624, 1,060 fungal, "intermediate" | 1 (P40526) | broad |

Why these three: each is the Pfam handle of a wall cluster that the first run could not reach (BIG1_C covers 17 percent of the query; the three
hits miss both limits of the rule), and each sits in a large fungal PANTHER family, so a working specific model would cover many proteins. The
order is the owner's. The names and E-values are read from `docs/reports/data/pr6_pfam_scan_2026-10-22/scan_per_cluster.tsv`; `read_extras` and
`check_against_scan` stop if the list, its order, a name, an accession, a cluster, an E-value or a coverage differs from the committed scan, or if a hit
meets the Pfam rule (then it belongs to the first run). The PANTHER columns come from `docs/reports/data/panther_family_scoping/family_summary.tsv`
and are metadata (`extras_metadata.tsv`): they are not a filter and not truth.

## The test (the same as the first run)

`hmmsearch --cut_ga` of the three models over the Kiss positives (132), the 161 shared hard negatives (137 N1, 24 N2; plain accession) and the four
kept proteomes (Afum_Af293_UniProt 9,647, Calb_SC5314 6,212, Cimm_RS 9,910, Scer_S288C 6,722; 32,491 proteins). The family's OWN evidenced members
are all members of the listed cluster (the first run's rule): **4 sequences in all (Big1 1, KNR4/SMI1 2, DFG10 1)**, and the committed
`family_members.tsv` shows it. The own-member recall is circular (`leakage=tuned_on_truth`); with 1 or 2 members it says almost nothing. This is why the
tier derived by `cellsurface_sorting_hat.family_tier` is `smoke` at most (`descriptive_shared_sets`). A0A1D8PE53 is a member of KS0006 although the scan's hit
sequence of the cluster is P38678; whether it hits `SMI1_KNR4` is a result of the run, not an assumption.

## Files

| File | Purpose |
|---|---|
| `pr6_family_specificity_extras_list.tsv` | the owner-approved list (accession, name, cluster, gene, PANTHER family, scan E-value and coverage, expected, note) |
| `pr6_family_specificity_extras.py` | `candidates` (list + scan -> the four candidate files of the first run's format + `extras_metadata.tsv`), `check-models`, `summarize`, `convert` (dry run) |
| `pr6_family_specificity_extras.sbatch` | HPCC job, DRAFT, NOT SUBMITTED: pin check, inputs check, `hmmfetch -f` of 3 names, `check-models`, `hmmsearch --cut_ga` over six FASTA files, `summarize` |
| `docs/reports/data/pr6_family_extras_candidates_2026-10-24/` | the committed candidate tables; the job makes them again and stops if one differs |

## Run order (all commands NOT run; each HPCC step waits for the owner's explicit go-ahead)

Sign-off: the owner approved the extras list by the decision of 2026-10-24. The clone update and the job need the owner's go-ahead (D21).

0. Locally, offline: `python3 -I analysis/wall_fusion_truth/pr6_family_specificity_extras.py candidates --out /new/empty/dir` and compare with the
   committed folder (the tests do this).
1. HPCC, update the clone (a write): `git -C $B/adhesionPred-kiss53 pull --ff-only origin main`; it must print the merge commit of this PR; the tree must
   be clean.
2. Submit (a write; the inputs are those of the first run):
   ```
   B=/bigdata/stajichlab/jstajich/projects
   W=$B/adhesionPred/_workdir/kiss_cluster_split
   sbatch --export=ALL,PROJ_ROOT=$B/adhesionPred-kiss53,CODE_COMMIT=<merge commit>,IN=$W/in,OUT=$W/family_specificity_extras_out,PFAM_RELEASE=38.2 \
     -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log \
     $B/adhesionPred-kiss53/analysis/wall_fusion_truth/pr6_family_specificity_extras.sbatch
   ```
   Partition `short`, 4 cores, 8 GB, 1 hour; writes only to a new `$OUT`. It stops with a message if a name (`BIG1_C`, `SMI1_KNR4`, `Steroid_dh`) is missing
   from Pfam-A 38.2, if a model has no GA line, if the release directory is not `*Pfam38.2`, or if `$OUT` exists.
3. Read the outputs (read only): `scp` the small tables back.
4. Record PR (results and a short report). Only with the owner's sign-off per family, and as a dry run first:
   `PYTHONPATH=src python3 analysis/wall_fusion_truth/pr6_family_specificity_extras.py convert --specificity <family_specificity.tsv> --source <report> --date <date>`.

## Expected outputs (`$OUT`)

The tables of the first run for three rows: `family_specificity.tsv`, `positives_hits.tsv`, `negatives_hits.tsv`, `proteome_hits.tsv`, `provenance.tsv` (it names
the list and its sha256 instead of the Pfam rule, and the tier cap), `provenance.txt`, six `domtbl_*.txt`. Accessions, counts and FASTA header text; no sequence.

## What a result would mean

- A quiet family (no N1 hits, a few proteome hits that read as the right protein) supports a specific model for Big1 or KNR4/SMI1 on this narrow test, but
  with 1 or 2 own members the tier is `smoke`, and a hit on a partial domain (BIG1_C covers 17 percent) says nothing about the full protein.
- `Steroid_dh` is a general steroid-dehydrogenase domain: **many hits on non-members (negatives and proteomes) are expected, and that is a valid
  finding**: it would show that the Pfam handle of DFG10 is not specific and that a narrower model is needed (DFG10 itself is in a PANTHER family that
  the scoping table calls intermediate).
- Hits on the other Kiss positives are descriptive cross-reaction counts, not errors. No result enters the panel without a `family_measure.tsv` row and
  the owner's sign-off in its own PR.

## Not confirmed

- HPCC was not touched: the proteome files, `hmmfetch -f` by name for `BIG1_C`, `SMI1_KNR4`, `Steroid_dh` in Pfam-A 38.2 (names come from the scan table; a
  missing name stops the job), and that `hmmfetch -f` exits non-zero on a missing name (the count and `check-models` checks back it up).
