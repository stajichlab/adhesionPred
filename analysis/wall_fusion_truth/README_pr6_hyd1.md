# PR 6 for hyd1_like: build and leave-cluster-out (DRAFT, nothing run on HPCC)

Design, decisions needed and the staged plan: `docs/superpowers/specs/2026-10-23-pr6-hyd1-like-design.md`. Target: KS0017, A0ACE4KFX7 (H. capsulatum, adhesion slice, Eurotiomycetes). Custom models are report-only; every cell is `descriptive` (the Kiss-unit row is labelled `partial`); no call, status, tier, `narrow_ok`, `family_measure.tsv` or `family_table.tsv` row is written. The ecm33 files are not edited; the ecm33 code is imported by path.

## Files

| File | Purpose |
|---|---|
| `pr6_hyd1_members.py` | stdlib; members (query + drift-supported hits, query coverage >= 0.5, named clade) and clade units from the committed drift record |
| `pr6_hyd1_build.py` | `cluster-units`, `manifest` (clade, cluster or `all`), `tasks`, `check-alignment`, `cutoff`; reuses `pr6_ecm33_build.py` for the independence helpers, alignment checks, cutoff rule and freeze guard (only the unit kinds and the two training-rule readings are new) |
| `pr6_hyd1_report.py` | `recall_by_unit.tsv` with independent comparators, `u_indep` and `u_with_ceiling`, Clopper-Pearson intervals in every row (function imported from `pr6_ecm33_eval_report.py`), the external line, the panel-counted numbers; refuses an unfrozen freeze or cutoffs made with another one |
| `pr6_hyd1_negatives_precheck.py` / `.sbatch` | pre-freeze check of the 161 shared negatives and the four proteomes (family rule and hydrophobin Pfam panel) plus a relatedness table |
| `pr6_hyd1_extract.sbatch` | members, UniRef50 sequences, MMseqs2 30 percent clusters, cluster units |
| `pr6_hyd1_build.sbatch` | the array: manifest, MAFFT, `hmmbuild`, cutoff, then the held-out searches |
| `pr6_hyd1_freeze_PROPOSED.json` | the freeze proposal (`frozen: false`); the owner copies it to `pr6_hyd1_freeze.json` |
| `tests/wall_fusion_truth/test_pr6_hyd1.py` | fixtures and the committed drift record |

Reuse note: no ecm33 file was refactored. `pr6_ecm33_report.py` is imported for `read_tsv`, `write_tsv`, `sha256_file`, `recovered_by_score`; `pr6_ecm33_negatives_precheck.py` for the tblout parser and `build`. What was copied or rewritten: `unit_row` (new comparator logic), `make_manifest` (units that are clades or clusters), `unit_table`/`select_members` (different input table), and the three sbatch scripts (adapted, trimmed).

## Run order (exact commands, NOT run; each HPCC step needs the owner's go-ahead)

Owner decisions of 2026-10-23 are recorded in `pr6_hyd1_freeze_PROPOSED.json` (`decided_by_owner_2026-10-23`), spec section 9. The pre-check job (30321490) is done.

Local check first (no HPCC):
```
python3 -I analysis/wall_fusion_truth/pr6_hyd1_members.py \
  --drift-hits docs/reports/data/pr6_drift_check_2026-10-22/drift_hits.tsv \
  --hits-manifest docs/reports/data/pr6_drift_check_2026-10-22/hits_manifest.tsv --out /path/new_dir
```
Expected: 30 members (query + 29) in 4 clades; 3 excluded; units Sordariomycetes `unit`, Dothideomycetes `unit`, Leotiomycetes `too_few_members`, Eurotiomycetes `kiss_single` (reading `total`).

1. **Sign the freeze** (the owner only): `cp pr6_hyd1_freeze_PROPOSED.json pr6_hyd1_freeze.json`, set `"frozen": true` and `owner_signature`, commit. Nothing is built before it. `panel_ignored` (19 ids), `negatives_ignored` ([]) and `heldout_truth.external_set` (6 ids) are already filled.
2. **Update the HPCC clone** to the merge commit of the PR that holds the freeze and these files (go-ahead needed).
3. **Extract job** (go-ahead needed). It reads the signed freeze for the training-rule reading and the external set:
   ```
   sbatch --export=ALL,PROJ_ROOT=...,CODE_COMMIT=<sha>,OUT=$W/hyd1_extract \
     -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log .../pr6_hyd1_extract.sbatch
   ```
   Expected: `members/` (members.tsv 30 rows, excluded.tsv 3, units.tsv, member_ids.txt 29 ids), `hyd1_members.faa` (30 sequences) and `.sha256`, `clu_cluster.tsv`, `units_cluster.tsv`, `cysteines.tsv` (descriptive), `external.faa` (6 sequences), `ext_cluster.tsv` (clusters over members + external), `provenance.txt`. Checks: 29 ids found in UniRef50; no duplicate id; no external id is a member.
4. **Build array** (go-ahead needed): compute N from `pr6_hyd1_build.py tasks` (see the header of `pr6_hyd1_build.sbatch`; the list now also holds one `external` task), then `sbatch --array=0-$N%4 --export=ALL,PROJ_ROOT=...,CODE_COMMIT=<sha>,EXTRACT=$W/hyd1_extract,OUT=$W/hyd1_build,PFAM_RELEASE=38.2 ... pr6_hyd1_build.sbatch`. Results in `<out>/<unit>/<mode>/`. The job stops without a frozen freeze and when the extract job used another reading than the freeze.
5. **Report** (local, no HPCC): `python3 -I analysis/wall_fusion_truth/pr6_hyd1_report.py --results <hyd1_build> --freeze analysis/wall_fusion_truth/pr6_hyd1_freeze.json --out <dir>`. Every row prints the exact 95 percent Clopper-Pearson interval for each recall and for recovered-of-u (`*_ci95`); a wide interval is not a result. The `line` column marks headline, descriptive, external and Kiss-unit rows; the extra calls and hard-negative rate are printed with the `panel_ignored` ids left out and with them counted.

## Checks built in

Held-out ids absent from every manifest and alignment (also the excluded 30 percent neighbours); no manifest, task list or cutoff without a frozen freeze; the cutoff function takes no positive score; the held-out searches start only after `cutoff.json` exists; seed comparators are skipped in the Kiss unit; the report refuses cutoffs made with a different freeze.

## Unconfirmed

HPCC resources and run times (guesses); `hmmfetch -f` by model name for the eight panel models (names from `data/sorting_hat/family_table.tsv`; HsbA name unverified against Pfam-A 38.2); `phmmer`/`jackhmmer` `-Z` behavior; MMseqs2 option set; that `Cimm_RS` etc. are in `$PROTEOME_DIR` (the ecm33 default path is reused); the inside-clade Eurotiomycetes hits stay out (no script drafted for their drift check); the extract/build steps for the external line (`external.faa`, `ext_cluster.tsv`, the `external` task) are new and untested on HPCC. Nothing in this folder has been run on HPCC.

## Result of the build (2026-10-23)

Pre-check 30321490, extract 30325384 and build array 30335767 (10 tasks, code `9f4843b`, freeze sha256 `65ed9bb6...552ff2`) all ran; the report ran locally. No unit is testable (u_indep at most 2 against `min_u` 3; Sordariomycetes under cluster30 is not buildable, 4 training sequences). The relaxed hydrophobin comparator already recovers most held-out members (for example 19 of 19 Sordariomycetes under `none`), the Kiss unit (n = 1) is recovered, and every cell is descriptive. Record: `docs/reports/data/pr6_hyd1_like_2026-10-23/`; report: `docs/reports/2026-10-23-pr6-hyd1-like-build.md`.
