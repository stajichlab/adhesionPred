# Hydrophobin validation: implementation plan

*2026-10-08. Revision 2 (after plan review 1: 1 blocker, 8 major, 6 minor, all addressed). Spec: `docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md` (revision 3 plus the
owner answers D6, D7 and D8 of 2026-10-08). Reviews of the spec: `...-design-review-1.md`, `...-review-2.md`.
Author: Claude Code (claude-sonnet-5-5). The plan review is `2026-10-08-hydrophobin-validation-review-1.md`.*

## Global constraints

1. Python is `/usr/bin/python3.12` for Mycelium scripts. Project code and tests use the project environment
   (`PYTHONPATH=$PWD/src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat tests/calibration_truth tests/hydrophobin_truth -q -p no:cacheprovider`).
   Analysis tests live in `tests/<analysis dir>/` and load the script with `importlib`, as `tests/calibration_truth` does.
   There is no `tests/analysis/` directory. The baseline on this branch before any work is 701 passed, 7 skipped.
2. Compute runs on SLURM (`short` or `short_gpu` partition), not on the login node, when a step takes more than
   a few minutes. Use `$SCRATCH` (`${SCRATCH:?}`) for transient files and copy results to `/bigdata`. Never use
   `BASH_SOURCE[0]` in scripts that run on SLURM.
3. Large text outputs are compressed (`.gz` for portable files, `.zst` for internal intermediates). New
   `bin/` or `analysis/` Python scripts that read large text accept `.gz`.
3a. Job sizing: group work into jobs of about 1 to 1.5 hours. Earlier per-proteome jobs ran 0.3 to 10 minutes, so
   each phase submits one job per step over all proteomes (a loop inside the job), not one job per proteome.
4. Commits: ruff and ruff-format run as hooks. After a failed commit, fix, re-add and re-commit, and check
   `git log`. Commit messages end with the attribution line given for this session.
5. No push, PR or merge without the owner. Work stays on local branches. This includes H0: the note fix is committed
   on the PR #75 branch locally, and the owner decides when it is pushed.
6. No result is reported unless a command produced it. Unknowns are listed as unknown.
7. Tests are written before the code they test (red, then green), as in the earlier plans.

## Branch and dependency map

Facts checked in review 1: `per-call-status` branched from `9a3202c` and does not contain PR #75's commits. PRs #75 and #76
are open. `per-call-status` merges into `hydrophobin-validation` with no conflict (`git merge-tree`). `modules/` is
identical on both branches.

| Phase | Tasks | Depends on | Branch |
|---|---|---|---|
| A (now) | H0, H3, H3b, H2, H2b, H2c, H1c | nothing merged (H1c touches only `modules/`) | H0 on `signoff-cfem-hydrophobin`; others on `hydrophobin-validation` |
| B | H4 (after H3 and H2); H1 | H4 touches only `modules/`. H1 touches `categories.yaml`, `outputs.py` and golden files that differ on `per-call-status`, so H1 waits until PR #76 and #75 are merged, or until `per-call-status` is merged into this branch locally | `hydrophobin-validation` |
| C | H5, H5b, H6, H1b, H7 | phases A and B | `hydrophobin-validation` |

If the owner prefers not to wait, the work branch can merge `per-call-status` locally (conflict-free). The H1 commit
message states which route was used.

## Phase A

### H0. Correct the PF01185 specificity note (PR #75 branch)

Files: `data/sorting_hat/family_table.tsv`.
Steps:
1. Count hits per proteome for PF01185 from `_workdir/sorting_hat/calibration/hydrophobin_discovery/*.domtbl`
   (`awk '!/^#/ && $4=="Hydrophobin"{print $1}' FILE | sort -u | wc -l` per file). Expected from the review: 6, 5, 6, 1 and 2 in Af293, A1163, W72310,
   RS and ER3.
2. Edit the note to the counted values.
3. Run the family-table test. Commit on `signoff-cfem-hydrophobin`. Acceptance: tests pass, the note matches the
   counts.

### H3. Read the sources and freeze the pattern

Files: `docs/paper/05-literature-verification.md` (new section 6), `data/sorting_hat/cys8_spacing.yaml` (new).
Steps:
1. Read Wessels 1994, Linder 2005 and Sunde 2008 for the class I and class II cysteine spacing. Use PubMed
   (`get_full_text_article` where open) and the project's own paper folder if the PDFs are there. If a source cannot be read in full, record
   "not read" and do not use its numbers.
2. For each spacing number record paper, page or figure, and the quoted sentence.
3. Read the full text of Kubicek 2008 (PMC2253510) and Yang 2006 (PMC1780129). Record any spacing, any
   pattern, and what Yang 2006 compared with.
   Source PMIDs are looked up first and written into the file (Wessels 1994, Linder 2005, Sunde 2008, Kubicek 2008 PMID 18186925,
   Seidl-Seiboth 2011 PMID 21424760, Yang 2006 PMID 17217508, Peñas 1998 PMID 9758836). If no source with a spacing can be read in full, **stop and
   report to the owner**. Do not write a spacing from memory.
4. Write `cys8_spacing.yaml`: for class I, class II and `other`, the min and max of each gap, each with a source
   field. Add a sentence on which spacing the module uses first.
5. Commit with message "freeze literature spacing". **No truth protein is used in this task.** Acceptance (checkable): `grep -c "source:" data/sorting_hat/cys8_spacing.yaml` equals the number of spacing entries;
   `git log --format=%h -1 -- data/sorting_hat/cys8_spacing.yaml` returns a hash that is an ancestor of the first H2 commit
   (`git merge-base --is-ancestor`). The hash is recorded in H2's provenance file.

### H3b. Wider prior-art search

Files: `docs/paper/05` (section 5 extended).
Steps: search PubMed for pattern-based or HMM-based hydrophobin identification and surface-active protein
prediction; search bioRxiv; read the InterPro entries for the Pfam models and the Pfam release notes. Also read Seidl-Seiboth 2011 (find a route; the PubMed converter gave no PMCID) and Peñas 1998 (PMC106595, open). Record each source with PMID or URL, what it says, and its bearing on the rescue rule. Acceptance: the novelty note in `docs/paper/04` is updated with what was found, or says "no further prior art found in these sources" with the sources named; the count of "PMID" or "doi" strings in the new section of `docs/paper/05` is at least the number of sources listed.

### H2. Truth set (after H3)

Files: `analysis/hydrophobin_truth/` (`build_truth.py`, `make_split.py`, tables, `provenance.json`, `README.md`).
Steps:
1. **Test first.** `tests/hydrophobin_truth/test_hydrophobin_truth.py` (new, loads `analysis/hydrophobin_truth/build_truth.py` with `importlib`): tier assignment on a small fixture of UniProt JSON
   entries (the T2 rule applies only to entries whose protein name or function text mentions hydrophobin or rodlet; the 15 non-matching query rows, for example BRLA_EMENI, are never positives; an entry with ECO:0000269 on function is T2, an entry with only ECO:0000255 is T3, an entry with a
   PMID-backed manual row is T1); a fixture where a protein is named hydrophobin by rule only is T3; the split
   function puts every member of a cluster in one part and is the same for the same seed.
2. `build_truth.py`: fetch UniProt JSON for the accessions in `sp_hydrophobin_query.tsv` (with evidence codes),
   assign tiers per annotation, add T1 entries from PubMed searches (each with a PMID in `manual_entries.tsv`),
   restrict T1 phenotypes to hydrophobin-specific ones, and write `truth_all.tsv` with columns: accession, species, tier,
   pmid, evidence_code, sequence_sha, source.
3. Cluster with MMseqs2 (30% identity, coverage 0.5, as in `analysis/calibration_truth/repeat_call_truth/`).
   Write `clusters.tsv`.
4. `make_split.py`: 30% development and 70% test by cluster, seed recorded, written once to `split.tsv`. The file is
   committed before any measurement is run.
5. Write `seen.tsv`: the 9 no-model entries, RodA to RodG, PSH_FLAVE, CFTH1, FBH1, the three Af293 GPI proteins
   (Q4WUX0, Q4WBD4, Q4WBC1), the proteins named in review 2 (B0YAH7, B0Y4E8, KAK9640788.1, F00FD2C2_006465-T1,
   F00FD2C2_006466-T1), and every protein of the 20 screening candidates from `analyse.py` output in the in-scope proteomes
   (regenerate the list and store it).
5a. Write `curated_entries.tsv` for hydrophobins already in `data/curated/` (grep the curated tables; report the count, even if 0).
5b. Write `hard_negatives_swissprot.tsv` (cross-species hard negatives of spec 4.3 item 2: CFEM, cerato-platanin, killer toxins, small
   cysteine-rich secreted proteins, PIR and Ccw12-like, ATP8/ATP9) for the sequence-level report table only.
6. Split use: the development part is the only part on which the pattern or a threshold may change. The status files are
   written from the **test** part. If the test part has fewer than 10 positive clusters, the report says the split gives no
   usable test and the status is `smoke` with leakage `partial` (spec 4.4).
7. Print counts: T1, T2, T3 entries; clusters per tier; per species. **Report the counts as measured.** If the T1+T2 clusters are fewer than 20 per
   group, say so in the README.
Tests are in `tests/hydrophobin_truth/test_hydrophobin_truth.py`.
Acceptance: tests pass; counts printed and written to the README; the split file exists and is committed; no
rescue-rule result has been computed.

### H2b. Map truth proteins to the proteomes

Files: `analysis/hydrophobin_truth/map_to_proteomes.py`, `proteome_map.tsv`.
Proteomes: the in-scope eight in `_workdir/sorting_hat/*.faa` and the E1 files in
`/bigdata/stajichlab/shared/projects/Fungi_5k/input/` (*B. bassiana* ARSEF_2860, *F. velutipes* 6-3, *F. fulva* Race5_Kim, *F. graminearum* PH-1,
*P. expansum* MD-8, *P. ostreatus* PC9, *T. asperellum* FT101, *T. virens* Gv29-8).
Steps:
1. Test first: a fixture proteome and two truth proteins, one exact match, one at 96% identity and 95% coverage (mapped), one at
   80% (not mapped).
2. Map each truth protein to its proteome by exact sequence; otherwise by BLAST (ncbi-blast 2.14.0+, at least 95% identity
   over at least 90% of the length, best hit, one protein per hit). Record the identity.
3. Write `proteome_map.tsv` and a list of truth proteins with no match, per species.
4. Build per-species negative populations as in the repeat work: curated proteins and reviewed, secreted UniProt proteins of the
   species mapped by sequence, minus any protein with a hydrophobin label of any tier. Write `negatives_<species>.tsv`, labelled "assumed".
5. Proteome IDs: Fungi_5k files use `F…-T1` style IDs, so mapping is by sequence, never by name.
Acceptance: tests pass; the unmatched list is stored; every mapping has an identity value.

### H2c. Inputs for `calibrate truth` (new, from review 1 blocker)

Files: `analysis/hydrophobin_truth/make_calibration_inputs.py`, `calibration/<species>/truth.tsv`, `taxa.tsv`.
`calibrate truth` needs a `--truth` table with columns `id, label, cluster` where `id` is the proteome protein ID, and a taxon ID
per proteome. Read the `truth` subcommand in `calibration/cli.py` and the repeat-call inputs under
`analysis/calibration_truth/repeat_call_truth/` first, and copy their column names and rules.
Steps:
1. Test first: a fixture proteome with two positives and three negatives gives a truth table with the three columns, a cluster for every row,
   and no id outside the proteome.
2. Cluster all rows of one species (positives **and** negatives) with MMseqs2 at 30% identity and 0.5 coverage, so every negative has a cluster.
   Do this for the test part (status table) and for all rows (full table).
3. Write one `truth.tsv` per species (label values as the CLI expects) and `taxa.tsv` with the NCBI taxon ID of each proteome. Resolve
   taxon IDs from the NCBI taxonomy dump in `/srv/projects/db/taxonomy`, not from memory. Record any proteome whose strain has no node.
Acceptance: the CLI's own input check accepts each `truth.tsv`; tests pass.

## Phase B

### H1. Separate Pfam module and `hydrophobin_domain`

Files: `src/cellsurface_sorting_hat/modules/pfam.py`, `modules/cli.py`, `categories.yaml`, `outputs.py`,
`data/sorting_hat/family_table.tsv`, tests `test_data_and_scripts.py`, `test_cli.py`, golden files in
`tests/cellsurface_sorting_hat/golden/`.
Steps:
1. Tests first (red):
   - the module set in `test_data_and_scripts.py` includes `pfam_hydrophobin` and `pfam_hsba`;
   - a Pfam fixture with a hydrophobin-model hit gives `pfam_hydrophobin.hit = 1` and `pfam_adhesion.hit = 0` for the protein;
   - `wall_family_domain` is `not_called` for that protein; `hydrophobin_domain` is `called`;
   - a surface protein (R0 called) with only a hydrophobin hit is not in `other_surface_no_mechanism`;   - the report text lists the new calls and shows the category 2c line as `hydrophobin_domain` OR `hsba_domain` (display only);
   - the toy fixture in `test_cli.py` has `pfam_hydrophobin` and `pfam_hsba` tables, so the `other_basis` assertions hold.
2. Move rows. PF01185, PF06766, PF22354, PF28987, PF29785, PF29802, PF29465 to module `pfam_hydrophobin`; PF12296 to `pfam_hsba`.
   Add `pfam_hydrophobin` and `pfam_hsba` to `MODULES`. Keep `active` values as they are.
3. `categories.yaml`: add `hydrophobin_domain` and `hsba_domain`, both placed before the `other_*` calls; add both to the mechanism
   lists of both `other_*` calls.
4. Update `outputs.py` report text and the golden files. Regenerate the goldens only after reading each diff and
   confirming that the only changes are the new calls and the changed mechanism lists.
5. Run the full test suite. Commit. Acceptance: all tests pass; the diff of the goldens contains only intended changes.

Caveat recorded in the commit message: `config_sha256` changes, so existing call files are stale until H1b.

### H1c. Per-module family digest (phase A: needs only `modules/`)

Files: `modules/cli.py` (`_family_digest`), tests.
Steps: test first: changing a row of module A does not change the artefact digest of module B; changing a row of module A changes the digest of A. Then
change `_family_digest` to hash only the rows of the module. Check the effect on existing status files: module
identities of `pfam_adhesion` change once; record this in the commit message. Acceptance: tests pass.

### H4. `cys8_pattern` module (after H3, H2)

Files: `src/cellsurface_sorting_hat/modules/cys8.py` (new), `modules/cli.py`, `data/sorting_hat/cys8_spacing.yaml`, tests.
Steps:
1. Tests first, with **synthetic** sequences built from the frozen spacing (one class I-like, one class II-like, one with a
   gap out of range, one with seven cysteines) and a protein with no matching pattern; plus R0 `called` and `not_called` and a
   missing R0 row. Expectations: `hit=1` only when pattern and R0 `called`; missing R0 gives empty `hit` and state `ok`;
   invalid sequence gives the `na_invalid` row from `invalid_row` (as other wrappers), so the module state stays `ok`. **The 9 entries and any seen protein are not used in tests.**
2. Implement the wrapper. It reads the R0 module table through `_condition_table`, writes the columns `hit`, `spacing_class`, `n_cys`,
   `length`, and puts the pattern, the spacing source (file hash), and the `conditions` entry (R0 module and identity) into
   `params`.
3. Test that a changed spacing file changes `params_hash`, and that a new R0 identity changes it.
4. Register the module and add it to the module CLI. Do **not** edit `categories.yaml` here.
Acceptance: tests pass; no change to `categories.yaml`.

## Phase C

### H5. Runs on the in-scope proteomes

Steps:
1. Reuse `scripts/sorting_hat/pfam_hmmsearch.sbatch` (it takes the family table and writes `provenance.json`). Add a short driver
   `analysis/hydrophobin_truth/run/hydrophobin_runs.sbatch` only if a loop over proteomes is needed; test it with a dry run on one small FASTA. The driver uses
   `$SCRATCH` and absolute paths. One job for all proteomes of this step (sizing rule 3a).
2. Run `hmmsearch --cut_ga` of the full family table on the eight in-scope proteomes (old domain tables predate the new families).
3. Run the module wrappers: `pfam_hydrophobin`, `pfam_hsba`, `pfam_adhesion`, `cys8_pattern`, plus the TMHMM module (CFEM is active with `no_tm`, so
   the `pfam` command refuses to run without `--tm-module`).
4. Run the core `cellsurface_sorting_hat` run per proteome to write `run.json` with the current `config_sha256`. `calibrate truth --call-status` needs it.
5. Check each module state is `ok` and count `not_assessable` rows.
Acceptance: module tables and `run.json` exist; `config_sha256` in each `run.json` equals the hash of the current `categories.yaml`; states are `ok`.

### H5b. Runs on the extension proteomes (owner D7: yes, local files)

Steps: for the eight E1 proteomes in `Fungi_5k/input`, run SignalP 6 (R0) once as a single job on `short_gpu` (pass `-p short_gpu` at submit; the sbatch
header says `exfab`), then TMHMM, Pfam, the module wrappers and the core run, as in H5. Record the SignalP run record and the taxon ID per species. Acceptance: module tables, R0 tables
and `run.json` exist for every E1 proteome, or the missing species is listed with the reason.

### H6. Measurement

Files: `analysis/hydrophobin_truth/measure.py`, `docs/reports/2026-10-xx-hydrophobin-validation.md`, status files.
Steps:
1. Test first: sensitivity, specificity, Wilson interval and cluster bootstrap on a fixture with known answers, and the rule 6.3 function.
   Expected values come from an independent calculation (for example `statsmodels` or `scipy`), stored in the test. Cases: 0/0 gives "no evidence";
   4/4 gives "no evidence" (below the 5-cluster floor; Wilson lower bound 0.51); 5/5 gives "keep" (lower bound 0.57); 4/5 gives "no evidence" (lower bound 0.38).
   Cluster rule: rescue-only calls are clustered (MMseqs2, 30% and 0.5) across proteomes; a cluster whose members get different outcomes is `unresolved`.
2. Measure `pfam_hydrophobin` (module status) and, for the rescue, the sequence-level table and per-proteome rescue-only calls.
   Metrics as in spec section 6.1. T1+T2 and T3 in separate tables. Report the result with and without the "seen" proteins.
3. Review by hand every miss, every candidate false positive (including rescue-only calls and Pfam hits with no label), with BLAST top hit.
   Record outcomes (`hydrophobin`, `not_hydrophobin`, `unresolved`; T4 where no paper).
4. Apply rule 6.3. If the rescue is kept: add `hydrophobin_protein` to `categories.yaml` (after `hydrophobin_domain`, before the `other_*` calls;
   update the mechanism lists, tests and goldens as in H1) and measure the call. If not: record the result and leave `categories.yaml` unchanged.
4a. Produce the sequence-level tables: per-model hits, the cross-species hard-negative table, and the leave-one-species-out table on the Swiss-Prot set.
5. Write status files with `calibrate truth`, from the test part of the split (inputs from H2c) (leakage `partial`, with the note of spec section 5.3). Expect `smoke`. Do not write a status that the
   rules refuse.
6. Write the report. State every limit: counts, assumed negatives, partial leakage, the circularity of pattern-named labels.
Acceptance: report and status files exist and `calibrate truth` accepted them; every number in the report comes from a stored command output.

### H1b. Re-measure the repeat call after the last config edit

Steps: after the last edit of `categories.yaml`, **first re-run the core run for S288C and *C. albicans*** (the `run.json` `config_sha256` must equal the new config, or `--call-status` refuses), then re-run `calibrate truth --call-status` for `tandem_repeat_protein` in S288C and *C. albicans*. Compare with the earlier values
(sensitivity 0.435 and 0.462; specificity 0.987 and 0.976). If they differ, stop and report. Acceptance: the numbers match or the difference is explained.

### H7. Paper notes

Files: `docs/paper/02` (ledger rows), `docs/paper/04`, `docs/paper/05`, `docs/HANDOFF-2026-10-xx.md`.
Steps: add ledger rows for the hydrophobin truth set and measurement; update open questions; update the novelty note with H3b; hand-off section for the
next session. Acceptance: every number in the docs also appears in the report (`grep` each number); the hand-off names the commit hashes.

## Review and stop points

- This plan goes to an independent review before any Phase A code or analysis starts.
- After H3 the owner can inspect the frozen spacing. After H2 the counts are reported. If the counts say fewer than 20 T1+T2 clusters,
  continue (the result is `smoke` either way) and say so in the report.
- Stop and report if: a test cannot be made red first; a command output disagrees with this plan; the repeat-call numbers in H1b differ.
