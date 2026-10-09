# Hydrophobin extended level (relaxed Pfam first, custom HMM tested against it): implementation plan

*2026-10-08. Revision 2 (after plan review 1: 2 blockers, 6 major, 11 minor findings, all addressed; review in
`2026-10-08-hydrophobin-extended-level-review-1.md`). Spec: `docs/superpowers/specs/2026-10-08-hydrophobin-custom-hmm-design.md` (revision 4, decision E11).
Author: Claude Code (claude-sonnet-5-5). This revision has not been reviewed again.*

## Global constraints

1. Python is `/usr/bin/python3.12`. Test command A (no tools on PATH): `PYTHONPATH=$PWD/src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat tests/hydrophobin_truth tests/calibration_truth -q -p no:cacheprovider`.
   Baseline before this plan: 838 passed, 7 skipped. Test command B (tools on PATH): `module load hmmer/3.4 mafft/7.505 mmseqs2/17-b804f ncbi-blast/2.14.0+` then the same line. Tests that call HMMER, mafft or MMseqs2 are split into a pure-function
   test (fixture tables, runs under A) and a tool test marked `skipif(shutil.which(...) is None)` (runs under B). Each task states the expected counts under A and B when it is committed.
2. SLURM for anything over a few minutes. `short` and `short_gpu` cap at 2 hours. Use `$SCRATCH` (`${SCRATCH:?}`), absolute paths, never `BASH_SOURCE[0]`. Group work into jobs of about 1 to 1.5 hours. Drivers go in `analysis/hydrophobin_truth/run/`.
   `tests/hydrophobin_truth/test_run_scripts.py` (new, task L0b) applies the same checks as the `scripts/sorting_hat` lint tests to that directory (`SCRATCH:?`, no `BASH_SOURCE`, `set -euo pipefail`, no `/tmp`, `/usr/bin/python3.12`).
3. Large text outputs are compressed (`.gz`). Scripts that read large text accept `.gz`.
4. Tests first (red, then green). Commits: ruff and ruff-format run as hooks; after a failed commit, fix, re-add, re-commit, check `git log`. Commit messages end with `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
5. No push, PR or merge without the owner. Work stays on `hydrophobin-validation`.
6. Report only numbers that a stored command produced. State unknowns as unknown. Every status is `smoke`.
7. **Pre-registration (spec 6.1).** Nothing outside the *tuning* data is scored before the freeze commit of L4. Before the freeze: tuning proteomes, tuning hard-negative parts, per-fold alignments and HMM builds. After the freeze: held-out clusters,
   LP proteins, test hard-negative parts, test proteomes and all of L5 to L7. The sha256 of every pre-freeze score file goes into `freeze.json`. Every later script checks `git merge-base --is-ancestor <freeze commit> HEAD` and that `freeze.json` is unchanged.
   Numbers that the plan reviewer computed before the freeze are recorded in `freeze.json` as "seen by the author".
8. Stop and report if: a test cannot be made red first; a command output disagrees with this plan or the spec in a way that changes a decision.

## Order

L0b, L1, L1b, L2, L3, L4a, L4 (freeze) | owner stop 1 | L5 | owner stop 2 (HMM decision) | L6a, L7a, L8 | owner stop 3 (decisions) | L8b, L6c, L6b, L7b, L9.

## L0b. Lint test for run drivers

File `tests/hydrophobin_truth/test_run_scripts.py`. Test first: a fixture script without `SCRATCH:?` fails the check; the real drivers in `analysis/hydrophobin_truth/run/` pass. Acceptance: passes under A.

## L1. LP clusters and v2 reserve

Files: `analysis/hydrophobin_truth/literature/lp_clusters.py`, `lp_table.tsv`, `v2_reserve.tsv`, `tests/hydrophobin_truth/test_lp_clusters.py`.
1. Test first (pure function): `reserve(lp_clusters, t12_clusters, unverified)` returns LP clusters with no T1/T2 member and no unverified protein; fixtures for each case.
2. Cluster `jensen_sequences.faa` (50) with the 131 T2 sequences by MMseqs2 17 (30%, 0.5, neutral IDs). **This joint clustering is used only for the reserve.** Folds use `clusters_positives.tsv` (L2).
3. For the 8 proteins whose stated pattern is not found (An07g03340, An08g09880, AN6401.2, ACLA_001890, ACLA_048810, ACLA_072820, ACLA_018290, ACLA_007980): print the actual gaps, the cysteine count and length beside the stated pattern. Mark `matches_after_check` or `unresolved`. Nothing is changed from memory.
4. Find the Jensen proteins in the Af293 and A1163 proteomes by exact sequence match, so LP proteins can be excluded from tuning data.
5. Print counts (LP clusters, with a T1/T2 member, reserved; the reviewer found 12 clusters and 4 reservable with 5 proteins). Look for T1 candidates in the supplements and `NOTES.md`; report "none" if none.
Acceptance: pure tests pass under A; counts stored; if reserved clusters < 3, the table header and `README.md` say the v2 test uses post-freeze T1 only.

## L1b. T2 label review

Files: `analysis/hydrophobin_truth/review_t2_text.py`, `t2_text_review.tsv`, `tests/hydrophobin_truth/test_review_t2_text.py`.
1. Test first: a function classifying a FUNCTION or SUBCELLULAR LOCATION text as `hydrophobin_related` (hydrophobin, rodlet, surface hydrophobicity, spore wall, aerial hyphae, amphipathic, self-assembl) or `other`, on fixture strings.
2. Run over `uniprot_entries.json.gz`; one row per T2 entry. Print the `other` entries for the owner.
3. **Relabel policy:** flagged entries stay T2 for v1 unless the owner decides before L2. The choice is recorded in `freeze.json`.
Acceptance: tests pass; list stored; count of `other` printed.

## L2. Folds and the proteome split

Files: `analysis/hydrophobin_truth/make_folds.py`, `folds.tsv`, `proteome_split.tsv`, `tests/hydrophobin_truth/test_make_folds.py`.
1. Test first: `loco_folds(clusters)` gives one fold per T2 cluster; `split_species_groups(groups, seed, min_truth_groups)` is reproducible, gives 5 and 5, and puts at least two groups with truth positives in each part; the three *A. fumigatus* strains are one group.
2. Folds come from `clusters_positives.tsv` restricted to T2 (28 clusters; sha256 recorded). The 6 Pfam-missed clusters are flagged. T3 members of clusters are not training proteins.
3. Proteome split: the 12 proteomes of `run_list.tsv` form 10 species groups. Seed 20261008. Write `proteome_split.tsv` (proteome, group, part).
4. The fixed 1,000-protein sample of the earlier draft is dropped (it never limits a cutoff). Tuning uses whole tuning proteomes.
Acceptance: tests pass under A; 28 folds, 6 flagged; the split file has 12 rows and 5 and 5 groups.

## L3. Hard-negative sets

Files: `analysis/hydrophobin_truth/make_hard_negatives.py`, `hard_negatives.tsv.gz`, `hard_negative_queries.md`, `hard_negative_thresholds.json`, `tests/hydrophobin_truth/test_make_hard_negatives.py`.
Groups: CFEM (PF05730), cerato-platanin (PF07249), HsbA (PF12296), PIR (PF00399) and Ccw12-like, small secreted cysteine-rich proteins.
1. Write the exact UniProt queries and Pfam accessions to `hard_negative_queries.md` before fetching. Small secreted cysteine-rich group: reviewed fungal entries, keyword Toxin or Secreted, length 40 to 250, at least 6 cysteines, no hydrophobin name, no hydrophobin-class Pfam hit.
2. Test first: a protein with a hydrophobin-class hit is removed from every group and listed in `removed_hydrophobin_like.tsv`; group assignment from a hit table.
3. Fetch, hmmsearch the group Pfam models, cluster **within each group** (30%, 0.5), split clusters into tuning and test parts by seed. (Deviation from spec 6.1, which clustered the negatives together: groups are clustered separately.)
4. Print counts per group and part before any scoring. Fix the per-group threshold in `hard_negative_thresholds.json` (default: call rate at most 5% of the test part). **The HsbA group has no pass threshold (E11)**; its rate is reported.
Acceptance: tests pass; counts printed; thresholds file committed before L4.

## L4a. R0 SignalP on the reference sequences

One job on `short_gpu` (`-p short_gpu` at submit), reusing `scripts/sorting_hat/signalp_gpu.sbatch` inline (`bash -l`) on one FASTA: T1, T2, T3, LP, hard negatives. The proteomes already have R0 (`step1_rule@R0`). Convert with `cellsurface_sorting_hat_module signalp`. Record the SignalP version and mode.
Acceptance: R0 table for every sequence; identity recorded; counts of `called` per tier printed (this is a count of R0 calls, not a method result).

## L4. Relaxed procedure, per-fold HMM build, freeze

Files: `analysis/hydrophobin_truth/relaxed.py`, `lco.py`, `make_relaxed_hmm.py`, `ship_inputs.py`, tests `test_relaxed.py`, `test_make_relaxed_hmm.py`, `test_lco.py`, and (at the freeze) `freeze.json`.
1. **Scores, tuning data only.** `hmmsearch` of the seven hydrophobin-class Pfam models (derived template, step 4) over the **tuning proteomes** and the **tuning hard-negative parts**, both option sets (default filters, `--nobias`), reporting everything at or above 0 bits (`-T 0 --domT -1000`), storing `--tblout` and `--domtblout`. One job. Cutoffs are applied offline.
2. Pure tests first (`test_relaxed.py`): `apply_cutoff(scores, cutoff)` decides a hit by the full-sequence score; `lowest_cutoff(tuning_scores, per_proteome_limit, hard_neg_rates, floor)` returns the lowest score at or above the floor meeting both conditions, ties to the higher cutoff, and returns the **floor when no negative scores above it**;
   no positive score is an input; the option rule picks the option with more recovered training Pfam-missed clusters on a fixture.
3. Tool test (`test_make_relaxed_hmm.py`, under B): the derived file has the seven models and `GA` lines `cutoff -1000.00`; `hmmsearch --cut_ga` on a synthetic FASTA gives the same hits as offline thresholding at the cutoff, **including a fixture protein whose domains all score below 0**; the domain table header has `--cut_ga`, the `# Query file:` line and `# [ok]`; `pfam.parse_domtblout` accepts it.
4. Build a **template** `.hmm` with `hmmfetch` from `/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm`. Record the real path that `current` points to and the sha256 of the source models. The final GA lines are written at step 8.
5. HMM per fold (`lco.py`): for each of the 28 folds, align the training proteins with `mafft --auto` (7.505; the command line is frozen), `hmmbuild`, and run the cysteine-column check (the eight conserved columns present in at least 80% of training sequences; checked on the alignment, not on scores). A fold that fails the check **stops the run and is reported**.
   Test first (`test_lco.py`): pure fold logic and the column check on a synthetic alignment; tool test under B on a synthetic family. Score only the tuning proteomes and tuning hard-negative parts with each fold's HMM.
6. Cutoffs by `lowest_cutoff` for relaxed Pfam (one cutoff per option; no positive score used) and per fold for the HMM. The shipped HMM cutoff uses all 28 clusters.
7. The relaxed Pfam procedure is not nested: every fold gives the same relaxed cutoff (spec 6.1). The report says so.
8. **Freeze commit.** Write the final `hydrophobin_relaxed.hmm` (GA lines carry the cutoff for the option chosen by the option rule) only as a file in `analysis/`, not yet in `data/`. Write `freeze.json`: aligner and command line, cutoff rule, floor 0, the two options, thresholds of spec 6.2 to 6.5, hard-negative thresholds hash, `proteome_split.tsv` hash, `folds.tsv` and `clusters_positives.tsv` hashes, the relaxed Pfam source path and hashes, the L1b relabel choice, the sha256 of every pre-freeze score file, and the numbers the reviewers computed before the freeze (listed as "seen by the author").
Acceptance: tests pass under A and B (counts stated in the commit); `freeze.json` is committed; no file with a held-out, LP, test-part or test-proteome score exists (check by listing the score directory).

**Owner stop 1** after the freeze: the frozen cutoffs and options are shown to the owner.

## L5. Evaluation (after the freeze)

Files: `analysis/hydrophobin_truth/evaluate.py`, `tests/hydrophobin_truth/test_evaluate.py`, tables in `docs/reports/data/sorting_hat/hydrophobin_ext/`.
0. First: re-derive the spec section 1 numbers (relaxed Pfam at 8.5 and 15 bits, with and without `--nobias`, extra calls per proteome, which of the 9 and 6 clusters). If a number differs, **report it and keep going**; it cannot change the freeze.
1. Test first: recovery per cluster by the spec definition; `u`; the HMM decision rule (u of 0 or 1 gives "not testable"; u of 2 or more gives ceil(u/2)); fixtures for each branch.
2. Leave-one-cluster-out: Pfam GA (with and without the R0 and cysteine conditions), relaxed Pfam, `cys8_pattern`, HMM. Recall by cluster and by accession with Wilson intervals, labelled `partial`.
3. Hard negatives per group on the test part. HsbA rate in its own column.
4. LP secondary set: recall on the Jensen proteins that Pfam GA misses and that are not reserved. Reserved LP clusters are not scored.
5. Count T2 entries lost to the 8-cysteine condition (6 expected) and labelled proteins with no R0 call.
Acceptance: tests pass; tables stored; every number comes from a stored script run.

**Owner stop 2** after L5: the HMM decision (added, or "not testable", or failed) and the relaxed recall result.

## L6a. Module code (no `categories.yaml` edit)

Files: `src/cellsurface_sorting_hat/modules/hydrophobin_relaxed.py` (new), `modules/cli.py`, tests under `tests/cellsurface_sorting_hat/modules/`.
1. Tests first: the wrapper reads the domain table (`parse_domtblout`, `--cut_ga`, `# Query file:`, model names **and lengths** (`qlen`)), the R0 table (`_condition_table`), the full-sequence cysteine count; `hit` = a hit and at least 8 cysteines and R0 `called`; empty `hit` if R0 is missing; `na_invalid` for invalid sequences; `params` hold cutoff, search options, score type, minimum cysteines, R0 condition; the `.hmm` is in `ModuleSpec.artefacts` so a changed file changes the module identity; a domain row with a score below 0 is read.
2. Implement. Copy the frozen `.hmm` to `data/sorting_hat/hydrophobin_relaxed.hmm` with its provenance file. **The job script is not added to `submit_modules.sh` here** (so a failed level adds no cost to every run).
Acceptance: tests pass under A and B.

## L7a. Proteome runs and cost (module tables only)

One job per step over the 12 proteomes (hmmsearch with the frozen `.hmm` and options, the wrapper, no core run needed). Cost from the module tables: extra calls = relaxed `hit` with no `pfam_hydrophobin` and no `pfam_hsba` call and no T1/T2/LP label, per 10,000 proteins, reported for the **test proteomes** (out of sample) and the tuning proteomes (in sample), with labelled extra calls and `pfam_hsba` overlaps in their own columns (E11).
If u is 2 or more and the HMM passed the recall part, also run the all-cluster HMM on the 12 proteomes. If u is 0 or 1, this step is skipped and the report says why.
Acceptance: module tables `ok`; the cost table stored.

## L8. Evidence sheets and extra-call clusters

1. Cluster the extra calls across proteomes (MMseqs2, 30%, 0.5). The ship-rule script uses these clusters.
2. Evidence sheet per extra call and per unlabelled strict call: length, cysteine count and gaps, R0, TMHMM, Pfam, relaxed and HMM scores, BLAST top hits (the 174 known hydrophobins and Swiss-Prot 2023_03), any paper. Write `evidence_sheets.tsv`. The owner's choices go in `owner_decisions.tsv` (id, tier T5, decision, reason, date). I write no decision.
3. Note in the report that the *P. ostreatus* strain of the Xu 2021 gene list was not checked against PC9.
Acceptance: a sheet row for every extra call; file ready for the owner.

**Owner stop 3**: decisions on the sheets.

## L8b. Ship-rule script

`analysis/hydrophobin_truth/ship_rule.py`, `tests/hydrophobin_truth/test_ship_rule.py`. Test first on fixtures for each branch: recall rule (6.2), cost in every test proteome (6.3), per-group hard-negative rule with HsbA excluded (6.4), precision rule with at least 5 resolved unlabelled clusters and Wilson lower bound at least 0.5 (6.5), the combined ship rule. It reads L5, L7a, `owner_decisions.tsv` and the extra-call clusters and writes `ship_decision.json`.
Acceptance: tests pass; `ship_decision.json` states each condition and the result.

## L6c. Call and config (only if `ship_decision.json` says ship)

Files: `src/cellsurface_sorting_hat/categories.yaml`, `outputs.py`, `docs/superpowers/specs/2026-10-04-orchestrator-design.md` (limit 4 text, kept equal to the report text), `scripts/sorting_hat/submit_modules.sh` and the new job script, tests and golden files.
1. Tests first. Named tests that the change breaks and must be updated: in `test_engine.py` `test_packaged_config_loads_and_has_the_spec_calls`, `test_other_not_surface_leaves_out_unknown_calls_and_names_them`, `test_bad_config_is_refused` (replace `d["calls"][N]` by a name lookup), `test_more_bad_config_is_refused` (same); add `hydrophobin_relaxed` to `base_modules`. In `test_cli.py` `test_golden_calls`, `test_wide_table_and_run_json`, `test_report_has_the_known_limits_and_counts`, `test_absent_required_modules_are_listed`, `test_run_output_is_pinned` (add the table to `build()`; read every golden diff line). In `test_outputs.py` the known-limits test. In `test_data_and_scripts.py` `test_submit_modules_starts_every_job_and_each_has_a_script` (job list, `.sbatch` set) and `test_submit_modules_passes_absolute_paths_to_sbatch_and_supports_dry_run` (row and line counts).
2. `hydrophobin_extended` = `or: [{flag: pfam_hydrophobin.hit}, {flag: hydrophobin_relaxed.hit}]`, placed before the `other_*` calls, **added** to their mechanism lists. Add the job script to `submit_modules.sh` now. If the HMM passed, add `hydrophobin_ext` the same way (L6a and L6c steps for it).
3. If the ship decision says no: `categories.yaml` is not edited, the module stays in the tool as a reported module, and the report says why.
Acceptance: all tests pass; the golden diff contains only intended lines.

## L6b. Re-measure stale call files

After the last `categories.yaml` edit (or, if none, after confirming that none is needed): core run for S288C and *C. albicans* into new `out_*` directories, then `calibrate truth --call-status` for `tandem_repeat_protein` with the same inputs. Expect S288C 0.435 [0.067, 0.667] with specificity 0.987, and *C. albicans* 0.462 [0.000, 0.788] with 0.976. If they differ, stop and report.

## L7b. Core runs and status files

Only if the call was added: core runs on the 12 proteomes into `out_ext`, then `calibrate truth --call-status` for `hydrophobin_extended` in the species with truth, leakage `tuned_on_truth`, `smoke`. Report in `docs/reports/2026-10-xx-hydrophobin-extended-level.md` with the spec section 9 limits, the unlabelled-hydrophobin caveat, and that the cost limit binds only at low cutoffs with `--nobias` (reviewer numbers: *F. graminearum* 12.51, *B. bassiana* 11.61, PC9 6.68, *B. dermatitidis* 6.17 per 10,000 at 0 bits).

## L9. Docs

`docs/paper/02` ledger rows, `docs/paper/04` limits (every status `smoke`; the HMM may be "not testable"; T5 is judgement), `docs/paper/05` prior art, `docs/HANDOFF-2026-10-08-hydrophobin-extended.md`, memory entry. Acceptance: a script (`analysis/hydrophobin_truth/check_doc_numbers.py`, tested) finds each number of the report in a stored table.
