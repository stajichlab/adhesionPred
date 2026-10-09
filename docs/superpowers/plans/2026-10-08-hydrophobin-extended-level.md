# Hydrophobin extended level (relaxed Pfam first, custom HMM tested against it): implementation plan

*2026-10-08. Revision 1, not yet reviewed. Spec: `docs/superpowers/specs/2026-10-08-hydrophobin-custom-hmm-design.md` (revision 3 plus
decision E11). Reviews of the spec: `...-design-review-1.md`, `...-review-2.md`. Author: Claude Code (claude-sonnet-5-5).*

## Global constraints

1. Python is `/usr/bin/python3.12`. Tests: `PYTHONPATH=$PWD/src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat tests/hydrophobin_truth tests/calibration_truth -q -p no:cacheprovider`.
   Analysis tests live in `tests/hydrophobin_truth/` and load scripts with `importlib`. Baseline before this plan: 838 passed, 7 skipped.
2. SLURM for anything over a few minutes. Partitions `short` and `short_gpu` cap at 2 hours. Use `$SCRATCH` (`${SCRATCH:?}`), absolute paths, never `BASH_SOURCE[0]`.
   Group work into jobs of about 1 to 1.5 hours. One-off drivers go in `analysis/hydrophobin_truth/run/`, not `scripts/sorting_hat/`.
3. Large text outputs are compressed (`.gz`). Scripts that read large text accept `.gz`.
4. Tests first (red, then green). Commits: ruff and ruff-format run as hooks. After a failed commit, fix, re-add, re-commit, check `git log`. Commit messages end with
   `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
5. No push, PR or merge without the owner. Work stays on `hydrophobin-validation`.
6. Report only numbers that a stored command produced. State unknowns as unknown. Every status is `smoke`.
7. **Pre-registration.** Task L4 ends with a freeze commit. No held-out cluster result, no proteome result and no LP result is looked at before that commit. The freeze fixes: the aligner, the
   cutoff procedure (spec 6.1), the search options candidates, the thresholds of 6.2 to 6.5, the fixed negative sample, the hard-negative rates and the fold list.
8. Stop and report if: a test cannot be made red first; a command output disagrees with this plan or the spec; a number in the spec's section 1 does not re-derive.

## Dependencies

L0 is done. L1, L1b, L2, L3 do not need each other except where noted. L4a needs L1 to L3 (sequence lists). L4 needs L2, L3, L4a. L5 needs the L4 freeze. L6 needs L5. L6b after the last
`categories.yaml` edit. L7 needs L6. L8 needs L7. L9 last.

## L1. Literature check, LP clusters, v2 reserve

Files: `analysis/hydrophobin_truth/literature/lp_clusters.py` (new), `lp_table.tsv`, `v2_reserve.tsv`, tests `tests/hydrophobin_truth/test_lp_clusters.py`.
Steps:
1. Test first: a pure function `reserve(lp_clusters, t12_clusters, unverified)` returns LP clusters with no T1/T2 member and no unverified protein; fixtures with one cluster that contains a T2 protein, one with an unverified
   protein, one clean.
2. Take `jensen_sequences.faa` (50) and the T1/T2 sequences (131). Cluster together with MMseqs2 17-b804f (30% identity, coverage 0.5, neutral IDs because MMseqs2 rewrites `sp|ACC|NAME`).
3. For the 8 proteins whose stated cysteine pattern is not found (An07g03340, An08g09880, AN6401.2, ACLA_001890, ACLA_048810, ACLA_072820, ACLA_018290, ACLA_007980): print the actual gaps found in the resolved sequence, the
   cysteine count and length, and the stated pattern side by side. Mark each `matches_after_check` or `unresolved`. No pattern is changed from memory.
4. Write `lp_table.tsv` (id, species, cluster, in_t12_cluster, unverified, reserved). Print counts: LP clusters, clusters with a T1/T2 member, reserved clusters.
5. Look for T1 candidates: grep the Jensen and Xu supplements and `NOTES.md` for protein-level statements. Report "none" if none.
Acceptance: tests pass; counts printed and stored; if reserved clusters < 3, write that in the table header and in `README.md` (the v2 test then uses post-freeze T1 only).

## L1b. T2 label review

Files: `analysis/hydrophobin_truth/review_t2_text.py`, `t2_text_review.tsv`, test.
Steps:
1. Test first: a function that classifies a FUNCTION or SUBCELLULAR LOCATION text as `hydrophobin_related` (mentions hydrophobin, rodlet, surface hydrophobicity, spore wall, aerial hyphae, amphipathic, self-assembl)
   or `other` on fixture strings.
2. Run over `uniprot_entries.json.gz`. Write one row per T2 entry with the experimental text and its class.
3. Print the `other` entries. They are listed for the owner. Nothing is relabelled here.
Acceptance: tests pass; the list is stored; the count of `other` entries is printed.

## L2. Folds and the fixed negative sample

Files: `analysis/hydrophobin_truth/make_folds.py`, `folds.tsv`, `negatives_fixed.tsv.gz`, test.
Steps:
1. Test first: `loco_folds(clusters)` yields one fold per T1/T2 cluster (28 now) with the held-out cluster IDs; `sample_negatives(proteome_ids, exclude, n, seed)` is reproducible and removes excluded IDs;
   the cluster split of the sample puts whole clusters in one part.
2. Folds: one per T1/T2 cluster from `clusters_positives.tsv` restricted to T2 (28). Pfam-missed clusters are flagged (6).
3. Fixed sample: 1,000 random proteins per measured proteome (12 proteomes of `run_list.tsv`), seed 20261008, minus any T1/T2/T3/LP protein. Cluster all sampled proteins with MMseqs2 and split the clusters
   into a tuning part (50%) and a test part (50%) by seed.
4. Record all of it. Nothing is scored.
Acceptance: tests pass; `folds.tsv` has 28 folds and 6 flagged; the sample file has 12,000 rows or fewer with the removals stated.

## L3. Hard-negative sets

Files: `analysis/hydrophobin_truth/make_hard_negatives.py`, `hard_negatives.tsv.gz` (group, accession, species, cluster, part), `hard_negative_queries.md`, test.
Groups (spec E6): CFEM (PF05730), cerato-platanin (PF07249), HsbA (PF12296), PIR (PF00399) and Ccw12-like, small secreted cysteine-rich proteins (killer toxins, effector candidates).
Steps:
1. Write the exact UniProt queries and Pfam family accessions for each group to `hard_negative_queries.md` first (reviewed Swiss-Prot, Fungi). For the small secreted cysteine-rich group, define: reviewed fungal entries with
   keyword Toxin or Secreted, length 40 to 250, at least 6 cysteines, no hydrophobin name and no hydrophobin-class Pfam hit. State the rule before fetching.
2. Test first: group assignment from a Pfam hit table (a protein with a hydrophobin-class hit is removed from every group and listed in `removed_hydrophobin_like.tsv`).
3. Fetch, hmmsearch the group Pfam models for membership, cluster (30%, 0.5), split clusters into tuning and test parts by seed.
4. Print counts per group and part **before any scoring**. Fix the per-group threshold (default: call rate at most 5% of the test part) in `hard_negative_thresholds.json`.
Acceptance: tests pass; counts printed; thresholds file committed before L5.

## L4a. R0 SignalP on the reference sequences

Files: `analysis/hydrophobin_truth/run/r0_reference.sbatch` (one job), outputs under `analysis/hydrophobin_truth/r0_reference/`.
Steps:
1. Make one FASTA: T1, T2, T3, LP (the 50 Jensen), hard negatives, fixed negative sample. Neutral IDs.
2. One job on `short_gpu` (`-p short_gpu` at submit) reusing `scripts/sorting_hat/signalp_gpu.sbatch` inline (`bash -l`), same build and mode (6.0h fast, eukarya). Record the SignalP module and version.
3. Convert with `cellsurface_sorting_hat_module signalp` so the R0 call has the same identity logic.
Acceptance: R0 table exists for every sequence; the identity is recorded; counts of `called` per tier are printed.

## L4. Relaxed Pfam procedure, HMM per fold, freeze

Files: `analysis/hydrophobin_truth/relaxed.py`, `lco.py`, `make_relaxed_hmm.py`, `data/sorting_hat/hydrophobin_relaxed.hmm` (written at the freeze), tests.
Steps:
1. **Scores once.** Run `hmmsearch` of the seven hydrophobin-class Pfam models over all reference sequences, the 12 proteomes and the sample, in two option sets (default filters; `--nobias`), at a low reporting
   threshold, writing `--tblout` and `--domtblout`. One job per option set group (sizing rule). Store the full-sequence and domain scores. Cutoffs are then applied offline.
2. Test first (`test_relaxed.py`): `apply_cutoff(scores, cutoff)` with a hit decided by the full-sequence score; `lowest_cutoff(tuning_negative_scores, max_rate)` returns the lowest score with call rate at or below the rate and
   ties go to the higher cutoff; no positive score is an input to `lowest_cutoff`.
3. Test first (`test_make_relaxed_hmm.py`): the derived file has the seven models, `GA` lines `cutoff 0.0` (sequence cutoff, domain 0.0), and a run of `hmmsearch --cut_ga` on a small synthetic FASTA gives the same hits as
   offline thresholding at the cutoff; the domain table header has `--cut_ga`, the `# Query file:` line and `# [ok]`, and `pfam.parse_domtblout` accepts it. (Review 2 showed a file with `GA 8.50 8.50` works; the domain value is now 0.0.)
4. Build the derived `.hmm` with `hmmfetch` from `/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm`. Record the real path that `current` points to and the sha256 of the source models.
5. **Re-derive the review's numbers** (spec section 1, items 1 to 6): relaxed Pfam at 8.5 and 15 bits, with and without `--nobias`, extra calls per proteome with the R0 and 8-cysteine conditions, which of the 9 Pfam-missed proteins and 6
   clusters are reached. Print the table. If any number differs from the spec, stop and report.
6. HMM build (`lco.py`): per fold, align the training proteins with `mafft` (7.505; chosen over `famsa` unless a test shows a problem), `hmmbuild`, `hmmsearch` of the held-out cluster, the hard-negative tuning part and the fixed sample.
   A script checks the eight conserved cysteine columns (present in at least 80% of the training sequences, stated before running). Test first on a small synthetic family.
7. Nested cutoff per fold by `lowest_cutoff`, search option chosen the same way. The shipped cutoff over all 28 clusters by the same rule.
8. **Freeze commit.** Write `freeze.json` (aligner and version, cutoff procedure, option candidates, thresholds of spec 6.2 to 6.5, rate `r`, fixed sample hash, folds hash, hard-negative thresholds hash, the relaxed-Pfam source path and
   hashes). Commit. After this commit nothing in `freeze.json` changes.
Acceptance: tests pass; the review's numbers are re-derived and printed; `freeze.json` is committed before L5; no L5 output exists before it (check by file times and commit order).

## L5. Evaluation

Files: `analysis/hydrophobin_truth/evaluate.py`, `docs/reports/data/sorting_hat/hydrophobin_ext/` tables.
Steps:
1. Test first: recovery per cluster by the spec definition (one Pfam-missed member called at the fold cutoff), `u` (clusters not recovered by relaxed), and the HMM decision rule (u of 0 or 1 gives "not testable";
   u of 2 or more gives the ceil(u/2) rule); fixtures cover each branch.
2. Leave-one-cluster-out: Pfam GA (with and without the R0 and cysteine conditions), relaxed Pfam, `cys8_pattern`, HMM. Recall by cluster and by accession with Wilson intervals, labelled `partial`.
3. Hard negatives per group on the test part; HsbA in its own column (E11).
4. LP secondary test set: recall on the Jensen proteins that Pfam GA misses and that are not reserved, reported separately.
5. Count T2 entries lost to the 8-cysteine condition (6 expected) and labelled proteins with no R0 call.
Acceptance: tests pass; tables stored; every number comes from a stored script run.

## L6. Modules and call (only what the ship rule allows)

Files: `src/cellsurface_sorting_hat/modules/` (new wrapper for `hydrophobin_relaxed`), `modules/cli.py`, `data/sorting_hat/hydrophobin_relaxed.hmm`, `hydrophobin_relaxed.provenance.json`, `categories.yaml`, tests, golden files.
Steps:
1. Tests first: a wrapper `hydrophobin_relaxed` reads the domain table (parser needs `--cut_ga`, the `# Query file:` check, the model names), the R0 table (`_condition_table`), and the full-sequence cysteine count; `hit` = a hit and at least 8 cysteines and R0 `called`;
   empty `hit` if R0 is missing; `na_invalid` for invalid sequences; `params` hold cutoff, search options, score type, minimum cysteines and the R0 condition; the `.hmm` file is in `ModuleSpec.artefacts` so that a changed file changes the module identity.
2. Implement. Add the job script to the standard run set (`scripts/sorting_hat/`), and update the script lint tests (`SCRATCH:?`, trap, job set).
3. Engine tests first: `hydrophobin_extended` = `or` of `pfam_hydrophobin.hit` and `hydrophobin_relaxed.hit`, placed before the `other_*` calls, added to their mechanism lists; passes `call_eligible`; golden files regenerated only after reading every changed line.
4. **Ship rule gate (spec 6.5).** Only if L5 and L7 pass the ship rule, edit `categories.yaml`. If not, the module stays in the tool, `hydrophobin_extended` is not added, and the report says why. If the HMM passes, repeat steps 1 to 3 for `hydrophobin_ext`.
Acceptance: tests pass; the golden diff contains only intended lines.

## L6b. Re-measure stale call files

After the last `categories.yaml` edit: run the core command for S288C and *C. albicans* into a new `out_*` directory, then `calibrate truth --call-status` for `tandem_repeat_protein` with the same inputs as before. Expect the same numbers
(S288C 0.435 [0.067, 0.667], 0.987; *C. albicans* 0.462 [0.000, 0.788], 0.976). If they differ, stop and report.

## L7. Proteome runs and cost

Steps: one job per step over the 12 proteomes (hmmsearch with the shipped `.hmm`, wrappers, core run into `out_ext`), using the drivers in `analysis/hydrophobin_truth/run/` as models. Count extra calls per 10,000 proteins for every proteome, with labelled calls in
a separate column and HsbA hits in their own column (E11). Write status files for `hydrophobin_extended` where the truth tables allow (`calibrate truth --call-status`, leakage `tuned_on_truth`). Write the report
`docs/reports/2026-10-xx-hydrophobin-extended-level.md` with the limits of the spec section 9.
Acceptance: module states `ok`; `config_sha256` in each `run.json` equals the hash of the current `categories.yaml`; the report lists every number's source.

## L8. Evidence sheets and owner decisions

Steps: build the evidence sheet per extra call and per unlabelled strict call (length, cysteine count and gaps, R0, TMHMM, Pfam, relaxed and HMM scores, BLAST top hits against the 174 known hydrophobins and Swiss-Prot, any paper). Write
`analysis/hydrophobin_truth/evidence_sheets.tsv`. The owner's choices go in `owner_decisions.tsv` (id, tier T5, decision, reason, date). Apply the precision rule of spec 6.5 once decisions exist.
Acceptance: a sheet row for every extra call; the file is ready for the owner. No decision is written by me.

## L9. Docs

`docs/paper/02` ledger rows, `docs/paper/04` limits (state: every status `smoke`; the HMM test may be "not testable"; T5 is judgement), `docs/paper/05` prior art update, `docs/HANDOFF-2026-10-08-hydrophobin-extended.md`, memory entry.
Acceptance: every number in the docs appears in a report or table.

## Review and stop points

This plan goes to an independent review before any task starts. Stop points for the owner: after L4 (the freeze commit and the re-derived numbers), after L5 (the HMM decision), after L8 (decisions needed).
