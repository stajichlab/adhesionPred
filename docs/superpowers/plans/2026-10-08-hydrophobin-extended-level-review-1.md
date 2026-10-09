# Review 1 of the hydrophobin extended-level implementation plan

*2026-10-08. Independent review of `2026-10-08-hydrophobin-extended-level.md` (revision 1, commit 2690846), branch
`hydrophobin-validation`. Reviewer: Claude Code (claude-opus-5-5). The plan and the spec were not edited. Nothing is committed.
Every number below was computed by the reviewer with the commands in section 0, unless it is marked as quoted.*

**Pre-registration note.** Sections 2 and 3 contain proteome results (relaxed-Pfam calls per proteome, a random proteome
sample) computed before the freeze commit. The author will read them before the freeze. `freeze.json` and the leakage text
(spec 6.8) should record that these numbers were seen.

## 0. Method

- Code read: `modules/base.py` (`ModuleSpec`, `invalid_row`, `write_module`), `modules/pfam.py` (`parse_domtblout`),
  `modules/cli.py` (`_condition_table`, `pfam`, `cys8`), `modules/cys8.py`, `engine.py` (`call_eligible`), `categories.yaml`,
  `outputs.py` (`_known_limits`), `calibration/cli.py` (`_check_call_run`), `scripts/sorting_hat/*.sbatch`, `submit_modules.sh`,
  `tests/cellsurface_sorting_hat/test_data_and_scripts.py`, `test_engine.py`, `test_cli.py`, `test_outputs.py`.
- Test command of global constraint 1, run as written: **838 passed, 7 skipped** (the 7 are `shellcheck` skips). Matches the plan.
- Scratch copy of `src`, `tests`, `data/sorting_hat`, `scripts/sorting_hat`: `hydrophobin_extended` added after `hsba_domain` and to
  both `other_*` mechanism lists, then `pytest tests/cellsurface_sorting_hat`.
- HMMER 3.4. The seven hydrophobin-class models fetched with `hmmfetch -f` from
  `/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm` (resolves to `.../2026-01-27-Pfam38.2`).
  GA lines rewritten with `sed`.
- `hmmsearch -T 0 --domT 0 --tblout --domtblout`, default filters and `--nobias`, on the 12 proteomes of `run_list.tsv` and on
  `nopfam9.faa`. Joined to `step1_rule@R0.call`, `cys8_pattern.n_cys`, `pfam_hydrophobin.hit`, `pfam_hsba.hit` in each work
  directory and to `proteome_map.tsv`. A conditioned relaxed-only call = best full-sequence score >= cutoff, R0 `called`,
  `n_cys` >= 8, no `pfam_hydrophobin` and no `pfam_hsba` hit.
- A random sample of 1,000 unlabelled proteins per proteome (Python `random.Random(20261008)`, my own draw; it is not the plan's L2
  sample).
- MMseqs2 17-b804f `easy-cluster --min-seq-id 0.3 -c 0.5` on the 131 T2 sequences alone, and on those 131 plus the 50
  `jensen_sequences.faa` (neutral IDs). These were compared with `clusters_positives.tsv`.
- mafft 7.505 (`--auto`) and `hmmbuild` on 5 T2 sequences, on all 131 T2 sequences, and on 1 sequence.
- Scratch outputs are in the session scratchpad, not in the repository.

## 1. Spec requirements and review 2 findings against plan tasks

| Item | Task | State |
|---|---|---|
| N1 recovery definition, `u` rule | L5 step 1 | mapped |
| N1 tie-break | L4 step 2 (ties go to the higher cutoff) | mapped |
| N2 search options frozen | L4 steps 1, 7, 8 | partial: the rule that chooses the option is not defined (B2) |
| N3 full-sequence cutoff, domain GA 0.0 | L4 step 3 | partial: domain GA 0.0 drops some full-sequence hits (M2) |
| N4 cutoff from negatives only | L4 steps 2, 7 | mapped |
| N5 one ship rule | L6 step 4 | ordered before its inputs exist (B1) |
| N6 tuning proteins in their own column | L7 | partial: the clustering of extra calls for 6.5 has no step (m8) |
| N7 reserve by a rule without v1 output | L1 | mapped |
| N8 cost limit does not bind | none | not addressed; it binds at low cutoffs with `--nobias` (B2) |
| N9 fixed sample defined | L2 | mapped; overlaps the cost proteomes (m7) |
| N10 integration files | L6 | partial: several tests and `outputs.py` are not named (M3) |
| N11 provenance, seven GA lines, base name and `qlen` | L4 steps 3, 4; L6 | partial: `qlen` check not named (m9) |
| N12 `tuned_on_truth` for relaxed | L7 | mapped |
| N13 L1b flagged entries; *P. ostreatus* strain | L1b | partial (m4, m9) |
| Spec 6.2 HMM cost <= relaxed + 2 per 10,000 | none | missing (M5) |
| Spec 6.5 precision for HMM-only calls | none | missing (M5) |
| Spec 6.1 sample "cluster-split with the hard negatives" | L2, L3 | clustered separately (M6) |
| Spec 3.3 one shared cluster definition | L1, L2 | two different clusterings (M6) |
| E11 HsbA | L3, L5 step 3, L7 | partial: HsbA threshold and "HsbA protein" not defined (m3) |
| Spec 6.5 "`hydrophobin_relaxed` stays as a reported module" | L6 step 4 | how it is reported is not stated (m6) |
| Spec 3.2 rule 5, 6.7, 6.8, 7, 9.4 | L1, L5, L7, L8 | mapped |

## 2. Facts checked

| Plan or spec claim | Re-derived | Verdict |
|---|---|---|
| Baseline 838 passed, 7 skipped | 838 passed, 7 skipped | correct |
| 50 Jensen sequences, 131 T1/T2 | 50; 131 T2, 0 T1 | correct |
| L2: 28 folds, 6 flagged (`clusters_positives.tsv`, T2 rows) | 131 T2 rows in 28 clusters; Pfam-missed clusters {HFBE_PENEN} (2), {HFB3, HYD1_TRIAP, HYD2_GIBZE} (3), {HYD1_GIBZE} (1), {HCF1, HCF2} (3), {PSH_FLAVE} (1), {HCF4} (1) | correct |
| 8 unverified Jensen IDs exist in `jensen_sequences.faa` | all 8 present once (as `jensen|ID`) | correct |
| L6b: S288C 0.435 [0.067, 0.667], 0.987; *C. albicans* 0.462 [0.000, 0.788], 0.976 | `_workdir/sorting_hat/*/status/calls/tandem_repeat_protein.json`: same | correct |
| `short`, `short_gpu` cap at 2 hours | `sinfo`: 2:00:00 both | correct |
| `current` Pfam link | `2026-01-27-Pfam38.2`; the seven GA lines are all `NN.NN NN.NN` | correct |
| Spec 1.5: 8.5 bits, default filters, conditions: 0 to 3 extra calls; 10 in total | 10 (Calb 1, Bder 1, Fful 3, Fgra 3, Pexp 1, Post 1); 6 unlabelled, 4 labelled | correct |
| Spec 1.5: 15 bits: 0 or 1 extra calls | 0 or 1 in every proteome | correct |
| Spec 1.3: bit scores of the 9 (default filters) | 23.9, 23.6, 18.6, 15.9, 8.6; four with no row at `-T 0` | correct |
| Spec 1.6: HCF1 30.4, HCF4 16.6 with `--nobias` | 30.4, 16.6 | correct |
| `--nobias` scores not in the spec | HFBE_PENEN 18.6 (DewD), HYD1_TRIAP 4.3, HFB3_HYPVG 2.4; all 9 have a score >= 0 | new |
| `call_eligible(hydrophobin_extended)` | `(True, "")` in the scratch config | correct |
| mafft 7.505, MMseqs2 17-b804f installed | yes (`module avail`) | correct |

Feasibility checks of item 3 of the request:

- **Offline thresholding vs `--cut_ga`.** Three sequences (PSH_FLAVE, HFBE_PENEN, HYD2_GIBZE), GA rewritten to `C 0.00;` for C = 8.5,
  15.0, 16.0, 23.95, default and `--nobias`. The `--cut_ga` tblout, the `--cut_ga` domtblout and the offline threshold on the
  `-T 0` tblout gave the same proteins in all 8 runs. The header has `# Query file:`, the option line with `--cut_ga`, and `# [ok]`.
  The equivalence fails for proteins whose domains all score below 0 (M2).
- **`--nobias` with `--tblout`/`--domtblout` at a low threshold** gives what L4 step 1 needs (full-sequence and domain scores). All 24
  proteome searches plus the 9-sequence searches ran in one short interactive call on 2 CPUs. Output: 368 KB in total.
- **MMseqs2 on LP plus T2** runs (181 sequences). Results in M6 and in the L1 numbers below.
- **mafft and hmmbuild** run on 5 sequences, on 131 sequences, and `hmmbuild` runs on 1 sequence. On all 131 T2 (`mafft --auto`): 1,735
  columns, `LENG 134`, `EFFN 4.5`. Eight columns have Cys in >= 95% of sequences; no other column reaches 50%. The 80% check passes on
  the full set.
- **L1 numbers (my run).** 12 LP clusters; 8 contain a T2 protein; 4 contain no T2 and no unverified protein. Reserved: {An09g05530,
  JGI43184}, {ATEG_10285}, {AFLA_060780}, {AFLA_063080} (5 proteins). ACLA_001890 and ATEG_08089 cluster with RODF_ASPFU (as review 2).
  So the reserve is >= 3 clusters. Three of the 7 Jensen Pfam misses are reserved. The 6.7 secondary set is then AO090012000143,
  AFLA_014260 (identical sequences) and ATEG_08089, after ACLA_001890 (unverified) is removed: 3 proteins, 2 unique sequences.

## 3. Findings

### B1. BLOCKER. The ship-rule gate (L6 step 4) comes before its inputs, and no task applies the rule after L8

Evidence.
- L6 step 4: "Only if L5 and L7 pass the ship rule, edit `categories.yaml`". L7 "needs L6". So L7 has not run when L6 decides.
- The ship rule (spec 6.5) needs the precision term. That needs owner T5 decisions on evidence sheets (L8). L8 comes after L7.
- L7 writes status files for `hydrophobin_extended` with `calibrate truth --call-status`. `_check_call_run` refuses a call that is
  not in the config (`call_eligible` returns "unknown call") and refuses a run whose `config_sha256` differs from `--config`.
  So L7 cannot write them unless the call is already in a config.
- L8 says "Apply the precision rule of spec 6.5 once decisions exist". No task then edits `categories.yaml`, regenerates golden
  files, runs L6b, or writes status files.

Fix. Re-order:
1. L6a: module `hydrophobin_relaxed`, converter, job script. No `categories.yaml` edit.
2. L7a: proteome runs. Cost from the module tables (relaxed hit, no `pfam_hydrophobin`, no `pfam_hsba`, no label). This is what the
   reviewer did, so no call is needed.
3. L8: evidence sheets. Owner stop.
4. L8b: `ship_rule.py` (tested, red first) reads L5, L7a and `owner_decisions.tsv` and writes `ship_decision.json`.
5. L6c (only if L8b passes): `categories.yaml`, engine tests, golden files. Then L6b. Then L7b: core runs and status files.

`cellsurface_sorting_hat --config` and `calibrate --config` accept another file. A candidate config can be used if status
files are wanted before the decision. State which.

### B2. BLOCKER. The values that decide the frozen cutoff are not set by any task, and on the data the proteome sample does not bind

Evidence.
- L4 step 8 freezes "rate `r`". No task chooses `r`. Spec 6.1 points to 6.4, which is the per-group pass rate for the **test** part
  (default 5%). That is a different quantity.
- `lowest_cutoff` has no floor. The lowest possible cutoff is the reporting threshold of L4 step 1 ("a low reporting
  threshold"). Its value is not named.
- L4 step 7: "search option chosen the same way". With negatives only, `lowest_cutoff` gives a cutoff per option. It gives no
  rule to choose between options. Choosing by recall would use positive scores, against L4 step 2.
- My sample of 12,000 unlabelled proteins (1,000 per proteome): conditioned relaxed-only calls at >= 0 bits: **0** (default), **7**
  (`--nobias`); at 5 bits 0 and 2; at 8.5 bits 0 and 2; at 15 bits 0 and 0. A rate of 7 in 12,000 is 0.06%. With a 50% tuning part
  (about 6,000) the smallest non-zero rate is about 0.017%. So for any `r` above about 0.06%, the proteome sample returns the floor.
  The hard-negative tuning part might bind. It is not built yet, so that is unknown.
- At 0 bits the outcome is set. Default filters: 5 of 9 Pfam-missed proteins in 5 of 6 clusters (HCF4 has no row), so u = 1.
  `--nobias`: all 9 in 6 of 6 (lowest HFB3_HYPVG 2.4), so u = 0. In both cases the HMM is "not testable on the primary set".
  R0 has not been run on the Swiss-Prot sequences, so these are upper bounds.
- At 0 bits with `--nobias`, unlabelled extra calls per 10,000 proteins exceed the 5 limit of 6.3 in 4 proteomes: *F. graminearum*
  12.51, *B. bassiana* 11.61, *P. ostreatus* PC9 6.68, *B. dermatitidis* 6.17. Default filters at 0 bits: maximum 2.86 (PC9).
  So N8 ("the limit does not bind") holds for default filters only.

Fix. Before L4 step 7, the owner sets and the plan writes into `freeze.json`:
- `r` for the hard-negative tuning part and for the sample tuning part, as a rate or a count, and whether they are pooled.
- A floor (a minimum cutoff in bits), and the reporting threshold of L4 step 1 at or below it.
- The candidate cutoff grid (for example every observed negative score plus the floor).
- The option rule. Either fix one option now (for example `--nobias`, with default filters as a reported sensitivity column). Or
  choose by recall on the training Pfam-missed proteins inside each fold, and say that positives are used for the option only.
- Add tests: `lowest_cutoff` returns the floor when no negative scores above it; the option rule on a fixture.

### M1. MAJOR. Results are produced, and in L4 step 5 printed, before the freeze

Evidence.
- Constraint 7: "No held-out cluster result, no proteome result and no LP result is looked at before that commit."
- L4 step 5 prints "extra calls per proteome with the R0 and 8-cysteine conditions" at 8.5 and 15 bits, with and without
  `--nobias`, before step 8. The `--nobias` proteome numbers are not in the spec, so this is new proteome information.
- L4 step 1 scores "all reference sequences, the 12 proteomes and the sample". That includes LP proteins, the hard-negative **test**
  part and the sample **test** part. L4 step 6 runs `hmmsearch` of the held-out cluster for each fold. All exist on disk before the
  freeze.
- Constraint 7 does not name the hard-negative test part or the sample test part.
- `lowest_cutoff` uses no positive score (L4 step 2). So nothing positive needs to be scored before the freeze. The cysteine-column
  check (L4 step 6) uses the alignment, not scores.

Fix (cleanest).
- Before the freeze: score only the hard-negative tuning part and the sample tuning part (both options), per fold for the HMM.
  Build per-fold alignments and HMMs and run the cysteine check.
- Move L4 step 5 to the first step after the freeze (L5 step 0). The freeze does not depend on those numbers. The "stop if a
  number differs" rule still applies there.
- Move every search of held-out clusters, LP, test parts and full proteomes to L5 and L7.
- Add the test parts to constraint 7. Record the sha256 of every pre-freeze score file in `freeze.json`. L5 scripts check that the
  freeze commit is an ancestor of `HEAD` (`git merge-base --is-ancestor`) and that `freeze.json` is unchanged.

### M2. MAJOR. Domain GA 0.0 drops proteins that the offline full-sequence threshold counts

Evidence.
- With `--cut_ga`, a domain row is written only if the domain score is at least the domain GA. `parse_domtblout` reads domain rows only.
- `--nobias -T 0 --domT 0` on the 12 proteomes: 5 proteins have a full-sequence score of 0.2 to 4.3 bits and no domain at >= 0
  (YPR064W 4.3, F0349401_006848-T1 2.6, KAK9646131.1 2.0, XP_001245331.2 0.7, F1BB8A46_006295-T1 0.2). At model level, RODE_ASPFU
  scores 6.2 to Eas with best domain -3.2.
- With B2, a cutoff in this range is likely.
- `GA    6.00 -100.00;` is accepted by `hmmsearch --cut_ga` (exit 0). The domtblout then has the RODE_ASPFU Eas row at -3.2.
  `parse_domtblout` read the 28 rows.

Fix. Set the domain GA to a large negative value (for example -1000.00, test that HMMER accepts it), or let the wrapper read the
`--tblout` full-sequence score. Add a test fixture with a protein whose domains all score below 0. The L4 step 3 equivalence test
must include such a protein.

### M3. MAJOR. L6 does not name the tests and code that the call change breaks

Evidence (scratch copy, call after `hsba_domain` and in both `other_*` mechanism lists): 9 failures.
- `test_engine.py`: `test_packaged_config_loads_and_has_the_spec_calls` (fixed call list);
  `test_other_not_surface_leaves_out_unknown_calls_and_names_them` (`other_basis` gains `hydrophobin_extended`, because `base_modules`
  has no `hydrophobin_relaxed`); `test_bad_config_is_refused` (two cases use `d["calls"][6]`, which moves).
- `test_cli.py`: `test_golden_calls`, `test_wide_table_and_run_json`, `test_report_has_the_known_limits_and_counts`,
  `test_absent_required_modules_are_listed`, `test_run_output_is_pinned` (golden files).
- `outputs._known_limits` describes `hydrophobin_domain`. `test_outputs.py::test_the_spec_and_the_report_list_the_same_known_limits`
  requires it to equal the list in `docs/superpowers/specs/2026-10-04-orchestrator-design.md` (9 items). A text change needs both.
- `test_data_and_scripts.py`: `test_submit_modules_starts_every_job_and_each_has_a_script` (exact job list, and the set of all
  `.sbatch` files must equal it); `test_submit_modules_passes_absolute_paths_to_sbatch_and_supports_dry_run` (5 rows, 5 dry-run lines).
  The plan names only "SCRATCH:?, trap, job set".
- `test_engine.py::test_more_bad_config_is_refused` uses `d["calls"][-2]` and `[1]`. They pass now, but depend on the position.

Fix. List these files in L6 step 3. Add `hydrophobin_relaxed` to `base_modules` and to `build()` in `test_cli.py`. Replace the
index-based mutations with name lookups. State the new known-limits text and the orchestrator-spec edit.

### M4. MAJOR. Tests that call HMMER, mafft or MMseqs2 cannot run under the plan's test command

Evidence. In the shell that runs the test command, `hmmsearch`, `hmmbuild`, `mafft` and `mmseqs` are not on `PATH` (`which`).
L4 step 3 ("a run of `hmmsearch --cut_ga` on a small synthetic FASTA") and L4 step 6 ("Test first on a small synthetic family")
need them. They will fail, or, with `skipif`, skip, so red-first cannot be shown and the skip count changes.

Fix. Split each such test into a pure-function test (fixture tables) and a tool test marked `skipif(shutil.which(...) is None)`.
Give a second command: `module load hmmer/3.4 mafft/7.505 mmseqs2/17-b804f` then the same `pytest` line, with the expected count.

### M5. MAJOR. The HMM's cost and precision are never measured

Evidence. Spec 6.2 adds the HMM only at "a cost no higher than the relaxed baseline plus 2 calls per 10,000 proteins" and if it
passes 6.4 and 6.5. L7 runs only "hmmsearch with the shipped `.hmm`". L8 makes sheets for "every extra call" of the shipped level.
No task runs the all-cluster HMM on the 12 proteomes or makes sheets for HMM-only calls.

Fix. If u >= 2 and the HMM passes the recall part, L7 also runs the all-cluster HMM with its frozen cutoff on the 12 proteomes, and
L8 includes HMM-only calls. If u is 0 or 1, say that this step is skipped and why. B2 shows u is likely 0 or 1.

### M6. MAJOR. Three different cluster definitions; one changes the Pfam-missed clusters

Evidence.
- `clusters_positives.tsv` (T2 plus T3), T2 alone, and the L1 joint run (T2 plus LP) all give 28 T2 clusters. The memberships differ.
  T2 plus T3 vs T2 alone: 4 clusters differ (HUM2_MYCMD, HYD20_PLEO1, BHP3_BOTFB and HYD2_BIOOC move between clusters).
- Joint T2 plus LP: **5** Pfam-missed clusters, not 6. HYD1_GIBZE joins {HFB3, HYD1_TRIAP, HYD2_GIBZE}. HFBE_PENEN separates from
  RODF_ASPFU.
- Spec 3.3: "One shared definition of cluster is used for all tests." Spec 6.1: the sample is "cluster-split with the hard
  negatives". L2 and L3 cluster the sample and the hard negatives separately.

Fix. State in L2 that folds come from `clusters_positives.tsv` (record its sha256 in `freeze.json`). State in L1 that the joint
clustering is used only for the reserve. Either cluster the sample and the hard negatives together (spec 6.1), or record the
deviation in the plan.

### m1. MINOR. L4 builds the shipped `.hmm` (step 4) before the cutoff exists (step 7)

Fix. Step 4 builds a template and records the source hashes. Step 8 writes the final GA lines.

### m2. MINOR. The aligner mode and the action on a failed cysteine check are not stated

The plan names mafft 7.505 but no mode. `--auto` worked on 131 T2 sequences (section 2). The smallest leave-one-cluster-out
training set has 91 proteins (131 minus the 40-protein cluster), and `hmmbuild` runs on 1 sequence. So a fold with "few sequences" is
not a risk. Fix: freeze the mafft command line. State what happens when a fold fails the 80% check (stop and report, no silent
skip).

### m3. MINOR. E11 is not complete in L3 and L7

L3 step 4 sets "the per-group threshold" for every group. E11 excludes HsbA from the hard-negative rate. Fix: state that the HsbA
group has no pass threshold. Define an "HsbA protein" in proteomes (for example `pfam_hsba` hit at GA). Proteins near HsbA but below
its GA then count as cost; say so.

### m4. MINOR. L1b flagged entries (N13)

"Nothing is relabelled here" does not say whether an owner relabel after L2 is allowed. A relabel changes folds after the freeze.
Fix: flagged entries stay T2 for v1, unless the owner decides before L2. Record the choice in `freeze.json`.

### m5. MINOR. Test files and acceptance checks not named

L1b, L2, L3 say "test" without a file. `lco.py` has no test file. L8 and L9 have no command that checks acceptance. Fix: name the
files (`tests/hydrophobin_truth/test_review_t2_text.py`, `test_make_folds.py`, `test_make_hard_negatives.py`, `test_lco.py`).
For L9, a script that checks that each number in the docs appears in a stored table.

### m6. MINOR. A failed relaxed level still adds a job to every run, and "reported" is not defined

L6 step 2 adds the job to `submit_modules.sh` before the ship decision. If the level fails, every run pays for a search that no
call reads. Spec 6.5 says the module "stays as a reported module". No call or `evidence:` entry reads it, so the core report does not
show it. Fix: add the job to the standard run only with L6c (B1), or add `hydrophobin_relaxed.hit` to the `evidence:` list and say so.

### m7. MINOR. The tuning sample overlaps the cost proteomes, and the counted calls are not defined

The sample is drawn from the 12 proteomes that 6.3 measures. My draw holds 4 `pfam_hydrophobin` GA hits in 12,000 proteins
(default filters). L2 removes T1/T2/T3/LP proteins but not GA hits. Fix: the tuning rate counts conditioned relaxed-only calls (no
`pfam_hydrophobin`, no `pfam_hsba`), as 6.3. Report cost with and without the sample tuning part. L2 must say how LP proteins are
found in Af293 and A1163 (exact sequence match).

### m8. MINOR. The clustering of extra calls for 6.5 has no step

Spec 6.5 clusters unlabelled extra calls across proteomes. L8 has no clustering step or parameters. Fix: add it to L8 (MMseqs2,
30%, 0.5) and to the ship-rule script.

### m9. MINOR. Remaining review 2 items

- N11: the wrapper should also compare model lengths (`qlen`), not only names and the query-file base name.
- N8: write in the plan that the cost limit binds only at low cutoffs with `--nobias` (B2 numbers).
- N13: the *P. ostreatus* strain of Xu 2021 vs PC9 is not checked. Say so in the L7 report.

### m10. MINOR. The nested procedure is not nested for relaxed Pfam

No positive is an input, and the tuning negatives are the same in every fold. So every fold cutoff for relaxed Pfam equals the
shipped cutoff. The relaxed leave-one-cluster-out numbers are then the shipped-cutoff numbers on the training proteins. Fix: state
this in L5 and in the report.

### m11. MINOR. Drivers in `analysis/hydrophobin_truth/run/` are not linted

The lint tests read `scripts/sorting_hat/` only. Fix: run the same checks (`SCRATCH:?`, no `BASH_SOURCE`, `set -euo pipefail`,
no `/tmp`, `/usr/bin/python3.12`) on the new drivers, or add the directory to a test.

## 4. Compliance with the plan's constraints and the owner's rules

- SLURM, `${SCRATCH:?}`, no `BASH_SOURCE[0]`, absolute paths: stated in constraint 2. The `scripts/sorting_hat` lint tests enforce it
  there. Not enforced for `analysis/hydrophobin_truth/run/` (m11).
- 2-hour partitions: true. The L4 step 1 searches are small (24 proteome searches in one short call on 2 CPUs; 368 KB output). One
  job is enough.
- Compression: the score tables are small. `negatives_fixed.tsv.gz` and `hard_negatives.tsv.gz` are compressed. No issue.
- No push, PR or merge: stated. Owner stop points: after L4, L5, L8. B1 adds one: the ship decision after L8, before any config edit.
- Python 3.12 path and test command: correct.

## 5. Verdict

**Not ready to execute.** The paths, function names, scripts and the test command are real. The GA-rewrite mechanism and offline
thresholding agree on the tested proteins. mafft, hmmbuild and MMseqs2 run on the real sequences. Two blockers remain. B1: the ship
rule is applied before its inputs exist, and nothing applies it afterwards. B2: the freeze cannot be written, because the rate, the
floor and the option rule are not set; on the data the proteome sample does not bind, so these values decide the cutoff and the HMM
outcome. M1 (pre-freeze results) and M2 (domain GA 0.0) change what the freeze contains and should be fixed in the same revision.
M3 to M6 need a few sentences each. After a revision that fixes B1, B2 and M1 to M6, the plan can start at L1.
