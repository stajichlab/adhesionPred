# Review 1: hydrophobin validation implementation plan (revision 1)

*2026-10-08. Independent review by Claude Code (claude-opus-5-5). Reviewed:
`docs/superpowers/plans/2026-10-08-hydrophobin-validation.md` (commit 1b4ac58, branch `hydrophobin-validation`)
against the spec (revision 3 plus section 9 owner answers) and spec reviews 1 and 2. The plan and spec were
not edited. Every statement below was checked against a file, code or a command output. Where I did not
check something, I say so.*

## How I checked

- Git: `git merge-base`, `git merge-base --is-ancestor`, `git log` on `hydrophobin-validation`,
  `per-call-status`, `signoff-cfem-hydrophobin`, `main`, `origin/main`; `git merge-tree --write-tree
  per-call-status HEAD` (clean); `gh pr list --state all`.
- Code on `per-call-status` with `git show`: `modules/cli.py` (`run`, `_write`, `_family_digest`,
  `_condition_table`), `calibration/cli.py` (`truth`, `_read_truth`, `_check_call_run`), `engine.py`
  (`_eval_other`), `status.py` (`status_from_measure`), `call_status.py` (`LEAKAGE_VALUES`),
  `tests/cellsurface_sorting_hat/test_cli.py` (`test_golden_calls`, `test_run_output_is_pinned`),
  `test_data_and_scripts.py`. On this branch: `modules/pfam.py`, `modules/base.py`, `categories.yaml`,
  `data/sorting_hat/family_table.tsv`, `scripts/sorting_hat/*.sbatch`, `pyproject.toml`.
- `git diff --stat HEAD per-call-status -- src tests`: `modules/` and `categories.yaml` are identical on both
  branches; `engine.py`, `outputs.py`, `call_status.py`, `calibration/cli.py`, `test_cli.py` and the golden
  files differ or exist only on `per-call-status`.
- Ran the test command of the plan's Global constraint 1 on this branch, with `tests/calibration_truth`
  added: 701 passed, 7 skipped.
- Files: `_workdir/sorting_hat/` (FASTA, module tables, run records, status files, `sacct.txt`), the
  discovery domain tables, `/bigdata/stajichlab/shared/projects/Fungi_5k/input/`, `analysis/hydrophobin_truth/`.
- `sinfo -s`; `module avail` for ncbi-blast, MMseqs2, signalp.
- Wilson lower bounds with z = 1.959964 (the value in `calibration/intervals.py`).

## Facts that check out

- H0 counts. `awk '!/^#/ && $4=="Hydrophobin"{print $1}' FILE | sort -u | wc -l`: Af293 6, A1163 5,
  W72310 6, Cimm_RS 1, Bder_ER3 2, Calb 0, Cneo 0, Scer 0. These equal the plan's expected values. The
  discovery run used `--cut_ga` (`hydrophobin_discovery/run.sh`).
- Q4WUX0, Q4WBD4 and Q4WBC1 are each in `_workdir/sorting_hat/Afum_Af293_UniProt.faa` once. All three
  headers say "GPI anchored protein".
- All eight E1 files exist in `Fungi_5k/input/` as `<name>.proteins.fa`: B. bassiana ARSEF_2860 (9478
  proteins), F. velutipes 6-3 (13341), F. fulva Race5_Kim (13560), F. graminearum PH-1 (11193), P. expansum
  MD-8 (10624), P. ostreatus PC9 (10483), T. asperellum FT101 (10088), T. virens Gv29-8 (10835). No `*` in
  any sequence. IDs are of the form `F1BB8A46_000001-T1`, not UniProt or NCBI accessions.
- `ncbi-blast/2.14.0+`, `MMseqs2/17-b804f` and `signalp/6-gpu` exist as modules. Earlier SignalP jobs ran on
  `short_gpu` (job 29551710, 4 min 12 s; job 29551711, 2 min 29 s; `Afum_Af293_UniProt/sacct.txt`).
- `calibrate truth --module pfam_hydrophobin --call hydrophobin_domain` is allowed by the CLI: the call
  reads one module (`calibration/cli.py` on `per-call-status`, the `truth` branch of `run`). `--leakage
  partial` is a valid choice (`call_status.py:37`). With leakage other than `none` the cap is `smoke`.
- `_family_digest`, `_condition_table`, `pfam.MODULES` (`modules/pfam.py:22`) and `call_eligible`
  (`engine.py` on `per-call-status`) exist with the names the plan uses.
- The repeat values for H1b (0.435, 0.462, 0.987, 0.976) equal the call files in
  `_workdir/sorting_hat/{Scer_S288C,Calb_SC5314}/status/calls/tandem_repeat_protein.json` and the repeat
  report section 8.
- No script in `scripts/sorting_hat/` uses `BASH_SOURCE`. All `.sbatch` files use `${SCRATCH:?}`.
- `per-call-status` and this branch merge without conflict (`git merge-tree`).

## Findings

### 1. BLOCKER: no step produces the inputs that `calibrate truth` needs

Evidence:
- `calibrate truth` reads `--calls-long` (a core run output) and refuses a run whose module identity is
  not the one in the work directory (`_check_run_identity`). With `--call-status` it also needs `run.json`
  next to it with `config_sha256` equal to the config, and every read module `ok` (`_check_call_run`).
- No task runs the core command (`cellsurface_sorting_hat ... --workdir --out --taxon`). H5 runs hmmsearch
  and wrappers only. H5b says "the full chain". H1b says only "re-run `calibrate truth --call-status`".
  After H1 the config hash changes, so H1b fails with "the run used another config" unless the core
  command runs first on S288C and *C. albicans*.
- `--truth` must be a TSV with columns `id`, `label` (1 or 0) and `cluster` (`_read_truth`). The plan's
  truth file `truth_all.tsv` has `accession, species, tier, pmid, evidence_code, sequence_sha, source`.
  `id` must be the proteome ID (for E1, IDs such as `F1BB8A46_000001-T1`). No task writes a per-proteome
  file in this format. Negatives (H2b step 4) get no cluster; `cluster` must not be empty.
- The E1 core runs need `--taxon` or `--taxon-map` (`cli.py:67-68` on `per-call-status`). The plan names no
  taxon ID for any E1 proteome.

Fix:
- Add a step after H5 and H5b: run the core command on each work directory (in-scope and E1) with the
  config of that moment, with the taxon ID named per proteome.
- In H1b: run the core command on S288C and *C. albicans* first, then `calibrate truth --call-status`
  with the same truth files (`truth_v3.*.tsv`), calibration sets and leakage (`tuned_on_truth`) as the
  existing files. Write the full commands in the plan.
- Add a step (end of H2b) that writes `truth.<proteome>.tsv` with `id, label, cluster` per work
  directory, for the test part of the split (finding 9). Cluster positives and negatives together per
  proteome, as `build_truth.py` of the repeat work does.
- In H6 step 5, write the exact commands: `--module pfam_hydrophobin --call hydrophobin_domain` for the
  module status, and `--call-status --call hydrophobin_protein` only if the rescue is kept.

### 2. MAJOR: the branch facts in the dependency map are wrong

Evidence:
- Plan line 32: `per-call-status` "already contains PR #75's commits through the shared base `484247e`".
  `git merge-base HEAD per-call-status` is `9a3202c` (PR #73 merge). `484247e` (PR #74 merge) is not an
  ancestor of `per-call-status`. `signoff-cfem-hydrophobin` (PR #75) is not an ancestor of
  `per-call-status`. It is an ancestor of this branch.
- Local `main` is at `6d67077` (PR #71). `origin/main` is at `484247e`. "Rebased on `main`" with the
  local branch would rebase on an old commit.
- H0 commits on `signoff-cfem-hydrophobin`. PR #75 changes only if that branch is pushed. Global constraint
  5 forbids a push without the owner. The H0 commit is also not on this branch, which already contains
  the old PR #75 commits, and H1 moves the PF01185 row.
- `gh pr list`: #75 and #76 are OPEN.

Fix: correct line 32. Say "rebase on `origin/main` after `git fetch`". Or say "merge `per-call-status`
into this branch" (`git merge-tree` shows no conflict). For H0, say that the owner pushes, and that
this branch takes the H0 commit (merge or cherry-pick) before H1.

### 3. MAJOR: `cys8_pattern` with state `error` for invalid sequences makes the module `partial`

Evidence:
- Plan H4 step 1: "invalid sequence gives state `error`".
- `_write` (`modules/cli.py`) sets the run state to `partial` if any row is `error`.
  `_check_call_run` refuses a read module that is not `ok`.
- Module convention: an invalid protein gets `invalid_row(p)`, state `na_invalid` (`modules/base.py:46-47`).
  `_write` does not count `na_invalid` as an error.
- `Scer_S288C/modules/step1_rule@R0.tsv.gz` has 8 rows `na_invalid`. With the plan's rule, `cys8_pattern`
  would be `partial` in S288C.

Fix: use `invalid_row` (`na_invalid`) for invalid sequences. Add a test: a FASTA with one invalid protein
gives run state `ok`. Review 2 finding 5 asked for this.

### 4. MAJOR: rule 6.3 counts rescue-only calls "once per cluster", but no step clusters them

Evidence:
- Spec 6.3: each rescue-only call is "counted once per cluster (so three *A. fumigatus* strains count
  once)".
- H2 step 3 clusters truth proteins only (`clusters.tsv`). Rescue-only calls are proteome proteins, mostly
  not in the truth set (spec F7: none of the in-scope regex candidates has a label).
- No step defines the cluster of a rescue-only call, or the outcome of a cluster whose members get
  different review outcomes.

Fix: add to H6 a step that clusters all rescue-only calls of all measured proteomes with the H2 MMseqs2
settings, writes `rescue_clusters.tsv`, and gives one outcome per cluster. State the rule for mixed
clusters (for example: `unresolved` unless all resolved members agree).

### 5. MAJOR: the rule 6.3 test cases do not test the rule

Evidence (Wilson lower bound, z = 1.959964):
- 4/4: 0.510. 5/5: 0.566. 5 with 4 correct: 0.376. 6 with 5: 0.437. 8 with 7: 0.529.
- Plan H6 step 1: "4 resolved clusters gives 'no evidence'". Only 4/4 tests the 5-cluster floor: with
  3/4 the Wilson bound alone already rejects.
- "5 clusters with 4 correct gives 'keep' only if the Wilson lower bound is at least 0.5" is not a fixed
  expectation. The bound is 0.376, so the answer is "no evidence".
- No case gives "keep".

Fix: fixed cases with fixed answers: 4/4 gives "no evidence" (floor); 5/5 gives "keep"; 4/5 gives "no
evidence"; 0 resolved gives "no evidence"; `unresolved` outcomes do not enter n. Use
`calibration/intervals.wilson`.

### 6. MAJOR: the "seen" list is incomplete

Evidence:
- Spec 5.3 item 2: the list "holds every protein named in the review files and in this spec".
- Plan H2 step 5 lists the 9 entries, RodA to RodG, PSH_FLAVE, CFTH1, FBH1 and Q4WUX0, Q4WBD4, Q4WBC1.
- Review 2 (finding 1 and 11) names more: B0YAH7, B0Y4E8 (A1163), KAK9640788.1 (W72310),
  `F00FD2C2_006465-T1` and `F00FD2C2_006466-T1` (Bder ER3). Review 2 counts 20 screening-regex candidates in
  the in-scope proteomes. All were seen. They are the likely rescue-only calls of rule 6.3.

Fix: build `seen.tsv` by script from (a) the IDs in both spec reviews and the spec, and (b) every
candidate in the output of `analysis/sorting_hat_run/hydrophobin_discovery/analyse.py`. Store the
command. Add "every seen protein" to the H6 step 3 review list (spec 6.2).

### 7. MAJOR: the T2 rule can label a non-hydrophobin as positive

Evidence:
- `sp_hydrophobin_query.tsv` has 189 rows. 15 do not have "hydrophobin" or "rodlet" in the name (spec F2).
  Row 2 is BRLA_EMENI, the conidiation regulator BrlA.
- Plan H2 step 2 fetches UniProt JSON "for the accessions in `sp_hydrophobin_query.tsv`". The H2 test says
  "an entry with ECO:0000269 on function is T2". An experimentally supported function of BrlA would then
  make BrlA a T2 positive.

Fix: restrict H2 to the 174 named entries (or state the filter). Require that the ECO:0000269 evidence is
on a statement that the protein is a hydrophobin (name, family or function text). Add a test: an entry with
ECO:0000269 on an unrelated function is not a positive.

### 8. MAJOR: HPC steps are not sized, and some cannot run as written

Evidence:
- Per-proteome job times in `Afum_Af293_UniProt/sacct.txt` and others: hmmsearch 18 to 47 s, repeats 1.6 to
  2.7 min, BLAST 41 to 60 s, TMHMM 9 to 10 min, SignalP 2.5 to 4.2 min. `submit_modules.sh` submits 5 jobs
  per proteome. For 8 E1 plus 8 in-scope re-runs, that is many jobs of a few minutes. The user rule is
  about 1 to 1.5 hours per job.
- `signalp_gpu.sbatch` has `#SBATCH -p exfab`. The plan says `short_gpu`. Earlier runs used `short_gpu`,
  so the partition was given at submit time. The plan does not say so.
- The `pfam` module command refuses to run without `--tm-module` while an active family has
  `second_condition=no_tm` (`modules/cli.py`, `run`). CFEM (PF05730) is active with `no_tm`
  (`family_table.tsv` line 2). E1 needs TMHMM. H5b names SignalP only.
- H5 step 1 writes a new `hydrophobin_runs.sbatch`. `pfam_hmmsearch.sbatch` already searches the whole
  family table with `--cut_ga`, uses `$SCRATCH`, checks the model names and writes `provenance.json`.
- H5 step 2 says "the wrappers `pfam_hydrophobin`, `pfam_hsba`". There is one `pfam` subcommand. It writes
  every module in `pfam.MODULES`, including `pfam_adhesion` and `pfam_allergen`.

Fix: one GPU job for SignalP over all 8 E1 proteomes (about 8 x 4 min), one CPU job for TMHMM over all 8
(about 80 min), one CPU job for hmmsearch over all 16 proteomes. Reuse the existing `.sbatch` logic in a
loop, or write one loop script with a test in `test_data_and_scripts.py` like the existing script tests.
Name the module commands for E1: `signalp`, `tm`, `pfam` (with `--sp-module` and `--tm-module`),
`cys8_pattern`. Add proteome provenance records for E1 as for the in-scope runs.

### 9. MAJOR: the split is made but never used, and the spec's test-part rule is missing

Evidence:
- Spec 4.4 item 2: "If the test part has fewer than 10 positive clusters, the report says the split gives
  no usable test". Spec 6.3: "The status files record the measurement on the test part."
- Plan H2 step 6 uses a different rule ("fewer than 20 per group"; "group" is not defined).
- Plan H6 does not say that status files use the test part only. No task uses the development part.

Fix: in H2, apply the spec rule (fewer than 10 positive test clusters). In H6, measure status files on the
test part; report all-data numbers only in report tables. State that the development part is unused if
the pattern is not tuned (H3 freezes it from the literature).

### 10. MINOR: test locations and test commands

Evidence:
- There is no `tests/analysis/` directory. Analysis tests live in `tests/<analysis directory>/` and load the
  script with `importlib` (`tests/calibration_truth/test_adjudicate_mechanism_labels.py` lines 8-14).
  `pyproject.toml` sets `testpaths = ["tests"]`.
- The Global constraint 1 command runs `tests/cellsurface_sorting_hat` only. It does not run new analysis
  tests.
- The golden files and `test_run_output_is_pinned` exist only on `per-call-status`. The regenerate command
  is `UPDATE_GOLDEN=1 pytest tests/cellsurface_sorting_hat/test_cli.py -k pinned` (`test_cli.py` docstring).
- (Line numbers in this item are from `per-call-status`.) The toy fixture of `test_cli.py` writes `pfam_adhesion` and `pfam_allergen` tables only (lines 42-43). If
  `hydrophobin_domain` and `hsba_domain` become mechanisms and the fixture has no table for them, they are
  `not_assessable`, and `_eval_other` adds them to `other_basis`. The assertions for ENZ1 and STAR1
  (`other_basis == "cocci_specificity_rank_top15"`, lines 143 and 146) then fail.
- The Pfam routing test (hit goes to `pfam_hydrophobin`, not `pfam_adhesion`) belongs in
  `tests/cellsurface_sorting_hat/modules/test_pfam.py` or `test_module_cli.py`. H1 does not list them.

Fix: use `tests/hydrophobin_truth/`. Add it to the test command. Name the regenerate command. Say whether
the toy fixture gets `pfam_hydrophobin` and `pfam_hsba` tables (then `other_basis` stays) or the
assertions change. List the module test files.

### 11. MINOR: H1c and H4 do not need PR #76

Evidence: `modules/cli.py`, `modules/pfam.py` and `modules/base.py` are identical on this branch and on
`per-call-status`. `_condition_table` is on this branch. H1c and H4 touch only these files and new files.

Fix: optional. H1c and H4 can move to Phase A. H1 still needs the golden files and `outputs.py` of #76.

### 12. MINOR: spec items with no task

- Cross-species hard-negative table (spec 4.3 item 2).
- Leave-one-species-out table at sequence level (spec 7 item 6).
- Curated entries from `data/curated/` as a truth source (spec 4.2 item 3); H2 uses them only for negatives.
- Where the 2c display combination (`hydrophobin_domain` OR `hsba_domain`) is made (spec 3.2; review 2
  finding 7).
- Seidl-Seiboth 2011 (spec 5.1) is not in H3. Peñas 1998 full text (spec 8) is not in H3 or H3b.
- `analysis/hydrophobin_truth/README.md` lines 6 to 8 are still pasted HTTP headers (review 2 finding 10).
- Sensitivity of the pattern on the 165 Pfam-positive entries (spec 6.1) needs the pattern without the R0
  gate. H4 does not say that the pattern function is callable on its own.

Fix: add each to a task, or list it as dropped with the reason.

### 13. MINOR: H3 has no source plan for the class I and II spacing

Evidence: `docs/paper/` holds no PDF. No file under the repository (outside `_workdir`) matches
Wessels, Linder, Sunde, Kubicek or Yang by name. The plan gives no PMID for Wessels 1994, Linder 2005 or
Sunde 2008. I did not check their access.

Fix: give the PMIDs. Add a stop point: if no class I or class II spacing has a read source, stop before
H4 and ask the owner.

### 14. MINOR: details in H2b and H5

- H2b step 1 says "two truth proteins" and then lists three.
- "95% identity over at least 90% of the length": the plan does not say query or subject length. The E1
  proteomes are re-annotations (IDs `F...-T1`), so model lengths can differ from UniProt.
- `_workdir/sorting_hat/*.faa` matches 11 files, including `Afum_Af293_Fungi5k.faa`,
  `Scer_S288C_nodubious.faa` and `iuis_fungal_allergens.faa`. Name the eight files.
- The three *A. fumigatus* strains are three work directories. The plan does not say which get a
  `pfam_hydrophobin` status entry, and the 7 Rod proteins are in each.

Fix: correct the count, name the length, list the eight files, and state the strain rule for status
entries.

### 15. MINOR: acceptance checks that a command cannot run

- H3 ("every number has a source line"), H3b and H7 ("docs match the report") have no command.
- H5's new script has no test. `test_data_and_scripts.py` has tests for the other job scripts (for
  example lines 190-225).
- H6 should reuse `calibration/intervals.wilson` and `cluster_bootstrap`, not new code.

Fix: for H3, a small schema check of `cys8_spacing.yaml` (every gap has `min`, `max`, `source`). For H5,
a script test. For H7, a grep that the report numbers appear in the docs.

## Mapping of review 2 findings to the plan

| Review 2 finding | In the plan |
|---|---|
| 1 rules 6.2 and 6.3; strains once | Partly. T4 outcome in H6 step 3. Cluster counting has no clustering step (finding 4). |
| 2 R0 identity in params | Yes. H4 steps 2 and 3. |
| 3 order of config edits | Yes. H4 does not edit `categories.yaml`; H1b is last. H1b lacks the core re-run (finding 1). |
| 4 E1 coverage and strain mapping | Yes. F. velutipes is in E1; H2b maps by BLAST. |
| 5 row states | No. The plan uses `error` (finding 3). |
| 6 mechanism placement | Yes. H1 step 3. |
| 7 `hsba_domain` role and 2c display | Role yes. Display location no (finding 12). |
| 8 sources for H3 | Partly (finding 13). |
| 9 H1 dependency wording | Yes, but with a wrong branch fact (finding 2). |
| 10 H0 commit item, README lines | H0 commit item removed. README lines not fixed (finding 12). |
| 11 "seen" list | No (finding 6). |
| 12 `estimated` rule | Spec fixed. Not needed in the plan. |
| 13 keep decision on all data | Spec fixed. H6 step 5 cites spec 5.3. |

## Not checked

- Whether the H3 papers are open access, and what spacing they give.
- Whether any truth sequence maps to an E1 proteome.
- The count of reviewed secreted UniProt proteins per E1 species (the negative population size).
- The full test suite on `per-call-status` (the handoff reports 788 passed, 7 skipped; I did not run it).

## Verdict

**Not ready to execute.** One blocker: no step makes the core-run output and the `id, label, cluster`
truth files that `calibrate truth` needs, so H1b and H6 cannot run as written (finding 1). Eight major
findings need plan edits: branch facts (2), invalid-row state (3), rescue clustering (4), the rule 6.3
test cases (5), the "seen" list (6), the T2 filter (7), HPC sizing and missing TMHMM for E1 (8), and use
of the split (9). Phase A tasks H0 (after the push question in finding 2), H3 and H3b can start once
findings 2 and 13 are fixed. These are text changes to the plan.
