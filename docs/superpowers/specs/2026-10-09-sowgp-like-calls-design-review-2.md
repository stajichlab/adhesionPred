# Review 2 of the M5 spec (SOWgp-like calls, revision 2)

*2026-10-09. Independent second-round review of `docs/superpowers/specs/2026-10-09-sowgp-like-calls-design.md` (rev 2) against `2026-10-09-sowgp-like-calls-design-review-1.md`. Reviewer: Claude Code (claude-opus-5-5). No spec or code file was edited. Each statement below was checked by reading the named file or by a recount with `/usr/bin/python3.12` on the named table.*

## 1. Verdict

**PASS** (no blocker). Both first-round blockers are resolved. Three new MAJOR findings must be fixed in the spec before the plan is written. All three are text fixes.

Counts: 0 BLOCKER, 3 MAJOR, 8 MINOR.

## 2. Findings

### M1. MAJOR: entries without positives give status `unvalidated`, not `smoke`

Evidence.
- `status.py` lines 120 to 122: `if not sens or measure.get("n_pos", 0) < 1: return UNVALIDATED`. `cluster_bootstrap` returns `sensitivity: None` when there are no positives (`calibration/intervals.py` lines 104 to 108).
- Spec section 4 has five entries with no positives: `sowgp_ortholog` in S288C, *C. albicans* and *B. dermatitidis* ER-3, and `pro_cys_array_protein` in S288C and *C. albicans*. Each of these gets `unvalidated`.
- Spec section 4 says "both calls are capped at `smoke`". Work plan M5c has the check "status files `smoke`". That check fails for these five entries.

Fix. State that specificity-only entries get `unvalidated` and that their specificity is kept in the `measure` block only. Change the M5c check to "`smoke` for RS and *C. posadasii*; `unvalidated` with a recorded specificity for the specificity-only entries". Or drop the specificity-only entries and give their specificity in the report.

### M2. MAJOR: the `categories.yaml` edit makes every existing call status file stale; the plan has no re-measurement

Evidence.
- `call_status.py` lines 177 to 180: `stale_reason` returns "config differs" when `source.config_sha256` differs from the run's config hash. `call_status.py` line 8: "A stale file is never used."
- M5b2 adds two calls, two thresholds and mechanism-list entries to `categories.yaml`. The config hash changes. The `tandem_repeat_protein` call status files (`docs/reports/data/sorting_hat/call_status/{Scer_S288C,Calb_SC5314}/tandem_repeat_protein.owner_definition.config-with-surface-attachment.json`, the current ones per the M2 spec line 41) become stale.
- The M2 spec line 131 says "M5 ... require[s] the repeat re-measurement again". The M5 spec line 48 says a threshold change "changes `config_sha256` only", and M5b2's check is "config hash change noted". Neither says that the existing statuses stop applying or that they must be measured again.

Fix. In M5b2 (or a new task after it) re-run `calibrate truth --call-status` for `tandem_repeat_protein` in S288C and *C. albicans* with the new config, using `analysis/calibration_truth/repeat_call_truth/truth_v3.{Scer_S288C,Calb_SC5314}.tsv`. Rewrite line 48: "any `categories.yaml` change changes `config_sha256`, which makes every call status file stale".

### M3. MAJOR: the `sowgp_unit` module as written makes every non-hit protein `not_assessable`; the `sowgp_ortholog` expression is not given

Evidence.
- Spec section 3.1: "one row per protein with at least one row in the `hmmsearch --cut_ga` table: `hit` (1), ...". Read literally, proteins without a hit get no row.
- `modules/base.py` `write_module` (lines 50 to 81) writes only the rows it is given; it does not fill missing proteins.
- `engine.py` `_row` (lines 181 to 188) returns `None` for a missing row. A `flag` or `call` leaf then gives `NOT_ASSESSABLE` (lines 196 to 201). So `sowgp_ortholog` would never be `not_called`, and section 4's negatives ("all other RS proteins") would all be `not_assessable`. `calibrate truth` counts them in "not assessable", not as negatives (`calibration/cli.py` lines 475 to 483).
- The existing pattern writes one row per FASTA protein with `hit` 0 when there is no table row (`modules/hydrophobin_relaxed.py` `relaxed_rows`, lines 72 to 100).
- Section 3.2 gives exact YAML for `pro_cys_array_protein`. Section 3.1 gives none for `sowgp_ortholog`.

Fix. Write "one row per protein in the FASTA; `hit` 1 and `n_units` the row count when the protein has at least one `--cut_ga` row, else `hit` 0 and `n_units` 0; invalid sequences get `invalid_row`". Give the YAML, for example `- name: sowgp_ortholog` / `expr: {flag: sowgp_unit.hit}`, and confirm that it reads one module so `calibrate truth --module sowgp_unit` applies. In the M5a check, write "0 proteins with `hit` 1" for S288C and ER-3, not "0 rows".

### m1. MINOR: the BAD1 control row repeats "not verified", but review 1 verified that it was not run

Evidence. `analysis/cocci_repeats/unit_bad1_like.tsv` does not exist. `analysis/cocci_repeats/logs/` has no log for script 38 or 39 (checked again). Spec section 2 row "BAD1 control" still says "whether it was run is not verified".

Fix. Write "not run: `unit_bad1_like.tsv` absent, no log".

### m2. MINOR: the roadmap still expects the *B. dermatitidis* run to call BAD1 under (b)

Evidence. `docs/ROADMAP-2026-10-08-DRAFT.md` line 335 (M5 row): "*B. dermatitidis* run calls BAD1 under (b) only". Spec 3.2 and O3 default: BAD1 fails the rule and is not covered.

Fix. Say in the spec that the roadmap M5 check changes with O3, or update the roadmap in M5e.

### m3. MINOR: `36_unit_calibrate.py` has no bit-score floor

Evidence. The script has only `--floor-e` (lines 148 to 152); the found test is `e <= floor_e` (line 237). It also reads `fullscore` from `f[5]` (line 109), which is `qlen`. Spec 3.1 and M5a say "rerun `36_unit_calibrate.py` with a bit-score floor of 26".

Fix. M5a: add a `--floor-bits` option (domain score) and fix or remove the `fullscore` field. Then rerun.

### m4. MINOR: the label tier of the ortholog positives is not stated, and the RS entry can hold one protein only

Evidence. Section 4, row 1: positives are "PMID 12065484 for SOWgp itself; other orthologs from `05_sowgp_pangenome.py`". An RS entry can only match RS proteins; other strains' orthologs are "truth rows without a call" (`calibration/cli.py` lines 475 to 478). Orthologs found by `05_sowgp_pangenome.py` are labels from sequence similarity. Under the hydrophobin tier table cited in section 4 (`2026-10-08-hydrophobin-validation-design.md` line 145), such a label is T3 ("secondary set ..., never in the primary result"). Whether PMID 12065484 names the *C. posadasii* gene directly was not checked in this review.

Fix. RS row: positive = `XP_001245172.2` only. Silveira row: state the tier of `QVM09276.1` (T1 only if a paper names that gene; else T3) and the source.

### m5. MINOR: the `tuned_on_truth` label for the RS positive is verified, not only plausible

Evidence. `docs/reports/data/sorting_hat/panel.tsv` line 2: `XP_001245172.2`, `tandem_repeat_protein`, tuning `tuned`. The `pro_cys_array_protein` call reads the same `repeat02`/`repeat14` modules. The source of the 15% and 4% thresholds is still unrecorded, as the spec says.

Fix. Section 4 row 4: "`tuned_on_truth` (SOWgp tuned the repeat modules, panel.tsv; the source of the composition thresholds is not recorded)".

### m6. MINOR: name the truth files for the specificity-only `pro_cys_array_protein` entries

Evidence. Section 4 row 5 says "repeat-truth proteins that are not Pro/Cys-rich arrays". The files are `analysis/calibration_truth/repeat_call_truth/truth_v3.Scer_S288C.tsv` and `truth_v3.Calb_SC5314.tsv` (commit 5b27195). The spec does not name them or say how a label 1 there becomes label 0 here.

Fix. Name the files and the relabel rule.

### m7. MINOR: YAML order constraints are not stated

Evidence. `engine.py` lines 117 to 126: a `mechanism` entry must be an earlier ungated call; a `basis_calls` entry must be an earlier ungated call. `sowgp_ortholog` must precede `other_not_surface` (`categories.yaml` line 95); `pro_cys_array_protein` must precede `surface_attachment_candidate` (line 69) if it joins its `basis_calls`.

Fix. One sentence in 3.1/3.2 or in M5b2.

### m8. MINOR: line references and run list

Evidence. Section 4 cites `cli.py` lines 461 to 495 for the one-taxon rule. The one-taxon check is at `calibration/cli.py` lines 449 to 453; the truth-taxon check call is at line 489. Section 4 lists *C. albicans* for `sowgp_ortholog`, but M5a runs only RS, *C. posadasii*, *B. dermatitidis* and S288C.

Fix. Cite `calibration/cli.py` lines 442 to 495. Add *C. albicans* (`_workdir/sorting_hat/Calb_SC5314.faa` exists) to the M5a run list, or drop it from section 4.

## 3. First-round findings: status

| # | Review 1 finding | Status in rev 2 |
|---|---|---|
| 1 | BLOCKER BAD1 fails rule | Resolved. 3.2 "Not covered"; O3 default drop BAD1. Roadmap still disagrees (m2). |
| 2 | BLOCKER positive set | Resolved. Positives are external labels only; the 20 and 74 are T4 in 3.3 with cluster counts. Tier of ortholog positives still open (m4). |
| 3 | Sensitivity floor | Resolved in text (section 2, 3.1). The rerun needs a script change (m3). |
| 4 | Domain vs sequence score | Resolved (section 2 row 2, risk 6). |
| 5 | Call YAML | Resolved. YAML matches the node grammar (`categories.yaml` lines 4 to 9; `engine.py` lines 151 to 163, 192 to 215); no `ref`; reads three modules; `call_eligible` (lines 452 to 475) accepts it. |
| 6 | Detector choice | Resolved (O2b, M5b). |
| 7 | Name vs rule | Resolved by stating the meaning; rename left to O2. Accepted option of review 1. |
| 8 | Truth per species | Mostly resolved: one row per species, ID mapping, *C. posadasii* proteome, Fungi_5k moved to report. New issues M1, m4, m6. |
| 9 | GA design | Resolved (`GA -1000.00 26.00;`, decision from `--cut_ga` rows, sequence score reported only). Module row coverage is new issue M3. |
| 10 | Cutoff inside a gap | Resolved (3.1, risk 2). |
| 11 | Onygenales set | Resolved (section 2 row 4). |
| 12 | Panel path | Resolved. |
| 13 | Name and gate vs M2 | Resolved. M2 spec line 55 now names `sowgp_ortholog` and `pro_cys_array_protein`; M5 uses one category. |
| 14 | `other_*` and split models | Resolved (3.1 Composites, O1, risk 4). |
| 15 | Adhesion claim | Resolved (section 1 and 3 cite PMID 12065484). |
| 16 | Task order | Resolved (stop point; M5b2 separate). |
| 17 | Leakage labels | Resolved (section 2 last row, section 4). See m5. |

## 4. Items checked and correct

Recounts on the analysis tables (all match the spec):
- `unit_genomes.tsv`: 6,313 rows (fungi5k 5,813; cocci_pangenome 493; cocci_longread 7). Label `Coccidioides_immitis_RS` occurs twice.
- `unit_hmm_hits.tsv.gz`: 2,055 rows; `full_score` is 47 in every row. *Coccidioides* 1,748 domain hits; 1,293 at 26 or more in 487 proteins and 455 proteomes; lowest of these 43.9. Others: 307 rows, all `dom_of` 1, max 19.7. *Coccidioides* below 26: 455 hits, 18.3 to 25.8. 487 *Coccidioides* proteins with any hit.
- `unit_calibration_sensitivity.tsv` (natural flank): 0.6 all found; 0.5 1.000/0.953; 0.4 0.923/0.057, median best 26.0; 0.3 0.013. `36_unit_calibrate.py` default `--copies 4`, `--floor-e 1.0`.
- `unit_architecture_onygenales.tsv.gz`: 644,399 rows, 78 strains (RS and Silveira each twice); 1,457 repeat calls; 83 pass Pro+Cys > 15 and Cys >= 4, in 50 proteomes; 9 at period 47 in 9 proteomes; 74 remaining in 27 periods. `37_unit_search.sh` line 36 selects Fungi_5k Onygenales plus `cocci_longread`.
- `class2a_candidates.tsv`: 41 rows; 20 / 10 / 11; 12 of 20 pass the rule; 4 period-47 SOWgp orthologs as named; Ser/Thr rows in CiB10637 1, CiB10992 1, VFC140 2, Cpos1038 4, Cpos3700 1, Silveira 1.
- RS raw repeat tables (`_workdir/sorting_hat/Cimm_RS/raw/repeats/`): the rule gives 6 calls with repeat14 and 7 with either; `XP_001239868.2` is repeat02 only (Pro 22.3, Cys 5.0); `XP_001248770.2` (period 27, Pro 3.5, Cys 12.1) and `XP_001248220.2` (period 43, Pro 9.2, Cys 16.1) are called.
- `14_repeat_detect_general.py` `composition()` (line 448) gives rounded `pct_ser_thr`, `pct_pro`, `pct_cys` on the sequence passed in.
- `modules/repeats.py`: call = period > 0, coverage >= 0.25 (line 11), copies >= 2.5 (line 12), rule at line 88.

Code and docs:
- `engine.py` line 192 region: `{call: X}` reads module X's `call` column. `call_eligible` refuses `ref` (line 464) and needs two or more modules (line 473).
- `calibrate truth`: one taxon only; `--module` path refuses a call that reads more than one module; `--leakage` values include `in_reference` and `tuned_on_truth` (`call_status.py` line 37); any value except `none` caps at `smoke` (`calibration/cli.py` line 520).
- `status.py` lines 112 to 133: `estimated` needs 20 positives, 20 negatives and 20 clusters per class; the spec's "neither can reach `estimated`" holds.
- `hydrophobin_relaxed.model_info` and `check_table` accept any single GA pair and check `--cut_ga`, query file name and model length; they can be reused as the spec says.
- `2026-10-08-hydrophobin-validation-design.md` line 147: T4 "Never added to the truth set, never used for sensitivity."
- `docs/reports/data/sorting_hat/panel.tsv` lines 2 and 3: SOWgp rows as stated.
- `REPORT_2026-10-01_sowgp_depth_vs_repeats.md` line 88 (39 of 488) and `REPORT_2026-10-01_sowgp_longread_annotation.md` line 80 (3 of 6 loci split): risk 4 numbers hold.
- `_workdir/sorting_hat/` has no *C. posadasii* proteome (checked).
- Not re-run in this review (taken from review 1): the HMMER GA-line tests, the 1,000 synthetic arrays (264 proteins), the sequence-versus-domain counts, and the RS `XP_001245172.2` domain scores.
