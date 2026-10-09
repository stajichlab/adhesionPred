# Review 2 of the custom hydrophobin HMM design

*2026-10-08. Second independent review of `2026-10-08-hydrophobin-custom-hmm-design.md` (revision 2, commit 33718c5),
branch `hydrophobin-validation`. Reviewer: Claude Code (claude-opus-5-5). The spec was not edited. Nothing here is
committed. Every number below was computed by the reviewer with the commands in section 0, unless it is marked as quoted.*

## 0. Method

- HMMER 3.4 (`module load hmmer/3.4`), MMseqs2 17-b804f, Python 3.12. Scratch outputs are in the session scratchpad, not in the repository.
- (a) GA rewrite: `hmmfetch` of Hydrophobin and Eas from `/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm`
  (`current` is a symlink to `2026-01-27-Pfam38.2`). `GA` lines rewritten with `sed` to `8.50 8.50;` and to `8.50 0.00;`.
  `hmmsearch --cut_ga --domtblout --tblout` on `analysis/hydrophobin_truth/nopfam9.faa`. The table was read with
  `cellsurface_sorting_hat.modules.pfam.parse_domtblout`.
- (b) Clusters: the 131 T2 sequences of `truth_all.tsv`, `mmseqs easy-cluster --min-seq-id 0.3 -c 0.5`. The seven
  hydrophobin-class models from `_workdir/sorting_hat/calibration/hydrophobin_discovery/models.hmm`, searched with `--cut_ga`, with
  rewritten GA lines, with `-T 0 --domT 0`, with `--nobias`, and with `--max`.
- Cost: the seven models with `-T 0 --domT 0` (default filters, and again with `--nobias`) on the 12 proteomes of `run_list.tsv`.
  Joined to `pfam_hydrophobin`, `pfam_hsba`, `step1_rule@R0` and `cys8_pattern.n_cys` in each work directory, and to
  `proteome_map.tsv` for labels. Extra call = hit at the cutoff, not hit by `pfam_hydrophobin` or `pfam_hsba`, at least 8 Cys,
  R0 `called`.
- (c) LP: `jensen_resolution.tsv`, `jensen_sequences.faa`, `protein_lists.tsv` (`cys_pattern_stated`, `N{x}` read as `.{x}`).
  Pfam `--cut_ga` and `-T 0` on the 50 sequences. Exact-sequence match against T2/T3. A joint MMseqs2 clustering of the 7 LP
  Pfam misses with the 131 T2 sequences.
- (e) Code read: `engine.py` (`_validate`, `_eval_other`, `call_eligible`, `modules_of_call`), `categories.yaml`,
  `modules/pfam.py`, `modules/base.py`, `modules/cli.py`, `cli.py` (`absent_modules`), `scripts/sorting_hat/pfam_hmmsearch.sbatch`,
  status file locations.

## 1. Status of the review 1 findings

| # | Status | Evidence in revision 2 |
|---|---|---|
| B1 | RESOLVED (intent) | Relaxed Pfam is the baseline (E10, 4.1). The HMM criterion is relative (6.2). The new criterion has its own defects: N1, N2, N3. |
| M1 | RESOLVED | 6.2 counts clusters, gives the Wilson interval for 3 of 6 ([0.19, 0.81], reviewer-checked) and calls it a smoke test. The cluster list is only by reference to review 1 (N13). |
| M2 | RESOLVED (design) | Nested leave-one-cluster-out (6.1); alignment by script (4.2). The HMM cutoff is still tuned on in-sample scores (N4). |
| M3 | PARTIAL | Fold negatives are named (6.1). The proteome sample has no source, size or cluster rule, and no criterion uses its fold scores (N9). |
| M4 | RESOLVED | Cutoff chosen on training clusters and frozen (6.1). Cost measured once (6.3). False-positive loop moved to v2 (5, 8). |
| M5 | RESOLVED | 6.3 defines the count (no `pfam_hydrophobin`, no `pfam_hsba`, no T1/T2 label), per 10,000 proteins, in each of the 12 proteomes; labelled calls in a separate column; risk 9.5. The limit does not bind (N8). |
| M6 | PARTIAL | GA rewrite, job script, `ModuleSpec.artefacts`, query-file check: present. The score is not a full-sequence score as stated (N3). `params` does not list the minimum cysteine count (N10). |
| M7 | PARTIAL | Per-group rule, tuning/test split, three methods on the same negatives, HsbA question to the owner: present (5, 6.4). What a failed group does to the ship decision is not stated (N5). |
| M8 | RESOLVED | Task L4a. |
| M9 | NOT RESOLVED | The reserve rule (3.2.5) depends on v1 output, is used by v1 evaluation (6.7), and leaves about one protein (N7). The post-freeze T1 route in section 8 is fine. |
| M10 | PARTIAL | Rule 6.3 kept as 6.5. It counts the tuning proteins as precision evidence (N6). |
| M11 | RESOLVED | Section 1.2 and 6.8. |
| M12 | RESOLVED | Grouped models dropped (4.2). |
| m1 | RESOLVED | `partial` label (6.8); build by script (4.2). |
| m2 | RESOLVED | 28 clusters (2, 3.1, 9.1). Reviewer re-clustered the 131 T2: 28. |
| m3 | RESOLVED | Bit scores and the 9-sequence search are stated (1.3). The "no hit" for 4 proteins depends on the bias filter (N2). |
| m4 | RESOLVED | Six T2 entries listed (6.2). Reviewer-checked: 6, 6, 7, 7, 7, 7 Cys; all six have a Pfam GA domain hit. |
| m5 | RESOLVED | 6.6. |
| m6 | RESOLVED | L6b; "added", not "replaces" (4.1). The only call status files on disk are `tandem_repeat_protein` (`_workdir/sorting_hat/Calb_SC5314/status/calls`, `docs/reports/data/sorting_hat/call_status/`), so L6b covers them. |
| m7 | PARTIAL | L1b reviews T2 text. What happens to a flagged entry (dropped before L2 or kept) is not said (N13). |
| m8 | RESOLVED | Section 8. |
| m9 | RESOLVED | 3.2 rule 2. |

## 2. Facts checked

| Spec claim | Re-derived | Verdict |
|---|---|---|
| 3.1: T2 131 in 28 clusters | 131; 28 | correct |
| 6.2: 9 missed in 6 clusters (3, 2, 1, three singletons) | {HYD1_TRIAP, HFB3_HYPVG, HYD2_GIBZE} 3 of 3; {HCF1, HCF2 of HCF1-3} 2 of 3; {HFBE_PENEN with RODF_ASPFU} 1 of 2; HYD1_GIBZE, PSH_FLAVE, HCF4_FULFL singletons | correct |
| 1.3: bit scores 23.9, 23.6, 18.6, 15.9, 8.6; four with no hit | same, with default filters | correct with default filters; see N2 |
| 1.5: 8.5 bits reaches 5 of 9 | full-sequence score: 5 of 9, in 5 of 6 clusters. GA lines `8.50 8.50` as 4.1 prescribes: 4 of 9, in 4 of 6 clusters (HYD2_GIBZE domain score 7.1) | correct for full-sequence only; see N3 |
| 1.5: extra calls 0 to 3; *F. fulva* and *F. graminearum* 3 | same totals at 8.5 bits, full-sequence. Of the 10 calls, 4 are labelled (HCF2, HYD1_GIBZE, HYD2_GIBZE, HFBE_PENEN) and 6 are unlabelled | correct |
| 1.5: at 15 bits 0 or 1 extra; 3 of 6 proteome misses back | same | correct |
| 3.2: 50 resolved | 50 of 50 `resolved` | correct |
| 3.2: stated pattern found in 41 of 49 | 41 of 49 | correct |
| 3.2: 43 of 50 Pfam GA hits; the 7 named | 43; same 7 names | correct |
| 4.2: mafft 7.505, famsa 2.4.1, no clustalo | `module avail`: mafft 7.490, 7.505; famsa 2.4.1 | correct |

Check (a): the GA rewrite works. `hmmsearch --cut_ga` on a file with `GA 8.50 8.50;` exits 0. The header has
`# Option settings: ... --cut_ga ...`, `# Query file: relaxed.hmm` and `# [ok]`. `parse_domtblout` accepts the table and
returns 5 domain rows (PSH_FLAVE, HFBE_PENEN x2, HCF2_FULFL, HYD1_GIBZE). The parser reads domain rows and field 14 (domain score) only.

Check (e): `hydrophobin_extended` as written reads two modules, has no `ref`, is not `kind: other`, and names no step 1 module.
`call_eligible` returns true. `_validate` accepts it as a mechanism if it is placed before the `other_*` calls. The R0 dependency
sits inside the module (as `cys8_pattern`), so `modules_of_call` does not list R0. That is consistent with the earlier spec 3.4.

## 3. New findings

### N1. BLOCKER. The HMM decision rule is not defined, and in the readings that can be computed it is impossible or needs 100%

Evidence.
- Pfam relaxed (default filters) can score 5 of the 9 proteins at any cutoff. HYD1_TRIAP, HFB3_HYPVG, HCF1_FULFL and HCF4_FULFL have no
  domain row even at `-T 0 --domT 0`.
- "A cluster is recovered" is not defined. Two of the 6 clusters have both reachable and unreachable missed members
  ({HCF1, HCF2}; {HYD1_TRIAP, HFB3_HYPVG, HYD2_GIBZE}).
- The fold cutoff "maximises recall of the Pfam-missed training proteins". Many cutoffs give the same recall. No tie-break is given.
- Leave-one-cluster-out results for relaxed Pfam (default filters, full-sequence score, R0 assumed `called`, all 9 have >= 8 Cys):
  - "highest cutoff with maximal training recall", "any missed member": 4 of 6 (in the HYD1_TRIAP fold the training cutoff is 15.9, so HYD2_GIBZE at 8.6 is lost). 2 clusters remain.
  - "lowest cutoff the negative rate allows", "any member": up to 5 of 6. 1 cluster remains (HCF4).
  - "all missed members": 3 of 6. 3 clusters remain.
- The HMM must recover "at least 2 more" clusters. With 1 left it cannot pass. With 2 left it must recover both. The rule does not say
  whether "more" is a difference of counts or the clusters the HMM recovers that relaxed does not (the union gain that matters to an `or` call).
- R0 has not been run on the Swiss-Prot sequences (L4a), so these counts are upper bounds.

Fix. Define "cluster recovered" (for example: at least one Pfam-missed member called). Define the tie-break (for example: the highest
cutoff that reaches the maximal training recall). Define the gain as the Pfam-missed clusters called by the HMM in a fold and not by relaxed
in the same fold. Write the relaxed numbers above into the spec and set the HMM threshold against the clusters that relaxed leaves
(for example "at least 1 of the remaining clusters, and no cluster lost"). Or state that the HMM is expected to fail and why it is still built.
Decide N2 first, because it changes the remaining clusters.

### N2. MAJOR. The HMMER bias filter hides two "Pfam-missed" proteins. The search options are not part of the freeze

Evidence.
- `hmmsearch --cut_ga --nobias` (also `--max`, or `--F1 1 --F2 1 --F3 1`) on `nopfam9.faa`: HCF1_FULFL hits Eas at 30.4 bits
  (domain 30.0). The Eas GA is 27.0. With the default filters there is no row.
- With `--nobias`, HCF4_FULFL scores 16.6 (domain 14.7) to Eas. HYD1_TRIAP (4.3) and HFB3_HYPVG (2.4, domain 5.2) stay low.
- On the 131 T2 sequences, `--cut_ga --nobias` hits 123 (default 122).
- In the 12 proteomes, `--cut_ga --nobias` adds one GA hit to `pfam_hydrophobin`: FEBB419A_005025-T1 in *F. fulva*, which
  `proteome_map.tsv` maps to Q00367 (HCF1_FULFL). No other proteome changes.
- With `--nobias`, relaxed Pfam can reach all 6 clusters (HCF4 at <= 14.7 bits on the domain score, HYD2_GIBZE at <= 7.1). At 14.5 bits it reaches 5 of 6.
- `scripts/sorting_hat/pfam_hmmsearch.sbatch` uses default filters. The spec does not name the filter options for relaxed Pfam or for the HMM.
  6.1 freezes aligner, cutoff procedure, thresholds and negatives, not search options.

Fix. Name the `hmmsearch` options of each method in 4.1 and 4.2 and include them in the freeze, `params` and provenance. Report the
relaxed baseline with and without `--nobias`. Say in section 1 that one of the 9 is a filter artefact at GA. Hydrophobin
sequences are cysteine- and glycine-rich, so a composition filter is a plausible cause, but this review did not test the cause.
Decide whether the HMM comparison uses the same options as the baseline (it should).

### N3. MAJOR. "One frozen full-sequence bit score" with equal sequence and domain GA is not a full-sequence cutoff

Evidence.
- With `--cut_ga`, a sequence is kept if its full-sequence score is at least GA1 and a domain row is written only if the domain score
  is at least GA2. `parse_domtblout` reads domain rows only.
- GA `8.50 8.50`: HYD2_GIBZE (sequence 8.6, domain 7.1) is in the `--tblout` but not in the `--domtblout`. 4 of 9 proteins, 4 of 6 clusters.
- GA `8.50 0.00`: 5 of 9, as in section 1.5.
- In the proteomes, relaxed-only unlabelled calls at 7 bits with `--nobias` differ by method: *F. graminearum* 5 (equal GA) and 6 (sequence only);
  *P. ostreatus* 3 and 4.

Fix. Choose one rule and use it in every number. Either set the domain GA to a fixed low value and call the cutoff "full-sequence",
or keep equal values and call the rule "sequence and best domain both at least T". Re-state 1.5 with that rule.

### N4. MAJOR. The HMM cutoff is chosen on scores of its own training proteins

Evidence. In a fold, the HMM is built from all training clusters. These include the Pfam-missed training proteins. "Maximise recall of
the Pfam-missed training proteins" then reads in-sample scores. Recall will be near the maximum over a wide range of cutoffs, so the
negative-rate limit or the tie-break chooses the cutoff. Relaxed Pfam does not have this problem (its models are fixed). The two methods
are tuned in different ways.

Fix. Score the training proteins out of sample inside each fold (an inner leave-one-cluster-out over the 27 training clusters;
`hmmbuild` is fast). Or set the HMM cutoff from the training negatives only (for example the highest negative score plus a fixed
margin), and use the same rule for relaxed Pfam.

### N5. MAJOR. No single ship rule joins 6.2 to 6.5

Evidence.
- 6.2 (recall), 6.3 (cost), 6.4 (per-group hard negatives) and 6.5 (precision) are separate tests. The spec does not say which must all
  pass for relaxed Pfam or for the HMM to enter `hydrophobin_extended`.
- 6.4 says when a group "passes". It does not say what a failed group does (M7 partial).
- 6.2 says relaxed passes "at the cost limit of 6.3" in the leave-one-cluster-out runs. 6.3 measures cost once, at the frozen cutoff.
  The fold cutoffs have no cost number.
- If relaxed Pfam fails, E10 says "the relaxed Pfam level stays" only for the case that the HMM adds nothing. The case "relaxed fails" has no outcome.

Fix. Add one table: relaxed ships if A and B and C and D; the HMM is added if relaxed ships and E and F. Give the outcome when relaxed fails
(no `hydrophobin_extended` call; or the HMM is tested against Pfam GA).

### N6. MAJOR. The precision term counts the tuning proteins, and leans on T5

Evidence.
- At 8.5 bits (default filters, full-sequence) there are 10 relaxed-only proteome calls with the conditions. 4 are T2-labelled:
  HCF2_FULFL, HYD1_GIBZE, HYD2_GIBZE, HFBE_PENEN (4 clusters). These are the Pfam-missed proteins that the cutoff is chosen to recover (6.1).
  6 are unlabelled: C6_00820W_A (*C. albicans*, 11 Cys), F00FD2C2_006465-T1 (*B. dermatitidis*), FEBB419A_011564-T1 and FEBB419A_012043-T1
  (*F. fulva*), F0349401_010222-T1 (*F. graminearum*), F4DD442B_006382-T1 (*P. ostreatus*).
- With only the 4 labelled clusters, n = 4 < 5, so 6.5 fails. With 10 resolved clusters, the Wilson lower bound is 0.49 for 8 of 10
  and 0.60 for 9 of 10. So the gate needs owner (T5) decisions, and 9 of 10 must be hydrophobins.
- E7 says T5 is "reported beside the experimental tiers and never merged". 6.5 lets T5 resolve clusters for the keep rule. That is a merge
  into a decision. Rev 3 rule 6.3 used T4 in the same way, so this is a known choice, but the spec does not say it.

Fix. Exclude proteins used to choose the cutoff from 6.5, or report the gate with and without them. State that the gate then rests on
T5, and that E7 allows T5 in a decision rule but not in a status.

### N7. MAJOR. The v2 literature reserve is circular and almost empty (M9 not resolved)

Evidence.
- 3.2.5 reserves the set "before v1 training" but defines it as the Jensen proteins "not recovered by v1". That needs v1 first.
- The Jensen Pfam misses are 7. ACLA_001890 is one of the 8 unverified, so 6 remain. AO090012000143 and AFLA_014260 are identical
  sequences, so 5 are unique.
- Pfam scores (default filters, domain): AO090012000143 19.7, ATEG_08089 18.3, AFLA_063080 17.3, AFLA_060780 15.8. ATEG_10285 has no row
  (6.2 bits with `--max`). Relaxed Pfam at <= 15.8 bits recovers all except ATEG_10285. The reserve is then 1 protein.
- 6.7 scores the same LP set in v1. So v1 "sees" it in evaluation, against E4 ("proteins that v1 never saw").
- In a joint clustering (30%, 0.5) ATEG_08089 and ACLA_001890 fall in one cluster with RODF_ASPFU, a T2 training protein. The LP test is
  not out of cluster for every protein.
- 14 of the 50 Jensen sequences are identical to T2 Swiss-Prot entries (for example AFUA_5G09580 = RODA_ASPFU). Rule 3.2.3 handles
  duplicates. The secondary test should name which LP proteins remain after it.

Fix. Fix the reserve now by a rule that does not use v1 output (for example, Jensen proteins whose cluster has no T1/T2 member,
listed in L1). Do not score the reserve in v1 (move it out of 6.7). Or drop the LP reserve and use only T1 proteins added after the
freeze date (section 8).

### N8. MINOR. The cost limit does not bind

Evidence. Unlabelled relaxed-only calls with the conditions, maximum per proteome, per 10,000 proteins:
default filters 1.61 at 8.5 bits, 2.68 at 5 bits (*F. graminearum*); `--nobias` 4.47 at 7 bits (equal GA, *F. graminearum*).
All are below 5. Without the R0 and 8-Cys conditions there are 46 relaxed-only rows at 8.5 bits in the 12 proteomes. With them, 10. So
the cysteine condition does the filtering, and the cost limit does not choose between cutoffs.

Fix. State this. Either keep the limit as a sanity bound, or set it from these numbers before the freeze.

### N9. MINOR. The fold proteome sample is not defined

6.1 names "a fixed, cluster-split proteome sample" for scoring folds. Source proteomes, size, the clustering it is split by, and
any criterion that uses its fold results are not given. If it comes from the 12 cost proteomes, say so.

### N10. MINOR. Integration details for L6

- A run without `hydrophobin_relaxed` is tolerated (`cli.py` lists it in `absent_modules`). Then `hydrophobin_extended` is
  `not_assessable` for every protein without a Pfam hit, and `_eval_other` lists it as left out in `other_basis` for those proteins.
  L6 must add the relaxed `hmmsearch` job to the standard run.
- L6 should name the test and output files that change, as H1 did: `tests/cellsurface_sorting_hat/test_data_and_scripts.py`
  (module-set assertion), golden `calls.long.expected.tsv` and `report.expected.md`, `outputs.py`.
- `params` should hold the GA pair, the score rule (N3), the search options (N2), the minimum cysteine count (8) and the R0 condition.
- 4.2 writes the three-branch call without `{flag: ...}` nodes. Write it in the config syntax.

### N11. MINOR. Provenance of the derived file

`/bigdata/operations/pkgadmin/srv/projects/db/pfam/current` is a symlink to `2026-01-27-Pfam38.2`. Record the resolved path.
`Pfam-A.hmm` has GA lines in two formats (`23.00 23.00;` and `27 27;`). The build script should check that exactly seven GA lines
were rewritten. The `# Query file:` line holds the path as given on the command line, so the wrapper should compare the base name
and also the model names and lengths (`qlen`).

### N12. MINOR. Labels for relaxed Pfam

6.8 says "the shipped model is trained on all clusters ... `tuned_on_truth`". The relaxed cutoff is also chosen on all 28 clusters,
which include the 6 proteome misses. Say that 6.8 applies to relaxed Pfam too. Section 8 names the catalog `predicted_by_hmm_v1`.
If v1 is relaxed Pfam only, use a neutral label.

### N13. MINOR. Small gaps

- 6.2 refers to review 1 for the cluster list. Put the list in the spec (section 2 above has it).
- L1b flags T2 entries. Say whether a flagged entry leaves T2 before L2 (that can change the 28 clusters and the 6 Pfam-missed clusters).
- 6.3 cites 40 *P. ostreatus* hydrophobin genes. The strain of Xu 2021 was not checked against PC9. At 8.5 bits relaxed Pfam makes 1
  unlabelled extra call in PC9, so the effect is small at present.

## 4. Contradictions with the earlier spec and owner decisions

- Earlier spec 3.4 (`hydrophobin_protein` = Pfam OR `cys8_pattern`): the report applied rule 6.3 and did not add it. Revision 2 keeps
  `cys8_pattern` as a side column. No contradiction.
- Earlier spec risk 7 ("the cutoff is not changed here") and section 11 of revision 2: `hydrophobin_domain` stays at GA. Relaxed Pfam is a
  second level. No contradiction.
- E7 versus 6.5: see N6.
- E4 ("a v2 model is tested on proteins that v1 never saw") versus 6.7 and 3.2.5: see N7.
- E10 ("the custom HMM ... is added only if it recovers more"): 6.2 makes "more" mean "at least 2 more clusters". With the numbers in N1
  that is at most the full remainder. The owner should confirm the threshold with those numbers in view.

## 5. Verdict

**Not ready for a plan.** The module, the GA-rewrite mechanism and the call fit the existing code (checks a and e). The facts in
sections 1, 3.2 and 6.2 are correct as stated for default HMMER filters and a full-sequence score. The blocker is N1: the HMM decision
rule is not defined well enough to compute, and in the computable readings it cannot be met or needs every remaining cluster. N2 and N3
change the numbers that N1 depends on and must be fixed first. N4 to N7 need a sentence or a rule each. After a revision 3 that fixes
N1 to N7, a plan can follow.
