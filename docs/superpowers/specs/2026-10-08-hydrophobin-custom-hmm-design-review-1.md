# Review 1 of the custom hydrophobin HMM design

*2026-10-08. Independent review of `2026-10-08-hydrophobin-custom-hmm-design.md` (rev 1, commit 092d26c), branch
`hydrophobin-validation`. Reviewer: Claude Code (claude-opus-5-5). The spec was not edited. Every number below was
re-derived by the reviewer with the commands named in section 0. Nothing here is committed.*

## 0. Method

- Truth and clusters: `analysis/hydrophobin_truth/truth_all.tsv`, `sp_hydrophobin_query.tsv`, `clusters_positives.tsv`,
  `nopfam9.domtbl`, `h6/summary.tsv`, `proteome_map.tsv`, `run_list.tsv`. Read with `/usr/bin/python3.12`.
- T2-only clustering (the spec's L2 rule): `mmseqs easy-cluster --min-seq-id 0.3 -c 0.5` (mmseqs2/17-b804f) on the 131 T2 sequences.
- Pfam checks: the seven hydrophobin-class models taken with `hmmfetch` from
  `_workdir/sorting_hat/calibration/hydrophobin_discovery/models.hmm` (Pfam 38.2), HMMER 3.4.
  (a) `hmmsearch --cut_ga --domtblout` on the 131 T2 sequences. (b) `hmmsearch -T 8 --domT 0 --tblout` on the 12
  proteomes of `run_list.tsv`, full-sequence bit score, joined to the existing `pfam_hydrophobin` and `step1_rule@R0`
  module tables in each work directory, and to a cysteine count of the full sequence.
- Code read: `categories.yaml`, `modules/cys8.py`, `modules/pfam.py`, `modules/cli.py`, `modules/base.py`,
  `engine.py` (`_validate`, `_eval`, `call_eligible`, `call_hash`), `calibration/cli.py` (`truth`, `_check_call_run`).
- HPC modules: `module avail` for hmmer, mafft, muscle, famsa, clustalo, clustal-omega, clustalw, mmseqs2, and a grep of
  the full list for clustal, hhsuite and t-coffee.
- Scratch outputs are in the session scratchpad, not in the repository.

## 1. Facts in sections 1 and 9

| Spec claim | Re-derived | Verdict |
|---|---|---|
| 9 of 174 named entries (5.2%) have no hydrophobin-class Pfam cross-reference; all 9 are T2 | 9 of 174 = 0.0517; all 9 T2 | correct |
| Pfam finds 36 of 42 labelled hydrophobins in 8 proteomes | `h6/summary.tsv` pfam_only TP sum 36, positives 42 | correct, see M11 |
| 5 of the 9 have a weaker Eas or Hydrophobin hit, E 2e-8 to 1e-3 | `nopfam9.domtbl`: HCF2 2.1e-8, PSH 2.4e-8, HYD1_GIBZE 7.6e-7, HFBE 7.4e-6, HYD2_GIBZE 9.9e-4 | correct, see m3 |
| PSH_FLAVE is 0.9 bits below the domain GA | domain score 22.1, PF01185 GA 23.00 23.00 | correct. Full-sequence score is 23.9, above the sequence GA. See M6 |
| T2 131, T3 43 | 131 and 43 | correct |
| 122 of 131 T2 (93%) have a Pfam cross-reference | 122, 0.931 | correct |
| Spacing rule recovers 1 of 6 Pfam misses; matches 60% of T2 | report: HYD2_GIBZE only; 78 of 131 = 59.5% | correct |
| HYD2B_BEAB2 has 6 Cys, HFBA_PENEN 7 | 6 and 7 | correct, but incomplete, see m4 |
| "About 130 T2 proteins in 30 clusters" (2, 9.1) | 30 clusters is the count over all 174 entries. T2 alone: 28 clusters (both in `clusters_positives.tsv` and in a T2-only re-clustering) | minor error, m2 |

An extra check supports the spec's Pfam-missed definition. On the 131 T2 sequences, `hmmsearch --cut_ga --domtblout` hits
exactly the 122 entries with a cross-reference (symmetric difference 0). So "no cross-reference" and "no Pfam 38.2 GA domain
hit" pick the same 9 proteins.

## 2. Findings

### B1. BLOCKER. The success criterion can be met by lowering the Pfam cutoff, and no such baseline is in section 6.8

Evidence.
- Full-sequence bit scores of the 9 against the Pfam models: PSH_FLAVE 23.9, HCF2_FULFL 23.6, HYD1_GIBZE 18.6, HFBE_PENEN 15.9,
  HYD2_GIBZE 8.6. The other 4 have no hit at E <= 1e-3 (`nopfam9.domtbl`).
- A full-sequence cutoff of 8.5 bits on the seven Pfam models, with the same conditions as the HMM branch (R0 `called`, at
  least 8 Cys), reaches 5 of 9 by score. (R0 was not run on the Swiss-Prot sequences; see M8.)
- The same rule in the 12 proteomes gives these calls that the current `pfam_hydrophobin` does not make: Af293 0, A1163 0,
  W72310 0, S288C 0, *C. albicans* 1, *C. immitis* RS 0, *B. dermatitidis* ER3 1, *B. bassiana* 0, *F. fulva* 3,
  *F. graminearum* 3, *P. expansum* 1, *P. ostreatus* PC9 1. The maximum is 3, below the limit of 5.
- At 15 bits the extra calls are 0 or 1 per proteome. They include 3 of the 6 labelled proteome misses
  (HCF2_FULFL 23.6, HYD1_GIBZE 18.6, HFBE_PENEN 15.9).

So the criterion of section 6.2 and 6.3 is met by Pfam with a lower cutoff, without a new HMM. A pass would not show that
the custom HMM adds anything. Section 11 excludes "changing the Pfam cutoffs". That is about the tool. A baseline is a
measurement, not a change.

Fix. Add a baseline "Pfam hydrophobin models, full-sequence bit cutoff chosen by the same procedure on the same development
data, same R0 and cysteine conditions" to 6.8. Make the criterion relative. For example: the HMM recovers at least k more
Pfam-missed clusters than that baseline, at an equal or lower count of extra calls. Fix k before L4. Say which score
(sequence or domain) the baseline uses (M6).

### M1. MAJOR. The recall target counts proteins, but the 9 Pfam-missed proteins are 6 clusters

Evidence (T2-only clusters, 30% identity, coverage 0.5):

| Cluster (representative) | Members | Pfam-missed |
|---|---|---|
| A0A6V8R0V1 | HYD1_TRIAP, HFB3_HYPVG, HYD2_GIBZE | 3 of 3 |
| Q7Z9L5 | HCF1_FULFL, HCF2_FULFL, HCF3_FULFL | 2 of 3 |
| Q4WEK0 | HFBE_PENEN, RODF_ASPFU | 1 of 2 |
| I1RDP9, Q0KKA0, Q7Z9L4 | HYD1_GIBZE, PSH_FLAVE, HCF4_FULFL | singletons |

"At least half" of 9 proteins is 5. Two clusters (3 + 2) reach it. Counted by cluster, half is 3 of 6. The Wilson 95%
interval for 3 of 6 is about [0.19, 0.81]. The target cannot tell a useful model from a weak one. Literature T1 proteins
may add clusters, but section 3.2.4 says the count is unknown.

Fix. State the target in clusters. Report the interval with it. Pre-register the number of Pfam-missed clusters that must
be reached. Say in the spec that with 6 clusters the test is a smoke test of recall, not an estimate.

### M2. MAJOR. The development split and the leave-cluster-out (LCO) test overlap, and the development part is too small

Evidence.
- 28 T2 clusters. A 30% development part is about 8 clusters.
- With 8 of 28 clusters in development, the chance (hypergeometric) that 0 of the 6 Pfam-missed clusters fall in
  development is 0.10, and that 3 or more fall in it is 0.21.
- The class II (PF06766 cross-reference) proteins are in 4 clusters of sizes 33, 3, 1 and 1. The chance that no class II
  cluster falls in development is about 0.24. Then architecture (b) cannot be built on development data.
- Two clusters hold 40 and 33 of 131 T2 proteins (56%).
- Section 6.1 holds out "every cluster once". That includes the development clusters, which chose the cutoff and the
  architecture. Their LCO results are in-sample.
- Section 4.1.1 checks the alignment "by eye". That cannot be repeated per fold. If one alignment of all proteins is cut per
  fold, the held-out proteins shaped the alignment.

Fix. Choose one of two designs and write it in the spec. (a) Nested: in each fold, choose the cutoff (and the architecture)
on the training clusters only. (b) Split: put all Pfam-missed clusters in the test part on purpose (a stratified split),
choose the cutoff on the Pfam-positive development clusters and development hard negatives, and report LCO only on test
clusters. In both, build each fold's alignment with a script (for example `mafft` or `famsa`, then `hmmbuild`), not by eye.
`hmmbuild` is fast, so plain leave-one-cluster-out (28 folds) is feasible. It removes the need for the grouping rule.

### M3. MAJOR. "Every group has at least one positive cluster" says nothing, and the negatives in LCO are not defined

Evidence. Only positives are clustered (3.3). Every LCO group therefore has a positive cluster. The LCO test needs
negatives to set and check the cutoff. Section 6.1 does not say which negatives are used in a fold, or whether the hard
negatives (L3) are split.

Fix. Say which negatives each fold uses (hard negatives split by cluster into development and test, or a fixed negative
set never used for tuning). Replace the grouping sentence with the fold rule from M2.

### M4. MAJOR. The cutoff is tuned on the evaluation data in two places

Evidence.
- Section 5: "A cutoff choice is judged against the criterion in section 6". The criterion includes the LCO recall and the
  proteome cost. That is choosing the cutoff on the test.
- Section 5: after the model is built, families from its false positives in the measured proteomes join the hard
  negatives, "and the cutoff is re-checked". The same proteomes then measure the cost limit (6.3). The cost limit becomes
  in-sample.
- Section 6.7 says the cutoff is fixed before any proteome result is seen. The two sentences in section 5 contradict it.

Fix. For v1, choose the cutoff once on development data (4.1.3) and freeze it (6.7). Measure the proteome cost once. Move the
false-positive loop to v2 (section 8), with new proteomes for its test.

### M5. MAJOR. The cost limit is not defined well enough to measure

Evidence and open questions.
- "Calls that Pfam does not make" may include labelled true positives. The 6 labelled proteome misses are T2 training
  proteins of the shipped model, so the HMM will call them. *F. fulva* has 3 of them, which uses 3 of the 5.
- It is not said whether `pfam_hsba` counts as "Pfam". HsbA had 9 of 16 hits in Af293 (rev 3, F6), and section 5 lists
  HsbA as a hard negative.
- Proteome sizes are 6,212 to 13,560 proteins (`grep -c '>'` on the 12 FASTA files). "About 10,000" does not say whether
  the limit scales.
- "Each measured proteome" is not listed. There are 8 with truth and 4 without.
- An unlabelled true hydrophobin counts as cost. Rev 3 F8 cites 40 *P. ostreatus* hydrophobin genes (PMID 33636611).
  `pfam_hydrophobin` calls 10 in PC9. The strain of that paper was not checked by this review, so the size of this effect is
  unknown.

Fix. Define the count as HMM-only calls with no T1/T2 label, after removing proteins called by `pfam_hydrophobin` or
`pfam_hsba` (or say otherwise). Give the limit per 10,000 proteins. List the proteomes. Report labelled HMM-only calls in a
separate column. Say in the risks that true unlabelled hydrophobins count as cost.

### M6. MAJOR. `hydrophobin_ext` cannot reuse the Pfam path as written, and the score is not defined

Evidence.
- Wrappers do not run tools. A job script runs `hmmsearch` (`scripts/sorting_hat/pfam_hmmsearch.sbatch`). The `pfam`
  subcommand of `cellsurface_sorting_hat_module` reads its `--domtbl`. E5 says the module "runs `hmmsearch`".
- `pfam.parse_domtblout` refuses a table without `--cut_ga` in its header. A run with `-T <bits>` would be refused.
- `parse_domtblout` reads the domain score (field 14). "Best HMM score" in 4.2 does not say sequence or domain.
  PSH_FLAVE shows the two differ at the cutoff: sequence 23.9, domain 22.1.
- The current pattern for module inputs works. `cys8` reads R0 with `_condition_table` and puts R0's
  `params_hash` and `artefact_hash` in `params["conditions"]`. `hydrophobin_ext` can do the same.

Fix.
- Add a job script that runs `hmmsearch` with `data/sorting_hat/hydrophobin_ext.hmm`. Add a wrapper subcommand that reads the
  table.
- Write the frozen cutoff into the HMM file as `GA` lines. `--cut_ga` then works, the existing parser check holds, and the
  cutoff travels with the file. Otherwise add a parser that accepts `-T`/`--domT` and checks the header value against
  `params`.
- Name the score (sequence or domain) in the spec and in `params`.
- Pass the `.hmm` file in `ModuleSpec.artefacts`. `artefact_hash` is a hash of `name:sha256`, so a changed HMM changes the
  module identity. Then the call file goes stale (03 section 6a). Put in `params`: cutoff, score type, minimum cysteine count,
  model version, and the R0 condition. The provenance JSON does not need to be in the identity if the HMM hash is.
- Check in the wrapper that the domain table came from that HMM (the `# Query file:` header line, or the model names and
  lengths).

### M7. MAJOR. Hard negatives have no pass rule and conflict with the HsbA decision

Evidence.
- Section 5 refers to section 6 for the judgement. Section 6 has no hard-negative criterion. Only 6.5 asks to report the
  rate.
- HsbA is a hard negative here (E6). Rev 3 D6 says the owner counts HsbA as a hydrophobin for cataloging, and category 2c
  displays `hydrophobin_domain` with `hsba_domain`. `hsba_domain` is a mechanism of the `other_*` calls in `categories.yaml`.
- The negative lists of CFEM, cerato-platanin, HsbA and PIR come from Pfam membership. The Pfam hydrophobin models are not
  scored on them in the plan.
- No count of these sets exists yet. Task L3 makes them.

Fix. Add a per-group limit to section 6 (for example, a call rate at or below the rate of the relaxed-Pfam baseline of B1).
Split hard negatives into development and test if the cutoff uses them. Report Pfam GA, relaxed Pfam and the HMM on the same
hard negatives. Ask the owner whether an HMM hit on an HsbA protein is a false positive or an expected overlap.

### M8. MAJOR. No task runs SignalP on the truth and hard-negative sequences

Evidence. Section 6.2 measures recall "with the signal-peptide and cysteine conditions". R0 exists only for the 12
proteomes (`analysis/hydrophobin_truth/h5_signalp.*.log` lists proteomes). The LCO test runs on Swiss-Prot and literature
sequences. UniProt has no signal-peptide annotation for 2 T2 entries (RODD_ASPFU, RODE_ASPFU). That is not R0.

Fix. Add a task before L5: run the R0 SignalP job (same build and mode) on the T1, T2, T3 and hard-negative sequences, and
record the identity.

### M9. MAJOR. The v2 test set of E4 is not reserved

Evidence. E4 tests v2 on "held-out literature proteins that v1 never saw". E2 and L1 put every verified literature T1 protein
into v1 training. Section 8 builds v2 from promoted catalog proteins. No literature protein is held back.

Fix. In L1 or L2, reserve a literature set by cluster, before v1 training. Or define the v2 test as T1 proteins added after
the v1 freeze date, and record that date.

### M10. MAJOR. The new criterion is weaker than rule 6.3 of rev 3, and no reason is given

Evidence. Rev 3 rule 6.3 kept a rescue only with at least 5 resolved rescue-only clusters and a Wilson lower bound of
precision of at least 0.5. The new criterion has no precision term. It caps the count of extra calls. Owner review (section
7) is collected but does not enter the decision.

Fix. Either keep a precision term on resolved HMM-only calls by cluster, as in rule 6.3, or state why the count cap replaces
it.

### M11. MAJOR. The proteome numbers and the Swiss-Prot numbers are the same proteins

Evidence. The 42 positives are 42 proteome-protein pairs of 31 unique Swiss-Prot accessions (the 3 *A. fumigatus* strains
repeat RodA to RodG). The 6 proteome misses are HCF1, HCF2, HCF4 (*F. fulva*), HYD1_GIBZE, HYD2_GIBZE (*F. graminearum*)
and HFBE_PENEN (*P. expansum*). They are 6 of the 9 Swiss-Prot misses (report 2026-10-08 section 4; README H2b). The shipped
model is trained on them.

Fix. Say in sections 1 and 6 that the proteome recall of the shipped model on these 6 is a training-set result. Report
Pfam recall as 25 of 31 unique accessions beside 36 of 42 pairs.

### M12. MAJOR. Architecture (b) depends on class labels that come from the rules under test

Evidence. The only class labels at hand are the Pfam family (PF01185 class I, PF06766 class II) and the `cys8` spacing class.
The "divergent/unclassified" group would be mostly the Pfam-missed clusters (4 of the 6 clusters have only Pfam-missed
members). A group model of those proteins is then built from the recall test set. In LCO it has 5 clusters left. Each extra
model also adds a cutoff to tune on about 8 development clusters.

Fix. Say how a training protein gets its group, from which source. Use one shared cutoff or fix the per-model cutoffs by one
rule. Report the recall of the "divergent" model on Pfam-missed clusters only from folds that hold them out.

### m1. MINOR. The leakage label of the LCO numbers

Evidence. The author has seen the 9 (`nopfam9.domtbl`, `nopfam9_regex.tsv`, the review lists). Section 6.6 calls the LCO
numbers "the honest estimate". Rev 3 5.3 labelled the same knowledge `partial`.

Fix. Label the LCO numbers `partial` too. Keep the 9 out of any manual alignment edit.

### m2. MINOR. Cluster count

"About 130 T2 proteins in 30 clusters" (2 and 9.1). T2 proteins are in 28 clusters. 30 is the count over all 174 entries.

### m3. MINOR. E-values in section 1 come from a 9-sequence search

E-values scale with database size. In a proteome search at E <= 1e-3, only HCF2_FULFL and HYD1_GIBZE are found among the
labelled misses (this review, run (b)). Give bit scores (23.9, 23.6, 18.6, 15.9, 8.6) and say the database size.

### m4. MINOR. The cysteine condition loses 6 T2 proteins, not 2

T2 with fewer than 8 Cys in the full sequence: HYD2B_BEAB2 6, HYD2_CORMI 6, HFBA_PENEN 7, QID3_TRIHA 7, HYD3_BIOOC 7,
HFBD_PENEN 7. T3: HYD3_COPC7 6. All 6 T2 have a Pfam GA domain hit, so `hydrophobin_extended` keeps them through the Pfam
branch. The HMM branch alone misses them. The LCO table of the HMM alone should show this. Also say whether "8 cysteines"
is the full sequence or the aligned region.

### m5. MINOR. The two branches of `hydrophobin_extended` use different conditions

The seven Pfam rows in `data/sorting_hat/family_table.tsv` have no `second_condition`. The HMM branch needs R0 `called` and at
least 8 Cys. The call is valid, but the baselines in 6.8 compare unlike rules. Report Pfam with and without the same
conditions.

### m6. MINOR. The config change makes existing call files stale

Adding `hydrophobin_extended` to `categories.yaml` (and to the mechanism lists) changes `config_sha256`. The
`tandem_repeat_protein` call files go stale again (ledger row H4 re-measured them after the last change). Add a task like
H1b after the last edit. Say whether `hydrophobin_extended` replaces `hydrophobin_domain` in the mechanism lists or is added.

### m7. MINOR. T2 labels are not reviewed

`analysis/hydrophobin_truth/README.md` says the T2 rule does not check that the experimental statement is
hydrophobin-specific, and asks for a review before a label is final. The spec trains on T2 with no review task.

### m8. MINOR. "Independent evidence" in section 8 includes orthology

Orthology to a T1/T2 protein "by a method that is not the HMM" is still sequence similarity to the training set. Say which
methods count, and that a BLAST or phmmer hit alone does not.

### m9. MINOR. Literature tiers

Rev 3 F8 says Xu 2021 found the 40 *P. ostreatus* genes by the cysteine pattern. Most literature proteins from genome surveys
will be "literature, predicted". They must never enter the Pfam-missed test set, because their selection used the pattern or
a domain search. State this in 3.2.

## 3. Implementability summary (item 2 of the request)

- `hydrophobin_extended` as written (`or` of two `flag` nodes on two modules, no `ref`, no literal step 1 module, not
  `kind: other`) passes `call_eligible`. `call_hash` covers the expression and `ENGINE_SEMANTICS`. Module identities are
  checked separately by `_check_call_run` through `run.json`.
- A mechanism must be an earlier ungated call (`engine._validate`). The call is ungated, so it is valid if it is placed
  before `other_not_surface`.
- With Pfam `hit` = 0 and HMM `hit` empty (missing R0), the `or` gives `not_assessable`. That is the `cys8` behaviour and is
  correct.
- `calibrate truth --call-status` can measure the call on a proteome run if both modules are `ok` in `run.json` and the
  config hash matches. It cannot produce LCO numbers. LCO needs an analysis script, and gives no status file. No call reads
  `hydrophobin_ext` alone, so that module gets no module status. This is consistent with decision D3.
- Tools on HPCC (`module avail`): hmmer 1.8.5, 2.3.2, 3.3.2, 3.4 (default; `hmmbuild`, `hmmalign`, `hmmsearch` present),
  mafft 7.490 and 7.505, muscle 3.8.31, 5.1 and 5.3, famsa 2.4.1, tcoffee 13.45.65, mmseqs2 13, 15, 17-b804f. No clustalo,
  clustal-omega or clustalw module was found.

## 4. Task order

Missing or misplaced tasks:
1. R0 SignalP on truth and hard-negative sequences (M8). Before L4.
2. Reserve a v2 literature test set (M9). In L1 or L2.
3. The relaxed-Pfam baseline (B1). In L4 (same tuning procedure) and L5.
4. Fold rule and negatives per fold written before L4 (M2, M3).
5. Hard-negative development/test split (M7). In L3.
6. Re-measure stale call files after the `categories.yaml` edit (m6). After L6.
7. The false-positive loop of section 5 moves to v2 (M4).

## 5. Verdict

**Not ready for a plan.** The module and call can be built with the existing code patterns (section 3). The blocker is
the success criterion: a lower Pfam cutoff already meets it in the data at hand (B1). The evaluation also needs a fold rule
that keeps tuning data out of the test (M2 to M4), a measurable cost limit (M5), and a hard-negative rule (M7). After
these are fixed in a revision 2, a plan can follow.
