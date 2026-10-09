# Extended hydrophobin evidence level: relaxed Pfam first, custom HMM tested against it

*2026-10-08. Revision 4 (after independent reviews 1 and 2 of this spec and review 1 of the plan, which found that the proteome sample never limits the cutoff; the cutoff rule and the proteome split are new in this revision; review 1: 1 blocker, 12 major, 9 minor; review 2: 1 blocker, 6 major, 6 minor; all addressed below; reviews are
`2026-10-08-hydrophobin-custom-hmm-design-review-1.md` and `-review-2.md`). Author: Claude Code (claude-sonnet-5-5). Builds on
`2026-10-08-hydrophobin-validation-design.md` (rev 3) and `docs/reports/2026-10-08-hydrophobin-validation.md`. The owner agreed
the decisions in section 2 in a question-by-question interview on 2026-10-08. Status: draft. No code is written.*

## 1. Why

1. The Pfam hydrophobin models miss known hydrophobins. In Swiss-Prot, 9 of 174 named entries (5.2%) have no hydrophobin-class Pfam
   cross-reference. All 9 have experimental (T2) evidence. On the 131 T2 sequences, `hmmsearch --cut_ga` hits exactly the 122
   entries with a cross-reference, so "no cross-reference" and "Pfam miss" are the same 9 proteins (reviewer-checked).
2. In the 8 measured proteomes Pfam finds 36 of 42 labelled hydrophobin-proteome pairs. By unique Swiss-Prot accession it is 25 of 31
   (the three *A. fumigatus* strains repeat RodA to RodG). The 6 proteome misses are 6 of the 9 Swiss-Prot misses, so they are not an
   independent result (HCF1, HCF2, HCF4 *F. fulva*; HYD1_GIBZE, HYD2_GIBZE *F. graminearum*; HFBE_PENEN *P. expansum*).
3. Bit scores of the 9 against the seven Pfam models (full sequence, searching 9 sequences): PSH_FLAVE 23.9 (its domain score is 22.1,
   0.9 bits below the domain gathering cutoff of 23.0), HCF2_FULFL 23.6, HYD1_GIBZE 18.6, HFBE_PENEN 15.9, HYD2_GIBZE 8.6. Four
   (HYD1_TRIAP, HFB3_HYPVG, HCF1_FULFL, HCF4_FULFL) have no hit at E <= 1e-3. E-values scale with database size, so bit scores are used from here on.
4. The 8-cysteine spacing rule built from published spacings (`cys8_pattern`) recovers 1 of the 6 proteome misses and matches 60% of
   T2 hydrophobins. It is not a usable rescue.
5. The reviewer ran the seven Pfam models on the 12 proteomes at a full-sequence cutoff of 8.5 bits with the conditions R0 `called` and at least 8
   cysteines. Extra calls per proteome were 0 to 3 (Af293, A1163, W72310, S288C, *C. immitis*, *B. bassiana* 0; *C. albicans*, *B. dermatitidis*,
   *P. expansum*, *P. ostreatus* PC9 1; *F. fulva* and *F. graminearum* 3). That cutoff reaches 5 of the 9 by score (R0 was not run on the Swiss-Prot
   sequences, so the signal-peptide condition is unchecked there). At 15 bits there are 0 or 1 extra calls and 3 of the 6 proteome misses come back.
   These numbers are the reviewer's. They are re-derived in task L4.

6. **The HMMER bias composition filter hides at least one of the 9 (review 2, to be re-derived in L4).** With `--nobias`, HCF1_FULFL scores 30.4 bits against Eas, above the Eas gathering cutoff
   of 27, so Pfam GA would recover it (123 T2 hits instead of 122, and the *F. fulva* proteome protein FEBB419A_005025). HCF4_FULFL scores 16.6 bits with `--nobias` (domain 14.7). The strict level
   keeps the standard options. The search options are therefore part of the relaxed level and are frozen with the cutoff (section 4.1).

Consequence: a relaxed Pfam cutoff is the cheapest extended level and must be the baseline. A custom HMM is justified only if it adds
recall beyond it.

Prior art exists for pattern-plus-domain searches (Yang 2006, PMID 17217508; Jensen 2010, PMID 21182770). This work does not claim novelty
of the approach. It adds a measured call with a status inside the tool.

## 2. Decisions agreed with the owner (2026-10-08)

| # | Decision |
|---|---|
| E1 | Optimise for recall with **two evidence levels**. `hydrophobin_domain` (Pfam, gathering cutoff) stays unchanged. `hydrophobin_extended` is the second level. |
| E2 | Train any custom model on T1 and T2 proteins only. T3 is held out as a separate check. |
| E3 | Small numbers are accepted. All results are `smoke`. Leave-cluster-out testing shows whether a model is consistent without memorising. |
| E4 | Order: freeze v1 and measure it; scan Fungi_5k with v1 as a **separate project**; predicted members are never truth without independent evidence; a v2 model is tested on proteins that v1 never saw. |
| E5 | Shipped models are versioned data files in `data/sorting_hat/` with a provenance file. |
| E6 | Hard negatives: CFEM, cerato-platanin, HsbA, small secreted cysteine-rich proteins (killer toxins, effector candidates), cell wall proteins with cysteines (PIR, Ccw12-like). More are added from the model's false positives, in v2 only. |
| E7 | The owner's adjudication is a separate tier (T5, "owner-reviewed"), recorded with a reason per protein, reported beside the experimental tiers and never merged. |
| E8 | The owner supplied Linder 2005, Sunde 2008, Wessels 1994 and the supplements of Jensen 2010 and Xu 2021 (`to_import/hydrophobins/`). Merge order of PRs #76, #75 and this branch stays with the owner. |
| E9 | Literature proteins without protein-level evidence form a "literature, predicted" tier: held out from training, used as an extra test set, reported separately. |
| E10 | **The first extended level is a relaxed Pfam cutoff.** The custom HMM is built and tested against that baseline, and is added only if it recovers more (section 6). If it adds nothing, that is reported and the relaxed Pfam level stays. |
| E11 | An extended-level or HMM hit on an HsbA protein is an **expected overlap**, not a false positive (owner, 2026-10-08). It is reported in its own column and excluded from the cost count and from the HsbA hard-negative rate. |

A machine-learning classifier was considered and is not part of this design. With about 130 experimentally supported proteins in 28 clusters
it is likely to overfit or to relearn the cysteine pattern. This is an expectation, not a result.

## 3. Data

### 3.1 Tiers

| Tier | Meaning | Training | Measurement |
|---|---|---|---|
| T1 | Hydrophobin identity or function shown in a paper at protein level (purification, mass spectrometry, rodlet or film formation, contact angle, or a hydrophobin-specific deletion phenotype). PMID recorded. | yes | primary |
| T2 | Swiss-Prot reviewed, experimental evidence (ECO:0000269) on a FUNCTION, SUBCELLULAR LOCATION or SUBUNIT comment | yes | primary |
| T3 | Named hydrophobin by rule or similarity only | no | separate check |
| LP | Literature, predicted: listed in a paper as a hydrophobin on the strength of a genome screen, EST or expression only | no | separate test set, secondary |
| T4 | Review judgement from sequence and BLAST, no paper | no | reported only |
| T5 | Owner-reviewed (E7) | no (v1) | separate column |

Counts now: T2 131 (in 28 clusters; 30 clusters is the count over all 174 named entries), T3 43, no T1. T2 labels have not been reviewed for hydrophobin
specificity (`analysis/hydrophobin_truth/README.md`). Task L1b spot-checks every T2 FUNCTION or location text for hydrophobin-related content
and lists entries whose experimental statement is about something else.

### 3.2 Literature (task L0 done)

The supplied files give **no protein-level evidence** and no accessions for 90 listed proteins (`analysis/hydrophobin_truth/literature/NOTES.md`):
- Jensen 2010 (PMID 21182770): 50 proteins from nine *Aspergillus* genomes. They are genome predictions made by a pattern, size and signal-sequence screen.
  All 50 gene IDs were resolved to sequences (FungiDB-31 local files, the local A1163 proteome, UniProt for Af293). The stated cysteine pattern is found
  exactly in 41 of 49 resolved sequences. For 8 it is not (current annotation or a transcription error). These 8 are not used until checked.
  With the seven hydrophobin-class Pfam models at the gathering cutoff, 43 of the 50 have a hit and 7 do not
  (AO090012000143, ATEG_10285, ATEG_08089, AFLA_060780, AFLA_014260, AFLA_063080, ACLA_001890).
- Xu 2021 (PMID 33636611): 40 gene names and primers, no class, accession or sequence. Not usable without the main text.

Rules:
1. A literature protein is T1 only with protein-level evidence in the paper. None of these 90 qualifies. They are LP.
2. LP proteins never enter training, and **never enter the primary recall criterion** (section 6). Their selection used a pattern or a domain search.
   They form a secondary test set: recall of each method on the Jensen proteins that Pfam misses, reported separately.
3. A protein already in Swiss-Prot is not added twice. Strain orthologs (the same protein in several genomes) are one cluster.
4. T1 entries found later (papers with protein-level evidence) carry the PMID and the table or figure. No accession is written from memory.
5. **v2 reserve.** Before any v1 result is seen, the LP proteins are clustered together with the T1/T2 proteins (30% identity, 0.5 coverage). LP clusters that contain no T1/T2 protein and none of the 8 unverified proteins are
   reserved as the v2 test set, and are excluded from the v1 secondary test (6.7). The reservation does not depend on any method's result. LP proteins that cluster with a T2 protein (review 2: ATEG_08089 and ACLA_001890 cluster with RODF_ASPFU)
   are not eligible. If fewer than 3 clusters are reserved, the v2 test uses only T1 proteins added after the v1 freeze date.

### 3.3 Clusters and folds

MMseqs2 at 30% identity and 0.5 coverage. The folds and the Pfam-missed clusters come from `analysis/hydrophobin_truth/clusters_positives.tsv` (28 T2 clusters, the 9 Pfam-missed proteins in 6; sha256 recorded in the freeze). The joint clustering of LP and T2 proteins is used only to choose the v2 reserve (a different clustering gives 5 Pfam-missed clusters, not 6, because HYD1_GIBZE joins another cluster).

## 4. Models

### 4.1 Level 1 (strict) and level 2 (extended)

- `hydrophobin_domain`: unchanged (Pfam 38.2, `--cut_ga`, module `pfam_hydrophobin`).
- **Relaxed Pfam (E10, v1 extended level).** A derived file `data/sorting_hat/hydrophobin_relaxed.hmm`: the same seven hydrophobin-class Pfam models,
  fetched with `hmmfetch`, with the `GA` lines rewritten to one frozen full-sequence bit score. `hmmsearch --cut_ga` then applies it, the existing
  parser check (`# --cut_ga` in the header) holds, and the cutoff travels with the file. The sequence GA is the frozen full-sequence cutoff. The domain GA is set to **-1000.00**, so the domain test never binds and a hit is decided by the full-sequence score (the plan reviewer found that a domain GA of 0.0 drops proteins whose domains all score below 0, for example
  with `--nobias`; `GA 6.00 -100.00;` is accepted by HMMER). Task L4 tests this with `hmmsearch --cut_ga`, the parser, and a fixture protein whose domains all score below 0. The search options (default filters, or `--nobias`) are chosen with the cutoff by the nested procedure (section 6.1), frozen in the job script, and named in `params`.
  A provenance file records the Pfam release, the source model hashes, the cutoff, how it was chosen, and the file hash.
- **Module `hydrophobin_relaxed`** (one row per protein): `hit` = the protein has a hit to any of the seven models at the frozen cutoff **and** at least 8 cysteines
  in the full sequence **and** the R0 call is `called`. The R0 module identity goes into `params["conditions"]`, as in `cys8_pattern`. The `.hmm`
  file goes into `ModuleSpec.artefacts` so that its hash is in the module identity. A row has an empty `hit` when it would be a hit except that the R0 row is
  missing. The wrapper checks that the domain table came from this file (query file name and model names). A job script runs `hmmsearch`.
- **Call:**

```yaml
- name: hydrophobin_extended
  expr:
    or: [{flag: pfam_hydrophobin.hit}, {flag: hydrophobin_relaxed.hit}]
```

  The call is ungated, reads two modules, has no `ref`, and passes `call_eligible`. It is placed before the `other_*` calls and **added** to their mechanism
  lists (it does not replace `hydrophobin_domain`). Because it is a superset, every Pfam call is also an extended call. The report shows whether a protein is
  called by Pfam, by the relaxed level only, or by both.

### 4.2 Custom HMM (candidate, tested against the relaxed baseline)

One family-wide HMM. Grouped models (class I, class II, divergent) are dropped for v1: the class labels come from Pfam or from the spacing rule under test, and
the divergent group would be built from the Pfam-missed proteins that make up the recall test.

Build per fold, by script, never by eye: align the training proteins with `mafft` (7.505) or `famsa` (2.4.1) (the tool is chosen in the plan; `clustalo` is
not installed), run `hmmbuild`, then `hmmsearch` of the held-out proteins. The alignment is checked on the cysteine columns by a script (the eight conserved
columns must be present in at least a stated fraction of sequences).

If the HMM passes section 6, it ships as `data/sorting_hat/hydrophobin_ext.hmm` (GA lines carry the frozen cutoff), with a provenance file (training list with
cluster IDs and tiers, aligner and HMMER versions, commands, file hash), and a module `hydrophobin_ext` built like `hydrophobin_relaxed`. The call then becomes
`or: [pfam_hydrophobin.hit, hydrophobin_relaxed.hit, hydrophobin_ext.hit]`. If it fails, the model and its numbers are reported and nothing else changes.

`cys8_pattern` and its spacing file stay as a reported side column. They are not part of any default call.

## 5. Hard negatives

Groups (E6): CFEM; cerato-platanin; HsbA; small secreted cysteine-rich proteins (killer toxins, effector candidates); cell wall proteins with cysteines
(PIR, Ccw12-like). Sources: Swiss-Prot reviewed fungal entries per group, plus members of the corresponding Pfam families. Each group is clustered. Groups are split by
cluster into a tuning part and a test part before any model is scored. The tuning part and a fixed sample of proteome negatives are the only negatives
that may be used to choose a cutoff. R0 SignalP is run on all of them (task L4a). Counts per group are reported before scoring.

**HsbA** (E11, owner decision): a hit on an HsbA protein is an expected overlap, reported in its own column, and excluded from the cost count (6.3) and from the HsbA hard-negative rate. HsbA proteins stay in the hard-negative table so the overlap is visible.

Other look-alikes (LysM effectors, fungal defensins and antifungal proteins, expansin-like proteins) are not verified. They are added in v2 from the model's
false positives in new proteomes, not in v1 (the same proteomes cannot both add a negative group and then test the cost).

## 6. Evaluation and success criteria

### 6.1 Folds

Leave-one-cluster-out over the 28 T1/T2 clusters (`hmmbuild` is fast). In each fold, all choices (alignment, cutoff for the HMM and for the relaxed Pfam) are made on the
training clusters only (nested): the cutoff is the **lowest bit score at or above a floor of 0 bits** at which both conditions hold:
(i) on the **tuning proteomes** the unlabelled extra calls (conditioned on R0 and at least 8 cysteines, excluding proteins that `pfam_hydrophobin` or `pfam_hsba` call and excluding T1/T2/LP proteins) are at most
5 per 10,000 proteins in each tuning proteome (the same limit as 6.3), and (ii) the call rate on the tuning part of each hard-negative group (HsbA excluded, E11) is at most the per-group rate fixed in L3.
No positive protein's score sets a cutoff, so a model is never tuned on the scores of its own training proteins. Ties go to the higher cutoff. Review of the plan showed that a 1,000-protein sample per proteome never limits the
cutoff (0 or 7 conditioned calls in 12,000 proteins at 0 bits), so whole proteomes are used for tuning, split from the proteomes used for the cost test.

**Proteome split (seeded, once, before any scoring).** The 12 proteomes form 10 species groups (the three *A. fumigatus* strains are one group). The groups are split 5 and 5, with the constraint that each part holds at least two groups with truth positives.
The tuning proteomes tune cutoffs. The test proteomes measure cost (6.3). The split is stored in `proteome_split.tsv`.

**Search option.** For the relaxed Pfam level and the HMM the default filters and `--nobias` each get a cutoff by the rule above. The option that recovers more Pfam-missed **training** clusters in the fold is chosen. Positive scores are used only
to choose between these two options, never to set a cutoff. The option is frozen with the cutoff.

The shipped cutoff is chosen by the same rule over all 28 clusters and frozen. The relaxed Pfam cutoff does not depend on the held-out cluster (no positive is an input), so every fold gives the same relaxed cutoff and the relaxed
leave-one-cluster-out numbers are the shipped-cutoff numbers on the training proteins. This is stated in the report. Only the HMM is rebuilt per fold.

The v1 freeze is a commit that fixes the aligner and its command line, the cutoff rule, the floor, the two search options, the thresholds of 6.2 to 6.5, the hard-negative rates, the proteome split, the folds and the cluster file hash, **before** any held-out, test-part,
LP or test-proteome result is looked at. Before the freeze only the tuning parts (tuning proteomes, tuning hard negatives) and the per-fold alignments are scored or built. The hashes of the pre-freeze score files go in the freeze file.

### 6.2 Recall (primary: T1 and T2 only)

Counted **in clusters**. The Pfam-missed set is the 9 Swiss-Prot entries, in 6 clusters (3, 2, 1 and three singletons; the clusters are listed in the review). Half of 6 is 3, and the Wilson 95% interval for 3 of
6 is about [0.19, 0.81]. This is a smoke test of recall, not an estimate. The same is reported for unique accessions.

- **Definition.** A Pfam-missed cluster is *recovered* by a method in a fold when at least one Pfam-missed member of the held-out cluster is called by the method (with its R0 and cysteine conditions) at the
  fold's cutoff. Recovery is counted over the 6 Pfam-missed clusters, once per cluster.
- **Relaxed Pfam passes** if it recovers at least 3 of the 6 clusters at the limits of 6.3 and 6.4. Its numbers are reported even if it fails.
- **Custom HMM.** Let u be the number of the 6 clusters that the relaxed baseline (with its frozen options) does not recover in the same runs.
  - If u is 0 or 1, the HMM cannot add recall that the primary test can detect. The result is "not testable on the primary set". The HMM is not shipped in v1. Its leave-one-cluster-out numbers on all 28 clusters and on the LP set are reported,
    and it stays a v2 candidate.
  - If u is 2 or more, the HMM is added only if it recovers at least ceil(u/2) of those u clusters, at a cost no higher than the relaxed baseline plus 2 calls per 10,000 proteins, and it passes 6.4 and 6.5.
  - These rules are fixed now. They can change only before the freeze.
- The recall of the HMM alone, and the count of T2 entries lost to the 8-cysteine condition, are reported. Six T2 entries have fewer than 8 cysteines in the full
  sequence (HYD2B_BEAB2 6, HYD2_CORMI 6, HFBA_PENEN 7, QID3_TRIHA 7, HYD3_BIOOC 7, HFBD_PENEN 7). All six have a Pfam hit, so `hydrophobin_extended` keeps them through the Pfam branch. The HMM branch alone does not.

### 6.3 Cost (proteomes)

Measured on the **test proteomes** (out of sample for the cutoff). The tuning proteomes are reported too and are in-sample. Calls made by the extended level that `pfam_hydrophobin` and `pfam_hsba` do not make, and that have no T1/T2 label, counted per proteome and per 10,000 proteins, in each of the 12
proteomes of `analysis/hydrophobin_truth/run_list.tsv` (8 with truth, 4 without; proteome sizes 6,212 to 13,560). The limit is **5 per 10,000 proteins** in each test proteome. `pfam_hsba` hits are the expected overlap (E11). Proteins that are near HsbA but below its gathering cutoff count as cost.
Labelled extra calls are reported in a separate column. Unlabelled true hydrophobins count as cost. Probably hydrophobins (for example in *P. ostreatus*, which has 40 published hydrophobin
genes) inflate the count. The cost is measured once on the frozen cutoff. The false-positive review loop is not run in v1.

### 6.4 Hard negatives

Per group, the call rate of Pfam GA, relaxed Pfam and the HMM on the **test** part. A group passes if the call rate is at most the rate fixed in L3 before scoring (default 5% of
the group). The pass rule is per group, not pooled. The HsbA group has no pass threshold (E11); its call rate is reported.

### 6.5 Precision term and the ship rule

Extra calls (relaxed-only or HMM-only) that are **not** T1/T2 training or tuning proteins are clustered across proteomes (the T1/T2 proteins that a cutoff is tuned to recover are labelled and reported in their own column; at 8.5 bits, 4 of the
10 relaxed-only proteome calls are such proteins, per review 2). Each unlabelled cluster gets an outcome from the evidence sheet: `hydrophobin`, `not_hydrophobin`, `unresolved`. Outcomes come from owner decisions (T5) or from a paper (T1).
T5 is used **only** for this ship decision and for the owner's own table column. It is never merged into recall, specificity or any training set (E7).
A level needs at least 5 resolved unlabelled clusters and a Wilson 95% lower bound for precision of at least 0.5. The result is reported with T5 only, with T1 only, and the number of `unresolved` clusters. If fewer than 5 are resolved, the result is
"no evidence that it helps" and the level is not shipped.

**Ship rule (one rule for both levels).** A level (relaxed Pfam, or the HMM) is added to `categories.yaml` only if all of the following hold: the recall rule of 6.2; the cost limit of 6.3 in every one of the 12 proteomes;
the per-group hard-negative rule of 6.4; the precision rule of this section. If the relaxed Pfam level fails any of them, `hydrophobin_extended` is not added, `hydrophobin_relaxed` stays as a reported module, and the numbers are
published with the reason.

### 6.6 Baselines and conditions

The same table shows Pfam GA with and without the R0 and cysteine conditions (the family table rows for the seven models have no second condition), relaxed Pfam, the frozen-spacing rule, and the HMM, on the same proteins.

### 6.7 Secondary test set (LP)

Recall of each method on the Jensen proteins that Pfam GA misses (7 now, minus any of the 8 unverified), reported separately and never mixed into 6.2. These are predictions, so this is a check on consistency, not a measurement.

### 6.8 Leakage and status

The author has seen the 9 Pfam-missed entries, the regex results, the review lists and the baseline numbers. Leave-cluster-out numbers are labelled `partial`, not "honest estimate".
The shipped cutoff and the shipped HMM (if any) are chosen on all clusters, so a status measured on the same proteins is `tuned_on_truth`. This holds for the relaxed Pfam level too. Every status file is `smoke`. No `estimated`.
`calibrate truth --call-status` measures the call on a proteome run. It cannot produce leave-cluster-out numbers. Those come from an analysis script and give no status file.

## 7. Owner adjudication (E7)

For every extra call with no T1/T2 label (relaxed-only, HMM-only, and the 22 unlabelled Pfam calls) the owner receives an evidence sheet: length, cysteine count and spacing, signal peptide,
TMHMM, Pfam, relaxed and HMM scores, BLAST top hits (the 174 known hydrophobins, Swiss-Prot), any paper. The owner chooses `hydrophobin`, `not hydrophobin` or `unsure` and gives a reason.
Decisions are stored in `analysis/hydrophobin_truth/owner_decisions.tsv` (id, tier T5, decision, reason, date). T5 appears in separate columns.

## 8. Expansion (separate project, E4)

After v1 is frozen and measured, scan the Fungi_5k proteomes with v1 and write hits to a catalog labelled `predicted_by_hmm_v1`. A catalog entry becomes a positive only with independent
evidence: a paper with protein-level evidence, an experimental UniProt annotation, or a structural or functional assay. Orthology to a T1/T2 protein does not count by itself, and a BLAST or phmmer hit alone does not.
A v2 model is built only from promoted proteins and is tested on the reserved LP set and on T1 proteins added after the v1 freeze date, which is recorded. The false-positive review loop and any new hard-negative groups
belong to v2, with new proteomes for the cost test. Phylogeny and gene-gain-and-loss analysis use the catalog and are a research analysis, not a tool feature. The scan size is estimated in that project.

## 9. Risks

1. **Small numbers.** 28 T2 clusters; the Pfam-missed proteins are in 6; 3 to 7 positives per measured proteome. Every pass is `smoke`.
2. **Circularity.** Hydrophobins are named for the cysteine pattern. A model measures conformity to known members, not discovery.
3. **T2 is not fully independent of Pfam.** 122 of 131 T2 entries (93%) have a Pfam cross-reference. The model may learn what Pfam already knows. The 9 Pfam-missed entries are the real test and they are 6 clusters.
4. **Conditions.** The R0 signal-peptide condition may exclude real hydrophobins (Lovett 2022, bioRxiv, reports candidates without a predicted signal peptide). The cysteine condition excludes 6 T2 entries with fewer than 8 cysteines (kept by the Pfam branch).
   Both numbers are reported.
5. **Unlabelled true hydrophobins count as cost** (6.3).
6. **Owner labels are judgement** (T5).
7. **Few clusters for a custom model.** About 28 clusters for one family-wide alignment of diverse sequences. A weak HMM is a likely result. That outcome is accepted and reported.

## 10. Work plan

| Task | Content | Output |
|---|---|---|
| L0 | Literature extraction and Jensen ID resolution (done) | `analysis/hydrophobin_truth/literature/` |
| L1 | Verify the 8 Jensen proteins whose stated pattern is not found; check for any T1 candidates; reserve the LP v2 set | tables |
| L1b | Review every T2 FUNCTION or location text for hydrophobin-specific content | list of flagged entries |
| L2 | Define the 28 clusters, the leave-one-cluster-out folds and the fixed negative sample (seeded); record all of it before any scoring | `folds.tsv`, `negatives_fixed.tsv` |
| L3 | Hard-negative sets per group, cluster split into tuning and test; counts; the per-group threshold fixed before scoring | `hard_negatives_*.tsv` |
| L4a | R0 SignalP (same build and mode) on T1, T2, T3, LP and hard-negative sequences; record the identity | R0 table |
| L4 | Re-derive the relaxed-Pfam numbers; implement the nested cutoff procedure; build the HMM per fold; **freeze commit** (procedure, thresholds, negatives) | scripts, frozen commit |
| L5 | Leave-one-cluster-out evaluation: Pfam GA, relaxed Pfam, spacing rule, HMM; hard-negative rates; LP set | tables |
| L6a/L6c | (order changed in the plan) L6a: module, converter and job script, no `categories.yaml` edit. L6c, only after the ship decision: `categories.yaml` edit. Module `hydrophobin_relaxed` and call `hydrophobin_extended` (tests first); the HMM module only if it passes; `categories.yaml` edit only if the ship rule (6.5) passes. Includes: the module-set test, golden files read line by line, the job script added to the standard run set, `params` (cutoff, search options, score type, minimum cysteines, R0 condition), provenance recording the real path of the `current` Pfam link | code, tests, shipped `.hmm` files, provenance |
| L6b | Re-measure the stale `tandem_repeat_protein` call files after the last `categories.yaml` edit | status files |
| L7 | Proteome runs and cost measurement (12 proteomes); status files for `hydrophobin_extended` | report |
| L8 | Evidence sheets; owner T5 decisions | sheets, `owner_decisions.tsv` |
| L9 | Paper notes, hand-off | docs |

Each task is its own commit. The plan and an independent review follow this spec.

## 11. Not in this design

Fungi_5k scan and the evolution study (section 8); a machine-learning classifier; grouped HMMs; changing the Pfam gathering cutoff of `hydrophobin_domain`; merging or pushing PRs; reading the alignment pictures of
Linder 2005 and Sunde 2008 for spacing numbers (they are pictures and were not transcribed).
