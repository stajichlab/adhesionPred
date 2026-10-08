# Custom hydrophobin HMM and the `hydrophobin_extended` call: design

*2026-10-08. Revision 1. Author: Claude Code (claude-sonnet-5-5). Builds on
`2026-10-08-hydrophobin-validation-design.md` (rev 3) and the measurement in
`docs/reports/2026-10-08-hydrophobin-validation.md`. The owner agreed the eight decisions in section 2 in a
question-by-question interview on 2026-10-08. Status: draft for independent review. No code is written.*

## 1. Why

1. The Pfam hydrophobin models miss known hydrophobins. In Swiss-Prot, 9 of 174 named entries (5.2%) have no hydrophobin-class
   Pfam cross-reference, and all 9 have experimental (T2) evidence. In the 8 measured proteomes, Pfam finds 36 of 42 labelled
   hydrophobins. Five of the 9 Swiss-Prot misses have a weaker hit to the Eas or Hydrophobin models at E-values from 2e-8 to 1e-3,
   and one of them, PSH_FLAVE, is 0.9 bits below the domain gathering cutoff.
2. The 8-cysteine spacing rule built from published spacings recovers 1 of the 6 labelled proteins that Pfam misses in those
   proteomes. It matches 60% of T2 hydrophobins. It is too narrow to use as a rescue.
3. A profile built from the known hydrophobins can capture family-level detail that a spacing range cannot. This is a recall gain
   on known members. It is not evidence of discovery, because hydrophobins are partly defined by the cysteine pattern (spec rev 3,
   section 4.1).

Prior art exists for pattern-plus-domain searches (Yang 2006, PMID 17217508; Jensen 2010, PMID 21182770). This work does not claim
novelty of the approach. It adds a measured call with a status inside the tool.

## 2. Decisions agreed with the owner (2026-10-08)

| # | Decision |
|---|---|
| E1 | Optimise for recall, with **two evidence levels**. `hydrophobin_domain` (Pfam) stays unchanged. A separate `hydrophobin_extended` call adds the custom HMM. |
| E2 | Train on T1 and T2 proteins only. Add hydrophobins from the literature as T1, after each is verified in the paper, with the PMID, when they are not curated in Swiss-Prot. T3 entries are held out and used only as a separate check. |
| E3 | Success criterion (section 6). All results are `smoke`. Leave-cluster-out testing shows that the model is consistent without memorising. |
| E4 | Order: freeze v1 and measure it; scan Fungi_5k with v1 as a **separate project**; predicted members are never truth until independent evidence; a v2 model is tested on held-out literature proteins that v1 never saw. |
| E5 | The HMM ships as a versioned data file `data/sorting_hat/hydrophobin_ext.hmm` with a provenance file. A module `hydrophobin_ext` runs `hmmsearch` with a frozen score cutoff. |
| E6 | Hard negatives: CFEM proteins, cerato-platanin proteins, HsbA proteins, small secreted cysteine-rich proteins (killer toxins, effector candidates), and cell wall proteins with cysteines (PIR, Ccw12-like). Other look-alikes are added from the model's own false positives. |
| E7 | The owner's adjudication counts as a separate tier "owner-reviewed" (T5), recorded with a reason per protein. It is reported beside the experimental tiers, never merged. |
| E8 | The owner supplied Linder 2005, Sunde 2008, Wessels 1994, the supplements of Jensen 2010 and Xu 2021 (`to_import/hydrophobins/`). Merge order of PRs #76, #75 and this branch stays with the owner. |

A machine-learning classifier was considered. With about 130 experimentally supported proteins in about 30 clusters, it is likely to
overfit or to relearn the cysteine pattern. It is not part of this design. This is an expectation, not a result.

## 3. Training data

### 3.1 Tiers

| Tier | Meaning | Training | Measurement |
|---|---|---|---|
| T1 | Hydrophobin identity or function shown in a paper at protein level (purification, mass spectrometry, rodlet or film formation, contact angle, or a hydrophobin-specific deletion phenotype). PMID recorded. | yes | primary |
| T2 | Swiss-Prot reviewed, experimental evidence (ECO:0000269) on a FUNCTION, SUBCELLULAR LOCATION or SUBUNIT comment | yes | primary |
| T3 | Named hydrophobin by rule or similarity only | **no** | separate check |
| T4 | A review judgement from sequence and BLAST with no paper | no | reported only |
| T5 | Owner-reviewed (E7) | not in v1 training | separate column, never merged |

Counts now: T2 131, T3 43 (`analysis/hydrophobin_truth/truth_all.tsv`). No T1 entries yet.

### 3.2 Literature expansion

Sources supplied: Jensen 2010 supplements (hydrophobins of nine *Aspergillus* genomes) and Xu 2021 supplements (40 *P. ostreatus*
genes). The extraction is in `analysis/hydrophobin_truth/literature/` (task L0). Rules:
1. A literature protein becomes T1 only if the paper gives protein-level evidence for that protein. A gene prediction, an EST or an
   expression profile is not T1. Such proteins are listed as "literature, predicted" and held out like T3.
2. Each T1 entry records the PMID, the table or figure, and the accession or the sequence. No accession is written from memory.
3. A protein already in Swiss-Prot is not added twice.
4. The number of proteins the literature adds, and how many clusters they open, is measured and reported before training. It is unknown now.

### 3.3 Clusters and splits

MMseqs2 at 30% identity and 0.5 coverage on all T1 and T2 proteins, as before. The development/test split of clusters is made once,
with a seed, before any model is built, and stored. Leave-cluster-out groups (section 6) are made from the same clusters.

## 4. The model

### 4.1 Build

1. Align the training proteins per cluster, then align cluster representatives (a profile-profile or progressive alignment;
   the tool is chosen in the plan after checking what is installed). Alignment quality is checked by eye on the cysteine columns.
2. `hmmbuild` from the alignment. Two architectures are compared on the **development** clusters only: (a) one family-wide HMM,
   (b) a small set of HMMs by group (class I, class II, divergent/unclassified), called as "any". The winner is chosen on the
   development part and frozen.
3. The cutoff is a **bit score**, not an E-value, so it does not depend on the proteome size. It is chosen on the development
   clusters and frozen.
4. The profile, the alignment, the training list with cluster IDs and tiers, the HMMER version, the build command and the file hash
   go in `data/sorting_hat/hydrophobin_ext.provenance.json`.

### 4.2 Module and call

`hydrophobin_ext` (module, one row per protein): `hit` = best HMM score at or above the frozen cutoff **and** at least 8 cysteines
**and** the R0 call is `called`. The R0 identity is in the module params, as for `cys8_pattern`. Fields also include `score`,
`n_cys`, `length` and `pattern_match` (from the frozen spacing, for comparison only). A row has an empty `hit` when it would be a hit
except that the R0 row is missing.

```yaml
- name: hydrophobin_extended
  expr:
    or: [{flag: pfam_hydrophobin.hit}, {flag: hydrophobin_ext.hit}]
```

The strict call `hydrophobin_domain` is unchanged. The extended call is a superset: every Pfam call is also an extended call. The report shows
which protein is called by Pfam, by the HMM only, or by both. `hydrophobin_extended` is ungated, reads two modules and has no `ref`, so it can get a
call status file and can be a mechanism of the `other_*` calls (placed before them).

`cys8_pattern` and its spacing file stay in the tool as a reported side column. They are not part of any default call.

## 5. Hard negatives and look-alikes

Groups (E6): CFEM; cerato-platanin; HsbA; small secreted cysteine-rich proteins (killer toxins, effector candidates); cell wall
proteins with cysteines (PIR, Ccw12-like). Sources: Swiss-Prot reviewed fungal entries per group, and Pfam family membership for
the first three and the last (the Pfam models of those families are used only to build the negative list). For each group the report gives
the HMM call rate. A cutoff choice is judged against the criterion in section 6, not against one pooled number.

Other look-alikes (LysM effectors, fungal defensins and antifungal proteins, expansin-like proteins) were named from general knowledge and are
not verified. The data-driven route: after the model is built, review what it calls wrongly in the measured proteomes. A family that
keeps appearing joins the hard-negative set, and the cutoff is re-checked on the development clusters.

## 6. Evaluation and success criterion

1. **Leave-cluster-out (LCO).** Train on all clusters but one group, score the held-out group, repeat until every cluster has been
   held out once. Groups are formed so that every group has at least one positive cluster. The number of groups is set in the plan from the cluster count.
2. **Recall target.** In the LCO runs the HMM (at the frozen cutoff, with the signal-peptide and cysteine conditions) recovers at least
   half of the labelled hydrophobins that Pfam misses. The Pfam-missed set is the 9 Swiss-Prot entries plus any T1 literature proteins
   that Pfam misses (counted when the literature is added).
3. **Cost limit.** At the frozen cutoff, the calls that Pfam does not make number at most 5 per proteome of about 10,000 proteins, in each measured proteome.
   Each is reviewed by hand with the evidence sheet (section 7).
4. **If the criterion fails**, `hydrophobin_extended` is not added to `categories.yaml`. The module and the profile remain, reported only.
5. **Reported with every number:** the counts of positives and clusters, the interval, that negatives are assumed, that the numbers are small,
   and per hard-negative group the call rate.
6. **Leakage.** The shipped model is trained on all clusters, so a status measured on the same proteins is `tuned_on_truth`. The LCO numbers are the honest
   estimate of generalisation and are reported as such. Every status file is `smoke`. No `estimated`.
7. **Pre-registration.** The cutoff, the architecture and this criterion are written to a commit before any held-out or proteome result is looked at.
8. **Baselines in the same table:** Pfam alone, the frozen-spacing rule, and the HMM, on the same proteins.

## 7. Owner adjudication (E7)

For every called protein with no T1/T2 label (strict Pfam calls, HMM-only calls) the owner receives an evidence sheet: sequence length, cysteine count and spacing, signal peptide, TMHMM,
Pfam and HMM scores, BLAST top hits (the 174 known hydrophobins and Swiss-Prot), any paper found. The owner chooses `hydrophobin`, `not hydrophobin` or `unsure`
and gives a reason. The decisions are stored in `analysis/hydrophobin_truth/owner_decisions.tsv` (id, tier T5, decision, reason, date), as in the repeat-call work.
T5 labels appear in separate columns of every table. They are not merged with T1/T2 and are not used to train v1.

## 8. Expansion (separate project, E4)

After v1 is frozen and measured, scan the Fungi_5k proteomes with v1 and write the hits to a catalog labelled `predicted_by_hmm_v1`. A catalog entry becomes
a positive only with independent evidence (a paper, an experimental UniProt annotation, or orthology to a T1/T2 protein by a method that is not the HMM). A v2 model
is built only from promoted proteins and is tested on held-out literature proteins that v1 never saw. Phylogeny and gene-gain-and-loss analysis use the catalog and are a research
analysis, not a tool feature. Size of the scan is not estimated here. It is estimated in that project.

## 9. Risks

1. **Small numbers.** About 130 T2 proteins in 30 clusters, 3 to 7 per measured proteome. The LCO test has coarse steps. A pass is `smoke`.
2. **Circularity.** Hydrophobins are named for the cysteine pattern. The model measures conformity to known members. It says little about unknown hydrophobins.
3. **Tiers are not fully independent of Pfam.** T2 is based on UniProt experimental evidence, but the entries were collected by name. Some are in Swiss-Prot because a Pfam-based rule fired. T3 is
   held out for this reason, and the 122 T2 entries with a Pfam cross-reference are 93% of the T2 set (122 of 131). The model may learn what Pfam already knows. The 9 no-Pfam entries and any literature proteins
   Pfam misses are the real test, and there are few of them.
4. **Alignment of a diverse family.** Class I and class II share little sequence beyond the cysteines. A single alignment may be poor. The architecture comparison (4.1) addresses this.
5. **Signal-peptide condition.** Lovett 2022 (bioRxiv) reports hydrophobin candidates without a predicted signal peptide. The condition is kept (D10 of the earlier spec). The number of labelled
   hydrophobins without an R0 call is reported.
6. **Cysteine condition.** At least 8 cysteines excludes hydrophobins with 6 or 7 (HYD2B_BEAB2 has 6, HFBA_PENEN has 7 in the full sequence). The number of labelled proteins lost to this condition is reported.
7. **Owner labels are judgement.** T5 is reported separately for that reason.

## 10. Work plan

| Task | Content | Output |
|---|---|---|
| L0 | Extract spacing statements and protein lists from the supplied papers (running) | `analysis/hydrophobin_truth/literature/` |
| L1 | Verify literature proteins; add T1 entries with PMID; report counts and clusters; merge into the truth table | `manual_entries.tsv`, updated `truth_all.tsv` |
| L2 | Re-cluster; make and store the development/test split and the LCO groups before any model | `clusters.tsv`, `split_hmm.tsv.gz`, `lco_groups.tsv` |
| L3 | Hard-negative sets per group (E6) | `hard_negatives_*.tsv` |
| L4 | Alignment and `hmmbuild`; architecture comparison and cutoff on development clusters; write the choices to a commit (pre-registration) | `models/`, a frozen commit |
| L5 | LCO evaluation against the criterion; baselines | table |
| L6 | Module `hydrophobin_ext` and call `hydrophobin_extended` (tests first); `categories.yaml` edit only if the criterion passes | code, tests |
| L7 | Runs on the measured proteomes; measurement; per-group call rates | status files, report |
| L8 | Evidence sheets for the owner; T5 decisions | sheets, `owner_decisions.tsv` |
| L9 | Paper notes, hand-off | docs |

Each task is its own commit. The plan and an independent review follow this spec.

## 11. Not in this design

Fungi_5k scan and evolution study (section 8); a machine-learning classifier; changing the Pfam cutoffs; merging or pushing PRs; Wessels 1994, Linder 2005 and Sunde 2008 beyond the statements the extraction records.
