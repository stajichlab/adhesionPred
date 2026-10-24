# PR 6 decision note: what to do next with custom models and PANTHER families

*2026-10-24. Author: Claude Code (claude-sonnet-5-5). For the project owner. This note RECOMMENDS and DECIDES NOTHING. No HPCC or SLURM step was run, no script, freeze or result table was changed. Every number was re-read from a committed file (named in brackets). The Ecm33 independent evaluation (merged as #135) is `docs/reports/2026-10-24-pr6-ecm33-independent-eval.md`.*

## Words used in this note (each defined once)

| Word | Plain meaning |
|---|---|
| HMM, profile model | A statistical description of a protein family, built from an alignment of its members. It scores how well a new protein fits the family. |
| Pfam | A public library of family HMMs. Our current "panel" is 39 Pfam models. |
| GA (gathering threshold) | The score cut-off the Pfam curators set for each family. "Pfam GA" means: use the model exactly as Pfam ships it. |
| PANTHER | A second public library of protein families, built independently of Pfam. It is also an HMM library, so agreement with it is agreement between two profile models, not proof of function. |
| report-only | The model is run and its hits are listed in a separate section of the report. It changes no call, tier, status or panel row. |
| Held-out | Proteins kept out of the model's training so they can be used to test it. |
| cluster30 | A way to make a held-out test fairer: remove from training every protein that is at least 30 percent identical (over 50 percent of its length) to a held-out protein, so the model cannot just recognise close relatives. |
| u | In our tests, the number of held-out proteins that every cheap comparison method missed. It is the only place where a custom model could show it adds something. If u is below 3 the test is declared "not testable", because 1 or 2 proteins cannot show anything. |
| jackhmmer, phmmer | Cheap search programs. phmmer searches with one protein. jackhmmer repeats the search three times, adding what it finds to its query. "Seed-only" means it starts from the single Kiss protein and no alignment. |
| Specificity test | Counts how many proteins a family model hits that should NOT be hits (shared negatives and four proteomes). A family that hits everything cannot go into the panel. |
| Clade | A fungal group (for example Sordariomycetes). |

## 1. Where PR 6 stands, in one paragraph

The 39-model Pfam panel misses 44 of the Kiss control clusters (12 adhesion, 28 wall, 4 fusion). We scanned all of Pfam: 2 clusters hit a panel model, 38 hit another Pfam family, 4 hit nothing [pr6_pfam_scan_2026-10-22/scan_per_cluster.tsv]. We then built two custom models for the gaps. Both reached the same place: a custom model works (it finds held-out proteins and beats Pfam), but it does not beat a cheap comparator. Ecm33: the custom model finds 469 of 469 held-out members, but so does an alignment-seeded jackhmmer, and seed-only jackhmmer ties it on the PANTHER test (14 of 15, 185 of 187). hyd1_like: the existing relaxed hydrophobin model already finds most held-out members (19/19, 10/10, 8/10, 5/7, 5/6 external). PANTHER gave an independent set of groupings, and it also shows a few extra families that the Pfam candidate rule left out. The central tension is: what is a custom model worth when a cheap search ties it? The honest answer from the data is "a fixed, fast, versioned artifact" for Ecm33 and "little" for hyd1_like. Meanwhile the biggest untested piece, the specificity test of the 30 candidate Pfam families, has not been run.

## 2. Decided versus open

Decided by the owner (all from the committed records):

| Date | Decision | Source |
|---|---|---|
| 2026-10-22 | Whole-Pfam scan first; custom HMM only where no Pfam hit | spec section 7, decision 1 |
| 2026-10-22 | A custom model is report-only (decision 6); wall tested leave-clade-out, adhesion leave-cluster-out `partial` | spec section 7, decisions 4 and 6 |
| 2026-10-22 | Candidates = Pfam hit covering at least 50 percent OR E-value 1e-30 or better: 26 clusters, 30 families; each gets a specificity test before any panel addition | spec section 7b, S1 and S2 |
| 2026-10-22 | Coverage rule for supported hits = query coverage at least 0.5; aga1 and ece1 left out of the HMM work; PR 6 targets are ecm33 and hyd1_like | drift-check report section 4b |
| 2026-10-23 | The Ecm33 `all` model is integrated report-only, unvalidated | ecm33 build report section 6 |
| 2026-10-24 | KEEP the custom Ecm33 model, report-only, and record that it may need additional development | this task's brief; recorded in `data/custom_models/ecm33_all_2026-10-23/README.md` |

Open:

1. What "additional development" of Ecm33 means and what to do first.
2. hyd1_like: ship, do not ship, or run the inside-clade drift check first.
3. Which of the extra PANTHER families to pursue, and how.
4. Running the specificity test of the 30 candidate Pfam families (not yet run), and whether to add the extra families to it.
5. aga1 and ece1 (issue #121).
6. Whether to stop at Pfam plus the Ecm33 model.

## 3. Every family in play

"Evidence" is how much we actually know: strong, moderate, weak, none. "Cost" is what the next useful step costs.

| Family or cluster | What we know | What a model would add | Evidence | Cost of next step |
|---|---|---|---|---|
| Ecm33 (KS0086) | Pfam PF12454 reaches GA on 28 of 187 PANTHER clusters (48 relaxed); custom model 185 of 187; seed-only jackhmmer also 185 [eval report, table in section 3] | A fixed fast model that covers what Pfam misses | Moderate for "finds family members"; none for "separates adhesion/wall from other proteins beyond the family itself" | Low (see 4.1) |
| hyd1_like (KS0017) | No test unit usable (u at most 2); relaxed hydrophobin model finds most held-out members [pr6-hyd1-like-build.md section 5] | Little beyond the relaxed hydrophobin model | Weak (tiny n, wide intervals) | 0 to 1 small HPCC job |
| aga1 (KS0080) | 11 of 64 jackhmmer hits passed the old looser rule, 0 under the query-coverage rule; matches cover 15 to 31 percent of a protein that is 64 percent S/T/N/Q/G/P [drift-check report] | A non-homology (sequence-feature) approach, not an HMM | Weak | High, new design |
| ece1 (KS0107) | 0 of 126 jackhmmer hits supported: drift confirmed [drift-check report] | Nothing buildable | None | Not worth it now |
| 26 candidate clusters, 30 Pfam families | Pfam already hits them (14 by at least 50 percent cover; 12 by E-value only) [candidate_clusters.tsv]. Specificity not measured | Possible panel additions if specific | Moderate for "hits the member"; none for specificity | 1 HPCC job (4.4) |
| PANTHER versus those 26 | 18 of the 26 clusters have a PANTHER family; the PANTHER family has 1.0 to 1.1 times the fungal proteins of the cluster's Pfam for 16 of them; 1.19 (GH63 cluster), 1.38 (TPT) are broader; 8 clusters (kre1, neg1, iff4, cps1, cot-5, ncl1, ncl2, ncl3) have no PANTHER family [missed44_cluster_panther.tsv, family_summary.tsv] | PANTHER adds almost no new members for these | Moderate | None needed |
| Big1 PTHR28285 (big1, KS0087) | 1,349 fungal proteins, 14 clades, specific. Pfam BIG1_C PF27678 covers 17 percent of the query, E 8e-30, just under the 1e-30 candidate cut; BIG1_N PF27676 covers 45 percent, E 2.5e-16 [scan_per_cluster.tsv] | A wall-synthesis family the panel lacks | Moderate | Low (4.3) |
| KNR4/SMI1 PTHR47432 (cot-2/smi1, KS0006) | 1,432, 16 clades, specific. Pfam SMI1_KNR4 covers 32 percent, E 5e-19 | Same | Moderate | Low |
| DFG10 PTHR14624 (KS0090) | 1,060, 15 clades, labelled "intermediate". Pfam Steroid_dh covers 42 percent, E 7e-09 | Same, but Steroid_dh is a broad enzyme domain | Weak to moderate | Low |
| Chitin Synthase Regulator PTHR46430 (chr4 KS0065, chs4_skt5 KS0084) | 4,935 fungal proteins, 31 clades, labelled specific, but 246 bacterial proteins; Pfam hits are Sel1 repeats (E 6.7e-05 and 4e-08), a repeat that occurs in many unrelated proteins. Needs a manual check | A chitin-synthase regulation family | Weak until checked by hand | Medium (manual) |
| PTHR35192 (cfl1, KS0056) | 1,521 proteins labelled "PROTEIN, PUTATIVE-RELATED", only 5 clades with at least 5 proteins | An adhesion singleton, but the family name says nothing | Weak | Low to medium |
| PMT PTHR10050 (AGOS_AAL050W, KS0137) | 5,658 fungal proteins, labelled intermediate; the same family already enters through PMT PF02366 via two candidate clusters (KS0016, KS0035) | Nothing new | Moderate | None extra |
| irs4 (KS0123, PTHR11216) and opy2 (KS0126, PTHR14596) | irs4 is labelled broad (4,694 proteins, 28 clades), not usable. opy2 has 203 proteins, 10 clades | Nothing now | Weak | Not worth it now |
| Other weak partial clusters: dfi1 (KS0128), ncl4 (KS0148), pga13 (KS0125) | Pfam hit by E-value only above the cut; no PANTHER family | Unknown | None | Not worth it now |

Check on the 12 "weak partial" clusters left out by the candidate rule (38 minus 26): KS0056, KS0006, KS0065, KS0084, KS0087, KS0090, KS0123, KS0125, KS0126, KS0128, KS0137, KS0148. For 9 of them (all except KS0125, KS0128, KS0148) PANTHER gives a family. The panel Pfam does not carry Big1, KNR4/SMI1 or DFG10 (the fungal count of proteins carrying a panel Pfam is 0 for each) [family_summary.tsv].

**Answer to "are you sure nothing else was identified by PANTHER?"** Not quite nothing. For the 26 candidate clusters PANTHER adds essentially nothing (same protein sets). For the clusters the Pfam rule left out, PANTHER does show five specific-looking families (Big1, KNR4/SMI1, DFG10, Chitin Synthase Regulator, PTHR35192) that no panel model carries. Two of these (Chitin Synthase Regulator, DFG10) and PTHR35192 need a manual check before anyone trusts them. One small difference from the brief: the GH30 cluster (neg1) has no PANTHER family in the committed table, so it is not among the "1.2 to 1.9 times bigger" families; the ratios I could verify are 1.19 (GH63) and 1.38 (TPT), and 1.6 for the trehalose phosphatase domain alone (PTHR10788 5,703 versus PF02358 3,535).

## 4. Options per open question

### 4.1 Ecm33: the owner decided to keep it. What could "additional development" mean?

| Option | Pros | Cons | Cost | Risk |
|---|---|---|---|---|
| A. Independent held-out set from curated orthologs (Swiss-Prot or literature, not from a search from Q5AGC4) | Removes the circularity of both tests so far; the only test of "is this Ecm33" | Few curated Ecm33s exist (13 reviewed fungal proteins in PTHR31018 [family_summary.tsv]); tiny n | Mostly local (no HPCC to assemble); one small HPCC job to score | Result may again be "not testable"; wide intervals |
| B. Calibrated cutoff (positive-score-informed, or a second cutoff) | Cutoff 29.6 is set from the false-hit side only, so the margin is untested | The margin is already large (lowest held-out score 86.9 against cutoffs of at most 45.9 [build report]); little to gain | Local | Low |
| C. GA-style provenance (write down exactly how the cut-off was chosen, as Pfam does) and versioning rules | Cheap; makes the model citable and reproducible | Does not add evidence | Local, owner step to approve wording | None |
| D. A fixed second comparator (the jackhmmer recipe written down and pinned, run next to the model) | Shows the custom model is no worse than the cheap method on every new proteome | The two tie, so it adds no new information | One HPCC job | Low |
| E. Do nothing more | No cost | The model stays "unvalidated" in every report | None | Reports keep carrying the label |

Cheap first step (recommended): C now (no cost); then A as a local list-building step so the owner can see how many curated orthologs exist before any job is requested. Defer B and D.

### 4.2 hyd1_like

| Option | Pros | Cons | Cost | Risk |
|---|---|---|---|---|
| Ship report-only like Ecm33 | Consistent with Ecm33; the custom model hits 5 of 6 external curated proteins and the Kiss protein itself (1 of 1) | The relaxed hydrophobin model already recovers most held-out members; two cutoffs sit at the floor (0.1 and 0.2 bits); Sordariomycetes cluster30 is not buildable | Local plus owner sign-off | Adds a second model that duplicates the first; confuses readers |
| Do not ship; keep the record | Honest about the data; no new model to maintain | Any member the relaxed model misses stays missed (it missed 2 of 7 Dothideomycetes, 2 of 10 in one cluster, 1 of 6 external) | None | We might later want it |
| Inside-clade drift check first (16 jackhmmer inside-Eurotiomycetes hits plus 3 no-clade hits) | Settles whether those 19 are real | Does not change the main finding that the relaxed model already works | One small HPCC job (pairwise phmmer, needs sequences that live on HPCC) | Low value |

Recommendation: do not ship now, keep the record, and skip the drift check unless the owner wants to ship. In the Kiss unit the relaxed panel also hits the Kiss protein (1 of 1) [hyd1-like build report, table in section 4].

### 4.3 The extra PANTHER families (Big1, KNR4/SMI1, DFG10, Chitin Synthase Regulator, PTHR35192)

The circularity caveat: PANTHER is itself an HMM, so a PANTHER-based member set tests agreement of two profile models. It cannot prove function. For every family below the specificity test (counting wrong hits) is the part that is independent of PANTHER.

| Option | Pros | Cons | Cost | Risk |
|---|---|---|---|---|
| A. Test their existing Pfam models (BIG1_N PF27676 and BIG1_C PF27678, SMI1_KNR4 PF09346, Steroid_dh PF02544) in the same specificity test | No new model; the Pfam GA lines already exist; Big1, KNR4/SMI1 are labelled specific by PANTHER | Pfam covers only 17 to 45 percent of each query; the test script pins its list of 30 families, so adding needs a change to the pinned list (an owner-approved edit, not made here) | One extra small HPCC job (or fold into the main job if the owner approves the list change first) | Partial coverage may mean a Pfam hit that is a shared domain |
| B. Build a custom HMM from the PANTHER member set | Covers the whole protein | Circular; heavier; more models to maintain | Several HPCC jobs, owner sign-offs | High; ties a cheap method as with Ecm33 |
| C. Manual check first for Chitin Synthase Regulator (Sel1 repeats; bacterial members) and PTHR35192 | Cheap; prevents testing a repeat family | Requires the owner or a curator reading | Owner step | Low |
| D. Leave them out | No cost | Misses Big1/KNR4 which look specific | None | Gap stays |

Recommended order: Big1 first (most specific, 14 clades, its Pfam hit just missed the E-value cut: 8e-30 against 1e-30), KNR4/SMI1 second, DFG10 third, then the manual check of PTHR46430 and PTHR35192. Do A (not B).

### 4.4 The 30 candidate Pfam families: run the specificity test

| Option | Pros | Cons | Cost | Risk |
|---|---|---|---|---|
| A. Run the pinned test as it is (30 families) | Already written, tested and drafted; this was decided on 2026-10-22; nothing can enter the panel without it | Own-member recall is circular (the families were chosen because they hit those members) | One HPCC job, owner go-ahead | Several broad families (for example Glyco_hydro_30, Glyco_tranf_2_3, PMT) may show many wrong hits; that is the point |
| B. Run it after adding the extra PANTHER-linked Pfams | One job instead of two | Needs an owner-approved change to the pinned list and its pin check | One HPCC job plus a small edit | Delay if the edit is contested |
| C. Skip it | None | Breaks the owner rule that no family enters the panel without it | None | No panel additions possible |

### 4.5 aga1 and ece1 (issue #121)

| Option | Pros | Cons | Cost | Risk |
|---|---|---|---|---|
| Start a non-homology (sequence-feature) approach now | Only way to cover Ser/Thr-rich adhesin anchors | New design, overlaps the existing repeat and glycosylation modules, cannot be tested leave-clade-out (family leaves training) | High (new module, new tests) | High; open-ended |
| Rerun the pairwise step to save alignment coordinates, then try a domain-level model of one aga1 segment | Cheap rerun; may find a conserved segment | "Same segment" must be defined before looking; evidence so far is short local matches | One small HPCC job | Segment may be low-complexity |
| Leave as a known gap | No cost | Gap stays | None | None |

Recommendation: leave as a known gap now; revisit after 4.4 (the specificity result shows how much the panel can cover).

### 4.6 The cost of doing nothing more

Stopping at Pfam plus the Ecm33 model leaves: 4 clusters with no Pfam hit (Ecm33 is covered report-only by the custom model; hyd1_like only by the relaxed hydrophobin model; aga1 and ece1 are not covered); 26 candidate clusters never specificity-tested, so none of the 30 candidate families can be added; the five extra PANTHER families unused; Ecm33 labelled "unvalidated". The panel's known misses (44 clusters, 0 of 12 seed-independent adhesion clusters hit [issue #121]) stay. No risk of harm, but the central finding of PR 6 would be "a gap map with no panel change".

## 5. Recommended order of work

| # | Step | HPCC job? | Owner step | Why this order |
|---|---|---|---|---|
| 1 | Write GA-style provenance and versioning for the Ecm33 model (4.1 C); list how many curated Ecm33 orthologs exist (4.1 A) | No | Approve wording | Free; closes the 2026-10-24 decision |
| 2 | Run the specificity test of the 30 candidate Pfam families as pinned (4.4 A) | Yes, 1 job | Go-ahead for the job | The one step that unlocks panel additions; already decided 2026-10-22 |
| 3 | Read the measurement table and sign (or not) the panel additions | No | Sign-off per family | Owner rule: no panel addition without sign-off |
| 4 | Manual check of PTHR46430 and PTHR35192; if the owner agrees, add Big1, KNR4/SMI1 and DFG10 Pfams to a second small specificity run (4.3 A) | Yes, 1 small job | List approval, go-ahead | Cheap, uses existing Pfam models |
| 5 | hyd1_like: keep the record, do not ship (4.2) | No | Confirm | Duplicates an existing model |
| 6 | aga1/ece1: keep as known gap; reopen after step 3 (4.5) | No | None | Needs a new design |
| 7 | Ecm33 development beyond step 1 (curated-ortholog test) only if step 1 shows enough orthologs | Maybe 1 small job | Go-ahead | Avoids a "not testable" third time |

What would change the recommendation: (a) if the specificity test shows most candidate families hit many wrong proteins, steps 4 and 7 become more valuable (a custom, narrower model may be needed); (b) if step 1 finds almost no curated Ecm33 orthologs, drop the independent-set plan and rely on provenance and versioning; (c) if the owner wants a hyd1_like result in the Kiss gap table, ship it report-only; (d) if the manual check shows PTHR46430 is mostly Sel1 repeats, drop it.

## 6. What we do NOT know

- Whether any of the 30 candidate Pfam families is specific (the test has not been run on HPCC) [README_pr6_family_specificity.md].
- Whether the custom Ecm33 model finds Ecm33 beyond what a similarity search finds; both tests so far are circular or PANTHER-based. The PANTHER test is `not testable (u < 3): u = 1`, with only 15 headline clusters (14/15 custom, 95 percent interval 0.681 to 0.998).
- Whether live UniRef50 representatives of 2026-10-23 behave like the 2025_03 release used for training (136 of the 187 ids are not in the 2025_03 file) [eval report section 1].
- Whether PTHR46430, PTHR35192 and DFG10 are wall functions or broader families.
- Whether the 16 inside-clade hyd1_like hits and the 3 no-clade hits are real.
- Whether any of this matters for adhesion: PANTHER agreement and Pfam hits measure family membership, not adhesion.
