# PR 6: hyd1_like build and leave-cluster-out design (DRAFT for the owner)

*2026-10-23. Author: Claude Code (claude-sonnet-5-5). Nothing here was run on HPCC or submitted to SLURM. Binding decisions: `2026-10-22-pr6-custom-hmms-design.md` sections 4 to 7b; `docs/reports/2026-10-22-pr6-drift-check.md` section 4b (the rule is query coverage >= 0.5; hyd1_like stays a target); `docs/reports/2026-10-22-pr6-gap-search-uniref50.md`. Precedent: the ecm33 run (`docs/reports/2026-10-23-pr6-ecm33-build.md`, `analysis/wall_fusion_truth/README_pr6_ecm33.md`) and the hydrophobin run (`analysis/hydrophobin_truth/`). Everything marked PROPOSAL is for the owner; this document decides nothing.*

## 1. What the ecm33 run taught, and what this design does about it

| Lesson | Design response |
|---|---|
| Every comparator tied the custom model; u = 0 in all 11 units, so `u >= 3` made 0 units testable. The iteration comparator was seeded by the same training alignment. | Comparators independent of the training alignment (section 5). The one comparator that reaches the members by construction (jackhmmer from the Kiss protein) is a reported ceiling and never enters u. |
| The held-out members came from the same search that supplied training. | Section 6: an independent held-out set is proposed from the committed hydrophobin truth. It is an owner decision (Q4). Without it every cell stays circular and `descriptive`. |
| A true family member was among the shared negatives (Q4WNS8), found only by a pre-check. | A pre-freeze negatives job (section 8, step 1) runs before the freeze, over the 161 negatives and the four kept proteomes, with the family rule and the hydrophobin-class Pfam panel. |
| Fixed E-values over a tiny target file look better than in the original search. | Seed comparators use `-Z` = UniRef50 sequence count (70,198,728, round 2 provenance). |

## 2. The member set (recounted from the committed TSVs)

`pr6_hyd1_members.py` reads `docs/reports/data/pr6_drift_check_2026-10-22/drift_hits.tsv` and `hits_manifest.tsv`. A hit is a member when its class is `supported_homolog`, its **query** coverage is >= 0.5 (owner decision 4b) and it has a named fine clade. Counts (the unit tests recompute them from the files):

| Item | n |
|---|---|
| Drift hits of A0ACE4KFX7 (jackhmmer outside fungal hits) | 43 (32 `supported_homolog`, 11 `unsupported`) |
| Supported with query coverage >= 0.5 | 32 (lowest query coverage 0.531) |
| ... in a named fine clade | **29**: Sordariomycetes 19, Dothideomycetes 7, Leotiomycetes 3 |
| ... without a fine clade (`(no clade)`, listed in `excluded.tsv`, never a clade, not counted) | **3**: UniRef50_A0A1J9P1Q5 (E 2.3e-32, query cover 0.99), UniRef50_A0A2J6TR40 (1.4e-21, 0.98), UniRef50_A0A9P4ULB3 (8.4e-19, 0.857); all taxid 716546 |
| Query itself (A0ACE4KFX7, 98 residues, Eurotiomycetes, `locus_confirmed`) | 1 |
| **Members** | **30** = query + 29 hits, in 4 named clades |
| Reached by one-pass phmmer (`in_phmmer = yes`) | 3 of the 29 named members, 2 of the 3 no-clade hits; the other 27 of the 32 only by iteration |

Two things to see. (1) The strongest outside homologs (E 1e-32 and 1e-21, near full length) are the no-clade hits, which the rule excludes. (2) Inside-clade (Eurotiomycetes) hits exist (16 UniRef50 ids in the round 2 table: 8 by phmmer at E 4.4e-47 to 9.3e-6, 8 more by jackhmmer) but the drift check covered only outside hits, so they have no query-coverage check and are **not** members (Q3).

A member is a search hit that passed the drift check. It is not a curated hydrophobin and not a proven hyd1_like family member. The Kiss annotation of KS0017 is "hydrophobin"; no hydrophobin Pfam model hits the query at GA (`scan_per_cluster.tsv`: `no_pfam_hit`; PF22354 Eas appears only in the weak table, E up to 1e-3).

## 3. Units and what n means

**n** is a count of UniRef50 representative sequences in a unit (members), never a count of independent clusters. The Kiss cluster KS0017 holds one protein, so the number of held-out Kiss clusters in this experiment is **1**. Members of one clade are related. Therefore **every cell is `descriptive`**; the only label the committed rules allow above that is `partial` (decision 4, adhesion leave-cluster-out), which I attach to the Kiss-unit row only. No cell is `estimated`; there is no interval.

Units (`pr6_hyd1_members.py`, `pr6_hyd1_build.py cluster-units`):

1. **The Kiss unit** `clade:Eurotiomycetes` = hold out A0ACE4KFX7 (the one Kiss protein; its clade holds no other member). Train on the 29 outside homologs. This is the leave-cluster-out of the one Kiss cluster, the number the adhesion rule asks for. n = 1: a single-protein line (`kiss_single`), never a rule unit (u <= 1 < min_u).
2. **Clade units**: each named clade of the members in turn. Status `unit` needs >= 5 members and the training rule met.
3. **Cluster units**: each MMseqs2 30 percent cluster (30 percent identity, 50 percent coverage) of the 30 members in turn; the held-out set is the whole cluster. These cannot be listed before the clustering (needs the sequences; step 2 of the plan). One column only: the held-out cluster is already removed, so `cluster30` equals `none`. Most clusters will be small; units need >= 5 members, so few or none may be testable. That is a result.

**Training rule** (decision 7, provisional): at least 5 sequences in at least 2 named clades outside the held-out unit, applied AFTER the independence rule. Two readings, both computed (Q2):

| Held-out unit | members | outside seq / clades | `total` (the ecm33 code) | `per_clade` (>= 2 clades each with >= 5) | Status, `total` | Status, `per_clade` |
|---|---|---|---|---|---|---|
| Sordariomycetes | 19 | 11 / 3 (Doth 7, Leo 3, query 1) | met | not met (only Doth >= 5) | unit | not_buildable |
| Dothideomycetes | 7 | 23 / 3 (Sord 19, Leo 3, query 1) | met | not met (only Sord >= 5) | unit | not_buildable |
| Leotiomycetes | 3 | 27 / 3 | met | met | too_few_members | too_few_members |
| Eurotiomycetes (Kiss unit) | 1 | 29 / 3 (Sord 19, Doth 7, Leo 3) | met | met (Sord, Doth) | kiss_single | kiss_single |

These are the counts before the 30 percent rule. The independence rule can only lower them; whether Sordariomycetes (training of 11) survives cluster30 is unknown until the clustering runs. Under `total`, the Sordariomycetes training set is 11 sequences, of which 1 is the query: small for an HMM.

The Kiss unit is the only unit that is not circular in the sense of ecm33: the held-out protein was not found by the search from other clades; the training members were found *from it* by search, so a model of them recovering it is still tied to that search. This is stated, not fixed (Q4).

## 4. Independence rule

For each unit: no held-out member, and in the headline mode no training member in the same 30 percent cluster as a held-out member, is in the build alignment (spec 4.4; `check-alignment` enforces it on the alignment ids, and a test checks the manifest). Headline `cluster30`; second descriptive column `none` (clade units and the Kiss unit; Q5). The Kiss unit is where this bites: the five one-pass hits are close homologs of the query and may well share a 30 percent cluster with it, so cluster30 may leave only the most distant ones.

## 5. Comparators (independent of the training alignment)

Computed on the held-out members of each unit with the same `-Z` for the seed searches. `u_indep` = held-out members missed by every comparator in the first three rows. The custom model must recover ceil(u/2) of u, and a unit is testable only if u >= `min_u` (3, proposed, carried over from ecm33).

| Comparator | What | Independent of training? |
|---|---|---|
| `pfam_ga` | the eight hydrophobin-class Pfam models in the panel (Hydrophobin PF01185, Hydrophobin_2 PF06766, Eas PF22354, DewD PF28987, Hyd1F PF29802, Hydrophobin_D PF29785, Hydrophobin_like PF29465, HsbA PF12296), `hmmsearch --cut_ga` | yes |
| `relaxed_panel` | `data/sorting_hat/hydrophobin_relaxed.hmm` (the seven hydrophobin-class models at the frozen 2.6-bit cutoff, `--nobias --cut_ga`, committed, chosen before this work) | yes |
| `phmmer_seed` | one-pass phmmer of the single Kiss protein over the held-out members | yes; not run in the Kiss unit (seed = held-out member) |
| `jackhmmer_seed` (ceiling) | `jackhmmer -N 3` from the single Kiss protein | **no**: the members were found by that iteration (27 of 32 only by it), so it recovers them by construction. Reported, never in u |

What u would be. I did not run any search. Reasoning from the record: the query has no Pfam GA hit, so `pfam_ga` is expected to recover few or none of the members; `phmmer_seed` by construction misses the members that only iteration reached (27 of the 32 supported hits; among the 29 named members 26), so for Sordariomycetes (16 of 19 iteration-only) u_indep could be large, up to 16, and for Dothideomycetes (7 of 7 iteration-only) up to 7 (minus whatever the Pfam panel recovers). If so these two units would be testable under `min_u` 3, unlike all 11 ecm33 units. This is an expectation to be tested, not a result. The cost: u large partly because `phmmer_seed` is a weak comparator (single pass from a 98-residue seed); the freeze names it as a cheap alternative, not as the best possible search, and the report says so. With the ceiling included u would be 0 (jackhmmer recovers the members), which is the ecm33 outcome; the report prints both `u_indep` and `u_with_ceiling` so the owner sees the gap.

## 6. The hold-out truth question (OWNER DECISION, Q4)

Within the member set the truth is search-defined (circular). The repo has independent hydrophobin material:

| Candidate set | Where | n | Comment |
|---|---|---|---|
| Swiss-Prot hydrophobins without a hydrophobin Pfam cross-reference | `analysis/hydrophobin_truth/nopfam9.faa`, `truth_all.tsv` (tier T2) | 9 | curated (UniProt evidence codes), Pfam-missed: the same blind spot as KS0017. Small. |
| Hydrophobin truth proteins T2/T3 | `analysis/hydrophobin_truth/positives_t2t3.faa`, `truth_all.tsv` | 174 named (131 T2, 43 T3) | curated, but almost all PF01185-type, so Pfam already covers them; informative for the model only if they are related to the Kiss protein |
| The 8 other Kiss hydrophobin positives | `data/controls/wall-fusion-kiss2019/controls.tsv` (A0A0D1BZL3, Q7S3P5, P79073, I1RDW3, A0ACE4JU10, P41746, E9QT94, Q04571) | 8 | Kiss-curated, other clusters; Pfam-covered |

The honest question is not "does the model hit hydrophobins" but "does a model of the Kiss protein's outside homologs recover curated hydrophobins that no Pfam hydrophobin model sees". Whether any curated hydrophobin is related to the Kiss protein is unknown; the pre-check job writes a relatedness table (phmmer of A0ACE4KFX7 over these sets, E <= 1e-3, descriptive) so the owner can decide with numbers. **Decided 2026-10-23 (section 9, item 4):** the set is the related truth proteins minus training members (n = 6). The original proposal was: (i) the nopfam9 set plus whichever truth proteins the relatedness table shows within E 1e-3 become an `external` descriptive line, (ii) every training member in a 30 percent cluster with an external protein is excluded from that evaluation's model (the same independence rule), (iii) the set is named and frozen before scoring. Implemented as the unit `external`; the freeze holds the ids.

## 7. Freeze contents (PROPOSALS; file `analysis/wall_fusion_truth/pr6_hyd1_freeze_PROPOSED.json`, `frozen: false`)

Decided by the owner already (not reopened): 5 extra calls per 10,000 (16.2455 over the 32,491 proteins of the four kept proteomes), the kept proteomes, adhesion leave-cluster-out reported `partial`, custom models report-only. Everything else is a proposal: `cutoff_floor` 0.0 (hydrophobin floor); `hard_negative_rate` 0.05 (hydrophobin value); `min_unit_members` 5; `training_rule.reading` `total`; independence `cluster30` headline and `none` second column; `min_u` 3; `negatives_ignored` empty until the pre-check job runs (candidate by reading controls.tsv only: B6QUR1, MP1, Pfam HsbA; the job decides); `proteome_other_hydrophobin` `ignore`; alignment check `min_columns` 40 (not the ecm33 100: the query is 98 residues and 31 of 32 hits are 61 to 197 residues, 100 occupied columns is out of reach) with occupancy 0.5 and median coverage 0.6; no cysteine-column check proposed (KS0017 has not been checked for the 8 cysteines); `evalue_db_size` 70198728; the comparator list of section 5; `heldout_truth.external_set` null.

## 8. Staged HPCC plan (nothing run; EACH step needs the owner's go-ahead, D21)

| Step | Job | Needs | Output | Owner reads |
|---|---|---|---|---|
| 0 | this PR merged; HPCC clone updated to the merge commit | go-ahead | | |
| 1 | DONE 2026-10-23 (job 30321490): `pr6_hyd1_negatives_precheck.sbatch` | | family and panel tables, relatedness | results in section 9 |
| 2 | `pr6_hyd1_extract.sbatch` (after the freeze is signed; it reads the freeze for the reading and the external set): members, FASTA, 30 percent clusters, `units_cluster.tsv`, `cysteines.tsv`, `external.faa`, `ext_cluster.tsv`. Scores nothing | go-ahead | members, units, cluster units | the unit list, how small the clusters are |
| 2b | optional: an inside-clade drift check (query coverage of the 16 Eurotiomycetes hits). A script is not drafted | Q3 | | |
| 3 | the freeze commit (the owner only) | decisions of section 9 | `pr6_hyd1_freeze.json`, `frozen: true` | |
| 4 | `pr6_hyd1_build.sbatch` array (tasks from `pr6_hyd1_build.py tasks`): manifest, MAFFT, `hmmbuild`, cutoff, only then the held-out searches | go-ahead | `<unit>/<mode>/` tables | |
| 5 | `pr6_hyd1_report.py` locally; a record PR; a report | | `recall_by_unit.tsv` | Q-after: ship or not |

No call, status, tier, `narrow_ok`, `family_measure.tsv` or `family_table.tsv` row is written by any of this.

## 9. Owner decisions of 2026-10-23 (all recommendations accepted) and open items

Pre-check result (job 30321490, code 5fb6b48, read-only outputs in `adhesionPred-ecm33-results/hyd1_negcheck/`): 0 family members among the 161 shared negatives (phmmer of A0ACE4KFX7, E <= 1e-5); 0 proteome proteins called family by phmmer; 1 hydrophobin-class Pfam hit among the negatives (B6QUR1, MP1, *T. marneffei*, stratum laccase_etc, HsbA GA score 168.7, not called family by phmmer); 19 hit rows among the four kept proteomes (18 distinct proteins: Afum_Af293 16, Cimm_RS 2; E9QT94 hits two models), none called family by phmmer; 7 curated hydrophobin truth proteins within E 1e-3 of the Kiss protein (I1RDP9 2.3e-10, G9MKB2 2.7e-10, I1RDQ0 3.7e-9, A0A6V8R0V1 2e-8, Q8X1T7 1.3e-5, A0A067P943 7.3e-4, Q9Y8F0 8.1e-4). The query has no Pfam GA hit.

Decided (recorded in `pr6_hyd1_freeze_PROPOSED.json` under `decided_by_owner_2026-10-23`; the file stays `frozen: false` until the owner signs):

1. **Build hyd1_like at all:** yes, report-only (as recommended; not separately asked, unchanged).
2. **Training rule reading = `total`** (the ecm33 code), stated as the weaker reading; the `per_clade` counts are still printed in `units.tsv`.
3. **Members** = the query + the 29 drift-supported hits with a named clade = 30 members in 4 named clades. The 16 inside-clade Eurotiomycetes hits and the 3 supported no-clade hits stay OUT for now.
4. **External held-out set (descriptive, tiny):** the curated hydrophobin truth proteins within E <= 1e-3 of the Kiss protein, MINUS training members. Recomputed from the committed drift table: A0A6V8R0V1 is UniRef50_A0A6V8R0V1 (Sordariomycetes, supported, query cover 0.551), so it is a member, stays in the training set and is removed from the external set. **External set n = 6:** I1RDP9, G9MKB2, I1RDQ0, Q8X1T7, A0A067P943, Q9Y8F0. Scored as the unit `external` (one column, cluster30: training members in a 30 percent cluster with an external protein are removed). n = 6 is a count of proteins, not of independent clusters; no estimate.
5. **Independence columns:** cluster30 headline, `none` second descriptive column.
6. **Comparators:** the alignment-independent set of section 5 (eight hydrophobin-class Pfam models at GA, the relaxed panel, one-pass phmmer of the Kiss protein skipped in the Kiss unit); jackhmmer from the single protein only as a ceiling outside u; `min_u` 3; the Kiss unit (n = 1) is never a rule unit.
7. **Hard-negative rate 0.05, `negatives_ignored` = []** (no family member among the shared negatives).
8. **`ignore`:** hydrophobin-class panel hits that phmmer does not call family are left out of the extra-call count and the hard-negative rate, **listed by id in the freeze** (`panel_ignored`: 18 proteome proteins + B6QUR1 = 19 distinct ids). The report also prints the extra-call count and the hard-negative rate with those ids counted (columns `*_panel_counted`, at the headline cutoff, and `cutoff_if_panel_counted`).
9. **Alignment check:** `min_columns` 40 (the query is 98 residues); no cysteine check as a stop condition; the cysteine count of the 30 members is written to `cysteines.tsv` by the extract job as a descriptive note.
10. **`-Z` 70198728** on the seed searches.

Reporting rule added with the decisions: every row of the report (headline units, descriptive column, external line, Kiss unit) prints the exact two-sided 95 percent Clopper-Pearson interval for each recall and for recovered-of-u (the function of `pr6_ecm33_eval_report.py`, imported). A wide interval is not a result.

Still open: (a) the inside-clade drift check of the 16 Eurotiomycetes hits (step 2b, no script drafted); (b) whether the 3 no-clade hits are ever used as training only; (c) `cutoff_floor` 0.0, `min_unit_members` 5, alignment `occupancy` 0.5 and `min_median_coverage` 0.6 remain proposals; (d) the freeze signature itself.

## 10. What this does not show

The members are homologs by iterated search, not curated hyd1_like. Clades and clusters are related sequences. One Kiss cluster is held out, so no estimate exists. A custom model that beats the independent comparators on these units would still not show adhesion function, and it changes no call. The scripts are untested on HPCC; `esl-reformat`-free (the jackhmmer comparator here takes the single query sequence), but the options of `hmmfetch -f`, `mmseqs easy-cluster` and `jackhmmer -Z` are unconfirmed on the cluster.
