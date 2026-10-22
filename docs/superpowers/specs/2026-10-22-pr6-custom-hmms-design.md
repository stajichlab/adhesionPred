# PR 6: custom profile HMMs for what Pfam misses, with leave-clade-out (DRAFT for owner review)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Status: DRAFT. Every decision below is a PROPOSAL with a recommended answer; none is decided. This is the own spec that the plan (section "PR 6") requires before work starts. Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`. Decisions already made by the owner: D8, D9, D19, D20 (2026-10-17 spec) and the Kiss decisions of 2026-10-22 (`docs/reports/2026-10-18-kiss-feasibility.md`, section 10). Numbers below were recomputed from `candidate_groups.tsv` (this folder's `2026-10-22-pr6-custom-hmms-data/`), the Task 5.4 outputs (job 30233470, `recall_by_slice.tsv`, `clusters_missed.tsv`) and the cluster record (`docs/reports/data/kiss_cluster_split_2026-10-22/`).*

## 1. What the Task 5.4 baseline showed

`hmmsearch --cut_ga` of the 39 Pfam models (21 wall and fusion, 18 active family rows; Pfam 38.2) on the 76 positive clusters. A cluster counts as hit when at least one member hits a model of its own slice (`cluster_recall_any`).

| Slice | Seed class (slice models) | Clusters | Hit | 95% interval |
|---|---|---|---|---|
| adhesion | seed-like | 16 | 16 | 0.81 to 1.00 |
| adhesion | seed-independent | 12 | 0 | 0.00 to 0.24 |
| wall | seed-like | 14 | 14 | 0.79 to 1.00 |
| wall | seed-independent | 30 | 2 | 0.02 to 0.21 |
| fusion | seed-independent | 4 | 0 | 0.00 to 0.49 (descriptive) |

Pfam finds what is close to a Pfam seed and almost nothing else. The 44 clusters it misses (adhesion 12, wall 28, fusion 4) are the candidate list. Three of the 28 wall misses (pir1, ecm33, cfmA) are `seed-like` for all 39 models (close to an adhesion-family seed) but missed by the wall models.

**Caveats that carry into this design.**
- These are family coverage numbers, not the sorting-hat call's sensitivity.
- `seed-like` is close to circular (a cluster at 30% identity to a seed is found by that model almost by construction). Only the `seed-independent` rows measure discovery.
- The interval is a Wilson interval on clusters taken as independent. It is a proxy.
- A "miss" is a miss of the 39-model panel. It is not yet a miss of Pfam (section 2, step 6.0).

## 2. A miss of the panel is not yet a gap in Pfam

The panel holds 39 chosen models. Several missed groups are families that may have well-known Pfam domains (O-mannosyltransferases, kexin-like proteases, alpha-glucan synthases, glycosidases). If so, the fix is to add a Pfam model to the panel, not to build a custom HMM. That is cheaper, has a Pfam GA cutoff already, and needs no leave-clade-out (the Pfam model is not trained on our data). I have not checked this: no search against the whole of Pfam-A was run.

**Step 6.0 (proposal, first):** `hmmsearch --cut_ga` of the 44 missed clusters' members against all of Pfam-A 38.2, and report per cluster which families hit (none, one, several). Only the clusters with no Pfam hit at all stay custom-HMM candidates. This could shrink the list a lot before any HMM is built.

## 3. The 44 missed clusters grouped as candidate targets

The groups are **my grouping by the Kiss annotation text** (`candidate_groups.tsv`, column `kiss_annotation`), not by homology. A group of unrelated proteins cannot share one HMM. Section 4 turns groups into homologous families (step 6.1). Columns: clusters = missed MMseqs2 clusters; shared = clusters whose Kiss family also has members in another clade (the others are `clade_only`, decision 3); max/clade = the most missed clusters of the group in one fine clade.

| Group | Slice | Clusters | Clades (fine) | Shared | Max/clade | Genes |
|---|---|---|---|---|---|---|
| A1 hydrophobin-like | adhesion | 2 | Eurotiomycetes, Ustilaginomycetes | 0 | 1 | hyd1_like, rep1 |
| A2 yeast surface adhesin-like | adhesion | 4 | Saccharomycotina | 0 | 4 | aga1, dfi1, ece1, iff4 |
| A3 Neurospora lectins | adhesion | 4 | Sordariomycetes | 0 | 4 | ncl1 to ncl4 |
| A4 adhesion singletons | adhesion | 2 | Eurotiomycetes, Tremellomycetes | 0 | 1 | cfl1, glfB |
| F1 pheromone receptors | fusion | 3 | Eurotiomycetes, Saccharomycotina, Ustilaginomycetes | 2 | 2 | pra1, preA, ste2, ste3 |
| F2 so (fusion) | fusion | 1 | Sordariomycetes | 0 | 1 | so_soft |
| W1 glycosylation and glycan supply | wall | 7 | 6 clades (Chytridiomycota to Ustilaginomycetes) | 4 | 5 | pmt1, pmt2, pmt4, cot-5, cwh41, dfg10, gfa1 |
| W2 chitin synthase regulation | wall | 5 | Saccharomycotina, Sordariomycetes, Taphrinomycotina | 1 | 4 | chs4/skt5, chs5, chs6, chs7, chr4 |
| W3 beta-1,6-glucan (Kre family) | wall | 5 | Saccharomycotina, Sordariomycetes | 0 | 4 | kre1, kre5, kre6, big1, neg1 |
| W4 alpha-1,3-glucan and other synthases | wall | 2 | Eurotiomycetes, Sordariomycetes | 0 | 1 | ags1, ags2, cps1 |
| W5 GPI wall proteins and signalling | wall | 6 | Eurotiomycetes, Saccharomycotina | 0 | 5 | pir1, ecm33, cfmA, irs4, pga13, opy2 |
| W6 processing and other | wall | 3 | Eurotiomycetes, Saccharomycotina, Sordariomycetes | 3 | 3 | kex2, tps2/orlA, cot-2/smi1 |

Missed clusters by stratum: wall 8 shared and 20 clade-only; adhesion 0 shared and 12 clade-only; fusion 2 shared and 2 clade-only.

**What is feasible under the owner's floor (5 clusters, fine level primary).** The floor was set on all positive clusters. Here the test units are only the *missed* clusters, so the counts are smaller:

| Held-out clade | Slice | Missed clusters | Shared | Clade-only | Reaches floor 5 |
|---|---|---|---|---|---|
| Saccharomycotina | wall | 21 | 8 | 13 | all, shared and clade-only |
| Eurotiomycetes | wall | 8 | 5 | 3 | all and shared; clade-only descriptive |
| Sordariomycetes | wall | 7 | 4 | 3 | all only; both strata descriptive |
| Saccharomycotina | adhesion | 4 | 0 | 4 | no |
| Sordariomycetes | adhesion | 4 | 0 | 4 | no |

(A cluster with members in two clades counts once per clade, so the wall rows add to more than 28.) Findings:

1. **Per group, no held-out clade reaches 5 clusters except where a group is large and spans one clade** (W1 and W5 in Saccharomycotina have 5). A per-group leave-clade-out is not an estimate anywhere. Per-group results are `descriptive`.
2. **Per slice (all custom models of the slice together), only wall is testable:** Saccharomycotina, Eurotiomycetes and (all-clusters only) Sordariomycetes. Adhesion has at most 4 missed clusters in any clade, so it is `descriptive` at the clade level. Fusion is `descriptive` (decision 5).
3. **Clade-only clusters cannot be recovered by transfer.** If every Kiss member of a family is in the held-out clade (all 12 adhesion misses, 20 of 28 wall misses), removing that clade removes the whole family from training. Recall 0 is then by design; it measures discovery of a family never seen, and only a homolog from outside Kiss (section 4) can change it. For clade-restricted genes (A2, A3: the Candida and *Saccharomyces* adhesin-like genes; the *Neurospora* lectins) there may be no homolog outside the clade at all. These stay out of reach of a homology model, and need a different kind of module (composition, repeats, domain architecture), which is outside this PR.
4. **The shared wall clusters are where leave-clade-out is informative:** 8 clusters (4 in W1, 1 in W2, 3 in W6). These are the families that Kiss shows in two or more clades, such as the PMT O-mannosyltransferase cluster (9 members in 5 clades) and the kexin and gfa1 clusters.

## 4. Design (all PROPOSAL)

**4.1 Unit of model.** A homologous family, not a function group. Step 6.1 clusters the members of the 44 clusters (and later their outside homologs) by a sensitive search (jackhmmer or MMseqs2 profile search, E-value stated before looking) and writes `families.tsv`. For example, whether ncl1 to ncl4 are one family is unknown now: they are four MMseqs2 clusters at 30% identity.

**4.2 Unit of evaluation.** The held-out clade, per slice, with every custom model of the slice searched together. Reported per stratum (shared, clade-only), with n and status by decisions 3 and 4. Per-group and per-family numbers are descriptive columns. Where a slice has no clade at the floor (adhesion), leave-cluster-out is reported, labelled `partial`, never `estimated`.

**4.3 Training alignments.** Three sources, in order of trust:
1. Kiss members of the family outside the held-out unit (the control set, accessions from `data/controls/wall-fusion-kiss2019/`).
2. Outside orthologs found by sequence search (reciprocal best hit or a profile search) in UniProt reviewed fungal proteins (the second source of the feasibility report, section 7; Swiss-Prot gave wall_remodeling 55 entries in Saccharomycotina and 11 in Basidiomycota, wall_synthesis 18 in Basidiomycota, and 0 in Chytridiomycota). That source is queried by Pfam family, so it does not reach a family without a Pfam model; a sequence search is needed. UniProt FASTA staged to HPCC would be a new input (public, hash recorded).
3. The proteomes of the refresh (`scripts/sorting_hat/refresh_3_10/inputs.tsv`: nine run, five staged but not run), as FASTA files under the older HPCC work directories. The 2026-10-17 staging directory `$NEW` was deleted on 2026-10-22, so anything staged from it would be rebuilt.

Family representatives (the 35 *A. fumigatus* orthologs) may be used as training seeds (decision 6). They are never test members and never counted toward a floor.

**4.4 Independence rule.** For each held-out unit (a clade, or a cluster for leave-cluster-out) the build alignment contains no sequence that is the unit's member, and none that clusters with it at 30% identity and 50% coverage (the same MMseqs2 rule as Task 5.3), and no Kiss family-mate in the held-out clade. The build writes a manifest (training accessions with cluster ids and clades, aligner and HMMER versions, commands, file hashes). A test checks that the held-out accessions and their 30% neighbours are absent from the manifest (the precedent: `analysis/hydrophobin_truth/lco.py`).

**4.5 Build and cutoff.** Align with MAFFT L-INS-i (the precedent used it after `--auto` failed a pre-registered column check; availability on HPCC to be checked), `hmmbuild`, then a gathering cutoff set from the training clades only. Proposed rule (hydrophobin precedent, section 6.1): the lowest bit score at or above a floor at which, on proteomes and shared negatives that contain no held-out data, the unlabelled extra calls are at most a stated rate per 10,000 proteins (5 in the precedent), and the call rate on the shared hard negatives (161 rows; 6 `#neg` copies of overlap accessions are in mixed clusters) is at most a stated rate. No positive's score sets a cutoff. A freeze commit fixes the rule, the aligner command line, the cutoff floor, the held-out units and the thresholds before any held-out score is looked at.

**4.6 Reporting per held-out unit.** For each unit: recovered or not per missed cluster, the score and the cutoff, the nearest panel model, and the Pfam baseline on the same unit. The comparison is paired: a custom model counts as adding recall only for clusters that the best baseline misses (let u be that number). **Lesson of the hydrophobin run:** u was 0, so the custom HMM was "not testable on the primary set". Here u is large by construction (all 44), but after step 6.0 and a relaxed-GA Pfam baseline it may shrink; the minimum u for a testable result and the rule "at least ceil(u/2) recovered" are fixed before the freeze.

**4.7 False hits.** Each custom model is scored on the shared negatives (per group, never pooled) and on the wall false-hit panel for its slice, descriptively (wall rows are never negatives for the adhesion call, D20). A model that hits an adhesion family's negatives above the rate limit is reported and not shipped.

## 5. Status and tier consequences

- A custom model gets at most `smoke` until its leave-clade-out (or leave-cluster-out) measurement is committed (D9 step 3). `estimated` needs the floor of decision 4 on its own stratum.
- The tier derivation (`src/cellsurface_sorting_hat/family_tier.py`) knows two bases today (`descriptive_shared_sets`, `per_family_reviewed`) and has no hold-out field. The plan's test "the derivation refuses `estimated` for a custom model without the hold-out field" therefore needs a code change first: a new basis (for example `held_out`) and columns (`holdout_unit`, `holdout_kind`, `n_clusters_held_out`) in `family_measure.tsv`. PROPOSAL: PR 6a adds these with tests first.
- `narrow_ok` is untouched. Any narrow-call use of a custom model is an owner decision (D4), made per model after the numbers exist.
- The wall slice has no call: wall models are report-only (`wall_fusion_models.tsv`). A custom wall model would join that file or a sibling, not `family_table.tsv`, unless the owner adds a call for it. Where a custom model lives is an open question (section 7, Q6).

## 6. Staged plan (HPCC jobs need the owner's go-ahead each time)

Conventions: SLURM partition `short`, `${SCRATCH:?}`, no `/tmp`, pin check on the code commit and inputs, scripts committed unedited after they run, results in their own PR, tests first, findings in an issue rather than edits of a run record. Heavy work never on the login node.

| Step | What | Where | Output | Stop point |
|---|---|---|---|---|
| 6.0 | Whole-Pfam scan of the members of the 44 missed clusters | HPCC job (one, minutes) | per-cluster Pfam hits | owner reads it; decides which clusters stay candidates and whether to add Pfam models to the panel |
| 6.1 | Homology grouping of the remaining members, and a Pfam-free family list | HPCC job (MMseqs2/jackhmmer) | `families.tsv` | owner reads the families |
| 6.2 | Assemble training sets (Kiss, UniProt reviewed fungal proteins, proteomes) and the independence manifests | local + HPCC | manifests, hashes | review of the manifests; freeze commit |
| 6.3 | Per held-out unit: MAFFT, `hmmbuild`, cutoff on training data (SLURM array) | HPCC | models, cutoffs | none until the freeze is committed |
| 6.4 | Score held-out units; baselines on the same units; false hits | HPCC | tables | owner reads them |
| 6.5 | Report and `family_measure.tsv` rows with the new basis | repo | PR | owner decides `narrow_ok` and any call |

PR order: **6a** this spec reviewed (grill) and merged; **6b** tests first, the tier-derivation change (new basis and columns), step 6.0 and 6.1 scripts; **6c** results of 6.0 and 6.1 (record PR); **6d** the freeze commit; **6e** builds and held-out scoring; **6f** report. If step 6.0 empties the list, the PR ends at 6c with a proposal to add Pfam models to the panel.

## 7. Open questions for the owner (ordered by dependency; each with my recommendation)

1. **Run the whole-Pfam scan (6.0) before anything is built?** Recommend yes. It is one cheap job and may remove candidates.
2. **What is the unit of a model?** (a) homologous family from step 6.1, (b) my function groups, (c) one HMM per missed cluster. Recommend (a); the groups of section 3 are not homology.
3. **What is the unit of evaluation?** (a) held-out clade per slice with all models together, (b) per model. Recommend (a), with per-model rows as descriptive columns.
4. **Which slices are attempted?** Recommend wall (3 clades reach the floor: Saccharomycotina, Eurotiomycetes, Sordariomycetes on all clusters); adhesion as leave-cluster-out `partial` only (at most 4 missed clusters per clade); fusion out (decision 5 stands; 4 clusters).
5. **Clade-only clusters (adhesion 12, wall 20) have no Kiss relative outside their clade.** Do we try to find outside homologs, or accept that homology models cannot reach them? Recommend try with the UniProt search in 6.2, report how many have any outside homolog, and leave the rest to a non-homology module outside this PR.
6. **Where does a custom model live?** (a) report-only (wall models file or a new sibling), (b) a family row in `family_table.tsv` with tier at most `smoke` and `narrow_ok=no`, (c) a call. Recommend (a) for wall, and (b) only for an adhesion family once there is a leave-cluster-out number; (c) is a later decision per model.
7. **Minimum training set per model.** Provisional proposal: at least 5 training sequences in at least 2 clades outside the held-out unit (as the floor of 5 is provisional too). Recommend adopt as provisional.
8. **Cutoff rule and limits.** Adopt the hydrophobin precedent (at most 5 unlabelled extra calls per 10,000 proteins; hard-negative rate fixed before scoring)? Recommend yes, stated numbers fixed at the freeze.
9. **Comparators beyond Pfam GA.** Add (a) a relaxed-GA Pfam baseline (the `hydrophobin_relaxed` precedent) and (b) a jackhmmer or phmmer iteration seeded by the held-in members. Recommend both; a custom HMM must beat the cheap alternatives to be worth shipping.
10. **Leave-species-out inside a clade** for clade-only clusters (for example hold out *Candida* and train on *Saccharomyces* within Saccharomycotina). Recommend include as `descriptive` only, never `estimated` (species are not independent clades).
11. **ESM models stay out** (issue #15 gate of 300 clusters; the control set has 76 positive clusters). Recommend unchanged.
12. **Rows to re-check in the Kiss slice.** `glfB` (a UDP-galactofuranose transporter) is in the adhesion category; the `ncl` lectins and `cot-2/smi1` rest on the Kiss annotation text only; seven of the 44 clusters are `portal_validated` and four contain a `locus_confirmed` member (all accepted by decision 0). Recommend keep them and keep the flags; the freeze lists them so a result can be re-run without them.

## 8. Not in this design

No ESM or fine-tuned model. No change to any call, `narrow_ok`, tier or status file by this document. No HPCC job has been run for PR 6. The grouping table is a proposal made from annotation text; it is not evidence of homology.

## 9. Files

`2026-10-22-pr6-custom-hmms-data/candidate_groups.tsv`: the 44 missed clusters with slice, my group, genes, clades, stratum, Kiss family number, both seed classes, `jgi_accept` and the Kiss annotation (ids and text only; no sequences). Test: `tests/wall_fusion_truth/test_pr6_candidate_groups.py`.
