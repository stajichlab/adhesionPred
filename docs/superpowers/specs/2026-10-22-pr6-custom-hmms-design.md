# PR 6: custom profile HMMs for what Pfam misses, with leave-clade-out

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Status: owner decisions recorded (section 7, final 2026-10-22); the body follows them. This is the own spec that the plan (section "PR 6") requires before work starts. Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`. Decisions already made by the owner: D8, D9, D19, D20 (2026-10-17 spec) and the Kiss decisions of 2026-10-22 (`docs/reports/2026-10-18-kiss-feasibility.md`, section 10). Numbers below were recomputed from `candidate_groups.tsv` (this folder's `2026-10-22-pr6-custom-hmms-data/`), the Task 5.4 outputs (job 30233470, `recall_by_slice.tsv`, `clusters_missed.tsv`) and the cluster record (`docs/reports/data/kiss_cluster_split_2026-10-22/`).*

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

**Step 6.0 (decision 1, first):** `hmmsearch --cut_ga` of the members of the 44 missed clusters against all of Pfam-A 38.2, reporting per cluster which families hit (none, one, several). Only a cluster with no Pfam hit at all becomes a custom-HMM target. A cluster that hits is dropped from the custom list and flagged: the proposal for it is to add that Pfam model to the panel (a separate, later decision per model). The fusion clusters are scanned too (the scan is one job and the result is descriptive); no custom fusion model is built (decision 4).

## 3. The 44 missed clusters grouped as candidate targets

The groups are **my grouping by the Kiss annotation text** (`candidate_groups.tsv`, column `kiss_annotation`), not by homology. By decision 2 they are reporting labels only: a model is built per homologous family found by sequence search (step 6.1), never per group. Columns: clusters = missed MMseqs2 clusters; shared = clusters whose Kiss family also has members in another clade (the others are `clade_only`, decision 3); max/clade = the most missed clusters of the group in one fine clade.

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

## 4. Design (follows the owner decisions of section 7)

**4.1 Unit of model (decision 2).** One model per homologous family. Step 6.1 clusters the members of the clusters that survive step 6.0 by a sensitive sequence search (jackhmmer or MMseqs2 profile search; the E-value is written before the results are read) and writes `families.tsv`. Function groups are labels in the report. Whether ncl1 to ncl4 are one family is unknown now: they are four MMseqs2 clusters at 30% identity.

**4.2 Unit of evaluation (decisions 3, 4, 10).** The held-out clade, per slice, with every custom model of the slice searched together. Per-model and per-family numbers are descriptive columns.
- **Wall:** leave-clade-out at the fine level (primary), with the coarse level as the second test where eligible. Reported per stratum (shared, clade-only), with n and status by decisions 3 and 4 of the Kiss decisions.
- **Adhesion:** leave-cluster-out, reported as `partial`, never `estimated` (at most 4 missed clusters in any clade).
- **Fusion:** out (no custom model). Fusion clusters appear only in the step 6.0 scan, descriptively.
- **Leave-species-out inside a clade** (for example hold out *Candida* and train on *Saccharomyces*): `descriptive` only.

**4.3 Training sequences (decisions 5, 7).** In order of trust:
1. Kiss members of the family outside the held-out unit (the control set, `data/controls/wall-fusion-kiss2019/`).
2. Outside homologs found by `jackhmmer` or `phmmer` against UniProt (reviewed fungal proteins are the first choice; the exact set, its download date and hash are fixed at step 6.2). The search is by sequence, so it also reaches families without a Pfam model, which the Pfam-keyed second source of the feasibility report (section 7 there) does not. Step 6.2 reports **how many clade-only clusters (adhesion 12, wall 20 among the 44) have any outside homolog**. The ones with none are left to a non-homology module outside this PR.
3. Family representatives (the 35 *A. fumigatus* orthologs) may be training seeds (Kiss decision 6); never test members, never counted toward a floor.

**Minimum training set (decision 7, provisional like the floor):** at least 5 sequences in at least 2 clades outside the held-out unit. A unit that cannot meet it is recorded as `not_buildable` and reported, not built.

**4.4 Independence rule.** For each held-out unit (a clade, or a cluster for leave-cluster-out) the build alignment contains no member of the unit, none that clusters with it at 30% identity and 50% coverage (the Task 5.3 rule), and no Kiss family-mate in the held-out clade. The build writes a manifest (training accessions with cluster ids and clades, aligner and HMMER versions, commands, file hashes). A test checks that the held-out accessions and their 30% neighbours are absent from the manifest (precedent: `analysis/hydrophobin_truth/lco.py`).

**4.5 Build and cutoff (decision 8).** Align with MAFFT L-INS-i (the precedent used it after `--auto` failed a pre-registered column check; availability on HPCC to be checked), `hmmbuild`, then a cutoff set from the training clades only, by the hydrophobin rule: the lowest bit score at or above a floor at which the unlabelled extra calls on the false-hit proteomes (4.8) are at most **5 per 10,000 proteins** and the call rate on the shared hard negatives (161 rows; 6 `#neg` copies of overlap accessions sit in mixed clusters) is at most a rate fixed in advance. Over the four kept proteomes (32,491 proteins) the extra-call limit is 16.2, so **about 16 extra calls in all**. Both limits are written at the freeze commit, before any held-out score is looked at. No positive's score sets a cutoff. The freeze also fixes the aligner command line, the cutoff floor, the held-out units, the minimum training set and the thresholds of 4.6, and lists the weak-evidence rows by name (decision 12).

**4.6 Reporting per held-out unit (decision 9).** For each unit: recovered or not per missed cluster, the score and the cutoff, the nearest panel model, and the comparators on the same unit:
- Pfam GA (the Task 5.4 baseline);
- a relaxed-GA Pfam baseline (the `hydrophobin_relaxed` precedent: the same models with a frozen lower cutoff);
- a `jackhmmer` or `phmmer` iteration seeded by the held-in members.
The comparison is paired: a custom model counts as adding recall only for clusters the best comparator misses (let u be that number). **Lesson of the hydrophobin run:** u was 0 and the custom HMM was "not testable on the primary set". The minimum u for a testable result, and the rule "at least ceil(u/2) recovered", are fixed at the freeze.

**4.7 False hits.** Each custom model is scored on the shared negatives (per group, never pooled) and on the wall false-hit panel for its slice, descriptively (wall rows are never negatives for the adhesion call, D20). A model that hits above the fixed rate is reported and not shipped.

**4.8 Proteome set (decision 9b).** One proteome per species in any analysis or testing, from now on. Proteomes with internal-ID headers are removed (`Afum_Af293_Fungi5k`, `Bder_ER3`). The kept proteomes are `Afum_Af293_UniProt` (9,647 proteins), `Calb_SC5314` (6,212), `Cimm_RS` (9,910) and `Scer_S288C` (6,722): 32,491 in all (`scripts/sorting_hat/refresh_3_10/inputs.tsv`). The other proteomes of that file (`Afum_A1163`, `Afum_W72310`, `Scer_S288C_nodubious`, and the two removed) are not used. The 2026-10-17 run records and the 2026-10-22 re-run records are history and are not changed; the rule applies to what is run from now on, including the PR 6 false-hit checks. The proteome split used for tuning and testing the cutoff (hydrophobin precedent: tuning and test halves, seeded and stored before scoring) has to fit four species; how many tuning proteomes remain is fixed at the freeze.

## 5. Status and tier consequences

- A custom model gets at most `smoke` until its leave-clade-out (wall) or leave-cluster-out (adhesion) measurement is committed (D9 step 3). `estimated` needs the floor of Kiss decision 4 on its own stratum; an adhesion number is `partial`.
- The tier derivation (`src/cellsurface_sorting_hat/family_tier.py`) knows two bases today (`descriptive_shared_sets`, `per_family_reviewed`) and has no hold-out field. The plan's test "the derivation refuses `estimated` for a custom model without the hold-out field" therefore needs a code change first: a new basis (for example `held_out`) and columns (`holdout_unit`, `holdout_kind`, `n_clusters_held_out`) in `family_measure.tsv`. PR 6b does this, tests first.
- **Where a model lives (decision 6):** a custom model is **report-only** (the wall models file or a sibling). A row in `family_table.tsv` is possible only for an adhesion family that has a leave-cluster-out number, and that is a separate owner decision per family. No call is added by this PR.
- `narrow_ok` is untouched. Any narrow-call use is an owner decision per model (D4).

## 6. Staged plan (each HPCC job needs the owner's go-ahead)

Conventions: SLURM partition `short`, `${SCRATCH:?}`, no `/tmp`, a pin check on the code commit and the inputs, scripts committed unedited after they run, results in their own PR, tests first, findings in an issue and not in edits of a run record. Heavy work is never on the login node. All false-hit and tuning runs use the four kept proteomes (4.8).

| Step | What | Where | Output | Stop point |
|---|---|---|---|---|
| 6.0 | One job: `hmmsearch --cut_ga` of the members of the 44 missed clusters against all of Pfam-A 38.2 (`/srv/projects/db/pfam/current`, release check as in Task 5.4) | HPCC job (minutes) | per-cluster Pfam hits | owner reads it; clusters that hit are dropped and flagged for the panel; the rest are the custom targets |
| 6.1 | Homology grouping of the surviving clusters' members by sequence search | HPCC job | `families.tsv` | owner reads the families |
| 6.2 | Outside homologs by `jackhmmer` or `phmmer` against UniProt (set and hash fixed here); report how many clade-only clusters have any; apply the minimum training set (5 sequences in 2 clades); assemble the manifests; mark `not_buildable` units | HPCC job plus local staging of the UniProt set | manifests, counts | review; the freeze commit (limits written, weak-evidence rows listed by name) |
| 6.3 | Per held-out unit: MAFFT, `hmmbuild`, cutoff from training data and the four kept proteomes (SLURM array) | HPCC | models, cutoffs | none before the freeze is committed |
| 6.4 | Score the held-out units; the three comparators (Pfam GA, relaxed-GA Pfam, `jackhmmer` or `phmmer` iteration) on the same units; false hits on the four proteomes and the shared negatives | HPCC | tables | owner reads them |
| 6.5 | Report, and `family_measure.tsv` rows with the new basis | repo | PR | owner decides any family row, `narrow_ok` or call |

PR order: **6a** this spec merged; **6b** tests first, the tier-derivation change, the step 6.0 and 6.1 scripts; **6c** results of 6.0 and 6.1 (record PR); **6d** the freeze commit; **6e** builds and held-out scoring; **6f** report. If step 6.0 leaves no cluster without a Pfam hit, the PR ends at 6c with the list of Pfam models proposed for the panel.

## 7. Owner decisions (final, 2026-10-22)

| # | Decision | Where it is in the body |
|---|---|---|
| 1 | Whole-Pfam scan first (step 6.0). Only clusters with no Pfam hit become custom-HMM targets; clusters that hit are dropped or flagged for adding that Pfam model to the panel. | 2, 6 |
| 2 | One model per homologous family found by sequence search (step 6.1). Function groups are reporting labels only. | 3, 4.1 |
| 3 | Evaluate per held-out clade per slice with all models of the slice together; per-model rows are descriptive columns. | 4.2 |
| 4 | Wall: leave-clade-out. Adhesion: leave-cluster-out, reported as `partial`. Fusion stays out. | 4.2, 5 |
| 5 | Search for outside homologs in step 6.2 (`jackhmmer` or `phmmer` against UniProt); report how many clade-only clusters have any; the rest are left to a non-homology module outside this PR. | 4.3, 6 |
| 6 | A custom model is report-only (wall models file or a sibling). A family row only for an adhesion family with a leave-cluster-out number, as a separate decision. | 5 |
| 7 | Minimum training set: at least 5 sequences in at least 2 clades outside the held-out unit (provisional). A unit that cannot meet it is `not_buildable`. | 4.3 |
| 8 | Hydrophobin precedent: at most 5 unlabelled extra calls per 10,000 proteins and a fixed hard-negative rate, both written at the freeze. Over the kept proteomes (32,491 proteins) the extra-call limit is about 16. | 4.5 |
| 9 | Comparators: a relaxed-GA Pfam baseline (`hydrophobin_relaxed` precedent) and a `jackhmmer` or `phmmer` iteration seeded by the held-in members. | 4.6 |
| 9b | One proteome per species in any analysis or testing. Proteomes with internal-ID headers are removed (`Afum_Af293_Fungi5k`, `Bder_ER3`); kept: `Afum_Af293_UniProt`, `Calb_SC5314`, `Cimm_RS`, `Scer_S288C`. The 2026-10-17 run records are history and unchanged; the rule applies from now on, including the PR 6 false-hit checks. | 4.8 |
| 10 | Leave-species-out inside a clade: `descriptive` only. | 4.2 |
| 11 | No ESM models (issue #15, gate of 300 clusters, stands). | 8 |
| 12 | Keep `glfB`, the `ncl` lectins, `cot-2/smi1` and the weak-evidence rows with their flags. The freeze lists them by name so that a result can be re-run without them. | 4.5 |

## 7b. Owner decisions after the step 6.0 scan (final, 2026-10-22)

Step 6.0 (job 30235946, record in `docs/reports/data/pr6_pfam_scan_2026-10-22/`, note `docs/reports/2026-10-22-pr6-pfam-scan.md`) found: of the 44 missed clusters, 2 hit a panel model, 38 hit another Pfam family and 4 have no Pfam hit. The owner decided how to use this:

| # | Decision |
|---|---|
| S1 | The candidates for the next step are the union of a Pfam hit covering at least 50 percent of the sequence OR a best hit with an E-value of 1e-30 or better (26 clusters: the 14 with at least 50 percent coverage plus 12 partial hits; the 12 weaker partial hits stay descriptive). Coverage alone is a poor proxy, because a short domain in a long protein is normal. |
| S2 | Every candidate family gets a specificity test before any panel addition: hits over the shared negatives and over the four kept proteomes, recorded in a measurement table with a script-derived tier (at most `smoke`). Nothing enters the active panel without a `family_measure.tsv` row and the owner's sign-off; `narrow_ok` stays `no` until signed. |
| S3 | The outside-homolog search (steps 6.1 and 6.2) runs for the 4 gaps (hyd1_like KS0017, aga1 KS0080, ece1 KS0107, ecm33 KS0086), starting with the ecm33 check (Ecm33 PF12454 is in the panel yet gives no GA hit on Q5AGC4). A model is built only where the training rule of decision 7 is met; the rest are reported as `not_buildable` known gaps. |
| S4 | The two pieces are drafted in parallel; the HPCC jobs run one at a time with the owner's go-ahead for each: the small gap search first, the specificity tests after. |

The 14 "add to panel" proposals of the scan come from the same candidate list; they need no separate decision. The two wall clusters hit by adhesion family models (pir1 KS0082, cfmA KS0116) stay descriptive.

## 8. Not in this design

No ESM or fine-tuned model (decision 11). No fusion model (decision 4). No change to any call, `narrow_ok`, tier or status file by this document. No non-homology module (composition, repeats, domain architecture) for the clade-restricted genes; that is a separate piece of work. No HPCC job has been run for PR 6. The grouping table is a reporting label made from annotation text; it is not evidence of homology.

## 9. Files

`2026-10-22-pr6-custom-hmms-data/candidate_groups.tsv`: the 44 missed clusters with slice, my group, genes, clades, stratum, Kiss family number, both seed classes, `jgi_accept` and the Kiss annotation (ids and text only; no sequences). Test: `tests/wall_fusion_truth/test_pr6_candidate_groups.py`.
