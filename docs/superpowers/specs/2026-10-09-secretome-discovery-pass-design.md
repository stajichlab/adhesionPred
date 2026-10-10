# B: Discovery pass on the secretome (Onygenales and Aspergillaceae)

*2026-10-09. Revision 2 (after `...-secretome-discovery-pass-design-review-1.md`: 1 blocker, 9 major, 9 minor; all addressed below). Author: Claude Code (claude-sonnet-5-5). Part of the tier restructuring agreed with the owner on 2026-10-09 (order B, then C, then A). Status: draft for independent review. No code is written. Numbers come from files named in the text or from the review's recomputation.*

## 1. The question

In the scan of 56 Fungi_5k proteomes (`docs/reports/2026-10-09-fungi-scan-descriptive.md`), 97% of the proteins with a signal peptide have no mechanism evidence (sum of `other_surface_no_mechanism[R0]` over sum of `signal_peptide_protein[R0]` = 0.967 in `proteome_summary.tsv`; per-proteome medians: Onygenales 412.5 of 427, *Aspergillus* 784 of 808.5, medians computed from `proteome_summary.tsv`; the report table shows upper medians for the orders). The owner wants to know whether the Onygenales and the Aspergillaceae carry patterns in this pool that suggest additional candidates for (1) interaction with the host, (2) protection from immune recognition, or (3) general cell surface functions.

This is a descriptive discovery pass. It makes no call and does not change `categories.yaml`. The output is a table of families and features that differ between the target clades and a background, with candidate lists. **Clade enrichment cannot tell host interaction apart from other lineage traits.** The pass states associations, not selection and not function.

## 2. Owner decisions that fix the design (2026-10-09)

1. The unit of analysis is the **genus**. A pattern is reported only if it appears in **at least 3 independent genera** of the target clade.
2. The background is larger than the 56 scan proteomes. The owner accepts that phylogeny confounds the result.
3. For *Aspergillus*: Aspergillaceae is the main target; a within-*Aspergillus* comparison is a descriptive supplement.
4. Tiers: broad (cell wall and surface), narrow (adhesion), host interaction. This pass feeds the host-interaction tier and the broad tier.

Changes in rev 2 that the owner has not yet seen (decisions O1 to O3 in section 8): the background becomes all passing genera (about 600) instead of a draw of 150; a tree-based count of origins is reported beside each p-value; precomputed Fungi_5k annotations are reused where a check shows they agree.

## 3. Data

Counts from `/bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv` (5,813 proteomes). Genus counts below are **named genera**: rows with an empty GENUS are excluded (Onygenales has 3 such proteomes, "Onygenales sp. PD_10", "PD_12", "PD_40"; Aspergillaceae has 2).

| Group | Proteomes | Distinct GENUS values (empty string included) | Named genera |
|---|---|---|---|
| Onygenales (target 1) | 71 | 28 | 27 |
| Aspergillaceae (target 2) | 406 | 8 | 7 |
| Eurotiales (includes Aspergillaceae) | 483 | 15 | 14 |
| Pezizomycotina | 2,802 | 654 | not recomputed |

The number of named genera that **pass the selection rule** is 26 or fewer for Onygenales and 7 or fewer for Aspergillaceae (review recomputation); the exact numbers come from B1 and are listed with the failing genera.

**Selection rule** (three parts, all required): (1) the proteome FASTA exists in `input/`; (2) at least 4,000 proteins (`MIN_PROTEINS` in `select_proteomes.py`); (3) `NCBI_TAXONID` is in `_workdir/sorting_hat/taxdump/nodes.dmp` (the tool needs it). *Paracoccidioides* (taxids 502780 and 502779; both proteomes have more than 4,000 proteins) fails part 3. B1 recovers it through `merged.dmp` or the species taxid and records the mapping. Deterministic choice within a genus: the proteome with the highest BUSCO completeness if the file in Fungi_5k gives it, otherwise the highest protein count; BUSCO or protein count is recorded for every chosen proteome and shown beside the results (a poor gene model removes signal peptides).

- **Targets:** every named Onygenales genus and every named Aspergillaceae genus that passes the rule.
- **Background:** **every** named genus that passes the rule in Pezizomycotina outside Onygenales and Aspergillaceae (the review counts 603), including the other Eurotiales genera. Using all of them removes the seed dependence. Stability is checked in the analysis by 5 random subsamples of 150 genera, stratified by class with allocation proportional to class size and capped at 40% of the draw.
- **Supplement (5.3):** up to 40 *Aspergillus* species proteomes, spaced across the species tree (section 5.3).

## 4. Existing results in Fungi_5k, and what is run

`/bigdata/stajichlab/shared/projects/Fungi_5k/results/function/` holds, for all 5,813 proteomes: `pfam/*.domtblout.gz` (full Pfam-A, `hmmsearch --cut_ga`; release 37.x, not 38.2), `predgpi/*.predgpi.gff3.gz`, `signalp/` (SignalP 6), `tmhmm/`, `cazy/`, `merops/`, `targetP/`, `wolfpsort/`, `aiupred/`, `repeatmasker/`. `Phylogeny/fungi_tree/final_tree.nw` is a species tree. `results/Fungi5K_cluster.tsv.gz`, `orthogroups_mmseqs.tsv.gz` and `mmseqs_nopfam_genedump.*` hold MMseqs2 clusters, including clusters with no Pfam hit.

Decisions:

- **Family models needed by the tool:** the 24 models of `family_table.tsv` are searched again with Pfam 38.2 (several are new in 38: PF28987, PF29785, PF29802, PF29465, PF30910). The 56 scan proteomes already have these tables.
- **SignalP 6, TMHMM:** reuse the Fungi_5k results only after B2a shows that they give the same `signal_peptide_protein[R0]` and TM calls as the scan run on the 56 scan proteomes (same mode and organism setting). If they differ, the pass runs its own `signalp_gpu.sbatch` and `tmhmm` steps.
- **Full Pfam-A and PredGPI:** reuse the Fungi_5k results. The report states that the release is Pfam 37.x for these tables and that the 24 family models come from 38.2.
- **New computations:** repeat02, repeat14, allergen BLAST and the family models for the new proteomes (the scan scripts, `run_list.tsv`), the core command (`calls.long.tsv.gz`), and composition from the sequence.

**Job count (stated here because it fixes the design).** Per-proteome scan rates from `sacct` (jobs 29681841 to 29681858): TMHMM 8 to 11 min, repeats about 2 min, allergen about 1.5 min, SignalP about 0.9 min on GPU, family models about 1.5 min. With SignalP and TMHMM reused, one proteome costs about 5 min (repeats, allergen, family models, converters). About 660 proteomes (the 603 background genera plus the targets plus the supplement, less the reused scan ones) cost about 55 hours, about 40 jobs at 1.25 to 1.5 hours each. If B2a fails and SignalP and TMHMM must be run, the job count exceeds the limit of 40 to 60 and the plan returns to the owner. The plan splits the full list of proteomes into one-job groups and does not exceed 60 submitted jobs at any time.

**The 56 scan proteomes** were run on 2026-10-09 with 10 active families. Today 18 are active (commits `cbfdc45` and `d94a8ab`). Their module tables and calls are **not reused as they are**. B2 re-runs the converters and the core command on them with the current family table (the raw Pfam tables hold hits for all 24 models, so no new `hmmsearch` is needed). B2 checks that `run.json` and the `pfam_*.json` of every proteome list the same active set, and stops if they do not.

All new jobs use `${SCRATCH:?}` for node-local I/O, absolute paths and never `BASH_SOURCE[0]`. The Pfam path is the resolved `2026-01-27-Pfam38.2/Pfam-A.hmm`, not `current`. Outputs are written as `.tsv.gz`.

## 5. Analysis

### 5.1 Pools

Two pools per proteome, because the "unexplained" pool removes the classes that matter most for the three purposes (proteins with a tandem repeat call, for example the Ser/Thr-rich repeat glycoproteins):
- **Pool U (unexplained):** `signal_peptide_protein[R0]` called and `other_surface_no_mechanism[R0]` called, and no mature TM helix and no C-terminal HDEL or KDEL (those are reported separately as ER/Golgi/membrane proteins).
- **Pool S (secretome):** `signal_peptide_protein[R0]` called, same TM and HDEL/KDEL exclusion. It contains the proteins with mechanism evidence. A family that is enriched only in pool S is then visible.

### 5.2 Families (Pfam-A) and Pfam-free clusters

For each family, each pool, and each target (Onygenales, Aspergillaceae):
- **Presence:** the genus's chosen proteome has at least one pool protein with a hit.
- **Copy number:** hits per 1,000 pool proteins in the genus's proteome.
- **Tests (all families with at least one present genus, so BH is not applied after a label-based filter):**
  - Fisher exact test on genus presence, **two-sided**, and the odds ratio with both genus fractions.
  - Mann-Whitney test on copies per 1,000 pool proteins across genera, so lineage-specific expansions are found.
  - BH over all families tested for each pool and target. The number of families tested is taken from the data and stated.
- **Reporting rule:** a family is reported only with at least 3 target genera present (for presence) or at least 3 target genera above the background median (for copy number). The rule limits what is reported. It has no statistical meaning.
- **Fisher p and q are a ranking statistic under an independence assumption that does not hold.** Every output table says this next to the column. In the review's computation, a family present in all 7 Aspergillaceae genera and in 10% of 150 background genera gives p = 4.2e-7, and one gain on the family's ancestral branch would give the same result.
- **Origins on the tree:** for each reported family, the number of independent origins on `final_tree.nw` by Fitch parsimony on presence/absence over all genera. A family with one origin is labelled "clade-specific, one origin". The count is printed beside every q. No other tree-based test is run.
- **Pfam-free units:** the same test on MMseqs2 clusters that have no Pfam hit (the Fungi_5k cluster tables), at the genus level. These clusters can enter the candidate lists, so proteins like SOWgp (no Pfam domain) are not excluded by design. B4 first checks how the cluster tables are organised and records whether they can be joined to the chosen proteomes.
- **Sister comparison** (descriptive, no test): the target against the 7 other Eurotiales genera, and against the rest of Pezizomycotina.

### 5.3 Composition, GPI and the supplement

- **Features per genus (pool U and pool S):** the fraction of pool proteins with Cys above 5%, Pro above 8%, Ser+Thr above 25%, N-glycosylation sequons per 100 aa above 1.0, and the GPI-predicted fraction (PredGPI results of Fungi_5k). **The four levels are fixed here before any data is seen**; they are chosen from the SOWgp, FLO11 and Pir examples in the reports and are not tuned. The test is Mann-Whitney across genera with the effect size and the genus counts, BH corrected. No candidates come from one feature alone.
- **Supplement within *Aspergillus*:** up to 40 species proteomes, chosen by spacing across the species tree (the subtree of *Aspergillus*; farthest-point sampling from `final_tree.nw`). `samples.csv` has no section column, so groups are the clades obtained by cutting that subtree at a depth that gives about 8 groups (the depth is written in the plan). For the families found for Aspergillaceae, the fraction of species with the family by clade. No test.
- **Pathogen contrast inside Onygenales (descriptive):** a table for each reported family, with genera that include mammal-pathogenic species against the other Onygenales genera. The list of genera and the source for each is written in the plan before the run. The n is small and no test is run.

### 5.4 Candidate lists and the role table

For each reported family (and each reported Pfam-free cluster): the protein, species, signal peptide, GPI call, composition, repeat call, CAZy and MEROPS class if any, and the antigen call (**C. immitis RS only**; the antigen module is absent in every other run).

The role table (host interaction, immune evasion, general surface):
- The default value for every family is "no stated role".
- Every role row cites its basis: the quoted Pfam or InterPro text, or a PMID. A row that rests on assistant reasoning only carries the label T4 in the row.
- Roles for the families are written from the sorted accession list **before** the statistics are attached.
- The first page of the report says: the role column is assistant judgement, not evidence.

## 6. Limits

1. Phylogeny confounds the result. Genus counting reduces repeated sampling of a genus. It does not remove deeper shared ancestry. The origin count (5.2) shows how many independent gains support a family; it is a parsimony count, not a test.
2. The Onygenales target includes saprotrophs and keratinophiles (*Onygena*, *Amauroascus*, *Malbranchea*, *Uncinocarpus*), a bee pathogen (*Ascosphaera*), dermatophytes and the dimorphic pathogens. The Aspergillaceae target includes genera that are not known as mammal pathogens (*Monascus*, *Xeromyces*, *Hamigera*, *Penicilliopsis*). A family enriched in either clade can be a lineage marker.
3. The pools depend on the active family set (18 of 24 today). Each proteome's `run.json` records it, and B2 checks that all proteomes agree.
4. Pfam-A covers known domains. Full Pfam-A for the pass is release 37.x; the 24 family models are 38.2.
5. SignalP 6 does not separate wall proteins from secreted enzymes. `other_surface_no_mechanism` is a residual and not a negative set.
6. PredGPI and SignalP 6 were not measured here against truth. GPI is a description.
7. Genus sampling: one proteome per genus; gene-model quality varies (BUSCO or protein count is recorded). Where a second proteome per genus exists, B4 repeats the presence test with it as a check.
8. About 3,000 families are expected; the actual number comes from the data.

## 7. Work plan

| Task | Content | Check |
|---|---|---|
| B1 | `select_genera.py`: targets, the full background, the supplement; `selection_discovery.tsv`; the named-genus counts; the failing genera and the *Paracoccidioides* mapping | counts per group in section 3 are filled in with the real numbers |
| B2a | Compare Fungi_5k SignalP 6 and TMHMM results with the scan run on the 56 proteomes | agreement table; the decision to reuse or run |
| B2 | Re-run converters and the core command on the 56 scan proteomes with the current family table; run the new proteomes in one-job groups (job plan first) | all jobs COMPLETED; every `run.json` lists the same active set |
| B3 | Pools U and S; composition; join with Fungi_5k Pfam-A, PredGPI, CAZy, MEROPS; tests first with a small fixture | one proteome checked by hand |
| B4 | The tests of section 5 with tests on synthetic genus tables (a planted family is found; a family in 2 genera is not reported; the origin count of a one-gain family is 1); the 5-subsample stability table; the cluster-table check | tests pass |
| B5 | Pathogen-genus list with sources, the role table (blind to the statistics), the report and candidate lists | numbers copied from the output files |

Each task is its own commit. An independent review of this spec (rev 2) precedes B1. Stop points: an owner decision (O1 to O3), a failed review gate, B2a failing, or a job that fails twice.

## 8. Decisions for the owner

| # | Question | Default |
|---|---|---|
| O1 | Use all passing Pezizomycotina genera (about 600) as the background instead of a draw of 150? | yes (removes seed dependence; the cost is about 40 jobs) |
| O2 | Report the Fitch origin count beside every q-value? | yes |
| O3 | Reuse the Fungi_5k SignalP 6, TMHMM, Pfam-A 37.x and PredGPI results if B2a shows agreement? | yes |
