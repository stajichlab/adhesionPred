# B: Discovery pass on the unexplained secretome (Onygenales and Aspergillaceae)

*2026-10-09. Revision 1. Author: Claude Code (claude-sonnet-5-5). Part of the tier restructuring agreed with the owner on 2026-10-09 (order B, then C, then A). Status: draft for independent review. No code is written. Numbers come from files named in the text.*

## 1. The question

The scan of 56 Fungi_5k proteomes (`docs/reports/2026-10-09-fungi-scan-descriptive.md`) shows that 97% of proteins with a signal peptide have no mechanism evidence (medians per proteome: Onygenales 413 of 427, *Aspergillus* 784 of 809). The owner wants to know whether the Onygenales and the Aspergillaceae carry patterns in this pool that suggest additional candidates for (1) interaction with the host, (2) protection from immune recognition, or (3) general cell surface functions.

This is a descriptive discovery pass. It makes no call and it does not change `categories.yaml`. Its output is a table of families and features that are more common in the target clades than in a background, with a candidate list. It states associations, not selection and not function.

## 2. Owner decisions that fix the design (2026-10-09)

1. The unit of analysis is the **genus**. A pattern is reported only if it appears in **at least 3 independent genera** of the target clade.
2. The background is larger than the 56 scan proteomes: one proteome per genus from Fungi_5k, about 150 to 200 genera. The owner accepts that phylogeny confounds the result.
3. Aspergillaceae (8 genera, 406 proteomes) is the main target for the *Aspergillus* question (the genus *Aspergillus* alone has 228 proteomes in 1 genus). A within-*Aspergillus* comparison is a descriptive supplement (section 5.3).
4. The tiers are broad (cell wall and surface), narrow (adhesion) and host interaction. This pass feeds the third tier and the broad tier.

## 3. Data

Counts from `/bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv` (5,813 proteomes):

| Group | Proteomes | Genera |
|---|---|---|
| Onygenales (target 1) | 71 | 28 |
| Aspergillaceae (target 2) | 406 | 8 |
| Eurotiales (includes Aspergillaceae) | 483 | 15 |
| Pezizomycotina | 2,802 | 654 |
| Other Ascomycota | 1,255 | 125 |
| Basidiomycota | 1,294 | 387 |

- **Selection:** one proteome per genus by the rule of `analysis/fungi_scan/select_proteomes.py` (skip fragmentary proteomes under 4,000 proteins; deterministic by species and strain). Target genera: all 28 Onygenales genera and all 8 Aspergillaceae genera, where each has a proteome that passes the 4,000-protein rule (genera without one are listed). Background: a stratified random sample of Pezizomycotina genera not in the targets (seed recorded), about 150, stratified by class so that no class dominates. Eurotiales genera outside Aspergillaceae (7) go to the background and are also reported as the sister comparison.
- **Supplement (5.3):** up to 40 *Aspergillus* species proteomes, chosen across the species in `samples.csv`.
- The 56 proteomes of the earlier scan are reused where they are the chosen proteome of a genus. Their outputs stay as they are.

## 4. Pipeline

Per selected proteome, with the existing scan scripts (`analysis/fungi_scan/run/`, `scan_step.sbatch`) extended with the new proteomes in `run_list.tsv`:

1. The seven module steps already run in the scan (SignalP 6, TMHMM, the Pfam family models, repeat02, repeat14, allergen BLAST), then the core command. This gives `calls.long.tsv.gz`.
2. **Unexplained secretome** = proteins with `signal_peptide_protein[R0]` called and `other_surface_no_mechanism[R0]` called (the definition in `docs/CLASSES.md`). The definition does not use the broad and narrow tiers (workstream A is later). Whatever changes in the config before this pass runs is recorded in `run.json` of each proteome.
3. **New steps, on the unexplained secretome only:**
   - A search of the full Pfam-A 38.2 (`hmmsearch --cut_ga`, the file the family scan used).
   - PredGPI (`predgpi/202001`) for GPI-anchor prediction. DeepLoc 2.1 is not used in this pass (its leakage to the training data is not known).
   - Composition from the sequence: Pro, Cys, Ser+Thr, Asp+Glu, N-glycosylation sequons (N-X-S/T, X not P), length, and the share of residues in repeat arrays from repeat14.
4. Output per proteome: a table of unexplained-secretome proteins with all Pfam hits, the GPI call and the composition.

Job sizing follows the HPCC rule (1 to 1.5 hours of real runtime per job). The scan timings (TMHMM 0.9 to 1.5 hours for 8 proteomes; SignalP 0.4 hours for about 28 proteomes) set the groups. The estimate of job count goes into the plan and is checked against the 40 to 60 job limit given on 2026-10-08.

## 5. Analysis

### 5.1 Families (Pfam-A) in the unexplained secretome

For each Pfam family, and each target clade (Onygenales, Aspergillaceae):
- Presence in a genus = at least one unexplained-secretome protein in the genus's selected proteome has a hit.
- Count genera with the family in the target and in the background.
- Test: Fisher exact test on genus counts (target present, target absent, background present, background absent); Benjamini-Hochberg over all families tested; report the odds ratio and both genus fractions.
- **Minimum:** the family is reported only with at least 3 target genera present. With 28 Onygenales genera and 8 Aspergillaceae genera, the minimum applies to both.
- Also report the sister comparison (Onygenales or Aspergillaceae against the other Eurotiales and Pezizomycotina separately) so a reader sees whether the pattern depends on the background.

### 5.2 Composition and GPI

For each feature (fraction Cys, Pro, Ser+Thr, sequons per 100 aa, GPI-predicted fraction): the per-genus value is the median over the unexplained-secretome proteins of the genus's proteome, or the fraction of proteins above a fixed level stated in the plan. Target against background by the Mann-Whitney test on genera, BH corrected, with the effect size. No candidates come from a single feature alone.

### 5.3 Supplement within *Aspergillus*

Descriptive only. For families found in 5.1 for Aspergillaceae, the fraction of the 40 *Aspergillus* species proteomes that have the family, by species group as given in `samples.csv`. No test.

### 5.4 Candidate lists

For each family reported in 5.1 (and for each protein that carries it): the protein, species, signal peptide, GPI call, composition, repeat call, and any hit to the antigen or allergen modules. The list of "plausible role" is a table in this spec's plan that maps Pfam descriptions to the three purposes (host interaction, immune evasion, general surface). It is written by the assistant from the Pfam descriptions and the literature, and it is a tier T4 review (assistant judgement, no paper per family unless stated). It is not a truth set.

## 6. Limits

1. The background is not independent of the targets: phylogeny. Genus-level counting reduces repeated sampling of a genus. It does not remove deeper shared ancestry, for example Aspergillaceae genera that share a family from a common ancestor. The result is an association.
2. The unexplained secretome depends on the family table: 18 of 24 families are active today. A family that is active in one run and not in another changes the pool. Record the active set in each run.
3. Pfam-A covers known domains. A family that is new or unknown does not appear. Proteins with no Pfam hit are reported as a count per proteome, and their composition is part of 5.2.
4. SignalP 6 predicts signal peptides. It does not separate wall proteins from secreted enzymes. The tool's `other_surface_no_mechanism` is a residual and is not a negative set.
5. PredGPI and SignalP 6 were not measured against truth here. GPI is a description.
6. Multiple testing: about 3,000 families are tested; with BH and a 3-genus minimum the number of reported families is small. The minimum of 3 genera has no statistical meaning.
7. Genus sampling: 28 Onygenales genera are all genera with a passing proteome; the background is a sample, so the result depends on the seed. The plan repeats the background draw with 5 seeds and reports the families found in all 5.

## 7. Work plan

| Task | Content | Check |
|---|---|---|
| B1 | `select_genera.py` (extends the scan selection): the target genera, the stratified background, the supplement; `selection_discovery.tsv` | counts per group match section 3; genera without a passing proteome listed |
| B2 | Run the module steps and the core command on new proteomes (job plan first) | all jobs COMPLETED; run.json per proteome |
| B3 | New steps: full Pfam-A, PredGPI, composition on the unexplained secretome; tests first with a small fixture | one proteome run by hand agrees with the script |
| B4 | The statistics of section 5 with tests on synthetic genus tables (a planted family is found; a family in 2 genera is not reported) | tests pass; 5-seed stability table |
| B5 | Report and candidate lists; the role table | numbers copied from the output files |

Each task is its own commit. An independent review of this spec precedes B1. Stop points: an owner decision, a failed review gate, or a job that fails twice.
