# Adhesion prediction, round 2: benchmark, cruft removal and a fresh layout study (plan, 2026-11-26)

Owner decisions of 2026-11-26 are marked **D1–D4**. Nothing here has run. Every HPCC step and SLURM job needs the owner's go-ahead,
with the exact command shown first.

## Why

An inventory of every component (survey of 2026-11-26, summarised below) found that none has been measured outside the clade it was
built in with independent truth and enough positives. Every trained model so far reduces to amino-acid composition plus length:

| Model | Score | Composition + length | Source |
|---|---|---|---|
| `surface_glyco` ESM-2 + LR | ROC-AUC 0.999 leave-family-out; ~12% adhesin precision on S288C | 0.995 | `docs/model-review/2026-09-27-review-and-framework-plan.md` §4.1, §4.3 |
| Stage-2 ESM C 300M + LR | PR-AUC 0.939 leave-genome-out; **ROC-AUC 0.598 leave-clade-out** | 0.907 | same, §4.4–4.5 |
| AGA1 feature classifier (GB/RF/LR) | AUC 0.97 cluster-grouped CV | 0.973 (composition alone); no other feature group matters | `docs/reports/data/aga1_feature_classifier_2026-10-25/ablation.tsv` |

So the models learned "long, Ser/Thr-rich, repetitive glycoprotein", not adhesion. Two consequences shape this plan:

1. Any new signal must beat **composition + length**, or it is the same signal again.
2. A contrast against the whole proteome (or even the whole surface proteome) is largely won by composition. The decisive contrast is
   adhesins against **hard negatives**: wall and surface proteins with a known, non-adhesion job and the same composition.

What has held up: the signal-peptide gate (R0; sens 0.85, spec 0.97, held out to *A. nidulans*), PredGPI (yeasts), Pfam family calls
(precise on known families, 0/12 on seed-independent Kiss clusters), the hydrophobin domain call, and the calibration tooling.

## Parts

### 2a. Run the cross-clade structure search as built

`analysis/crossclade_domain_search/` README steps 5–8 (AFDB models, Reseek + Foldseek, 2-of-3 rule, gene-model check). Unchanged; it
also seeds the domain library in 2b.

### 2b. Benchmark, domain library and ablation harness

**Truth (D3: both contrasts, tiered; D4: curate, then freeze).**

1. *Curation pass* (dated, one PR): add experimentally shown adhesins and hard negatives from under-covered clades (Basidiomycota,
   Mucoromycota, Pezizomycotina beyond *Aspergillus*, chytrids if any). Each entry keyed by UniProt or NCBI accession (no JGI IDs),
   with an evidence code: `ko_assay` (deletion changes adhesion in an assay), `het_expr` (heterologous expression confers adhesion),
   `binding` (purified protein or domain binds a ligand or host cell), plus PMIDs. Moonlighting cytoplasmic proteins (GAPDH, HSP60,
   enolase) are kept but labelled `mechanism=moonlighting` and left out of the layout study (they have no secretory layout).
   Starting files: `data/curated/adhesins/adhesins.tsv` (E1), `data/controls/repeat-adhesins-other-species/`,
   `data/controls/wall-fusion-kiss2019/`, `data/controls/hard-negatives-shared/`.
2. *Hard negatives*: Pir, Cwp, Gas, Crh, Ecm33/Sps2, Dfg5/Dcw1, GPI glucanases, Yps aspartyl proteases, Ccw/Tip/Srp families and their
   orthologs with direct evidence of a non-adhesion role. A protein with mixed evidence is excluded, not guessed.
3. *Owner sign-off*, then the set is **frozen and versioned** (`data/benchmark/adhesion_v1/`, `SHA256SUMS`, `VERSION`); later results
   cannot change it. A v2 is a new folder.

**Contrasts.** Every component is reported against both: (C1) adhesins vs hard negatives, which decides keep or remove; (C2) adhesins
vs all other secreted or GPI proteins in the same proteomes, descriptive.

**Splits.** Leave-one-clade-out at class level (Saccharomycetes, Schizosaccharomycetes, Eurotiomycetes, Sordariomycetes,
Agaricomycotina/Pucciniomycotina, Mucoromycota), plus cluster-grouped CV (MMseqs2 30% identity clusters never split). A clade is
scored only if it holds at least 5 positive clusters; otherwise it is listed as "not scored" with its counts.

**Baseline every component must beat.** Logistic regression on 20 amino-acid frequencies + log length, refit inside each split.

**Domain library.** Cropped, folded adhesion head domains with structure clusters: PA14 (Epa/Flo5), Flo11 domain, Als N-terminal
Ig-like, CFEM, CPL1 (PF21671), LysM, fasciclin, Hyr/Iff N-terminal, Bys1/thaumatin, hydrophobin, BAD1/WI-1 repeat, plus round 2a
candidates that pass. Each domain carries its own empirical null (background split as in round 2a).

**Ablation rule.** A held-out component returns only if, on C1, it improves on the baseline (or on baseline + kept components) in the
pooled leave-clade-out score with a 95% bootstrap interval (by cluster) that excludes zero. Metrics: PR-AUC and recall at a fixed 1%
false-positive rate on C1, with intervals.

### 2c. Fresh layout study (D-scale: one proteome per species, broad)

**Sequences.** One reference proteome per species across all fungal classes, species list and accessions from PopSeq (per the
metadata rule; missing fields are fixed there by issue, not queried here). Secreted proteins only (R0 gate), plus GPI status.

**Per-protein description.**

- *Layout tokens*: signal peptide → folded N-terminal domain (structure cluster id from the domain library or a de novo Foldseek
  cluster) → stalk → anchor (GPI, Pir repeat, Kre9/Cwp-type, transmembrane, none). Distances and lengths of each segment.
- *Stalk profiles* (windowed, position-normalised): hydrophobicity (Kyte–Doolittle and a hydrophobic-moment variant), net charge,
  β-aggregation / amyloid-core propensity (as in the Als T-region), N- and O-glycosylation site density, Ser/Thr fraction, proline
  spacing, repeat period, copy number and copy-to-copy divergence.
- *Label-free evolutionary signals*: gene-family size and expansion per lineage, distance to chromosome ends where assemblies allow,
  repeat copy-number variation between strains of a species.

**Tests.** For each feature: effect size and 95% interval for C1 and C2, and whether it adds to the composition + length baseline in
leave-clade-out (same rule as 2b). Unsupervised: cluster layouts across the broad set and report which clusters are enriched for truth
positives and whether they recur across clades. Results are hypotheses unless they pass the 2b rule; the report says which is which.

**Multiple testing.** Feature tests are pre-listed in the 2c README before results are read; FDR (Benjamini–Hochberg) across them.

## Cruft (D2: remove now, tag first)

Tag the current `main` commit `legacy-2026-11-26` (owner pushes the tag), then one PR removes from `main`:

- `src/surface_glyco/` (ESM-2 classifier; no shipped model; equals the composition baseline) and its entry points and tests;
- analyses built on it: `analysis/adhesion_properties`, `analysis/kingdom_survey`, `analysis/embedding_clustering`,
  `analysis/model_review/01–03`;
- the `cys8` module subcommand (finds 1 extra hydrophobin; no call reads it); its code and spacing move to
  `analysis/hydrophobin_truth/` so the hydrophobin validation still reproduces;
- ~~PR 6 ecm33 custom HMMs~~: kept as report data after all. No code in `src/` reads them, and three test files and the PR 6
  evaluation record pin them by checksum, so removing them would break the record of why they don't help.

**Held out** (kept in the repo, off in every call and aggregate until they pass the 2b rule): `hydrophobin_relaxed`, lectin models,
`cys_rich` and `expression` lookups, the AGA1 feature classifier and `model_surface_candidate`, `disorder`, `repeats` (repeat02,
repeat14) and the calls that read them (`gpi_long_disorder_protein`, `tandem_repeat_protein`, `cell_wall_adhesion_candidate`,
`surface_attachment_candidate`). Status files of held-out calls stay but are marked `held_out`.

**Kept**: signal-peptide gate (R0), `gpi`, Pfam family calls (as domain labels, not adhesion claims), hydrophobin domain call,
calibration and status tooling, the allergen and antigen modules (not adhesion components; out of scope here).

Two catches the removal PR must handle: the distribution in `pyproject.toml` is itself named `surface_glyco` (with the three
`surface_glyco_*` scripts and `package-data`), so the project needs renaming (proposal: `adhesionpred`, keeping
`cellsurface_sorting_hat` as the package); and the Phase B jobs of `analysis/step1_compare/` (the folder that produced R0, kept)
import the ESM-2 embedding helper `surface_glyco.embeddings`, so that helper moves to `analysis/step1_compare/jobs/esm_embed.py`
rather than disappearing (owner choice 2026-11-26).

The removal PR lists every dependent file it touches; the sorting-hat acceptance harness must still give 0 differences for kept calls.

## Order

1. Round 2a (needs the HPCC back; go-ahead per step).
2. Cruft PR (tag first) — independent of the HPCC.
3. Curation pass → owner sign-off → freeze `adhesion_v1`.
4. Benchmark harness and baseline (code + fixture tests), then ablation of kept and held-out components.
5. Layout study: feature extraction on the broad set (HPCC), pre-listed tests, report.

## Open questions for later

- Which assemblies are good enough for the telomere-distance feature.
- Whether de novo structure clusters in 2c need AFDB coverage for proteins over 2,700 residues (ESMFold or crop-and-predict).
