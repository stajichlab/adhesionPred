# Adhesion Protein Embedding-Space Clustering — Design

Date: 2026-09-04
Author: Jason Stajich (with Claude)

## Purpose

Investigate whether the ~750,000 adhesion-predicted fungal proteins
(`analysis/adhesion_properties/`) fall into distinct structural/functional
sub-types (e.g. FLO11-like mucin/GPI-anchored adhesins vs. Als-like vs.
other architectures) by clustering them in protein-language-model
embedding space, then validating any clusters found against the
length/composition/domain data already computed.

Two embedding spaces are used, answering two different questions:

1. **The classifier's own embedding space** (`esm2_t12_35M_UR50D`, layer 6,
   mean-pooled — the exact representation `adhesion_predict` was trained
   and predicts on) — reveals what the *current classifier* is
   discriminating on internally. Directly relevant to the prior analysis's
   caveat that the CAZy AA1/laccase enrichment might be a classifier
   artifact: if AA1-flagged proteins form a tight, separate cluster here,
   that's evidence the classifier is keying on a composition quirk rather
   than a real adhesin feature.
2. **ESM Cambrian (ESM C) 300M** (`esmc-300m-2024-12`, EvolutionaryScale,
   released Dec 2024, successor to ESM2, permissively licensed) — a
   substantially stronger, independent representation not shaped by our
   small 35M-parameter classifier, better positioned to reveal genuine
   biological sub-architectures regardless of what the classifier itself
   can resolve.

## Non-goals

- No retraining or modification of `adhesion_predict`'s classifier.
- No new sequence generation (ESM3-style) — this is representation/
  clustering only.
- No changes to `analysis/kingdom_survey/` or `analysis/adhesion_properties/`
  beyond reading their already-committed/on-disk tables.
- Not a rigorous hypothesis-testing pass like the prior two analyses —
  this is exploratory pattern-finding; cluster "significance" is judged by
  domain/property coherence, not a formal statistical test.

## Prior-art check (verified, not assumed)

Grepped the full repo (files + git history, including diff content via
`git log -S`) for any prior use of ESM C, ESM Cambrian, or ESM3 — no
matches. `Changes.md`'s 1.0.0 release notes confirm only `esm2_t6_8M_UR50D`
and `esm2_t12_35M_UR50D` have ever been used in this project. ESM C is
entirely new to this codebase.

## Key technical facts established during design (not assumed)

- The classifier's actual embedding recipe (`src/adhesion_predict/embeddings.py`):
  `esm2_t12_35M_UR50D`, `repr_layers=[6]` (the *middle* layer of a 12-layer
  model, not the final layer), mean-pooled over the token dimension,
  truncated to 1022 residues. This exact recipe must be replicated for the
  "classifier's own embedding space" — using a different layer or the final
  layer would not answer "what does the classifier see."
- Embedding dimensionality: `esm2_t12_35M_UR50D` layer 6 → 480-dim.
  `esmc-300m-2024-12` → 960-dim (30 layers), mean-pooled over residue
  tokens excluding BOS/EOS.
- 749,697 proteins × 480 or 960 floats (float32) is ~1.4GB / ~2.9GB per
  embedding matrix — small enough to fit in memory for PCA/clustering on a
  single node; no distributed-computing infrastructure needed.
- `fair-esm` (already a dependency) loads ESM2 variants; ESM C requires a
  **new** dependency — either the `esm` PyPI package's
  `.from_pretrained_esm('esmc_300m')`, or HuggingFace `transformers`'
  `AutoModelForMaskedLM.from_pretrained("Biohub/esmc-300m-2024-12")`. Not
  yet verified to actually load in this environment (network access,
  potential gating) — this is the single biggest unknown in this design
  and is checked first, before building the rest of the pipeline around it.
- GPU resources are available on this cluster: `gpu`/`preempt_gpu`
  partitions include A100/H100/ada6000 nodes (`sinfo` confirmed). Embedding
  extraction for ~750K short protein sequences through a ≤300M-parameter
  transformer is expected to take well under an hour on a modern GPU, but
  this is an estimate, not measured — the first real-data extraction run
  is itself a check on this.

## Protein set

All 749,697 adhesion-predicted protein IDs from
`analysis/adhesion_properties/tables/protein_universe.csv` (`group="adhesion"`
rows) — same set already analyzed for composition/domains. Peptide
sequences are re-fetched from `functionalDB/function.duckdb`'s
`gene_proteins` table via the already-built, already-reviewed
`analysis/adhesion_properties/db.py`'s `fetch_lengths` (reused, not
duplicated — it already handles the `.protein` suffix quirk and the
DuckDB memory-limit fix from the prior project's incident).

## Pipeline

1. **Feasibility check** (first task, gates everything else): load both
   models on a handful of real sequences, confirm both produce embeddings
   of the expected dimensionality with no errors. If ESM C cannot be
   loaded in this environment (dependency/network/gating issue), stop and
   report back before building anything further — do not silently fall
   back to a different model.
2. **Extract embeddings** for both models over the full 749,697-protein
   set, as two separate GPU SLURM batch jobs (same operational pattern
   established for the prior project's memory-constrained CPU job:
   explicit resource request, `$SCRATCH` for any spill/temp data, chunked/
   incremental writes so a mid-run failure doesn't lose all progress).
   Output: `embeddings_esm2_classifier.npy` (749697×480) +
   `embeddings_esmc300m.npy` (749697×960), each with a parallel
   `protein_ids_*.txt` (row order).
3. **Reduce + cluster**, independently per embedding space: PCA to ~50
   components (captures most variance, denoises, makes clustering
   tractable), then HDBSCAN (density-based — does not require
   pre-specifying a cluster count, appropriate since the number of
   architectures is unknown) for cluster labels; UMAP 2D projection for
   visualization (on the full set if tractable, otherwise a fixed random
   subsample sized for plot legibility — decided at implementation based
   on actual UMAP runtime at this N).
4. **Validate against existing data**: for each embedding space's
   clustering, join cluster labels against
   `analysis/adhesion_properties/tables/protein_sequence_properties.csv`
   and `protein_domain_flags.csv` (both already on disk, regenerable via
   that project's own scripts if absent) to compute, per cluster: N,
   median length/Ser-Thr-Pro%/hydrophobicity, fraction with each domain
   flag, and top enriched Pfam/CAZy/MEROPS domains (reusing
   `analysis/adhesion_properties/domains.py`'s already-fixed,
   deterministic `top_domain_table` logic, applied per-cluster instead of
   adhesion-vs-background).
5. **Compare the two clusterings**: a contingency table / cross-tabulation
   between ESM2-classifier-space cluster labels and ESM-C-space cluster
   labels for the same proteins (same protein, two different embedding
   spaces, do they agree on groupings?), plus a targeted check: do the
   CAZy AA1-family-flagged proteins concentrate in one cluster (either
   embedding space), or are they scattered — directly informing the prior
   analysis's "classifier artifact vs. real biology" caveat.
6. **Figures**: UMAP scatter colored by cluster (both embedding spaces);
   UMAP scatter colored by key property/domain features (Ser-Thr-Pro%,
   has_signal_peptide, has_cazy) to visually check whether clusters align
   with known features; a per-cluster domain-enrichment summary
   (ranked table, not a heatmap, per the earlier project's established
   preference against min-max-normalized heatmaps over few rows).
7. **Report**: `analysis/embedding_clustering/REPORT.md` — cluster counts
   found per embedding space, per-cluster property/domain signatures,
   whether recognizable architectures emerge (and whether they resemble
   FLO11-like vs. Als-like vs. other patterns, stated as a hypothesis the
   data suggests, not a confirmed identification), the two-embedding-space
   agreement analysis, the AA1/laccase concentration check, and caveats
   (exploratory/unsupervised — no significance testing; HDBSCAN's
   density-based "noise" label for un-clustered points is expected and
   should be reported as its own bucket, not hidden; embedding models
   themselves have their own biases from pretraining data).

## File layout

```
analysis/embedding_clustering/
  00_feasibility_check.py        # gates everything else
  embed.py                       # shared embedding-extraction helpers
  01_extract_embeddings_esm2_classifier.py   # GPU SLURM job
  02_extract_embeddings_esmc300m.py          # GPU SLURM job
  cluster.py                     # shared PCA/HDBSCAN/UMAP functions
  03_cluster_esm2_classifier.py
  04_cluster_esmc300m.py
  validate.py                    # per-cluster property/domain summaries
  05_validate_and_compare.py      # cross-tab + AA1 concentration check
  06_make_figures.py
  07_generate_report.py
  tables/
  figures/
  REPORT.md
```

## Dependencies

- New: a package to load ESM C (`esm` PyPI package or `transformers` +
  the HuggingFace mirror weights — decided during the Task 1 feasibility
  check based on what actually works in this environment), `hdbscan`,
  `umap-learn`, `scikit-learn` (PCA — check if already present).
- Reused, not duplicated: `analysis/adhesion_properties/db.py` (peptide
  fetch), `domains.py` (top-domain logic), `src/adhesion_predict/embeddings.py`
  (the classifier's exact embedding recipe, called with its existing
  `get_esm_embeddings` function rather than reimplemented).

## Verification plan

- Task 1's feasibility check must pass (both models load, produce
  correctly-shaped embeddings on ≥5 real sequences) before any
  full-scale extraction is attempted.
- Row-count reconciliation: both embedding matrices' row count must equal
  749,697 (or the exact count of successfully-embedded proteins, with any
  shortfall — e.g. from a sequence exceeding a model's max length —
  explicitly logged, not silently dropped).
- Sanity bounds: no NaN/Inf in extracted embeddings; cluster labels
  include the expected HDBSCAN noise label (-1) as its own reported
  category, not discarded.
