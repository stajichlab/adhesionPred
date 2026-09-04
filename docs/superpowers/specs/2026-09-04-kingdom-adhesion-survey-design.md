# Kingdom-wide Fungal Adhesion Protein Survey — Design

Date: 2026-09-04
Author: Jason Stajich (with Claude)

## Purpose

The `adhesion_predict` tool has already been run genome-wide over the Fungi_5k
collection, producing one CSV per species in `results/` listing every protein
called "Adhesion" with a probability score. This design covers a new,
one-shot analysis that aggregates those existing results against the
taxonomic metadata in `samples.csv` to ask: **do specific fungal clades
(phylum/class/order, drilling to family/genus where interesting) carry a
higher composition of adhesion-predicted proteins than others**, both in
raw counts and normalized by proteome size? Output is a citable Markdown
report with figures and tables, laying groundwork for a later
phylogeny-aware re-analysis.

## Non-goals

- No re-running of the model or changes to `src/adhesion_predict`.
- No formal phylogenetic comparative methods (PGLS, phylogenetic ANOVA) in
  this pass — flagged as a known limitation and left for a follow-up once a
  tree is imposed.
- No interactive dashboard — deliverable is a static, versionable report.

## Inputs

- `results/*.adhesion_predict.csv` (5,805 files) — columns `id,prediction,probability_adhesion`,
  already filtered to "Adhesion" calls only (default `predict.py` behavior).
- `/bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv` (5,813 rows) — columns
  `ASMID,SPECIESIN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,
  CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,LOCUSTAG`.
- `/bigdata/stajichlab/shared/projects/Fungi_5k/input/*.proteins.fa.fai` (5,813 files) —
  used only for line counts (= total protein count per proteome), not sequence content.

## Join logic (the crux of this design)

Result filenames are derived from `species_strain` (e.g.
`Absidia_glauca_CBS_101.48_substr._RVII-324_met-.adhesion_predict.csv`), which
does **not** match `samples.csv` directly (`samples.csv` keys on `ASMID` and
`LOCUSTAG`, not this display name). The reliable join key is the **LOCUSTAG
prefix embedded in every protein ID**, e.g. `F07B100A_000481-T1` → prefix
`F07B100A`, which equals `samples.csv`'s `LOCUSTAG` column exactly.

Procedure per result file:
1. Read the first data row's `id` field; split on `_`, take token 0 as the
   candidate LOCUSTAG.
2. Look up that LOCUSTAG in a `samples.csv`-derived dict to get taxonomy +
   `ASMID`.
3. If the result file has **zero** adhesion calls (0 data rows, i.e. no
   proteins predicted positive), fall back to matching by filename stem
   against the `.proteins.fa` basenames in the shared input dir to still
   recover a LOCUSTAG (via that file's own header prefixes) — necessary
   because such files carry no `id` to extract from. This is expected to be
   a small number of species.
4. Get total protein count for that species: `wc -l` equivalent (line count)
   of the matching `<same-basename>.proteins.fa.fai` in the shared input dir.
5. Any species where step 2 or 4 fails is written to
   `tables/unmatched_species.csv` with a reason column, and excluded from
   downstream stats/figures — not silently dropped.

Reconciliation check: `matched_count + unmatched_count` must equal 5,805
(the number of result files present), and a manual spot-check of 3 known
species' LOCUSTAG joins against `samples.csv` is done before trusting the
full table.

## Output: master species table

`analysis/kingdom_survey/tables/species_adhesion_summary.csv`, one row per
matched species:

```
locustag, asmid, species_name, phylum, subphylum, class, subclass, order,
family, genus, total_proteins, adhesion_count, adhesion_fraction,
mean_adhesion_prob, median_adhesion_prob
```

- `adhesion_count` = number of data rows in the result CSV (0 if the file
  has header only).
- `adhesion_fraction` = `adhesion_count / total_proteins`.
- `mean_adhesion_prob` / `median_adhesion_prob` computed over the
  `probability_adhesion` column of adhesion-called proteins only (undefined/NaN
  when `adhesion_count == 0`, not treated as zero).

## Taxonomic aggregation & statistics

Primary ranks: **phylum → class → order**. Family/genus reserved as an
ad hoc drill-down into whichever order(s) look most interesting after the
primary pass (not run uniformly across every family/genus up front — that
would produce hundreds of near-empty groups).

Per rank, per group:
- N species, median + IQR of `adhesion_fraction`, median + IQR of
  `adhesion_count`, median of `mean_adhesion_prob`.
- Groups with **N < 5** species are shown in tables (clearly marked
  "small-N, descriptive only") but excluded from the formal hypothesis test
  to avoid single-species noise driving conclusions.
- Kruskal-Wallis test across groups with N ≥ 5 at that rank, separately for
  `adhesion_fraction` and `adhesion_count`.
- Where Kruskal-Wallis is significant (p < 0.05): Dunn's post-hoc pairwise
  comparisons, Benjamini-Hochberg corrected, plus epsilon-squared effect
  size for the omnibus test.
- Output: `tables/stats_by_phylum.csv`, `tables/stats_by_class.csv`,
  `tables/stats_by_order.csv`, each with group summary + omnibus test result;
  pairwise post-hoc results in a companion `*_posthoc.csv` per rank.

## Figures

All figures follow the dataviz skill's guidance (loaded at implementation
time, before writing plotting code) for palette/consistency/theme-awareness
where applicable (these are static PNG/SVG for a Markdown report, so theme
handling applies loosely — light-mode print-friendly is the target).

1. **Ordered box/violin of `adhesion_fraction` by phylum** (all matched
   species), groups ordered by median descending.
2. **Ordered box/violin of `adhesion_fraction` by order**, restricted to
   the top ~20 orders by species count (to stay legible), ordered by median.
3. **Scatter: `total_proteins` vs `adhesion_count`**, colored by phylum,
   with a diagonal reference line for "expected count under the
   kingdom-wide median fraction" — the key plot for visually separating
   "just has a bigger genome" from real compositional enrichment.
4. **Ranked summary heatmap**: rows = phylum (or top orders), columns =
   {median fraction, median count, median probability}, values
   normalized/colored per column for at-a-glance comparison.
5. **Distribution of `mean_adhesion_prob` by phylum** — secondary signal:
   even at similar fractions, is one clade's adhesion calls more
   confidently positive?

Figures saved to `analysis/kingdom_survey/figures/*.png` (and `.svg` for
anything destined for print).

## Report

`analysis/kingdom_survey/REPORT.md`: narrative structure —
1. Summary of scope (N species, N phyla/classes/orders covered, matched vs
   unmatched).
2. Headline finding(s) with the two or three most relevant figures inlined.
3. Full per-rank breakdown (remaining figures + links to full stats tables).
4. **Explicit caveats section**: (a) species within a genus/family are not
   statistically independent — shared ancestry inflates apparent
   significance for closely related, over-sampled clades; Kruskal-Wallis
   here treats them as i.i.d., which this pass does not correct for; (b)
   proteome annotation completeness/quality varies by genome and could
   confound raw counts (mitigated somewhat by using fraction, not just raw
   count); (c) classifier was trained on a small positive set (FLO/ALS1-like
   proteins from *Saccharomyces*/*Candida*) so may under-call adhesion-like
   proteins in very divergent lineages (e.g. early-diverging fungi) whose
   domains don't resemble the training set — flag phyla with unusually low
   fractions as a training-bias hypothesis to test, not just a biological
   conclusion.
5. Pointer to the follow-up: imposing phylogeny (from `nf_phyling` output or
   known taxonomy-as-tree) for a phylogenetically-corrected re-analysis.

## File layout (final)

```
analysis/kingdom_survey/
  01_build_species_table.py
  02_taxonomic_stats.py
  03_make_figures.py
  tables/
    species_adhesion_summary.csv
    unmatched_species.csv
    stats_by_phylum.csv (+ _posthoc.csv)
    stats_by_class.csv (+ _posthoc.csv)
    stats_by_order.csv (+ _posthoc.csv)
  figures/
    *.png, *.svg
  REPORT.md
```

## Verification plan

- Reconciliation: matched + unmatched species counts sum to 5,805.
- Manual spot-check of 3 known species' taxonomy joins against `samples.csv`.
- Sanity bounds: `adhesion_fraction` ∈ [0, 1] for every row; no negative
  counts; `mean_adhesion_prob` ∈ [0, 1] or NaN.
- Visual sanity check: figures render without errors and axis ranges/labels
  are sane before finalizing the report.
