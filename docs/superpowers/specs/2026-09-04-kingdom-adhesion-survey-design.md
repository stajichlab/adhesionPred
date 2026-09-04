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

## Dependencies

- `statsmodels` — already present in the environment (0.14.4), used for the
  mixed-effects model.
- `scikit-posthocs` — **not currently installed**; needed for Dunn's
  post-hoc test. Add to `pyproject.toml` before implementation.
- `scipy` for Kruskal-Wallis (already a project dependency).

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

## Taxonomy-based pseudo-phylogenetic correction

Taxonomy is a nested proxy tree (phylum ⊃ class ⊃ order ⊃ family ⊃ genus),
just without branch lengths. Treating all species as independent samples in
the naive tests above inflates significance whenever a clade is
oversampled (e.g. one genus with 40 sequenced species dominating an order's
median). Two corrections are run alongside the naive tests, not instead of
them, and the report explicitly compares all three:

1. **Genus-averaged (clade-collapsed) tests.** Before running
   Kruskal-Wallis/Dunn's at phylum/class/order, first collapse to one row
   per genus (mean `adhesion_fraction`, mean `adhesion_count`, mean
   `mean_adhesion_prob` across that genus's species), then run the same
   tests on genus-level means. This removes the dominant pseudoreplication
   source (many congeneric species) with no new dependencies. Output:
   `tables/stats_by_phylum_genus_avg.csv` (and class/order equivalents),
   same shape as the naive tables.
2. **Nested mixed-effects model (stretch, still in scope per user request).**
   Fit `statsmodels.formula.api.mixedlm("adhesion_fraction ~ C(<rank>)", data,
   groups=data["genus"])` for each of phylum/class/order as the fixed effect,
   with **genus as a single random intercept**. This is a pragmatic
   simplification: `MixedLM` doesn't cleanly support 3-4 levels of nested
   random effects (order/family/genus all nested under a phylum fixed
   effect) without heavy variance-component (`vc_formula`) machinery: genus
   is chosen as the one random-effects level because it's the finest
   per-species grouping and captures the largest share of non-independence.
   Report the fixed-effect Wald test p-value and the estimated genus
   variance component (how much of the residual variance is "shared
   ancestry within genus" vs. residual noise). Requires adding
   `statsmodels` as a dependency if not already present. Output:
   `tables/mixedlm_by_phylum.csv`, `_by_class.csv`, `_by_order.csv` with
   fixed-effect estimate, p-value, and genus variance component per rank.
   New dependency check: confirm `statsmodels` is available in the project
   environment before implementation; if not, add it to `pyproject.toml`.
3. **Reporting the comparison.** For each rank, the report shows naive vs.
   genus-averaged vs. mixed-model results side by side. Where they largely
   agree, that's a more trustworthy signal. Where the naive test is
   significant but the genus-averaged/mixed-model result is not, that's
   flagged explicitly as "likely an oversampling artifact, not a genuine
   clade effect" rather than silently favoring one number.
4. **Remaining limitation, still stated in the caveats section:** this
   corrects for oversampling within genus but does not use real branch
   lengths or account for uneven divergence times between clades (a true
   phylogenetic comparative method would) — that remains the scope of the
   future tree-based follow-up.

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
   significance for closely related, over-sampled clades. This pass
   partially corrects for it via genus-averaging and a genus-random-effect
   mixed model (see "Taxonomy-based pseudo-phylogenetic correction" above),
   but neither uses real branch lengths or divergence times, so residual
   non-independence beyond the genus level remains uncorrected; (b)
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
    stats_by_phylum_genus_avg.csv (+ _posthoc.csv)
    stats_by_class_genus_avg.csv (+ _posthoc.csv)
    stats_by_order_genus_avg.csv (+ _posthoc.csv)
    mixedlm_by_phylum.csv
    mixedlm_by_class.csv
    mixedlm_by_order.csv
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
