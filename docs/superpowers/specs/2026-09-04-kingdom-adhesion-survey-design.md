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

- `results/*.adhesion_predict.csv` (**5,803 files** — confirmed by
  `ls results/*.adhesion_predict.csv | wc -l`; the earlier count of 5,805
  included `names.txt` and the `models/` symlink present in that directory,
  not actual result files) — columns `id,prediction,probability_adhesion`,
  already filtered to "Adhesion" calls only (default `predict.py` behavior).
- `/bigdata/stajichlab/shared/projects/Fungi_5k/samples.csv` (5,813 rows) — columns
  `ASMID,SPECIESIN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,
  CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,LOCUSTAG`.
- `/bigdata/stajichlab/shared/projects/Fungi_5k/input/*.proteins.fa.fai` (5,813 files) —
  used both for line counts (= total protein count per proteome) and, per the
  revised join logic below, as the join anchor itself.
- Provenance: results were generated with `esm2_t12_35M_UR50D` (per
  `run_fungi5k*.sh`); the minimum observed `probability_adhesion` is 0.5027,
  consistent with the model's internal 0.5 decision threshold. Since
  `predict.py` writes only "Adhesion" calls by default, **no threshold
  sweep is possible from these files** — `adhesion_fraction` is fixed to
  whatever the model's built-in 0.5 cutoff calls positive. This is stated
  as a permanent limitation in the report, not something this analysis can
  probe further without rerunning prediction with `--show-all`.

## Join logic (the crux of this design) — revised

Result filenames are derived from `species_strain` (e.g.
`Absidia_glauca_CBS_101.48_substr._RVII-324_met-.adhesion_predict.csv`), which
does **not** match `samples.csv` directly (`samples.csv` keys on `ASMID` and
`LOCUSTAG`, not this display name).

The original design proposed extracting the LOCUSTAG from the *result* CSV's
first data row. The Opus-reviewed, verified approach is **more reliable and
also doubles as a QC check**, so it replaces that:

1. **Primary join key: filename stem → `.fai` → LOCUSTAG.** Result filename
   stems match `input/*.proteins.fa` basenames exactly, 5,803/5,803 (verified
   by set-diff). Take the *first* protein ID in that species'
   `<stem>.proteins.fa.fai`, split on `_`, take token 0 as the LOCUSTAG. This
   was verified to succeed for all 5,813 `.fai` files with zero exceptions —
   i.e. it is the input side of the join, independent of whatever is (or
   isn't) in the result CSV.
2. Look up that LOCUSTAG in the `samples.csv`-derived dict to get taxonomy +
   `ASMID`.
3. Get total protein count for that species: line count of the same `.fai`
   file used in step 1.
4. **QC cross-check (catches mismatched annotation sets):** also read the
   result CSV's own first data row and extract its LOCUSTAG prefix the same
   way. If it disagrees with the `.fai`-derived LOCUSTAG, the species is
   **not** silently joined — it's written to `tables/mismatched_locustag.csv`
   with both prefixes, for manual diagnosis, and excluded from stats. This is
   a real, observed case, not a hypothetical: `Neohortaea_acidophila_CBS_113389`
   has result IDs that are NCBI `XP_…` accessions while its `.fai` starts
   `F009D8E7_000001-T1` — its results appear to come from a different
   protein annotation than the one used for the proteome-size denominator,
   which would silently corrupt its `adhesion_fraction` if joined naively.
5. **The "zero adhesion calls" fallback in the original design is removed —
   it was dead code.** Verified across all 5,803 result files: none are
   header-only (minimum call count is 1, median 113, max 1,148), so every
   result file has at least one data row to read for step 4's cross-check.
6. Any species where step 2 or 3 fails (no LOCUSTAG match, or no `.fai`) is
   written to `tables/unmatched_species.csv` with a reason column.
7. **Reverse reconciliation (missing from the original design):** also
   diff `samples.csv` LOCUSTAGs against the set of successfully-joined
   LOCUSTAGs to find species with **no result file at all**. There are
   confirmed to be exactly **10 such species, all in genus *Colletotrichum***
   (samples.csv has 101 *Colletotrichum* total — a genus/family drill-down
   there would silently under-count by 10 without this check). Written to
   `tables/missing_results.csv`.

Reconciliation check: `matched_count + unmatched_count (step 6) +
mismatched_count (step 4) = 5,803` (result files present), **and**
`matched_count + missing_count (step 7) = ` the count of samples.csv rows
whose LOCUSTAG appears in any `.fai` file. A manual spot-check of 3 known
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
  `probability_adhesion` column of adhesion-called proteins only. (No species
  has `adhesion_count == 0` in practice — verified, min is 1 — so this is not
  expected to produce NaN, but the column is still not treated as zero if it
  ever does.)

**Blank taxonomy handling.** `samples.csv` has real gaps: PHYLUM 1 blank,
CLASS 27, ORDER 66, GENUS 147, FAMILY 176, and **SUBCLASS 2,786 blank**
(more than half the rows — SUBCLASS is dropped from all aggregation/figures
entirely, kept in the master table only for reference). For every other
rank, a blank value means that species is **excluded from that rank's
grouping** (not pooled into a fake `""`/`"unknown"` group) — it still
appears normally in coarser ranks that aren't blank for it.

## Dependencies

- `statsmodels` — already present in the environment (0.14.4), used for the
  mixed-effects model.
- `scikit-posthocs` — **not currently installed**; needed for Dunn's
  post-hoc test. Add to `pyproject.toml` before implementation.
- `scipy` for Kruskal-Wallis (already a project dependency).

## Taxonomic aggregation & statistics

Primary ranks: **phylum → class → order**. Family/genus drill-down is
**pre-registered, not chosen post-hoc**: after the primary pass, drill down
into the top 3 orders by epsilon-squared effect size (not by p-value, given
the sample-size issue below), and label that section of the report
explicitly as "exploratory / hypothesis-generating" rather than confirmatory.

Phylum representation is **wildly unbalanced** (Ascomycota 4,057,
Basidiomycota 1,294, ... down to phyla with N=1, e.g. Sanchytriomycota,
Cryptomycota). At n≈5,800, Kruskal-Wallis will read as significant on
essentially any real-world clade difference, however small — **epsilon-
squared effect size is the headline statistic, p-values are secondary and
always reported alongside it, never alone.**

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
  size for the omnibus test (reported regardless of significance).
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
   species), groups ordered by median descending. **Every group is
   annotated with its N** (e.g. "Ascomycota (n=4,057)"), and any group with
   N < 5 is rendered in grey/desaturated rather than full color, so the
   reader can't mistake a 1-species phylum for a real distribution — sample
   size is the main confound here, so it has to be visible in the figure
   itself, not just in a table the reader may not open.
2. **Ordered box/violin of `adhesion_fraction` by order**, restricted to
   the top ~20 orders by species count (to stay legible), ordered by median,
   same N-annotation/greying convention as (1).
3. **Scatter: `total_proteins` vs `adhesion_count`**, colored by phylum,
   with a diagonal reference line for "expected count under the
   kingdom-wide median fraction" — the key plot for visually separating
   "just has a bigger genome" from real compositional enrichment.
4. **Ranked summary table (not a heatmap)**: phylum (or top orders) ×
   {median fraction, N, median count, median probability}, sorted by median
   fraction. A per-column min-max-normalized heatmap over only ~9 phyla
   was considered and rejected — with that few rows, min-max scaling
   visually inflates trivial differences into apparent hot/cold extremes,
   which would mislead more than a plain sorted table.
5. **Distribution of `mean_adhesion_prob` by phylum** — secondary signal:
   even at similar fractions, is one clade's adhesion calls more
   confidently positive? Same N-annotation/greying convention as (1).

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
   conclusion; (d) taxonomic sampling is wildly unbalanced (phylum N ranges
   from 4,057 down to 1) — at this sample size Kruskal-Wallis reads as
   significant on almost any real difference, so p-values are reported
   alongside epsilon-squared effect sizes, never alone, and single-species
   or near-empty groups are visually flagged rather than presented as if
   comparable to well-sampled ones; (e) `adhesion_fraction` is defined
   entirely by the model's fixed internal 0.5 probability threshold —
   results files contain only positive calls, so no threshold sensitivity
   analysis is possible from existing data; a different cutoff could shift
   which clades look "enriched."
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
    missing_results.csv
    mismatched_locustag.csv
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

- Reconciliation: matched + unmatched + mismatched species counts sum to
  5,803 (result files present); matched + missing-results counts reconcile
  against the samples.csv-side LOCUSTAG set (expect 10 missing, all
  *Colletotrichum*, confirmed during design review).
- Manual spot-check of 3 known species' taxonomy joins against `samples.csv`,
  plus confirmation that `Neohortaea_acidophila_CBS_113389` lands in
  `mismatched_locustag.csv` rather than being silently joined.
- Sanity bounds: `adhesion_fraction` ∈ [0, 1] for every row; no negative
  counts; `mean_adhesion_prob` ∈ [0, 1].
- Visual sanity check: figures render without errors, N-annotations and
  small-N greying are visibly correct, and axis ranges/labels are sane
  before finalizing the report.

## Design review

Reviewed by an independent Opus pass against the actual data before
implementation (2026-09-04). Findings incorporated above: corrected file
count (5,803, not 5,805), reverse reconciliation catching 10 missing
*Colletotrichum* species, a mismatched-annotation case
(`Neohortaea_acidophila_CBS_113389`) requiring a dedicated QC check rather
than silent exclusion, removal of dead-code fallback logic (verified no
result file is header-only), a stronger primary join key
(`.fai`-anchored instead of result-CSV-anchored), blank-taxonomy handling
including dropping SUBCLASS, replacing a misleading small-N heatmap with a
sorted table, N-annotation/greying on all group figures, pre-registering
the family/genus drill-down rule, and two additional caveats (sampling
imbalance, threshold fixed at 0.5 with no sweep possible).
