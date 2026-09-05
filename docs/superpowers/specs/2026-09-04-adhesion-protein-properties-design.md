# Adhesion Protein Molecular Properties & Functional Domain Survey — Design

Date: 2026-09-04
Author: Jason Stajich (with Claude)

## Purpose

The kingdom-wide survey (`analysis/kingdom_survey/`) established *how much*
of each fungal proteome is predicted adhesion-like and whether that varies
by clade. This design covers a follow-on analysis asking *what these
proteins are like*: their length, amino acid composition, and functional
domain annotations (Pfam, CAZy, MEROPS, signal peptide, transmembrane
topology, subcellular targeting) — compared against a background sample of
non-adhesion proteins from the same species, broken down by taxonomic
clade, to look for patterns that suggest function.

## Non-goals

- No ESM-2 embedding/clustering analysis in this pass — that question
  (whether "encoder/transformer processing" reveals sub-types within the
  adhesion-predicted set) is deferred to a follow-up phase, informed by
  what this analysis finds, per the agreed decomposition.
- No re-running of the model or changes to `src/adhesion_predict/`.
- No changes to `analysis/kingdom_survey/` beyond importing its already-
  generic `stats.py` functions (see "Reuse" below) — this is a new,
  independent analysis directory.

## Inputs

- `results/*.adhesion_predict.csv` (5,802 matched species, per
  `kingdom_survey`) — `id, prediction, probability_adhesion` per adhesion-
  called protein.
- `analysis/kingdom_survey/tables/species_adhesion_summary.csv` — per-
  species taxonomy (phylum/class/order/family/genus) and LOCUSTAG, already
  matched and reconciled. Used as the taxonomy lookup so this analysis
  doesn't re-derive the join from scratch.
- `/bigdata/stajichlab/shared/projects/Fungi_5k/input/*.proteins.fa.fai` —
  full per-species protein ID lists, used to sample the non-adhesion
  background.
- `/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb`
  (278GB, read-only) — tables used: `gene_proteins` (length, peptide
  sequence), `pfam`, `cazy_overview`, `merops`, `signalp`, `tmhmm`,
  `targetp`. All confirmed to cover the full 5,813-species set and to use
  the same `protein_id` convention as `results/` (e.g.
  `F07B100A_000481-T1`), except `gene_proteins.protein_id` which carries a
  `.protein` suffix (`F07B100A_000481-T1.protein`) — confirmed by direct
  query during design review, not assumed.

## Scope addition (flagged, approved)

`signalp`/`tmhmm`/`targetp` are pulled in addition to the three domain
databases the user named (Pfam/CAZy/MEROPS), because fungal adhesins
(FLO11/ALS-type) are classically GPI-anchored surface proteins with a
signal peptide and are *not* enzymes — secretion/topology features are
plausibly more informative about adhesin function specifically than
enzyme-family databases.

## Data-source facts established during design review (not assumptions)

- `pfam`: a protein can have many rows (avg 2.1/protein, max 270) — one row
  per domain hit, so "has a Pfam hit" requires `COUNT(DISTINCT protein_id)`
  and identity tallies must count **distinct proteins per pfam_id**, not
  raw hit rows (a single repeat-domain protein would otherwise dominate a
  frequency count).
- `tmhmm`, `targetp`, `signalp`: exactly one row per protein that has any
  hit at all — presence is just "does a row exist for this protein_id",
  no dedup needed.
- `targetp.prediction` has only two observed values: `SP` (secreted) and
  `mTP` (mitochondrial-targeted) — the table contains **no row** for a
  protein with no predicted targeting signal (absence of a row = "no
  target signal predicted", not a third category value).
- `cazy_overview.cazyme_fam` embeds the hit's alignment range in the string
  (e.g. `GH114(59-261)`) — must strip the `(...)` suffix to get the clean
  family code (`GH114`) before tallying.
- `merops.merops_id` is a per-sequence MEROPS accession (e.g.
  `MER0080922`), not a peptidase family/clan code — there is no
  family-level lookup table in this database, so MEROPS results are
  reported as raw accession-level hit presence/frequency, with this
  limitation stated in the report rather than presented as if it were
  family-level (avoids implying more specificity than the data supports).
- `gene_proteins.length` and `.peptide` are exact (verified against a known
  protein: `F07B100A_000481-T1.protein`, length 386, matches expectation).

## Protein universe: adhesion set + background comparator

1. **Adhesion set**: every `id` from every matched species' result CSV
   (5,802 species, ~750K proteins per the earlier kingdom-wide count),
   annotated with its species' LOCUSTAG/taxonomy (via
   `species_adhesion_summary.csv`) and its `probability_adhesion`.
2. **Background set**: for each species, a random sample of protein IDs
   from that species' `.fai` file, **excluding** IDs already in the
   adhesion set, sized to match that species' adhesion count 1:1 (or all
   remaining non-adhesion proteins if fewer are available). This keeps the
   background roughly the same size as the adhesion set (~750K) and
   balanced per-species, so clade-level comparisons aren't skewed by which
   species happen to be oversampled in one group vs. the other.
3. Both sets get the same downstream feature extraction (length, AA
   composition, domain flags) so every reported feature has a same-species,
   same-scale comparison.
4. Sampling is seeded (fixed random seed) for reproducibility, recorded in
   the report.

## Per-protein features

**Sequence properties** (from `gene_proteins.peptide`, computed in Python
after pulling the adhesion+background subset — not the full 56M-row
table):
- `length` (from `gene_proteins.length` directly, exact)
- `pct_ser`, `pct_thr`, `pct_pro` (individual), `pct_ser_thr_pro` (combined
  — the classic fungal-adhesin mucin-like/Ser-Thr-Pro-rich signature)
- `pct_cys`
- `aromaticity` (fraction Phe+Trp+Tyr)
- `mean_hydrophobicity` (Kyte-Doolittle scale, per-residue mean)
- `net_charge_ph7` (simple count of Lys+Arg minus Asp+Glu, not a full pKa
  model — stated as an approximation in the report)

**Domain/topology flags** (boolean presence, from the DuckDB tables):
- `has_pfam`, `has_cazy`, `has_merops`, `has_signal_peptide` (any `signalp`
  row), `has_tm_helix` (`tmhmm.PredHel >= 1`), `targetp_category` (`SP`,
  `mTP`, or `none` if no row)

**Top-domain identity tables** (separate from the per-protein flag table):
for Pfam, CAZy, and MEROPS independently: count of **distinct adhesion
proteins** hit by each specific domain ID/family/accession, alongside the
same count for the background set, so the report can show enrichment
(adhesion rate vs. background rate) for the top ~20 most common hits, not
just raw adhesion-side frequency.

## Aggregation & statistics

Reuses `analysis/kingdom_survey/stats.py`'s `aggregate_summary` and
`kruskal_wallis_by_rank` directly (imported via the same `sys.path.insert`
pattern `kingdom_survey`'s own CLI scripts already use) — both are generic
over `(df, rank_col, value_col)` with no kingdom-survey-specific
assumptions, so duplicating them here would be pure repetition. `dunn_posthoc`
is reused the same way where a per-rank omnibus test is significant.

Per continuous feature (length, each composition metric), per rank
(phylum/class/order, same `MIN_GROUP_N = 5` convention): naive
Kruskal-Wallis + epsilon-squared comparing adhesion vs. background groups
*within* each clade is not the right test shape here (the grouping factor
of interest is clade, with adhesion/background as a second factor) — instead:
- Compute `aggregate_summary`-style per-clade median/IQR **separately** for
  the adhesion set and the background set (two calls, same function).
- Run `kruskal_wallis_by_rank` **within the adhesion set** across clades
  (does the feature vary by clade among adhesion proteins — same shape as
  the kingdom survey's existing tests) and **separately within the
  background set** (does the same feature vary by clade among ordinary
  proteins) — comparing the two clade-level patterns side by side shows
  whether an apparent clade effect on, say, length is specific to adhesion
  proteins or just reflects that clade's overall protein-length tendency.
- A simple two-group (adhesion vs. background) Mann-Whitney U test **within
  each clade** answers "is this feature actually different for adhesion
  proteins vs. the rest of the proteome, in this clade" — run per rank per
  clade with the same `MIN_GROUP_N` gating on both groups' sample sizes,
  BH-corrected across clades within a rank.

For domain-presence flags (binary), per clade: fraction-with-hit for
adhesion vs. background, no formal test in this pass (descriptive
comparison only, consistent with how the kingdom survey already treats
some binary/categorical breakdowns) — the tables and grouped bar charts
carry the comparison, since with many small-N clades a formal proportions
test would be as unreliable as the kingdom survey's own p-value caveat
already explains for continuous features.

## Figures

Dataviz skill invoked at implementation time (as it was for
`kingdom_survey`), same palette/N-annotation conventions:
1. Boxplot of `length` by phylum, adhesion vs. background as paired boxes
   per phylum group (N-annotated, small-N greyed, same convention as
   `kingdom_survey`'s `boxplot_by_rank`).
2. Boxplot of `pct_ser_thr_pro` by phylum, same paired adhesion/background
   layout — the figure most directly testing the "adhesins are
   Ser/Thr/Pro-rich" hypothesis across clades.
3. Grouped bar chart: fraction with a Pfam hit / CAZy hit / MEROPS hit /
   signal peptide / TM helix, adhesion vs. background, by phylum (one
   panel or facet per domain type, or one bar chart per domain type —
   decided at implementation based on legibility with 9 phyla).
4. Top-domain enrichment bar chart: for Pfam (and separately CAZy), the
   top ~15-20 most common domains among adhesion proteins, with adhesion
   rate and background rate as paired bars, sorted by enrichment ratio
   (adhesion rate / background rate) rather than raw frequency — a domain
   that's common everywhere isn't as informative as one that's
   specifically enriched among adhesion proteins.
5. Scatter: `length` vs. `probability_adhesion` (does the model's
   confidence correlate with protein length?).
6. Scatter: `pct_ser_thr_pro` vs. `probability_adhesion` (does the
   Ser/Thr/Pro signature correlate with model confidence — a sanity check
   on whether the classifier is actually picking up this known biological
   signature, or something else entirely).
7. Scatter: per-species mean adhesion-protein `length` vs. that species'
   `adhesion_fraction` (from `kingdom_survey`'s summary table) — does
   "typical adhesin size" relate to "how much of the proteome is
   adhesion-like" at the species level?

## Report

`analysis/adhesion_properties/REPORT.md`, same structure convention as
`kingdom_survey/REPORT.md`: Scope (adhesion set size, background sample
size and method, per-species matching), Headline findings (hand-written,
sourced from real numbers the same way `kingdom_survey`'s Task 8 was),
per-feature breakdowns with embedded figures and tables, a
top-domain-enrichment section, caveats (MEROPS accession-vs-family
limitation, net-charge approximation, background-sampling randomness,
non-independence across clades — same taxonomy caveat as the kingdom
survey, restated here since it applies equally), and a pointer to the
deferred embedding-clustering follow-up.

## File layout

```
analysis/adhesion_properties/
  db.py                          # DuckDB connection + query helpers
  properties.py                  # AA composition / sequence property functions
  01_build_protein_universe.py   # CLI: adhesion set + background sample -> tables/protein_universe.csv
  02_query_functional_domains.py # CLI: domain flags + top-domain tables
  03_compute_sequence_properties.py  # CLI: length + AA composition per protein
  04_aggregate_and_test.py       # CLI: merge + per-clade stats
  05_make_figures.py             # CLI: all figures
  06_generate_report.py          # CLI: REPORT.md assembly
  tables/
  figures/
  REPORT.md
```

## Verification plan

- Reconciliation: protein_universe.csv's adhesion-set row count matches
  the sum of adhesion_count across all matched species in
  `species_adhesion_summary.csv`.
- Spot-check 3 known protein IDs' length/peptide against a direct DuckDB
  query, matching the design-review verification method already
  established for `gene_proteins`.
- Sanity bounds: all percentage features in [0, 100]; `has_*` flags are
  boolean; background sample size per species ≤ that species' total
  protein count.
- Visual sanity check of all figures before finalizing the report, same
  as `kingdom_survey`.
