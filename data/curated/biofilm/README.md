# Biofilm label table (DRAFT)

`biofilm.tsv`, built by `analysis/curation/build_biofilm.py` from QuickGO experimental
(ECO:0000269) annotations in Fungi to *biofilm formation* (GO:0042710) and *regulation of
single-species biofilm formation* (GO:1900190), with descendants.

Kept separate from adhesion on purpose: most biofilm genes are transcriptional regulators
or matrix enzymes, not adhesins. Of 207 rows, about half are regulators.

| class | n | meaning |
|---|---|---|
| `biofilm_surface` | 39 | surface/secreted: candidate matrix, adhesin or cell-wall effector |
| `biofilm_regulator` | 102 | regulation term, or a transcription factor / kinase / chromatin regulator |
| `biofilm_other` | 66 | intracellular protein annotated to the biofilm process itself |

`direction` (`promotes` / `restrains` / `unknown`) comes from the positive- and
negative-regulation terms. It is `unknown` for nearly every row, because the bare process
term does not record direction — a known gap that expert review should fill.

`adhesin_table` cross-references `../adhesins/adhesins.tsv`: 17 of these proteins are
curated adhesins, so the two labels genuinely overlap and should not be merged.

**Coverage caveat:** 179 of 207 rows are *Candida albicans*, reflecting where biofilm
genetics has been done. This table cannot support a kingdom-wide biofilm predictor yet.
