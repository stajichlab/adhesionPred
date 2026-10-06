# Draft curation table for task 03 (repeat mechanism label)

*Generated 2026-10-06 from `data/curated/adhesins/adhesins.tsv` and the MMseqs2 clusters of
`analysis/calibration_truth/c1_truth_count.py` (30% identity, 50% coverage). Draft input for
`docs/agent-tasks/03-repeat-mechanism-controls.md`. It is not a control set.*

`curation_table.tsv` holds the 186 rows with an accession: `adhesin`, `hard_negative` and
`surface_other_adhesion_phenotype`. The last five columns are **empty on purpose**:

| Column | To fill |
|---|---|
| `mechanism_label` | `repeat_avidity`, `single_interface`, `other` or `unknown` (adhesins only) |
| `mechanism_source_pmid` | The paper that gives the evidence for the label |
| `mechanism_quote` | The supporting sentence, short |
| `reviewer`, `review_date` | Who and when |

Rules: label from the paper, not from a repeat detector. Do not edit `adhesins.tsv`; record changes
in `manual_overrides.tsv`. `cluster_has_both_adhesin_and_hard_negative = yes` marks the 2 clusters
that need a decision.
