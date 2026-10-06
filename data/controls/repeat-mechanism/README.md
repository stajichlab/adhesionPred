# Draft curation table for task 03 (repeat mechanism label)

*Generated 2026-10-06 from `data/curated/adhesins/adhesins.tsv` and the MMseqs2 clusters of
`analysis/calibration_truth/c1_truth_count.py` (30% identity, 50% coverage). Draft input for
`docs/agent-tasks/03-repeat-mechanism-controls.md`. It is not a control set.*

`curation_table.tsv` holds the 186 rows with an accession: `adhesin`, `hard_negative` and
`surface_other_adhesion_phenotype`. The last six columns are **empty on purpose**:

| Column | To fill |
|---|---|
| `mechanism_label` | One class code of `docs/TOOL-ARCHITECTURE.md` section 2: `2a` repeat/avidity, `2b-i` CFEM, `2b-ii` small Cys-knot, `2b-iii` Bys1, `2c` hydrophobin, `2d` moonlighting; or `other`, or `unknown` (adhesins only) |
| `label_confidence` | `high` (the paper states the mechanism), `medium` (the paper implies it), `low` (inferred from the family) |
| `mechanism_source_pmid` | The paper that gives the evidence for the label |
| `mechanism_quote` | The supporting sentence, short |
| `reviewer`, `review_date` | Who and when |

Rules: label from the paper, not from a repeat detector. Do not edit `adhesins.tsv`; record changes
in `manual_overrides.tsv`. `cluster_has_both_adhesin_and_hard_negative = yes` marks the 2 clusters
that need a decision.

## Comparing two independent fills

```
python3.12 analysis/calibration_truth/compare_mechanism_labels.py \
  --agent-a data/controls/repeat-mechanism/curation_table.agent-A.tsv \
  --agent-b data/controls/repeat-mechanism/curation_table.agent-B.tsv \
  --out data/controls/repeat-mechanism/comparison
```

It writes `comparison.tsv`, `review_queue.tsv` and `summary.md`. The quote check is a keyword check.
A fail means that a person reads the row. A pass does not prove the label.
