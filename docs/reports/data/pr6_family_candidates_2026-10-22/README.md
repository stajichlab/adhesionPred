# PR 6: candidate Pfam families for the specificity tests (2026-10-22)

Derived offline and deterministically from the committed step 6.0 record
(`docs/reports/data/pr6_pfam_scan_2026-10-22/scan_per_cluster.tsv`) and `docs/reports/data/kiss_baseline_2026-10-22/clusters_missed.tsv`:

    python3 -I analysis/wall_fusion_truth/pr6_family_specificity.py candidates --out <new dir>

Rule (owner decision 1 of 2026-10-22; both limits inclusive): a `hit_other_pfam` cluster is a candidate when one of its
hits covers at least 0.5 of the sequence OR has an E-value of 1e-30 or better. A family is a candidate when one of ITS
hits meets a rule. The clusters of class `hit_panel_model` and `no_pfam_hit` are not candidates.

| File | What it is |
|---|---|
| `candidate_clusters.tsv` | the 26 candidate clusters, their candidate families and the rule met |
| `candidate_families.tsv` | the 30 candidate families, their clusters, the rule met, and whether the family is already a row of `family_table.tsv` |
| `family_members.tsv` | the family's own evidenced members (all members of its candidate clusters), the denominator of the specificity test (D19) |
| `families.txt` | the model NAMES for `hmmfetch -f` |

The job `analysis/wall_fusion_truth/pr6_family_specificity.sbatch` makes these tables again and stops if one differs from this
directory. This is a list of families to TEST. It is not a panel addition, and it changes no tier or `narrow_ok`.
