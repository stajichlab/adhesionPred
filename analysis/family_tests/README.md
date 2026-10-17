# Family tests (M4, descriptive)

Tests each family of `data/sorting_hat/family_table.tsv` (24 families; the 14 inactive ones are the target) against:

- positives: 108 proteins with `cls = adhesin` in `data/curated/adhesins/adhesins.tsv` (54 E1, 41 E2, 14 E3 after the owner moved two PA14 rows to E3 on 2026-10-09; sequences from UniProt);
- negatives: the shared hard-negative set (`data/controls/hard-negatives-shared/`, 137 N1 and 24 N2 rows; 76 N1 clusters);
- the 56 scan proteomes (raw Pfam tables of `_workdir/fungi_scan/`).

Method: `hmmfetch` the 24 models by name from Pfam 38.2 (`.../2026-01-27-Pfam38.2/Pfam-A.hmm`), `hmmsearch --cut_ga` on each FASTA, then `family_tests.py`.
Output: `docs/reports/data/sorting_hat/family_tests/family_tests.tsv`.

Limits. The E2 labels were made partly from the same Pfam families (family-defined, README of the curated table), so E2 hits are circular and are not sensitivity. The positives are not clustered, so E1 hit counts are protein counts. There are no estimated sensitivity or specificity numbers here.
