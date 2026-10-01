# Columns of the D1 outputs

`01_extract_go_truth.py` writes `truth_set.tsv.gz`, `counts.tsv` and `extract_log.json`.

Note: the `nohom_*` columns of `counts.tsv` reproduce `d1_count.py` (labels recomputed without
homology codes). They are NOT used as truth. The `direct_*` columns are the direct-evidence truth.

## truth_set.tsv.gz (one row per gene)

| Column | Meaning |
|---|---|
| source_id | Identifier of the source in `species.tsv`. |
| species | Species name from `species.tsv`. |
| taxon_id | NCBI taxon id from `species.tsv`. |
| in_clade | Clade label from `species.tsv`. |
| role | Role of the source (train, test_clade, alternate_file, undecided). |
| gene_id | Gene identifier in GAF column 2. |
| symbol | Gene symbol in GAF column 3. |
| synonym1 | First synonym in GAF column 11. |
| label | Truth label under the `non_iea` policy (P-ext, N-int, N-sec or ambiguous). |
| subset | For P-ext: `wall` or `extracellular-only`; empty otherwise. |
| stratum | `subset` if set, else `label`. |
| tier | Evidence tier of the source (T-a for GO). |
| label_no_homology | Label recomputed without homology evidence codes (traceability only). |
| label_experimental | Label recomputed with experimental evidence codes only. |
| homology_only | `yes` if `label` differs from `label_no_homology`, else `no`. |
| pm_candidate | `yes` if the gene is a plasma-membrane candidate. |
| evidence_codes | Sorted evidence codes of the gene's cellular component rows. |
| surface_evidence | Non-IEA evidence codes on surface terms. |
| internal_evidence | Non-IEA evidence codes on internal terms. |
| internal_evidence_htp_only | `yes` if all internal evidence codes are high-throughput codes. |
| secretory_evidence | Non-IEA evidence codes on secretory-pathway terms. |
| source_file | GAF file name. |
| source_sha256 | SHA-256 of the GAF file. |
| source_date | `date-generated` header of the GAF file. |
| obo_sha256 | SHA-256 of the GO OBO file. |

## counts.tsv (one row per source)

| Column | Meaning |
|---|---|
| source_id | Identifier of the source. |
| primary_db | Most frequent database prefix in GAF column 1 (rows of other databases are dropped). |
| genes_cc | Genes with at least one kept cellular component row. |
| genes_noniea_cc | Genes whose evidence codes are not only IEA. |
| p_ext | Genes labelled P-ext. |
| p_ext_wall | P-ext genes in the wall subset. |
| p_ext_extonly | P-ext genes in the extracellular-only subset. |
| n_int | Genes labelled N-int. |
| n_sec | Genes labelled N-sec. |
| ambiguous | Genes labelled ambiguous. |
| pm_candidates | Genes flagged as plasma-membrane candidates. |
| cc_iea_triples | Distinct (gene, term, evidence) cellular component triples with IEA. |
| cc_triples | All distinct cellular component triples. |
| cc_iea_frac | `cc_iea_triples / cc_triples`. |
| all_aspects_iea_frac | IEA fraction of distinct triples over all three GO aspects. |
| obsolete_rows | GAF rows dropped for an obsolete term. |
| unknown_term_rows | Kept GAF rows whose term is absent from the OBO file. |
| exp_p_ext, exp_n_int, exp_n_sec, exp_ambiguous | Label counts with experimental evidence only. |
| nohom_p_ext, nohom_n_int, nohom_n_sec, nohom_ambiguous | Label counts without homology codes. Reproduce `d1_count.py`. NOT truth. |
| direct_p_ext, direct_n_int, direct_n_sec, direct_ambiguous | Genes whose `label` is X and `homology_only` is `no`. This is the direct-evidence truth. |
| ambiguous_htp_only | Ambiguous genes whose internal evidence is high-throughput only. |
