# Columns of the step 1 outputs

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

## extract_log.json (from 01_extract_go_truth.py)

`all_sources` is `true` when the run processed every source in `species.tsv`. It is `false` after a run with `--sources` that selected a subset. `sources` lists the processed sources (one object each, with `source_id`).

## truth_sequences.tsv.gz (02_attach_sequences.py, one row per matched gene)

| Column | Meaning |
|---|---|
| source_id, gene_id, label | As in `truth_set.tsv.gz`. |
| fasta_id | First token of the FASTA header. |
| length | Length of the cleaned sequence. |
| seq_sha256 | SHA-256 of the cleaned sequence (upper case, no whitespace, `J` to `L`, no `*`). |
| sequence | The cleaned sequence. |

The file does not exist after a STOP. A STOP happens when a source matches no gene, when a P-ext or ambiguous gene has no usable sequence, or when a sequence has a non-ASCII character.

## unmatched_ids.tsv (02_attach_sequences.py)

| Column | Meaning |
|---|---|
| source_id, gene_id, symbol, synonym1, label | As in `truth_set.tsv.gz`. |
| reason | `no_fasta_record` (no FASTA record has the gene key) or `empty_sequence` (the record has no residue after cleaning). |

## sequence_counts.tsv (02_attach_sequences.py, one row per source)

| Column | Meaning |
|---|---|
| source_id | Identifier of the source. |
| truth_set_sha256 | SHA-256 of the `truth_set.tsv.gz` that was read. |
| fasta_file, fasta_sha256 | FASTA file name and its SHA-256. |
| fasta_duplicate_records | Identical duplicate FASTA records that were merged. |
| genes, matched, unmatched | Gene counts. |
| unmatched_p_ext, unmatched_ambiguous, unmatched_n_int, unmatched_n_sec | Unmatched genes by label. |

## sequence_run.json (02_attach_sequences.py)

Keys: `sources` (selected source ids), `all_sources` (`true` when all sources of `species.tsv` were selected), `truth_set_sha256`, `fasta_sha256` (per source), `git_commit` (or `unknown`), `python`, `arguments`. It has no time stamp. A run with `--sources` replaces the full outputs with partial ones. Check `all_sources` before you use `truth_sequences.tsv.gz`. 02 refuses a `truth_set.tsv.gz` unless `extract_log.json` says `all_sources: true`; `--allow-partial-truth-set` overrides this.

## d8_triage.tsv (03_triage_pm.py, one row per PM candidate)

A PM candidate is a gene with `pm_candidate == "yes"` in `truth_set.tsv.gz`. Matched UniProt
entries are listed in one order in all list columns.

| Column | Meaning |
|---|---|
| source_id, gene_id, symbol | As in `truth_set.tsv.gz`. |
| d8_class | `P-gpi`, `PM-TM` or `pm-unresolved` (rules in `d8_triage.py`). |
| d8_reason | Why the class was given (for example the accession and the ECO codes). |
| uniprot_accessions | Accessions of the matched UniProt entries, separated by `,`. Empty if no entry matched. |
| reviewed | `yes` or `no` for each matched entry, separated by `,`. |
| gpi_eco | For each entry, the ECO codes of its GPI-anchor Lipidation features (joined by `,`); entries are separated by `;`. |
| tm_count | For each entry, the number of Transmembrane features; entries are separated by `;`. |
| tm_eco | For each entry, the ECO codes of its Transmembrane features (joined by `,`); entries are separated by `;`. |

## d8_gpi_outside_pext.tsv (03_triage_pm.py, one row per gene and UniProt entry)

Genes whose `label` is not `P-ext` but that match a reviewed UniProt entry of the organism with a
GPI-anchor Lipidation feature that has evidence ECO:0000269. The script does not relabel them.

| Column | Meaning |
|---|---|
| source_id, gene_id, symbol, label | As in `truth_set.tsv.gz`. |
| uniprot_accession | The matched UniProt entry. |
| gpi_eco | ECO codes of the GPI-anchor features of that entry, separated by `,`. |

## d8_counts.tsv (03_triage_pm.py, one row per source with at least one PM candidate)

| Column | Meaning |
|---|---|
| source_id | Identifier of the source. |
| pm_candidates | PM candidates of the source. |
| p_gpi, pm_tm, pm_unresolved | PM candidates in each D8 class. |
| organism_curated_gpi_entries | Entries from the organism GPI query that are reviewed and have a GPI-anchor feature with ECO:0000269. |
| curated_gpi_outside_pext | Rows in `d8_gpi_outside_pext.tsv` for the source. |
| uniprot_release | `X-UniProt-Release` header of the first candidate query response (of the organism query if that is empty). |

## truth_set_triaged.tsv.gz (03_triage_pm.py, one row per gene)

All columns of `truth_set.tsv.gz`, plus `d8_class` and `d8_reason` (as in `d8_triage.tsv`; empty
for genes that are not PM candidates).

- For PM candidates, `stratum` is overwritten with `d8_class`. `label` stays `P-ext`.
- A later phase must therefore filter on `d8_class` (`P-gpi`, `PM-TM`, `pm-unresolved`) to use
  the triage, not on `label`.
- `truth_set.tsv.gz` is not changed. It keeps the GO-only `stratum`.

## d8_run.json (03_triage_pm.py)

Keys: `sources` (selected source ids), `all_sources` (`true` when all sources of `species.tsv`
were selected), `truth_set_sha256` (SHA-256 of the `truth_set.tsv.gz` that was read; nothing
compares it later), `uniprot_release` (per source with PM candidates), `git_commit`, `python`,
`arguments`. It has no time stamp.

## keyword_tier.tsv.gz (04_build_keyword_tier.py, one row per kept T-c protein)

| Column | Meaning |
|---|---|
| accession, gene, genome, taxon_id | Copied from `data/curated/surface/surface.tsv`. |
| length | Length of the cleaned sequence. |
| seq_sha256 | SHA-256 of the cleaned sequence (`seqhash.seq_sha256`). |
| tier | Always `T-c`. |

## keyword_tier_removed.tsv (04_build_keyword_tier.py, one row per removed T-c protein)

`reason` is the first matching rule: `literature_accession`, `spombe_taxon`, `heldout_accession`, `no_sequence`, `literature_hash`, `heldout_hash`. `matched` names the matched protein (seed gene, taxon id, or `source_id:gene_id`). Held-out proteins are the labelled genes (`ambiguous` included, `unlabelled` excluded) of every source whose role is not `train`.

Limitation: removal is by accession and by exact cleaned-sequence hash only. Near-identical orthologs or paralogs under other accessions remain in the training table. Homology clustering in the later dataset plan must remove them.

## keyword_sequences.fasta.gz (04_build_keyword_tier.py)

UniProtKB FASTA records of the accessions in `surface.tsv` and the literature seeds, as returned
by UniProt (gzip). 04 downloads it with `--fetch` only when it is missing. Its SHA-256 and
UniProt release are in `keyword_sequences.json`.

## keyword_sequences.json (04_build_keyword_tier.py)

Keys: `sha256` of `keyword_sequences.fasta.gz`, `uniprot_release`, `accessions` (sorted list of the requested accessions), `accessions_sha256` (SHA-256 of that list joined by newlines), `requested`, `returned`, `fetched_on`. A cache is refused if the list differs from the current `surface.tsv` and seeds file.

## keyword_tier_run.json (04_build_keyword_tier.py)

Keys: `all_sources` (`true` only when `extract_log.json` and `sequence_run.json` both say `all_sources: true`), `truth_set_sha256` (SHA-256 of the current `truth_set.tsv.gz`; 04 stops unless it equals the value in `sequence_run.json`), `inputs_sha256`, `uniprot_release`, `kept`, `removed`, `missing_sequences`, `seeds` (per literature-seed accession: `gene`, `in_keyword_tier`, `has_sequence`), `git_commit`, `python`, `arguments`. It has no time stamp. The run STOPs when a requested accession or a seed has no sequence, unless `--allow-missing-sequences` or `--allow-missing-seed-sequences` is given (the second disables the hash rule for those seeds). A seed that is not a T-c row needs no hash rule for removal of itself, but its sequence still defines the literature hash.
