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
| no_uniprot_entry | PM candidates with no matching UniProt entry (a subset of `pm_unresolved`). |
| gpi_feature_no_evidence | UniProt entries returned for the PM candidates that have at least one GPI-anchor Lipidation feature without an evidence code. Such a feature cannot make the gene P-gpi. |
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

## d8_uniprot/ (03_triage_pm.py)

Raw UniProt JSON pages, one file per page (`<source_id>_candidates_<n>_<page>.json` and
`<source_id>_gpi_<page>.json`). A run overwrites the pages that it fetches. It does not delete
pages from older runs, so stale pages can remain. Use the pages whose names match the current
run, or delete the folder before a run.

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

# Phase B outputs (features and embeddings)

Scripts 05, 06 and 07 and the jobs in `jobs/` write to `$STEP1_WORKDIR/phaseb/`. Later phases
join on `(source_id, gene_id)` and on `seq_sha256`. `emb_row` and `emb_cterm_row` give the rows
of the embedding matrices.

## sequence_members.tsv.gz (05_prepare_sequences.py, one row per protein of an input set)

| Column | Meaning |
|---|---|
| set_id | Input set in `sequence_sets.tsv` (`truth`, `uniprot_kw`, `<species>_proteome`, ...). |
| source_id | For `truth`: the `source_id` of `truth_sequences.tsv.gz`. Else equal to `set_id`. |
| gene_id | For `truth`: the GAF gene ID. For `uniprot_kw`: the UniProt accession. Else the first FASTA header token. |
| seq_sha256 | SHA-256 of the cleaned sequence (`seqhash.seq_sha256`). |
| length | Length of the cleaned sequence. |

## unique_sequences.tsv.gz (05_prepare_sequences.py, one row per unique sequence)

| Column | Meaning |
|---|---|
| row | Position in seq_sha256 order (0-based). Row of `emb/<model>.nterm.npy`. |
| seq_sha256 | SHA-256 of the cleaned sequence. Sort key. |
| length | Length of the cleaned sequence. |
| cterm_row | For sequences longer than 1,022 aa: row of `emb/<model>.cterm.npy` (0-based, in seq_sha256 order). Empty otherwise. |
| sequence | The cleaned sequence. Every character is in the ESM-2 alphabet. |

## unique_sequences.fasta.gz (05_prepare_sequences.py)

The unique sequences as FASTA, in `row` order. The header is the `seq_sha256`. J1 reads it.

## prepare_run.json (05_prepare_sequences.py)

**Keys:**

- `inputs` per set, see below
- `members`
- `unique_sequences`
- `unique_residues`
- `over_max_residues` unique sequences longer than 1,022 aa
- `all_sources` `true` when `sequence_run.json` and `keyword_tier_run.json` both say `all_sources: true`
- `truth_set_sha256` copied from `sequence_run.json` (empty if the log cannot be read); 07 stops if it differs from `d8_run.json`
- `git_commit`
- `python`
- `arguments`

**Keys of each `inputs` object:**

- `file`
- `sha256`
- `members`
- `empty_records`
- `dropped_duplicate_records` records whose `gene_id` occurred before in the same set with the same sequence (kept once; 0 for the truth set)

No time stamp.

## chunk_plan.tsv (06_plan_embedding.py, one row per embedding chunk)

| Column | Meaning |
|---|---|
| chunk_id | `<window>_<NNNN>`, numbered in plan order. |
| window | `nterm` (first 1,022 residues) or `cterm` (last 1,022 residues; sequences > 1,022 aa only). |
| n_seqs | Number of chunk members. |
| residues | Sum of the window lengths. |
| members_sha256 | SHA-256 of the member seq_sha256 values joined by newlines, in chunk order. |

## chunk_members.tsv.gz (06_plan_embedding.py, one row per chunk member)

| Column | Meaning |
|---|---|
| chunk_id | As in `chunk_plan.tsv`. |
| position | Row in the chunk array (0-based). |
| row | `row` of the sequence in `unique_sequences.tsv.gz`. |
| seq_sha256 | SHA-256 of the sequence. |

## job_plan.json (06_plan_embedding.py)

**Keys:**

- `models`
- `rates_residues_per_s` per model
- `model_load_s` per model
- `batch_size` per model
- `chunk_residues`
- `chunks`
- `residues_per_model`
- `unique_sequences`
- `total_seconds`
- `n_jobs`
- `seconds_per_job`
- `longest_job_seconds`
- `time_minutes`
- `rate_source` `J0` or `assumed`
- `throughput_sha256` empty when `rate_source` is `assumed`
- `unique_sequences_sha256`
- `git_commit`
- `python`
- `arguments`

**Keys added when `rate_source` is `J0`:**

- `j0_sample_sha256`
- `j0_device`
- `j0_gpu_name`

The formula is in the docstring of `chunk_plan.py`.

## j0/throughput.json (jobs/throughput_pilot.py, J0)

**Keys:**

- `schema` `step1-phaseb-j0-throughput/2`
- `device`
- `gpu_name`
- `torch`
- `torch_cuda`
- `esm`
- `n_proteins`
- `residues`
- `seed`
- `sample_sha256`
- `truth_sequences_sha256`
- `batch_sizes`
- `model_load_s` per model
- `runs` one object per model and batch size
- `best` per model
- `git_commit`
- `python`

**Keys of each `runs` object with status `ok`:**

- `model`
- `batch_size`
- `batch_failures` batches that `get_esm_embeddings` retried one sequence at a time, for example after out of memory
- `status` `ok`, `batch_failures`, or `skipped_<n>`
- `seconds`
- `proteins_per_s`
- `residues_per_s`
- `peak_mem_bytes` `null` on CPU
- `dim`

A run with status `skipped_<n>` has only `model`, `batch_size`, `batch_failures` and `status`.

**Keys of each `best` object:**

- `batch_size`
- `residues_per_s`
- `proteins_per_s`
- `seconds`

`best` holds, per model, the `ok` run with the highest `residues_per_s`.

## j0/gpu_cpu_diff.json (jobs/gpu_cpu_diff.py, J0)

**Keys:**

- `schema` `step1-phaseb-gpu-cpu-diff/1`
- `device_a`
- `device_b`
- `gpu_name`
- `torch`
- `n_windows`
- `windows_sha256`
- `batch_size`
- `models` per model, see below
- `git_commit`

**Keys of each `models` object:**

- `max_abs_diff`
- `mean_abs_diff`
- `max_rel_diff`
- `min_cosine`
- `repeat_identical_on_a`
- `repeat_max_abs_diff_on_a`

No threshold is applied.

## j0/nvidia_smi.csv (jobs/j0_pilot.sh, J0)

Output of `nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv` on the J0 node:
one header line and one line per GPU. It records the GPU and driver that J0 measured.

## signalp/part_NNN/ (jobs/j1_features.sh, J1)

SignalP 6 output for one part of `unique_sequences.fasta.gz`: `prediction_results.txt.gz`
(the done marker of the part), `output.gff3.gz`, `region_output.gff3.gz` (gzip `-n`), and
`input.sha256` (SHA-256 of the part FASTA). The format is described in `feature_parsers.py`.

## predgpi/part_NNN.tsv.gz (jobs/j1_features.sh, J1)

`jobs/predgpi_scores.py` output for one part, one row per sequence.

| Column | Meaning |
|---|---|
| id | The `seq_sha256`. |
| length | Sequence length. |
| gpi_call | `highly_probable`, `probable`, `weakly`, `none`, or `too_short` (40 aa or less). |
| gpi_prob | The PredGPI CLI score: 1.0, 0.70, 0.55, or 0. |
| omega | Omega site (GPI calls only). |
| fpr | PredGPI estimated false positive rate. Lower is more GPI-like. |
| svm | SVM output. |

`part_NNN.input.sha256` holds the SHA-256 of the part FASTA.

## features_unique.tsv.gz (07_build_features.py, one row per unique sequence)

| Column | Meaning |
|---|---|
| row, seq_sha256, length, cterm_row | As in `unique_sequences.tsv.gz`. |
| ser_thr_frac | (count of S + count of T) / length, 6 decimals. |
| sp_prediction | SignalP 6: `SP` or `OTHER`. Empty if SignalP gave no call. |
| sp_prob | SignalP 6 `SP(Sec/SPI)` probability. |
| sp_other_prob | SignalP 6 `OTHER` probability. |
| sp_cs_end | Last residue of the signal peptide (`CS pos: X-Y`, X). SP only. |
| sp_cs_prob | Probability of the cleavage site. SP only. |
| gpi_call, gpi_prob, gpi_omega, gpi_fpr, gpi_svm | PredGPI columns `gpi_call`, `gpi_prob`, `omega`, `fpr`, `svm`. |

## features.tsv.gz (07_build_features.py, one row per member)

| Column | Meaning |
|---|---|
| set_id, source_id, gene_id, seq_sha256, length | As in `sequence_members.tsv.gz`. |
| label, subset, stratum, d8_class, homology_only, role | From `truth_set_triaged.tsv.gz` for `truth` members. Empty for other sets. |
| ser_thr_frac, sp_prediction, sp_prob, sp_other_prob, sp_cs_end, sp_cs_prob, gpi_call, gpi_prob, gpi_omega, gpi_fpr, gpi_svm | As in `features_unique.tsv.gz`. |
| emb_row | Row of `emb/<model>.nterm.npy` (the `row` of the sequence). |
| emb_cterm_row | Row of `emb/<model>.cterm.npy`. Empty for sequences of 1,022 aa or less: use `emb_row` (their C-terminal window is the whole sequence). |

## feature_coverage.tsv (07_build_features.py, one row per set_id)

| Column | Meaning |
|---|---|
| set_id | Input set. |
| members | Members of the set. |
| unique_sequences | Distinct seq_sha256 values in the set. |
| no_signalp | Unique sequences of the set without a SignalP call. |
| no_predgpi | Unique sequences of the set without a PredGPI call. |
| over_1022 | Unique sequences of the set longer than 1,022 aa. |
| gpi_too_short | Unique sequences of the set of 40 aa or less (PredGPI gives no score). |

## features_run.json (07_build_features.py)

**Keys:**

- `tool_outputs_sha256` per J1 file, key = path under `phaseb/`
- `truth_set_sha256` the hash that `d8_run.json` and `sequence_run.json` share
- `input_sha256` see below
- `unique_sequences`
- `members`
- `missing_signalp`
- `missing_predgpi`
- `sp_predictions` counts
- `gpi_calls` counts
- `all_sources` `true` when `d8_run.json` and `prepare_run.json` say `all_sources: true`
- `git_commit`
- `python`
- `arguments`

**Keys of `input_sha256`:**

- `truth_set_triaged.tsv.gz`
- `sequence_members.tsv.gz`
- `unique_sequences.tsv.gz`

Script 07 reads `sequence_run.json`, `d8_run.json` and `phaseb/prepare_run.json`. It stops if the `truth_set_sha256` values of the three differ. It also stops if the `sha256` that `prepare_run.json` records for the truth sequences file (`inputs.<truth set id>.sha256`) differs from the file now. The truth set ids are the sets with kind `truth` in `sequence_sets.tsv` (`--sets`). Script 07 does not need J2 or the assembly step.

## emb/<model>/<chunk_id>.npy (jobs/embed_chunks.py, J2)

float32 array, one row per chunk member in `chunk_members.tsv.gz` order, ESM-2 layer 6,
mean over residue tokens (`surface_glyco.embeddings.get_esm_embeddings`). The sidecar
`<chunk_id>.json` is the done marker.

**Keys of the sidecar:**

- `window`
- `members_sha256`
- `repr_layer`
- `batch_size`
- `device`
- `seconds`
- `batch_failures` number of "Warning: batch failed" lines of `get_esm_embeddings` for this chunk (0 unless `--allow-batch-failures`)
- `chunk_id`
- `model`
- `shape`
- `dtype`
- `array_sha256` SHA-256 of the raw array bytes

## emb/<model>.nterm.npy (jobs/assemble_embeddings.py)

float32, shape (unique sequences, dim): row `row` is the N-terminal window embedding of that
sequence. `emb/<model>.cterm.npy` has shape (sequences > 1,022 aa, dim) and row `cterm_row`.
`assemble_embeddings.window_matrix` builds the full C-terminal matrix for M8-C and M35-C.
dim is 320 for `esm2_t6_8M_UR50D` and 480 for `esm2_t12_35M_UR50D`.

## emb/chunk_manifest.tsv (jobs/assemble_embeddings.py, one row per model and chunk)

| Column | Meaning |
|---|---|
| model | ESM-2 model name. |
| chunk_id, window, n_seqs, members_sha256 | As in `chunk_plan.tsv`. |
| array_sha256 | SHA-256 of the raw bytes of `emb/<model>/<chunk_id>.npy` (as in its JSON). |

## emb/embedding_run.json (jobs/assemble_embeddings.py)

**Keys:**

- `models` per model and window, see below
- `unique_sequences`
- `chunks`
- `unique_sequences_sha256`
- `chunk_plan_sha256`
- `git_commit`
- `python`

Evaluation scripts must check that `features_run.json` `input_sha256["unique_sequences.tsv.gz"]` equals `embedding_run.json` `unique_sequences_sha256`.

**Keys of each window object in `models`:**

- `shape`
- `dtype`
- `array_sha256`

# Phase C outputs (evaluation)

Scripts 08 to 12 in `phasec/` write to `$STEP1_WORKDIR/phasec/`. Each run JSON records
`outputs_sha256` for its other outputs; the next script stops when a file differs from it.
Tables join on `seq_sha256`. List columns hold sorted unique values, comma separated.
No Phase C job has been run on the real data; the tests run these scripts on a small fixture.

## phasec/eval_table.tsv.gz (08_build_eval_tables.py, one row per unique sequence)

| Column | Meaning |
|---|---|
| seq_sha256 | Hash of the sequence (Phase A `seqhash.seq_sha256`). |
| origin | `go` (truth genes) or `tc` (keyword tier T-c rows that survive ruling C-6). |
| class | `pos`, `neg` or `excluded` (`labelmap.class_of`). A hash with a `pos` or `neg` member and an excluded member is `excluded` (ruling C-14). |
| label | Truth labels of the members (list); `T-c` for T-c rows. |
| subset | `wall` or `extracellular-only` of the members (list); empty for T-c rows. |
| stratum | Reporting strata of the members (list, `labelmap.stratum_of`): `wall`, `extracellular-only`, `PM-TM`, `pm-unresolved`, `N-int`, `N-sec`, `ambiguous`; `T-c` for T-c rows. |
| d8_class | D8 classes of the members (list). |
| homology_only | `no` when any member has the label without homology codes (direct evidence); `yes` otherwise; empty for T-c rows. |
| internal_evidence_htp_only | `yes` when any member says so, else `no`; empty for T-c rows. |
| source_ids | Truth sources of the members (list); `T-c` for T-c rows. |
| gene_ids | Gene IDs of the members, or the UniProt accessions of the T-c rows (list). |
| species | Species (truth) or genome (T-c) names (list). |
| roles | Source roles from `species.tsv` (list); `tc` for T-c rows. |
| clades | `in_clade` of the sources, or the clade from `phasec/tc_taxon_clades.tsv` (list). |
| taxon_ids | NCBI taxon IDs (list). |
| length | Sequence length. |
| emb_row, emb_cterm_row | As in `features.tsv.gz`. |

## phasec/eval_literature.tsv (08_build_eval_tables.py, one row per literature seed with an accession)

| Column | Meaning |
|---|---|
| accession | UniProt accession from `uniprot_query` of `eurotiomycetes_seeds.tsv`. |
| gene, moonlighting, species, order | From `eurotiomycetes_seeds.tsv`. |
| lit_class | The seed file's `class` (`adhesin`, `hard_negative`). It describes adhesion, not location. |
| seq_sha256, length, emb_row, emb_cterm_row | From the `uniprot_kw` member; empty without a sequence. |
| literature_positive | `yes` when the row has a sequence and `moonlighting` is not `YES` (spec 3.2; `hard_negative` rows count, ruling C-9 amended); else `no`. |

## phasec/eval_dedupe_log.tsv (08_build_eval_tables.py, one row per dropped or reclassified member)

| Column | Meaning |
|---|---|
| origin | `go` or `tc`. |
| source_id, gene_id | Truth source and gene, or `T-c` and the accession. |
| seq_sha256 | Hash of the member. |
| class | Class of the member. |
| reason | `both_classes` (hash dropped), `class_and_excluded` (row becomes excluded), `go_label_wins` (T-c row dropped, ruling C-6). |
| detail | The classes of the GO members with this hash. |

## phasec/eval_sequences.fasta.gz (08_build_eval_tables.py)

One record per hash of `eval_table.tsv.gz` and `eval_literature.tsv`, sorted by hash. The header
is the `seq_sha256`. The MMseqs2 input of 09.

## phasec/build_run.json (08_build_eval_tables.py)

**Keys:**

- `all_sources` `true` (08 stops otherwise)
- `truth_set_sha256` shared by `d8_run.json`, `features_run.json`, `keyword_tier_run.json`
- `input_sha256` see below
- `go_members_by_source_class` GO members per source and class, before dedupe
- `labelled_genes_without_sequence` per source
- `table_rows_by_origin_class` key `<origin>:<class>`
- `log_rows_by_reason`
- `tc_rows` rows of `keyword_tier.tsv.gz`
- `tc_rows_dropped_by_go_class` T-c rows dropped by ruling C-6, by the classes of the GO members
- `go_members_sharing_a_tc_hash` key `<class>:<label>:<d8_class>`
- `literature_rows`
- `literature_positives`
- `literature_no_accession` seed genes without an accession
- `sequences` records in `eval_sequences.fasta.gz`
- `git_commit`
- `library_versions` Python, numpy, scipy, scikit-learn
- `arguments`
- `outputs_sha256` see below

**Keys of `input_sha256`:**

- `truth_set_triaged.tsv.gz`
- `features.tsv.gz`
- `features_unique.tsv.gz`
- `unique_sequences.tsv.gz`
- `keyword_tier.tsv.gz`
- `species.tsv`
- `eurotiomycetes_seeds.tsv`
- `tc_taxon_clades.tsv`

**Keys of `outputs_sha256`:**

- `eval_table.tsv.gz`
- `eval_literature.tsv`
- `eval_dedupe_log.tsv`
- `eval_sequences.fasta.gz`

## phasec/clusters.tsv.gz (09_make_splits.py, one row per sequence)

| Column | Meaning |
|---|---|
| seq_sha256 | Hash of a sequence of `eval_sequences.fasta.gz`. |
| cluster_id | Hash of the MMseqs2 cluster representative (`easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0`). |

## phasec/split_members.tsv.gz (09_make_splits.py, one row per split, fold, sequence and part)

A sequence can have two rows in one split with different parts (for example `test` and `test_lit`).

| Column | Meaning |
|---|---|
| split_id | `S1`, `S2-<source>`, `S3-<clade>` or `FULL` (`splits.py`). |
| fold | S1 fold 0 to 4; `0` for the other splits. |
| seq_sha256 | Hash of the sequence. |
| part | `train`, `train_tc` (V-kw only), `test`, `test_tc` (S1, scored, never truth), `test_lit`. |
| origin | `go`, `tc` or `lit`. |
| class | `pos`, `neg` or `excluded`. |
| cluster_id | As in `clusters.tsv.gz`. |

## phasec/tc_removed.tsv (09_make_splits.py, one row per T-c row removed from a split)

| Column | Meaning |
|---|---|
| split_id, fold | As in `split_members.tsv.gz`. |
| seq_sha256 | Hash of the T-c row. |
| gene_ids, taxon_ids | As in `eval_table.tsv.gz`. |
| rule | `a_test_protein`, `b_cluster_mate` (ruling C-5), `c_test_taxon` (ruling C-7). |

## phasec/max_identity.tsv.gz (09_make_splits.py, one row per S2 or S3 split and test sequence)

| Column | Meaning |
|---|---|
| split_id | An S2 or S3 split. |
| seq_sha256 | A test or literature sequence. |
| max_identity | Highest `fident` against the GO training proteins of the split (`easy-search -s 7.5 -c 0.5 --cov-mode 0`); empty when there is no hit. Self hits are skipped. |
| below_0.3 | `yes` when there is no hit or the identity is below 0.3 (ruling C-4). |

## phasec/logs/mmseqs.log (09_make_splits.py)

The output of the MMseqs2 calls. 09 writes it to `phasec/logs/` (on shared storage), also after
a STOP. It is not a hashed output and no later script reads it.

## phasec/splits_run.json (09_make_splits.py)

**Keys:**

- `all_sources`
- `truth_set_sha256`
- `input_sha256` see below
- `mmseqs_version` output of `mmseqs version`
- `mmseqs_commands` the two MMseqs2 commands without paths and threads
- `cluster_tsv_sha256` the raw MMseqs2 cluster table
- `seed` 20261001
- `sequences`
- `clusters`
- `splits` split ids in order
- `members_by_split_fold_part` key `<split>|<fold>|<part>`
- `tc_removed_by_split_fold_rule` key `<split>|<fold>|<rule>`
- `identity_below_0.3_by_split`
- `git_commit`
- `library_versions`
- `arguments`
- `outputs_sha256` see below

**Keys of `input_sha256`:**

- `eval_table.tsv.gz`
- `eval_literature.tsv`
- `eval_sequences.fasta.gz`
- `species.tsv`

**Keys of `outputs_sha256`:**

- `clusters.tsv.gz`
- `split_members.tsv.gz`
- `tc_removed.tsv`
- `max_identity.tsv.gz`

## phasec/scores.tsv.gz (10_fit_and_score.py, one row per unit, candidate and scored sequence)

A sequence that is in two scored parts of one unit is scored once. Its `part` lists the parts,
sorted and comma joined (for example `test,test_lit`).

| Column | Meaning |
|---|---|
| split_id, fold | The outer training set (unit). |
| variant | `V-go` or `V-kw`. |
| candidate | `B0` (log length), `B1` (amino-acid composition and log length), `R0`, `R1`, `R2`, `M8`, `M35`, `M8-C`, `M35-C`, `H`. |
| seq_sha256 | The scored sequence. |
| part | Its part in the split; `all` for FULL (every Phase B unique sequence). |
| score | Logistic-regression decision value; empty for rules. |
| prob | Platt probability (ruling C-10); empty for rules. |
| call | `1` when score >= the unit's threshold (Youden's J), or the rule call; else `0`. |

## phasec/scores_run.json (10_fit_and_score.py)

10 stops when a file that `build_run.json` or `splits_run.json` lists differs from its recorded
hash, or when `splits_run.json` `input_sha256` differs from `build_run.json` `outputs_sha256`
(the hash chain 08, 09, 10).

**Keys:**

- `all_sources`
- `truth_set_sha256`
- `input_sha256` see below
- `embedding_array_sha256` key `<model>.<window>`
- `seed` 20261001
- `candidates`
- `c_grid` values of C tried: 0.001, 0.003, 0.01, 0.1, 1, 10 (ruling C-12)
- `inner_folds`
- `units` key `<split>|<fold>|<variant>`, one object per candidate
- `git_commit`
- `library_versions`
- `arguments`
- `outputs_sha256` see below

**Keys of `input_sha256`:**

- `split_members.tsv.gz`
- `clusters.tsv.gz`
- `features_unique.tsv.gz`
- `unique_sequences.tsv.gz`

**Keys of `outputs_sha256`:**

- `scores.tsv.gz`

**Keys of a rule object in `units`:**

- `n_train`
- `n_train_pos`
- `g` PredGPI class cut (null for R0)
- `t` Ser+Thr cut, 0.20 to 0.40 in steps of 0.05 (null for R0 and R1)
- `j` Youden's J on the training rows
- `rule_grid` Youden's J on the training rows (no test data) for every grid cell: R1 `{g: J}` (3 numbers), R2 `{g: {t: J}}` with t keys `0.20` to `0.40` (15 numbers); null for R0

`rule_grid` (final review I-2) shows how J changes over the grid. The fitted cell is the cell
with the highest J (ties: the stricter cell). R2 calls a subset of the proteins that R0 calls,
so its recall and FPR cannot exceed those of R0 on the same rows; its J can.

**Keys of a logistic-regression object in `units`:**

- `n_train`
- `n_train_pos`
- `C`
- `h_variant` the ESM variant of H; null for the other candidates
- `threshold` decision value cut (Youden's J on the inner out-of-fold values)
- `platt_a`
- `platt_b`
- `inner_pr_auc` inner out-of-fold PR-AUC per setting (key `<C>`, or `<variant>:<C>` for H)
- `convergence_warnings`

## phasec/metrics.json (11_evaluate.py)

Values are objects `{value, lo, hi, n_defined}`: the point value, the 2.5th and 97.5th
percentiles over the resamples where the metric is defined, and their number. The default is
2,000 cluster-bootstrap resamples. `null` means not defined (for example recall without
positives). A truth block without negatives (the literature set) has `null` precision, PR-AUC,
precision at recall and precision at the rule's recall (spec 3.2: recall only). A stratum
without positives or without negatives has `null` for the metrics that need them.

**Keys:**

- `schema` `step1-phasec-metrics/1`
- `settings` resamples, seed, cut-offs, candidates, variants; see below
- `test_sets` key = test set name, see below
- `agreement` R2 against each ML candidate per proteome set and per truth class; see below
- `score_sources` per proteome set: counts of `oof`, `in_sample`, `final`
- `named_panel` one object per panel protein (parent spec section 6); see below
- `context` T-c rows of C. immitis and C. posadasii in V-kw training, per split and fold; see below
- `tc_removed` as `tc_removed_by_split_fold_rule` in `splits_run.json`
- `fitted_settings` key `<split>|<fold>|<variant>` (the keys of `units` in `scores_run.json`), one object per candidate; see below
- `grid_edge_counts` per grid parameter: `at_edge` (fitted values at the first or last grid value) and `n` (fitted values); see below

**Keys of `settings`:**

- `calibration_bins`
- `candidates`
- `ci_level`
- `estimate_half_width` recall half-width limit of the label `estimate` (ruling C-8)
- `estimate_min_direct_positives` count floor of 20 direct-evidence positives (ruling C-8)
- `fpr_level`
- `identity_cutoff`
- `long_cutoff`
- `n_resamples`
- `prevalences`
- `recall_levels`
- `saturation_auc`
- `seed` 20261001
- `variants`

**Keys of a test set object:**

- `split`
- `kind` `pooled`, `species`, `clade` or `literature`
- `truth` key `direct` or `all`, see below
- `recall_half_width` per label candidate (R2 and ML; V-go; direct truth)
- `fpr_half_width`
- `n_direct_positives` positives of the direct truth (stratum `all`)
- `floor_met` `true` when `n_direct_positives` is at least 20 (ruling C-8)
- `label` `estimate` (every recall half-width at most 0.10 and `floor_met`) or `smoke test` (ruling C-8)
- `max_recall_half_width`
- `max_fpr_half_width`
- `zero_width_recall_interval` label candidates whose recall interval has lo equal to hi
- `ambiguous` score distribution of the ambiguous genes
- `ambiguous_htp_only` the same for the high-throughput-only sub-stratum
- `lists` `pm-unresolved` and `P-gpi` genes with their calls
- `calibration` S2 test sets only: Brier score and 10 reliability bins per ML candidate

**Keys of a truth object:**

- `n` positives and negatives per stratum
- `metrics` per stratum, variant and candidate
- `variant_effect` V-kw minus V-go per stratum, candidate and metric (paired)
- `vs_rule` precision at the recall of R2 and recall at the FPR of R2 per ML candidate
- `prevalence` precision at assumed prevalence 0.01, 0.03, 0.05, 0.1 (assumed, not measured)
- `nsec_fpr_at_rule_recall` N-sec FPR at the recall of R2 and the paired differences

**Strata (keys of `n` and of `metrics` of a truth object; `STRATA` in `11_evaluate.py`):**

- `all` every row of the truth
- `wall` positives with subset `wall` (positives only)
- `extracellular-only` positives with subset `extracellular-only` (positives only)
- `N-int` negatives of stratum `N-int` (negatives only)
- `N-sec` negatives of stratum `N-sec` (negatives only)
- `PM-TM` negatives of stratum `PM-TM` (negatives only)
- `long` proteins longer than `long_cutoff` (1,022 aa)
- `identity_below_0.3` proteins whose maximum identity to the GO training proteins is below 0.3

Every test set has the first seven strata. `identity_below_0.3` exists only for the S2 and S3
test sets (also `S3-Eurotiomycetes:literature`), because 09 writes `max_identity.tsv.gz` rows
for S2 and S3 splits only; the S1 test sets do not have it.

**Keys of `context`:**

- `onygenales_tc_rows_in_vkw_training` key `<split>|<fold>`: number of T-c rows with taxon 246410 (C. immitis RS) or 443226 (C. posadasii) in the V-kw training of that unit; a unit without such rows is absent

**Keys of a rule object in `fitted_settings`:**

- `n_train`
- `n_train_pos`
- `g`
- `t`
- `j`
- `at_grid_edge`

**Keys of a logistic-regression object in `fitted_settings`:**

- `n_train`
- `n_train_pos`
- `C`
- `h_variant`
- `threshold`
- `platt_a`
- `platt_b`
- `inner_pr_auc`
- `convergence_warnings`
- `at_grid_edge`

The values are those of the unit object in `scores_run.json` (its `rule_grid` is not copied).

**Keys of `at_grid_edge` (only the parameters that are not null for the candidate):**

- `C` `true` when C is 0.001 or 10 (the ends of `c_grid`)
- `g` `true` when g is `highly_probable` or `weakly` (the ends of the g grid)
- `t` `true` when t is 0.2 or 0.4 (the ends of the t grid)

**Keys of `grid_edge_counts`:**

- `C`
- `g`
- `t`

**Keys of `agreement`:**

- `proteomes` per proteome set, variant and ML candidate: agreement counts with R2
- `truth` per class (`pos`, `neg`), variant and ML candidate: agreement counts with R2
- `score_source_counts` counts of `oof`, `in_sample`, `final` for `proteomes` and `truth`

**Keys of a named panel object:**

- `name`
- `id`
- `why`
- `found` `true` when the protein is in a proteome set

**Keys added when `found` is `true`:**

- `seq_sha256`
- `set_ids`
- `class`
- `label`
- `stratum`
- `literature_class`
- `score_source`
- `calls` key `<candidate>|<variant>`
- `scores` key `<candidate>|<variant>`

## phasec/findings.json (11_evaluate.py)

11 stops unless exactly 3 S2 test sets exist (`findings.N_S2`); finding (c) records `n_s2` and
`expected_s2`.

**Keys:**

- `schema` `step1-phasec-findings/1`
- `a_b1_not_saturated` B1 ROC-AUC below 0.99 on S1:all and every S2 test set (direct, V-go); see below
- `b_ml_beats_b1_and_r2_on_nsec_s1` per ML candidate: B1 and R2 N-sec FPR minus the ML N-sec FPR at the recall of R2, interval above 0; see below
- `c_same_under_s2` the statement of (b) on every S2 test set; see below

**Keys of `a_b1_not_saturated`:**

- `threshold`
- `b1_roc_auc` per test set
- `holds`

**Keys of `b_ml_beats_b1_and_r2_on_nsec_s1`:**

- `test_set` `S1:all`
- `candidates` one object per ML candidate, see below
- `holds_for` the ML candidates whose `holds` is `true`
- `n_resamples` the number of bootstrap resamples
- `min_defined_fraction` 0.95: a difference counts only when it is defined in at least this fraction of `n_resamples` (final review M-3)

**Keys of a candidate object in `b_ml_beats_b1_and_r2_on_nsec_s1`:**

- `B1` the difference as `{value, lo, hi, n_defined}`
- `R2` the same
- `n_defined` `{B1: n, R2: n}`: the resamples in which each difference is defined
- `beats_B1` `true` when `lo` of `B1` is above 0 and `n_defined` of `B1` is at least `min_defined_fraction` x `n_resamples`
- `beats_R2` the same for `R2`
- `holds` both `true`

**Keys of `c_same_under_s2`:**

- `test_sets` one object of the form of (b) per S2 test set
- `n_s2` the number of S2 test sets found
- `expected_s2` 3 S2 test sets
- `holds_for` ML candidates in `holds_for` of every S2 test set (empty unless `n_s2` equals `expected_s2`)

## phasec/proteome_calls.tsv.gz (11_evaluate.py, one row per protein of every proteome set)

| Column | Meaning |
|---|---|
| set_id, source_id, gene_id, seq_sha256 | As in `sequence_members.tsv.gz`; sets with kind `download` or `site` in `sequence_sets.tsv`. |
| class, label, origin | From `eval_table.tsv.gz` when the hash is there; else empty. |
| score_source | `oof` (hash in an S1 fold), `in_sample` (in a FULL training table but in no S1 fold), `final` (in no training table). |
| call_<candidate>_<variant> | `1` or `0`. |
| score_<candidate>_<variant> | Decision value; empty for rules. |

## phasec/evaluate_run.json (11_evaluate.py)

11 stops when a file that `build_run.json`, `splits_run.json` or `scores_run.json` lists differs
from its recorded hash, or when the hash chain 08, 09 (`input_sha256` against `outputs_sha256`)
or 09, 10 is broken.

**Keys:**

- `all_sources`
- `truth_set_sha256`
- `input_sha256` see below
- `seed` 20261001
- `n_resamples`
- `git_commit`
- `library_versions`
- `arguments`
- `outputs_sha256` see below

**Keys of `input_sha256`:**

- `scores.tsv.gz`
- `split_members.tsv.gz`
- `clusters.tsv.gz`
- `max_identity.tsv.gz`
- `eval_table.tsv.gz`
- `eval_literature.tsv`
- `sequence_members.tsv.gz`
- `sequence_sets.tsv`
- `species.tsv`

**Keys of `outputs_sha256`:**

- `metrics.json`
- `findings.json`
- `proteome_calls.tsv.gz`

## phasec/report.md (12_report.py)

The owner's report. Every number in it is a value of `metrics.json` or `findings.json`.
12 checks this before it writes. The check is a guard against typing errors and has cell-level
tests; it cannot detect a value that comes from another cell. The report has a decision table of
all candidates per test set, a glossary, fixed caveat lines, the finding tables and the table
`Fitted settings per unit` with a line that counts the settings at a grid edge; it prints `n/a`
for a `null` value. A finding-table interval that rests on fewer resamples than `n_resamples`
shows `(n_defined k of B)`.

## phasec/report_run.json (12_report.py)

**Keys:**

- `all_sources`
- `truth_set_sha256`
- `input_sha256` see below
- `git_commit`
- `library_versions`
- `arguments`
- `outputs_sha256` see below

**Keys of `input_sha256`:**

- `metrics.json`
- `findings.json`

**Keys of `outputs_sha256`:**

- `report.md`

## phasec/logs/wall.<job>.txt (phasec/c1_evaluate.sh)

Two kinds of line per step of job C1: `C1 step <step> start seconds=<seconds>` before the step
and `C1 step <step> wall_seconds=<seconds>` after it. A step that the time limit kills has a
`start` line and no `wall_seconds` line. The file is copied on every exit of the job.
