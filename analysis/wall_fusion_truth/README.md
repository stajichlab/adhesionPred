# Kiss 2019 truth set: feasibility stage (plan Task 5.1)

This folder holds the feasibility stage for the wall, adhesion and fusion truth set built from Kiss et al. 2019 (Nat Commun, PMID 31501435). No model is run. The report is `docs/reports/2026-10-18-kiss-feasibility.md`.

## Files

| File | Purpose |
|---|---|
| `map_kiss.py` | Resolves organism codes, maps rows to UniProt sequences, counts rows, clades and clusters, counts the second source. Uses only the Python standard library. |
| `results/query_cache.json` | Every NCBI and UniProt answer that the script used (URL to text). The script runs offline from this file. |
| `results/kiss_species_72.tsv` | The 72-genome table of Kiss Supplementary Data 1 (MOESM5), reduced to six columns. |
| `results/organism_resolution.tsv` | The 16 organism codes: species, NCBI taxon, lookup scope, lineage, clade. |
| `results/row_mapping.tsv` | One line per row of the 181-row table: slice, organism, clade, mapping status, reason, UniProt accession. |
| `results/counts_per_slice.tsv`, `loss_reasons.tsv`, `rows_per_organism.tsv` | Row counts, loss reasons and rows per organism. |
| `results/clusters_per_clade.tsv`, `floors_per_slice.tsv` | Clusters per clade and the clades that a floor of 3, 5 or 8 clusters admits, in three scenarios. |
| `results/reconciliation.tsv`, `data3_full_counts.tsv`, `family_rows_cluster_overlap.tsv` | Counts against report 06 and the full Data 3. |
| `results/second_source_uniprot_reviewed.tsv`, `second_source_summary.tsv` | Count of UniProt reviewed entries of the wall and fusion Pfam families, per clade. |
| `results/mapped_sequences.faa` | UniProt sequences of the mapped rows (small, derived). These are not JGI sequences. |

## How to run

```
python3 -I analysis/wall_fusion_truth/map_kiss.py                 # offline, from the cache
python3 -I analysis/wall_fusion_truth/map_kiss.py --online        # fill the cache (public NCBI and UniProt APIs only)
uv run --no-project --with openpyxl python3 -I analysis/wall_fusion_truth/map_kiss.py \
    --species-xlsx MOESM5.xlsx --data3-xlsx MOESM7.xlsx           # rebuild the two derived tables
```

Output order is fixed (sorted). A second offline run gives the same bytes. The supplement files are untrusted data: download them to a scratch folder, read them with `python3 -I`, never execute them. They are not in the repository.

## Rules of the mapping

- A JGI protein ID cannot be turned into a sequence by a public route (see the report, section 3). A row counts as MAPPED only when an identifier in the row itself (a locus tag such as Afu3g12700, NCU09332, BATDEDRAFT_20500) matches exactly one UniProt entry inside the genome strain of the paper (or the species, where the table names no strain).
- A gene name that matches one UniProt entry is only a CANDIDATE. It is never counted as mapped.
- A row without a JGI ID (Data 4 and Data 8 rows) is mapped to the A. fumigatus ortholog locus tag of its own `ortholog_Afu` column. It stands for the family. It is not the Kiss protein.
- The sequence must be 20 to 15000 residues long, use amino-acid letters, and have no run of 5 or more X.
- The independent unit is the Kiss `cluster_number`. It is a proxy. It is not a 30 percent identity cluster. The MMseqs2 clustering of D8 is Task 5.3 (HPCC) and is not done here.

## Tests

`tests/wall_fusion_truth/test_map_kiss.py` uses a 3-row fixture and a stub lookup. It needs no network.
