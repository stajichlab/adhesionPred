# Family-matched non-allergen negatives (all-kingdom allergens)

Built 2026-10-10. Plant and animal groups only. Fungal sets are not used.

## Query

For each (Pfam family, group) pair, the query sent to UniProt REST (`/uniprotkb/search`, fields `accession,length,organism_name,protein_name,reviewed`):

```
reviewed:true AND xref:pfam-{PF} AND NOT keyword:KW-0020 AND NOT cc_allergen:allergic
AND NOT antigen AND NOT cancer AND NOT taxonomy_id:9606 AND NOT taxonomy_id:10239
AND length:[10 TO 850] AND {TAX}
```

- plant: `TAX` = `taxonomy_id:33090` (Viridiplantae)
- animal: `TAX` = `(taxonomy_id:33208 AND NOT taxonomy_id:9606)` (Metazoa, not human)
- `antigen` and `cancer` are free-text terms, as written in the task.
- `taxonomy_id:10239` removes viruses.

## Groups

Each positive comes from `allergenonline_all/allergenonline_new.tsv`. Its AllergenOnline types come from `allergenonline_all/allergenonline_fungi.tsv` (column `uniprot`).
- plant: any type in {Aero Plant, Food Plant, Gliadin}. 26 positives have two types; the plant rule applies if one of them is plant.
- animal: all other positives. This includes Bacteria, Protozoan and Vaccine types, which get Metazoa sampling by this rule.

Positive families come from UniProt `xref_pfam` for each positive (batches of 40, OR query). A family used by both groups is sampled separately for each group.

## Filters applied after the query

- Remove every positive accession (1679).
- Remove every accession in `allkingdom/kingdom_negatives.tsv` (1669 accessions in that file).
- Remove any negative already picked by an earlier family (sorted order).

## Sampling

- Cap: 10 negatives per family per group. A pool of 10 or fewer gives all of its members.
- Seed: 20261010, one `random.Random` stream.
- Order: (family, group) pairs sorted by Pfam accession, then group. Each pool is sorted by accession before sampling.

## Outputs

- `family_negatives.faa`: sequences, IDs `K_<accession>`.
- `family_negatives.tsv`: one row per negative.
- `counts.md`: counts per group and family, cluster counts, positives without Pfam.
- `cluster/`: MMseqs2 `easy-cluster` output on `all_input.faa` (positives with IDs `P_<accession>`, then negatives `K_<accession>`).
- `raw/`: Pfam per positive, candidates before exclusion, per-pair stats and the step log.
- `scripts/`: the code that produced these files (`common.py`, `step1_pfam.py`, `step2_select.py`, `step5_checks.py`, `step6_counts.py`).

Family names come from the InterPro REST entry for each Pfam accession.

## Clustering

`mmseqs easy-cluster all_input.faa clu <tmp> --min-seq-id 0.3 -c 0.5 --cov-mode 0 -s 7.5 --threads 4` with MMseqs2 build 17-b804f (non-AVX2).
