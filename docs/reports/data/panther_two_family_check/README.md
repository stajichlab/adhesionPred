# PANTHER two-family hand check (2026-10-24)

Ids and metadata only (no sequences). Public read-only APIs, no login; raw answers cached outside the repo
(`~/.cache/adhesionPred-panther-twofamily`). No model was run. Report: `docs/reports/2026-10-24-panther-two-family-check.md`.
Scripts (stdlib only, `python3 -I`; `common.py` is the cached GET, copied from `panther_family_scoping/scripts` with a new cache folder and User-Agent):
`s1_family_tables.py REPO OUT`, `s2_pfam_overlap.py OUT_TSV`, `s3_kiss_mapping.py REPO OUT`. Run in that order. Re-running on another day can give other numbers (new UniProt or InterPro release).

## Tables

| File | Content |
|---|---|
| `family_totals.tsv` | Per family: proteins by domain of life and kingdom (UniProt lineage text). |
| `subfamily_sizes.tsv` | Per PANTHER subfamily: name, all/fungal/reviewed counts, fungal share, fungal fine clades with at least 5, bacteria, archaea, non-fungal eukaryotes. `(family only)` = in the family but without a subfamily. |
| `fungal_clades_by_subfamily.tsv` | Fungal proteins per fine clade and subfamily (clade rule of `analysis/wall_fusion_truth/map_kiss.py`). |
| `domain_counts.tsv` | Pfam and InterPro entries on the members (fungal, per subfamily, and non-fungal): count, scope size, fraction. Entries seen on at least 2 percent or 20 proteins. |
| `reviewed_members.tsv` | All Swiss-Prot members (all taxa): subfamily, gene, name, organism, length, Pfam, function and location comments (cut). |
| `kiss_in_family.tsv`, `kiss_mapping.tsv` | Every Kiss accession (data/controls clusters.tsv, controls.tsv, clusters_missed.tsv) found in the two families, with subfamily and the step 6.0 Pfam scan hits. |
| `pfam_vs_family_counts.tsv` | Fungal UniProt proteins carrying a Pfam model, and how many of them are in the family. |

## API calls (all 2026-10-24; User-Agent `adhesionPred-twofamily/1.0`; 0.25 s pause; retry with backoff)

- InterPro: `GET https://www.ebi.ac.uk/interpro/api/entry/panther/{FAM}` (name, description, counters),
  `GET .../entry/panther/{FAM}/?subfamilies` (subfamily names), `GET .../entry/{pfam|interpro}/{id}` (entry names, descriptions).
- UniProt: `GET https://rest.uniprot.org/uniprotkb/stream?query=xref:panther-{FAM}&format=tsv&fields=accession,reviewed,organism_name,organism_id,lineage,lineage_ids,xref_pfam,xref_interpro,protein_name,gene_names,length,xref_panther`
  (all taxa); the same with `(reviewed:true)` plus `cc_function,cc_subcellular_location`;
  `.../uniprotkb/search?query=(xref:pfam-{PF}) AND (taxonomy_id:4751) [AND (xref:panther-{FAM})]&size=1` (count from `x-total-results`);
  `.../taxonomy/search?query=tax_id:{id}` (clade names); `GET .../uniprotkb/J9VHS9.json` (one entry, read by hand).
- PubMed (by hand, titles and abstracts only): PMIDs 22737071, 28039134, 24901238, 35896747.

## Caveats

- UniProt PANTHER cross-references lag InterPro (PTHR46430: UniProt 5,247 proteins, InterPro 5,476; PTHR35192: 1,521 and 1,598). Counts here are UniProt counts. The search count (4,924 fungal) and stream count (4,935) of one query differ by under 1 percent.
- Domain fractions are UniProt cross-references of InterPro matches, not a new search.
