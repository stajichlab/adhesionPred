# PANTHER family scoping (2026-10-23)

Ids and metadata only (no sequences). Public read-only APIs, no login, cached raw answers kept outside the repo
(`~/.cache/adhesionPred-panther`). No model was run and nothing here says how well any model performs.
Report: `docs/reports/2026-10-23-panther-family-scoping.md`.

## Tables

| File | Rows | Content |
|---|---|---|
| `accession_to_panther.tsv` | 112 | One row per UniProt accession of (a) the 44 clusters missed by the Pfam panel (`clusters_missed.tsv`) and (c) the 112 control rows (`controls.tsv`); the two sets overlap (`source_set`). Columns: PANTHER family, subfamily, family name, PANTHER score (an E-value as reported by InterPro; 0 means underflow), all InterPro entries, Pfam entries. |
| `family_summary.tsv` | 61 | Per PANTHER family: name, InterPro entry, size in all taxa, fungal proteins (UniProt stream), reviewed and unreviewed fungal, fine clades with at least 5 proteins, top clades, Metazoa, Viridiplantae, Bacteria and Archaea counts, fungal share, `broad_or_specific` (count rule, see report section 6), fungal proteins with a panel Pfam, with a shipped Pfam, number of our accessions in the family. |
| `family_clade_counts.tsv` | 1,418 | Fungal proteins per fine clade per family (long format). |
| `family_member_overlap.tsv` | 89 | Our accessions (a)(c) that are in each family, with fine clade and UniProt status. |
| `pfam_to_panther.tsv` | 296 | For 69 Pfam families of the panel (30 candidates of the specificity test, 18 active rows of `family_table.tsv`, 21 wall/fusion models): fungal protein count, and the PANTHER families found over a sample (all reviewed fungal proteins plus the first 500 unreviewed ones), with fractions. |
| `missed44_cluster_panther.tsv` | 44 | The 44 missed clusters with the Pfam scan class of step 6.0 and the PANTHER family of their accessions (JGI-only clusters are listed with 0 UniProt accessions). |
| `controls_panther_coverage.tsv` | 3 | Controls with a PANTHER family, by slice. |
| `test_set_screen.tsv` | 60 | Mechanical screen of the streamed families (5 criteria, report section 6). |

`scripts/`: the scripts that were run, in order `s1` to `s4` (stdlib only, `python3 -I`; `common.py` is the cached GET).
Usage is in each docstring. Re-running gives different numbers when UniProt or InterPro have a new release.

## API calls (all 2026-10-23; User-Agent `adhesionPred-scoping/1.0`; 0.25 s pause between calls; retry with backoff)

1. `s1`: per accession, three InterPro calls:
   `GET https://www.ebi.ac.uk/interpro/api/entry/panther/protein/uniprot/{acc}` (family, subfamily, score),
   `GET https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{acc}`,
   `GET https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}`.
2. `s2`: per Pfam family, UniProt REST:
   `GET https://rest.uniprot.org/uniprotkb/search?query=(xref:pfam-{PF}) AND (taxonomy_id:4751)&fields=accession&size=1` (count from the
   `x-total-results` header; also with `AND (reviewed:true)`),
   `GET .../uniprotkb/stream?query=(xref:pfam-{PF}) AND (taxonomy_id:4751) AND (reviewed:true)&format=tsv&fields=accession,reviewed,organism_id,xref_pfam,xref_panther`,
   `GET .../uniprotkb/search?query=... AND (reviewed:false)&format=tsv&fields=...&size=500` (the sample of unreviewed ones).
3. `s3`: per PANTHER family with at least 10 percent of some Pfam sample, or named by an accession of (a)(c):
   `GET https://www.ebi.ac.uk/interpro/api/entry/panther/{PTHR}` (name, integrated InterPro entry, all-taxa counters);
   UniProt counts `xref:panther-{PTHR}` for all taxa, fungi (4751), NOT fungi, Metazoa (33208), Viridiplantae (33090), Bacteria (2), Archaea (2157);
   `GET .../uniprotkb/stream?query=(xref:panther-{PTHR}) AND (taxonomy_id:4751)&format=tsv&fields=accession,reviewed,organism_id,lineage_ids,xref_pfam`
   (not for PTHR33048, 66,799 fungal proteins: counts only); clade names `GET https://rest.uniprot.org/taxonomy/search?query=tax_id:...`.
   Fine clade rule: `analysis/wall_fusion_truth/map_kiss.py::clades_from_lineage` (Pezizomycotina and Basidiomycota to class, other Ascomycota to subphylum, other phyla to phylum), applied to the UniProt lineage ids.
4. `s4`: joins only (no calls). Also done once by hand: `GET https://www.pantherdb.org/services/oai/pantherdb/supportedgenomes` (144 reference genomes)
   and `.../geneinfo?geneInputList=P38248&organism=559292` (PANTHER v19 returned PTHR31018:SF3; no E-value field in that answer).

## Caveats

- The UniProt PANTHER cross-references lag the InterPro release (PTHR31018: InterPro 1,929 fungal proteins, UniProt stream 1,757). Family counts here are UniProt counts.
- UniProt search count and stream count of one query differ by under 1 percent; the stream count is used where a stream was fetched.
- The Pfam sample of `pfam_to_panther.tsv` is not random for large families (all reviewed plus the first 500 returned by UniProt): fractions are indications, not estimates.
- `broad_or_specific` is a count rule and does not read the function; see the report.
