# Curated Ecm33-family orthologs (count, 2026-10-24)

Ids and metadata only (no sequences). Report: `docs/reports/2026-10-24-ecm33-curated-ortholog-count.md`.
Nothing here scores a model, and no call, tier, status or panel row depends on it.

## Files

- `curated_orthologs.tsv`: 21 rows = 13 Swiss-Prot reviewed entries of PANTHER PTHR31018 (set `swissprot_PTHR31018`) and 8 literature-named
  proteins with an accession (set `literature_named`; five of them are also Swiss-Prot entries, so the table lists them twice on purpose,
  once per set). Columns: `accession, set, organism, clade_fine, reviewed, pmid, how_resolved, uniref50, uniref50_is_training_member,
  uniref50_is_independent_test_cluster, separated_30pct, independent`.
  - `uniref50_is_training_member`: the UniRef50 cluster id is a `hit_id` of `docs/reports/data/pr6_ecm33_build_2026-10-23/members/members.tsv`
    (a protein in that cluster is at least 50 percent identical to a training member, so it is far above 30 percent).
  - `uniref50_is_independent_test_cluster` and `separated_30pct`: the cluster is one of the 187 test clusters of
    `docs/reports/data/pr6_ecm33_independent_eval_2026-10-24/extract/test_clusters.tsv`, and whether `separation.tsv` says it shares no
    30 percent identity / 50 percent coverage MMseqs2 cluster with a training member. `not_applicable` when the cluster is a training member.
  - `independent`: `yes` only if the cluster is not a training member and is separated at 30 percent. Result: 0 of 21 rows.
- `scripts/build_table.py`: builds the table from the committed tables (stdlib; `python3 -I scripts/build_table.py REPO_ROOT OUT_TSV`).

## Exact API calls (all read-only, no login, 2026-10-24, UniProt release 2026_03 of 02-September-2026; at least 1.5 s between calls)

1. Swiss-Prot entries of the family (answer: 13 entries, header `x-total-results: 13`):
   `GET https://rest.uniprot.org/uniprotkb/search?query=(taxonomy_id:4751) AND (reviewed:true) AND (xref:panther-PTHR31018)&fields=accession,organism_id,organism_name,protein_name,gene_names,length,xref_pfam,xref_panther&format=tsv&size=500`
2. Cross-check by name and by Pfam/InterPro (answer: 15 entries = the same 13 plus two namesakes that are not family: Q09750 pst1 of S. pombe, a SIN3 homolog, and A0A1D8PHR5 PST1 of C. albicans, a quinone oxidoreductase):
   `GET https://rest.uniprot.org/uniprotkb/search?query=(taxonomy_id:4751) AND (reviewed:true) AND (xref:pfam-PF12454 OR xref:interpro-IPR024238 OR xref:interpro-IPR051648 OR gene:ECM33 OR gene:SPS2 OR gene:PST1 OR gene:SPS22 OR protein_name:ecm33)&fields=accession,organism_name,protein_name,gene_names,xref_panther,xref_pfam&format=tsv&size=500`
3. Literature-named proteins: `accession:J4UGA7`, `accession:G9MJ35`, `accession:N1RSC8` (single-accession lookups taken from the 2026-10-23 scoping), plus name searches
   (`organism_name:"Metarhizium robertsii" AND protein_name:ecm33`, `organism_name:"Aspergillus flavus" AND (protein_name:ecm33 OR gene:AFLA_113120 OR protein_name:"cell wall protein")`,
   `organism_name:"Podosphaera xanthii" AND protein_name:ecm33`, `organism_name:"Nakaseomyces glabratus" AND protein_name:ecm33`), all against `https://rest.uniprot.org/uniprotkb/search` with `format=tsv`. The three unresolved organisms returned no usable Ecm33 entry (see the report).
4. UniRef50 cluster of N1RSC8: `GET https://rest.uniprot.org/uniref/search?query=(identity:0.5) AND (uniprot_id:N1RSC8)&format=tsv&fields=id` (answer `UniRef50_A0A428P255`).
5. PMIDs were re-read on 2026-10-24 from PubMed (title and organism checked): 24549768, 30840111, 35183338, 15583168, 16804168, 20032302, 29696338, 42325904.

Raw answers are cached outside the repo in `/Users/jstajich/projects/adhesionPred-ecm33-curated-cache/` (`q1.tsv`, `q2.tsv`, `n1.tsv`, and the response headers `h1.txt`, `h2.txt`).
