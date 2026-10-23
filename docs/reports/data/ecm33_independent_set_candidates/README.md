# Ecm33 independent-set candidates (scoping, 2026-10-23)

Ids and metadata only (no sequences). Derived from public read-only APIs on 2026-10-23. No model was run and nothing here
says anything about the performance of a model. Report: `docs/reports/2026-10-23-ecm33-independent-set-scoping.md`.

## Files

- `candidates.tsv`: 1,929 fungal UniProt accessions that InterPro assigns to the PANTHER family PTHR31018 (InterPro
  IPR051648, "Cell Wall Integrity and Assembly Regulator": the Ecm33 / Pst1 / Sps2 / Sps22 family). Columns:
  `accession, organism, taxid, clade_fine, source, evidence, tier, panther_subfamily, panther_evalue, uniref50,
  overlap_with_training`.
  - `clade_fine`: the fine level of `analysis/wall_fusion_truth/map_kiss.py::clades_from_lineage` (the subphylum for Ascomycota outside Pezizomycotina,
    e.g. Saccharomycotina and Taphrinomycotina; the class for Pezizomycotina and Basidiomycota; the phylum otherwise), computed from the UniProt taxonomy lineage. The same rule reproduces `clade_fine` of all 481 rows of
    `docs/reports/data/pr6_ecm33_build_2026-10-23/members/members.tsv` (0 mismatches).
  - `tier`: `T1_swissprot` (UniProt reviewed, 13), `T2_named` (UniProt protein name contains Ecm33, Sps2, Sps22, Pst1,
    "GPI-anchored cell wall organization" or "cell surface GPI-anchored"; 892; unreviewed names, many set by automatic
    annotation), `T3_panther_only` (1,024; usually "Uncharacterized protein" or "Receptor L-domain domain-containing
    protein").
  - `panther_subfamily`: SF3 (Ecm33 / Pst1 type, 1,831), SF12 (Sps2 / Sps22 type, 93), SF11 (S. pombe ecm33 type, 5).
  - `panther_evalue`: the PANTHER HMM E-value reported by InterPro for the family match.
  - `uniref50`: UniRef50 cluster of the accession. `overlap_with_training`: `in_training_cluster` when that cluster id is
    a `hit_id` of `members.tsv` (the 481 training members, including the query Q5AGC4), otherwise `new_cluster`.
- `stats.txt`: counts per clade and overlap (output of `scripts/build_candidates.py`).
- `scripts/`: the scripts that were run (stdlib only, run with `python3 -I`).

## How it was obtained (all 2026-10-23, no login)

1. Family members (InterPro 1,929 fungal proteins, paged with the `next` cursor, `scripts/fetch_pan.py`):
   `GET https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/panther/PTHR31018/taxonomy/uniprot/4751?page_size=100`
2. Swiss-Prot and UniProt annotation (cross-check, 1,747 fungal accessions with a PANTHER cross-reference, the UniProt
   release is older than the InterPro one; `scripts/fetch_up.py`):
   `GET https://rest.uniprot.org/uniprotkb/search?query=(taxonomy_id:4751) AND (xref:panther-PTHR31018)&fields=accession,reviewed,organism_id,organism_name,protein_name,gene_names,length,xref_pfam,xref_panther&format=tsv&size=500`
   and the Swiss-Prot query `(taxonomy_id:4751) AND (reviewed:true) AND (gene:ECM33 OR gene:SPS2 OR gene:PST1 OR gene:SPS22 OR protein_name:Ecm33 OR protein_name:"cell wall protein ECM33" OR protein_name:Sps2 OR protein_name:Pst1 OR xref:pfam-PF12454 OR xref:interpro-IPR024238)`.
3. UniRef50 cluster of each accession (`scripts/fetch_idmap.py`): `POST https://rest.uniprot.org/idmapping/run`
   (`from=UniProtKB_AC-ID`, `to=UniRef50`), then `GET .../idmapping/uniref/results/stream/{jobId}?format=json`. It
   returned 1,841 of 1,929 clusters; the remaining 88 were filled one by one (`scripts/fetch_uniref_fill.py`):
   `GET https://rest.uniprot.org/uniref/search?query=(identity:0.5)+AND+(uniprot_id:{accession})&format=tsv&fields=id`.
4. Clade: `GET https://rest.uniprot.org/taxonomy/search?query=tax_id:A OR tax_id:B ...&fields=id,rank,scientific_name,lineage&format=json`
   (`scripts/fetch_tax.py`; 40 taxids per call).
5. Table and counts: `python -I scripts/build_candidates.py (Python 3.10 or newer) WORKDIR members.tsv OUTDIR`.
6. OrthoDB, InterPro entry metadata and PubMed searches that are only described in the report:
   `https://data.orthodb.org/v12/search?query=ecm33&level=4751`, `.../v12/tab?id=1159006at4751`,
   `https://www.ebi.ac.uk/interpro/api/entry/...`, `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=...`.

## Caveats

- `new_cluster` means only that the UniRef50 cluster is not a training member. It is not proof of less than 30 percent
  identity to every training sequence: the freeze's `cluster30` separation has to be computed (MMseqs2, HPCC, not run).
- InterPro and UniProt are live resources: counts change with each release. Re-fetching gives a different table.
- The family call is made by an HMM (PANTHER), so this is independent of the phmmer search from Q5AGC4 but not of
  sequence similarity. T3 rows have no name or functional evidence beyond that call.
