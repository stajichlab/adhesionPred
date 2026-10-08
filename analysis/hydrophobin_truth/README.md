# Hydrophobin truth work (2026-10-08)

- `sp_hydrophobin_query.tsv`: UniProt Swiss-Prot (reviewed), Fungi (taxon 4751).
  Query: `(protein_name:hydrophobin OR gene:hyd* OR keyword:KW-0291 OR cc_function:hydrophobin) AND taxonomy_id:4751 AND reviewed:true`
  via https://rest.uniprot.org/uniprotkb/stream, TSV fields accession,id,protein_name,gene_names,organism_name,length,xref_pfam,xref_interpro,ft_signal,cc_function. Fetched 2026-10-08. UniProt release 2026_03 (2026-09-02).
access-control-expose-headers: Link, X-Total-Results, X-UniProt-Release, X-UniProt-Release-Date, X-API-Deployment-Date
x-uniprot-release: 2026_03
x-uniprot-release-date: 02-September-2026.
  The Pfam column is the UniProt cross-reference, not an hmmsearch at Pfam 38.2 GA.
- `nopfam9.faa`, `nopfam9.domtbl`: the 9 named entries with no hydrophobin-class cross-reference; hmmsearch -E 1e-3 against the 18 models in
  `_workdir/sorting_hat/calibration/hydrophobin_discovery/models.hmm`.
- `nopfam9_regex.tsv`: screening regex result per protein (regex in `analysis/sorting_hat_run/hydrophobin_discovery/analyse.py`). Screening only, not the rule.
