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

## H2 results (2026-10-08, `build_truth.py`, UniProt release 2026_03)

- `uniprot_entries.json.gz`: JSON of the 189 query entries (with evidence codes). `truth_all.tsv`: tier per entry.
- Tier counts: T2 131, T3 43, none 15 (the 15 query rows whose name has no hydrophobin or rodlet). No T1 manual entries yet (`manual_entries.tsv` absent).
  Tiers use UniProt evidence codes only (experimental ECO:0000269 from PubMed on a FUNCTION, SUBCELLULAR LOCATION or SUBUNIT comment). The Pfam cross-reference is not used.
  The T2 rule does not check that the experimental statement is hydrophobin-specific. Review before any label is final.
- MMseqs2 17-b804f, 30% identity, coverage 0.5, on the 174 named entries (`clusters_positives.tsv`): **30 clusters**, 28 with a T2 member. Two clusters hold 71 and 40 proteins.
- All 9 entries without a Pfam hydrophobin cross-reference are T2.
- T2 proteins and clusters per species (species with at least 5 T2 proteins): *P. ostreatus* 14 in 2 clusters; *B. bassiana* 9 in 6; *A. fumigatus* 7 in 5; *P. expansum* 7 in 5;
  *E. nidulans* 6 in 5; *V. dahliae* 5 in 2; *G. zeae* 5 in 5; *F. fulva* 5 in 3; *G. moniliformis* 5 in 3.
  Clusters are across all 174 entries, so a species' clusters are shared with other species.
- The development/test split is not made here. It needs clusters over positives **and** assumed negatives per species, so it is made in H2c, before any measurement.

## H2b results (2026-10-08): truth proteins mapped to proteomes

`map_to_proteomes.py` (BLASTP 2.14.0+, at least 95% identity and 90% query coverage, best bitscore, one proteome protein per truth protein), job 29636016.
Outputs: `proteome_map.tsv`, `proteome_map_dropped.tsv` (hits below the thresholds or lost to a conflict), `taxa.tsv` (taxon IDs from `/srv/projects/db/taxonomy` or the run records).

| Proteome | Truth proteins mapped |
|---|---|
| Afum_Af293 / A1163 / W72310 | 7 / 6 / 5 |
| Scer_S288C, Calb_SC5314, Cimm_RS, Bder_ER3, Cneo_H99 | 0 each |
| Bbas_ARSEF2860 | 7 |
| Fvel_6-3 | 5 |
| Fful_Race5 | 4 |
| Fgra_PH-1 | 5 |
| Pexp_MD-8 | 5 |
| Post_PC9 | 4 |
| Tasp_FT101 | 6 |
| Tvir_Gv29-8 | 0 |

The 9 entries without a Pfam hydrophobin cross-reference that map to a measured proteome: HFBE_PENEN (P. expansum), HYD1_GIBZE and HYD2_GIBZE (G. graminearum PH-1),
HCF2_FULFL, HCF1_FULFL and HCF4_FULFL (F. fulva Race5). The other three do not map: HYD1_TRIAP (best hit in T. asperellum FT101 92.6% identity at 100% coverage),
HFB3_HYPVG (best hit in T. virens Gv29-8 67.6%), and PSH_FLAVE (98.4% identity but 85% query coverage in F. velutipes 6-3, below the 90% rule; left unmapped, not relaxed).
So the rescue gain can be measured in three species (P. expansum, G. graminearum, F. fulva) with 1, 2 and 3 entries. Those are small numbers.
Cneo_H99 has no run record and no strain node in names.dmp, so it is not used for status.
