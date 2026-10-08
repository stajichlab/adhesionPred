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

## H2c results (2026-10-08): inputs for `calibrate truth`, split made once

`make_calibration_inputs.py` (job 29636315, MMseqs2 17-b804f, 30% identity, coverage 0.5, seed 20261008, 30% of clusters to development).
Label 1 = mapped T1/T2 truth protein; mapped T3 proteins are excluded; every other proteome protein is label 0 and **assumed** negative (absent annotation).
Every row has a cluster (all proteome proteins were clustered). `split.tsv.gz` holds the dev/test part of every row and is the record made before any measurement.
The per-species `calibration/<proteome>/truth_{all,test,dev}.tsv` tables (13 MB) are not committed. They come from the script and the stored split.

| Proteome | Positives (T1/T2) | Assumed negatives | Clusters | Positives in dev / test |
|---|---|---|---|---|
| Afum_Af293 | 7 | 9640 | 7697 | 1 / 6 |
| Afum_A1163 | 6 | 9936 | 7961 | 2 / 4 |
| Afum_W72310 | 5 | 10551 | 8105 | 0 / 5 |
| Bbas_ARSEF2860 | 7 | 9471 | 7720 | 0 / 7 |
| Fful_Race5 | 4 | 13556 | 10823 | 1 / 3 |
| Fgra_PH-1 | 5 | 11188 | 8750 | 3 / 2 |
| Pexp_MD-8 | 5 | 10619 | 8009 | 2 / 3 |
| Post_PC9 | 3 | 10479 | 7617 | 0 / 3 |

Skipped: Scer_S288C, Calb_SC5314, Cimm_RS, Bder_ER3, Cneo_H99, Tvir_Gv29-8 (0 positives), Fvel_6-3 (0 T1/T2 positives mapped; its mapped proteins are T3), Tasp_FT101 (2 positives, below the minimum of 3).

Reading of the counts: the test part has 2 to 7 positive proteins in each species, and positives per proteome are 3 to 7 in total. The spec rule (at least 10 positive clusters in the test part) is not met in any species.
So no species has a usable test split. The status files will be `smoke`, leakage `partial`, from the whole proteome table, and the report will say so. The development part still serves the rule
that the pattern changes only on development data.
