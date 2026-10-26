# Aggregate report with allergen and antigen modules: 10 proteomes (2026-10-26)

Third aggregate report. Adds `allergen_homology` and `antigen_lookup` to the work directories of
`docs/reports/data/aga1_aggregate_domains_2026-10-26/`. Inputs stay outside the repository; no sequences are stored.

Made with `python3 -I analysis/sorting_hat_run/aga1_allergen_antigen_modules.py --runs-root <runs> --taxdump nodes.dmp --taxa analysis/sorting_hat_run/aga1_acceptance_taxa.tsv --cimm-uniprot-fasta <UP000001261.faa>`
(local, no HPCC job), then `cellsurface_sorting_hat_aggregate --manifest manifest.tsv --taxdump nodes.dmp --names names.dmp --rank class --out <new folder>`.

* **Allergen:** BLASTP (local blastp 2.17.0+; `-evalue 1 -max_target_seqs 200 -seg no`, the options of `scripts/sorting_hat/blast_allergen.sbatch`) of each proteome against the IUIS fungal allergen sequences (111 sequences, 5 skipped; database sha256 `86a559cf...a16a7`).
* **Antigen:** module `antigen`, applicable taxon 246410 only. Every other proteome is `not_applicable`. The ranking is `analysis/cocci_antigens/cocci_antigen_ranking.tsv` (sha256 `cf04b9d4...318bc`), keyed on *C. immitis* RS gene IDs. The UniProt proteome of *C. immitis* RS has accessions, so the accession-to-gene map was made from the `GN=CIMG_` field of the UniProt headers: 9,442 of 9,712 accessions map (sha256 of the map `452177f4...e977`). The `protein_map.tsv` that made the earlier Cimm record is on HPCC only and was not used; this is a different identifier route, not a re-run of that record.
* Engine output: `<species>/out_allergen_antigen`; `categories.yaml` sha256 unchanged.

| sha256 | file |
|---|---|
| 29baf95dfcf50db3c75d3ae7db91c4e31330a2dc0e112943958444c52f3e67e2 | counts_by_group.tsv |
| 3b7fa046fe1a79a5c480adb4b16351ad9513dd05f321c91b517bf30f0fe80b06 | counts_by_proteome.tsv |
| e31a36e20cfe4158056ef231809891e6f273755c8c9d75b43ded732c366c17dd | proteomes.tsv |
| 16fec78d7b5f5b56e39709c9b0da91e69cdc921b35441256fd20adcc1f717f6d | report.md |

## Reading the numbers

* `cocci_specificity_rank_top15`: *C. immitis* RS 1,334 of 8,679 assessable proteins (153.7 per 1,000, as expected for a top 15% rule). Every other proteome is `not_applicable`, which the report shows as NA (nothing is assessable there), not as 0. `serodiagnostic_marker_candidate`: 135 in *C. immitis*, none elsewhere.
* `iuis_allergen_similarity` (identity at least 35% over at least 80 residues to an IUIS fungal allergen): 8.8 to 13.9 per 1,000 in every proteome; a similarity call, not an allergen call, and unvalidated.
* `iuis_allergen_homolog` is still 100% not assessable: `pfam_allergen` is `unavailable` because both allergen families are inactive in the family table (owner decision of 2026-10-17, issue #19).
* About 11% of *C. immitis* proteins are not assessable for the antigen rank (accessions with no `GN=CIMG_` or no ranking row).
