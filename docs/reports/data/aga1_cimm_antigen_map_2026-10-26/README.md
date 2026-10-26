# C. immitis RS: accession-to-gene map for the antigen ranking (2026-10-26)

**Problem.** The third aggregate report (`docs/reports/data/aga1_aggregate_allergen_antigen_2026-10-26/`) mapped UniProt accessions to ranking genes with the `GN=CIMG_` field of the FASTA header: 9,442 of
9,712 accessions mapped and 1,033 proteins were `not_in_reference` for the antigen rank. A review of that report flagged the 1,033 and asked to reconcile the route with the older HPCC `protein_map.tsv` (RefSeq `XP_` to gene).

**Finding.** The gene symbol in the header (for example `GN=COQ3`) hides the locus name. The UniProt table of the proteome has the locus names (`Gene Names (ORF)`) for all 9,712 accessions:

| route | mapped accessions |
|---|---|
| A. UniProt locus name (`Gene Names (ORF)`, first `CIMG_` name) | 9,712 of 9,712 |
| B. UniProt RefSeq cross-reference through the HPCC `protein_map.tsv` | 9,712 of 9,712 |
| A and B give a different gene | 4 (Q1E9B7, Q1DRC2, I7MN50, Q1DVI3: older locus tags of renamed genes; A is current) |
| header `GN=CIMG_` (the earlier map) | 9,442 |

**Decision recorded here (route A).** The map is route A; the HPCC protein map is not needed. With it the antigen module gives 8,947 `ok` and 765 `not_in_reference` rows (was 8,679 and 1,033). The 765 are genes that are **not in the ranking**
(8,986 genes; 8,822 of them are covered by the proteome), so the remaining gap is a property of the ranking, not of the identifiers. Proteins with a percentile of 15 or less: 1,344 (was 1,334). 108 genes have more than one accession in the proteome (up to 3);
each accession gets the rank of its gene (the module's rule), so a duplicated gene counts once per accession.

**Inputs (outside the repository).** UniProt REST table of `proteome:UP000001261` (`fields=accession,gene_orf,gene_names,xref_refseq`, TSV, 9,712 rows), sha256 `fb97967b18e0730026da814d1c031c05436d59d0eb60ca44bf02b1f709597e13`;
HPCC `protein_map.tsv` sha256 `ca0615211e865559207f98f00e3c676769552331871625ef3dbcea7483c83806`; the map made by route A: sha256 `1dd9d5c7bb9769c3c6c41439b1a584745485322bf8eba86cfec9bb5a41e24332`; ranking sha256 `cf04b9d4d26f25b9c78cd073276aebf4bd19966693232a1e3c5903b5369318bc`.

**Code.** `analysis/sorting_hat_run/aga1_allergen_antigen_modules.py` now takes `--cimm-uniprot-table` (the REST table) instead of the FASTA. The engine outputs of the third report were not re-made; they use the older map.
The *C. immitis* work directory now holds the new `antigen_lookup` module, so the next engine run of that proteome gives 765 and not 1,033 unassessable proteins for `cocci_specificity_rank_top15`.
