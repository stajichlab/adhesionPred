# Kiss 2019 truth set: feasibility report (Task 5.1)

Date: 2026-10-18. Branch: `kiss2019-feasibility`. Plan: `docs/superpowers/plans/2026-10-17-sorting-hat-review-fixes-and-wall-fusion-testing.md`, PR 5, Task 5.1. Spec points: D8, D9 (step 0), D20; finding FD-A8 in `docs/reviews/2026-10-17-decisions-d1-d18-review-fable.md`.

No model was run. No row was added to a truth file. All numbers come from files in `analysis/wall_fusion_truth/results/` (made by `analysis/wall_fusion_truth/map_kiss.py`) or from the queries named below.

**Update 2026-10-22.** The owner downloaded the MycoCosm gene-catalog protein files of 14 of the 16 genomes (they are archived outside git). The rows with a JGI ID were then mapped by exact ID against these files (section 4b). Sections 1, 3 and 4 to 6 keep the numbers of the first, UniProt-only stage as measured; section 4b and the added scenarios in section 6 give the result with the files; section 4c (a second update of the same day) adds a UniProt lookup by the locus tag of the JGI header; sections 8 to 10 are revised.

## 1. Answer in short

- The table has 181 rows. 144 rows have a JGI protein ID (143 distinct IDs). 37 rows have none. The 37 are all family-level rows from Data 4 and Data 8. All 651 rows of Data 3 have a JGI ID.
- **No public route turns a JGI protein ID into a sequence.** The JGI pages are behind a challenge page (HTTP 403). UniProt has no JGI cross-reference. The HPCC proteomes use their own IDs. Section 3 gives the evidence. **The owner's downloaded JGI protein files close this gap (section 4b).**
- I mapped a row only when an identifier in the row itself (a locus tag) matches one UniProt entry in the genome strain. **34 of 144 rows (23.6 percent) map this way.** 26 of the 34 are A. fumigatus. 77 more rows have a unique gene-name hit in UniProt. I keep these as candidates and do not count them as mapped. 33 rows are lost.
- With the UniProt-only strict mapping (34 rows), leave-clade-out was not possible for any slice at any floor. **With the JGI files, 133 of 144 JGI-ID rows are strictly mapped (92.4 percent): 34 by locus tag and 99 by exact ID in the file (74 verified by a name or sequence check of the row, 5 by a UniProt lookup of the header locus tag, 20 through a validated portal). Leave-clade-out becomes possible for adhesion (Saccharomycotina and Sordariomycetes at floor 5, coarse level Saccharomycotina and Pezizomycotina) and for wall (Eurotiomycetes, Saccharomycotina, Sordariomycetes), and stays impossible for fusion (at most 2 clusters per clade).** Section 6 has the tables and the caveat on clusters that also occur in other clades.
- **Owner decisions of 2026-10-22 are final** (section 10): accept all 99 file-mapped rows, keep the name-conflict rule, fine clade level primary, floor 5 for `estimated` with two cluster strata (shared and clade-only), fusion descriptive only, family representatives never test controls.
- The Kiss `cluster_number` is a proxy for the independent unit. It is not a 30 percent identity cluster. The MMseqs2 clustering of D8 is Task 5.3 (HPCC) and is not done here.
- The counts of the spec agree with the data: the paper's 23 families and 45 genes equal Data 3 (adhesion 23 clusters, 45 rows). The 51 rows and 26 clusters of the 181-row table include 6 family-level rows from Data 4 and Data 8 (section 5).

## 2. Organisms and clades (step 1)

Source of each name: Kiss Supplementary Data 1 (MOESM5, "List of 72 genomes", columns Species ID and Species Name). I checked every name in NCBI Taxonomy (E-utilities `esearch` and `efetch`). The clade comes from the NCBI lineage, not from `node_of_emergence`. Coarse clade: phylum, and for Ascomycota the subphylum. Fine clade: as coarse, but class for Pezizomycotina and Basidiomycota.

| Code | Species (Kiss Data 1) | NCBI name | Taxon ID | Lookup scope (UniProt) | Phylum | Subphylum | Class | Coarse clade | Fine clade |
|---|---|---|---|---|---|---|---|---|---|
| Ashgo11 | Eremothecium gossypii | Eremothecium gossypii | 33169 | 284811 Eremothecium gossypii ATCC 10895 | Ascomycota | Saccharomycotina | Saccharomycetes | Saccharomycotina | Saccharomycotina |
| Aspfu1 | Aspergillus fumigatus | Aspergillus fumigatus | 746128 | 330879 Aspergillus fumigatus Af293 | Ascomycota | Pezizomycotina | Eurotiomycetes | Pezizomycotina | Eurotiomycetes |
| Batde5 | Batrachochytrium dendrobatidis | Batrachochytrium dendrobatidis | 109871 | 684364 Batrachochytrium dendrobatidis JAM81 | Chytridiomycota | Chytridiomycotina | Chytridiomycetes | Chytridiomycota | Chytridiomycota |
| Botci1 | Botrytis cinerea | Botrytis cinerea | 40559 | 40559 Botrytis cinerea | Ascomycota | Pezizomycotina | Leotiomycetes | Pezizomycotina | Leotiomycetes |
| Canalb1 | Candida albicans | Candida albicans | 5476 | 237561 Candida albicans SC5314 | Ascomycota | Saccharomycotina | Pichiomycetes | Saccharomycotina | Saccharomycotina |
| Copci | Coprinopsis cinerea | Coprinopsis cinerea | 5346 | 5346 Coprinopsis cinerea | Basidiomycota | Agaricomycotina | Agaricomycetes | Basidiomycota | Agaricomycetes |
| CryneH991 | Cryptococcus neoformans | Cryptococcus neoformans | 5207 | 5207 Cryptococcus neoformans | Basidiomycota | Agaricomycotina | Tremellomycetes | Basidiomycota | Tremellomycetes |
| Crypa2 | Cryphonectria parasitica | Cryphonectria parasitica | 5116 | 660469 Cryphonectria parasitica EP155 | Ascomycota | Pezizomycotina | Sordariomycetes | Pezizomycotina | Sordariomycetes |
| Fusgr1 | Fusarium graminearum | Fusarium graminearum | 5518 | 5518 Fusarium graminearum | Ascomycota | Pezizomycotina | Sordariomycetes | Pezizomycotina | Sordariomycetes |
| Hisca1 | Histoplasma capsulatum | Histoplasma capsulatum | 5037 | 5037 Histoplasma capsulatum | Ascomycota | Pezizomycotina | Eurotiomycetes | Pezizomycotina | Eurotiomycetes |
| Neucr2 | Neurospora crassa | Neurospora crassa | 5141 | 367110 Neurospora crassa OR74A | Ascomycota | Pezizomycotina | Sordariomycetes | Pezizomycotina | Sordariomycetes |
| Sacce1 | Saccharomyces cerevisiae | Saccharomyces cerevisiae | 4932 | 559292 Saccharomyces cerevisiae S288C | Ascomycota | Saccharomycotina | Saccharomycetes | Saccharomycotina | Saccharomycotina |
| Schco3 | Schizophyllum commune | Schizophyllum commune | 5334 | 578458 Schizophyllum commune H4-8 | Basidiomycota | Agaricomycotina | Agaricomycetes | Basidiomycota | Agaricomycetes |
| Schpo1 | Schizosaccharomyces pombe | Schizosaccharomyces pombe | 4896 | 4896 Schizosaccharomyces pombe | Ascomycota | Taphrinomycotina | Schizosaccharomycetes | Taphrinomycotina | Taphrinomycotina |
| Trire2 | Trichoderma reesei | Trichoderma reesei | 51453 | 51453 Trichoderma reesei | Ascomycota | Pezizomycotina | Sordariomycetes | Pezizomycotina | Sordariomycetes |
| Ustma22 | Ustilago maydis | Mycosarcoma maydis (differs) | 5270 | 5270 Mycosarcoma maydis | Basidiomycota | Ustilaginomycotina | Ustilaginomycetes | Basidiomycota | Ustilaginomycetes |

Notes.
- The table has 16 organism codes (the task text counts 17 because the 37 rows without ID have an empty code). 15 of 16 names equal the NCBI scientific name. Ustilago maydis is now Mycosarcoma maydis in NCBI (taxon 5270, same taxon).
- UniProt lookups use the genome strain when Data 1 names it and NCBI Taxonomy has a taxon with exactly that name (8 of 16 organisms use a strain taxon: Ashgo11, Aspfu1, Batde5, Canalb1, Crypa2, Neucr2, Sacce1, Schco3). For the other 8 organisms I used the species taxon. I did not guess a strain.
- `Cryptococcus neoformans` (CryneH991) uses the species taxon 5207, because NCBI has no scientific name equal to "Cryptococcus neoformans var. grubii H99". Locus tags CNJ01930 and CND01240 are JEC21 tags and found no entry (section 4).

## 3. How I tried to get sequences (step 2)

| Order | Route | Result |
|---|---|---|
| a | Data deposit of the paper | The Data availability statement names only Supplementary Data 1-4 and 8. They hold accession numbers of genomes, not protein sequences. Data 1 gives the JGI portal URL and the NCBI assembly ID of each genome. There is no deposit of the sequences. |
| a2 | JGI Mycocosm public pages (`mycocosm.jgi.doe.gov/cgi-bin/dispGeneModel?db=Aspfu1&id=4407` and `/Aspfu1/Aspfu1.home.html`) | HTTP 403, a "Just a moment..." challenge page from the site. I did not try to pass it. JGI downloads also need a login. I did not log in. |
| b | UniProt cross-reference | No entry of the UniProt cross-reference database list (`rest.uniprot.org/configure/uniprotkb/allDatabases`) contains "jgi" (case-insensitive search of the downloaded list). A query `xref:jgi-4407` returned 0 entries. |
| b2 | UniProt by gene name and organism | Used as the second step (sections 4 and 6). It gives gene names and locus tags, not JGI IDs. |
| c | NCBI protein | Exact-string queries `4407_Aspfu1`, `20500_Batde5` and `5432_Sacce1` returned 0 records. The probe was small. I found no route. |
| d | HPCC proteomes (read-only) | See below. No JGI ID in the headers. |

HPCC commands (read-only, nothing written there). Directory: `/bigdata/stajichlab/shared/projects/Fungi_5k`.

| Command | Result |
|---|---|
| `ls input \| grep -c "proteins.fa$"` | 5813 proteomes |
| `head -n1 input/Aspergillus_fumigatus_Af293.proteins.fa` | `>F85C5601_000001-T1 F85C5601_000001` (renamed IDs) |
| `grep -m1 -c -i "jgi" input/Aspergillus_fumigatus_Af293.proteins.fa` | 0 |
| `grep -m1 "AFUA_3G12700\|Afu3g12700" gff3/Aspergillus_fumigatus_Af293.gff3` | no line |
| `head -n 200000 gff3/Aspergillus_fumigatus_Af293.gff3 \| grep -c -i "proteinId\|jgi"` | 0 |
| `ls input \| grep -E "^(...)\.proteins.fa$"` for 10 names | 8 names found (Batde5 JAM81, Botci B05.10, Canalb SC5314, Copci okayama7#130, Cryptococcus JEC21, Schco H4-8, Trire QM6a, Ustma 521). Strains differ from Kiss for Copci (AmutBmut), Cryptococcus (H99) and Botci. No S288C or Hisca NAm1 proteome under that name. Not used. |
| `ls .../adhesionPred/_workdir/step1`, `.../fungi_scan` | no Kiss proteome. Not used. |

Headers of the other HPCC proteomes (Neurospora crassa OR74A, Candida albicans SC5314, Ustilago maydis 521) have the same renamed form `<prefix>_NNNNNN-T1`. A JGI ID cannot be matched to them.

Consequence at that stage: every row with a JGI ID is counted as having **no ID-verified sequence**. I did not match an ID to a sequence by guess. (The owner later supplied the files; section 4b.)

## 4. What I did instead

1. **Locus tag in the row.** The field `gene_or_family_id` of a row can hold a locus tag: `Afu3g12700` (turned into `AFUA_3G12700` by one fixed rule), `NCU09332`, `AGOS_AAL050W`, `BATDEDRAFT_20500`, and the forms `CNJ01930`, `CND01240`, `UM05433`, `ACR290W`. I queried UniProt with `(gene_exact:<tag>) AND (taxonomy_id:<scope>)`. A row maps when exactly one entry (reviewed first) is found. For the 2 Batde5 rows the number in the tag equals the JGI ID (`20500`, `35524`): this is a second check, and it agrees. Result: 34 rows.
2. **Family rows without a JGI ID.** Data 4 and Data 8 rows give an A. fumigatus ortholog locus tag. 35 of the 37 map to the A. fumigatus protein. This protein stands for the family. It is not the Kiss protein of another species. These rows add Eurotiomycetes only. They are reported apart (column "Family rep").
3. **Gene name.** A gene name that gives exactly one UniProt entry in the strain (or species) is a candidate: 77 rows. I do not count them as mapped. A name from a paper is not the JGI protein model, and a name can hold two genes.
4. **Plausibility check.** A sequence must have 20 to 15000 residues, only amino-acid letters and no run of 5 or more X. All 34 + 35 accepted sequences passed. `results/mapped_sequences.faa` holds the 52 distinct UniProt sequences (UniProt sources, not JGI). Its SHA-256 is in the section "Files" below.

### Counts per slice

Slices follow D20 and are never pooled: adhesion (`adhesion`), wall (`cell wall`, `cell wall (GPI anchor assembly)`), fusion (`cell fusion*`, `pheromone receptor`).

| Slice | Rows | Kiss clusters | With JGI ID | Without JGI ID | Mapped by locus tag | Family rep (no ID) | Name candidate | Lost |
|---|---|---|---|---|---|---|---|---|
| adhesion | 51 | 26 | 45 | 6 | 4 | 4 | 27 | 16 |
| wall | 122 | 44 | 94 | 28 | 29 | 28 | 46 | 19 |
| fusion | 8 | 6 | 5 | 3 | 1 | 3 | 4 | 0 |
| all | 181 | 76 | 144 | 37 | 34 | 35 | 77 | 35 |

Rows lost, with reason (all slices):

| Slice | Reason | Rows |
|---|---|---|
| adhesion | jgi_id_has_no_public_sequence_route; gene_name_ambiguous | 2 |
| adhesion | jgi_id_has_no_public_sequence_route; gene_name_no_hit | 12 |
| adhesion | no_jgi_id_and_no_afu_locus | 2 |
| wall | jgi_id_has_no_public_sequence_route; gene_name_no_hit | 15 |
| wall | jgi_id_has_no_public_sequence_route; locus_tag_no_hit_in_genome_strain | 4 |

The 33 lost rows with a JGI ID: 27 have a gene name with no hit (mostly N. crassa names such as chs1 and gel4; also 8 Basidiomycota rows with Schco3 and CryneH991 names), 2 have a name with two hits (als3_als1, pga62_flo1), 4 carry a locus tag of another strain or an old form (ACR290W, CND01240, CNJ01930, UM05433). The 2 lost rows without JGI ID have no A. fumigatus locus (`NA`, and one C. albicans ORF).

### Rows per organism code

| Organism code | Clade | Slice | Rows | Locus tag | Family rep | Name candidate | Lost |
|---|---|---|---|---|---|---|---|
| (no JGI ID) |  | adhesion | 6 | 0 | 4 | 0 | 2 |
| (no JGI ID) |  | wall | 28 | 0 | 28 | 0 | 0 |
| (no JGI ID) |  | fusion | 3 | 0 | 3 | 0 | 0 |
| Ashgo11 | Saccharomycotina | wall | 4 | 3 | 0 | 0 | 1 |
| Aspfu1 | Pezizomycotina | adhesion | 3 | 3 | 0 | 0 | 0 |
| Aspfu1 | Pezizomycotina | wall | 22 | 22 | 0 | 0 | 0 |
| Aspfu1 | Pezizomycotina | fusion | 1 | 1 | 0 | 0 | 0 |
| Batde5 | Chytridiomycota | wall | 2 | 2 | 0 | 0 | 0 |
| Botci1 | Pezizomycotina | adhesion | 2 | 0 | 0 | 2 | 0 |
| Canalb1 | Saccharomycotina | adhesion | 11 | 0 | 0 | 9 | 2 |
| Canalb1 | Saccharomycotina | wall | 25 | 0 | 0 | 25 | 0 |
| Copci | Basidiomycota | adhesion | 2 | 0 | 0 | 2 | 0 |
| CryneH991 | Basidiomycota | adhesion | 1 | 0 | 0 | 0 | 1 |
| CryneH991 | Basidiomycota | wall | 2 | 0 | 0 | 0 | 2 |
| Crypa2 | Pezizomycotina | adhesion | 1 | 0 | 0 | 0 | 1 |
| Fusgr1 | Pezizomycotina | adhesion | 1 | 0 | 0 | 1 | 0 |
| Hisca1 | Pezizomycotina | adhesion | 2 | 0 | 0 | 0 | 2 |
| Hisca1 | Pezizomycotina | wall | 2 | 0 | 0 | 0 | 2 |
| Neucr2 | Pezizomycotina | adhesion | 6 | 1 | 0 | 1 | 4 |
| Neucr2 | Pezizomycotina | wall | 19 | 2 | 0 | 4 | 13 |
| Neucr2 | Pezizomycotina | fusion | 1 | 0 | 0 | 1 | 0 |
| Sacce1 | Saccharomycotina | adhesion | 8 | 0 | 0 | 8 | 0 |
| Sacce1 | Saccharomycotina | wall | 14 | 0 | 0 | 14 | 0 |
| Sacce1 | Saccharomycotina | fusion | 2 | 0 | 0 | 2 | 0 |
| Schco3 | Basidiomycota | adhesion | 3 | 0 | 0 | 0 | 3 |
| Schpo1 | Taphrinomycotina | wall | 3 | 0 | 0 | 3 | 0 |
| Trire2 | Pezizomycotina | adhesion | 2 | 0 | 0 | 2 | 0 |
| Ustma22 | Basidiomycota | adhesion | 3 | 0 | 0 | 2 | 1 |
| Ustma22 | Basidiomycota | wall | 1 | 0 | 0 | 0 | 1 |
| Ustma22 | Basidiomycota | fusion | 1 | 0 | 0 | 1 | 0 |

By clade (144 rows with a JGI ID): Saccharomycotina 64 (Canalb1 36, Sacce1 24, Ashgo11 4), Pezizomycotina 62, Basidiomycota 13 (Ustma22 5, Schco3 3, CryneH991 3, Copci 2), Taphrinomycotina 3, Chytridiomycota 2. None of the 13 Basidiomycota rows maps by locus tag. 5 are candidates and 8 are lost.

## 4b. JGI-file step (2026-10-22)

**Files.** The owner downloaded the gene-catalog protein FASTA file of 14 portals from MycoCosm (not Aspfu1 and not Batde5). `results/jgi_files.tsv` lists all 23 downloaded files with sha256, md5, protein count and the check IDs (no sequence, no local path). The md5 of every file equals the md5 in the JGI download manifests (`Download_2426889_File_Manifest.csv`, `Download_2426926_File_Manifest.csv`). The script checks the sha256 of each used file against this table before it reads the file, and it refuses to run on a mismatch. The files stay in `to_import/kiss2019_jgi/` (outside git). Portal names that differ from the Kiss codes: Ustma22 is Ustma2_2, CryneH991 is Cryne_H99_1, Ashgo11 is Ashgo1_1, Botci1 is Bcine1, Copci is Copci_AmutBmut1. Canalb1 has two versions; only `20161025` contains the Kiss IDs (`20160604` has another ID space). The 28 rows of Aspfu1 (26) and Batde5 (2) have no file; all 28 are mapped by locus tag already.

**Rule.** A row is looked up by its exact protein ID in `jgi|<portal>|<id>|<name>`. The found protein is checked against the row, in this order (column `jgi_check_how`):

1. a locus tag of the row, if the header uses the same naming scheme (`NCU08192` against `NCU08192T0`; `UM05433` against `UMAG_05433` by number); another scheme (`CND01240` against `CNAG_00996`) is not comparable;
2. a gene name of the row equal to the header name (Schpo1: `bgs1`);
3. the UniProt candidate of the row: its gene names carry the header locus, or its sequence is identical. A similar sequence that is not identical (difflib ratio 0.90 or more) is not evidence. A different sequence counts as a mismatch.

`agree` maps the row (status `mapped_jgi_file`, `jgi_accept = row_verified`). A mismatch does not map it. A row with no evidence (`unverifiable`) maps only when its whole portal is validated: at least 3 rows of the portal agree and at most 10 percent of the rows with evidence disagree (`jgi_accept = portal_validated`). Rows mapped by locus tag keep their status. The sequence must also pass the plausibility check of section 4. The tables hold the sequence length and the sha256 of the sequence, never the sequence (the script scans every written table for the found sequences and stops if one is there; the committed UniProt cache and `mapped_sequences.faa` hold public UniProt sequences and are not scanned).

**Result of the file step for the 110 rows that the UniProt stage did not map by locus tag** (before the locus-tag lookup of section 4c; `results/jgi_portal_summary.tsv`, `results/row_mapping.tsv`):

| Portal | Rows with ID | Found in file | Check: agree | differ | no evidence | Not mapped by locus tag | Mapped, row verified | Mapped, portal validated | Not mapped |
|---|---|---|---|---|---|---|---|---|---|
| Canalb1 | 36 | 36 | 33 | 1 | 2 | 36 | 33 | 2 | 1 |
| Sacce1 | 24 | 24 | 24 | 0 | 0 | 24 | 24 | 0 | 0 |
| Neucr2 | 26 | 26 | 9 | 0 | 17 | 23 | 6 | 17 | 0 |
| Ustma2_2 | 5 | 5 | 4 | 0 | 1 | 5 | 4 | 1 | 0 |
| Hisca1 | 4 | 4 | 0 | 0 | 4 | 4 | 0 | 0 | 4 |
| Schco3 | 3 | 3 | 0 | 0 | 3 | 3 | 0 | 0 | 3 |
| Cryne_H99_1 | 3 | 3 | 0 | 0 | 3 | 3 | 0 | 0 | 3 |
| Schpo1 | 3 | 3 | 3 | 0 | 0 | 3 | 3 | 0 | 0 |
| Copci_AmutBmut1 | 2 | 2 | 0 | 1 | 1 | 2 | 0 | 0 | 2 |
| Bcine1 | 2 | 2 | 0 | 2 | 0 | 2 | 0 | 0 | 2 |
| Trire2 | 2 | 2 | 2 | 0 | 0 | 2 | 2 | 0 | 0 |
| Crypa2 | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 0 | 1 |
| Fusgr1 | 1 | 1 | 1 | 0 | 0 | 1 | 1 | 0 | 0 |
| Ashgo1_1 | 4 | 4 | 4 | 0 | 0 | 1 | 1 | 0 | 0 |
| Aspfu1, Batde5 | 28 | 0 | n/a | n/a | n/a | 0 | 0 | 0 | 0 |
| all | 144 | 116 | 80 | 4 | 32 | 110 | 74 | 20 | 16 |

All 116 IDs were found in the files. Of the 110 rows, **94 were mapped by file at this step and 16 were not**. The 94 rows hold 93 distinct sequences (cot-2 and gs-1 are two Kiss rows of one protein, 8184_Neucr2). By slice (status after the file step; section 4c has the final table):

| Slice | Rows | Kiss clusters | With JGI ID | Without JGI ID | Mapped by locus tag | Mapped by file | Family rep (no ID) | Name candidate | Lost |
|---|---|---|---|---|---|---|---|---|---|
| adhesion | 51 | 26 | 45 | 6 | 4 | 29 | 4 | 1 | 13 |
| wall | 122 | 44 | 94 | 28 | 29 | 61 | 28 | 0 | 4 |
| fusion | 8 | 6 | 5 | 3 | 1 | 4 | 3 | 0 | 0 |
| all | 181 | 76 | 144 | 37 | 34 | 94 | 35 | 1 | 17 |

Rows with a JGI ID that were strictly mapped after this step (locus tag or file): adhesion 33 of 45, wall 90 of 94, fusion 5 of 5; together 128 of 144. The 17 lost rows are the 2 rows without ID and A. fumigatus locus (section 4) and 15 rows with an ID; the 1 candidate is coh2.

**The gene-name candidates of section 4 are validated by the files.** Of the 77 candidates, 71 have a sequence identical to the JGI protein, 1 (ste2, Sacce1) is similar and agrees by locus, 1 (coh2) is similar and has no other evidence, and 4 differ (coh1, hyr1, bhp2, bhp3). So 72 of 77 (93.5 percent) are confirmed and 71 of 77 (92.2 percent) are identical. This is a measured error rate for the cheapest route (a gene name searched in UniProt inside the species or strain taxon).

**The 16 rows that were not mapped at this step:**

| Rows | Portal | Why |
|---|---|---|
| coh1 | Copci_AmutBmut1 | the JGI protein (120 aa) differs from the UniProt candidate (113 aa, ratio 0.58); the candidate may come from another strain (UniProt has Okayama 7 for this species) or another gene model |
| coh2 | Copci_AmutBmut1 | similar to the candidate (110 aa, ratio at least 0.90) but not identical; no other evidence; stays a candidate |
| hyr1 | Canalb1 | JGI protein 1249 aa, candidate 919 aa; the same portal agrees for 33 other rows, so this is a model or name difference, not a version problem |
| bhp2, bhp3 | Bcine1 | JGI proteins 691 and 151 aa; candidates 98 aa. Bcine1 is the 2025 annotation of B05.10; Kiss used Botci1 v1.0. 0 of 2 rows agree, so this file does not share the Kiss ID space |
| bad1, hyd1_like, pmt1, ags1 | Hisca1 | no name evidence (names do not occur in the header; no UniProt candidate); the portal has no agreeing row, so it is not validated |
| sc1, sc3, sc6 | Schco3 | same |
| cfl1, CND01240, CNJ01930 | Cryne_H99_1 | same (CND and CNJ are JEC21 tags; the header uses CNAG tags) |
| cryp | Crypa2 | same |

Section 4c runs the UniProt query by the header locus tags for these rows.

## 4c. UniProt lookup by the locus tag of the JGI header (2026-10-22, second update)

**Rule (owner approved).** For each row that is still unmapped after section 4b, has a found JGI protein and a locus tag in the header name (`HCAG_00594T0` gives `HCAG_00594`, `BCIN_02g00390.1` gives `BCIN_02g00390`; gene-model names such as `fgenesh1_kg.2_#_1532_...` have none), I queried UniProt with `(gene_exact:<tag>) AND (taxonomy_id:<scope>)` (public API, `--online`, answers cached in `results/query_cache.json`). Scope: the genome strain of Kiss Data 1 when NCBI Taxonomy resolves that term to one taxon, else the scope of section 2. This changed one scope: "Histoplasma capsulatum NAm1" resolves to taxon 2059318, now named *Histoplasma mississippiense*. With the species taxon 5037 all four Hisca1 tags gave no hit; with 2059318 three give a hit. One exact entry (reviewed first) is required, and its sequence must equal the JGI protein. Then the row is `mapped_jgi_file` with `jgi_accept = locus_confirmed`. A hit with another sequence is reported as disagree. Rows whose ID is absent from the file would stay unmapped (none: all 116 IDs were found).

One more rule that I added: when the row has a gene-name candidate (section 4) and the locus-tag entry is a different UniProt entry, the row's name contradicts the protein at this ID. The row is not mapped (`confirmed_name_conflict`).

**What a recovery proves.** The locus-tag hit shows that the file protein equals that UniProt entry. It does not show that the Kiss row means this protein: for the 5 recovered rows the UniProt gene names are only the locus tag, so there is no name evidence. `locus_confirmed` is therefore weaker than `row_verified`, and the row-verified scenario does not include it.

Ten rows have a header locus tag; six have none (coh1, coh2, sc1, sc3, sc6, cryp). Outcome of the ten:

| Row | Portal | Header locus tag | UniProt result | Outcome |
|---|---|---|---|---|
| bad1 | Hisca1 | HCAG_00594 | A0ACE4JU10, 137 aa, identical | recovered |
| hyd1_like | Hisca1 | HCAG_06031 | A0ACE4KFX7, 98 aa, identical | recovered |
| pmt1 | Hisca1 | HCAG_02731 | A0ACE4K441, 776 aa, identical | recovered |
| ags1 | Hisca1 | HCAG_11199 | no hit | stays lost |
| cfl1 | Cryne_H99_1 | CNAG_00795 | J9VHS9, 309 aa, identical | recovered |
| CND01240 | Cryne_H99_1 | CNAG_00996 | J9VNY5, 767 aa, identical | recovered |
| CNJ01930 | Cryne_H99_1 | CNAG_04763 | J9VXX5, 807 aa; JGI protein 810 aa | disagree, stays lost |
| hyr1 | Canalb1 | CAALFM_C500730WA | Q59XA7 (reviewed, gene names HYR3 IFF2), 1249 aa, identical | name conflict, not mapped |
| bhp2 | Bcine1 | BCIN_02g00390 | A0A384J844, 691 aa, identical | name conflict, not mapped |
| bhp3 | Bcine1 | BCIN_03g05320 | A0A384JCJ6, 151 aa, identical | name conflict, not mapped |

**Bcine1 (the re-annotation issue).** `Bcine1` is the 2025 annotation of B05.10 (11,409 proteins in `FilteredModels2`); Kiss used Botci1 v1.0. The two locus tags resolve in UniProt and give exactly the JGI sequences, so the file agrees with UniProt. But the Kiss names bhp2 and bhp3 point to hydrophobin candidates of 98 aa (A0A384JW35, A0A384JJM8), while the proteins at IDs 6405 and 8060 are 691 and 151 aa. The IDs were renumbered between the two annotations. So the locus tags **do not help** for Bcine1: they show that the file is consistent with UniProt, and that its IDs are not the Kiss IDs. Both rows stay unmapped and the portal cannot be used for Kiss IDs. A file of the original Botci1 v1.0 annotation would be needed. The same holds for hyr1: the protein at ID 56210 is Hyr3/Iff2 in UniProt, not Hyr1; the other 33 Canalb1 rows with evidence agree, so I read this as an ID or label difference of one row.

**Result.** 5 rows are recovered: Hisca1 3 (adhesion 2, wall 1) and Cryne_H99_1 2 (adhesion 1, wall 1). Final status of the 144 rows with a JGI ID (`results/counts_per_slice.tsv`):

| Slice | Rows | Kiss clusters | With JGI ID | Without JGI ID | Mapped by locus tag | Mapped by file | Family rep (no ID) | Name candidate | Lost |
|---|---|---|---|---|---|---|---|---|---|
| adhesion | 51 | 26 | 45 | 6 | 4 | 32 | 4 | 1 | 10 |
| wall | 122 | 44 | 94 | 28 | 29 | 63 | 28 | 0 | 2 |
| fusion | 8 | 6 | 5 | 3 | 1 | 4 | 3 | 0 | 0 |
| all | 181 | 76 | 144 | 37 | 34 | 99 | 35 | 1 | 12 |

The 99 rows mapped by file: 74 row verified, 5 locus confirmed, 20 portal validated. Strictly mapped rows with a JGI ID (locus tag or file): adhesion 36 of 45, wall 92 of 94, fusion 5 of 5; **together 133 of 144 (92.4 percent)**. Row-verified rows: 34 + 74 = 108 (unchanged).

**The 11 rows that stay unmapped:**

| Rows | Portal | Why |
|---|---|---|
| coh1, coh2 | Copci_AmutBmut1 | no locus tag in the header (gene-model names); coh1 differs from its candidate, coh2 is similar but not identical (stays a candidate) |
| sc1, sc3, sc6 | Schco3 | gene-model names, no locus tag, no evidence |
| cryp | Crypa2 | same |
| hyr1 | Canalb1 | locus-tag entry is Hyr3/Iff2: name conflict |
| bhp2, bhp3 | Bcine1 | locus-tag entries are other proteins than the named hydrophobins: IDs renumbered |
| CNJ01930 | Cryne_H99_1 | UniProt entry 3 aa shorter than the JGI protein: disagree |
| ags1 | Hisca1 | no UniProt entry for HCAG_11199 |

## 5. Reconciliation with report 06 and the paper (step 4)

Source of the full numbers: Data 3 (`S7.xlsx`, present, sheet `gene list`), re-derived by `map_kiss.py --data3-xlsx` (`results/data3_full_counts.tsv`). All 144 Data 3 rows of the 181-row table were found in the xlsx by (gene, JGI ID, cluster).

| Slice | Rows (181-row table) | Clusters (181-row table) | Rows (Data 3 only) | Clusters (Data 3 only) | Paper or report 06 |
|---|---|---|---|---|---|
| adhesion | 51 | 26 | 45 | 23 | "23 families" and "45 genes": **agree with Data 3 only** |
| wall | 122 | 44 | 94 | 40 | report 06: cell wall biogenesis 94 rows, 40 clusters: **agrees with Data 3 only** |
| fusion | 8 | 6 | 5 | 3 | no fusion set in the paper (zero hits for "fusion" in the full text, report 06 row 22) |

- Adhesion: Data 3 has 45 rows in 23 clusters. This equals the paper (45 genes, 23 families). The 181-row table adds 6 family-level rows: 3 in clusters that Data 3 already has (417, 4291, 7971) and 3 in new clusters (338, 3660, 4240). That gives 51 rows and 26 clusters. Both pairs of numbers are right. They count different sets.
- Wall: Data 3 gives 94 rows and 40 clusters. The table adds 28 family-level rows and 4 new clusters (109, 1348, 179, 2065).
- Fusion: 5 rows of Data 3 (so_soft, ste3, ste2, pra1, preA) are picked by keyword, not by a Kiss category. The other 3 come from Data 8 (clusters 1764, 2511, 3355). The fusion slice is a hand-picked list. Kiss did not define it.
- Data 3 as a whole has 651 rows, 362 clusters, and a JGI ID in every row (`data3_full_counts.tsv`).
- The 10 clusters that are new in the family-level rows (adhesion 3, wall 4, fusion 3) exist only as A. fumigatus representatives. They have no JGI ID and no other organism.

## 6. Floors for leave-clade-out (steps 3 and 5)

Independent unit: the Kiss `cluster_number` (a proxy). A clade is eligible for hold-out at floor k when it has at least k distinct clusters among the rows that pass. Three scenarios:

- **upper bound**: all 144 rows with a JGI ID, as if the owner gets the JGI protein files.
- **strict**: rows mapped by locus tag (34 rows). This is the only scenario that stands on identifiers now.
- **strict plus candidates**: 34 plus the 77 gene-name candidates (not accepted; shown to size the gain).
- **strict plus JGI file** (added 2026-10-22): the 34 locus-tag rows plus the 99 rows mapped by file (133 rows). This is the scenario that stands on identifiers now.
- **strict plus JGI file, row verified**: as before, but without the 20 portal-validated and the 5 locus-confirmed rows (108 rows).

The first three scenarios are the UniProt-only stage of the first run and keep their numbers. The tables below give, for the two added scenarios, the numbers of eligible clades.

Clades: "coarse" and "fine" as in section 2. Cells give the number of eligible clades. The names are in `results/floors_per_slice.tsv`.

Upper bound:

| Slice | Level | Clades present | Floor 3 | Floor 5 | Floor 8 |
|---|---|---|---|---|---|
| adhesion | coarse | 3 | 3 | 2 | 2 |
| adhesion | fine | 7 | 4 | 2 | 1 |
| wall | coarse | 5 | 2 | 2 | 2 |
| wall | fine | 7 | 3 | 3 | 3 |
| fusion | coarse | 3 | 0 | 0 | 0 |
| fusion | fine | 4 | 0 | 0 | 0 |

Strict:

| Slice | Level | Clades present | Floor 3 | Floor 5 | Floor 8 |
|---|---|---|---|---|---|
| adhesion | coarse | 1 | 1 | 0 | 0 |
| adhesion | fine | 2 | 0 | 0 | 0 |
| wall | coarse | 3 | 1 | 1 | 1 |
| wall | fine | 4 | 1 | 1 | 1 |
| fusion | coarse | 1 | 0 | 0 | 0 |
| fusion | fine | 1 | 0 | 0 | 0 |

Strict plus candidates:

| Slice | Level | Clades present | Floor 3 | Floor 5 | Floor 8 |
|---|---|---|---|---|---|
| adhesion | coarse | 3 | 2 | 1 | 1 |
| adhesion | fine | 6 | 1 | 1 | 1 |
| wall | coarse | 4 | 2 | 2 | 2 |
| wall | fine | 5 | 3 | 2 | 2 |
| fusion | coarse | 3 | 0 | 0 | 0 |
| fusion | fine | 4 | 0 | 0 | 0 |

Strict plus JGI file (133 rows):

| Slice | Level | Clades present | Floor 3 | Floor 5 | Floor 8 |
|---|---|---|---|---|---|
| adhesion | coarse | 3 | 3 | 2 | 2 |
| adhesion | fine | 5 | 4 | 2 | 1 |
| wall | coarse | 5 | 2 | 2 | 2 |
| wall | fine | 7 | 3 | 3 | 3 |
| fusion | coarse | 3 | 0 | 0 | 0 |
| fusion | fine | 4 | 0 | 0 | 0 |

Strict plus JGI file, row verified (108 rows):

| Slice | Level | Clades present | Floor 3 | Floor 5 | Floor 8 |
|---|---|---|---|---|---|
| adhesion | coarse | 3 | 2 | 1 | 1 |
| adhesion | fine | 4 | 1 | 1 | 1 |
| wall | coarse | 5 | 2 | 2 | 2 |
| wall | fine | 6 | 3 | 2 | 2 |
| fusion | coarse | 3 | 0 | 0 | 0 |
| fusion | fine | 4 | 0 | 0 | 0 |

Eligible clades, strict plus JGI file (clusters, rows; floors 5 and 8 unless noted):

- adhesion, coarse: Pezizomycotina (9 clusters, 14 rows) and Saccharomycotina (10, 18) at floors 3, 5 and 8; Basidiomycota (4, 4) at floor 3 only.
- adhesion, fine: Saccharomycotina (10, 18) at all floors; Sordariomycetes (6, 9) at floors 3 and 5; Eurotiomycetes (3, 5) and Ustilaginomycetes (3, 3) at floor 3 only.
- wall, coarse: Pezizomycotina (23, 42) and Saccharomycotina (29, 43) at all floors.
- wall, fine: Eurotiomycetes (16, 23), Saccharomycotina (29, 43) and Sordariomycetes (14, 19) at all floors.
- fusion: no clade, at any floor or level (at most 2 clusters, 1 to 2 rows per clade).
- Row verified only: adhesion coarse floor 5 and 8 keep only Saccharomycotina (9 clusters, 16 rows); adhesion fine keeps only Saccharomycotina; wall fine Sordariomycetes drops to 4 clusters (6 rows), so floors 5 and 8 keep Eurotiomycetes and Saccharomycotina.

**Clusters that also occur in another clade.** The Kiss cluster is a family. When a held-out clade has a cluster that also has members in a training clade, the held-out proteins have homologs in training, and the hold-out is weaker. In the strict plus JGI file scenario (`results/clusters_per_clade.tsv`, column "also in other clade"):

| Slice | Level | Clade | Clusters | Also in another clade | Only in this clade |
|---|---|---|---|---|---|
| adhesion | coarse | Pezizomycotina | 9 | 0 | 9 |
| adhesion | coarse | Saccharomycotina | 10 | 0 | 10 |
| adhesion | fine | Saccharomycotina | 10 | 0 | 10 |
| adhesion | fine | Sordariomycetes | 6 | 0 | 6 |
| wall | coarse | Pezizomycotina | 23 | 13 | 10 |
| wall | coarse | Saccharomycotina | 29 | 13 | 16 |
| wall | fine | Eurotiomycetes | 16 | 12 | 4 |
| wall | fine | Saccharomycotina | 29 | 13 | 16 |
| wall | fine | Sordariomycetes | 14 | 10 | 4 |

For adhesion no cluster crosses a clade, so the hold-out is clean. For wall, 13 of the 29 Saccharomycotina clusters and 12 of the 16 Eurotiomycetes clusters also have members in another clade. If the unit is the cluster that occurs only in the held-out clade, then at floor 5 wall keeps Saccharomycotina (16) and, at the coarse level, Pezizomycotina (10), and drops Eurotiomycetes and Sordariomycetes (4 each). I give this as a size, not as a decision on the rule.

Clusters and rows per clade of the three UniProt-stage scenarios (fine level, plus the coarse clades Basidiomycota and Pezizomycotina where they exist; the per-clade rows of the two added scenarios are in `results/clusters_per_clade.tsv`):

Upper bound:

| Slice | Level | Clade | Rows | Clusters | Clusters also in another clade |
|---|---|---|---|---|---|
| adhesion | coarse | Basidiomycota | 9 | 4 | 0 |
| adhesion | coarse | Pezizomycotina | 17 | 9 | 0 |
| adhesion | fine | Agaricomycetes | 5 | 1 | 1 |
| adhesion | fine | Eurotiomycetes | 5 | 3 | 0 |
| adhesion | fine | Leotiomycetes | 2 | 1 | 1 |
| adhesion | fine | Saccharomycotina | 19 | 10 | 0 |
| adhesion | fine | Sordariomycetes | 10 | 6 | 1 |
| adhesion | fine | Tremellomycetes | 1 | 1 | 0 |
| adhesion | fine | Ustilaginomycetes | 3 | 3 | 1 |
| fusion | coarse | Basidiomycota | 1 | 1 | 1 |
| fusion | coarse | Pezizomycotina | 2 | 2 | 1 |
| fusion | fine | Eurotiomycetes | 1 | 1 | 1 |
| fusion | fine | Saccharomycotina | 2 | 2 | 1 |
| fusion | fine | Sordariomycetes | 1 | 1 | 0 |
| fusion | fine | Ustilaginomycetes | 1 | 1 | 1 |
| wall | coarse | Basidiomycota | 3 | 1 | 1 |
| wall | coarse | Pezizomycotina | 43 | 23 | 13 |
| wall | fine | Chytridiomycota | 2 | 1 | 1 |
| wall | fine | Eurotiomycetes | 24 | 16 | 12 |
| wall | fine | Saccharomycotina | 43 | 29 | 13 |
| wall | fine | Sordariomycetes | 19 | 14 | 10 |
| wall | fine | Taphrinomycotina | 3 | 2 | 1 |
| wall | fine | Tremellomycetes | 2 | 1 | 1 |
| wall | fine | Ustilaginomycetes | 1 | 1 | 1 |

Strict:

| Slice | Level | Clade | Rows | Clusters | Clusters also in another clade |
|---|---|---|---|---|---|
| adhesion | coarse | Pezizomycotina | 4 | 3 | 0 |
| adhesion | fine | Eurotiomycetes | 3 | 2 | 0 |
| adhesion | fine | Sordariomycetes | 1 | 1 | 0 |
| fusion | coarse | Pezizomycotina | 1 | 1 | 0 |
| fusion | fine | Eurotiomycetes | 1 | 1 | 0 |
| wall | coarse | Pezizomycotina | 24 | 16 | 1 |
| wall | fine | Chytridiomycota | 2 | 1 | 1 |
| wall | fine | Eurotiomycetes | 22 | 16 | 1 |
| wall | fine | Saccharomycotina | 3 | 1 | 1 |
| wall | fine | Sordariomycetes | 2 | 1 | 1 |

Strict plus candidates:

| Slice | Level | Clade | Rows | Clusters | Clusters also in another clade |
|---|---|---|---|---|---|
| adhesion | coarse | Basidiomycota | 4 | 2 | 0 |
| adhesion | coarse | Pezizomycotina | 10 | 4 | 0 |
| adhesion | fine | Agaricomycetes | 2 | 1 | 1 |
| adhesion | fine | Eurotiomycetes | 3 | 2 | 0 |
| adhesion | fine | Leotiomycetes | 2 | 1 | 1 |
| adhesion | fine | Saccharomycotina | 17 | 9 | 0 |
| adhesion | fine | Sordariomycetes | 5 | 2 | 1 |
| adhesion | fine | Ustilaginomycetes | 2 | 2 | 1 |
| fusion | coarse | Basidiomycota | 1 | 1 | 1 |
| fusion | coarse | Pezizomycotina | 2 | 2 | 1 |
| fusion | fine | Eurotiomycetes | 1 | 1 | 1 |
| fusion | fine | Saccharomycotina | 2 | 2 | 1 |
| fusion | fine | Sordariomycetes | 1 | 1 | 0 |
| fusion | fine | Ustilaginomycetes | 1 | 1 | 1 |
| wall | coarse | Pezizomycotina | 28 | 19 | 13 |
| wall | fine | Chytridiomycota | 2 | 1 | 1 |
| wall | fine | Eurotiomycetes | 22 | 16 | 10 |
| wall | fine | Saccharomycotina | 42 | 29 | 13 |
| wall | fine | Sordariomycetes | 6 | 4 | 4 |
| wall | fine | Taphrinomycotina | 3 | 2 | 1 |

What the tables say.

- **Fusion**: at most 2 clusters in any clade, in every scenario. No clade passes even floor 3. Fusion has 8 rows and 6 clusters (3 clusters in Data 3). It can only be descriptive recall. Leave-clade-out and leave-cluster-out are both too thin for a statement.
- **Adhesion**: upper bound, floor 5: Saccharomycotina (10 clusters, 19 rows) and Sordariomycetes (6 clusters, 10 rows). At floor 8 only Saccharomycotina. With strict mapping, the 4 rows are A. fumigatus and N. crassa: no clade passes floor 5. With candidates, only Saccharomycotina passes floor 5. In the upper bound, the 10 Saccharomycotina clusters have no member in another clade. A held-out clade therefore removes whole clusters, and nothing of the family stays in the training clades.
- **Wall**: upper bound, floors 3, 5 and 8: Eurotiomycetes (16 clusters), Saccharomycotina (29) and Sordariomycetes (14). Strict mapping: only Eurotiomycetes (16 clusters, 22 of the 24 rows are the A. fumigatus proteins). The other clades keep 1 cluster. Hold-out of Eurotiomycetes would leave 7 strict rows (3 clades, 1 cluster each) for training.
- **Basidiomycota** and the non-hyphal clades are under every floor. Basidiomycota has 4 clusters in adhesion (upper bound, 9 rows), 1 in wall and 1 in fusion. Chytridiomycota has 1 cluster (2 rows). Taphrinomycotina has 2 clusters (3 rows).

### What this means for D9 step 3 (leave-clade-out)

With the JGI files (strict plus JGI file, which is the scenario that stands on identifiers):

- **Fusion: leave-clade-out is not possible.** Use descriptive recall only. Say so with the n (8 rows, 6 clusters; 5 of the 5 rows with ID are mapped, but at most 2 clusters per clade).
- **Adhesion: possible at floor 5** for Saccharomycotina and Sordariomycetes (fine level), or Saccharomycotina and Pezizomycotina (coarse level). No adhesion cluster crosses clades, so the hold-out is clean. The training side after holding out Saccharomycotina is Sordariomycetes, Eurotiomycetes, Ustilaginomycetes and Basidiomycota rows (18 rows); after holding out Sordariomycetes it is small too. Only Saccharomycotina passes floor 8.
- **Wall: possible** for Eurotiomycetes, Saccharomycotina and Sordariomycetes at all floors (fine level), but 10 to 13 of their clusters also have members in another clade, so the clean part is 16 clusters for Saccharomycotina and 4 for each of the other two. The result should be reported with this number.
- The locus-tag lookup (section 4c) changed no floor of 5 or 8: it adds 5 rows and one clade at floor 3 (adhesion, fine: Eurotiomycetes reaches 3 clusters).
- In every slice a leave-clade-out result rests on Saccharomycotina, Sordariomycetes, Eurotiomycetes. A conclusion about Basidiomycota (3 clusters in adhesion, 1 in wall), Chytridiomycota or Taphrinomycotina cannot come from Kiss rows (section 7 gives the second source).
- The UniProt-only stage said the opposite for strict mapping (not possible for any slice). That statement is superseded.
- Report every number with its n. Do not pool slices.

## 7. Second source (step 6)

For Basidiomycota and the non-hyphal clades the second source is **UniProt reviewed (Swiss-Prot) entries of the wall and fusion Pfam families**, plus literature rows (not counted here: I ran no literature search and cannot give a number). I used the 23 wall and fusion models of `analysis/wall_fusion_report/wall_fusion_models.tsv` (21 wall models and 2 fusion models).

Query for each model and clade (public UniProt REST, `rest.uniprot.org/uniprotkb/search`):

`(xref:pfam-<PFID>) AND (reviewed:true) AND (taxonomy_id:<clade taxon>)`

Clade taxa (from the NCBI lineage of the resolved organisms): Basidiomycota 5204, Saccharomycotina 147537, Chytridiomycota 4761. I did not count Mucoromycota or Zoopagomycota: no organism of the table is in them. The cap of 500 results was never reached. Per-model counts are in `results/second_source_uniprot_reviewed.tsv`, with the exact query on each line.

Distinct reviewed entries per model set (an entry that matches two models counts once):

| Model set | Clade | Distinct reviewed entries | Distinct organisms |
|---|---|---|---|
| cell_fusion_candidate_models | Basidiomycota | 5 | 5 |
| cell_fusion_candidate_models | Chytridiomycota | 0 | 0 |
| cell_fusion_candidate_models | Saccharomycotina | 13 | 12 |
| wall_remodeling | Basidiomycota | 11 | 7 |
| wall_remodeling | Chytridiomycota | 0 | 0 |
| wall_remodeling | Saccharomycotina | 55 | 13 |
| wall_structural_other | Basidiomycota | 0 | 0 |
| wall_structural_other | Chytridiomycota | 0 | 0 |
| wall_structural_other | Saccharomycotina | 26 | 8 |
| wall_synthesis | Basidiomycota | 18 | 3 |
| wall_synthesis | Chytridiomycota | 0 | 0 |
| wall_synthesis | Saccharomycotina | 9 | 2 |

- Basidiomycota: wall_remodeling 11 entries in 7 organisms; wall_synthesis 18 entries in 3 organisms; wall_structural_other 0; fusion 5 entries in 5 organisms (all PRM1, PF27758; Fig1 gives 0).
- Chytridiomycota: 0 in every set. Swiss-Prot cannot supply a chytrid row for these families.
- Saccharomycotina: 55, 9, 26 and 13 entries. This clade is already the best covered in Kiss.
- Limits: these are family-level hits, not checked for wall function (for example GH5 and the CE4 family also hold non-wall proteins). The entries are not clustered, so they are not independent units (Task 5.3 does that). For the adhesion slice this source gives no count, because no Pfam model list exists for it here.
- I did not add any row to a truth file.

## 8. What I could not do

- With the JGI files and the locus-tag lookup, 133 of 144 JGI-ID rows are strictly mapped; 11 are not (sections 4b and 4c). 28 of the 133 (Aspfu1 and Batde5) map by locus tag only and have no JGI file. Of the 99 rows mapped by file, 20 rely on a validated portal and 5 on a locus-tag lookup without name evidence, not on a check of the row itself.
- Whether the downloaded annotation is the one Kiss used cannot be proven from the files. Evidence per portal is the share of agreeing rows (Canalb1 33 of 34 rows with evidence, Sacce1 24 of 24, Schpo1 3 of 3, Trire2 2 of 2, Ustma2_2 4 of 4, Neucr2 9 of 9). For Bcine1 0 of 2 agree and the locus tags show that its IDs were renumbered (section 4c), so it is not used. For Hisca1, Schco3, Cryne_H99_1 and Crypa2 there is no name evidence; the 5 rows of Hisca1 and Cryne_H99_1 that were recovered rest on the locus-tag lookup only.
- Future truth sets and controls should be keyed by NCBI or UniProt accessions, not by JGI IDs (owner's rule of 2026-10-21). The JGI files are an archived one-off, outside git (`to_import/kiss2019_jgi/`, with `FILE_TABLE.tsv`); only their hashes are in `results/jgi_files.tsv`.
- I did not run the 30 percent identity clustering (Task 5.3). Kiss families are not identity clusters. A Kiss cluster can hold paralogs, and two Kiss clusters can be related. Cluster counts here are an upper bound of independent units.
- I did not count literature rows.
- In this task I read only Data 3 (S7) and Data 1 (MOESM5, downloaded to a scratch folder) of the supplement. I did not read the figure legends.
- Strain scope: for 8 organisms I used the species taxon. A gene-name candidate can come from another strain of the species.

## 9. Files

All in `analysis/wall_fusion_truth/` (README there). Tests: `tests/wall_fusion_truth/test_map_kiss.py` (3-row fixture, no network) and `test_map_kiss_jgi.py` (fixture FASTA files, no JGI data). The script runs offline from `results/query_cache.json`; the JGI step needs `--jgi-dir` (the committed results were made with it; a run without it writes the UniProt-only columns). New tables: `results/jgi_files.tsv`, `results/jgi_portal_summary.tsv`, `results/name_conflict_review.tsv`, `results/family_rep_clusters.tsv`. Hash of `results/row_mapping.tsv` (with the JGI columns and the locus-tag lookup): b4459f3a3d877dea9a817c916df5dc5209e0654fb82af5aee3c2227c3ad38b0d. The cache `results/query_cache.json` holds the 25 added answers (14 UniProt: 10 locus-tag queries and the 4 Hisca1 queries with the species taxon; 11 NCBI Taxonomy queries for the strain terms), so an offline rerun gives the same bytes. The UniProt-only run gave bc35a684d65495e1e132e6444dba5c705543b88607765137ff51228250ffbc14. Hash of `results/mapped_sequences.faa` (unchanged): 9a261a58b54307c01b1d3db2ba6924b0e30bf95403cee9752d10597082fee19a.

## 10. Owner decisions (final, 2026-10-22)

These replace the open questions of the earlier version. The script and tables follow them.

| # | Decision | How the tables follow it |
|---|---|---|
| 0 | Accept all 99 file-mapped rows (74 `row_verified`, 5 `locus_confirmed`, 20 `portal_validated`). Filter later. | `row_mapping.tsv` keeps `jgi_accept`; the scenario `strict_plus_jgi_file` uses all 133 strict rows, `strict_plus_jgi_file_row_verified` uses 108. |
| 1 | Keep the name-conflict rule: hyr1, bhp2 and bhp3 stay unmapped. | `results/name_conflict_review.tsv`: row, JGI ID, header name, file protein accession, gene names and length, named UniProt candidate accession and length, why. 3 rows. |
| 2 | The fine clade level is the primary hold-out. The coarse level is a second, harder test wherever eligible. | Column `role` (`primary` or `second_harder_test`) in `floors_per_slice.tsv` and `clusters_per_clade.tsv`. |
| 3 | The floor counts all clusters present. Every result is reported in two strata: shared (the cluster has members in another clade) and clade-only. A stratum is `estimated` only with at least 5 clusters of its own, else `descriptive`. | `clusters_per_clade.tsv`: `clusters_also_in_other_clade` (shared), `clusters_clade_only`, `status_all_clusters`, `status_shared`, `status_clade_only` (`none` when the stratum has 0 clusters). |
| 4 | Floor 5 for `estimated`. Clades with 3 or 4 clusters are shown as `descriptive`. | Same columns; `ESTIMATED_FLOOR = 5` in the script. |
| 5 | Fusion is a descriptive-only slice: the Pfam baseline is reported on its 5 rows, with no `estimated` claim and no custom model yet. | `stratum_status` returns `descriptive` for every fusion stratum. |
| 6 | The 35 A. fumigatus family representatives are flagged `family_rep`. They may be training seeds and appear in the baseline recall report. They are never test controls and never counted toward a floor. | Column `control_role` in `row_mapping.tsv` (`family_rep` 35, `kiss_row` 144, `unmapped_family_row` 2). They carry no clade, so no clade table counts them. `results/family_rep_clusters.tsv` lists the clusters that carry one. |

**Family representatives and clusters.** 27 clusters carry a family representative (adhesion 4, fusion 3, wall 20). 9 of them have no strictly mapped Kiss row, so the representative is the only mapped member (adhesion 2: clusters 338 and 3660; fusion 3; wall 4). The other 18 also have a Kiss row with a clade.

**Result under these decisions** (scenario `strict_plus_jgi_file`, floor 5 on all clusters; `results/clusters_per_clade.tsv`, `results/floors_per_slice.tsv`). Clusters: all, shared, clade-only; status of each stratum.

Primary (fine level):

| Slice | Clade | Rows | Clusters | Shared | Clade-only | Shared status | Clade-only status |
|---|---|---|---|---|---|---|---|
| adhesion | Saccharomycotina | 18 | 10 | 0 | 10 | none | estimated |
| adhesion | Sordariomycetes | 9 | 6 | 0 | 6 | none | estimated |
| wall | Saccharomycotina | 43 | 29 | 13 | 16 | estimated | estimated |
| wall | Eurotiomycetes | 23 | 16 | 12 | 4 | estimated | descriptive |
| wall | Sordariomycetes | 19 | 14 | 10 | 4 | estimated | descriptive |

Second, harder test (coarse level):

| Slice | Clade | Rows | Clusters | Shared | Clade-only | Shared status | Clade-only status |
|---|---|---|---|---|---|---|---|
| adhesion | Saccharomycotina | 18 | 10 | 0 | 10 | none | estimated |
| adhesion | Pezizomycotina | 14 | 9 | 0 | 9 | none | estimated |
| wall | Saccharomycotina | 43 | 29 | 13 | 16 | estimated | estimated |
| wall | Pezizomycotina | 42 | 23 | 13 | 10 | estimated | estimated |

Descriptive only: adhesion fine Eurotiomycetes (3 clusters, 5 rows) and Ustilaginomycetes (3 clusters, 3 rows), adhesion coarse Basidiomycota (4 clusters, 4 rows), every other clade with fewer than 5 clusters, and the whole fusion slice (5 mapped rows with a clade; at most 2 clusters per clade). Wall Eurotiomycetes and Sordariomycetes are eligible at floor 5 by all clusters, but their clade-only stratum has 4 clusters, so it is `descriptive`.
