# CGD biofilm literature links: leads for curation

*2026-10-07. Source: `to_import/biofilm_literature_topic_search_results.csv` (export from the Candida Genome Database, topic "Biofilms"; the file is in `to_import/` and is not tracked). This note records what the file holds and which leads it gives. Nothing here is curated evidence.*

## What the file is

- One row per gene and paper pair. Columns: Topic, PubMed ID, Year, Reference, Gene Name, Systematic Name, Species.
- 4,128 rows, 925 papers (1994 to 2026), 1,599 systematic names, 1,170 named genes.
- 3,395 rows are *C. albicans* SC5314. Others: *C. glabrata* 193, *C. auris* 95, *C. parapsilosis* 42, *C. tropicalis* 32. 355 rows have no species.
- 664 rows have no gene name. 171 of them belong to one paper, PMID 21414038 (Bonhomme 2011).
- A link means CGD filed the paper under the topic for that gene. It does not show that the gene has a role in biofilms. Large expression screens link many genes (PMID 21414038: 420 gene links, 19322777: 101, 16151249: 101, 15964282: 95).

## Overlap with our curated tables

- 345 curated gene names exist in `data/curated/` (adhesins, biofilm, hard negatives, seeds).
- 216 of the 1,170 named genes in the CGD file are already curated.
- Genes with most papers: ALS3 (145), HWP1 (127), EFG1 (100), ALS1 (69), BCR1 (62), ECE1 (53), ERG11 (53), CDR1 (50), TEC1 (47).
- A name filter (my choice, not a CGD class) picked 120 genes from cell wall and adhesin families (ALS, HWP, EAP, RBT, PGA, YWP, SIM, IFF, HYR, ECM, CRH, SAP, CHT, KRE, AWP, EPA and others). 70 of them are not in our curated tables. List: `data/curated/biofilm/cgd_biofilm_literature_leads.tsv` (gene, systematic names, species, paper count, PMIDs, in-curated flag).
- Top uncurated ones by paper count: SAP4 (14), CHT2 (10), FLO8 (8), PHR2 (8), KRE1 (7), RBT4 (7), CHT3 (4), KRE6 (4), PGA45 (4).

## Where each lead could be used

| Agent task | Use |
|---|---|
| 03 repeat-mechanism controls | RBT1 (Q59TP1) has 13 CGD biofilm papers: 17586721, 19151323, 19541532, 19837954, 20150241, 20709785, 23194472, 23766273, 26772652, 29062088, 29085811, 33580263, 41165448. None is checked here for a repeat statement. Braun 2000 (PMID 10978273) says "alleles of two different sizes" per CGD. The full text is still to be read. HWP1 and ALS3 have the most papers and may supply repeat-region statements. |
| 04 wall-family domain controls | PGA, AWP, EPA, RBT and HWP-family genes are candidate members or non-members for the Pfam families (Flo11, GLEYA, Hyphal_reg_CWP, Candida_ALS). |
| 07 hard negatives | SAP (secreted proteases), CHT (chitinases), KRE, CHS, MNN, GAS (wall synthesis and remodeling enzymes) are surface-located or secreted and are not adhesins by function. They are candidate hard negatives, to be checked one by one. |
| Biofilm label table (`data/curated/biofilm/`) | Direction of effect (promotes or restrains) is `unknown` for nearly all rows. The papers linked here can fill it. |

## Limits

- The CGD topic is a literature link, not an annotation of function.
- 179 of the 207 rows of `biofilm.tsv` are already *C. albicans*. This file adds more of the same bias.
- The family filter will miss genes without a standard name and may include genes that are not surface proteins.
- Open each paper before citing it. The CSV gives a citation line only.

## Second file: cell wall and adherence

`to_import/cellwall_CGD_literature.csv` (CGD export, same columns). 3,195 rows: topic "Cell wall properties/components" 2,086 and "Adherence" 1,109. 763 papers, 910 gene names, 1,166 systematic names. Species: *C. albicans* 2,247, *C. glabrata* 551, *C. auris* 71, *C. parapsilosis* 58, 255 without species. Look up here first for tasks 03 (repeat statements), 04 (family members) and 07 (hard negatives).

Cross-reference (2026-10-07): 909 named genes, 206 already in our curated tables, 703 not. 441 genes have at least one "Adherence" row; 278 of them are not curated. List: `data/curated/biofilm/cgd_adherence_literature_leads.tsv` (columns: gene, systematic names, species, paper count, adherence rows, cell wall rows, PMIDs, in-curated flag). Uncurated genes with the most papers: PHR2 (17), GSC1 (16), CDR1 (12), SAP4 (12), INT1 (12), ERG11 (11), ECE1 (11), SOD5 (11), GSL1 (9), HSP90 (9). Many of these are not surface adhesins (drug pumps, signalling, wall synthesis). A paper link under "Adherence" shows that CGD filed it there, not that the protein binds a surface. RBT1 has 3 adherence rows and 9 cell wall rows. ALS7 has 20 papers in the curated group (see the label conflict in item 13 of the worklist).
