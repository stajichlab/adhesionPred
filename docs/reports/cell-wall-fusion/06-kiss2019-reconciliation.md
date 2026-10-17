# 06. Kiss et al. 2019: reconciliation of reports 02 and 03

Paper: Kiss E, et al. Comparative genomics reveals the origin of fungal hyphae and multicellularity. Nat Commun 2019; PMID 31501435; PMC6733946; DOI 10.1038/s41467-019-12085-w.

## 1. Sources and method

- Full text: PubMed Central via the PubMed MCP full-text tool (PMC6733946, per PubMed ID conversion of PMID 31501435). Table 1 was read from the Europe PMC full-text XML (`<table-wrap id="Tab1">`) because nature.com and pmc.ncbi.nlm.nih.gov pages required a browser check and could not be fetched. Figure legends were not read (not retrievable); figure panel contents are therefore UNVERIFIABLE here.
- Supplements: downloaded `41467_2019_12085_MOESM{7,8,12}_ESM.xlsx` (HTTP 200, valid Excel files) and parsed with openpyxl (`python -I`), treated as data only.
- Mapping verified from the sheet title cell A1: MOESM7 = "Supplementary Data3: Complete list of collected 651 hyphal morphogenesis-related genes."; MOESM8 = "Supplementary Data4: List of gene families showing duplication events at the time of the emergence of filamentous fungi (BCZ nodes), based on COMPARE analysis."; MOESM12 = "Supplementary Data8: List of 414 gene families identified by ANOVA, with p < 0,05 considered as significant."
- A wrong PMC identifier was tried once first (it returned an unrelated lipid-raft paper) and discarded.

## 2. Answer: the two numbers measure different things, and both are right

Report 02's "23 adhesion-related families" and report 03's "45 adhesion genes" are the same set counted at two levels.

| Quantity | Value | Level | Where |
|---|---|---|---|
| Adhesion-related families | **23** | families (MCL clusters) | Main text, Results, section "No expansion of kinase, receptor and adhesive repertoires": "adhesion (23 clusters)" and "We identified 23 families of putative adhesion-related proteins in fungi, including adhesins, flocculins, hydrophobins, various lectins, and glycosylphosphatidylinositol-anchored cell wall proteins." Also the Fig. 2 repertoire analysis (figure legend not read). |
| Adhesion genes | **45** | genes | Methods, "Analyses of key multicellularity-related genes": "Families of adhesive proteins were identified based on experimentally characterized genes collected from the literature. We identified 45 genes, which mostly grouped into flocculins, lectins, hydrophobins, and other (GPI)-modified cell wall adhesins." |
| Adhesion rows in Data 3 | 45 rows, 23 distinct `cluster number` values | genes and families | MOESM7, sheet `gene list` |

So the 45 literature genes fall into 23 clusters. The paper treats adhesion as a repertoire-evolution test set (kinases 954 clusters, histidine kinases 96, receptors 183, adhesion 23), not as one of the "functional groups" of the hypha-gene analysis.

Both reports are therefore consistent with the primary source on this point. Report 02's phrase "defines 23 families" is accurate. Report 03's 45 is accurate but counts genes. The two numbers should not be compared. Neither is a "correlated set": the 23/45 adhesion set is curated from literature; the 414-family set is statistically derived.

The paper's text also says the adhesion families were not expanded at hyphal origin (17 to 14 copies in BCZ nodes). The 17 here is a copy number, unrelated to the 17 families of Table 1. Beware of this coincidence of numbers.

### Four distinct sets in the paper (do not mix)

| Set | Size | Nature | Source in file |
|---|---|---|---|
| Curated hypha-morphogenesis genes | 651 genes in 362 families (clusters) | curated from 519 publications; includes the 45 adhesion genes | Data 3 |
| Of those 362 families, BCZ-origin and conserved | 17 | families; Table 1 | main text and Table 1 |
| Of the 362, with duplications in BCZ nodes | 93 (25.7%) | families | text; Data 4 has 93 rows (computed) |
| Genome-wide ANOVA screen | 414 families (114 originate in BCZ nodes per text) | statistically correlated with hyphae, from all families | Data 8 |

## 3. Counts computed from the xlsx

### Data 3 (MOESM7), sheet `gene list`

Header row (row 3): `functional category | gene name | protein ID (JGI) | protein function | ortholog in Aspergillus fumigatus | cluster number | node of emergence | Source`. Rows 1 and 2 are title and blank. Data rows: 651 (rows with non-empty `functional category` and `gene name`; no other non-empty rows). Distinct `cluster number`: 362. Distinct `protein ID (JGI)`: 642 (a few IDs are shared by more than one row).

| `functional category` | Rows (genes) | Distinct `cluster number` |
|---|---|---|
| polarity establishment | 107 | 63 |
| vesicle transport | 107 | 57 |
| cell wall biogenesis | **94** | 40 |
| signaling | 89 | 46 |
| septation | 58 | 33 |
| actin cytoskeleton | 55 | 35 |
| transcriptional regulation | 53 | 41 |
| adhesion | **45** | **23** |
| microtubular transport | 32 | 19 |
| cell cycle regulation | 11 | 5 |
| Total | 651 | 362 |

Note: the per-category cluster counts sum to 362, so the categories do not share clusters.

Discrepancy inside the paper: the text says nine functional groups (no adhesion group) with 651 genes in 362 families, and gives 107 / 107 / 11 as the largest/smallest. The file has ten categories, and 651 and 362 include the adhesion rows (nine non-adhesion categories sum to 606 genes and 339 clusters). So the paper's "nine groups" sentence omits adhesion while the 651 and 362 totals include it. Counting only the nine groups is not what the totals in the text mean.

### Data 4 (MOESM8), sheet `duplications in BCZ nodes`

Header (row 3): `number of duplications in BCZ nodes | functional category | cluster number | A.fumigatus ortholog | S. cerevisiae ortholog | ortholog from other species`. 93 data rows, 93 distinct `cluster number` values; matches the text's "Ninety-three (25.7%) of the 362". By category: signaling 18, cell wall biogenesis 16, polarity establishment 14, transcriptional regulation 11, vesicle transport 9, actin cytoskeleton 7, septation 7, microtubular transport 6, adhesion 3, cell cycle regulation 2. Adhesion rows: clusters 417 (glfB), 4291 (iff4/C. albicans), 7971 (Afu1g01540), 1 duplication each.

### Data 8 (MOESM12), sheet `ANOVA`

Header (row 7): `cluster number | p value | F value | A. fumigatus ortholog | A. fumigatus ortholog function | S. cerevisiae ortholog | S. cerevisiae ortholog function | C. albicans ortholog | UC fungi | MC fungi`. 414 data rows, 414 distinct `cluster number`. There is no functional category column in this sheet, so "cell wall" or "adhesion" for these rows can only come from annotation text.

### Cell wall genes in the hypha set

94 genes (Data 3, `functional category` = `cell wall biogenesis`), in 40 clusters. By `node of emergence`, the 17 genes at BCZ-type nodes are the rows of clusters 293, 374, 435, 843, 3398 and 9065 (see section 5).

## 4. Table 1 of the paper: the 17 families (Europe PMC XML)

Caption: "List of the 17 gene families whose emergence coincides with the evolution of hyphae". Columns: Emergence of gene family; A. fumigatus ortholog; S. cerevisiae ortholog; Functional category and putative function. Text of the Results: two transcriptional regulator families (StuA, MedA), six cell wall, three actin, three polarity, two signaling and one cell cycle (sum 17).

| Emergence node | A. fumigatus | S. cerevisiae | Category (Table 1) |
|---|---|---|---|
| MRCA Dikarya + Mucoro + Zoopago + Chytrid + Blastocladio | Afu7g03870 | PAN1 | Actin cytoskeleton |
| same | crh3 | UTR2 | **Cell wall biogenesis** |
| same | gel7 | GAS1 | **Cell wall biogenesis** |
| same | Afu6g04940 | BNR1 | Polarity establishment |
| same | Afu4g04120 | BEM1 | Polarity establishment |
| same | stuA | PHD1 | Transcriptional regulation |
| same | medA | NA | Transcriptional regulation |
| MRCA Dikarya + Mucoro + Zoopago + Chytrid | Afu6g07910 | SLM1 | Actin cytoskeleton |
| same | Afu8g04520 | SLA1 | Actin cytoskeleton |
| same | Afu4g06130 | WHI2 | Cell cycle regulation |
| same | Afu4g00620 | DFG5 | **Cell wall biogenesis** |
| same | Afu8g02320 | NA | **Cell wall biogenesis** (ortholog of N. crassa cps1) |
| same | chsD | NA | **Cell wall biogenesis** (class VI chitin synthase) |
| same | rgsB | RAX1 | Polarity establishment |
| same | Afu2g08800 | SSY1 | Signaling |
| same | ricA | NA | Signaling |
| MRCA Dikarya + Mucoro + Zoopago | kre6 | KRE6 | **Cell wall biogenesis** |

Count: 7 + 9 + 1 = 17. Cell wall: 6. Table 1 contains no adhesion-category, hydrophobin or fusion family.

## 5. Cross-check of Table 1 cell wall families against Data 3 (computed)

All six are in Data 3 under `cell wall biogenesis`, with matching `node of emergence`:

| Family | Data 3 cluster | `node of emergence` (Data 3) | Data 3 gene rows |
|---|---|---|---|
| crh3 / UTR2 | 374 | lca of Dikarya, Mucoromycota, Zoopagomycota, Chytridiomycota, Blastocladiomycota | crh1, crh2_utr2, csf4 (A. fumigatus ortholog column: crh3/Afu3g09250) |
| gel7 / GAS1 | 435 | same five-clade node | gas1, phr1, gel4, gel7_Afu6g12410 |
| DFG5 | 843 | lca of Dikarya, Mucoromycota, Zoopagomycota, Chytridiomycota | dfg5 (2 rows) |
| cps1 (Afu8g02320) | 3398 | same four-clade node | cps1 |
| chsD | 9065 | same four-clade node | chs6, chsD_Afu1g12600 |
| kre6 | 293 | lca of Dikarya, Mucoromycota, Zoopagomycota | kre6 |

The Table 1 grouping (A: five clades; B: four clades; C: three clades) matches these nodes. Data 4 lists all six clusters with these duplication counts: crh3 374 (6), gel7 435 (4), kre6 293 (2), dfg5 843 (1), cps1 3398 (3), chsD 9065 (1).

## 6. Claim-by-claim check of report 03

Statuses: CONFIRMED / CORRECTED / UNVERIFIABLE (not checkable from the three xlsx files and the full text I could read).

| # | Claim in report 03 | Status | Detail |
|---|---|---|---|
| 1 | 651 hypha genes in 362 families, 519 publications, from 9 functional groups (section 1.1) | CORRECTED | 651 and 362 confirmed (file and text). "519 publications": stated in text, not checked in the file. The file has 10 categories (adhesion included in the 651); the text calls them nine groups. |
| 2 | Data 3 category counts (polarity 107, vesicle 107, cell wall 94, signaling 89, septation 58, actin 55, transcription 53, adhesion 45, microtubule 32, cell cycle 11) | CONFIRMED | All ten counts match. |
| 3 | Cell wall: 94 genes | CONFIRMED | 94 rows, 40 clusters. |
| 4 | Adhesion: 45 genes | CONFIRMED | 45 rows. |
| 5 | Adhesion "24 distinct clusters in the file" (section 2.3) | CORRECTED | **23** distinct `cluster number` values. This equals the paper's "23 clusters"; the report's 24 is an off-by-one and also made the 23 look like a separate quantity. |
| 6 | "The paper text gives 9 groups... The 9-group totals in the text are probably the file minus adhesion. I did not reconcile row counts" | CORRECTED | Not minus adhesion: the 651/362 totals include adhesion (nine groups alone sum to 606 genes). The text lists nine groups and omits adhesion; the file lists ten. The 107 / 107 / 11 figures in the text match the file. |
| 7 | "This [adhesion category] is the set behind the 23 adhesion clusters analysis" | CONFIRMED | 45 genes in 23 clusters; Methods says 45 genes. |
| 8 | 17 families in Table 1; six cell wall: crh3, gel7/GAS1, DFG5, cps1, chsD, kre6 | CONFIRMED | Matches Table 1 exactly (section 4). |
| 9 | Table 1 group membership (A: 7 families, B: 9, C: 1) and per-family annotations | CONFIRMED | Matches the XML; "polarity" is "Polarity establishment" in the table (text calls it "polarity maintenance"). |
| 10 | Table 1 "surface or intracellular" column and "GPI-anchored" for gel7/DFG5 | UNVERIFIABLE (partly) | The column is the report's own inference. File text does say "GPI anchored protein" for gas1 and "GPI-anchored protein"/"GPI-linked cell surface protein" for dfg5. cps1 membrane location is not in the files. |
| 11 | 414 families; 114 BCZ-origin; 20 transporters; 84 uncharacterized | CONFIRMED (text); 414 CONFIRMED (file). | Data 8 has 414 rows. The 114/20/84 figures are stated in text; the sheet has no column to reproduce them. |
| 12 | 93 families with duplications in BCZ nodes | CONFIRMED | Text and Data 4 (93 rows). |
| 13 | Cell wall cluster table, clusters 293 (kre6), 374 (crh1/crh2/csf4), 435 (gas1/phr1/gel4/gel7), 843 (dfg5), 166, 171, 209, 234, 3228, 4470, 4935, 635, 636, 3602, 12180 and their nodes | CONFIRMED | Checked in Data 3. Minor: cluster 209 node is "lca of Opisthokonta, E.histolytica, T.trahens" (the report says "Opisthokonta + amoebozoa"; T. trahens is not an amoebozoan). Clusters 3398 (cps1) and 9065 (chsD) are not in the report's table but are in Data 3. |
| 14 | csf4 "GPI anchored cell wall glycosidase"; Gel4 "CAZy family 17" | CONFIRMED | Quoted text in `protein function`. |
| 15 | Cluster 4166: bgt1/bgt2 with JGI IDs 3834_Neucr2 and 5869_Neucr2, source PMID 26541633 | UNVERIFIABLE | Rows exist (bgt2_eglc, bgt1; node lca of Opisthokonta); JGI IDs and source citation were not re-read. |
| 16 | Adhesion cluster table (220, 417, 1012/2024/7971/48004 ncl1-4, 4291, 5266, 5817, 6453, 7118, 7187, 7276, 8398, 8498, 8594/16358/14985, 16148, 44817/50327/126186/140788) | CONFIRMED | All 23 clusters present with the listed genes; hfb1 is labelled "galectin" in cluster 7118. |
| 17 | Adhesion set "mostly yeast and Candida" (section 1.1 caveats) | CORRECTED | By JGI suffix: Canalb1 11, Sacce1 8 (19 of 45, 42 percent); the rest are N. crassa 6, hydrophobins from U. maydis/S. commune/C. cinerea/others. Largest share, not a majority. |
| 18 | Fusion-adjacent Data 3 rows: so_soft (cluster 13739, Pezizomycotina, "promote cell-cell fusion"), ste3/preA (668), ste2 (9439), adn3 (8029) | CONFIRMED (so_soft, 668, 9439); adn3 not rechecked | so_soft: `septation`, node "lca of Pezizomycotina"; JGI ID 2491_Neucr2; source citations include Fleissner & Glass 2007 (SO in hyphal fusion). |
| 19 | Data 8 hits (435 gel2 0/46; 374 crh3 2/51; 179 eng2 1/38; 2065 bgt3 0/41; 635 0/51; 171 0/51; 3660 0/40; 3355 PRM1 2/48; 2511 0/41; 1764 3/49; 1348 GPI12 3/52; 338 EMP70 3/52; 209/234/368 3-4 / 43-53) | CONFIRMED | UC/MC values match (`UC fungi`, `MC fungi`). 4240: the A. fumigatus and S. cerevisiae columns are NA; the C. albicans column reads "C2_08770C Putative adhesin-like protein"; 0/38. |
| 20 | Data 3 cluster 3359 "glycosyl-transfer domain" | UNVERIFIABLE as surface-related | Annotation reads "predicted transferase activity"; it does not state cell wall. Not included in the TSV. |
| 21 | Download list, Data N = MOESM(N+4) | CONFIRMED for Data 3, 4, 8 | Title cells match (section 1). Others unchecked. |
| 22 | "Cell fusion not a focus of Kiss et al." | CONFIRMED | The words "fusion", "anastomosis" and "fuse" occur zero times in the Europe PMC full-text XML (including references). Fusion appears only in supplement annotations. |

## 7. Report 02 check

| Claim in report 02 | Status |
|---|---|
| "defines 23 families of putative adhesion-related proteins" | CONFIRMED (Results, "23 families of putative adhesion-related proteins"). |
| "reports no expansion of adhesive, kinase or receptor families when filamentous fungi arose" | CONFIRMED (abstract and Results). |
| "The full text has no occurrence of 'fusion' or 'anastomosis'" | CONFIRMED for the XML full text. The supplement does contain fusion annotations (so_soft; Data 8 clusters 2511, 3355, 1764). |
| "Hydrophobins are in the 23 adhesion-related families" | CONFIRMED. Hydrophobin-annotated rows are in clusters 220, 5817, 7118, 16148, 50327, 126186 (hfb1 is labelled "galectin"). |

The wording in report 02 that the paper "defines" the 23 families is slightly generous: they are the families containing 45 literature-curated adhesion genes (Methods), not families defined by a sequence profile. Per-family Pfam/InterPro IDs are not given in the files read.

## 8. Extracted rows

TSV: `docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` (181 rows plus header; columns `gene_or_family_id, cluster_number, category, source_sheet, source_category_value, protein_id_jgi, ortholog_Afu, node_of_emergence, annotation, selection_rule`).

| Source | Category | Rows | Rule |
|---|---|---|---|
| Data 3 `gene list` | cell wall | 94 | `functional category` = cell wall biogenesis |
| Data 3 `gene list` | adhesion | 45 | `functional category` = adhesion |
| Data 3 `gene list` | cell fusion | 1 | so_soft; protein function says "promote cell-cell fusion" |
| Data 3 `gene list` | pheromone receptor | 4 | ste3, pra1, ste2, preA; keyword, labelled fusion-adjacent |
| Data 4 `duplications in BCZ nodes` | cell wall / adhesion | 16 / 3 | `functional category` column (family-level; duplicates Data 3 clusters) |
| Data 8 `ANOVA` | cell wall / adhesion / fusion | 12 / 3 / 3 | hand-picked by annotation keywords (cell wall, glucan, chitin, GPI anchor, adhesin, adhesion, fusion, syncytium); the sheet has no category column |

There is no "cell fusion" functional category in the files. Vesicle fusion (SNAREs, exocyst, vacuole fusion) rows were not included. The MAPK and G protein pheromone-pathway rows (fus3, mpkB, ste50, gpaA/gpgA) were also left out, as they are not fusion factors by annotation. Rows overlap by design: a gene may appear in Data 3 and (as a family) in Data 4 and Data 8; cross-reference by `cluster_number`.

## 9. Open points

- Table 1 gives A. fumigatus and S. cerevisiae ortholog names only; no Pfam/InterPro IDs. Mapping to our Pfam families needs a separate step (report 05).
- Figure legends (Fig. 2 repertoire, Fig. 3 category counts) were not read, so the 17 to 14 copy-number change and the Fig. 2 panels are quoted from text only.
- `uv run` created `.venv` and `uv.lock` in this worktree during the first openpyxl attempt. An attempt to delete them was denied by the permission system, so they remain untracked and should be removed by the user (`.venv` is probably gitignored; `uv.lock` shows in `git status`).
