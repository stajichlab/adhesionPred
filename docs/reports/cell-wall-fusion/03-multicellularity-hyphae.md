Errata applied 2026-10-17 from reports 05 to 07; items still unverified are marked.

# 03. Multicellularity and hyphal-formation genes vs cell wall remodeling, fusion and adhesion

Date: 2026-10-17. Literature research only; no code changed.
Sources: PubMed metadata (via the PubMed MCP), PMC full text of PMC6733946 through Europe PMC XML, and the paper's Supplementary Data xlsx files, which were downloaded and parsed. Anything not read directly is marked "unverified".

## 0. What the two papers are

The two PMIDs are different kinds of paper. Only the first is about multicellularity.

| | PMID 31501435 | PMID 26541633 |
|---|---|---|
| Title | Comparative genomics reveals the origin of fungal hyphae and multicellularity | Role of BGT-1 and BGT-2, two predicted GPI-anchored glycoside hydrolases/glycosyltransferases, in cell wall remodeling in Neurospora crassa |
| Authors | Kiss E, Hegedus B, Viragh M, Varga T, Merenyi Z, Koszo T, Balint B, Prasanna AN, Krizsan K, Kocsube S, Riquelme M, Takeshita N, Nagy LG | Martinez-Nunez L, Riquelme M |
| Year, journal | 2019, Nat Commun 10:4080 | 2015 (online Nov 2015), Fungal Genet Biol 85:58-70 |
| DOI | 10.1038/s41467-019-12085-w | 10.1016/j.fgb.2015.11.001 |
| Type | Comparative genomics: gene-family reconstruction, ancestral state reconstruction, enrichment tests | Single-gene-pair cell biology: GFP tagging, live-cell imaging, deletion mutants |
| Full text read | Yes (PMC6733946, all sections and Table 1) | No. Abstract only. No PMC copy exists and the publisher page returned HTTP 403. |

PMID 26541633 is not a multicellularity paper. It is a localization and phenotype study of two N. crassa genes. It still connects to PMID 31501435 directly: both BGT genes are entries in the Kiss et al. hypha-gene dataset (section 2.3). Riquelme M is a co-author of both.

---

## 1. Summaries

### 1.1 Kiss et al. 2019 (PMID 31501435)

- **Question.** How was the genetic toolkit for hyphal multicellularity assembled? Did fungi expand kinases, receptors and adhesins, as other multicellular lineages did?
- **Data.** 72 proteomes. Per the text: 4 plesiomorphically unicellular fungi, 41 hyphal, 13 secondarily simplified (yeast-like) fungi, and 14 non-fungal relatives. All fungal phyla are covered. Microsporidia were left out.
- **Method.**
  1. Species tree from 455 single-copy orthologs (RAxML, IQ-TREE, PhyloBayes).
  2. Ancestral state reconstruction of hyphae (threshold model, phytools). The origin falls on three early nodes, the "BCZ nodes": the Blastocladiomycota, Chytridiomycota and Zoopagomycota splits. Posterior probability for the hyphal state is 0.53, then 0.68, then 0.92.
  3. MCL protein families, gene-tree/species-tree reconciliation (Notung), and Dollo mapping of duplications and losses (COMPARE pipeline).
  4. A manually curated set of 651 hypha-morphogenesis genes in 362 families, from 519 publications, in 9 functional groups.
  5. A genome-wide permutation ANOVA for families whose duplication rate peaks in BCZ nodes, giving 414 families.
  6. Gene-structure, domain-architecture and yeast gene-loss analyses.
- **Main findings.**
  - Hyphae arose early, probably gradually across the BCZ nodes.
  - There was no expansion of Ser/Thr kinases, receptors or adhesion families at the origin of hyphae. Adhesion families (23 clusters) contracted slightly from 17 to 14 copies in BCZ nodes. They expanded later in the Agaricomycotina (class 1 hydrophobins and Cfl1 homologs) and in the Saccharomycotina (yeast adhesins and lectin-like cell wall proteins).
  - 50% of hypha families (181) predate hyphae. 45% (164) arose after. Only 17 (4.7%) originated in BCZ nodes and stayed conserved (Table 1).
  - Duplications in BCZ nodes were significantly enriched in two functional groups: cell wall biogenesis and transcriptional regulation. No single family was individually significant.
  - Phagocytosis-related genes were exapted for hyphae: Arp2/3 retained, ELMO and DOCK retained, WASP retained, WAVE only in early-diverging fungi, WASH lost.
  - 414 families show correlated evolution with hyphae. 114 originated in BCZ nodes. 20 are transporters (possible osmotrophy confound) and 84 are uncharacterized.
  - Yeasts lose hypha genes less than expected from genome-wide loss. Adhesion and microtubule transport are among the six groups lost more than expected.
- **Caveats.**
  - The authors say that traits co-distributed on the tree (osmotrophy, cell wall changes) can produce spurious correlation with hyphae.
  - The adhesion set rests on 45 literature genes, mostly yeast and Candida [corrected 2026-10-17, see 06 section 6 row 17: by JGI suffix Candida albicans 11 and S. cerevisiae 8 = 19 of 45 (42 percent); the largest share, not a majority; the rest are N. crassa 6 and hydrophobins from U. maydis, S. commune, C. cinerea and others]. The authors note that adhesins in early-diverging fungi are poorly annotated.
  - Hyphae are scored as a binary trait.
  - Family-level "origin" depends on homology detection. Divergent homologs can look like de novo families.
  - The paper counts duplications at nodes, not expression or function.

### 1.2 Martinez-Nunez and Riquelme 2015 (PMID 26541633)

Abstract only; full text not accessed.

- **Question.** Where do N. crassa BGT-1 (NCU06381) and BGT-2 (NCU09175) act, and what do they do in cell wall remodeling?
- **Data and method.** N. crassa. Each protein was tagged with GFP: after the signal peptide for BGT-1, and before the predicted GPI attachment site for BGT-2. Live imaging, single and double deletion mutants, calcofluor white (CW) and Congo red (CR) sensitivity, and conidiation phenotypes.
- **Findings (from the abstract).**
  - Both are predicted GPI-anchored glycoside hydrolases with additional glycosyltransferase activity.
  - Both sit at the plasma membrane of the hyphal apical dome, excluding the very tip and the Spitzenkorper (Spk), where chitin and beta-1,3-glucan synthases occur.
  - Both also localize at the leading edge of new septa, at unreleased interconidial junctions, at conidial poles, at germling and hyphal fusion sites, and at branch emergence sites.
  - The bgt-2 deletion and the double mutant are more resistant to CW and CR, which suggests an altered wall. Growth rate is not strongly reduced. The mutants have conidiation defects, with conidia not separating.
  - The authors propose a division of labor at the tip. Glucans are synthesized at the tip, and hydrolyzed in the subapical zone to make ends for branching and crosslinking.
- **Caveats.**
  - The predictions are in silico (GPI site, GH and GT annotation). No biochemical activity is shown in the abstract.
  - Localization is by GFP fusion. The tag position was chosen to avoid the GPI signal.
  - Family assignment (for example GH17 or GH72) is not stated in the abstract: "unverified".

---

## 2. Gene sets and families named

All rows are from Kiss et al. unless the PMID column says otherwise. Per-family Pfam/InterPro IDs are not given in the paper text for these families. They are "unverified" here except where the text names a domain. "Function" is the paper's own annotation.

### 2.1 Table 1 of Kiss et al. (the 17 families originating in BCZ nodes)

Source: Table 1, main text, PMID 31501435. Parsed from the Europe PMC XML. Entries without a gene name in the table are marked "NA" there. The table has three emergence groups.

Group A: MRCA of Dikarya, Mucoromycota, Zoopagomycota, Chytridiomycota and Blastocladiomycota.

| Family (A. fumigatus / S. cerevisiae ortholog) | Group | Function (paper) | Surface or intracellular (my inference) |
|---|---|---|---|
| PAN1 (Afu7g03870 / PAN1) | actin | endocytic adaptor; recruits Arp2/3 | intracellular |
| crh3 / UTR2 | cell wall | chitin transglycosylase at polarized-growth sites; transfers chitin to beta-1,6 and beta-1,3 glucans | cell wall; the family also contains GPI-anchored members (csf4 in Data 3) |
| gel7 / GAS1 | cell wall | beta-1,3-glucanosyltransferase; cell wall remodeling at germination and branching | GPI-anchored (Data 3 says "GPI anchored protein" for Gas1) |
| BNR1 (Afu6g04940) | polarity | formin; actin cable assembly | intracellular |
| BEM1 (Afu4g04120) | polarity | scaffold for actin reorganization | intracellular |
| stuA / PHD1 | transcription | APSES; yeast-filament transition | nuclear |
| medA / NA | transcription | APSES-like regulator | nuclear |

Group B: MRCA of Dikarya, Mucoromycota, Zoopagomycota and Chytridiomycota.

| Family | Group | Function (paper) | Surface or intracellular |
|---|---|---|---|
| SLM1 (Afu6g07910) | actin | PtdIns(4,5)P2 effector | intracellular |
| SLA1 (Afu8g04520) | actin | endocytic patch protein | intracellular |
| WHI2 (Afu4g06130) | cell cycle | stress response; filamentous growth | intracellular |
| DFG5 (Afu4g00620) | cell wall | mannosidase; bud formation and filamentous growth | GPI-anchored (Data 3: "GPI-anchored protein", alpha-1,6-mannanase) |
| cps1 (Afu8g02320) | cell wall | polysaccharide synthase ortholog of N. crassa cps1 | integral membrane enzyme (unverified) |
| chsD (class VI chitin synthase) | cell wall | chitin biosynthesis | plasma membrane enzyme |
| rgsB / RAX1 | polarity | bud site selection | membrane-associated |
| SSY1 (Afu2g08800) | signaling | SPS amino-acid sensor | plasma membrane |
| ricA | signaling | G-protein GEF | intracellular |

Group C: MRCA of Dikarya, Mucoromycota and Zoopagomycota.

| Family | Group | Function (paper) | Surface or intracellular |
|---|---|---|---|
| kre6 / KRE6 | cell wall | beta-1,6-glucan biosynthesis | Golgi/plasma membrane enzyme |

So 6 of 17 are cell wall biogenesis. None are adhesins, hydrophobins, lectins or cell-fusion factors.

### 2.2 Other families named in the main text

| Family or gene | How identified | Pattern | Where in paper |
|---|---|---|---|
| Ser/Thr kinases (954 clusters), hybrid HKs (96), receptors (183), adhesion (23) | curated and InterPro-based sets | no expansion at BCZ nodes; kinases contracted 572 to 529; adhesion 17 to 14 | Fig. 2 and "No expansion of kinase, receptor and adhesive repertoires" |
| Hybrid histidine kinases (class III, X) | domain-based | expanded in MRCA of Chytridiomycota + higher fungi (24 dup, 10 loss); second wave in Mucoromycota + Dikarya | same section; Supplementary Fig. |
| Pheromone receptors | GPCR scan | the only GPCR family conserved across fungi (19 GPCR families in fungi vs 135-583 genes in animals) | same section |
| Class 1 hydrophobins and Cfl1 homologs | adhesion set | expanded in Agaricomycotina, tied to fruiting bodies | same section |
| Yeast adhesins and lectin-like cell wall proteins (Candida) | adhesion set | expanded in Saccharomycotina | same section |
| Cell wall biogenesis, transcription regulation | enrichment of BCZ duplications across the 9 groups | the only two groups significantly enriched (Fisher, BH-corrected) | "The evolution of hypha morphogenesis genes" |
| Chitin synthases (class V and VII gain a myosin motor domain) | domain-architecture test | more domains in multicellular fungi; 4 families in total | Supplementary Data 7 |
| Arp2/3, ELMO, DOCK, WASP retained; WAVE early-diverging only; WASH lost | phagocytosis gene list (106 genes) | exaptation | Fig. on phagocytosis genes |
| Bgt3, RgsB, Gel2 (A. fumigatus), Bem1, Rax1 | 414-family screen | originate or duplicate at BCZ nodes | "Genome-wide screen" section |
| Tropomyosin TPM1, Dip1, Dpp1 (farnesol synthesis) | 414-family screen | BCZ origins | same section |
| NADPH oxidases | yeast loss analysis | lost in all yeast clades; complete loss in Saccharomycotina | yeast section |
| Dynactin, gamma-tubulin complex, kinesins, dynein | yeast loss analysis | contracted or lost in yeasts | yeast section |

### 2.3 Gene lists from the supplement (parsed from the xlsx files)

Only the supplement carries gene-level detail, so I parsed it. The counts below are from the files, not the text.

**Supplementary Data 3 (651 genes).** Rows per functional category:

| Category | Rows |
|---|---|
| polarity establishment | 107 |
| vesicle transport | 107 |
| cell wall biogenesis | 94 |
| signaling | 89 |
| septation | 58 |
| actin cytoskeleton | 55 |
| transcriptional regulation | 53 |
| adhesion | 45 |
| microtubular transport | 32 |
| cell cycle regulation | 11 |

Two observations:
- The paper text gives 9 groups and does not count "adhesion" as one of them. The file has 10 categories, including "adhesion" (45 rows). This is the set behind the "23 adhesion clusters" analysis. The 9-group totals in the text are probably the file minus "adhesion". This is an inference; I did not reconcile row counts. [corrected 2026-10-17, see 06 section 3 and 6 row 6: this is wrong. The 651 genes and 362 families in the text include the adhesion rows; the nine non-adhesion categories alone sum to 606 genes and 339 clusters. The text lists nine groups and omits adhesion, while the file has ten categories. The 107 / 107 / 11 figures in the text match the file.]
- Columns: functional category, gene name, JGI protein ID, protein function, A. fumigatus ortholog, cluster number, node of emergence, source citation.

Cell wall category examples, cluster number and node of emergence from Data 3:

| Cluster | Genes (examples) | Emergence node | Note |
|---|---|---|---|
| 166 | PMT1/2/4 O-mannosyltransferases | LCA of Opisthokonta | O-mannosylation of surface proteins |
| 171 | exg1 | Opisthokonta | exo-1,3-beta-glucanase |
| 209 | chs1/2/3/8, chsA/B/C/G (classes I, II, III) | Opisthokonta + amoebozoa [corrected 2026-10-17, see 06 section 6 row 13: Data 3 node is "lca of Opisthokonta, E.histolytica, T.trahens"; T. trahens is not an amoebozoan] | chitin synthases |
| 234 | chsE, csmB (class V, VI-myosin type) | LCA of Fungi | myosin-motor chitin synthase |
| 293 | kre6 | BCZ-type node (Dikarya + Muco + Zoopago) | beta-1,6-glucan |
| 374 | crh1, crh2 (utr2), csf4 | Dikarya + Muco + Zoopago + Chytrid + Blasto | chitin transglycosylase; csf4 described as "GPI anchored cell wall glycosidase" |
| 435 | gas1, phr1, gel4, gel7 | same BCZ node | GPI-anchored beta-1,3-glucanosyltransferases, Gas/Phr/Gel; "CAZy family 17" for Gel4 as stated in the file [unresolved conflict, see 05] (reports 01 and 04 put GEL proteins in GH72; 05 could not resolve this; no side is chosen here) |
| 4166 | N. crassa bgt2 (NCU09175?), bgt1 (NCU06381?) | LCA of Opisthokonta | source citation is PMID 26541633; see below |
| 843 | dfg5 | BCZ-type node | GPI-anchored alpha-1,6-mannanase |
| 3228 | ecm33 | Opisthokonta | GPI-linked cell wall protein |
| 4470 | sun41 / sun1 | LCA of Ascomycota | beta-1,3-glucan modifying glycosidase; Sun-family |
| 4935 | pir1 | Ascomycota | cell wall protein crosslinking beta-1,3 glucan |
| 635 | fks1, bgs1, bgs4, gsc1 | Opisthokonta | beta-1,3-glucan synthase |
| 636 | cts1, cht2, chiA1 | Eukaryotes | chitinase; cht2 and chiA1 are GPI-linked |
| 3602 | ags1, ags2 | Dikarya | alpha-1,3-glucan synthase |
| 12180 | cfmA | Dikarya | GPI-anchored CFEM-domain protein |

Note on Data 3 row for cluster 4166: the two N. crassa genes are listed under JGI IDs 3834_Neucr2 (bgt2_eglc) and 5869_Neucr2 (bgt1). The file does not state the NCU numbers. I have not mapped JGI IDs to NCU06381 and NCU09175 myself, so that mapping is "unverified". The A. fumigatus ortholog column gives Afu1g11460 (named bgt1 there) for bgt2, and NA for bgt1.

Adhesion category (45 rows, 24 distinct clusters in the file [corrected 2026-10-17, see 06 section 6 row 5: 23 distinct cluster numbers, equal to the paper's "23 clusters"; 24 was an off-by-one]; full list in Data 3):

| Cluster | Genes | Emergence node |
|---|---|---|
| 220 | hum2, sc1, sc3, sc6, coh1, coh2 (class I hydrophobins) | Dikarya |
| 417 | glfB (UDP-galactofuranose transporter) | Opisthokonta + T. trahens |
| 1012, 2024, 7971, 48004 | ncl1-ncl4 (mannose or GlcNAc binding lectins) | various |
| 4291 | iff4, hyr1 (GPI-modified) | Saccharomycotina |
| 5266 | aga1 (adhesin) | Opisthokonta |
| 5817 | hum3 | Opisthokonta |
| 6453 | cfl1 (hypha-specific matricellular protein) | Basidiomycota |
| 7118 | cryp, hfb2, hyd5, NCU08192, bhp2, bhp3 (class II hydrophobins); hfb1 is labelled "galectin" in the file | Dikarya + Mucoromycota |
| 7187 | pga62_flo1, fig2 | Eukaryotes |
| 7276 | flo1, flo5, flo9, flo10 | S. cerevisiae clade |
| 8398 | rep1 (repellents) | Opisthokonta + T. trahens |
| 8498 | als3_als1 | Saccharomycotina + Pezizomycotina |
| 8594, 16358, 14985 | wap1, rbt5, rbt1, hwp1, eap1, flo11, bsc1 | Candida clade or Saccharomycotina |
| 16148 | bad1, rodA, rodB (class I hydrophobins) | Eurotiomycetes, Leotiomycetes, Sordariomycetes |
| 44817, 50327, 126186, 140788 | dfi1, eas, hyd1_like, ece1 | various |

The file lists the "node of emergence" per cluster as the paper's inferred origin. I did not verify the labels against the trees. Where the file shows hfb1 as "galectin" in a hydrophobin cluster I report it as written. It is not verified.

Further rows in Data 3 that concern fusion or surface contact, which fall outside the cell wall category:
- septation: `so_soft` (cluster 13739, Pezizomycotina; "promote cell-cell fusion, reversibly associate with pores"), Woronin-body genes (hex, wsc, lah, spa1), septins (aspA/C), peroxisomal and RAM-network genes.
- signaling: pheromone receptors ste3/preA (cluster 668), ste2/preB (cluster 9439), mpkB, ubc4, ptc1/3/ptcG.
- transcription: sfl1, flo8, adn3 (cluster 8029, "adhesion, transcriptional regulator of flocculation"), snf5, hgrA (sebA).

**Supplementary Data 8 (414 families).** Columns: cluster number, p-value, F-value, A. fumigatus ortholog and function, S. cerevisiae ortholog and function, C. albicans ortholog, number of unicellular fungi with the family (UC), number of multicellular fungi (MC). I searched the 414 rows for surface and fusion keywords. Hits that matter here (cluster / ortholog / annotation as in the file; UC = unicellular fungi with the family, MC = multicellular):

| Cluster | Ortholog | File annotation | UC / MC |
|---|---|---|---|
| 435 | gel2 (Afu6g11390) | GPI-anchored 1,3-beta-glucanosyltransferase | 0 / 46 |
| 374 | crh3 | chitin transglycosylase | 2 / 51 |
| 179 | eng2 (Afu2g14360) | GPI-anchored endo-beta-1,3-glucanase | 1 / 38 |
| 2065 | bgt3 (Afu5g08780) | putative cell wall glucanase | 0 / 41 |
| 635 | fks1 | beta-1,3-glucan synthase | 0 / 51 |
| 209, 234, 368 | chsB, chsE, chsF | chitin synthases | 3 to 4 / 43 to 53 |
| 636, 109 | chi4, chiB1 | chitinases | |
| 171 | exg1 | exo-1,3-beta-glucanase | 0 / 51 |
| 3359 | Afu8g01730 | glycosyl-transfer domain | 0 / 38 |
| 3660 | Afu7g05340 | "conserved protein predicted adhesin; expression induced in biofilm" | 0 / 40 |
| 4240 | C. albicans C2_08770C | "putative adhesin-like protein" (only a C. albicans ortholog listed) | 0 / 38 |
| 3355 | Afu4g11210 | ortholog of S. cerevisiae PRM1 (pheromone-regulated membrane fusion protein) | 2 / 48 |
| 2511 | Afu4g08570 | "role in conjugation with cellular fusion" (positive regulation) | 0 / 41 |
| 1764 | Afu7g02350 | role in conidium formation, syncytium formation | 3 / 49 |
| 84, 1204 | gpaA, gpgA (G-protein alpha and gamma, pheromone response) | | |
| 1348 | Afu5g12550 | GPI12, GPI-anchor assembly (second step) | 3 / 52 |
| 338 | Afu6g14250 | EMP70 ortholog, "cellular adhesion and filamentous growth" | 3 / 52 |

So the screen finds GPI anchor synthesis (GPI12), several GPI-anchored glucan remodelers, a predicted adhesin, and PRM1-like fusion genes. I read these off the file. Their significance is mine: the paper text does not discuss them individually, except Bgt3 and Gel2 by name.

### 2.4 Downloadable lists (all links checked, HTTP 200)

Base: `https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12085-w/MediaObjects/`

| Supplementary Data (paper numbering) | File | Content |
|---|---|---|
| (descriptions) | `41467_2019_12085_MOESM4_ESM.pdf` | file list |
| 1 | `..._MOESM5_ESM.xlsx` | 72 genomes |
| 2 | `..._MOESM6_ESM.xlsx` | copy numbers of hypha families per species |
| 3 | `..._MOESM7_ESM.xlsx` | 651 hypha genes |
| 4 | `..._MOESM8_ESM.xlsx` | families with BCZ duplications |
| 5 | `..._MOESM9_ESM.xlsx` | Fisher tests |
| 6 | `..._MOESM10_ESM.xlsx` | Welch t-tests |
| 7 | `..._MOESM11_ESM.xlsx` | domain architecture changes |
| 8 | `..._MOESM12_ESM.xlsx` | 414 ANOVA families |
| 9 | `..._MOESM13_ESM.xlsx` | yeast loss tests |
| 10-17 | `..._MOESM14_ESM.txt` to `..._MOESM21_ESM.zip` | MCL clusters, alignment, trees, scripts |

Naming: Data N is MOESM(N+4). I verified the mapping for Data 3, 4 and 8 from file contents. Data 1, 2, 5, 6, 7 and 9 I only matched through the description PDF. The adhesion gene set (45 genes) is a sheet category inside Data 3, not a separate file. The 1 to 17 numbering in the PDF is in the order listed, and MOESM1 to MOESM3 are the main Supplementary Information, peer review file and reporting summary (inferred from the Europe PMC supplement list, not opened). Paper also says Data 10 contains the MCL clusters; I did not open it.

For PMID 26541633: no supplementary data or gene list was accessed. The only gene IDs are the two in the title (NCU06381, NCU09175).

---

## 3. Overlap analysis

Labels: **[E]** evidence read in the sources, **[I]** my inference.

### 3.1 Cell wall synthesis and remodeling

- [E] Cell wall biogenesis is one of two functional groups with significantly enriched duplications in BCZ nodes. Six of the 17 BCZ-origin families are cell wall genes (crh3, gel7, DFG5, cps1, chsD, kre6).
- [E] The Data 3 cell wall category is 94 genes. It includes chitin synthases (classes I to VI), beta-1,3 and alpha-1,3 glucan synthases, PMT O-mannosyltransferases, chitinases, Gas/Gel glucanosyltransferases, Crh transglycosylases, Dfg5/Dcw1, Sun proteins and BGT1/2.
- [E] Cell wall and transcription genes are the two groups without significant gene-structure change.
- [E, PMID 26541633] BGT-1 and BGT-2 are present at hyphal apical dome plasma membrane, septa and hyphal fusion sites. Their deletion alters CW and CR sensitivity.
- [I] The wall-remodeling families that duplicated at hyphal origin (Gas/Gel, Crh, Dfg5, chitin synthase classes) are the same enzyme classes that Riquelme's tip model places at the apical dome. This is consistent but the two papers do not test each other.

### 3.2 Cell-cell adhesion

- [E] Adhesion families did not expand at the origin of hyphae. They contracted from 17 to 14 copies, and expanded later in Agaricomycotina and Saccharomycotina.
- [E] The authors say adhesion "might not be key in vegetative hyphae" and propose that communication may run through volatiles such as farnesol, not canonical receptors.
- [E] Genes lost more than expected in yeast include adhesion-related ones.
- [I] The paper therefore gives little support for "adhesin" gene families driving hyphal formation. The hyphal origin signal sits in cell wall enzymes and polarity genes.
- [I] Class 1 and class 2 hydrophobins and Cfl1 are the adhesion-related families that do track complex multicellularity (fruiting bodies), not hyphae.

### 3.3 Cell fusion

- [E] The curated set contains so_soft (SO/SOFT), pheromone receptors, MAPK components and Woronin-body genes. The 414-family screen contains PRM1-like (Afu4g11210) and "conjugation with cellular fusion" (Afu4g08570) families.
- [E, PMID 26541633] BGT-1-GFP and BGT-2-GFP localize at germling and hyphal fusion sites.
- [I] Cell fusion is not a focus of Kiss et al. Fusion genes appear only as entries in the lists, and the main text does not analyse them as a class. The paper gives no test of their evolutionary pattern.

### 3.4 GPI-anchored or secreted surface proteins vs intracellular

Annotation of GPI anchoring below is from the Data 3 and Data 8 text fields (curated by the authors), not from my own prediction. I did not run any GPI predictor.

| Class | Families | Basis |
|---|---|---|
| Stated GPI-anchored or cell-wall-resident | Gas/Gel/Phr (cluster 435), Dfg5, Ecm33, csf4 (crh cluster), cht2, chiA1, cfmA, eng2, Sun/Pir (cell wall proteins), hydrophobins, Flo/Als/Hwp/Eap/Rbt, Hyr1/Iff4, Cfl1, Rep1, lectins (ncl) | [E] file annotations |
| Membrane-embedded enzyme (not GPI) | chitin synthases, Fks1/Bgs, Ags, Kre6, PMTs, Cps1 | [E] function; membrane location [I] |
| Intracellular or cytosolic | Pan1, Sla1, Slm1, Bem1, Bnr1, Whi2, RicA, exocyst, SNAREs, septins, Woronin-body proteins, transcription factors | [E] function |
| Plasma membrane receptors and sensors | pheromone receptors Ste2/Ste3, Ssy1, Sho1, Gpr1 | [E] |
| Secreted | Kex2/KexB protease, chitinases (some) | [I] |

### 3.5 What a surface-glycoprotein detector would likely count as "adhesion"

Context: the repo's own review concluded the classifier behaves as a cell-surface glycoprotein detector with about 12% adhesin precision (AGENTS.md; memory note). This section is inference throughout [I], since I did not run the model on any of these genes.

- **Likely to score as adhesion-like.** GPI-anchored, Ser/Thr-rich or glycosylated wall proteins: Gas/Gel/Phr family members, Dfg5, Ecm33, Sun and Pir, CFEM proteins (cfmA), Cfl1, Flo/Als/Hwp-type adhesins, hydrophobins (small, cysteine-rich, secreted; composition may differ from adhesins), and BGT-1/BGT-2 (predicted GPI-anchored). The Kiss et al. "adhesion" category itself contains cell wall proteins that are not adhesins in a strict sense.
- **Unlikely.** Intracellular polarity, actin, septation, vesicle and transcription families, and multi-pass membrane enzymes (chitin and glucan synthases) without a glycosylated ectodomain.
- **Consequence for the sorting hat.** [I] Cell wall remodelers such as Gas/Gel and BGT would count as false positives for "cell_wall_adhesion_candidate" unless a Pfam or functional rule separates enzymes from adhesins. The enzymes are the part that tracks hyphal origin. The adhesin families do not.
- Pfam for these classes: I did not look them up. For example, Gel/Gas is glycoside hydrolase family 72 and BGT-1/2 may fall in GH17, but I did not verify either, and no Pfam ID is given in the sources: "unverified".

---

## 4. Follow-up reading (PMIDs resolved against PubMed metadata)

| PMID | Title | Journal, year | Why |
|---|---|---|---|
| 30902897 | Transcriptomic atlas of mushroom development reveals conserved genes behind complex multicellularity in fungi (Krizsan et al.) | PNAS 2019 | Cited by Kiss et al. for kinase diversity in complex multicellular Basidiomycota; expression side of complex multicellularity |
| 32191325 | Unmatched Level of Molecular Convergence among Deeply Divergent Complex Multicellular Fungi (Merenyi et al.) | Mol Biol Evol 2020 | Same group; convergent gene families in fruiting-body fungi, relevant to the hydrophobin and Cfl1 expansion |
| 32743136 | Off the wall: The rhyme and reason of hyphal morphogenesis (Verdin et al.) | Cell Surf 2019 | Review of cell wall biosynthesis and the hyphal tip, from Riquelme's group; PubMed shows no PMC id (access not checked) |
| 36786575 | AoSte12 Is Required for Mycelial Development, Conidiation, Trap Morphogenesis, and Secondary Metabolism by Regulating Hyphal Fusion in Nematode-Trapping Fungus (Bai et al.) | Microbiol Spectr 2023 | Peripheral: hyphal fusion tied to adhesive trap formation |

I read only the abstracts of these. I did not verify whether 30902897 and 32191325 provide gene lists. A discarded search hit (PMID 33634751, autophagy assay guidelines) was irrelevant and is not used.

Not found in PubMed with my queries, so not cited: a review of N. crassa cell fusion genes; a GH72 family review. I did not cite them from memory.

---

## 5. Open items

1. Full text of PMID 26541633 (publisher 403, no PMC). Needed for GH family, signal-peptide and GPI-site evidence, and Pfam or domain information.
2. Map Data 3 JGI IDs (3834_Neucr2, 5869_Neucr2) to NCU numbers.
3. Open Data 2 and Data 4 to get species-level copy numbers for the cell wall families, if the sorting hat wants clade patterns.
4. If wanted, run the shipped classifier on the cell wall and adhesion categories of Data 3 (A. fumigatus and N. crassa orthologs) to test section 3.5 empirically. Not done here.

## Errata log

Applied 2026-10-17 from reports 05, 06, 07. Line numbers are those of the file before this errata (the added header line shifts later lines by 2).

1. Line 46. Old: `mostly yeast and Candida.`. New: same text with an added marker. Source: 06 s6 #17.
2. Line 149. Old: `The 9-group totals in the text are probably the file minus "adhesion". This is an inference; I did not reconci`. New: same text with an added marker. Source: 06 s3, s6 #6.
3. Line 158. Old: `| Opisthokonta + amoebozoa |`. New: same text with an added marker. Source: 06 s6 #13.
4. Line 162. Old: `"CAZy family 17" for Gel4 as stated in the file |`. New: same text with an added marker. Source: 05 Summary item 4, Table B GH17 row.
5. Line 175. Old: `45 rows, 24 distinct clusters in the file`. New: same text with an added marker. Source: 06 s6 #5.

Items still unverified (not changed): Gel4 family (conflict marked above); cluster 4166 JGI IDs and source citation (06 row 15); cluster 3359 (06 row 20); 519 publications (06 row 1).
