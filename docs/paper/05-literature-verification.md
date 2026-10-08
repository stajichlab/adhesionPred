# 05. Verification of literature claims written from memory

*2026-10-06. Checked against PubMed (abstracts), the NCBI taxonomy dump in `/srv/projects/db/taxonomy`, and the
UniProt proteome records (live queries on this date). Only abstracts were read. No full text.*

The Plan 2 review notes (`docs/superpowers/plans/2026-10-05-cellsurface-sorting-hat-modules-and-calibration-review-1.md`)
name these papers but do not record the sentence each one supported. This file records what each source says.
Before a paper cites one of them, compare the sentence you want to write with the "What the source shows" column.

## 1. Citations

| Cited as | Resolved to | What the source shows (from the abstract) | Status |
|---|---|---|---|
| Vaknin 2014 | Vaknin Y, ..., Osherov N. *Fungal Genet Biol* 63:55-64 (2014; online 2013-12-20). PMID 24361821, doi 10.1016/j.fgb.2013.12.005 | *A. fumigatus* has three CFEM-domain GPI-anchored proteins (CfmA-C). Deletion mutants are more sensitive to Congo Red and Calcofluor White. No growth or germination defect. No role in heme uptake or biofilm. No change in virulence in insect or mouse models. | Exists. Year is 2014 in the volume, 2013 online. |
| Liu 2016 | Liu H, ..., Filler SG. *Nat Microbiol* 2:16211 (2016). PMID 27841851, doi 10.1038/nmicrobiol.2016.211 | The thaumatin-like protein CalA is on the *A. fumigatus* cell surface, binds integrin alpha5beta1 and mediates epithelial and endothelial invasion. The deletion mutant has reduced virulence in mice. | Exists. This is my match by topic (CalA appears in the project). The review notes do not name the title. |
| Balajee 2007 | Balajee SA, ..., Rooney AP. *Eukaryot Cell* 6:1392-9 (2007). PMID 17557880, doi 10.1128/EC.00164-07 | CSP typing: sequence of a putative cell surface protein gene (Afu3g08990) with tandem repeats and point mutations separates outbreak isolates. Af293 has a CSP gene structure "substantially different" from the other isolates. | Exists. Topic match is my inference. Other 2007-era Balajee papers exist, so confirm this is the one meant. |
| Fedorova 2008 | Fedorova ND, ..., Nierman WC. *PLoS Genet* 4:e1000046 (2008). PMID 18404212, doi 10.1371/journal.pgen.1000046 | Genome of the clinical isolate A1163 and comparison with Af293. Core genes are 99.8% identical. Up to 2% of genes are unique to each isolate. Variable genes can be as low as 40% identical. Strain-specific genes cluster in subtelomeric islands. | Exists. Supports "annotation and gene content differ between Af293 and A1163". |
| Gravelat 2013 | Gravelat FN, ..., Sheppard DC. *PLoS Pathog* 9:e1003575 (2013). PMID 23990787, doi 10.1371/journal.ppat.1003575 | Galactosaminogalactan (a polysaccharide, made via uge3) is the dominant adhesin of *A. fumigatus*. It mediates adherence to plastic, fibronectin and epithelial cells and masks beta-glucan. | Exists. **It is a polysaccharide adhesin. A protein-based adhesion predictor cannot detect it.** Do not cite it as a protein adhesin. |

## 2. Other claims

| Claim | Result | Evidence |
|---|---|---|
| A1163 is a CEA10 derivative (CBS 144.89 = FGSC A1163 = CEA10) | **Correct for the UniProt and NCBI naming.** | UniProt UP000001699: "strain CBS 144.89 / FGSC A1163 / CEA10", taxon 451804, 9,942 proteins, assembly GCA_000150145.1. NCBI has separate strain nodes for A1163 (451804) and CEA10 (505235). Both have the same parent (746128). Whether A1163 and CEA10 are the same isolate in the lab sense is not tested here. |
| Taxon 330879 = Af293 | **Correct.** | NCBI: "Aspergillus fumigatus Af293", strain, parent 746128. UniProt UP000002530: taxon 330879, 9,647 proteins. |
| Taxon 451804 = A1163 | **Correct.** | NCBI: "Aspergillus fumigatus A1163", strain, parent 746128. |
| Taxon 746128 = *A. fumigatus* (species) | **Correct.** | NCBI: species, scientific name "Aspergillus fumigatus", parent 2720872. W72310 has no strain node. Its taxon is the species ID. |
| gp43 (Paracoccidioides) binds laminin and fibronectin | **Supported.** | PMID 16698299 (Mendes-Giannini 2006, *Microbes Infect* 8:1550-9): "gp43 bound both fibronectin and laminin" in ligand affinity assays. PMID 10962270 (Hanna 2000, *Microbes Infect* 2:877-84): anti-gp43 serum abolished 85% of binding to Vero cells; gp43-deficient mutant 113M was less adherent. gp43 is also named a laminin receptor in that abstract. |
| gp43 is a moonlighting glucanase on the cell surface | **Not checked.** | The abstracts above say nothing about glucanase activity. The UniProt entry Q9HDL9 is annotated as a glucan 1,3-beta-glucosidase. |

## 3. Corrections for the other documents

- `docs/paper/04`, sections 2 and 4: updated to point to this file.
- The A1163 UniProt count of 9,942 is correct as of 2026-10-06. The Plan 2 plausibility range 9,000 to 10,500 holds.
- Still unverified: the YPS3 accession (A0ACF1AYZ2), the BcLysM1 adhesion sentence, the three cell wall reviews
  (Gow 2017, Riquelme 2020, Gow 2023), and the WHO/IUIS counts.

## 4. Further checks, 2026-10-07

| Item | Result | Evidence |
|---|---|---|
| YPS3 accession `A0ACF1AYZ2` (task 08, agent-C) | **Matches by name.** UniProt TrEMBL entry A0ACF1AYZ2: "Yeast-phase specific protein yps-3", gene YPS3, ORF `I7I48_01120`, *Histoplasma ohiense* (taxon 2902605), 137 aa. Not checked: that the sequence is the YPS3 of the 2005 paper, or that the strain is G217B. | UniProt REST, 2026-10-07 |
| BcLysM1 adhesion sentence (PMID 39655398) | **Supported by the abstract.** "contribution of BcLysM1 in infection initiation and in adhesion to bean leaf surfaces were demonstrated" and "a dual role in mycelial adhesion and suppression of chitin-triggered host immunity". The words "hydrophobic surfaces" are not in the abstract. | PubMed abstract, [doi](https://doi.org/10.1002/jobm.202400552) |
| WHO/IUIS counts (120 molecules, 31 species) | **Internally consistent.** The stored table `analysis/allergen_scoping/iuis_fungal_allergens.tsv` has 120 rows and 31 species. The isoallergen table has 116 rows with a protein sequence (111 written to the FASTA, 5 skipped). allergen.org was not queried again. | recount, 2026-10-07 |
| Table B3 of `docs/paper/02` | **Matches** `metrics.json` (24 values). | recomputed, 2026-10-07 |
| RBT1 CGD statements | **Verified** in CGD's own notes, with PMIDs 10978273, 19837954, 21414038. Braun 2000 full text: the two RBT1 alleles differ by a segment (aa 612 to 640), not by repeat number. | CGD API, `to_import/genetics0031.pdf` |

Still not verified: the three cell wall reviews (Gow 2017, Riquelme 2020, Gow 2023), the other "copied" rows of `docs/paper/02` (the A rows, B1 to B2, B4 to B6, C and D rows), the strain identity of YPS3, and the allergen.org live counts.

## 5. The 8-cysteine hydrophobin pattern: prior art (2026-10-08)

Searched PubMed only. Abstracts only. No full text. Not a complete novelty search.

| Source | What the abstract says | Bearing on the rescue rule |
|---|---|---|
| Yang 2006, *BMC Bioinformatics* 7(Suppl 4):S16. PMID 17217508, doi 10.1186/1471-2105-7-S4-S16 | Primary-structure analysis of hydrophobins: BLAST, MEME motifs, and MAST search of nr with the motifs and the "C-CC-C-C-CC-C" pattern. 9 new candidates after filtering by pattern, domain and length. | A pattern-based search for new hydrophobins exists. Read the full text before any novelty claim. |
| Kubicek 2008, *BMC Evol Biol* 8:4. PMID 18186925 | Hydrophobins have eight conserved cysteines. Class I and II are separated by hydropathy and solubility. Class II was found only in ascomycetes. *Trichoderma* has up to 10 class II genes. | Defines the family and the two classes. |
| Seidl-Seiboth 2011, *J Mol Evol* 72:339-51. PMID 21424760 | *Trichoderma* hydrophobins that deviate from the two-class scheme in hydropathy, cysteine spacing and surface pattern; they form separate clades inside ascomycete class I. | Cysteine spacing varies. A fixed spacing pattern can miss such proteins. |
| Xu 2021, *Microbiol Res* 247:126723. PMID 33636611 | 40 hydrophobin genes in *P. ostreatus*; all contain eight cysteines with a conserved spacing pattern; 33 are class I. | The pattern is used to identify family members at genome scale. |
| De Vries 1999, *Eur J Biochem* 262:377-85. PMID 10336622 | CFTH1 from *Claviceps fusiformis* has three class II hydrophobin domains in one protein, each preceded by a Gly/Asn-rich region. | A length cap on the rescue rule would miss this protein. |
| Peñas 1998, *Appl Environ Microbiol* 64:4028-34. PMID 9758836 | Hydrophobins are small (about 100 +/- 25 residues), cysteine-rich proteins in the cell wall rodlet layer. | Source of the usual length range. |
| Pitocchi 2026, *Int J Biol Macromol* 378:153896. PMID 42546942, doi 10.1016/j.ijbiomac.2026.153896 | PAC3, an 83-residue surface-active protein of *Acremonium sclerotigenum*, lacks the eight-cysteine motif. The authors propose a new family of fungal protein biosurfactants. | A surface-active protein can exist without the motif. The rule has a stated limit. |

Not read: Wessels 1994, Linder 2005, Sunde 2008 (the sources for class I and II cysteine spacing). Task H3 of
`docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md` reads them.

## 6. Cysteine spacing and wider prior art for the 8-cysteine rule (2026-10-08, tasks H3 and H3b)

Rules for this section. No truth protein, no sequence analysis and no regex were used. No spacing is from memory.
Every quote below was copied from text that I read in this task. Unit of every spacing number: residues between
two cysteines. The numbers are also in `data/sorting_hat/cys8_spacing.yaml`. Sources that disagree are listed
side by side. They are not merged.

### 6.1 Read status of each source

| Source | PMID / doi | Status |
|---|---|---|
| Wessels 1994, *Annu Rev Phytopathol* 32:413-437 | no PMID in PubMed; doi 10.1146/annurev.py.32.090194.002213 | NOT READ. Not in PubMed (author and title queries gave no hit). Annual Reviews page returned 403. Citing papers attribute numbers to it. Those are not used as read. |
| Linder 2005, *FEMS Microbiol Rev* 29:877-96 | PMID 16219510, doi 10.1016/j.femsre.2005.01.004 | NOT READ. Abstract only. No PMC copy. Publisher, OUP and doi pages returned 403. VTT page has the abstract only. |
| Sunde 2008, *Micron* 39:773-84 (PubMed title "Structural analysis of hydrophobins") | PMID 17875392, doi 10.1016/j.micron.2007.08.003 | NOT READ. Abstract only. No PMC copy. |
| Kubicek 2008, *BMC Evol Biol* 8:4 | PMID 18186925, PMC2253510, doi 10.1186/1471-2148-8-4 | READ IN FULL. The PubMed tool text drops the subscripts. The consensus in Figure 2 is a picture. I rendered page 5 of the BMC PDF and read the consensus from it. |
| Seidl-Seiboth 2011, *J Mol Evol* 72:339-51 | PMID 21424760, doi 10.1007/s00239-011-9438-3 | READ IN FULL. No PMCID. The Springer page shows the full text and Table 3 (page `.../tables/3`). |
| Yang 2006, *BMC Bioinformatics* 7(Suppl 4):S16 | PMID 17217508, PMC1780129, doi 10.1186/1471-2105-7-S4-S16 | READ IN FULL (Europe PMC XML). It gives no spacing numbers. |
| Peñas 1998, *Appl Environ Microbiol* 64:4028-34 | PMID 9758836, PMC106595, doi 10.1128/aem.64.10.4028-4034.1998 | READ IN FULL. The PubMed tool returned an empty full text. I read the PMC web page. |
| Mgbeahuruike 2013, *BMC Evol Biol* 13:240 (found by citation chain) | PMID 24188142, PMC3879219, doi 10.1186/1471-2148-13-240 | READ IN FULL. |
| Jensen 2010, *BMC Res Notes* 3:344 (found by search) | PMID 21182770, PMC3020181, doi 10.1186/1756-0500-3-344 | READ IN FULL. |
| Li 2021, *Int J Mol Sci* 22:643 (found by search) | PMID 33440688, PMC7827705, doi 10.3390/ijms22020643 | READ IN FULL. |

None of Wessels 1994, Linder 2005 and Sunde 2008 could be read. Seven other sources with spacing numbers were read
in full. So the stop rule of task H3 did not apply.

### 6.2 Spacing statements, verbatim

| Source and place | Quote |
|---|---|
| Seidl-Seiboth 2011, Abstract | "Hydrophobins are conventionally grouped into two classes (class I and II) according to their solubility in solvents, hydropathy profiles and spacing between the conserved cysteines. Here we describe a novel set of hydrophobins from Trichoderma spp. that deviate from this classification in their hydropathy, cysteine spacing and protein surface pattern." |
| Seidl-Seiboth 2011, Table 3 (title and footnote) | "Spacing between conserved cysteines in class I, class II hydrophobins and the HFBs from this study". Footnote: "Numbers indicate the numbers of amino acids occurring between the two cysteine residues (consecutively numbered from N to C terminus), and were taken from all ascomycete HFBs investigated in this study (Tables 2, 4). Values in italics specify minimum and maximum numbers" |
| Seidl-Seiboth 2011, Table 3 values (columns 1_2, 3_4, 4_5, 5_6, 7_8) | Row I: 6, 33-39, 19-25, 5, 15-17. Row I new: 6-7, 9-13, 5-10, 5, 8-12. Row II: 10, 11, 16, 8, 10. The page prints the header as "Class, Cysteine-#, 1_2, 3_4, 4_5, 5_6, 7_8". I read the five numbers after the class as the five columns. Row II agrees with Mgbeahuruike 2013 and with the Kubicek consensus (6.3), which supports that reading. |
| Seidl-Seiboth 2011, Results | "A characteristic feature of all detected Trichoderma class I HFBs hydrophobins was an amino acid stretch of at least 4 N residues that immediately followed the C6/C7 pair." and "The second clade, termed clade B, did not exhibit these features, and also contained 1-2 aa less in the spacers between the 8 Cs." and "a 45-65 aa core structure containing the eight Cs" |
| Seidl-Seiboth 2011, Discussion | "their other properties (hydropathy plots, cysteine spacing) are clearly different from the known class I and II HFBs and suggest that the concept of grouping into class I and II needs to be expanded." |
| Seidl-Seiboth 2011, Methods | "We used BLASTP, Psi-BLAST and HMM profiles (PF01185, IPR 001338) to identify ascomycetous class I HFBs in other fungal genera" |
| Kubicek 2008, Background | "In the primary sequence, the most important feature common to all hydrophobins is the characteristic pattern of eight Cys-residues, which gives rise to a common disulfide network" |
| Kubicek 2008, Results (protein structure) | "most of the predicted HFBs had the expected structure of 90 - 110 amino acids, which includes a 15-20 aa signal peptide, the 65 aa core structure displaying the eight cysteines which are predicted to have four 4 beta-strands and a single helix." |
| Kubicek 2008, Methods | "class I hydrophobins (identified according to the criteria described by Linder et al. []; e.g. by the difference in the number of amino acids between the conserved cysteins and their hydropathy profile) removed." |
| Kubicek 2008, Figure 2 (consensus line, read from the page image) | C1 P X G L X(3-4) P Q C2 C3 X(3) V L G V al X L D C4 X(2) P X(9) F X(3) C5 X(3) G X(4) C6 C7 V V P I X(4) al L C8. Caption: "X denotes any amino acid, and the subscript the number of them; "al" denotes any aliphatic, hydrophobic amino acid". Gaps by counting this string by hand: C1-C2 9 to 10, C2-C3 0, C3-C4 11, C4-C5 16, C5-C6 8, C6-C7 0, C7-C8 10. This is the class II consensus for Trichoderma/Hypocrea only. |
| Mgbeahuruike 2013, Results | "A comparison was made between the aligned sequences and already published sequence consensus of class I, C-X5-7-C-C-X19-39-C-X-8-23-C-X5-C-C-X6-18-C-X2-13 [33] and class II, C-X9-C-C-X11-C-X14-16-C-X8-C-C-X10-C-X6-7 [14] hydrophobins." (ref 33 is Kershaw and Talbot 1998, ref 14 is Kubicek 2008. Both not read here as primary sources.) |
| Mgbeahuruike 2013, Table 1 (C1/C2, C3/C4, C4/C5, C5/C6, C7/C8) | Class I (basidiomycetes): 6, 26-33, 12-13, 6, 13. Class I (ascomycetes): 6-7, 26-39, 18-21, 6-8, 15-17. Class II: 9-10, 11, 15-16, 2-7, 10. T. terrestris protein 159967: 7, 5, 8, 5, 12. U. maydis protein 5010: 6, 49, 17, 5, 16. Footnote: "Cysteine residues C2/C3 and C6/C7 are adjacent in all characterized hydrophobins." |
| Mgbeahuruike 2013, Results | "The hydrophobins from the thermophilic fungus Thielavia terrestris and the corn smut fungus Ustilago maydis deviated from the remaining analyzed hydrophobins in the length of the region between cysteine residues C3 and C4" |
| Jensen 2010, Background | "considerable variation is seen in the cysteine spacing of class I hydrophobins, while less variation is seen for class II hydrophobins [7]" (ref 7 is Kershaw and Talbot 1998) |
| Jensen 2010, Results | "Forty-four of the identified hydrophobins displayed class I cysteine spacing pattern, but only twenty-four had a characteristic class I hydropathy plot" and "A common feature in 44 of the 50 hydrophobins is a conserved spacing of five amino acids between the fifth and sixth cysteines, while the remaining six hydrophobins contain either seven or eight amino acids." |
| Jensen 2010, Results (A. terreus) | "ATEG_06492 displayed a characteristic class I hydrophobin cysteine spacing pattern (CN{7}CCN{40}CN{16}CN{5}CCN{17}C), whereas a class II hydrophobin spacing pattern was observed for ATEG_04730 (CN{10}CCN{11}CN{16}CN{8}CCN{10}C)." |
| Jensen 2010, Results | "They have a similar cysteine pattern of CN{5-13}CCN{17}CN{7-12}CN{7}CCN{8-12}C (where N signifies any other amino acid than cysteine) ... patterns that differ from both class I and class II hydrophobins and can therefore theoretically not be placed in either class." |
| Jensen 2010, Abstract | "twenty-six of the identified hydrophobins were intermediate forms." |
| Li 2021, Results | "CmHYD1 contained Pfam06766 ... It had a cysteine pattern of CX9-CCX8-CX19-CX8-CCX10-C" and "Four hydrophobin units had the same cysteine pattern of CX9-CCX11-CX16-CX8-CCX10-C, which was a little different from CmHYD1 and CmHYD2, but corresponded to the consensus defined for the fungal class II hydrophobins" and "the cysteine pattern was displayed as CX7-CCX37 or 25-CX17-CX5-CCX10-C" (class I, CmHYD3 and CmHYD4). |
| Li 2021, Methods 4.4 | "All protein sequences bearing the hydrophobin type I (C-X5-7-C-C-X19-39-C-X8-23-C-X5-C-C-X6-18-C-X2-13) or type II (C-X9-CC-X11-C-X14-16-C-X8-C-C-X10-C-X6-7) signature sequences were retrieved [37]." (ref 37 is Mgbeahuruike 2013.) |
| Peñas 1998, Results | "The predicted protein contains two clusters of cysteines spaced like the clusters described previously for hydrophobins (33) (C-X6-CC-X31-C, C-X5-CC-X12-C)." The paper does not say which pattern is class I or II. |
| Peñas 1998, Introduction and Discussion | "characterized by the conserved pattern of spacing of the eight cysteine residues present in their sequences". "Class I hydrophobins form very stable complexes that are insoluble in SDS and contain cysteine doublets followed by stretches of hydrophilic amino acids. Class II hydrophobins are soluble in SDS, and their cysteine doublets are immediately followed by hydrophobic residues." |
| Yang 2006, Abstract and Results | "Based on the newly found motifs and the well-known C-CC-C-C-CC-C pattern we used MAST to search the entire nr database." and "the sequences that don't have the eight-cysteine residues pattern (C-CC-C-C-CC-C) by the perl program". No numbers of residues between cysteines are given. |

### 6.3 Where the sources disagree

All differences are in the table in 6.2 and in the yaml. Listed gap by gap, first for class II.

| Gap | Seidl-Seiboth 2011 | Mgbeahuruike 2013 | Kubicek 2008 (counted) | Li 2021 signature | Jensen 2010 example |
|---|---|---|---|---|---|
| C1-C2 | 10 | 9-10 | 9-10 | 9 | 10 |
| C3-C4 | 11 | 11 | 11 | 11 | 11 |
| C4-C5 | 16 | 15-16 | 16 | 14-16 | 16 |
| C5-C6 | 8 | 2-7 | 8 | 8 | 8 |
| C7-C8 | 10 | 10 | 10 | 10 | 10 |

The one clear conflict is class II C5-C6. Mgbeahuruike Table 1 gives 2-7. All four other places give 8. I did not check which is right. The Mgbeahuruike text does not explain the 2-7.

Class I (ascomycetes unless stated):

| Gap | Seidl-Seiboth 2011 | Mgbeahuruike 2013 asco | Mgbeahuruike 2013 basidio | Li 2021 signature |
|---|---|---|---|---|
| C1-C2 | 6 | 6-7 | 6 | 5-7 |
| C3-C4 | 33-39 | 26-39 | 26-33 | 19-39 |
| C4-C5 | 19-25 | 18-21 | 12-13 | 8-23 |
| C5-C6 | 5 | 6-8 | 6 | 5 |
| C7-C8 | 15-17 | 15-17 | 13 | 6-18 |

Class I ranges differ between sources and between basidiomycetes and ascomycetes. Two sources say class I spacing varies more than class II spacing: Jensen 2010 (Background) and the Lovett 2022 preprint (Box 1: "The spacing between the cysteines is highly variable." for class I, "The spacing between the cysteines is highly conserved." for class II). Three sources say proteins exist that fit neither class: Seidl-Seiboth 2011, Jensen 2010 (the four unclassified proteins and "26 intermediate forms") and Mgbeahuruike 2013 (two proteins).

Yang 2006 compared with: BLAST of 183 UniProt "hydrophobin" entries against nr, MEME motifs, MAST search, a perl filter for the eight-cysteine pattern "C-CC-C-C-CC-C", and a SMART domain check. It reports that the nine candidates "contain the hydrophobin domain stored in the Pfam database" (Figure 3 caption). It did not compare against Pfam HMM hits as a baseline and did not use HMMs for search. It did not use a signal peptide.

### 6.4 Wider prior art for the rescue rule (task H3b)

Rule under test: add proteins that have the 8-cysteine pattern plus a signal peptide to a Pfam-based hydrophobin call.
Searches done 2026-10-08: 7 PubMed field-tag queries (titles with "hydrophobin" and prediction, HMM, genome, Pfam, cysteine, secretome, classification, evolution; titles of about 140 returned hits checked; 9 papers read in full or in part); two web searches limited to bioRxiv; InterPro API for 20 hydrophobin entries; Europe PMC.
The bioRxiv MCP tool is date-only, so I did not use it.

| Source | What it says | Bearing on the rescue rule |
|---|---|---|
| Jensen 2010, PMID 21182770, doi 10.1186/1756-0500-3-344 (full text) | Genome screen of 9 Aspergillus genomes for hydrophobins by "the criteria of minimum eight cysteines, two cysteine pairs, a size of app. 100 AA and the cysteine pattern". "All identified hydrophobins had theoretical signal sequences". "Forty-five of the identified proteins contained domains classifying them as hydrophobins by Pfam. The remaining five hydrophobins could not be classified." Authors write that BLAST-and-motif methods "may be missed as hydrophobins have high sequence diversity" and criticise Yang 2006 for that. | Closest prior art. It uses the eight-cysteine pattern plus signal sequence to find hydrophobins without a Pfam requirement. It reports 5 of 50 with no Pfam hit. Its data cannot show the false-positive cost, because it did not screen non-hydrophobins. |
| Lovett, Kasson, Gandier 2022, bioRxiv doi 10.1101/2022.08.19.504535 (preprint; read the abstract, Introduction, Box 1, Methods 4.3, Results 5.1-5.2 and Discussion from the PDF; not every section) | Mines 45 proteomes with a "6-Cys pattern" (C2C3-X-C4-X-C5-X-C6C7) "No constraints were imposed on the distance between cysteine residues", plus hmmsearch with PF01185 and PF06766 (Pfam 34.0) and SignalP 5.0. "While only 15 hydrophobin candidates were pulled exclusively by the permissive 6-Cys motif, it is important to note that these sequences would otherwise not have been associated with the hydrophobin family." Also "996" other 6-Cys proteins were called hydrophobin-like and "anticipated to contain a wide range of sequences including those that are functionally and structurally unrelated to hydrophobins". "NRC_EAS was not captured by the Pfam HMMs for hydrophobins, as permissively applied in the study". Some candidates have no secretion signal (Group 10: 48.7% with a predicted signal, 51.3% without). Preprint, not peer reviewed. I did not check for a journal version (the bioRxiv API shows no published version). | Prior art for the combination of a cysteine pattern, Pfam and a signal peptide. It reports that a loose cysteine pattern adds few hydrophobins beyond Pfam (15 in the whole survey) and many other proteins. It also reports hydrophobin candidates without a signal peptide, which a signal-peptide condition would miss. This is a count from one preprint, not a test of our rule. |
| Li 2021, PMID 33440688, doi 10.3390/ijms22020643 (full text) | Finds hydrophobins in ascomycete genomes by BLAST, "All the proteins predicted to contain Pfam06766 and Pfam01185", and "type I ... or type II ... signature sequences", then takes the union; "The classification was checked by Pfam carefully"; SignalP 5.0 for signal peptides. | Prior art for a union of Pfam and a cysteine-spacing signature in a genome survey. It gives no count of proteins found by the signature only. |
| Mgbeahuruike 2013, PMID 24188142, doi 10.1186/1471-2148-13-240 (full text) | Genome survey of two wood-degrading basidiomycetes (335 sequences from 41 species in the tree). Keeps sequences with a hydrophobin domain from InterProScan. Predicts signal peptides with SignalP 3.0. Classes assigned with Kyte-Doolittle plots and the published signatures. | Uses a domain call, then spacing for class. Not a rescue rule. |
| Yang 2006, PMID 17217508 (full text) | See 6.3. Pattern "C-CC-C-C-CC-C" used as a filter after BLAST and MEME/MAST. Every candidate also had to have the SMART hydrophobin domain. | Pattern-based search exists. Their filter needs the domain, the opposite of a rescue. Jensen 2010 states that this approach can miss divergent proteins. |
| Seidl-Seiboth 2011, PMID 21424760 (full text) | Trichoderma hydrophobins with cysteine spacing that deviates from class I and II. Found with BLAST, PSI-BLAST and HMM profiles (PF01185, IPR001338). | The spacing in Table 3 for "I new" differs from class I. A class I spacing range from other sources would not match them (compare 6.3). |
| Bouqellah and Farag 2023, PMID 38004644, PMC10672791 (keyword search of the text, not read in full) | In silico study of 45 class II proteins. "all the proteins carried signal peptides". Uses SignalP 6.0. | Supports a signal peptide in class II. Not a rule test. |
| AlphaFold hydrophobin study 2025, PMID 40852847, doi 10.1002/pro.70279 (abstract and keyword search of text) | AlphaFold-based classification of more than 7,000 UniProt class I and II hydrophobins. Defines "non-canonical hydrophobins with extended disordered N-terminal tails to have at least 70 residues preceding the first Cys". | A length or position cap before C1 would exclude some of them. This agrees with the design decision of no cap. |
| De Vries 1999, PMID 10336622 (from section 5, abstract only) | Three class II domains in one protein. | No cap on length. |
| InterPro entries, read through the InterPro API (https://www.ebi.ac.uk/interpro/api/entry/) | IPR001338 "Class I Hydrophobin" (member databases Pfam and SMART): "This entry represents class I hydrophobins found in fungi." "characterised by eight cysteine residues arranged in a strictly conserved motif". 5,597 proteins. IPR010636 "Class II hydrophobin" (Pfam, CDD, PANTHER): "restricted to ascomycetes" in Pfam PF06766; 1,622 proteins. IPR019778 "Class I Hydrophobin, conserved site": "a cysteine-rich conserved site found at the C-terminal end of class I hydrophobin sequences". Further hydrophobin entries not integrated into those two: PF22354 (Eas), PF28987 and IPR062132 (dewD/F), PF29465 and IPR063347 (hydrophobin-like protein 1), PF29785 (hydrophobin D), PF29802 and IPR062192 (class I hydrophobin F), plus SM00075 and CDD cd23505, cd23507, cd23508, cd23516. IPR049715 (BslB) is bacterial. | The Pfam families that the project already activates (family_table.tsv, Pfam 38.2) include PF22354, PF28987, PF29465, PF29785 and PF29802. The InterPro pages give no spacing numbers. The InterPro text confirms the family is defined by eight cysteines but does not define a spacing. |
| Pfam PF01185 and PF06766 pages (through InterPro) | PF06766: "a family of fungal hydrophobins that seems to be restricted to ascomycetes ... eight cysteine residues arranged in a strictly conserved motif. ... Note that some family members contain multiple copies." PF01185 has no description text in the API. | Same as above. |

Not read: Kershaw and Talbot 1998 (PMID 9501475) and Wessels 1997 (PMID 8922117), which other papers cite for the class I and II spacing; Wessels 1994.

Summary. Pattern-based and signal-peptide-based hydrophobin identification is published (Jensen 2010, Li 2021, Lovett 2022 preprint). A Pfam-plus-pattern union is published too. I found no source that measures sensitivity and false-positive cost of a rescue of the form "8-cysteine spacing pattern plus signal peptide added to a Pfam call" against a truth set. This statement covers only the sources in this section. It is not a complete novelty search. Google Scholar, Scopus and the full text of the unread papers were not searched.
