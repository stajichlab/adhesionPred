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
