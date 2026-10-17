# Cell wall biosynthesis and remodeling: gene sets from three reviews

Date: 2026-10-17. Literature research only. No code was changed.

## 0. Important limitation (read first)

Only ONE of the three papers could be read in full text.

| Requested ID | Resolved to (PubMed) | Full text obtained? |
|---|---|---|
| PMID 28429675 | Steinberg G, Penalva MA, Riquelme M, Wosten HA, Harris SD. "Cell Biology of Hyphal Growth." Microbiol Spectr 2017 Apr;5(2). doi 10.1128/microbiolspec.FUNK-0034-2016. PMC11687463 | No. PMC returned the abstract only. |
| PMC5968459 (= PMID 29643171) | Riquelme M, Aguirre J, Bartnicki-Garcia S, Braus GH, Feldbruegge M, Fleig U, Hansberg W, Herrera-Estrella A, Kaemper J, Kueck U, Mourino-Perez RR, Takeshita N, Fischer R. "Fungal Morphogenesis, from the Polarized Growth of Hyphae to Complex Reproduction and Infection Structures." Microbiol Mol Biol Rev 2018 Apr 11;82(2):e00068-17. doi 10.1128/MMBR.00068-17 | No. Publisher blocks XML; PMC web page is behind a browser check; Europe PMC and BioC endpoints failed. Abstract only. |
| PMID 32743136 (PMC7389182) | Verdin J, Sanchez-Leon E, Rico-Ramirez AM, Martinez-Nunez L, Fajardo-Somera RA, Riquelme M. "Off the wall: The rhyme and reason of Neurospora crassa hyphal morphogenesis." Cell Surf 2019 (epub 2019-03-08);5:100020. doi 10.1016/j.tcsw.2019.100020 | Yes (about 53 KB of text). |

Consequences:
- Every gene row below is sourced from PMID 32743136 (full text read). Species names are italic in the source and were stripped by the text extraction. I inferred N. crassa from the NCU locus IDs and the title, and I say so where it matters.
- For PMID 28429675 and PMID 29643171 I can only state what the abstracts say. I did NOT see their gene tables. Any claim about what they list is marked "not seen". The task asked for coverage of S. cerevisiae, C. albicans, A. fumigatus and C. neoformans. Those species are NOT supported by any text I could read, apart from the S. cerevisiae comparisons inside 32743136. The owner should obtain the two reviews as PDFs if species-wide tables are needed.
- PubMed MCP search results for this topic mostly returned unrelated items. Follow-up PMIDs in section 5 were each resolved to a title by get_article_metadata.

## 1. What each review is

**PMID 28429675, Steinberg et al. 2017, Microbiol Spectr (Fungal Kingdom series).** Review of the cell biology of hyphal tip growth in filamentous fungi. Per the abstract: hyphae grow by polarized exocytosis at the apex; motors and the cytoskeleton deliver membranes and proteins; the delivered enzymes include cell wall synthases that are exocytosed for local wall synthesis; endocytosis opposes exocytosis; early endosomes act as spatial organizers; septal pores allow cytoplasmic communication. Focus is trafficking and cytoskeleton, not wall chemistry. Species focus is not stated in the abstract (filamentous fungi in general). Wall gene content: not seen.

**PMC5968459 / PMID 29643171, Riquelme et al. 2018, Microbiol Mol Biol Rev.** Review of fungal morphogenesis from polarized hyphal growth to reproductive and infection structures. Per the abstract: tip elongation needs membrane insertion and extension of the "rigid chitin-containing cell wall"; secretory vesicles carry cell wall-synthesizing enzymes and accumulate in the Spitzenkoerper; cell integrity, branching and morphogenesis depend on vesicle and cytoskeleton dynamics; differentiation of reproduction or interaction structures reprograms the basic machinery. Species focus: filamentous fungi, many model species (12 authors; the abstract does not name species). Wall gene content: not seen. Keywords include "cell wall".

**PMID 32743136, Verdin et al. 2019, Cell Surface.** Review of cell wall composition, assembly, remodeling and the secretory logistics of wall enzymes, centred on Neurospora crassa (NCU IDs throughout; N. crassa named in title). Main claims, as read in the full text:
- The hyphal wall is layered: inner alkali-insoluble fibrillar skeleton (beta-1,3-glucan with beta-1,6 branches, small chitin fraction) and an outer alkali-soluble amorphous layer carrying cell wall proteins (CWPs) covalently bound to beta-1,3-glucan via GPI remnants, or directly via an alpha-1,6 bond to N-linked galactan/mannan (section 2.1).
- Wall growth needs a balance of synthesizing enzymes, loosening (chitinase/glucanase) enzymes and rigidifying (cross-linking) enzymes (section 2.2, "integrated model").
- Beta-1,3-glucan (FKS-1 / RHO-1 / GS-1) and chitin (seven CHS) synthases are delivered by distinct vesicle populations (macrovesicles and chitosomes) to the Spitzenkoerper; Rabs, the exocyst and the cytoskeleton direct delivery (section 5).
- Function of most GH-family wall enzymes in hyphae is unproven; the authors repeatedly say "no evidence" or "inconclusive" (sections 4.1, 4.2).

## 2. Gene-set table

Source column: "32743136 s.X" = section X of Verdin et al. (read in full text). Family IDs are quoted as the paper states them. Anything I did not read is "unverified". Species column uses N. crassa locus IDs (NCU) because the paper uses them; S. cerevisiae or other species appear only where the paper makes the comparison, and species names were partly lost in extraction, so those species assignments are marked "inferred".

| Functional group | Representative genes | Species | Family / domain as stated | Role in wall remodeling | Source |
|---|---|---|---|---|---|
| beta-1,3-glucan synthase complex | FKS-1 (NCU06871, catalytic); RHO-1 (NCU01484, regulatory GTPase); GS-1 (NCU04189, putative GSC part) | N. crassa | Fks1 present across Ascomycota, Basidiomycota and several other phyla; GS-1 orthologs only in Asco/Basidiomycota; GS-1 ortholog Knr4/Smi1 (S. cerevisiae, inferred). No Pfam stated | Synthesis of the main structural polysaccharide at secretion sites | 32743136 s.3.1 |
| Rho1 control of glucan synthesis and CWI | RHO-1, RGF-1 (NCU00668, GEF), LRG-1 (GAP), PKC-1, MAK-1 MAPK | N. crassa | Rho GTPase; no family ID | RHO-1 signals to PKC-1 and MAK-1 (cell wall integrity pathway); S. cerevisiae Knr4/Smi1 is linked to PKC1-SLT2 cascade | 32743136 s.3.1 |
| Chitin synthases (CHS) | chs-1 NCU03611, chs-2 NCU05239, chs-3 NCU04251, chs-4 NCU09324, chs-5 NCU04352, chs-6 NCU05268, chs-7 NCU04350 (gene-name-to-NCU pairing is by order of listing in the paper; names were stripped, so treat the pairing as "likely") | N. crassa | Catalytic domain PF03142 (stated). Division 1 (classes III, II, I = CHS-1,-2,-3) also PF08407 and PF01644. Division 2 (classes IV, V, VII) have cytochrome b5-like domain; CHS-5 and CHS-7 have N-terminal myosin motor-like domain (MMD) and C-terminal DEK. CHS-6 class VI, signal peptide. Mucoromycota have classes I, IV, V, VII plus class VIII | Chitin synthesis at tip and septa; zymogen activation by proteolysis (proteases not identified) | 32743136 s.3.2 |
| CHS trafficking and regulation | CSE-7 (NCU05720, Chs7 ortholog, ER export of CHS-4); CBS-5 NCU07435 (Chs5 ortholog); BUD-7 NCU04511; CSA-1 NCU09322, CSA-2 NCU02592, CSR-3 NCU02351 (Chs4 orthologs); PFA-4 NCU02118, PFA-3 NCU01267; PP1 NCU00043; SEC-14 NCU02263 | N. crassa; S. cerevisiae Chs3 pathway (Chs3, Chs4, Chs5, Chs7, Pfa4, Bni4) described as the comparison | no family IDs | Getting CHS to the surface in an active state | 32743136 s.5.1 |
| Vesicle and exocytosis machinery | Rabs YPT-1, YPT-31, SEC-4; exocyst SEC-3,-5,-6,-8,-15, EXO-70, EXO-84; KIN-1; actin/microtubule motors | N. crassa | no family IDs | Delivery of synthases (not wall remodeling itself) | 32743136 s.5 |
| Chitinases | CHIT-1 (NCU02184, signal peptide plus predicted GPI); 12 putative GH-18 genes in genome | N. crassa; Cts1 and ChiA comparisons (species lost in extraction) | GH-18 (stated). Beta-N-acetylglucosaminidases GH-20 (stated) | Loosening of chitin; 10 of the 12 are non-essential; role in hyphal remodeling "inconclusive" | 32743136 s.4.1 |
| Alpha-1,3-glucanase / glucanases | NCU06010 (GH-17, mutanase, expressed during conidiation); NCU07076 (GH-81 putative beta-1,3-glucanase); NCU03914 (GH-5 beta-1,3-exoglucanase; orthologs Exg1,2,3); "at least 38" glucan-modifying proteins in the genome | N. crassa; Agn1 comparison (species lost in extraction) | GH-17, GH-81, GH-5 (stated) | Loosening; separation of cells and septal material; alpha-1,3-glucan only in spores | 32743136 s.4.1 |
| GPI-anchored beta-1,3-endoglucanases / glucanosyltransferases (BGT) | BGT-1 (NCU06381), BGT-2 (NCU09175); 64% and 47% identity to Bgt2 (S. cerevisiae, inferred) | N. crassa | GH-17 (stated) | Branching and remodeling of beta-1,3-glucan; deletions raise resistance to Calcofluor White and Congo Red; conidiation defects | 32743136 s.4.1, 4.2 |
| GEL / Gas family (beta-1,3-glucanosyltransferases) | gel-1 (NCU01162), gel-2 (NCU06781), gel-3 (NCU08909), gel-4 (NCU07253), plus NCU06850 (name stripped; 5 genes in total). Gas1, Phr family named as homologs in other fungi (exact species lost in extraction) | N. crassa | GH-72 (stated), "GEL"; predicted GPI signals | Elongation and cross-linking of beta-1,3-glucan; creates anchoring sites for mannoproteins, galactomannans, chitin; may matter for glycoprotein incorporation | 32743136 s.4.2 |
| CRH family transglycosidases | Crh1p, Crh2p, Crr1p (S. cerevisiae, inferred); N. crassa MWG-1, GH-16-7 NCU05974, CRF-1, GH-16-11 NCU09117, GH-16-14 NCU09672 (46-48% identity to Crh1); 13 GH-16 proteins in genome | S. cerevisiae (inferred), N. crassa | GH-16 (stated); GPI anchored | Transfer chitin to beta-1,6- and beta-1,3-glucan (cross-links control morphogenesis); function in N. crassa unproven | 32743136 s.4.2 |
| GPI-anchored wall proteins (CWPs) | HAM-7 NCU00881 (anastomosis receptor); ACW-2 NCU00957, ACW-3 NCU05667, ACW-7 NCU09133 (Kre9 domain, beta-glucan assembly); ACW-5 NCU07776, ACW-6 NCU03530 (CFEM domain) | N. crassa | Kre9 domain, CFEM domain (named, no Pfam ID stated) | Structural and receptor CWPs; CWP set also includes chitinases, chitosanases, glucanases, glucanosyltransferases | 32743136 s.2.1 |
| Beta-1,6-glucan synthesis genes (KRE family), Dfg5/Dcw1, mannoprotein biosynthesis | not named in the text I read, apart from "Kre9 domain" above | all | unverified | unverified | not seen in any of the three papers (Dfg5, Dcw1 do not appear in the text read) |
| Calcineurin, Slt2/Mpk1, Cek1 and HOG regulation | PKC-1, MAK-1 (CWI) named only. Calcineurin: not seen | N. crassa | unverified | MAK-1 pathway maintains integrity | 32743136 s.3.1; others not seen |
| Species other than N. crassa and S. cerevisiae (C. albicans, A. fumigatus, C. neoformans) | not seen in any text read | | unverified | | none |

Note on S. pombe: the text states that S. pombe lacks wall chitin and that Bgs1p and Bgs4p are FKS-1 orthologs at septa. The species name was stripped in extraction; I inferred it from the Bgs1/Bgs4 names. Treat as inferred.

## 3. Adhesion relevance

### 3a. Evidence stated in the papers (32743136 only)

- CWPs are covalently attached to beta-1,3-glucan through GPI remnants or to the galactan/mannan core (s.2.1). That is the attachment route for GPI-anchored wall proteins, and adhesins that use it sit in the outer amorphous layer.
- GEL/GH-72 enzymes create anchoring sites for mannoproteins and may be needed for glycoprotein incorporation (s.4.2). This is a paper statement about N. crassa GEL, with no biochemical evidence in N. crassa ("no evidence of the biochemical activity").
- alpha-1,3-glucan "agglutinates" chitin and beta-1,3-glucan in other fungi to protect their exposure (s.2.1).
- Chitinase, glucanase and Crh roles in hyphal remodeling are called inconclusive or unproven (s.4.1, 4.2).
- Calcofluor White and Congo Red resistance changes in bgt mutants show altered wall architecture (s.4.1).
- The review contains no adhesin, no Als/Epa/Flo discussion, and no beta-glucan masking or immune-exposure discussion (nothing seen; N. crassa focus).

### 3b. My inference (not stated in the papers)

Plausibly adhesion-relevant (could change what is displayed or how it is retained):
1. GPI-anchored CWP cross-linkers (GEL/GH-72, CRH/GH-16, Dfg5/Dcw1-like GH76): they set how many GPI-CWPs, which can include adhesins, are retained in the wall. Dfg5/Dcw1 are not in the text I read; the GH76 link is from my own background knowledge and is unverified.
2. Kre9-domain GPI proteins (ACW-2,-3,-7) and beta-1,6-glucan machinery: beta-1,6-glucan is a common attachment point for GPI-CWPs in yeast (background knowledge, unverified here).
3. alpha-1,3-glucan and chitin/glucan exposure (chitin synthases, chitinases, glucanases): they change surface masking and thickness. Section 5 of this report lists papers on beta-glucan masking in Candida.
4. CWI signalling (RHO-1, PKC-1, MAK-1; Cek1 in Candida per follow-up papers) regulates wall remodeling and, in Candida, unmasking.

Background housekeeping (likely to be false positives for a cell-surface glycoprotein detector, mostly because of signal peptide, GPI anchor, S/T-rich or N-glycosylated sequence rather than adhesive function):
1. GPI-anchored glucanosyltransferases and glucanases (BGT, GEL, GH-17, GH-72, GH-16). Many carry signal peptide plus GPI, which a surface-glycoprotein classifier will score high. They are enzymes with catalytic domains (stated GH families) and no adhesion claim in the paper.
2. GH-18 chitinases with signal peptide and GPI (CHIT-1).
3. CFEM-domain and Kre9-domain GPI proteins (ACW family). CFEM has a pathogenesis link per the paper, but no adhesion claim is made.
4. Intracellular machinery (FKS, CHS, Rho, Rab, exocyst, motors): multi-pass membrane or cytosolic, unlikely to score as surface glycoproteins, but CHS and FKS are polytopic membrane proteins that a hydrophobicity-based model may confuse with surface proteins.

Practical suggestions (inference): (a) treat PF03142 (CHS catalytic domain, stated in 32743136), GH-16, GH-17, GH-18, GH-72 as a "wall enzyme" exclusion or demotion set for adhesin calls, and verify each family ID against Pfam/CAZy/dbCAN before coding it; (b) keep GPI-CWP crosslinkers visible as a separate "cell wall remodeling" class so they are not called adhesins; (c) verify the claims about Dfg5/Dcw1 (GH76), KRE genes and Candida masking against the follow-up papers below before coding.

## 4. Gaps

- No gene tables from PMID 28429675 or PMID 29643171 were read.
- No wall gene sets for C. albicans, A. fumigatus or C. neoformans were obtained from the three papers.
- Dfg5, Dcw1, KRE (beta-1,6-glucan), calcineurin/Crz1, Slt2/Mpk1, Gas1 details, Pfam IDs other than PF03142, PF08407, PF01644: no support from texts read.
- S. cerevisiae and S. pombe attributions inside 32743136 are inferred because species names were lost in extraction.

## 5. Follow-up papers (each PMID resolved to the stated title in PubMed)

| PMID | Citation | Why |
|---|---|---|
| 39694372 | Wen TT et al. "Fungal beta-1,3-glucanosyltransferases: A comprehensive review on classification, catalytic mechanism and functional role." Int J Biol Macromol 2024 | GH72/GH17 transglycosylase families, domains, mechanism; fills the GEL/Gas/BGT gap |
| 36211971 | Ortiz-Ramirez JA, Cuellar-Cruz M, Lopez-Romero E. "Cell compensatory responses of fungi to damage of the cell wall induced by Calcofluor White and Congo Red ... A review." Front Cell Infect Microbiol 2022 | Wall compensation (CWI, chitin up-regulation) after wall damage |
| 42805532 | Soni A et al. "Environmental and genetic regulation of mannan architecture in Candidozyma auris (formerly Candida auris) ...: A review." Int J Biol Macromol 2026 | Mannan remodeling that masks beta-glucan; OCH1, MNT5; HOG and PKC/CWI. Candida auris, not C. albicans |
| 41891702 | King AE et al. "Loss of Fgr41 in [Candida albicans] attenuates virulence and increases proinflammatory immune responses ... dependent on beta(1,3)-glucan but not dectin-1." Infect Immun 2026 | C. albicans beta-glucan unmasking; Cek1 MAPK; a cell wall secreted protein linked to masking |
| 40531863 | Yang G et al. "Mitochondrial anchor protein Num11 is key to pathogenicity of [Candida albicans] by affecting mitochondrial function and cell wall masking." Virulence 2025 | Cdc42-Cek1 MAPK and wall remodeling with beta-glucan/chitin exposure |
| 42813837 | Kamasah JS et al. "Caspofungin-induced beta-(1,3)-glucan unmasking is associated with enhanced innate immune clearance of [Candida auris]." Microbiol Spectr 2026 | Echinocandin-driven wall remodeling and beta-glucan exposure |

Species names in some of these titles were dropped by the PubMed tool output; I filled them in brackets from the abstracts (Candida auris for 42813837 from keywords; Candida albicans for 41891702 and 40531863 from MeSH terms). I did not read these papers beyond abstracts.

Not found by my searches: classic S. cerevisiae reviews on GPI-CWP cross-linking (Dfg5/Dcw1/Crh) and beta-1,6-glucan. Recommend a targeted PubMed search by author (for example Klis, Kapteyn, Lesage, Gow) and verifying the PMIDs before use. No PMIDs are given here for them because I did not verify any.
