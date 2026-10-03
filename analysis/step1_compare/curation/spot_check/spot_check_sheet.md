# Spot-check sheet (15 of 59 rows, seed 20261003)

For each row decide: **ok**, **wrong** (say why) or **unsure**. Reviewer verdicts are hidden. Pool: all draft rows from H99 pass 2 and the U. maydis pilot, minus 3 rows the reviewer rejected (SOD1, CIG1 plasma membrane, VCX1).

## H99:14 LAC1 (Cryptococcus neoformans var. grubii H99) -> P-ext

- Accession J9VY90; GO GO:0005618; code IDA; selected_by_predictor no; PMID 11500433
- Evidence: A monoclonal antibody to the C. neoformans laccase was generated and used to show localization in the cell walls of representative serotype A (H99) and serotype D (B-3501) strains by immunoelectron microscopy. | strain/experiment: strain: H99 stated in the quote (plus B-3501); text names the construct CNLAC1-GFP; the antibody sentence says 'the C. neoformans laccase', mapped to LAC1 by that CNLAC1 mention, not by an ID in the sentence; text checked: PMC98673 PMC HTML page; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## H99:15 APH1 (Cryptococcus neoformans var. grubii H99) -> P-ext

- Accession J9VHR6; GO GO:0005576; code HDA; selected_by_predictor yes; PMID 25227465
- Evidence: We performed a proteomic analysis on secretions obtained from the more-virulent encapsulated serotype A strain of C. neoformans, H99, and identified CNAG_02944 (designated Aph1) as the acid phosphatase responsible for the extracellular acid phosphatase activity previously observed by others | strain/experiment: strain: H99, from the quote; H99 culture secretions by mass spectrometry; selected_by_predictor yes: paper says only Aph1 is predicted secreted among candidate acid phosphatases; text checked: PMC4172073 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## H99:17 CNAG_05312 (Cryptococcus neoformans var. grubii H99) -> P-ext

- Accession J9VPN6; GO GO:0005576; code HDA; selected_by_predictor no; PMID 26453029
- Evidence: This approach led to the unique identification of the novel secreted protein (CNAG_05312) that was specifically associated with modulation of Pka1 activity and not found in other proteomic studies. | strain/experiment: strain: H99, Methods: 'The C. neoformans var. grubii wild-type strain H99 (WT) and the PGAL7::PKA1 strain with galactose-inducible/glucose repressible expression of PKA1 were used for this study'; culture supernatant LC-MS/MS; text checked: PMC4600298 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## H99:19 CRZ1 (Cryptococcus neoformans var. grubii H99) -> N-int

- Accession J9VE33; GO GO:0005634; code IDA; selected_by_predictor no; PMID 23251520
- Evidence: Exposure to elevated temperature (30–37°C vs 25°C) and extracellular calcium caused calcineurin-dependent nuclear accumulation of Crz1-GFP. | strain/experiment: strain: H99, Methods: 'The strains Crz1-GFP, Δcrz1, Crz1-GFP:Pab-dsRed were created from H99'; location is condition dependent (cytosol and nuclei at 25 C; puncta under salt and heat shock); text checked: PMC3520850 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## H99:29 UGG1 (Cryptococcus neoformans var. grubii H99) -> N-sec

- Accession J9VG93; GO GO:0005783; code IDA; selected_by_predictor no; PMID 40424315
- Evidence: Furthermore, green fluorescence protein (GFP)-tagged Ugg1, Mns1, and Mns101 proteins colocalized with the ER marker, supporting their role as confirming them to be functional ERQC components based on their subcellular localization in the ER (Figure 3—figure supplement 2). | strain/experiment: strain: H99; Methods of the same paper build deletions in 'the C. neoformans serotype A strain H99'; the GFP strains were made in 'the WT strain', not named H99 in that sentence; GFP fusion colocalised with ER-Tracker; ID: CNAG_03648 (Methods); text checked: PMC12113280 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## H99:30 MNS1 (Cryptococcus neoformans var. grubii H99) -> N-sec

- Accession J9VN97; GO GO:0005783; code IDA; selected_by_predictor no; PMID 40424315
- Evidence: Furthermore, green fluorescence protein (GFP)-tagged Ugg1, Mns1, and Mns101 proteins colocalized with the ER marker, supporting their role as confirming them to be functional ERQC components based on their subcellular localization in the ER (Figure 3—figure supplement 2). | strain/experiment: strain: H99; Methods of the same paper build deletions in 'the C. neoformans serotype A strain H99'; the GFP strains were made in 'the WT strain', not named H99 in that sentence; GFP fusion colocalised with ER-Tracker; ID: CNAG_02081 (Methods); text checked: PMC12113280 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## H99:6 CDA1 (Cryptococcus neoformans var. grubii H99) -> P-ext

- Accession J9VPD7; GO GO:0016020; code IDA; selected_by_predictor no; PMID 40424315
- Evidence: Cda1 was primarily detected in the insoluble cellular protein fraction, encompassing the cell wall and membrane compartments (Figure 7G, left), consistent with its cell surface localization via a GPI anchor. | strain/experiment: strain: paper builds mutants in H99 (Methods: 'introduced into the C. neoformans serotype A strain H99'); the WT strain used for this blot is not named as H99 in the sentence; anti-Cda1 western of total, soluble, insoluble fractions; wall and membrane not separated; secreted fraction also blotted (Fig 7G right) but WT secreted level not stated in text; text checked: PMC12113280 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## UM:13 Pit2 (Ustilago maydis 521) -> P-ext

- Accession A0A0D1EAR7; GO GO:0005576; code IDA; selected_by_predictor no; PMID 36224193
- Evidence: The SG200 strain expressing Pit2-mCherry, an effector that was previously shown to localize in the maize apoplast20, was used as a positive control for secretion. ... Confocal microscopy confirmed that both Erc1-mCherry and Pit2-mCherry showed fluorescent signals on the outer surface of fungal hyphal tips, while the Int.mCherry showed fluorescent signals only inside the fungal cell (Fig. 3a). | strain/experiment: strain: SG200 named; Pit2-mCherry used as the positive control for secretion in this paper; cytosolic mCherry as negative control; accession: mapped by exact gene name PIT2 (paper text does not give the UMAG ID); proteome header GN=PIT2 UMAG_01375; text checked: PMC9556619 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## UM:19 Llp1 (Ustilago maydis 521) -> P-ext

- Accession A0A0D1E6R1; GO GO:0005618; code IDA; selected_by_predictor unknown; PMID 39997458
- Evidence: To assess whether Llp1, expressed under control of its native promoter, attaches to fungal hyphae during colonization, leaf samples infected with SG200 Llp1-HA were harvested 4 dpi and were partially macerated and then subjected to anti-HA immunostaining as described before [28,47]. Tight association of Llp1 to the fungal cell wall was confirmed in this experiment (Figure 3b). | strain/experiment: strain: SG200 Llp1-HA named; HA tag at the native locus (C-terminal); GPI status not stated in the text read; partially macerated infected leaves, anti-HA immunostain; accession: accession by UMAG_00027; text checked: PMC11857070 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## UM:2 Rsp3 (Ustilago maydis 521) -> P-ext

- Accession A0A0D1DYI3; GO GO:0005618; code IDA; selected_by_predictor unknown; PMID 29703884
- Evidence: To investigate whether Rsp3 when expressed from its native promoter attaches to fungal hyphae during colonization, leaf samples infected by SG200Δrsp3-rsp3-HA or SG200 Pcmu1-mCherry-AvitagHA30 expressing cytosolic mCherry-AvitagHA from the cmu1 promoter were subjected to anti-HA immunostaining after partially macerating the infected plant tissue. Rsp3-HA was detected around the outside of fungal hyphae (Fig. 4c, upper panel) while non-secreted mCherry-AvitagHA could not be detected (Fig. 4c, lower panel). ... Secreted Rsp3-HA binds to the U. maydis cell wall. | strain/experiment: strain: SG200 derivatives named in the quote (SG200Δrsp3-rsp3-HA, native promoter, non-permeabilised anti-HA immunostain of infected leaves); wall term rests on the paper's figure title; C-terminal HA on a protein the paper says has no GPI anchor; accession: accession by UMAG_03274 stated in paper; checked in UniProt REST and proteome FASTA; text checked: PMC5923269 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## UM:21 Xyn11A (Ustilago maydis 521) -> P-ext

- Accession Q4P0L3; GO GO:0005576; code IDA; selected_by_predictor no; PMID 34947062
- Evidence: As expected, given the presence of signal peptides in their N-terminal regions, Xyn2 and Xyn11A, but not Xyn3, were found to be secreted (Figure 1B,C). ... (C) Secretion of Xyn2, Xyn11A, and Xyn3 (tagged with GFP) was assayed in a colony secretion assay. Wild-type SG200 strain was used as control for proper colony washing, and SG200 cells expressing cytoplasmic GFP under control of constitutive otef promoter served as cell lysis control (SG200-GFP). | strain/experiment: strain: SG200 (as for Xyn2); GFP-tagged colony secretion assay; accession: accession by UMAG_06350; text checked: PMC8706147 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## UM:24 Sir2 (Ustilago maydis 521) -> N-int

- Accession A0A0D1CD62; GO GO:0005634; code IDA; selected_by_predictor no; PMID 37113216
- Evidence: Consistent with their localization motifs, we observed that Sir2 and Hst4 displayed nuclear localization (Figure 1B). ... (B) Subcellular localization of the indicated U. maydis sirtuins tagged with eGFP in its endogenous loci. | strain/experiment: strain: NOT NAMED for the eGFP-tagged strains in the text read; the deletions were made in SG200 and the controls are called wild-type, so SG200 is probable but not stated (strain_unstated); eGFP at the endogenous locus; the sentence names two sirtuins and is shared with the other protein of the pair; accession: accession by UMAG_00963; text checked: PMC10126416 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## UM:29 Tay1 (Ustilago maydis 521) -> P-ext

- Accession A0A0D1C393; GO GO:0005576; code IDA; selected_by_predictor yes; PMID 34166468
- Evidence: c. Secretion of Tay1 during maize infection. Left: Plants infected with the U. maydis strain SG200ΔC carrying the construct Ptay1:SPTay1-mCherry-Tay1 in the ip locus show accumulation of the mCherry signal at the hyphal edges, tips and cell-to-cell crosses. Right: plants infected with the same strain and plasmolyzed with 1 M mannitol show accumulation of the mCherry signal in the hyphae as well as the apoplastic space. | strain/experiment: strain: SG200ΔC (SG200 derivative) named; N-terminal mCherry behind the Tay1 signal peptide from the tay1 promoter, plasmolysis; selected as SignalP-predicted secreted Pleiades; accession: accession by UMAG_03752; text checked: PMC8224859 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## UM:7 Stp4 (Ustilago maydis 521) -> P-ext

- Accession A0A0D1C5E3; GO GO:0005576; code IDA; selected_by_predictor yes; PMID 33941900
- Evidence: To test whether Stp complex members are secreted, we expressed all HA-tagged Stp proteins from a constitutive promoter. In the supernatants of their respective cultures, Stp1, Stp2, Stp3 and Stp4 proteins could be detected by western blot ... a-d, HA-specific western blot to reveal secretion of Stp proteins after constitutive expression from the otef promoter in AB33 (a), SG200Δkex2 (b) or SG200 (c, d). | strain/experiment: strain: AB33, SG200Δkex2 or SG200 depending on panel (panel for this protein not checked; Stp1 needed the kex2 deletion strain); HA-tagged, otef promoter (overexpressed); selected as predicted effector genes; supernatant western with tubulin lysis control; accession: accession by UMAG_12197 from the paper's data availability statement; text checked: PMC8159752 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____

## UM:9 Stp3 (Ustilago maydis 521) -> P-ext

- Accession A0A0D1EAN2; GO GO:0005576; code IDA; selected_by_predictor yes; PMID 33941900
- Evidence: To localize Stp1 and Stp3 at higher resolution, immunoelectron microscopy was performed on leaf sections infected with strains expressing HA fusion proteins from their native promoters. Antibody labelling was carried out prior to fixation and thin sectioning of maize leaves, and the only fungal hyphae accessible to the label were those inside plant cells that had been cut open so that plant cytoplasm leaked out. Here, the immunolabel was found outside the hyphal cell wall and in contact with protrusions extending into the plant cell for both Stp1-HA and Stp3-HA. ... In an untagged SG200 strain, protrusions were seen, but no specific labelling was detected (Fig. 2d). | strain/experiment: strain: SG200 derivatives (HA fusions from native promoters, untagged SG200 control named); immuno-EM, label outside the hyphal wall; the paper also describes the Stp complex as membrane anchored; accession: accession by UMAG_00715; text checked: PMC8159752 EuropePMC XML; retrieved 2026-10-02
- Your verdict: ____  Comment: ____
