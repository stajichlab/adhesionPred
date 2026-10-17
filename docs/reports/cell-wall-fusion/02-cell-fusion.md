# Fungal cell fusion and its relation to adhesion and cell wall remodeling

Literature research only. No code was changed. Date: 2026-10-17.
Sources: PubMed metadata (all PMIDs below were resolved by the PubMed tool and the title matches), and PMC full text through the PubMed tool or Europe PMC XML. Where only an abstract was available this is stated.

Conventions: "Stated" means a sentence in the cited paper says it. "Inference" is mine and is labelled. "Unverified" means I did not find it stated in a paper I read. Pfam IDs are given only when they come from this repo's own `data/sorting_hat/family_table.tsv` (and are marked as such). None of the fusion papers read name a Pfam ID, so the Pfam column is "unverified (no ID stated)" unless a domain name is stated, which is then written as the domain name.

---

## 1. Scope and papers read

Scope: mating cell fusion (S. cerevisiae, S. pombe, Candida, Cryptococcus), vegetative germling/hyphal fusion (Neurospora crassa, Fusarium oxysporum), and the adhesion step that precedes fusion. Basidiomycete clamp connections are only touched (see open questions).

### The two papers the owner asked about (what they actually are)

| PMID | Title (year) | What it is | Fusion or hyphal fusion content |
|---|---|---|---|
| 31501435 | Comparative genomics reveals the origin of fungal hyphae and multicellularity (Nat Commun 2019, Kiss et al.) | Comparative genomics of 72 fungal genomes on the origin of hyphae. Not a cell fusion paper. | The full text has no occurrence of "fusion" or "anastomosis". It is useful only for adhesion: it defines 23 families of putative adhesion-related proteins and reports no expansion of adhesive, kinase or receptor families when filamentous fungi arose (stated, abstract and results). |
| 26541633 | Role of BGT-1 and BGT-2, two predicted GPI-anchored glycoside hydrolases/glycosyltransferases, in cell wall remodeling in Neurospora crassa (Fungal Genet Biol 2015, Martinez-Nunez and Riquelme) | Primary paper, GFP localisation and mutants in N. crassa. Abstract only was available to me (not in PMC). | Stated in the abstract: BGT-1-GFP and BGT-2-GFP were seen at germling and hyphal fusion sites (also septa, conidial poles, branch sites, hyphal apical dome subapex). Mutants show altered resistance to calcofluor white and Congo red and conidiation defects. No fusion defect is reported in the abstract. |

### Reviews and primary papers used

| PMID | Title (year) | Depth read | Why used |
|---|---|---|---|
| 35138122 | Chemotropism and Cell-Cell Fusion in Fungi (Microbiol Mol Biol Rev 2022, Clark-Cotton, Jacobs, Lew) | Abstract only (PMC text returned empty) | Best recent review across yeast and filamentous fungi; chemotropism, exploratory polarization, local wall degradation. |
| 35148940 | A focus on yeast mating: From pheromone signaling to cell-cell fusion (Semin Cell Dev Biol 2022, Sieber, Coronas-Serna, Martin) | Abstract only | S. cerevisiae and S. pombe pheromone signalling and fusion; secretion concentrated at contact leads to local wall digestion. |
| 31152053 | Molecular mechanisms of chemotropism and cell fusion in unicellular fungi (J Cell Sci 2019, Martin) | Abstract only | Wall expansion to wall digestion transition; membrane fusion mechanism still open. |
| 23466674 | Mate and fuse: how yeast cells do it (Open Biol 2013, Merlini, Dudin, Martin) | Full text (PMC3718343) | Gene-level source for Fus1, Fus2, Prm1, Fig1, Kex2, Dni1, Fus3, Ste proteins (Table 1 of that paper). |
| 31001214 | Communicate and Fuse: How Filamentous Fungi Establish and Maintain an Interconnected Mycelial Network (Front Microbiol 2019, Fischer and Glass) | Full text (PMC6455062, via Europe PMC XML) | Gene-level source for N. crassa MAK-1/MAK-2/SOFT/STRIPAK/NOX, membrane fusion genes, HAM genes. Over 70 genes involved. |
| 25329559 | A model for cell wall dissolution in mating yeast cells (PLoS One 2014, Huberman and Murray) | Full text (PMC4199604) | Scw4/Scw11 putative glucanases; wall dissolution is contact driven; agglutinin-mediated adhesion raises local enzyme concentration (model). |
| 30586320 | Cell fusion in yeast is negatively regulated by components of the cell wall integrity pathway (Mol Biol Cell 2018, Hall and Rose) | Full text (PMC6594448) | Mid2 and Pkc1 hold fusion back; links wall integrity sensing to Fus2/Cdc42. |
| 33397995 | Single-cell fluidic force microscopy reveals stress-dependent molecular interactions in yeast mating (Commun Biol 2021, Mathelie-Guinlet et al.) | Full text (PMC7782832) | Sexual agglutinin biology: GPI, Als-like architecture, binding force, pheromone induction. Key for adhesion overlap. |
| 24612417 | Sexual biofilm formation in Candida tropicalis opaque cells (Mol Microbiol 2014, Jones, Hirakawa, Bennett) | Abstract only | Sexual agglutinin-like FGR23, mating plus biofilm. |
| 26438189 | The Candida albicans agglutinin-like sequence family of adhesins (Future Microbiol 2015, Cota and Hoyer) | Abstract only | Als family adhesion mechanism (binds flexible C termini). |
| 29176784 | PRM1 and KAR5 function in cell-cell fusion and karyogamy ... Cryptococcus (PLoS Genet 2017, Fu and Heitman) | Abstract plus matching sentences from full text | Prm1 conservation in a basidiomycete yeast; four predicted TM domains. |
| 24514900 | Membrane organization and cell fusion during mating in fission yeast requires multipass membrane protein Prm1 (Genetics 2014, Curto et al.) | Abstract only | S. pombe Prm1 and Dni (Fig1-related claudin-like) proteins. |
| 27531696 | Suppressor of fusion, a Fusarium oxysporum homolog of Ndt80, is required for nutrient-dependent regulation of anastomosis (Fungal Genet Biol 2016, Shahi et al.) | Abstract only | Ndt80 family and anastomosis (Suf deletion gives hyper-fusion). |
| 31501435, 26541633 | see above | | |

Not found or not retrievable: no dedicated review of basidiomycete clamp-connection fusion since 2012 turned up in my PubMed queries (queries were AND-combined and returned only Schizophyllum papers that I did not read). Two items of the owner's list were not found as a dedicated paper: a primary paper on Scw/Fus-related glucanases for filamentous fungal fusion, and a Pfam-annotated list of fusion proteins.

---

## 2. Gene-set table by process step

Species abbreviations: Sc = S. cerevisiae, Sp = S. pombe, Nc = N. crassa, Ca = C. albicans, Ct = C. tropicalis, Cn/Cd = Cryptococcus neoformans/deneoformans, Fo = F. oxysporum, Pa = Podospora anserina.
"Family or domain" is only what the paper states; otherwise unverified.

### 2a. Chemotropic signalling

| Step | Genes (species) | Ortholog statement in the papers | Family or domain (stated) | Evidence | PMID |
|---|---|---|---|---|---|
| Pheromone, receptor | a-factor, alpha-factor; Ste2 (alpha-factor receptor), Ste3 (a-factor receptor) (Sc); Mam2, Map3 (Sp); Ste2 and Ste3 also required for sexual biofilm in Ct | Sc and Sp receptors are both GPCRs with 7 TM | GPCR | Stated | 23466674, 24612417 |
| Heterotrimeric G protein | Gpa1, Ste4, Ste18 (Sc); Gpa1, Gnr1 (putative) (Sp) | Sp Ggamma not identified | G protein | Stated | 23466674 |
| MAPK scaffold and cascade | Ste5 scaffold, Ste20 (PAK), Ste11/Ste7/Fus3 and Kss1 (Sc); Byr2/Byr1/Spk1 (Sp); Ste12 TF (Sc), Ste11 (Sp) | Nc MAK-2 pathway (NRC-1, MEK-2, MAK-2) is orthologous to the Sc pheromone response pathway; PP-1 is the Ste12 ortholog; HAM-5 scaffold is not an ortholog of Ste5 and is restricted to Pezizomycotina | Kinases; Ste12-like TF with homeobox-like STE domain and C2H2 Zn motifs | Stated | 23466674, 31001214 |
| Polarity | Cdc42, Cdc24 (GEF), Bem1, Far1 (Sc); Cdc42, Scd1, Scd2 (Sp); Nc CDC-42, RAC-1, CDC-24 | Orthologs named in Nc review | GTPase, GEF | Stated | 23466674, 31001214 |
| Nc somatic (non-pheromone) signalling | MAK-2 complex: HAM-5, MAK-2, MEK-2, NRC-1, STE-50, STE-20, RAS-2; HYM-1 (NDR kinase scaffold); oscillates with ~4 min period at germ tube tips | Pheromone GPCRs and G proteins are dispensable for somatic germling and hyphal fusion, but needed for mating fusion | Kinases, scaffolds | Stated (cell biology, mutants, AP-MS) | 31001214 |
| Nc CWI branch | MIK-1, MEK-1, MAK-1; scaffold SOFT (also named HAM-1); oscillates out of phase with MAK-2 ("ping-pong" model, SOFT proposed for signal secretion, MAK-2 for reception - a model) | SOFT is Pro40 in Sordaria; Sc has orthologs of mak-1 pathway kinases | Kinases, scaffold | Stated | 31001214 |
| STRIPAK | PPG-1, PP2A-A, HAM-3 (striatin), HAM-2, HAM-4 (membrane anchors), MOB-3; SCI1 | Mouse striatin complements Nc ham-3 mutant | Phosphatase complex; striatin has calmodulin binding motif | Stated | 31001214 |
| Transcription factors | PP-1 (Ste12-like), ADV-1 (Zn(II)2Cys6) | ADV-1 clear homologs restricted to Pezizomycotina | Zn cluster, C2H2 | Stated | 31001214 |
| Non-self recognition before and after fusion | doc loci (greenbeard genes, 5 communication groups); het loci sec-9, plp-1 mediate post-fusion death | | Unverified | Stated | 31001214 |
| Fo nutrient regulation | Suf (Ndt80-family homolog; hyper-fusion when deleted); Nc Vib-1 is the homologous regulator of heterokaryon incompatibility | Suf is the Fo homolog of Nc Vib-1 | Ndt80 (p53-like) family TF | Stated (abstract) | 27531696 |

### 2b. Adhesion and agglutination (before fusion)

| Genes (species) | Ortholog statement | Family or domain (stated) | Evidence | PMID |
|---|---|---|---|---|
| alpha-agglutinin (MATalpha, Sc); gene symbol not stated in that paper, so "SAG1" is unverified here | Structure "highly reminiscent of the Als adhesin family of C. albicans" | 650 aa; signal peptide (19 aa); GPI anchor addition signal; glycosylated C-terminal stalk (~300 aa); N-terminal beta-sheet-rich Ig-like region homologous to Als N-terminal domains; contains the a-agglutinin binding site | Stated | 33397995 |
| a-agglutinin: Aga1 (725 aa, wall anchor subunit) plus Aga2 (69 aa, binding subunit) (MATa, Sc) | | Aga2 joined to Aga1 by two disulfide bonds; Aga2 C-terminal peptide inserts in the alpha-agglutinin binding cleft; both agglutinins are covalently linked to wall polysaccharide through the GPI remnant | Stated | 33397995 |
| Agglutinin regulation | Pheromone raises agglutinin expression; agglutinins constitutively low; single a-alpha bond ~100 pN, requires pheromone induction; tension strengthens binding | | | Stated | 33397995 |
| FGR23 (Ct) | "agglutinin-like protein"; enhances mating and sexual biofilm | Suggests conserved role of sexual agglutinins in mating, cell cohesion and biofilm (authors' wording) | Agglutinin-like (Pfam not stated) | Stated (abstract) | 24612417 |
| Als family (Ca) | Ca Als binds flexible C termini of many ligands | Als binding domain structure solved | Als domain | Stated (abstract) | 26438189 |
| Nc | No adhesin, agglutinin or flocculin is named in the Nc communication review (the text has no occurrence of "adhes" or "agglutin") | | | Stated by absence in the text I searched | 31001214 |

### 2c. Cell wall remodeling at contact

| Genes (species) | Role stated | Family or domain (stated) | Evidence | PMID |
|---|---|---|---|---|
| Fus1 (Sc) | O-glycosylated 1-pass plasma membrane protein; localised by Cdc42/Cdc24/Chs5; needed for vesicle clustering and wall removal; double fus1 fus2 mutants arrest with wall between cells | 1-pass TM, O-glycosylated | Stated | 23466674 |
| Fus2 (Sc) | Binds GTP-Cdc42; with Rvs161 (BAR domain amphiphysin); expanding ring proposed to remove wall remnants | Putative Rho-GEF domain; BAR (Rvs161) | Stated | 23466674, 30586320 |
| Kel1, Rvs161, Spa2, Bni1, Tpm1, Myo2 (Sc) | Actin cables and vesicle delivery; Myo2 moves chitin synthase 3 and Ste5 to the shmoo tip | | Stated | 23466674, 30586320 |
| Scw4, Scw11 (Sc) | Putative cell wall glucanases implicated in mating; Scw11 promoter has a Ste12 site; deleting SCW4 and SCW11 reduces pheromone-induced death by 40% in the confinement assay. Enzymes are described as released at the shmoo tip (hypothesis); "putative" in the paper, no activity shown there | Glucanase homology (stated as "homology to known cell wall glucanases") | Stated, partly hypothesis | 25329559 |
| Mid2, Wsc1 (Sc CWI sensors), Pkc1, Rho1 | Mid2 and Pkc1 negatively regulate fusion via Cdc42 and Fus2; mid2 deletion causes "mating-induced death" by premature fusion. Sensors have highly O-mannosylated extracellular domains | TM CWI sensors | Stated | 30586320 |
| Fks1 (Sc beta-1,3-glucan synthase, Rho1 target) | Wall deposition | | Stated | 30586320 |
| HAM-7, WSC-1 (Nc) | Both encode GPI-anchored wall proteins required for wall maintenance and MAK-1 phosphorylation; only Delta ham-7 is defective in germling and hyphal fusion | GPI-anchored (stated) | Stated | 31001214 |
| BGT-1 (NCU06381), BGT-2 (NCU09175) (Nc) | Predicted GPI-anchored glycoside hydrolase/glycosyltransferase; localised to germling and hyphal fusion sites among others; single and double mutants alter stress sensitivity; no fusion defect reported in abstract | GH, glycosyltransferase, GPI anchor site (predicted) | Stated (abstract) | 26541633 |
| Chitin synthase 3 (Sc) | Delivered by Myo2 to shmoo tip | | Stated | 23466674 |
| NCU05836 (Mnl2p-like ER mannosidase), NCU05916 (putative alpha-1,3-mannosyltransferase) (Nc) | ADV-1 targets; mutants reduced in communication | Glycan-processing enzymes | Stated | 31001214 |
| LAO-1 (Nc) | L-ascorbate-oxidase-like laccase; mutant completely defective in germling fusion | Multiple copper-oxidase domains | Stated | 31001214 |
| Sp Dni1 | Claudin-like Fig1-related; roles in plasma membrane and wall remodeling at fusion | Fig1-related | Stated | 23466674, 24514900 |

### 2d. Membrane fusion

| Genes (species) | Role stated | Family or domain (stated) | Evidence | PMID |
|---|---|---|---|---|
| Prm1 (Sc, Sp, Nc PRM-1, Cn/Cd PRM1) | Plasma membrane fusion after wall breakdown; 4-pass TM; two extracellular loops; not the fusogen (Merlini reading). Sp prm1 deletion nearly blocks fusion. | 4-pass TM (Phobius for Cryptococcus) | Stated | 23466674, 24514900, 29176784, 31001214 |
| Fig1 (Sc), Dni1 (Sp), FIG-1 (Nc) | 4-pass TM Ca2+ influx regulator; Nc FIG-1 is a low-affinity Ca2+ channel | Claudin-like (Sp/Sc) | Stated | 23466674, 31001214 |
| Kex2 (Sc) | Golgi protease; processes alpha-factor; acts additively with Prm1 | Protease | Stated | 23466674 |
| LFD-1, LFD-2, NCU01697 (Nc) | Mutants chemotrop normally but fail fusion with lysis; LFD-1 is Pezizomycotina-specific 1-TM; LFD-2 is Golgi/ER; NCU01697 has predicted acyltransferase domain | 1-TM (LFD-1); acyltransferase (NCU01697) | Stated | 31001214 |
| Ergosterol: erg-2 (Nc) and ergosterol/sphingolipid/PIP2 enrichment at fusion site (Sc) | Sterol composition affects the switch from communication to fusion | Lipid biosynthesis | Stated | 31001214, 23466674 |
| Sp Fus1 formin (unrelated to Sc Fus1), Cdc8, Myo51 | Formin essential for Sp fusion | Formin | Stated | 23466674 |

### 2e. Other fusion-required Nc genes (membrane or secretory)

| Genes (Nc) | Stated feature | PMID |
|---|---|---|
| NOX-1 complex: nox-1, nor-1 (NoxR), HAM-6 (TM protein, Pezizomycotina-specific), BemA, Cdc24, Cdc42, RacA; IQG-1 scaffold links Nox1 and Mak1 in a relative | NADPH oxidase; superoxide cannot cross membranes so direct ROS signalling is "highly unlikely"; plasma membrane localisation of NOX-1 needed for germling fusion | 31001214 |
| ham-8 (4-TM MARVEL-like), ham-9 (pleckstrin domain), ham-10 (putative C2 domain), ham-11 (Sordariaceae-only PM protein), ham-13, ham-14, cch-1 (Ca2+ channel), cse-1/pik-1/nfh-1 (Ca-dependent exocytosis), Arp2/3 subunits, BAR-protein amphiphysin, exocyst Sec15 | Domains named as stated in the review; Pfam IDs not given | 31001214 |

### 2f. Post-fusion

| Item | Statement | PMID |
|---|---|---|
| Karyogamy: Kar5 (Cn/Cd nuclear membrane fusion), Kar3 and microtubules (Sc) | Nuclear congression and fusion after cell fusion | 29176784, 23715006 (title: When yeast cells meet, karyogamy!; not read beyond abstract) |
| Heterokaryon incompatibility (het) cell death | Post-fusion non-self checkpoint in Nc | 31001214 |
| Cryptococcus unisexual cycle | Cell fusion is dispensable for solo unisexual reproduction in C. deneoformans | 29176784 |

---

## 3. Which proteins a cell-surface classifier or Pfam call could pick up

This section is mostly inference built on the stated features above.

### Likely visible (signal peptide, GPI anchor and/or heavily glycosylated secreted or wall proteins)

| Protein | Stated feature supporting visibility | Comment (inference) |
|---|---|---|
| Sc alpha-agglutinin | Signal peptide, GPI signal, glycosylated stalk, Als-like Ig N-terminal domain (33397995) | The repo already lists Pfam Candida_ALS_N (PF11766) and Candida_ALS (PF05792) as inactive in `family_table.tsv`. Whether Sc alpha-agglutinin hits these is untested by me; the Als-like description only supports it being a candidate. |
| Sc Aga1 (GPI wall anchor) | Wall-anchored through GPI remnant (33397995) | Aga2 is a 69 aa peptide; very short proteins may fall below a length filter (inference). |
| Ct FGR23 and Candida agglutinin-like proteins | Agglutinin-like (24612417) | Same Als-type family logic. |
| Nc HAM-7 and WSC-1 | GPI-anchored wall proteins (31001214) | Fusion relevance is for HAM-7 only. |
| Nc BGT-1 and BGT-2 | Predicted GPI-anchored GH (26541633) | Enzymes, not adhesins. A wall-protein class call (GPI, GH) fits; an adhesin call would be wrong. |
| Sc Scw4 and Scw11 putative glucanases | Putative wall glucanases (25329559) | Secreted-wall-enzyme class. Wall-degrading enzymes being free in the wall is stated generally for yeast; I did not find a statement that these two are GPI anchored. |
| Sc Fus1 (O-glycosylated 1-pass TM) | O-glycosylated plasma membrane protein (23466674) | TM filter would remove it in a no-TM rule (inference). |

### Probably not visible to a surface classifier or Pfam adhesion-family call

- Signalling kinases and scaffolds (Ste5, Ste20, Ste11, Ste7, Fus3, Kss1, MAK-1, MAK-2, NRC-1, MEK-1/2, SOFT, HAM-5), STRIPAK, GTPases (Cdc42, Rho1), PP2A and transcription factors (Ste12/PP-1, ADV-1, Suf). Intracellular; no signal peptide expected.
- Multi-pass or TM plasma membrane proteins (Prm1, Fig1, Dni1, FIG-1, LFD-1, HAM-6, HAM-8, HAM-11, Ste2, Ste3, Mid2 and Wsc1 sensors). They sit in the membrane. They would not be wall proteins; with a no-TM rule they would be excluded, and with a signal-peptide-only rule they might be misjudged either way (inference).
- Golgi or ER proteins (Kex2, LFD-2, mannosidase, mannosyltransferase).
- Cytoskeleton and secretory proteins (formins, myosins, exocyst, BAR proteins).
- NOX complex proteins (membrane enzymes; stated to act from the plasma membrane).

The practical consequence (inference): a Pfam or surface classifier would catch the adhesion step (agglutinins) and some wall enzymes, and would essentially miss the signalling and membrane fusion machinery. The fusion core is mostly signalling plus membrane proteins, and cell wall remodeling at the contact site is done by few, mostly "putative" enzymes in the papers read.

---

## 4. Overlap of fusion with adhesion biology: evidence versus inference

### Evidence (stated in papers read)

1. Sexual agglutinins are cell-surface adhesins that mediate cell-cell contact as "the first initial step in the cell mating process that leads to fusion" in S. cerevisiae (33397995). They are GPI-anchored glycoproteins, and alpha-agglutinin shares architecture with Als adhesins (33397995).
2. Pheromone induces agglutinin expression; agglutinins are constitutively low (33397995). Binding strengthens with prolonged contact and mechanical tension (33397995).
3. Mating cells "initially bind to each other at their shmoo tips using mating agglutinins" before wall dissolution (25329559). Contact via agglutinins is proposed to raise the local concentration of wall remodeling enzymes, restricting wall dissolution to touching surfaces (25329559, a model supported by confinement experiments).
4. In C. tropicalis, an agglutinin-like gene (FGR23) enhances both mating and sexual biofilm formation, and the biofilm depends on pheromone receptors Ste2 and Ste3 (24612417). The authors suggest a conserved role for sexual agglutinins in mating, cell cohesion and biofilm formation.
5. In the filamentous Ascomycete Nc, the reviews I read describe somatic fusion as chemotropic signalling (MAK-2, SOFT, STRIPAK, NOX) and membrane fusion genes, and do not mention adhesins, agglutinins or flocculins at all (31001214, text search). The pheromone GPCRs are dispensable for somatic fusion (31001214).
6. Across 72 fungal genomes, adhesion-related families (flocculins, lectins, hydrophobins, GPI wall proteins, 23 families) did not expand when filamentous fungi arose, and adhesion "becomes key only in fruiting bodies or attachment to host surfaces" (31501435).
7. Wall integrity machinery is a brake on fusion: Mid2 and Pkc1 negatively regulate Sc fusion (30586320). In Nc, HAM-7 (GPI wall protein) feeds MAK-1 and is needed for fusion (31001214). Both link wall-resident or wall-sensing proteins to fusion control.
8. Wall remodeling enzymes appear at fusion sites in Nc (BGT-1, BGT-2 localisation, 26541633), and putative glucanases Scw4 and Scw11 contribute in Sc (25329559).

### Inference (mine, not stated)

- Sexual agglutinin biology is the cleanest overlap: the same GPI-anchored, Als-like surface protein class that adhesion callers target is the one that precedes fusion. A call of "cell wall adhesion candidate" on an agglutinin is therefore a call of both adhesion and mating-contact function. Which of the two is true for an unknown Als-like or flocculin-like protein cannot be separated by sequence alone.
- Flo-type flocculins (FLO1, FLO5, FLO9, FLO10, FLO11 in this repo's table) were not discussed in the fusion papers I read. I found no statement linking them to fusion. Flocculation as lectin-mediated cell-cell adhesion is a different phenotype from sexual agglutination. Any overlap is by shared surface architecture (GPI, repeats), not shown function.
- Cryptic or non-agglutinin adhesion in hyphal fusion (Nc, Fusarium) is not supported by the papers read; fusion tubes may rely on chemotropic contact rather than a dedicated adhesin (hypothesis, evidence is absence of mention).
- The biofilm link is supported only for Candida sexual biofilms (item 4). I do not infer that fusion machinery is part of ordinary biofilm.
- Hydrophobins are in the 23 adhesion-related families of 31501435 but are not tied to fusion in any paper read.

---

## 5. Open questions and a short follow-up reading list

### Open questions

1. Do Sc alpha-agglutinin, Aga1 and Candida sexual agglutinin-like proteins hit the Als Pfam models already in `family_table.tsv` (PF11766, PF05792, both inactive)? This is testable on proteomes in the repo; I did not run it.
2. Are Scw4 and Scw11 GPI anchored or free secreted proteins, and what are their Pfam domains? Not stated in the papers read; check UniProt/Pfam directly (unverified).
3. Which Nc or Fusarium wall enzymes (glucanases, chitinases) are induced at fusion? Only BGT-1/BGT-2 localisation was found; no knockout fusion phenotype in the abstract. Search with Neurospora gene names (e.g. NCU numbers) should be done.
4. Basidiomycete mating and clamp fusion genes: no review found in this search; Fischer and Glass mention clamp connections only briefly (31001214). Candidates retrieved but not read: PMID 31416583 (Schizophyllum commune septa), 29129696 (mushroom mitosis and nuclear migration), 30062773 (Thn1 and pheromone response in S. commune).
5. Candida albicans mating fusion genes, Fus1/Fus2 equivalents in filamentous fungi: not covered; only C. tropicalis sexual biofilm and Cryptococcus PRM1/KAR5 were read.
6. Sc Fus1 versus Sp Fus1 are unrelated proteins with the same name (stated in 23466674). Any ortholog-based mapping across species must use the mapping in that review, not names.
7. Whether PRM-1, LFD-1 and HAM-8-like proteins have recognisable Pfam domains is not stated in these papers (MARVEL, pleckstrin and C2 are named as domain types only).
8. Full text of 35138122, 35148940 and 31152053 was not read (only abstracts; PMC returned empty text for the first). Their gene-level details should be checked before this report is used for a gene list.

### Follow-up reading list (all PMIDs verified through PubMed metadata; not all read in full)

| PMID | Title (year) | Why |
|---|---|---|
| 35138122 | Chemotropism and Cell-Cell Fusion in Fungi (2022) | Read full text for cell wall degradation in filamentous fungi. |
| 31152053 | Molecular mechanisms of chemotropism and cell fusion in unicellular fungi (2019) | Wall digestion and membrane fusion model. |
| 35148940 | A focus on yeast mating: From pheromone signaling to cell-cell fusion (2022) | Updated Fus1/Fus2/Prm1 summary. |
| 34382996 | Ultrastructural plasma membrane asymmetries in tension and curvature promote yeast cell fusion (J Cell Biol 2021) | Membrane fusion mechanism; I only grepped it. |
| 26438189 | The Candida albicans agglutinin-like sequence family of adhesins (2015) | Als structure, binding peptides. |
| 24612417 | Sexual biofilm formation in Candida tropicalis opaque cells (2014) | FGR23 and agglutinin-like proteins. |
| 27531696 | Suppressor of fusion, a Fusarium oxysporum homolog of Ndt80 (2016) | Ndt80-family anastomosis regulator. |
| 33397995 | Single-cell fluidic force microscopy ... yeast mating (2021) | Re-read for the agglutinin binding mechanism; sequence details (H273, disulfides). |

---

## Method notes

- PubMed searches treated all words as AND, so long queries returned nothing; short queries were used.
- Gene names lost in plain-text extraction of PMC6455062 were recovered from Europe PMC full-text XML.
- No claim in this report is taken from my own memory except where marked "unverified"; e.g. the S. cerevisiae SAG1 symbol and any Pfam IDs for fusion proteins are not stated in the papers and are not asserted.
