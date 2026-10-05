# Fungal cell wall and surface gene classes against the tools in this repository

**Working note, 2026-10-04.** Input to `docs/superpowers/specs/2026-10-04-orchestrator-design.md`.

> **What this is.** A map from the gene classes that three fungal cell wall publications discuss to
> the tools this repository has today. **What this is not.** It is not a test. No class below was
> run through any tool for this note. "Detectable" means the tool design fits the class. It does
> not mean the tool has been measured on it.

**Scope decision, 2026-10-04.** The owner set the synthase, glucanase and degradation-enzyme classes (sections 2 and 3) as later work. The tool's main goal is to find cell surface proteins for adhesion, antigen, allergen and surface glycoprotein questions. See the orchestrator spec, section 3.7.

## 1. Sources and how I read them

| Source | What I read | Limit |
|---|---|---|
| Gow, Latgé, Munro 2017, *The Fungal Cell Wall: Structure, Biosynthesis, and Function* (PMID 28513415, [DOI](https://doi.org/10.1128/microbiolspec.FUNK-0035-2016)) | Abstract (PubMed) and a machine-made list of section headings and named molecules from the PMC page | The class list below is from that heading summary. I did not read the full text. Check names before citing. |
| Riquelme, Munro, Gow 2020, editorial for the special issue *The fungal cell wall: biology, biosynthesis and biotechnology* (PMID 32743149, [DOI](https://doi.org/10.1016/j.tcsw.2020.100037)) | Full text (PubMed Central) | It is an editorial. It lists nine papers and gives no gene list of its own. |
| Gow, Casadevall, Fang 2023, *Top five unanswered questions in fungal cell surface research*, The Cell Surface (PMC10654581) | Machine-made summary of the PMC page | Same limit as the 2017 row. |

Pfam accession numbers below were looked up in the local `Pfam-A.hmm`
(`/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm`, 30,134 models). I confirmed that each
accession exists and has the name shown. I did not check that the family covers every gene named in
the papers.

## 2. The main finding for design

The papers describe the cell wall by **function and process**. Most of their gene classes are not
the population that step 1 targets. Step 1 asks: is this protein secreted or GPI-anchored? Three
kinds of class fall outside that question:

1. **Membrane-embedded synthases.** Chitin synthases, β-1,3-glucan synthase, α-1,3-glucan synthases.
   They have transmembrane segments and no signal peptide in the step 1 sense.
2. **Cytosolic and signaling proteins.** Pkc1, Bck1, Rho1, MAPKs, calcineurin, Hog1.
3. **Non-protein components.** β-glucans, chitin, mannan, GXM, galactosaminoglycan, melanin.
   The tools work on proteins. They can flag the enzymes that make these.

Two kinds of class sit inside the step 1 population but are **not adhesins**. The architecture
documents treat them as hard negatives for the adhesin question. The papers treat them as central
cell wall classes: GPI-anchored remodeling enzymes (GEL, PHR, GAS; yapsins) and hydrolases
(chitinases, glucanases). A sorter that has "adhesin" as its only cell wall label would call these
"other". The orchestrator needs a separate label for them.

## 3. Class by class

Column "Today" says what the repository has. Column "Needs" says what is missing.

| Class (examples named in the sources) | Where it sits | Family signal (Pfam, local release) | Today | Needs |
|---|---|---|---|---|
| **GPI-anchored and secreted cell wall proteins** (Als, Iff, Epa, Hwp1, Eap1, Sun41; GPI protein counts: 66 predicted in *S. cerevisiae*, over 100 in *Candida*) | Step 1 population | none (generic) | Step 1 measured (rule R2 recall 0.42, FPR 0.006 on S1:all). Not decided. No model ships. | Decision on rule, ML or hybrid. GPI truth rows (none exist). |
| **Als adhesin family** | Step 2a/2b | Candida_ALS_N PF11766, Candida_ALS PF05792 | Repeat detectors, validated in Saccharomycotina only. No Als-specific scan. | HMM wrapper with a specificity test. Clade scope is *Candida*. |
| **CFEM proteins** (Rbt5, Pga10/Rbt51, Csa1/Wap1; also Ag2/PRA) | Step 2b-i | CFEM PF05730 | Fold work done (`docs/reports/2026-09-29-class2b-structure.md`). No scan wrapper. | Wrapper and specificity test. |
| **Hydrophobins** (RodA, RodB) | Step 2c | Hydrophobin PF01185, Hydrophobin_2 PF06766; also Eas PF22354, DewD PF28987 | HMMs exist. No wrapper. | Wrapper. Check which of the four HMMs cover the named genes. |
| **Transglycosidases** (GEL, PHR, GAS; Sps2, Dfg, Plb, Crh, Yps families) | Step 1 population (GPI-anchored enzymes) | Glyco_hydro_72 PF03198 (GEL/GAS family, to verify); Glyco_hydro_16 PF00722 (Crh-type, to verify) | Step 1 may call them. The 2026-09-27 review scored Gel1 (a *Coccidioides* enzyme) at 0 and called that arguably correct for an adhesin question. No enzyme-class label. | New label "remodeling enzyme" and a family table. |
| **Yapsins and other cell wall proteases** (Yps1, SAP9, SAP10) | Step 1 population | Asp PF00026 (generic, not yapsin-specific) | None. SAP9 is in the named panel as a negative. | A family table that separates yapsins from other aspartyl proteases. Needs curated examples. |
| **Chitinases, endo-β-1,3-glucanases, chitin deacetylases** | Secreted enzymes | Glyco_hydro_18 PF00704; Cellulase (GH5) PF00150; Polysacc_deac_1 PF01522 | None. | Family table. GH18 and GH5 are large families with intracellular members. |
| **Chitin synthases** (classes I to VII; Chs1, Chs2, Chs3, Chs6, Chs8; Mcs1) | Membrane synthase | Chitin_synth_1 PF01644, Chitin_synth_2 PF03142, Chitin_synth_1N PF08407 | None. | Family call is easy. The class I to VII split needs phylogenetic placement, not a Pfam hit. |
| **β-1,3-glucan synthase** (Fks, Gsc; regulator Rho1) | Membrane synthase | Glucan_synthase PF02364 | None. | Family call. Echinocandin resistance hot spots are a separate question. |
| **α-1,3-glucan synthases** (AGS1-3), **β-1,6-glucan** (KRE2, KRE5, KRE6, KRE9) | Membrane and Golgi enzymes | Kre9_KNH PF10342 for Kre9. Others not looked up. | None. | Pfam lookup for the others. |
| **Mannosyltransferases, galactofuran enzymes** (Ugm1, Gfsa), **GAG deacetylase** (ADG3) | Golgi and membrane enzymes | not looked up | None. | Lookup. Low priority. |
| **GPI anchor biosynthesis** (Gwt1, the target of fosmanogepix) | Endoplasmic reticulum membrane | GWT1 PF06423; Gpi1 PF05024; Gpi16 PF04113; Gaa1 PF04114; GPI2 PF06432 | None. | Family call. Relevant because a Gwt1 inhibitor changes the GPI-protein population that step 1 predicts. |
| **Melanin enzymes** (laccase; PKS alb1/pksP) | Mixed | Cu-oxidase PF00394 (generic multicopper oxidase) | None. Laccases are a planned hard-negative class (issue #14). | Family table. |
| **Wall integrity sensors and signaling** (Wsc, Mid2, Mtl1, Sho1, Sln1; Pkc1, Bck1, MAPKs, calcineurin, Hog1) | Membrane sensors; cytosolic kinases | WSC PF01822 | None. MSB2 and HKR1 are in the named panel as negatives. | Out of scope for v1 (cytosolic). The WSC domain also occurs in secreted enzymes, so a WSC hit alone is ambiguous. |
| **Septation and polarized growth** (Bni4, Chs4, septins, exocyst, polarisome) | Cytosolic | not looked up | None. | Out of scope. |
| **Immune and diagnostic targets** (Als1 and Als3 N-terminal regions, Mp65, Aspf3; GXM; anti-β-glucan antibodies) | Mixed | none | Antigen layer for *Coccidioides* only. Aspf3 is called an allergen in the 2017 summary. | Allergen layer (issue #19, parked). Antigen layer for other clades needs confounder genomes. |
| **Allergens** (Aspf3 named; Alternaria Alt a 1 and Aspergillus Asp f 4 have Pfam models AltA1 PF16541 and Allergen_Asp_f_4 PF25312) | Mixed (secreted and intracellular) | AltA1 PF16541; Allergen_Asp_f_4 PF25312 | None. No data directory. | Curation, then a homology-based module. Issue #19 lists the questions. |

## 4. What the tools can add to these questions

I list only what follows from the tool design and the measurements already made. None of it is
tested for these classes.

1. **A count and a list of the step 1 population per genome.** The papers cite counts of GPI
   proteins (66 in *S. cerevisiae*, over 100 in *Candida*). Step 1 produces such a list for any
   proteome. Its false-positive and false-negative rates are known only for yeast and Eurotiomycetes
   (`docs/model-review/STATUS.md`).
2. **A split of the surface population into adhesin-like, enzyme-like and other.** This needs the
   family table in section 3 and a rule for the order of calls. It does not exist.
3. **Family presence and copy number across the Fungi_5k proteomes.** A family scan by Pfam over
   488 or more proteomes is cheap. It answers "which clades have which wall-synthesis families".
   This repository has run the same kind of scan for PF28404 (`docs/reports/2026-10-02-pf28404-family.md`).
4. **Antigen and allergen candidates.** Only for *Coccidioides*, and only the antigen side.

## 5. What the tools cannot add

- Composition, cross-linking and masking of the wall. These are chemistry, not sequence.
- Environmental remodeling (the stress responses in the 2017 and 2023 pieces). The antigen layer
  uses spherule versus mycelium expression for *Coccidioides* only. No general expression module
  exists.
- Moonlighting surface proteins (Hsp60, gp43). The architecture documents place them out of scope.
- Sub-class assignment inside a large family (for example chitin synthase classes I to VII) from a
  single-sequence call.

## 6. Open points

1. I did not read the full text of the 2017 and 2023 papers. A reader should check the class list
   against the papers before it goes into a spec as a requirement.
2. Pfam mappings for KRE, AGS, mannosyltransferases and Ugm1 are not looked up.
3. Whether one Pfam model per class is specific enough is untested. Each family table row needs a
   specificity test with negative controls from the same genomes (compare the plan for the CFEM,
   Bys1 and hydrophobin scans).
