# Lectin class (host-interaction tier): design

*2026-10-10. Revision 1. Author: Claude Code (claude-sonnet-5-5). Task: `docs/agent-tasks/09-lectin-class-spec.md`. Owner decision of 2026-10-10: lectins become a new candidate class in the host-interaction tier, spec first. Status: draft for independent review (Opus). No code is written. `categories.yaml`, `family_table.tsv`, status files, `data/curated/adhesins/adhesins.tsv` and the agent task files are not changed. Every number below comes from a file or PMID named in the text. Items I could not verify are listed in section 10.*

## 1. The question

Which fungal proteins should the tool call "lectin" in the host-interaction tier, with which models, which truth, and which status? The answer in short:

1. The class can be defined by a small set of Pfam models that have fungal structure evidence (section 3.1). Most other lectin-like Pfam models are ubiquitous or have no fungal function evidence. They are shown, not called.
2. The truth set is `smoke` at best. Section 4.4 gives the cluster counts. No species reaches 20 independent positive clusters. In most species the status stays `unvalidated`.
3. The truth is partly circular for a family-defined class. Section 4.1 says how the spec limits this.
4. I found no fungal-side lectin paper for *Coccidioides*, *Histoplasma* or *Blastomyces* (section 2, F12). The discovery-pass families in Onygenales are associations, not functions.

## 2. Facts that this design depends on

Numbers were measured on 2026-10-10. "Local" files are in `_workdir/` of the main checkout and are not in git.

| # | Fact | Source |
|---|---|---|
| F1 | Active families that overlap lectin folds today: PA14 (PF07691) and GLEYA (PF10528), both in the narrow (adhesion) tier. PA14: 9 of 54 E1, 1 of 41 E2, 8 of 14 E3 adhesin positives; 0 N1 and 4 N2 hard negatives; 128 scan proteins in 32 proteomes. GLEYA: 12 of 54 E1, 1 of 41 E2, 3 of 14 E3; 0 N1, 1 N2; 42 scan proteins in 9 proteomes. | `docs/reports/data/sorting_hat/family_tests/family_tests.tsv` |
| F2 | The nine other lectin models tested in the family tests (Jacalin, Ricin_B_lectin, Fungal_lectin_2, RicinB_lectin_2, Laminin_G_3, GBD_Y3, Lectin_leg-like, LDL, Fungal_lectin) hit **none** of the 108 adhesin positives. The 34 hits of the 11-model file are GLEYA 16 and PA14 18. The 5 hits on hard negatives are PA14 on four *A. fumigatus* beta-glucosidases (Q4WU49, Q4WLY1, Q4WL79, Q4WA69) and GLEYA on the *V. dahliae* entry G2WSY1. | `_workdir/family_tests/lectin_adhesin_positives.domtbl`, `lectin_hard_negatives.domtbl`, `lectin_names.txt` (local); method in `analysis/family_tests/README.md` |
| F3 | In the shared hard-negative set the four PA14 glucosidases and G2WSY1 are **N2** rows (stratum `domain_non_member`). Per the owner decision of 2026-10-09 (Q5), only N1 rows count toward the 20-cluster floor. The four glucosidases are one cluster (HN0069); G2WSY1 is HN0018. | `data/controls/hard-negatives-shared/controls.tsv`; `docs/HANDOFF-2026-10-09-overnight.md` |
| F4 | Pfam 38.2 has 74 models whose name or description matches lectin, agglutinin, galectin or glycan binding (case-insensitive). The match is loose: it includes plectin repeats (substring "lectin"), viral and bacterial haemagglutinins and animal lectins. About 20 are plausible in fungi (section 3.1). | `_workdir/discovery/pfam_names.tsv` (local), my `grep` of it |
| F5 | Pool counts (653 proteomes, one per genus: 27 Onygenales, 7 Aspergillaceae, 619 background; pool S 436,120 proteins, pool U 422,071 proteins, so U is 96.8% of S). Proteins in pool S / pool U, and genera with the family (Onygenales / Aspergillaceae / background): Jacalin PF01419 308 / 307, 3 / 2 / 251; Ricin_B_lectin PF00652 270 / 270, 1 / 1 / 216; Fungal_lectin_2 PF18647 251 / 251, 15 / 3 / 152; RicinB_lectin_2 PF14200 146 / 146, 1 / 0 / 121; Laminin_G_3 PF13385 74 / 74, 1 / 0 / 53; GBD_Y3 PF22803 70 / 70, 4 / 0 / 47; Fungal_lectin PF07938 32 / 31, 6 / 0 / 22; LDL PF21691 34 / 30, 6 / 0 / 25; CVNH PF08881 560 / 558, 3 / 2 / 228; LysM PF01476 1,797 / 1,698, 26 / 7 / 478; Lectin_leg-like PF03388 45 / 45, 1 / 0 / 41; Calreticulin family PF00262 31 / 31, 1 / 0 / 30; Gal-bind_lectin PF00337 2 / 2; FB_lectin PF07367 0; CCL2-like PF21595 0. GLEYA 339 / 0 and PA14 175 / 0 (they are mechanisms, so they leave pool U by definition). I recounted these from the pool files; the first eight rows match the task file. | `_workdir/discovery/pools/*.pools.tsv.gz` (local), `analysis/secretome_discovery/pools.py`, `analysis/fungi_scan/selection_discovery.tsv`; recount script described in section 10 |
| F6 | Discovery pass, Onygenales pool U, q < 0.05: PF21691 6 of 27 genera, q 0.0135; PF07938 6 of 27, q 0.0163; PF18647 15 of 27, q 0.0261. Pfam-A in the pools is release 37.x; GA cutoffs here are from 38.2. The pass shows association with the clade. It makes no call and shows no function. | `docs/reports/data/sorting_hat/secretome_discovery/families_U_onygenales.tsv`; `docs/reports/2026-10-10-secretome-discovery-pass.md` |
| F7 | Pfam 38.2 gathering cutoffs (sequence = domain): Fungal_lectin PF07938 25.3; LDL PF21691 27; Fungal_lectin_2 PF18647 26.5; Jacalin PF01419 22; Ricin_B_lectin PF00652 29; RicinB_lectin_2 PF14200 25.6; Laminin_G_3 PF13385 27.6; GBD_Y3 PF22803 28; Lectin_leg-like PF03388 21.4; Gal-bind_lectin PF00337 22.1; FB_lectin PF07367 25; CCL2-like PF21595 27; CVNH PF08881 24.5; LysM PF01476 20.9 (domain 11.8); PA14 24; GLEYA 26.2. | `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm` (`GA` lines) |
| F8 | UniLectin: licence is Creative Commons Attribution 4.0 (CC BY 4.0) "to all copyrightable parts of UniLectin", with a disclaimer that some data may be covered by other rights. UniLectin3D home page (updated 2026-03-11): 3,005 3D structures (1,916 with interacting glycan), 794 distinct lectin sequences, 266 glycans, 1,246 articles. Access: web search and a public search API (POST `/api/search-module`, at most 100 documents per call). I found no bulk download link on the home page. | `https://unilectin.unige.ch/pages/license` and home page, fetched 2026-10-10 |
| F9 | UniLectin3D fungal entries (origin "Fungal and yeast lectins"). The phrase query reports **numFound 202** structures. Enumeration by class, species and family names gave **197 structures** (the API caps each call at 100 and ignores a start offset, so 5 are missing; probably in the capped classes galectin-like and Ricin-like). The 197 structures belong to **57 accessions** (56 UniProt, 1 UniParc ID without a UniProt FASTA), 19 UniLectin classes, 36 species names (35 fungal). One species (*Ralstonia solanacearum*, bacterial) is mislabelled as fungal origin. Published counts for comparison: 157 fungal structures of 1,740 (9%) as of 2018-08-29 (PMID 30239928); 194 fungal structures, 129 with a ligand (Lebreton 2021, abstract, PMID 34200153). | `analysis/lectin_class_spec/unilectin_fungal_sequences.tsv` (57 rows), `analysis/lectin_class_spec/unilectin_fungal_table.py` (method in its docstring) |
| F10 | Clusters of the 56 fungal UniLectin3D sequences (MMseqs2 17-b804f `easy-cluster --min-seq-id 0.3 -c 0.5`): **33 clusters**. Excluding the PA14 class, the yeast Emp L-type class and *Ralstonia* (the first is owned by an active call, the second is house-keeping, the third is bacterial): 43 sequences, **25 clusters**. Of these, 12 sequences in **8 clusters** are from Ascomycota (*Aleuria*, *Arthrobotrys*, *A. fumigatus*, *A. oryzae*, *Scedosporium*, *Cladosporium*, *Zymoseptoria*, *Sclerotinia*, *Gibberella*, *Magnaporthe*, *Neurospora*, *Tuber*). The most clusters in one species is 3 (*Coprinopsis*). *A. fumigatus* has 1 sequence (FleA, Q4WW81). | same table; `cluster30` column |
| F11 | Pfam GA hits on the 56 fungal UniLectin3D sequences (`hmmscan --cut_ga`, Pfam 38.2): 46 of 56 have at least one hit. Mapping of UniLectin3D class to Pfam model: AAL-like (FleA, AOL, AAL, SapL1, *Arthrobotrys* AOL, and one bacterial *Ralstonia* lectin) to Fungal_lectin PF07938 (6 of 6); Lyophyllum ginkbilobin-like to LDL PF21691; fungal fruit body lectin to FB_lectin PF07367 (5 of 5 with a sequence); galectin-like to Gal-bind_lectin PF00337 (4 of 4); Coprinus trefoil lectin-like to CCL2-like PF21595; cyanovirin-like to CVNH PF08881 (4 of 4); LysM-like to LysM PF01476; Ricin-like and the *Sclerotinia* trefoil lectin to RicinB_lectin_2 PF14200 (5 of 5); PA14 yeast adhesin to PA14 and GLEYA (10 of 10); yeast Emp L-type to Lectin_leg-like PF03388. **No Pfam model hits** the CML-like, CBM67-like, monocot-lectin-like and fungal prism classes, and 5 of 6 Laetiporus trefoil sequences (the sixth hits ETX_MTX2, a toxin model). Ten sequences have no hit at all. Tectonin hits PF19193 and the Fve (Ig-like) sequence hits PF09259. Hits to PVL-like propellers are FG-GAP_3, not a lectin model. **No fungal UniLectin3D sequence hits Jacalin, Ricin_B_lectin, Laminin_G_3, Fungal_lectin_2 or GBD_Y3.** | `analysis/lectin_class_spec/unilectin_fungal_sequences.tsv` (`pfam_ga_hits` column) |
| F12 | UniProt `Signal` feature on the 56 sequences (UniProt JSON, 2026-10-10): present in 9 of 10 PA14 sequences, 2 of 2 Emp L-type, 2 of 2 LysM, 1 of 1 CBM67-like; **absent** in all AAL-like sequences (including FleA Q4WW81), all fruit body, galectin, ricin-like, prism, trefoil and cyanovirin sequences. UniProt may omit a signal peptide that SignalP 6 would predict. I have not run SignalP on these sequences. | UniProt REST `accessions` JSON (`ft_signal`), fetched 2026-10-10 |
| F13 | FleA (*A. fumigatus*, UniLectin class AAL-like, family "AOL, AFL / FleA", Q4WW81) is already in `adhesins.tsv` with class `indirect_regulator`, evidence EXP, PMID 27058347. | `data/curated/adhesins/adhesins.tsv` (read only) |
| F14 | FleA literature (PubMed, abstracts read): Kerr 2016, PLoS Pathog 12(4):e1005555, PMID 27058347, [doi:10.1371/journal.ppat.1005555](https://doi.org/10.1371/journal.ppat.1005555), PMC4825926. FleA has twelve binding sites for fucosylated structures; protein and conidia bind lung mucin in a fucose-dependent way; macrophages take up ΔfleA conidia less well; mice infected with ΔfleA conidia had **more severe** pneumonia; the authors conclude FleA is not a virulence factor. Sakai 2019, Med Mycol 57(1):71-83, PMID 29370403, [doi:10.1093/mmy/myx163](https://doi.org/10.1093/mmy/myx163): fleA disruption gave earlier germination and **enhanced** pathogenicity in mice; six fucose sites contribute unequally; immune response in cultured cells differs for point mutants. Dussouy 2021, Org Biomol Chem 19(14):3234-3240, PMID 33885578, [doi:10.1039/d1ob00152c](https://doi.org/10.1039/d1ob00152c): calls the role of FleA "controversial"; reports a FleA co-crystal with a thiofucoside. The two knockout studies agree in direction. | PubMed |
| F15 | PA14 lectin-type adhesins: Veelders 2010, PNAS 107(52):22511-6, PMID 21149680, [doi:10.1073/pnas.1013210108](https://doi.org/10.1073/pnas.1013210108): Flo5A domain structures with high-mannose ligands; PA14-like fold with a DcisD calcium site. Ielasi 2012, Acta Cryst D 68:210-7, PMID 22349222, [doi:10.1107/S0907444911054898](https://doi.org/10.1107/S0907444911054898): Epa1p ligand domain with lactose, PA14 superfamily. Diderrich 2015, J Biol Chem 290(32):19597-613, PMID 26105055, [doi:10.1074/jbc.M115.655654](https://doi.org/10.1074/jbc.M115.655654): 17 Epa paralogs, most recognise glycans with terminal galactosides and confer epithelial cell adhesion. de Groot and Klis 2008, Mol Microbiol 68(3):535-7, PMID 18331471, [doi:10.1111/j.1365-2958.2008.06182.x](https://doi.org/10.1111/j.1365-2958.2008.06182.x): a commentary; a pentapeptide on a surface loop of the PA14 domain sets glycan specificity. | PubMed (abstracts) |
| F16 | Survey: Lebreton 2021, J Fungi 7(6):453, PMID 34200153, [doi:10.3390/jof7060453](https://doi.org/10.3390/jof7060453), PMC8227253. The abstract says 109 lectin sequence motifs from the UniLectin3D classification screened 1,223 MycoCosm species and gave **33,485** putative lectin sequences (MycoLec). The task file says 28,691; I could not reconcile the two (the abstract is the only part I re-read in this session). The predictions are not experimentally validated. Related: Bonnardel 2021, Nucleic Acids Res 49(D1):D1548-D1554, PMID 33174598, [doi:10.1093/nar/gkaa1019](https://doi.org/10.1093/nar/gkaa1019): UniLectin classification of 35 folds, 109 classes (at least 20% sequence similarity), 350 families (at least 70%). UniLectin3D itself: Bonnardel 2019, Nucleic Acids Res 47(D1):D1236-D1244, PMID 30239928, [doi:10.1093/nar/gky832](https://doi.org/10.1093/nar/gky832), PMC6323968 (full text read). | PubMed |
| F17 | Endemic dimorphic fungi. Searches "Coccidioides lectin" (32 hits), "Histoplasma lectin" (49 hits) and "Blastomyces BAD-1 adhesin" (4 hits) returned host-side work (Dectin-1, Dectin-2, mannose receptor; PMIDs 42171338, 32358020, 29610256, 24391211, 24379281, 35259948) and tools that use plant or animal lectins as probes (Grys 2016, J Proteome Res 15(10):3463-3472, PMID 27546806, [doi:10.1021/acs.jproteome.5b01054](https://doi.org/10.1021/acs.jproteome.5b01054): GSL-II and sWGA bind *C. posadasii* spherules and endospores). I read only the first 8 of 32 and first 8 of 49 hits, by abstract. **I found no paper that shows a fungal-encoded lectin in *Coccidioides*, *Histoplasma* or *Blastomyces*.** BAD-1 binds heparin (a glycosaminoglycan) through tandem repeats and WxxWxxW motifs, not through a lectin fold (Brandhorst 2013, PLoS Pathog 9(7):e1003464, PMID 23853587, [doi:10.1371/journal.ppat.1003464](https://doi.org/10.1371/journal.ppat.1003464); Beaussart 2015, mBio 6(5):e01403-15, PMID 26396244, [doi:10.1128/mBio.01403-15](https://doi.org/10.1128/mBio.01403-15)). BAD-1 belongs to the repeat call, not to this class. | PubMed |
| F18 | Mushroom lectin review: Erjavec 2012, Trends Biotechnol 30(5):259-73, PMID 22341093, [doi:10.1016/j.tibtech.2012.01.004](https://doi.org/10.1016/j.tibtech.2012.01.004). Abstract only (publisher page was not reachable in the earlier session). It covers higher fungi (Agaricomycetes). | PubMed |
| F19 | Tool structure. A Pfam family is a row of `data/sorting_hat/family_table.tsv` with a `module` name. A new module name (as `pfam_hydrophobin` and `pfam_hsba`) gives its own call and its own module status. A one-module call gets the module status; a call status needs two or more modules (hydrophobin spec decision D3). `docs/CLASSES.md`: calls are evidence calls, composites (`_candidate`) or `other_*`. Gated calls depend on R0. | `docs/CLASSES.md`; `docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md` section 3.2 |

## 3. Design

### 3.1 Which models define the class (task item 1)

Rule: a model joins the call (set M1) only if it has a fungal 3D structure in UniLectin3D or a paper-backed sugar-binding or host-binding result, **and** it is not ubiquitous in the secretome, **and** it does not belong to an active narrow-tier family or to a house-keeping class. Models that fail the rule go to a reported-only set (M2) or an exclusion set (M3). The cutoff for every model is the Pfam GA (F7). No custom cutoff is used in v1.

| Set | Model | UniLectin class (F11) | Evidence level | Why surface or host relevant | Pool S proteins (F5) |
|---|---|---|---|---|---|
| M1 | Fungal_lectin PF07938 | AAL-like (PropLec6A) | 3D structure (UniLectin3D: 5 fungal sequences from *A. fumigatus*, *A. oryzae*, *Scedosporium*, *Aleuria*, *Arthrobotrys*, plus one bacterial; ligand status not checked per entry; FleA has a thiofucoside co-crystal, PMID 33885578); host-binding and knockout papers for FleA (F14) | FleA binds host mucin and macrophage fucose (PMID 27058347). A host-recognition lectin, not shown to be an adhesin | 32 |
| M1 | LDL PF21691 | Lyophyllum ginkbilobin-like | 3D structure; the UniLectin record title names an alpha-galactosyl binding lectin (PMID 25504801 from the UniLectin record; not read by me), one mushroom species | Host relevance not shown. Included for the discovery-pass signal (F6) | 34 |
| M1 | Gal-bind_lectin PF00337 | galectin-like | 3D structure; mushroom (Coprinopsis, Agrocybe) | Pathogen relevance not shown | 2 |
| M1 | FB_lectin PF07367 | fungal fruit body lectin | 3D structure; mushroom | No pool S protein; rare in pathogens | 0 |
| M1 | CCL2-like PF21595 | Coprinus trefoil lectin-like | 3D structure; mushroom | No pool S protein | 0 |
| M2 (report only) | Fungal_lectin_2 PF18647 | none (no UniLectin3D fungal sequence hits, F11) | Pfam description "Putative glycan binding domain"; no fungal function evidence found | Enriched in Onygenales pool U (F6), 251 pool proteins, in 170 of 653 genera | 251 |
| M2 | GBD_Y3 PF22803 | none | Pfam description only | Pool 70, 51 genera | 70 |
| M2 | Jacalin PF01419 | none among fungal structures | Plant lectin family; fungal function not shown here | Ubiquitous: 256 of 653 genera | 308 |
| M2 | Ricin_B_lectin PF00652 | none among fungal structures | Same; the fold also occurs in enzymes (from memory, not checked) | Ubiquitous: 218 genera | 270 |
| M2 | RicinB_lectin_2 PF14200 | Ricin-like, *Sclerotinia* trefoil | 3D structure (ligand not checked) for mushroom and *Sclerotinia* | Ubiquitous: 122 genera; hits 5 of 5 UniLectin3D sequences of the classes | 146 |
| M2 | Laminin_G_3 PF13385 | none | "Concanavalin A-like lectin/glucanases superfamily" per Pfam description: mixed function | 54 genera | 74 |
| M2 | CVNH PF08881 | cyanovirin-like | 3D structure (4 ascomycete sequences, ligand not checked) | Ubiquitous: 233 genera; function in fungi not shown | 560 |
| M2 | LysM PF01476 | LysM-like | 3D structure (*Cladosporium* Ecp6, *Zymoseptoria* Mg1LysM); chitin binding, not a glycan-lectin in the narrow sense | Ubiquitous: 511 genera, 1,797 pool proteins. Would flood the call | 1,797 |
| M3 (exclude) | Lectin_leg-like PF03388 | yeast Emp L-type (Emp46p, Emp47p) | UniProt subcellular location for both: Golgi and ER membrane | House-keeping cargo receptors in the secretory pathway. Lebreton 2021 names ERGIC-53 and related ER lectins as house-keeping (task file reading of the full text; abstract does not state it) | 45 |
| M3 | Calreticulin family PF00262 | calnexin and calreticulin | Lebreton 2021 names them as house-keeping ER lectins; the *C. neoformans* calnexin CNE1 (J9VLH0) is already an N1 hard negative (PMID 34566931) | ER protein | 31 |
| M3 | Malectin PF11721 | not in the fungal UniLectin3D set | ER glycoprotein quality control (from memory, not checked) | 1 | 1 |
| Owned elsewhere | PA14 PF07691, GLEYA PF10528 | PA14 yeast adhesin | Active in `wall_family_domain`, narrow tier | The lectin class never re-calls them | 175, 339 |

Genus counts in the M2 rows are Onygenales + Aspergillaceae + background genera from F5 (for example Fungal_lectin_2: 15 + 3 + 152 = 170). The numbers of genera for Jacalin (3 + 2 + 251 = 256) and others are computed the same way.

**Ubiquitous classes.** Jacalin, Ricin-B, CVNH, LysM and Fungal_lectin_2 each occur in 26% to 78% of the 653 genera (LysM 511 of 653). A call on them would mark most fungi. The spec does not call them. It shows them as a count in the report (`lectin_like_families`) so that a reader can see the proteins. A discriminator that is not a model score (for example architecture: a protein that is mostly the lectin domain, or a lectin domain with a signal peptide and no enzyme domain) can promote an M2 model later. That needs its own spec and truth.

**M1 is small on purpose.** The sum of pool S proteins over the M1 models is 68 (overlap between models not checked) in 436,120 pool S proteins. This is a statement about the candidate class, not a finding about fungi: the models cover the structurally known lectins, which are mostly from mushrooms.

**Models without a Pfam entry.** Five UniLectin3D classes of fungal lectins have no Pfam hit (F11): CML-like, CBM67-like, monocot-lectin-like, fungal prism lectin (PhoSL-like) and Laetiporus trefoil lectin-like (5 of 6). UniLectin-derived HMMs for them would be new custom models (as the hydrophobin custom HMM work). v1 does not build any. Reason: all five classes are from mushrooms and none has pathogen or host evidence in what I read. The report states the gap.

### 3.2 Gate (task item 2)

Two forms, following the hydrophobin and repeat precedents:

- `lectin_domain` (evidence call, **ungated**): a hit to an M1 model at GA. It reads the module `pfam_lectin`. It does not need R0.
- `secreted_lectin_candidate` (composite, **gated**): `lectin_domain` and R0 signal peptide. It reads two modules, so it can have a call status.

Reason for no gate on the evidence call: F12. UniProt shows no signal peptide for FleA and for most structurally known fungal lectins. If SignalP 6 agrees, an R0 gate would remove FleA, the one protein with host-binding evidence. The pool design (signal peptide by construction) would never see it. Task L2 runs R0 on the 56 truth sequences. If R0 calls FleA, the gated form can become the default measured call. If not, the report says that the host-interaction lectin evidence is mostly in proteins without a predicted signal peptide. (FleA is described as expressed on conidia; whether it reaches the surface by a signal-peptide-independent route is not shown in the abstracts I read.)

No GPI or TM condition in v1. A mature TM helix is not a condition for lectin domains. Reuse the existing tool rules (R0 for the gate; TMHMM and GPI shown as evidence fields).

### 3.3 Relation to existing calls and double counting (task item 5)

| Existing call | Relation |
|---|---|
| `wall_family_domain` (PA14, GLEYA, CFEM, Flocculin, ...) | Disjoint by construction. M1 contains no model that is active in `family_table.tsv`. PA14 and GLEYA stay in the narrow tier. A protein with both a PA14 hit and an M1 hit is called by both calls, and the report shows both (this can happen only if a protein has two lectin domains; the F11 mapping shows none among the 10 PA14 structures). |
| `surface_attachment_candidate` | Not changed. It does not read `lectin_domain`. A lectin is not an attachment mechanism in the v1 definition unless the owner says so (decision D5). |
| `other_surface_no_mechanism`, `other_not_surface` | `lectin_domain` joins the mechanism lists (default, as for hydrophobin D2). A secreted protein with an M1 hit is then not "no mechanism". Effect on the discovery pools: pool U loses the M1 proteins (at most 31 + 30 + 2 = 63 of 422,071 from F5). The discovery tables stay valid because they used the run with 18 active families. |
| `tandem_repeat_protein` | Independent. A repeat protein with a lectin domain is called by both. |
| Antigen and allergen calls | Independent. |

**Narrow tier or host-interaction tier.** The split follows what the paper shows about the lectin, not about the domain:

- Adhesion-like (narrow): lectin domains that mediate cell-cell or cell-host attachment as a wall protein. PA14 and GLEYA (Flo, Epa; PMIDs 21149680, 22349222, 26105055). These stay where they are.
- Host-recognition (host interaction): soluble or conidial lectins that the host recognises or that bind host glycans without showing attachment as the function. FleA (F14). The class in this spec is this group. FleA is not called an adhesin, because both knockout studies found worse disease without it (F14).
- Unknown role: every other M1 model. The call reports "lectin domain". The tier label for them is "host interaction (candidate)", not a measured role.

### 3.4 Calls, modules and what changes in the config (stated, not done)

| Item | Content |
|---|---|
| Family table rows | M1 models get `module=pfam_lectin`, `class` text "host-interaction lectin", `source_pmid` from section 2, `specificity_note` with pool counts, `active=no` until the owner activates (decision D1). M2 models get `module=pfam_lectin_like`. M3 models are not rows. |
| Modules | `pfam_lectin` (M1), `pfam_lectin_like` (M2, report only, read by no call). Same wrapper as `pfam_hydrophobin`. Per-module family digest (hydrophobin task H1c) must exist first, or adding rows changes the identity of every Pfam module. |
| Calls | `lectin_domain` = `{flag: pfam_lectin.hit}`; `secreted_lectin_candidate` = `lectin_domain` and `signal_peptide_protein[R0]`. |
| Category | host interaction and immune evasion (tier from the owner decision of 2026-10-09). The current `categories.yaml` has five categories; M2/workstream A decides the category name. This spec does not pick a name. |
| Effect | `config_sha256` changes, so existing call status files are stale (hydrophobin spec 3.5). Batch with workstream A so the repeat call is re-measured once. |

## 4. Truth set (task item 3)

### 4.1 Circularity and tiers

The class is defined by Pfam models. Three problems follow:

1. **A positive that is "found by the model" is not evidence.** The rule of the task is kept: a protein is not truth because a Pfam model hits it.
2. **UniLectin3D classes and Pfam models share sequences.** Pfam seeds may contain the same structure sequences (not checked). Leakage for all structure-derived rows is `in_reference` until checked.
3. **A recall measured on structure-defined proteins measures how well the model covers the structure it was built from.** It is not recall on unknown proteins. The report says this beside every number.

Tiers (follow the hydrophobin tiers):

| Tier | Evidence | Use |
|---|---|---|
| T1 | A paper reports that the protein binds a glycan (structure with ligand, ITC, glycan array or inhibition by a sugar) **and** the sequence is a fungal protein with a UniProt accession. The PMID is recorded. | primary positives |
| T2 | UniLectin3D curated lectin structure without a recorded ligand, or a Swiss-Prot reviewed entry with experimental sugar-binding annotation (ECO:0000269) | primary positives |
| T3 | Family or name by sequence similarity only | secondary, reported apart, never primary |
| T4 | Review judgement from sequence, domain and architecture evidence without a paper | reported only, never truth, never used for sensitivity |
| H flag | A paper reports host binding or a host-relevant phenotype for the protein (example: FleA, PMIDs 27058347 and 29370403). The flag is an attribute of a T1 or T2 row. It does not change the tier. | host-interaction column in the report |

### 4.2 Sources

1. UniLectin3D (CC BY 4.0, F8). The enumeration used 57 accessions (F9). Licence allows reuse with attribution; the output table cites the three UniLectin papers. A full export would need the owner to ask the UniLectin team or to confirm that the API use is acceptable (I found no stated rate limit or bulk download).
2. UniProt reviewed entries with sugar-binding evidence (task L1 counts them; I did not).
3. PubMed: papers with glycan-binding or host-binding experiments for fungal proteins in the species of the status entries. Found so far: FleA (F14). None for the endemic fungi (F17).
4. `data/curated/adhesins/adhesins.tsv` (read only): FleA is there as `indirect_regulator` (F13). The lectin truth may reuse the row's PMID. It does not change the row.

### 4.3 Negatives

| Negative | Source | Status |
|---|---|---|
| ER house-keeping lectin-fold proteins | *C. neoformans* CNE1 (J9VLH0), already in the shared set as N1 (`secretory_non_surface`, PMID 34566931); yeast Emp46p and Emp47p (Q12396, P43555), localisation from UniProt record only, **no paper read**; calnexin and calreticulin orthologs in the six Phase C species (task L3 finds the IDs and sources) | N1 for CNE1. Emp46/47 need a paper before they are N1. |
| Enzymes with a lectin-like fold | The four PA14 *A. fumigatus* glucosidases (F3) are N2 and hit PA14, which is not in M1. They do not test specificity of `lectin_domain`. Other enzymes with a carbohydrate-binding module of the M2 models are not yet identified (L3). I do not cite any from memory. | none usable now |
| Per-species population | As in the hydrophobin spec 4.3: reviewed or curated secreted proteins of the species minus any lectin label, labelled **assumed** | needs owner confirmation (D4) |
| Shared hard negatives | `data/controls/hard-negatives-shared/`: 137 N1 rows in 76 clusters, none above 16 per species (F3 and handoff recount). No row is a lectin-domain look-alike except CNE1. | reuse, N1 only toward the floor |

### 4.4 Cluster counts and the honest answer (task item 3)

| Level | Positives | Clusters at 30% | Against the 20-cluster floor |
|---|---|---|---|
| All 56 fungal UniLectin3D sequences | 56 sequences | 33 (F10) | above 20, but this is a sequence table, not a species entry, and it includes bacterial, ER and PA14 rows |
| After removing PA14, Emp L-type and *Ralstonia* | 43 | 25 | above 20, mostly Agaricomycetes (mushroom) |
| Ascomycota only | 12 | 8 | below 20 |
| One species, best case | 1 to 3 | 1 to 3 (*Coprinopsis* 3; *A. fumigatus* 1) | far below 20 |
| Onygenales (*Coccidioides*, *Histoplasma*, *Blastomyces*) | 0 | 0 | none |

Statements:

1. **The realistic status is `smoke` at most.** `estimated` needs 20 positive and 20 negative clusters in one species entry and an interval half-width of at most 0.10 (`COMMON-RULES.md` rule 7). No species has more than 3 positive clusters. Pooled results across species are report tables, not status entries (hydrophobin decision D1).
2. **Most entries will be `unvalidated`.** Only *A. fumigatus* (FleA, 1 cluster) has a positive in a Phase C species that I can name. The status rules decide whether one cluster gives `smoke`; this spec does not assume it. Specificity-only entries (no positives) get `unvalidated` (sowgp spec section 4).
3. **N1 negatives per species are also short.** The largest N1 count in one species is 16 clusters (*S. cerevisiae*, handoff recount), below 20.
4. The pooled 25-cluster sequence table has a use: it measures how many fungal lectin classes the M1 models cover (F11). That is a coverage statement, not a sensitivity for a species.

## 5. Curation workflow (task item 4)

Pattern: `analysis/hydrophobin_truth/curation_rubric.md` (written before any decision; evidence lines; decision rules in order; outcomes `lectin (T4)`, `not lectin (T4)`, `unresolved`; T4 never in truth; the owner can overrule).

### 5.1 Evidence lines (written before any decision)

| Line | Meaning | Independent of the Pfam call? |
|---|---|---|
| A | A 3D structure with a glycan ligand exists for the protein or for an ortholog at 90% identity or more over 90% of length (UniLectin3D or PDB) | yes |
| B | Hit to an M1 model at GA | no |
| C | A paper reports glycan binding or host binding for this gene product (PMID, quote) | yes |
| D | Reciprocal best hit at 40% identity and 70% coverage of both to a T1 or T2 protein | partly |
| E | Other function: a non-lectin domain covers most of the sequence, or the protein is annotated ER or Golgi house-keeping with a source | yes |
| F | Signal peptide by R0 | yes (not a lectin line; used for the gated call and the report) |

### 5.2 Decision rules (order fixed)

1. `lectin (T4)` if C holds, or if A holds with B or D.
2. `not lectin (T4)` if E holds and C does not hold.
3. `unresolved` otherwise (only B, only D, or conflicting evidence).

### 5.3 What an agent can do and what needs the owner

| Task | Agent | Owner |
|---|---|---|
| Enumerate UniLectin3D fungal entries, map to UniProt, cluster, run Pfam and R0 (done in part, F9 to F11) | yes | none |
| Fetch UniProt evidence codes, find reviewed sugar-binding entries, map to the Phase C proteomes by sequence (at least 95% identity over 90% of length) | yes | none |
| PubMed search per species for glycan-binding and host-binding papers; quote the sentence | yes | none |
| Apply the rubric to unlabelled hits in the species proteomes (T4) | yes | may overrule |
| Choose M1, M2 and M3 membership | proposal only | decides (D1, D2) |
| Accept "assumed" negatives (D4) | no | decides |
| Decide tier label of each role (adhesion, host recognition, unknown) | proposal only | decides |
| Contact UniLectin about a bulk export, if wanted | no | decides |

## 6. Measurement and report (task item 6)

### 6.1 Calls and status

| Call | Reads | Status measured |
|---|---|---|
| `lectin_domain` | `pfam_lectin` (one module) | module status of `pfam_lectin` (a one-module call has no call status, as for `hydrophobin_domain`) |
| `secreted_lectin_candidate` | `pfam_lectin`, R0 | derived (composite) |
| M2 families | `pfam_lectin_like`, read by no call | none; counts only |

Status per species entry: `unvalidated` (no truth in the species), `smoke` (some truth, too little for a number), `estimated` (not reachable, section 4.4). A status file records leakage `in_reference` for structure-derived rows. Truth tables use the species and taxon rules of `COMMON-RULES.md` rule 5 and one species per entry.

### 6.2 What the report shows (per proteome and per run)

1. Count of `lectin_domain` proteins, and of `secreted_lectin_candidate`, with the family name of each protein.
2. A separate line "lectin-like families (not called)": count per M2 family. The reader sees PF18647, Jacalin and others without a call.
3. For each call the status line and the leakage statement, and in words: "the truth is structure-derived; recall measures model coverage of known structures".
4. A sequence-level table (pooled, no status): the 56 UniLectin3D sequences, their classes, M1 hits and R0 result; the five classes with no Pfam model (section 3.1).
5. T4 outcomes, in separate columns, never in sensitivity.
6. For *Coccidioides*, *Histoplasma* and *Blastomyces*: "no fungal lectin literature found by the searches in F17; calls here are domain evidence only".

### 6.3 Discovery-pass families

| Family | Handling |
|---|---|
| PF07938 | In M1. The enrichment (6 of 27 Onygenales genera, q 0.0163, F6) is shown beside the call as an association with the clade. No function is claimed for Onygenales members. FleA-orthology is a sequence statement (Pfam hit at GA). |
| PF21691 | In M1 for the same reason (6 of 27, q 0.0135). The structure is from a mushroom; the call says "alpha-galactosyl lectin domain". |
| PF18647 | In M2: 15 of 27 genera (q 0.0261) but no fungal structure and no function evidence (F11). Shown as a count and an enrichment. Promotion to M1 needs a structure or a paper (decision D2). |

## 7. Decisions for the owner (stop points)

D1 and D2 must be decided before L4 (family table rows). Defaults are my recommendation. None is data-driven where I say so.

| # | Question | Default |
|---|---|---|
| D1 | Is v1 M1 = PF07938, PF21691, PF00337, PF07367, PF21595 (all with fungal structure evidence), with the rows inactive until a review of the first run? | Yes. The choice rests on structure evidence only, not on any measured performance. |
| D2 | M2 models (PF18647, PF22803, PF01419, PF00652, PF14200, PF13385, PF08881, PF01476): reported only, or in the call? | Reported only. Reason: ubiquity (F5) and no fungal function evidence. |
| D3 | M3 (house-keeping: PF03388, PF00262, PF11721) excluded from the call and reported as a note? | Excluded. The Emp46/47 localisation needs a paper first (4.3). |
| D4 | Accept per-species "assumed" negatives for this class (as D8 of the hydrophobin spec)? | Yes, labelled "assumed". |
| D5 | Does a lectin domain count as a mechanism for `surface_attachment_candidate`? | No. Only for the `other_*` mechanism lists. |
| D6 | Gate: ungated evidence call plus gated composite (section 3.2)? | Yes. Revisit after L2 shows R0 on FleA. |
| D7 | Role label: "host recognition" for FleA-type lectins and "host interaction (candidate)" for other M1; PA14 and GLEYA stay narrow? | Yes. |
| D8 | Use the UniLectin3D API for a one-time enumeration (CC BY 4.0, attribution in the table), or ask for a bulk export? | API enumeration, with attribution. No bulk export request without the owner. |
| D9 | Custom HMMs for UniLectin classes with no Pfam model (five classes, 3.1)? | Not in v1. |
| D10 | Is a lectin domain shown in the report as an antigen or allergen angle? | No. |

## 8. Owner questions, one at a time

1. **Is M1 right for v1?** (D1, D2.) Recommended answer: yes, M1 = five structure-backed models; M2 reported only. Reason: the class stays small and every model has fungal structure evidence; the ubiquitous models would mark 26% to 78% of genera.
2. **Gate.** (D6.) Recommended: ungated evidence call plus a gated composite, and decide after L2 whether FleA has an R0 signal peptide.
3. **Do you accept that the status will be `smoke` at most and mostly `unvalidated`?** (Section 4.4.) Recommended: yes. More truth would need new experiments in the endemic fungi, not more curation of mushroom lectins.
4. **Role of FleA.** (D7.) Recommended: label FleA a host-recognition lectin, not an adhesin; both knockout studies found worse disease without it. `adhesins.tsv` keeps its `indirect_regulator` row; the lectin truth uses its PMID without editing it.
5. **Should I stop here, or fund a literature pass on lectins in the endemic fungi?** Recommended: one more targeted search (L3b) before any claim, because my search covered only the first 8 hits of two of the three queries and read abstracts only (F17).

## 9. Risks

1. **Truth is structure-derived and mostly from mushrooms.** A recall number says little about pathogen lectins. Stated in the report.
2. **Circularity with Pfam.** Leakage `in_reference`; status capped at `smoke`.
3. **Gate.** If SignalP agrees with UniProt (F12), the host-interaction lectins that matter are not in the secretome pools of the discovery pass. The discovery enrichments for PF07938 and PF21691 then describe the small subset with a signal peptide.
4. **Pfam release mismatch.** Pools use Pfam-A 37.x, GA values 38.2 (F6).
5. **UniLectin3D enumeration is a lower bound** (197 of 202, F9). One accession has no sequence. Counts can change.
6. **The mislabel** of *Ralstonia* (F9) shows that the origin field needs a check per row.
7. **Pool U is 96.8% of pool S** (F5): "no mechanism" is a weak filter, so pool U is not a test of whether a lectin family is unexplained.
8. **Role over-reading.** A lectin domain does not show the protein binds a host glycan or triggers inflammation. Lebreton 2021 states that its predictions are not experimentally validated (F16).

## 10. What I could not verify, and assumptions

- **Not read in full text:** Kerr 2016, Sakai 2019, Dussouy 2021, Veelders 2010, Ielasi 2012, Diderrich 2015, Brandhorst 2013, Beaussart 2015, Grys 2016, Erjavec 2012, Lebreton 2021 (only abstracts, in this session). UniLectin3D 2019 full text was read. PMC full text was not opened for the others.
- **Lebreton 2021 numbers:** abstract 33,485 putative lectins; task file 28,691. Not reconciled. The statement that Lebreton names calnexin, calreticulin and ERGIC-53 as house-keeping comes from the task file, not from my reading.
- **PMID 25504801** (Lyophyllum structure) comes from the UniLectin API record; I did not look it up in PubMed.
- **UniLectin3D completeness:** 197 of 202 structures enumerated. No bulk download found. Licence read from the site (CC BY 4.0); I did not read the terms-of-use page or check for an API rate limit.
- **FleA signal peptide:** not computed (F12). SignalP was not run.
- **Pfam seed overlap** with UniLectin3D structures: not checked.
- **From memory, not checked:** that the Ricin-B fold occurs in enzymes; that malectin is an ER glycoprotein-quality-control protein; Emp46/47 as cargo receptors (only their UniProt location is read). They appear in the text with that flag.
- **Pool recount:** I recounted F5 with a throwaway script (read `pools.tsv.gz` of the 653 proteomes, split `pfam` on `;`, count by pool and by group); the script is not committed. The first eight rows equal the task file. Local files are not in git, so a reviewer on another machine cannot repeat it without the main checkout.
- **Search limits:** PubMed queries for the endemic fungi were few and abstract-level. A negative result here is "not found by these searches", not "does not exist".
- **Assumption:** that "lectin" in this tool means a Pfam-defined glycan-binding domain. A protein that binds glycans without such a model (BAD-1 and heparin, F17) is not in the class.
- **No SLURM jobs were run.** I ran MMseqs2 and `hmmscan` on 56 sequences on the login node.

## 11. Work plan

| Task | Content | Check |
|---|---|---|
| L0 | Owner answers D1 to D3 and questions 1 to 5 | decisions recorded in this file |
| L1 | Commit the enumeration as a reproducible table with `unilectin_fungal_table.py`; add UniProt evidence codes and reviewed status for the 56 accessions; complete the 5 missing structures (query by class and species in smaller slices) | table with ECO codes; count of T1, T2, T3 rows by species and cluster |
| L2 | Run R0 (SignalP 6) and TMHMM on the 56 sequences with the tool's own wrappers (a small run; check the queue first) | FleA R0 result known; D6 revisited |
| L3 | Negatives: IDs and sources for calnexin, calreticulin and Emp46/47 orthologs in the Phase C species; look for enzymes with CBM13-type or other M2 domains with a paper; map to proteomes by sequence | N1 rows with PMID; counts of clusters |
| L3b | Targeted literature pass: fungal lectins in Onygenales, *Aspergillus*, *Candida*; read the full text of Kerr 2016 and Sakai 2019 | list with PMID and quote; `unverified_notes.md` |
| (stop) | Owner reviews `counts.md` and `unverified_notes.md` | |
| L4 | Family table rows (inactive), module `pfam_lectin` and `pfam_lectin_like`, calls; tests first; after workstream A decides category names; per-module family digest first | tests pass; `config_sha256` change noted |
| L5 | Run `hmmsearch --cut_ga` on the Phase C proteomes and the discovery pool proteomes; counts per proteome | module tables |
| L6 | Truth tables per species, status files; apply the rubric to unlabelled M1 hits (T4) | status files; report |
| L7 | Docs: `docs/CLASSES.md` entries (generated), paper notes | class-table test passes |

Each task is its own commit. An independent review of this spec (revision 1) precedes L1.

## 12. Sources and files

Files read: `docs/agent-tasks/README.md`, `COMMON-RULES.md`, `09-lectin-class-spec.md`; `docs/CLASSES.md`; `docs/HANDOFF-2026-10-09-overnight.md`; `docs/reports/2026-10-10-secretome-discovery-pass.md`; `analysis/family_tests/README.md`; `analysis/hydrophobin_truth/curation_rubric.md`; `docs/superpowers/specs/2026-10-09-sowgp-like-calls-design.md`; `docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md`; `data/sorting_hat/family_table.tsv`; `data/controls/hard-negatives-shared/controls.tsv`; `data/curated/adhesins/adhesins.tsv` (read only). Local, read only, in the main checkout: `_workdir/discovery/pools/`, `_workdir/discovery/pfam_names.tsv`, `_workdir/family_tests/`. New files from this task: `analysis/lectin_class_spec/unilectin_fungal_table.py` and `unilectin_fungal_sequences.tsv`. The intermediate files (API responses, sequences, MMseqs2 and hmmscan outputs) were in `$SCRATCH` and are not committed.
