# Lectin class (host-interaction tier): design

*2026-10-10. Revision 2. Author: Claude Code (claude-sonnet-5-5). Task: `docs/agent-tasks/09-lectin-class-spec.md`. Owner decision of 2026-10-10: lectins become a new candidate class in the host-interaction tier, spec first. Status: draft, revised after independent review 1 (`2026-10-10-lectin-class-design-review-1.md`, Opus). A second independent review of this revision is not done. No code is written. `categories.yaml`, `family_table.tsv`, status files, `data/curated/adhesins/adhesins.tsv` and the agent task files are not changed. Every number below comes from a file or PMID named in the text. Items I could not verify are listed in section 10.*

## Changes from revision 1

| Review finding | What I did |
|---|---|
| Required 1 (H1): FleA role argument invalid; papers missing | Rewrote F14 and section 3.3. A disease outcome does not show whether a protein mediates attachment. The spec now states the dispute (Houser 2013, Kerr 2016, Sakai 2019, Richard 2018, Duca 2022, SapL1 2021, Dussouy 2021). The host-interaction label rests on host-glycan binding and immune effects. "The one protein with host-binding evidence" is removed. D7 and owner question 4 rewritten. Houser 2013 read in full text; the others at abstract level. |
| Required 2 (H2): paracoccin missed | Added F20. Paracoccin is encoded by PADG_03347 (UniProt C1G842). I ran `hmmscan --cut_ga` on it: the only hit is Glyco_hydro_18 (PF00704). It hits no lectin model, so M1 does not see it. Section 1, 4.4 and 6.2 changed from "none" to "one, not reachable by M1". |
| Required 3 (H3): LysM and CVNH placement | Stated one rule for M1 and M2 (section 3.1): host-interaction evidence level H and ubiquity U. Ecp6 (LysM) and MoCVNH papers added (F21). LysM stays out of the call for ubiquity (78% of genera). The rule moves four revision-1 M1 models to M2, because they have no host evidence (section 3.1). This changes the default of D1. D2 gets the architecture-defined LysM sub-call option. |
| Required 4 (H4): M1 coverage not stated | Stated in F10, 4.4 and owner question 3. Five revision-1 M1 models: 16 of 43 sequences, 7 of 25 clusters; Ascomycota 5 of 12 sequences, 3 of 8 clusters. PF07938 alone: 5 of 43 sequences, 3 of 25 clusters. I reproduced both counts with the committed `cluster_check.py`. |
| Required 5 (M1): `smoke` from one in-reference positive | Stated in 4.4 statement 2. The code gives `smoke` for one positive with a sensitivity (`status.py` `status_from_measure`, read in this revision). Leakage other than `none` caps at `smoke` (`call_status.py`). The *A. fumigatus* number would be 1 of 1 on a protein that may be in the model's reference set. The option "no status file" is added to D4b. |
| Required 6 (M2): no negative tests M1 specificity | Stated in 4.3. I scanned the 137 N1 negatives and six Phase C proteomes (F23). 0 of 137 N1 hit an M1 model; the N1 set was not built to contain look-alikes, so this says little. Task L3 now says: build look-alike negatives. |
| Required 7 (M3): Emp46/47 paper | Cited Sato and Nakano 2002 (PMID 12134087, F22). Abstract read. "No paper read" removed. |
| Required 8 (M4): "FleA-orthology"; Aspergillaceae 0 of 7 | Wording removed. Section 6.3 now says family membership, reports 0 of 7 Aspergillaceae genera, and states that pool members have a signal peptide by construction and FleA does not. |
| Required 9 (M5): D5 precedent | Stated in 3.3 and D5. `surface_attachment_candidate` already reads hydrophobin and HsbA (`docs/CLASSES.md`). D5 is shown as an open question with both options. |
| Required 10 (M6): ungated call and surface location | Added to 3.2 and 6.2 as fixed report text. |
| Required 11 (L1 to L7) | L1: F1 now cites `data/sorting_hat/family_table.tsv`. L2: Lebreton count resolved in F16 (33,518 total, 28,691 filtered; abstract 33,485); house-keeping statement cites the paper. L3: the per-module family digest exists (`_family_digest` in `modules/cli.py`); dependency dropped. L4: recount script, FASTA, clusters and outputs committed (section 12). L5: the discovery-pass reason for LDL is removed. L6: P-type lectin (MRH domain) added to M3. L7: wording removed ("small on purpose", "every model has fungal structure evidence", the opinion on new experiments). |
| Work-plan notes | Per-module digest dependency dropped. Added N1 and six-proteome scan (F23, task L2b, done as a login-node run). L3b gets the new papers. D7 rewritten. |
| New searches (task) | Added F21 (Beauveria BbLec1, M. oryzae MoCVNH structure paper 2015, Ielasi 2016 Als/Epa, Kock 2018 PA14 yeasts). Searches for *Coccidioides*, *Histoplasma* and *Blastomyces* found host-side papers only. |
| Disagreements | None with a finding. One adjustment: the review proposes keeping M1 and stating the reason. I applied the stated rule to every model instead. The result is a smaller M1 (one model) and a different D1 default. The owner can restore the five-model M1 (D1 option B). |

## 1. The question

Which fungal proteins should the tool call "lectin" in the host-interaction tier, with which models, which truth, and which status? The answer in short:

1. The call can be defined by one Pfam model with host-interaction evidence and low ubiquity: Fungal_lectin PF07938 (AAL-like fucose-binding lectins, section 3.1). Other lectin-like models are shown, not called. This is a smaller class than in revision 1.
2. The truth set is `smoke` at most. Section 4.4 gives the cluster counts. No species reaches 20 independent positive clusters. In most species the status stays `unvalidated`.
3. The truth is partly circular for a family-defined class. Section 4.1 says how the spec limits this.
4. A published host-binding lectin exists in Onygenales: paracoccin of *Paracoccidioides* (F20). The Pfam lectin models do not detect it. I found no fungal-lectin paper for *Coccidioides*, *Histoplasma* or *Blastomyces* (F17, not found by these searches). The discovery-pass families in Onygenales are associations, not functions.
5. The role of FleA is disputed in the literature (F14). The spec does not decide it.

## 2. Facts that this design depends on

Numbers were measured on 2026-10-10. "Local" files are in `_workdir/` of the main checkout and are not in git.

| # | Fact | Source |
|---|---|---|
| F1 | Active families that overlap lectin folds today: PA14 (PF07691) and GLEYA (PF10528), both in the narrow (adhesion) tier. Both are `active=yes` in the family table (GLEYA since 2026-10-09; PA14 has `second_condition=signal_peptide`). The family-test file `family_tests.tsv` has `active=no` for GLEYA because it predates the activation. Family tests: PA14 9 of 54 E1, 1 of 41 E2, 8 of 14 E3 adhesin positives; 0 N1 and 4 N2 hard negatives; 128 scan proteins in 32 proteomes. GLEYA 12 of 54 E1, 1 of 41 E2, 3 of 14 E3; 0 N1, 1 N2; 42 scan proteins in 9 proteomes. | `data/sorting_hat/family_table.tsv` (active state); `docs/reports/data/sorting_hat/family_tests/family_tests.tsv` (counts) |
| F2 | The nine other lectin models tested in the family tests (Jacalin, Ricin_B_lectin, Fungal_lectin_2, RicinB_lectin_2, Laminin_G_3, GBD_Y3, Lectin_leg-like, LDL, Fungal_lectin) hit **none** of the 108 adhesin positives. The 34 hits of the 11-model file are GLEYA 16 and PA14 18. The 5 hits on hard negatives are PA14 on four *A. fumigatus* beta-glucosidases (Q4WU49, Q4WLY1, Q4WL79, Q4WA69) and GLEYA on the *V. dahliae* entry G2WSY1. | `_workdir/family_tests/lectin_adhesin_positives.domtbl`, `lectin_hard_negatives.domtbl`, `lectin_names.txt` (local); method in `analysis/family_tests/README.md` |
| F3 | In the shared hard-negative set the four PA14 glucosidases and G2WSY1 are **N2** rows (stratum `domain_non_member`). Per the owner decision of 2026-10-09 (Q5), only N1 rows count toward the 20-cluster floor. The four glucosidases are one cluster (HN0069); G2WSY1 is HN0018. N1: 137 rows in 76 clusters; most clusters in one species: *S. cerevisiae* 16, *A. fumigatus* 15. | `data/controls/hard-negatives-shared/controls.tsv`; `docs/HANDOFF-2026-10-09-overnight.md` |
| F4 | Pfam 38.2 has 74 models whose name or description matches lectin, agglutinin, galectin or glycan binding (case-insensitive). The match is loose: it includes plectin repeats (substring "lectin"), viral and bacterial haemagglutinins and animal lectins. About 20 are plausible in fungi (section 3.1). | `_workdir/discovery/pfam_names.tsv` (local) |
| F5 | Pool counts (653 proteomes, one per genus: 27 Onygenales, 7 Aspergillaceae, 619 background; pool S 436,120 proteins, pool U 422,071 proteins, so U is 96.8% of S). Recounted with a committed script. Per family: pool S proteins / pool U proteins; genera with the family in Onygenales / Aspergillaceae / background; total genera and share of 653. Jacalin PF01419 308 / 307; 3 / 2 / 251; 256 (39%). Ricin_B_lectin PF00652 270 / 270; 1 / 1 / 216; 218 (33%). Fungal_lectin_2 PF18647 251 / 251; 15 / 3 / 152; 170 (26%). RicinB_lectin_2 PF14200 146 / 146; 1 / 0 / 121; 122 (19%). Laminin_G_3 PF13385 74 / 74; 1 / 0 / 53; 54 (8%). GBD_Y3 PF22803 70 / 70; 4 / 0 / 47; 51 (8%). Fungal_lectin PF07938 32 / 31; 6 / 0 / 22; 28 (4%). LDL PF21691 34 / 30; 6 / 0 / 25; 31 (5%). CVNH PF08881 560 / 558; 3 / 2 / 228; 233 (36%). LysM PF01476 1,797 / 1,698; 26 / 7 / 478; 511 (78%). Lectin_leg-like PF03388 45 / 45; 1 / 0 / 41; 42. Calreticulin PF00262 31 / 31; 1 / 0 / 30. Gal-bind_lectin PF00337 2 / 2; 0 / 0 / 2. FB_lectin PF07367 0. CCL2-like PF21595 0. Malectin PF11721 1 / 1. GLEYA 339 / 0 and PA14 175 / 0 (mechanisms, so they leave pool U by definition). The first eight rows equal the task file. | `analysis/lectin_class_spec/pool_recount.py` and `pool_recount.tsv` (committed; the input `_workdir/discovery/pools/*.pools.tsv.gz` is local) |
| F6 | Discovery pass, Onygenales pool U, q < 0.05: PF21691 6 of 27 genera, q 0.0135; PF07938 6 of 27, q 0.0163; PF18647 15 of 27, q 0.0261. Pfam-A in the pools is release 37.x; GA cutoffs here are from 38.2. The pass shows association with the clade. It makes no call and shows no function. | `docs/reports/data/sorting_hat/secretome_discovery/families_U_onygenales.tsv`; `docs/reports/2026-10-10-secretome-discovery-pass.md` |
| F7 | Pfam 38.2 gathering cutoffs (sequence = domain): Fungal_lectin PF07938 25.3; LDL PF21691 27; Fungal_lectin_2 PF18647 26.5; Jacalin PF01419 22; Ricin_B_lectin PF00652 29; RicinB_lectin_2 PF14200 25.6; Laminin_G_3 PF13385 27.6; GBD_Y3 PF22803 28; Lectin_leg-like PF03388 21.4; Gal-bind_lectin PF00337 22.1; FB_lectin PF07367 25; CCL2-like PF21595 27; CVNH PF08881 24.5; LysM PF01476 20.9 (domain 11.8); PA14 24; GLEYA 26.2. PF07938 has 3 seed sequences (NSEQ 3); LDL 14, FB_lectin 12, CCL2-like 17 (review 1 check). | `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm` (`GA`, `NSEQ` lines) |
| F8 | UniLectin: licence is Creative Commons Attribution 4.0 (CC BY 4.0) "to all copyrightable parts of UniLectin", with a disclaimer that some data may be covered by other rights. UniLectin3D home page (updated 2026-03-11): 3,005 3D structures (1,916 with interacting glycan), 794 distinct lectin sequences, 266 glycans, 1,246 articles. Access: web search and a public search API (POST `/api/search-module`, at most 100 documents per call). I found no bulk download link on the home page. Review 1 did not re-check this. | `https://unilectin.unige.ch/pages/license` and home page, fetched 2026-10-10 |
| F9 | UniLectin3D fungal entries (origin "Fungal and yeast lectins"). The phrase query reports **numFound 202** structures. Enumeration by class, species and family names gave **197 structures** (the API caps each call at 100 and ignores a start offset, so 5 are missing; probably in the capped classes galectin-like and Ricin-like). The 197 structures belong to **57 accessions** (56 UniProt, 1 UniParc ID without a UniProt FASTA), 19 UniLectin classes, 36 species names (35 fungal). One species (*Ralstonia solanacearum*, bacterial) is mislabelled as fungal origin. Published counts for comparison: 157 fungal structures of 1,740 (9%) as of 2018-08-29 (PMID 30239928); 194 fungal structures, 129 with a ligand (Lebreton 2021, abstract, PMID 34200153). | `analysis/lectin_class_spec/unilectin_fungal_sequences.tsv` (57 rows), `unilectin_fungal_table.py` (method in its docstring) |
| F10 | Clusters of the 56 fungal UniLectin3D sequences (MMseqs2 17-b804f `easy-cluster --min-seq-id 0.3 -c 0.5`): **33 clusters**. Excluding the PA14 class, the yeast Emp L-type class and *Ralstonia* (owned by an active call, house-keeping, bacterial): 43 sequences, **25 clusters**. Of these, 12 sequences in **8 clusters** are from Ascomycota. The most clusters in one species is 3 (*Coprinopsis*). *A. fumigatus* has 1 sequence (FleA, Q4WW81). **M1 coverage of this table (Pfam GA hit to an M1 model):** the five revision-1 models (PF07938, PF21691, PF00337, PF07367, PF21595) hit **16 of 43 sequences, in 7 of 25 clusters**; in Ascomycota 5 of 12 sequences, 3 of 8 clusters (all AAL-like). PF07938 alone (the revision 2 M1) hits 5 of 43 sequences (*Aleuria*, *Arthrobotrys*, *A. fumigatus*, *A. oryzae*, *Scedosporium*), in **3 of 25 clusters**. Cluster-level recall of any M1 choice on this table is therefore 0.28 (five models) or 0.12 (PF07938), and every hit is a structure-derived sequence (leakage `in_reference` until checked). | `analysis/lectin_class_spec/unilectin_fungal_sequences.tsv`, `unilectin_fungal_sequences.faa`, `unilectin_fungal_clusters.tsv`; counts by `cluster_check.py` (reproduces the partition of the `cluster30` column) |
| F11 | Pfam GA hits on the 56 fungal UniLectin3D sequences (`hmmscan --cut_ga`, Pfam 38.2): 46 of 56 have at least one hit. Mapping of UniLectin3D class to Pfam model: AAL-like (FleA, AOL, AAL, SapL1, *Arthrobotrys* AOL, and one bacterial *Ralstonia* lectin) to Fungal_lectin PF07938 (6 of 6); Lyophyllum ginkbilobin-like to LDL PF21691; fungal fruit body lectin to FB_lectin PF07367 (5 of 5 with a sequence); galectin-like to Gal-bind_lectin PF00337 (4 of 4); Coprinus trefoil lectin-like to CCL2-like PF21595; cyanovirin-like to CVNH PF08881 (4 of 4); LysM-like to LysM PF01476; Ricin-like and the *Sclerotinia* trefoil lectin to RicinB_lectin_2 PF14200 (5 of 5); PA14 yeast adhesin to PA14 and GLEYA (10 of 10); yeast Emp L-type to Lectin_leg-like PF03388. **No Pfam model hits** the CML-like, CBM67-like, monocot-lectin-like and fungal prism classes, and 5 of 6 Laetiporus trefoil sequences (the sixth hits ETX_MTX2, a toxin model). Ten sequences have no hit at all. Tectonin hits PF19193 and the Fve (Ig-like) sequence hits PF09259. Hits to PVL-like propellers are FG-GAP_3, not a lectin model. **No fungal UniLectin3D sequence hits Jacalin, Ricin_B_lectin, Laminin_G_3, Fungal_lectin_2 or GBD_Y3.** | `analysis/lectin_class_spec/unilectin_fungal_sequences.tsv` (`pfam_ga_hits` column) |
| F12 | UniProt `Signal` feature on the 56 sequences (UniProt JSON, 2026-10-10): present in 9 of 10 PA14 sequences, 2 of 2 Emp L-type, 2 of 2 LysM, 1 of 1 CBM67-like; **absent** in all AAL-like sequences (including FleA Q4WW81), all fruit body, galectin, ricin-like, prism, trefoil and cyanovirin sequences. Houser 2013 (full text) states that AFL is "predicted not to contain a signal peptide or transmembrane region" (SignalP 4.0), and that no targeting signal is predicted although FleA was found in the secreted proteome of *A. fumigatus* and other lectins without a signal peptide have been found on the cell surface. I have not run SignalP 6 on these sequences. | UniProt REST `accessions` JSON (`ft_signal`), fetched 2026-10-10; PMID 24340081 (PMC3858362) |
| F13 | FleA (*A. fumigatus*, UniLectin class AAL-like, Q4WW81) is already in `adhesins.tsv` with class `indirect_regulator`, evidence EXP, PMID 27058347. | `data/curated/adhesins/adhesins.tsv` (read only) |
| F14 | **FleA and AAL-like lectins: the evidence and the dispute.** (a) Houser 2013, PLoS One 8(12):e83077, PMID 24340081, [doi:10.1371/journal.pone.0083077](https://doi.org/10.1371/journal.pone.0083077), PMC3858362, full text read. Recombinant AFL (= FleA) binds fucosylated glycans (glycan array: Lewis Y preferred; hemagglutination); structure with a selenofucoside (PDB 4agi); polyclonal antibodies stain *A. fumigatus* conidia (indirect: the peptide block did not inhibit staining); AFL induces IL-8 in BEAS-2B cells and L-fucose inhibits it. The authors write "We believe that AFL may be present on the conidial surface and it could be the protein responsible for conidia attachment to the human lung epithelium" (a belief, not a result) and that it "might be an important virulence factor". (b) Kerr 2016, PLoS Pathog 12(4):e1005555, PMID 27058347, [doi:10.1371/journal.ppat.1005555](https://doi.org/10.1371/journal.ppat.1005555), PMC4825926, abstract: twelve fucose-binding sites; protein and conidia bind lung mucin in a fucose-dependent way; macrophages take up ΔfleA conidia less well; mice infected with ΔfleA conidia had more severe pneumonia; the authors conclude FleA is not a virulence factor. (c) Sakai 2019, Med Mycol 57(1):71-83, PMID 29370403, [doi:10.1093/mmy/myx163](https://doi.org/10.1093/mmy/myx163), abstract: fleA disruption gave earlier germination and enhanced pathogenicity in mice. (d) Richard 2018, Sci Rep 8:15699, PMID 30356167, [doi:10.1038/s41598-018-33902-0](https://doi.org/10.1038/s41598-018-33902-0), abstract: contact with bronchial epithelial cells inhibits conidial germination, and the process involves PI3K on the host side "and the lectin FleA on the fungal side". (e) Duca 2022, ACS Chem Biol 17(12):3515-3526, PMID 36414265, [doi:10.1021/acschembio.2c00708](https://doi.org/10.1021/acschembio.2c00708), abstract: multivalent fucosides "reduced the adhesion of [*A. fumigatus*] spores to lung epithelial cells"; the paper frames FleA and SapL1 as anti-adhesion targets (an inhibitor result, so an indirect test of the lectin). (f) Martinez-Alarcon 2021, Sci Rep 11:16109, PMID 34373510, [doi:10.1038/s41598-021-95008-4](https://doi.org/10.1038/s41598-021-95008-4), abstract: SapL1 of *Scedosporium apiospermum* (UniProt A0A084FYP2) was found as a FleA homolog, binds fucosylated glycans (array, ITC), has a crystal structure with a fucoside, and "binds to bronchial epithelial cells in a fucose-dependent manner". (g) Dussouy 2021, Org Biomol Chem 19(14):3234-3240, PMID 33885578, [doi:10.1039/d1ob00152c](https://doi.org/10.1039/d1ob00152c), abstract: calls the role of FleA "controversial". (h) Balogh 2003, Fungal Genet Biol 39(2):128-35, PMID 12781671, [doi:10.1016/s1087-1845(03)00023-9](https://doi.org/10.1016/s1087-1845(03)00023-9), abstract: deletion of the AAL-like lectin AOL in *Arthrobotrys oligospora* (UniProt G1XA82) gave no significant change in germination, growth or pathogenicity. **What the papers show.** Two knockout studies (b, c) show that FleA is not needed for virulence in mice, and that loss of FleA does not reduce disease. A disease outcome does not show whether FleA mediates attachment. Papers (a), (b), (d), (f) show host-glycan or host-cell binding by FleA or SapL1, and (a), (d) an epithelial response. Paper (e) supports an attachment role by an inhibitor result. Papers (b), (c), (h) show that gene loss gives no loss of function in the readouts used. **The role (adhesin, host-recognition lectin, or neither) is disputed. The spec does not decide it.** | PubMed (all abstracts read; (a) full text read) |
| F15 | PA14 lectin-type adhesins: Veelders 2010, PNAS 107(52):22511-6, PMID 21149680, [doi:10.1073/pnas.1013210108](https://doi.org/10.1073/pnas.1013210108): Flo5A domain structures with high-mannose ligands. Ielasi 2012, Acta Cryst D 68:210-7, PMID 22349222, [doi:10.1107/S0907444911054898](https://doi.org/10.1107/S0907444911054898): Epa1p ligand domain with lactose. Diderrich 2015, J Biol Chem 290(32):19597-613, PMID 26105055, [doi:10.1074/jbc.M115.655654](https://doi.org/10.1074/jbc.M115.655654): 17 Epa paralogs, most recognise glycans with terminal galactosides and confer epithelial cell adhesion. de Groot and Klis 2008, Mol Microbiol 68(3):535-7, PMID 18331471, [doi:10.1111/j.1365-2958.2008.06182.x](https://doi.org/10.1111/j.1365-2958.2008.06182.x): commentary. Kock 2018, Front Microbiol 9:2581, PMID 30425696, [doi:10.3389/fmicb.2018.02581](https://doi.org/10.3389/fmicb.2018.02581): PA14/Flo5-like adhesins of methylotrophic yeasts (*Komagataella*; Flo1/Cea1 binds terminal GlcNAc and promotes adhesion to chitin when expressed in *S. cerevisiae*). Ielasi 2016, mBio 7(4):e00584-16, PMID 27406561, [doi:10.1128/mBio.00584-16](https://doi.org/10.1128/mBio.00584-16): glycan arrays of *C. albicans* Als1p, Als3p and *C. glabrata* Epa1, Epa6, Epa7 "with lectin activity" and predicted host receptors. Als1p and Als3p are not in a lectin model of this spec; I did not check which Pfam family covers their glycan-binding region. | PubMed (abstracts) |
| F16 | Survey: Lebreton 2021, J Fungi 7(6):453, PMID 34200153, [doi:10.3390/jof7060453](https://doi.org/10.3390/jof7060453), PMC8227253. Full text read in review 1: "33,518 putative lectin sequences (28,691 with filter at 0.25 quality score)"; the abstract gives 33,485. Section 3.4 of the paper names calnexin, calreticulin, the P-type lectin and p58/ERGIC-53 as house-keeping (review 1). The predictions are not experimentally validated. Related: Bonnardel 2021, Nucleic Acids Res 49(D1):D1548-D1554, PMID 33174598, [doi:10.1093/nar/gkaa1019](https://doi.org/10.1093/nar/gkaa1019): UniLectin classification (35 folds, 109 classes, 350 families). UniLectin3D: Bonnardel 2019, Nucleic Acids Res 47(D1):D1236-D1244, PMID 30239928, [doi:10.1093/nar/gky832](https://doi.org/10.1093/nar/gky832), PMC6323968. | PubMed; review 1 |
| F17 | Endemic dimorphic fungi (*Coccidioides*, *Histoplasma*, *Blastomyces*). PubMed searches: "Coccidioides lectin" (32 hits), "Histoplasma capsulatum secreted lectin-like protein galectin" (0), "Blastomyces lectin OR Histoplasma lectin fungal carbohydrate-binding protein" (25), "Histoplasma lectin" (49, revision 1). I read titles or abstracts for 8 of the 32 *Coccidioides* hits in revision 1 and 8 titles in revision 2 (I did not count the overlap, and I did not open the rest), and for 5 of 25 and 8 of 49 of the others. All read hits are host-side work (Dectin-1, Dectin-2, mannose receptor, MBL, surfactant protein) or tools that use plant or animal lectins as probes (Grys 2016, J Proteome Res 15(10):3463-3472, PMID 27546806, [doi:10.1021/acs.jproteome.5b01054](https://doi.org/10.1021/acs.jproteome.5b01054)). **I found no paper that shows a fungal-encoded lectin in *Coccidioides*, *Histoplasma* or *Blastomyces*.** This is "not found by these searches". BAD-1 binds heparin through tandem repeats and WxxWxxW motifs, not a lectin fold (Brandhorst 2013, PMID 23853587, [doi:10.1371/journal.ppat.1003464](https://doi.org/10.1371/journal.ppat.1003464); Beaussart 2015, PMID 26396244, [doi:10.1128/mBio.01403-15](https://doi.org/10.1128/mBio.01403-15)). BAD-1 belongs to the repeat call. | PubMed |
| F18 | Mushroom lectin review: Erjavec 2012, Trends Biotechnol 30(5):259-73, PMID 22341093, [doi:10.1016/j.tibtech.2012.01.004](https://doi.org/10.1016/j.tibtech.2012.01.004). Abstract only. Higher fungi (Agaricomycetes). | PubMed |
| F19 | Tool structure. A Pfam family is a row of `data/sorting_hat/family_table.tsv` with a `module` name. A new module name (as `pfam_hydrophobin` and `pfam_hsba`) gives its own call and its own module status. A one-module call gets the module status; a call status needs two or more modules (hydrophobin spec decision D3). A family row can carry a `second_condition` (PA14 has `signal_peptide`). The per-module family digest exists: `_family_digest(families, module)` in `src/cellsurface_sorting_hat/modules/cli.py` hashes the rows of one module only. `docs/CLASSES.md`: `surface_attachment_candidate` reads repeat, wall-family domain, hydrophobin and HsbA. | `docs/CLASSES.md`; `docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md` section 3.2; `modules/cli.py` |
| F20 | **Paracoccin (*Paracoccidioides brasiliensis*, Onygenales, Ajellomycetaceae).** Alegre-Maller 2014, PLoS Negl Trop Dis 8(12):e3317, PMID 25474158, [doi:10.1371/journal.pntd.0003317](https://doi.org/10.1371/journal.pntd.0003317): paracoccin (PCN) is a GlcNAc-binding lectin; recombinant PCN activates TLR2 and TLR4; the interaction with TLR2 "depended on carbohydrate recognition" (it was affected by mutation of TLR2 N-glycosylation sites); rPCN protected mice. Oliveira 2017, PLoS One 12(8):e0184010, PMID 28846733, [doi:10.1371/journal.pone.0184010](https://doi.org/10.1371/journal.pone.0184010): PCN is a GlcNAc-binding protein with N-acetyl-beta-D-glucosaminidase activity; found in hyphal tips and the cell wall. Pitangui 2021, Front Mol Biosci 8:700797, PMID 34532342, [doi:10.3389/fmolb.2021.700797](https://doi.org/10.3389/fmolb.2021.700797) (PMC8438136, full text read): PCN has one chitinase region and "four potential lectin sites" found by sequence comparison; the gene is PADG_03347 (Gene ID 22582669); signal peptide cleavage after residue 22, GH18 domain at 29-218; endogenous PCN overexpression gives more severe disease through the chitinase site, and soluble PCN gives protection through the carbohydrate-binding sites. Pitangui 2023, Virulence 14(1):2150455, PMID 36599817, [doi:10.1080/21505594.2022.2150455](https://doi.org/10.1080/21505594.2022.2150455): PCN level changes virulence in an insect larvae model. **Pfam check (this revision):** UniProt C1G842 (Pb18, "Chitinase", 218 aa, signal peptide 1..22 by SignalP in UniProt, gene PADG_03347) matches the paper's length, signal site and GH18 range. `hmmscan --cut_ga` against Pfam 38.2 gives one hit: Glyco_hydro_18 PF00704 (domain score 35.8, positions 31-176). At E < 1e-3 the same single hit. **No lectin model hits paracoccin. None of the 18 models of the scan hits it.** The lectin sites were predicted from sequence by the authors; I did not check whether GlcNAc binding is a property of the GH18 fold. | PubMed; UniProt C1G842; `hmmscan` on 2026-10-10 (log not committed; command in section 12) |
| F21 | **Host-interaction evidence for other fungal carbohydrate-binding families.** (a) LysM: de Jonge 2010, Science 329:953-5, PMID 20724636, [doi:10.1126/science.1190859](https://doi.org/10.1126/science.1190859): *Cladosporium fulvum* Ecp6 "sequesters chitin oligosaccharides ... to prevent elicitation of host immunity". Bolton 2008, Mol Microbiol 69(1):119-36, PMID 18452583, [doi:10.1111/j.1365-2958.2008.06270.x](https://doi.org/10.1111/j.1365-2958.2008.06270.x): Ecp6 is a virulence factor (RNAi in *C. fulvum*; expression in *F. oxysporum* raised virulence); LysM domains are "recognized as carbohydrate-binding modules"; chitin binding is by homology modelling in this paper. (b) CVNH-LysM: Koharudin 2011, Structure 19(5):662-74, PMID 21565701, [doi:10.1016/j.str.2011.03.004](https://doi.org/10.1016/j.str.2011.03.004): *M. oryzae* MGG_03307 (UniProt G4N991) has an NMR structure; both domains bind carbohydrate (NMR titrations); live-cell imaging shows it in the appressorium; "Our results suggest that MGG_03307 plays a role in the early stages of plant infection". Koharudin 2015, Structure 23(11):2143-54, PMID 26455798, [doi:10.1016/j.str.2015.07.023](https://doi.org/10.1016/j.str.2015.07.023): crystal structures with GlcNAc oligomers; the abstract states that the function "is associated with protecting M. oryzae from the host immune response in plants". I read only the abstracts. (c) Fruit body lectin domain: Peng 2021, Fungal Biol 125(11):914-922, PMID 34649678, [doi:10.1016/j.funbio.2021.06.005](https://doi.org/10.1016/j.funbio.2021.06.005): *Beauveria bassiana* BbLec1 (a protein with a "Fruit Body_domain") binds chitobiose and chitin; the knockout has an impaired cell wall and weaker virulence in insects. The role is cell-wall and eisosome function; a host target is not shown. (d) Host side, out of scope: Han 2019, Plant Physiol 179:1416-1430, PMID 30696749 (a *M. oryzae* chitinase binds a rice jacalin-related lectin); Wawra 2019, New Phytol 222:1493-1506, PMID 30688363 (*Serendipita indica*, Basidiomycota, WSC beta-glucan-binding lectins; not Ascomycota and not in this spec). | PubMed (abstracts) |
| F22 | Emp46p and Emp47p (*S. cerevisiae*, Q12396 and P43555, Lectin_leg-like PF03388). Sato and Nakano 2002, Mol Biol Cell 13(7):2518-32, PMID 12134087, [doi:10.1091/mbc.e02-01-0027](https://doi.org/10.1091/mbc.e02-01-0027), abstract read: Emp46p is a homolog of Emp47p with a single C-terminal transmembrane domain and a carbohydrate recognition domain signature; localised to Golgi membranes, cycles to the ER; the double disruption gives "a marked defect in the secretion of a subset of glycoproteins"; "required for the export of specific glycoprotein cargo from the endoplasmic reticulum". Both have a UniProt signal peptide and a TM helix (review 1), so they would pass an R0 gate. | PubMed |
| F23 | **Scan with the models (login node, `hmmsearch --cut_ga`, Pfam 38.2, whole proteomes, no signal peptide condition).** The 137 N1 shared negatives: **0 hit any of the five revision-1 M1 models or any M2 model**; 1 hit an M3 model (CNE1, Calreticulin). Six Phase C species (proteomes: *S. cerevisiae* S288C, *C. albicans* SC5314, *A. fumigatus* Af293 (UniProt), *C. neoformans* H99 from `_workdir/sorting_hat/*.faa`; *A. nidulans* FGSC A4 UP000000560 and *Mycosarcoma maydis* 521 UP000000561 downloaded from UniProt on 2026-10-10 because I found no local copy; the local copies of the first four are the files the tool was run on). Proteins with an M1-set hit: *A. fumigatus* 3 (FleA Q4WW81 and Q4WMH4 hit PF07938; Q4WAA0 hits FB_lectin PF07367); all other five species 0. Per-model counts (proteins): PF07938 2, PF07367 1 (*A. fumigatus*); Jacalin 3, 2, 1, 0, 0, 1 (Afum, Anid, Calb, Cneo, Mmay, Scer); Ricin_B 1 (Afum), 5 (Cneo); LysM 7 (Afum), 17 (Anid), 2 (Cneo), 2 (Mmay); CVNH 3 (Afum), 2 (Anid); Lectin_leg-like 1 to 3 in every species (Emp46, Emp47 and a third protein in *S. cerevisiae*); PA14 and GLEYA as in F1. Paracoccin hits none (F20). | `analysis/lectin_class_spec/m1_scan.sh`, `m1_summary.py`, `n1_fasta.py`, `m1_scan_hits.tsv`, `m1_scan_hits_per_model.tsv` |

## 3. Design

### 3.1 Which models define the class (task item 1)

**Rule (one rule for every model).** Two inputs are stated for every model. Both come from the facts table.

- **H, host-interaction evidence from a paper** for a fungal member of the model's family. H1: the paper shows binding to a host glycan or host cell, or an immune effect that depends on glycan binding. H2: the paper shows a virulence or infection-stage phenotype or a localisation, but not the host target. H0: no such paper found by my searches. H0 is not proof of absence.
- **U, ubiquity**: the share of the 653 genera with at least one pool S protein of the family (F5). I call U "rare" below 10% and "common" at 10% or above. **The 10% threshold is my choice. No data supports it.** The owner decides it (D2).

A model joins the call (set **M1**) if H is H1, U is rare, and the model is not owned by an active narrow-tier family and not a house-keeping class. All other models go to a reported-only set **M2** or an exclusion set **M3**. An M2 model with H1 or H2 and common U is a candidate for a later architecture-defined sub-call (below). The cutoff for every model is the Pfam GA (F7). No custom cutoff is used in v1.

| Set | Model | UniLectin class (F11) | H (papers) | U (genera of 653) | Pool S proteins (F5) |
|---|---|---|---|---|---|
| **M1** | Fungal_lectin PF07938 | AAL-like (PropLec6A) | **H1**: FleA and SapL1 bind fucosylated glycans, FleA binds mucin and induces IL-8, SapL1 binds epithelial cells (F14 a, b, f). Role disputed (F14) | rare (28, 4%) | 32 |
| M2 | LysM PF01476 | LysM-like | **H1**: Ecp6 sequesters chitin oligosaccharides and is a virulence factor (F21 a) | common (511, 78%) | 1,797 |
| M2 | CVNH PF08881 | cyanovirin-like | **H2**: MoCVNH carbohydrate binding in vitro, appressorium location, a stated role in immune protection (F21 b, abstracts) | common (233, 36%) | 560 |
| M2 | FB_lectin PF07367 | fungal fruit body lectin | **H2**: BbLec1 (F21 c) | rare (0 in pools; 1 protein in *A. fumigatus*, F23) | 0 |
| M2 | LDL PF21691 | Lyophyllum ginkbilobin-like | H0. Structure and sugar binding: *Lyophyllum* (PMID 25504801, van Eerde 2014, Glycobiology, per review 1) | rare (31, 5%) | 34 |
| M2 | Gal-bind_lectin PF00337 | galectin-like | H0 (mushroom structures) | rare (2, 0.3%) | 2 |
| M2 | CCL2-like PF21595 | Coprinus trefoil lectin-like | H0 (mushroom structure) | rare (0) | 0 |
| M2 | Fungal_lectin_2 PF18647 | none (F11) | H0; Pfam description "Putative glycan binding domain" | common (170, 26%) | 251 |
| M2 | GBD_Y3 PF22803 | none | H0; Pfam description only | rare (51, 8%) | 70 |
| M2 | Jacalin PF01419 | none among fungal structures | H0 | common (256, 39%) | 308 |
| M2 | Ricin_B_lectin PF00652 | none among fungal structures | H0 | common (218, 33%) | 270 |
| M2 | RicinB_lectin_2 PF14200 | Ricin-like, *Sclerotinia* trefoil | H0 (structures; ligand not checked) | common (122, 19%) | 146 |
| M2 | Laminin_G_3 PF13385 | none | H0; Pfam description "Concanavalin A-like lectin/glucanases superfamily" (mixed function) | rare (54, 8%) | 74 |
| M3 | Lectin_leg-like PF03388 | yeast Emp L-type (Emp46p, Emp47p) | House-keeping: ER-to-Golgi cargo receptors (F22) | n/a | 45 |
| M3 | Calreticulin PF00262 | calnexin, calreticulin | House-keeping ER lectins (Lebreton 2021 section 3.4, F16); CNE1 (J9VLH0) is an N1 negative (PMID 34566931) | n/a | 31 |
| M3 | P-type lectin (MRH domain) PF07915 (PRKCSH, "Glucosidase II beta subunit-like protein") and PF02157 (Man-6-P_recep) | not in the fungal UniLectin3D set | House-keeping by Lebreton 2021 (the P-type lectin, F16). That these two Pfam models are the MRH models of the yeast ERAD lectin is **from UniProt annotation of an OS-9 homolog (PF07915); not checked in a paper** | n/a | not counted |
| M3 | Malectin PF11721 | not in the fungal UniLectin3D set | ER glycoprotein quality control: **from memory, not checked** | n/a | 1 |
| Owned elsewhere | PA14 PF07691, GLEYA PF10528 | PA14 yeast adhesin | Active in `wall_family_domain`, narrow tier (F1, F15) | n/a | 175, 339 |

**Consequence of the rule.** Revision 1 had five models in M1 (PF07938, PF21691, PF00337, PF07367, PF21595). The rule keeps one. Four have no host evidence (H0 or H2) and would join a host-interaction call only because they are rare and have a structure. That is not a host-interaction criterion. They go to M2 (reported). This is a choice that the owner can reverse (D1 option B). The effect on the numbers is in F10 and F23: with PF07938 only, 3 of 25 structure clusters are covered, and 2 proteins are in the six Phase C proteomes (both *A. fumigatus*).

**Ubiquitous classes.** Jacalin, Ricin-B, CVNH, LysM, Fungal_lectin_2 and RicinB_lectin_2 each occur in 19% to 78% of the 653 genera. A call on them would mark many fungi. The spec does not call them. It shows them as a count (`lectin_like_families`) so that a reader can see the proteins. LysM is out of the call because of ubiquity. LysM is also a carbohydrate-binding module of fungal chitinases (Pitangui 2021, F20, lists CBM50 = LysM among them); I did not check the domain content of the 7 LysM proteins of the *A. fumigatus* scan (F23). It is not out because chitin is not a glycan: chitin is a glycan, UniLectin3D lists LysM-like as a lectin class, and the Ecp6 paper shows a host-immunity function (F21).

**Option for a later sub-call (D2).** A discriminator that is not a model score can promote an M2 model. Example for LysM: a protein with a signal peptide, at least one LysM domain, and no catalytic domain of a Pfam model (an "effector-like LysM protein"). This needs its own spec and its own truth. Ecp6 (B3VBK9) and Mg1LysM (F9XHX3) in the UniLectin3D table are the first positives. The same idea could apply to CVNH-LysM (MoCVNH).

**Models without a Pfam entry.** Five UniLectin3D classes of fungal lectins have no Pfam hit (F11): CML-like, CBM67-like, monocot-lectin-like, fungal prism lectin (PhoSL-like) and Laetiporus trefoil lectin-like (5 of 6). UniLectin-derived HMMs would be new custom models. v1 does not build any. All five classes are from mushrooms, and I found no pathogen or host evidence for them in what I read. **Paracoccin is a second kind of gap:** a host-binding protein whose lectin activity no Pfam lectin model sees (F20). The report states both gaps.

### 3.2 Gate (task item 2)

Two forms, following the hydrophobin and repeat precedents:

- `lectin_domain` (evidence call, **ungated**): a hit to an M1 model at GA. It reads the module `pfam_lectin`. It does not need R0.
- `secreted_lectin_candidate` (composite, **gated**): `lectin_domain` and R0 signal peptide. It reads two modules, so it can have a call status.

Reason for no gate on the evidence call: F12. UniProt shows no signal peptide for FleA and for most structurally known fungal lectins, and Houser 2013 reports no signal peptide by SignalP 4.0. If SignalP 6 agrees, an R0 gate (or a `second_condition=signal_peptide` on the family row, as PA14 has) would remove FleA, which is the protein with the most host-interaction evidence. The pool design (signal peptide by construction) never sees it. Task L2 runs R0 on the 56 truth sequences. If R0 calls FleA, the gated form can become the default measured call. If not, the report says that the host-interaction lectin evidence is mostly in proteins without a predicted signal peptide. Houser 2013 and Kerr 2016 report FleA on conidia and in the secreted proteome; how it reaches the surface without a predicted signal peptide is not shown in the texts I read.

**What the ungated call lets in.** A hit to PF07938 in any protein of the proteome, with or without a signal peptide. Cytoplasmic lectins are included. The call does **not** mean that the protein is on the cell surface or secreted. The report text for every `lectin_domain` row says so (section 6.2, item 3).

No GPI or TM condition in v1. Reuse the existing tool rules (R0 for the gate; TMHMM and GPI shown as evidence fields).

### 3.3 Relation to existing calls and double counting (task item 5)

| Existing call | Relation |
|---|---|
| `wall_family_domain` (PA14, GLEYA, CFEM, Flocculin, ...) | Disjoint by construction. M1 contains no model that is active in `family_table.tsv`. PA14 and GLEYA stay in the narrow tier. A protein with both a PA14 hit and an M1 hit is called by both calls, and the report shows both. No such protein is in the 10 PA14 structures (F11). |
| `surface_attachment_candidate` | Not changed. It does not read `lectin_domain`. See D5: this departs from the hydrophobin and HsbA precedent. |
| `other_surface_no_mechanism`, `other_not_surface` | `lectin_domain` joins the mechanism lists (default, as for hydrophobin D2). A secreted protein with an M1 hit is then not "no mechanism". Effect on the discovery pools: pool U loses the M1 proteins (PF07938 only: 31 of 422,071). The discovery tables stay valid because they used the run with 18 active families. |
| `tandem_repeat_protein` | Independent. A repeat protein with a lectin domain is called by both. |
| Antigen and allergen calls | Independent. |

**D5 and the precedent.** `docs/CLASSES.md` shows that `surface_attachment_candidate` already reads hydrophobin and HsbA. Hydrophobin and HsbA are surface-active proteins, not attachment mechanisms in a strict sense. A rule "only strict attachment mechanisms" would remove them. The rule "any surface-active or surface-binding mechanism that a paper links to host or substrate contact" would include a lectin. Revision 1 chose the opposite of the hydrophobin precedent without saying so. Both options are in D5. My default is No for v1, for one reason: the attachment role of FleA is the disputed point (F14). Adding the lectin to an attachment composite would assert the role that the literature does not settle. This reason is about FleA; it does not apply to other M1 members. The owner decides.

**Narrow tier or host-interaction tier.** The split follows what a paper shows about the protein, not about the domain:

- Narrow (adhesion): proteins for which a paper shows an adhesion function by gain or loss of function, for example Epa and Flo (F15: Epa paralogs "confer epithelial cell adhesion"; Cea1 promotes adhesion to chitin). PA14 and GLEYA stay there.
- Host interaction: proteins for which papers show host-glycan binding or an immune effect, and for which an adhesion function is not shown by gain or loss of function. FleA and SapL1 (F14). **Whether FleA is also an adhesin is disputed (F14). The tier label does not settle that.**
- Unknown role: other models. The report says "lectin domain", not a role.

### 3.4 Calls, modules and what changes in the config (stated, not done)

| Item | Content |
|---|---|
| Family table rows | M1 (PF07938) gets `module=pfam_lectin`, `class` text "host-interaction lectin", `source_pmid` from section 2, `specificity_note` with pool counts, `active=no` until the owner activates (decision D1). M2 models get `module=pfam_lectin_like`. M3 models are not rows. |
| Modules | `pfam_lectin` (M1), `pfam_lectin_like` (M2, report only, read by no call). Same wrapper as `pfam_hydrophobin`. The per-module family digest exists (F19), so adding rows to one module does not change the identity of the other Pfam modules. No dependency. |
| Calls | `lectin_domain` = `{flag: pfam_lectin.hit}`; `secreted_lectin_candidate` = `lectin_domain` and `signal_peptide_protein[R0]`. |
| Category | host interaction and immune evasion (tier from the owner decision of 2026-10-09). The current `categories.yaml` has five categories; workstream A decides the category name. This spec does not pick a name. |
| Effect | `config_sha256` changes, so existing call status files are stale (hydrophobin spec 3.5). Batch with workstream A so the repeat call is re-measured once. |

## 4. Truth set (task item 3)

### 4.1 Circularity and tiers

The class is defined by Pfam models. Three problems follow:

1. **A positive that is "found by the model" is not evidence.** The rule of the task is kept: a protein is not truth because a Pfam model hits it.
2. **UniLectin3D classes and Pfam models share sequences.** Pfam seeds may contain the same structure sequences (not checked; there is no `Pfam-A.seed` in the Pfam 38.2 directory). PF07938 has 3 seed sequences (F7). Leakage for all structure-derived rows is `in_reference` until checked.
3. **A recall measured on structure-defined proteins measures how well the model covers the structure it was built from.** It is not recall on unknown proteins. The report says this beside every number.

Tiers (follow the hydrophobin tiers):

| Tier | Evidence | Use |
|---|---|---|
| T1 | A paper reports that the protein **binds a glycan** (structure with ligand, ITC, glycan array or inhibition by a sugar) **and** the sequence is a fungal protein with a UniProt accession. The PMID is recorded. | primary positives |
| T2 | UniLectin3D curated lectin structure without a recorded ligand, or a Swiss-Prot reviewed entry with experimental sugar-binding annotation (ECO:0000269) | primary positives |
| T3 | Family or name by sequence similarity only | secondary, reported apart, never primary |
| T4 | Review judgement from sequence, domain and architecture evidence without a paper | reported only, never truth, never used for sensitivity |
| H flag | A paper reports host-glycan or host-cell binding or a host-relevant phenotype (H1 or H2 of section 3.1). Example: FleA (PMIDs 24340081, 27058347, 29370403), SapL1 (34373510). The flag is an attribute of a T1 or T2 row. It does not change the tier. | host-interaction column in the report |

### 4.2 Sources

1. UniLectin3D (CC BY 4.0, F8). The enumeration used 57 accessions (F9). Licence allows reuse with attribution; the output table cites the three UniLectin papers. A full export would need the owner to ask the UniLectin team or to confirm that the API use is acceptable.
2. UniProt reviewed entries with sugar-binding evidence (task L1 counts them; I did not).
3. PubMed: papers with glycan-binding or host-binding experiments for fungal proteins in the species of the status entries. Found so far: FleA and SapL1 (F14), paracoccin (F20; not in a Phase C species), Ecp6 and MoCVNH (F21; not in Phase C species). None for *Coccidioides*, *Histoplasma* or *Blastomyces* (F17).
4. `data/curated/adhesins/adhesins.tsv` (read only): FleA is there as `indirect_regulator` (F13). The lectin truth may reuse the row's PMID. It does not change the row.

### 4.3 Negatives

| Negative | Source | Status |
|---|---|---|
| ER and Golgi house-keeping lectin-fold proteins | *C. neoformans* CNE1 (J9VLH0), in the shared set as N1 (`secretory_non_surface`, PMID 34566931). Yeast Emp46p and Emp47p (Q12396, P43555): **PMID 12134087 (F22)** shows Golgi localisation and a role in ER export of glycoprotein cargo. Calnexin and calreticulin orthologs in the six Phase C species: the scan finds one Calreticulin-family protein per species (F23); sources to be found in L3. | CNE1: N1 in the shared set. Emp46/47: candidate N1 rows (stratum `secretory_non_surface`) with a source now. They are **not added** by this spec (the shared set is not edited). |
| Enzymes with a lectin-like fold | The four PA14 *A. fumigatus* glucosidases (F3) are N2 and hit PA14, which is not in M1. They do not test specificity of `lectin_domain`. Other enzymes with a carbohydrate-binding module of the M2 models are not yet identified (L3). I cite none from memory. | none usable now |
| Per-species population | As in the hydrophobin spec 4.3: reviewed or curated secreted proteins of the species minus any lectin label, labelled **assumed** | needs owner confirmation (D4) |
| Shared hard negatives | `data/controls/hard-negatives-shared/`: 137 N1 rows in 76 clusters, none above 16 per species (F3). | reuse, N1 only toward the floor |

**No current negative tests the specificity of M1.** CNE1 (calnexin), Emp46/47 (Lectin_leg-like) and the glucosidases (PA14) hit models that are not in M1. I scanned the 137 N1 negatives with the five revision-1 M1 models and all M2 models: **0 hit** (F23). That result is expected from how the negatives were chosen. They were not selected as look-alikes of PF07938. A specificity near 1.0 on them would carry no information about M1. The real look-alikes in the six Phase C proteomes are the unlabelled M1-set hits: for PF07938, *A. fumigatus* Q4WMH4 (weak: sequence score 30.2 against GA 25.3); for FB_lectin (if the owner restores option B), Q4WAA0. I did not search the literature for these two proteins. They are candidates for T4 review (section 5), not negatives.

### 4.4 Cluster counts and the honest answer (task item 3)

| Level | Positives | Clusters at 30% | Against the 20-cluster floor |
|---|---|---|---|
| All 56 fungal UniLectin3D sequences | 56 sequences | 33 (F10) | above 20, but this is a sequence table, not a species entry, and it includes bacterial, ER and PA14 rows |
| After removing PA14, Emp L-type and *Ralstonia* | 43 | 25 | above 20, mostly Agaricomycetes (mushroom) |
| Of the 43, covered by the five revision-1 M1 models | 16 | **7** | far below 20 |
| Of the 43, covered by PF07938 (revision 2 M1) | 5 | **3** | far below 20 |
| Ascomycota only (all 43-set rows) | 12 | 8 | below 20 |
| Ascomycota, covered by PF07938 | 5 | 3 | far below 20 |
| One species, best case | 1 to 3 | 1 to 3 (*Coprinopsis* 3; *A. fumigatus* 1) | far below 20 |
| Onygenales (*Coccidioides*, *Histoplasma*, *Blastomyces*) | 0 | 0 | none in a Pfam lectin model. Paracoccin (*Paracoccidioides*) is one published host-binding lectin that no M1 model detects (F20) |

Statements:

1. **The realistic status is `smoke` at most.** `estimated` needs 20 positive and 20 negative clusters in one species entry and an interval half-width of at most 0.10 (`COMMON-RULES.md` rule 7). No species has more than 3 positive clusters. Pooled results across species are report tables, not status entries (hydrophobin decision D1).
2. **What `smoke` means here.** The code gives `smoke` to a measure with at least one positive and a sensitivity (`status_from_measure` in `src/cellsurface_sorting_hat/status.py`; a measure with a positive but no specificity is at most `smoke`). Leakage other than `none` caps a status at `smoke` (`call_status.py`). So *A. fumigatus* would get `smoke` from FleA, 1 of 1. FleA is probably in the model's reference set: PF07938 has 3 seed sequences and its description names AAL (seed membership not checked). A sensitivity of 1 of 1 on a probable in-reference protein carries no information about sensitivity. **Option (D4b): write no status file for this call until a positive exists that is not structure-derived.** Without a status file the call is `unvalidated`.
3. **Most entries will be `unvalidated`.** Only *A. fumigatus* (FleA, 1 cluster) has a positive in a Phase C species that I can name. Specificity-only entries (no positives) get `unvalidated` (sowgp spec section 4).
4. **N1 negatives per species are also short.** The largest N1 count in one species is 16 clusters (*S. cerevisiae*), below 20. They do not test M1 (section 4.3).
5. **M1 recall on the pooled structure table.** The five revision-1 models cover 7 of 25 clusters (0.28) and PF07938 covers 3 of 25 (0.12). The other structure-derived fungal lectin clusters are not detected by the M1 choice. Part of the gap is by the rule of section 3.1 (models without host evidence are in M2). Part is real: five UniLectin classes have no Pfam model at all (F11).
6. The pooled 25-cluster table has a use: it measures how many fungal lectin classes the M1 models cover. That is a coverage statement, not a sensitivity for a species.

## 5. Curation workflow (task item 4)

Pattern: `analysis/hydrophobin_truth/curation_rubric.md` (written before any decision; evidence lines; decision rules in order; outcomes `lectin (T4)`, `not lectin (T4)`, `unresolved`; T4 never in truth; the owner can overrule).

### 5.1 Evidence lines (written before any decision)

| Line | Meaning | Independent of the Pfam call? |
|---|---|---|
| A | A 3D structure with a glycan ligand exists for the protein or for an ortholog at 90% identity or more over 90% of length (UniLectin3D or PDB) | yes |
| B | Hit to an M1 model at GA | no |
| C | A paper reports **glycan binding** for this gene product (PMID, quote) | yes |
| C-H | A paper reports host binding or a host-relevant phenotype without glycan binding. It sets the H flag. It does not make a lectin call. | yes |
| D | Reciprocal best hit at 40% identity and 70% coverage of both to a T1 or T2 protein | partly |
| E | Other function: a non-lectin domain covers most of the sequence, or the protein is annotated ER or Golgi house-keeping with a source | yes |
| F | Signal peptide by R0 | yes (not a lectin line; used for the gated call and the report) |

### 5.2 Decision rules (order fixed)

1. `lectin (T4)` if C holds, or if A holds with B or D.
2. `not lectin (T4)` if E holds and C does not hold.
3. `unresolved` otherwise (only B, only D, C-H only, or conflicting evidence).

### 5.3 What an agent can do and what needs the owner

| Task | Agent | Owner |
|---|---|---|
| Enumerate UniLectin3D fungal entries, map to UniProt, cluster, run Pfam and R0 (done in part, F9 to F11) | yes | none |
| Fetch UniProt evidence codes, find reviewed sugar-binding entries, map to the Phase C proteomes by sequence (at least 95% identity over 90% of length) | yes | none |
| PubMed search per species for glycan-binding and host-binding papers; quote the sentence | yes | none |
| Apply the rubric to unlabelled hits in the species proteomes (T4) | yes | may overrule |
| Choose M1, M2 and M3 membership, and the ubiquity threshold | proposal only | decides (D1, D2) |
| Accept "assumed" negatives (D4) | no | decides |
| Decide tier label of each role (adhesion, host interaction, unknown) | proposal only | decides |
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
2. A separate line "lectin-like families (not called)": count per M2 family, with the H and U values of section 3.1. The reader sees PF18647, Jacalin, LysM and others without a call.
3. **Fixed text beside every `lectin_domain` row:** "This call is a hit to a lectin domain model. It does not show that the protein is on the cell surface or secreted. It does not show that the protein binds a host glycan." For each call the status line and the leakage statement, and in words: "the truth is structure-derived; recall measures model coverage of known structures".
4. A sequence-level table (pooled, no status): the 56 UniLectin3D sequences, their classes, M1 hits and R0 result; the five classes with no Pfam model (section 3.1). With the coverage numbers: five-model 16 of 43 sequences and 7 of 25 clusters; PF07938 5 of 43 and 3 of 25.
5. T4 outcomes, in separate columns, never in sensitivity.
6. For *Coccidioides*, *Histoplasma* and *Blastomyces*: "no fungal lectin literature found by the searches in F17; calls here are domain evidence only". For *Paracoccidioides*: "paracoccin is a published host-binding lectin (F20); no lectin model hits it".

### 6.3 Discovery-pass families

| Family | Handling |
|---|---|
| PF07938 | In M1. The enrichment (6 of 27 Onygenales genera, q 0.0163, F6) is shown beside the call as an association with the clade. **A hit is family membership, not orthology to FleA.** The pools contain a signal peptide by construction, and FleA has none by UniProt and by SignalP 4.0 (F12). PF07938 is in 0 of 7 Aspergillaceae genera of the pools (F5). So the Onygenales PF07938 proteins in pool U may be a different subgroup from FleA and SapL1. No function is claimed for them. |
| PF21691 | In M2 (no host evidence; rare). The enrichment (6 of 27, q 0.0135) is shown as an association. The structure is from a mushroom. |
| PF18647 | In M2: 15 of 27 genera (q 0.0261) but no fungal structure and no function evidence (F11); common (26% of genera). Promotion needs a structure or a paper (decision D2). |

## 7. Decisions for the owner (stop points)

D1 and D2 must be decided before L4 (family table rows). Defaults are my recommendation. None is data-driven about performance. The defaults of D1, D2 and D7 follow the rule of section 3.1, and the rule has not been tested on any data.

| # | Question | Default |
|---|---|---|
| D1 | Which models are in the call (M1)? Option A: PF07938 only (rule of 3.1: host evidence H1 and rare). Option B: the five revision-1 models (PF07938, PF21691, PF00337, PF07367, PF21595; structure evidence, rare, mostly no host evidence). Rows inactive until a review of the first run. | A. Reason: the call sits in a host-interaction tier, and only PF07938 has H1 evidence. Cost: A covers 3 of 25 structure clusters; B covers 7 of 25. |
| D2 | M2 handling and the ubiquity threshold (10% of genera, my choice, no data). Are M2 models reported only, or can an architecture-defined sub-call (for example effector-like LysM: signal peptide, LysM, no catalytic domain) be specified later? | M2 reported only. Sub-call for LysM and CVNH-LysM: a separate spec, after the owner agrees. |
| D3 | M3 (house-keeping: PF03388, PF00262, PF07915, PF02157, PF11721) excluded from the call and reported as a note? | Excluded. Emp46/47 now have a paper (F22). |
| D4 | Accept per-species "assumed" negatives for this class (as D8 of the hydrophobin spec)? | Yes, labelled "assumed". |
| D4b | Write a status file for `lectin_domain` when the only positive is a structure-derived in-reference protein (*A. fumigatus* FleA, 1 of 1, `smoke`)? | No. Keep the call `unvalidated` until a non-structure-derived positive exists. |
| D5 | Does a lectin domain count as a mechanism for `surface_attachment_candidate`? The hydrophobin and HsbA precedent says surface-active proteins count (`docs/CLASSES.md`). | No for v1, because the attachment role of FleA is disputed (F14). The owner may choose Yes for consistency with the precedent. Only for the `other_*` mechanism lists, the answer is Yes. |
| D6 | Gate: ungated evidence call plus gated composite (section 3.2)? | Yes. Revisit after L2 shows R0 on FleA. |
| D7 | Role label. Option A: "host-glycan-binding lectin, role disputed" for PF07938 proteins (based on F14: binding and immune effects shown; adhesion disputed); no label of adhesin or non-adhesin. Option B: label FleA an adhesin. Option C: label FleA not an adhesin. PA14 and GLEYA stay narrow (gain-or-loss-of-function adhesion papers, F15). | A. The papers do not support B or C. B rests on inhibitor and hypothesis statements (Duca 2022, Houser 2013). C rests on disease outcome, which does not test attachment. |
| D8 | Use the UniLectin3D API for a one-time enumeration (CC BY 4.0, attribution in the table), or ask for a bulk export? | API enumeration, with attribution. No bulk export request without the owner. |
| D9 | Custom HMMs for UniLectin classes with no Pfam model (five classes, 3.1), and for lectins that Pfam does not see (paracoccin)? | Not in v1. |
| D10 | Is a lectin domain shown in the report as an antigen or allergen angle? | No. |

## 8. Owner questions, one at a time

1. **Which models are in the call?** (D1, D2.) Recommended answer: PF07938 only. Reason: the rule of section 3.1 asks for a paper that shows host-glycan binding (H1) and a rare model. Only PF07938 passes. Four mushroom-structure models have no host evidence and become reported-only. Host evidence exists for one revision-1 M1 class (AAL-like) and for two M2 classes (LysM, with chitin sequestration and virulence; CVNH, weaker). The cost is stated: the call covers 3 of 25 structure clusters. If the owner prefers coverage, option B covers 7 of 25.
2. **Gate.** (D6.) Recommended: ungated evidence call plus a gated composite, and decide after L2 whether FleA has an R0 signal peptide. The report says that the ungated call does not show surface location.
3. **Do you accept that the status will be `smoke` at most and mostly `unvalidated`, and that no status file is written for the one in-reference positive?** (Section 4.4, D4b.) Recommended: yes. M1 recall on the structure table is 7 of 25 clusters (five models) or 3 of 25 (PF07938), so any recall number says little about unknown lectins. I do not claim that more experiments are needed in the endemic fungi: the evidence I found shows that published host-binding lectins exist outside the six Phase C species (SapL1, paracoccin, Ecp6), and a literature-based truth set from them has not been built.
4. **Role of FleA.** (D7.) Recommended: label PF07938 proteins "host-glycan-binding lectin, role disputed". Do not decide adhesin or non-adhesin. `adhesins.tsv` keeps its `indirect_regulator` row; the lectin truth uses its PMID without editing it.
5. **Attachment composite.** (D5.) Recommended: No for v1. The hydrophobin and HsbA precedent points the other way. State your choice.
6. **Next literature step.** Recommended: L3b as listed in section 11 (it now includes the new papers and has a stop rule). Not more searching before the owner answers questions 1 to 5.

## 9. Risks

1. **Truth is structure-derived and mostly from mushrooms.** A recall number says little about pathogen lectins. Stated in the report.
2. **Circularity with Pfam.** Leakage `in_reference`; status capped at `smoke`.
3. **Gate.** If SignalP 6 agrees with UniProt and SignalP 4.0 (F12), the host-interaction lectins that matter are not in the secretome pools of the discovery pass. The discovery enrichments for PF07938 describe the small subset with a signal peptide.
4. **Pfam release mismatch.** Pools use Pfam-A 37.x, GA values 38.2 (F6).
5. **UniLectin3D enumeration is a lower bound** (197 of 202, F9). One accession has no sequence. Counts can change.
6. **The mislabel** of *Ralstonia* (F9) shows that the origin field needs a check per row.
7. **Pool U is 96.8% of pool S** (F5): "no mechanism" is a weak filter, so pool U is not a test of whether a lectin family is unexplained.
8. **Role over-reading.** A lectin domain does not show that the protein binds a host glycan or triggers inflammation. Lebreton 2021 states that its predictions are not experimentally validated (F16).
9. **A Pfam-defined class misses lectins that Pfam does not model.** Paracoccin is one (F20). The report lists this gap; it does not close it.
10. **The M1/M2 rule is untested.** The ubiquity threshold (10%) and the H levels are my choices. A different threshold changes M1: at 5% PF07938 and LDL stay rare, at 1% PF07938 becomes common (4% of genera).

## 10. What I could not verify, and assumptions

- **Not read in full text:** Kerr 2016, Sakai 2019, Dussouy 2021, Richard 2018, Duca 2022, Martinez-Alarcon 2021 (SapL1), Balogh 2003, Veelders 2010, Ielasi 2012 and 2016, Diderrich 2015, Kock 2018, Brandhorst 2013, Beaussart 2015, Grys 2016, Erjavec 2012, all paracoccin papers except Pitangui 2021, de Jonge 2010, Bolton 2008, Koharudin 2011 and 2015, Peng 2021, Sato and Nakano 2002. Read in full text: Houser 2013 (this revision), Pitangui 2021 (this revision), Lebreton 2021 and UniLectin3D 2019 (earlier). All other statements about papers come from abstracts. A sentence from an abstract may omit a limit that the full text states.
- **Paracoccin identity:** I matched gene PADG_03347 to UniProt C1G842 by gene name, length (218 aa), signal peptide site (22) and GH18 range. I did not align the UniProt sequence to the sequence in a paper. The Pfam result applies to C1G842. Whether GlcNAc binding by paracoccin is by a lectin site or by the GH18 site is not known to me.
- **Pfam seed overlap** with UniLectin3D structures: not checked (no seed file).
- **SignalP 6 was not run** on FleA or any truth sequence. The statements about signal peptides come from UniProt and from Houser 2013 (SignalP 4.0).
- **Houser 2013 conidia staining** is antibody based; the paper notes that the peptide block did not inhibit staining. I take the surface location of FleA as "suggested".
- **The UniLectin licence and API behaviour** (F8) were read in revision 1 and not re-checked by review 1 or by me.
- **M3 P-type lectin models:** that PF07915 and PF02157 cover the MRH domain of the yeast ERAD lectin comes from a UniProt annotation (PF07915 on an OS-9 homolog of *Paracoccidioides*) and from Pfam descriptions. I did not check a paper. Malectin as an ER quality-control protein is **from memory, not checked**.
- **Scan inputs.** *A. nidulans* and *M. maydis* proteomes are UniProt reference proteomes (UP000000560, UP000000561), not necessarily the annotation used for the Phase C truth sets. The scan is whole-proteome hits, not a pool S scan.
- **Search limits.** PubMed queries for the endemic fungi, *Candida*, *Magnaporthe*, *Aspergillus* and *Scedosporium* were few and mostly read at title or abstract level (F17, F21). A negative result is "not found by these searches", not "does not exist". I did not search *Candida albicans* beyond Als, Epa and "lectin-like"; I did not search for other *Aspergillus* lectins beyond FleA (the "Aspergillus fumigatus lectin conidia host glycan" search returned 18 hits of which I read 8). *Blastomyces* and *Histoplasma* searches returned 25 and 49 hits of which I read 5 and 8.
- **Assumption:** that "lectin" in this tool means a Pfam-defined glycan-binding domain. A protein that binds glycans without such a model (BAD-1 and heparin, F17; paracoccin, F20) is not in the class.
- **No SLURM jobs were run.** I ran MMseqs2 and `hmmsearch` or `hmmscan` on the login node on 56 sequences, 137 sequences, 1 sequence and six proteomes (about 45,000 proteins).

## 11. Work plan

| Task | Content | Check |
|---|---|---|
| L0 | Owner answers D1 to D5, D7 and questions 1 to 6 | decisions recorded in this file |
| L1 | Add UniProt evidence codes and reviewed status for the 56 accessions; complete the 5 missing structures (query by class and species in smaller slices). The enumeration table, FASTA and clusters are committed (section 12). | table with ECO codes; count of T1, T2, T3 rows by species and cluster |
| L2 | Run R0 (SignalP 6) and TMHMM on the 56 sequences with the tool's own wrappers (a small run; check the queue first) | FleA R0 result known; D6 revisited |
| L2b | **Done in this revision (login node):** the 137 N1 negatives and the six Phase C proteomes scanned with the five revision-1 M1 models and the M2 and M3 models (F23). Left: scan the N1 set and proteomes with the final M1 after D1. | `m1_scan_hits.tsv`; N1 hits 0; *A. fumigatus* 3 M1-set proteins |
| L3 | Negatives: IDs and sources for calnexin, calreticulin and Emp46/47 orthologs in the Phase C species (Emp46/47 source: PMID 12134087); **look-alike negatives for PF07938**: reviewed proteins with a beta-propeller or other lectin-like fold and a paper that gives another function; look at Q4WMH4 (*A. fumigatus*) and its orthologs; enzymes with CBM13-type or other M2 domains with a paper | N1 rows with PMID; counts of clusters |
| L3b | Targeted literature pass. **Must include** Houser 2013, Kerr 2016, Sakai 2019, Richard 2018, Duca 2022, Martinez-Alarcon 2021 (SapL1), the paracoccin papers (PMIDs 25474158, 28846733, 34532342, 36599817), Ecp6 (20724636, 18452583), MoCVNH (21565701, 26455798), BbLec1 (34649678), Emp46/47 (12134087). Read the full text of Kerr 2016 and Sakai 2019. New searches: other *Aspergillus*, *Candida* (Als lectin region), *Coccidioides*, *Histoplasma*, *Blastomyces* fungal lectins. **Stop rule:** if L3b finds a host-binding lectin family that no M1 model detects (paracoccin already meets this), the owner decides whether to build a custom model (D9) before L4. | list with PMID and quote; `unverified_notes.md` |
| (stop) | Owner reviews `counts.md` and `unverified_notes.md` | |
| L4 | Family table rows (inactive), module `pfam_lectin` and `pfam_lectin_like`, calls; tests first; after workstream A decides category names. No digest dependency (F19). | tests pass; `config_sha256` change noted |
| L5 | Run `hmmsearch --cut_ga` on the Phase C proteomes and the discovery pool proteomes; counts per proteome | module tables |
| L6 | Truth tables per species, status files (only if D4b allows); apply the rubric to unlabelled M1 hits (T4) | status files; report |
| L7 | Docs: `docs/CLASSES.md` entries (generated), paper notes | class-table test passes |

Each task is its own commit. An independent review of this revision precedes L1.

## 12. Sources and files

Files read: `docs/agent-tasks/README.md`, `COMMON-RULES.md`, `09-lectin-class-spec.md`; `docs/CLASSES.md`; `docs/HANDOFF-2026-10-09-overnight.md`; `docs/reports/2026-10-10-secretome-discovery-pass.md`; `analysis/family_tests/README.md`; `analysis/hydrophobin_truth/curation_rubric.md`; `docs/superpowers/specs/2026-10-09-sowgp-like-calls-design.md`; `docs/superpowers/specs/2026-10-08-hydrophobin-validation-design.md`; `data/sorting_hat/family_table.tsv`; `data/controls/hard-negatives-shared/controls.tsv`; `data/curated/adhesins/adhesins.tsv` (read only); `src/cellsurface_sorting_hat/status.py` and `call_status.py`; review 1. Local, read only, in the main checkout: `_workdir/discovery/pools/`, `_workdir/discovery/pfam_names.tsv`, `_workdir/family_tests/`, `_workdir/sorting_hat/*.faa`.

Committed in `analysis/lectin_class_spec/` (revision 2):

| File | Content |
|---|---|
| `unilectin_fungal_table.py`, `unilectin_fungal_sequences.tsv` | revision 1: enumeration script and table |
| `unilectin_fungal_sequences.faa` (new) | the 56 UniProt sequences, header = accession (from UniProt, CC BY 4.0, fetched 2026-10-10); made by `make_sequences.py` |
| `unilectin_fungal_clusters.tsv` (new) | MMseqs2 17-b804f `easy-cluster --min-seq-id 0.3 -c 0.5`; `cluster_check.py` shows the same partition as the `cluster30` column and prints the 16 of 43 and 7 of 25 coverage counts |
| `pool_recount.py`, `pool_recount.tsv` (new) | F5 recount (needs the local pools directory) |
| `m1_scan.sh`, `n1_fasta.py`, `m1_summary.py`, `m1_scan_hits.tsv`, `m1_scan_hits_per_model.tsv` (new) | F23 scan of the N1 set, paracoccin and six proteomes. Pfam 38.2, `hmmsearch --cut_ga`, hmmer 3.4 |

Paracoccin command (not scripted): `hmmscan --cut_ga --domtblout out.domtbl /srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm C1G842.faa` (C1G842 from `https://rest.uniprot.org/uniprotkb/C1G842.fasta`); and the same with `-E 1e-3` instead of `--cut_ga`. Intermediate files (proteome FASTAs, domtbl files, MMseqs2 working files) were in a scratch directory and are not committed.
