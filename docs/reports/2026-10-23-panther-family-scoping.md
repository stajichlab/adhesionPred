# PANTHER (via InterPro) as a source of comprehensive, independent family definitions: scoping

*2026-10-23. Author: Claude Code (claude-sonnet-5-5). Scoping only: no HPCC or SLURM step, no model scored, no sequence downloaded, and nothing here says how well any model performs. Tables, scripts and the exact API calls with dates: `docs/reports/data/panther_family_scoping/` (README there). Earlier work this builds on: `docs/reports/2026-10-23-ecm33-independent-set-scoping.md` (PTHR31018 only).*

## 1. Answer in short

- **PANTHER gives a family for most of the wall and fusion proteins and for almost none of the adhesins.** Of the 112 control rows (all have a UniProt accession; the 20 JGI-key rows of `controls_local_sequence.tsv` were skipped, they have none), 91 have a PANTHER family: wall 74 of 76, fusion 5 of 5, adhesion 12 of 31. Of the 57 UniProt accessions in the 44 clusters missed by the Pfam panel, 50 have one. The 21 control rows with none are 19 adhesion rows (class I hydrophobins hum2, bad1, rodA, rodB, eas, and hyd1_like; the flocculin / Flo11 / Hwp1 group flo1, flo5, flo9, flo10, flo11, fig2, rbt1, pga62_flo1, hwp1; aga1, ece1, iff4, dfi1) and two wall rows (kre1, pga13). hyd1_like (A0ACE4KFX7), aga1 (P32323) and ece1 (Q07730) have none, as expected.
- **PANTHER is a more comprehensive family than the Pfam model in one case that matters here: Ecm33.** PTHR31018 has 1,757 fungal proteins in UniProt (1,929 in InterPro), against 555 for PF12454; only 541 (31 percent) of the 1,757 carry PF12454. For most other panel families the PANTHER family and the Pfam family are nearly the same set of proteins (for 16 of the 21 wall/fusion models 92 to 100 percent of the Pfam sample falls in one PANTHER family and the two fungal counts are close; the other 5 split across families, section 4), so PANTHER adds few new members there.
- **Independent held-out set (question 1):** yes, for the families in section 5; the best is still Ecm33, then the chitin synthase export / trafficking families Chs7 (PTHR35329), Chs6 (PTHR31975), Chs5 (PTHR47351), Big1 (PTHR28285) and Kre6 (PTHR31361). **Better member set for building custom HMMs (question 2):** only Ecm33 (3.2 times more fungal proteins than the Pfam family); for the 13 of the 38 clusters whose PANTHER family is specific, the family is the same protein set as the Pfam family that already hits (fungal counts within 3 percent for Chs7, Chs5, Chs6, Kre6, Ste2 and Ste3, 13 percent for Ags; not compared for the others), and for the other 14 of the 27 clusters with a family it is broader than the wall function. **Nothing useful (question 3):** all adhesin groups listed above, hyd1_like, aga1, ece1, iff4, dfi1, kre1, pga13, and 7 of the 38 clusters that are JGI-only.
- Everything is **circular in method**: PANTHER is a profile-HMM family. Agreement between a custom HMM and a PANTHER family measures agreement between two profile models, and the PANTHER family may be wider than the function (section 7).

## 2. How the questions were answered

Per accession (InterPro REST, three calls each): PANTHER family and subfamily with the score InterPro reports, all InterPro entries, Pfam entries (`accession_to_panther.tsv`, 112 accessions: union of the 57 UniProt accessions of `clusters_missed.tsv` and the 112 `controls.tsv` rows). Per Pfam family of the panel (30 candidates of the specificity test, 18 active `family_table.tsv` rows, 21 wall/fusion models: 69 families, none in two sets): the PANTHER families that hold its fungal proteins (UniProt REST, all reviewed fungal proteins plus the first 500 unreviewed ones as a sample; `pfam_to_panther.tsv`). Per PANTHER family (61 families: the 41 named by an accession plus the families holding at least 10 percent of a Pfam sample): fungal and non-fungal counts, reviewed count, a fungal protein stream with UniProt lineage ids for the fine clade (rule of `analysis/wall_fusion_truth/map_kiss.py::clades_from_lineage`), Pfam cross-references (`family_summary.tsv`, `family_clade_counts.tsv`). The UniProt PANTHER cross-references lag InterPro (Ecm33: 1,757 against 1,929), so family counts here are UniProt counts.

PANTHER score field: InterPro gives one value per protein per family match, an E-value of the family HMM. Of the 91 mapped accessions, 58 report 0 (underflow), the largest is 6.4e-4 (opy2, PTHR14596:SF72), and 1 of 91 is above 1e-5. One oddity: in the all-taxa listing of PTHR31018 one protein (Paraphaeosphaeria minitans) has score 580, so the field is not always an E-value and needs a check before any cutoff is used. The PANTHER site itself (`geneinfo`, v19) returns family, subfamily and GO terms and no E-value; its reference proteomes are 144 genomes, of ours only S. cerevisiae, C. albicans, S. pombe, N. crassa, U. maydis, C. neoformans and B. dendrobatidis (no Aspergillus), so InterPro is the route for scoring all of UniProt.

## 3. (a) The 44 clusters the Pfam panel missed

`missed44_cluster_panther.tsv`. Classes are those of the step 6.0 scan (`docs/reports/data/pr6_pfam_scan_2026-10-22/`).

| Class (step 6.0) | Clusters | PANTHER family for at least one UniProt accession | JGI-only (not testable here) | UniProt accession, no PANTHER family |
|---|---|---|---|---|
| hit_other_pfam | 38 | 27 | 7 (KS0112 neg1, KS0146 cps1, KS0147 cot-5, KS0148 to KS0151 ncl) | 4 (iff4, dfi1, kre1, pga13) |
| hit_panel_model | 2 | 2 (pir1 PTHR47254, cfmA PTHR37928) | 0 | 0 |
| no_pfam_hit | 4 | 1 (ecm33 PTHR31018:SF3, E 7.9e-30) | 0 | 3 (hyd1_like, aga1, ece1) |

The 27 clusters with a family, by the breadth of the family (count rule of section 6): 13 specific, 9 intermediate, 5 broad.

- Specific: SMI1 (KS0006, PTHR47432), Kre6 (KS0010, PTHR31361), Chs5 (KS0063, PTHR47351), Chr4 and Chs4 (KS0065 and KS0084, both PTHR46430), Big1 (KS0087, PTHR28285), Chs7 (KS0088, PTHR35329), Chs6 (KS0091, PTHR31975), Ags (KS0115, PTHR47182), Ste2 (KS0030, PTHR28009), Ste3 and pra1 (KS0067 and KS0079, PTHR28097), cfl1 (KS0056, PTHR35192, 1,521 Basidiomycota proteins, named only "PROTEIN, PUTATIVE-RELATED").
- Intermediate: Kex2 (KS0007, PTHR42884, subtilisin-like proprotein convertases, 13,000 Metazoa), Kre5 (KS0009, PTHR11226), the three PMT clusters (KS0016, KS0035, KS0137, PTHR10050, subfamilies SF46, SF50, SF51 follow the PMT1/2/4 classes), Dfg10 (KS0090, PTHR14624), opy2 (KS0126, PTHR14596), so_soft (KS0138, PTHR46689), rep1 (KS0143, PTHR48128: only 1 fungal protein of 561 in all taxa, so it is not a fungal family).
- Broad (do not use as a positive set): glfB (KS0120, PTHR11132 triose-phosphate transporters, 9,969 fungal and 22,000 plant proteins), gfa1 (KS0034, PTHR10937), Cwh41 (KS0114, PTHR10412), orlA and tps2 (KS0119, PTHR10788), irs4 (KS0123, PTHR11216).

## 4. (c) Controls, family table and the 21 models

- Controls by slice (`controls_panther_coverage.tsv`): adhesion 12 of 31, fusion 5 of 5, wall 74 of 76. The 12 adhesion rows with a family sit in the class II hydrophobin family (hyd5, hfb1, hfb2, NCU08192: PTHR42341), CFEM (wap1_csa1, rbt5), ALS (als3_als1: PTHR33793), the multifunctional adhesion family PTHR10068 (eap1, 201 fungal proteins, 1,126 Metazoa), Keratinocyte-proline-rich (bsc1: PTHR48138, 41 fungal proteins), PTHR35192 (cfl1), PTHR48128 (rep1) and the transporter family of glfB.
- Pfam families of the panel and PANTHER (`pfam_to_panther.tsv`). **Wall synthesis and remodeling models** (GH72, GH16, GH17, GH76, GH18, Cellulase, chitin synthases, glucan synthase, Kre9, SUN, X8, Tos1, BtgE, polysaccharide deacetylase): for 16 of the 21 models 92 to 100 percent of the Pfam sample falls in one named PANTHER family (for example PF03198 to PTHR31468 97 percent, PF02364 to PTHR12741 100 percent, PF29664 to PTHR16631 99 percent, PF12454 to PTHR31018 98 percent; PTHR31468 has 6,166 fungal proteins against 6,163 for PF03198); 5 do not (Cellulase 57 percent to PTHR31297, GH18 54 percent, polysaccharide deacetylase 66 percent, Kre9_KNH 37 percent, Fig1 85 percent). Some Pfam families split across several PANTHER families (GH18: PTHR11177 54 percent, PTHR47700 18 percent; PF10342 Kre9_KNH: PTHR28154 37 percent, PTHR40633 25 percent; CFEM PF05730: PTHR37928 37 percent, PTHR33048 22 percent), and some PANTHER families bundle several Pfams (PTHR22914 chitin synthase holds PF01644, PF03142 and PF08407). **Adhesion Pfam families of the active table**: PANTHER matches only ALS (PF05792 78 percent and PF11766 94 percent to PTHR33793), class II hydrophobin (PF06766 83 percent to PTHR42341), CFEM, HsbA (PF12296 84 percent to PTHR38123), PIR (PF22799 97 percent to PTHR47254) and PA14 (PF07691 95 percent to a cellulose-degrading GH3 family, so that one is not a wall-adhesion family). Class I hydrophobin PF01185 (5,422 fungal proteins), Flo11 (PF10182), flocculin repeats (PF00624, PF13928), Eas (PF22354), GLEYA (PF10528), Hyr1 (PF15789), Hyphal_reg_CWP (PF11765) and Hydrophobin_like / D / F (PF29465, PF29785, PF29802) have no PANTHER family holding more than 9 percent of their proteins.
- Interpretation: for the panel, PANTHER is a second, independently built grouping that mostly reproduces what Pfam already groups. It does not widen the adhesion side.

## 5. Families by size, clade spread and specificity (`family_summary.tsv`, `family_clade_counts.tsv`, `test_set_screen.tsv`)

Fine clades are the fine level of the project rule (class for Pezizomycotina and Basidiomycota). "clades>=5" is the number of fine clades with at least 5 fungal proteins. Counts are UniProt (see section 2).

| PANTHER family | Name | Fungal | Reviewed | Clades>=5 | Non-fungal (Metazoa / plants / bacteria) | Fungal Pfam coverage | Biggest clades |
|---|---|---|---|---|---|---|---|
| PTHR31018 | Ecm33 / Pst1 / Sps2 | 1,757 | 13 | 15 | 178 / 47 / 879 | 541 with PF12454 (31 percent) | Sordariomycetes 423, Saccharomycotina 380, Eurotiomycetes 292, Dothideomycetes 264, Mucoromycota 98 |
| PTHR35329 | Chitin synthase export chaperone (Chs7) | 2,957 | 19 | 29 | 0 / 2 / 0 | PF12271 3,005 (same set) | Sordariomycetes 705, Agaricomycetes 464, Eurotiomycetes 457 |
| PTHR31975 | Bud site selection 7 related (Chs6 / ChAPs) | 1,749 | 5 | 24 | 0 / 25 / 1 | PF09295 1,732 (same) | Sordariomycetes 357, Saccharomycotina 299, Eurotiomycetes 240 |
| PTHR47351 | Chs5 (export / cell fusion regulator) | 1,521 | 3 | 21 | 0 / 0 / 0 | PF16893 1,506 (same) | Sordariomycetes 349, Eurotiomycetes 235, Agaricomycetes 213 |
| PTHR28285 | Big1 | 1,349 | 8 | 14 | 0 / 1 / 0 | PF27676 / PF27678 (not in the panel) | Sordariomycetes 356, Eurotiomycetes 237, Saccharomycotina 225 |
| PTHR31361 | Kre6 related | 3,325 | 6 | 22 | 0 / 49 / 0 | PF03935 3,414 (same) | Agaricomycetes 1,522 (46 percent), Dothideomycetes 341, Saccharomycotina 313 |
| PTHR46430 | Chitin synthase regulator (Chs4 / Chr4, Sel1 repeats) | 4,935 | 10 | 31 | 4 / 4 / 246 | Sel1 PF08238 (not in the panel) | Sordariomycetes 1,091, Eurotiomycetes 707, Agaricomycetes 626 |
| PTHR47182 | Ags (cell wall alpha-1,3-glucan synthase) | 1,794 | 5 | 11 | 0 / 20 / 3 | PF26108 etc. 1,580 | Eurotiomycetes 662, Sordariomycetes 337, Agaricomycetes 316 |
| PTHR28097 | Ste3 (a-factor receptor) | 3,188 | 9 | 17 | 0 / 0 / 0 | PF02076 3,272 (same) | Agaricomycetes 1,819 (57 percent), Sordariomycetes 386 |
| PTHR28009 | Ste2 (alpha-factor receptor) | 1,037 | 5 | 9 | 0 / 1 / 0 | PF02116 1,050 (same) | Sordariomycetes 348, Eurotiomycetes 229, Dothideomycetes 179 |
| PTHR31030 | Prm1 (plasma membrane fusion) | 1,544 | 32 | 21 | 0 / 1 / 0 | PF27758 1,553 (same; shipped) | Sordariomycetes 365, Eurotiomycetes 236, Agaricomycetes 219 |
| PTHR31468 | Gas1 (GH72) | 6,166 | 24 | 21 | 1 / 6 / 36 | PF03198 6,163 (same; shipped) | Sordariomycetes 1,686, Eurotiomycetes 1,365, Dothideomycetes 1,131 |
| PTHR42341 | Class II hydrophobin | 1,325 | 45 | 6 | 0 / 0 / 0 | PF06766 1,601 (shipped) | Sordariomycetes 896 (68 percent), Dothideomycetes 268 |
| PTHR33793 | Alpha-agglutinin (ALS) | 645 | 16 | 3 | 0 / 0 / 0 | PF05792 495, PF11766 319 (shipped) | Saccharomycotina 403, Sordariomycetes 223, Dothideomycetes 14 |
| PTHR10050 | PMT (Dolichyl-P-mannose protein mannosyltransferase) | 5,658 | 16 | 32 | 3,439 / 4 / 4,990 | PF02366 5,382 (not shipped) | broad enzyme family (intermediate) |

The family numbers are of the whole family; our own accessions cover few clades of it. For example the 3 control accessions in PTHR31468 are in 2 fine clades (Eurotiomycetes, Saccharomycotina) of a 21-clade family; PTHR35329 has 2 controls in 2 of 29 clades; PTHR31018 has 1 (Q5AGC4, Saccharomycotina) of 15; chitin synthase PTHR22914 has 11 controls in 2 of 35 clades. Only PMT (14 accessions, 5 clades) is broader. The training set of the Ecm33 custom model (481 members, `pr6_ecm33_build_2026-10-23`) has members in 19 of the 21 fine clades; the earlier scoping compared it with the PANTHER set (725 proteins in 187 non-training UniRef50 clusters).

Reviewed (Swiss-Prot) counts are small everywhere except chitin synthase (128), the class II hydrophobins (45), Prm1 (32) and the Gas1 family (24). Most PANTHER-called proteins are unreviewed, so a name-based tier is weak for most families (see the Ecm33 scoping for the T1 to T3 tiers).

## 6. Broad or specific

`broad_or_specific` in `family_summary.tsv` is a count rule and does not read the function: broad if the family has at least 20,000 proteins in all taxa, or at least 5,000 with at least 100 proteins in each of 3 non-fungal kingdoms; intermediate if at least 2 kingdoms hold 100 or more, or the fungal share is below 0.5; specific otherwise. It is a sorting aid. Read the names too: PTHR46430 (Chs4 / Chr4) is called specific by the count but its Sel1-repeat family has 246 bacterial members, so check that the fungal members are the Chs4 / Chr4 clade before using it; PTHR31018 is called intermediate (41 percent non-fungal, 879 bacterial) but its fungal members are one family (restrict to taxon 4751); PTHR12741 (GT48 glucan synthase) is called intermediate because of 8,700 plant callose synthases and is fungal-clean when restricted to fungi. A family that covers a large share of non-wall proteins (PTHR11132, PTHR10788, PTHR10937, PTHR10963 GH16 with 14,000 bacterial members, PTHR11216, PTHR43806) is not a usable positive set. InterPro descriptions of PANTHER families are model-written and unchecked (`is_llm: true`), and PANTHER subfamily names can be wrong (PTHR31018:SF3 is labelled "RECEPTOR PROTEIN-TYROSINE KINASE"); only ids are used here.

Screen (`test_set_screen.tsv`, 5 criteria: at least 500 fungal proteins, 8 fine clades, fungal share 0.90 or more, 5 reviewed, at least 1 of our accessions): 12 families meet all five (Gas1, DCW1 GH76, PTHR46430, CFEM PTHR37928, Ste3, Chs7, Ags, Chs6, Kre9/Knh1 PTHR28154, Big1, SUN, Ste2); 22 more meet four. The screen is only a sorting aid; it does not know whether the Pfam of the family is already in the shipped panel.

## 7. Ranking as independent test sets and staged steps (not run)

The score combines how many proteins the family adds beyond what the shipped Pfam panel already calls, clade spread, specificity and reviewed members. Rank 1 to 5 are the ones I would use; 6 to 10 are possible.

| Rank | Family | Why | Limits |
|---|---|---|---|
| 1 | PTHR31018 Ecm33 | The only family where PANTHER is much wider than the Pfam (1,216 of 1,757 fungal proteins have no PF12454); a custom model exists; 15 clades>=5; 13 reviewed; already staged | Intermediate breadth, 330 or so "Receptor L-domain" proteins not settled as Ecm33 proper; a wider-than-function family is possible |
| 2 | PTHR35329 Chs7 | 2,957 fungal, 29 clades, 19 reviewed, fungal share 0.998, no shipped Pfam | Same size as PF12271, so no new members, only an independent grouping |
| 3 | PTHR31975 Chs6 | 1,749, 24 clades | 5 reviewed |
| 4 | PTHR47351 Chs5 | 1,521, 21 clades | 3 reviewed; fn3 and BRCT domains are generic, check the Pfam split |
| 5 | PTHR28285 Big1 | 1,349, 14 clades, 8 reviewed, no shipped Pfam | one control accession |
| 6 | PTHR31361 Kre6 | 3,325, 22 clades | 46 percent Agaricomycetes; 6 reviewed |
| 7 | PTHR46430 Chs4 / Chr4 | 4,935, 31 clades | Sel1-repeat family, check specificity |
| 8 | PTHR47182 Ags | 1,794, 11 clades | 5 reviewed |
| 9 | PTHR28097 Ste3, PTHR28009 Ste2 | fusion receptors, 3,188 and 1,037 | one clade holds 57 percent (Ste3) |
| 10 | PTHR31468 Gas1 and PTHR42341 class II hydrophobin | large and well supported (6,166 with 24 reviewed; 1,325 with 45 reviewed) | the Pfam model already calls the same proteins, so they test the Pfam, not a gap |

Staged steps for ranks 1 to 5 (none run; all small):
1. Fetch the family's fungal accession list from InterPro (`/protein/UniProt/entry/panther/{PTHR}/taxonomy/uniprot/4751`, paged) with E-value, subfamily, UniProt status and the fine clade; keep ids and metadata only.
2. Map each accession to its UniRef50 cluster and drop the training members' clusters, then cluster the rest with the training members at 30 percent identity and 50 percent coverage on HPCC (the freeze's rule), one representative per cluster as the unit.
3. Fetch sequences at job time into the work directory (never committed).
4. For ranks 2 to 5 first decide with the owner whether a custom HMM is wanted at all, since the Pfam hit exists; then hmmsearch the frozen custom model and the comparators, report recall per fine clade and by tier.
5. Use the PANTHER split of a Pfam family as a hard-negative source: proteins of the same Pfam in a different PANTHER family (for example PF10342 Kre9_KNH across PTHR28154, PTHR40633 and PTHR35185; CFEM across PTHR37928, PTHR33048 and PTHR35607).
6. Report separately a stratum of PANTHER E-value above 1e-5 (weak calls) and clades with no training members.

## 8. Caveats

- PANTHER is a profile-HMM family built on phylogenetic trees; a custom HMM agreeing with it shows agreement of two profile models. A protein called by both may belong to a wider family than the wall function (Receptor L-domain, Sel1, GH16 including bacterial members).
- UniProt cross-references lag InterPro (1,757 against 1,929 fungal proteins for PTHR31018); counts change at each release (PANTHER 19.0 per the PANTHER service; InterPro and UniProt release numbers were not recorded). The UniProt search count and the stream count of one query differ by under 1 percent.
- The Pfam sample of `pfam_to_panther.tsv` is all reviewed plus the first 500 unreviewed proteins UniProt returns, so fractions are indications for large families.
- Licence (not re-checked, to confirm before any redistribution): UniProt data CC BY 4.0, InterPro data CC0 (as stated by EBI), PANTHER terms of the PANTHER site. This work commits ids and metadata only.
- Clades with few or no training members (Entomophthoromycota, Mortierellomycota, Chytridiomycota in Ecm33) can only test a model trained on everything.
- Fine-clade assignment uses the UniProt lineage and the project rule; Basidiomycota classes are as in `map_kiss.py`, so counts of Agaricomycetes can dominate (Kre6, Ste3).

## 9. Not done

No sequence was downloaded and no model or HMM was run; no SLURM or HPCC. JGI-key rows (9 keys in the 44 clusters, 20 rows of `controls_local_sequence.tsv`) were skipped, as they have no accession. Orthologs by name in SGD, CGD, OrthoDB and the literature were not repeated (see the Ecm33 scoping). No cutoffs, thresholds for E-value or subfamily choices were set.
