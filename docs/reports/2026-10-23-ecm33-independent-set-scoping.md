# Independent held-out set for the custom Ecm33 model: scoping

*2026-10-23. Author: Claude Code (claude-sonnet-5-5). Scoping only: no HPCC job was run, no model was scored, and nothing here says how well any model performs. Data: `docs/reports/data/ecm33_independent_set_candidates/` (README there lists the exact API calls). Context: `docs/reports/2026-10-23-pr6-ecm33-build.md`, question 3 (a held-out set that is not the product of a search from Q5AGC4).*

## 1. Answer in short

- Independent evidence of Ecm33-family membership exists and is large. The best source is the **PANTHER family PTHR31018** (InterPro IPR051648), which InterPro reports for **1,929 fungal UniProt proteins** in 21 fine clades. It was assigned by a different method (PANTHER family HMMs built on phylogenetic trees) than our training set (phmmer of Q5AGC4 against UniRef50, E <= 1e-5).
- Overlap with the 481 training members is by UniRef50 cluster: **1,204 of the 1,929 proteins (62 percent) sit in a cluster that is a training member (317 clusters); 725 proteins in 187 other clusters, in 18 fine clades, are not in a training cluster.** Of the 481 training members, 164 (34 percent) are not in any PANTHER-called cluster.
- Everything else is small or not independent: Swiss-Prot has 13 family entries and all 13 are in training clusters; the literature is a handful of papers that name Ecm33 orthologs by homology to ECM33; OrthoDB gives a group list with no UniProt accessions.
- Recommendation: use the 725 proteins of 187 non-training clusters (PANTHER set minus training clusters), report recall per fine clade over one representative per cluster, stratified by evidence tier, after a cluster30 separation from the training members (a staged HPCC step). Details in section 5.

## 2. What exists, by source

| # | Source | What | Count (fungi) | Clades | How membership was decided | Independent of our phmmer search? | Access |
|---|---|---|---|---|---|---|---|
| 1 | UniProt Swiss-Prot (reviewed) | Family entries by gene or name (ECM33, PST1, SPS2, SPS22, Ecm33 names) plus xref PF12454 / IPR024238 | 13 with a PTHR31018 call (12 plus the Meu10 entry O74411; the query also returned two namesakes that are not family: PST1_SCHPO Q09750, a SIN3 homolog, and PST1_CANAL A0A1D8PHR5, a quinone oxidoreductase) | Saccharomycotina 9 (4 S. cerevisiae strains), Taphrinomycotina 2, Eurotiomycetes 2 | Curated by UniProt annotators, in part by homology and sequence-rule annotation | Method yes; but not as a test set: all 13 are in training clusters | UniProt REST, CC BY 4.0 (as stated by UniProt) |
| 2a | InterPro IPR051648 / PANTHER PTHR31018 | Family entry (3,471 proteins in all taxa; "mannoprotein outer layer", Ecm33 / Pst1 / Sps2 / Sps22) | **1,929** (SF3 Ecm33/Pst1 type 1,831; SF12 Sps2/Sps22 type 93; SF11 S. pombe type 5); 13 reviewed, 1,916 TrEMBL | 21 fine clades | Family HMM from a phylogenetic tree of reference proteomes; each protein is scored against the HMM (E-value reported) | Yes (different method and data path). Still sequence similarity | InterPro REST; no login (InterPro data is released as CC0 per EBI, not re-checked here) |
| 2b | Pfam PF12454 (Ecm33) | The model we already use | 557 proteins in InterPro (all taxa); 541 of the 1,747 UniProt PANTHER proteins carry it | n/a | Pfam HMM | Not an independent standard here: it misses P38248 and scores Q5AGC4 at 5.5 against GA 23.2. Pfam PF12454 is not integrated into an InterPro entry | InterPro REST |
| 2c | InterPro IPR036941 (Receptor L-domain superfamily, SSF/CATH) | Structural superfamily seen in the Ecm33 proteins | 17,951 proteins in all taxa | n/a | Structure-based superfamily | Too broad (includes non-fungal receptors); not usable | InterPro REST |
| 3 | OrthoDB v12 | Fungal-level orthologous group `1159006at4751` ("cell wall") contains ECM33, PST1, SPS2 and SPS22 genes of S. cerevisiae S288C, C. albicans SC5314, S. pombe, A. fumigatus Af293 | 4,277 genes in 1,938 genomes (49 genes with no fine clade); smaller groups `331559at4751` (212 genes) and `425402at4751` (11) are fragments of the family | about 27 named clades | Best reciprocal hit clustering over genomes | Yes, but OrthoDB gives genome-specific gene ids, not UniProt accessions, so the overlap with `members.tsv` can only be stated by clade (below) | OrthoDB public API, no login; license not re-checked |
| 4 | Literature | Papers that characterise ecm33 genes (section 4) | about 15 organisms named; accessions found by name for only 3 of 8 spot checks | Saccharomycotina, Taphrinomycotina, Eurotiomycetes, Sordariomycetes, Leotiomycetes | Authors call orthologs by sequence homology to ECM33; function shown by deletion mutants | Method yes. Not a systematic set | PubMed (public) |
| 5 | SGD, CGD, Alliance, eggNOG | Orthologs of ECM33 | not obtained | | | | SGD `fungal_homology` endpoint returned an error; the Alliance orthologs call returned 0 records for the id tried (S000000282, not verified as the ECM33 id); eggNOG offers bulk downloads only (not fetched). Left for later |

Notes on what the PANTHER and UniProt calls are. All 13 Swiss-Prot hits have a PTHR31018 call. S. cerevisiae P38248 (Ecm33), Q12355 (Pst1), P08459 (Sps2), P25380 (Sps22) and Q4WNS8 (A. fumigatus ecm33) are all in the family, so the two true Ecm33s that sit in the shared hard negatives (P38248, Q4WNS8) are family members by this independent source too. PANTHER names the subfamilies of PTHR31018 with a wrong label ("receptor protein-tyrosine kinase"), so the subfamily ids are used only as ids. The UniProt names of 892 unreviewed rows contain Ecm33 / Sps2 / Pst1 or similar, but most of them are set by automatic annotation from the same domains, so they are not an extra line of evidence. The 330 or so "Receptor L-domain domain-containing protein" rows are in the family by PANTHER; whether they are Ecm33 proper is not settled by any source here.

## 3. Counts by clade and overlap (PANTHER set, `candidates.tsv`)

`new_clusters` is the number of distinct UniRef50 clusters not in the training members; `strong` is the number of those with a PANTHER E <= 1e-10 (the 563 proteins in 109 clusters at that level; 725 / 187 without the cut). Training members per clade are from `members.tsv`.

| Clade (fine) | Proteins | In a training cluster | In a new cluster | New clusters | New clusters, strong | Training members |
|---|---|---|---|---|---|---|
| Sordariomycetes | 457 | 204 | 253 | 39 | 23 | 57 |
| Saccharomycotina | 381 | 285 | 96 | 21 | 19 | 107 |
| Eurotiomycetes | 331 | 263 | 68 | 19 | 7 | 55 |
| Dothideomycetes | 294 | 209 | 85 | 40 | 27 | 106 |
| Leotiomycetes | 108 | 29 | 79 | 18 | 6 | 38 |
| Mucoromycota | 99 | 68 | 31 | 10 | 9 | 33 |
| Kickxellomycota | 71 | 65 | 6 | 3 | 1 | 37 |
| Orbiliomycetes | 58 | 19 | 39 | 6 | 4 | 12 |
| Entomophthoromycota | 40 | 0 | 40 | 21 | 2 | 0 |
| Taphrinomycotina | 18 | 14 | 4 | 1 | 1 | 11 |
| Glomeromycota | 16 | 15 | 1 | 1 | 1 | 8 |
| Lecanoromycetes | 14 | 14 | 0 | 0 | 0 | 4 |
| Mortierellomycota | 13 | 0 | 13 | 2 | 2 | 1 |
| Pezizomycetes | 11 | 9 | 2 | 1 | 1 | 5 |
| Chytridiomycota | 5 | 1 | 4 | 4 | 4 | 1 |
| Zoopagomycota | 3 | 2 | 1 | 1 | 1 | 1 |
| Tremellomycetes | 3 | 2 | 1 | 1 | 1 | 0 |
| Xylonomycetes | 2 | 2 | 0 | 0 | 0 | 1 |
| Basidiobolomycota | 2 | 1 | 1 | 1 | 1 | 1 |
| Agaricomycetes | 2 | 1 | 1 | 1 | 1 | 1 |
| Geoglossomycetes | 1 | 1 | 0 | 0 | 0 | 2 |

Reading it: eight clades have at least 5 new clusters at any E (Sordariomycetes, Dothideomycetes, Saccharomycotina, Entomophthoromycota, Eurotiomycetes, Leotiomycetes, Mucoromycota, Orbiliomycetes). With the strong cut, six clades have at least 5 (Dothideomycetes 27, Sordariomycetes 23, Saccharomycotina 19, Mucoromycota 9, Eurotiomycetes 7, Leotiomycetes 6). Entomophthoromycota (21 new clusters, 0 training members) and Mortierellomycota, Chytridiomycota are clades with essentially no training members: only a model trained on all members can be tested there, not a leave-clade-out one. Taphrinomycotina, Glomeromycota, Pezizomycetes and Geoglossomycetes have too few for a rate.

PANTHER E-values: median 3.4e-24 for proteins in training clusters, 4.1e-20 for proteins in new clusters; 76 of the 725 new-cluster proteins have E above 1e-5 (33 of 1,204 in training clusters). The new-cluster set is therefore somewhat weaker on average, as expected for proteins not reached by the phmmer search, but 563 of 725 are at E <= 1e-10.

Overlap in the other direction: 164 of the 481 training members (Dothideomycetes 47, Eurotiomycetes 27, Sordariomycetes 23, Saccharomycotina 22, Leotiomycetes 19, Mucoromycota 15, others 11) are in no cluster that contains a PANTHER PTHR31018 protein. Either PANTHER misses them or they are not family members; this scoping cannot tell. It bears on how much the training members should be trusted, and it is a second, separate use of the PANTHER set (a purity check of the members).

OrthoDB (source 3), by clade, group `1159006at4751`: Sordariomycetes 946, Saccharomycotina 883, Eurotiomycetes 825, Leotiomycetes 422, Dothideomycetes 376, Kickxellomycota 203, Mucoromycota 192, Entomophthoromycota 79, Lecanoromycetes 75, Orbiliomycetes 56, Mortierellomycota 52, Chytridiomycota 29, Pezizomycetes 24, Taphrinomycotina 22, Glomeromycota 11, and a few genes in smaller clades (Lichinomycetes, Tremellomycetes and others). The counts are genes, many strains of one species, and not comparable with unique proteins. The clade pattern agrees with PANTHER, which is the main thing it adds: a second, differently built grouping (best reciprocal hits) gives the same fungal span. Overlap with `members.tsv` by accession is not available (no UniProt accession in the table).

## 4. Literature (PMIDs found by PubMed search, 2026-10-23)

Search: `(Ecm33 OR ECM33 in title/abstract) AND (fungal OR yeast OR Candida OR Aspergillus OR Saccharomyces OR Schizosaccharomyces)`: 53 records. Papers that characterise an ecm33 gene or its product, by organism:

- Saccharomyces cerevisiae: 15583168 (PST1 and ECM33, GPI proteins and cell wall integrity), 12583915, 10734606 (Las21), 29029364 (Ecm33 and TORC1).
- Candida albicans: 16400176, 22665950, 26870022, 26290404, 18546157, 23243062.
- Schizosaccharomyces pombe: 20032302 (`ecm33+` target of Atf1 and Mbx1), 24454826.
- Aspergillus fumigatus: 16804168, 16672465.
- Aspergillus flavus: 29696338 (AFLA_113120).
- Beauveria bassiana and Metarhizium robertsii: 24549768 (Bbecm33, Mrecm33).
- Fusarium oxysporum f. sp. cubense: 35183338 (FocECM33).
- Trichoderma virens: 30840111 (a GPI-anchored cell wall protein).
- Podosphaera xanthii (powdery mildew; Leotiomycetes): 42325904.
- Proteomics and GPI surveys that name Ecm33 proteins: 12845604 (genome-wide GPI proteins), 29223801, 26878023.

Accession look-ups by name in UniProt: Beauveria bassiana J4UGA7 and Trichoderma virens G9MJ35 are PANTHER-called; Fusarium oxysporum f. sp. cubense N1RSC8 is named "Protein ecm33" but is not in the InterPro PTHR31018 list; A. flavus (AFLA_113120), Metarhizium robertsii, Podosphaera xanthii, C. glabrata and C. parapsilosis gave no entry by those queries. Mapping the loci named in each paper to accessions is a manual job and was not done. The papers decide the orthology by BLAST to ECM33, so they add functional evidence (a deletion phenotype in cell wall integrity) for a handful of proteins, not an independent membership call; all organisms are in clades already in training. Curated count: about 15 organisms, 3 mapped to accessions here.

## 5. Independence assessment and recommendation

What is independent of the phmmer search: the PANTHER family calls and OrthoDB groups were made without Q5AGC4 and without UniRef50 as the search space. What is not: Pfam PF12454 (our comparator, and built from related seeds), UniProt names assigned by automatic annotation (they follow the same domains), and anything that is in a training UniRef50 cluster (it is a near-copy, at 50 percent identity or more, of a training sequence). Both PANTHER and phmmer are sequence-similarity methods, so a protein that is found by both is "family by two methods", not "family by function". The Swiss-Prot and literature entries add function for a few proteins and no independent membership beyond that.

Recommended test set:

1. **Primary: PANTHER PTHR31018 proteins in UniRef50 clusters that are not training members** (725 proteins, 187 clusters, 18 fine clades; `candidates.tsv` with `overlap_with_training == new_cluster`). Use one representative accession per cluster (the cluster is the unit, as with the build's "units are clades, not independent clusters" caveat), and a second cut at PANTHER E <= 1e-10 (563 proteins, 109 clusters) as the headline, the full set as the descriptive column.
2. Report by tier: T1 Swiss-Prot (none of the 13 is in a new cluster, so T1 is a positive-control list and not part of the test), T2 named (321 of the 725), T3 PANTHER only (404 of the 725). Report SF12 (Sps2 / Sps22 type; 7 of the 725) apart from SF3 (Ecm33 / Pst1 type; 718), because it is an open question whether a custom Ecm33 model should call the sporulation paralogs.
3. **Separation to stage (not run):** cluster these accessions with the training members at 30 percent identity and 50 percent coverage (the freeze's `cluster30` rule) and drop any protein in the same cluster as a training member. UniRef50 separation is only a first filter.
4. **How to compute recall with the frozen model (staging, no HPCC step run):** fetch the sequences of the chosen accessions from UniProt at job time into the work directory (sequences are not committed), run the frozen `all` custom model (trained on the 481 members, cutoff 29.6 from the build) with the same `hmmsearch` call and the frozen cutoff, and count the proteins with a score at or above it. Recall per fine clade = hits / proteins over cluster representatives, with the same comparators (Pfam GA, relaxed Pfam 15.0, the iteration comparator) so that u (members missed by all comparators) can be computed for the first time on non-circular members. The per-clade leave-out models are not the right models here: the set is independent of training by construction, so the `all` model is used, with clades only as a stratification. The false-call side stays with the frozen hard-negative rate (P38248 and Q4WNS8 stay out of that rate, they are in this set as true members).
5. **Use the PANTHER set a second way, as a purity check of the 481 members:** the 164 members outside any PANTHER cluster should be inspected before anyone reads 100 percent recall as meaning anything.

Caveats, to carry into any use: (i) PANTHER is itself an HMM, so recall on this set measures agreement between two profile models, not truth; a protein that PANTHER and the custom model both call could still be a non-Ecm33 member of a wider family (the 330 "Receptor L-domain" rows); (ii) PANTHER E-values above 1e-5 (76 of the 725) are weak calls; (iii) T3 rows have no name; (iv) the counts are for the InterPro and UniProt releases of 2026-10-23 and change; (v) 88 of 1,929 accessions had no cluster in the bulk id mapping and were filled by single queries, 0 left unmapped; (vi) a cluster in a clade with no training members (Entomophthoromycota, Mortierellomycota) tests extrapolation to a clade the model never saw, and should be reported apart; (vii) no sequence identity to training was computed, so "new cluster" is a lower bound on independence.

What we would need to stage: the 725-protein (or 563) accession list as a work-directory input; an MMseqs2 `cluster30` job against the training members; sequence fetch by accession; an `hmmsearch` of the frozen `all` model, plus the three comparators; a recall table by clade and tier. All are small; none was run.

## 6. hyd1_like (hydrophobin-like): can the same approach give an independent set?

Not directly. The PR 6 query is A0ACE4KFX7 (98 aa, Histoplasma, Eurotiomycetes) and it carries no Pfam, InterPro or PANTHER signature in UniProt. InterPro does have large hydrophobin entries (IPR001338 Class I, PF01185 and SM00075: 5,621 fungal proteins; IPR010636 Class II, PF06766 and PANTHER PTHR42341: 1,620 fungal proteins), so a second grouping exists for classical hydrophobins. But the hyd1_like cluster was chosen because it is a gap (no Pfam hit), so classical hydrophobin entries are a different family, and a set built from them would test something else. A family-level independent set for hyd1_like would need either OrthoDB groups that contain A0ACE4KFX7 (not looked up here) or a PANTHER/InterPro signature for it, which none has; the same OrthoDB route as for Ecm33 (group lists by genome ids) is the only candidate and its fungal-level group would need checking for how broad it is. This was not done.

## 7. Not done

No sequence was downloaded, no HPCC or SLURM step run, no model scored, no change to the freeze, the calls or any status. SGD, CGD, Alliance and eggNOG ortholog lists were not obtained (section 2). Loci in the papers were not mapped to accessions beyond three spot checks.
