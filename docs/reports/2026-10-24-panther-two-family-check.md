# Hand check of two PANTHER families: PTHR46430 and PTHR35192 (2026-10-24)

Read-only desk research. Public InterPro and UniProt web services and four PubMed abstracts. No model was run, nothing here says how well any model performs, and no freeze, call, tier or panel was changed. Tables and exact API calls: `docs/reports/data/panther_two_family_check/` (README there).

Terms, once each. **PANTHER**: a database of protein families defined by hidden Markov models (HMMs, statistical profiles of a protein family), grouped into families (PTHRnnnnn) and subfamilies (PTHRnnnnn:SFnn). **Pfam**: another HMM database, of protein domains. **InterPro**: a portal that merges the member databases and lists which entries match a protein. **Swiss-Prot (reviewed)**: the hand-curated part of UniProt. **Sel1-like repeat (Sel1, PF08238)**: a short ~36 residue repeat that stacks into a curved helical surface; related to the tetratricopeptide repeat (TPR); found in many unrelated proteins. **Fine clade**: the clade rule of the project (class for Pezizomycotina and Basidiomycota, subphylum for other Ascomycota, phylum otherwise). **Kiss cluster**: a cluster of the Kiss 2019 gene list used as truth set.

## 1. Verdicts

| Family | Verdict | One-line reason |
|---|---|---|
| PTHR46430 Chitin Synthase Regulator | Do not use as one set. Broad as a whole. Usable only as a descriptive, subfamily-level list, and only with the Sel1 caveat. | 98 percent of fungal members are defined by the Sel1 repeat; the family splits into three subfamilies of which only SF1 (Skt5-related) and SF2 (Chs4 factor 4) are chitin or glucan synthase regulators by literature. |
| PTHR35192 "PROTEIN, PUTATIVE-RELATED" | A real, specific, basidiomycete-only family (name is a placeholder; it is the CPL1/PriA family). Usable as a candidate positive set for "CPL1-like secreted proteins", but not as an adhesion set. | 100 percent fungal, 5 clades, one InterPro family (IPR038955) and one Pfam domain (PF21671, 96.8 percent). Function in the literature is a virulence effector and a fruiting-related protein, plus (via Kiss) an adhesin. |

## 2. PTHR46430

### 2.1 Size and taxa (UniProt, 2026-10-24; InterPro counts in brackets where they differ)

| Group | Proteins | Share |
|---|---|---|
| All taxa | 5,247 (InterPro 5,476) | |
| Fungi | 4,935 | 94.1 % |
| Bacteria | 246 (229 Pseudomonadati, 15 Bacillati, 2 others) | 4.7 % |
| Archaea | 1 | |
| Other eukaryotes | 65 (57 with no kingdom, 4 Metazoa, 4 Viridiplantae) | 1.2 % |

Fungal reviewed (Swiss-Prot): 10. The earlier scoping said 246 bacterial proteins; this confirms it (246).

### 2.2 Subfamilies (InterPro `?subfamilies`: 3 subfamilies; UniProt also holds 269 proteins with the family but no subfamily)

| Subfamily | Name (PANTHER) | All | Fungal | Fungal share | Swiss-Prot | Bacteria | Fungal clades with >= 5 |
|---|---|---|---|---|---|---|---|
| PTHR46430:SF1 | CHITIN SYNTHASE REGULATOR SKT5-RELATED | 1,799 | 1,708 | 94.9 % | 6 | 77 | 23 |
| PTHR46430:SF2 | CHITIN SYNTHASE REGULATORY FACTOR 4 | 1,532 | 1,487 | 97.1 % | 1 | 32 | 22 |
| PTHR46430:SF3 | ACTIVATOR OF C KINASE PROTEIN 1 | 1,647 | 1,606 | 97.5 % | 3 | 35 | 24 |
| (family only, no subfamily) | | 269 | 134 | 49.8 % | 0 | 102 | 7 |

Most of the 246 bacterial proteins sit in SF1 to SF3 (144) or have no subfamily (102). Bacterial members are therefore not confined to a loose fringe: the family-only group is half bacterial, and the three subfamilies each carry 32 to 77 bacterial proteins.

### 2.3 The Kiss proteins

| Kiss cluster | Gene, protein | Subfamily | Step 6.0 Pfam hit (E value, share of the Pfam model covered) |
|---|---|---|---|
| KS0084 chs4_skt5 | P34226, S. cerevisiae Skt5 (Chs4) | SF1 Skt5-related (specific type) | Sel1 PF08238, 4e-08, 0.353 |
| KS0065 chr4 | O94614, S. pombe Chr4 | SF2 Chitin synthase regulatory factor 4 (specific type) | Sel1 PF08238 6.7e-05, 0.264; TPR_2 PF07719 and TPR_8 PF13181 |

So the two Kiss proteins fall in two different, narrower subfamilies, not in a "broad" bucket; but neither subfamily name is the same as the other, and S. pombe Chr4 (the functional homolog) does not group with S. cerevisiae Skt5 at subfamily level. No other Kiss accession (112 clustered accessions plus 44 missed clusters) is in PTHR46430 (`kiss_in_family.tsv`).

### 2.4 Domains of the fungal members (`domain_counts.tsv`)

| Entry | Fungal members with it |
|---|---|
| Sel1-like repeat IPR006597 | 4,883 of 4,935 (98.9 %) |
| Pfam Sel1 PF08238 | 4,843 (98.1 %) |
| TPR-like helical superfamily IPR011990 | 4,921 (99.7 %) |
| Tetratricopeptide repeat IPR019734 | 467 (9.5 %; 30.9 % in SF2, under 2 % in SF1 and SF3) |
| Other Pfam repeat models (TPR_1, TPR_6, TPR_12, TPR_16, TPR_8, TPR_2) | under 0.2 % |
| ChAPs (Chs5-Arf1-binding proteins) | not found as an entry on any member |

Per subfamily the fractions are 96 to 99 percent Sel1 in each. The members are Sel1-repeat proteins; the PANTHER family is a group of Sel1 proteins, and no second, family-specific domain exists to separate them. Bacterial and other non-fungal members are 88.8 percent Pfam Sel1 plus some MYND zinc finger (8.7 %) and peptidoglycan-binding domain (4.5 %) proteins.

How broad is the Sel1 model itself? Fungal proteins with PF08238: 17,527; of these only 4,832 (27.6 percent) are in PTHR46430 (`pfam_vs_family_counts.tsv`). A Sel1 hit therefore says little about Chs4-type membership.

### 2.5 Swiss-Prot members (`reviewed_members.tsv`, 10 fungal; 0 non-fungal)

| Subfamily | Protein | Organism | Curated function |
|---|---|---|---|
| SF1 | Skt5 / Chs4 P34226 | S. cerevisiae | Activator of chitin synthase CHS3 (experimental) |
| SF1 | Csr2 J9VSG5 | C. neoformans | Activator of CHS3 (experimental) |
| SF1 | Shc1 P39000, A6ZR53 | S. cerevisiae (two strains) | Needed for CHS3 activity in sporulation (experimental) |
| SF1 | Chr1 O59732 | S. pombe | "Involved in chitin biosynthesis" (by inference, ECO:0000305) |
| SF1 | Skt5 P9WEM4 | Malassezia restricta | CHS3 activator, by similarity |
| SF2 | Chr4 O94614 | S. pombe | Keeps Bgs1 (1,3-beta-glucan synthase) at poles and septum under wall stress (experimental) |
| SF3 | Ack1 Q07622 | S. cerevisiae | Cell wall integrity pathway; localizes the Rho1 exchange factor Rom2 (experimental) |
| SF3 | Cfh1 Q09897, Cfh2 O94486 | S. pombe | "May play a role in cell wall organization", by similarity to Chr4 |

All reviewed members are cell wall regulators of some kind, which is independent support for a wall-related grouping. But the support is by name, and only 3 or 4 of them are experimentally tied to chitin or glucan synthase activity; Ack1 acts on a signalling pathway (Rho1 localization), not on a synthase.

### 2.6 Is a Chs4 (Skt5) regulator subfamily separable from generic Sel1 proteins?

By PANTHER, yes: SF1 (1,708 fungal proteins, 6 Swiss-Prot with CHS3-activation functions) and SF2 (1,487 fungal, 1 Swiss-Prot). By domain, no: members all carry Sel1 repeats in the same share, so a domain model cannot separate them; only a family-level or subfamily-level HMM can. Whether SF1 and SF2 contain only true synthase regulators is not established: they are 94.9 and 97.1 percent fungal, with 77 and 32 bacterial proteins, and their sizes (1,708 and 1,487 fungal proteins, 3,195 together) look larger than one regulator per genome would give, so paralogs or distant Sel1 proteins are probably included. I did not count proteomes, so this is an expectation, not a measurement.

### 2.7 Verdict PTHR46430

Broad, with a possible narrow core. Do not use the family as one test set: it is defined by a promiscuous repeat, has a bacterial component, and its third subfamily (SF3, Ack1/Cfh) is signalling or "by similarity". If the owner wants a Chs4 test set, use SF1 only, plus SF2 for S. pombe type regulators, as a descriptive set, and never score a Sel1 hit as evidence. Unclear: whether SF1 members outside a few experimentally studied species actually activate chitin synthase.

## 3. PTHR35192

### 3.1 Size and taxa

| Group | Proteins |
|---|---|
| All taxa (UniProt) | 1,521 (InterPro 1,598) |
| Fungi | 1,521 (100 %) |
| Bacteria, Archaea, Metazoa, plants | 0 |
| Swiss-Prot | 2 |
| Subfamilies | 1 (PTHR35192:SF2, named "APPLE DOMAIN-CONTAINING PROTEIN", all 1,521 proteins) |

Fungal clades (>= 5 proteins): Agaricomycetes 601, Tremellomycetes 580, Pucciniomycetes 169, Microbotryomycetes 130, Dacrymycetes 38; Mixiomycetes 3. All Basidiomycota.

### 3.2 What the family is

The family name "PROTEIN, PUTATIVE-RELATED" is a PANTHER placeholder; it is not a catch-all of unrelated proteins. Evidence: all 1,521 proteins carry InterPro family IPR038955 "Protein PriA/CPL1"; 1,472 (96.8 percent) carry the Pfam domain PF21671 "CPL1-like"; none carries any other Pfam or InterPro entry at 2 percent or more (49 have no Pfam entry). InterPro describes it as a group of fungal proteins including PriA from Lentinula edodes and CPL1 from Cryptococcus neoformans. Fungal proteins with PF21671: 1,640; 89.8 percent of them are in the family. The subfamily name "Apple domain-containing" does not match the domain data; the Pfam entry named no apple domain. I did not resolve that discrepancy.

Swiss-Prot members: CPL1 (J9VQE5, C. neoformans, CNAG_02797, 199 aa): a secreted virulence effector that polarizes macrophages through Tlr4 (Dang et al. 2022, Nature, PMID 35896747, doi 10.1038/s41586-022-05005-4); priA (Q01200, Lentinula edodes): "may play a role during the beginning of fruiting".

### 3.3 The Kiss protein and Cfl1

Kiss cluster KS0056 gene cfl1, accession J9VHS9 (C. neoformans H99, CNAG_00795, 309 aa, unreviewed, UniProt name "Protein CPL1-like domain-containing protein") is in the family. It is not CPL1 itself (that is J9VQE5). It is the only Kiss accession in the family.

What Cfl1 is, per PubMed abstracts: a cell-surface and secreted protein, "the first identified adhesin in the Basidiomycota", downstream of the transcription factor Znf2, regulating morphology, cell adhesion, biofilm formation and virulence (Wang, Zhai, Lin 2012, PMID 22737071, doi 10.1371/journal.ppat.1002765). Cfl1 has four paralogs in C. neoformans (Dha1, Dha2 and two more), all secreted, sharing a conserved C-terminal "SIGC" domain and an N-terminal signal peptide (Gyawali et al. 2017, PMID 28039134, doi 10.1128/AEM.02967-16). The family is therefore functionally varied: an adhesin, an immune effector, and a fruiting-related protein. I did not establish that the SIGC domain is Pfam PF21671, nor that CNAG_00795 is the Cfl1 locus beyond the Kiss table's own "locus_confirmed" mapping.

### 3.4 Verdict PTHR35192

Specific by composition (fungal only, one domain), real family, but function is mixed and the clade spread is narrow (Basidiomycota only, 5 clades with >= 5). It cannot support a clade-spread truth test beyond basidiomycetes and the project clade floor would still be met only for those five clades. As an adhesion set it is not supported: the only adhesin evidence is one paper on one protein, and the family's best-characterized member is an immune effector. Fit as a set for "CPL1/Cfl1-like secreted basidiomycete proteins", with a descriptive label.

## 4. Circularity and independent evidence

PANTHER families are HMMs built from sequence clustering, and the Pfam models are HMMs too. Finding that a Pfam model or a new model recovers a PANTHER family measures how far two profile models agree, not whether the grouping is biologically right. Independent evidence here is:

- PTHR46430: Swiss-Prot curated functions for 10 members, all wall-related, 4 or 5 with experimental support (Skt5, Csr2, Shc1, Chr4, Ack1). This supports "some of the family regulate wall synthesis" but not that all 4,935 do.
- PTHR35192: the PubMed work above (Cfl1 adhesin, CPL1 effector), and InterPro's own family text. No experiment ties the 1,521 members to one function.

## 5. What a specificity test would look like (existing Pfam models, accessions from `scan_per_cluster.tsv`)

| Cluster | Model hit by the cluster | Expected informativeness |
|---|---|---|
| KS0065 chr4 | PF08238 Sel1 (E 6.7e-05, model share 0.264); PF07719 TPR_2; PF13181 TPR_8; weak: PF00515 TPR_1, PF13174 TPR_6, PF13424 TPR_12, PF13432 TPR_16 | Low. Only 27.6 percent of fungal Sel1 proteins are in the family; a test would count Sel1 proteins as hits and show mostly non-Chs4 proteins. The TPR models hit 0 to 3 family members. The test would confirm that the repeat is not specific; it would not test a Chs4 model. |
| KS0084 chs4_skt5 | PF08238 Sel1 (E 4e-08, 0.353) | Same: low. |
| KS0056 cfl1 | PF21671 CPL1-like (E 6.8e-20, model share 0.184) | Moderate. 89.8 percent of fungal PF21671 proteins sit in PTHR35192, so the model is specific to the family. But it hits only 18 percent of the model length for this protein, and the family is only Basidiomycota; the test measures the CPL1 family, not adhesion. |

A model that would be informative for the Chs4 group is a purpose-built one from SF1 and SF2 members (not one existing Pfam model); that is a later PR 6 option, not part of this check.

## 6. What I could not establish

- Whether the 1,708 and 1,487 fungal members of SF1 and SF2 are all true chitin or glucan synthase regulators (only 7 are curated).
- Why the PTHR35192 subfamily is named "Apple domain-containing protein" and whether PF21671 equals the SIGC domain.
- Whether CNAG_00795 is the published Cfl1 locus (Kiss table says confirmed; not re-checked here).
- UniProt lags InterPro by 229 (PTHR46430) and 77 (PTHR35192) proteins; I did not look at which are missing.
- Species coverage per proteome (how many Chs4-type genes per genome).
