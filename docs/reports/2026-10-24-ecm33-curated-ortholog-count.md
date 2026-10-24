# How many curated Ecm33 orthologs could make an independent test?

*2026-10-24. Author: Claude Code (claude-sonnet-5-5). Counting only: no model was run, no HPCC or SLURM step, no call, tier, panel or `narrow_ok` changed. Ids-only table and the exact API calls: `docs/reports/data/ecm33_curated_orthologs/`. Context: `2026-10-23-ecm33-independent-set-scoping.md` (section 2, row 1, said 13 Swiss-Prot entries) and `2026-10-24-pr6-ecm33-independent-eval.md`.*

## 1. Answer in short

- **Zero.** No curated Ecm33-family fungal protein is both outside the training set and separated from it by 30 percent identity. A genuinely independent, curated test set of Ecm33 does not exist today.
- Swiss-Prot (the hand-reviewed part of UniProt) has **13** reviewed entries in the PANTHER family PTHR31018. A fresh read-only UniProt query on 2026-10-24 (release 2026_03) returned the same 13 accessions as `candidates.tsv`. They are 9 distinct UniRef50 clusters, 5 species and 3 clades (Saccharomycotina, Taphrinomycotina, Eurotiomycetes). **All 13 sit in a UniRef50 cluster that is a training member**, so each is at least 50 percent identical to a training protein.
- Papers that name Ecm33 orthologs: 10 organisms in the 2026-10-23 scoping list. I resolved **7** to an accession (4 of them are Swiss-Prot entries of 4 organisms; 3 by name search in unreviewed UniProt); **3 could not be resolved** (Metarhizium robertsii, Aspergillus flavus, Podosphaera xanthii). Of the 7 resolved, 6 are in training clusters and the seventh (F. oxysporum N1RSC8) is in one of the 187 test clusters but shares a 30 percent cluster with 8 training members, so it is not separated either. **0 of 7 are independent.**
- The only truly separated Ecm33-like proteins are the 66 PANTHER test clusters that are separated at 30 percent (in 13 clades, 21 of them Entomophthoromycota). Their evidence is PANTHER agreement, not curation: 12 carry an Ecm33-like UniProt name (`T2_named`, unreviewed, many named automatically) and 54 are `T3_panther_only` (mostly "uncharacterized"). None is reviewed.

## 2. Words used (each defined once)

| Word | Plain meaning |
|---|---|
| Ecm33 family | Fungal cell-surface proteins anchored by a lipid tail (GPI anchor) that help build the cell wall; in yeast the family has Ecm33, Pst1, Sps2 and Sps22 (names of genes in baker's yeast). |
| Swiss-Prot, reviewed | Entries of UniProt that a curator read and checked. TrEMBL (unreviewed) entries are computer-made. |
| PANTHER PTHR31018 | A public family model (itself a profile HMM) that groups the Ecm33 proteins. It is a second method, not proof of function. |
| UniRef50 cluster | A group of proteins at least 50 percent identical to each other, represented by one protein. The training set is a list of such clusters. |
| 30 percent separated | No training protein shares a 30 percent identity / 50 percent coverage group (MMseqs2) with the test protein. Only then does a test show more than "finds close relatives". |

## 3. Swiss-Prot reviewed entries of PTHR31018 (13)

Re-verified: UniProt `reviewed:true AND xref:panther-PTHR31018`, fungi (taxonomy 4751), 2026-10-24: 13 results; the accessions equal the 13 `T1_swissprot` rows of `docs/reports/data/ecm33_independent_set_candidates/candidates.tsv` (checked by the test). A wider name/Pfam query returns 15, the 13 plus two namesakes that are not family (Q09750, a SIN3 homolog, and A0A1D8PHR5, a quinone oxidoreductase).

| Accession | Organism | Clade | Name | UniRef50 cluster | In training? |
|---|---|---|---|---|---|
| P38248, A6ZL22, B5VE42, C7GQJ1 | S. cerevisiae (S288c, YJM789, AWRI1631, JAY291) | Saccharomycotina | Ecm33 | UniRef50_P38248 | yes (4 accessions, 1 cluster) |
| Q12355, A6ZY20 | S. cerevisiae (S288c, YJM789) | Saccharomycotina | Pst1 | UniRef50_Q12355 | yes (2 accessions, 1 cluster) |
| P08459 | S. cerevisiae S288c | Saccharomycotina | Sps2 | UniRef50_P08459 | yes |
| P25380 | S. cerevisiae S288c | Saccharomycotina | Sps22 | UniRef50_P25380 | yes |
| Q5AGC4 | C. albicans SC5314 | Saccharomycotina | Ecm33 (the Kiss query) | UniRef50_Q5AGC4 | yes (the query itself) |
| O13960 | S. pombe | Taphrinomycotina | ecm33 | UniRef50_O13960 | yes |
| O74411 | S. pombe | Taphrinomycotina | Meu10 | UniRef50_O74411 | yes |
| Q4WNS8 | A. fumigatus Af293 | Eurotiomycetes | ecm33 | UniRef50_Q4WNS8 | yes |
| D4AZK9 | Arthroderma benhamiae | Eurotiomycetes | ARB_01627 | UniRef50_D4AZK9 | yes |

Counts: 13 accessions, 9 UniRef50 clusters (4 + 2 + 1 + 1 + 1 + 1 + 1 + 1 + 1 accessions), 5 species, 3 clades. Whether a cluster is a training member was read from `members.tsv` (481 members). Because the training set is the clusters themselves, a Swiss-Prot protein in such a cluster is at least 50 percent identical to a training member and cannot be separated at 30 percent. P38248 and Q4WNS8 are the two true Ecm33s among the shared negatives; they were left out of the false-hit count by the freeze, and they are training clusters too.

## 4. Literature-named proteins

PMIDs were re-checked in PubMed on 2026-10-24 (the title and organism match). Accession resolution by name is not the same as the locus in the paper; where I did it by name it is labelled so.

| Organism (PMID) | Accession | How resolved | In training? | Independent? |
|---|---|---|---|---|
| S. cerevisiae Pst1 / Ecm33 (15583168) | P38248, Q12355 | Swiss-Prot | yes | no |
| C. albicans Ecm33 | Q5AGC4 | Swiss-Prot (PMIDs of this organism are in the scoping, not re-read here) | yes | no |
| S. pombe ecm33 (20032302) | O13960 | Swiss-Prot | yes | no |
| A. fumigatus ecm33 (16804168) | Q4WNS8 | Swiss-Prot | yes | no |
| B. bassiana Bbecm33 (24549768) | J4UGA7 | UniProt name search ("Ecm33-like protein", PANTHER SF3); locus not checked | yes (UniRef50_A0A179IKX0) | no |
| T. virens Ecm33 (30840111) | G9MJ35 | UniProt name search; locus not checked | yes (UniRef50_G9MJ35) | no |
| F. oxysporum f. sp. cubense FocECM33 (35183338) | N1RSC8 | UniProt name search ("Protein ecm33"); no PANTHER call; locus not checked | no; test cluster UniRef50_A0A428P255, not separated (8 training members in its 30 percent cluster); the custom model scores it 256.7 | no |
| M. robertsii Mrecm33 (24549768) | none | UniProt name search for *M. robertsii* gave no entry; other Metarhizium species have "ECM33-like" entries (for example M. anisopliae C0LZ60), in a test cluster (UniRef50_A0A8H4PWU2) that is not separated (128 training members in its 30 percent cluster) | n/a | no |
| A. flavus AFLA_113120 (29696338) | none | no Ecm33 entry found by name or locus | n/a | n/a |
| P. xanthii PxECM33 (42325904) | none | no entry found by name | n/a | n/a |

So: 10 organisms named, 7 resolved (4 reviewed, 3 by name), 3 unresolved, 0 independent. All 7 resolved organisms are in clades that are already in training (Saccharomycotina, Taphrinomycotina, Eurotiomycetes, Sordariomycetes). Candida glabrata has two unreviewed "Cell wall protein ECM33" entries (Q6FK20, Q6FV55), both in training clusters (no paper was mapped, so they are not in the table above). In the ids-only table the five Swiss-Prot proteins that also have a paper appear in both sets, so the table has 21 rows, not 16 distinct proteins.

## 5. Conclusion in numbers

| Question | Answer |
|---|---|
| Swiss-Prot reviewed PTHR31018 entries | 13 (9 UniRef50 clusters, 5 species, 3 clades) |
| of them not in a training cluster | 0 |
| literature-named organisms resolved to an accession | 7 of 10 |
| of them separated from training at 30 percent | 0 |
| separated PANTHER test clusters | 66 of 187 (12 `T2_named`, 54 `T3_panther_only`), 13 clades; none reviewed |

A curated, independent test of 0 proteins is possible now. The family is conserved, so every curated protein known has a close relative in training. If the owner wants a test that is both independent and curated, the options are (a) to curate some of the 12 separated `T2_named` clusters by hand (read the papers, check the gene), (b) to ask a curator or database (Swiss-Prot, OrthoDB with accessions) to add reviewed entries from clades outside Saccharomycotina, Taphrinomycotina and Eurotiomycetes, or (c) to accept PANTHER agreement as the only independent line, which is what the 2026-10-24 evaluation did.

## 6. Caveats

- UniProt and InterPro counts change every release (this read: UniProt 2026_03). The 2026-10-23 scoping used the InterPro API, which lists 1,929 fungal proteins; UniProt's PTHR31018 cross-reference list is older and shorter, but the reviewed set is the same 13.
- A UniRef50 cluster id can change between UniRef releases; the training set (2025_03) and the live ids (2026) were compared by id, as in the evaluation record.
- Paper-to-locus matching was not done for the three name-resolved proteins; a wrong match would change nothing here because each is in a training cluster or a non-separated test cluster.

No tier, panel, call, `narrow_ok` or status change is made by this record.
