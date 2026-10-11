# Claims I could not check (lectin validation, 2026-10-11)

Each line is a claim in the tables or the report that I did not verify. Items that came from earlier work are marked "inherited".

## Positives

1. **UniLectin3D ligands.** For the 39 `structure_backed` rows, the sugar component ids come from the RCSB non-polymer entity names (UniLectin3D record and RCSB GraphQL, fetched 2026-10-11). I did not check in the structure papers that the sugar sits in the binding site. NAG and NDG entries can be glycosylation. The tier is T2 for this reason.
2. **UniLectin3D PMIDs.** The PMIDs in the `source` column that say "from the UniLectin3D record, abstract not read" were not read.
3. **A7UNK4 (LDL, *Lyophyllum decastes*).** The re-fetch of UniLectin3D did not return this accession (API cap). The structure count (4) comes from the committed table of the spec work. I did not see its PDB ids or ligand.
4. **Swiss-Prot evidence codes.** Where a row cites a UniProt FUNCTION statement with ECO:0000269 and PMIDs, I read the UniProt text. I did not read the cited papers, except where the PMID is in the `source` column with an abstract quote. Examples: AAL (P18891), Mg1LysM follow-up papers 32574207 and 33797163, Q06100 (the abstract of PMID 8999822 does not say which isolectin was expressed).
5. **Vd2LysM (G2X4U8).** That the lineage-specific effector in the abstract of PMID 27911046 is Vd2LysM comes from the UniProt record. I did not read the full text.
6. **T. rubrum LysM1 and LysM2 as H1a.** From the abstract of PMID 30947244 only. I did not read the full text. I did not check which human glycoproteins were used.
7. **Paracoccin quotes.** The abstracts of PMIDs 16476564, 24743161 and 25474158 were read in PubMed. The statements about the GlcNAc site and the size conflict (70 kDa versus 27 to 28 kDa) are inherited from the L3b report (full text read there).
8. **FleA role.** The dispute (Sakai 2019, PMID 29370403) is inherited. Sakai 2019 is not verified in the PDF.
9. **Species names.** The species is the species-rank NCBI name above the UniProt organism. For A0A1U7Q1Z0 the UniProt organism is the genus *Psathyrella* (taxon 71729), so the species field holds a genus.

## Negatives

10. **ER house-keeping label.** The reason is the spec decision (row M3 of spec 3.1) and Lebreton 2021 (PMID 34200153, inherited). Twelve of the 18 rows have FUNCTION text by similarity or without an evidence code. Six have ECO:0000269 (P43555, Q12396, P27825, Q99220, A0A0D1C6P2, J9VLH0).
11. **Enzyme class.** The 8 alpha-galactosidases have evidence by similarity or none. I did not check whether the ricin-B-like module of these enzymes binds sugar. UniProt gives them the keyword Lectin.
12. **N1 rows.** Inherited from the shared set. I did not re-read their sources. I did check the sequences hit no model (that is a result, not a check of the labels).
13. **Assumed rows.** No reason exists. UniProt Pfam cross-references come from UniProt release 2026_03, not from the Pfam 38.2 file used for the scan.

## Method

14. **Seed leakage.** The seeds are the current InterPro seed alignments. The version string of each of the 20 models equals the version in the local Pfam 38.2 file. `none` means no hit at 30% identity by mmseqs2 `-s 7.5`. Remote homology can be missed.
15. **Q8NJT4 and Q8TGE0 as duplicates.** I judged them to be the same protein as Q4WW81 and Q2UNX8 from the protein names, species and PF07938 scores (484.1 against 489.1; 518.4 against 525.1). I did not align them.
16. **D4AL26 and the Onygenales group of spec F27.** D4AL26 (365 aa) is 92.3% identical to F3F37487_006503-T1 (337 aa) of *Trichophyton benhamiae* CBS 112371 (mmseqs2). It carries the same strain as D4AL26 (CBS 112371). I did not check whether they are the same gene.
17. **Wilson intervals** treat each cluster as independent. Clusters at 30% identity are not a proof of independence.
18. **Truncated gene models.** I looked for missed positives that are truncated. I found none among the 17 misses, but I did not compare lengths to full-length family members.
