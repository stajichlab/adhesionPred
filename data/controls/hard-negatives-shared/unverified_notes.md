# Unverified notes and open questions (task 07)

*2026-10-08, task 07 agent (claude-opus-5-5). Nothing in this file is in `controls.tsv` as a fact.*

## Claims I could not check

1. **Aspergillus Wsc sensors.** *A. fumigatus* Wsc1, Wsc3 and MidA (PMID 22220813) and *A. nidulans*
   WscA and WscB (PMID 21926329) have characterised sensor functions. I could not map the gene
   names to accessions. UniProt has eight *A. fumigatus* and eighteen *A. nidulans* WSC or Mid2-domain
   entries without gene names. The full texts were not available to me. A guess such as
   "ANIA_05660 is wscA" is from memory, not checked. These rows are not in the table. The
   `mucin_sensor` stratum therefore has no non-yeast row.
2. **Coccidioides GEL1.** PMID 12761077 describes *C. posadasii* GEL1. I did not find its locus
   tag in the text that I read. Several GH72 entries exist (for example J3KBG5, J3KGG0, A0A0E1RUX2,
   A0ACG8DAS1). Not added.
3. **Histoplasma SOD3 accession (row A0ACF1AQZ4).** I took the accession from the UniProt gene name
   SOD3 in *Histoplasma ohiense* (G217B lineage). The paper (PMID 22615571) text that I searched does
   not state the strain or a locus tag. The length (232 aa) fits the C-terminal GPI signal near
   residue 205 that the paper describes. Check the mapping.
4. **Rsp3 accession (row A0A0D1DYI3).** PMC5923269 gives UMAG_03274. UniProt A0A0D1DYI3 has ORF name
   UMAG_03274. I did not check the repeat count against the paper.
5. **Possible adhesion phenotypes of homologs in *Candida*.** From memory, not checked: *C. albicans*
   ecm33, crh-family and phr/gas mutants may have altered adhesion in some papers. I did not search
   these. If true, the rule "function includes adhesion in another organism" may apply to the
   Ecm33, Crh (GH16) and Gel/Gas (GH72) rows. Only the cspA double-mutant result (PMID 20656913) for
   *A. fumigatus* ECM33 and GEL2 is checked; those two rows are flagged in `notes`.
6. **Systematic moonlighting check.** I applied the "adhesion in another organism" rule to the
   families where I found evidence (Msb2, Pra1/Aspf2, hydrophobins, Hsp60, Sub3, PLB1, Rep1,
   Cda/MP84, CalA). I did not search adhesion literature for every family in the table (for example
   laccases, aspartic proteases, chitinases, LPMOs, CFEM effectors). A row without a flag is not
   proven free of such reports.
7. **Task 08 ENG2 quote (row XP_045276491.1).** I read the abstract of PMID 41061037 and the quote
   matches it. I did not read PMID 40667097, which task 08 also cites.
8. **UniProt function text as the quote.** For 139 rows the quote is the UniProt FUNCTION comment,
   with ECO:0000269 PubMed evidence. I read 41 PubMed abstracts, not the abstract of every cited
   PMID. The UniProt curators made the link from paper to function.
9. **Pfam seed coverage.** I compared rows with the seed alignments of the 25 Pfam families in
   `family_table.tsv` (InterPro API, Pfam 38.2 seeds, 2,578 seed sequences). The leakage value uses
   the seed domain sequences. The full-length sequences of 74 seed members are no longer in UniProt.
   The search used the domain sequences, so these are covered for the domain only.
10. **SignalP 6 training data.** Overlap between these rows and the SignalP 6 training set is not
    measured (same limit as for R0 in `docs/paper/03-status-and-validation-rules.md`, decision 1).

## Open questions for the owner

1. **Msb2 family.** The task file names Msb2 as the example of the `mucin_sensor` stratum. The
   *A. nidulans* ortholog MsbA influences adhesion (PMID 25294314), and `manual_overrides.tsv` already
   classes msbA as `surface_other_adhesion_phenotype`. I rejected all Msb2 rows, including
   *S. cerevisiae* MSB2, to follow the task rule. Do you want Msb2 back as a flagged hard negative?
2. **What does `serves=step1` mean?** I read step 1 as "outside the plasma membrane" (docs/paper/01).
   A wall enzyme is then a step 1 positive. I marked `step1` only for signal-peptide or membrane
   proteins that are not at the surface (ER, vacuole, plasma membrane sensors; 9 rows). I added the
   value `adhesion_level` for the other rows. Please confirm.
3. **Strata that are not in the task list.** I added `secretory_non_surface` (5 rows, step 1
   look-alikes) and `wall_structural_other` (9 *S. cerevisiae* rows carried from `adhesins.tsv` that
   are not enzymes: CWP1, CWP2, DAN1, TIR1, TIR3, TIR4, EGT2, KRE9, FIT1). Keep, rename or drop?
4. **Mixed clusters.** EXG1 and *S. pombe* exg1 share a 30% cluster with gp43 (Q9HDL9, an E2
   adhesin that is also a GH5 glucanase). The four *A. fumigatus* PA14 beta-glucosidases share a
   cluster with A0A1D8PJB5 and A0AAW0V6D7, which `adhesins.tsv` labels E2 adhesin **only by the PA14
   domain**. The README of `data/curated/adhesins/` says PA14-only hits are E3 or hard negatives. I
   suggest that these two positives move to E3. Then the PA14 cluster is a clean negative cluster.
5. **N2 rows.** 24 rows are N2. They test architecture only. Should N2 rows count toward the
   20-cluster floor?
6. **Pooling strata in one species.** Only *S. cerevisiae* has 20 or more negative clusters. Do you
   accept a pooled-strata negative set for one call in one species?

## Weakest rows (for the PR text)

- The 4 *Drepanopeziza brunnea* CFEM rows (N2; UniProt says "Appears to function ... may play a role").
- The 4 *A. fumigatus* PA14 beta-glucosidases (N2; mixed cluster).
- EXG1 and exg1 (mixed cluster with gp43).
- *A. fumigatus* ECM33 and GEL2 (adhesion loss in double mutants with cspA).
- *H. ohiense* SOD3 (accession mapping not checked).
- *C. immitis* CTS1 (N2; homology to a purified antigen).
- Blastomyces ENG2 (N2; NCBI accession, no UniProt entry).
