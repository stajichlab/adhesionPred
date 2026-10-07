# Findings from agent-B curation pass (repeat-mechanism controls)

Reviewer: agent-B (opencode/glm-5)
Date: 2026-10-06
Output: `curation_table.agent-B.tsv` (23 columns = source 22 + `evidence_basis`)

## Method

- All 84 unique PMIDs cited by the 96 E1/E2 adhesin rows were fetched from PubMed and checked against the row's protein (identity checks via UniProt REST where needed).
- Full texts were additionally read for PMIDs 17337634 (PMC), 20656913 (PMC), 21283544 (PMC), 32286952 (PMC), 39455573 (PMC), 40595627 (PMC), 37769084 (PMC).
- Labels come only from the row's cited papers. Where a cited paper does not support the row's protein, the row is labeled from family evidence (`evidence_basis=family_inference`, confidence `low`) and listed below.
- Quote conventions: quotes are verbatim from the cited paper (abstract or full text of the same PMID); `[...]` marks an elision between two verbatim spans. For rows without PMIDs, the quote is the row's own `evidence_summary` text (the family evidence the label rests on); such rows have an empty `mechanism_source_pmid`.
- `evidence_basis` values: `stated_in_paper` = label grounded in a cited paper; `family_inference` = label rests on the row's family/domain evidence (no usable per-protein citation). `background_knowledge` was not used.
- Sheppard 2004 (15128742) full text: no open-access copy exists online (paywalled; no PMCID; Cloudflare-blocked publisher sites; nothing in Europe PMC/Unpaywall/OpenAlex/Wayback). The owner supplied the PDF (`to_import/PIIS0021925819710679.pdf`); text extracted locally with pypdf and used for the ALS6 and ALS9 rows (findings 7 and 23).
- One owner-supplied source was used: PMID 17600078 for the ALS9 row (not present in the row's `pmids` field). It was verified against PubMed before use.

## A. Rows where the cited paper does not support the row's adhesion claim

1. **A0A2H1A5L3** (E1, *Candidozyma auris*, PA14-family protein B9J08_02633, 2608 aa). Sole PMID 40595627 studies Als4112 (= B9J08_004112), not this protein. The full text does not mention B9J08_02633. The row's GO IMP annotation cannot be traced to this paper. Labeled `unknown/low/family_inference`.
2. **A0A1D8PJ68** (E2, *C. albicans* orf19.9293). Sole PMID 29726301 studies ORF19.1725 (and orf19.5557 etc.); orf19.9293 is not mentioned. Labeled `unknown/low/family_inference`.
3. **Q5A312 / ALS7** (E1, *C. albicans*). Sole PMID 17510860 studies Als1p and Als5p only; Als7 is not mentioned. Labeled `unknown/low/family_inference`.
4. **Q5AL03 / HYR1** (E1, *C. albicans*). PMID 21283544 assigns Hyr1 a role in resistance to neutrophil killing ("Overexpression of HYR1, but not HWP1, in the bcr1/bcr1 background, significantly rescued the higher susceptibility phenotype of this mutant to killing by leukocytes"). The paper contains no adhesion assay for Hyr1. The GO IMP adhesion annotation is not supported by this paper. Labeled `unknown/low/family_inference`.
5. **P41746 / rodA** (E2, *Aspergillus fumigatus*). PMID 9742203 is a study of the *Magnaporthe* conidiation regulator ACR1; RodA is not mentioned. Labeled `2c/low/family_inference` from the row family "class I hydrophobin".
6. **P52751 / MPG1** (E2, *Magnaporthe*). In PMID 9742203, MPG1 appears only as a hydrophobin gene whose expression fails to turn off in acr1(-) spores; no adhesion function is shown. Labeled `2c/low/family_inference`.
7. **A0A1D8PQ86 / ALS9** (E1, *C. albicans*) — resolved with owner-supplied evidence. The 15128742 abstract does not name Als9, but the owner supplied the paper PDF and the full text confirms Als9 was tested: "Als9p-expressing S. cerevisiae adhered above background levels only to laminin". The owner also supplied PMID 17600078 (Zhao, Oh, Hoyer 2007, Microbiology 153:3476-85, "Unequal contribution of ALS9 alleles to adhesion between Candida albicans and human vascular endothelial cells"), which is ALS9-specific: "Deletion of ALS9 significantly reduces C. albicans adhesion to human vascular endothelial cell monolayers", the ALS tandem-repeat architecture is stated, and a gene fusion carrying the ALS9-2 5' domain with the ALS9-1 tandem repeats and 3' domain restores wild-type adhesion. Labeled `2a/medium/stated_in_paper` from 17600078. 17600078 is not in the row's `pmids` field (`15128742;17510860`) - suggest adding it. 17510860 remains Als1p/Als5p-only (mismatch).
8. **O74670 / MSG** (E2, *Pneumocystis*). PMIDs 32363390 and 32561915 assert adhesion only as background ("...with importance in adhesion and immune recognition"; "MSGs are presumably responsible for antigenic variation and adhesion to host cells"). Neither paper contains an adhesion assay. Labeled `unknown/low/stated_in_paper` (assertion-level support).
9. **Q6FUW3 / EPA3** (E1, *C. glabrata*) — indirect support only. PMID 30348666 shows an Epa3 requirement for biofilm formation and azole accumulation; the abstract contains no direct adhesion assay. Labeled `unknown/low/stated_in_paper`.

## B. PMID mismatches inside otherwise-supported rows (minor)

10. **Q59L12 / ALS3**: among 11 cited PMIDs, 17510860 studies Als1p/Als5p only.
11. **Q5A8T4 / ALS1**: cited 15042589 (survey of potential cell-surface proteins; mentions Als1 only as a known adhesin) and 22359502 (adherence-regulator network screen) contain no ALS1-specific mechanism. The 2a label rests on 15116430.
12. **Q6FUW5 / EPA1**: cited 39349750 is about *C. glabrata* yapsin aspartyl proteases and is unrelated to Epa1. The `other/high` label rests on 10417386.

## C. Notes affecting mechanism labels

13. **A4D962 / BAD1**: cited abstracts (12023375, 15585870) demonstrate receptor binding ("BAD1 binds yeast to macrophages (Mphi) via CR3 and CD14...") but do not mention tandem repeats. The row's family field says "BAD4 tandem repeat (EGF-like)" (sic - likely BAD1). TOOL-ARCHITECTURE lists BAD1 as a 2a exemplar; the cited papers do not document repeat-mediated adhesion, so the row is labeled `other/high`.
14. **P32323 / AGA1**: cited papers document Aga1p as the anchorage subunit of a-agglutinin ("...consists of an anchorage subunit Aga1p and a receptor binding subunit Aga2p"), not a repeat mechanism. Labeled `other/high`. TOOL-ARCHITECTURE lists AGA1 as a 2a exemplar (repeat coverage 0.75); owner may want to re-check that exemplar choice against the cited papers.
15. **Q2LC49 / MAD1**: the cited paper (17337634, full text) documents "six tandem repeats comprised of 12 residues, GKETTPAQQTTP" in the Thr-rich middle domain - labeled `2a/medium`. The CFEM domain that UniProt annotates on MAD1 is not mentioned in the paper.
16. **Q9HDL9 / gp43**: `2d/medium` combines the row's glucanase annotation with ECM-receptor evidence in 10962270/16698299 ("gp43 bound both fibronectin and laminin"). The cited abstracts do not directly demonstrate moonlighting of an active enzyme.
17. **Q96V71 / SOWgp82**: PMID 12065484 predates the *Coccidioides immitis*/*C. posadasii* split; isolate C735 is called "C. immitis" in the paper.
18. **E9P9G2, C4R2D7, C4XZ24 / FLO11 rows**: 32286952 documents homotypic Flo11A-domain adhesion (labeled `other/high`) but contains no repeat evidence. The repeat-length paper 19160455 is cited only in the P08640 row. Suggest adding 19160455 to the other FLO11 rows if repeat evidence is wanted for them.
19. **A0A2H0ZI42 / IFF4109** (E2): the row cites no PMID. However, 37769084 - already cited by the SCF1 row - demonstrates Iff4109 adhesion (full text: "...discovered an uncharacterized adhesin, Surface Colonization Factor (SCF1), and a conserved adhesin, IFF4109, that are essential for colonization of inert surfaces and mammalian hosts"; "IFF4109, but not SCF1, mediates adhesion through cell surface hydrophobicity"). UniProt IFF09_CANAR also cites 37769084. Suggest adding 37769084 to this row; that would support an `other` label (hydrophobicity-mediated adhesion) instead of the current family-inference `unknown`.
20. **pfl rows (Q874R4, Q9URU4, Q8TFG9, P0CU05, Q7Z9I1, Q92344, Q9P5N1, Q9P7Q2)**: PMID 23236291 demonstrates flocculation for the pfl genes collectively ("Overexpression of the pfl(+) genes singly was sufficient to trigger flocculation"); individual pfl genes are not distinguished in the abstract. All labeled `unknown/low`.
21. **Hyr1-family E2 rows (23 rows, PF11765)**: the E1 anchors (AWP2, HYR1, IFF4) are themselves `unknown/low` - the family has no demonstrated mechanism. All labeled `unknown/low/family_inference`.
22. **P38894 / FLO5 and P39712 / FLO9**: support comes from 19420680, which expresses S. cerevisiae FLO1/FLO5/FLO9/FLO10 in K. marxianus and observes flocculation; 19087208 (also cited) mentions Flo10/Flo11 only. No FLO5/FLO9-specific mechanism is stated. Labeled `unknown/low`.
23. **Q5A2Z7 / ALS6**: upgraded `other/medium` to `other/high` using the owner-supplied Sheppard 2004 full text: "Als6p expression results in adherence to gelatin alone", and the Als6-NT/Als5-CT chimera "adhered only to gelatin, as did S. cerevisiae expressing Als6p" - Als6p-mediated, N-terminal-domain-determined substrate-specific adherence is stated and demonstrated. The tandem-repeat and C-terminal portions of Als5p/Als6p are "virtually identical", so the paper argues against a repeat-based mechanism for the Als5p/Als6p functional difference.
24. **Q5A8T7 / ALS5**: 12200964 refers to the protein as Ala1/Als5p (ALA1/ALS5 nomenclature).
25. **Candida_ALS (PF05792) E2 rows (A0AAW0V956, A0AAW0V9U7, P0CU38)**: PF05792 is the ALS tandem-repeat domain (UniProt Q5A8T4/ALS1: PF05792 present in 12 copies, matching the 12 UniProt REPEAT features). The `2a/low` labels assume the repeat mechanism documented for E1 anchors ALS1/ALS3/ALS4112; per-protein repeat numbers were not verified.
26. **Q6FTA2 (Candida_ALS_N, PF11766)**: PF11766 is the ALS N-terminal (substrate-binding) domain, not the repeat domain. Labeled `other/low` from the anchors' N-terminal-domain mechanism (ALS5 amino-acid-patch recognition; ALS1 N-terminal binding region).
27. **Q59TP1 / RBT1 (Flo11, PF10182)**: PF10182 is the Flo11A N-terminal domain (UniProt Q59TP1: 55..278), not a repeat domain. Labeled `other/low` from the Flo11A homotypic-adhesion mechanism documented for the FLO11 anchors; this is not a repeat-domain inference.
28. **Flocculin-repeat (PF00624) E2 rows (Q6FR82, Q6FQC7, Q6FUW7)**: PF00624 is a tandem-repeat domain shared with FLO1 (2a/high). Labeled `2a/low/family_inference`; per-protein repeats were not verified.

## Label distribution (96 rows)

| label | count | | confidence | count | | basis | count |
|---|---|---|---|---|---|---|
| 2a | 19 | | high | 25 | | stated_in_paper | 57 |
| 2b-iii | 1 | | medium | 10 | | family_inference | 39 |
| 2c | 3 | | low | 61 | | background_knowledge | 0 |
| 2d | 2 | | | | | | |
| other | 20 | | | | | | |
| unknown | 51 | | | | | | |
