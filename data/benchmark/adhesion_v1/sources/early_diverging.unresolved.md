# Early-diverging fungi: candidates NOT put in the table

Curated 2026-11-26; second pass 2026-11-26 after owner feedback. Scope: Mucoromycota,
Zoopagomycota, Chytridiomycota, Blastocladiomycota, Microsporidia, Cryptomycota. Every PMID below
was confirmed with the PubMed tools and every accession with UniProt REST or NCBI eutils; items are
here because one of those checks failed or because the evidence does not meet the inclusion rule.

Literature retrieved through PubMed. DOIs are given inline.

Owner's tiering rules applied on the second pass: moonlighting cytoplasmic proteins and proteins
with only indirect adhesion phenotypes belong in separate tiers, and a null adhesion result alone
does not make a hard negative.

## Resolved on the second pass (now in the table)

| Candidate | Accession | Label | How the accession was recovered |
|---|---|---|---|
| Bd CBM18-LL1 | A0A097F7X9 (EMBL KM099424 / AIT18659.1), 465 aa | hard_negative | Direct UniProt query on taxon 109871; the entry's gene name is literally `CBM18-LL1` and it comes from a Stajich & Liu deposit filed under the pre-publication title of PMID 25819009. The closed-access table was never needed. |
| NbSWP16 | R0MN98 (EOB14338.1), 383 aa | adhesin | Table 1 of the open-access review PMID 32132983 (PMC7040029) lists NbSWP16 as 44.0 kDa / GenBank EOB14338; verified live and corroborated by the sequence itself (signal peptide, basic HBM stretches, exactly three copies of the proline-rich repeat EPPHHHPKPCPPC). |
| AlocSWP2 | Q6E6F1 (EMBL KX255658 / API85832.1), 222 aa | hard_negative | Found while checking the same review table; the GenBank record carries /gene="SWP2" and /note="AlocSWP2" and was deposited by the authors of PMID 28277606. |
| NbHSWP11 (already in the table, flag cleared) | ABV48899.1, 446 aa | adhesin | No longer inferred: the GenBank CDS carries `/gene="HSWP11"`, and PMID 32132983 Table 1 lists NbSWP11 as 446 aa / EF683111 / Yang et al. 2014 (= PMID 24813020). |

## Still unverifiable after a second attempt

| Candidate | Organism | Why still unresolved |
|---|---|---|
| EcRBLL-1, AaRBLL-1, AaRBLL-2 | *Encephalitozoon cuniculi*, *Anncaliia algerae* | PMID 34808329 ([DOI](https://doi.org/10.1016/j.parint.2021.102518)) is still not in PMC and Europe PMC reports `isOpenAccess: N`, `inEPMC: N`, `hasSuppl: N`. Europe PMC advertises a "Free pdf" author-manuscript link at sciencedirect.com/science/article/am/pii/S1383576921002361, but that URL returns **HTTP 403 Forbidden**. No preprint exists: a Europe PMC search across MED and PPR for RBLL + Anncaliia/Encephalitozoon returns no preprint version. Per instructions, no paywall was paid for or circumvented, and the attempt stops here. The abstract names no locus, and guessing which of the many ricin-B-lectin genes is "RBLL-1" would be inventing an accession. **Only the full text or an author deposit can close this.** |
| Bd CBM18-LL2, and the tyrosinase-like and deacetylase-like CBM18 members | *B. dendrobatidis* | Only CBM18-LL1 was ever deposited. GenBank `lectin-like protein AND Batrachochytrium[orgn]` returns exactly one record (KM099424); the accession block KM099420-KM099435 contains no sibling Bd deposit; UniProt has one `CBM18*` gene name for taxon 109871. The other family members exist only as BDEG locus tags inside the closed-access table of PMID 25819009. |

## Positive but indirect — the owner's separate tier, not a row here

- **NcSWP8** (*Vairimorpha (Nosema) ceranae*, PMID 39858910,
  [DOI](https://doi.org/10.3390/microorganisms13010142)). Accession IS resolvable: locus NCER_100151
  = GenBank **EEQ83062.1**, 189 aa, which matches the paper's stated 21.37 kDa (UniProt
  A0A0F9WIV3 is a 172 aa entry from a different assembly and is NOT the paper's protein). Held back
  on purpose: the only host-interaction result is an antibody-blocking reduction in **infection
  rate**, not in adherence, and the direct results (co-IP and Far-Western with NcPTP2/NcPTP3, spore
  wall localisation) are structural. Calling it an adhesin overreaches the rule; calling it a hard
  negative risks burying a real adhesin, because its signal is positive-but-indirect rather than
  null. It belongs in the indirect tier. This is ready to promote the moment the owner decides which
  tier takes it — the accession is verified.
- **NbRBL** (*Nosema bombycis* ricin-B lectin, PMID 27649890,
  [DOI](https://doi.org/10.1093/abbs/gmw093)). Adding recombinant NbRBL *enhanced* spore adhesion to
  BmN cells rather than blocking it, and the paper's own localisation calls NbRBL intracellular. The
  wider NbRBL family (52 genes) is catalogued in the open-access PMID 35736035
  ([DOI](https://doi.org/10.3390/jof8060551)), where the one functionally tested member, NbRBL28, is
  a secreted nuclear effector that drives host cell-cycle genes — a secreted-effector function, not
  adhesion. Indirect tier at best.

## Fails the inclusion rule (evidence insufficient)

- **NbPTP1 / E. hellem PTP1 O-mannosylation** (PMIDs 15501763,
  [DOI](https://doi.org/10.1128/IAI.72.11.6341-6350.2004); 32194029,
  [DOI](https://doi.org/10.1016/j.jip.2020.107350)). Mannose pre-treatment of host cells reduced
  infection and ConA binds PTP1, but the binding measured is a plant lectin binding the parasite
  protein, not the parasite protein binding a host ligand. The Sf9 heterologous-expression paper
  measures parasite proliferation, not adhesion, so it is not `het_expr` for adhesion either.
- **Encephalitozoon intestinalis host integrins** (PMID 28903465,
  [DOI](https://doi.org/10.1093/femsle/fnx169)). Soluble α3β1 and α5β1 integrin reduced spore
  adherence, implicating ~100 RGD-containing *E. intestinalis* proteins, but no individual fungal
  protein was identified or tested.
- **Bombyx BmTLP / Bmtutl-519, shrimp LvTLP** (PMIDs 23542093, 30528488, 37506856). These are the
  *host* receptors; recorded in the SWP26 and EhSWP1 notes, not fungal proteins.
- **E. intestinalis SWP1 and SWP2** (PMID 11598081,
  [DOI](https://doi.org/10.1128/IAI.69.11.7057-7066.2001)). Developmental expression and exospore
  localisation only; no function assay in either direction, so not usable as hard negatives under
  `function_shown` either.
- **Nosema bombycis SWP30 / SWP32 / SWP25** (PMID 18563739,
  [DOI](https://doi.org/10.1002/pmic.200700584)). Proteomic identification and immunolocalisation
  only. Accessions are now known from PMID 32132983 Table 1 (EF683101, EF683103, EF683102) should
  anyone want them later, but there is still no functional assay.
- **Rhizopus FTR1 as an adhesin**: rejected as an adhesin (virulence attenuation with no adhesion
  assay); kept in the table as a hard negative on its demonstrated transport function.
- **CotH expression / miRNA studies** (PMID 40488471, [DOI](https://doi.org/10.1128/spectrum.02852-24))
  and reviews (PMIDs 37360854, 35342843, 38445820, 42198648, 24355916, 29522765, 32132983):
  expression or narrative only; used to find primary papers and accessions, never as evidence of
  function.

## Already present, deliberately not re-added

| Protein | Accession | Where |
|---|---|---|
| Rhizopus delemar CotH2 | I1C4E4 | `existing.tsv` row 33 (adhesin, E1) |
| Rhizopus delemar CotH3 | I1CFE1 | `existing.tsv` row 34 (adhesin, E1) |
| Rhizopus delemar CotH7 | I1BJN3 | `existing.tsv` row 35 (adhesin, E1) |
| Nosema bombycis SWP7 | B3STP1 | `existing.tsv` row 32 (adhesin, E1) |

The CotH–GRP78 work (PMID 24355926, [DOI](https://doi.org/10.1172/JCI71349); PMID 31206021,
[DOI](https://doi.org/10.1126/sciadv.aaw1327)) supports those three rows but yielded nothing new.
SWP7 is the partner of the added SWP9 and is the same experimental unit.

## Taxa searched that returned nothing usable

- **Zoopagomycota** (Entomophthorales, *Conidiobolus*, *Basidiobolus*, *Entomophthora*): no
  protein-level adhesion evidence. Conidial adhesion to insect cuticle is described only at
  whole-cell level.
- **Blastocladiomycota** (*Allomyces*, *Blastocladiella*, *Catenaria*): nothing; no named zoospore
  adhesion protein in PubMed.
- **Chytridiomycota** beyond Bd CBM18-LL1: *Spizellomyces* and *B. salamandrivorans* return no
  protein-level adhesion or binding study.
- **Cryptomycota** (*Rozella* and relatives): nothing; no experimental protein work exists.
- **Mucoromycota beyond the CotH family**: *Mucor circinelloides*, *Lichtheimia*, *Cunninghamella*
  and *Apophysomyces* searches return only genome surveys (PMID 28923009, PMID 27447865) that
  annotate CotH copies with no adhesion assay on a named protein.
- **Nosema/Vairimorpha ceranae**: only NcSWP8 (above), which is indirect.

## Known weaknesses of what was written

1. **Q8SWH2 (E. cuniculi Hsp70 ECU02_0100) should probably be dropped by the owner's new rules.**
   Its hard-negative status rests on a *null* adherence result (anti-C1 antibody did not block,
   where anti-EnP1 did) plus internal localisation; the Hsp70 chaperone "function" is inferred from
   homology, not assayed in PMID 20811483. Under "a null adhesion result alone does not make a hard
   negative" and "moonlighting cytoplasmic proteins are a separate tier", this row fails twice over.
   It was left in place only because the second-pass instruction was to keep the existing 12 rows
   unchanged. **Owner decision needed.**
2. Thirteen of the fifteen rows rest on antibody- or soluble-protein-blocking adherence assays in
   cell culture, the only evidence available for genetically intractable microsporidia. Effect sizes
   are small (10–20% reduction in adherence) and none of these organisms has a knockout.
3. Redundancy: the two EnP1 rows (Q8SWL3, A7TZU4) are orthologues from one paper and are not
   independent clusters; the same applies to SWP9 (new) versus SWP7 (already present); and NbSWP26,
   NbSWP16, EhSWP1 and AlocSWP2 all carry heparin-binding motifs, so the positive set and one of the
   negatives share the feature that most obviously separates them.
4. Clade balance after the second pass: Microsporidia 8 adhesins + 4 hard negatives,
   Mucoromycetes 0 + 2, Chytridiomycetes 0 + 1, and still nothing from Blastocladiomycota,
   Zoopagomycota or Cryptomycota. The set tests the heparin-binding-motif mechanism in one clade
   rather than the early-diverging fungi as a whole.
5. NbSWP12 (hard negative) and NbSWP9 (adhesin) both bind chitin spore coats, so the two sets are
   not separable on chitin binding. Note also that PMID 28031263
   ([DOI](https://doi.org/10.1128/IAI.00872-16)) shows NbSWP9 additionally interacts with NbPTP1 and
   NbPTP2 and tethers the polar tube, i.e. the same structural role that puts NbSWP5 and NcSWP8 on
   the negative/indirect side. NbSWP9 is a positive only because of the separate adherence-blocking
   result in PMID 25605761.
6. NbSWP5 appears under two accessions, both 186 aa and both legitimate: UniProt G0Z414 (EMBL
   HQ881497, the SWP5 paper's own deposit, used in the table) and EF683105/ABV48893 from the earlier
   Wu et al. HSWP series, which is the one PMID 32132983 Table 1 cites. Deduplicate on sequence, not
   accession.
