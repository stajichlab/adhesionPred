# Filamentous Ascomycota (non-yeast, non-A. fumigatus) — unresolved and rejected candidates

Date: 2026-11-26. Every PMID below was confirmed with the PubMed tools; accessions were checked
against `https://rest.uniprot.org/uniprotkb/<ACC>.json` or NCBI eutils esummary. Items are here
because they fail the inclusion rule, because no sequence could be pinned down, or because they are
a judgement call the owner should settle.

## A. Already present — not re-added

Checked against `existing.tsv` and `data/controls/repeat-adhesins-other-species/counts.md`:

- **Pyricularia oryzae MPG1** (P52751) and **MHP1** (G4MWK2) — already present (E2/E3).
- **Metarhizium anisopliae MAD1** (Q2LC49) and **MAD2** (Q2LC47) — already present (E1). The
  primary paper (PMID 17337634, "The MAD1 adhesin of Metarhizium anisopliae links adhesion with
  blastospore production and virulence to insects, and the MAD2 adhesin enables attachment to
  plants") was re-found and matches.
- **Coccidioides SOWgp** (Q8NK60/Q8NK61/Q96V71) — already present.
- **Paracoccidioides gp43** (Q9HDL9) — already present.
- **Pneumocystis jirovecii Msg** (O74670) and **PCINT1** — already present.
- **Blastomyces dermatitidis BAD1/WI-1** (A4D962) and **Histoplasma Yps3**, **Histoplasma Hsp60**
  (P50142) — already present.
- **Botrytis cinerea BcLysM1** — already present. The primary paper is PMID 39655398 (LysM effector,
  mycelial adhesion to bean leaf surfaces plus chitin-triggered immunity suppression); it matches.
- **Fusarium Hyd3**, **Talaromyces GAPDH**, **Penicillium HfbA** (negative), **Botrytis Bhp1/2/3**
  (negatives), **Blastomyces Eng2** (negative), **Pyricularia MoMsb2** (negative, mucin sensor),
  **Blumeria Lip1** (negative), **Zymoseptoria ZtGT2** (negative), **Coccidioides CTS1** (negative),
  **Histoplasma Cbp1** (negative), **Talaromyces HSP60** (negative) — already present.

## B. Verified but deliberately excluded (judgement calls)

### B1. Developmental / signalling mutants with an adhesion side-effect — NOT adhesins

- **Botrytis cinerea BcSUN1** (PMID 28163701) — **RESOLVED 2026-11-26, now row 24 of the table as a
  hard_negative.** Owner decision: a wall/surface protein whose deletion reduces adhesion but which
  has no adhesive domain, the loss being downstream of its structural or enzymatic role, is exactly
  the trap a predictor must avoid, so it belongs in the negatives rather than being dropped.
  Accession resolved from the paper text (Bcin06g06040, EnsemblFungi B05.10) to UniProt
  **A0A384JKR5** (B. cinerea B05.10, 471 aa, consistent with the reported 48 kDa glycoprotein),
  fetched live. The second SUN gene in the genome (Group-II, Bcin07g06600 = A0A384JNI4, 530 aa) was
  never knocked out and is not entered.
- **Beauveria bassiana Bbhog1** (PMID 19363067). Δhog1 has reduced attachment to insect cuticle but
  it is a HOG MAP kinase; the attachment defect is explained by down-regulation of hyd1/hyd2.
  Signalling mutant, excluded.
- **Beauveria bassiana Ssr4** (PMID 32317391). Chromatin-remodelling SWI/SNF co-subunit; conidial
  adhesion loss is a transcriptional consequence. Excluded.
- **Pyricularia oryzae Tra1p** (PMID 23733043) — transcription factor, spore adhesion phenotype.
  Excluded (not a surface protein).
- **Pyricularia oryzae MoPmt2** (PMID 27199956) — protein O-mannosyltransferase; "adhesive
  properties" defect is a glycosylation effect. Excluded.
- **Pyricularia oryzae EMP1** (PMID 15055545) — extracellular matrix protein required for
  appressorium formation and pathogenicity; the phenotype reported is appressorium development, not
  a quantitative adhesion assay. Excluded for lack of an adhesion measurement.
- **Botrytis cinerea Bcchs3a** (PMID 20672878) — chitin synthase; "altered adhesion" is downstream of
  a wall-synthesis defect and the protein is a membrane-embedded enzyme, not a wall adhesin.
- **Fusarium oxysporum FoFLP1/3/4/5** (PMID 35487139) — fasciclin-like proteins knocked down by
  RNAi; only conidiation and disease severity were measured. No adhesion assay → excluded.
- **Pyricularia oryzae MoFLP1** — the homologous case WAS included (ko_assay, conidial adhesion
  measured), but flagged in the table as a weak positive because the GFP fusion localised to the
  vacuole membrane and the mutant is pleiotropic.

### B2. "Adhesin" genes with no adhesion phenotype — cannot be hard negatives either

- **Beauveria bassiana Adh1/Mad1** (UniProt J4US32, BBA_02419, 619 aa) and **Adh3**
  (J4WB32, BBA_03909, 242 aa), PMID 34169325. Both disruptants were indistinguishable from wild type
  in every phenotype examined, including conidial adherence. They are named adhesins with a *null*
  adhesion result, which is informative but is not "a demonstrated non-adhesion function", so under
  the strict rule they qualify neither as adhesins nor as hard negatives.
  **Owner decision 2026-11-26: these stay OUT of the table — a null adhesion result is not evidence
  of a non-adhesion function.** Recorded here only, with that reasoning.
- **Beauveria bassiana hyd1A** in the 2011 study (PMID 21375591) likewise showed *no* adhesion
  effect; the 2024 study (PMID 39724799) reports the opposite. Entered as an adhesin on the strength
  of the 2024 result, with the conflict recorded in the notes column.
- **Clonostachys rosea Hyd1/Hyd3** (PMID 24483277, Hypocreales, in scope). ΔHyd3 had reduced
  Arabidopsis root colonization but ΔHyd1 did not, and the double mutant colonized *better* than
  wild type. Colonization, not an adhesion assay, and the direction is incoherent → excluded.

## C. Candidates that could not be pinned to a sequence

- **Colletotrichum** — no entry. The only protein-level adhesion evidence is the 110 kDa spore-coat
  glycoprotein of *C. graminicola* (mAb UB20 blocking, microsphere binding), which has no resolvable
  sequence; it is already in the project's `rejected.tsv`. PubMed searches for Colletotrichum
  adhesion genes with deletion mutants returned nothing usable. **Colletotrichum remains a gap.**
- **Magnaporthe spore tip mucilage** — repeatedly reported as the adhesive of the conidial apex, but
  no protein component has been identified and knocked out. PMID 18051333 (ECM and integrin-like
  proteins in conidial adhesion) and PMID 21732962 (enzymatic detachment of germlings by
  collagenase/gelatinase/MMPs) both demonstrate an adhesive ECM but name no gene or accession.
  PMID 27315248 attributes germling adhesion to chitosan, a polysaccharide, not a protein.
- **Blumeria graminis** — no adhesin found. No effector or surface protein of *Blumeria* has a
  quantitative adhesion assay with a defined gene; the taxon is not transformable, so only the
  existing Lip1 hard negative stands.
- **Trichoderma harzianum ThSWO and ThFFSL**, PMID 42664349 (Sci Adv 2026) — **RESOLVED
  2026-11-26, now rows 25 and 26 of the table.** The PMC full text (PMC13524027), "Gene and protein
  sequence information" section, states the accessions outright: ThSWO = NCBI **XP_024771705.1**
  (verified by eutils esummary: 499 aa, hypothetical protein M431DRAFT_484109, *T. harzianum*
  CBS 226.95, taxid 983964 — consistent with the paper's stated 52.13 kDa theoretical mass) and
  ThFFSL = UniProt **A0A2T3ZSN5** (verified live: 301 aa, "Fucose-specific lectin", same strain and
  taxid — consistent with the stated ~32 kDa). The experiments were done in strain NJAU4742
  (genome LVVK00000000.1), so the accessions are the reference-strain sequences the authors
  themselves used, not NJAU4742 deposits; this is recorded in the notes column. Both ΔThswo and
  Δ*Thffsl* knockouts show significantly reduced attached hyphal biomass on rice roots (Fig. 6D, 6G),
  ΔThswo is rescued by exogenous ThSWO protein (Fig. 6E), and ThSWO-ThFFSL binding was measured by
  ITC. Also noted while reading the full text: the six co-acting O-glycosidases are
  OPB36282.1, OPB37664.1, OPB38657.1, OPB41138.1 (this one is in fact cited as ThFFSL in the
  proteomics figure), OPB41830.1, OPB43763.1, and the proteomics are in PRIDE PXD068647. The
  OPB41138.1 / A0A2T3ZSN5 double citation of ThFFSL in the same paper is an inconsistency in the
  source; the table uses the accession the Methods section designates as the ThFFSL sequence.
- **Trichoderma/other hyphosphere saSSCPs (HFS family)**, PMID 33538393. Recombinant hydrophobins,
  cerato-platanins and the new HFS family condition glass and PET surfaces to a moderately
  hydrophilic state that "significantly favoured subsequent hyphal attachment". This is substrate
  conditioning by a soluble protein, not protein-mediated adhesion of the cell, and the effect on
  root attachment was mixed → excluded, but it is the reason Sm1/Epl1 is entered as a hard negative
  rather than ignored.
- **Fusarium oxysporum Msb2** (PMID 21441438), the Sordariomycete counterpart of the already-present
  MoMsb2 hard negative. Excluded only because no UniProt entry exists for the *F. oxysporum* f. sp.
  *lycopersici* 4287 Msb2 (searches of taxid 426428 returned nothing; the available Msb2 entries are
  from f. sp. *cubense* and f. sp. *raphani*). Easy to add if a strain-matched accession is
  acceptable.
- **Pneumocystis carinii MSG / gpA** — host-side binding partners are well documented (surfactant
  protein A and D carbohydrate recognition domains bind MSG, PMIDs 9201957 and 11306442), and
  vitronectin/fibronectin bind P. carinii beta-glucan rather than a protein (PMID 11509330). These
  show host lectins binding the fungal surface, not a fungal protein binding a host ligand, and the
  P. jirovecii Msg is already present, so no new row was added. Searching found no Pneumocystis
  adhesin with a knockout (the organism is not genetically tractable); PMID 26215665 is a Flo8-like
  *transcription factor*, excluded as a signalling/regulatory protein.
- **Talaromyces marneffei Mp1p** — PubMed query failed to return a clean result set in this session;
  not assessed. Worth a second pass.
- **Metarhizium Mcl1** (collagen-like protein). Searches returned only a promoter-engineering paper
  (PMID 30800586) and reviews; the original Wang & St Leger work describes Mcl1 as an *anti*-adhesive
  / immune-evasion coat on hemolymph blastospores, which would make it a hard negative, but the
  primary paper was not located and verified here. Not entered.

## D. Dothideomycetes — per-genus second pass (2026-11-26)

The first attempt used one combined query and was swamped. This is the redone pass, one genus at a
time. Net result: **one new row (Bipolaris maydis Hyp1, hard_negative) and no adhesin positives.**
The clade is not empty of adhesion biology, but it is almost empty of *named proteins with
sequences* — in three of the six genera the adhesive was characterised biochemically and never
identified. That is the clade-level finding, not an artefact of the search.

| genus | queries run (hit counts) | verdict |
|---|---|---|
| *Zymoseptoria* | `Zymoseptoria AND (adhesion OR adherence OR attachment)` (9); `Zymoseptoria AND (adhesin OR hydrophobin OR "cell wall protein")` (0) | nothing new; see Z1 |
| *Alternaria* (incl. *alternata*, *brassicicola*) | `Alternaria adhesion` (35, all titles screened); `Alternaria brassicicola AND (adhesion OR adherence OR attachment OR adhesin)` (7); `Alternaria spore adhesion appressorium attachment` (0) | no adhesin, no usable negative; see A1, A2 |
| *Cochliobolus* / *Bipolaris* | `Cochliobolus[TIAB] AND (adhesion OR attachment OR adhesin OR hydrophobin)` (5); `Bipolaris[TIAB] AND (adhesion OR attachment OR adhesin)` (6) | **1 hard negative added**; see C1 |
| *Parastagonospora* / *Stagonospora nodorum* | `(Parastagonospora OR "Stagonospora nodorum" OR "Phaeosphaeria nodorum") AND (adhesion OR attachment OR adhesin OR hydrophobin)` (3) | real adhesion biology, no sequence; see P1 |
| *Venturia inaequalis* | `"Venturia inaequalis" AND (adhesion OR attachment OR adhesin OR hydrophobin OR subcuticular)` (13, all screened); `Venturia inaequalis appressorium protein` (0) | real adhesion biology, no sequence; see V1 |
| *Pyrenophora* | `Pyrenophora[TIAB] AND (adhesion OR attachment OR adhesin OR hydrophobin)` (4) | genuinely nothing; see Y1 |

Note on PubMed behaviour for anyone repeating this: the bare term "Cochliobolus" is auto-mapped to
the Ascomycota MeSH heading and returns thousands of irrelevant records, so `[Title/Abstract]` has to
be forced. Likewise `Bipolaris` alone pulls in human keratitis/peritonitis case reports.

**C1. *Bipolaris maydis* (*Cochliobolus heterostrophus*) Hyp1–Hyp4**, PMID 41356074 (2025).
ADDED: Hyp1 as a hard negative (M2UQJ0). Single, triple and quadruple hydrophobin disruptants showed
no change in colony hydrophobicity, growth, conidiation, stress tolerance, virulence or sexual
reproduction, yet each gene restored hydrophobicity when expressed in an *A. nidulans* ΔdewA mutant.
The paper gives no accessions; it states the sequences came from the JGI *C. heterostrophus* C5 v2.0
gene catalog and reports Hyp1 = 117 aa, Hyp2 = 100 aa, Hyp3 = 157 aa, Hyp4 = 182 aa. Those four
lengths are each unique among the twelve C5 hydrophobin entries in UniProt, giving an unambiguous
mapping: **Hyp1 = M2UQJ0 (117 aa), Hyp2 = M2UUF3 (100 aa), Hyp3 = M2URC1 (157 aa),
Hyp4 = M2UIG9 (182 aa)**, all four fetched and verified live. Only Hyp1 (class I, the DewA/RodA
counterpart) was entered; the three class II paralogues are available on identical evidence if the
owner wants them. Caveat carried in the notes column: adhesion was never assayed in *B. maydis*, and
the non-adhesion function was demonstrated heterologously rather than in the native host.

**Z1. *Zymoseptoria tritici* ZtGT2 — back-annotation patch (2026-11-26).**
This pass found the primary paper (PMID 29020037), which supports the existing hard-negative label
for a stronger reason than was previously recorded: the authors state explicitly that **"adhesion to
leaf surfaces was unaffected in ΔZtGT2 mutants"**, while the type-2 glycosyltransferase is required
for hyphal extension on solid surfaces and for virulence. That turns an assumed negative into a
tested one — the same quality of negative as *M. oryzae* Cut2. No other *Zymoseptoria* surface
protein has an adhesion assay, and the genus has no adhesin positive.

*Where the row actually is.* **ZtGT2 is not in `filamentous_asco.tsv`** (that file has no
*Zymoseptoria* row; it was excluded from this batch as already present). The row lives in the
project's own curation at
`/Users/jstajich/projects/adhesionPred/data/controls/repeat-adhesins-other-species/hard_negatives_other_species.tsv`,
**line 11**, and that table uses a different schema (`accession, gene, protein_name, organism,
taxon_id, stratum, evidence_level, pmid, quote, serves, tuned, leakage, curator, date, notes`) with
a `quote` field rather than `assay_summary`. The project tree is read-only for this work, so the
patch is recorded here rather than applied.

*Accession verified live:* **F9WWD1** — *Zymoseptoria tritici* (strain CBS 115943 / IPO323),
taxid 336722, "Type 2 glycosyltransferase", gene GT2, 471 aa; the UniProt entry itself cites
PMID 29020037. This matches the existing row exactly, so no accession change is needed.

*Patch 1 — `pmid` field: no change required.* The row **already carries 29020037**; the PMID was
never missing. Nothing to add.

*Patch 2 — `quote` field (this is the real gap).* The stored quote is about virulence, not adhesion:

- OLD: `"The inability to extend surface hyphae resulted in a complete loss of virulence on plants."`
- NEW: `"Despite this, adhesion to leaf surfaces was unaffected in deltaZtGT2 mutants and global RNAseq-based gene expression profiling highlighted that surface-sensing and protein secretion was also largely unaffected."`

*Patch 3 — `notes` field.* Append: `Adhesion to leaf surfaces was measured and was UNAFFECTED in the deletion mutants (PMID 29020037), so this is a tested negative, not an assumed one; the demonstrated function is synthesis of an outer cell wall polysaccharide required for hyphal extension on solid surfaces. Back-annotated 2026-11-26.`

*Equivalent `assay_summary` in the `filamentous_asco.tsv` schema*, should the row ever be migrated
into this table (clade_class `Dothideomycetes`, label `hard_negative`, mechanism
`type 2 glycosyltransferase synthesising an outer cell wall polysaccharide required for hyphal
extension on solid surfaces`, evidence_code `function_shown`, pmids `29020037`, uniprot_length `471`):

> A T-DNA mutant unable to extend hyphae specifically on solid surfaces mapped to the predicted type 2 glycosyltransferase ZtGT2; deletion mutants lose surface hyphal extension and plant virulence, and electron microscopy shows the protein maintains the outermost cell wall surface, but the authors measured adhesion to leaf surfaces and found it unaffected in the deletion mutants, with surface-sensing and protein secretion also largely unchanged by RNAseq.

*Not added:* the three *Bipolaris maydis* class II hydrophobin paralogues (Hyp2 = M2UUF3,
Hyp3 = M2URC1, Hyp4 = M2UIG9). Owner decision 2026-11-26: one row from that evidence is enough, as
three near-identical paralogues from a single experiment would skew the Dothideomycete fold toward
one paper. Their accessions are recorded in C1 above if that is ever revisited.

**A1. *Alternaria* — verified absence of any adhesin.** All 35 hits for `Alternaria adhesion` were
screened by title: they are allergy/ILC2 immunology (integrins, ICAM-1, VCAM-1 — *host* adhesion
molecules, not fungal), nanomaterial and pesticide-formulation papers, keratitis case reports and
biocontrol screens. Not one concerns an *Alternaria* adhesion protein. *A. alternata* and
*A. brassicicola* searched separately return nothing either. **Alternaria is genuinely empty at the
protein level.**

**A2. *Alternaria brassicicola* spore-surface lipase**, PMID 10556710 (verified). An 80 kDa serine
lipase on the surface of ungerminated spores; anti-lipase antibodies reduced blackspot lesions by
90% on intact cauliflower leaves but not on dewaxed leaves, and the authors propose it "interact[s]
closely with epicuticular leaf waxes for adhesion and/or penetration". This is the best
Dothideomycete hard-negative candidate on function (a demonstrated enzyme at the cell surface), but
it fails verification: the protein was never sequenced (it was detected with a heterologous
anti-*Botrytis* lipase antibody), there is no accession, and its adhesion role is explicitly
speculative. Not entered.
Also screened and rejected: *A. brassicicola* **AbMak1** (PMID 28955470), a flavin-dependent
monooxygenase whose deletion alters conidial wall morphology and "hydrophilic adhesion" — an
intracellular enzyme acting through melanisation, so neither a surface protein nor an adhesin.

**P1. *Parastagonospora (Stagonospora) nodorum***, PMID 17998157 (verified). Conidia adhere to
polystyrene within minutes; SDS/chitinase/lyticase release, fluorescent-lectin labelling, and
protease and lectin treatments in a quantitative adhesion assay together show that **pre-formed
spore-surface glycoproteins mediate the rapid adhesion**. The evidence standard for adhesion is met,
but no individual protein was purified, named or sequenced — the paper works with SDS-PAGE bands.
Nothing can be entered. **Follow-up if the owner wants this genus: the *S. nodorum* genome is public
and the Green lab's later work may have identified the bands.**

**V1. *Venturia inaequalis* spore tip glue (STG)**, PMID 18943251 (verified). Non-germinated conidia
adhere to hydrophobic surfaces within minutes via a localised apical "spore tip glue" droplet;
histochemistry shows protein plus carbohydrate, lectin labelling shows β-galactose and
N-acetylglucosaminyl residues, and protein-synthesis inhibitors do not block adhesion, so the
adhesive is pre-formed. Again a demonstrated protein-containing adhesive with no identified protein.
Also screened: the *V. inaequalis* secretome (PMID 36678413), effector repertoire (PMID 36329441),
cell wall carbohydrate dynamics (PMID 37039647) and the melanin-deficient isolate (PMID 36675856) —
none assays adhesion for a named protein, and the appressorium work (PMID 18943680) implicates
melanin, not a protein.

**Y1. *Pyrenophora* — verified absence.** Only four hits across *P. teres* and *P. tritici-repentis*:
a bacterial-ghost pesticide-delivery paper, two on the Ptr ToxA host-selective toxin (a necrosis
inducer taken up by host cells, not an adhesin), and a 1983 study of barley mesophyll cells binding
to *P. teres* hyphae that names no fungal molecule. No adhesion protein and no surface-enzyme
negative.

## E. Other taxa searched with no qualifying result

Pezizomycetes and Eurotiomycetes other than *Aspergillus nidulans* and the dimorphic human pathogens
returned no protein with a direct adhesion assay. *Ophiostoma* sensu stricto returned nothing, but
its order-mate *Sporothrix* turned out to be the richest source in this work (six rows).
