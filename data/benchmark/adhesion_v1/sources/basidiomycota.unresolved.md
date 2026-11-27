# Basidiomycota curation — unresolved and rejected candidates

Curated 2026-11-26. Scope: Basidiomycota, all subphyla. Every PMID below was confirmed with the
PubMed tools and every accession with `rest.uniprot.org` or NCBI eutils; nothing here is from memory.
Candidates are listed with the reason they are **not** in `basidiomycota.tsv`.

Headline finding: the experimentally supported Basidiomycota adhesin literature is very small. After
removing what is already in the benchmark, only **two new adhesin rows** survive the inclusion rule,
both from *Schizophyllum commune*. Several famous candidates fail, and the reasons are informative.

---

## Already present — not re-added

| Protein | Organism | Where it already is |
|---|---|---|
| Cfl1 (cell flocculin 1) | *Cryptococcus neoformans* | `existing.tsv` (row `CFL1 … adhesin E1`) and `counts.md` Cryptococcus group |
| Cda3 / MP84 (mannoprotein 84) | *Cryptococcus neoformans* | `counts.md` Cryptococcus group (`2d` adhesin). UniProt P82476, 412 aa, verified |
| Cpl1 | *Cryptococcus neoformans* | `counts.md` Cryptococcus group, as a **negative** (`repeat_non_adhesin`) |
| Rep1 (repellent 1) | *Mycosarcoma maydis* / *Ustilago maydis* | `existing.tsv`, Q99109, 652 aa, verified. PMID 17159213 shows Δrep1 colonies are fully wettable and lose attachment of hyphae to Teflon — good supporting evidence for the existing row |
| Lep1 (late effector protein 1) | *Mycosarcoma maydis* | `counts.md` Mycosarcoma group. PMID 33843063 (Fukada 2021) verified: lep1 mutants show attenuated hyphal aggregation in tumours, constitutive expression induces cell aggregation |

---

## Failed the inclusion rule (adhesion never measured for the protein)

### Ustilago maydis Rsp3 — repetitive secreted protein 3
- PMID 29703884 verified (Ma et al. 2018, Nat Commun, DOI 10.1038/s41467-018-04149-0).
- What is shown: Rsp3 decorates the hyphal surface and binds two secreted maize DUF26 proteins
  (AFP1, AFP2); AFP1 is antifungal against the rsp3 mutant but not against an Rsp3-expressing strain;
  silencing AFP1/AFP2 partly rescues virulence of rsp3 mutants.
- Why rejected: this is **shielding of the hyphal surface from an antifungal lectin**, not adhesion.
  No adhesion, attachment, aggregation or binding-to-substrate assay was done. Protein–protein
  interaction with a host defence protein is not the `binding` evidence code in the sense used here
  (binding of an adhesin to a ligand/cell that it attaches to). A strong repeat-containing surface
  protein that a repeat-avidity predictor will call positive — recommend adding it as a *labelled
  uncertain* row later rather than as either positive or negative now.

### Ustilago maydis Rsp1 / Hum3 pathogenicity pair
- PMID 17917743 verified (Müller 2008). Δhum3 Δrsp1 loses pathogenicity after penetration, but
  yeast growth, mating, aerial hyphae and surface hydrophobicity are unaffected; **no adhesion
  assay**. Hum3 is in the table as a hard negative on the strength of PMID 17159213 instead.

### Cryptococcus neoformans Cda2 / "Mp98 adhesin"
- The "Mp98 is an adhesin" claim traces to T-cell-stimulation and mannoprotein-antigen work, not to
  an adhesion assay. Verified enzymatic evidence instead (PMIDs 30459196, 34712865), so Cda2 is in
  the table as a **hard negative**. Flagged because the user's candidate list named it as an adhesin.

### Cryptococcus neoformans Blp1 / Gat204 (capsule-independent antiphagocytic proteins)
- PMID 21402362 verified (Chun, Brown, Madhani 2011, DOI 10.1016/j.chom.2011.02.003).
- Shows capsule-independent **inhibition of phagocytosis**, scored as uptake by macrophages. There
  is no adhesion/attachment assay and the phenotype is anti-uptake, the opposite of adhesion.
  No "Bar1" gene exists in *C. neoformans* under that name in the verified literature — the user's
  "Bar1/Blp1" candidate resolves to BLP1 only.

### Cryptococcus neoformans Dha1, Dha2, Cfl105 (Cfl1/SIGC family)
- PMID 28039134 verified (Gyawali 2017). Shows cell-type-specific secretion and copper regulation
  of the Cfl1 homologues. No adhesion assay for any of them. Only Cfl1 (already present) has one.

### Cryptococcus neoformans 14-3-3 (CNAG_05235)
- Described as an adhesin in the Camacho & Casadevall review (PMID 30060601, verified). The primary
  work used a copper-repressible knockdown that also reduced vesicle secretion, capsule size, growth,
  morphology and cell division. Pleiotropic and cytoplasmic/vesicular rather than a wall adhesin;
  the adhesion change is not attributable to the protein acting as an adhesin. Rejected.

### Pisolithus tinctorius SRAP32 (symbiosis-regulated acidic polypeptides)
- PMID 10517026 verified. Up-regulated in ectomycorrhiza, localised to the hyphal surface by
  immunoelectron microscopy, contains an RGD motif. **No adhesion assay at all** — localisation plus
  an RGD motif is exactly the kind of inference the inclusion rule excludes.

### Puccinia graminis f. sp. tritici RGD-binding spore effectors (fibronectin-III/BRCA1-domain protein; VPS9/CUE protein)
- PMID 21873196 verified (Nirmala 2011, PNAS). RGD peptides block formation of the urediniospore
  adhesion structures, and the two proteins were purified on RGD-affinity resin — but what was
  measured for the proteins themselves is RPG1 phosphorylation/degradation and hypersensitive
  response, i.e. avirulence effector activity. Neither protein was tested in an adhesion assay, and
  no accession is given for either. Rejected on both counts.

### Trichosporon asahii cell-surface molecules (CSM) binding vitronectin / alpha-1-antitrypsin / albumin
- PMID 37258132 verified. Surface plasmon resonance and pull-downs were done on an undefined
  cell-surface **extract**; no individual fungal protein was identified. PMID 41627572 (verified,
  Mycopathologia 2026) is in-silico plus qRT-PCR expression of Cfl1/Cpl1/Dha1/Dha2 homologues — both
  fail the rule (no protein identity + no functional assay / expression only).

---

## RESOLVED — enzymes with an indirect adhesion phenotype (owner decision 2026-11-26: hard_negative)

These three were originally held back because they have quantitative adhesion phenotypes on
deletion yet no adhesive domain. **Owner decision of 2026-11-26: label them `hard_negative`**, on
the grounds that a surface/wall protein whose deletion reduces adhesion as a downstream consequence
of its enzymatic or biosynthetic function is precisely the trap a predictor must not fall into.
All three are now **rows 10-12 of `basidiomycota.tsv`** and are no longer unresolved. Each row's
`notes` field carries the sentence *"adhesion phenotype on deletion but indirect mechanism; owner
decision 2026-11-26 to label hard_negative"*, so they can never be mistaken for untested negatives.

| Protein | Accession (verified) | Length | Primary PMIDs (verified) | Mechanism recorded |
|---|---|---|---|---|
| Plb1 phospholipase B, *C. neoformans* H99 | `Q9P8P2` (PLB1, CNAG_06085, taxid 235443) | 637 | 16487740 | `phospholipase_B` |
| Mpr1 M36 fungalysin, *C. neoformans* H99 | `J9VXZ9` (CNAG_04735, taxid 235443) | 831 | 24895304 | `m36_protease` |
| Cps1 hyaluronic acid synthase, *C. neoformans* | `Q6YNJ4` (CPS1, EMBL AY063511, taxid 5207) | 449 | 16790766; 17545316 | `hyaluronic_acid_synthase` |

Verification details and caveats:

- **Plb1** — Ganendren 2006 (PMID 16487740, DOI 10.1016/j.micinf.2005.10.018) works in the H99
  background: Δplb1 adheres ~69% less to A549 lung epithelium; anti-PLB1 antibodies and
  PLB-selective (but not LPL- or LPTA-selective) inhibitors cut wild-type adhesion by 33–58%;
  **0.1–2 mM exogenous palmitic acid restores adhesion of the live Δplb1 mutant dose-dependently
  (to 67%)**, which is the clean demonstration that the fatty acid, not the protein, does the
  sticking. The same paper evidences all three enzyme activities.
- **Mpr1** — Vu 2014 (PMID 24895304, DOI 10.1128/mBio.01101-14) is in the H99 serotype A
  background. Proteolytic activity was measured directly on the recombinant protein (fluorescent
  protease assay, comparable to a trypsin control). Two caveats recorded in the row:
  (i) **this is the most contested of the three** — expression of CnMPR1 alone lets *S. cerevisiae*
  attach to and cross the brain endothelium, which in isolation would read as `het_expr`; it is
  labelled negative on mechanism, not on weak evidence, and should be excluded from any sensitivity
  estimate that treats negatives as untested; (ii) **the accession assignment is an inference** —
  the paper names no locus tag, and CNAG_04735 is the *only* peptidase M36 (Pfam PF02128) protein
  in the H99 genome (checked across all *Cryptococcus* entries: exactly one per genome).
- **Cps1** — synthase function shown directly by heterologous expression (CPS1 in *S. cerevisiae*
  and *E. coli* makes hyaluronic acid, PMID 17545316) plus hyaluronan ELISA and site-directed
  mutagenesis of the conserved glycosyltransferase domains (PMID 16790766). **Strain**: the
  deletion was made in a **serotype D** strain; `Q6YNJ4` is the cloned CPS1 product (449 aa). The
  H99 serotype A orthologue is `J9VUH0` / CNAG_04320 (451 aa) if a single-strain reference proteome
  is required. Cps1 is an integral membrane glycosyltransferase, so it is also a negative for the
  surface-localisation feature.

---

## No accession available

### Uromyces viciae-fabae spore-surface cutinase and two serine esterases
- PMID 12297669 verified (Deising, Nicholson, Haug, Howard, Mendgen 1992, Plant Cell,
  DOI 10.1105/tpc.4.9.1101). This is the **strongest functional adhesion evidence in the rust
  literature**: living spores form adhesion pads that stick to the broad-bean cuticle; autoclaved
  spores form pads that do not stick; the serine-esterase inhibitor DFP blocks adhesion of living
  spores; and **adding back the purified cutinase plus the two non-specific esterases to autoclaved
  spores restores adhesion to the host cuticle** (a reconstitution experiment).
- Why unresolved: the 1992 paper identifies the proteins only biochemically. No sequence was ever
  deposited — UniProt has **no** cutinase and in fact no entries at all for *Uromyces viciae-fabae*
  (checked by taxonomy query). Without an accession there is no row. If the owner wants rust
  coverage, the follow-up is to find a modern *U. viciae-fabae* / *U. appendiculatus* genome and
  identify the cutinase orthologue, which would be a prediction, not curation.

---

## Clades searched that returned nothing usable

Searches were run with the PubMed tools (several phrasings each; PubMed ANDs every term, so queries
were kept to 2–4 terms) plus one extended web search.

- **Malassezia** (Ustilaginomycotina): no protein with any adhesion assay. The *M. sympodialis*
  genome paper (PMID 23341551, verified) explicitly notes Malassezia has **few GPI cell-wall
  proteins and lacks the wall proteins that mediate adhesion and biofilm formation in other fungi** —
  a genuinely useful negative result for the benchmark. PMID 34379190 (MGL_3153 GPI protein) is
  homology modelling plus qPCR only. Remaining hits are host-cell immunology, not fungal adhesins.
- **Pucciniomycotina rust fungi** (*Puccinia*, *Uromyces*, *Phakopsora*): only the two entries above.
  No rust protein qualifies.
- **Rhodotorula, Sporobolomyces**: nothing, confirming the earlier finding recorded in
  `data/controls/repeat-adhesins-other-species/counts.md`. No defined protein adhesin exists.
- **Cryptococcus gattii**: no adhesin-level work distinct from *C. neoformans*.
- **Coprinopsis cinerea, Agaricus, Pleurotus, other Agaricomycetes**: hydrophobins are characterised
  for wettability/self-assembly, not for attachment, except in *Schizophyllum*. The *Coprinopsis*
  lectins (CCL2, CGL2) bind glycans directly but are cytoplasmic nematotoxins, not adhesins.
- **Rhizoctonia solani** (infection cushions): no protein-level attachment study.
- **Ectomycorrhizal Agaricomycetes** (*Pisolithus*, *Laccaria*): only SRAP32 above.

---

## Notes on the hard negatives that were included

- *U. maydis* **Hum2** and **Hum3** are the best negatives found anywhere in this clade: PMID
  17159213 tested attachment to Teflon in the hydrophobin deletion strains and found it unchanged,
  while Δrep1 lost it. Hum3 in particular (828 aa, 17 KEX2-separated Rep1-like repeats) is the same
  architecture as the Rep1 positive, so it is a direct test of the repeat-avidity call.
- *S. commune* **SC4** and *C. neoformans* **Cda1/Cda2** sit in the same protein families as
  positives already in the benchmark (SC3; Cda3/MP84). These within-family positive/negative pairs
  are intentional and are flagged in the `notes` column, but they will show up as homology conflicts
  in any clustering QC, exactly like the four conflicts already flagged in `counts.md`.
- **Cig1** and **Pit2** are included because their non-adhesion functions were shown directly
  (heme-binding spectrum of recombinant protein; protease inhibition by protein and by a synthetic
  motif peptide), not merely annotated.
