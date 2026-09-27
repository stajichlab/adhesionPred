# Basidiomycete yeast adhesion and biofilm: what is actually known

Literature survey, 2026-09-27. Companion to `data/curated/biofilm/basidiomycete_seeds.tsv`.
Written against the coverage gap flagged in `data/curated/biofilm/README.md`: 179 of 207
biofilm rows are *Candida albicans*, and the only curated basidiomycete adhesin was Cfl1,
whose accession could not be resolved.

Everything below traces to a PMID or a database record. Where I could not resolve something,
it says so instead of guessing.

---

## 1. Headline results

**1. The Cfl1/CPL1 family accessions are resolved.** All five *C. neoformans* H99 members,
mapped through the FungiDB community gene-name table
(`fungidb/cryptococcus` -> `Community_annotations/H99_names.csv`) and then to UniProt:

| gene | CNAG locus | UniProt | length | source of the locus |
|---|---|---|---|---|
| CFL1 | CNAG_00795 | **J9VHS9** | 309 | FungiDB name table, cites PMID 22737071 |
| CPL1 | CNAG_02797 | **J9VQE5** (reviewed) | 199 | FungiDB name table |
| DHA1 | CNAG_07422 | **J9VSS0** | 327 | FungiDB name table |
| DHA2 | CNAG_06082 | **J9VY16** | 317 | PMID 28039134 only |
| CFL105 | CNAG_03454 | **J9VVN7** | 409 | PMID 28039134 only |

These are exactly the five Pfam **PF21671** ("Protein CPL1-like", InterPro IPR048661,
PANTHER PTHR35192 PriA/Cpl1_fungi) proteins in the H99 proteome, which is an independent
consistency check on the mapping. DHA2 and CFL105 are *not* in the FungiDB name table; their
loci rest on a single paper (PMID 28039134) and a WebFetch of its PMC full text, so treat
them as one notch less certain than the other three.

**2. Basidiomycete adhesin architecture is not ascomycete adhesin architecture.**
Pfam counts from `rest.uniprot.org` (2026-09-27), Basidiomycota (taxon 5204) vs
Ascomycota (taxon 4890):

| Pfam | domain | Basidiomycota | Ascomycota |
|---|---|---|---|
| PF00624 | Flocculin repeat | **0** | 331 |
| PF10182 | Flocculin_t3 | 5 | 215 |
| PF05792 | Candida_ALS | 3 | 492 |
| PF11766 | Candida_ALS_N | 6 | 313 |
| PF10528 | Flo11 | 8 | 1674 |
| PF07691 | PA14 | 675 | 3335 |
| **PF21671** | **CPL1-like / SIGC** | **1640** | **0** |
| PF01185 | Hydrophobin_1 | 4340 | 1061 |

The training positives for the shipped model are FLO/ALS/FLO11 (see
`docs/model-review/2026-09-27-review-and-framework-plan.md` §5.2). Those domains are
effectively absent from the Basidiomycota. The one basidiomycete family with a demonstrated
adhesin in it is absent from the Ascomycota. There is no domain overlap to transfer on.

**3. The composition is wrong too, not just the domain.** Computed from the UniProt
sequences of the five H99 family members and two *C. albicans* Als proteins:

| protein | length | % Ser+Thr | % Cys |
|---|---|---|---|
| Cfl1 (J9VHS9) | 309 | 20.7 | **8.7** |
| Cpl1 (J9VQE5) | 199 | 14.6 | 5.0 |
| Dha1 (J9VSS0) | 327 | 21.7 | 6.7 |
| Dha2 (J9VY16) | 317 | 19.6 | 6.6 |
| Cfl105 (J9VVN7) | 409 | 21.3 | 10.0 |
| Als1 (Q5A8T4) | **1260** | **36.8** | 0.8 |
| Als3 (Q59L12) | 1155 | 34.5 | 1.0 |

Cfl1-family proteins are ~4x shorter, far less Ser/Thr-rich, and roughly **8x more
cysteine-rich** than Als adhesins. They are secreted and shed (PMID 23798436), not GPI-anchored.
Every feature the review doc identifies as the model's decision boundary — "secreted +
low-complexity Ser/Thr-rich + long" (§5.1) — points the *wrong way* for this family. The
model will not just be imprecise on basidiomycetes; it should be expected to have near-zero
recall on the one basidiomycete adhesin family we can name.

**4. Family membership does not imply adhesion, within this very family.** Cfl1 is an
adhesin (PMID 22737071). Its paralog **Cpl1 is a secreted immune effector** that drives M2
macrophage polarization through host TLR4 and is essential for virulence, with no reported
adhesion or biofilm role (PMID 35896747, *Nature*). Same domain, same secretion, opposite
label. This is the basidiomycete version of the YWP1 problem and is why every PF21671 row
in the seeds table is E2/E3, never E1.

---

## 2. Experimentally demonstrated vs. inferred vs. hypothesis

**Experimentally demonstrated (E1):**
- *C. neoformans* **Cfl1** — adhesion, biofilm, morphogenesis, paracrine signalling; the
  first adhesin described in the Basidiomycota (PMID 22737071, 23798436, 24567775).
- *C. neoformans* **Znf2** — regulator; necessary and sufficient for filamentation, acts on
  pathogenicity partly through cell adhesion; Cfl1 is its downstream adhesion target
  (PMID 22737071). Note the sign flip worth recording: ZNF2 overexpression *promotes*
  adhesion and filamentation and *abolishes* virulence.
- *C. neoformans* **Crz1** — deletion affects biofilm formation (PMID 23640031). Direction
  not stated in the abstract, so the table says `unknown`.
- *U. maydis* **Rep1** — deletion abolishes attachment to hydrophobic Teflon (PMID 17159213).
  Already in `data/curated/adhesins/literature_seeds.tsv`.
- *C. neoformans* **Cpl1**, **Blp1** — demonstrated surface/secreted function that is
  explicitly *not* adhesion (PMID 35896747, 21402362). Useful stage-2 negatives.

**Inferred from homology or from a weaker phenotype (E2):**
- **Cfl105**, **Dha1**, **Dha2** — secretion confirmed; overexpression adherence effects
  ranged from "strong" (Cfl105) to "modest" (Dha1/Dha2); no deletion adhesion phenotype
  (PMID 28039134).
- *U. maydis* **Hum3**, **Rsp1** — the double mutant is non-pathogenic, but surface
  hydrophobicity, mating and attachment were *unaffected*; the defect is host penetration
  (PMID 17917743). Calling these adhesins would be over-reading.

**Hypothesis only (E3):**
- All four *T. asahii* flocculins. The evidence is qRT-PCR (CFL1 and DHA1 strongly
  biofilm-induced) plus in-silico docking to albumin/collagen IV/fibronectin/hemoglobin/
  laminin (PMID 41627572). No deletion, no binding assay, and **no accession published**.
- Every `adhesin_candidate_family` row (Rhodotorula, Naganishia, Papiliotrema, Filobasidium,
  Cutaneotrichosporon, Apiotrichum, Kwoniella, Sporidiobolus). These are Pfam-membership
  pools, nothing more.

---

## 3. Which species have real molecular work, and which have none

| lineage | molecular genetics on adhesion/biofilm? | what exists |
|---|---|---|
| *Cryptococcus neoformans* | **Yes** | Deletion/overexpression genetics, secretion assays, ChIP, GPI-protein mutant libraries. The only basidiomycete yeast where this is true. |
| *Ustilago maydis* | **Yes**, but plant-pathogen, not a yeast biofilm | rep1/hum2/hum3/rsp1 knockouts (PMID 17159213, 17917743) |
| *Trichosporon asahii* | Expression only | qRT-PCR + in silico (PMID 41627572); a PDT transcriptome flags A1Q1_04029, A1Q1_01345, A1Q1_08069, A1Q1_01456 as possible morphogenesis/biofilm genes (PMID 40382755) — unnamed, undirected, not included in the table |
| *Malassezia* spp. | **None at gene level** | Phenotype only: hydrophobicity/adherence/biofilm on plastics and catheters (PMID 28340187, 17510859, 22682201); keratinocyte adhesion mediated by *host* glycosaminoglycans (PMID 34052141) |
| *Rhodotorula* spp. | **None** | PubMed "Rhodotorula adhesin gene" = **0 hits**. Phenotype only (PMID 29058139, 42296353, 23827647, 25843277) |
| *Sporobolomyces*/*Sporidiobolus* | **None** | PubMed "Sporobolomyces adhesion surface" = **0 hits** |
| *Naganishia*, *Papiliotrema*, *Filobasidium* | **None** | Genome-derived PF21671 pools only |

So the honest position is: **one species** (*C. neoformans*) carries essentially all the
basidiomycete evidence, and one protein in it (Cfl1) carries the only E1 adhesin label.

---

## 4. Capsule, polysaccharides, and the trap they set

This needs to be explicit because it is the single easiest way to mislabel *Cryptococcus*.

- **GXM (glucuronoxylomannan) is a polysaccharide, not a protein.** Cryptococcal biofilm
  formation *depends on the capsule* and correlates with the capsular polysaccharide's
  ability to bind the solid support (PMID 16177306, 17513597). EDTA blocks biofilm by
  reducing GXM shedding and vesicle secretion (PMID 22941091).
- **Direction is not obvious even for antibody.** Anti-GXM antibody in solution *prevents*
  biofilm (PMID 16177306); the same specificity *coated on the surface* immobilizes cells and
  *accelerates* biofilm (PMID 19251903). A naive "involved in biofilm" label would merge these.
- **Cps1**: cps1 loss and hyaluronidase both reduce binding to human brain endothelial cells,
  and the host receptor is CD44 (PMID 17545316, 18248627). But Cps1 is a membrane
  glycosyltransferase — the adhesive molecule is its hyaluronic acid product. And the
  *Neurospora* ortholog makes a polysaccharide that is **not** hyaluronic acid (PMID 24953997),
  so even the chemistry is not settled.
- **Pbx1/Pbx2**: deletion gives *clumpy cells* and dry colonies, but the demonstrated role is
  fidelity of GXM synthesis (PMID 17337638). **A cell-aggregation phenotype in Cryptococcus is
  not evidence of an adhesin.** These are in the table as hard negatives for exactly this reason.
- The FungiDB H99 name table also carries CAP10/CAP59/CAP60/CAP64/CAS3/CAS32-35/CAS41, all
  capsule-biosynthesis proteins. None are adhesins; several will look like surface
  glycoproteins. They are a ready-made hard-negative set and were not enumerated in the TSV
  only to keep it reviewable.

---

## 5. Does the GPI cell-wall-protein story even exist in basidiomycetes?

Partly, and unevenly.

- *C. neoformans* does have GPI-anchored, Ser/Thr-rich mannoproteins: **MP88** (CNAG_00776)
  has a C-terminal Ser/Thr-rich O-glycosylated region plus a predicted GPI site, and the
  original paper notes at least 11 further H99 genes share the architecture (PMID 12228274).
  A systematic screen of 47 predicted GPI-anchored proteins found phenotypes for MP88 and the
  unnamed CNAG_00137 (reduced phagocytosis, stress sensitivity) — but **no adhesion assay was
  run in that study** (PMID 35640861). So there is a cryptococcal stage-1 surface class; there
  is no cryptococcal stage-2 evidence inside it.
- *Malassezia* appears to have lost much of it. The genome paper states Malassezia "do not
  appear to have many cell wall-localized GPI proteins and lack other cell wall proteins
  previously identified in other fungi" (PMID 23341551), and UniProt returns **zero** PF21671
  proteins in the genus. Malassezia adheres and forms biofilm anyway. This is a real
  biological absence, and it means a protein-centric adhesin predictor may have nothing to
  find in the most abundant fungus on human skin.
- *Ustilago*/Ustilaginomycotina use a third solution entirely: **repellents** (Pfam PF29797),
  which UniProt shows are confined to Ustilaginomycotina (Ustilago, Sporisorium, Tilletia,
  Moesziomyces, Pseudozyma, Testicularia, Meira...) and are absent from Agaricomycotina and
  Pucciniomycotina. Conversely PF21671 is absent from Ustilaginomycotina. The two candidate
  basidiomycete adhesion systems are **phylogenetically complementary**.

Caveat on a number I did *not* put in the table: UniProt `keyword:KW-0336` (GPI-anchor) counts
are 8 for *C. neoformans* H99 vs 96 for *C. albicans* SC5314. That ratio mostly reflects
curation depth, not biology, and should not be quoted as evidence. If we want a real
comparison we have to run NetGPI/PredGPI over both proteomes ourselves.

---

## 6. What this means for the model

Mapping onto the tiers in `docs/model-review/2026-09-27-review-and-framework-plan.md` §6:

- **T3 (out-of-clade) is currently untestable for Basidiomycota**, because we have exactly one
  E1 positive (Cfl1) plus one plant-pathogen positive (Rep1). Two positives is not a clade test.
- The clean, cheap first experiment: **score Cfl1 (J9VHS9), Cfl105, Dha1, Dha2 and Cpl1 with
  the shipped model and with the stage-1/stage-2 models.** My prediction from §1.3 above is
  that all five score low, and that the model's calls on *C. neoformans* instead land on the
  Ser/Thr-rich GPI mannoproteins (MP88 and relatives) and the capsule machinery. If that is
  what happens, it is a concrete, publishable failure mode and it is already set up.
- Cpl1 is a near-ideal **paralog-level hard negative**: same family, same secretion, no
  adhesion. Very few hard negatives are this well matched.
- The `adhesin_candidate_family` rows should **not** be used as training positives. They are a
  pool to prioritize, and the Kwoniella expansion (278 PF21671 proteins, the largest in
  UniProt, in the sister genus to Cryptococcus) is the most interesting signal in that pool.

---

## 7. What to do next, in order

1. **Verify the five accessions by sequence, not by ID.** Pull J9VHS9 and confirm against the
   Cfl1 sequence in PMID 22737071 / 23798436 supplementary data. The CNAG->UniProt mapping is
   consistent from two directions (FungiDB name table, PF21671 membership) but has not been
   checked at the residue level, and DHA2/CFL105 rest on one paper.
2. **Score the family with the current model** (item above). Cheap, decisive.
3. **Resolve the *T. asahii* orthologs.** Blast/MMseqs the five H99 PF21671 proteins against
   the *T. asahii* var. *asahii* CBS 8904 proteome (the 8 PF21671 members are A1Q2_07842,
   A1Q2_03191, A1Q2_03930, A1Q2_00039, A1Q2_05189, A1Q2_00505, A1Q2_00709, A1Q2_01196) and
   fill in the blank `uniprot_query` cells. Do not assign them by gene name alone — the
   Trichosporon paper assigned names purely by homology and published no accession.
4. **Build a PF21671 HMM profile and a family alignment** across the genera in the table.
   This gives a basidiomycete family-held-out fold for T2 that does not exist today.
5. **Enumerate the Cryptococcus hard negatives properly**: capsule biosynthesis (CAP*/CAS*),
   GPI mannoproteins without adhesion evidence (MP88, MP98/Cda2, Cig1, CNAG_00137), Pbx1/Pbx2,
   Blp1, Cpl1. That is ~20 well-attested stage-1-positive / stage-2-negative proteins from one
   genome, which is more hard negatives than the whole current basidiomycete positive set.
6. **Do not extend the biofilm table from GO for these species.** The direction problem is
   worse here than in Candida: the same antibody prevents or promotes cryptococcal biofilm
   depending on whether it is free or surface-bound, and the dominant matrix component is a sugar.
7. **Accept, and write down, that Malassezia and Rhodotorula may not be solvable by a protein
   adhesin model.** For Rhodotorula there is not a single published adhesion gene. For
   Malassezia the genome argues the GPI cell-wall protein repertoire is largely absent.
   Kingdom-wide coverage claims should carry that caveat.

---

## 8. Known weaknesses of this survey

- **One person, one day, PubMed + UniProt.** No manual reading of full texts beyond abstracts
  and two PMC fetches (PMC5311391, PMC12864254). The CNAG assignments for DHA2 and CFL105 come
  from a WebFetch summary of the PMC full text and should be re-read directly.
- PubMed ANDs and auto-expands every term, so short queries were used throughout; that finds
  the obvious papers and will have missed things phrased unusually. Notably, `Cryptococcus
  Cfl1 adhesin` returns only **4** records in all of PubMed, which is itself the finding.
- Non-yeast Agaricomycetes (mushroom-forming fungi) were not surveyed even though they hold a
  large share of the PF21671 proteins (Mycena 72, Lentinula 41, Armillaria 14...). If the
  family turns out to be an adhesin family, that is where most of it lives.
- No structure work was done. AlphaFold models of the SIGC domain would settle quickly whether
  it is plausibly a ligand-binding module or purely a structural/cross-linking one.
- `direction` is `unknown` for more rows than I would like, and I left it unknown rather than
  inferring it — consistent with the warning in `data/curated/README.md`.
