# Protein language and structure models for mapping fungal adhesins across clades, with a focus on *Rhodotorula*

Date: 2026-11-25. Literature and tool survey; nothing was run for this report. Written by Claude Code (claude-opus-5-5) from a PubMed and web
search done by a sub-agent, then spot-checked (section 6 says which items were checked and which were not).

## Why this survey

The sorting hat was built and tuned on ascomycete adhesins. In the AGA1 acceptance run of 2026-10-26
(`docs/reports/data/aga1_sorting_hat_acceptance_2026-10-26/`), the two basidiomycetes got almost no calls:

| Species | Proteins | Broad (`ser_thr_gpi_wall_candidate`) | Narrow (`flocculin_like_candidate`) | Model (`model_surface_candidate`) |
|---|---|---|---|---|
| *R. toruloides* NP11 | 8,160 | 2 | 1 | 1 |
| *C. neoformans* H99 | 7,427 | 2 | 0 | 1 |
| The six yeasts (range) | 5,070 to 6,067 | 10 to 28 | 2 to 9 | 12 to 45 |

That might be biology (basidiomycete yeasts may not use Ser/Thr-stalk GPI adhesins much) or a blind spot of the rules. The question here is which
newer models could find adhesion and biofilm proteins in *Rhodotorula* (Pucciniomycotina, Microbotryomycetes), far from the
Saccharomycotina where nearly all adhesin truth data come from. The second question is a general classifier strategy for fungal adhesion and
biofilm proteins.

The 2026-09-27 review (`docs/model-review/2026-09-27-review-and-framework-plan.md`) found the ESM-2 classifier mostly detects composition and
tandem repeats, and that a bigger model is not the fix. #15 keeps an ESM-based classifier behind a gate of 300 independent positive clusters,
and #121 repeats that rule. So the main use of new models here is **finding remote homologs and carrying labels across clades**, not training a
better end-to-end classifier.

## 1. Tools, ranked for this use

| Rank | Tool | What it adds | Cost | Licence |
|---|---|---|---|---|
| 1 | **ProstT5 + Foldseek** | Predicts the Foldseek 3Di structure alphabet from sequence, so proteomes without structures can be searched at structure-level sensitivity. Built into Foldseek (`createdb --prostt5-model`). Heinzinger et al. 2024, NAR Genom Bioinform, doi:10.1093/nargab/lqae150. | About 3B parameters; one GPU with 16 GB or more (CPU works but is slow). | MIT |
| 2 | **Reseek** | Structure search more sensitive than Foldseek and DALI, with calibrated E-values (useful when a homology cut-off is needed). Edgar 2024, Bioinformatics, doi:10.1093/bioinformatics/btae687; https://drive5.com/reseek. | CPU only; needs 3D models (AFDB or ESMFold). | GPL-3.0 |
| 3 | **TEA (The Embedded Alphabet)** | A contrastively trained 20-letter alphabet made from pLM embeddings, searchable with MMseqs2; reported close to structure-search sensitivity with no structures needed. Pantolini et al., bioRxiv doi:10.1101/2025.11.27.690975 (ICLR 2026); https://huggingface.co/PickyBinders/tea. | One GPU embedding pass, then CPU search. | Not checked |
| 4 | **ESM C 300M / 600M** | Successor to ESM-2 for embeddings; 300M is reported to roughly match ESM-2 650M. The weights are already in the HPCC Hugging Face cache. | One GPU of 16 GB or more. Proteins over about 2,000 aa need windows. | Reported as MIT in the blog and as "Cambrian Open" (300M) or non-commercial (600M) on the model cards. Fine for academic use. |
| 5 | **ProtT5-XL-U50, Ankh3** | Relative position embeddings, so no hard length limit; memory still grows with the square of length. Matters for adhesins of 1,000 to 5,000 aa. ProtT5 is the encoder behind PLM-BLAST, EBA and DeepLoc 2. Ankh3: arXiv 2505.20052. ProtT5 is in the HPCC cache. | One GPU; estimated 40 GB in half precision for about 5,000 aa (not measured). | ProtT5 academic-free; Ankh3 not checked |
| 6 | Profluent E1 | Retrieval-augmented encoder (embeds a protein with its homologs as context); could help nearly orphan *Rhodotorula* proteins. https://www.profluent.bio/showcase/e1. Not benchmarked on adhesins. | 150M to 600M | Permissive click-through |
| 7 | PLM-BLAST / EBA | Residue-level alignment of embeddings: finds a binding domain inside a long protein, which mean-pooling cannot. | CPU after embedding | — |
| 8 | SaProt 650M | Structure-aware tokens; a second view when a model exists. Same 1,022-residue limit as ESM-2. | One GPU | MIT |

Localization and PTM predictors: keep NetGPI/PredGPI and SignalP 6 (already modules). DeepLoc 2.1 (membrane type) and AIUPred (already used for
disorder) are cross-checks. NetOGlyc 4.0 is trained on human mucin-type sites and fits fungal O-mannosylation poorly. FaaPred (2010) and
FungalRV (2011) are composition-based SVMs and probably have the same saturation problem as our first model. **No fungal adhesin predictor
built on a protein language model (2024 to 2026) was found.** The nearest template is Fungtion (ESM-1b SVM ensemble for fungal effectors,
https://step3.erc.monash.edu/Fungtion).

Not recommended now: ESM3-open (generative, non-commercial), DPLM/DPLM-2, Dayhoff, xTrimoPGLM and AMPLIFY (no evidence of an advantage for this
task; AMPLIFY is a cheap Apache-2.0 baseline), gLM2 (prokaryotic genomic context), Evo 2 (DNA, needs H100s, no protein-function record), TM-Vec
and DeepFRI (superseded for this use; GO terms too coarse).

## 2. Strategies

### A. Domain-anchored cross-clade search (first to try)

1. Cut known adhesins into the folded binding domain and the stalk or repeat region (repeat module plus AFDB/ESMFold pLDDT, keeping segments
   with pLDDT 70 or more as the domain; `analysis/class2b_structure/02_confidence_profile.py` already extracts such cores).
2. Search only the binding domains against *Rhodotorula* and other basidiomycete proteomes with ProstT5 + Foldseek, Reseek and TEA. A hit
   counts when at least 2 of the 3 agree.
3. Require surface context from the existing modules (signal peptide plus GPI or wall anchor, or disorder plus a Ser/Thr stalk).

Precedent: in *Pichia kudriavzevii*, structure modelling found 11 flocculin/adhesin-like proteins with no sequence similarity to *Candida*
adhesins (Alvarado et al. 2023, PLoS Pathog, PMID 37196016).

**Test without leakage.** Queries come from Saccharomycotina only. Score recovery of basidiomycete adhesins that never enter the query or
training sets (*C. neoformans* Cfl1/Dha1/Dha2/Cpl1/Cfl105; *Trichosporon asahii* homologs; *Ustilago maydis* Rep1/Hum2/Hum3). Decoys are
same-length surface glycoproteins (SED1, TIR4, MSB2, HKR1-like) and a shuffled-repeat control. Report precision at the scale of published
inventories (about 10 to 30 calls per genome), not ROC-AUC.

### B. Embedding nearest-neighbour labels by mechanism

Instead of one "adhesin" label, use mechanism classes: repeat stalk (avidity), folded binding domain (affinity), hydrophobin or repellent,
matrix, and wall-remodelling enzyme (GH16, GH72, GH76; Crh, Gas, Dfg5). Embed with ESM C 300M and ProtT5 using windows (for example 1,000 aa
with 250 overlap), keep the maximum over windows plus a separate N-terminal 500-aa vector, and mask padding (#10). Assign by nearest labelled
neighbours and abstain when the nearest neighbour is far.

**Test.** Leave one clade out (Saccharomycotina, Pezizomycotina, Basidiomycota) with grouping by orthogroup or family. Add a composition-only
baseline (amino-acid frequencies and length) that the embeddings must beat on the held-out clade; otherwise the embedding is relearning Ser/Thr
composition.

### C. A pan-*Rhodotorula* orthogroup atlas (no labels needed)

Build orthogroups (OrthoFinder, or MMseqs2 on TEA strings) across available *Rhodotorula* genomes (*R. toruloides* NP11, *R. mucilaginosa*,
*R. graminis* WP1, *R. diobovata*) and other Microbotryomycetes. Keep groups that are secreted or anchored, disordered or repeat-rich, and
lineage-specific, especially with copy-number expansions near chromosome ends (the *C. auris* Hil expansion is the precedent: Smoak et al. 2023,
Genetics, PMID 36794645). Rank them with A and B; check against biofilm versus planktonic expression where data exist.

### D. LoRA or contrastive fine-tuning of ESM C on binding-domain crops only

Training on binding domains (not whole proteins) stops the model from shortcutting on repeats or composition; hard negatives are surface
glycoproteins from the same organisms. Test: family hold-out plus clade hold-out as in B, and check that predictions barely change when the stalk
is deleted or shuffled. **This stays behind the #15 gate** (about 100 non-redundant adhesins exist now; the gate is 300 clusters).

## 3. *Rhodotorula* and basidiomycete literature

**No *Rhodotorula* adhesin gene has been characterized** (PubMed search, 2026-11-25). The nearest work is phenotypic:

- Giraldo et al. 2026, PNAS: *R. taiwanensis* sticks to insect cuticle through hydrophobic interactions; no genes. PMID 42296353.
- Yockey et al. 2019, mSphere: *R. mucilaginosa* envelope and capsule; chitin more exposed than in *C. neoformans*. PMID 31167944.
- Gattlen et al. 2011, Biofouling: *R. mucilaginosa* biofilm formation, repeatability and cell attachment. PMID 21981116.
- Antunes et al. 2024, Biotechnol Biofuels Bioprod: in an evolved *R. toruloides* strain, cell-surface biogenesis and remodelling genes change
  the most. PMID 38807231.

Genomes: *R. toruloides* (Zhu et al. 2012, PMID 23047670; chromosome-level CBS 14, PMID 34648882); *R. mucilaginosa* JGTA-S1 (PMID 30615101)
and clinical DH4148 (PMID 42788761); *R. graminis* WP1 (PMID 26441909); *R. diobovata* (PMID 33924147).

Basidiomycete adhesion and surface proteins (truth-set candidates; see #13):

- *C. neoformans* **Cfl1**, the first adhesin described in Basidiomycota (Wang, Zhai and Lin 2012, PLoS Pathog, PMID 22737071); shed into the
  matrix as a paracrine signal (PMID 23798436).
- The **Cfl1/Dha1/Dha2/Cpl1/Cfl105** family shares a C-terminal domain called SIGC and is basidiomycete-specific (PMID 28039134).
  Pfam **PF21671 "Protein CPL1-like"** (InterPro IPR048661) covers *C. neoformans* CPL1; InterPro lists 72 Pucciniomycotina proteins in 12
  proteomes with it (checked 2026-11-25). Whether PF21671 is the same region as the SIGC domain of PMID 28039134 was not checked.
- *Trichosporon asahii* homologs of the four are up-regulated in biofilm (Dos Santos et al. 2026, Mycopathologia, PMID 41627572).
- *U. maydis* **Rep1 repellents** (repeat peptides cut by Kex2) took over attachment to hydrophobic surfaces; Hum2 is a class I hydrophobin;
  Hum3 is a hydrophobin preceded by 17 Kex2-separated repeats (PMID 17159213). This "repeat precursor" class may be missed by the current modules.
- *Malassezia*: the *M. furfur* secreted protease Mfsap1 controls dispersal and adhesion (PMID 36442106); *M. sympodialis* genome PMID 23341551.
- *C. neoformans* vesicle proteome 2026: cell-wall-integrity mutants lose GPI-anchored mannoproteins from vesicles (PMID 42925063).

Families to search for first in *Rhodotorula* proteomes: CPL1-like (PF21671), hydrophobins, and Rep1-like Kex2 repeat precursors.

## 4. Proposed next step

Strategy A on the HPCC as a test: Saccharomycotina binding-domain queries against *R. toruloides* NP11 and *C. neoformans* H99, with the
*C. neoformans* Cfl1 family as the held-out positive check. The ProstT5 step runs on a GPU (`exfab`); the searches run on CPU. Plan and scripts in
a follow-up change.

## 5. Compute already on the HPCC

- Hugging Face cache (`~/.cache/huggingface/hub`): ESM C 300M and 600M, ESM-2 150M and 650M, ESMFold, ProtT5-XL half-precision encoder,
  Evo 2 1B.
- Foldseek: `module load foldseek/10-941cd33-gpu` (and a CPU build); static binary at `_workdir/bin/foldseek`; PDB database in
  `_workdir/foldseek_db/`. Foldseek needs AVX2 nodes (`epyc`, not `c01`).
- GPU partitions: `exfab` (one node, 2 Ada 6000, 30 days), `short_gpu` (2 h), `gpu`.

## 6. What was checked and what was not

Checked on 2026-11-25 (Crossref or InterPro API): the DOIs of ProstT5 (Heinzinger 2024), Reseek (Edgar 2024) and TEA (Pantolini, bioRxiv);
PF21671 name and taxonomic spread. The survey draft gave PF06493 for the SIGC domain; **that is wrong** (PF06493 is DUF1096).

All 20 PubMed IDs in this report were checked against NCBI esummary (first author, journal, year, title); three details of the draft were corrected
(first author of PMID 21981116, the species of Mfsap1, the year of PMID 41627572). Only titles were checked, not the papers' content.

Not checked: licences of Ankh3 and
TEA; ESM C's true training length; memory use of ProtT5 or Ankh3 at about 5,000 aa; AFDB coverage of *R. mucilaginosa* and *R. graminis*
(NP11, UP000016926, has models); citations for SaProt, PLM-BLAST, EBA, DeepLoc 2.1 and Evo 2; the claim that Ankh3 beats ESM-2 (the authors'
own result). bioRxiv was not searched directly, so 2025 to 2026 *Rhodotorula* preprints may be missing.

## Sources

https://huggingface.co/EvolutionaryScale/esmc-600m-2024-12 ; https://evolutionaryscale.ai/blog/esm-cambrian ; https://arxiv.org/abs/2505.20052 ;
https://www.biorxiv.org/content/10.1101/2025.11.27.690975 ; https://www.profluent.bio/showcase/e1 ; https://drive5.com/reseek ;
https://www.ebi.ac.uk/interpro/entry/pfam/PF21671/ ; https://step3.erc.monash.edu/Fungtion ; PubMed records by PMID above.
