# Literature import log

This project has no Mycelium living-repo (`.living/`) structure, so imports are tracked here
instead: one entry per paper dropped into `to_import/`, recording the citation, what was
extracted, and where it landed. Source files stay in `to_import/` for provenance; they are not
deleted after processing.

**Copyright note: this is a public GitHub repo (`stajichlab/adhesionPred`).** Publisher PDFs
(Elsevier, ASM, JSTOR/OUP, etc.) are NOT committed — `to_import/*.pdf` is gitignored. Only the
extracted facts/sequences/citations go into git, which is fair use; the PDFs themselves stay
local-only in this working directory. An openly-licensed source (CC-BY, like the MDPI paper
below) is fine to commit in full and is the one exception already in the repo.

---

## molecules-23-03145.xml

**Citation:** Duarte Escalante E, Frías De León MG, Martínez García LG, Herrera J, Acosta
Altamirano G, Cabello C, Palma G, Reyes Montes MdR. "Selection of Specific Peptides for
*Coccidioides* spp. Obtained from Antigenic Fractions through SDS-PAGE and Western Blot Methods
by the Recognition of Sera from Patients with Coccidioidomycosis." *Molecules* 2018, 23(12),
3145. DOI: [10.3390/molecules23123145](https://doi.org/10.3390/molecules23123145). PMID:
30513599, PMC6321320.

**Read:** 2026-09-28. Full JATS XML (`to_import/molecules-23-03145.xml`), Table 3
("Identification of proteins that recognize *Coccidioides* spp.") specifically.

**What Table 3 contains:** 3 proteins identified by mass spec from SDS-PAGE/western-blot bands
(100 kDa ×2, 50 kDa ×1) that were most frequently recognized by antibodies in sera from 27
patients with confirmed coccidioidomycosis, with BLAST cross-reactivity against related fungi
(*C. immitis*/*C. posadasii*, *H. capsulatum*, *A. fumigatus*, *P. brasiliensis*, *P. lutzi*).

**Extracted into:**
- `analysis/cocci_antigens/literature_antigens.fa` — full-length sequences for the 3
  accessions (A0A0J6YAG6, A0A0J8R1I4, A0A0J6F383). All three are now **deleted from UniProtKB**
  ("not part of a reference proteome") and were recovered from **UniParc** instead
  (UPI0000D87081, UPI00065EBA54, UPI0001A7E586) — UniProt's archive keeps sequences after
  UniProtKB deletion.
- `analysis/cocci_antigens/literature_antigens_meta.tsv` — per-protein band, source species,
  and the paper's own cross-reactivity BLAST results (species, % identity), transcribed
  directly from Table 3.

**Not yet done:** these 3 proteins are not yet wired into `02_score_antigens.py` as a scored
input (distinct from the IEDB confounder panel in `00_download_iedb_antigens.py` — these are
serum-validated *Coccidioides* antigens, not cross-reactivity confounders from other fungi).

**Why this matters beyond antigen scoring:** A0A0J6YAG6 ("GPI anchored serine/threonine-rich
protein") is architecturally similar to the tandem-repeat surface proteins (class 2a) in the
queued retrain-test task (`docs/HANDOFF-HPCC.md` §4) — worth checking against
`analysis/model_review/repeat_structure_transfer.py` as a possible *Coccidioides*-native
positive example, not just an antigen candidate.

---

## 1-s2.0-S0378111996004866-main.pdf

**Citation:** Zhu Y, Yang C, Magee DM, Cox RA. "*Coccidioides immitis* Antigen 2: analysis of
gene and protein." *Gene* 1996, 181:121-125. DOI:
[10.1016/S0378-1119(96)00486-6](https://doi.org/10.1016/S0378-1119(96)00486-6).

**Read:** 2026-09-28. Full text (5 pages, text-extracted with `pypdf`; no poppler/pdftotext
available on this HPCC node, so the Read tool's PDF page-rendering is unavailable here — see
note at the bottom of this log).

**What it contains:** the original cloning of the Ag2/PRA gene (`pra`) from *C. immitis*
mycelial-phase genomic DNA by PCR: 582 bp ORF split by two introns, 194 aa translation product,
N-terminal signal peptide, C-terminal GPI-anchor signal, and an 11-repeat Pro-Thr tetrapeptide
(TXX'P) region (aa 89–141) that is the dominant predicted antigenic/surface-exposed epitope.

**Extracted into:** `analysis/cocci_antigens/literature_antigens.fa` and
`..._meta.tsv` — full 194 aa sequence. **Not reconstructed from this paper's own (OCR-noisy)
sequence figure** — instead fetched directly from UniProt accession **Q6QJA6** (`pra` gene,
"Proline-rich antigen", *C. immitis*, 194 aa), confirmed to match this paper's reported ORF
length and translation length exactly.

**Why this matters:** Ag2/PRA is already one of the acceptance-test anchor genes in
`analysis/cocci_antigens/NOTES.md` (matched by gene name against the pangenome), but the
pipeline had no primary-literature reference sequence or citation for it until now.

---

## 30099284.pdf

**Citation:** Pappagianis D, Smith CE, Kobayashi GS, Saito MT. "Studies of Antigens from Young
Mycelia of *Coccidioides immitis*." *J Infect Dis* 1961, 108(1):35-44. JSTOR:
[30099284](http://www.jstor.com/stable/30099284).

**Read:** 2026-09-28. Full text (11 pages, `pypdf`).

**What it contains:** a methods paper — autolysis of young (rather than aged) mycelial cultures
to speed up and improve the yield/potency of coccidioidin antigen preparations for the
precipitin and complement-fixation (CF) serologic tests. Foundational for the "coccidioidin"
crude-antigen preparations used clinically for decades, but predates protein purification or
sequencing of any specific antigen.

**Extracted:** nothing — no named protein, accession, or sequence to extract. Kept for
historical/methods citation only (why coccidioidin is prepared the way it is).

---

## calhoun1986.pdf

**Citation:** Calhoun DL, Osir EO, Dugger KO, Galgiani JN, Law JH. "Humoral Antibody Responses
to Specific Antigens of *Coccidioides immitis*." *J Infect Dis* 1986, 154(2):
[full citation in PDF; no DOI, pre-DOI era].

**Read:** 2026-09-28. Full text (8 pages, `pypdf`).

**What it contains:** immunoblot survey of patient and rabbit-immune sera against 3
coccidioidal extracts, identifying antigens at **100, 60, 45, and 70 kDa** recognized by
patient sera. The **100 kDa antigen co-migrates with the conventional tube-precipitin (TP)
antigen** by agar diffusion (TP antigen is the historical name for what later molecular work
identified as SOWgp). The 70 kDa antigen distinguished pulmonary from disseminated disease
(reactive in 5/6 pulmonary sera, absent in disseminated-disease sera) — a specificity signal,
not just a presence/absence one.

**Extracted:** nothing sequence-level — this is 1986 immunoblot work, pre-cloning of any of
these antigens. Kept for citation: it is the primary source tying the "100 kDa"/TP-antigen
identity to SOWgp, and for the pulmonary-vs-disseminated 70 kDa finding, which is not
represented anywhere else in this repo's antigen scoring.

---

## cole-et-al-1989-isolation-of-antigens-with-proteolytic-activity-from-coccidioides-immitis.pdf

**Citation:** Cole GT, Zhu S, Pan S, Yuan L, Kruse D, Sun SH. "Isolation of Antigens with
Proteolytic Activity from *Coccidioides immitis*." *Infect Immun* 1989, 57(5):1524-1534.

**Read:** 2026-09-28. Full text (11 pages, `pypdf`).

**What it contains:** biochemical isolation of two proteinase antigens from the saprobic
(mycelial) phase soluble conidial wall fraction: **AgII**, a 60 kDa serine proteinase (weak/no
reactivity with 21 patient sera by ELISA — the authors suggest it is poorly presented to the
host during infection), and **AgCS**, separated into 39 and 19 kDa fractions with most of the
proteolytic activity.

**Extracted:** nothing sequence-level — pre-cloning biochemistry. Kept for citation: AgII's
*negative* serology result is a useful negative control/caveat for any future proteinase-based
antigen candidate (having proteolytic activity or being immunodominant in culture does not
imply patient-antibody reactivity).

---

## Note on PDF handling in this environment

This HPCC node has neither `poppler-utils` (`pdftoppm`/`pdfinfo`) nor a PDF Python library
installed, so the Read tool's page-rendering path fails outright, and figures/images inside
PDFs (e.g. sequence-trace figures) are not recoverable as text — only what the PDF's embedded
text layer contains. `pypdf` was installed user-locally (`pip3 install --user pypdf`, resolves
under `/usr/bin/python3.12` specifically — the default HPCC `python3` is 3.9 and does not have
it) to extract that text layer directly. Where a paper's key sequence exists only inside a
figure image (e.g. the Ivey 2003 Vaccine paper below), it is not recoverable this way, and is
reported as unavailable rather than reconstructed by hand from partial OCR fragments.

---

## 1-s2.0-S0264410X03004857-main.pdf

**Citation:** Ivey FD, Magee DM, Woitaske MD, Johnston SA, Cox RA. "Identification of a
protective antigen of *Coccidioides immitis* by expression library immunization." *Vaccine*
2003, 21:4359-4367.

**Read:** 2026-09-28. Full text (9 pages, `pypdf`).

**What it contains:** expression library immunization (ELI) screen of a *C. immitis*
spherule-phase cDNA library (800–1000 genes) in BALB/c mice, narrowed by successive
fractionation to a single protective clone, **7-3-5-5 = ELI-Ag1**: 672 bp ORF, 224 aa, 19 aa
N-terminal signal peptide, 15 aa C-terminal GPI-anchor site, partial homology to a
*Neurospora crassa* hypothetical protein. Table 1 lists predicted antigenic determinant and
MHC-binding peptide fragments (partial sequence coverage, ~12–24 aa each) with position
numbers.

**Extracted:** nothing full-length. **The full ELI-Ag1 sequence is not retrievable**: (1) it
was checked against UniProt/NCBI under the names "ELI-Ag1", "ELI1", and by length — no
*Coccidioides* entry exists (only a Histoplasma ortholog, **Q1HRW5**, "ELI-ag1-like protein",
218 aa — a homolog, not this paper's sequence); (2) the paper's own Figure 5, which shows the
full nucleotide+translated sequence, is a rendered image in the PDF, not text, so it can't be
extracted without the poppler/OCR tooling this node lacks (see note above). Table 1's short
antigenic-determinant fragments were deliberately **not** stitched into a reconstructed
"full-length" sequence — partial fragments with position numbers are not the same as a verified
sequence, and presenting a hand-assembled guess as ELI-Ag1's sequence would be worse than
having none.

**If the full sequence is needed later:** options are (a) request `poppler-utils` be installed
on this node to render Fig. 5 for visual/OCR reading, (b) contact the authors or check GenBank
under author name (Cox RA / Magee DM) rather than protein name, or (c) treat the Histoplasma
ortholog Q1HRW5 as a stand-in for homology search only, never as ELI-Ag1 itself.

---

## PIIS0021925819710679.pdf

**Citation:** Sheppard DC, Yeaman MR, Welch WH, Phan QT, Fu Y, Ibrahim AS, Filler SG, Zhang M,
Waring AJ, Edwards JE Jr. "Functional and structural diversity in the Als protein family of
*Candida albicans*." *J Biol Chem* 2004, 279(29):30480-30489. DOI:
[10.1074/jbc.M401929200](https://doi.org/10.1074/jbc.M401929200). PMID: 15128742 (checked against
PubMed on 2026-10-06: the DOI and title match this PDF).

**Read:** 2026-10-06. Full text (10 pages, text-extracted with `pypdf`; no poppler on this node).
The PDF is gitignored (Elsevier). Only the extracted facts below are committed.

**Why it was imported:** PMID 15128742 was already listed in the `pmids` of the `ALS9`, `ALS3`,
`ALS6`, `ALS5` and `ALS1` rows of `data/curated/adhesins/adhesins.tsv`, but nothing had been
extracted from the paper. The repeat-mechanism curation cited a different paper (PMID 15116430, the
ALS1 paper) as the mechanism source for ALS9.

**What the paper contains:** `ALS1`, `ALS3`, `ALS5`, `ALS6`, `ALS7` and `ALS9` were cloned and
expressed in *S. cerevisiae* S150-2B (ADH1 promoter). `ALS2` and `ALS4` could not be amplified.
Surface expression was confirmed by flow cytometry. Adherence was tested on gelatin, fibronectin,
laminin, FaDu epithelial cells and endothelial cells. N-terminal domain swaps (Als5p/Als6p) put the
substrate specificity in the N-terminal domain. Homology models place the N-terminal domains in the
immunoglobulin superfamily.

**ALS9 (A0A1D8PQ86), as extracted:**
- Adhered above background to **laminin only**. Not to gelatin, fibronectin, epithelial cells or
  endothelial cells (page 4).
- Surface expression was detected: 11.4% of cells above background (4-fold) with antiserum A and
  33.9% (13-fold) with antiserum B (Table II). Antiserum A gave the lowest signal of the six.
- Not named as invasive in the text (Als3p, Als1p and Als5p are). Fig. 3 was not read.
- N-terminal model: collagen-binding protein of *S. aureus* (PDB 1d2p) as the primary homolog, and
  Als2p and Als9p share the same primary, secondary and tertiary homolog (page 6). Als2p, Als4p and
  Als9p form a third structural group (group C; page 7).
- **No statement about tandem repeats of ALS9.** The paper says only that ALS genes in general have
  tandem repeats (page 3).

**Other family members extracted:** Als1p, Als3p and Als5p bound all substrates tested; Als6p bound
gelatin only; **Als7p bound none of the substrates tested**.

**Extracted into:** `data/controls/repeat-mechanism/paper_extractions.tsv` (12 rows; page, method,
result, a short quote, what it supports, limits). No sequence was taken from the paper.

**Consequences for the curation tables (not applied to the agent tables):**
- ALS9 `evidence_level` E1 is supported: direct, heterologous, narrow (laminin).
- The ALS9 mechanism label `2a` remains `family_inference`. This paper does not give ALS9 repeats and
  puts the binding specificity in the folded N-terminal domain. The cited source 15116430 is about
  ALS1.
- The ALS7 row (Q5A312, *C. albicans*) is listed as an E1 adhesin from PMID 17510860, but this paper
  found no adherence for Als7p in five substrates. The two results conflict. A person should read both.
- ALS2 and ALS4 get no support from this paper.

**Not done:** Fig. 1 percentages and Fig. 3 (invasion) were not extracted from the figures. The
Als9p sequence variant used (plasmid from a library of unstated strain background) was not checked
against UniProt A0A1D8PQ86.
