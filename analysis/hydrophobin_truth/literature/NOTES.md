# Notes on the literature extraction (2026-10-08)

Files here: `spacing_quotes.md`, `protein_lists.tsv`, `tools/` (the three scripts I used). There is no `protein_sequences.faa`. No supplement contains an amino-acid sequence (see section 3).
Nothing was committed to git.

## 1. What could be read, and how

| Input | Result |
|---|---|
| Read tool on PDFs | Failed. `pdftoppm` is not installed. `pdfinfo`, `pdftotext`, PyPDF and PyMuPDF are also absent. |
| What worked for PDFs | Ghostscript 9.27 (`gs -sDEVICE=txtwrite`, and `-sDEVICE=png16m` to render pages that I then looked at). `tools/colx.py` rebuilds two-column text from the `txtwrite` XML output. |
| Linder 2005 (20 pages) | Read in full for sequence-related content (journal pp. 877-885). Clean text layer. Pages 886-896 were searched by keyword only (self-assembly, applications, references). Fig. 3 (alignment) is not visible in my render. I saw the caption only. |
| Sunde 2008 (12 pages) | Read the text for pp. 773-778 in full. Later pages (polymer structure, TEM, AFM) searched by keyword only. Fig. 1 (alignment) is a picture. I looked at it and did not transcribe it. |
| Wessels 1994 (25 pages) | Hydrophobin section (pp. 422-430) read in full. The rest of the review (cell wall enzymes, chitin) was not read. The PDF is a scan with an OCR text layer that has many errors. Quotes are OCR-corrected and marked so. |
| Jensen 2010 Additional file 2 (1 page) | The text layer is scrambled. I rendered the page at 200 dpi and read the table from the image. See the caution in section 3. |
| Jensen 2010 Additional file 3 (3 pages) | A rotated alignment picture (Jalview). No text layer. At 60 dpi it is not readable. No sequences were extracted. |
| Xu 2021 mmc1 to mmc4 (.doc) | Converted. antiword, catdoc, libreoffice, soffice and python olefile are not installed. I wrote a small pure-Python reader for the Word 97 format (OLE container, piece table) in `tools/doc2txt.py`. It gave clean text for all four files. Table layout is kept as tab-separated cells. |
| Wösten 1994, PMID 8005099 | PubMed record only (abstract). Not full text. Recorded in `spacing_quotes.md` section 4. |

The earlier file `docs/paper/05-literature-verification.md` section 6 lists Wessels 1994, Linder 2005 and Sunde 2008 as NOT READ. These three were read in this task. Result: they give no numeric spacing (see below).

## 2. Findings that matter for the spacing task (H3)

1. None of the three reviews (Wessels 1994, Linder 2005, Sunde 2008) gives a numeric Cys-to-Cys spacing in its text. Each gives only the qualitative pattern: eight Cys, the 2nd and 3rd adjacent, the 6th and 7th adjacent (Linder p. 878).
2. The numbers exist only as alignment pictures (Linder Fig. 3, Sunde Fig. 1, Wessels Fig. 1 is hydropathy, not sequence). If the owner wants numbers from these papers, someone must read the pictures. I did not do that.
3. Sunde 2008 p. 774 and Linder 2005 p. 883 both say that class II has more uniform spacing and loop lengths than class I. They give no numbers. I did not compare these statements with the numbers in `data/sorting_hat/cys8_spacing.yaml`.
4. Wessels 1994 defines class I and II by hydropathy pattern (aligned on the Cys), not by Cys spacing. Wessels says cerato-ulmin has homology "particularly with respect to the spacing of the eight cysteine residues". So spacing is shared across both classes at the level of the 1994 data.
5. The Jensen 2010 Additional file 2 has per-protein spacing in the form `CN{a}CCN{b}CN{c}CN{d}CCN{e}C`. That gives residue counts for the five gaps for 49 proteins (1 is a fragment with no pattern). These are numbers from the paper itself. They are in the column `cys_pattern_stated` in `protein_lists.tsv`. Gap order: C1-C2 = a, C3-C4 = b, C4-C5 = c, C5-C6 = d, C7-C8 = e. In this notation the adjacent pairs C2-C3 and C6-C7 are the "CC". (I take N{x} as x residues. The file does not define the notation. For example ATEG_04730 has 10, 11, 16, 8, 10. That is the same set of numbers that 05-literature-verification.md 6.3 lists for class II in other papers. I did not check the others.)
6. One protein in the Jensen table is class II by all three criteria: ATEG_04730 (A. terreus). Jensen's title says that Aspergillus hydrophobins "cannot be clearly divided into two classes". The table shows that: 26 of 50 are "Intermediate" in the theoretical class.

## 3. Protein lists

### Jensen 2010 (PMID 21182770), Additional file 2

- 50 rows. 9 genomes of 8 species (A. fumigatus has two strains). Per genome: A. oryzae RIB40 2; A. niger CBS 513.88 8; A. niger ATCC 1015 7; E. nidulans FGSC A4 6; A. fumigatus AF293 5; A. fumigatus A1163 4; A. terreus NIH 2624 5; A. flavus NRRL 3357 7; A. clavatus NRRL 1 6.
- The table has three class columns. Counts over 50 proteins:

| Column | I | II | Intermediate | "-" (no class) |
|---|---|---|---|---|
| Class based on cysteine pattern | 45 | 1 | 0 | 4 |
| Class based on hydropathy plot | 24 | 11 | 15 | 0 |
| Theoretical class | 23 | 1 | 26 | 0 |

- In `protein_lists.tsv`, `class_as_stated` is the "Theoretical class" column. "Intermediate" is written as `other (Intermediate)`. "-" is written as `unclassified`. The other two columns are kept as extra columns. The owner should decide which column to use (see section 5).
- `accession_if_given`: none. The table gives gene or locus IDs only (AO0900..., An..g..., JGI IDs for A. niger ATCC 1015, AN....2 for A. nidulans, AFUA_/AFUB_, ATEG_, AFLA_, ACLA_). No UniProt or GenBank accession is in the file. The ID is in the `protein_id` column.
- Evidence: the file gives none. The proteins are predicted from genome sequences. I did not read the Jensen main text in this task. Earlier work in the repo read it in full (05-literature-verification.md 6.1). I do not use that here.
- Row `JGI128530` is a fragment. The table says "Fragment (similar to An07g03340)". It has no pattern.
- Some A. niger CBS 513.88 and ATCC 1015 rows, and the A. fumigatus AF293 and A1163 rows, are strain orthologs. For example the pattern for An07g03340, JGI194815 is the same family. I did not map orthologs. If the truth set needs one entry per protein family, the owner must decide how to merge strain copies. A. niger ATCC 1015 JGI IDs are not mapped to CBS 513.88 IDs in the file.
- Transcription caution. I read the table from a 200 dpi render, because the text layer is scrambled. I did not run a second independent read. 36 of the 49 patterns also appear unbroken in the scrambled text layer and match my transcription. The other 13 patterns and 31 of the 50 gene IDs could not be matched in the scrambled text and were checked by eye only. Please spot-check before using the numbers as frozen parameters.
- No sequences. Additional file 3 is a picture of an alignment.

### Xu 2021 (PMID 33636611), mmc1 to mmc4

| File | Content | Proteins |
|---|---|---|
| mmc1.doc | Table S1. Primers for hydrophobin gene amplification (forward and reverse primer, 40 rows). | 40 gene names |
| mmc2.doc | Table S2. Primers for qRT-PCR (18 hydrophobin genes plus a Tubulin reference). | 18 of the 40 (all appear in S1). Tubulin is not a hydrophobin and is not in the list. |
| mmc3.doc | Table S3. Primers for the FBH1 RNAi vector and transformant detection. | none new |
| mmc4.doc | Figure S1 caption only: "The hydropathy patterns of P. ostreatus hydrophobins. The red pentagram represents the position of cysteine doublets." The figure itself did not come out of the .doc (no picture stream; I did not find the image). | none |

- 40 gene names: Vmh1, Vmh2, Vmh3, FBH1, Hydph5 to Hydph21 (17), POH1 to POH3 (3), Po.hyd1 to Po.hyd16 (16). 4 + 17 + 3 + 16 = 40. This agrees with "40 genes" in the title/abstract of the paper as given in the task.
- `class_as_stated`: not stated in any of the four files. All 40 are "not stated".
- `accession_if_given`: none. No UniProt, GenBank, JGI or gene ID is in the files. Names only.
- Evidence: none stated in these files. Primers exist for all 40 genes. qRT-PCR primers exist for 18. A primer list is not evidence of expression or of protein. The results (expression data, which genes were detected) are in the main paper, which I did not have.
- No sequences. Primers only. I did not use the primers as sequence evidence.
- The numbering has a gap. Names Hydph1 to Hydph4 do not appear. Vmh1-3 and FBH1 may be the same genes as Hydph1-4, but the supplement does not say so. I did not assume it.
- Possible duplicates by identical forward primer (computed, not assumed): Vmh2 and Po.hyd9 have the same forward primer `AGTCAATCTACTCGCAATGTTCTCC` (reverse primers differ). Hydph15 and Hydph16 have the same forward primer `CTCCCCTGCTCACTGAC` (reverse primers differ). This may mean the same gene under two names, or primer reuse for different genes. Not resolved. All 40 names are kept as separate rows.

## 4. Counts

| Source | Proteins | Class I | Class II | Intermediate / other | Unclassified | Not stated |
|---|---|---|---|---|---|---|
| Jensen 2010 (theoretical class) | 50 | 23 | 1 | 26 | 0 | 0 |
| Jensen 2010 (cysteine-pattern class) | 50 | 45 | 1 | 0 | 4 | 0 |
| Jensen 2010 (hydropathy class) | 50 | 24 | 11 | 15 | 0 | 0 |
| Xu 2021 | 40 | 0 | 0 | 0 | 0 | 40 |
| Total rows in `protein_lists.tsv` | 90 | | | | | |

### Experimental evidence level (what the files show)

| Source | Protein-level evidence | Gene-level evidence | Prediction only |
|---|---|---|---|
| Jensen 2010 Additional file 2 | none stated | none stated | all 50 (genome predictions, as far as the file shows) |
| Xu 2021 supplements | none stated | primers for all 40 and qRT-PCR primers for 18 (no results in the files) | all 40 by default, because the files do not state anything else |

Neither supplement documents a purified protein, a mass spectrum or an antibody. None of the 90 rows can be used as a protein-level positive from these files alone. For the truth set the owner would need to take the evidence from the main text of each paper.

Statements about protein-level evidence in the three reviews (not in the lists): Wessels 1994 p. 425: only Sc3 and Sc4 of the class I proteins in his Fig. 1 had been identified as proteins by 1994. Linder 2005 p. 880: HFBII crystal structure (T. reesei). Sunde 2008: structures of HFBI, HFBII (crystal) and EAS (NMR). These are not part of the lists above.

## 5. Decisions the owner must make

1. Which Jensen class column to use as the label: cysteine-pattern (45 I, 1 II, 4 none), hydropathy (24 I, 11 II, 15 Intermediate), or theoretical (23 I, 1 II, 26 Intermediate). The first two columns disagree for 9 proteins labelled I by pattern and II by hydropathy. Possible circularity: the Jensen cysteine-pattern class comes from Cys spacing, and the rescue rule also uses Cys spacing. The design spec has a rule against circularity (section 4.1). I did not read that section closely and did not apply it. I only flag it.
2. Whether Jensen proteins with no experimental evidence enter the truth set at all. The file states no evidence. The design spec section 4 covers the truth set. I did not check it for this question.
3. Whether to take the Cys spacing numbers from the Jensen table (a paper that gives per-protein numbers) as a source for H3. They are model-free numbers from a paper, but each class is defined by the authors. The numbers for the five gaps are not a consensus; they are 49 individual proteins.
4. The Xu 2021 main text is needed to get gene IDs, classes and expression evidence. The supplements give none. A request to the owner for the PDF of the main text, or its Table 1.
5. Whether the pictures (Sunde Fig. 1, Linder Fig. 3, Jensen Additional file 3) should be read for numbers. I judged that reading residues from pictures is too error-prone to do silently.
6. The two papers disagree on the size of the HFBII hydrophobic patch (Linder p. 880: about 12% of the surface; Sunde p. 777: about 20% of the accessible area). Not relevant to the spacing rule. Recorded for completeness.

## 6. Things not checked

- I did not compare any number here with `data/sorting_hat/cys8_spacing.yaml` or with Kubicek, Seidl-Seiboth or Mgbeahuruike.
- I did not read the Jensen or Xu main texts.
- I did not check the references cited inside the reviews.
- I did not run any sequence analysis or regex on any protein.

## Follow-up, 2026-10-08 (resolution of the Jensen 2010 gene IDs)

`tools/resolve_jensen_ids.py` resolved all 50 Jensen 2010 gene IDs to protein sequences (`jensen_sequences.faa`, `jensen_resolution.tsv`) from FungiDB-31 (local), the local A1163 UniProt proteome and UniProt REST (Af293 AFUA_ IDs). A. niger ATCC 1015 JGI IDs were mapped as ASPNIDRAFT_<number>. These are current annotations, not the 2010 gene models.

Checks (commands in this session):
- The cysteine pattern stated in the Jensen supplement is found exactly in the resolved sequence for 41 of 49 proteins. It is not found for 8 (An07g03340, An08g09880, AN6401.2, ACLA_001890, ACLA_048810, ACLA_072820, ACLA_018290, ACLA_007980). Possible causes: a different gene model, or an error in the transcription. Not resolved. Those 8 are not used until checked.
- `hmmsearch --cut_ga` of the seven hydrophobin-class Pfam models (Pfam 38.2): 43 of 50 have a hit. Seven do not: AO090012000143, ATEG_10285, ATEG_08089, AFLA_060780, AFLA_014260, AFLA_063080, ACLA_001890. Jensen reported 5 without a Pfam domain in 2010. The difference is the Pfam release and model set.
- Evidence level: gene prediction only (Jensen's pattern, size and signal-sequence screen). Tier "literature, predicted", held out from training.

## L1 result (2026-10-08, `lp_clusters.py`; commands and tables in this folder)

- Joint clustering of the 50 LP proteins with the 131 T2 proteins (MMseqs2 17, 30% identity, coverage 0.5): 12 LP clusters; 8 have a T2 member; **4 clusters (5 proteins) are reserved for the v2 test** (`v2_reserve.tsv`). This clustering is used only for the reserve. Folds use `clusters_positives.tsv`.
- Gap check of the 8 proteins whose stated pattern was not found (`unverified_gap_report.tsv`):
  - **Five ACLA proteins are a row offset in the transcription, not a sequence problem.** The actual gaps of each ACLA sequence equal the stated gaps of the next row, in a cycle:
    ACLA_048810 actual (7,39,21,5,17) = stated for ACLA_072820; ACLA_072820 actual (5,32,6,5,13) = stated for ACLA_018290; ACLA_018290 actual (7,36,18,5,17) = stated for ACLA_007980;
    ACLA_007980 actual (7,16,6,5,26) = stated for ACLA_001890; ACLA_001890 actual (7,33,11,5,15) = stated for ACLA_048810. The `cys_pattern_stated` column of these five rows in `protein_lists.tsv` is shifted by one row. Not corrected in the table. The 2010 supplement image would confirm.
  - An07g03340 has 13 cysteines (258 aa). The stated pattern matches only its first five gaps. Unresolved.
  - An08g09880: last gap 9 in the sequence, 10 stated. Unresolved.
  - AN6401.2: last gap 16 in the sequence, 35 stated. Unresolved.
- All 8 stay excluded from the reserve and from any use until checked.
