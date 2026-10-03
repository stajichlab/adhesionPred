# H99 curation pilot, pass 2: report

Curator: Claude Code (claude-sonnet-5-5), 2026-10-02. All counts below were computed from the TSV files in this directory.
I opened no score, prediction or blinded file. I did not edit `h99_pilot/`. I ran no `git commit`.
I did not open the H99 GOA file. Gene-level `expected_label` values come from my own rows (and, for CDA1, from the reviewer's statement about GOA ISS terms). They are not checked against GOA.

## Files
- `candidates.tsv` (43 entries) and `candidate_list.sha256`. The hash was written before `draft_rows.tsv`. `sha256sum -c` passes.
- `draft_rows.tsv` (32 rows), `row_provenance.tsv` (origin, text file and quote check for each row), `gpi_leads.tsv` (3 rows), `query_log.tsv` (36 entries), `texts/` (saved XML, HTML and abstract files).
- Helper scripts: `fetch.py`, `ctx.py`, `qc.py`, `build_rows.py`, `es.py`, `up.py`, `lg.py`.

## Counts (draft_rows.tsv)
- Rows: 32. Distinct genes (UniProt accession): 20.
- Rows by expected_label: P-ext 23, N-int 5, N-sec 4. Genes by label: P-ext 12, N-int 4, N-sec 4.
- Distinct genes with expected_label P-ext and evidence_level direct: **12** (APH1, BIM1, CDA1, CDA2, CEL1, CIG1, CNAG_05312, LAC1, PLB1, PQP1, QSP1, cnap1). All 32 rows are `direct`. There are no `transfer` rows.
- Rows by evidence_code: IDA 29, HDA 3.
- selected_by_predictor: no 28, yes 3 (APH1 x2, cnap1), unknown 1 (PLB1 raft row).
- Text source: 31 rows from full text, 1 row from an abstract only (PLB1, PMID 17947228, no PMC copy exists).
- Machine quote check: 32 of 32 quotes found in the saved text file (whitespace normalised). Quotes were re-copied from Europe PMC XML where it exists. For 4 papers only a PMC HTML page was available (PMC1828480, PMC1398056, PMC98673, and PMC3318296 for the skipped RAC1). For PMC13220255 the XML came from NCBI efetch.
- New vs carried over: 18 rows come from pass 1 (12 carried with only a date added or a strain note fixed, 6 corrected: CDA1, APH1 x2, cnap1, SOD2 x2). 14 rows are new.
- New genes (9, not in pass 1): CIG1, SOD1, BIM1, CEL1, UGG1, MNS1, MNS101, NOP1, VCX1. New P-ext genes: CIG1, BIM1, CEL1 (3). Pass 1 had 8 P-ext genes (APH1 was `ambiguous`). Now 12.
- GPI leads: 3 rows (CDA2 x2, PLB1 17947228). Quotes re-copied. The CDA2 line 2 quote now ends "in C. neoformans." (the pass-1 copy lost the species name).
- Candidates: 43 entries. Status: row_written 20, no_usable_evidence 7, skipped_with_reason 15, not_screened 1.

## Pass-1 corrections applied
- APH1: expected_label P-ext, selected_by_predictor yes. Quote re-copied with "C. neoformans, H99".
- cnap1 (May1): selected_by_predictor yes.
- SOD2: code IDA. The pass-1 quote was an abstract sentence and did not support the strain or the assay. I replaced it with two full-text sentences (Fig. 6 fractionation, mitochondrion and cytosol). The strain is H99 only by inference from two sentences in PMC7961099: the original "H99 clinically isolated strain" is called the "stud" strain, and "the current studies presented here used the fully virulent 'stud' strain". No sentence says "we used H99" for this assay. This is stated in the row note.
- CDA1: expected_label P-ext, GO:0016020, full sentence with "Cda1CS" re-copied. A second CDA1 row was added from PMID 40424315.
- Actin (CNAG_00483): dropped. H99 GOA has an IEA actin cortical patch term whose ancestors include cell periphery. The derived label is `unlabelled`, so actin cannot be a negative control (reviewer finding; I did not re-check GOA).
- Strain notes fixed with full-text sentences: CNAG_05312 (H99 stated), CRZ1 ("created from H99"), PLB1-GFP (H99 stated in Methods; "WT cells expressing Plb1-GFP" in the legend).
- PLB1 raft row (16524904) and LAC1 (17101662): the strain was "not stated" in the abstracts. I got the full text from PMC HTML pages. PLB1: "(H99)". LAC1: "Strain H99 lac1 ura5 ... recipient strain". Both are now resolved.
- GPI lead PLB1 15826239 dropped. The anchor-motif experiment was done in S. cerevisiae (heterologous), and the pass-1 quote dropped "(GPI-)".

## New candidates from the reviewer
- CIG1 (J9VUB8, CNAG_01653): 2 rows. HDA from the 26453029 secretome text. IDA from PMID 41313167 (Cig1:mCherry mostly at the plasma membrane, H99, overexpressed from the ACT1 promoter). The secretome sentence names "Cig1"; the CNAG ID comes from the paper's table and from UniProt.
- SOD1 (J9VLJ9): 1 row, cytosol (N-int), PMC7961099, same H99 inference as SOD2. Open conflict: the 16524904 abstract says SOD1 is enriched in raft membrane fractions. If the owner adds that as a membrane row, the gene is no longer N-int.
- RAC1 (22327008): skipped. The PMC HTML page says the parental strain was JEC21 (serotype D). It is not H99. I found no CNAG ID.
- ECM33, GAS1, PDA1: screened. PubMed gave 0 hits for each. For Europe PMC I read 25 titles per query (reviews or unrelated work). I found no H99 localisation experiment. Status no_usable_evidence. The pass-1 query with 619 hits was replaced by gene-specific queries. About 100 remaining hits on the combined query were not title-screened.

## Other new rows (extended search)
- BIM1 (J9VHN6; PMID 31932719): plasma membrane, wall, and partial supernatant (3 rows, one sentence). The wild-type control is H99 (DTY758). The parent strain of the bim1 deletion and the Bim1-HA strains is not named in the text I read. Treat H99 as unverified for these 3 rows.
- CEL1 (J9VH79; PMID 37099613): wall association by Zymolyase. H99 background stated. Accession matched by gene name CEL1 (the paper gives CNAG_00601; the UniProt header gives no ORF).
- LAC1 second paper (PMID 11500433): immunoelectron microscopy, cell walls of H99 and B-3501. Mapped to LAC1 through the CNLAC1-GFP mention elsewhere in the text.
- UGG1, MNS1, MNS101 (PMID 40424315): GFP fusions colocalise with ER-Tracker (N-sec, IDs from Methods). The GFP strains were made in "the WT strain". The paper builds its deletions in H99, but the WT strain is not named H99 in that sentence.
- NOP1 (J9VR32) and VCX1 (J9VDQ4) from the organelle marker toolkit (PMID 42085564; H99 alpha mCherry, KN99a GFP). These are overexpressed H3-promoter fusions. NOP1 is mapped to the single fibrillarin entry by protein name. The VCX1 row is weak: the sentence reports vacuole morphology, not Vcx1 localisation itself.
- CDA1 second row (PMID 40424315): anti-Cda1 western, insoluble fraction (wall plus membrane, not separated). I did not write a row for secretion because the text does not state the wild-type secreted level.

## Rows with a strain caveat
The spec requires H99 or a stated H99-congenic strain. These rows do not fully meet this:
- BIM1 x3 (parent strain of the tagged strains not named).
- SOD2 x2 and SOD1 (H99 inferred from the "stud" sentences).
- UGG1, MNS1, MNS101 (WT strain not named H99 in the sentence).
- CDA1 row from 40424315 (WT strain for the blot not named).
- NOP1, VCX1 (H99 stated for the strain set; the tag per marker not identified).
KN99 and KN99alpha rows (CDA2 x3, CDA1, QSP1, PQP1) are accepted as H99-congenic. The QSP1 Methods sentence states it.

## Papers screened
- Full text read for some part: 22 papers (22354955, 30459196, 27212659, 25227465, 26453029, 27977806, 23251520, 33567338, 41313167, 40777840 (partial), 31932719, 37099613, 40424315, 42085564, 11500433, 17101662, 16524904, 22327008, 33114434, 36670130, 42016327, 42367114). Plus grep-only screens of 38742885, 34566931, 38953322, 33158259, 35737716, 39119141, 28039134, 25841021.
- Titles screened in Europe PMC lists: about 25 per query for 18 queries. I did not count unique titles. Many repeat across queries.
- Pass-1 papers are not re-counted here.

## Skipped and why
- RAC1 (JEC21), PMA1 (two ATPase entries, no ID in text), DNJ1, DCP1, CAP6, KTR3, RAB5, MJR1, ANT1, ATG8 (accession not unique or not found by name; CNAG IDs are only in supplementary tables), CMP1 (no CNAG ID in text), OLP1, WOS2, YPT7 (no unique accession or no location sentence resolved in the time available), STR1/STR3 (S. pombe), Ricin B-like proteins (no localisation experiment), Cfl1 family 28039134 (CFL1 has no UniProt entry found by GN; strain of the secretion assay not resolved), strain secretomes 25841021 (KN99alpha; protein lists have no CNAG IDs in the text).
- Pass-1 skips carried over: CPL1, YOR1, LAC2, CDA3, URE1, MPK1, KRP1, CFO1, MP88, Cap64.
- Observation, not acted on: 28039134 gives CPL1 as CNAG_02797, and UniProt J9VHS9 has that ORF ("CPL1-like"). The pass-1 seed accession for CPL1 is J9VQE5. I did not resolve this.
- not_screened: 28039134 abstract-level screen for H99 and location was not done beyond the text grep.

## What I could not verify
- Expected labels were not checked against H99 GOA. Gene-level labels use my own rows only (CDA1 uses the reviewer's note).
- Whole-paper reading was not done for most papers. I searched full texts with regular expressions for localisation sentences. A paper with a relevant sentence phrased differently could have been missed.
- 40777840 (Cig1 liposome paper, H99) and 38742885 (May1 and cell wall) were only grepped.
- Search stopped after about 85 minutes, not because new searches stopped giving proteins. The last searches (secretion, wall release) still produced leads that I could not resolve to accessions. I have no data to say how many more H99 proteins exist.
- Europe PMC supplementary tables (CNAG IDs for the marker toolkit, secretome lists) were not read.

## Tool problems
- Europe PMC returns full text only for open-access records. PMC1828480, PMC1398056, PMC4452572, PMC5311391, PMC98673 and PMC3318296 were not available as XML. The PMC web page (`pmc.ncbi.nlm.nih.gov/articles/...`) with a browser user agent worked, and I converted it to text. Spaces around inline italics may differ from the print version, although I stripped inline tags.
- NCBI efetch gave XML for PMC13220255 only. It gave HTTP 429 for PMC12277061 (40543703) and an empty body for PMC11228865 and PMC4452572.
- PubMed `get_full_text_article` returned empty text for PMC3318296. Gene-name searches in PubMed returned 0 for Ecm33, Gas1 and Pda1.
- UniProt REST gave HTTP 400 for a long OR query. I used the local proteome FASTA headers (`UP000010091.fasta.gz`, sequence file, no scores) to check accessions.
- Writing to `/scratch/jstajich` directly was refused. I used the session scratchpad.
