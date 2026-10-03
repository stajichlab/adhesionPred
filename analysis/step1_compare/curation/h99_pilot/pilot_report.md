# H99 literature curation pilot: report

All counts below were computed from the TSV files in this directory, except "papers screened" (my own tally, see below).

## Counts
- Log lines in query_log.tsv: 67. PubMed search_articles calls: 31. PubMed full-text calls: 9. PubMed metadata calls: 7. Europe PMC curl calls: 15. UniProt REST calls: 4. One local FASTA grep line.
- Papers screened (title and abstract read through metadata): about 64 by my tally. Many were off topic (reviews, drug studies, other organisms).
- Papers that gave at least one row: 11 (PMIDs 16524904, 17101662, 17947228, 22354955, 23251520, 25227465, 26453029, 27212659, 27977806, 30459196, 33567338).
- Candidate proteins in candidates.tsv: 26 entries. Status: row_written 12, skipped 9, no_usable_evidence 5.
- Rows in draft_rows.tsv: 19, for 12 genes.
- Rows by expected_label: P-ext 12, N-int 4, ambiguous 2, N-sec 1.
- Rows by evidence_code: IDA 15, HDA 2, EXP 2.
- Quote source: 13 rows from full text, 6 rows from abstract only (marked [abstract]).
- gpi_leads.tsv: 4 rows, covering CDA2 (2 rows) and PLB1 (2 rows). One PLB1 row is from heterologous expression in S. cerevisiae.
- I stopped before 40 candidates because new searches gave few new H99 proteins.

## Quote checks
- 12 rows (11 full-text rows and 1 abstract row) were matched by script against text I saved (exact string match after whitespace normalisation). One mismatch was found and fixed (an abstract spelled MbetaCD in one retrieval and with the Greek letter in another).
- 7 rows were not machine-checked, because the text came inline and was not saved: QSP1, PQP1 (PMC5186401 full text), PLB1 from 17947228 (abstract), LAC1 (abstract), CRZ1 (abstract), SOD2 (2 rows, abstract). I copied them from the tool output. The reviewer should check these first.
- Some quotes carry text artifacts from the PubMed tool, for example "Cda1" with gene names missing in italics, "(, lanes 5 and 6)", "of, H99". I kept them as retrieved.

## Tool problems
- PubMed MCP: get_full_text_article returned an empty full_text for PMC1180731, PMC1398056 and PMC2224146, so those rows use abstracts. Tool output for large papers was written to files.
- Europe PMC: no full text for 9 of 15 requests (150 byte replies), including PMC5186401. The PubMed tool had that one.
- PubMed query translation is strict: queries with several specific terms often returned 0 results. Simple queries worked better.
- UniProt REST: gene-name lookups for CIG1, MPK1, CFO1, CFO2 in taxon 235443 returned nothing. I wrote no row for them because I will not guess an accession.
- Blinding: I did not open any file in the blinded list, and I ran no predictor.

## Limits and uncertainties
- Strain: several abstracts do not name the strain (LAC1, CRZ1, SOD2, PLB1 raft paper, CNAG_05312 secretome). I used evidence_level direct and said so in the note. The reviewer should confirm the strain is H99 or KN99 and not another serotype A or D isolate.
- Cda2 paper: it also saw Cda2 in a cytosolic fraction. The authors call it shearing from the wall. I wrote no cytosol row. The three Cda2 rows each carry the gene-level label P-ext, not the label for that single term.
- APH1: extracellular plus vacuole. The spec says N-sec means no wall or extracellular term, and ambiguous means surface plus internal evidence. Vacuole is not "internal" in that text, so P-ext may be the correct label. I used ambiguous. Please decide.
- CDA1: only a membrane-fraction statement, drawn mostly from a catalytic mutant against wild type. Weakest of the localisation rows.
- cnap1 mapping to UniProt J9VS02: made by protein name (Major aspartyl peptidase 1 = May1), not by CNAG ID.
- CRZ1 is nuclear only under some conditions; it forms puncta under salt and heat shock. SOD2 has a mitochondrial and a cytosolic isoform. Both give N-int.
- Selected_by_predictor: I set "no" when the paper did not pick proteins with a predictor. I did not read methods sections for every paper, so "no" for proteomic hits means "chosen by mass spectrometry", not "checked".
- PLB1 GFP rows (PMID 25227465): the paper states H99 for its proteomics. I did not confirm the strain of the GFP lines in the text I read.
- Seed proteins with no row: CPL1 (only var. neoformans paper found), YOR1 (paper says localisation failed), CDA3 (no experiment found), LAC2 (abstract gives no location). Seed proteins QSP1, PQP1, CDA1, CDA2, LAC1, cnap1, PLB1 have rows. I did not find literature for GO-labelled YOR1/CDA3/LAC2 surface location in this run, so their existing GO labels have no literature row here.
- Not pursued: the EV and secretome papers beyond 26453029 and 25227465 (many were reviews). Full-text EV protein lists (for example PMID 18039940 supplement) were not available to me, so no EV proteomics rows were written except those from secretome papers.

## Estimate for a longer run
I have no measured basis for a firm number. From this run: 11 papers gave 19 rows (12 genes) from about 64 screened papers. If a longer run could read supplementary protein tables (secretome, wall and EV proteomics) it could add many HDA rows, but I could not access any table here. For text-only evidence, I guess 10 to 25 more genes could be found with focused searches by gene name (for example CDA3, MP88, URE1, SOD1, CFO1, GFP-tagged intracellular controls), and Europe PMC or publisher full text. This is a guess, not a calculated figure.
