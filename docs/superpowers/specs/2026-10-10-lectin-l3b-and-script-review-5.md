# Review 5: lectin class L3b edits, L3b report and report script

*2026-10-10. Reviewer: Claude Code (claude-opus-5-5), independent review. Branch `lectin-review5` = `lectin-class-spec-r5` + merge of `lectin-report-script` (tip 28541c8) and `lectin-l3b` (tip fbe470d). No conflicts. I did not edit the spec, the script or the tests.*

## Verdicts

| Part | Verdict |
|---|---|
| A. L3b edits to the spec and the L3b report | **Pass with changes.** The quotes I checked are in the sources. Two places say more than the sources (A1, A2). Some text is stale (A3). |
| B. Report script `analysis/lectin_class_report/lectin_report.py` | **Pass with changes.** The counts are correct on real data. Two must-fix items: a note text that the L3b result made false (B1), and silent misreading of signal peptide files in other formats (B2). |

## What I ran

| Check | Result |
|---|---|
| `/usr/bin/python3.12 -m pytest tests/lectin_class_report` | 14 passed |
| Script on the 25 proteomes of `review3/rescan_sample.txt` (HMMER 3.4, Pfam 38.2, login node, 4 CPU, 15 s) and `compare_committed.py` against `review2/whole_proteome_hits.tsv` | **MATCH.** PF07938 36 = 36 proteins, 14 = 14 proteomes, max score difference 0.0. FB_lectin 2 = 2, LDL 1 = 1, Gal-bind 3 = 3, CCL2-like 0 = 0. |
| Same run against `r3/cbm_whole_proteome_hits.tsv` (F24 models; `compare_committed.py` does not cover them) | Chitin_bind_1: the same 62 proteins and scores (`diff` empty). CBM_14: 0 = 0. |
| SignalP 6 fast mode (binary of `review3/signalp6_check.sh`, reports 6.0h) on `f27_six/six.faa` | Data rows identical to `f27_six/signalp6_prediction_results.txt` (only the time stamp differs). |
| BLASTP (ncbi-blast 2.17.0+) of `six.faa` against `fleA_Q4WW81.faa` | Identical to `blastp_six_vs_FleA.tsv` (identity 38.39 to 44.70%, aligned 311 to 323 aa, E 4.49e-64 to 7.84e-79). One extra HSP for *O. corvina* (24.3%, 214 aa, E 8.4e-9) is not in the table. This does not change F27. |
| The six sequences against their source FASTA (`six.tsv` paths); FleA against current UniProt Q4WW81 | All six identical (lengths 327 to 356). FleA identical (MD5 of the sequence). |
| Script on *C. neoformans* H99, *C. immitis* RS, *A. fumigatus* Af293 (`_workdir/sorting_hat/*.faa`) | The method note "RicinB_lectin_2 in C. neoformans H99" is true: 3 proteins at sequence GA, 1 at domain GA. |
| Archive hashes `docs/reports/data/lectin_l3b/archive_sha256.txt` | 39 of 39 files OK. |
| `bash -n` on the sbatch file; `set -euo pipefail` with `module load hmmer/3.4` | Both OK. No `BASH_SOURCE`, no path from the script location. |

## A. L3b report and spec edits

### Claims checked against the sources (20)

Archived full texts are in `_workdir/curation_texts_archive/lectin_l3b/`. I searched them with a context script. Abstracts came from the PubMed tool.

| # | Claim (spec or L3b report) | Source checked | Result |
|---|---|---|---|
| 1 | Kerr 2016: internalisation count, "Internalized conidia were counted ... calcofluor white" | 27058347.txt | Found |
| 2 | Kerr: "50% less" (RAW264.7), "reduction in phagocytosis of 40%" (human AM); 6 replicates, 3 experiments or 3 donors | 27058347.txt | Found |
| 3 | Kerr: "PFA-fixed WT or ΔfleA"; mucin binding of three deletion strains "much weaker"; 2EHex and fucose inhibit WT conidia | 27058347.txt | Found |
| 4 | Kerr: no epithelial cell assay | 27058347.txt | "epithelial" occurs 3 times, all in the Introduction or as background. No assay. Agrees. |
| 5 | Kerr: no complement ("Despite repeated attempts"); "high on resting conidia"; "largely absent in hyphae"; no signal peptide by SignalP 4.0; "shed significant amounts" | 27058347.txt | Found |
| 6 | Kerr: dermatophyte genus names in the Results | 27058347.txt | Found: "dermatophytes (Trichophyton, Microsporum and Arthroderma)". The correction note in the spec is right. |
| 7 | Sakai 2019: no obvious difference in adhesion to A549, J774, MH-S; phagocytosis rate not different | PubMed abstract; a third WebFetch of the Oxford page (by me) | The third fetch returned the same sentences. **I could not get the PDF or a raw page**: curl with a browser agent gets HTTP 403; Unpaywall `oa_status` closed; Europe PMC `inEPMC` N. The text is still from a summarising tool. |
| 8 | Richard 2018: no fleA mutant; exogenous FleA reverses inhibition; αMF restores activity; adhesion 97.5 and 99.2% vs 91.9 and 98% | 30356167.txt | Found. No mutant, deletion or knockout word in the text. |
| 9 | Fernandes 2017: knockdown to 60, 40, 60% of WT; A549 adherence by CFU at 2 h did not differ; phagocytosis by peritoneal macrophages at 4 h similar; 27-kDa band | 28720727.txt | Found (see A6 on the knockdown numbers). |
| 10 | Alegre 2014: PADG-3347 by proteomics; five exons, four introns; "single 28-kDa band"; 70 kDa "assumption"; rPCNexon4 "only residual" NAGase, four times lower; both rPCNs bind laminin | 24743161.txt | Found (see A5 on the 28 kDa product). |
| 11 | Pitangui 2021 is a review | 34532342.xml (`article-type="review-article"`); PubMed type "Journal Article; Review" | Agrees |
| 12 | Pitangui: GH18 29-218; sites 30-37, 55-62, 98-107, 120-126; Prosite 98-99 and 123-126; template 2XTK, 47.47%; silencing sentence cites Fernandes 2017 | 34532342.txt | Found |
| 13 | Tronchin 2002: 32-kDa fucose-specific lectin (Ishimaru 1996, Letter) cited as separate; protease-deficient mutant shows distinct protein | 12438366.txt | Found (see A8) |
| 14 | Powers-Fletcher 2011 is a calnexin paper; FleA and fucose do not occur | 22163332.txt | "calnexin" 109 times; "FleA", "fucose", " AFL" 0 times. Agrees. |
| 15 | Houser 2013: formaldehyde fixation, "PBS + 0.2% Triton" before antibody | 24340081.txt | Found |
| 16 | *Histoplasma*: Ramos-Martinez 2026, 10 isolates, beta-galactose, partial purification, sequence not known | PubMed abstract; 41805922.txt ("identify its genetic sequence") | Found |
| 17 | *Histoplasma*: Taylor 2004 beta-galactose inhibits attachment to erythrocytes | PubMed abstract | Found. The abstract also names an ionic mechanism with sialic acid and heparan sulfate; the spec omits it. |
| 18 | *Coccidioides* CFP28: 231 aa, 19 aa signal peptide, "single DOMON-like type 9 CBM (CBM9-like-3)", rises with GlcNAc, binding not tested | PubMed abstract 25023485 | Agrees. The abstract has no binding assay. |
| 19 | *Blastomyces*: Sugar and Picard 1990, sugars did not inhibit; "non-lectin-like attachment" | PubMed abstract | Agrees. The study used an avirulent strain; the spec does not say so. |
| 20 | *Talaromyces*: Hamilton 1998, NANA abolished adherence; "pH shift" caveat | 9826390.txt | Found |

### Findings, ranked

**A1 (medium, change before the spec is final). The reason given for the Kerr/Sakai phagocytosis difference is stronger than the sources.** F14 (c) and "What the papers show" say the papers differ on phagocytosis "with different conidia treatment (fixed against agar-overlay colony counting, which needs living conidia)". The colony-count method is the Sakai **adhesion** and **8 h survival** method. The method of the Sakai **phagocytosis-rate** assay (Fig S1) was not seen (L3b report section 5 says so). So "different conidia treatment" is not shown for the phagocytosis comparison. The L3b report says more: section 3.1 item 2 writes "phagocytosis of live conidia" and proposed edit 1 writes "live ΔfleA conidia", with no "my inference" mark. In my third fetch, the Sakai Discussion gives its own reasons: a different strain (Kerr used TGJF5.3), a different stage of fleA expression, Kerr's lack of a complemented strain, and "different cell lines used in the phagocytosis and adhesion assays". Proposed text: "The papers differ on phagocytosis. Possible reasons named by Sakai are a different strain, a different cell line and the lack of a complemented strain in Kerr. The conidia treatment differs for the adhesion assay (Kerr used PFA-fixed conidia; Sakai counted colonies). For the phagocytosis-rate assay the Sakai method was not read." Mark the live-conidia point as an inference wherever it is kept.

**A2 (medium). "Only Sakai tested adhesion to cultured cells" is too narrow.** Richard 2018 measured conidial adhesion to BEAS-2B cells (97.5 and 99.2%) and the L3b report quotes it. Sakai is the only **knockout** adhesion test. F14 and the L3b report should say "the only gene-loss adhesion test to cultured cells" in both places.

**A3 (medium). Stale text left in the spec after the L3b edits.** (1) Section 3.1 table, PF07938 row, still reads "Sakai 2019 found probably no change in conidial adhesion after fleA loss (F14 c, not verified in the full text)". L3b replaced "probably" everywhere else. (2) Section 6.2 item 6 holds the new endemic-fungi wording, but the script note `ENDEMIC_NOTE` does not (see B1). (3) Section 10 "Revision 5 items" item 1 still says the Sakai phagocytosis text comes from review 4's fetch; L3b fetched it again. These are text errors, not new claims.

**A4 (low). The 70 kDa question.** The spec says the conflict is "partly resolved". Alegre 2014 states the re-assignment to the 28-kDa product and calls the 70 kDa a prior assumption. That supports the spec wording. Neither paper explains the 70 kDa band, and the spec says so. No change needed.

**A5 (low). 28 kDa is a recombinant band.** Alegre 2014 reports 28 kDa for refolded recombinant full-length protein (rPCNfull) on SDS-PAGE, and 22 kDa as the **predicted** size of rPCNexon4. The spec text could be read as a native size. Add "recombinant".

**A6 (low). Knockdown wording.** Fernandes 2017 reports expression at 60, 40 and 60% of WT, that is, reductions of about 40 to 60%. The L3b report writes "knockdown (40 to 60% of wild-type expression)" in section 3.3 and "to 60%, 40%, and 60%" in the table. Both appear, so a reader can mix "40 to 60% of WT" with "a 40 to 60% reduction". Use one form.

**A7 (low). Two small omissions in the endemic-fungi table.** Sugar and Picard 1990 used an avirulent *Blastomyces* strain (abstract). Taylor 2004 also names a sialic-acid and heparan-sulfate ionic mechanism. Neither changes the conclusion.

**A8 (low). Tronchin 2002.** The claim "the authors treat the sialic-acid lectin as a different protein" is supported by the Introduction and Discussion wording ("further studies are needed to characterize this lectin"; the fucose lectin is cited as a separate identification). The paper does not state in one sentence that the two are different proteins. The spec's "treat ... as distinct" is a fair reading; keep "not counted as FleA evidence".

**A9 (no change, information).** Sakai 2019 is still not verified in the PDF. I tried a browser user agent (403), Unpaywall (closed) and Europe PMC (not in PMC). A third WebFetch of the Oxford page returned the same sentences as the two earlier ones, plus the Discussion reasons in A1. Three fetches of one summarising route are not the same as the PDF. The spec's caveat is correct and should stay.

### Contradictions with D7, D9 or the owner-decisions block

I found none. The F14, F17 and F20 rewrites keep the family-only label of D7 and the "no custom model" default of D9. The D7 cell now cites the Sakai result and its limit, which is consistent with F14 (c). F17 and section 1 point 4 agree ("no fungal-lectin gene"). The one inconsistency with the new text is in the **script** (B1), not in the spec.

### F27 numbers

Re-verified in this review: the six sequences match their source proteomes (lengths 327, 336, 337, 337, 337, 356); SignalP 6.0h fast mode gives SP for all six with P(SP) 0.999666 to 0.999810 and the cleavage sites of the committed table; BLASTP identities 38.4 to 44.7% with E 4.5e-64 to 7.8e-79. The spec reports these as a similarity check and not orthology, which matches what the data show.

One point for the record: a signal peptide here is not a property of the Onygenales group alone. In my 25-proteome run, 8 proteins reach PF07938 score 100 or more and SignalP 6 gives 2 of them a signal peptide (both Onygenales: *A. uncinatum* FA49FE15_006535-T1 and *O. corvina* F1215351_005656-T1). The second *A. uncinatum* protein (FA49FE15_005258-T1, score 213.6, named in F27 as not in pool S) is "OTHER", which fits the spec text.

## B. Report script

### Correctness

The hit parsing is right. `hmmsearch --cut_ga --tblout` is used, and the fields read are target (1), query accession (4), sequence E-value (5), sequence score (6) and best-domain score (9). The script states in the report that its counts come from the sequence table, while the tool's Pfam module reads `--domtblout`; it prints `domain_ge_ga` per hit for that reason. I confirmed the named example: *C. neoformans* H99 has 3 RicinB_lectin_2 proteins at sequence GA and 1 at domain GA. GA values are read per model from the extracted HMM block, and a model with no `GA` line is a fatal error.

Other things I checked in the code and in runs: `.gz` and `.zst` are decompressed into the temp directory and removed after the search (hmmsearch cannot rewind a stream); the temp directory is `--tmpdir` or `$SCRATCH`, and the script stops when neither is set; no path is derived from the script location; the sbatch wrapper takes every path from the environment, uses `${SCRATCH:?}`, has a `trap` to clean up, and copies output back with a temporary name. Output files are written with a temporary name and `os.replace`. Two inputs whose base names collide are a fatal error.

Wording: the fixed text, the no-call text and the PF07938 label are character-for-character the D7 and 6.2 texts of the spec. I found no function claim, no call and no status output. `PF07938_LABEL` is exactly "PF07938 (AAL-like / FleA-family) domain hit".

Hand-copied values: the `MODELS` table copies H and U from spec 3.1, and the report says so in the text. I checked all 15 rows against the 3.1 table and against F24 and F25. They agree, including "not scanned; at least ..." for the pool S values. The 15 rows are the M1 and M2 rows of the 3.1 table; M3 and the narrow-tier models are notes only, as the spec asks.

### Findings, ranked

**B1 (must fix). `ENDEMIC_NOTE` contradicts the L3b result.** The script prints "Coccidioides, Histoplasma and Blastomyces: the literature searches of spec F17 found no fungal lectin literature." After L3b, F17 and section 6.2 item 6 say the opposite in part: *Histoplasma* has a published beta-galactose-binding activity (no gene), and *Coccidioides* CFP28 has a CBM9-like domain (no binding data). The script was written on the pre-L3b text. Replace the note with the 6.2 item 6 wording, which also covers *Paracoccidioides*. The script omits the paracoccin sentence of 6.2 item 6 in this note, although `PARACOCCIN_NOTE` carries it; keeping both is fine.

**B2 (must fix, or document). `read_signalp` accepts files it cannot read and answers "no" for every protein.** The rule is: any second word that is not `OTHER` or `NO_SP` means a signal peptide. Tested:
- SignalP 4.1 short format: column 2 is a number (`0.110`), so every protein is read as **having** a signal peptide. In my test the protein whose `?` column is `N` came out "yes".
- A `yes/no` table (`no`, `N`, `0`): every value is read as a signal peptide.
- SignalP 5.0 format: `SP(Sec/SPI)` and `OTHER` are in column 2, so this format works by luck.
- A SignalP 6 file with CRLF line ends: the last column keeps `\r`; the prediction column is not last, so it still works. A **one-column** CRLF list would store IDs with a trailing `\r` and match nothing.

The failure is silent and it flips the (a2) section. Two cheap fixes: accept only a known set of prediction words (`SP`, `LIPO`, `TAT`, `TATLIPO`, `PILIN`, `OTHER`, `NO_SP`, `Y`, `N`, `yes`, `no`) and stop on anything else; and `line.rstrip()` instead of `rstrip("\n")`. Until then, the docstring should say that only SignalP 5/6 output and a bare ID list are supported.

**B3 (medium). A duplicate ID inside one proteome is not reported and is counted twice.** The check compares ID sets between proteomes. With one file holding the same ID twice, the hits table gets two identical rows and the counts table says 3 hit proteins for 3 records, 2 of which are the same ID. The spec's option (c) report is a count of proteins, so a repeated ID inflates it. Count distinct IDs per proteome, or warn when `len(ids) != n_records`. `scan_fasta` already has both numbers.

**B4 (medium). A protein with two domain hits is one row per model, which is right, but the counts column name can mislead.** `n_hit_domain_ge_ga` counts hit **rows** whose best-domain score reaches the domain GA, not domains. With one protein and two models (my test protein G4N991 hits LysM and CVNH) the per-model counts are correct. No data bug; rename or explain the column.

**B5 (low). Empty FASTA is a fatal error for the whole run.** `hmmsearch` exits 1 on an empty file and the script stops, so one empty proteome in a list of many loses the run. A header-only file (no residues) works and gives 0 proteins. Decide whether an empty file should be a warning with a 0 row.

**B6 (low). ID matching for the name list.** `lookup_sp` splits `sp|ACC|NAME` to find the accession, but `named_section` matches the name list on the full ID only. A name list that holds `Q4WW81` finds nothing in a UniProt-header proteome, while `sp|Q4WW81|LECF_ASPFU` works. The docstring promises accession matching for signal peptide data only, so this is a documentation gap more than a bug; make both paths the same.

**B7 (low). An ID with spaces.** The first word is the ID, as documented. A header `>my protein one` gives ID `my`. Correct but easy to misread in the report; the report could print the full header for named rows.

**B8 (low). Lowercase and CRLF FASTA.** Both work (HMMER handles them). Scores are the same as for the uppercase, LF version in my tests.

**B9 (low). No `README` or run record in `analysis/lectin_class_report/`.** The committed example report (`f27_six/example_report_six_proteomes.md`) does not say which six proteomes were used or give the command. Add the command line and input list, as the other `analysis/` folders do.

**B10 (information).** `blastp_six_vs_FleA.tsv` keeps one HSP per query. My re-run gives a second, weak HSP for *O. corvina* (24.3%, E 8.4e-9). State the filter in the file header.

## C. For the owner, not asked by the spec

1. **The (a2) section needs signal peptide data that the report does not produce.** Under (c-i) the user must run SignalP 6 separately and pass the file. Is that part of the v1 report, and who runs it? Without `--signalp` the report has no (a2) section at all.
2. **Which proteomes the v1 report covers** (14 Phase C proteomes, the 653 discovery proteomes, or the 6 of F27) and where the output is stored. The spec does not say.
3. **How a protein count is defined when IDs repeat** (B3), which matters if Fungi_5k files ever hold duplicates.
4. **Whether the report must name the six F27 proteins even in runs that do not contain them.** Today a name list row that is absent prints "not in the input proteomes, or no PF07938 hit", which mixes two cases.
5. **Whether CFP28 of *Coccidioides* gets a Pfam scan** (the L3b report lists it as an open item, and it is a one-sequence login-node run).
6. **Whether the Sakai 2019 PDF should be obtained through the library.** Every route I tried is closed. The adhesion direction in the spec rests on one summarising tool.

## Must-fix versus minor

- **Must fix before the report is used:** B1, B2.
- **Should fix:** A1, A2, A3, B3.
- **Minor:** A4 to A8, B4 to B10.
- **No change:** A9 (the caveat is right as written).