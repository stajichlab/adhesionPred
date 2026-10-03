# U. maydis curation pilot: independent review

Reviewer: Claude Code (claude-opus-5-5), 2026-10-02/03. Per-row verdicts are in `review_report.tsv`.
I opened no score, prediction or blinded file. I ran no predictor. I did not read `_workdir/step1_compare/phasec/`.

Method:
- I re-fetched all 12 cited papers myself (Europe PMC `fullTextXML`). I did not use the curator's `texts/`.
- I checked all 12 PMIDs in Europe PMC. Each one resolves to the cited paper about the named *U. maydis* protein.
- A script matched each quote segment against my copies. 21 of the 30 rows matched after whitespace and quote normalisation. The other 9 matched after I removed all non-alphanumeric characters. The differences come from superscript reference numbers and subscripts (for example "AvitagHA30", "(Fig. 2a)25", "Cmu11-21"). No word differs. Result: 30 of 30 rows are verbatim.
- Labels: I wrote my own script around the repo's `labels.classify`, `go_obo.parse_obo` and `gaf.filter_gaf`. Inputs: `MYCMD-uniprot.gaf.gz` (taxon 5270) and `go-basic.obo`. I computed each gene's label from GOA alone, from the curated rows alone, and from both merged, under all 3 policies.
- Accessions: I checked all 23 against the `UP000000561.fasta.gz` headers. I checked the UMAG IDs with UniProt REST.

## Row verdicts (30)

| Outcome | Rows | n |
|---|---|---|
| confirmed, no change needed | 7, 8, 12, 14, 15, 16, 18, 27, 28, 29 | 10 |
| confirmed, keep with note (strong promoter, strain or field fix) | 1-6, 9, 10, 13, 17, 19-22, 30 | 16 |
| confirmed, but the curator's strain or selected_by_predictor field is wrong | 11 (Pep1 selected_by_predictor should be yes); 23-26 (strain IS stated) | 5 |
| rejected | none | 0 |
| cannot_verify | none | 0 |

The three groups overlap. Rows 19 and 20 also carry a wrong strain note. All 30 rows have the verdict `confirmed`. I rejected no row. Many rows need a note.

Field errors to fix before merge:
- **Row 11, Pep1: selected_by_predictor must be `yes`.** The paper says: "systematic analysis of effector genes in U. maydis which is solely based on two criteria: the protein should carry a secretion signal and the predicted product should be novel". Fig. 1A cites SignalP.
- **Rows 23-26, sirtuins: the strain is stated.** Methods: "All the strains used in this study are derived from the haploid pathogenic SG200 strain". Remove `strain_unstated`.
- **Rows 19-20, Xyn2 and Xyn11A: the strain is stated.** Table 1 lists "SG200 xyn2:GFP ... P otef: xyn2:GFP" and "SG200 xyn11A:GFP ... P otef". These fusions are overexpressed from the otef promoter. The curator's note says neither fact.
- **Rows 3-6, Stp1-Stp4: replace "panel not checked"** with the panel mapping (question 2).
- **Rows 12, 13 and 16: accession wording.** The proteome headers carry only GN=PIT2 and GN=CMU1. The UMAG IDs come from UniProt REST, not from the header.

Checks that passed for all 30 rows:
- GO terms: 6 terms are used (GO:0005576, GO:0005618, GO:0005886, GO:0005783, GO:0005634, GO:0005739). All are known and not obsolete.
- Evidence code: IDA on all 30 rows. This is right, because each row is a single-protein experiment. No row should be HDA.
- evidence_level: `direct` on all rows. This is consistent.
- Format: all 30 rows have "retrieved 2026-10-02". `reviewer` and `review_date` are empty.

## Answers to the questions

**1. Strong promoter or overexpression.** My rule:
- Keep: native promoter or native locus.
- Keep with note: strong promoter, if the assay has a lysis control or a non-secreted control and no evidence in the paper contradicts it. Such an assay shows that the protein CAN be secreted. It does not show where the native protein is.
- Drop: strong promoter, if the assay cannot separate an intracellular signal from an extracellular one.

| Row | Gene | Decision | Reason |
|---|---|---|---|
| 2 | Rsp3 supernatant | keep with note | otef; Rsp3(d24-60) is retained (non-secreted control). Redundant with native row 1 |
| 3-6 | Stp1-4 supernatant | keep with note | otef; tubulin lysis control. Stp1, Stp3 and Stp4 also have native-promoter rows (7, 8, 29). **Stp2 rests only on this row** |
| 9-10 | Stp5, Stp6 PM fraction | keep with note | otef; Pit1 PM, Cmu1 secreted and mCherry cytosol controls. If overexpression sends the protein to the ER instead, the label is still N-sec |
| 17 | Sts2 | keep with note | pit2 promoter; signal on the hyphal "edge"; no plasmolysis; mCherry-alone control is cytoplasmic. Native Sts2-3xHA is found in the host nuclear fraction, which also implies secretion |
| 19-20 | Xyn2, Xyn11A colony | keep with note | otef; cytoplasmic GFP lysis control; Xyn3-GFP is not secreted (internal negative) |
| 21 | Xyn11A apoplast | keep with note | pit2 promoter; full-length band in apoplastic fluid; free mCherry is also present |
| 30 | Cpl1 wall | keep with note, **weakest row** | otef; the text never says whether the filaments were permeabilised; only "faint" Cpl1-HA in the supernatant and "prominent amounts" in cells. The wall claim is the authors' conclusion. This row is Cpl1's only evidence. The owner may drop it |

Genes whose P-ext label rests only on strong-promoter rows: Stp2, Sts2, Xyn2, Xyn11A, Cpl1 (5 genes).

**2. Stp1-Stp4 panels (Extended Data Fig. 4; methods of 33941900).**
- Panel a, AB33-derived strains, otef: Stp2, Stp3, Stp4. Methods: "To visualize secretion of Stp2-HA, Stp3-HA and Stp4-HA, AB33-derived strains expressing the respective genes constitutively from the otef promoter were generated".
- Panel b, SG200dkex2, otef: Stp1 (legend and methods).
- Panels c and d, SG200 (d with protease inhibitors): this is consistent with the methods for SG200Potef-HA-Stp5 and SG200Potef-Stp6-HA (cOmplete was added). The legend does not name the protein per panel. That mapping is my inference.
- Native-promoter immuno-EM (rows 7 and 8): Stp1 and Stp3, SG200 derivatives.

**3. Sirtuin strain: named.** See the field errors above. SG200-derived is stated in the methods.

**4. Rsp3 wall term.** The figure title plus "around the outside of fungal hyphae" (native promoter, non-permeabilised) is enough for P-ext. For the wall subset it is weak.
- The paper states that the attachment mechanism is unknown.
- A clearer sentence is the Fig. 4b legend: "Secreted Rsp3-HA binds to the fungal cell wall in filaments grown on artificial surface ... without prior permeabilization". That panel uses otef.
- Fig. 4a immunogold reports Rsp3-HA "primarily detected inside the fungal cytosol and in the biotrophic interface". This is probably an in-transit signal. It is not recorded, and it has no label effect.

**5. Pit2 and Cmu1 accessions: confirmed by UMAG ID.**
- UniProt REST: A0A0D1EAR7 = PIT2, ORF UMAG_01375. A0A0D1DWQ2 = CMU1, ORF UMAG_05731.
- Searches by UMAG ID return only these entries. Each gene name occurs once in UP000000561.
- I did not reopen the original Pit2 and Cmu1 papers to read the IDs there.
- Both genes are already P-ext in GOA (EXP), so these rows change no label.

**6. Labels for Stp5, Stp6 and Pdi1.**
- **Stp5 and Stp6 as N-sec: correct under `labels.py`.** PM is in SECRETORY, and GOA has no wall or extracellular term for them. Stp6 has native-promoter "speckles" that the paper calls "on the surface of fungal hyphae". That is GO:0009986 (cell surface). Its ancestors do not meet any `labels.py` set (I checked), so the label does not change. Owner question: both are predicted TM members of a "cell surface-exposed" complex. They are N-sec by rule, but they are debatable negatives for a surface-protein predictor.
- **Pdi1 as N-sec: correct.** GO:0005783 has GO:0012505 as an ancestor (I checked). GOA has only IBA ER and IEA ER lumen. Risk: the Cpl1 paper says "a significant amount is also localized to the cell wall (Marín-Menguiano et al., 2019)". I found no such sentence in the Pdi1 full text. I did not read its S1 Table (DIGE cell-wall extract, 6 proteins). If Pdi1 is in that table, an HDA wall row would make Pdi1 P-ext, and row 22 would become a conflict.

**7. Addendum (Mer1, Tay1, Stp4 second paper, Cpl1).**
- `sha256sum -c` passes for both hash files.
- The addendum does **not** predate the rows. `candidates_addendum.tsv` mtime is 23:20:44. `draft_rows.tsv` mtime is 23:20:36. The status column already says `row_written`.
- The report discloses this, and the addendum keeps `candidates.tsv` unchanged. So the handling is honest. But the addendum is a record after the fact, not a pre-registration.
- `candidates.tsv` (23:02:16, statuses `row_planned`) and its hash (23:02:34) are older than the last write of `draft_rows.tsv`. I cannot prove that they are older than its first write.
- All curator files went into one commit (a28aaa4, 23:28:40). Spec section 6 requires the hash to be committed before any score join. That is still possible, because no join has happened.

**8. Skipped candidates (sample of 8 decisions, from my own fetches).**

| Candidate | Curator reason | My check |
|---|---|---|
| Cda1/2/4/7 (36946727) | periphery only | **Right.** The text says "cell periphery", "puncta on the cell periphery" and "septa". GO:0071944 is not in SURFACE or SECRETORY |
| PR-1La (37716995) | GO cell surface only | **Wrong reason.** Fig. 4c is an otef supernatant western ("Secretion of PR-1La ... Tubulin was used as an internal control for a non-secreted cytosolic protein"). By the curator's own standard (rows 2-6, 19-20) this supports an extracellular IDA row for UMAG_01204 (selected_by_predictor would need checking) |
| Vp1 (34436129) | no UMAG ID | No UMAG ID in the text: right. The evidence is mixed anyway (C-terminally tagged Vp1-HA is absent from the supernatant and the apoplast). Skip is OK |
| Lip3 (42093215) | no UMAG ID | No lip3 UMAG ID in the text: right. Evidence is lipase activity in deletion supernatants, not localisation of Lip3. Skip is OK |
| Abc1 (39487654) | surface punctate | Not checked in depth. UMAG_02796 is named as an ABC transporter |
| Lep1 (33843063) | abstract only | Abstract: "Lep1 is bound to the cell wall of biotrophic hyphae". UniProt has lep1 = UMAG_11940. The H99 review accepted one abstract-only lead (PLB1). Under that precedent a `strain_unstated`, abstract-only wall row is possible. The skip is strict, not wrong |
| Sta1 (32112567) | abstract only | Abstract: "attached to the cell wall of filamentous hyphae". No UniProt gene-name hit for sta1, so the accession is unresolved. Skip is OK |

**Proteomics papers.**
- **18456523 and 18508396 are not proteomics papers.** 18456523 is a review that classifies the predicted secretome. 18508396 is an in-silico analysis of the genome. Neither can give HDA rows. The report's statement that they "could give HDA rows" is wrong.
- **36354345** (PMC9746322) has no Europe PMC XML. It has a PMC HTML page, which I read. It re-annotates an older secretome and makes no new MS measurement. Its HDA source is ref. 10.
- **That source is 22300648** (Couturier et al. 2012, BMC Genomics, PMC3298532, open access). It is LC-MS/MS of the culture secretome of strain "UM521 FGSC 9021": "86 proteins were detected". The protein list is in Additional file 2 (I did not read it). This is a real HDA (extracellular region) source for the reference strain. Proteins were found by MS, so they were not selected by a predictor. **The curator retrieved this PMID in query 1 and did not use it.**
- query_log line 34 says 36354345 was fetched as "full text XML". The curator's `texts/` has only `PMID36354345_abstract.txt`. The report's "abstract only" is right. The log line is wrong.

**gpi_leads.tsv (0 rows): I agree.** Europe PMC query: "Ustilago maydis" AND (PI-PLC OR "phosphatidylinositol-specific phospholipase" OR "omega site" OR "GPI anchor attachment"). It returned 60 hits. None is a *U. maydis* GPI-anchor experiment. The Cda papers (36946727, 33653886) call CDAs only "putatively" or "predicted" GPI anchored. The 0 rows are correct. The header has the spec section 5 columns plus `override_tm`.

## GOA merge (my recomputation)

The GAF filter kept 13,034 CC rows. It dropped 45 `other_db` rows and 1 NOT row. Taxon filtering dropped no rows.

| Group | Genes | GO-only label | Merged label |
|---|---|---|---|
| P-ext from curated rows only (all policies) | Rsp3, Stp1, Stp2, Stp3, Stp4, Pep1, Sts2, Llp1, Xyn2, Xyn11A, Mer1, Tay1, Cpl1 (13) | unlabelled | P-ext x3 |
| P-ext already in GOA | Pit2, Erc1, Cmu1 (3) | P-ext x3 | P-ext x3 |
| Negatives from curated rows only, all policies | Stp5, Stp6 (N-sec), Hst5 (N-int) | unlabelled | expected |
| Negatives from curated rows under no_homology and experimental only | Pdi1 (N-sec), Sir2, Hst4, Hst6 (N-int) | non_iea: expected; others: unlabelled | expected x3 |

**Conflicts: 0.** This confirms the curator's claim. The merged label equals `expected_label` for all 23 genes under all 3 policies. The curator's `goa_label_check.tsv` agrees with my output in every cell I compared (all 23 genes).

## Report counts

Checked against the files, and correct: 30 rows; 23 genes; rows by label P-ext 23 / N-int 4 / N-sec 3; genes 16 / 4 / 3; IDA 30; selected_by_predictor yes 13 / no 11 / unknown 6; candidates 50 (row_planned 20, skipped_with_reason 28, no_usable_evidence 2); addendum 4; query_log 34; both hashes pass.

Not correct:
- "9 P-ext genes without yes": this becomes 8 after the Pep1 fix.
- The strain caveats for the sirtuins and Xyn2/Xyn11A are wrong (strains are stated).
- The description of 18456523 and 18508396 as proteomics papers is wrong.
- query_log line 34 (36354345 "full text XML") is wrong.
- "Proteome header GN=PIT2 UMAG_01375" is inexact (the header has no UMAG ID).

Not verified: the number of titles screened, and the yield-per-hour figures.

## Final counts

| Count | Genes | n |
|---|---|---|
| P-ext, direct, at least one confirmed row | Rsp3, Stp1, Stp2, Stp3, Stp4, Pep1, Pit2, Erc1, Cmu1, Sts2, Llp1, Xyn2, Xyn11A, Mer1, Tay1, Cpl1 | **16** |
| minus genes whose only rows are selected_by_predictor=yes (Erc1, Mer1, Tay1, Stp1-4, Pep1 after the fix) | Rsp3, Pit2, Cmu1, Sts2, Llp1, Xyn2, Xyn11A, Cpl1 | **8** |
| of the 16: P-ext only through curated rows | the 16 minus Pit2, Erc1, Cmu1 | 13 |
| of the 16: evidence only from strong-promoter rows | Stp2, Sts2, Xyn2, Xyn11A, Cpl1 | 5 |
| of the 8 non-yes genes: native or unstated promoter, not overexpression-only | Rsp3, Llp1, Pit2, Cmu1 (Pit2 and Cmu1 promoter not stated; both already P-ext by GOA EXP) | 4 |

Negatives that survive: N-int Sir2, Hst4, Hst5, Hst6 (4). N-sec Stp5, Stp6, Pdi1 (3), with the owner questions in question 6.

My recommended headline is **16 P-ext genes, 8 without predictor-selected genes**. If the owner drops strong-promoter-only evidence, the counts are 11 and 4. If the owner drops only Cpl1, they are 15 and 7.

Caveat on selected_by_predictor: almost every *U. maydis* effector here comes from the SignalP-predicted secretome at some point. The `no`/`unknown` genes are controls (Pit2, Cmu1), genes found by family search (Xyn2, Xyn11A), genes chosen by expression (Llp1, Cpl1), or genes whose selection the paper does not explain (Rsp3, Sts2). The "without yes" count depends on that field and is fragile.

## What I could not verify

- Supplementary material: Pdi1 S1 Table, Couturier 2012 Additional file 2, Stp Extended Data blots per panel, the Cce1 complementation construct promoter, and Supplementary Table S1 of the sirtuin paper.
- Full texts of Lep1 (33843063), Sta1 (32112567), 18456523 and 18508396 (no open full text). I read their abstracts only.
- Abc1 (39487654): I did not read its localisation passage in depth.
- Whether the first write of `draft_rows.tsv` was after the hash of `candidates.tsv`. Only last-modified times exist.
- I read the passages around the hits, not each paper in full.
