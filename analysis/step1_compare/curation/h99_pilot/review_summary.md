# H99 curation pilot: independent review

Reviewer: Claude Code (claude-opus-5-5), 2026-10-02. Per-row verdicts are in `review_report.tsv`.
All counts below come from `review_report.tsv`, `draft_rows.tsv` and `gpi_leads.tsv`.
I opened no blinded file and ran no predictor.

## Draft rows (19)

| Outcome | Rows | Count |
|---|---|---|
| confirmed | 2, 3, 4, 6, 7, 8, 10, 11, 15, 16, 17 | 11 |
| confirmed with a note (`strain_unstated` only) | 9 (PLB1 raft), 12 (LAC1) | 2 |
| rejected | 5, 13, 14, 18, 19, 20 | 6 |
| cannot_verify | none | 0 |

Rejected rows by reason:
- `wrong_label`: 4 (row 5 CDA1, rows 13 and 14 APH1, row 20 actin).
- `wrong_code`: 2 (rows 18 and 19 SOD2: the assay is fractionation with markers, so IDA, not EXP).
- No `wrong_pmid`, `quote_not_found`, `claim_not_supported`, `strain_not_h99`, `wrong_term` or `wrong_accession`.

Text source: 16 rows were checked against full text and 3 against the abstract only (rows 8, 9, 12). Full text is not
available for these 3 papers from the PubMed tool, Europe PMC or NCBI efetch.

Other checks:
- PMIDs: all 11 PMIDs exist and are about the named protein in *C. neoformans*.
- GO terms: all 8 terms exist in `go-basic.obo` as cellular-component terms and are not obsolete.
- Accessions: all 12 accessions are in UP000010091 with the right protein.
  - cnap1 = J9VS02. I confirmed this by ID, not only by name. The paper lists "CNAG_05872 ... MAY1". UniProt J9VS02 has ORF CNAG_05872 and the synonym MAY1, which cites PMID 27977806.
  - Actin: the file holds P48465 ("Actin", GN=CNAG_00483). The curator's summary ("actin P48465") and the file name the same entry. The paper gives no gene ID, so this mapping is by protein name.
- Strain: the curator's notes say several strains were "not stated". The full texts resolve 4 of them to H99:
  - CNAG_05312: "The C. neoformans var. grubii wild-type strain H99 (WT) ... were used for this study".
  - CRZ1: "The strains Crz1-GFP, Δcrz1, Crz1-GFP:Pab-dsRed were created from H99".
  - SOD2: "the current studies presented here used the fully virulent 'stud' strain".
  - PLB1-GFP: "Wild-type C. neoformans var. grubii strain H99 ... were used in this study".
  - CDA1, CDA2, QSP1 and PQP1 use KN99 or KN99alpha. These are H99-congenic serotype A strains, so I accepted them.
- `selected_by_predictor`: 2 rows are wrong. They are not counted as rejections, because the verdict list has no code for this.
  - Row 16 (May1): the authors chose deletion candidates among the MS-identified "peptidases with predicted signal sequences" (SignalP 4.0). The value should be `yes`.
  - Row 13 (APH1): the paper says "Only one of these enzymes, Aph1, is predicted to be secreted and was identified in our proteomic analysis". The value should be `yes`, or `unknown` at most. The curator's own note says this, but the field says `no`.
- Format: the spec (section 4) requires a retrieval date in `evidence_note`. No row has one.

## gpi_leads.tsv (4 lines)

| Line | Gene | PMID | Verdict |
|---|---|---|---|
| 2 | CDA2 | 22354955 | confirmed. The quote ends "Cda2 in." The paper says "Cda2 in C. neoformans." (the italic species name was dropped) |
| 3 | CDA2 | 22354955 | confirmed (exact; omega-site deletion, KN99) |
| 4 | PLB1 | 17947228 | confirmed (exact abstract; H99, PI-PLC release, YW3548) |
| 5 | PLB1 | 15826239 | claim_not_supported. The anchor-motif experiment is in *S. cerevisiae* (pYES2). The quote also drops "(GPI-)": the abstract reads "create PLB1(GPI-) resulted". The only *C. neoformans* data in the abstract is beta-glucanase release from "wild-type C. neoformans" (strain unstated). That shows wall location, not the anchor. |

Totals: 3 confirmed, 1 rejected, 0 cannot_verify.

## Answers to the curator's questions

- **APH1:** P-ext is right, and `ambiguous` is wrong.
  - In `analysis/step1_compare/labels.py`, vacuole (GO:0005773) is in SECRETORY, not in INTERNAL. INTERNAL holds only cytosol, nucleus and mitochondrion.
  - Extracellular region plus vacuole therefore gives P-ext.
  - The brief's rule leads to the same answer, because a vacuole is not inside evidence.
- **CDA1:** this is the weakest row. The quote shows a crude "cell membrane fraction" (Fig. 5 legend: "M, cell membrane fraction"), mainly for the catalytic mutant. The sentence that the curator cut off names a CW lane, but the text does not say whether Cda1 is in it.
  - For the row alone, N-sec is right.
  - The gene does not get N-sec. H99 GOA already has non-IEA ISS terms for CDA1: fungal-type cell wall and extracellular region. The parent rule gives P-ext under the `non_iea` policy, and `unlabelled` under `no_homology`. So `expected_label` should be P-ext, or the row will go to `curated_conflicts.tsv`.
  - The paper does not separate plasma membrane from internal membranes. GO:0016020 (membrane) would be the more faithful term.
- **CDA2 (rows 2-4):** P-ext is right for each row as a gene-level label.
  - The paper sees Cda2 in the cytosolic fraction. The authors attribute this to shearing: "The appearance of Cda2 in this fraction is likely a consequence of shearing from the wall during sample preparation". I agree that no cytosol row should be written.
  - All 3 terms are already in H99 GOA with code IDA, so these rows add a PMID and a quote but no new term.
- **Actin (row 20):** N-int fails under the parent rule. H99 GOA has an IEA term, actin cortical patch (GO:0030479). Its is_a/part_of ancestors include cell periphery (GO:0071944), which is in ANY_SECRETORY. The derived label is therefore `unlabelled`. A cytosol-only negative control needs a protein with no IEA periphery or membrane terms.
- **PQP1 (row 7):** I confirm the row. However, GO:0009986 (cell surface) has no wall or extracellular ancestor. The row therefore adds no surface term, and PQP1 is P-ext only through the GOA extracellular-region EXP row.

## Skipped and no_usable_evidence candidates: possible misjudgements

1. **CIG1: too strict.** The gene name lookup failed. However, the proteome has J9VUB8, "Cytokine inducing-glycoprotein", GN=CNAG_01653. The 26453029 secretome table lists "CNAG_01653 Cytokine-inducing glycoprotein" in strain H99. That is the same evidence standard as the CNAG_05312 HDA row.
2. **SOD1: too strict.** The curator looked only at the abstract of 33567338. The full text (PMC7961099, H99) has H99 fractionation with markers: "High levels of cytosolic Sod1 were detected during Cu sufficiency, with a small fraction associated with the mitochondria". That supports a cytosol IDA row for J9VLJ9. The raft-membrane finding (16524904, abstract only) would add a membrane term and block N-int, so the owner should decide on that row.
3. **RAC1: too strict, and inconsistent with the curator's own practice.** The abstract of 22327008 says Rac1 "depended on Wsp1 for its vacuolar membrane localization". That could be an N-sec row. The curator wrote rows for LAC1 and CRZ1 on unstated-strain abstracts, but not for this one. However, UniProt has no RAC1 gene name in taxon 235443, so the CNAG ID must come from the full text (PMC3318296, which Europe PMC did not return).
4. **Cap64: inconsistent.** The curator skipped it for "strain not stated", which is the same state as LAC1. The skip is still right for other reasons. The data is an overexpressed mCherry fusion in "patch-like" acidic compartments. The accession is also unresolved: the proteome has only "CAP64 gene product-related" entries.
5. **ECM33/GAS1/PDA1: too loose.** The search returned 619 hits and none were screened. `no_usable_evidence` is the wrong status, because nothing was evaluated. Use a status such as `not_screened`.

The other skips are right:
- CPL1: only var. neoformans data.
- YOR1: the paper says localisation failed.
- KRP1: the data is from *C. gattii*.
- MPK1: the Cda2 paper gives it as "data not shown". A MAP kinase CNAG_04514 is J9VSI6, but the paper gives no ID.
- URE1 and LAC2: no H99 location statement in the text that I could read.

## Quality of quote handling

- Draft rows: 17 of 19 quotes matched exactly, allowing whitespace and typography only.
  - 15 matched by script against saved full text or abstracts.
  - Rows 6 and 7 matched by eye against the PubMed-tool full text, because the tool output was not saved to a file.
- 2 quotes needed correction. Both come from formatted text that the PubMed tool drops:
  - Row 5: "(Cda1)" should be "(Cda1CS)". This changes which protein the sentence describes (the catalytic mutant). The sentence was also cut before "as can be seen on lanes marked M for membrane and CW for cell wall fractions in Fig. 5".
  - Row 13: "strain of, H99" should be "strain of C. neoformans, H99".
- gpi_leads: 2 of 4 quotes were exact and 2 were corrected (lines 2 and 5).
- I did not reject any row for these formatting losses. However, the row 5 loss changes the meaning. A curator should re-copy quotes from Europe PMC XML when it exists.
- The curator said 7 rows were not machine-checked. All 7 match: rows 6, 7, 8 and 12 by my check, and rows 17, 18 and 19 against both the abstract and the full text.

## Time and tool problems

- PubMed `get_full_text_article` returned an empty `full_text` for PMC1828480 (LAC1) and PMC1398056 (PLB1 raft). Europe PMC `fullTextXML` returned a 150-byte "not available" reply for PMC5186401, PMC1398056, PMC1828480 and PMC1180731. NCBI efetch returned no body for these, because the publisher does not allow XML download.
- PMID 17947228 has no PMC copy. The curator's log maps PMC2224146 to it, but that PMC ID belongs to PMID 18039940 (the vesicle paper).
- The full text for CRZ1 (PMC3520850) and SOD2 (PMC7961099) was available from Europe PMC. The curator used only the abstracts, and so missed the strain statements and the fractionation assay for SOD2.
- Review cost: about 25 tool calls (one PubMed metadata batch, 3 PubMed full-text calls, 11 Europe PMC fetches, 4 NCBI efetch calls, a few UniProt REST calls, and local GOA/OBO/FASTA checks). I did not measure wall-clock time.
