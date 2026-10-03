# H99 curation pilot, pass 2: independent review

Reviewer: Claude Code (claude-opus-5-5), 2026-10-02. Per-row verdicts are in `review_report.tsv`.
All counts come from `review_report.tsv`, `draft_rows.tsv`, `gpi_leads.tsv` and `candidates.tsv`.
I opened no score, prediction or blinded file. I ran no predictor. I did not read `_workdir/step1_compare/phasec/`.

Method:
- I re-fetched every cited paper myself. I did not use the curator's `texts/`, except one grep of `PMC4747415.txt` and `PMC5311391.txt` (see the candidates section).
  - Europe PMC `fullTextXML`: 11 papers.
  - PMC HTML page: 5 papers (PMC5186401, PMC1398056, PMC1828480, PMC98673, PMC13220255).
  - NCBI efetch: PMC13220255.
  - Europe PMC abstract: PMID 17947228.
- A script checked each quote against my copies after normalisation (whitespace, quotes, dashes).
- I checked labels with the repo's own `labels.classify`, `go_obo` and `gaf.filter_gaf`. The inputs were the H99 GOA file `_workdir/step1_compare/downloads/313589.C_neoformans_var_grubii_H99.goa` (taxon 235443) and `go-basic.obo`. I computed each gene's label from GOA alone, from the curated rows alone, and from both merged, under all 3 policies.
- I checked accessions against the `UP000010091.fasta.gz` headers, and checked ORF names with UniProt REST.

## Draft rows (32)

| Outcome | Rows | Count |
|---|---|---|
| confirmed | 2-5, 7-21, 23, 28, 29, 32 | 23 |
| confirmed with a note (`strain_unstated` only) | 6 (CDA1), 25-27 (BIM1), 30 (MNS1), 31 (MNS101) | 6 |
| rejected | 22 (SOD1), 24 (CIG1 PM), 33 (VCX1) | 3 |
| cannot_verify | none | 0 |

Rejected rows:
- **Row 22, SOD1, `wrong_label`.** The row itself (cytosol, IDA, PMC7961099) is sound. However, H99 GOA already has 2 rows for SOD1 from PMID 16524904: GO:0005886 plasma membrane (EXP) and GO:0044853 plasma membrane raft (IDA). The parent rule therefore gives N-sec under all 3 policies, with or without this row. The GO-only label is already N-sec. As written, the row goes to `curated_conflicts.tsv`.
- **Row 24, CIG1 plasma membrane, `claim_not_supported`.** See question 1.
- **Row 33, VCX1, `claim_not_supported`.** The sentence reports vacuole shape. It uses Vcx1-mCherry as an assumed marker. No sentence says that the Vcx1 location was tested, and no co-stain with a vacuole dye is named. The VCX1 label does not change: GOA has GO:0005774 vacuolar membrane IDA (PMID 20889719), so VCX1 is N-sec from GO alone.

Notes on confirmed rows. None of these is a rejection.
- **CDA1 (rows 5 and 6).** The rows add only GO:0016020 (membrane). This term is in ANY_SECRETORY, but it is not in SECRETORY or SURFACE. After the merge, CDA1 is P-ext only under `non_iea`, and only through the ISS (homology) terms for wall and extracellular region. Under `no_homology` and `experimental` it is `unlabelled`. `expected_label` P-ext matches the default label, so no conflict row is written. However, **CDA1 is not a direct (homology_only=no) P-ext gene.**
- **PQP1 (row 8).** GO:0009986 (cell surface) has no ancestor in any `labels.py` set, so the row changes no label. PQP1 is P-ext through the GOA EXP extracellular row (same PMID), which is not homology evidence. The quote does not name Pqp1. The next sentence (residual cleavage in pqp1Δ) supports the attribution.
- **PLB1 raft (row 10).** `selected_by_predictor` should be `no`, not `unknown`. Nothing in the text says a predictor chose PLB1.
- **LAC1 (row 14).** The 2001 antibody is "to the C. neoformans laccase". LAC2 was described later, so the text does not exclude cross-reaction.
- **SOD2 (rows 20 and 21).** Fig. 6A uses the SOD2C strain (sod2Δ complemented with WT SOD2 at the safe haven locus). I accept "the current studies presented here used the fully virulent 'stud' strain" (stud = original H99) as a strain statement, as pass 1 did.
- **HDA-only surface evidence.** APH1, CNAG_05312 and CIG1 (after row 24 is dropped) have only HDA surface evidence, so they fall in `surface_evidence_htp_only`.

Checks that passed for all 32 rows:
- PMIDs: all 17 PMIDs exist in Europe PMC/PubMed. Each concerns the named protein in *C. neoformans*.
- Quotes: 32 of 32 are verbatim in my own re-fetched copies (31 full text, 1 abstract). The gpi_leads quotes are verbatim too (3 of 3).
- GO terms: all 11 terms are known and not obsolete.
- Evidence codes: IDA 29, HDA 3. All are appropriate.
- Accessions: all 20 are in UP000010091 with the right protein.
  - Checked by ORF ID: BIM1 (paper gives J9VHN6), CEL1 (UniProt J9VH79 ORF CNAG_00601), UGG1, MNS1 and MNS101 (CNAG IDs in the GN field), CIG1 (CNAG_01653) and cnap1 (CNAG_05872).
  - Checked by name only: NOP1 (the only fibrillarin entry, which GOA already annotates with IDA nucleus, PMID 24520056) and VCX1.
- `evidence_level`: all rows are `direct` with non-homology codes. This is consistent.
- `selected_by_predictor`: right, except row 10 (should be `no`).
- Format: all 32 `evidence_note` values carry "retrieved 2026-10-02". `reviewer` and `review_date` are empty, as expected before review.

## Answers to the questions

1. **CIG1 Cig1:mCherry PM row (41313167): not usable as direct evidence.** Reject it.
   - The fusion is "C-terminus mCherry-tagged" on a protein that the authors call GPI-anchored. GPI attachment removes the C-terminal signal, so the mCherry signal need not show where mature Cig1 is.
   - The fusion is overexpressed from the ACT1 promoter, and the paper says "In the wild-type parental strain H99, Cig1 is lowly expressed".
   - No membrane marker is named.
   - The imaged strain is "a C. neoformans strain expressing Fbp1:FLAG". That sentence does not name its background.
   - Effect: CIG1 stays P-ext through the HDA row 23. It loses its non-IEA PM term, so it is no longer an `is_pm_candidate` (D8 input).
2. **SOD1: the N-int label cannot stand.** This does not depend on the raft paper's abstract. H99 GOA already carries the raft finding as PM EXP and raft IDA rows, so the derived label is N-sec under every policy. A curated row cannot override GOA (spec 3.4). The options are (a) set `expected_label` to N-sec, or (b) drop SOD1 as a negative. Whether a KCN-sensitive SOD activity in raft fractions is good PM evidence is a GOA question for the owner. Curation cannot change it.
3. **NOP1: keep** (confirmed). The mapping by name is supported, because GOA already puts an IDA nucleus row from another paper on the same accession J9VR32. The fusion is overexpressed, but it is co-stained with DAPI in H99α. **VCX1: drop** (see above). Dropping it costs nothing, because GOA already makes VCX1 N-sec.
4. **BIM1 (unnamed parent strain).** The strain table (Supplementary Dataset 1) is the only place that names the parent of bim1Δ (DTY1000). I did not read it. The text says "a series of isogenic single or combined mutants were generated", and every figure uses "wild-type H99 (DTY758)" as the WT. H99 is probable, but the text does not state it. I count BIM1 as confirmed with a note. A strict "H99 stated" rule would drop it.
5. **UGG1, MNS1, MNS101 as N-sec: correct.** GO:0005783 (ER) has GO:0012505 (endomembrane system) as an ancestor, and GO:0012505 is in SECRETORY. GOA has no surface or internal term for these genes (only IEA ER, ER lumen and membrane). So the label is N-sec under all 3 policies, with GOA merged.
   - UGG1 strain: GFP-Ugg1 was integrated "into the CNAG_03648 locus of the ugg1Δ strain", and ugg1Δ was made in "serotype A strain H99". This makes the H99 background stated.
   - MNS1 and MNS101: the fusions were integrated into "the WT strain", and the paper does not name H99 in that sentence. So these 2 rows get `strain_unstated`.
6. **PMIDs 41313167 and 40424315 exist and are the cited papers.**
   - 41313167 is Avina SL et al. 2026, *Infect Immun*, "Mannoprotein Cig1 contributes to the immunogenicity of a heat-killed F-box protein Fbp1 *Cryptococcus neoformans* vaccine model", PMC12797982.
   - 40424315 is Mota C et al. 2025, *eLife*, "Evolutionary unique N-glycan-dependent protein quality control system ...", PMC12113280.
   - Both quotes are verbatim.

## GOA cross-check (the curator did not do this)

| Gene | GO-only label (non_iea / no_homology / experimental) | Merged with curated rows | Conflict with `expected_label`? |
|---|---|---|---|
| SOD1 | N-sec / N-sec / N-sec | N-sec under all 3 | **yes** (expected N-int) |
| CDA1 | P-ext / unlabelled / unlabelled | same | no under non_iea, but P-ext rests on ISS only |
| APH1, CNAG_05312, CIG1, BIM1, CEL1 | unlabelled | P-ext under all 3 | no. These 5 genes are P-ext only because of curated rows |
| UGG1, MNS1, MNS101 | unlabelled | N-sec under all 3 | no |
| CDA2, QSP1, PQP1, PLB1, LAC1, cnap1 | P-ext under all 3 | same | no |
| CRZ1, SOD2, NOP1 | N-int under all 3 | same | no |
| VCX1 | N-sec under all 3 | same | no |

The rejection of rows 22, 24 and 33 does not change any merged label in this table. SOD1 is N-sec whether or not row 22 is kept.

## gpi_leads.tsv (3 rows)

| Line | Gene | PMID | Verdict |
|---|---|---|---|
| 2 | CDA2 | 22354955 | confirmed (verbatim incl. "in C. neoformans."; KN99; PI-PLC and Triton X-114) |
| 3 | CDA2 | 22354955 | confirmed (verbatim; omega-site deletion) |
| 4 | PLB1 | 17947228 | confirmed (verbatim abstract; H99; PI-PLC release, YW3548; abstract only) |

Format gap: spec section 5 says `curated_gpi.tsv` keeps `source_id`, `gene_id` and `note`, and gains `species`, `uniprot_accession`, `evidence_level`, `evidence_note`, `reviewer`, `review_date` and `override_tm`. `gpi_leads.tsv` has only `symbol`, `uniprot_accession`, `pmid`, `evidence_note` and `note`. These columns must be added before the merge.

## Candidates: checks of the skip reasons

I checked 6 skips that were made on accession grounds. Two reasons are factually wrong.
- **DNJ1: wrong.** The local proteome header has `sp|J9VKM5|DNJ1_CRYN9 ... DNJ1`. This is a unique entry (UniProt ORF CNAG_01347). The toolkit paper reports "ER-mCherry overlapped strongly with ER-Tracker staining" for Dnj1-mCherry. The Dnj1 paper (34566931, PMC8461255) used Dnj1-GFP to confirm ER localisation. A row is possible. It would not change the label, because GOA already has GO:0005788 ER lumen EXP from 34566931, so DNJ1 is N-sec from GO alone. The paper also uses one strain ID (YSB11832, "SH:PH3-DNJ1-mCherry") for both the Dnj1 and the Mjr1 images. This is an error in the paper.
- **CFL1 (28039134): the reason is wrong.** The curator said that CNAG_00795 has no entry. The proteome has `tr|J9VHS9| ... GN=CNAG_00795` ("CPL1-like domain-containing protein"). The skip may still be right for another reason: I did not establish the strain for the secretion assay, and the paper also uses serotype D genomes.
- **CPL1 observation in pass2_report.md: wrong.** J9VHS9 is CNAG_00795, not CNAG_02797. UniProt REST gives J9VQE5 = CPL1 = ORF CNAG_02797. The pass-1 seed accession is therefore consistent with 28039134. No conflict exists.
- **PMA1: right.** There are 2 "Plasma membrane ATPase" entries (J9VWK6 CNAG_06400, J9VVZ0 CNAG_03565), and the text gives no ID. The paper also says "The plasma membrane marker Pma1 also labeled vacuolar structures".
- **ATG8, RAB5: right.** No entry has these gene names. ATG8 matches several unnamed "Autophagy-related protein" entries. The toolkit text has no CNAG IDs except the safe-haven flanks (CNAG_00777/00778).
- **RAC1 (22327008): skip confirmed.** The PMC page says "(serotype D) strain JEC21 was used as the parental strain".

**ECM33/GAS1/PDA1.** I listed all 124 hits of the curator's query 7, not a sample. This includes the 99 that the curator did not screen.
- None of the 124 is an H99 localisation experiment on Ecm33, Gas1 or Pda1. Most are reviews, or papers on *Candida*, *Aspergillus* or plant pathogens.
- Hit 36 (PMID 26878023, PMC4747415) includes *C. neoformans* var. *grubii* H99 cell-wall and secretome MS. It states: "the Ecm33 homologs, which were abundant in the filamentous fungi studied here, were not detected in the fractions from yeast". This is negative for ECM33. The paper is a possible HDA source for other H99 cell-wall proteins (in supplementary tables, not read). It filtered by predicted signal peptides, so any rows from it need `selected_by_predictor` = yes.
- The `no_usable_evidence` status for these 3 genes is acceptable for query 7. The other queries (8: 442, 9: 119 unscreened hits) stay unscreened.

## pass2_report.md counts

Correct when I recount them from the files:
- 32 rows and 20 genes.
- Rows by label: P-ext 23, N-int 5, N-sec 4. Genes by label: 12, 4, 4.
- IDA 29, HDA 3.
- `selected_by_predictor`: no 28, yes 3, unknown 1.
- 43 candidates: row_written 20, no_usable_evidence 7, skipped_with_reason 15, not_screened 1.
- 18 carried rows (12 carried, 6 corrected) and 14 new rows. 9 new genes.
- `sha256sum -c candidate_list.sha256` passes.

Not correct or not verifiable:
- **"12 distinct P-ext direct genes" is too high.** See the final counts below.
- The DNJ1, CFL1 and CPL1 statements are wrong (see above).
- "Papers screened: 22" and the title-screen totals: I did not verify these.
- Process: `h99_pass2/` is untracked in git. Spec section 6 requires `candidate_list.sha256` to be committed before any join to scores. The `candidates.tsv` statuses record outcomes (`row_written`), so the hash fixes a list that was written after the curation was done. This is acceptable only because no score join has happened.

## Final surviving counts

The rule from the brief: distinct genes with `expected_label` P-ext, `evidence_level` direct, and at least one row with verdict confirmed or confirmed with a note only.

| Count | Genes | n |
|---|---|---|
| By the brief's rule | APH1, BIM1, CDA1, CDA2, CEL1, CIG1, CNAG_05312, LAC1, PLB1, PQP1, QSP1, cnap1 | **12** |
| Also requires the gene's P-ext label to rest on non-homology evidence after the GOA merge (spec 3.3 headline: homology_only=no) | the 12 above minus CDA1 | **11** |
| Also requires H99 to be stated, not probable (drops BIM1) | 11 minus BIM1 | **10** |
| Of the 11: surface evidence only HDA | APH1, CNAG_05312, CIG1 | 3 |
| Of the 11: P-ext only because of curated rows (GO-only label unlabelled) | APH1, CNAG_05312, CIG1, BIM1, CEL1 | 5 |

Negatives that survive:
- N-int: CRZ1, SOD2, NOP1 (3 genes).
- N-sec: UGG1, MNS1, MNS101 (3 genes).
- SOD1 derives N-sec, but its row was rejected for the wrong label. VCX1 is N-sec from GO, but its row was rejected.

My recommended headline count is **11**. It is 10 if the owner requires H99 to be named in the text for the tagged strain.

## What I could not verify

- Supplementary tables: BIM1 Dataset 1 (strain parents), the toolkit paper's marker CNAG IDs, and the 40424315 strain list (Supplementary file 1A).
- Whether the GOA SOD1 PM and raft rows are good evidence. They come from KCN-sensitive activity in raft fractions.
- The total of titles screened in pass 2 for queries other than query 7.
- The full text of PMID 17947228 (no PMC copy exists).
- Time: about 45 tool calls (17 Europe PMC metadata, 16 full-text fetches, about 10 UniProt REST calls, plus local GOA/OBO/FASTA checks). I did not measure wall-clock time.
