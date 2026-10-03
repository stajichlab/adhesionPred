# Ustilago maydis curation pilot: report

Curator: Claude Code (claude-sonnet-5-5), 2026-10-02. Counts below were computed from the TSV files in this directory (scripts: `build_rows.py`, `goa_check.py`, and the inline count script whose output I copied).
I opened no score, prediction or blinded file. I did not edit any existing file. I ran no `git commit`.
Reference proteome: `UP000000561.fasta.gz` (Umay_MYCMD, taxon 5270). GOA file: `MYCMD-uniprot.gaf.gz`. `gene_id` equals the UniProt accession, as spec section 5 says for this source.

## Files
- `candidates.tsv` (50 entries) and `candidate_list.sha256`. The hash was written before `draft_rows.tsv`. `sha256sum -c` passes.
- `candidates_addendum.tsv` (4 entries) and `candidate_addendum.sha256`. I added 4 genes (Mer1, Tay1, Stp4 second paper, Cpl1) after the first hash, because I found them while reading papers that were already fetched. They are in a second hashed file so that `candidates.tsv` stays unchanged.
- `draft_rows.tsv` (30 rows), `row_provenance.tsv`, `gpi_leads.tsv` (header only, 0 rows), `goa_label_check.tsv`, `query_log.tsv` (34 entries), `texts/` (XML, PMC HTML text and abstract files).
- `candidates.tsv` has 6 columns. It has no `origin` column (the H99 file has one).
- Helper scripts: `fetch.py`, `qc.py`, `es.py`, `up.py`, `pmchtml.py`, `ql.py`, `build_rows.py`, `goa_check.py` (run with `/usr/bin/python3.12`, because the repo code needs Python 3.10 or newer).

## Counts (draft_rows.tsv)
- Rows: 30. Distinct genes (UniProt accession): 23.
- Rows by expected_label: P-ext 23, N-int 4, N-sec 3. Genes by label: P-ext 16, N-int 4, N-sec 3.
- Rows by evidence_code: IDA 30. There are no HDA rows and no `transfer` rows. All 30 rows are `direct`.
- Text source: 30 rows from full text, 0 rows from an abstract only. All 12 cited papers had Europe PMC XML, so no row rests on PMC HTML text.
- selected_by_predictor: yes 13 rows, no 11 rows, unknown 6 rows. Genes where every row is `yes`: Erc1, Mer1, Stp1, Stp2, Stp3, Stp4, Tay1 (7 P-ext genes). Reports that drop `yes` genes keep 9 P-ext genes (Cmu1, Cpl1, Llp1, Pep1, Pit2, Rsp3, Sts2, Xyn11A, Xyn2).
- Machine quote check: 30 of 30 rows pass. I re-read `draft_rows.tsv`, split each `evidence_note` into its quoted segments, and found every segment in the saved text (whitespace normalised). Quotes were cut from the saved Europe PMC XML text by anchor strings, so they are verbatim. A row has 1 to 3 segments joined by " ... ".
- All 30 `evidence_note` values contain "retrieved 2026-10-02". `reviewer` and `review_date` are empty.
- Accessions: 23 of 23 are in the UP000000561 FASTA. Check method: UniProt REST ordered-locus lookup by the UMAG ID that the paper states, then membership in the FASTA. Two exceptions, mapped by exact gene name against the FASTA header: Pit2 (GN=PIT2, UMAG_01375) and Cmu1 (GN=CMU1, UMAG_05731). The papers used for those two rows do not give the UMAG ID.
- Distinct genes with expected_label P-ext and evidence_level direct: **16** (Rsp3, Stp1, Stp2, Stp3, Stp4, Pep1, Pit2, Erc1, Cmu1, Sts2, Llp1, Xyn2, Xyn11A, Mer1, Tay1, Cpl1).
- Genes whose label comes only from the curated rows (GOA alone does not give the expected label), from `goa_label_check.tsv`:
  - `non_iea` policy: 16 genes. 13 are P-ext (Rsp3, Stp1, Stp2, Stp3, Stp4, Pep1, Sts2, Llp1, Xyn2, Xyn11A, Mer1, Tay1, Cpl1). 3 are N-sec or N-int (Stp5, Stp6, Hst5).
  - `no_homology` and `experimental` policies: 20 genes. 13 are P-ext (same list). 7 are negatives (Stp5, Stp6, Pdi1, Sir2, Hst4, Hst5, Hst6). Pdi1, Sir2, Hst4 and Hst6 have only IBA or IEA rows in GOA, so the GOA label under these policies is `unlabelled`.
  - Pit2, Erc1 and Cmu1 are already P-ext in GOA under all 3 policies.
- Conflicts with GOA: **0 genes**. I merged GOA and curated rows and ran `labels.classify` under all 3 policies. The merged label equals `expected_label` for all 23 genes. All GO terms are known and not obsolete.
  - Observation, not a conflict: GOA has Erc1 as IDA GO:0030446 (hyphal cell wall) and GO:0044158. My Erc1 rows claim the extracellular region only. The paper says TEM showed no specific Erc1 accumulation in the fungal wall.
  - Observation: GOA has Hst6 as IBA nucleus, mitochondrion and cytosol. The label is N-int with or without my row.

## Papers
- Papers with at least one row: 12 (PMIDs 19197359, 29703884, 29745456, 31730668, 33941900, 34166468, 34947062, 36224193, 37113216, 37171083, 37872143, 39997458).
- Rows and genes per paper: 33941900 (Stp paper) 9 rows, 7 genes (Stp1 to Stp6 and the Pit2 control). 37113216 (sirtuins) 4 rows, 4 genes. Stp1 to Stp6 and the 4 sirtuins are 10 of the 23 genes. All other papers gave 1 or 2 rows each.
- Papers screened: 52 PMIDs have a saved abstract or text file in `texts/`. Full text was fetched for 41 of them (38 XML, 3 PMC HTML). I searched the full texts with regular expressions for localisation and secretion sentences. I read in full only the passages around the hits. I did not read the whole of any paper.
- Abstract only (no full text found): 17042749, 27563844, 20663961, 16314447, 17405809, 18456523, 18508396, 33843063, 38742361, 32112567, 36354345, 20587773, 30510169. None gave a row.
- Fetched but not searched in depth: 35083222, 33802393, 37590419, 42429560, 40141077, 33653886, 22589719.
- Titles screened in Europe PMC lists: about 12 to 40 per query for 20 queries (`query_log.tsv`). I did not count unique titles.

## Time and yield
- Start 22:38 and end about 23:25 (system clock), so about 47 minutes. I stopped before the 90-minute limit because the later searches returned few new papers that I could read as full text with a direct localisation sentence. I have no data on how many more U. maydis rows exist.
- U. maydis: 30 rows, 23 genes, 16 P-ext genes from 12 papers in about 47 minutes. Rows per paper 2.5. Genes per paper 1.9. About 29 genes per hour and 20 P-ext genes per hour.
- H99 pass 2 (from `pass2_report.md` and `review_summary.md`): 32 rows, 20 genes, 12 P-ext genes (11 after review) from 17 PMIDs, in about 85 minutes, after a pass 1 whose time I do not know. Rows per PMID 1.9. Genes per PMID 1.2. About 14 genes per hour and 8 P-ext genes per hour, counting pass 2 time only.
- The comparison is not like for like. Reasons:
  - H99 pass 2 carried 18 rows from pass 1, so its hours understate the total effort.
  - U. maydis papers in this pilot come from effector biology. Two papers gave 10 genes. Without them the yield is 13 genes from 10 papers (1.3 per paper).
  - My U. maydis rows have had no independent review. The H99 review rejected 3 of 32 rows and found 2 counts too high. I expect some of my rows to fail review.
  - U. maydis has a larger set of papers that use HA-tagged secretion assays and mCherry plasmolysis. This is why the P-ext yield is higher. It is not a measure of curator speed.

## Row caveats the reviewer should check
- Strain not named for the tagged strain in the text I read (`strain_unstated`): Sir2, Hst4, Hst5, Hst6 (4 rows; deletions were made in SG200 and the controls are called wild-type). Xyn2 and Xyn11A colony-secretion rows (SG200 is named only for the control in the legend; the tagged strains are not named). Stp1 to Stp4 supernatant rows (the strain is AB33, SG200Δkex2 or SG200 depending on the panel; I did not map panel to protein).
- Overexpression or strong promoter: Rsp3 supernatant row, Stp1 to Stp4 supernatant rows, Stp5 and Stp6 (otef), Cpl1 (otef), Sts2 (pit2 promoter), Xyn11A apoplast row (pit2 promoter). None of these proteins is described as GPI anchored in the text I read. Llp1 GPI status is not stated in the text I read.
- Pdi1: C-terminal GFP on a protein with an HDEL motif. The tag may mask the motif. The signal still co-localised with the ER marker mRFP:HDEL.
- Stp5 and Stp6 (N-sec): plasma membrane fraction with Pit1 as marker. The same paper shows Stp6-mCherry in speckles outside the hyphae. I did not count that for Stp6.
- Cmu1 and Pit2: accession by exact gene name. Both are already P-ext in GOA.
- Wall term for Rsp3 rests on the paper's figure title ("Secreted Rsp3-HA binds to the U. maydis cell wall.") and the sentence "detected around the outside of fungal hyphae". Llp1 and Cpl1 use sentences that state the cell wall.
- selected_by_predictor `yes` rests on a sentence in the paper that says the genes were predicted effectors or SignalP-predicted (Stp family, Erc1, Pleiades). `unknown` was used when I did not find such a sentence. I did not read each introduction in full.

## GPI leads
- `gpi_leads.tsv` has a header and 0 rows. The header has the `curated_gpi.tsv` columns from spec section 5 plus `override_tm`.
- I found no U. maydis experiment on a GPI anchor (PI-PLC release, Triton X-114 partition, omega-site mutation). The CDA papers (36946727, 33653886) call 5 CDAs "putatively" or "predicted" GPI anchored. The PR-1 paper (37716995) names GPI-anchored yeast proteins. I did not make a lead from predictions (spec section 4: predictors are not evidence).

## Skipped candidates and why (28 in `candidates.tsv`)
- Surface or periphery only, with no wall or extracellular GO term, so the label cannot be derived: Cda1, Cda2, Cda4, Cda7 (36946727), PR-1La (37716995), Abc1 (39487654).
- Reporter or activity evidence only: Cts1 (21808052, 32733418, 40348193). The Rrm4 paper says wall association "probably".
- Abstract only, no full text: Chs5, Chs6, Chs7, Mcs1, Chs3/4/6/8 antibodies (16314447, 17042749, 20663961, 27563844), Lep1 (33843063), Sta1 (32112567), Msb2 and Sho1 (20587773), Tin2 (30510169).
- No UMAG ID in the paper, so the accession cannot be verified: Vp1 (34436129), Lip3 (42093215).
- Location is cited from another paper, not tested: Jps1, Pit1, Mmf1, Pdi1 wall statement in the Cpl1 paper (37171083, which cites the Pdi1 paper; I found no wall claim in the Pdi1 paper itself).
- Hedged statement: ApB73 (27279632, "seems to stay attached"). I did not write a row.
- Label not derivable: peroxisome, lipid droplet, endosome proteins, Xyn3 (not secreted, no cytosol claim). Peroxisome, lipid droplet and endosome are not in `INTERNAL` or `SECRETORY` in `labels.py`.
- Controls used only as lysis controls (actin, tubulin): the pellet is not the cytosol, and the mapping would be by protein name.
- Xyn1: band size in apoplast western is unexpected; secretion is cited from an earlier paper.
- Hst2: nuclear only in mitotic cells.
- No usable evidence found: Hum2, Hum3, Rep1 (no hydrophobin localisation paper among 20 titles); the seed list members Mig1, Mig2, Um05505, Rbf1, Mfa1, Pra1, Gas1, Ecm33, Ssr1, Cel1, Cbh, Egl1, Pel1, Sod1 were not searched one by one.
- Supplementary tables not read: the Pdi1 DIGE paper lists 6 cell wall and 11 secreted glycoproteins in S1 Table; the Pleiades paper and others may have more.

## What I could not verify
- No independent review of any row has been done.
- GOA check uses the repo code (`labels.classify`, `go_obo.parse_obo`, `gaf.filter_gaf` with taxon 5270). I did not check whether the 5270 taxon filter drops strain-level GOA rows. The filter dropped 45 `other_db` rows and 1 NOT row and no taxon rows.
- Strain of the tagged strains for the 4 sirtuin rows and the Xyn rows (see caveats).
- Whether the Stp1 to Stp4 supernatant panels use SG200, AB33 or SG200Δkex2 for each protein.
- Whether Llp1 is GPI anchored.
- Whole-paper reading was not done. A paper with a relevant sentence phrased differently could have been missed.
- Unread: supplementary tables (strain lists, secretome lists) and the PDF figures. All claims rest on text and legends.
- The U. maydis secretome and wall proteomics papers (18456523, 18508396, 36354345) have no full text in Europe PMC. I read only their abstracts. They could give HDA rows.

## Tool problems
- PubMed `search_articles` returned only PMIDs. I used Europe PMC (curl) for titles and full text.
- Gene-name PubMed searches returned 0 or few hits for several seed genes. I did not retry them with other names.
- `fetch.py` for 8 PMIDs took longer than 120 seconds once and ran in the background. It finished without error. PMC3237082, PMC99612 and PMC305698 had no Europe PMC XML. I got PMC HTML text with a browser user agent (`pmchtml.py`) and read only PMC3237082 (Cts1).
- Writing to `/scratch/jstajich` directly was refused. I wrote only inside `umay_pilot/`.
