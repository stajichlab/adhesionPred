# 05 Pfam and CAZy family ID verification

Date: 2026-10-17. Scope: every PF\d{5} in reports 01 to 04 in this directory plus `data/sorting_hat/family_table.tsv` (all 24 rows, taken in full; every row is an adhesion, wall, hydrophobin or allergen row; the two allergen rows PF16541 and PF25312 are included). 53 distinct Pfam IDs.

## Method and source

- Pfam source: InterPro API `https://www.ebi.ac.uk/interpro/api/entry/pfam/PFxxxxx` (JSON, HTTP 200 for all 53; fetched with curl on 2026-10-17). Fields used: `metadata.name.short`, `metadata.name.name`, `metadata.type`, `metadata.description`. No fallback endpoint was needed.
- Every ID returned a live Pfam entry. The API record carries no obsolete/merged flag, so "no deprecation seen" means only that the entry resolves. Clan membership was not checked.
- CAZy source: `https://www.cazy.org/GHnn.html` (WebFetch). CAZy family pages list activities only; none of them mention GPI anchoring or fungal protein names (Gas1, Gel, Phr, Dfg5, Dcw1, Crh, Bgl2, Utr2). So "fungal GPI-anchored members described" is NO for every family below.
- Names are marked MATCH only when the claimed name equals the fetched `short` or `name`. Claims about protein examples (for example "Rbt5", "Dfg5") are not in InterPro records and are marked separately.
- Report 04 line numbers refer to `04-repo-coverage-and-gaps.md`; FT = `family_table.tsv` row number (header is line 1).

## Chitin synthase question (the brief versus the repo)

- PF01644 = **Chitin_synth_1** ("Chitin synthase-like glycosyltransferase domain", type family). Record text: GT domain of chitin synthases classes I, II and III; "also describes the GT domain found in hyaluronan synthase found in animals".
- PF03142 = **Chitin_synth_2** ("Chitin synthase", type family). Record text: fungal chitin synthase enzymes, EC 2.4.1.16.
- PF08407 = Chitin_synth_1N, "N-terminal domain of Chitin synthase (PF01644)".
- Verdict: the **repo is right** (04:14, `docs/reports/2026-10-04-cell-wall-gene-classes-vs-tools.md:56`). The orchestrator brief, which gave PF03142 = Chitin_synth_1 and PF01644 = Chitin_synth_2, was wrong (swapped).
- Source: https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01644 and .../PF03142 and .../PF08407.

## Table A: Pfam IDs

| ID | Claimed name/function (report, file:line) | Verified name | Type | Result | Source |
|---|---|---|---|---|---|
| PF01644 | Chitin_synth_1 (04:14, 04:57, 04:126); Division 1 domain (01:40) | Chitin_synth_1, Chitin synthase-like glycosyltransferase domain | family | MATCH (repo names). Brief's "Chitin_synth_2" MISMATCH | [PF01644](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01644) |
| PF03142 | Chitin_synth_2 (04:14, 04:126); "CHS catalytic domain" (01:40, 01:80) | Chitin_synth_2, Chitin synthase | family | MATCH on name. Brief's "Chitin_synth_1" MISMATCH. 01 calls it a "domain" (paper's wording); Pfam type is family | [PF03142](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF03142) |
| PF08407 | Chitin_synth_1N (04:14, 01:40) | Chitin_synth_1N, Chitin synthase N-terminal | family | MATCH | [PF08407](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF08407) |
| PF02364 | Glucan_synthase, beta-1,3-glucan synthase (04:44, 04:58, 04:127) | Glucan_synthase, 1,3-beta-glucan synthase component | family | MATCH | [PF02364](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF02364) |
| PF03198 | Glyco_hydro_72, GH72 GEL/GAS/PHR (04:44, 04:54, 04:119; "to verify" 04:230) | Glyco_hydro_72, Glucanosyltransferase | domain | MATCH. Description: GPI-anchored beta(1-3)glucanosyltransferases | [PF03198](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF03198) |
| PF00722 | Glyco_hydro_16, GH16 Crh/Utr2 (04:54, 04:120; "to verify" 04:230) | Glyco_hydro_16, Glycosyl hydrolases family 16 | domain | MATCH (name). Crh/Utr2 link not in record | [PF00722](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00722) |
| PF00332 | GH17, Bgl2 (04:70, 04:122; "name unverified" 04:231) | Glyco_hydro_17, Glycosyl hydrolases family 17 | domain | MATCH (name). Bgl2 link not in record | [PF00332](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00332) |
| PF03663 | GH76 Dfg5/Dcw1 (04:70, 04:121; "my recall" 04:231) | Glyco_hydro_76, Glycosyl hydrolase family 76; description "Family of alpha-1,6-mannanases" | repeat | MATCH on family name. Dfg5/Dcw1 link NOT in record. Pfam type is "repeat" | [PF03663](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF03663) |
| PF00704 | Glyco_hydro_18, GH18 chitinases (04:44, 04:56, 04:69, 04:124) | Glyco_hydro_18, Glycosyl hydrolases family 18 | domain | MATCH | [PF00704](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00704) |
| PF00150 | Cellulase (GH5), Exg1 (04:56, 04:123) | Cellulase, Cellulase (glycosyl hydrolase family 5) | domain | MATCH (name). Exg1 link not in record | [PF00150](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00150) |
| PF01522 | Polysacc_deac_1, chitin deacetylase CE4 (04:56, 04:70, 04:125) | Polysacc_deac_1, Polysaccharide deacetylase | family | MATCH on Pfam name. CE4 mapping not fetched: UNVERIFIABLE | [PF01522](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01522) |
| PF00026 | Asp, generic yapsin/Sap family (04:55, 04:135) | Asp, Eukaryotic aspartyl protease | domain | MATCH | [PF00026](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00026) |
| PF10342 | Kre9_KNH, Kre9/Knh1 (04:59, 04:130, 04:71) | Kre9_KNH, Kre9/KNH-like N-terminal Ig-like domain | domain | MATCH | [PF10342](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF10342) |
| PF03856 | SUN (Sun4, Sim1, Uth1, Nca3) (04:71, 04:131) | SUN, Beta-glucosidase (SUN family); lists Nca3, Sun4, Sim1 | family | MATCH (Uth1 not in record) | [PF03856](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF03856) |
| PF12454 | Ecm33 / Sps2 (04:71, 04:132) | Ecm33, GPI-anchored cell wall organization protein | family | MATCH (Sps2 not in record) | [PF12454](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF12454) |
| PF10287 | Tos1 (04:71, 04:133) | YJL171C_Tos1_C, Cell wall protein YJL171C/Tos1, C-terminal | domain | MATCH | [PF10287](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF10287) |
| PF10290 | Tos1 (04:71, 04:133) | YJL171C_Tos1_N, N-terminal | domain | MATCH | [PF10290](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF10290) |
| PF07983 | X8 domain (04:70, 04:71, 04:119) | X8, X8 domain | domain | MATCH | [PF07983](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF07983) |
| PF06423 | Gwt1 (04:60, 04:138) | GWT1 | family | MATCH | [PF06423](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF06423) |
| PF05024 | Gpi1 (04:60, 04:138) | Gpi1, N-acetylglucosaminyltransferase component (Gpi1) | family | MATCH | [PF05024](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF05024) |
| PF04113 | Gpi16 (04:60, 04:138) | Gpi16, GPI transamidase component | domain | MATCH | [PF04113](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF04113) |
| PF04114 | Gaa1 (04:60, 04:138) | Gaa1, Gaa1-like GPI transamidase component | family | MATCH | [PF04114](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF04114) |
| PF06432 | GPI2 (04:60, 04:138) | GPI2, Phosphatidylinositol N-acetylglucosaminyltransferase | family | MATCH | [PF06432](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF06432) |
| PF10333 | "Pig-X annotation" (04:138) | Pga1, GPI-Mannosyltransferase II co-activator | family | Record name is Pga1; the Pig-X equivalence is not in the record: UNVERIFIABLE for the Pig-X claim | [PF10333](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF10333) |
| PF01822 | WSC, wall integrity sensors (04:61, 04:136) | WSC, WSC domain | domain | MATCH | [PF01822](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01822) |
| PF00187 | Chitin_bind_1, CBM18 (04:68, 04:72, 04:139) | Chitin_bind_1, Chitin recognition protein | domain | MATCH on Pfam name. CBM18 mapping: CAZy CBM18 page says CBM18 binds chitin but gives no Pfam link, no fungal rows: UNVERIFIABLE | [PF00187](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00187) |
| PF01476 | LysM (04:70, 04:72, 04:139) | LysM, LysM domain | domain | MATCH | [PF01476](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01476) |
| PF01607 | CBM_14, Avr4 (04:72, 04:139) | CBM_14, Chitin binding Peritrophin-A domain | domain | MATCH on name. CAZy CBM14 page says chitin binding; no Avr4 mention. Avr4 link UNVERIFIABLE | [PF01607](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01607) |
| PF07938 | `PF07938_LABEL` code identifier only (04:197), no name claimed | Fungal_lectin, Fungal fucose-specific lectin | repeat | MATCH (no claim made; recorded for completeness) | [PF07938](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF07938) |
| PF05730 | CFEM; "GPI wall protein family (Rbt5, Pga10, Ag2/PRA)" (04:26, FT:2) | CFEM, CFEM domain | domain | MATCH on name. Rbt5/Pga10/Ag2 examples not in the record we read | [PF05730](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF05730) |
| PF04681 | Bys1, yeast-phase surface protein (CalA) (04:30, FT:3) | Bys1, Blastomyces yeast-phase-specific protein | domain | MATCH (name). CalA link not in record | [PF04681](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF04681) |
| PF01185 | Hydrophobin (04:28, FT:4) | Hydrophobin, Fungal hydrophobin | family | MATCH | [PF01185](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF01185) |
| PF06766 | Hydrophobin_2 (04:28, FT:5) | Hydrophobin_2, Fungal hydrophobin | family | MATCH | [PF06766](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF06766) |
| PF28987 | DewD (04:28, FT:6) | DewD, Unclassified hydrophobin dewD | family | MATCH | [PF28987](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF28987) |
| PF22354 | Eas (04:28, FT:7) | Eas, Hydrophobin | family | MATCH | [PF22354](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF22354) |
| PF29785 | Hydrophobin_D (04:28, FT:17) | Hydrophobin_D, Unclassified hydrophobin D | domain | MATCH | [PF29785](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF29785) |
| PF29802 | Hyd1F (04:28, FT:18) | Hyd1F, Class I hydrophobin F | domain | MATCH | [PF29802](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF29802) |
| PF29465 | Hydrophobin_like (04:28, FT:19) | Hydrophobin_like, Hydrophobin-like domain | domain | MATCH | [PF29465](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF29465) |
| PF12296 | HsbA (04:29, FT:20) | HsbA, Hydrophobic surface binding protein A | family | MATCH | [PF12296](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF12296) |
| PF07691 | PA14, flocculin/Epa lectin domain (04:27, FT:16) | PA14, PA14 domain | domain | MATCH on name. "Also in beta-glucosidases" not checked against the record text | [PF07691](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF07691) |
| PF00399 | PIR (04:31, 04:140, FT:15) | PIR, Yeast PIR protein repeat | repeat | MATCH | [PF00399](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00399) |
| PF22799 | PIR1-like_C, mannoprotein (04:32, FT:24) | PIR1-like_C, Cell wall mannoprotein PIR1-like, C-terminal domain | domain | MATCH | [PF22799](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF22799) |
| PF15789 | Hyr1 (04:33, FT:23) | Hyr1, Hyphally regulated cell wall GPI-anchored protein 1 | repeat | MATCH | [PF15789](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF15789) |
| PF11765 | Hyphal_reg_CWP (04:33, FT:14) | Hyphal_reg_CWP, Hyphally regulated cell wall protein N-terminal | repeat | MATCH | [PF11765](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF11765) |
| PF00624 | Flocculin (04:34, FT:21) | Flocculin, Flocculin repeat | repeat | MATCH | [PF00624](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00624) |
| PF13928 | Flocculin_t3 (04:34, FT:22) | Flocculin_t3, Flocculin type 3 repeat | repeat | MATCH | [PF13928](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF13928) |
| PF10182 | Flo11 (04:34, 02-none, FT:12) | Flo11, Flo11 domain | domain | MATCH | [PF10182](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF10182) |
| PF10528 | GLEYA (04:34, FT:13) | GLEYA, GLEYA domain | domain | MATCH | [PF10528](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF10528) |
| PF11766 | Candida_ALS_N (04:35, 04:79, 02:128, 02:175, FT:8) | Candida_ALS_N, Agglutinin-like protein, N-terminal domain | domain | MATCH | [PF11766](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF11766) |
| PF05792 | Candida_ALS (04:35, 02:128, 02:175, FT:9) | Candida_ALS, Candida agglutinin-like (ALS) | repeat | MATCH | [PF05792](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF05792) |
| PF30910 | ALS_M (04:35, FT:25) | ALS_M, Agglutinin-like protein middle domain | domain | MATCH | [PF30910](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF30910) |
| PF16541 | AltA1 allergen (04:36, FT:10) | AltA1, Alternaria alternata allergen 1 | domain | MATCH | [PF16541](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF16541) |
| PF25312 | Allergen_Asp_f_4 (04:36, FT:11) | Allergen_Asp_f_4, Allergen Asp f 4-like domain | domain | MATCH | [PF25312](https://www.ebi.ac.uk/interpro/api/entry/pfam/PF25312) |

Family-table `name` column check: all 24 FT names equal the InterPro `short` names above (PF05730 CFEM, PF04681 Bys1, PF01185 Hydrophobin, PF06766 Hydrophobin_2, PF28987 DewD, PF22354 Eas, PF11766, PF05792, PF16541, PF25312, PF10182, PF10528, PF11765, PF00399, PF07691, PF29785, PF29802, PF29465, PF12296, PF00624, PF13928, PF15789, PF22799, PF30910). No mismatch.

Report 01 and 03 contain no Pfam IDs other than PF03142, PF08407, PF01644 (01:40, 01:80, 01:86). Report 02 has only PF11766 and PF05792. Report 03 has none (03:66, 03:295 say Pfam "unverified").

## Table B: GH / CAZy claims (cazy.org)

| Family | Claim (report, file:line) | CAZy activity (fetched) | GPI-anchored fungal members described? | Result | Source |
|---|---|---|---|---|---|
| GH72 | GEL/GAS/PHR GPI beta-1,3-glucanosyltransferases (01:69, 01:80, 03:162, 03:295, 04:119) | beta-1,3-glucanosyltransglucosidase (EC 2.4.1.-), retaining | No (Pfam PF03198 text does say GPI-anchored) | MATCH on activity | https://www.cazy.org/GH72.html |
| GH17 | BGT-1/2 "may fall in GH17" (03:295); Bgl2/Exg-like (04:122); Gel4 "CAZy family 17 as stated in the file" (03:162) | endo/exo-1,3-glucosidases plus beta-1,3-glucanosyltransglucosidase and laminarin beta-1,6-glucanosyltransglucosidase | No | Family plausible for transglycosylases (MATCH on activity). The Gel4 = GH17 label (03:162) is UNVERIFIABLE here and sits against 01:69 and 04:119, which put GEL in GH72. Check Gel4 in CAZy/UniProt | https://www.cazy.org/GH17.html |
| GH16 | CRH/Utr2 (01:69, 04:120) | includes chitin beta-1,3/1,6-glucanosyltransferase, xyloglucan:xyloglucosyl transferase, plus hydrolases | No | MATCH on activity (Crh-type transferase listed). Crh1/Utr2 names not on page | https://www.cazy.org/GH16.html |
| GH76 | Dfg5/Dcw1 alpha-1,6-mannanases (01:69, 04:121) | alpha-glucosidase (EC 3.2.1.20) and mannan endo-alpha-1,6-mannosidase (EC 3.2.1.101) | No | Family activity MATCH. Dfg5/Dcw1 link UNVERIFIABLE from these records (still "recall") | https://www.cazy.org/GH76.html |
| GH18 | Cts1/Cts2 chitinases (01:80, 04:124) | chitin endo-beta-1,4-N-acetylglucosaminidase and others (10 activities) | No | MATCH on activity | https://www.cazy.org/GH18.html |
| GH5 | Exg1 glucanase (04:56, 04:123) | includes glucan exo-beta-1,3-glucosidase (EC 3.2.1.58) | No | MATCH on activity (Exg1 not named) | https://www.cazy.org/GH5.html |
| CBM18 | Chitin_bind_1 (04:68) | binds chitin; no fungal rows | n/a | UNVERIFIABLE Pfam link; binding MATCH | https://www.cazy.org/CBM18.html |
| CBM14 | CBM_14 / Avr4 (04:72) | binds chitin; no Avr4 mention | n/a | UNVERIFIABLE Avr4 link | https://www.cazy.org/CBM14.html |
| CE4 | chitin deacetylase (04:125) | not fetched | n/a | UNVERIFIABLE | none |

## Summary of mismatches and consequences

1. **Chitin synthase names, brief only.** Correct: PF01644 = Chitin_synth_1 (GT domain of classes I-III chitin synthases, also hyaluronan synthase), PF03142 = Chitin_synth_2 (fungal chitin synthase family), PF08407 = Chitin_synth_1N. The repo (04:14 and the 2026-10-04 report) is right; the brief had them swapped. No Pfam name in reports 01 to 04 or in `family_table.tsv` is wrong.
2. **Consequence of the swap if it had been used.** A rule written as "PF03142 = class I-III GT domain" would have tied the wrong model to the Division 1 / class I-III statement (01:40). Also, PF01644 matches animal hyaluronan synthase per its record, so it is not chitin-synthase-specific. Neither 01:80 (PF03142 as a broad CHS demotion model) nor 04:57/04:126 depends on the swapped names; they use the repo's names, which are correct. 01:40 describes PF03142 as "catalytic domain (stated)"; the record says family, but this is a paper-sourced wording, not an ID error.
3. **Unverified links that rest on correct IDs.** The IDs PF03663 (Glyco_hydro_76), PF00332 (Glyco_hydro_17), PF00722 and PF03198 are right, which settles the "my recall" and "to verify" ID readings in 04:230-231. What remains unverified is which proteins sit in them: Dfg5/Dcw1 in GH76, Bgl2 in GH17, Crh/Utr2 in GH16, and Exg1 in GH5 (CAZy gives families and activities, not these gene names).
4. **Open conflict.** 03:162 says Gel4 is "CAZy family 17 as stated in the file", while 01:69 and 04:119 put GEL proteins in GH72. This was not resolved from the fetched records and needs a check of Gel4 itself.
5. **Not verifiable from the fetched records:** PF10333 as Pig-X (record says Pga1), CE4 for PF01522, Pfam-to-CBM links for PF00187 and PF01607, Avr4, CalA, and the CFEM protein examples.
6. **No NOT FOUND, no obsolete or merged entries observed** (all 53 resolved). Clan and merge history were not checked.
