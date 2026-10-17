# Lectin class design: review 3

*2026-10-10. Independent review of `docs/superpowers/specs/2026-10-10-lectin-class-design.md`, revision 3 (commit `2cc718c`, branch `lectin-class-spec-r3`). Reviewer: Claude Code (claude-opus-5-5). The spec is not edited. Review files are in `analysis/lectin_class_spec/review3/`. Login-node runs only. No SLURM jobs. Section numbers and fact numbers (F1 to F26) refer to revision 3.*

## Verdict

**Pass with changes. The changes are not minor. The spec needs a short revision 4 before it goes to the owner.** A full fourth review is not needed. The author can close the four required changes (section 4) and state how each was checked.

The work of revision 3 is sound. Every number that I re-ran is correct (section 2). All review 2 findings are fixed in the files, not only in words (section 1). The recommendation (option (c)) is defensible, and the spec labels the ranking as a judgement.

Four problems remain. Two of them change what the owner reads about the main decision:

1. **The recommended option (c) cannot be built as section 3.4 describes it** (H1). Inactive family rows produce no output in the current code, the Pfam module writes no score, and the module list is a fixed tuple. The report that (c) promises (section 6.2 items 1 and 2) has no implementation path in the spec.
2. **The spec's own rule does not give an empty M1** (H2). CBM_14 (Avr4) has H1b, whole-proteome U 4.7%, and is not owned or house-keeping. By the rule of section 3.1 it joins M1. D11 removes it with a reason ("the rule is about H1a") that the rule text does not say.
3. **The option (a) population is not described, and it is the most specific signal in the data** (M1). The 32 pool S PF07938 proteins are 6 full-length, high-score proteins, all in Onygenales (dermatophytes and relatives), plus 26 weak short hits in other genera. The spec calls them "a different subgroup" without data.
4. **The FleA role summary omits a knockout adhesion assay** (M2). The Sakai 2019 abstract says that adhesion of ΔfleA conidia to cultured cells was measured. D7 says the "not an adhesin" view "rests on disease outcome".

## 1. Review 2 findings: fixed in the files?

| Review 2 finding | Result | Evidence |
|---|---|---|
| H1. U on the wrong population | **Fixed.** U is now defined on the population the call reads (3.1). Whole-proteome U is in F25 and the 3.1 table. Options (a), (b), (c) are presented. D1, section 1, 4.4, owner question 1 and Risk 10 are rewritten. But see H2: the statement "the ungated M1 is empty" does not follow from the rule. | spec 3.1, F25; my recount (section 2) |
| H1 (plan). SignalP 6 before D1 | **Fixed.** SignalP 6 was run (F26). L2 is first in section 11. | `signalp6/`; my re-run (section 2) |
| H2. Paracoccin | **Fixed.** Coltri 2006 added and correct. The length claim is removed. The 70 kDa versus 218 aa conflict is stated. The partial GH18 hit, the lectin sites inside GH18, the silencing statement and the D9 consequence are all stated. One new point: Pitangui 2021 is a review (L1 below). | F20, D9, 6.2 item 6 |
| H3. D7 contradicts the report text | **Fixed.** The label names the family. 6.2 item 3, D7 and owner question 4 agree. Small remainder: the option (c) column of 3.4 uses "lectin domain" in the class text, the (a)/(b) column uses "domain". Use one form. | 3.4, 6.2, D7 |
| M1. Name-search scope; Avr4 | **Fixed in words and data.** F4 states the limit. Avr4 was scanned (F24), and CBM_14 and Chitin_bind_1 were scanned on 653 proteomes. But the result creates H2 below. | F24; `r3/` files; my recount |
| M2. H bias; H1 mixed mechanisms | **Fixed.** H1a and H1b are separate. The bias is stated in 3.1 and Risk 10. The CVNH row gives the reason for H2. | 3.1, Risk 10 |
| M3. One-model call | **Fixed.** Stated; the call is named for the family. | 3.1, 3.4 |
| M4. Five papers | **Fixed.** All five are cited. I checked all five PMIDs (section 2). | F14, F16 |
| M5. Look-alikes in F23 | **Fixed.** | 4.3, L3 |
| L1. `cluster_check.py` | **Fixed.** I re-ran it (section 2). | commit `4d8bdee` |
| L2. Proteome release date | **Not done; stated** in section 10. Acceptable. | section 10 |
| L3. Heading | **Fixed.** | 4.4 |
| L4. "by SignalP in UniProt" | **Fixed.** | F20 |
| Q4. D5 overlap in the report | **Fixed.** | D5 |

**The author's one disagreement (SignalP 6 lowers option (a)).** I agree with the fact and in part with the conclusion. SignalP 6 gives no signal peptide to any of the 6 AAL-like structure sequences (I reproduced this, section 2). So the H1a evidence (FleA, SapL1) does not come from proteins that a gated call would call. That is a valid reason not to rank (a) above (b). But the spec does not say what (a) **does** call. My check (M1) shows that (a) calls a clade-specific group of full-length secreted AAL-like proteins in Onygenales. That is a fact the owner needs in order to judge (a). The ranking "(a) third" can stay as a judgement, but the option table must describe the population.

## 2. Checks I ran

| Claim | Check | Result |
|---|---|---|
| F25 whole-proteome U | `review3/recount_u.py` on `review2/whole_proteome_hits.tsv` and `r3/cbm_whole_proteome_hits.tsv`, against `analysis/fungi_scan/selection_discovery.tsv` (653 proteomes, 653 distinct genera, 6,847,098 proteins) | PF07938: 997 proteins, 392 genera (60.0%); Onygenales 11, Aspergillaceae 4. Score 50 or more 188 (28.8%); 100 or more 104 (15.9%); 200 or more 52 (8.0%); 300 or more 40 (6.1%). FB_lectin 75 / 68 (10.4%). LDL 39 / 36 (5.5%). Gal-bind_lectin 41 / 31 (4.7%). CCL2-like 2 / 2 (0.3%). CBM_14 54 / 31 (4.7%). Chitin_bind_1 2,222 / 572 (87.6%). **All match.** Output: `review3/recount_u.txt`. |
| Independent hmmsearch | `review3/rescan_pf07938.sh`: `hmmsearch --cut_ga` (HMMER 3.4), PF07938 from Pfam 38.2 (`hmmfetch`; GA 25.3, NSEQ 3), 25 proteomes drawn at random (seeded `shuf`) from the 653 | 36 proteins in 14 of 25 proteomes. **Identical protein set and identical scores** to the committed table (maximum difference 0.0). Output: `review3/compare_rescan.txt`, `rescan_sample.txt`. |
| F26 SignalP 6 | Read `signalp6/prediction_results.txt` (57 rows: 15 SP, 42 OTHER). Re-ran `signalp6 --organism eukarya --mode fast` on the 6 AAL-like sequences plus Ecp6 (B3VBK9) and Emp46 (Q12396) (`review3/signalp6_check.sh`) | All 6 AAL-like OTHER; P(SP) FleA 0.000002, SapL1 0.000000, AAL 0.000420, AOL 0.000017, *A. oryzae* 0.000001, *Ralstonia* 0.002805. Ecp6 SP (cleavage 18-19), Emp46 SP 0.711751 (cleavage 46-47). **Identical to the committed values.** Note: the module `signalp/6` points to the 6.0i package; its binary reports "6.0h". F26 says "6.0h ... module `signalp/6`". Both are true; say so (L4). |
| F25 truth scores, F24 Avr4 | Read `r3/pf07938_truth_scores.tbl`, `r3/avr4_hmmscan.tbl`; CBM_14 GA from Pfam 38.2 | Q2UNX8 525.1, FleA 489.1, SapL1 300.7, AAL 158.6, AOL 132.6, *Ralstonia* 42.3. Avr4 CBM_14 28.7, E 1.8e-6, GA 22.1. **Match.** |
| F10 clusters | MMseqs2 17-b804f `easy-cluster --min-seq-id 0.3 -c 0.5` on the committed FASTA; `cluster_check.py` | 56 sequences, 33 clusters, same partition; 43 / 25 after exclusions; five models 16 / 7; PF07938 5 / 3. My cluster file is **byte-identical** to `unilectin_fungal_clusters.tsv`. Extra set (five models plus CVNH and RicinB_lectin_2): 24 sequences, 12 clusters. |
| Pool S PF07938 population (new) | `review3/pools_pf07938_scores.py`: pool S proteins with PF07938 in the local pools (`_workdir/discovery/pools/`), joined to their Pfam 38.2 whole-proteome scores | See M1. Output: `review3/pools_pf07938_scores.txt`. |
| Pfam seed leakage (new) | Downloaded the PF07938 seed alignment (InterPro API, `annotation=alignment:seed`); `mmseqs easy-search` of the AAL-like truth sequences against the 3 ungapped seed sequences | See L2. Output: `review3/pf07938_seed.faa`, `pf07938_seed_vs_truth.m8`. |
| Code path for option (c) (new) | Read `src/cellsurface_sorting_hat/modules/pfam.py` and `modules/cli.py` | See H1. |

**PMIDs.** I fetched PubMed metadata for 12 PMIDs: the 7 new ones (16476564, 23695231, 25760594, 30277619, 12438366, 22163332, 17153926) and 5 older ones with load-bearing claims (34373510, 27058347, 29370403, 34532342, 26455798). All exist. Title, journal, volume and pages match the spec for all 12. I compared each claim with the abstract:

- Coltri 2006: 70-kDa minor protein; GlcNAc binding; peptides "did not match any known protein"; surface of yeasts at budding regions; laminin binding inhibited by GlcNAc; TNF-alpha; "can favor yeast adhesion" as a suggestion. **Match.**
- Kuboi 2013: 32-kDa L-fucose hemagglutinin from mycelium; `AfufleA`, 314 aa; role "not clear"; "may enhance attachment". **Match.**
- Houser 2015: six non-equivalent sites; micromolar affinity for oligosaccharides on human epithelia; "might play an important role". **Match.**
- Lehot 2018: nanomolar multivalent fucosides; "highly efficient in preventing fungal spore adhesion to bronchoepithelial cells". **Match.** The abstract also ascribes part of the effect to conidium aggregates. The spec could add this; it weakens the inhibitor line a little.
- Tronchin 2002: 32-kDa sialic-acid-specific lectin from conidia; distinct from alkaline protease; no gene named. **Match.**
- Powers-Fletcher 2011: ClxA is an ER lectin chaperone; ΔclxA stress defects; virulent in two models. **Match.**
- van den Burg 2006 (Avr4): "chitin-binding lectin containing an invertebrate chitin-binding domain (CBM14)"; protects against plant chitinases; "counter-defensive virulence factor". **Match.**
- Martinez-Alarcon 2021, Kerr 2016, Koharudin 2015: **match.**
- **Sakai 2019: claim incomplete.** See M2.
- **Pitangui 2021: PubMed article type is "Review".** See L1.

## 3. Findings, ranked by severity

### High

**H1. Option (c), the recommended option, cannot be built as written, and the spec does not say what must be built.**

- Section 3.4, option (c) column: "All lectin-like models get `module=pfam_lectin_like`, `active=no` or report-only". The "Effect" row says "None on `config_sha256` if no call is added and rows are inactive (not checked in code)".
- The code: `pfam_rows` uses only active families (`wanted = {... if f.module == module and f.active}`, `modules/pfam.py` line 135). A module with no active family writes every protein as `unavailable` with the note "no active family" (`modules/cli.py` lines 292-305). So inactive rows give no hit and no count.
- `MODULES = ("pfam_adhesion", "pfam_allergen", "pfam_hydrophobin", "pfam_hsba")` (`modules/pfam.py` line 22). A family row with `module=pfam_lectin_like` raises `FamilyTableError`. A new module name is a code change.
- The Pfam module rows hold `hit` and `families`. They hold no score. Section 6.2 item 1 under (c) promises "the count of PF07938 hits per proteome, with the score of each". The module cannot supply the score.
- An active row needs `active_by` and `active_date` (`pfam.py` lines 68-69), that is, an owner activation. The spec does not ask the owner for it.

Consequence: under (c) as written, the report of 6.2 items 1 and 2 would be empty. The value of (c) to the owner rests on that report ("nothing is lost that the report does not still show", 3.1 reason 5). Revision 4 must say what (c) builds: active rows in a new module that no call reads (and confirm that the engine accepts a module that no call reads), a code change to `MODULES`, a score column or a separate report table, and the owner activation. Add these to L4 and to D1 or a new decision. Check the `config_sha256` claim in code, or keep "not checked" visible in the option table.

**H2. Under the spec's own rule, the ungated M1 is {CBM_14}, not empty. D11 overrides the rule with a reason the rule does not contain.**

- Rule (3.1): "A model joins a call if H is H1a or H1b, U is rare, and the model is not owned by an active narrow-tier family and not a house-keeping class."
- CBM_14 row of the table: H1b (Avr4), whole-proteome U 31 of 653 genera (4.7%, rare), not owned, not house-keeping. It meets every condition.
- 3.1 says: "Only these two models have H1 (H1a for PF07938, H1b for LysM). So the ungated M1 is empty." Section 1 point 1 and owner question 1 repeat the conclusion. This is false for the table as revised.
- D11 reason (1): "The rule and the tier are about lectins that bind host glycans (H1a)." The rule text admits H1b. Reason (2) (one paper; the only known protein scores 28.7 against GA 22.1) is a valid separate reason.

Revision 4 must choose one: (i) change the rule to H1a only for this tier, and say that H1b models (LysM, CBM_14) go to a later chitin-binding effector spec by rule; or (ii) keep the rule, state that it gives M1 = {CBM_14}, and state that D11 overrides it, with reason (2). Option (i) is simpler and matches the tier text of 3.3. Either way, correct section 1, 3.1 and owner question 1.

### Medium

**M1. The option (a) population is not described. In the data it is a clade-specific group of secreted, full-length AAL-like proteins in Onygenales, plus a tail of weak hits.**

My join of the 32 pool S PF07938 proteins to their Pfam 38.2 scores (`review3/pools_pf07938_scores.txt`):

| Group | Proteins (genera) | Sequence score | Length (pool table) |
|---|---|---|---|
| Onygenales: *Myotisia* sp. PD_48, *Microsporum canis*, *Arthroderma uncinatum*, *Onygena corvina*, *Nannizzia gypsea*, *Trichophyton benhamiae* | 6 (6) | 288.3 to 330.6 | 327 to 337 aa (3 checked) |
| Other genera (background) | 26 (22) | 26.2 to 66.6 | 119 aa in the one I checked |

- At score 100 or more, pool S holds 6 proteins in 6 genera, **all Onygenales**: 6 of 27 Onygenales genera against 0 of 626 other genera. FleA scores 489.1 and is 315 aa, so the 6 are in the FleA score range and length range.
- All 6 are in pool U (no current mechanism). They are the discovery-pass PF07938 enrichment (F6: 6 of 27, q 0.0163).
- In Aspergillaceae, the high-scoring PF07938 proteins (*Penicillium rubens* 561.6, *Hamigera avellanea* 404.6, *Evansstolkia* 241.2) are **not** in pool S. That fits FleA (no signal peptide).
- PubMed search "(Trichophyton OR Microsporum OR Arthroderma OR Nannizzia OR dermatophyte) AND (fucose-binding OR fucose-specific OR AAL OR Aleuria aurantia lectin) AND (gene OR protein OR secreted)": 0 hits. A broader search found activity-level work only: Esquenazi 2003 (PMID 12628546, [doi:10.1016/S0928-8244(03)00007-5](https://doi.org/10.1016/S0928-8244(03)00007-5)) and Esquenazi 2004 (PMID 15059626, [doi:10.1016/j.resmic.2003.12.002](https://doi.org/10.1016/j.resmic.2003.12.002)) report mannose- and galactose-binding adhesin activity on *Trichophyton* microconidia, with no gene named. That specificity is not fucose. I do not claim a link to the PF07938 proteins.

What this changes. The spec's statement (option (a) row) that the 32 pool S proteins "are a different subgroup" is now partly measured: 6 of them are full-length AAL-like proteins with a predicted signal peptide, restricted to one clade. The other 26 are weak, short hits. A gate plus a score floor (for example R0 and score 100) would call exactly the 6 Onygenales proteins. This variant is not in the option table. It still has no H evidence: no paper shows binding or function for a dermatophyte AAL-like protein. So it does not by itself overturn the recommendation of (c). But it is the most specific pattern in the data, it is in the clade the owner cares about, and the owner should see it when choosing among (a), (b) and (c). Revision 4 should add this table (or the counts) to F25 or section 6.3, describe the population in the option (a) row, list the gate-plus-score variant (or say why not), and add the 6 proteins to L3/L6 as T4 review targets. Under (c), 6.2 should report them by name.

**M2. Sakai 2019 tested adhesion of ΔfleA conidia. F14 and D7 do not say so.**

- Sakai 2019 abstract (PMID 29370403, [doi:10.1093/mmy/myx163](https://doi.org/10.1093/mmy/myx163)): "The pattern of conidial phagocytosis and adhesion to cultured cells did not explain this enhanced pathogenicity."
- So the knockout study measured adhesion to cultured cells. The abstract does not give the direction or size of the effect.
- F14 lists Sakai 2019 only for "earlier germination and enhanced pathogenicity". The summary says "Papers (b), (c), (h) show that gene loss gives no loss of function in the readouts used" and "A disease outcome does not show whether FleA mediates attachment". D7: "C rests on disease outcome, which does not test attachment."
- Kerr 2016 also has a knockout readout near attachment: macrophages "bind and phagocytose" ΔfleA conidia less well.

Revision 4 must add the Sakai adhesion assay to F14 (as an abstract-level statement without direction) and rewrite the D7 reason for not choosing C. The full text of Sakai 2019 is already in L3b. It should be read before the owner answers question 4, because it is the one knockout adhesion test of FleA that I found.

### Low

**L1. Pitangui 2021 is a review.** PubMed article type for PMID 34532342 is "Journal Article; Review". F20 uses it for the gene ID, the signal peptide, the GH18 range, the lectin-site coordinates and "PCN-silencing does not affect the adherence or phagocytosis". Say that it is a review, and that the silencing result comes from primary work it cites (not identified; a PubMed search "paracoccin silencing adherence Paracoccidioides" gave 0 hits).

**L2. Pfam seed leakage can now be stated, not "not checked".** The PF07938 seed (InterPro, Pfam 38.2 entry PF07938.18) has 3 sequences: K1XQV1 (*Marssonina brunnea*), B6HN96 (*Penicillium rubens*), B8NIA0 (*Aspergillus flavus*). Best identity to a seed: *A. oryzae* Q2UNX8 100% (302 aligned residues, to B8NIA0); FleA 81.7%; SapL1 43.7%; AOL 36.6%; AAL 32.8%. All five fungal truth sequences are at or above 30% identity to a seed, so under COMMON-RULES rule 6 each is `in_reference`. This confirms the spec's leakage label and the statement that the truth sequences cannot set a cutoff for option (b). Replace "not checked" in 4.1 item 2, 4.4 statement 2 and section 10.

**L3. Als N-terminal domain is a lectin-type model owned elsewhere and is not named.** Candida_ALS_N PF11766 ("Agglutinin-like protein, N-terminal domain") contains "agglutinin", so it is in the F4 name list. It is active in `pfam_adhesion` (narrow tier) in `family_table.tsv`. Ielasi 2016 (F15) reports glycan arrays of Als1p and Als3p. F15 says "I did not check which Pfam family covers their glycan-binding region". Add PF11766 to the "Owned elsewhere" row of 3.1 and to F1, and close the F15 sentence.

**L4. Small corrections.**
- F26: the module `signalp/6` loads the 6.0i package; its binary reports "6.0h". State both, so that the R0 wrapper run (L2) can match the version string (the tool's `check_signalp` compares version strings, `tests/cellsurface_sorting_hat/calibration/test_phasec.py`).
- `r3/cbm_whole_proteome_summary.txt` prints "seq score >= 25.3" for CBM_14 and Chitin_bind_1. 25.3 is the PF07938 GA, not theirs (CBM_14 GA 22.1). At 25.3, CBM_14 is 53 proteins in 30 genera. The spec quotes the GA counts correctly; label or remove the 25.3 lines.
- Lehot 2018: the abstract ascribes the anti-adhesion effect partly to "the formation of conidium aggregates". Add one clause to F14 (k).

**L5. Literature additions (no change to the outcome).**
- Mentlak 2012, Plant Cell 24(1):322-35, PMID 22267486, [doi:10.1105/tpc.111.092957](https://doi.org/10.1105/tpc.111.092957), abstract: *Magnaporthe oryzae* Slp1, a secreted LysM protein, binds chitin, competes with rice CEBiP and "is required by M. oryzae for full virulence". This is an Ascomycota H1b case with a gene-loss phenotype. It belongs beside Ecp6 in the LysM row and in the future chitin-binding effector spec.
- Liu and Stajich 2015, Fungal Genet Biol 77:31-9, PMID 25819009, [doi:10.1016/j.fgb.2015.03.003](https://doi.org/10.1016/j.fgb.2015.03.003), abstract: CBM18 (Chitin_bind_1) lectin-like proteins of *Batrachochytrium dendrobatidis* bind chitin and protect *Trichoderma* against chitinase in vitro. The Chitin_bind_1 row says "H0 (no paper read)". The paper is chytrid, not Ascomycota, and has no host target, so H is H1b-like at most. Chitin_bind_1 stays out (87.6% of genera).
- My other searches for fungal-encoded lectins in Ascomycota pathogens ("Aspergillus fumigatus"[Title] AND lectin[Title] NOT FleA: 14 hits, all FleA papers already cited, calnexin, or host-side; *Sporothrix*, *Fusarium*, *Botrytis*, *Magnaporthe*, *Pneumocystis*, *Candida* with lectin and host terms: 20 of 90 titles read, all host-side, plant lectins as probes, or Epa work already cited) found no other fungal lectin with host-binding data. This is "not found by these searches".

**Claims that rest on abstracts only.** All of F14 except Houser 2013; all of F21; Coltri 2006; van den Burg 2006; all five review 2 papers. The spec says so in section 10. The claims I compared match their abstracts. The two that carry the most weight and should be read in full before D7 and D11 are Sakai 2019 (M2) and van den Burg 2006 (D11 relies on "one paper").

## 4. Answers to the specific questions

**Q3. Is option (c) sound, or an avoidance of the owner's request?** Sound, with conditions. The owner asked for help with models, truth data, literature and curation. Under (c) the spec still gives: a model table with H and U for 16 models, truth tiers, a curation rubric, negatives, a literature base with checked PMIDs, and a work plan to build the truth that a later call needs (L1, L3, L3b). That is the requested help. The data do not support a call now: the ungated GA call fires in 60% of genera, the gated call has no H evidence for what it calls, and a score cutoff has no truth to set it (L2 shows that every candidate truth sequence is in the model's reference). Choosing no call is not avoidance. But (c) is only useful if its report exists (H1) and shows the Onygenales group (M1). Without those, (c) delivers nothing to a user of the tool.

**Fairness of (a) and (b).** (b) is fairly presented, including that the cutoff decides the rule result. (a) is fairly presented on SignalP but not on what it calls (M1). The ranking is stated as reasoning only, in 3.1, D1 and owner question 1. That is honest.

**Is the U rule principled, and is the 10% threshold still decisive?** U is now measured on the right population. The 10% threshold is moot for the ungated GA call (60%) and for (a) (4%; any threshold above 4.3% gives "rare"). It is decisive only for (b), between a cutoff of 100 (15.9%) and 200 (8.0%). The spec says this. One limit the spec does not state: U counts genera with at least one hit. It does not measure how many proteins per proteome a call would flag, which is what "flooding" means for a user. For (c) this does not matter.

**Q4. Validation design.**
- Circularity and leakage: correct, and now confirmed by the seed check (L2).
- 20-cluster floor: arithmetic correct (F10 reproduced). "smoke at most" is honest.
- D4b (no status file for one in-reference positive): consistent with `status.py` and `call_status.py` as review 2 checked. I agree.
- D5 double counting: consistent for both answers.
- Biology: FleA wording is fair, with the M2 gap. Paracoccin handling is now correct (L1 is a labelling point). LysM and CVNH placement is correct by ubiquity; add Slp1 (L5). D11 needs the H2 fix. ER house-keeping (M3): PF07915 is supported by the OS-9 hits (review 2); calnexin has a paper; Malectin stays "from memory" with 1 pool S protein, an acceptable risk.

**Q6. Work plan and owner questions.**
- The plan is realistic in size. L4 must change for H1 (code change to `MODULES`, score output, owner activation of report-only rows). L2 (R0 wrapper run) matters only for (a) and (b); under (c) it can follow D1.
- Owner question 1 is answerable but dense. After M1 it should name the Onygenales group in one sentence. After H2 it must not say "M1 is empty" without the D11 qualifier.
- Owner question 3 asks two things (accept `smoke`/`unvalidated`, and no status file). Split it or say "both".
- A question is missing under (c): "Do you agree to activate report-only lectin rows in a new module that no call reads?" (H1). The owner must sign `active_by`.
- Questions 2, 4 to 8 are answerable without extra context. Question 4 should wait for the Sakai full text (M2).

## 5. Required changes for revision 4

1. **(H1)** Describe what option (c) builds: a new module name in `MODULES`, active report-only rows (owner activation), a score output or a separate report table, and confirmation that a module read by no call is accepted by the engine. Update 3.4, 6.2, L4 and add an owner question. Check or flag the `config_sha256` claim.
2. **(H2)** Make the rule and D11 agree: H1a-only rule for this tier, or M1 = {CBM_14} with an explicit D11 override. Fix section 1 point 1, 3.1 and owner question 1.
3. **(M1)** Describe the pool S PF07938 population (6 full-length Onygenales proteins at 288-331, 26 weak hits at 26-67). Describe it in the option (a) row and in 6.3. State the gate-plus-score variant or the reason to leave it out. Name the 6 proteins as T4 review targets and as report rows under (c).
4. **(M2)** Add the Sakai 2019 adhesion assay to F14 and rewrite the D7 reason against option C. Read Sakai 2019 in full before owner question 4, or say that question 4 waits for it.

Minor (can go with the above): L1 (Pitangui 2021 is a review), L2 (seed leakage confirmed), L3 (PF11766 owned elsewhere), L4 (SignalP version text, 25.3 lines, Lehot clause), L5 (Slp1, CBM18 paper), one form of the PF07938 class text in 3.4, split owner question 3.

## 6. Ready for the owner?

**No, not yet. Revision 4 is needed.** Items 1 to 3 change what the owner reads about the main decision. Item 4 changes the reasoning of D7. All four are text changes plus one small data table that already exists (`review3/pools_pf07938_scores.txt`). None needs a new compute run. After revision 4, the author can send the spec to the owner without a fourth full review, provided the change table shows how each of the four items was checked.

## 7. What I did not do

- I did not read any paper in full text in this review. All literature checks are at abstract level.
- I did not re-run the whole-proteome scan on all 653 proteomes. I re-ran PF07938 on 25 of them (exact match) and recounted all models from the committed tables.
- I did not run SignalP on the 6 Onygenales pool S proteins myself; their signal peptide comes from the pool definition (R0 in the discovery pass). I checked the length of 3 of the 6 in the pool table, not all 6.
- I did not test whether the engine accepts a module that no call reads (H1). I read the Pfam module code only.
- I did not check whether the 6 Onygenales proteins are orthologs of each other or of FleA. Their scores and lengths are the only evidence that they are full-length AAL-like proteins.
- Searches beyond those listed in L5 and M1 were not made. Negative results are "not found by these searches".
- Scratch files (`_review3_scratch/` in the worktree) are not committed. The scripts and small outputs are in `analysis/lectin_class_spec/review3/`.
