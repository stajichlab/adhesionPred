# Lectin class design: review 4

*2026-10-10. Independent review of `docs/superpowers/specs/2026-10-10-lectin-class-design.md`, revision 4 (commit `730fab2`, branch `lectin-class-spec-r4`). Reviewer: Claude Code (claude-opus-5-5). The spec is not edited. Scope: the four required changes of review 3 (H1, H2, M1, M2), a PMID and data spot check, and readiness for the owner. Login-node runs only. No SLURM jobs. New review files are in `analysis/lectin_class_spec/review4/`.*

## Verdict

**Pass with changes. The changes are small. The spec can go to the owner after two short text fixes (section 4, "must fix").** No fifth review is needed.

All four review 3 items are closed in substance:

- H1: every statement about the code in section 3.4 is correct. I read the code myself. Option (c-i) needs no change to `src/`.
- H2: the H1a-only rule is stated with reasons, applied to the models in the table, and the CBM_14 outcome is shown as an owner option (D11, question 7). One set label does not follow the rule (finding M-B).
- M1: F27 and the whole-proteome count are correct. The option (a2) text and its circularity are fair. One sentence in section 1 and owner question 1 overstates the clade restriction (finding M-A).
- M2: my own fetch of the publisher page returns the same Sakai 2019 result as the spec, plus Discussion sentences. The D7 wording is fair. F14 can now say more (finding L-A).

## 1. Checks I ran

| Claim | Check | Result |
|---|---|---|
| 3.4 (1) `MODULES` fixed tuple; other names rejected | Read `modules/pfam.py` lines 22, 60-61 | **Correct.** `MODULES = ("pfam_adhesion", "pfam_allergen", "pfam_hydrophobin", "pfam_hsba")`; `load_family_table` raises `FamilyTableError` for any other module. |
| 3.4 (2) only active families; no active family gives `unavailable` | `pfam.py` lines 135-137; `modules/cli.py` lines 276-307 | **Correct.** `pfam_rows` filters `f.active`; the command loops over `pfam.MODULES` and writes all proteins `unavailable`, note "no active family", when none is active. |
| 3.4 (3) no score in Pfam rows | `pfam.py` line 24 and lines 168-170 | **Correct.** `COLUMNS = ["hit", "families"]`. `parse_domtblout` reads `seq_score` (line 109) and the domain score, but the rows do not write them. |
| 3.4 (4) active row needs `active_by`, `active_date` | `pfam.py` lines 68-69 | **Correct.** |
| 3.4 (5) `config_sha256` = SHA-256 of `categories.yaml` bytes | `engine.py` `load_config` lines 82-95 | **Correct.** `sha256=hashlib.sha256(raw).hexdigest()` on `path.read_bytes()`. The family table does not enter it. The per-module digest (`_family_digest`, `modules/cli.py` 225-232) hashes only that module's rows. |
| 3.4 (6) engine loads every table in `modules/` | `cli.py` `load_modules` lines 183-191; `engine.py` `referenced_modules` 403-409; `cli.py` 394-399 | **Correct.** `load_modules` globs `*.tsv.gz` and `*.json`. `referenced_modules` is used only for `absent_modules`. Evidence fields come only from `cfg.evidence` (`collect_evidence`). |
| Call status stays valid when a module that no call reads is added | `call_status.py` `stale_reason` lines 177-197 | **Yes** (not stated in the spec; it supports (c-ii)). Staleness checks `config_sha256` and the modules that the call reads. A new module read by no call does not make a call status file stale. |
| (c-i) buildable with no `src/` change | Spec 3.4; review 2 and r3 scan scripts exist | **Yes.** It is an `analysis/` script plus a report table. The committed scripts `review2/scan653.sh` and `r3/scan_cbm.sh` already do this kind of scan. |
| F27 pool S population | Read `review3/pools_pf07938_scores.txt`; counted genera; read lengths of all 32 from `_workdir/discovery/pools/` | **Counts correct:** 6 Onygenales proteins at 288.3 to 330.6; 26 others at 26.2 to 66.6 in **22** genera; at 100 or more 6 proteins, 6 genera. **Two small errors** (finding L-B): Onygenales lengths are 327, 336, 337, 337, 337, **356** (not "327 to 337" or "about 330"); 5 of the 26 weak hits are 400 to 454 aa, so they are not all "short". |
| Whole-proteome Onygenales at score 100 or more | Ran the `awk` command of spec section 12 | **Correct:** 7 proteins in 6 genera; *A. uncinatum* FA49FE15_005258-T1 at 213.6. FA49FE15_005258-T1 is not in the pool table (so not pool S). |
| Whole-proteome genera at score 100 or more by group | Same table, grouped | Onygenales 6, Aspergillaceae 3, background 95 (104 total, matches F25). See finding M-A. |
| Pool S under Pfam 38.2 (new) | `review4/pool_vs_382.py`: every Pfam 38.2 GA PF07938 hit joined to pool membership | **36** pool S proteins hit PF07938 at Pfam 38.2 GA, not 32. The 4 extra have no Pfam 37.x PF07938 label in the pool table; all are background, scores 27.3 to 37.6. The score-100 set is unchanged (the same 6). Finding L-C. Output: `review4/pool_vs_382.txt`. |
| F28 seed leakage | Read `review3/pf07938_seed_vs_truth.m8` | **Correct:** Q2UNX8 1.000 (302 aa), FleA 0.817, SapL1 0.437, AOL 0.366, AAL 0.328, *Ralstonia* 0.354 (82 aa). Rule 6 in `COMMON-RULES.md` uses 30% identity. |
| F10 clusters | MMseqs2 b804fbe (17-b804f) `easy-cluster --min-seq-id 0.3 -c 0.5`; `cluster_check.py` | **Reproduced:** 56 / 33, same partition; 43 / 25; five models 16 / 7; PF07938 5 / 3. Output **byte-identical** to `unilectin_fungal_clusters.tsv`. |
| F24 CBM_14 | Counted `r3/cbm_whole_proteome_hits.tsv` | 54 proteins, 31 genera. **Correct.** |

**PMIDs (PubMed metadata and abstracts, this review).** I checked 10: 29370403 (Sakai 2019), 22267486 (Mentlak 2012, Slp1), 25819009 (Liu and Stajich 2015), 27058347 (Kerr 2016), 34373510 (Martinez-Alarcon 2021), 16476564 (Coltri 2006), 17153926 (van den Burg 2006), 24340081 (Houser 2013), 30277619 (Lehot 2018), 34532342 (Pitangui 2021). Title, journal, volume and pages match for all 10. The claims match the abstracts:

- Slp1: binds chitin, suppresses chitin-induced immunity, competes with CEBiP, "required by M. oryzae for full virulence". **Match.** H1b is the right level: the bound glycan is fungal chitin.
- Liu and Stajich 2015: CBM18 lectin-like genes; chitin binding in vitro; co-localisation with chitinous cell boundaries when expressed in *S. cerevisiae*; protection of *T. reesei* against chitinase. **Match.** No host target, so "H1b-like at most" is fair. (Note for transparency: the owner is a co-author of this paper. This does not change the classification.)
- Pitangui 2021: article type "Journal Article; Review". **Match** with F20 as revised.
- Lehot 2018: the "conidium aggregates" clause is in the abstract. **Match** with F14 (k).

## 2. Findings, ranked by severity

No high-severity finding.

### Medium

**M-A. Section 1 and owner question 1 overstate the clade restriction of PF07938.**

- Section 1 point 1: "Option (a2) calls 6 full-length PF07938 proteins that occur only in Onygenales (6 of 27 genera in pool S, 0 of 626 other genera, F27)."
- Owner question 1: "PF07938 has 6 full-length, secreted-type proteins in 6 Onygenales genera ... and none in 626 other genera at score 100 or more".
- The data: on whole proteomes, PF07938 at score 100 or more is in **98 non-Onygenales genera** (3 Aspergillaceae, 95 background) and in 6 of 27 Onygenales genera. As a share: 22% of Onygenales genera, 15.7% of the other genera. The clade restriction exists only after the pool S filter (signal peptide, no mature TM helix, no HDEL/KDEL; discovery-pass report). The proteins themselves are not Onygenales-only. Only the pool S form of them is.
- F27 and the (a2) row state the whole-proteome number (104 genera). But section 1 and question 1 are what the owner reads first, and "occur only in Onygenales" is false as written.
- Fix: in section 1 and question 1, say "6 full-length PF07938 proteins **in pool S** (signal peptide, no TM, no ER retention) occur only in Onygenales; without the pool S filter, PF07938 at score 100 or more is in 98 other genera". One sentence each.

**M-B. CBM_14 and Chitin_bind_1 are in a set that the rule does not define, and the report does not show them.**

- Rule (3.1): "Models that do not join go to a reported-only set M2 or an exclusion set M3." M3 is for house-keeping models.
- CBM_14 (H1b) and Chitin_bind_1 are put in a third set, "out of v1". LysM, also H1b, is in M2 and is reported (6.2 item 2). So two H1b models are reported and treated by the rule, and one is not. The rule as amended would put CBM_14 in M2.
- Under option (c) the report shows "all lectin-like models". It does not say whether CBM_14 and Chitin_bind_1 appear.
- Fix: put CBM_14 and Chitin_bind_1 in M2 (reported, not called), or state in 3.1 why they are outside both M2 and the report (for example: outside the F4 name-search scope, deferred to the chitin-binding effector spec). Either is fine. The table and the rule must agree. One sentence.

### Low

**L-A. Sakai 2019: the direction can now be stated with more support, and there is a second conflict with Kerr 2016.**

- I fetched the Oxford Academic page (`academic.oup.com/mmy/article/57/1/71/4820940`) with a web-fetch summary tool. It is the same kind of tool as the author used, so this is a second fetch, not a PDF reading. It returned the Results sentences that the spec quotes, and also:
  - Results: "the rates of phagocytosis of the disruptant and parental strains were not different (Fig. S1)."
  - Discussion: "The conidial phagocytosis by and adhesion to cultured cells were not significantly different ... suggesting that these factors are not responsible for the enhanced pathogenicity."
  - Discussion: "Our data did not support the FleA-dependence of conidial binding to the cultured cells."
  - Discussion: the authors note that Kerr 2016 reported decreased phagocytosis and adhesion of ΔfleA conidia, and that "different cell lines used in the phagocytosis and adhesion assays might also be responsible".
  - The tool reported that the Figure 2 legend names a Student t test and that the paper also shows binding of purified FleA protein to host cells (Fig. 6), in a fucose-related way. I did not see the figure.
- Consequences for the spec: (1) F14 (c) and D7 are fair. "Probably no difference" is now supported by two fetches and by the authors' own Discussion. It is still "not read in the PDF". (2) F14 says "Kerr 2016 (b) has a nearby readout: macrophages take up ΔfleA conidia less well" but does not say that Sakai 2019 found **no** difference in phagocytosis. The two knockout papers disagree on uptake and on adhesion. Add this. It strengthens "the role is disputed" and does not change D7 option A.
- Question 4 can be answered now. I agree with the spec on this. Review 3 asked to wait for the full text; the Discussion sentences above remove most of that need.

**L-B. F27 lengths and "short".** Onygenales group lengths are 327 to 356 aa (*Onygena corvina* F1215351_005656-T1 is 356). The option (a) row says "about 330 aa". Of the 26 weak hits, 20 are 109 to 187 aa, 1 is 269 aa, and 5 are 400 to 454 aa (*Aphanocladium* x2, *Immersiella*, *Pseudomonodictys*, *Caloplaca*). "Weak, short hits" should be "weak hits, most of them short (20 of 26 under 190 aa)". This does not change any option.

**L-C. Pool S counts use Pfam 37.x labels; a Pfam 38.2 call would read 36 proteins.** F5 and F27 count 32 pool S proteins from the pool table's `pfam` column (Pfam 37.x). At Pfam 38.2 GA, 36 pool S proteins hit PF07938 (4 extra, all background, scores 27.3 to 37.6; 3 new genera). Option (a) would call about 36 proteins in about 31 genera (4.7%, still rare). Option (a2) is unchanged (the same 6). Add one line to F27 or Risk 4.

**L-D. The (a2) score floor is also fitted to the same data.** The spec says "any floor between 67 and 288 gives the same 6". So the floor sits in a gap that was seen in the same table. This is a second, smaller form of circularity beside the discovery-pass one. One clause in the (a2) row is enough.

**L-E. Small text points.**
- Owner question 9 ends with "Question 4 (label and role of FleA) can be answered now; the Sakai 2019 direction is still not read in the full PDF." This belongs in question 4.
- L0 lists "D1 to D5, D7, D9, D11, D12". D4b, D6, D8 and D10 have defaults but are not in L0 or in a question. State that their defaults stand unless the owner objects.
- Section 11 last line: "An independent review of this revision (items 1 to 3 of review 2) precedes L1." This is out of date.
- Header: "Revision 4 ... has no fourth review." Update to name this review.
- D6 says (a) and (a2) are gated by definition and asks whether an ungated (b) call should exist "beside them". Under a single-option answer to D1, D6 has nothing to decide. Say "only if the owner wants both a gated and an ungated call".
- 3.3 table: "pool U loses the M1 proteins (PF07938 only: 31 of 422,071)" is the option (a) number. Under (a2) it is 6.
- (c-ii) plan, two items not stated: (1) `scripts/sorting_hat/pfam_hmmsearch.sbatch` fetches every model named in the family table, active or not, so new rows need a new `hmmsearch` run on every proteome; (2) a new score column does not change module identity (`COLUMNS` is not hashed), so old and new tables of a module would share an identity. Neither matters for the recommended (c-i).

## 3. Answers to the review questions

**1. H1.** The spec's statements about the code are correct, line by line (section 1). (c-i) is buildable with no change to `src/`, `family_table.tsv` or `categories.yaml`. (c-ii) needs exactly what the spec lists, plus the two items in L-E. Risk 14 ("read, not run") is honest. I also did not run the code.

**2. H2.** The H1a-only rule is principled for a host-interaction lectin call: H1a is the only level whose evidence involves a host glycan or host cell. The reasons are stated, not asserted. The rule is applied correctly to PF07938 (H1a, common, so no ungated model), LysM, CVNH, FB_lectin and the H0 models. The CBM_14 alternative is presented fairly as D11 and question 7, with the statement that the rule alone would admit it under the revision 3 wording. The one gap is the "out of v1" set label (M-B).

**3. M1.** F27 and the whole-proteome claim are correct. The (a2) row describes the population and the circularity fairly and lists the limits (signal peptide from the pool definition, orthology not checked, no truth, no paper). The ranking (c), (a2), (b), (a) is labelled as reasoning only in 3.1, D1 and question 1. The overstatement in M-A and the fitted floor in L-D are the only fairness points.

**4. M2.** Sakai 2019 found no clear difference in conidial adhesion to A549, J774 and MH-S cells, and no difference in phagocytosis. The authors write that their data "did not support the FleA-dependence of conidial binding to the cultured cells". The spec's D7 and F14 wording ("probably no difference; not verified in the full text") is accurate and appropriately limited. The spec should add the phagocytosis result, which conflicts with Kerr 2016 (L-A).

**5. Spot checks.** 10 PMIDs match. F28 seed numbers match the committed file. `cluster_check.py` reproduces F10 exactly. F27 counts match; lengths need a small correction (L-B).

**6. New contradictions; owner questions.** Revision 4 introduced one misleading sentence (M-A) and one set label outside the rule (M-B). The rest is consistent. The owner questions (1, 2, 3a, 3b, 4 to 9) are answerable without extra context, except that question 1 needs the M-A fix to be read correctly. The defaults are sensible and agree with the data: no call in v1, the report outside the tool, H1a-only rule, no paracoccin model, label by family. None of the defaults rests on performance data, and the spec says so.

## 4. Required changes

**Must fix before the owner (text only, about four sentences):**

1. **(M-A)** Section 1 point 1 and owner question 1: say that the Onygenales restriction holds for pool S proteins only; without the pool S filter, PF07938 at score 100 or more is in 98 other genera.
2. **(M-B)** Put CBM_14 and Chitin_bind_1 in M2, or state why they are in neither M2 nor the report. Make the rule text and the 3.1 table agree.

**Minor, can fix later (before L4):**

3. (L-A) F14 (c): add the Sakai phagocytosis result and the Discussion sentence; state the conflict with Kerr 2016 on uptake and adhesion.
4. (L-B) F27 and the option (a) row: lengths 327 to 356; "most" weak hits short (20 of 26).
5. (L-C) F27 or Risk 4: 36 pool S proteins at Pfam 38.2 GA versus 32 by 37.x labels; (a2) unchanged.
6. (L-D) (a2) row: the floor is placed in a gap seen in the same data.
7. (L-E) Question 9 stray sentence; L0 list; section 11 last line; header; D6 scope; 3.3 pool U number; two (c-ii) plan items.

## 5. Ready for the owner?

**Yes, after the two must-fix items.** They are text changes and need no compute and no new review. The author can make them and send the spec. The main decision (D1) is presented with its options, its limits and an honest statement that no performance data supports any ranking.

## 6. What I did not do

- I did not read the Sakai 2019 PDF. The Discussion and Results sentences above come from a web-fetch summary of the publisher page, as for the author. Quotes are partial sentences.
- I did not run the sorting hat code. I read it.
- I did not re-run `hmmsearch` or SignalP. I used committed hit tables and the local pool tables.
- I did not check orthology of the 6 Onygenales proteins, their seed identity, or their signal peptides.
- I did not read the tests of the Pfam module beyond `grep` for module names (`tests/cellsurface_sorting_hat/test_data_and_scripts.py` line 57 lists the module names; a new module name needs a test change, as the spec says).
- Scratch MMseqs2 files are in `/tmp/r4clu_review4/` and are not committed. The new script and its output are in `analysis/lectin_class_spec/review4/`.
