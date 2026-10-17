# Review 1: lectin validation (branch `lectin-validation`, tip 036c2b2)

*2026-10-11. Reviewer: Claude Code (claude-opus-5-5), independent review. Scope: `docs/reports/2026-10-11-lectin-validation.md`, `analysis/lectin_validation/` (code, tables, `LOG.md`, `unverified_notes.md`), `analysis/lectin_validation/controls/`, `tests/lectin_validation/`. Governing documents: lectin spec (owner decisions block: option (c), report only), `docs/agent-tasks/COMMON-RULES.md`, hydrophobin validation spec. I did not edit the validation files, the spec or the scripts. All my runs were on the login node. No SLURM job was used. My scratch scripts are in `/tmp/lvrev_tools/` and are not committed.*

## Verdict

**Pass with changes.**

The numbers reproduce exactly. The truth rows are real proteins with correct sequences, species and taxon IDs. The report is mostly careful about circularity. But four statements in the report claim more than the data show, and two selection steps break the spirit of COMMON-RULES rule 3 without saying so. The text must change before the table is used as a report-grade table. No row needs to be removed. Section 7 gives the plain statement on use.

## 1. What I re-derived

| Check | Method | Result |
|---|---|---|
| Re-run of the pipeline | `run_all.sh` on a copy of the tree in `/tmp`, with the two Swiss-Prot raw files copied from the author's worktree (they are not in git) | Exit 0. `diff -rq` of `results/`, `controls/`, `out/clusters.tsv`, `out/truth_final.tsv`, `out/scan.hits.tsv` against the commit: **no difference**. Report script check: `mismatches 0`. |
| Sequence, species, taxon ID | UniProt REST for all 66 positives and all 166 reasoned negatives (232 rows); species rank from `nodes.dmp`; names from `names.dmp` | 232 of 232 sequences equal UniProt. 231 of 232 taxon IDs are the species-rank ancestor of the UniProt organism and match the `names.dmp` scientific name. The one exception is A0A1U7Q1Z0 (*Psathyrella*, genus 71729). It is disclosed (`unverified_notes.md` item 9). |
| Duplicate sequences | Exact-sequence grouping of `controls/sequences.faa` | 0 identical sequences under two accessions. |
| Q8NJT4 and Q8TGE0 | Pairwise comparison with Q4WW81 and Q2UNX8 | 313 of 314 and 307 of 310 residues match. They are the same genes (FleA, AOL). Rejection is right. `unverified_notes.md` item 15 ("I did not align them") can now be closed. |
| Leakage labels | Independent `mmseqs easy-search -s 7.5 -e 10` of all 729 sequences against the 1,627 seed sequences, then the same rule (identity 30% or more, E 1e-3 or less, target coverage 0.5 or more) | 0 label differences in 729 rows. The seed file is complete: the number of seed sequences per model equals `NSEQ` in the Pfam 38.2 HMM file for all 20 models. |
| PubMed sources | Abstracts of 13 PMIDs: 30947244, 34373510, 19100814, 8999822, 27911046, 28873459, 17153926, 22267486, 20724636, 21467214, 27174033, 16476564, 23695231 | Each quote in `curated_positives.tsv` is in the abstract. The T1 rows I read are honest T1 by spec 4.1 (glycan or chitin binding shown for the named protein). |
| Sub-GA beta-trefoil scores | `hmmsearch -T 0` with PF14200.13, PF00652.29 and PF07938.18 on the six missed beta-trefoil lectins | R4GRU4 24.9, R4GRU5 24.7, R4GRU6 24.6, B0D650 22.0, R4GRU9 21.2, Q7Z8V1 15.3 (GA 25.6). Same as `results/miss_subga.tsv`. No PF00652 hit at T 0. |
| Wilson interval | `lectin_metrics.wilson` on toy and real values | 0/10 gives 0.000 to 0.278. 10/10 gives 0.722 to 1.000. 19/30 gives 0.455 to 0.781. 0/430 gives 0.000 to 0.009. All correct. |
| Cluster-hit rule | All primary clusters checked for partial hits | No positive cluster is partly hit. "Any member hits" and "all members hit" give the same 19 of 30. |
| Status claim | Read `src/cellsurface_sorting_hat/status.py` (`status_from_measure`, lines 97 to 130) and `call_status.py` (`_entries`, lines 96 to 101) | Confirmed. See finding m9. |
| Tests | `/usr/bin/python3.12 -m pytest -q tests/lectin_validation` | 23 passed. |

## 2. Findings that must be fixed (text and disclosure; no row change needed)

### M1. The "seed-independent" recall counts the same clusters twice

Report lines 9, 92, 128 and 164 say that positives without a seed hit give 5 of 22 sequences and **3 of 14 clusters** (0.08 to 0.48), and that this is "the only seed-independent test of recall".

The leakage groups share clusters. In-reference positives are in 18 clusters and `none` positives are in 14 clusters, but there are only 30 primary clusters. Two clusters hold both kinds of row:

- Cluster `F2ST24`: J4URT3 and F2ST24 (`in_reference`) and F2SIK4 (`none`).
- Cluster `Q9P4R8`: Q9P4R8 and Q06100 (`in_reference`) and A0A2D0TCI3, Q206Z5, Q6WY08 (`none`).

These two clusters are two of the three "hit" clusters in the `none` group. Only one cluster that has **no** seed-similar member is hit: CNL (B2ZRS9). The correct cluster-level number is **1 of 12 clusters** (Wilson 0.01 to 0.35). The same error is in `results/family_recall.tsv`: the galectin `not_seed_similar` cluster (Q9P4R8) and the LysM `not_seed_similar` cluster (F2ST24) are the same clusters as in the `seed_similar` split.

`lectin_metrics.summarize` is correct for each group. The error is in the report: it treats per-group cluster counts as independent when the groups split clusters. Fix: define the seed-independent group at the cluster level (no member `in_reference`), and state 1 of 12. Add a test for it (section 6).

### M2. Most literature positives entered the candidate pool through the model under test (rule 3)

`make_acc_list.py:41` builds the candidate pool from Swiss-Prot entries with the keyword Lectin **or** a UniProt Pfam cross-reference to one of the 20 models. UniProt Pfam cross-references are InterPro matches of the same Pfam models at GA. `data/candidate_screen.tsv` shows that eight of the 13 added positives have no Lectin keyword and only a PF01476 cross-reference: F9XMT4, G2X4U8, H1UVH3, H1VSQ9, J4URT3, J5J7Y8, F2ST24, F2SIK4. G4N906 and Q00363 were added by name. B2ZRS9 and Q06100 have the keyword.

So for these eight LysM rows the recall (8 of 8) is fixed by how they were found. A LysM effector that the PF01476 model misses could not enter the pool. This affects:

- LysM family recall 11 of 11 (report line 90);
- `chitin_masking_h1b` 5 of 5 clusters and the T1 result 15 of 16 sequences (line 88);
- the claim that the two *T. rubrum* proteins are hit (lines 149, 161), since they entered because UniProt maps PF01476 to them. F2SIK4 is `none` for leakage but was selected by the model.

The report says that the F11 class mapping used Pfam hits (line 90). It does not say that candidate selection did. COMMON-RULES rule 3: "Do not pick positives because the module calls them." The rows are still true lectins, so they may stay. The report must state the selection route per row, and the recall of the literature rows must be labelled "by construction".

### M3. The `assumed` set is a construction check, not a measurement

`select_assumed.py:66` drops every protein with a UniProt Pfam cross-reference to any of the 20 models. The scan then uses the same models at GA. A hit is possible only where the UniProt Pfam release and Pfam 38.2 differ. The report says "partly by construction" (line 55). That understates it. It is almost entirely by construction.

The set also has no power. In all 38,451 reviewed fungal Swiss-Prot proteins, the 15 models give 97 protein-model hits (`results/swissprot_background.summary.tsv`), about 0.25%. A random sample of 497 with no filter at all would give about 1 hit. So a 0 of 497 result could not show anything, with or without the filter.

The Swiss-Prot background scan is the better form of this information and is already in the report. Fix: remove `assumed` from the false-hit tables (lines 10, 104, 159), or label it "construction check, not a false-hit rate". It also breaks the spirit of rule 3, for negatives.

### M4. The "no model" mapping of the beta-trefoil lectins hides a threshold miss

`measure.py:56` maps the UniLectin3D class "Laetiporus trefoil lectin-like" to `none`. The six rows of that class (R4GRU4, R4GRU5, R4GRU6, R4GRU9, B0D650, Q7Z8V1) are ricin-B-type beta-trefoil lectins. Five of them score 21.2 to 24.9 on RicinB_lectin_2 (GA 25.6), confirmed above. With the `none` mapping:

- the report says "every family that has a model among the 15 is hit in full", with "ricin-B-like 6/6" (line 90);
- the report says all 17 misses are in classes "without a lectin model among the 15" (line 147).

If these six rows are counted as the PF14200 family, PF14200 family recall is **6 of 12 sequences and 5 of 8 clusters**, not 6 of 6. The "model gap" wording in owner decision 4 (line 203) is not the right frame. The model exists. The GA cutoff excludes these lectins. Also, four of the five near-GA rows are *Boletus edulis* paralogs in one cluster (R4GRU5). "Five *Laetiporus*-type lectins" is **2 clusters** (R4GRU5, B0D650). Fix: state the mapping choice, give the PF14200 recall with and without these rows, and state the cluster count.

### M5. The T2 rule is applied in two ways

The report's T2 rule (line 24): "Swiss-Prot entry with an ECO:0000269 sugar-binding statement". G1X5Y2 (*Arthrobotrys oligospora* g276) is a Swiss-Prot entry with FUNCTION "Lectin that specifically binds to L-fucose (PubMed:35142828)", evidence ECO:0000269 (UniProt record, fetched 2026-10-11; also `data/candidate_screen.tsv`). It meets the rule. It was rejected as "T3" because the abstract of PMID 35142828 shows only a deletion phenotype (`build_truth.py:569`). I read the abstract. That is true.

But Q06100 (Cgl1) was accepted as T2 on its ECO:0000269 statement, while the abstract of PMID 8999822 does not say which of the two isolectins was expressed (`unverified_notes.md` item 4). The two rows are treated differently under the same rule. Either rule is defensible. The report must pick one and apply it to both rows. If G1X5Y2 enters, PF07938 has 6 positives in 4 clusters.

Also, report line 110 says D4AL26 "has no experimental evidence". UniProt gives D4AL26 PE 1 (protein level): it was found by mass spectrometry in the *A. benhamiae* secretome (PMID 21919205, Sriranganadane 2011). It has no **binding** evidence. The text should say that. It matters for the host-interaction question: a PF07938 protein of a dermatophyte is secreted in culture.

### M6. A false sentence about the UniParc row

Report line 24: "one UniParc ID, UPI000150AC21, has a sequence from the earlier spec work". `controls/rejected.tsv` rejects UPI000150AC21 ("UniParc ID without a sequence in the committed sequence file"), and it is not in `controls.tsv`. The sentence must go. (UPI000150AC21 is the *Sclerotium rolfsii* SRL entry; SRL is represented by A0A0M3KL30 and A0A0M3KL31.)

## 3. Minor findings

**m1. "in_reference" means domain similarity, not "in the seed".** Report line 164 says "The truth is mostly in the seeds". Only four positives are 98% or more identical to a seed sequence: I1RM03 (100%), Q2UNX8 (100%), Q00022 (99.2%), L8WGI4 (98.7%). Many LysM rows are `in_reference` because one 44-residue LysM domain is 35 to 59% identical to a bacterial, phage or algal seed domain (Ecp6 54.5% to LYTN_STAA8; Vd2LysM 40.0% to the phage Q9ZXH8; ELP1, ELP2, Slp1 to *Volvox* O65213). The galectin split is at the edge of both cutoffs: Q06100 is 31.3% at target coverage 0.508; Q206Z5, in the same 30% cluster, has no hit. Say "most positives have a domain at 30% identity or more to a seed domain". The rule is the one in COMMON-RULES rule 6, so the label is correct; the plain-English gloss is not.

**m2. Strata share clusters.** Stratum cluster counts add to 34, not 30. Clusters F2ST24 (h1a and h1b), P18891 (h1a and structure), Q2UNX8 (h1a and structure) and Q9P4R8 (structure and literature) span two strata. The `host_binding_h1a` 3 of 3 clusters is 3 clusters that each also hold seed-similar rows of other strata. Say so in section 5.1.

**m3. The negative-cluster floor uses the wrong rows.** Report line 177 and `controls/counts.md` give "most: 19" reasoned negative clusters per species. The owner decision of 2026-10-09 (spec F3) counts only N1 rows toward the 20-cluster floor. N1-only: *S. cerevisiae* 16, *A. fumigatus* 15. The conclusion (floor not met) does not change.

**m4. M3 hits on ER proteins are by construction.** The 18 `er_housekeeping` rows were selected by a UniProt M3 Pfam cross-reference (`build_truth.py`, M3 block). "The five M3 models hit all 18 ... They do what the spec wants of them" (line 107) should say this is expected from the selection. The same is true for 0 of 18 by the 15 models: M3 and M1/M2 are different Pfam families.

**m5. The enzyme false hits are hits on a seed protein.** All 8 alpha-galactosidases are one cluster (A1CBW8). P28351 is 100% identical to the PF00652 seed sequence AGALA_ASPNG/426-536; the other seven are 47 to 96%. The alpha-galactosidases were selected by name from the Pfam-or-keyword pool (`build_truth.py:429`). So "8 of 12 (0.39 to 0.86)" is one cluster with a known CBM13-type domain. The cluster-level 1 of 2 is the only fair number. The line "RicinB_lectin_2 PF14200 (2, at 28.8, below 30)" (line 106) is unclear; the two hits are above GA 25.6.

**m6. *T. rubrum* LysM as H1a needs a caveat.** PMID 30947244 shows binding to N-linked oligosaccharides of human skin glycoproteins **and** of fungal glycoproteins. LysM domains bind GlcNAc oligomers. The paper does not show a host-specific glycan. "That meets the H1a definition" (line 39) is a reading of the abstract. Add the fungal-glycoprotein binding to the quote, and leave the label to the owner (decision 1).

**m7. The report-script check is a plumbing check.** Only 16 of the 150 cells are non-empty (`results/report_script_check.tsv`). Both sides run HMMER 3.4 `--cut_ga` with the same Pfam 38.2 file. The check shows that the script parses and groups hits correctly. It is not an independent test of the models. The report says it covers hit sets only (line 195); add "16 non-empty cells".

**m8. Reproducibility depends on files outside git.** `run_all.sh` needs `data/raw/swissprot_fungi.tsv.gz` (not in git; the README says how to fetch it). `out/scan_swissprot.hits.tsv` is not made by `run_all.sh`; I did not find which command made it. The re-run works only with the author's worktree files. State the command that makes `scan_swissprot.hits.tsv`.

**m9. The `smoke` statement is correct.** The author did not re-read the code (line 180). I did. `status_from_measure` returns `smoke` unless there are 20 or more positives and negatives, 20 or more clusters of each (when recorded) and a half-width of 0.10 or less for both rates. `call_status._entries` caps any entry with leakage other than `none` at `smoke`. An entry may list several taxa, so a pooled entry is possible in code. A pooled entry would still be `smoke`: the cluster recall half-width is 0.16. A species without a positive is `unvalidated` (`n_pos < 1`). So "smoke at best" holds for any entry, pooled or not, from the size rule alone. The leakage cap is a second reason.

**m10. Output location.** COMMON-RULES says control sets go to `data/controls/<task-name>/`. This set is in `analysis/lectin_validation/controls/`. Move it or record the exception.

**m11. F28 values confirmed.** Report line 64 says the leakage method reproduces spec F28. My own search gives the same five values (FleA 81.7, SapL1 43.7, AOL 36.6, AAL 32.8, *A. oryzae* AOL 100.0).

## 4. Fungal lectins that the truth set missed

PubMed searches (2026-10-11): "(Avr4 OR Ecp6 OR LysM effector) AND chitin binding AND (fijiensis OR Dothistroma OR Rhizoctonia OR Trichoderma OR Penicillium OR Colletotrichum gloeosporioides OR Fusarium)" (13 hits); "fungal lectin AND (Aspergillus OR ...) AND (crystal structure OR glycan array OR recombinant) AND (fucose OR ...)" (67 hits, 30 read by title); a *Cordyceps* / *Aspergillus nidulans* / *Penicillium* / *Fusarium* / *Metarhizium* lectin search (113 hits, 38 read by title). Only abstracts were read. "Accession" means a UniProt entry I found by gene name; I did not check it against the paper's full text.

| Candidate | Species (phylum) | Evidence in the abstract | PMID | Accession | Note |
|---|---|---|---|---|---|
| CMRBL | *Cordyceps militaris* (Ascomycota) | Gene CCM_03832; recombinant; glycan array, type II LacNAc | 36632960 | G3JGU5 (gene CCM_03832; UniProt Pfam PF00652) | Would be the **only true positive for PF00652**, which now has 0 positives and 8 enzyme hits. High priority. |
| BbLec1 | *Beauveria bassiana* (Ascomycota) | "BbLec1 could bind to chitobiose and chitin" | 34649678 | Probably J4UPT0 (BBA_03816, 147 aa, the only PF07367 protein of strain ARSEF 2860) | Rejected for "no verifiable accession". J4UPT0 is a candidate; the gene ID must be checked in the full text. |
| CmLec4 | *C. militaris* | Binds host (silkworm) storage protein 2; knockout phenotype | 35718153 | not found by gene name | Insect-host binding; H-flag question. |
| CmRlec | *C. militaris* | Ricin B-type; "novel chitin-recognizing lectin" | 37073122 | not checked | |
| RsLysM | *Rhizoctonia solani* (Basidiomycota) | "displayed chitin-binding affinity" | 31076860 | not found by gene name | First Basidiomycota LysM effector. |
| Cg2LysM | *Colletotrichum gloeosporioides* | "showed chitin-binding activity" | 36876102 | not found by gene name | |
| TAL6 | *Trichoderma atroviride* | Binds polymeric chitin, not oligomers | 23289754; 31608044 | not found by gene name | Weak binder; not a pathogen. |
| MfAvr4 | *Pseudocercospora fijiensis* | Chitin binding, protects cell walls | 20368413 | not checked | CBM_14. |
| DsAvr4, CaAvr4, CbAvr4 | *Dothistroma septosporum*, *Cercospora apii*, *C. beticola* | "All, however, bind chitin" | 25845605 | not checked | CBM_14. |

Most new candidates are LysM or CBM_14. They would add Ascomycota clusters, but most would be `in_reference` and would also enter through a Pfam cross-reference. They do not reduce the circularity. CMRBL and BbLec1 add tests of PF00652 and PF07367 on non-structure rows.

Also seen, not a positive: J4KR73 (*B. bassiana*, "Peptidase M64, IgA", UniProt Pfam PF07938 only). If the name is right, it could be the first PF07938 look-alike negative that spec L3 asks for. It needs a paper before it can be a reasoned negative.

## 5. The five owner decisions in the report

| # | Report question | Right question? | Comment |
|---|---|---|---|
| 1 | Add *T. rubrum* LysM1/2 as H1a | Yes | Add caveat m6 (binding to GlcNAc cores of host and fungal glycoproteins). Under owner decision D11 (H1a only) this matters. |
| 2 | Accept `enzyme_lectin_like_domain` as an apart class | Yes | Better frame: these are true CBM13 domain hits on an enzyme. The question is whether the report labels proteins or domains. |
| 3 | Tectonin PF19193 and Fve PF09259 as M2 | Yes | FIP-Fve is an immunomodulatory protein in the Ig-like class; its lectin status is weaker than that of Tectonin 2. Decide separately. |
| 4 | Near-GA beta-trefoil family | Partly | It is a GA threshold question for an existing model (M4), not a model gap. Two clusters, not five lectins. |
| 5 | Use of the tables later | Yes | |

Questions the report does not raise:

6. **Candidate pool gated by the models (M2).** Accept with disclosure, or rebuild the literature positives from a search that does not use Pfam cross-references?
7. **`assumed` negatives (M3).** Drop them from the false-hit tables?
8. **T2 rule (M5).** Include G1X5Y2 (Swiss-Prot ECO:0000269), or tighten T2 and re-check Q06100?
9. **PF07938 look-alike negatives (spec L3).** Build them before any statement about PF07938 specificity? J4KR73 and Q4WMH4 are starting points.
10. **Leakage for short domains (m1).** Keep the COMMON-RULES 30% rule for LysM-sized domains, or add a near-identity tier (for example 90%) to separate "in the seed" from "similar to a seed domain"?
11. **New candidates (section 4).** Add CMRBL (G3JGU5) and BbLec1 (after the gene check)?

## 6. Tests

All 23 tests pass. The truth-table tests are useful: required columns, unique accessions, source and quote on every row, valid leakage values, a reason on every reasoned negative, sequence and cluster coverage, taxon names against `names.dmp`, disjoint positive and negative sets. The metric tests check Wilson at 0/n, n/n, 5/10 and n = 0, floors, and the "any member" cluster rule.

Gaps:

1. No test of `measure.py`. All grouped numbers in the report come from it, and M1 is an aggregation error that no test sees. Add a toy truth table where one cluster holds an `in_reference` and a `none` row, and assert that the seed-independent group is defined at cluster level.
2. `test_taxon_ids_and_species_match_names_dmp` checks the name only, not the rank. A genus (A0A1U7Q1Z0, *Psathyrella*) passes. Add a rank check with an explicit allow-list.
3. No test that every positive's selection route is recorded (M2).
4. No test that identical sequences do not appear under two accessions. (None exist now; I checked.)
5. `test_curator_tuned_date` hard-codes the curator name. It will fail when a person adds a row.

I found no bug in `wilson`, `hit_accessions` or `summarize` with toy cases.

## 7. Can the set be used as a report-grade table?

**Yes, after M1 to M6 are fixed in the text.** The rows, sequences, taxa, clusters, leakage labels and hit tables are correct and reproduce exactly. The fixes change sentences and some headline numbers (seed-independent recall 1 of 12 clusters; PF14200 family recall 6 of 12 with the beta-trefoil rows), not the rows.

It is **not** a calibration set. It cannot give a specificity for PF07938 (no relative was chosen as a negative), and the recall is mostly a measure of how the positives were selected and of seed similarity. The owner review gate of COMMON-RULES still applies.
