# M5: SOWgp-like calls (an ortholog call and a repeat-architecture call)

*2026-10-09. Revision 2 (after `2026-10-09-sowgp-like-calls-design-review-1.md`: 2 blockers, 9 major, 6 minor, all addressed below). Author: Claude Code (claude-sonnet-5-5). Milestone M5 of `docs/ROADMAP-2026-10-08-DRAFT.md`. Owner decision of 2026-10-08: SOWgp-like means **both** an ortholog call from the unit HMM and a repeat-architecture call. Status: draft for independent review. No code is written. Every number below is copied from a file named in the text or from the review's checks.*

## 1. The question

The project asked for antigenic candidates in *Coccidioides* "building on SOWgp and PRA1/PRA3" (`docs/reports/2026-09-27-coccidioides-antigen-findings.md`), and for proteins in other species that look like SOWgp. SOWgp (gene CIMG_04613, spherule outer wall glycoprotein) is a secreted protein with an unrepeated N-terminus and a tandem array of a 47 aa unit (`analysis/cocci_repeats/REPORT_2026-09-29_sowgp_repeat_structure.md`). It is a tandem-repeat protein like FLO11 but with a different composition: proline and aspartate rich, with cysteine at 6.5% (RS SOWgp: Pro 15.7%, Asp 13.6%, Cys 6.5%), not Ser/Thr rich. The SOWgp paper (Hung et al. 2002, PMID 12065484, doi:10.1128/IAI.70.7.3443-3456.2002) reports recombinant SOWgp that binds laminin, fibronectin and collagen IV, and a deletion strain that loses part of its ECM binding and is less virulent. So adhesion evidence exists for SOWgp itself. It does not exist for the look-alikes.

Today the tool flags SOWgp in *C. immitis* RS only because SOWgp tuned the repeat and antigen modules (`docs/reports/data/sorting_hat/panel.tsv`, rows `XP_001245172.2`, calls `tandem_repeat_protein` and `cocci_specificity_rank_top15`, tuning `tuned`). It has no call that means "SOWgp-like". Outside RS, a SOWgp-type protein appears as a generic `tandem_repeat_protein` plus a signal peptide. Nothing separates a Pro/Cys-rich array from a Ser/Thr array.

## 2. What exists (facts)

| Item | Where | Numbers |
|---|---|---|
| Unit profile `sowgp_unit.hmm` (34 aligned units, model length 47; no GA line) | `analysis/cocci_repeats/sowgp_unit.hmm`, `34_sowgp_unit_hmm.py` | Search over 6,313 proteomes (7 long-read *Coccidioides*, 493 pangenome strains, 5,813 Fungi_5k): **1,293 domain hits (487 proteins, 455 proteomes) at domain score 26 or more; all are *Coccidioides*.** The best hit elsewhere has domain score 19.7; no hit between 26 and 40 in either group. All 307 non-*Coccidioides* rows are single-domain hits. Null (shuffled proteomes): 0 hits. Every *Coccidioides* protein with any hit has a domain at 26 or more (487 of 487); 455 further *Coccidioides* domain hits score 18.3 to 26 (partial units). |
| Column bug in the 2026-09-30 search | `37_unit_search.sh` (awk prints `$6`), `36_unit_calibrate.py` (reads `f[5]`) | In `--domtblout`, field 6 is `qlen`, not the sequence score. The column `full_score` is 47 in all 2,055 rows. **Every score in the 2026-09-30 report is a domain score. The sequence score was never saved.** Do not use `full_score`. |
| Sensitivity floor of the profile | `REPORT_2026-09-30_sowgp_architecture_relatives.md` section 2.1; `36_unit_calibrate.py` | A unit counts as found at domain E-value 1.0 or less (`--floor-e`, default 1.0; domZ 6e7; about 21 bits with the model's Forward statistics). Arrays of four units at 60% identity: all found; 50%: 100% (95% all copies); 40%: 92% (6% all copies); 30%: 1%. The table's `median_best_score` at 40% identity is 26.0, so **at a 26-bit cutoff about half of the four-unit arrays at 40% identity have no unit called**. No indels simulated. |
| Architecture look-alikes in Onygenales | same report, section 3; `37_unit_search.sh` phase C runs `14_repeat_detect_general.py` | The scan ran on **78 Onygenales proteomes** (Fungi_5k Onygenales plus the 7 `cocci_longread` entries; two genomes appear twice: RS as `Coccidioides_immitis_RS` and `CimmitisRS_FungiDB`, Silveira likewise). Filter Pro+Cys over 15% and Cys at least 4%: 83 proteins in 50 proteomes (repeat14 only); 74 after removing 9 SOWgp proteins (period 47, 9 proteomes of 7 strains); 27 repeat periods, none at 47 aa. The 74 include duplicates from the twice-annotated genomes. Unreviewed. |
| Composition of known cases | 2026-09-30 report section 4; review checks | SOWgp RS: Pro 15.7, Asp 13.6, Cys 6.5. BAD1 (A4D962, 1,146 aa): Pro 3.8, Cys 8.0, sum 11.8%: **fails the rule**. The BAD1-like period-43 family in *Coccidioides* (Pro+Cys 10.4 to 12.1%) also fails. In RS the rule calls `XP_001248770.2` (period 27; Pro 3.5, Cys 12.1) and `XP_001248220.2` (period 43; Pro 9.2, Cys 16.1): Pro-poor, Cys-rich. |
| Curated class 2a repeat candidates | `docs/reports/2026-09-27-cocci-repeat-surface-proteins.md`, `class2a_candidates.tsv` (41), `class2a_candidates_general.tsv` (58) | 20 Pro/Cys-rich, 10 Ser/Thr-rich, 11 other (reported). **Only 12 of the 20 pass Pro+Cys over 15% and Cys at least 4%**; 8 have Cys 0 to 3.6% (labelled by `03_repeat_surface_candidates.py`, which has no cysteine floor). 4 of the 20 are SOWgp orthologs (CIMG_04613, CIB10637_003943, CIB10992_003451, QVM09276.1). The 20 fall in about 6 homology clusters (5 RS loci and 1 family without an RS gene model). |
| Repeat detectors | `02_repeat_profile.py`, `14_repeat_detect_general.py`; modules `repeat02`, `repeat14` | call = period above 0, coverage at least 0.25, copies at least 2.5; `composition()` in script 14 gives `pct_pro`, `pct_cys`, `pct_ser_thr`. In RS the 2026-09-30 rule gives 7 calls with either detector and 6 with repeat14 alone; `XP_001239868.2` (242 aa; Pro 22.3, Cys 5.0) is called by repeat02 only. |
| Spherule expression (RS) | `analysis/cocci_spherule/spherule_surface_table.tsv.gz` | 9,757 genes; 564 up at 48 h; SOWgp +9.09 log2 fold change (top) |
| Antigen ranking (RS) | `analysis/cocci_antigens/cocci_antigen_ranking.tsv` | module `antigen_lookup`; NOT CALIBRATED |
| BAD1 control | `38_unit_bad1_control.py` | whether it was run is not verified (roadmap section G) |
| Source of the 15% and 4% thresholds | none | No file records how they were chosen. |

## 3. Design

Two calls with different jobs. The names say what each call tests (homology; composition plus repeat). Adhesion evidence exists for SOWgp only (PMID 12065484); neither call is an adhesion claim for other proteins.

### 3.1 `sowgp_ortholog` (evidence call): SOWgp unit homologue

- **Meaning:** the protein has at least one domain hit to the SOWgp unit profile at domain score 26 or more. It is a family-membership call, not an architecture call. Outside *Coccidioides* it is expected to be `not_called` (section 2: best hit elsewhere 19.7). A `called` outside *Coccidioides* would be a finding.
- **Module `sowgp_unit`:** one row per protein with at least one row in the `hmmsearch --cut_ga` table: `hit` (1), `n_units` (number of domain rows in that table), `best_domain_score`, and `seq_score` (reported, not used in the decision). It is built like `hydrophobin_relaxed` (derived `.hmm`, `hmmsearch --cut_ga`, a wrapper that checks the table came from that file, params and file hash in the module identity). `model_info()` and `check_table()` are reused. The row builder is new: `relaxed_rows()` keys on the sequence score and keeps one row per protein, while this module needs domain rows (count and best score).
- **Derived `.hmm` file:** `sowgp_unit.hmm` plus a line `GA -1000.00 26.00;` before `STATS` (sequence cutoff -1000, domain cutoff 26). The domain score alone decides. This matches the calibration, which was done on domain scores. HMMER 3.4 accepts this line and writes the `--cut_ga` option line and `# [ok]`, so `pfam.parse_domtblout` and `check_table` accept the table. (`GA 26.00 26.00;` is an alternative; on 1,000 synthetic arrays both give the same 264 proteins. It can in principle drop a protein with one domain at 26 or more whose sequence score is lower, which is why `-1000` is chosen.)
- **Hit decision:** taken from the rows `--cut_ga` kept, never from the printed score (the printed value is rounded: two synthetic proteins print as 26.0 and were dropped).
- **Why the sequence score is not used:** on synthetic arrays a sequence-score decision calls arrays of weak units that no single domain supports (36% identity, 8 copies: sequence score at least 26 calls 34 of 40, domain score at least 26 calls 7; 33% and 12 copies: 20 versus 8). It has no calibration, because the sequence score was never saved (section 2).
- **Cutoff:** 26 bits is a **choice inside an empty interval, not a calibrated value**. Any cutoff between 19.7 and 40 gives the same decisions on the 6,313 proteomes. The data do not pick a value inside it. The reason for 26 is the median best score at 40% unit identity (26.0). Both statements go in the provenance file, with the column bug of section 2.
- **Reach:** at 26 bits about half of four-unit arrays at 40% identity are missed, so a relative below about 50% unit identity can be invisible. M5a reruns `36_unit_calibrate.py` with a bit-score floor of 26 and the report quotes that table instead of the E-value-floor table.
- **Gate:** none. The call does not require a signal peptide. (SOWgp has one; the call is about homology.) The report shows R0 beside it.
- **Category:** one category: cell wall and adhesion. The antigen angle is shown in the candidate report (section 3.3), not as a second category. M2 decides whether a call can have two categories; until it does, this spec uses one.
- **Composites:** `sowgp_ortholog` joins the mechanism lists of `other_not_surface` and `other_surface_no_mechanism`, so a split gene-model fragment that carries units but no repeat call does not fall into "other, no mechanism" (decision O1). It does not feed `cell_wall_adhesion_candidate` or `surface_attachment_candidate`.

### 3.2 `pro_cys_array_protein` (evidence call): Pro+Cys-rich tandem array

- **Meaning:** a protein that either repeat detector calls as a tandem repeat, with Pro+Cys over 15% and Cys at least 4% (whole protein). **Proline can be low.** The sum lets cysteine carry the call: a protein with Pro 3.5% and Cys 12.1% is called. The call therefore means "cysteine-rich array, with proline-rich ones included". It does not mean "SOWgp-like composition". Decision O2 asks whether to add a proline floor.
- **Not covered:** BAD1-type arrays (Pro 3.8%, Cys 8.0%, sum 11.8%) fail the rule. The call does not cover them. BAD1 is **not** a positive (decision O3).
- **Module `composition`:** one row per protein with `pct_pro`, `pct_cys`, `pct_ser_thr`, `n_residues`, and `pro_cys_pct` (Pro+Cys from residue counts, unrounded), for the **whole protein**. The 83 count in the 2026-09-30 report used rounded values; the strict and the non-strict form both give 83 on the Onygenales file. A repeat-region version of the columns is added only if M5b shows that the whole-protein version fails to reproduce the 83 proteins. Thresholds are in the call, not the module, so a threshold change changes `config_sha256` only.
- **Call (exact YAML):**

```yaml
- name: pro_cys_array_protein
  expr:
    and:
      - or: [{call: repeat02}, {call: repeat14}]
      - test: {module: composition, field: pro_cys_pct, op: ">", value: "$pro_cys_pct_min"}
      - test: {module: composition, field: pct_cys, op: ">=", value: "$cys_pct_min"}
```

  with `pro_cys_pct_min: 15` and `cys_pct_min: 4` added to `thresholds:`. `{call: X}` reads the `call` column of **module** X (`engine.py` line 192), so `repeat02` and `repeat14` are written as modules. The expression reads three modules and has no `ref` node, so `call_eligible()` accepts it and `calibrate truth --call-status` works.
- **Detector choice:** the 83 come from `repeat14` alone. M5b reproduces the 83 with `repeat14` only, and then reports the extra proteins that `repeat02` adds (RS: `XP_001239868.2`). If the extra calls are not wanted the call changes to `repeat14` only. Decision O2b.
- **Gate:** none (a signal peptide is shown beside it). A Ser/Thr-rich array such as FLO11 is `not_called`.
- **Category:** cell wall and adhesion.
- **Feeds:** `surface_attachment_candidate` already reads `tandem_repeat_protein`, so a Pro/Cys array is already inside the broader call. The sub-label of that call gains `pro_cys_array_protein` as a refinement (decision O1).

### 3.3 Candidate report (not a call)

A table per run: for each proteome, proteins with `pro_cys_array_protein` and a signal peptide, with period, units, composition, and (RS only) spherule fold change, antigen percentile and `sowgp_ortholog`. It answers "other plausible proteins like it in *Coccidioides*" and "anything that looks like it in other species". It is a table of candidates, not a validated list. The 20 class 2a Pro/Cys-rich candidates and the 74 look-alikes appear here as **T4 review outcomes** (assistant judgement, no paper), beside the measurement and never in it, with cluster counts (about 6 clusters for the 20), not protein counts.

## 4. Measurement

Rule from the hydrophobin work (`2026-10-08-hydrophobin-validation-design.md` line 147): T4 labels are never added to the truth set and never used for sensitivity. `calibrate truth` takes **one species per status entry** and refuses truth proteins outside it (`cli.py` lines 461 to 495). Each row below is its own entry and its own run.

| Call | Entry (taxon, proteome, run directory) | Positives (IDs in the run's namespace, with cluster) | Negatives | Leakage |
|---|---|---|---|---|
| `sowgp_ortholog` | *C. immitis* RS (`_workdir/sorting_hat/`, RefSeq IDs; `XP_001245172.2` is SOWgp) | SOWgp orthologs with an external label: PMID 12065484 for SOWgp itself; other orthologs from `05_sowgp_pangenome.py`; one homology cluster | all other RS proteins | `in_reference` (units come from all full-length pangenome SOWgp copies, `34_sowgp_unit_hmm.py`) |
| `sowgp_ortholog` | *C. posadasii* Silveira (a proteome file named in M5a and added to the run set; none exists in `_workdir/sorting_hat/` now) | its SOWgp ortholog | all other proteins | `in_reference` |
| `sowgp_ortholog` | S288C, *C. albicans*, *B. dermatitidis* ER-3 | none | all proteins (specificity only; `build_measure` allows entries without positives) | `in_reference` |
| `pro_cys_array_protein` | RS | SOWgp (external label only). BAD1 is not a positive. | proteins with an external label that is not a Pro/Cys array: none exist yet in RS; the 10 Ser/Thr class 2a negatives are outside RS, so they cannot be used here | `tuned_on_truth` is plausible but unverified: the source of the thresholds is not recorded |
| `pro_cys_array_protein` | S288C and *C. albicans* | none | repeat-truth proteins that are not Pro/Cys-rich arrays; specificity only | as above |

- **Truth IDs:** the analysis tables use FungiDB IDs (`CIMG_04613-t26_1-p1`); the tool's RS proteome uses RefSeq IDs. The truth tables need a mapped ID column. M5c builds it from the protein FASTA by exact sequence match and lists unmatched IDs.
- **The Ser/Thr negatives (10, in six non-RS proteomes: CiB10637 1, CiB10992 1, VFC140 2, Cpos1038 4, Cpos3700 1, Silveira 1)** need their proteomes in the run set. None is in `_workdir/sorting_hat/`. If M5c does not add them, they are reported as T4 only.
- **Fungi_5k specificity bound** is a report-level statement only: 0 domain hits at 26 or more in 5,813 Fungi_5k proteomes. It is not a status entry, because the tool was not run on those proteomes as separate species entries.
- **Status:** both calls are capped at `smoke` (leakage `in_reference` or `tuned_on_truth`), and SOWgp is one homology cluster, so neither can reach `estimated` (20 clusters per class needed). A status from these positives measures self-consistency only. The report says so and states no sensitivity or specificity for new proteins.
- `sowgp_ortholog` reads one module and gets a module status (module `sowgp_unit`). `pro_cys_array_protein` reads three modules and gets a call status.

## 5. Decisions for the owner (stop points)

O2, O3 and O4 must be decided before M5b starts.

| # | Question | Default |
|---|---|---|
| O1 | Does `pro_cys_array_protein` appear as a sub-label of `surface_attachment_candidate`? Does `sowgp_ortholog` join the mechanism lists of the `other_*` calls (and feed no other composite)? | sub-label yes; `sowgp_ortholog` joins the mechanism lists only |
| O2 | Keep Pro+Cys over 15% and Cys at least 4% as is (Pro can be low; the meaning in 3.2 says so), or add a proline floor and recount? | keep, with the stated meaning; the source of the thresholds is not recorded |
| O2b | Detector: `repeat02 or repeat14` (default) or `repeat14` only (reproduces the 83 exactly)? | either; M5b reports both |
| O3 | BAD1 is not a positive and BAD1-type arrays are outside the call (fails the rule, 11.8%). Or change the rule (for example Cys at least 4% and (Pro+Cys over 15% or Cys over 7%)) and record it as a new rule; then the 83 and 74 counts no longer apply. | drop BAD1; state the limit |
| O4 | The 20 class 2a candidates and the 74 look-alikes: reported as T4 outcomes beside the measurement, never as positives. | reported only |

## 6. Risks

1. Positives and thresholds come from the same proteins: any status is `smoke` with leakage, and says nothing about proteins unlike SOWgp.
2. The unit HMM has a sharp floor. A diverged relative gives `not_called` and the report must not imply absence. The 26-bit cutoff sits in an empty interval; the data do not choose a value in it.
3. A composition rule does not see layout. Pro/Cys-rich or Cys-rich proteins with a repeat but no relation to SOWgp are called (the 74 look-alikes are such cases, unreviewed).
4. SOWgp's gene models vary: in the pangenome 39 of 488 strains have no gene model but normal read depth, and 3 of 6 long-read loci split into two gene models (`REPORT_2026-10-01_sowgp_*`). A protein-level call misses split models; the mechanism-list membership (O1) keeps a fragment out of "other, no mechanism" but does not recover it as a SOWgp call.
5. Antigen and spherule columns exist for RS only.
6. The sequence score of the 6,313-proteome search was not saved; no sequence-score calibration exists.

## 7. Work plan

| Task | Content | Check |
|---|---|---|
| M5a | `sowgp_unit` module and derived `.hmm` (GA `-1000 26`); tests first, including one synthetic multi-unit fixture (RS-like array: 4 units, 22.5-bit unit not counted) and the derived file against `--cut_ga`; rerun `36_unit_calibrate.py` with a bit-score floor of 26; add a *C. posadasii* proteome to the run set; run on RS, *C. posadasii*, *B. dermatitidis*, S288C | RS calls CIMG_04613 with `n_units` 4 (`XP_001245172.2`: sequence score 324.8; domains 89.1, 95.6, 94.4, 66.7, 22.5, -3.0, -6.5); S288C and ER-3: 0 rows |
| (stop) | Owner decides O1 to O4 | |
| M5b | `composition` module; reproduce the 83 with `repeat14` alone on the 78 proteomes, then the 74 (removing 9 SOWgp), then report what `repeat02` adds; golden files read line by line | the counts reproduce or the difference is explained |
| M5b2 | the `categories.yaml` edit (calls, thresholds, mechanism lists), in its own task after the M5b counts are accepted | tests pass; config hash change noted |
| M5c | Truth tables with mapped IDs and status files (one entry per species, section 4) | status files `smoke`; unmatched IDs listed |
| M5d | Candidate report on the scan proteomes (section 3.3) | table per proteome; SOWgp row in RS |
| M5e | Docs: `docs/CLASSES.md` entries, paper notes | test passes |

Each task is its own commit. An independent review of this spec (rev 2) and of the plan precedes M5a.
