# M5: SOWgp-like calls (an ortholog call and a repeat-architecture call)

*2026-10-09. Revision 1. Author: Claude Code (claude-sonnet-5-5). Milestone M5 of `docs/ROADMAP-2026-10-08-DRAFT.md`. Owner decision of 2026-10-08: SOWgp-like means **both** an ortholog call
from the unit HMM and a Pro/Cys repeat-architecture call. Status: draft for independent review. No code is written. Every number below is copied from a file named in the text.*

## 1. The question

The project asked for antigenic candidates in *Coccidioides* "building on SOWgp and PRA1/PRA3" (`docs/reports/2026-09-27-coccidioides-antigen-findings.md`), and for proteins in other species that
look like SOWgp. SOWgp (gene CIMG_04613, spherule outer wall glycoprotein) is a secreted protein with an unrepeated N-terminus and a Pro/Cys-rich tandem array of a 47 aa unit
(`analysis/cocci_repeats/REPORT_2026-09-29_sowgp_repeat_structure.md`). It is a tandem-repeat protein like FLO11 but with inverted composition: proline and cysteine rich, not Ser/Thr rich
(`docs/TOOL-ARCHITECTURE.md`, section 2).

Today the tool flags SOWgp in *C. immitis* RS only because SOWgp tuned the repeat and antigen modules (`data/sorting_hat/panel.tsv`, rows `tuned`). It has no call that means "SOWgp-like". Outside RS,
a SOWgp-type protein appears as a generic `tandem_repeat_protein` plus a signal peptide. Nothing separates a Pro/Cys array from a Ser/Thr array.

## 2. What exists (facts)

| Item | Where | Numbers |
|---|---|---|
| Unit profile `sowgp_unit.hmm` (34 aligned units, model length 47) | `analysis/cocci_repeats/sowgp_unit.hmm`, `34_sowgp_unit_hmm.py` | Search over 6,313 proteomes (7 long-read *Coccidioides*, 493 pangenome strains, 5,813 Fungi_5k): all 1,293 hits at score 26 or more are *Coccidioides*; the best hit elsewhere scores 19.7; no hit between 26 and 40 in either group. Null (shuffled proteomes): 0 hits. |
| Sensitivity floor of the profile | `REPORT_2026-09-30_sowgp_architecture_relatives.md` section 2.1 | arrays of four units at 60% identity: all found; 50%: 100% (95% all copies); 40%: 92% (6% all copies); 30%: 1%. No indels simulated. |
| Architecture look-alikes in Onygenales | same report, section 3 | composition filter Pro+Cys over 15% and Cys at least 4%: 83 proteins in 50 proteomes; 74 after removing SOWgp; 27 repeat periods, none at 47 aa; none shares the SOWgp unit. Unreviewed. |
| Curated class 2a repeat candidates | `docs/reports/2026-09-27-cocci-repeat-surface-proteins.md`, `class2a_candidates.tsv` (41), `class2a_candidates_general.tsv` (58) | 20 Pro/Cys-rich, 10 Ser/Thr-rich, 11 other (reported). Only CIMG_04613 of the Pro/Cys-rich RS loci is spherule-induced. |
| Repeat detectors | `02_repeat_profile.py`, `14_repeat_detect_general.py`; modules `repeat02`, `repeat14` | call = period above 0, coverage at least 0.25, copies at least 2.5; the table also carries `pct_pro`, `pct_cys`, `pct_ser_thr` (computed by `composition()` in script 14; whether over the whole protein or the array is checked in task M5b) |
| Spherule expression (RS) | `analysis/cocci_spherule/spherule_surface_table.tsv.gz` | 9,757 genes; 564 up at 48 h; SOWgp +9.09 log2 fold change (top) |
| Antigen ranking (RS) | `analysis/cocci_antigens/cocci_antigen_ranking.tsv` | module `antigen_lookup`; NOT CALIBRATED |
| BAD1 control | `38_unit_bad1_control.py` | whether it was run is not verified (roadmap section G) |

## 3. Design

Two calls with different jobs. Neither is an adhesion claim.

### 3.1 `sowgp_ortholog` (evidence call): SOWgp unit homologue

- **Meaning:** the protein has at least one hit to the SOWgp unit profile at a frozen score. It is a family-membership call, not an architecture call. Outside *Coccidioides* it is expected to be `not_called` (section 2: no hit above 19.7). A `called` outside *Coccidioides* would be a finding.
- **Module `sowgp_unit`:** one row per protein: `hit` (1 when the best domain score is at least the frozen cutoff), `n_units` (domain hits at or above the cutoff), `best_score`. Built like `hydrophobin_relaxed`: a derived `.hmm` file with the GA line set to the frozen cutoff (domain cutoff equal to the sequence cutoff; to be checked in M5a that the sequence score of a multi-unit array does not change the decision, because arrays add up in the sequence score and not in the domain score), `hmmsearch --cut_ga`, a wrapper that checks the table came from that file, params and file hash in the module identity.
- **Cutoff:** 26 bits, from the calibration above (all 1,293 hits at 26 or more are *Coccidioides*; no hit between 26 and 40). The choice and its evidence are recorded in the provenance file. The reach is stated: a relative below about 40% unit identity is invisible.
- **Gate:** none. The call does not require a signal peptide. (SOWgp has one; the call is about homology.) The report shows R0 beside it.
- **Category:** cell wall and adhesion (listed) and antigen (cross-listed, because SOWgp is a top antigen). Not part of `cell_wall_adhesion_candidate`.

### 3.2 `pro_cys_array_protein` (evidence call): Pro/Cys-rich tandem array

- **Meaning:** a tandem-repeat protein (either detector calls it) whose composition is Pro+Cys over 15% and Cys at least 4%. These are the thresholds of the 2026-09-30 analysis
  (`REPORT_2026-09-30_sowgp_architecture_relatives.md` section 3). The cysteine floor is deliberate: a sum without a floor lets cysteine-free proteins through.
- **Module `composition`:** one row per protein with `pct_pro`, `pct_cys`, `pct_ser_thr`, `n_residues`, computed from the FASTA, for the **whole protein** (as in the 2026-09-30 analysis, so the 83 and 74 counts can be reproduced). A second set of columns for the repeat region is added only if M5b shows
  that the whole-protein version fails to reproduce the 83 proteins. Module identity: thresholds are in the call, not the module, so a threshold change changes `config_sha256` only.
- **Call:** `{call: tandem_repeat_protein}` and `{test: composition.pct_pro + pct_cys > 15}` and `{test: composition.pct_cys >= 4}`. The engine has `test` nodes on one field. A derived field `pro_cys_pct` (Pro+Cys) is written by the module so that a single `test` reads it.
- **Gate:** none (a signal peptide is shown beside it). A Ser/Thr-rich array such as FLO11 is `not_called`, which is the point of the call.
- **Category:** cell wall and adhesion.
- **Feeds:** `surface_attachment_candidate` already reads `tandem_repeat_protein`, so a Pro/Cys array is already inside the broader call. The sub-label (basis) of that call gains `pro_cys_array_protein` as a refinement (decision O1).

### 3.3 Candidate report (not a call)

A table per run: for each proteome, proteins with `pro_cys_array_protein` and a signal peptide, with period, units, composition, and (RS only) spherule fold change, antigen percentile and `sowgp_ortholog`. This answers "other plausible proteins like it in *Coccidioides*" and "anything that looks like it in other species". It is a table of candidates, not a validated list.

## 4. Measurement

| Call | Truth | Leakage | Expected status |
|---|---|---|---|
| `sowgp_ortholog` | positives: SOWgp orthologs (from `05_sowgp_pangenome.py` and the anchored family search; the pangenome table lists the member proteins per strain); negatives: all other proteins of *C. immitis* RS and *C. posadasii* Silveira (assumed), and the non-*Coccidioides* proteomes of the 6,313 set (specificity bound: 0 hits at 26 or more in 5,813 Fungi_5k proteomes) | `in_reference` (the profile was built from SOWgp units) | `smoke` at best |
| `pro_cys_array_protein` | positives: SOWgp orthologs, BAD1 and its orthologs, the 20 Pro/Cys-rich class 2a candidates after review; negatives: the 10 Ser/Thr-rich class 2a candidates, the S288C and *C. albicans* repeat-truth proteins (arrays that are not Pro/Cys-rich) | `tuned_on_truth` (the thresholds came from these proteins) | `smoke` at best |

No independent truth exists for the architecture call (roadmap section B). The report says so, and does not state a sensitivity or specificity for new proteins. `calibrate truth` measures a module or a call that reads two or more modules. `sowgp_ortholog` reads one module and gets a module status (module `sowgp_unit`). `pro_cys_array_protein` reads three modules (`repeat02`, `repeat14`, `composition`) and gets a call status.

## 5. Decisions for the owner (stop points)

| # | Question | Default |
|---|---|---|
| O1 | Does `pro_cys_array_protein` appear as a sub-label of `surface_attachment_candidate`? Does `sowgp_ortholog` feed any composite? | sub-label yes; `sowgp_ortholog` feeds nothing |
| O2 | Are the composition thresholds (Pro+Cys over 15%, Cys at least 4%) accepted as the definition of Pro/Cys-rich? | accepted, as used 2026-09-30 |
| O3 | Is BAD1 a positive for `pro_cys_array_protein`? (BAD1 is a Pro/Cys-rich repeat protein of Blastomyces; the 2026-09-30 report lists a period-43 BAD1-like family.) | yes, with its orthologs |
| O4 | Review of the 20 Pro/Cys-rich class 2a candidates as positives. | assistant review (tier T4, not owner-reviewed) under a rubric fixed first |

## 6. Risks

1. Positives and thresholds come from the same proteins: any status is `smoke` with leakage, and says nothing about proteins unlike SOWgp.
2. The unit HMM has a sharp floor (about 40% identity). A diverged relative gives `not_called` and the report must not imply absence.
3. A composition rule does not see layout. Pro/Cys-rich proteins with a repeat but no relation to SOWgp are called (the 74 look-alikes are such cases, unreviewed).
4. SOWgp's gene models vary: in the pangenome 39 of 488 strains have no gene model but normal read depth, and 3 of 6 long-read loci split into two gene models (`REPORT_2026-10-01_sowgp_*`). A protein-level call misses split models.
5. Antigen and spherule columns exist for RS only.

## 7. Work plan

| Task | Content | Check |
|---|---|---|
| M5a | `sowgp_unit` module and derived `.hmm` with frozen GA 26; tests first (synthetic arrays; the derived file against `--cut_ga`); run on RS, *C. posadasii*, *B. dermatitidis*, S288C | RS calls CIMG_04613; S288C 0 calls |
| M5b | `composition` module; check script 14's `composition()` scope; reproduce the 83 and 74 counts on the Onygenales set; call `pro_cys_array_protein`; tests, golden files read line by line | the counts reproduce or the difference is explained |
| M5c | Truth tables and status files for both calls (leakage as in section 4); the owner-independent review of the 20 candidates (O4) | status files `smoke`; reviewed list stored |
| M5d | Candidate report on the scan proteomes (section 3.3) | table per proteome; SOWgp row in RS |
| M5e | Docs: `docs/CLASSES.md` entries, paper notes | test passes |

Each task is its own commit. An independent review of this spec and of the plan precedes M5a.
