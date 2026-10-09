# Review 1 of the M5 spec (SOWgp-like calls, revision 1)

*2026-10-09. Independent review of `docs/superpowers/specs/2026-10-09-sowgp-like-calls-design.md`. Reviewer: Claude Code (claude-opus-5-5). No spec or code file was edited. Every statement below was
checked with a command or by reading the named file. Scratch outputs (synthetic arrays, derived `.hmm` files, domain tables) are in the session scratchpad and are not committed.*

## 0. Checks run

- `unit_hmm_hits.tsv.gz` (2,055 domain rows) joined to `unit_genomes.tsv` by proteome label and genus.
- `unit_architecture_onygenales.tsv.gz` (644,399 rows) filtered with the section 3 rule.
- `class2a_candidates.tsv` and `class2a_candidates_general.tsv` recounted and filtered with the M5 rule.
- `hmmsearch` (HMMER 3.4) of `sowgp_unit.hmm` on `_workdir/sorting_hat/Cimm_RS.faa`, `Scer_S288C.faa` and `Bder_ER3.faa` (`-E 1000 --domE 1000`), and with derived GA lines `26 26`, `26 -1000` and `-1000 26` (`--cut_ga`).
- 1,000 synthetic arrays (modal unit, 30 to 45% identity, 1 to 12 copies, RS-composition flanks) searched with the same files.
- The proposed `pro_cys_array_protein` rule applied to the raw `repeat02.tsv` and `repeat14.tsv` of RS, S288C, *C. albicans* SC5314 and *B. dermatitidis* ER-3 in `_workdir/sorting_hat/*/raw/repeats/`.
- `engine.py`, `categories.yaml`, `modules/hydrophobin_relaxed.py`, `modules/pfam.py`, `modules/repeats.py`, `calibration/cli.py`, `calibration/measure.py`, `call_status.py`, `docs/paper/03-status-and-validation-rules.md`, the M2 spec and the hydrophobin specs read.
- PubMed record of PMID 12065484 (the source the panel cites for SOWgp).

## 1. Facts of section 2: what holds

| Claim | Result |
|---|---|
| 6,313 proteomes (7 long-read, 493 pangenome, 5,813 Fungi_5k) | Holds (`unit_genomes.tsv`, column `source`). |
| All 1,293 hits at 26 or more are *Coccidioides*; best other 19.7; none between 26 and 40 | Holds. *Coccidioides*: 1,748 domain hits, 1,293 at 26 or more, all of these at 40 or more. Others: 307 hits, max 19.7. |
| Null: 0 hits | Holds (`unit_calibration_null.tsv` has a header only). |
| Sensitivity table | Numbers hold (0.923/0.057 at 40%, 1.000/0.953 at 50%, 0.013 at 30%). But see finding 3: the floor is not 26 bits. |
| 83 proteins in 50 proteomes; 74 after SOWgp; 27 periods, none at 47 | Holds on the final file. 1,457 repeat calls. `>15` and `>=15` both give 83. Report section 3 was marked provisional; the recount on the final file agrees. |
| `pct_pro`, `pct_cys` scope | Whole protein. `composition(seq)` in script 14 (line 486) and script 02 (line 146) receives the full sequence (trailing `*` stripped). Script 02 and 14 give identical composition for every protein in the four proteomes checked. M5b's scope question is answered. |
| 41 / 58 candidates; 20 / 10 / 11 | Holds for the 41-row table. The 58-row table has 21 / 11 / 26. |
| Only CIMG_04613 of the Pro/Cys RS loci is spherule-induced | Holds (report section 4.5). |
| Spherule table: 9,757 genes, 564 up, SOWgp +9.09 (top) | Holds. |
| BAD1 control not run | Holds. `38_unit_bad1_control.py` writes `unit_bad1_like.tsv`; the file does not exist; `logs/` has no 38 or 39 log. |
| `sowgp_unit.hmm` | Exists. `NAME SOWgp_unit`, `LENG 47`, `NSEQ 34`. It has **no GA, TC or NC line**. |

## 2. Findings

### 1. BLOCKER: BAD1 cannot be a positive for `pro_cys_array_protein` under the proposed rule (O3)

Evidence. BAD1 (`bad1_A4D962.fa`, 1,146 aa) has Pro 3.8%, Cys 8.0%: Pro+Cys = 11.8%, below 15%. The 2026-09-30 report (section 4) gives the same values (8.5, 3.8, 8.0). In
`Bder_ER3.faa` BAD1 is `F00FD2C2_006066-T1` (354 shared 12-mers with A4D962). Both detectors call it a repeat (period 24, 27 copies, coverage 0.55). Composition: Pro 3.5%, Cys 8.0%. The proposed call is
`not_called`. The BAD1-like period-43 family in *Coccidioides* also fails (Pro+Cys 10.4 to 12.1%). The O3 default makes every BAD1 ortholog a false negative by construction.

Fix. Choose one. (a) Drop BAD1 from the positives and state that the call does not cover BAD1-type (Cys-rich, Pro-poor) arrays. (b) Change the definition (for example Cys at least 4% and (Pro+Cys over 15% or
Cys over 7%)) and record it as a new rule, not the 2026-09-30 rule; the 83/74 counts then no longer apply. Decide this before M5b.

### 2. BLOCKER: the positive set for `pro_cys_array_protein` contradicts the rule and the tier precedent (O4, section 4)

Evidence.
- Of the 20 "Pro/Cys-rich" class 2a candidates, **only 12 pass** Pro+Cys over 15% and Cys at least 4%. The other 8 have Cys 0 to 3.6% (for example `CIMG_07912` 1.0%, `CPOS1038_006679` 0.0%,
  `CPOS1038_008584` 3.6%). They were labelled by `03_repeat_surface_candidates.py`, which has no cysteine floor (2026-09-30 report, section 3). As positives they are false negatives by design.
- 4 of the 20 are SOWgp orthologs (CIMG_04613, CIB10637_003943, CIB10992_003451, QVM09276.1). They are counted twice with the "SOWgp orthologs" positives.
- The 20 collapse to 5 RS loci and 1 family with no RS gene model (report 2026-09-27, section 4.5): about 6 homology clusters, not 20.
- Hydrophobin precedent: T4 (review judgement, no paper) is "never added to the truth set, never used for sensitivity; reported only" (`2026-10-08-hydrophobin-validation-design.md` line 147;
  custom-HMM spec tier table: T4 "reported only"). O4 uses assistant-reviewed candidates as positives. This is not consistent.

Fix. Positives = proteins with an external label only (SOWgp orthologs; BAD1 only if finding 1 is resolved). Report the 20 (and the 74) as T4 review outcomes beside the measurement, never in it. Give
cluster counts, not protein counts. State that a status from these positives measures self-consistency only.

### 3. MAJOR: the sensitivity floor in section 2 and the "about 40%" reach in 3.1 are not at the 26-bit cutoff

Evidence. `36_unit_calibrate.py` counts a unit as found when its domain i-E-value is at most 1.0 with `domZ` 6e7 (`--floor-e`, default 1.0; header line 222). With the model's Forward stats
(tau -3.8654, lambda 0.71923) this is about 21.0 bits; the hit table agrees (score 22.0 has E 0.52). The table's `median_best_score` at 40% identity is 26.0. So at a 26-bit cutoff about half of the 40%
arrays have no unit called, not 8%.

Fix. In section 2 say "floor = domain E at most 1 (about 21 bits)". In 3.1 state the reach at 26 bits: about half of four-unit arrays at 40% identity are missed (from the median), and give 50% as the
identity where detection is near complete. Better: rerun `36_unit_calibrate.py` with a bit-score floor of 26 in M5a and quote that table.

### 4. MAJOR: "hits" are domain hits; the sequence score was never recorded

Evidence. In `37_unit_search.sh` the awk prints `$6` into the column `full_score`. In `--domtblout`, field 6 is `qlen`. All 2,055 rows have `full_score` 47. `36_unit_calibrate.py` reads
`fullscore` from `f[5]` (also `qlen`). Every number in the 2026-09-30 report is a **domain** score. The 1,293 hits are 1,293 domains in **487 proteins in 455 proteomes**. All 307 non-*Coccidioides*
rows are single-domain hits.

Fix. Section 2: write "1,293 domain hits (487 proteins, 455 proteomes) at domain score 26 or more". Name the score type in 3.1 and in the provenance file. Record the column bug in the provenance
file so nobody uses `full_score`.

### 5. MAJOR: the call expression in 3.2 cannot be written as stated, and the written form decides call-status eligibility

Evidence.
- `{call: X}` reads the `call` column of **module** X (`engine.py` line 192). `tandem_repeat_protein` is a call, not a module. A call is read with `{ref: tandem_repeat_protein}`.
- `call_eligible()` refuses any expression with a `ref` node ("contains ref", line 464). With `ref`, `calibrate truth --call-status` refuses the call. Section 4 says the call gets a call status.
- `{test: composition.pct_pro + pct_cys > 15}` is not a node. A `test` node has `module`, `field`, `op`, `value` and reads one field (lines 151 to 157, 203 to 215). The spec then says a derived field
  `pro_cys_pct` solves this. The two sentences disagree.

Fix. Write the exact YAML:
```yaml
- name: pro_cys_array_protein
  expr:
    and:
      - or: [{call: repeat02}, {call: repeat14}]
      - test: {module: composition, field: pro_cys_pct, op: ">", value: "$pro_cys_pct_min"}
      - test: {module: composition, field: pct_cys, op: ">=", value: "$cys_pct_min"}
```
This reads three modules and has no `ref`, so it is eligible. Add the two thresholds to `thresholds:`. Define `pro_cys_pct` from residue counts (unrounded), and say how it relates to the
rounded sum used for the 83 (both strict and non-strict give 83 on the Onygenales file).

### 6. MAJOR: "either detector" does not reproduce the 83 count

Evidence. The 83 come from script 14 only (`37_unit_search.sh` phase C runs `14_repeat_detect_general.py`). The proposed call uses `repeat02 OR repeat14`. In RS the rule gives 7 calls with either
detector and 6 with repeat14 alone; `XP_001239868.2` (242 aa, Pro 22.3, Cys 5.0) is called by repeat02 only.

Fix. In M5b reproduce the 83 with repeat14 alone, then report the extra calls that repeat02 adds. Or define the call on repeat14 only. State which.

### 7. MAJOR: the call name says Pro/Cys, but the rule passes Pro-poor, Cys-rich proteins

Evidence. In RS, `XP_001248770.2` (period 27) is called with Pro 3.5% and Cys 12.1%; `XP_001248220.2` (period 43) with Pro 9.2% and Cys 16.1%. In the 2026-09-30 table, periods 27, 41, 80 and 25 have
median Pro 3.5 to 7.3%. The sum lets cysteine alone carry the rule. Also, the SOWgp paper (PMID 12065484) describes the repeat as enriched in proline (20.4 mol%) and aspartate (19.7%), not cysteine;
RS SOWgp is Pro 15.7%, Asp 13.6%, Cys 6.5%.

Fix. Either add a proline floor (and recount), or rename the call to say what it tests (for example "Pro+Cys-rich array, Cys at least 4%") and state in the meaning that Pro can be low. Owner decision
(extend O2).

### 8. MAJOR: the truth plan does not fit `calibrate truth` (one taxon per entry, IDs, runs)

Evidence.
- `calibrate truth` takes exactly one species taxon and refuses truth proteins outside it (`cli.py` lines 461 to 495; paper 03 section 4). Section 4 mixes *C. immitis*, *C. posadasii*, 5,813 Fungi_5k
  proteomes, S288C and *C. albicans* in one row. Each species is its own entry and its own run.
- "Negatives: non-*Coccidioides* proteomes of the 6,313 set" cannot be a status entry: there is no run of the tool on those proteomes and they span thousands of species. It is a report-level
  specificity bound only.
- In S288C and *C. albicans* there are no positives. An entry there measures specificity only (`build_measure` allows this; notes say "specificity not measured" only when there are no negatives).
- The 10 Ser/Thr negatives are all outside RS (CiB10637 1, CiB10992 1, VFC140 2, Cpos1038 4, Cpos3700 1, Silveira 1). None of these proteomes is in `_workdir/sorting_hat/`.
- IDs differ. The analysis tables use FungiDB IDs (`CIMG_04613-t26_1-p1`). The tool's RS proteome uses RefSeq IDs (`XP_001245172.2`). Truth tables need a mapped ID column; the spec does not say how.
- `_workdir/sorting_hat/` has no *C. posadasii* proteome, but M5a runs on it.

Fix. Rewrite section 4 as one row per species and call: taxon, proteome file, run directory, positives (IDs in the run's namespace, with cluster), negatives, leakage. Move the Fungi_5k bound to the report.
Name the *C. posadasii* proteome file and add it to the run set in M5a.

### 9. MAJOR: GA design for `sowgp_unit` is under-specified; the M5a check can be answered now

Evidence.
- `sowgp_unit.hmm` has no GA line. A derived file must add one before `STATS`. I tested `GA 26.00 26.00;`, `GA 26.00 -1000.00;` and `GA -1000.00 26.00;`; HMMER 3.4 accepts all three and `--cut_ga` writes
  the option line and `# [ok]`, so `pfam.parse_domtblout` and `hydrophobin_relaxed.check_table` accept the tables.
- RS: one protein hits even at E 1000: `XP_001245172.2`, sequence score 324.8, seven domains, scores 89.1, 95.6, 94.4, 66.7, 22.5, -3.0, -6.5. With `GA 26 26` or `GA -1000 26` the table has 4 rows
  (n_units 4). The 22.5 unit is not counted. S288C and *B. dermatitidis* ER-3: zero rows at E 1000 and with every GA file.
- Sequence versus domain on synthetic arrays (40 per cell): at 36% identity and 8 copies, sequence score at least 26 calls 34, domain score at least 26 calls 7. At 33% and 12 copies: 20 versus 8. At 30%
  and 2 copies: 4 versus 0. A sequence-score decision (the hydrophobin pattern, domain GA -1000) calls arrays of weak units that no single domain supports. It has no calibration: the sequence score of
  the 6,313-proteome search was not saved (finding 4).
- `GA 26 26` gave the same proteins as `GA -1000 26` on the 1,000 synthetic arrays (264 each). In 67 of 800 reported synthetic proteins the sequence score is below the best domain score, so `GA 26 26`
  can in principle drop a protein with one domain at 26 or more; it did not happen here.
- `relaxed_rows()` keys on `seq_score` and keeps the best row per protein. `sowgp_unit` needs domain rows: count and best domain score. `model_info()` and `check_table()` can be reused; the row
  builder cannot.
- Two synthetic proteins print as 26.0 but were dropped by `--cut_ga`: the printed value is rounded. The wrapper must take the decision from the `--cut_ga` rows, not re-test the printed score.

Fix. Specify `GA -1000.00 26.00;` (domain score alone binds; matches the calibration) or `GA 26.00 26.00;` with the reason. Say that `hit` = at least one row in the `--cut_ga` table and `n_units` =
row count. Say that the sequence score is reported, not used. Replace the M5a "to be checked" sentence with these results and keep one synthetic multi-unit fixture as a test.

### 10. MAJOR: the cutoff 26 is a choice inside a gap, not a calibrated value

Evidence. Any cutoff between 19.7 and 40 gives the same decisions on the 6,313 proteomes. *Coccidioides* also has 455 domain hits between 18.3 and 26 (partial units such as the 22.5 unit above). Every
*Coccidioides* protein with any hit has at least one domain at 26 or more (487 of 487).

Fix. In 3.1 say: "26 is inside the empty interval 19.7 to 40; the data do not choose a value within it." Give the reason for 26 (the median best score at 40% identity) and put both in the provenance file.

### 11. MAJOR: section 2 Onygenales set is described too broadly and has duplicates

Evidence. The architecture scan ran on 78 proteomes (Fungi_5k Onygenales plus the 7 `cocci_longread` entries; `37_unit_search.sh` line 36), not on the 571 Onygenales. RS appears twice
(`Coccidioides_immitis_RS`, `CimmitisRS_FungiDB`) and Silveira twice (`Coccidioides_posadasii_Silveira`, `CposadasiiSilveira2022_FungiDB`). "Removing SOWgp" removes 9 proteins (period 47) from 9
proteomes of 7 strains. The 74 include duplicate proteins from the twice-annotated genomes.

Fix. Say "78 Onygenales proteomes (two genomes present twice)". Say that "74 after SOWgp" removed 9 SOWgp proteins. Give a deduplicated count if it is used later.

### 12. MINOR: wrong path for the panel

Evidence. `data/sorting_hat/panel.tsv` does not exist. The file is `docs/reports/data/sorting_hat/panel.tsv`. Its SOWgp rows (`XP_001245172.2`, `tandem_repeat_protein` and
`cocci_specificity_rank_top15`, tuning `tuned`) match the spec's statement.

Fix. Correct the path.

### 13. MINOR: name and gate disagree with the M2 spec

Evidence. `2026-10-09-v1-class-list-design.md` line 46 names the M5 calls `sowgp_ortholog` and `secreted_pro_cys_array`. M5 names the second `pro_cys_array_protein` with no gate. M2 rule 5.1
gives each call "the category"; M5 cross-lists `sowgp_ortholog` under two categories. `docs/CLASSES.md` does not exist yet (M5e edits it).

Fix. Align the name and gate in one of the two specs. State whether a call may have two categories, in M2.

### 14. MINOR: `other_*` calls and split gene models

Evidence. `other_not_surface` and `other_surface_no_mechanism` list mechanism calls (`categories.yaml` lines 97 and 103). A protein called only by `sowgp_ortholog` (for example a split gene model
fragment that carries units but no repeat call) would still land in `other_surface_no_mechanism`. Risk 4 names split models but the design does not use the ortholog call to recover them.

Fix. State whether `sowgp_ortholog` joins the mechanism lists (O1).

### 15. MINOR: "Neither is an adhesion claim" next to "Category: cell wall and adhesion"

Evidence. PubMed, PMID 12065484 (Hung et al. 2002, [doi:10.1128/IAI.70.7.3443-3456.2002](https://doi.org/10.1128/IAI.70.7.3443-3456.2002)): recombinant SOWgp binds laminin, fibronectin and
collagen IV; the deletion strain partly loses ECM binding and is less virulent. So SOWgp itself has published adhesin evidence. The look-alikes have none. The spec does not cite this. The M2 spec (O4)
says SOWgp's adhesion role is not part of M2.

Fix. Cite PMID 12065484 for SOWgp. Say the call names homology or composition, and that adhesion evidence exists for SOWgp only.

### 16. MINOR: task order and ship rule

Evidence. The hydrophobin plan split module work (L6a) from the `categories.yaml` edit (L6c, after a ship rule). M5b adds the call in the same task as the module, with no ship rule. M5c depends on O3
and O4 being decided, but the spec does not mark them as stop points before M5b.

Fix. Mark O2, O3 and O4 as decisions needed before M5b. Put the `categories.yaml` edit in its own task after the M5b counts are accepted.

### 17. MINOR: leakage labels

Evidence. `in_reference` for `sowgp_ortholog` is right: the 34 units come from all full-length pangenome SOWgp copies (`34_sowgp_unit_hmm.py` docstring). For `pro_cys_array_protein`,
`tuned_on_truth` is plausible but unverified: no file says how 15% and 4% were chosen. Both labels cap at `smoke` (paper 03 section 5), and SOWgp is one homology cluster, so neither call can reach
`estimated` anyway (20 clusters per class needed).

Fix. Keep the labels. Add one sentence: the source of the 15% and 4% thresholds is not recorded.

## 3. Items checked and correct

- `calibrate truth --module sowgp_unit --call sowgp_ortholog` is allowed: the call reads one module (`cli.py` lines 457 to 470).
- `calibrate truth --call-status` is allowed for `pro_cys_array_protein` if finding 5 is applied (three modules, no `ref`).
- A threshold change changes `config_sha256` and `call_hash` (thresholds are resolved into the hash), so an old call status becomes stale. The spec's sentence on module identity is right.
- The derived `.hmm` file can travel in `ModuleSpec` artefacts as `hydrophobin_relaxed.hmm` does.
- `basis_calls` accepts `pro_cys_array_protein` (an earlier ungated call), so O1's default can be implemented.

## 4. Verdict

**Not ready for a plan.** Two blockers: the BAD1 positive fails the proposed rule (1), and the positive set mixes rule failures, duplicates and T4 review labels (2). Fix those, rewrite the call YAML
(5), restate the truth per species (8), and correct the floor and score-type facts (3, 4). Findings 9 and 10 can close the M5a "to be checked" item now. After that, a second review of the revised spec
is needed before the plan.
