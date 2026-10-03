# Design spec: Basidiomycota truth and curated GPI rows for step 1 validation (issue #50)

*Drafted 2026-10-02. Revision 3 on 2026-10-02, after two independent reviews by different models
(revision 1: 20 findings; revision 2: 16 findings; both "needs rework"). DRAFT. Revision 3 needs a third
independent review before any plan or code (owner instruction, 2026-10-02). No code, data or job exists
for this spec. Owner decisions are in section 10.*

Inputs: the step 1 spec `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md`
(cited "step 1 spec", with decisions Q1 to Q10 and rulings R-A to R-C); the Phase C spec
`docs/superpowers/specs/2026-10-01-step1-phase-c-evaluation-design.md` (rulings C-4, C-8);
`docs/model-review/STATUS.md`; `docs/HANDOFF-2026-10-02.md`; `analysis/step1_compare/` (README,
`species.tsv`, `curated_gpi.tsv`, `COLUMNS.md`, `labels.py`, `truth_table.py`, `02_attach_sequences.py`,
`03_triage_pm.py`, `d8_triage.py`, `07_build_features.py`, `phasec/`); the run outputs in
`_workdir/step1_compare/` (2026-10-02 run).

## 1. Purpose and non-goals

**Purpose.** Step 1 predicts surface localisation (cell wall, outer face of the plasma membrane, or
secreted). Two parts of its validation are weak or missing:

1. **Basidiomycota.** The Phase C Basidiomycota clade block has 16 direct-evidence positives
   (*Cryptococcus* H99 7, *Ustilago* 9) and 60 negatives. It is a smoke test. Decision Q8 says: curate
   Basidiomycota truth before the first release.
2. **GPI-anchored proteins.** `analysis/step1_compare/curated_gpi.tsv` has a header and no rows. The
   P-gpi label is a list, not a scored stratum (ruling R-B).

This spec defines how to add curated truth for both, and how that truth enters the pipeline.

**Non-goals.**
- It does not change a label definition (P-ext, P-gpi, N-int, N-sec, PM-TM, ambiguous; step 1 spec 2.2).
- It does not change ruling R-B (P-gpi follows the D8 rule). Owner decision, 2026-10-02: path (b).
- It does not train or choose a model. Gate values stay unset until after measurement.
- It does not decide if a smoke-test Basidiomycota result allows a release (decision D-A, section 10).
- It does not cover step 2, step 3, or the #14 hard-negative panel (step 2 work).
- It does not change `surface.tsv` or the UniProt-keyword tier T-c.

## 2. Facts checked

"Re-checked" means I re-ran the check on the repo or the run outputs. "Reviewer" means an independent
reviewer reported it and I did not re-check it.

| Fact | Value | Status |
|---|---|---|
| Rule for "estimate" (ruling C-8, `phasec/findings.py`) | recall half-width at most 0.10 for **every** one of R2, M8, M35, M8-C, M35-C, H, **and** at least 20 direct-evidence positives; otherwise "smoke test". The rule applies per test set | re-checked |
| Basidiomycota roles in `species.tsv` | H99 `test_clade`; JEC21 and CRYD1 `alternate_file`; *U. maydis* MYCMD `undecided` | re-checked |
| `splits.py` treats `test_clade` and `undecided` alike (S3 split) | yes | re-checked |
| `08_build_eval_tables.py` drops `alternate_file` rows; it stops when a truth row's `role` differs from `species.tsv` | yes (lines 179, 205; 165-169) | re-checked |
| `labels.classify` builds the label from GO **terms** (a surface term and an internal term give `ambiguous`; a surface term alone gives P-ext; secretory terms without surface terms give N-sec) | yes (`labels.py`) | re-checked |
| `labels.is_pm_candidate`: P-ext with a **non-IEA** plasma-membrane term | yes | re-checked |
| `truth_table.TRUTH_COLUMNS` has 25 columns | yes | re-checked |
| Step 02 stops, and deletes `truth_sequences.tsv.gz`, when a P-ext or ambiguous gene has no sequence | yes (`02_attach_sequences.py:118-128, 145-153`) | re-checked |
| Step 07 stops when a sequence has no SignalP or PredGPI call | yes (`07_build_features.py:179-184`) | re-checked |
| `d8_triage.classify_pm` returns P-gpi for a literature row **before** it checks the TM feature | yes (`d8_triage.py:108-110`) | re-checked |
| D8 uses `curated_gpi.tsv` only as a set of `(source_id, gene_id)` keys (`03_triage_pm.py:199-201`); a row that matches no candidate is ignored without an error | yes | re-checked |
| `d8_gpi_outside_pext.tsv` comes from reviewed UniProt entries with ECO:0000269 only (`03_triage_pm.py:113-129`) | yes | re-checked |
| D8 classes per source (non-alternate): PM-TM Scer 3, Calb 22, Spom 4, Afum 1, Anid 2, Umay 1 (33 genes); `pm-unresolved` Scer 9, Calb 37, Spom 6, Afum 8, Anid 3, H99 4, Umay 2; P-gpi Anid 1, Calb 1 | | re-checked |
| Direct-evidence P-ext genes that are plasma-membrane candidates, non-alternate sources | Saccharomycotina 67 (Scer 9, Calb 58); Taphrinomycotina 9; Eurotiomycetes 8; Basidiomycota 3 (H99 2, Umay 1) | re-checked |
| Phase C test sets include S1:all (pooled folds), S2-Scer_SGD, S2-Calb_CGD, S2-Spom_PomBase, and S3 sets for the clades whose sources are `test_clade` or `undecided` (Eurotiomycetes, Basidiomycota). There is no Saccharomycotina S3 set | | re-checked (`splits.py`, `11_evaluate.py`) |
| Tier `T-a` is hard-coded (`truth_table.py:141`); `phasec/` has no tier dimension; `dedupe.merge_group` has no `tier` column | yes | re-checked |
| Literature rows enter Phase C through `eval_literature.tsv` and part `test_lit`, which is fixed to Eurotiomycetes and takes no negatives | | re-checked (files); clade and negatives: reviewer |
| H99 P-ext genes | 11 (PQP1, CDA2, CDA1, QSP1, CPL1, YOR1, LAC2, cnap1, LAC1, CDA3, PLB1); direct evidence 9; CDA1 and CDA3 are `homology_only=yes` (ISS only) | re-checked |
| H99 genes that D8 marks `pm-unresolved` | CDA1, CDA2, CDA3, PLB1. D8 excludes them from scoring. The scored H99 set has 7 positives | re-checked |
| *U. maydis* GO P-ext genes | 62; 10 with direct evidence; 52 `homology_only=yes` | re-checked |
| *U. maydis* direct P-ext genes | lep1, See1, him3, rsp1, CMU1, afu1, ROW1, UMAG_00792, PIT2, rep1. ROW1 is a PM-TM gene, so 9 are scored | reviewer |
| Recall half-widths, Basidiomycota (`metrics.json`, V-go, direct) | H99: 7 positives, worst candidate 0.43. *U. maydis*: 9 positives, worst 0.33. Clade block: 16 positives, H 0.107, R2 0.158, M35 0.219, M8 0.222. All "smoke test"; floor not met | re-checked |
| Recall of R2 and M8 on Basidiomycota today | 0.125 and 0.750 | STATUS.md |
| Cluster-bootstrap half-width divided by binomial half-width: 1.18 to 1.58 (S1, S2 yeast sets); 0.88 to 1.14 (S3 sets) | | reviewer |
| Two tests named `test_tc_excludes_test_proteins` | `tests/step1_compare/test_keyword_tier.py:29` (D10) and `test_phasec_splits.py:144` (Phase C) | re-checked |
| `test_test_truth_sources` (`test_phasec_splits.py:133`) checks only that T-c rows are not test truth | | re-checked (name, line); content: reviewer |
| Phase C wall times for C1 steps 09, 10, 11 (`phasec/logs/wall.*.txt`): 61 s; 62 to 91 s; up to 723 s. Steps 08 and 12 were not timed | | re-checked |
| H99 YOR1 (J9VQH1, an ABC transporter) is P-ext on extracellular-region IDA (PMID:36247839) and extracellular-vesicle HDA (PMID:34377375). Its membrane terms are IEA only | | reviewer; re-checked by the second reviewer in the GOA file |
| UniProt 2026_03 reviewed entries with ECO:0000269 GPI evidence | S288C 3, *C. albicans* 3, *S. pombe* 0, *A. fumigatus* 0, *A. nidulans* 1 | step 1 spec 2.2; `d8_counts.tsv` (reviewer) |
| FungiDB downloads | HTTP 401 on 2026-09-27 | issue #14 |
| Curation order of the handoff on `main` | "#50 first, then #14" (`docs/HANDOFF-2026-10-02.md:99`, PR #52) | re-checked |

Not known: how many Basidiomycota cell-wall or secreted proteins have experimental evidence in the
literature; how many curated negatives exist; whether FungiDB gives H99 GO annotations that differ from
the GOA file; the effort of the literature work.

## 3. Design

Two workstreams. They share four genes (section 3.3).

### 3.1 Workstream A: Basidiomycota truth

**Step A1: FungiDB H99 check.** Compare H99 GO cellular-component annotations in FungiDB with
`313589.C_neoformans_var_grubii_H99.goa`. Report the genes with wall or extracellular terms that the GOA
file lacks, with their evidence codes. If FungiDB is unreachable, record that and go to A2. A2 does not
depend on A1. The result is a count. It changes no label.

**Step A2: Literature curation of H99 and *U. maydis*.**
- *Positives.* Proteins with experimental evidence of wall, surface or secreted location. Start from the
  11 H99 GO P-ext genes and the 10 *U. maydis* direct-evidence GO P-ext genes. Add proteins from the
  literature. Every row has a PMID and an evidence statement.
- *Negatives.* A negative row names a compartment through the `go_term` column (section 5). The label
  then comes from the parent rule (`labels.classify`), applied to the GO terms of the GOA file **plus**
  the curated term:
  - N-int needs a non-IEA internal term, and no endomembrane, plasma-membrane, vacuole, membrane,
    cell-periphery, wall or extracellular term at any evidence level, IEA included.
  - N-sec needs a non-IEA endomembrane, plasma-membrane or vacuole term, and no wall or extracellular
    term at any evidence level.
  A row whose derived label differs from its `expected_label` goes to `curated_conflicts.tsv`. The
  curator does not choose the label by hand.
- *Weak evidence.* A row that rests only on extracellular-vesicle proteomics or only on
  high-throughput codes (HDA, HMP, HEP, HGI, HTP) gets `htp_only=yes`. It stays out of headline metrics,
  as in ruling R-A. The GO seeds get the same re-check. YOR1 is the first case.
- *Ambiguous.* Surface and internal evidence together: decision Q9 applies.

### 3.2 Workstream B: curated GPI rows (path b, owner decision 2026-10-02)

**What a curated row does today.** `classify_pm` returns P-gpi for a literature row before it looks at
the TM feature. A row therefore changes two groups of genes in the D8 set:
- `pm-unresolved` genes become P-gpi.
- PM-TM genes become P-gpi (33 genes today, section 2). A training negative then turns into a positive.

**Owner decision, 2026-10-02: a curated row does not override a UniProt TM feature by default.** The
plan changes `classify_pm` so that a literature row on a gene with a TM feature keeps the class PM-TM and
writes the gene to `curated_conflicts.tsv`. The owner reviews that file. After review, the owner sets
`override_tm=yes` on the row in `curated_gpi.tsv`, and only then does the row give P-gpi. The default
value is `no`. The behaviour for reviewed UniProt ECO:0000269 entries does not change. No literature
rows exist today, so the change has no effect on current results. The plan checks the D8 tests.

**Bounds.** The direct-evidence bound for P-gpi positives (section 2) is: S2-Scer_SGD 9, S2-Calb_CGD 58,
S1:all 67 (cross-validation, both training sources pooled), and 9 (Taphrinomycotina), 8
(Eurotiomycetes), 3 (Basidiomycota). Ruling C-8 applies per test set, so P-gpi stays a **list** in every
reported test set. A scored P-gpi stratum needs a new ruling from the owner (section 7).

**Purpose of workstream B.**
1. Fill the P-gpi list with curated evidence, by clade.
2. Record curated GPI evidence for the model card and for error analysis.

Workstream B is not a GPI performance measure. The model card says "GPI-anchored protein performance is
not validated".

**Step B1: Collect literature rows.** One row per protein with experimental evidence of GPI anchoring.
Do not add a row for a protein that already has a reviewed UniProt ECO:0000269 GPI entry, because D8
reads those (7 entries, section 2). The six proteins that the step 1 spec names behave as follows
(re-checked): YPS1 is a D8 candidate; CWP1, SAG1 and CCW12 are already P-ext positives; GAS1 and TIP1 are
ambiguous. Rows for CWP1, SAG1 and CCW12 add literature evidence for the list and the card only.

**Step B2: Rules.**
- Predictor output never counts as evidence (decision Q2).
- A protein enters P-gpi only when it is P-ext, has a non-IEA plasma-membrane term, and has a curated
  row or a UniProt ECO:0000269 entry (D8 rule, ruling R-B), subject to the TM rule above.
- `gene_id` is the native identifier of the source (for example SGD `S000004924` for GAS1), because D8
  matches on `(source_id, gene_id)`. The UniProt accession goes in its own column.
- A new check (and test) requires that every row in `curated_gpi.tsv` matches a truth gene. The check
  writes `curated_gpi_unmatched.tsv` for rows that match no gene, and lists rows whose gene is outside
  P-ext. The current code ignores both cases silently.

### 3.3 Coupling of A and B

CDA1, CDA2, CDA3 and PLB1 are H99 P-ext genes that D8 marks `pm-unresolved`. They are excluded from
scoring. All three H99 wall genes are in this set, so the scored H99 wall stratum is empty. A literature
check can move them to P-gpi, or leave them unresolved. Treat these four genes as a joint A and B target.
Decision 1 (start B1 in parallel with A2) stays, because the literature search does not depend on A.

### 3.4 Entry of curated truth into the pipeline

I have not read all of `01_extract_go_truth.py` and the Phase C code in detail. The plan must confirm
every item below against the code.

**New file.** `analysis/step1_compare/curated_basidiomycota.tsv` (columns in section 5). Tier `T-b`.

**New merge step.** A step `01b` runs after `01_extract_go_truth.py` and before `02_attach_sequences.py`.
It works on GO **terms**, not on labels:
1. Index the source FASTA (`sequences.index_fasta`, mapping `uniprot`). Write `curated_unmatched.tsv` for
   each curated row with no sequence. Unmatched rows do not enter truth. Without this check, step 02
   stops the whole chain on one stale accession.
2. For each remaining row, add the pair `(go_term, evidence_code)` to the gene's `GeneRecord.rows`, or
   create the record for a gene that has no GO row.
3. Re-run `truth_table.build_truth_rows`. The label, subset, `homology_only`, `pm_candidate` (a non-IEA
   plasma-membrane term) and the evidence columns then follow from the parent rule.
4. Write `curated_conflicts.tsv`: genes whose label differs from `expected_label`, genes whose label
   changed against the GO-only label, genes that became ambiguous, and genes in the TM override case
   (section 3.2). The owner reviews this file.

A curated protein outside the reference proteome has no sequence and no SignalP or PredGPI call. Step 07
stops without those calls. Such a protein cannot enter Phase C. It stays in `curated_unmatched.tsv`.

**Tier.** A gene gets `T-a`, `T-b` or `T-a+T-b`. `dedupe.merge_group` joins the tiers of the members of a
hash group (sorted, unique, comma separated), as it does for `d8_class`. `11_evaluate.py` treats a test
set row as T-b when its tier list contains `T-b`.

**Columns of a merged row.**

| Column | Source |
|---|---|
| `source_id`, `species`, `taxon_id`, `in_clade`, `role` | `species.tsv` |
| `gene_id`, `symbol`, `synonym1` | the GO row, or the curated row (`synonym1` empty) |
| `label`, `subset`, `stratum` | `build_truth_rows` |
| `tier` | as above |
| `label_no_homology`, `label_experimental`, `homology_only` | `build_truth_rows`, from the added pairs |
| `pm_candidate` | `labels.is_pm_candidate` |
| `evidence_codes`, `surface_evidence`, `internal_evidence`, `internal_evidence_htp_only`, `secretory_evidence` | `build_truth_rows`; the curated pair adds its `evidence_code`, so step 08 does not stop on empty evidence |
| `source_file`, `source_sha256`, `source_date`, `obo_sha256` | GAF values when the gene has a GO row; otherwise the curated file name, its SHA-256, the latest `review_date`, and the pinned OBO hash |

**No label override.** A curated row cannot overwrite a GO label. It adds terms, and the parent rule
decides. Whether a `direct` curated row may outweigh GO `transfer` evidence is an evidence-policy
question for the owner (decision D-C).

**Phase C changes.**
- `08_build_eval_tables.py` carries `tier` into `eval_table.tsv.gz`.
- `11_evaluate.py` reports each Basidiomycota test set by tier and for the T-b-only subset.
- `12_report.py` prints the tier tables.
T-b rows enter through the GO test path (`eval_table`), not through `test_lit`. The `test_lit` path is
fixed to Eurotiomycetes and has no negatives.

**Role changes in `species.tsv`.** `Umay_MYCMD` changes from `undecided` to `test_clade` (decision 2).
No split changes, because `splits.py` treats both roles alike. Step 01 must still rerun, because `role`
is a truth column that step 08 compares with `species.tsv`. Decision 3 (JEC21 IBA truth as a separate
"homology-transfer" stratum) needs extra work, because step 08 drops `alternate_file` rows. The 52
homology-only *U. maydis* P-ext rows are a second homology-transfer set. The plan chooses how to produce
the stratum and whether it includes those rows. I do not propose a design here.

## 4. Evidence rules (both workstreams)

| Rule | Detail |
|---|---|
| Source of each row | A PMID. The PMID must resolve in PubMed. A row with none is not accepted |
| Quote | `evidence_note` holds the sentence from the paper that states the evidence, and the retrieval date |
| Evidence level | `direct`: the paper measures the location or the anchor in this protein. `transfer`: the paper measures it in an ortholog. Only `direct` rows count for headline metrics (step 1 spec 3.3) |
| Predictor-selected candidates | `selected_by_predictor` is `yes`, `no` or `unknown`. Some secreted-protein papers may choose candidates with SignalP (assumption, not checked). R0 would then recall them by construction. Reports show recall with and without `yes` rows |
| Second check | A second reviewer checks each row. A model may do the first pass. A model can repeat a wrong PMID. So the reviewer opens each PMID. The reviewer finds the quoted sentence. The reviewer records the date. The owner spot-checks a random sample (decision 5) |
| Accession | `gene_id` is the native identifier of the source in `species.tsv`. `uniprot_accession` is a separate column. A protein with no sequence in the reference proteome goes to `curated_unmatched.tsv` and is not in truth |
| No predictor labels | SignalP, PredGPI and TM predictions never set a label |
| Moonlighting | A protein with surface and internal evidence is `ambiguous` (Q9) |

## 5. Data format

`curated_basidiomycota.tsv` (new), one row per (gene, GO term, evidence):

| Column | Meaning |
|---|---|
| source_id | Source in `species.tsv` (`Cneo_H99_GOA` or `Umay_MYCMD`) |
| species | Species name (for readability; `source_id` is the key) |
| gene_id | Native identifier of the source (a UniProt accession for these two sources) |
| uniprot_accession | UniProt accession |
| symbol | Gene symbol, if any |
| go_term | GO cellular-component accession that the paper supports (for example GO:0005618 for the cell wall) |
| evidence_code | GO evidence code of the experiment (for example IDA, EXP) |
| expected_label | The label the curator expects after the parent rule: `P-ext`, `N-int`, `N-sec` or `ambiguous`. A check, not a label |
| evidence_level | `direct` or `transfer` |
| htp_only | `yes` when the evidence is vesicle proteomics or high-throughput only |
| selected_by_predictor | `yes`, `no` or `unknown` |
| pmid | PMID (several separated by `;`) |
| evidence_note | The quoted sentence, and the retrieval date |
| reviewer | Initials or model name of the second check |
| review_date | `YYYY-MM-DD` |

`curated_gpi.tsv` keeps `source_id`, `gene_id`, `symbol`, `pmid`, `note`. It gains `species`,
`uniprot_accession`, `evidence_level` (decision 4), `evidence_note`, `reviewer`, `review_date` and
`override_tm` (`no` by default; the owner sets `yes` after review, section 3.2). `gene_id` is the native
identifier of the source.

Both files are small plain text and stay uncompressed in git.

## 6. Leakage controls

| Risk | Control | Test or check |
|---|---|---|
| Curator sees model output. `proteome_calls.tsv.gz` holds calls and scores of every candidate for all H99 and *U. maydis* proteins | The curator (person or model) gets no read access to `_workdir/step1_compare/phasec/`. The curator logs each search query. The candidate list is hashed and committed before any join to scores. A protein found by the search is not dropped after scoring | process rule, recorded in the PR; `candidate_list.sha256` committed before the join |
| Residual exposure: aggregate Basidiomycota results are in `docs/model-review/STATUS.md` (lines 51, 93, 108) and in the report. A curator with repo access sees them | State this in the PR. The curator sees aggregates only. Gene-level calls stay in `_workdir` | process rule |
| Literature chose candidates by a predictor | `selected_by_predictor` column; report with and without | report format |
| Curated test proteins in T-c | D10 removes every curated accession and exact-sequence hash from T-c. Phase C rules `a` (accession and hash), `b_cluster_mate` and `c_test_taxon` apply | extend **both** `test_tc_excludes_test_proteins` (`test_keyword_tier.py:29`, `test_phasec_splits.py:144`) |
| Basidiomycota curated row in training | Basidiomycota stays a test clade. `splits.build` raises `SplitError` when a protein belongs to a training source and a test source | add a test on the new rows; keep the role in `species.tsv` |
| Homology across clades. Phase C allows cluster-mates of test proteins in training by design (ruling C-4). 3 of the 16 current positives have 0.3 or more identity to training (reviewer) | Report T-b rows by maximum-identity stratum, as for GO rows | report format |
| T-b and T-a are not independent, because A2 starts from the T-a genes | Report the T-b-only subset | report format |
| Positive and negative with the same sequence | Existing dedupe drops hashes in both classes | existing dedupe test |
| GPI rows bias toward yeast | Report P-gpi by clade. Do not pool | report format |

## 7. Acceptance and statistics

**The rule (ruling C-8).** A test set is an "estimate" when the recall half-width is at most 0.10 for
every one of R2, M8, M35, M8-C, M35-C and H, and the set has at least 20 direct positives. Otherwise it
is a "smoke test". The rule applies per test set.

**Size, binomial approximation.** n = 1.96² x p(1-p) / 0.10², rounded up. The candidate with recall
nearest 0.5 needs the most positives.

| Recall p | Positives needed |
|---|---|
| 0.125 (R2 on Basidiomycota today) | 43 |
| 0.750 (M8 on Basidiomycota today) | 73 |
| 0.5 (worst case) | 97 |

At today's recalls of R2 and M8, the binding value is **73**. The recalls of M35, M8-C, M35-C and H move
as new truth arrives, so the binding value can reach 97. The measured ratio of the cluster-bootstrap
half-width to the binomial half-width is 0.88 to 1.58 (reviewer). That ratio can move the number in
either direction. The number is an approximation, not a measured Basidiomycota value.

**Today (re-checked):** clade block 16 positives, half-widths H 0.107, R2 0.158, M35 0.219, M8 0.222; the
20-positive floor is not met. H99 alone has 7 positives (worst half-width 0.43). *U. maydis* alone has 9
(worst 0.33).

**Acceptance for workstream A.**
- Curated Basidiomycota truth is merged by step `01b`, with a source for every row, and
  `curated_unmatched.tsv` and `curated_conflicts.tsv` are written.
- Phase C reports each Basidiomycota test set by tier and for the T-b-only subset, with cluster-bootstrap
  intervals.
- The report states "estimate" or "smoke test" by the rule above.
- If the set is a smoke test, the model card says Basidiomycota is not validated. Decision D-A says if
  that outcome allows a release.

**Acceptance for workstream B.**
- `curated_gpi.tsv` has rows, each with a resolvable PMID, and the new match check passes.
- P-gpi is a list in every reported test set. A scored P-gpi stratum is outside this spec, because ruling
  C-8 is defined per test set. If the owner wants one, the owner makes a new ruling first.
- The report states the P-gpi counts by clade.

**Both.**
- All real-data metrics will move. Curated GPI rows can change training labels in *S. cerevisiae* and
  *C. albicans*, and so can change every candidate. The 2026-10-02 run is the baseline. The report shows
  before and after values for each candidate. The golden metrics file (`PHASEC_WRITE_GOLDEN=1`) covers the
  test fixture only. Only intended fixture paths change.
- No accuracy statement enters a README or card unless it comes from `metrics.json`.

## 8. Compute and effort

**Compute (measured).** C1 steps 09, 10 and 11 took 61 s, 62 to 91 s and up to 723 s. Steps 08 and 12
were not timed. Step 09 was timed once. The 1 to 1.5 hour sizing rule is for fan-out jobs and does not
apply to one job.

**Rerun chain.** A change in `truth_set.tsv.gz` or in D8 changes `truth_set_sha256`. The hash checks then
force the later steps: 01b, 02, 03, 04, 05, 07, 08, C1 and 12 (reviewer). Step 01 reruns when its inputs
or `species.tsv` change. The role change of `Umay_MYCMD` does change `species.tsv`. Step 03 needs the
UniProt network. The J1 and J2 jobs resume by hash. The proteome FASTA files of H99 and *U. maydis* are
already embedded, so a curated protein in those files has an embedding, unless the FASTA release changes
(reviewer).

**Effort.** I have no estimate of the effort. The size of the Basidiomycota and GPI literature is not
measured.

## 9. Risks

| Risk | How it is detected |
|---|---|
| FungiDB stays unreachable | A1 reports it. A2 does not depend on A1 |
| Too few experimental Basidiomycota positives | the counts in section 7; the smoke-test label stays |
| The 20-positive floor is not met per source | the clade block is the label unit today (decision D-B) |
| *Cryptococcus* capsule and secreted proteins blur "wall" and "extracellular" | `subset` column; recall per subset |
| Vesicle-proteome-only evidence labels a protein positive | `htp_only` column; excluded from headline |
| Curated negatives are rare, so FPR stays imprecise | report the interval; no gate from it |
| One curator introduces bias | second review (section 4); `reviewer` column |
| Curated rows duplicate GO rows | tier `T-a+T-b`; report the T-b-only subset |
| Training labels change and metrics move | baseline run; before and after table |
| A curated GPI row matches no gene and is ignored | the new check and `curated_gpi_unmatched.tsv` |
| A stale accession halts step 02 | step 01b writes `curated_unmatched.tsv` before the merge |
| A curated GPI row turns a TM negative into a positive | `override_tm=no` by default; `curated_conflicts.tsv` for owner review |

## 10. Decisions for the owner

**Decided on 2026-10-02** (the owner accepted each recommendation):
1. Start with A1, then B1 in parallel with A2. This matches the curation order on `main` (#50 before
   #14, `docs/HANDOFF-2026-10-02.md:99`).
2. Keep *U. maydis* as a second Basidiomycota set, reported beside H99.
3. Keep JEC21 IBA truth as a separate "homology-transfer" stratum with no weight in headline metrics.
4. Add `evidence_level` to `curated_gpi.tsv`.
5. An independent model run does the first review pass. The owner spot-checks a random sample.
6. Stop curation when the clade block reaches 73 direct positives (section 7), or at a time box that
   the owner sets later.
7. Workstream B follows path (b): P-gpi stays a list (ruling R-B unchanged).
8. A curated GPI row does not override a UniProt TM feature by default. The owner reviews each case and
   sets `override_tm=yes` to allow it.
9. Revision 3 of this spec gets a third independent review before any plan or code.

**Open (from the reviews).**
- **D-A. Does a smoke-test Basidiomycota result allow a release?** Q8 rejected "not validated flag
  only", and the step 1 spec section 11 says no model ships before Basidiomycota truth is in the
  evaluation. Decision 6 allows a time box. Recommendation: keep Q8 as written. If the time box ends
  below 73 positives, report to the owner for a ruling. The card flag alone does not satisfy Q8.
- **D-B. Pooled clade block versus "not pooled".** Phase C already reports a pooled
  `S3-Basidiomycota:clade` block, and the estimate label uses it. Recommendation: keep the clade block as
  the label unit, show each source beside it, and mark the clade block as pooled in the report.
- **D-C. Evidence policy: may a `direct` curated row outweigh GO `transfer` evidence?** Today a curated
  row only adds terms, and the parent rule decides (section 3.4). Recommendation: keep it that way. A
  conflict goes to `curated_conflicts.tsv`.
- **D-D. The JEC21 and *U. maydis* "homology-transfer" strata (decision 3).** They need a pipeline change
  because step 08 drops `alternate_file` rows. Recommendation: put the design choice in the plan, and size
  it after the plan names the change.

## 11. Deliverables

| ID | Deliverable |
|---|---|
| E1 | A1 report: FungiDB versus GOA for H99 (counts, codes), or a record that FungiDB was unreachable |
| E2 | `curated_basidiomycota.tsv` with PMIDs, quoted evidence, and second review; `candidate_list.sha256` and the query log |
| E3 | `curated_gpi.tsv` rows with PMIDs; `curated_gpi_unmatched.tsv` and the new match check with its test |
| E4 | Code: step `01b` merge (with `curated_unmatched.tsv` and `curated_conflicts.tsv`); `classify_pm` TM rule and `override_tm`; `truth_table.py` tier handling; `dedupe.merge_group` tier join; `08_build_eval_tables.py`, `11_evaluate.py`, `12_report.py` tier dimension; `species.tsv` roles; D10 removal; the tests in section 6 |
| E5 | Phase C rerun, baseline comparison table, report sections per tier, per source and per clade |
| E6 | `STATUS.md` update with new Basidiomycota and P-gpi numbers and intervals |

Done means: E1 to E6 exist; the owner has the report; the model card states the clade scope and the
validation state of Basidiomycota and GPI proteins.
