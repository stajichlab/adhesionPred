# Design spec: Basidiomycota truth and curated GPI rows for step 1 validation (issue #50)

*Drafted 2026-10-02. Revision 2 on 2026-10-02, after an independent review by a different model
(20 findings, verdict "needs rework"). DRAFT. Revision 2 needs a second independent review before any
plan or code. No code, data or job exists for this spec. Owner decisions are in section 10.*

Inputs: the step 1 spec `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md`
(cited "step 1 spec", with decisions Q1 to Q10 and rulings R-A to R-C); the Phase C spec
`docs/superpowers/specs/2026-10-01-step1-phase-c-evaluation-design.md` (rulings C-4, C-8);
`docs/model-review/STATUS.md`; `docs/HANDOFF-2026-10-02.md`; `analysis/step1_compare/` (README,
`species.tsv`, `curated_gpi.tsv`, `COLUMNS.md`, `03_triage_pm.py`, `truth_table.py`,
`phasec/08_build_eval_tables.py`, `phasec/findings.py`); the run outputs in
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
- It does not decide if a smoke-test Basidiomycota result satisfies Q8 (decision D-A, section 10).
- It does not cover step 2, step 3, or the #14 hard-negative panel (step 2 work).
- It does not change `surface.tsv` or the UniProt-keyword tier T-c.

## 2. Facts checked

"Re-checked" means I re-ran the check on the repo or run outputs after the review. "Reviewer" means
the independent reviewer reported it and I did not re-check it.

| Fact | Value | Status |
|---|---|---|
| Rule for "estimate" (ruling C-8, `phasec/findings.py`) | recall half-width at most 0.10 for **every** one of R2, M8, M35, M8-C, M35-C, H, **and** at least 20 direct-evidence positives; otherwise "smoke test" | re-checked |
| Basidiomycota roles in `species.tsv` | H99 `test_clade`; JEC21 and CRYD1 `alternate_file`; *U. maydis* MYCMD `undecided` | re-checked |
| `08_build_eval_tables.py` drops `alternate_file` rows | yes (lines 179, 205) | re-checked |
| H99 P-ext genes | 11 (PQP1, CDA2, CDA1, QSP1, CPL1, YOR1, LAC2, cnap1, LAC1, CDA3, PLB1); direct evidence 9; CDA1 and CDA3 are `homology_only=yes` (ISS only) | re-checked |
| H99 genes that D8 marks `pm-unresolved` | CDA1, CDA2, CDA3, PLB1. D8 excludes them from scoring. The scored H99 set has 7 positives | re-checked |
| Per-source and clade recall half-widths, Basidiomycota (`metrics.json`, V-go) | H99: 7 positives, worst candidate 0.43. *U. maydis*: 9 positives, worst 0.33. Clade block: 16 positives, H 0.107, R2 0.158, M35 0.219, M8 0.222. All "smoke test"; floor not met | re-checked |
| Direct-evidence P-ext genes that are plasma-membrane candidates (the only genes where a `curated_gpi.tsv` row can act), non-alternate sources | Saccharomycotina 67 (Scer 9, Calb 58); Taphrinomycotina 9; Eurotiomycetes 8; Basidiomycota 3 (H99 2, Umay 1) | re-checked |
| D8 uses `curated_gpi.tsv` only as a set of `(source_id, gene_id)` keys (`03_triage_pm.py:199-200`); a row that matches no candidate is ignored without an error | yes | re-checked |
| `d8_gpi_outside_pext.tsv` comes from reviewed UniProt entries with ECO:0000269 only (`03_triage_pm.py:114-129`); it does not read `curated_gpi.tsv` | yes | re-checked |
| Tier `T-a` is hard-coded (`truth_table.py:141`); `phasec/` has no tier dimension | yes | re-checked |
| Literature rows enter Phase C through `eval_literature.tsv` and part `test_lit`, which is fixed to Eurotiomycetes and takes no negatives | yes | re-checked (file names); clade and negatives: reviewer |
| Two tests named `test_tc_excludes_test_proteins` | `tests/step1_compare/test_keyword_tier.py:29` (D10) and `test_phasec_splits.py:144` (Phase C) | re-checked |
| `test_test_truth_sources` (`test_phasec_splits.py:133`) checks that T-c rows are not test truth only | | reviewer |
| Cluster-bootstrap half-width divided by binomial half-width: 1.18 to 1.58 (S1, S2 yeast sets); 0.88 to 1.14 (S3 sets) | | reviewer |
| Phase C wall times (`phasec/logs/wall.*.txt`): step 09 61 s; step 10 62 to 91 s; step 11 up to 723 s | | re-checked |
| Rerun chain after a change in `truth_set.tsv.gz` or D8: 01, 02, 03, 04, 05, 07, 08, C1, 12; step 03 needs the UniProt network; J1 and J2 resume by hash | | reviewer |
| *U. maydis* direct P-ext genes | 10 (lep1, See1, him3, rsp1, CMU1, afu1, ROW1, UMAG_00792, PIT2, rep1). ROW1 is a PM-TM gene, so 9 are scored | reviewer |
| H99 YOR1 (J9VQH1, an ABC transporter) is P-ext on extracellular-region IDA (PMID:36247839) and extracellular-vesicle HDA (PMID:34377375). Its membrane terms are IEA only | | reviewer |
| UniProt 2026_03 reviewed entries with ECO:0000269 GPI evidence | S288C 3, *C. albicans* 3, *S. pombe* 0, *A. fumigatus* 0, *A. nidulans* 1 | step 1 spec 2.2; counts in `d8_counts.tsv` (reviewer) |
| FungiDB downloads | HTTP 401 on 2026-09-27 | issue #14 |

Not known: how many Basidiomycota cell-wall or secreted proteins have experimental evidence in the
literature; how many curated negatives exist; whether FungiDB gives H99 GO annotations that differ
from the GOA file; the effort of the literature work.

## 3. Design

Two workstreams. They share four genes (section 3.3).

### 3.1 Workstream A: Basidiomycota truth

**Step A1: FungiDB H99 check.** Compare H99 GO cellular-component annotations in FungiDB with
`313589.C_neoformans_var_grubii_H99.goa`. Report the genes with wall or extracellular terms that the
GOA file lacks, with their evidence codes. If FungiDB is unreachable, record that and go to A2. A2 does
not depend on A1. The result is a count. It changes no label.

**Step A2: Literature curation of H99 and *U. maydis*.**
- *Positives.* Proteins with experimental evidence of wall, surface or secreted location. Start from the
  11 H99 and 10 *U. maydis* GO P-ext genes. Add proteins from the literature. Every row has a PMID and
  an evidence statement.
- *Negatives.* Every curated negative must pass the parent rule against the GOA file **at any evidence
  level**, IEA included:
  - N-int: no endomembrane, plasma-membrane, vacuole, membrane, cell-periphery, wall or extracellular
    term (step 1 spec 2.2).
  - N-sec: no wall or extracellular term (step 1 spec 2.2).
  A negative that fails the rule is `ambiguous` or is not used. The spec does not loosen the label
  definitions.
- *Weak evidence.* A row that rests only on extracellular-vesicle proteomics or only on
  high-throughput codes (HDA, HMP, HEP, HGI, HTP) gets `htp_only=yes`. It stays out of headline metrics,
  as in ruling R-A. The seeds from GO get the same re-check. YOR1 is the first case.
- *Ambiguous.* Surface and internal evidence together: decision Q9 applies.

### 3.2 Workstream B: curated GPI rows (path b, owner decision 2026-10-02)

A `curated_gpi.tsv` row changes a label only for a P-ext gene that is a plasma-membrane candidate in D8.
It moves that gene from `pm-unresolved` to P-gpi. The direct-evidence bound for P-gpi positives is in
section 2: Saccharomycotina 67, Taphrinomycotina 9, Eurotiomycetes 8, Basidiomycota 3. Clades are not
pooled. So P-gpi stays a **list** in every clade, except possibly Saccharomycotina (section 7).

**The purpose of workstream B is therefore:**
1. Fill the P-gpi list with curated evidence, by clade.
2. Record curated GPI evidence for the model card and for error analysis.
3. It is not a GPI performance measure. The model card says "GPI-anchored protein performance is not
   validated" unless section 7 shows otherwise.

**Step B1: Collect literature rows.** One row per protein with experimental evidence of GPI anchoring.
Do not add a row for a protein that already has a reviewed UniProt ECO:0000269 GPI entry. D8 reads those
(7 entries, section 2). The six proteins that the step 1 spec names behave as follows (re-checked):
YPS1 is a D8 candidate; CWP1, SAG1 and CCW12 are already P-ext positives; GAS1 and TIP1 are ambiguous.
So rows for CWP1, SAG1 and CCW12 add literature evidence for the list and the card only.

**Step B2: Rules.**
- Predictor output never counts as evidence (decision Q2).
- A protein enters P-gpi only when it is P-ext, has a non-IEA plasma-membrane term, and has a curated
  row or a UniProt ECO:0000269 entry (D8 rule, ruling R-B).
- `gene_id` is the native identifier of the source (for example SGD `S000004924` for GAS1), because D8
  matches on `(source_id, gene_id)`. The UniProt accession goes in its own column.
- A new check (and test) requires that every row in `curated_gpi.tsv` matches a truth gene. The check
  writes `curated_gpi_unmatched.tsv` for rows that match no gene, and lists rows whose gene is outside
  P-ext. The current code ignores both cases silently.

### 3.3 Coupling of A and B

CDA1, CDA2, CDA3 and PLB1 are H99 P-ext genes that D8 marks `pm-unresolved`. They are excluded from
scoring. All three H99 wall genes are in this set, so the scored H99 wall stratum is empty. A literature
check can move them to P-gpi or PM-TM. Treat these four genes as a joint A and B target. Decision 1
(start B1 in parallel with A2) stays, because the literature search does not depend on A.

### 3.4 Entry of curated truth into the pipeline

I have not read all of `01_extract_go_truth.py`, `truth_table.py` and the Phase C code in detail. The
plan must confirm every item below against the code.

**New file.** `analysis/step1_compare/curated_basidiomycota.tsv` (columns in section 5). Tier `T-b`.

**New merge step.** A step `01b` runs after `01_extract_go_truth.py` and before `02_attach_sequences.py`.
It merges curated rows into `truth_set.tsv.gz`, so that D8 (`03`) and all later steps see them. A curated
row passes through D8 like a GO row.

**Column mapping for a T-b row.**

| Column | Value |
|---|---|
| `tier` | `T-b`, or `T-a+T-b` when GO and curation agree on the gene |
| `homology_only` | `no` when `evidence_level=direct`; `yes` when `transfer` |
| `evidence_codes`, `surface_evidence`, `internal_evidence` | from the row's `evidence_code` (a GO evidence code, for example IDA or EXP), so that step 08 does not stop on empty evidence (`08_build_eval_tables.py:119-130`) |
| `pm_candidate` | from the GO file when the gene has a plasma-membrane term, else `no` |
| `subset`, `stratum`, `label` | from the row |
| `internal_evidence_htp_only` | `yes` when `htp_only=yes` |

**One row per gene, with a precedence rule (proposal; decision D-C).** A gene has one row per
`(source_id, gene_id)`.
- GO and curation agree: keep one row, tier `T-a+T-b`.
- They disagree: if the curated row is `direct` and the GO evidence is `transfer` or high-throughput
  only, the curated row wins. In every other case the gene becomes `ambiguous`.
- Every disagreement goes to `curated_conflicts.tsv`. The owner reviews that file.

**Phase C changes.** Phase C has no tier dimension today. Add one:
- `08_build_eval_tables.py` carries `tier` into `eval_table.tsv.gz`.
- `11_evaluate.py` reports each Basidiomycota test set by tier (T-a, T-b, T-a+T-b) and for the
  T-b-only subset.
- `12_report.py` prints the tier tables.
T-b rows enter through the GO test path (`eval_table`), not through `test_lit`. The `test_lit` path is
fixed to Eurotiomycetes and has no negatives.

**Role changes in `species.tsv`.** `Umay_MYCMD` changes from `undecided` to `test_clade`
(decision 2). Decision 3 (JEC21 IBA truth as a separate "homology-transfer" stratum) needs extra work,
because step 08 drops `alternate_file` rows. The plan must choose how to produce that stratum. I do
not propose a design here.

## 4. Evidence rules (both workstreams)

| Rule | Detail |
|---|---|
| Source of each row | A PMID. The PMID must resolve in PubMed. A row with none is not accepted |
| Quote | `evidence_note` holds the sentence from the paper that states the evidence, and the retrieval date |
| Evidence level | `direct`: the paper measures the location or the anchor in this protein. `transfer`: the paper measures it in an ortholog. Only `direct` rows count for headline metrics (step 1 spec 3.3) |
| Predictor-selected candidates | `selected_by_predictor` is `yes`, `no` or `unknown`. Many secreted-protein papers choose candidates with SignalP, so R0 recalls them by construction. Reports show recall with and without `yes` rows |
| Second check | A second reviewer checks each row. A model may do the first pass. A model can repeat a wrong PMID, so the reviewer opens the PMID, finds the quoted sentence, and records the date. The owner spot-checks a random sample (decision 5) |
| Accession | `gene_id` is the native identifier of the source in `species.tsv`. `uniprot_accession` is a separate column. A protein with no match in the reference proteome goes to `unmatched` and is not in truth |
| No predictor labels | SignalP, PredGPI and TM predictions never set a label |
| Moonlighting | A protein with surface and internal evidence is `ambiguous` (Q9) |

## 5. Data format

`curated_basidiomycota.tsv` (new), one row per gene:

| Column | Meaning |
|---|---|
| source_id | Source in `species.tsv` (`Cneo_H99_GOA` or `Umay_MYCMD`) |
| species | Species name (for readability; `source_id` is the key) |
| gene_id | Native identifier of the source (a UniProt accession for these two sources) |
| uniprot_accession | UniProt accession |
| symbol | Gene symbol, if any |
| label | `P-ext`, `N-int`, `N-sec` or `ambiguous` |
| subset | `wall` or `extracellular-only` for P-ext; empty otherwise |
| evidence_level | `direct` or `transfer` |
| evidence_code | GO evidence code of the experiment (for example IDA, EXP) |
| htp_only | `yes` when the evidence is vesicle proteomics or high-throughput only |
| selected_by_predictor | `yes`, `no` or `unknown` |
| pmid | PMID (several separated by `;`) |
| evidence_note | The quoted sentence, and the retrieval date |
| reviewer | Initials or model name of the second check |
| review_date | `YYYY-MM-DD` |

`curated_gpi.tsv` keeps `source_id`, `gene_id`, `symbol`, `pmid`, `note`. It gains `species`,
`uniprot_accession`, `evidence_level` (decision 4), `evidence_note`, `reviewer` and `review_date`.
`gene_id` is the native identifier of the source.

Both files are small plain text and stay uncompressed in git.

## 6. Leakage controls

| Risk | Control | Test or check |
|---|---|---|
| Curator sees model output. `proteome_calls.tsv.gz` already holds calls and scores of every candidate for all H99 and *U. maydis* proteins, and `report.md` lists the 16/60 results | The curator (person or model) gets no read access to `_workdir/step1_compare/phasec/`. The curator logs each search query. The candidate list is hashed and committed before any join to scores. A protein found by the search is not dropped after scoring | process rule, recorded in the PR; `candidate_list.sha256` committed before the join |
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
is a "smoke test". The candidate with recall nearest 0.5 needs the most positives.

**Size, binomial approximation.** n = 1.96² x p(1-p) / 0.10².

| Recall p | Positives needed |
|---|---|
| 0.125 (R2 on Basidiomycota today) | 42 |
| 0.75 (M8 on Basidiomycota today) | 72 |
| 0.5 (worst case) | 97 |

The cluster-bootstrap interval is not always wider than the binomial one. The reviewer measured ratios of
0.88 to 1.58 (section 2). So the number of positives needed is between about 42 and 97 per reported test
set, and may be lower or higher. It is not a measured Basidiomycota value.

**Today (re-checked):** clade block 16 positives, half-widths H 0.107, R2 0.158, M35 0.219, M8 0.222;
the 20-positive floor is not met. H99 alone has 7 positives (worst half-width 0.43). *U. maydis* alone
has 9 (worst 0.33).

**Acceptance for workstream A.**
- Curated Basidiomycota truth is merged by step `01b` with a source for every row.
- Phase C reports each Basidiomycota test set by tier and for the T-b-only subset, with cluster-bootstrap
  intervals.
- The report states "estimate" or "smoke test" by the rule above.
- If the set is a smoke test, the model card says Basidiomycota is not validated. Whether that outcome
  allows a release is decision D-A.

**Acceptance for workstream B.**
- `curated_gpi.tsv` has rows, each with a resolvable PMID, and the new check passes.
- P-gpi is a list in every clade. It becomes a scored stratum only if a clade has enough direct P-gpi
  positives to meet the rule above. Saccharomycotina has a bound of 67 and needs 42 to 97, so it is the
  only possible case. The report states the count.

**Both.**
- All real-data metrics will move. Curated GPI rows change training labels in *S. cerevisiae* and
  *C. albicans* (`pm-unresolved` becomes P-gpi), and so change every candidate. The 2026-10-02 run is the
  baseline. The report shows before and after values for each candidate. The golden metrics file
  (`PHASEC_WRITE_GOLDEN=1`) covers the test fixture only. Only intended fixture paths change.
- No accuracy statement enters a README or card unless it comes from `metrics.json`.

## 8. Compute and effort

**Compute (measured).** A Phase C run takes: step 09 61 s, step 10 62 to 91 s, step 11 up to 723 s. This
is under 15 minutes per run. The 1 to 1.5 hour sizing rule is for fan-out jobs and does not apply to one
job.

**Rerun chain.** A change in `truth_set.tsv.gz` or in D8 changes `truth_set_sha256`. The hash checks then
force 01, 01b, 02, 03, 04, 05, 07, 08, C1 and 12 (reviewer). Step 03 needs the UniProt network. J1 and J2
resume by hash. New H99 and *U. maydis* sequences are already embedded as whole proteomes, unless the
FASTA release changes (reviewer).

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

## 10. Decisions for the owner

**Decided on 2026-10-02** (the owner accepted each recommendation):
1. Start with A1, then B1 in parallel with A2. This matches the revised order in the handoff (#50
   before #14).
2. Keep *U. maydis* as a second Basidiomycota set, reported beside H99.
3. Keep JEC21 IBA truth as a separate "homology-transfer" stratum with no weight in headline metrics.
4. Add `evidence_level` to `curated_gpi.tsv`.
5. An independent model run does the first review pass. The owner spot-checks a random sample.
6. Stop curation when the section 7 size is reached, or at a time box that the owner sets later.
7. Workstream B follows path (b): P-gpi stays a list (ruling R-B unchanged).

**Open (new, from the review).**
- **D-A. Does a smoke-test Basidiomycota result allow a release?** Q8 rejected "not validated flag
  only", and the step 1 spec section 11 says no model ships before Basidiomycota truth is in the
  evaluation. Decision 6 allows a time box. Recommendation: keep Q8 as written. If the time box ends
  below the section 7 size, report to the owner for a ruling. The card flag alone does not satisfy Q8.
- **D-B. Pooled clade block versus "not pooled".** Phase C already reports a pooled
  `S3-Basidiomycota:clade` block, and the estimate label uses it. Recommendation: keep the clade block as
  the label unit, show each source beside it, and mark the clade block as pooled in the report.
- **D-C. Precedence rule for a curated row and a GO row (section 3.4).** Recommendation: the proposal in
  section 3.4, with the owner reviewing `curated_conflicts.tsv`.
- **D-D. The JEC21 "homology-transfer" stratum (decision 3).** It needs a pipeline change because step 08
  drops `alternate_file` rows. Recommendation: put the design choice in the plan, and size it after the
  plan names the change.

## 11. Deliverables

| ID | Deliverable |
|---|---|
| E1 | A1 report: FungiDB versus GOA for H99 (counts, codes), or a record that FungiDB was unreachable |
| E2 | `curated_basidiomycota.tsv` with PMIDs, quoted evidence, and second review; `candidate_list.sha256` and the query log |
| E3 | `curated_gpi.tsv` rows with PMIDs; `curated_gpi_unmatched.tsv` and the new match check with its test |
| E4 | Code: step `01b` merge; `truth_table.py` tier handling; `08_build_eval_tables.py`, `11_evaluate.py`, `12_report.py` tier dimension; `species.tsv` roles; D10 removal; the tests in section 6 |
| E5 | Phase C rerun, baseline comparison table, report sections per tier, per source and per clade |
| E6 | `STATUS.md` update with new Basidiomycota and P-gpi numbers and intervals |

Done means: E1 to E6 exist; the owner has the report; the model card states the clade scope and the
validation state of Basidiomycota and GPI proteins.
