# Design spec: Basidiomycota truth and curated GPI rows for step 1 validation (issue #50)

*Drafted 2026-10-02. Revision 4 on 2026-10-02, after three independent reviews by different models
(revision 1: 20 findings, "needs rework"; revision 2: 16 findings, "needs rework"; revision 3: 19
findings, "ready after fixes"). DRAFT. Revision 4 applies the third review. The owner decides if a fourth
review is needed before the plan. No code, data or job exists for this spec. Owner decisions are in
section 10.*

Inputs: the step 1 spec `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md`
(cited "step 1 spec", with decisions Q1 to Q10 and rulings R-A to R-C); the Phase C spec
`docs/superpowers/specs/2026-10-01-step1-phase-c-evaluation-design.md` (rulings C-4, C-8);
`docs/model-review/STATUS.md`; `docs/HANDOFF-2026-10-02.md`; `analysis/step1_compare/` (README,
`species.tsv`, `curated_gpi.tsv`, `COLUMNS.md`, `gaf.py`, `labels.py`, `truth_table.py`,
`02_attach_sequences.py`, `03_triage_pm.py`, `d8_triage.py`, `07_build_features.py`, `phasec/`); the run
outputs in `_workdir/step1_compare/` (2026-10-02 run).

## 1. Purpose and non-goals

**Purpose.** Step 1 predicts surface localisation (cell wall, outer face of the plasma membrane, or
secreted). Two parts of its validation are weak or missing:

1. **Basidiomycota.** The Phase C Basidiomycota clade block has 16 direct-evidence positives
   (*Cryptococcus* H99 7, *Ustilago* 9) and 60 negatives (`metrics.json`, `S3-Basidiomycota:clade`). It
   is a smoke test. Decision Q8 says: curate Basidiomycota truth before the first release.
2. **GPI-anchored proteins.** `analysis/step1_compare/curated_gpi.tsv` has a header and no rows. The
   P-gpi label is a list, not a scored stratum (ruling R-B).

This spec defines how to add curated truth for both, and how that truth enters the pipeline.

**Non-goals.**
- It does not change a label definition (P-ext, P-gpi, N-int, N-sec, PM-TM, ambiguous; step 1 spec 2.2).
- It amends one ruling only: the "until" clause of R-B (decision 12, section 10).
- It does not train or choose a model. Gate values stay unset until after measurement.
- It does not cover step 2, step 3, or the #14 hard-negative panel (step 2 work).
- It does not change `surface.tsv` or the UniProt-keyword tier T-c.

## 2. Facts checked

"Re-checked" means I re-ran the check on the repo or the run outputs. "Reviewer" means an independent
reviewer reported it and I did not re-check it.

| Fact | Value | Status |
|---|---|---|
| Rule for "estimate" (ruling C-8, `phasec/findings.py`) | recall half-width at most 0.10 for **every** one of R2, M8, M35, M8-C, M35-C, H, **and** at least 20 direct-evidence positives; otherwise "smoke test". The rule applies per test set | re-checked |
| Ruling R-A | keeps the D1 rule. It only reports `internal_evidence_htp_only` as a sub-stratum of ambiguous genes. It excludes nothing from the headline (step 1 spec 76-80, 445-446) | re-checked (reviewer) |
| Evidence-code sets | `HOMOLOGY_CODES` in `gaf.py:17`; `HIGH_THROUGHPUT_CODES` (HDA, HMP, HEP, HGI, HTP) in `labels.py:16`; the high-throughput codes are also in `EXPERIMENTAL_CODES` (reviewer), so an HDA-only positive counts as direct today | re-checked (sets); membership: reviewer |
| Basidiomycota roles in `species.tsv` | H99 `test_clade`; JEC21 and CRYD1 `alternate_file`; *U. maydis* MYCMD `undecided` | re-checked |
| `splits.py` treats `test_clade` and `undecided` alike (S3 split) | yes | re-checked |
| `08_build_eval_tables.py` drops `alternate_file` rows; it stops when a truth row's `role` differs from `species.tsv` | yes (lines 179, 205; 165-169) | re-checked |
| `labels.classify` builds the label from GO **terms** | yes (`labels.py`) | re-checked |
| `labels.is_pm_candidate`: P-ext with a **non-IEA** plasma-membrane term | yes | re-checked |
| `truth_table.build_truth_rows(filtered, ontology, info)` takes a parsed GAF object and builds the gene records itself with `collect_genes(filtered)` (from `filtered.cc_rows`); `count_rows` also reads `filtered.cc_rows` | yes (`truth_table.py:94-124, 177-178`) | re-checked |
| `truth_set.tsv.gz` keeps joined code lists per gene, not `(term, evidence)` pairs | | reviewer |
| `truth_table.TRUTH_COLUMNS` has 25 columns | yes | re-checked |
| Step 02 stops, and deletes `truth_sequences.tsv.gz`, when a P-ext or ambiguous gene has no sequence | yes (`02_attach_sequences.py:122-128, 149-153`) | re-checked |
| Step 07 stops when a sequence has no SignalP or PredGPI call | yes (`07_build_features.py:179-184`) | re-checked |
| `d8_triage.classify_pm` returns P-gpi for a literature row **before** it checks the TM feature | yes (`d8_triage.py:108-110`) | re-checked |
| D8 uses `curated_gpi.tsv` only as a set of `(source_id, gene_id)` keys (`03_triage_pm.py:199-201`); a row that matches no candidate is ignored without an error | yes | re-checked |
| `d8_gpi_outside_pext.tsv` comes from reviewed UniProt entries with ECO:0000269 only (`03_triage_pm.py:113-129`) | yes | re-checked |
| `truth_set_sha256` in `d8_run.json` is the SHA-256 of `truth_set.tsv.gz` (`03_triage_pm.py:216`). Step 08 also checks that `features_run.json` holds the SHA-256 of `truth_set_triaged.tsv.gz` (`08_build_eval_tables.py:105-110`) | yes | re-checked |
| `tests/step1_compare/test_d8.py:334-341` uses MSB2 (a TM gene) as a literature row and asserts P-gpi; lines 142 and 338 write the five-column `curated_gpi.tsv` header | yes | re-checked |
| D8 classes per source (non-alternate): PM-TM Scer 3, Calb 22, Spom 4, Afum 1, Anid 2, Umay 1 (33 genes); `pm-unresolved` Scer 9, Calb 37, Spom 6, Afum 8, Anid 3, H99 4, Umay 2; P-gpi Anid 1, Calb 1 | | re-checked |
| TM evidence of the PM-TM genes: all 39 PM-TM rows (33 non-alternate) have TM evidence ECO:0000255 or ECO:0000256 only (predicted); 14 of the 33 have one TM segment | | re-checked (39 rows, codes); 14 of 33: reviewer |
| Direct-evidence P-ext genes that are plasma-membrane candidates, by D8 class: with PM-TM genes Scer 9, Calb 58, Spom 9, Eurotiomycetes 8, Basidiomycota 3. **Without** PM-TM genes (P-gpi or `pm-unresolved` only): Scer 7, Calb 37, Spom 5, Afum 4, Anid 3, H99 2, Umay 0 | | re-checked |
| Phase C test sets: `S1:all` (pooled folds), `S1:Scer_SGD`, `S1:Calb_CGD`, `S2-Scer_SGD:Scer_SGD`, `S2-Calb_CGD:Calb_CGD`, `S2-Spom_PomBase:...`, and S3 sets for the clades whose sources are `test_clade` or `undecided` (Eurotiomycetes, Basidiomycota). There is no Saccharomycotina S3 set | | re-checked (`splits.py`, `11_evaluate.py`, `metrics.json` names) |
| Tier `T-a` is hard-coded (`truth_table.py:141`); `phasec/` has no tier dimension; `dedupe.merge_group` has no `tier` column; `11_evaluate.py` has a fixed copy list of columns (lines 186-196) | yes | re-checked |
| `dedupe._join` joins the values of a hash group with "," | yes | reviewer |
| Literature rows enter Phase C through `eval_literature.tsv` and part `test_lit`, which is fixed to Eurotiomycetes and takes no negatives | | re-checked (files); clade and negatives: reviewer |
| H99 P-ext genes | 11 (PQP1, CDA2, CDA1, QSP1, CPL1, YOR1, LAC2, cnap1, LAC1, CDA3, PLB1); direct evidence 9; CDA1 and CDA3 are `homology_only=yes` (ISS only) | re-checked |
| H99 genes that D8 marks `pm-unresolved` | CDA1, CDA2, CDA3, PLB1. D8 excludes them from scoring. The scored H99 set has 7 positives | re-checked |
| *U. maydis* GO P-ext genes | 62; 10 with direct evidence; 52 `homology_only=yes` | re-checked |
| *U. maydis* direct P-ext genes | lep1, See1, him3, rsp1, CMU1, afu1, ROW1, UMAG_00792, PIT2, rep1. ROW1 is a PM-TM gene, so 9 are scored | re-checked (second reviewer) |
| CWP1, SAG1, CCW12 | P-ext in Scer, `pm_candidate=no` | re-checked |
| Recall half-widths, Basidiomycota (`metrics.json`, V-go, direct) | H99: 7 positives, worst candidate 0.43. *U. maydis*: 9 positives, worst 0.33. Clade block: 16 positives, H 0.107, R2 0.158, M35 0.219, M8 0.222. All "smoke test"; floor not met | re-checked |
| Recall of R2 and M8 on Basidiomycota today | 0.125 and 0.750 | STATUS.md |
| Cluster-bootstrap half-width divided by binomial half-width: 0.88 to 1.14 for S3 sets; 1.18 to 1.58 for S1 and S2 yeast sets. Clade block: H 0.90, R2 0.98, M8 1.05, M35 1.14 | | reviewer |
| Two tests named `test_tc_excludes_test_proteins` | `tests/step1_compare/test_keyword_tier.py:29` (D10) and `test_phasec_splits.py:144` (Phase C) | re-checked |
| `test_test_truth_sources` (`test_phasec_splits.py:133`) checks only that T-c rows are not test truth | | re-checked (name, line); content: reviewer |
| Phase C wall times for C1 steps 09, 10, 11 (`phasec/logs/wall.*.txt`): 61 s; 62 to 91 s; up to 723 s. Steps 08 and 12 were not timed | | re-checked |
| H99 YOR1 (J9VQH1, an ABC transporter): P-ext on extracellular-region IDA (PMID:36247839) **and** extracellular-vesicle HDA (PMID:34377375). Its membrane terms are IEA only. It is not high-throughput-only | | re-checked (GOA file, second and third reviewer) |
| UniProt 2026_03 reviewed entries with ECO:0000269 GPI evidence | S288C 3, *C. albicans* 3, *S. pombe* 0, *A. fumigatus* 0, *A. nidulans* 1 | step 1 spec 2.2; `d8_counts.tsv` (reviewer) |
| FungiDB downloads | HTTP 401 on 2026-09-27 | issue #14 |
| Curation order of the handoff | "#50 first, then #14" on `origin/main` at commit `9eaf7ef` (`docs/HANDOFF-2026-10-02.md:99`) | re-checked |

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

**Step A2: Literature curation of H99 and *U. maydis*.** A curator (a person or a model) finds the rows
and writes them. A second reviewer checks them (section 4).
- *Positives.* Proteins with experimental evidence of wall, surface or secreted location. Start from the
  11 H99 GO P-ext genes and the 10 *U. maydis* direct-evidence GO P-ext genes. Add proteins from the
  literature. Every row has a PMID and an evidence statement.
- *Negatives.* A negative row names a compartment through the `go_term` column (section 5). The parent
  rule (`labels.classify`) then gives the label from the GO terms of the GOA file **plus** the curated
  term:
  - N-int needs a non-IEA internal term, and no endomembrane, plasma-membrane, vacuole, membrane,
    cell-periphery, wall or extracellular term at any evidence level, IEA included.
  - N-sec needs a non-IEA endomembrane, plasma-membrane or vacuole term, and no wall or extracellular
    term at any evidence level.
  A row whose derived label differs from its `expected_label` goes to `curated_conflicts.tsv`. The
  curator does not choose the label by hand.
- *High-throughput evidence (decision 11).* The truth set gets a derived column `surface_evidence_htp_only`.
  It is `yes` when every non-IEA surface evidence code of the gene, from GO and from curation together,
  is in `HIGH_THROUGHPUT_CODES`. This follows the same logic as `internal_evidence_htp_only` (ruling R-A).
  Genes with `yes` stay in the headline. Reports show them as a separate sub-stratum. YOR1 has IDA, so
  it is not high-throughput-only.
- *Ambiguous.* Surface and internal evidence together: decision Q9 applies.

### 3.2 Workstream B: curated GPI rows

**What a curated row does today.** `classify_pm` returns P-gpi for a literature row before it looks at
the TM feature. A row therefore changes two groups of genes in the D8 set:
- `pm-unresolved` genes become P-gpi.
- PM-TM genes become P-gpi (33 genes today, section 2). A training negative then turns into a positive.

**Owner decision 8 (2026-10-02): a curated row does not override a UniProt TM feature by default.** The
plan changes the code as follows:
- `literature_ids` in `03_triage_pm.py` becomes a dictionary from `(source_id, gene_id)` to the row's
  `override_tm` value.
- `classify_pm` takes that value. A literature row on a gene with a TM feature keeps the class PM-TM
  unless `override_tm=yes`.
- Step 03 writes each such gene to `d8_curated_conflicts.tsv` (a file that only step 03 writes). The
  owner reviews it, then sets `override_tm=yes` on the row if the owner agrees.
- A `curated_gpi.tsv` without the new columns stops step 03 with an error.
- The behaviour for reviewed UniProt ECO:0000269 entries does not change.

No literature rows exist today, so the change has no effect on current results. It breaks one existing
test: `test_d8.py:334-341` expects P-gpi for MSB2. The plan changes that test to expect PM-TM, and adds a
case with `override_tm=yes`. The plan also changes the headers at `test_d8.py:142` and `:338`.

**Asymmetry for the owner.** All PM-TM genes have predicted TM evidence only (section 2), and 14 of 33
have a single TM segment, which a C-terminal GPI signal can mimic (reviewer). Decision 8 lets a predicted
TM feature block curated experimental evidence. A reviewed UniProt ECO:0000269 GPI entry still wins over
the TM feature. The two kinds of evidence are treated differently. The review file exists to catch the
cases where the TM prediction is wrong.

**Bounds, per test set (section 2).** Direct-evidence P-ext genes that can become P-gpi:

| Test set | Default (no override) | With override (PM-TM genes included) |
|---|---|---|
| `S2-Scer_SGD:Scer_SGD` | 7 | 9 |
| `S2-Calb_CGD:Calb_CGD` | 37 | 58 |
| `S1:all` (both training sources) | 44 | 67 |
| Taphrinomycotina (*S. pombe*) | 5 | 9 |
| Eurotiomycetes | 7 | 8 |
| Basidiomycota | 2 | 3 |

**Decision 12 (amends ruling R-B).** P-gpi stays a **list** in a test set until that test set has at
least 20 direct-evidence P-gpi positives. Then P-gpi is scored in that set, and ruling C-8 gives it an
"estimate" or "smoke test" label. Only `S2-Calb_CGD` and `S1:all` can reach 20 under the default bounds,
and only if curation finds enough rows. The 20 is the C-8 floor, so no new number is introduced.

**Purpose of workstream B.**
1. Fill the P-gpi list with curated evidence, by clade.
2. Record curated GPI evidence for the model card and for error analysis.
3. Score P-gpi in a test set only when the 20-positive condition holds.

Workstream B does not measure GPI performance in general. The model card says "GPI-anchored protein
performance is not validated" unless a test set meets the condition.

**Step B1: Collect literature rows.** One row per protein with experimental evidence of GPI anchoring.
Do not add a row for a protein that already has a reviewed UniProt ECO:0000269 GPI entry, because D8
reads those (7 entries, section 2). The six proteins that the step 1 spec names behave as follows
(re-checked): YPS1 is a D8 candidate. CWP1, SAG1 and CCW12 are P-ext with `pm_candidate=no`, so D8 never
reads them and they cannot enter the P-gpi list. Rows for these three reach the model card only. GAS1 and
TIP1 are ambiguous.

**Step B2: Rules.**
- Predictor output never counts as evidence for a curated row (decision Q2).
- A protein enters P-gpi only when it is P-ext, has a non-IEA plasma-membrane term, and has a curated
  row or a UniProt ECO:0000269 entry (D8 rule), subject to the TM rule above.
- `gene_id` is the native identifier of the source (for example SGD `S000004924` for GAS1), because D8
  matches on `(source_id, gene_id)`. The UniProt accession goes in its own column.
- A new check (and test) requires that every row in `curated_gpi.tsv` matches a truth gene. The check
  writes `curated_gpi_unmatched.tsv` for rows that match no gene. It also lists rows whose gene is outside
  P-ext, and rows whose gene is P-ext but not a plasma-membrane candidate. The current code ignores all
  three cases silently.

### 3.3 Coupling of A and B

CDA1, CDA2, CDA3 and PLB1 are H99 P-ext genes that D8 marks `pm-unresolved`. They are excluded from
scoring. All three H99 wall genes are in this set, so the scored H99 wall stratum is empty. A literature
check can move them to P-gpi, or leave them unresolved. CDA1 and CDA3 stay homology-only even if they
become P-gpi. Treat these four genes as a joint A and B target. Decision 1 (start B1 in parallel with A2)
stays, because the literature search does not depend on A.

### 3.4 Entry of curated truth into the pipeline

I have not read `01_extract_go_truth.py` in full. The plan must confirm every item below against the
code.

**The merge happens inside step 01.** A separate step after 01 cannot work, because `truth_set.tsv.gz`
does not keep `(term, evidence)` pairs, and `build_truth_rows` builds the gene records from the parsed
GAF. Step 01 gets a `--curated` input and changes as follows:
1. Refactor `build_truth_rows` to take the gene records (`build_truth_rows(genes, ontology, info)`). The
   existing call passes `collect_genes(filtered)`.
2. Merge the curated `(go_term, evidence_code)` pairs into the `collect_genes` output. Create the record
   for a gene that has no GO row.
3. Compute `count_rows` from the merged records, so that `counts.tsv` matches `truth_set.tsv.gz`.
4. Record the SHA-256 of the curated file in `extract_log.json`.
5. Write the result to `truth_set.tsv.gz` (steps 02, 03 and 04 read that fixed name).

**Checks on curated input (step 01 stops on each).**
- `evidence_code` is not in the GO evidence-code list.
- `go_term` is unknown, obsolete, an `alt_id`, or not a cellular-component term.
- `evidence_level=direct` with a code in `HOMOLOGY_CODES`, or `evidence_level=transfer` with a code
  outside `HOMOLOGY_CODES`. The column `evidence_level` is a checked field. The truth column
  `homology_only` still comes from the code.
- `gene_id` does not match a source sequence (below).
Step 01 normalises `gene_id` with `sequences.normalize_id`, so that an isoform suffix does not create a
second record that step 02 rejects.

**Sequence check.** Step 01 indexes the source FASTA (`sequences.index_fasta`, the mapping and file from
`species.tsv`). It writes `curated_unmatched.tsv` for each curated row with no sequence. Unmatched rows do
not enter truth. Without this check, step 02 stops the whole chain on one stale accession. A curated
protein outside the reference proteome has no SignalP or PredGPI call, and step 07 stops without those
calls. So such a protein cannot enter Phase C.

**Gene-level columns from curated rows.** One gene can have several curated rows. The merged truth set
carries:

| Column | Rule |
|---|---|
| `tier` | a comma list: `T-a`, `T-b` or `T-a,T-b` |
| `surface_evidence_htp_only` | derived from the codes (section 3.1) |
| `selected_by_predictor` | `yes` when every curated row that supports the gene's label has `yes`; `no` when at least one has `no`; otherwise `unknown`. Empty for genes without curated rows |
| `curated_pmid` | the sorted PMIDs of all curated rows of the gene, joined with ";" |

These columns go through steps 01, 03 and 08 into `eval_table.tsv.gz`, and into the copy list of
`11_evaluate.py`. `dedupe.merge_group` joins them for a hash group like `d8_class`. The plan checks that
the join works for each.

**Columns of a merged row.**

| Column | Source |
|---|---|
| `source_id`, `species`, `taxon_id`, `in_clade`, `role` | `species.tsv` |
| `gene_id`, `symbol`, `synonym1` | the GO row, or the curated row (`synonym1` empty) |
| `label`, `subset`, `stratum` | `build_truth_rows` |
| `label_no_homology`, `label_experimental`, `homology_only` | `build_truth_rows`, from the merged pairs |
| `pm_candidate` | `labels.is_pm_candidate` |
| `evidence_codes`, `surface_evidence`, `internal_evidence`, `internal_evidence_htp_only`, `secretory_evidence` | `build_truth_rows`; the curated pair adds its `evidence_code`, so step 08 does not stop on empty evidence |
| `source_file`, `source_sha256`, `source_date`, `obo_sha256` | GAF values when the gene has a GO row; otherwise the curated file name, its SHA-256, the latest `review_date`, and the pinned OBO hash. `SourceInfo` holds one provenance per source today, so the plan adds per-gene provenance |
| new columns | `tier`, `surface_evidence_htp_only`, `selected_by_predictor`, `curated_overexpressed` (same rule as `selected_by_predictor`: `yes` when every supporting row is `yes`), `curated_pmid` |

**No label override.** A curated row cannot overwrite a GO label. It adds terms, and the parent rule
decides. Step 01 writes `curated_conflicts.tsv` for these label cases: a derived label that differs from
`expected_label`; a label that changed against the GO-only label; a gene that became ambiguous. The TM
override case is a separate file written by step 03 (section 3.2). The owner reviews both files. Whether
a `direct` curated row may outweigh GO `transfer` evidence is decision D-C.

**Tiers in Phase C.**
- `08_build_eval_tables.py` carries `tier` into `eval_table.tsv.gz`.
- `11_evaluate.py` reports each Basidiomycota test set by three tier strata: genes whose tier list is
  `T-a` only, `T-b` only, and `T-a,T-b`. It also reports the `surface_evidence_htp_only` sub-stratum.
- The C-8 label (estimate or smoke test) applies to the whole test set. Tier strata show metrics with
  intervals and no label.
- `12_report.py` prints the tier tables.
T-b rows enter through the GO test path (`eval_table`), not through `test_lit`. The `test_lit` path is
fixed to Eurotiomycetes and has no negatives.

**Role changes in `species.tsv`.** `Umay_MYCMD` changes from `undecided` to `test_clade` (decision 2). No
split changes, because `splits.py` treats both roles alike. Step 01 must still rerun, because `role` is a
truth column that step 08 compares with `species.tsv`. Decision 3 (JEC21 IBA truth as a separate
"homology-transfer" stratum) needs extra work, because step 08 drops `alternate_file` rows. The 52
homology-only *U. maydis* P-ext rows are a second homology-transfer set. The plan chooses how to produce
the stratum and whether it includes those rows. I do not propose a design here.

## 4. Evidence rules (both workstreams)

| Rule | Detail |
|---|---|
| Source of each row | A PMID. The PMID must resolve in PubMed. A row with none is not accepted |
| Quote | `evidence_note` holds the sentence from the paper that states the evidence, and the retrieval date |
| Evidence level | `direct`: the paper measures the location or the anchor in this protein. `transfer`: the paper measures it in an ortholog. Step 01 checks that the level fits the evidence code (section 3.4). Headline metrics count a gene as direct through `homology_only=no` (step 1 spec 3.3) |
| Predictor-selected candidates | `selected_by_predictor` is `yes`, `no` or `unknown`. Some secreted-protein papers may choose candidates with SignalP (assumption, not checked). R0 would then recall them by construction. The pilots confirmed this (many *U. maydis* effector papers chose genes by SignalP). Decision 13: genes whose every supporting row is `yes` form a report-only sub-stratum and stay out of the headline |
| Overexpressed fusions | `overexpressed` is `yes` (strong constitutive or heterologous promoter such as otef or ACT1), `no` (native promoter and locus stated) or `unknown`. Such a row shows that the protein can reach that place, not where the native protein sits at normal levels. Decision 14 |
| Curator and reviewer | The curator finds and writes the rows. A second reviewer checks each row. An independent model run does the first review pass, and the owner spot-checks a random sample (decision 5) |
| Reviewer actions | The reviewer opens each PMID. The reviewer finds the quoted sentence. The reviewer records the date. A model can repeat a wrong PMID, so a model alone is not enough |
| Accession | `gene_id` is the native identifier of the source in `species.tsv`. `uniprot_accession` is a separate column. A protein with no sequence in the reference proteome goes to `curated_unmatched.tsv` and is not in truth |
| Predictors as evidence | A SignalP, PredGPI or TM prediction never creates a curated row or a curated label. D8 already uses UniProt TM features to set the class PM-TM (section 3.2) |
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
| evidence_code | GO evidence code of the experiment (for example IDA, EXP). Vesicle proteomics uses HDA |
| expected_label | The label the curator expects after the parent rule: `P-ext`, `N-int`, `N-sec` or `ambiguous`. A check, not a label |
| evidence_level | `direct` or `transfer` (checked against `evidence_code`) |
| selected_by_predictor | `yes`, `no` or `unknown` |
| overexpressed | `yes`, `no` or `unknown` (section 4). Added by decision 14; the file has 15 columns |
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
| Residual exposure: aggregate Basidiomycota results are in `docs/model-review/STATUS.md` (lines 43, 51, 93, 107, 108, 115) and in the report. A curator with repo access sees them | State this in the PR. The curator sees aggregates only. Gene-level calls stay in `_workdir` | process rule |
| Literature chose candidates by a predictor | `selected_by_predictor` column; sub-stratum, outside the headline (decision 13) | report format |
| Curated test proteins in T-c | D10 removes every curated accession and exact-sequence hash from T-c. Phase C rules `a` (accession and hash), `b_cluster_mate` and `c_test_taxon` apply | extend **both** `test_tc_excludes_test_proteins` (`test_keyword_tier.py:29`, `test_phasec_splits.py:144`) |
| Basidiomycota curated row in training | Basidiomycota stays a test clade. `splits.build` raises `SplitError` when a protein belongs to a training source and a test source | add a test on the new rows; keep the role in `species.tsv` |
| Homology across clades. Phase C allows cluster-mates of test proteins in training by design (ruling C-4). 3 of the 16 current positives have 0.3 or more identity to training (reviewer) | Report T-b rows by maximum-identity stratum, as for GO rows | report format |
| T-b and T-a are not independent, because A2 starts from the T-a genes | Report the `T-b` only stratum | report format |
| Positive and negative with the same sequence | Existing dedupe drops hashes in both classes | existing dedupe test |
| GPI rows bias toward yeast | Report P-gpi by clade. Do not pool | report format |

## 7. Acceptance and statistics

**The rule (ruling C-8).** A test set is an "estimate" when the recall half-width is at most 0.10 for
every one of R2, M8, M35, M8-C, M35-C and H, and the set has at least 20 direct positives. Otherwise it
is a "smoke test". The rule applies per test set. "Direct positives" means `n_direct_positives` in
`metrics.json`. That count excludes `pm-unresolved` genes (for example CDA1, CDA2, CDA3 and PLB1).

**Size, binomial approximation.** n = 1.96² x p(1-p) / 0.10², rounded up. The candidate with recall
nearest 0.5 needs the most positives.

| Recall p | Positives needed |
|---|---|
| 0.125 (R2 on Basidiomycota today) | 43 |
| 0.750 (M8 on Basidiomycota today) | 73 |
| 0.5 (worst case) | 97 |

At today's recalls of R2 and M8, the binding value is **73**. The recalls of the other candidates move as
new truth arrives, so the binding value can reach 97. The Basidiomycota sets are S3 sets. The ratio of
the cluster-bootstrap half-width to the binomial half-width is 0.88 to 1.14 for S3 sets (reviewer). The
required number scales with the square of that ratio. At p = 0.75 it ranges from about 57 (0.88² x 73) to
about 95 (1.14² x 73). The number is an approximation, not a measured Basidiomycota value.

**Today (re-checked):** clade block 16 positives, half-widths H 0.107, R2 0.158, M35 0.219, M8 0.222; the
20-positive floor is not met. H99 alone has 7 positives (worst half-width 0.43). *U. maydis* alone has 9
(worst 0.33).

**Acceptance for workstream A.**
- Step 01 merges the curated Basidiomycota truth, with a source for every row. It writes
  `curated_unmatched.tsv` and `curated_conflicts.tsv`.
- Phase C reports each Basidiomycota test set by tier stratum and by the `surface_evidence_htp_only`
  sub-stratum, with cluster-bootstrap intervals.
- The report states "estimate" or "smoke test" by the rule above, for the whole test set.
- If the clade block is a smoke test, the model card says Basidiomycota is not validated. Under decision
  10, this outcome does not allow a release (section 10).

**Acceptance for workstream B.**
- `curated_gpi.tsv` has rows, each with a resolvable PMID, and the new match check passes.
- P-gpi is a list in every test set that has fewer than 20 direct P-gpi positives (decision 12). In a set
  that has 20 or more, P-gpi is scored and carries a C-8 label.
- `d8_curated_conflicts.tsv` exists, and the owner has reviewed it.
- The report states the P-gpi counts by test set and clade.

**Both.**
- Metrics can move. Basidiomycota rows are test-only and do not change training. Curated GPI rows change
  training labels only in the training species. By default they move a *S. cerevisiae* or *C. albicans*
  gene from `pm-unresolved` (excluded) to P-gpi (positive). With `override_tm=yes` they also move PM-TM
  genes. The 2026-10-02 run is the baseline. The report shows before and after values for each candidate.
  The golden metrics file (`PHASEC_WRITE_GOLDEN=1`) covers the test fixture only. Only intended fixture
  paths change.
- No accuracy statement enters a README or card unless it comes from `metrics.json`.

## 8. Compute and effort

**Compute (measured).** C1 steps 09, 10 and 11 took 61 s, 62 to 91 s and up to 723 s. Steps 08 and 12
were not timed. Step 09 was timed once. The 1 to 1.5 hour sizing rule is for fan-out jobs and does not
apply to one job.

**Rerun chains.** The plan confirms both chains against the hash checks.
- *Curated Basidiomycota rows (or a `species.tsv` role change).* Step 01 reruns, because its inputs
  change. The chain then runs 01, 02, 03, 04, 05, 07, 08, C1 and 12 (reviewer). Step 03 needs the UniProt
  network. The J1 and J2 jobs resume by hash. The proteome FASTA files of H99 and *U. maydis* are already
  embedded. So a curated protein in those files has an embedding, unless the FASTA release changes
  (reviewer).
- *Only `curated_gpi.tsv` changes.* `truth_set.tsv.gz` does not change, so `truth_set_sha256` does not
  change (re-checked). Step 03 produces a new `truth_set_triaged.tsv.gz`. Step 08 checks that step 07
  read that file, so the chain runs 03, 07, 08, C1 and 12 (reviewer). Steps 02, 04 and 05 are not needed.

**Effort.** I have no estimate of the effort. The size of the Basidiomycota and GPI literature is not
measured.

## 9. Risks

| Risk | How it is detected |
|---|---|
| FungiDB stays unreachable | A1 reports it. A2 does not depend on A1 |
| Too few experimental Basidiomycota positives | the counts in section 7; the smoke-test label stays |
| The 20-positive floor is not met per source | the clade block is the label unit today (decision D-B) |
| *Cryptococcus* capsule and secreted proteins blur "wall" and "extracellular" | `subset` column; recall per subset |
| High-throughput-only evidence labels a protein positive | `surface_evidence_htp_only` sub-stratum in the report |
| Curated negatives are rare, so FPR stays imprecise | report the interval; no gate from it |
| One curator introduces bias | second review (section 4); `reviewer` column |
| Curated rows duplicate GO rows | tier strata `T-a`, `T-b`, `T-a,T-b`; report the `T-b` only stratum |
| Training labels change and metrics move | baseline run; before and after table |
| A curated GPI row matches no gene and is ignored | the new match check and `curated_gpi_unmatched.tsv` |
| A stale accession halts step 02 | step 01 writes `curated_unmatched.tsv` before the merge |
| A curated GPI row turns a TM negative into a positive | `override_tm=no` by default; `d8_curated_conflicts.tsv` for owner review |
| A predicted TM feature blocks correct curated GPI evidence | the owner reviews `d8_curated_conflicts.tsv` and sets `override_tm` |
| A typo in an evidence code or GO term passes silently | step 01 stops on invalid codes, terms and level mismatches (section 3.4) |

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
7. P-gpi stays a list (path b). Decision 12 sets the condition.
8. A curated GPI row does not override a UniProt TM feature by default. The owner reviews each case and
   sets `override_tm=yes` to allow it.
9. Revision 3 of this spec got a third independent review. Revision 4 applies it.
10. **Keep Q8 as written (D-A).** A release needs an estimate-grade Basidiomycota set. If the time box
    ends below 73 direct positives, the owner makes a new ruling. The card flag alone does not satisfy
    Q8.
11. **High-throughput-only rows (was a wrong R-A claim).** Genes with `surface_evidence_htp_only=yes` stay
    in the headline. Reports show them as a separate sub-stratum.
12. **Amend ruling R-B.** P-gpi stays a list in a test set until that set has at least 20 direct P-gpi
    positives. Then P-gpi is scored there and gets a C-8 label.
13. **Predictor-selected genes are a report-only sub-stratum (owner, 2026-10-03).** A gene whose every
    supporting curated row has `selected_by_predictor=yes` is left out of the headline recall and reported
    beside it. Trade-off: the headline stays free of the circularity with SignalP-based gates (R0, H), and
    the headline loses positives (in the pilots, 2 of 12 H99 P-ext genes and 8 of 15 *U. maydis* P-ext
    genes). The pilot counts ignore GOA overlap and `pm-unresolved` exclusions. The count toward 73 uses the
    headline positives only.
14. **Overexpressed fusion rows stay, with a flag (owner, 2026-10-03, option B).** Add the `overexpressed`
    column. Reports show recall with and without genes whose every supporting row is `yes`.
    The Cpl1 row was dropped (permeabilisation not stated, faint supernatant signal). Whether such genes
    leave the headline is open (D-E).

**Open (from the reviews).**
- **D-E. Do genes supported only by `overexpressed=yes` rows leave the headline?** Today they stay, with a
  with/without report. In the pilots 3 *U. maydis* genes depend on such rows (Sts2, Xyn2, Xyn11A)
  once predictor-selected genes are removed. Recommendation: keep them in, flag them, and decide after the
  with/without report shows how much they move recall.
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
| E3 | `curated_gpi.tsv` rows with PMIDs; `curated_gpi_unmatched.tsv` and the new match check with its test; `d8_curated_conflicts.tsv` reviewed by the owner |
| E4 | Code, step 01: `--curated` input, refactored `build_truth_rows`, `count_rows` on merged records, input checks, sequence check, `curated_unmatched.tsv`, `curated_conflicts.tsv`, new truth columns, curated file hash in `extract_log.json`. Code, step 03: `literature_ids` dictionary, `classify_pm` TM rule and `override_tm`, `d8_curated_conflicts.tsv`, header check. Phase C: `dedupe.merge_group` joins, `08_build_eval_tables.py`, `11_evaluate.py` (copy list, tier strata, htp sub-stratum), `12_report.py`. `species.tsv` roles. D10 removal. Tests: section 6, plus changes to `test_d8.py:334-341`, `:142` and `:338` |
| E5 | Phase C rerun, baseline comparison table, report sections per tier stratum, per source and per clade |
| E6 | `STATUS.md` update with new Basidiomycota and P-gpi numbers and intervals |

Done means: E1 to E6 exist; the owner has the report; the model card states the clade scope and the
validation state of Basidiomycota and GPI proteins.
