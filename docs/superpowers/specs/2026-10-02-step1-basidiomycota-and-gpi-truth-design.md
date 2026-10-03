# Design spec: Basidiomycota truth and curated GPI rows for step 1 validation (issue #50)

*Drafted 2026-10-02. DRAFT, owner decisions 1 to 6 recorded (section 10). It needs an independent review (a different model) before any plan or
code. No code, data or job exists for this spec. Owner decisions are in section 10.*

Inputs: `docs/superpowers/specs/2026-09-30-surface-glycoprotein-model-design.md` ("step 1 spec",
cited by section and by decision Q1 to Q10, rulings R-A to R-C); `docs/model-review/STATUS.md`;
`docs/HANDOFF-2026-10-02.md`; `analysis/step1_compare/` (README, `species.tsv`, `curated_gpi.tsv`,
`COLUMNS.md`).

## 1. Purpose and non-goals

**Purpose.** Step 1 predicts surface localisation (cell wall, outer face of the plasma membrane, or
secreted). Two parts of its validation are weak or missing:

1. **Basidiomycota.** The test set has 16 positives (*Cryptococcus* 7, *Ustilago* 9) and 60 negatives.
   It is a smoke test (step 1 spec section 6 rule: a set is an estimate only when the recall interval
   half-width is 0.10 or less). Decision Q8 says: curate Basidiomycota truth before the first release.
2. **GPI-anchored proteins.** `analysis/step1_compare/curated_gpi.tsv` has a header and no rows. The
   P-gpi label is a list, not a scored stratum (ruling R-B). GPI-anchored protein performance is not
   measured.

This spec defines how to add curated truth for both, and how that truth enters the existing pipeline.

**Non-goals.**
- It does not change the label definitions (P-ext, P-gpi, N-int, N-sec, PM-TM, ambiguous; step 1 spec
  section 2.2).
- It does not train or choose a model. Gate values stay unset until after measurement.
- It does not cover step 2 or step 3 (adhesion mechanism, purpose layers), or the #14 hard-negative
  panel for step 2.
- It does not change `surface.tsv` or the UniProt-keyword tier T-c.

## 2. Facts checked

| Fact | Value | Source |
|---|---|---|
| Basidiomycota sources in the repo | H99 (`test_clade`), JEC21 and CRYD1 (`alternate_file`), *U. maydis* MYCMD (`undecided`) | `species.tsv` |
| H99 non-IEA cellular-component labels | 77 genes; 11 P-ext (CDA1-3, LAC1, LAC2, PLB1, CPL1, QSP1, PQP1, YOR1, cnap1); direct-evidence P-ext 9 | step 1 spec 3.3, 3.4 |
| JEC21 | 32 P-ext, all surface evidence IBA; no experimental P-ext gene; IEA fraction 0.520 | step 1 spec 3.3, 3.3 text |
| *U. maydis* direct truth | 10 P-ext / 6 N-int / 21 N-sec | step 1 spec 3.3 |
| Current Basidiomycota test set | 16 positives / 60 negatives; R0 recall 0.938, FPR 0.083 | STATUS.md |
| UniProt 2026_03 reviewed entries with experimental GPI evidence (ECO:0000269) | S288C 3 (GAS1, TIP1, YNL190W); *C. albicans* 3 (HWP1, YWP1, ECM331); *S. pombe* 0; *A. fumigatus* 0; *A. nidulans* 1 (chiA) | step 1 spec 2.2 |
| `curated_gpi.tsv` columns | `source_id`, `gene_id`, `symbol`, `pmid`, `note`; no data rows | the file |
| FungiDB downloads | returned HTTP 401 on 2026-09-27 | issue #14 body |
| Read depth of this spec | I read the step 1 spec in full, and the headers of `species.tsv`, `curated_gpi.tsv` and `COLUMNS.md`. I did not read `d8_triage.py`, `labelmap.py` or the Phase C evaluation code | - |

Not known: how many Basidiomycota cell-wall or secreted proteins have experimental evidence in the
literature; how many curated negatives exist; whether FungiDB still gives H99 GO annotations that
differ from the GOA file.

## 3. Design

Two workstreams. Each ends in a tracked table that the existing extraction and evaluation code reads.

### 3.1 Workstream A: Basidiomycota truth

**Step A1: FungiDB H99 check (first, small).** Compare the H99 GO cellular-component annotations in
FungiDB with `313589.C_neoformans_var_grubii_H99.goa`. Report: genes with wall or extracellular
terms that the GOA file lacks, and the evidence codes. If FungiDB is not reachable (HTTP 401 on
2026-09-27), record that, and go to A2. The result is a count, not a label change.

**Step A2: Literature curation of H99 (and optionally *U. maydis*).**
- Positives: proteins with experimental evidence of wall, surface or secreted location. Start from the
  11 H99 P-ext genes and the 9 *U. maydis* direct-evidence genes. Add proteins found in the
  literature (for example capsule-associated and wall proteins). Each row has a PMID and an evidence
  statement. I did not search the literature, and I give no expected count.
- Negatives: proteins with experimental evidence of an internal location (cytosol, nucleus,
  mitochondrion) and no surface evidence (N-int), and proteins in the secretory pathway with no
  wall or extracellular evidence (N-sec). Every negative has a PMID.
- Ambiguous cases (surface and internal evidence together) follow decision Q9: excluded from
  precision, reported as a stratum.

**Step A3: Entry into the pipeline.** Curated rows go in a new tracked file
`analysis/step1_compare/curated_basidiomycota.tsv` (columns in section 5). D9 reads it and merges it
with the GO-derived Basidiomycota truth. Each row keeps a `tier` of `T-b` (literature rows), as in the
step 1 spec (evidence tiers for positives). The GO-derived rows stay `T-a`. Reports show the two tiers
beside each other.

### 3.2 Workstream B: curated GPI rows

**Step B1: Collect literature rows.** One row per protein with experimental evidence of GPI
anchoring (for example biochemical evidence of the GPI anchor, or a mutant of the GPI attachment
signal that releases the protein). Start with the proteins that the step 1 spec names: GAS1, TIP1,
YPS1, SAG1, CWP1, CCW12. Add proteins from *C. albicans*, *S. pombe*, *A. fumigatus* and
Basidiomycota where the literature has direct evidence. I do not have PMIDs for these proteins in the
repo. They must come from the literature search.

**Step B2: Rules.**
- Predictor output never counts as evidence (decision Q2).
- A reviewed UniProt entry with ECO:0000269 for the GPI-anchor lipidation counts. This is already
  used by D8.
- A protein enters P-gpi only when it is also P-ext and has a non-IEA plasma-membrane term (D8 rule,
  ruling R-B). Proteins with curated GPI evidence that are outside P-ext (GAS1 and TIP1 are
  ambiguous) stay in `d8_gpi_outside_pext.tsv` and are not relabelled.
- Rows go in `analysis/step1_compare/curated_gpi.tsv` (existing columns).

**Step B3: Scoring.** P-gpi becomes a scored stratum only when the number of P-gpi rows reaches the
threshold in section 7. Before that it stays a list.

## 4. Evidence rules (both workstreams)

| Rule | Detail |
|---|---|
| Source of each row | a PMID, or a reviewed UniProt entry with an experimental evidence code. A row with neither is not accepted |
| Evidence level | `direct` (the paper measures the location or the anchor in this protein), or `transfer` (the paper measures it in an ortholog). Only `direct` rows count for headline metrics, as in step 1 spec 3.3 |
| Reviewer | each row is checked by a second person or an independent model run before it enters truth. Unreviewed rows stay in a draft file |
| Accession | a UniProt accession from the reference proteome in `species.tsv`. A protein with no accession in that proteome is listed as unmatched and is not in truth (as in `unmatched_ids.tsv`) |
| No predictor labels | SignalP, PredGPI and TM predictions never set a label. They are features of the rule |
| Moonlighting | a protein with surface and internal evidence is ambiguous (Q9). It is not a positive |

## 5. Data format

`curated_basidiomycota.tsv` (new), one row per gene:

| Column | Meaning |
|---|---|
| source_id | Source in `species.tsv` (`Cneo_H99_GOA` or `Umay_MYCMD`) |
| gene_id | UniProt accession |
| symbol | Gene symbol, if any |
| label | `P-ext`, `N-int`, `N-sec` or `ambiguous` |
| subset | `wall` or `extracellular-only` for P-ext; empty otherwise |
| evidence_level | `direct` or `transfer` |
| pmid | PMID (several separated by `;`) |
| evidence_note | One sentence: what the paper measured |
| reviewer | Initials or model name of the second check |
| review_date | `YYYY-MM-DD` |

`curated_gpi.tsv` keeps its five columns. The spec adds no column. If the owner wants
`evidence_level` there too, that is an open decision (section 10).

Both files are plain text and small, so they stay uncompressed in git.

## 6. Leakage controls

| Risk | Control | Test |
|---|---|---|
| Curated test proteins are in the UniProt-keyword tier T-c | D10 removes every curated accession and exact-sequence hash from T-c (step 1 spec 2.3) | extend `test_tc_excludes_test_proteins` to the new files |
| A curated Basidiomycota protein is in training | Basidiomycota is a test clade only (Q6, Q7). No curated row enters training | `test_test_truth_sources` (extend) |
| A positive and a negative have the same sequence | existing dedupe drops hashes in both classes | existing dedupe test |
| Curation looks at model output | the curator does not see model scores. The curation list is fixed before the Phase C re-run | process rule; recorded in the PR |
| Labels copy the rule | no predictor output as evidence (section 4) | review checklist |
| GPI rows bias toward yeast | report P-gpi by clade; do not pool the clades into one number | report format |

## 7. Acceptance

The test sets are "estimate" or "smoke test" by the rule in step 1 spec section 6.

**Rough size needed (an approximation, not a measurement).** With a simple binomial interval and an
assumed recall of 0.8, a half-width of 0.10 needs about 62 positives (1.96 x sqrt(0.8 x 0.2 / n) <=
0.10). Cluster bootstrap on homologous families widens the interval, so the real number is higher.
The current 16 positives give a half-width of about 0.20 under the same assumptions. The assumed recall
is not a measured value for Basidiomycota.

Acceptance for workstream A:
- Curated Basidiomycota truth is merged by D9 with a source for every row.
- The Phase C evaluation reports Basidiomycota per tier (T-a, T-b) with cluster-bootstrap intervals.
- The report states if the set is an estimate or a smoke test, by the rule above.
- If the set is still a smoke test, the model card states that Basidiomycota is not validated.

Acceptance for workstream B:
- `curated_gpi.tsv` has rows, each with a PMID.
- P-gpi is scored as a stratum only if the row count allows an interval with half-width of 0.10 or
  less for recall. Otherwise it stays a list and the report says so.

Both workstreams:
- The golden metrics file is regenerated with `PHASEC_WRITE_GOLDEN=1`. Only intended paths change.
- No accuracy statement enters a README or card unless it comes from `metrics.json`.

## 8. Effort and unknowns

I have no estimate of the effort. The size of the literature for Basidiomycota and GPI evidence is not
measured. A first pass of A1 (the FungiDB check) is small. A2 and B1 are literature work, as the
handoff says. The Phase C re-run after the truth changes needs compute that I have not sized. Run it
as a normal CPU job and size it to one to 1.5 hours, as the HPCC rules require.

## 9. Risks

| Risk | How it is detected |
|---|---|
| FungiDB stays unreachable | A1 reports it. A2 does not depend on A1 |
| Too few experimental Basidiomycota positives | the count in section 7; the smoke-test label stays |
| *Cryptococcus* capsule and secreted proteins blur "wall" and "extracellular" | `subset` column; recall per subset |
| Curated negatives are rare, so FPR stays imprecise | report the interval; do not state a gate from it |
| One curator introduces bias | second review (section 4); `reviewer` column |
| Curated rows duplicate GO rows | dedupe by accession and hash; report the overlap with T-a |

## 10. Decisions for the owner

**Decided 2026-10-02: the owner accepted the recommendation for items 1 to 6.** The text of each item
is the decision.

1. **Which workstream first?** Recommendation: A1, then B1 in parallel with A2. A1 is small and sets
   the A2 scope. B1 does not depend on A.
2. ***U. maydis* in or out?** Its 10 / 6 / 21 direct-evidence truth already exists. Recommendation:
   keep it as a second Basidiomycota set and report it beside H99, not pooled.
3. **JEC21 IBA truth:** keep as a separate "homology-transfer" stratum (step 1 spec Q8 option 2), or
   drop it? Recommendation: keep it, report it separately, and give it no weight in headline metrics.
4. **Add `evidence_level` to `curated_gpi.tsv`?** Recommendation: yes, to match section 4.
5. **Who is the second reviewer?** The step 1 workflow used an independent model for review.
   Recommendation: an independent model run for the first pass, then the owner spot-checks a random
   sample. I have no data on the error rate of either.
6. **Stop rule.** Recommendation: stop curation when the section 7 size is reached, or after a
   time box that the owner sets. I do not propose a time box without effort data.

## 11. Deliverables

| ID | Deliverable |
|---|---|
| E1 | A1 report: FungiDB versus GOA for H99 (counts, codes), or a record that FungiDB was unreachable |
| E2 | `curated_basidiomycota.tsv` with sources and second review |
| E3 | `curated_gpi.tsv` rows with PMIDs |
| E4 | D9 and D10 changes that read the new files, with the tests in section 6 |
| E5 | Phase C re-run, regenerated golden metrics, report section per tier and per clade |
| E6 | `STATUS.md` update with the new Basidiomycota and P-gpi numbers, with intervals |

Done means: E1 to E6 exist; the owner has the report; the model card states the clade scope.
