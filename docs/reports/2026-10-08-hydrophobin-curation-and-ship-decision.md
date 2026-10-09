# Hydrophobin unlabelled calls: evidence-based curation (tier T4) and the ship decision

*2026-10-08. Branch `hydrophobin-validation` (PR #77, draft). The owner could not curate `evidence_sheets.tsv`, so the 61 calls were curated from evidence and literature under `analysis/hydrophobin_truth/curation_rubric.md`, which was committed before any decision (`a9b7a32`).
Decisions are tier **T4: curator judgement from evidence, not owner-reviewed, not experimental**. They are kept in `curator_decisions.tsv`, apart from the owner columns, and never merged with T1/T2. The owner can audit and overrule any row.*

## 1. Method

- Evidence lines A to E as in the rubric, gathered by a separate agent (`analysis/hydrophobin_truth/curation/`: `curation_evidence.tsv`, `literature.tsv` (21 papers, 10 species), `literature_gene_map.tsv` (98 paper gene IDs mapped by sequence), `NOTES.md`, `scripts/`).
- Decisions applied by `apply_rubric.py` (tested, 7 tests): hydrophobin if C, or A with B or D; not hydrophobin if E and not C; a conflict between hydrophobin evidence and E is unresolved; otherwise unresolved.
- Line A (doublet architecture: 8 to 10 cysteines with the 2nd and 3rd and the 6th and 7th adjacent, R0 called, at most 300 aa) was tested on the known hydrophobins: it holds for 166 of 174 (161 yes; 5 fail only on length above 300; 8 have no such motif, among them RodD and the 7-cysteine HfbA).

## 2. Evidence counts (61 proteins)

| Line | Yes | Note |
|---|---|---|
| A doublet architecture | 24 | |
| B strict Pfam | 21 | |
| C named a hydrophobin in a paper or curated record | **0** | 2 partial mappings (below) |
| D reciprocal best hit with a known hydrophobin (40% identity, 70% coverage of both) | 6 | |
| E another function | 8 | 6 from a strong Swiss-Prot non-hydrophobin hit, 2 from non-hydrophobin Pfam coverage above 50% |

UniProt: 49 of 61 map at 95% identity and 90% coverage of both sequences. Twelve mapped TrEMBL entries are named "Hydrophobin", but all 12 names are rule or similarity based (RuleBase, ProtNLM, submitter), so they do not count as line C. FungiDB: 19 map to a product; none says hydrophobin. Sixteen proteins carry a non-hydrophobin Pfam domain.
No hydrophobin paper was found for *C. immitis*, *B. dermatitidis*, *C. albicans* or *S. cerevisiae*. Full text was read for 4 papers (Quarantin 2019, Valsecchi 2018, Luciano-Rosario 2022, Izumi 2026), the rest are abstract only.

## 3. Decisions (T4, not owner-reviewed)

| Group | hydrophobin | not hydrophobin | unresolved |
|---|---|---|---|
| strict and relaxed, unlabelled (18) | 16 | 0 | 2 |
| strict only, unlabelled (3) | 0 | 0 | 3 |
| **relaxed only, unlabelled (40)** | **0** | 8 | 32 |

- The 16 strict-and-relaxed proteins decided as hydrophobin hold line A (doublet architecture) together with the strict Pfam hit or an ortholog. A is independent of Pfam; B is not.
- The 8 relaxed-only proteins decided as not hydrophobin have a strong hit to an enzyme or other protein: cutinase, endoglucanase, PR5-like receptor kinase, two versatile peroxidases, SUN41 beta-glucosidase, an antimicrobial peptide, and one protein with an ARB_05566_N domain over 66% of its length (E rests on the Pfam coverage alone; PF28391 is a family named for an uncharacterised protein, so this is the weakest of the eight).
- **Unresolved is not a false positive and not a hydrophobin.** Under this rubric a relaxed-only call can become a hydrophobin only through C (none found) or A with D (none found), so 0 hydrophobins among the 40 means "no supporting evidence found", not "all wrong". The 32 unresolved include small cysteine-rich proteins with no hit at all.
- Two protein pairs are close to a named hydrophobin and fail the stated rule only on one-way coverage (the rule needs 95% identity and 90% coverage of both sequences): W72310 KAK9636619.1 (54 aa, 100% identical to RodD but 27% of RodD) and *P. expansum* FA144B67_000394-T1 (99.3% identical to HfbA, 90% of the query but 65% of HfbA). Both look like gene-model differences. They stay as the rubric decides (not line C); the owner may decide that one-way coverage is enough.
- Near misses to paper genes, not counted: W72310 KAK9641117.1 vs RodE 84.6%; *B. bassiana* F1BB8A46_000504-T1 vs Hyd2C 75.4%; PC9 F4DD442B_010567-T1 vs Hydph8 92.8%; PC9 F4DD442B_008436-T1 vs Hydph18 98.9% over 80% and 64%. The *P. ostreatus* UniProt entries that cite Xu 2021 are mostly from strain PC15, and the proteome here is PC9.
- Thresholds the evidence agent set itself (stated in `curation/NOTES.md`): line E "most of the sequence" is non-hydrophobin Pfam envelopes above 50% of the length; the one-to-one-ortholog clause of line D was not checked.

## 4. Ship decision (`ship_rule.py`, `ship_decision.json`)

| Condition | Result |
|---|---|
| Recall (at least 3 of 6 Pfam-missed clusters) | pass (6 of 6) |
| Cost in every test proteome (5 per 10,000) | **fail** (*B. bassiana* 10.55, *F. graminearum* 10.72; five others 1.0 to 1.9) |
| Hard negatives per group | pass (0 calls) |
| Precision of the extra calls (at least 5 resolved clusters, Wilson lower bound at least 0.5) | **fail**: the 40 relaxed-only calls form 36 clusters; 0 hydrophobin, 7 not hydrophobin, 29 unresolved; 7 resolved, lower bound 0.00 |

**The relaxed level is not shipped.** `hydrophobin_extended` is not added to `categories.yaml`. `hydrophobin_relaxed` stays as a module and a reported column.

## 5. Reading and options for the owner

- The relaxed level does what it was built for on labelled proteins: it recovers 8 of the 9 Swiss-Prot hydrophobins that the strict Pfam call misses (6 of 6 clusters). What fails is that, in whole proteomes, most of its extra calls cannot be shown to be hydrophobins, and a few are clearly other proteins.
- The strict level looks good on its unlabelled calls: 16 of its 21 unlabelled calls are supported by the independent doublet architecture, none has evidence of another function, 5 are unresolved. This is T4, and line A and the strict Pfam hit together carry most of that.
- Options (each would be a new version with its own pre-registration, not a reinterpretation): (1) keep the strict call only and report the relaxed module as an experimental column; (2) a stricter relaxed cutoff chosen to meet the cost limit, tested on new data; (3) require line A (doublet architecture) as an added condition to the relaxed level; the evidence suggests it separates the supported calls from the enzymes, but it was seen on these data and would need its own test.
- Not done: owner adjudication; the full-text read of the papers that are abstract-only; the T2 label review beyond keywords.
