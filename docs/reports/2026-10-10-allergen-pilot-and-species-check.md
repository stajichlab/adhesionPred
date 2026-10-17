# Allergen predictor: pilot, full MMseqs2 runs, and species check (2026-10-10)

Status: pipeline stages T0.2 to T0.5 have run on HPCC. The T0.6 build gate is not yet done. No model is trained.
Research use only. Nothing here is an allergenicity call.

## 1. Inputs

| Set | Rows | Build | Notes |
|---|---|---|---|
| Positives | 147 | `build/02_positives.py` | 99 IUIS, 33 KW-0020 curated, 15 similarity-only. 4 inactive UniProt entries rejected. 3 merged accessions mapped. |
| Candidates | 6,624 | `build/01_candidates.py --exclude positives.tsv` | Reviewed Swiss-Prot fungi, spec 3.2 filters. 3 accessions removed because they are positives (see section 4). |

## 2. Pilot (T0.4)

- The first two pilot jobs failed (29797130, 29797276). Causes: a SIGPIPE bug in the sampler, and the AVX2 MMseqs2 build crashing with SIGILL on node c03.
- Fixed: the sampler writes to a file before `head`, and all MMseqs2 jobs use the non-AVX2 build (17-b804f).
- The pilot then completed (29797481). On 1,000 sampled candidates: cluster 8.8 s, search 4.7 s, 6,213 pair rows.
- Projection to 6,624 candidates: search about 31 s, cluster about 117 s with 2x margin.

## 3. Full runs

| Job | Input | Result |
|---|---|---|
| Cluster (29798025) | 6,771 sequences | 4,011 clusters. Integrity PASS: one row per sequence, cluster IDs are members. 42 clusters hold both a positive and a candidate (expected under the pair-in-cluster exception). |
| Pairs (29798026) | 6,624 candidates vs 147 positives | 40,263 pair rows. Pairs at fident >= 0.6 with query and target coverage >= 0.5: 156 rows. At >= 0.7: 88 rows. |

## 4. Structural test that caught a build defect

`tests/species_check.py` check S1 (no accession is both a candidate and a positive) failed on the first run for O00089, P12547 and P41748.
- O00089 and P41748 are current primaries of merged accessions, which the candidate query did not map.
- P12547 is a WHO/IUIS allergen that UniProt does not tag as an allergen, so the candidate query returned it.
- Fix: `01_candidates.py` now excludes every positive accession (`--exclude`). The corrected set is 6,624. Clustering and pairs were rerun on it. S1 passes.

## 5. Species check (focus genera)

Positives, candidates, and candidates that pair with a positive (coverage rule applied):

| Genus | Positives | Candidates | Any pair | Pair at T_C 0.6 | Pair at T_C 0.7 |
|---|---|---|---|---|---|
| Aspergillus | 34 | 373 | 36 | 14 | 8 |
| Candida | 5 | 252 | 37 | 4 | 4 |
| Penicillium | 22 | 90 | 12 | 3 | 1 |
| Alternaria | 17 | 1 | 0 | 0 | 0 |
| Arthroderma | 15 | 7 | 3 | 1 | 0 |
| Trichophyton | 5 | 6 | 0 | 0 | 0 |
| Malassezia | 14 | 16 | 2 | 0 | 0 |
| Coccidioides | 0 | 2 | 0 | 0 | 0 |
| Paracoccidioides | 0 | 3 | 0 | 0 | 0 |
| Blastomyces | 0 | 1 | 0 | 0 | 0 |
| Ajellomyces (*Histoplasma*) | 0 | 2 | 0 | 0 | 0 |
| Emmonsia | 0 | 0 | 0 | 0 | 0 |
| Histoplasma | 0 | 0 | 0 | 0 | 0 |

Full table: `analysis/allergen_model/tests/species_out/species_check.md`.

## 6. Answers to the questions asked

- **Coccidioides or other Onygenales candidates?** Yes, but few. Coccidioides has 2 candidates (*C. posadasii* E9CXX8; *C. immitis* RS P0CI02). *Paracoccidioides* has 3, *Trichophyton* 6, *Arthroderma* 7, *Blastomyces* 1, *Ajellomyces* (*Histoplasma*) 2. None of the Onygenales candidates pair with a positive at the coverage rule. *Emmonsia* and *Histoplasma* have none. Because Onygenales are poorly represented in IUIS (task 05 notes only four Trichophyton proteases), zero pairs here does not show absence of allergens.
- **Aspergillus fumigatus and relatives:** 373 candidates from *Aspergillus*. 14 pair at T_C 0.6. One candidate (P67876) is identical (fident 1.000) to *Asp f 1* (P67875), a positive. That is the label-noise case the spec predicts: an unannotated allergen among the weak negatives. P34755 (0.989) pairs with P34754. The model must not treat these as clean negatives.
- **Alternaria:** 1 candidate, no pairs. The 17 Alternaria positives are the source of the signal, not the candidates.
- **Candida albicans:** 252 candidates, 37 pair with some positive, 4 at T_C 0.6. The 5 Candida positives include Q5AF03 and Q6YK78 at 0.98, so these are the closest cases.

## 7. What this does and does not show

- The pipeline runs end to end on HPCC for the clustering and the pairing search, with the pilot working.
- It shows that the candidate set is not clean: at least 14 Aspergillus and 4 Candida candidates pair with known positives at the T_C 0.6 rule, and one matches an allergen exactly.
- It does not show model performance. No model is trained. The label-noise limit in the spec applies.

## 8. Open items

- Decision on candidates that match a known positive at high identity (for example fident >= 0.9). Options: remove them from the weak-negative pool, or keep them and report the count as a limit. The spec currently keeps them and reports the limit. This needs an owner decision.
- The T0.6 build gate counts pairable candidates per arm (500 minimum). The counts above are the inputs to that gate.
