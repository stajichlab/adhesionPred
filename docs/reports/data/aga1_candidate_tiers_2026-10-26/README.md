# AGA1-class tiered candidate list and hand-check sample, 2026-10-26

Built from the per-protein flags of ten proteomes (`docs/reports/data/aga1_proteome_test_2026-10-25/`, `aga1_nonyeast_test_2026-10-26/`, `aga1_cpar_test_2026-10-26/`,
file `all_proteins_flags.tsv` in each) by `analysis/aga1_class/candidate_tiers.py`, with public UniProt annotations from
`analysis/aga1_class/fetch_uniprot_annotations.py` (UniProt release 2026_03). Nothing here is a verdict: the tiers are rules, the sample is
the set to be reviewed, and the annotation columns are there for the reviewer. Decisions: memory note of 2026-10-26 (owner grilling).

## Tier rules (frozen)

* **Tier 1**: narrow rule, or (broad rule and at least 2 of 3 models). **Tier 2**: not tier 1, and (broad rule, or at least 2 of 3 models with a
  signal peptide). **Tier 3**: not tier 1 or 2, a signal peptide, and among the 25 highest mean model ranks of its species.
* Tiers 1 and 2 together are capped at 30 novel proteins per species; tier 1 is never cut, tier 2 is cut by mean rank and the cut proteins
  move to tier 3 with the reason "tier 2 overflow".
* The mean rank is the mean of the three model scores' percentiles **within the species**.
* "Novel" means not among the 850 labelled proteins; those proteins keep their tier for reference (`review_eligible` 0) and are not sampled.
* **Hand-check sample**: all novel tier 1 proteins, plus 25 novel tier 2 and 15 novel tier 3 proteins, taken species in turn in sha256 order of
  `seed:tier:accession` (seed 20261026). 97 proteins: 57 + 25 + 15 (`sampling_record.txt` has the counts per species).

## Commands

```
python3 -I analysis/aga1_class/candidate_tiers.py assign --flags <three all_proteins_flags.tsv> --out tiers.tsv \
  --accessions-out accessions.txt --record-out sampling_record.txt
python3 -I analysis/aga1_class/fetch_uniprot_annotations.py --accessions accessions.txt --out uniprot_annotations.tsv
python3 -I analysis/aga1_class/candidate_tiers.py annotate --tiers tiers.tsv --annotations uniprot_annotations.tsv --out candidates_tiered.tsv
# hand_check_sample.tsv = the rows of candidates_tiered.tsv whose `sample` column is not empty (same column order)
```

## Files

| File | Content |
|---|---|
| `tiers.tsv` | tier, reason, review eligibility and sample membership per protein, before annotation (467 proteins) |
| `uniprot_annotations.tsv` | public UniProt annotation per protein: Pfam, InterPro, domain names, keywords, subcellular location (no sequence) |
| `candidates_tiered.tsv` | the two joined: the tiered candidate list |
| `hand_check_sample.tsv` | the 97 proteins to review |
| `sampling_record.txt` | seed, cap, quotas and counts per tier and species |

## sha256 of the committed files (a test checks every line)

```
98cfa9befcb344e3bf58fd72f2fe9168d52c5cf2bb38060d72be5cd7cfc9d598  candidates_tiered.tsv
e95df3eb80d364586c60c2466ca7dadb04cbabbb5e4242ff801f331f4948a811  hand_check_sample.tsv
744919e1e353159fd9473141ac302dcaada45ef6a6ffaccaa42d3d44ebeb0cf7  sampling_record.txt
62333588b482c5a17dd885d79eca254b9d9bb0a40ad1a6daebb41f25c7cb1d88  tiers.tsv
9b8e1f7fbc6dd99b1f1dd2e42501420e5e6ab4e281cd429d11c1385cac1aaffd  uniprot_annotations.tsv
```
