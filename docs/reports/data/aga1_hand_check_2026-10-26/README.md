# Hand-check of the AGA1-class candidate sample, 2026-10-26

The 97 proteins of `docs/reports/data/aga1_candidate_tiers_2026-10-26/hand_check_sample.tsv` (57 tier 1, 25 tier 2, 15 tier 3; ten proteomes) are
to be judged by the owner. This folder holds the evidence and the table for the verdicts. Nothing here is a verdict until the owner fills it in.

## What to edit and what to read

* **Edit `verdicts.tsv`** (97 rows). Fill `verdict`, `confidence`, `notes`, `reviewer`, `date` for each protein; do not change the other columns.
* **Read `evidence_sheets.md`**, one section per protein (automatic numbers, UniProt annotation, literature lead). `evidence_leads.tsv` and
  `evidence_auto.tsv` hold the same information as tables.
* **Check your file** at any time (it reports how many verdicts are still empty and stops on a value outside the vocabulary or an unknown accession):

```
python3 -I analysis/aga1_class/build_hand_check_sheet.py check --verdicts docs/reports/data/aga1_hand_check_2026-10-26/verdicts.tsv \
  --sample docs/reports/data/aga1_candidate_tiers_2026-10-26/hand_check_sample.tsv
```

## Verdict vocabulary (proposed; change it before you start if you want different meanings)

| verdict | meaning |
|---|---|
| `supported` | direct evidence that it is a cell-wall or surface adhesin-like protein of the flocculin/ALS/EPA/Hyr/Flo11/GPI-mannoprotein type |
| `plausible` | the family, domains and features fit, but no direct evidence for this protein |
| `not_wall_protein` | evidence that it is not (cytoplasmic or nuclear, a transporter, a transmembrane protein with no wall role ...) |
| `enzyme` | a catalytic enzyme (for example a GH18 chitinase or GH16/GH17 glucanase), whether or not it sits in the wall |
| `unknown` | cannot be judged from the evidence |

`confidence` is free text or high / medium / low as you prefer. The verdicts become the label table for the second version of the labels and for
recalibrating the thresholds per clade; the `enzyme` verdicts are the starting point for a signed list of enzyme families.

## How the evidence was made

* `evidence_auto.tsv` (`analysis/aga1_class/build_hand_check_sheet.py build`): the sample table plus numbers from the feature tables (composition,
  tail hydrophobicity, repeat z-score), nothing else.
* `evidence_leads.tsv` (`... leads`): leads drafted by four agents from public sources only, one batch per species group: UniProt (function, similarity,
  GO, keywords) and PubMed. They give no verdicts. **The leads are leads:** most UniProt entries here are unreviewed (TrEMBL), so function and GO text is
  often transferred from other species, and a "no literature found" for an uncharacterized protein mostly means that no gene-specific search turned
  anything up, not that nothing exists (the agents ran at most about three PubMed searches per protein).
* **Citations were checked in one way only.** All 23 PMIDs cited in the leads exist in PubMed, and their first author, year and topic match the way the leads
  use them. Whether each paper says what the lead claims about the specific protein was not independently checked; read the paper for any
  protein that matters to your decision. Seven rows have gene-specific lead confidence "high" (for example *C. parapsilosis* ALS3, ALS11 and ALS6 through locus-based
  papers; *C. glabrata* AWP1; *C. albicans* SCW11).

## Inputs not in the repository (sha256)

```
7f1eb913d176824eb55eeab8a80dda4200c6019ebbe817f6f3ff1facd4988c68  features_proteomes.tsv
bc3417c0c08e36e293822a3f2bcf26046b8e570aa478f894a0bfc70fd4bb71c6  features_nonyeast.tsv
63329781f124f9431640a9e9f34d944e9e9bb32a52bcb574935258372ad55749  features_cpar.tsv
78f991e3b30da18f6e4bba0aefaa3443d46da70436d736d32056a893f58a6ef8  leads_A_Cpar.tsv
bbb9affaa4d9a676952b32c0d3120a0ea01adf07ffe11aff5139334b67b0194b  leads_B_Cgla.tsv
297b91a8e8e705720e900f00a24e1a5010a3ff52bc8cc302d6a6eb72d9df492b  leads_C_yeasts.tsv
774ade25aaec39977b28d020b1ada772ade5966eee11bd3ee2180a201bfc4b11  leads_D_nonyeast.tsv
```

## sha256 of the committed files except `verdicts.tsv`, which the owner edits (a test checks every line)

```
dbf591d1b7c9fe65c211bbd2b4a202e7362744cdfe934de8ee172f9755ea0f9c  evidence_auto.tsv
142c35993ed06d4a759377b15e05911dc98d9164092eb454295d72b5aee64909  evidence_leads.tsv
8752f260c03d1772187fd16fa449213c0bd82be767eafa653e328edbfa9b7252  evidence_sheets.md
```
