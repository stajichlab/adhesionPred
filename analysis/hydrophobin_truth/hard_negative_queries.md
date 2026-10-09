# Hard-negative groups: queries and rules (task L3, written before any fetch)

2026-10-08. UniProt REST (`https://rest.uniprot.org/uniprotkb/search`), taxon 4751 (Fungi).

| Group | Source rule | Reviewed count (2026-10-08) |
|---|---|---|
| CFEM | `xref:pfam-PF05730` | 29 reviewed, 13,578 unreviewed |
| cerato-platanin | `xref:pfam-PF07249` | 18 reviewed, 2,178 unreviewed |
| HsbA | `xref:pfam-PF12296` | 3 reviewed, 3,881 unreviewed |
| cell wall, cysteine-containing | `xref:pfam-PF00399` (PIR) OR `gene:CCW12` OR `gene:CCW14` | 22 + 4 reviewed |
| small secreted cysteine-rich | reviewed, `keyword:KW-0964` (Secreted) OR `keyword:KW-0800` (Toxin), length 40 to 250, then at least 6 cysteines, no hydrophobin or rodlet in the name | counted after the cysteine filter |

Rules:
1. All reviewed entries of a group are used. For the four Pfam-defined groups, up to 300 unreviewed entries per group are added (seeded random sample, seed 20261008), so that a call rate has a usable denominator. Unreviewed entries are members by Pfam assignment, not by curation. The report says so.
2. A protein with a hit to any hydrophobin-class model (the seven Pfam models) at the gathering cutoff, or with hydrophobin or rodlet in its name, is removed from every group and listed in `removed_hydrophobin_like.tsv`.
3. Each group is clustered by itself (MMseqs2, 30% identity, coverage 0.5). Clusters are split 50/50 into a tuning part and a test part by seed. (Deviation from spec 6.1, which clustered the negatives together.)
4. The per-group pass threshold is the call rate on the **test** part, at most 5%. The HsbA group has no pass threshold (E11); its rate is reported.
5. Counts per group and part are printed before any scoring.

## Counts (2026-10-08, `make_hard_negatives.py`, `hard_negative_counts.json`)

| Group | Reviewed | Unreviewed sample | Clusters | Tuning part | Test part |
|---|---|---|---|---|---|
| CFEM | 29 | 300 | 91 | 201 | 128 |
| cerato-platanin | 18 | 300 | 13 | 30 | 288 |
| HsbA | 3 | 300 | 92 | 138 | 165 |
| cell wall, cysteine-containing | 26 | 300 | 49 | 88 | 238 |
| small secreted cysteine-rich | 132 | 0 | 67 | 58 | 74 |

No protein was removed as hydrophobin-like. Cerato-platanin has only 13 clusters, and one large cluster puts most of its proteins in the test part (30 in the tuning part). The tuning rate for that group rests on 30 proteins.
`hard_negative_thresholds.json` fixes the per-group pass rate (5% of the test part; none for HsbA).
