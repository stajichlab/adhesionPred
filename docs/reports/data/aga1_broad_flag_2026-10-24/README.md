# AGA1 broad flag candidate with the Ser/Thr stalk rule, 2026-10-24

Made by `analysis/aga1_class/stalk_features.py`, `stalk_grid.py` and `broad_flag.py` from the AIUPred record table
`docs/reports/data/aga1_disorder_2026-10-24/aga1_signalp_predgpi_aiupred.tsv` and the raw AIUPred file (which carries residue letters
and is not in the repository; its sha256 is in the disorder record README). Report: `docs/reports/2026-10-24-aga1-broad-flag-candidate.md`.

Flag = signal peptide yes AND PredGPI GPI yes AND a disordered run (AIUPred >= 0.5) of at least N residues that holds a stalk: a
30-residue window, wholly inside the run, that is at least 50% Ser/Thr (owner decisions of 2026-10-24: stalk rule 30 / 50%, N = 200
by default, 100 as the pre-declared sensitivity value). Every table carries its settings in the columns `stalk_rule`, `broad_min_run`
and `tail_rescue_min_hydro`. `nostalk` tables use the longest AIUPred run instead (the first candidate of PR 146). `rescue` tables add the
exploratory C-terminal tail rescue (see the report; off by default).

```
python3 -I analysis/aga1_class/stalk_features.py --raw <aiupred raw> --manifest analysis/aga1_class/manifest_accessions.txt --out stalk_features_w30_f0.5.tsv
python3 -I analysis/aga1_class/stalk_grid.py --raw <aiupred raw> --record <record table> --manifest <manifest> --out stalk_grid.tsv
python3 -I analysis/aga1_class/broad_flag.py --record <record table> --stalk-table stalk_features_w30_f0.5.tsv --out broad_flag_stalk_run200.tsv --summary-out summary_stalk_run200.tsv
# the others: --min-disordered-run 100 | --no-stalk-rule | --tail-rescue-min-hydro 1.4 or 2.0
```

| File | Content |
|---|---|
| `stalk_features_w30_f0.5.tsv` | per protein: longest stalk-holding run, best window Ser/Thr fraction, C-terminal tail hydrophobicity (no residues) |
| `stalk_grid.tsv` | counts per group over window 20/30/50 x Ser/Thr fraction 0.4/0.5/0.6 x run 100/200, plus the no-stalk setting |
| `broad_flag_*.tsv`, `summary_*.tsv` | the flag per protein and the yes / no / NA counts per group, for each setting |

## sha256

```
f48e66f340b8ba17852a0be8e9931226f13a0881ad8de19e2fa15d3559af9be2  broad_flag_nostalk_run100.tsv
a72ba5bfb758d70f2d7704c804de07838f7f40c75b8fb2b24c73de7fd12fa77a  broad_flag_nostalk_run200.tsv
7db66200055361e817eeb5cb1a999cbf9402c044907f4afed62b0ffb914732e6  broad_flag_stalk_run100.tsv
4c20fdae8128d37abb8e647fd84e4e4b540a58490341e0dd9cb725ea569eb0c7  broad_flag_stalk_run200.tsv
06c73f1fed69494df9a13094fc1eac46234aa4b119719a4c791602bd3c4f0621  broad_flag_stalk_run200_rescue1.4.tsv
df057c2bf9d98ec9fa490fb9949c9851a7c6340d5b8b990d130c0bb847f33738  broad_flag_stalk_run200_rescue2.0.tsv
4c14412a9726405e483f8d3d29ecd711c6b113ecbe6682ca5fe3dfab5a3fe7a0  stalk_features_w30_f0.5.tsv
789f0ecc48740fb9c58c42f9ec5586892e9dda9f5ff2cdc12d1fc904448fb80f  stalk_grid.tsv
66c3fbe2ddf0908c76d0a10c2016ba36ca7072240e4ea44ab2246ce48f9e89e7  summary_nostalk_run100.tsv
c0ba80487b48b7eedbbb6cca294af527572af7504adb20bd86a009fa73b347c6  summary_nostalk_run200.tsv
03eefe577ce14952d0a6ab4534246d9edf1971ed1c0ae139d94ac42283b03a1e  summary_stalk_run100.tsv
216fcd2f78d8e710b22a2ef6ff1d27ade10f9ef296b98ff6ceb2cf373ec44479  summary_stalk_run200.tsv
5143fffb132df309087e10a085935bd4c04e1ac29d279c1f67b3be66aaa807e4  summary_stalk_run200_rescue1.4.tsv
4de24bc09aa30b6e3e65ed8ef9d0384eb2aba815014dd9f40ec40ce7720d6abe  summary_stalk_run200_rescue2.0.tsv
```
