# AGA1 broad flag candidate, minimum disordered run 200 (and 100 as sensitivity), 2026-10-24

Made by `analysis/aga1_class/broad_flag.py` from `docs/reports/data/aga1_disorder_2026-10-24/aga1_signalp_predgpi_aiupred.tsv`
(report: `docs/reports/2026-10-24-aga1-broad-flag-candidate.md`). The flag is yes when SignalP 6 signal peptide, PredGPI GPI call
and an AIUPred disordered run of at least N residues are all yes (three-valued; the Ser/Thr stalk rule is not included). N is the
script's `--min-disordered-run`, default 200 (owner decision of 2026-10-24); 100 is the pre-declared sensitivity value. The value is
in every row (`broad_min_run`).

```
python3 -I analysis/aga1_class/broad_flag.py --record <the AIUPred record table> --min-disordered-run 200 \
  --out broad_flag_run200.tsv --summary-out summary_run200.tsv      # and the same with 100
```

| File | Content |
|---|---|
| `broad_flag_run200.tsv`, `broad_flag_run100.tsv` | one row per protein: parts, run length, `broad_min_run`, `broad_flag` |
| `summary_run200.tsv`, `summary_run100.tsv` | yes / no / NA counts per `role_group` |

## sha256

```
2192177e388f65fe18a7ea2e182c03c7ecaf6372f4aa02727427885af86e1492  broad_flag_run200.tsv
427e158df3d7b57dbd549056778cf33b40811cba9b288d656ee35fa077d6bb77  broad_flag_run100.tsv
c0ba80487b48b7eedbbb6cca294af527572af7504adb20bd86a009fa73b347c6  summary_run200.tsv
66c3fbe2ddf0908c76d0a10c2016ba36ca7072240e4ea44ab2246ce48f9e89e7  summary_run100.tsv
```
