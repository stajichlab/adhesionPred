# Ablation power calibration (2026-11-27)

Benchmark `adhesion_v1`, training pool `all`, bar `composition_length` (best of the declared family by inner cluster-grouped CV).

## What the benchmark actually weighs

- 50 positives and 75 hard negatives in 84 clusters across the scorable folds.
- **Effective sample size 16.4 positives** (33 clusters), not 50: paralogue families are one observation, not many.

## Does the decision rule mean what it says?

Rejection rate under a pure-noise component, nominal alpha = 0.05:

| component dims | rejection rate | mean delta under the null |
|---|---|---|
| 1 | 0.006 | -0.0067 |
| 5 | 0.003 | -0.0206 |
| 20 | 0.000 | -0.0582 |

A rejection rate above alpha means the rule is anti-conservative and every later component test is easier to pass than it claims. A non-zero mean delta means a component is rewarded or penalised for its feature count alone, which is why the null for a real component must be a permutation of that component, never the number zero.

## Minimum detectable effect at 80% power

| signal | effect | delta in stratified concordance |
|---|---|---|
| Gaussian, same in every clade (optimistic bound) | 1.48 | +0.041 |
| Gaussian, full in Saccharomycetes, 0.33x elsewhere (realistic) | 1.62 | +0.029 |
| domain call at 5% FPR (sensitivity needed) | 0.66 | +0.048 |

The binary row is the actionable one: it says what sensitivity a domain model would need, at the stated false-positive rate on hard negatives, before this benchmark could see it.

## How to read a null result

A component that fails to clear this bar has not been shown to be useless. It has been shown to be undetectable at this benchmark's size. Report the interval's upper bound, not the phrase "no improvement".
