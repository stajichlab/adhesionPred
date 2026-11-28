# Ablation power calibration (PR 0)

Answers one question, before any component of the project is tested:

> **How large an improvement would a component have to make before `adhesion_v1` could tell?**

If that minimum detectable effect is large, the finding is that the benchmark is too small, and
curating more truth outranks building an ablation harness. Deciding that cheaply, before building
the harness, is the entire purpose of this directory. **No real component is tested here.**

## Why this exists before the harness

An independent review implemented the original 2b harness design and measured it. Ten findings were
blocking; four were re-verified independently. The three that shaped this code:

| Finding | Measured | What it changed |
|---|---|---|
| Pooling out-of-fold scores across folds is incoherent | length-only AUC is 0.912 in Saccharomycetes and 0.327 in Eurotiomycetes — two folds learn "long ⇒ adhesin", the third the opposite; 47% of pooled pairs compare contradictory scores | primary metric is **within-fold stratified concordance**, which never compares across folds |
| The specified baseline is beaten by one of its own inputs | log length alone AP 0.690 vs 0.551 for the 21-feature model | the bar is the **best of a declared family**, with C tuned by inner cluster-grouped CV |
| ΔAP = 0 is the wrong null, and the error scales with feature count | pure noise is worth +0.007 at 1 dimension and −0.025 at 20 | the null is a **permutation of that component**, and the dimensionality bias is measured and reported |

## What it does

1. Loads the frozen benchmark, refusing to run unless it matches its own `SHA256SUMS`.
2. Builds leave-one-clade-out splits under the **strict cluster rule** (D-A): when a fold is held
   out, every protein whose MMseqs cluster touches that fold is dropped from training. Without it,
   five clusters that span folds leak homologues at ≥30% identity into training.
3. Chooses **the bar** — the best of length-only / composition-only / composition+length — by inner
   cluster-grouped CV inside each training part, and records the C it settled on.
4. Measures the **type-I rate** of the decision rule under pure-noise components of 1, 5 and 20
   dimensions. A rule whose nominal 5% is really 9.5% makes six later component tests a coin flip.
5. Measures **power curves** for three planted-signal shapes and reports the minimum detectable
   effect at 80% power.

## The three signal shapes, and why

- **homogeneous** — same effect in every clade. The most favourable case that could occur, so its
  curve is an optimistic bound, not an expectation.
- **heterogeneous** — full effect in the largest fold, attenuated elsewhere. Leave-one-clade-out
  cannot see a signal that does not transfer, because the model scoring a held-out clade never saw
  it. This is the realistic case for a domain family with narrow taxonomic range.
- **binary** — a domain call with a stated sensitivity on positives and false-positive rate on hard
  negatives. **This is the only row a biologist can act on**: it says what sensitivity a domain
  model would need before the benchmark could see it at all.

## Declared before running

Fixed in this README, in the code, and in the plan, so that reading the result cannot change them:

- **Metric** — within-fold stratified concordance. Average precision is reported descriptively,
  always beside its fold's prevalence, because it is prevalence-dependent and the folds sit at
  0.36 / 0.65 / 0.31.
- **Decision** — a one-sided paired permutation test by cluster, α = 0.05. Clusters are the
  exchangeable unit because a paralogue family is one observation, not nine.
- **Target power** — 80%.
- **Training pool** — stated with `--train-pool`, never implied. It changes the leakage drop counts
  (Pezizomycotina 5 with `all`, 2 with `scored`) and brings in ~26 positives from clades that are
  never tested. Default `all`; run both as a sensitivity.
- **Effective sample size** — reported as Kish ESS. The benchmark has 50 positives but an effective
  **16.4**, because ALS alone is 9 of the 33 positive clusters. Quote the effective number.

## Running

```bash
python3 analysis/ablation_power/run_power.py \
    --benchmark data/benchmark/adhesion_v1 --out /path/to/NEW_DIR
```

**Not `python3 -I`**, unlike the stdlib-only scripts elsewhere in this repo. Isolated mode implies
`-s`, which drops user site-packages, and on the HPCC numpy, scipy and scikit-learn live in
`~/.local/lib/python3.12/site-packages` -- with `-I` this dies on `import numpy`. Checking that the
imports work is not enough; they must be checked *with the flag the run will use*. Because the
result depends on packages outside the repo, their versions are recorded in `baseline.json`.

A few minutes on a login node, one core. The output directory must not exist. Writes `baseline.json` (the bar, the splits,
the per-fold report, the effective sample sizes, every setting), `type_i.tsv`, `power.tsv` and
`summary.md`.

## How to read a null result

A component that fails to clear this bar has **not** been shown to be useless. It has been shown to
be undetectable at this benchmark's size. Report the interval's upper bound — "rules out an
improvement above X" — never the bare phrase "no improvement". Treating "no evidence against" as
"evidence for keeping" is how the cruft this project just removed accumulated in the first place.

## Known gap

One bug found during development is worth remembering, because it is easy to reintroduce: the bar
and the component-augmented model **must be fitted at the same regularisation**. Fitting the bar
with a tuned C while components ran at a fixed C=1.0 made a pure-noise component appear to cost
0.15 concordance. `choose_bar` returns its per-fold C for exactly this reason, and two tests pin it.
