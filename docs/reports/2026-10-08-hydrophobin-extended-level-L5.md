# Hydrophobin extended level: L5 evaluation and aligner comparison

*2026-10-08. Branch `hydrophobin-validation` (local, not pushed). Freeze commit `f1986e2` (`analysis/hydrophobin_truth/freeze.json`). Spec revision 4, plan revision 2.
Tables: `docs/reports/data/sorting_hat/hydrophobin_ext/` (`l5_results.json`, `pfam_missed_clusters.tsv`, `aligner_comparison.tsv`). Leave-cluster-out numbers are labelled
`partial`: the author saw the Pfam-missed entries before the freeze. Negatives are assumed. Every number comes from the scripts in `analysis/hydrophobin_truth/`.*

## 1. Relaxed Pfam level (frozen: seven hydrophobin-class Pfam models, `hmmsearch --nobias`, 2.6 bits, R0 called, at least 8 cysteines)

| Result | Value |
|---|---|
| Pfam-missed T2 proteins called (of 9) | 8 (all but HFB3_HYPVG) |
| Pfam-missed clusters recovered (of 6) | 6, Wilson 95% [0.61, 1.00] |
| Recall rule (at least 3 of 6) | passed |
| Clusters left unrecovered (u) | 0 |
| T2 called by the strict level / by the extended level | 122 of 131 / 130 of 131 |

The strict level is the domain-row Pfam call of the tool (122 T2). The sequence-level `--cut_ga` table also lists PSH_FLAVE (123), whose sequence score passes and whose domain score does not.
Seven T2 entries are lost to the conditions: six have fewer than 8 cysteines (A0A0A2K253, A0A0A2K7M3, A0A2H4SHM4, J5JTB1, P52755, W8NQ71) and RodD (Q4WE22) has no R0 call. The strict branch keeps six of them.

## 2. HMM decision (spec 6.2)

u is 0, so the custom HMM cannot add recall that the primary test can detect: **not testable on the primary set**. Reported for information: the leave-cluster-out HMM calls 76 of 131 T2 proteins with default
filters (the frozen option for the HMM; 58%, Wilson [0.49, 0.66]) and recovers 5 of 6 Pfam-missed clusters; with `--nobias` it calls 118 of 131 (90%) and recovers 6 of 6. The HMM is not part of the v1 call unless the owner decides otherwise.

## 3. Hard negatives (test parts, assumed negatives)

Relaxed Pfam and the all-cluster HMM call 0 proteins in every group: CFEM 0 of 128, cerato-platanin 0 of 288, cell wall cysteine proteins 0 of 238, small secreted cysteine-rich 0 of 74, HsbA 0 of 165 (no threshold for HsbA, E11).
All groups with a threshold pass the 5% rule. Strict Pfam GA also calls 0.

## 4. Jensen 2010 proteins (secondary; predictions)

Of 37 scored (reserved clusters and 8 unverified proteins excluded), Pfam GA misses 3: AFLA_014260, AO090012000143, ATEG_08089. Relaxed Pfam and the all-cluster HMM recover all three. The frozen-spacing rule recovers none.
This is a consistency check on genome-screen predictions, not a measurement.

## 5. Aligner comparison (exploratory; does not change the freeze)

The question: would `famsa` or `muscle` align this family better than `mafft`? The same pipeline was run for `mafft` L-INS-i (frozen), `famsa` 2.4.1 and `muscle` 5.1: 28 leave-cluster-out HMMs plus the all-cluster HMM, cutoffs by the frozen rule on the tuning data, held-out scoring.

| Aligner | Option | T2 recall (of 131) | Pfam-missed clusters (of 6) | Pfam-missed proteins (of 9) | Jensen missed recovered (of 3) | Hard-negative max rate | Folds passing slots / ordinal / same-column (of 29) |
|---|---|---|---|---|---|---|---|
| mafft L-INS-i | default | 76 (0.58) | 5 | 6 | 3 | 0 | 29 / 28 / 27 |
| mafft L-INS-i | `--nobias` | 118 (0.90) | 6 | 8 | 3 | 0 | 29 / 28 / 27 |
| famsa | default | 77 (0.59) | 6 | 7 | 3 | 0 | 29 / 28 / 27 |
| famsa | `--nobias` | 117 (0.89) | 6 | 8 | 3 | 0 | 29 / 28 / 27 |
| muscle 5.1 | default | 74 (0.57) | 5 | 6 | 3 | 0 | 29 / 28 / 27 |
| muscle 5.1 | `--nobias` | 117 (0.89) | 6 | 8 | 3 | 0 | 29 / 28 / 27 |

Reading: the three aligners give the same result within one or two proteins out of 131, and the same alignment diagnostics (the same two folds fail the stricter checks with every aligner). The aligner is not what limits the HMM here.
The search option is: `--nobias` raises leave-cluster-out recall from about 0.57-0.59 to about 0.89-0.90 with every aligner. The HMMER composition-bias filter appears to hide cysteine-rich hydrophobins. This is a result on held-out members of the same 131 proteins, so it is exploratory. Changing the HMM option or the aligner would be a new freeze, tested on data not used here (the reserved Jensen clusters and any new experimentally supported hydrophobins).

## 6. Limits

All numbers use assumed negatives and T2 labels that were not reviewed beyond a keyword check. The Pfam-missed set is 9 proteins in 6 clusters, so the recall test is coarse. The proteome cost on the test proteomes is not measured yet (task L7a).
