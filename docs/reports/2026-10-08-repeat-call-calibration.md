# The repeat call has a measured status: per-call status entries and the first calibration

*2026-10-08. Branch `per-call-status` (local, not pushed). Research use only. Everything here is exploratory and the leakage cap applies: the best status any entry can reach is `smoke`.*

## 1. What was built

`tandem_repeat_protein` reads two modules (`repeat02` or `repeat14`). A status belonged to one module, so this call could never be calibrated (decision 3 of 2026-10-06). Per-call status entries remove that limit. Spec: `docs/superpowers/specs/2026-10-07-per-call-status-design.md`; plan and independent reviews: `docs/superpowers/plans/2026-10-08-per-call-status*.md`; rules: `docs/paper/03`, section 6a.

| piece | where | tests |
|---|---|---|
| call helpers (`reads_of_call`, `call_eligible`, `call_hash`) | `engine.py` | `test_call_helpers.py` (17) |
| call status file loader, validity, resolver | `call_status.py` | `test_call_status.py` (36) |
| writer with merge, lock, validation before write | `calibration/call_files.py` | `test_call_files.py` (14) |
| engine hook (leaf-call grouping, stale fallback) | `engine.evaluate` | `test_engine.py` (+9) |
| core command, report section, `run.json` field | `cli.py`, `outputs.py` | `test_cli.py` (+5, and a pinned full-output golden test) |
| `calibrate truth --call-status` | `calibration/cli.py` | `test_truth_call_status.py` (17) |

The suite went from 683 to 782 tests passed (7 skipped: shellcheck is not installed). The output of a run without call files is pinned and unchanged. Each task has mutation checks (about 40 in all) that fail named tests. A call file is never used when the config, the modules read, their identity or run state, or the call definition differ from the measurement.

## 2. The truth set

Per species, one row per protein, in `analysis/calibration_truth/repeat_call_truth/` (`truth_v2.<proteome>.tsv`, annotated versions, and the scripts that build them).

- **Positive:** a paper statement of a repeat region (adjudicated curation table), or two or more UniProt `Repeat` features.
- **Negative:** no UniProt `Repeat` feature and no Region, Compositional bias or Domain feature that mentions a repeat. **This is an assumed negative: an absent annotation is not proof.** UniProt lacks features for some repeat proteins (SOWgp, BAD1).
- **Population:** the 155 curated adhesin-table proteins of the two species, plus the reviewed UniProt proteins that have a signal peptide (344 for S288C, 162 for *C. albicans*), mapped to the run proteome by exact sequence. Family inference never counts as a positive (owner decision 5).
- **Clusters:** MMseqs2 at 30% identity and coverage 0.5, as in the C1 count.
- UniProt release 2026_03, fetched 2026-10-08.

## 3. Measured status of the repeat call

Cluster bootstrap, 95% intervals, leakage `tuned_on_truth` (FLO1, FLO5, FLO9, FLO10 and the Als proteins tuned the detectors).

| species (taxon) | positives (clusters) | negatives (clusters) | sensitivity [95%] | specificity [95%] | status |
|---|---|---|---|---|---|
| *S. cerevisiae* (4932) | 31 (22) | 297 (227) | 0.323 [0.043, 0.535] | 0.987 [0.962, 1.000] | `smoke` |
| *C. albicans* (5476) | 16 (9) | 210 (155) | 0.375 [0.000, 0.680] | 0.976 [0.938, 0.995] | `smoke` |

An earlier version on the 155 curated proteins only gave S288C 0.474 [0.111, 0.724] / 0.960 [0.800, 1.000] (19 positives, 14 clusters) and *C. albicans* 0.385 [0.000, 0.737] / 1.000 [0.952, 1.000] (13 positives, 6 clusters). Files: `docs/reports/data/sorting_hat/call_status/`.

How to read it:
- **Specificity is high** (at least 0.976 in both species, upper bound 1.0), but the negatives are assumed.
- **Sensitivity is low and the intervals are wide.** The call finds about a third of the positives.
- **The status cannot be `estimated`.** Sensitivity half-widths are about 0.25 to 0.34. The leakage cap would hold it at `smoke` in any case.
- *C. albicans* has 9 positive clusters (the Als family forms one or two), against a floor of 20. S288C has 22, past the cluster floor, but not the half-width rule.

## 4. Why sensitivity is low

Of the 31 S288C positives, 18 have no periodicity in either detector (period 0). The same holds for 8 of 16 in *C. albicans*. I looked up the UniProt `Repeat` features of these proteins (feature names and counts, fetched 2026-10-08).

| group | proteins | what UniProt records |
|---|---|---|
| **Numbered tandem arrays that the detectors miss** (S288C) | AGA1 (20 repeats, median unit 7 aa), HPF1 (18, 13 aa), SED1 (9, 14 aa), EGT2 (9, 35 aa), MSB2 (7, 17 aa), CNE1 (8) | repeats of 7 to 35 aa; three of them cover 0.6 to 0.84 of the protein |
| ... with only two repeats | KRE1, SAG1, CCW12 | 2 numbered repeats each |
| **Numbered tandem arrays that the detectors miss** (*C. albicans*) | EAP1 (25 repeats, median unit 6 aa, 0.60 of the protein), PGA18 (24, 8 aa, 0.41), ALS5, ALS6, ALS7 (Als repeats, 4 to 6 copies of about 32 aa, 0.16 of the protein) | |
| **Named repeat domains** | PEP1, VTH1, VTH2 (BNR), HRD3 (Sel), SCJ1 (CXXCXGXG), FMP27 and SPS22 (LRR), PGU1 (PbH), PIR5 (PIR); DSE1, ASC1, TUP1 (WD) | domain-type repeats |

So the misses are about half arrays and half repeat domains. The detectors are tandem-array detectors, so the domain-type misses are by design. **The array-type misses are not**: AGA1, HPF1, EAP1 and PGA18 have 18 to 25 short repeats and get no period. Detector scores explain part of it (`repeat14`, `rep_z_seq` against its gate `Z_MIN = 4.0`): EAP1 3.4, AGA1 3.7, SED1 3.7 (below the gate); ALS5 5.1, ALS6 5.4, EGT2 5.4 and MSB2 4.8 pass the gate but no repeat region is extracted (region score 0), and even with a region ALS5 would cover 0.16 of the protein, below the call's coverage cutoff of 0.25.

Two exploratory scans (nothing in the module was changed):
- **Coverage cutoff** (`threshold_scan_v2.tsv`). Lowering 0.25 to 0.15 gains a few positives (S288C 0.323 to 0.387; *C. albicans* 0.375 to 0.438) and costs 1 to 2 false positives. Lowering it further adds little.
- **An added sequence-level rule** (`zseq_rule_scan.tsv`): also call a protein when `rep_z_seq` is at least Z. At Z = 5.0 it gains EGT2 in S288C and ALS5 and ALS6 in *C. albicans*, and adds 1 and 4 false positives (S288C Q8TGR2; *C. albicans* PGA62, IFF5, RBR3, IFF9). CGD describes IFF5, RBR3 and IFF9 as GPI-anchored adhesin-like proteins, so several of these "false positives" are probably real repeat proteins that UniProt has no features for. That is the assumed-negative problem, and it means the specificity numbers above are lower bounds on the true specificity in an unknown amount. I do not propose a rule change on this evidence: a new rule needs a held-out truth set and the owner's decision.

Leads for the detector (not done): the short-unit arrays (EAP1, PGA18, AGA1; 6 to 8 aa units) and the Als repeats of ALS5, ALS6 and ALS7; FLO11 (Ser/Thr-rich, weak period-15 signal, coverage 0.04) is a known miss.

## 4a. *A. fumigatus*: no status written

The same method on the 167 reviewed, secreted UniProt proteins of *A. fumigatus* Af293 (taxon 330879; `truth_v2.Afum_Af293_UniProt.tsv`) gives 8 positives in 4 clusters and 159 negatives in 117 clusters. Result: sensitivity 0.000 [0.000, 0.490], specificity 1.000 [0.968, 1.000], no negative called. **I did not write a status entry for it.** All 8 positives are enzymes with repeat domains: seven PbH (pectate-lyase-like beta-helix) repeats in polygalacturonases and a xylogalacturonan hydrolase (pgaA, pgaB, pgaX, pgxB, pgxC, xghA, AFUA_1G17), and one BNR repeat set in Vps10. They are not surface arrays, so "sensitivity 0.0" would say nothing about the arrays the call is for, and the report would show a misleading number.

What this does show: the repeat call called none of 159 reviewed secreted *A. fumigatus* proteins. Whether it finds tandem arrays in filamentous fungi (SOWgp, BAD1 and CspA are the known cases) is **not measured**. The repeat call has a measurement only in two Saccharomycotina species. UniProt has no reviewed *A. fumigatus* array protein to test it on. The curated adhesin tables hold a few (CspA among them); a truth set for filamentous fungi needs curation (agent tasks 03 and 08).

## 5. Effect on the other calls

With the call files in place, a re-run gives statuses (S288C, 6,722 proteins; the run has the R0 status source):
- `tandem_repeat_protein`: `smoke` for 6,714 proteins (the 8 invalid sequences stay `unvalidated`).
- `cell_wall_adhesion_candidate[R0]`: 9 called proteins are `smoke` (decided by the repeat leg and R0), 4 called proteins are `unvalidated` because the PA14 leg also decided them and its module has no status. Not-called records are `smoke` only when every deciding leg is measured.
- *C. albicans* (6,212 proteins): 15 called candidates `smoke`.

These composite statuses are derived values (weakest measured leaf), not measurements of the composite.

## 6. What this does not show

- It does not show that the call separates adhesins from other wall proteins. The earlier curated-label check (`2026-10-07-classifier-runs-and-scale.md`) found that it does not on untuned proteins; this measurement counts repeat proteins of any kind.
- Negatives are assumed. Positives include repeat domains. Both inflate the apparent difficulty or ease in ways I cannot size.
- The truth proteins are the curated adhesin set plus reviewed UniProt proteins. Reviewed proteins are the better-studied ones.
- Leakage: the detectors were tuned on Als, Flo and SOWgp. The cap holds the status at `smoke`.

## 7. To reach `estimated`

Independent positive clusters (at least 20, and about 60 for a half-width of 0.10 at sensitivity 0.3, by the binomial approximation that ignores cluster structure), a truth set that no detector saw, and a decision on what "tandem repeat" means for the call (arrays only, or any repeat). The curation tasks 03 and 08 supply positives. A narrower truth definition (tandem arrays at a stated period) would give a sharper measure.
