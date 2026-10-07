# Final decisions report (2026-10-06)

Commits: 37afd51 (code and tests, plus spec 3.4 limit 9), a5c00c9 (plan text), 969f067 (README, CHANGELOG, AGENTS).

Cluster count reproduction: eval_table.tsv.gz rows with origin=go, homology_only=no, class pos/neg, source_ids split on comma
(a protein can be in two sources: Afum_ASPFU,Anid_EMENI), joined on seq_sha256 to clusters.tsv.gz, distinct cluster_id per
source and class. Result (pos/neg clusters): Anid 100/151, Afum 17/38, Scer 58/3156, Calb 113/410, Cneo 6/31, Umay 9/24.
Protein counts equal metrics.json (109/164, 19/45, 79/3785, 153/459, 7/32, 9/28).

Real-file phasec run (scratch workdir, fixture R0 record signalp 6.0h-gpu, mode fast): Anid estimated (half-width sens 0.0935,
spec 0.0227); Scer smoke (sens 0.1040, spec 0.0072); Calb smoke (0.1177, 0.0278); Afum smoke (0.1308, 0.0459);
Cneo smoke (0.2708, 0.1022); Umay smoke (0.1496, 0.1442). Time 10.9 s, 2.1 GB.

## Judgment calls (not stated in the instructions)
1. SignalP check: first version was module-name-through-major-version-and-gpu-tag; now a `module=` token in `tools.signalp` is compared exactly, and the bare-version rule is the fallback (strict boundary, whole token gpu).
2. `phasec` notes also carry `variant=R0` and a sentence on the widest-of interval; the record version text goes into the notes with `=` and `;` replaced.
3. Spec 3.4 limit 9 was committed with the code (commit 1), so the spec-equals-code test passes at every commit.
4. `truth --variant X` leaves out step 1 modules of other variants when counting the modules a call reads.
5. `evaluate` has an optional `measured_call_of`; the call is the LAST `call=` match in the notes.
6. Contributors are (module, call) pairs; a module measured on a composite call does not count for its own leaf call.
7. The calibration table got a "measured on call" column (the notes were not printed before) and a `measured_call` key in run.json calibration rows.
8. `module_identities` is a sorted list of {name, version, params_hash, artefact_hash} for every module with a record or table, including error states.
9. Cluster counting refuses a protein with no cluster and one in two clusters; strata are not widened; `entries_from_phasec` signature gained `cluster_counts`, `extra_notes`, `where_counts`.
10. Plan: three Global Constraints bullets, measured cluster counts in the R0 row, `--variant` and one-module text in Task 13, supersede notes on old listings.
11. `truth` check order: module read by call, one-module rule, run.json identity, then the rest; all before any write.
