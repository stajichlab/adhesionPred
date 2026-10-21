# Task 3.10 re-run after issue #93

*2026-10-21. Author: Claude Code (claude-sonnet-5-5). Follows `2026-10-18-task-3-10-refresh-results.md`. Each job was submitted with the owner's go-ahead.*

PR #96 added `species_limited` modules and raised `ENGINE_SEMANTICS` to 3. That changed `config_sha256` and every call hash, so the call status files of 2026-10-17 were stale. This re-run measured whether the change moved any result.

## What was run

Code commit `7ed2065` (clean clone on HPCC), `categories.yaml` sha256 `f50e1684…`, family table unchanged. The old `out_pre`, `out_final` and the two repeat status files were moved aside with the suffix `_before_93`.

| job | what | result |
|---|---|---|
| 30057968 | `out_pre`, 9 proteomes | all COMPLETED; 0 differing cells against `out_pre_before_93` |
| 30058018 | repeat call status, S288C and *C. albicans* | COMPLETED; sensitivity and specificity match the expected values; status `smoke`; no `.rejected` file |
| 30219753 | `out_final`, 9 proteomes | all COMPLETED; 0 differing lines against `out_final_before_93` in `calls.wide`, `calls.long`, `evidence`, `proteins` |

## Findings

- No count changed. Every proteome has an `antigen_lookup` table (job 2b), so the new rule has nothing to act on. This is the expected result.
- The repeat call status is `smoke` with the same measure, now with `config_sha256` `f50e1684` and `call_hash` `8644ca87`. The status is valid for the current config.
- Not tested by this re-run: the case the rule exists for, a proteome with no antigen table. The unit tests of PR #96 cover it. A run without job 2b was not made.

## Files

`data/sorting_hat/runs_2026-10-21_post93/` (run records) and the two `call_status/*/tandem_repeat_protein.owner_definition.config-of-2026-10-21-post93.json` files. On HPCC the `*_before_93` items, the clone, the wheelhouse and `$NEW` are kept until the owner approves a cleanup.
