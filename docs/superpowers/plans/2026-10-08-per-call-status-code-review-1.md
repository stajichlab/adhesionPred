# Code review 1 of the per-call status implementation, with dispositions

*2026-10-08. One independent review (Opus, read-only; it applied 30 mutations to a copy of the code: 27 killed, 3 survived). Verdict: **ready, with minor fixes**. No blocker or major finding; no path found by which a stale, hand-edited or malformed call file raises a status or applies to a wrong protein, taxon, variant or run.*

| id | finding | disposition |
|---|---|---|
| F1 | A name check swapped for a length check survived | Fixed: test `test_a_renamed_read_with_the_same_count_is_stale_not_an_error` |
| F2 | Plan T2.5 said a read-back failure leaves the previous file; the read-back runs after the write | Fixed in the plan and the code comment: validation before the write is what protects the old file |
| F3 | The lock file was made before the variant was checked (`../x` made a lock file) | Fixed: `check_variant` runs first; test parametrized over `R0`, `../x`, `a/b` |
| F4 | A per-variant eligible call without `--variant` gave a misleading refusal | Fixed: `check_variant` first in `_check_call_run` (no eligible per-variant call is packaged, so unit-tested only) |
| F5 | A refusal named the path twice | Fixed, with a test |
| F6 | The report column `file` held `valid` or `stale: ...` | Fixed: renamed `validity`, asserted in `test_cli.py` |
| F7 | Spec R6 listed clusters and leakage that the report omits | Spec corrected |
| F8 | The per-variant loop in `call_eligible` cannot fail today | Comment added: guard if the literal-name rule is relaxed |
| F9 | Module identity is checked against `run.json` outside the lock, `reads` re-read inside it | **Accepted, not changed.** One writer per work directory is the supported use (plan D8). The window is a module record changing between two lines of one command. Noted in the handoff. |
| F10 | Leakage values defined twice | Fixed: one list, `LEAKAGE_VALUES` |
| F11 | Formatter split an `if` around a comment | Fixed |
| F12 | The spec's "unrelated threshold change does not make the file stale" never holds, because `config_sha256` is compared first | Spec corrected |

Deviations from the plan that the reviewer found harmless: the hook returns a pair and not a triple; the stale basis uses `; ` and not `:`; `write_call_status` has no `reads` parameter (see F9); `run.json` `call_status_sources` is a superset of the planned fields; the alias `_reads_of_call` was removed; the derived-status note appears only in the Call calibration section (decision D5).

After the fixes: 788 tests pass, 7 skipped.
