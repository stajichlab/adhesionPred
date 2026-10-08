# Review 1 of the per-call status plan, with dispositions

*2026-10-08. One independent review (Opus, read-only; it ran the baseline tests: 683 passed, 7 skipped). Verdict: needs rework before implementation. 1 blocker, 7 major, 10 minor. I checked P1 against `cli.py` (`load_modules` records an identity before it checks the run state). All findings are applied in revision 2 of the plan (decisions D1 to D15) or decided there.*

| id | severity | finding | disposition |
|---|---|---|---|
| P1 | blocker | The reader checked names and identities but not run states, so a call status measured on two modules could apply when one module is `unavailable` | Applied (D10, T2.2): module states are passed to the resolver; reason `module state not ok`. Test and mutation added. |
| P2 | major | `call_eligible` was ambiguous for per-variant calls | Applied (D15, T1.2, T1.5) |
| P3 | major | A tenth known limit breaks `test_the_spec_and_the_report_list_the_same_known_limits`; no golden report text exists | Applied (D5): `_known_limits` unchanged; T4.6 |
| P4 | major | A static limits line changes every `report.md` (R5) | Decided (D5): the derived note is printed only when call files exist; the paper docs state it for all readers. A `status_basis` marker is left to the owner. |
| P5 | major | The golden comparison would be taken after Task 3 and compared gzip bytes; no fixture | Applied (Task 3a): golden files from the base commit, decompressed text, committed fixture |
| P6 | major | A hand-edited file could raise a status through taxon 0 or 1 or by editing the leakage cap | Applied (D11): taxa below 2 refused; `leakage` required in the notes and the cap applied at load |
| P7 | major | The lock test would rarely catch a missing lock | Applied (T2.5): deterministic probe through a patched `write_atomic` |
| P8 | major | No new-path test used `measured_call_of`, so a wrong leaf call would pass | Applied (T3.4) |
| P9 | minor | The core run would import numpy through `calibration.measure` | Applied (D1): `status_from_measure` moves to `status.py`; writer to `calibration/call_files.py` |
| P10 | minor | `resolve_entry` needs one identity; a call file has none | Applied (D1): `best_entry` helper |
| P11 | minor | `call_hash` over thresholds is redundant with `config_sha256` | Kept: the hash gives a clearer reason and covers `ENGINE_SEMANTICS`; documented |
| P12 | minor | `_reads_of_call` drops literal step 1 names | Applied (D15): calls that name a step 1 module literally are not eligible |
| P13 | minor | Identities must be written as strings | Applied (D11) |
| P14 | minor | No per-variant eligible call is packaged | Applied (T1.4): custom config |
| P15 | minor | `TOY_NODES` has no strain; the `write_run_json` helper lacks two keys | Applied (T2.3, Task 5 note) |
| P16 | minor | `report.md` contains the installed version | Applied (Task 3a) |
| P17 | minor | Only `unavailable` was tested | Applied (T2.2, Task 5) |
| P18 | minor | `fcntl` on a network file system | Applied (D8) |

Open items the reviewer listed are decided in D12 to D14.
