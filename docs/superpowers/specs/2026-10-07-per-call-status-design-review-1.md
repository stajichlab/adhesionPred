# Review 1 of the per-call status spec, with dispositions

*2026-10-07. One independent review (Opus, read-only, code and spec only, no tests run). Verdict: **needs rework before a plan**. The overall design is sound. Two blockers, eight major findings, six minor findings. I checked F1, F2, F3 and F4 against the code and the stored run files (engine.py:296-308, `modules_of_call`, `modules/pfam_allergen.json`, `out/run.json`). They hold. The other findings are the reviewer's and were not re-checked one by one.*

Dispositions: **A** = applied in revision 2 of the spec. **O** = needs an owner decision (the revision shows a proposal).

## Findings

| id | severity | finding | disposition |
|---|---|---|---|
| F1 | blocker | A record with no contributors (a `not_assessable` value) gets `unvalidated` today. Step 1 of the resolution would give it the call entry's status. | A. The call status applies only when `res.contributors` is not empty. Test and mutation added. |
| F2 | blocker | `modules_of_call` skips `kind: other`, so `reads` would be empty and the validity check always passes. Composite calls (`ref`) were allowed while "composite measured as a whole" was out of scope. | A. `--call-status` and the loader accept only calls with an `expr`, no `ref`, and not `kind: other`. |
| F3 | major | "A read module has no identity" is wrong. `load_modules` records an identity before it checks the state; `_check_run_identity` would accept an `unavailable` module. | A. Check `run.json` `module_states` for every read module. Refuse any state other than `ok`. |
| F4 | major | `truth` uses the packaged config. The run may have used `--config`. | A. `truth` gets `--config`. It refuses when `run.json` `config_sha256` differs. `call_hash` comes from that config. |
| F5 | major | "The call's expression" is not defined for the hash. | A. Section 3.1 defines the canonical object. |
| F6 | major | Contributors are `(module, leaf_call)` with no variant. A per-variant leaf file needs the variant. | A. Key `(leaf_call, variant)`; derivation stated. Test that an R0 file does not apply to R1. |
| F7 | major | A call entry at species S would beat a module entry at a strain below S. "Most specific wins" breaks across sources. | A, through Q4: one source per call. |
| F8 | major | `load_status_source` does not check status against measure. A hand-edited call file could raise a status. | A. The call-file loader requires `measure`, checks `status <= status_from_measure`, the file name, and the variant rule. |
| F9 | major | A bad status file stops the run (`InputError`). A file left after a call is renamed would stop every run. | O. Revision 2 proposes: refuse, with a message that names the file. See Q6. |
| F10 | major | The stored R0 files are under `_workdir/`, which is git-ignored. The golden test cannot run in CI. | A. Commit a fixture copy and pin the expected `calls.long`. |
| F11 | major | `outputs.py` is missing from the file list. `_reads_of_call` should move to `engine.py`. `--module` is `required=True`. | A. Section 5 updated. |
| F12 | minor | Basis string format not given. | A. Exact strings in 3.2. |
| F13 | minor | Concurrent `truth` runs can lose an update. | A. `fcntl` lock on the file; documented. |
| F14 | minor | No identity covers the engine semantics. | A. `ENGINE_SEMANTICS` constant in `call_hash`, bumped by hand. |
| F15 | minor | Unclear whether the call table appears without call files. | A. Only when call files exist (keeps `report.md` unchanged). |
| F16 | minor | A read module absent from the workdir is not a listed case. | A. Stale: "no module run record". |

## Answers to the spec's questions

- **Q1 (composite).** Yes, a leaf call status stands in for its modules. The engine already works this way for single-module leaves: contributors carry the leaf call name, so a module measured on `signal_peptide_protein` counts inside `cell_wall_adhesion_candidate`. Decision 5 as implemented is "weakest of the leaf-call statuses that decided it". The spec now says so, and extends it to multi-module leaves by grouping contributors by `(leaf_call, variant)`. Record: `docs/paper/03` decision 3 wording ("a call with no measurement of its own reports unvalidated") does not match the code for composite calls. **O**: fix the paper text. See Q7.
- **Q2 (unavailable module).** Refuse, based on `run_state`. Measuring with `pfam_allergen` unavailable measures a different rule, and the status would be wrong as soon as a family is activated.
- **Q3 (hash).** Not enough as written. Revision 2 defines it (F5), takes the config from the run (F4), hashes only the thresholds the call references, and adds `ENGINE_SEMANTICS` (F14).
- **Q4 (both sources).** A conflict is possible: `signal_peptide_protein.R0` can have a module file and a call file. Proposal: `--call-status` only for calls that read two or more modules and have no `ref` and no `other`; single-module calls keep `--module`. One source per call remains. **O**: see Q8.
- **Q5 (failure table).** Rows added in revision 2.

## Missing tests the reviewer named (now in the spec)

Status applied to `not_assessable` records; variant ignored in the lookup; `kind: other` or composite accepted; `call_hash` not following `ref` targets; no `run_state` check; no config hash check; a loader that accepts a status stronger than its measure; `_merge_entries` keeping old entries when only `call_hash` changed; a call entry overriding a more specific module entry. "Ignore `reads`" is split into "ignore names" and "ignore identities".

## Backward compatibility

The reviewer found no existing test or stored file that changes behaviour, provided: `evaluate` gets a new optional keyword argument; `status/calls/` is a subdirectory (the resolver loads only `status/<module>.json`); the report table appears only when call files exist. The argparse change to `--module` changes only the message for a missing `--module`.

The reviewer could not find a record of owner decision 8 in the code. It is in `docs/HANDOFF-2026-10-07.md`.
