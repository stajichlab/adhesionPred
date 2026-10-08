# Per-call status entries: implementation plan

*2026-10-08. **Revision 2**, after independent review 1 (`2026-10-08-per-call-status-review-1.md`: 1 blocker, 7 major, 10 minor; all applied or decided below). Plan for `docs/superpowers/specs/2026-10-07-per-call-status-design.md` (revision 2; owner decisions Q6 to Q8 recorded). Branch `per-call-status`. Test first in every task. One commit per task. Nothing is pushed until the owner says so. Revision 2 has not been re-reviewed.*

Python 3.12. Run tests with `PYTHONPATH=$PWD/src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q -p no:cacheprovider`. The baseline is 683 passed, 7 skipped. The pre-commit hook runs ruff 0.3.5 and ruff-format and aborts a commit that it reformats: re-add and commit again.

## Decisions taken in this plan (the spec left them open or implicit)

D1. **Module layout (review P9, P10).** `status.py` gains `status_from_measure` and its constants (`MIN_POSITIVES`, `MIN_NEGATIVES`, `MIN_CLUSTERS`, `MAX_HALF_WIDTH`), moved from `calibration/measure.py`, which re-exports them (so no numpy import enters the core run), and a helper `best_entry(entries, taxon, lineage)` that `resolve_entry` and the call resolver both use. A new `src/cellsurface_sorting_hat/call_status.py` holds the loader, the validity check and the resolver; it imports `status` and `engine`, and nothing from `calibration`. A new `src/cellsurface_sorting_hat/calibration/call_files.py` holds the writer, merge and lock; it imports `call_status`, `status` and `calibration.measure`. `engine.py` gets the call helpers (`reads_of_call`, `call_eligible`, `call_hash`, `ENGINE_SEMANTICS`) and does not import `call_status`.
D2. **Hook shape.** `evaluate(..., call_status_of=None)`. The hook is `call_status_of(call, variant, taxon)`. It returns `None` (no file for this call, or the taxon is not tested: module logic, no prefix), `(status, basis)` (applies), or `(None, reason)` (a file exists and is stale: module logic with the prefix `call status stale (<reason>):` on that leaf's basis items).
D3. **Byte identity (R5).** When no hook is given, or when it returns `None` for every leaf of a record, `evaluate` takes the old code path unchanged. The new path runs only when at least one leaf has a call status or a stale reason.
D4. **Grouping.** Contributors are `(module, leaf_call)`. The leaf variant is the record's label if the leaf call is `per_variant`, else the empty string. A leaf with an applying call status contributes one status and one basis item (`call:<leaf>:taxon:<tested>`), placed where the leaf's first item would stand in the old sorted order. Every other pair goes through `status_for` as today. The record status is the weakest over all items.
D5. **Derived-status note (review P3, P4).** `_known_limits` is not changed (a test pins exactly nine limits against the orchestrator spec, and a static line would change `report.md` for every run, against R5). The note "statuses of composite calls are derived: the weakest status of the measured leaf calls that decided them (owner decision Q7)" is printed only in the "Call calibration" section, so only when call files exist. The paper docs state it for every reader (Task 6). The owner may later ask for a marker in `status_basis`; that is a separate change.
D6. **Loader strictness (review F8).** The call-file loader requires `measure` in every entry, checks `status` against `status_from_measure`, rejects an unknown call, an ineligible call, a name that does not match, a variant error, and a taxon in two entries. `load_status_source` (module files) is unchanged.
D7. **Eager load, hard errors (Q6).** `CallStatusResolver.__init__` loads every `status/calls/*.json` and raises `InputError` naming the file for any structural problem (unknown or ineligible call, bad file, wrong name). A well-formed but stale file never raises.
D8. **Lock (review F13, P18).** `write_call_status` takes `fcntl.flock` on `status/calls/<file>.lock` and holds it through read, merge, write and read-back. `/bigdata` is a network file system on which `fcntl` locks may be unreliable (`cache.py` says so); the lock is documented as best effort, with one writer per work directory as the supported use.
D9. **Migration.** None. No change to module files or `load_status_source`.
D10. **Run state in the reader (review P1, blocker).** `CallStatusResolver` receives `loaded.states`. A file is stale with the reason `module state not ok: <m> (<state>)` when any read module's state is not `ok`, in addition to the identity checks. Order of reasons: `config differs`, `reads differ`, `no module run record: <m>`, `module state not ok: <m> (<state>)`, `identity differs: <m>`, `call_hash differs`.
D11. **Loader limits (review P6).** The call-file loader refuses a taxon below 2 (a status for the root or taxon 1 would apply to every protein). It requires `leakage: <value>` in `measure.notes` (the writer always writes it) and caps the status: any value other than `none` allows at most `smoke`. Identity values in `reads` are strings (review P13).
D12. **Stale basis text (open item).** For a stale leaf the basis item is `<module>:call status stale (<reason>); <module logic basis>`. A hook result is `(status, basis)` (the tested taxon is inside `basis`), not a triple.
D13. **Report rows (open item).** `rows()` gives one row per entry and tested taxon, with the tested taxon in the `taxon` column, plus one row per stale file with `valid=false`.
D14. **Flags (open items).** `--call-status` is a flag; the call comes from `--call`. `--config` applies to both modes, but only the `--call-status` path checks `config_sha256`. The call path writes `notes` as `call=<call>; variant=<v>; reads=<m1,m2>; leakage: <value>` plus the user's text. The legacy case (a module file with `call=` notes for a call that also has a call file) is tested in Task 4: the call file is used and the report shows both the module calibration row and the call row.
D15. **Eligibility (review P2, P12).** For a `per_variant` call the rule `reads >= 2` is checked for every variant label in `step1_variants`; for a call that is not per-variant it is checked once. A call that names a step 1 module literally (not through `{step1}`) is not eligible (reason `literal step1 module`).

## Task 1: engine helpers

Files: `src/cellsurface_sorting_hat/engine.py`, `tests/cellsurface_sorting_hat/test_call_helpers.py`.

Produces:
- `ENGINE_SEMANTICS = "1"` (bumped by hand when the evaluation rules change).
- `reads_of_call(cfg, call, variant="")`: `sorted(modules_of_call(cfg, call))`; with a variant, a step 1 module of another variant is removed (moved from `calibration/cli._reads_of_call`; that function becomes a thin alias).
- `call_eligible(cfg, call)` -> `(bool, reason)`, per D15. Eligible: the call exists, has `expr`, is not `kind: other`, its expression has no `ref` node (walk), names no step 1 module literally, and reads at least two modules for every variant it has. Reasons: `unknown call`, `kind other`, `contains ref`, `literal step1 module`, `reads fewer than two modules`.
- `call_hash(cfg, call, variant="")`: SHA-256 of canonical JSON (sorted keys, separators `(",", ":")`) of `{name, per_variant, variant, variant_module, expr, semantics}`. In `expr`, every `test` value `"$name"` is replaced by the number in `cfg.thresholds`, and `{step1}` in a module name by the variant module. Only thresholds the call references appear.

Tests (write first):
1. `reads_of_call` for `tandem_repeat_protein` is `["repeat02", "repeat14"]`; for `iuis_allergen_homolog` it is `["allergen_homology", "pfam_allergen"]`; for `signal_peptide_protein` with variant `R0` it is `["step1_rule@R0"]`.
2. `call_eligible`: true for the two calls above; false with the exact reason for `cell_wall_adhesion_candidate` (`contains ref`), `other_not_surface` (`kind other`), `signal_peptide_protein` (`reads fewer than two modules`), and an unknown name.
3. `call_hash`: stable; changes when a referenced threshold changes (copy the config to `tmp_path`, edit, reload); does not change when an unreferenced threshold changes; changes when `ENGINE_SEMANTICS` is monkeypatched; differs between variants of a per-variant call.
4. A config with a `ref` inside an otherwise two-module call is refused by `call_eligible`. A custom config with a per-variant call reading `{step1}` and one more module is eligible; one that names `step1_rule@R1` literally is refused (review P14).
5. `signal_peptide_protein` is not eligible although `reads_of_call` without a variant returns four step 1 modules (review P2).

Mutations (each makes a named test fail): hash ignores thresholds; hash includes unreferenced thresholds; `call_eligible` ignores `ref`.

## Task 2: call-file loader, validity, resolver

Files: new `src/cellsurface_sorting_hat/call_status.py`, `tests/cellsurface_sorting_hat/test_call_status.py`.

Produces:
- `CallSource(call, variant, config_sha256, call_hash, reads, entries, path)`.
- `load_call_source(path, cfg)`: structural checks in D6. Raises `ValueError` with the path in the message.
- `stale_reason(source, cfg, identities, states, config_sha256)` -> `""` or one reason of D10, checked in that order.
- `CallStatusResolver(workdir, cfg, identities, states, config_sha256, lineage)`:
  - loads `status/calls/*.json` (`InputError` from `cli` is not imported here; raise `CallStatusError(ValueError)` and let `cli` convert it);
  - `__call__(call, variant, taxon)` per D2, using `status.best_entry` (lineage, most specific wins);
  - `rows()` -> list of dicts for the report and `run.json` (file, call, variant, taxon, status, calibration set, counts, rates, valid, reason).
- (in `calibration/call_files.py`) `write_call_status(workdir, cfg, call, variant, config_sha256, reads, entries)`: validates entries with `_check_entries` and the D11 rules, takes the lock (D8), merges per the spec (other calibration sets stay, the same set is replaced, everything is dropped with a stderr message when `reads`, `call_hash` or `config_sha256` differ), writes atomically with `write_atomic`, reads back with `load_call_source`.

Tests (write first), all with a tiny config in `tmp_path` and the toy taxonomy of `conftest.py`:
1. A valid file loads. Each structural refusal (also: a taxon below 2; `leakage` missing from the notes; a status above `smoke` with leakage other than `none`; the same file edited so that a `smoke` entry says `estimated`): no `measure`; status stronger than the measure allows; name differs from `call`/`variant`; unknown call; ineligible call (`ref`, `other`, one module); variant given for a non-per-variant call; duplicate taxon.
2. `stale_reason`: one test per reason and one for the valid case. Identity tests change `version`, `params_hash` and `artefact_hash` separately. `module state not ok` is parametrized over `unavailable`, `error`, `not_run`, `partial` and a missing state (review P1, P17). A test where `repeat14` is `unavailable` and `repeat02` is called shows the call status does not apply (the OR is called through `repeat02` alone).
3. Resolver lineage (a `Lineage` built in the test with a strain below a species: `TOY_NODES` has none, review P15): a tested species applies to its strain; a sibling clade does not; the most specific tested taxon wins; a not-tested taxon returns `None`.
4. A stale file returns `(None, reason)`; a missing file returns `None`; a malformed file raises `CallStatusError` naming the file; an orphan file for a call that is no longer in the config raises.
5. `write_call_status`: merge rules (other set stays, same set replaced, `call_hash` change drops old entries with a message). The lock is tested deterministically (review P7): `write_atomic` is monkeypatched to try `fcntl.flock(LOCK_EX | LOCK_NB)` on the lock file from a second file descriptor and record that it fails while the writer is inside the write, and that it succeeds afterwards. A read-back failure leaves the previous valid file in place.

Mutations: ignore `reads` names; ignore identities; ignore `call_hash`; ignore `config_sha256`; ignore module states; refuse only `unavailable` and `error`; accept a status stronger than its measure; accept taxon 1; ignore the leakage cap; keep old entries when only `call_hash` changed; drop the lock.

## Task 3: engine hook

Files: `src/cellsurface_sorting_hat/engine.py`, `tests/cellsurface_sorting_hat/test_engine.py`.

Task 3a (before any engine change): add the golden test. A fixture directory under `tests/cellsurface_sorting_hat/golden/` holds a small workdir with R0, `repeat02`, `repeat14`, `pfam_adhesion`, allergen and antigen modules, a module status file for R0 and a taxon map. `test_golden_run` runs `cli.run` on it and compares the decompressed `calls.long.tsv.gz` text with a stored `calls.long.expected.tsv` and the `report.md` with the line `- version: ...` removed (review P5, P16) with a stored `report.expected.md`. Generate the expected files from the base commit (07882e0) and commit them before Task 3b. This is the fixture review F10 asked for.

Task 3b: `evaluate(cfg, protein_ids, taxa, modules, status_of, measured_call_of=None, call_status_of=None)` implementing D3 and D4.

Tests (write first, using the existing `run` helper extended with `call_status_of`):
1. Equivalence: for a multi-protein toy with several calls, `evaluate` with a hook that always returns `None` equals `evaluate` without the argument, record for record (value, status, basis).
2. A call status for `tandem_repeat_protein` applies to its records: status and basis `call:tandem_repeat_protein:taxon:<t>`; the module statuses of `repeat02` and `repeat14` are ignored for that call.
3. `not_assessable` records stay `unvalidated` with an empty basis whatever the hook returns (review F1).
4. Composite: `cell_wall_adhesion_candidate` takes the weakest over groups: the repeat leaf with a call status, the R0 leaf with a module status measured on `signal_peptide_protein` (supply `measured_call_of`, review P8). Assert the full basis string and its order. A leaf without a call status uses module logic.
5. False OR and false AND: the status comes from the false inputs, grouped by leaf; a call status on the false leaf applies.
6. `other_surface_no_mechanism`: a called mechanism with a call status contributes it; the surface contributors do not enter, as today.
7. Per-variant: a call status keyed `(signal_peptide_protein, "R0")` is looked up with variant `R0` for a record at label `R0` and is not used at label `R2`; a non-per-variant leaf is looked up with the empty variant even at label `R0`.
8. A stale hook result `(None, reason)` gives the module-logic status and the basis item prefix `call status stale (<reason>):`.
9. A taxon not tested returns `None`: module logic, no prefix.

Mutations: apply the call status to `not_assessable` records; ignore the variant; use the call status for a sibling clade (hook mutation); drop the weakest-over-groups rule; drop the stale prefix; pass the record's call instead of the leaf call to `status_for` (caught by test 4).

## Task 4: core command, report, run.json

Files: `src/cellsurface_sorting_hat/cli.py`, `outputs.py`, `tests/cellsurface_sorting_hat/test_cli.py`, `test_outputs.py`.

Change:
- `run()` builds `CallStatusResolver` after `load_modules` (needs `loaded.identities`, `cfg.sha256`, the lineage) and passes it as `call_status_of`. A `CallStatusError` becomes `InputError`.
- `RunInfo.call_status_sources` (default `[]`), written to `run.json` as `call_status_sources` (list of `{file, call, variant, valid, reason}`) and shown in the report as a "Call calibration" table only when the list is not empty (review F15). Stale files appear in it with the reason.
- The derived-status note of D5 is printed in the "Call calibration" section only.

Tests (write first):
1. End to end on a toy workdir: a valid `status/calls/tandem_repeat_protein.json` changes the status of that call in `calls.long.tsv.gz`; `run.json` lists it; the report has the table.
2. With no `status/calls/`, the golden test of Task 3a still passes unchanged (`calls.long` and `report.md` identical; `run.json` gains `call_status_sources: []`, the one change R5 does not cover).
3. A stale file (change a module version) leaves the module-logic statuses and lists the file with its reason in the report and `run.json`.
4. A file for an unknown call stops the run with exit code 2 and a message naming the file (Q6).
5. Legacy case (D14): a module file with `call=tandem_repeat_protein` notes and a call file for that call: the call file is used, and the report shows both the module row and the call row.
6. `test_the_spec_and_the_report_list_the_same_known_limits` passes unchanged.

## Task 5: `calibrate truth --call-status`

Files: `src/cellsurface_sorting_hat/calibration/cli.py`, `tests/cellsurface_sorting_hat/calibration/test_truth_call_status.py`.

Change: `--module` and `--call-status` are mutually exclusive and one is required; add `--config` (default: the packaged config). With `--call-status`:
1. `load_config(args.config)`; the call must pass `call_eligible`.
2. `run.json` next to `--calls-long` must have `config_sha256` equal to `cfg.sha256`, and `module_states` with state `ok` for every read module (any other state, or a missing key or file, is refused with the module and state named).
3. For each read module, `_check_run_identity` (existing) must pass. `reads` is built from the workdir's module records.
4. Truth reading, matching, `build_measure` and the leakage cap are unchanged (shared with the module path; factor the common part into a helper rather than copying).
5. The entry is written by `write_call_status`.
`--module` keeps its behaviour, including the refusal of a call that reads more than one module.

Note (review P15): the calibration tests' `write_run_json` helper writes neither `config_sha256` nor `module_states`; extend it in this task.

Tests (write first): the success path for `tandem_repeat_protein` on a toy truth table (sensitivity and specificity values checked against a hand count, leakage `tuned_on_truth` caps the status at `smoke`); then a second `cellsurface_sorting_hat` run uses the file. Refusals: a read module in each of `unavailable`, `partial`, `error` and `not_run` (parametrized, review P17); missing `module_states`; module identity mismatch; config hash mismatch; genus taxon; `kind: other`; composite; single-module call; both or neither of `--module` and `--call-status`.

## Task 6: documentation and rules

Files: `docs/paper/03-status-and-validation-rules.md` (section on call files: scope, validity, order, what a reader can rely on), `docs/paper/04-limits-and-open-questions.md`.

## Task 7: real data

Not code. After Tasks 1 to 6 pass, build truth tables from curated data (a separate analysis task) and run `truth --call-status` for `tandem_repeat_protein` on S288C and *C. albicans*. Report in `docs/reports/`. Use leakage `tuned_on_truth` where a truth protein tuned a detector. Do not raise any status above what the measure and the cap allow.

## Out of scope

A call file for a single-module call; for `ref` or `other` calls; a status for a composite measured as a whole; module-file migration; changes to the interval or `estimated` rules.

## Self-review against the spec

Every row of the spec's failure table has a test: not assessable (T3.3), `other`/`ref`/one module (T1.2, T2.1, T5), reads (T2.2), identity (T2.2), no run record (T2.2), `call_hash` (T1.3, T2.2), config (T2.2, T5), `module_states` (T5), no `measure` and status strength (T2.1), duplicate taxon (T2.1), file name (T2.1), variant rules (T2.1, T3.7), unknown call (T2.4, T4.4), legacy module file with `call=` notes (T3.2), concurrent writers (T2.5), leakage cap (T5).
