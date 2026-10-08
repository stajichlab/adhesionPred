# Per-call status entries: implementation plan

*2026-10-08. Plan for `docs/superpowers/specs/2026-10-07-per-call-status-design.md` (revision 2; owner decisions Q6 to Q8 recorded). Branch `per-call-status`. Test first in every task. One commit per task. Nothing is pushed until the owner says so. Not yet reviewed.*

Python 3.12. Run tests with `PYTHONPATH=$PWD/src /usr/bin/python3.12 -m pytest tests/cellsurface_sorting_hat -q -p no:cacheprovider`. The baseline is 683 passed, 7 skipped. The pre-commit hook runs ruff 0.3.5 and ruff-format and aborts a commit that it reformats: re-add and commit again.

## Decisions taken in this plan (the spec left them open or implicit)

D1. **Module layout.** A new module `src/cellsurface_sorting_hat/call_status.py` holds the loader, the validity check, the resolver and the writer. It imports `status` (`ENTRY_KEYS`, `_check_measure`, `weakest`, `resolve_entry`) and `calibration.measure` (`status_from_measure`, `_check_entries`, `make_entry`). `calibration.measure` never imports `call_status`, so there is no cycle. `engine.py` gets the call helpers (`reads_of_call`, `call_eligible`, `call_hash`, `ENGINE_SEMANTICS`) because they read the config. `call_status.py` imports `engine`. `engine.py` does not import `call_status`.
D2. **Hook shape.** `evaluate(..., call_status_of=None)`. The hook is `call_status_of(call, variant, taxon)`. It returns `None` (no file for this call, or the taxon is not tested: module logic, no prefix), `(status, basis)` (applies), or `(None, reason)` (a file exists and is stale: module logic with the prefix `call status stale (<reason>):` on that leaf's basis items).
D3. **Byte identity (R5).** When no hook is given, or when it returns `None` for every leaf of a record, `evaluate` takes the old code path unchanged. The new path runs only when at least one leaf has a call status or a stale reason.
D4. **Grouping.** Contributors are `(module, leaf_call)`. The leaf variant is the record's label if the leaf call is `per_variant`, else the empty string. A leaf with an applying call status contributes one status and one basis item (`call:<leaf>:taxon:<tested>`), placed where the leaf's first item would stand in the old sorted order. Every other pair goes through `status_for` as today. The record status is the weakest over all items.
D5. **Static report line.** `_known_limits` gets one line: composite call statuses are derived (owner decision Q7). That changes `report.md` for every run, so the golden report text in `tests/` is updated in Task 4. `calls.long.tsv.gz` stays byte-identical.
D6. **Loader strictness (review F8).** The call-file loader requires `measure` in every entry, checks `status` against `status_from_measure`, rejects an unknown call, an ineligible call, a name that does not match, a variant error, and a taxon in two entries. `load_status_source` (module files) is unchanged.
D7. **Eager load, hard errors (Q6).** `CallStatusResolver.__init__` loads every `status/calls/*.json` and raises `InputError` naming the file for any structural problem (unknown or ineligible call, bad file, wrong name). A well-formed but stale file never raises.
D8. **Lock (review F13).** `write_call_status` takes `fcntl.flock` on `status/calls/<file>.lock` around read, merge and write.
D9. **Migration.** None. No change to module files or `load_status_source`.

## Task 1: engine helpers

Files: `src/cellsurface_sorting_hat/engine.py`, `tests/cellsurface_sorting_hat/test_call_helpers.py`.

Produces:
- `ENGINE_SEMANTICS = "1"` (bumped by hand when the evaluation rules change).
- `reads_of_call(cfg, call, variant="")`: `sorted(modules_of_call(cfg, call))`; with a variant, a step 1 module of another variant is removed (moved from `calibration/cli._reads_of_call`; that function becomes a thin alias).
- `call_eligible(cfg, call)` -> `(bool, reason)`. Eligible: the call exists, has `expr`, is not `kind: other`, its expression has no `ref` node (walk), and `len(reads_of_call) >= 2` for at least the default variant. Reasons: `unknown call`, `kind other`, `contains ref`, `reads fewer than two modules`.
- `call_hash(cfg, call, variant="")`: SHA-256 of canonical JSON (sorted keys, separators `(",", ":")`) of `{name, per_variant, variant, variant_module, expr, semantics}`. In `expr`, every `test` value `"$name"` is replaced by the number in `cfg.thresholds`, and `{step1}` in a module name by the variant module. Only thresholds the call references appear.

Tests (write first):
1. `reads_of_call` for `tandem_repeat_protein` is `["repeat02", "repeat14"]`; for `iuis_allergen_homolog` it is `["allergen_homology", "pfam_allergen"]`; for `signal_peptide_protein` with variant `R0` it is `["step1_rule@R0"]`.
2. `call_eligible`: true for the two calls above; false with the exact reason for `cell_wall_adhesion_candidate` (`contains ref`), `other_not_surface` (`kind other`), `signal_peptide_protein` (`reads fewer than two modules`), and an unknown name.
3. `call_hash`: stable; changes when a referenced threshold changes (copy the config to `tmp_path`, edit, reload); does not change when an unreferenced threshold changes; changes when `ENGINE_SEMANTICS` is monkeypatched; differs between variants of a per-variant call.
4. A config with a `ref` inside an otherwise two-module call is refused by `call_eligible`.

Mutations (each makes a named test fail): hash ignores thresholds; hash includes unreferenced thresholds; `call_eligible` ignores `ref`.

## Task 2: call-file loader, validity, resolver

Files: new `src/cellsurface_sorting_hat/call_status.py`, `tests/cellsurface_sorting_hat/test_call_status.py`.

Produces:
- `CallSource(call, variant, config_sha256, call_hash, reads, entries, path)`.
- `load_call_source(path, cfg)`: structural checks in D6. Raises `ValueError` with the path in the message.
- `stale_reason(source, cfg, identities, config_sha256)` -> `""` or one of `config differs`, `reads differ`, `no module run record: <m>`, `identity differs: <m>`, `call_hash differs`. Checked in that order.
- `CallStatusResolver(workdir, cfg, identities, config_sha256, lineage)`:
  - loads `status/calls/*.json` (`InputError` from `cli` is not imported here; raise `CallStatusError(ValueError)` and let `cli` convert it);
  - `__call__(call, variant, taxon)` per D2, using `resolve_entry` semantics (lineage, most specific wins) on the entries;
  - `rows()` -> list of dicts for the report and `run.json` (file, call, variant, taxon, status, calibration set, counts, rates, valid, reason).
- `write_call_status(workdir, cfg, call, variant, config_sha256, reads, entries)`: validates entries with `_check_entries`, takes the lock (D8), merges per the spec (other calibration sets stay, the same set is replaced, everything is dropped with a stderr message when `reads`, `call_hash` or `config_sha256` differ), writes atomically with `write_atomic`, reads back with `load_call_source`.

Tests (write first), all with a tiny config in `tmp_path` and the toy taxonomy of `conftest.py`:
1. A valid file loads. Each structural refusal: no `measure`; status stronger than the measure allows; name differs from `call`/`variant`; unknown call; ineligible call (`ref`, `other`, one module); variant given for a non-per-variant call; duplicate taxon.
2. `stale_reason`: one test per reason and one for the valid case. Identity tests change `version`, `params_hash` and `artefact_hash` separately.
3. Resolver lineage: a tested species applies to its strain; a sibling clade does not; the most specific tested taxon wins; a not-tested taxon returns `None`.
4. A stale file returns `(None, reason)`; a missing file returns `None`; a malformed file raises `CallStatusError` naming the file; an orphan file for a call that is no longer in the config raises.
5. `write_call_status`: merge rules (other set stays, same set replaced, `call_hash` change drops old entries with a message); two processes writing the same file keep both entries (use `multiprocessing` with two different calibration sets).

Mutations: ignore `reads` names; ignore identities; ignore `call_hash`; ignore `config_sha256`; accept a status stronger than its measure; keep old entries when only `call_hash` changed; drop the lock.

## Task 3: engine hook

Files: `src/cellsurface_sorting_hat/engine.py`, `tests/cellsurface_sorting_hat/test_engine.py`.

Change: `evaluate(cfg, protein_ids, taxa, modules, status_of, measured_call_of=None, call_status_of=None)` implementing D3 and D4.

Tests (write first, using the existing `run` helper extended with `call_status_of`):
1. Equivalence: for a multi-protein toy with several calls, `evaluate` with a hook that always returns `None` equals `evaluate` without the argument, record for record (value, status, basis).
2. A call status for `tandem_repeat_protein` applies to its records: status and basis `call:tandem_repeat_protein:taxon:<t>`; the module statuses of `repeat02` and `repeat14` are ignored for that call.
3. `not_assessable` records stay `unvalidated` with an empty basis whatever the hook returns (review F1).
4. Composite: `cell_wall_adhesion_candidate` takes the weakest over groups: the repeat leaf with a call status, the R0 leaf with a module status. A leaf without a call status uses module logic.
5. False OR and false AND: the status comes from the false inputs, grouped by leaf; a call status on the false leaf applies.
6. `other_surface_no_mechanism`: a called mechanism with a call status contributes it; the surface contributors do not enter, as today.
7. Per-variant: a call status keyed `(signal_peptide_protein, "R0")` is looked up with variant `R0` for a record at label `R0` and is not used at label `R2`; a non-per-variant leaf is looked up with the empty variant even at label `R0`.
8. A stale hook result `(None, reason)` gives the module-logic status and the basis item prefix `call status stale (<reason>):`.
9. A taxon not tested returns `None`: module logic, no prefix.

Mutations: apply the call status to `not_assessable` records; ignore the variant; use the call status for a sibling clade (hook mutation); drop the weakest-over-groups rule; drop the stale prefix.

## Task 4: core command, report, run.json

Files: `src/cellsurface_sorting_hat/cli.py`, `outputs.py`, `tests/cellsurface_sorting_hat/test_cli.py`, `test_outputs.py`.

Change:
- `run()` builds `CallStatusResolver` after `load_modules` (needs `loaded.identities`, `cfg.sha256`, the lineage) and passes it as `call_status_of`. A `CallStatusError` becomes `InputError`.
- `RunInfo.call_status_sources` (default `[]`), written to `run.json` as `call_status_sources` (list of `{file, call, variant, valid, reason}`) and shown in the report as a "Call calibration" table only when the list is not empty (review F15). Stale files appear in it with the reason.
- `_known_limits` gets the derived-status line (D5).

Tests (write first):
1. End to end on a toy workdir: a valid `status/calls/tandem_repeat_protein.json` changes the status of that call in `calls.long.tsv.gz`; `run.json` lists it; the report has the table.
2. With no `status/calls/`, `calls.long.tsv.gz` is byte-identical to a run of the same fixture on the commit before this task (store the expected bytes' SHA-256 in the test) and `report.md` differs only by the new limits line.
3. A stale file (change a module version) leaves the module-logic statuses and lists the file with its reason in the report and `run.json`.
4. A file for an unknown call stops the run with exit code 2 and a message naming the file (Q6).

## Task 5: `calibrate truth --call-status`

Files: `src/cellsurface_sorting_hat/calibration/cli.py`, `tests/cellsurface_sorting_hat/calibration/test_truth_call_status.py`.

Change: `--module` and `--call-status` are mutually exclusive and one is required; add `--config` (default: the packaged config). With `--call-status`:
1. `load_config(args.config)`; the call must pass `call_eligible`.
2. `run.json` next to `--calls-long` must have `config_sha256` equal to `cfg.sha256`, and `module_states` with state `ok` for every read module (any other state, or a missing key or file, is refused with the module and state named).
3. For each read module, `_check_run_identity` (existing) must pass. `reads` is built from the workdir's module records.
4. Truth reading, matching, `build_measure` and the leakage cap are unchanged (shared with the module path; factor the common part into a helper rather than copying).
5. The entry is written by `write_call_status`.
`--module` keeps its behaviour, including the refusal of a call that reads more than one module.

Tests (write first): the success path for `tandem_repeat_protein` on a toy truth table (sensitivity and specificity values checked against a hand count, leakage `tuned_on_truth` caps the status at `smoke`); then a second `cellsurface_sorting_hat` run uses the file. Refusals: `unavailable` read module; missing `module_states`; module identity mismatch; config hash mismatch; genus taxon; `kind: other`; composite; single-module call; both or neither of `--module` and `--call-status`.

## Task 6: documentation and rules

Files: `docs/paper/03-status-and-validation-rules.md` (section on call files: scope, validity, order, what a reader can rely on), `docs/paper/04-limits-and-open-questions.md`.

## Task 7: real data

Not code. After Tasks 1 to 6 pass, build truth tables from curated data (a separate analysis task) and run `truth --call-status` for `tandem_repeat_protein` on S288C and *C. albicans*. Report in `docs/reports/`. Use leakage `tuned_on_truth` where a truth protein tuned a detector. Do not raise any status above what the measure and the cap allow.

## Out of scope

A call file for a single-module call; for `ref` or `other` calls; a status for a composite measured as a whole; module-file migration; changes to the interval or `estimated` rules.

## Self-review against the spec

Every row of the spec's failure table has a test: not assessable (T3.3), `other`/`ref`/one module (T1.2, T2.1, T5), reads (T2.2), identity (T2.2), no run record (T2.2), `call_hash` (T1.3, T2.2), config (T2.2, T5), `module_states` (T5), no `measure` and status strength (T2.1), duplicate taxon (T2.1), file name (T2.1), variant rules (T2.1, T3.7), unknown call (T2.4, T4.4), legacy module file with `call=` notes (T3.2), concurrent writers (T2.5), leakage cap (T5).
