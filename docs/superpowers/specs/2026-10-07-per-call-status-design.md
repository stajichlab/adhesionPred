# Per-call status entries: design

*2026-10-07. **Revision 2**, after review 1 (`2026-10-07-per-call-status-design-review-1.md`, verdict: needs rework; all findings applied or marked for the owner). Not yet re-reviewed. No code has been written. Owner decision 8 (2026-10-07): build this now, with a spec, a plan and an independent review before code.*

## 1. Problem

A status says how far a measurement supports a call's output in a taxon. Today a status belongs to a **module**: `<workdir>/status/<module>.json`, one entry per tested taxon. The entry notes name the call it was measured on (`call=signal_peptide_protein`). The engine counts a module status only for that call (`measured_call_of`). A composite call takes the weakest status of the contributors that decided it. Contributors are `(module, leaf_call)` pairs. So a composite call already takes its status from measurements of its parts. A module measured on `signal_peptide_protein` counts inside `cell_wall_adhesion_candidate`.

This works when one module decides one leaf call. It fails when the measured leaf call reads several modules:

- `tandem_repeat_protein` = `repeat02` OR `repeat14`. A truth table measures the OR. No single module can carry the measurement.
- `iuis_allergen_homolog` reads `allergen_homology` and `pfam_allergen`.

`calibrate truth` refuses such a call. These calls stay `unvalidated` whatever truth exists.

## 2. Requirements

R1. A truth table that measures a call with two or more modules writes a status for that **call** and taxon.
R2. A call status is valid only for the exact modules (name, version, parameters, artefact), the exact config, the call definition and the thresholds that the call read when it was measured.
R3. A stale or invalid call status is never used. It must never raise a status.
R4. Existing module status files stay valid and unchanged. The R0 files need no migration.
R5. A call without a valid call status behaves as today, byte for byte in `calls.long.tsv.gz` and `report.md`.
R6. The report shows call statuses with the measure (positives, negatives, clusters, sensitivity, specificity, leakage).
R7. Lineage rules stay as today: a status applies to the tested taxon and its descendants, the most specific tested taxon wins, a sibling clade does not inherit.
R8. One source per call (owner question Q8): a call has either a module file path or a call file path, never both.

## 3. Design

### 3.1 Scope and file

A call file is allowed only for a call that **has an `expr`, has no `ref` node, is not `kind: other`, and reads two or more modules** (after `{step1}` expansion for the variant). Today that is `tandem_repeat_protein` and `iuis_allergen_homolog`. The writer and the loader both refuse any other call. Single-module calls keep module files (`--module`).

File: `<workdir>/status/calls/<call>[.<variant>].json`. The file name must equal the content.

```json
{
  "call": "tandem_repeat_protein",
  "variant": "",
  "config_sha256": "<sha256 of the config file the run used>",
  "call_hash": "<see 3.1.1>",
  "reads": [
    {"name": "repeat02", "version": "...", "params_hash": "...", "artefact_hash": "..."},
    {"name": "repeat14", "version": "...", "params_hash": "...", "artefact_hash": "..."}
  ],
  "entries": [ {"taxa": [4932], "status": "smoke", "source": "...", "measure": { ... }} ]
}
```

`entries` use the existing entry and `measure` checks. `reads` lists every module the call reads for the variant, sorted by name. `variant` is given only for a `per_variant` call.

Loader rules (stricter than `load_status_source`): `measure` is required in every entry; `status` must not be stronger than `status_from_measure(measure)`; the file name must match call and variant; the call must exist in the config and be eligible; a taxon in two entries is refused.

#### 3.1.1 `call_hash`

SHA-256 of the canonical JSON (sorted keys, no spaces) of an object with:

- the call's `expr`, with every threshold reference replaced by the numeric value in force (only the thresholds the call references, so an unrelated threshold change does not make the file stale);
- `per_variant`;
- for `{step1}` calls, the variant label and the module name it maps to;
- `ENGINE_SEMANTICS`, a constant in `engine.py`, bumped by hand when the evaluation rules change (Kleene tables, NA handling, `_bad_value`).

Because call files are limited to calls without `ref`, no transitive hash is needed. If the rule on `ref` is ever relaxed, the hash must follow `ref` targets.

### 3.2 Resolution (engine)

`evaluate` gets a new optional keyword argument, `call_status_of(call, variant, taxon)`. It returns `(status, basis, tested_taxon)` or `None`. Without it, behaviour is unchanged (R5).

For a record (protein, call, variant, taxon T):

1. If `res.contributors` is empty (a `not_assessable` value), the record is `unvalidated` with an empty basis, as today. **The call status never applies here** (review F1).
2. Otherwise group the contributors by `(leaf_call, variant)`. The variant of a leaf is the composite's label if the leaf call is `per_variant`, else the empty string (review F6).
3. For each group, if a valid call status exists for `(leaf_call, variant)` and an entry applies to T by lineage, the group contributes that one status, with the basis `call:<leaf_call>:taxon:<tested>`.
4. Every other group contributes through the existing module logic (`status_for`).
5. The record's status is the weakest over the groups. The basis joins the group bases with `;`.

The same grouping serves composite calls (decision 5 as implemented), false AND and OR nodes, and the `called |= res.contributors` rule of `other_*` calls.

Stale or invalid sources: if a call file exists and is not valid, the group uses the module logic and its basis gets the prefix `call status stale (<reason>):`. The reasons are fixed strings: `reads differ`, `identity differs: <module>`, `no module run record: <module>`, `call_hash differs`, `config differs`. A call file for an unknown or ineligible call is a load error (Q6).

Valid means: every `reads` name equals the set of modules the call reads now; every identity equals the running identity; `call_hash` equals the hash computed now; `config_sha256` equals the run's config.

### 3.3 Calibration command

`calibrate truth` gets `--call-status`, mutually exclusive with `--module` (today `--module` is required). It also gets `--config` (default: the packaged file), used to read the call definitions.

With `--call-status`:

- The call must be eligible (3.1). Otherwise refuse.
- `run.json` must carry `config_sha256` equal to the hash of the loaded config, and `module_states` for every module the call reads. Each state must be `ok`. A state of `unavailable`, `partial`, `error` or `not_run`, or a missing state, is refused (review F3). A read module with `unavailable` state is refused even when a module record with an identity exists.
- Each read module's identity in `run.json` must equal the module record now in the work directory (the check that exists, `_check_run_identity`, applied to each module).
- The measure code is unchanged. It works on the call's final value in `calls.long`.
- Merge: entries of other calibration sets stay, entries of the same set are replaced, and all old entries are dropped (with a message) when `reads`, `call_hash` or `config_sha256` differ.
- The update takes an `fcntl` lock on `status/calls/<file>.lock`. The write is atomic as today.
- One call, one variant, one species-level taxon per command.

`--module` keeps its behaviour. It is still refused for a call that reads more than one module.

### 3.4 Report

The core report adds "Call calibration" (call, variant, taxon, status, calibration set, positives, negatives, sensitivity, specificity) **only when call files exist**, so `report.md` is unchanged otherwise (R5). Stale files are listed with the reason. `calls.long.tsv.gz` keeps its columns. `run.json` gets `call_status_sources` (file, call, variant, valid or the reason).

### 3.5 Failure modes

| case | behaviour |
|---|---|
| value is `not_assessable` (no contributors) | call status not applied (3.2 step 1) |
| call is `kind: other`, has `ref`, or reads fewer than two modules | writer refuses; loader refuses the file |
| `reads` names differ from the modules read now | stale: `reads differ` |
| identity differs (version, params, artefact) | stale: `identity differs: <module>` |
| a read module has no run record in the work directory | stale: `no module run record: <module>` |
| `call_hash` differs (expression, threshold, variant mapping, `ENGINE_SEMANTICS`) | stale |
| `config_sha256` differs from the run's config | stale (reader) or refused (writer) |
| a read module is `unavailable`, `partial`, `error`, `not_run`, or has no `module_states` entry in `run.json` | writer refuses |
| `run.json` has no `module_states` | writer refuses |
| entry has no `measure`, or status is stronger than the measure | loader refuses |
| same taxon in two entries | loader refuses |
| file name differs from `call`/`variant` | loader refuses |
| `variant` given for a call that is not `per_variant`, or missing for one that is | loader and writer refuse |
| call name not in the config, or call no longer eligible | load error naming the file (Q6) |
| a module file (`status/<module>.json`) with `call=` notes for a call that has a call file | the call file is used; the report lists both (legacy case; Q8) |
| two `truth` runs update the same file | serialised by the lock |
| status cap by leakage | as today (`smoke`) |

## 4. Tests (written before code)

1. Valid call entry: status and basis for a contributor group; a descendant taxon inherits; a sibling clade does not; the most specific tested taxon wins.
2. OR call with two modules: the call entry applies; module entries are ignored for that call.
3. `not_assessable` records stay `unvalidated` with an empty basis, whatever call files exist.
4. Per-variant lookup: a `signal_peptide_protein.R0` style file does not apply to R1.
5. Composite: `cell_wall_adhesion_candidate` takes the weakest over groups; a leaf with a call entry contributes one status; a leaf without one uses module logic. A false OR and a false AND take their statuses from the false inputs, grouped by leaf.
6. Each stale case in 3.5 gives the module logic with the reason prefix in the basis.
7. Loader refusals: no `measure`; status stronger than measure; wrong file name; ineligible call; duplicate taxon; variant rules.
8. Backward compatibility (golden): with no `status/calls/` directory, `calls.long.tsv.gz` and `report.md` are byte-identical to the stored output. The R0 module files are committed as a fixture under `tests/` (the `_workdir/` files are git-ignored).
9. `calibrate truth --call-status` on a fixture for `tandem_repeat_protein`: file written; a second engine run uses it; a changed config, changed module identity, and changed threshold each make it stale.
10. Writer refusals: `unavailable` read module; missing `module_states`; wrong run identity; config hash mismatch; genus taxon; `kind: other`; composite call; single-module call.
11. Merge: only `call_hash` changed, so old entries are dropped; same calibration set is replaced; other sets stay.
12. Concurrency: two writers on one file keep both entries (lock test).
13. Report: the call table appears only when call files exist.

Mutation checks (each must make a named test fail): apply the status to `not_assessable` records; ignore the variant in the lookup; ignore `reads` names; ignore `reads` identities (separate mutations); ignore `call_hash`; ignore `config_sha256`; accept `kind: other` or `ref` calls; skip the `module_states` check; let the loader accept a status stronger than its measure; keep old entries when only `call_hash` changed; let a call entry override a more specific module entry; use the call entry for a sibling clade.

## 5. Files that change

`src/cellsurface_sorting_hat/status.py` (call-file loader, validity, resolution), `engine.py` (`call_status_of` hook, grouping, `ENGINE_SEMANTICS`, `reads_of_call` moved here from `calibration/cli.py` so the engine and the writer share it), `cli.py` (load call sources, pass the hook, `run.json` field), `outputs.py` (`RunInfo`, `render_report` calibration block, `write_run_json`), `calibration/cli.py` (`--call-status`, `--config`, `--module` no longer `required=True`, lock), `calibration/measure.py` (the writer for call files), tests under `tests/cellsurface_sorting_hat/` plus the committed fixture, and `docs/paper/03-status-and-validation-rules.md` (rules; and the decision 3 wording, Q7).

Out of scope: a status for a composite call measured as a whole; a call file for single-module calls; migration of module files; any change to the interval or `estimated` rules; new truth tables.

## 6. Open items for the owner

Q6. A call file for a call that no longer exists (or is no longer eligible) stops the run, as every other bad status file does. The alternative is to skip it with a warning. Proposal: stop, with a message naming the file. A stale file for an existing call never stops the run.
Q7. `docs/paper/03` decision 3 says a call with no measurement of its own reports `unvalidated`. The code already gives a composite call the weakest status of its measured leaves. Proposal: fix the paper text to match the code.
Q8. One source per call: call files only for calls with two or more modules, module files for the rest. The R0 files stay as they are. Proposal: yes.
