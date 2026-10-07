# Per-call status entries: design

*2026-10-07. Draft for independent review. No code has been written. Owner decision 8 (2026-10-07): build this now, with a spec, a plan and an independent review before code.*

## 1. Problem

A status says how far a measurement supports a module's output in a taxon. Today a status belongs to a **module**: `_workdir/.../status/<module>.json`, one entry per tested taxon. The entry notes name the call it was measured on (`call=signal_peptide_protein`). The engine counts a module status only for that call (`measured_call_of`). A composite call takes the weakest status of the modules that decided it (decision 5).

This works when one module decides one call. It fails when the measured call reads several modules:

- `tandem_repeat_protein` = `repeat02` OR `repeat14`. A truth table measures the OR. The measurement cannot be attributed to one module.
- `iuis_allergen_homolog` reads `allergen_homology` and `pfam_allergen`.

`calibrate truth` refuses such a call (decision 3). These calls stay `unvalidated` for ever, whatever truth exists.

## 2. Requirements

R1. A truth table that measures a call writes a status for that **call** and taxon, whatever number of modules the call reads.
R2. A call status is valid only for the exact modules, module versions, parameters and artefacts that the call read when it was measured, and for the call definition and thresholds in force.
R3. A stale call status is never used. It must not raise a status silently.
R4. Existing module status files stay valid and unchanged. The R0 files need no migration.
R5. A call without a valid call status behaves as today.
R6. The report shows call statuses next to module statuses, with the measure (positives, negatives, clusters, sensitivity, specificity, leakage).
R7. Lineage rules stay as today: a status applies to the tested taxon and its descendants, the most specific tested taxon wins, a sibling clade does not inherit.

## 3. Design

### 3.1 File

`<workdir>/status/calls/<call>[.<variant>].json`, for example `tandem_repeat_protein.json` or `signal_peptide_protein.R0.json`.

```json
{
  "call": "tandem_repeat_protein",
  "variant": "",
  "call_hash": "<sha256 of the canonical JSON of the call's expression and the thresholds it references>",
  "reads": [
    {"name": "repeat02", "version": "...", "params_hash": "...", "artefact_hash": "..."},
    {"name": "repeat14", "version": "...", "params_hash": "...", "artefact_hash": "..."}
  ],
  "entries": [ {"taxa": [4932], "status": "smoke", "source": "...", "measure": { ... same as today ... }} ]
}
```

`entries` use the loader and the `measure` check that exist (`ENTRY_KEYS`, `MEASURE_KEYS`). `reads` lists every module in the call's expression for that variant, including a module in an OR branch, sorted by name.

### 3.2 Resolution (engine)

For a record (protein, call, variant) with taxon T:

1. If a call source exists for (call, variant) and is **valid**, find its entry by lineage (as `resolve_entry`). If one applies, the record's status is that entry's status and the basis is `call:<call>:taxon:<tested>`.
2. If the file exists but is not valid, record the reason in the basis (`call status stale: reads differ`) and use the existing module logic.
3. If no entry applies (taxon not tested) or no file exists, use the existing module logic (R5).

Valid means: the `reads` set equals the set of modules the call reads now (names), every identity (version, params hash, artefact hash) equals the running identity, and `call_hash` equals the current hash. Any difference is stale (R3).

Precedence: a valid call entry replaces the weakest-of-modules value for that call only.

### 3.3 Composite calls

Decision 5 stays: a composite call (for example `cell_wall_adhesion_candidate`) takes the weakest status of the modules that decided it. One change: where a deciding module belongs to a leaf call that has a valid call status for the taxon, that leaf's call status stands in for the status of the modules of that leaf. This is a question for review (section 6, Q1).

### 3.4 Calibration command

`calibrate truth` gets a mode for calls:

- `--module M` as today: writes `status/M.json`. Refused for a call that reads more than one module (unchanged).
- `--call-status` (no `--module`): writes `status/calls/<call>[.<variant>].json`. Allowed for any call. `reads` come from `run.json` `module_identities` for every module the call reads. Each identity must equal the module record now in the work directory (the check that exists, `_check_run_identity`, applied to each module).
- The measure code is unchanged: it already works on the call's final value in `calls.long`. Leakage and the `smoke` cap work as today. The merge works as `_merge_entries` does now: entries of other calibration sets stay, entries of the same set are replaced, and all old entries are dropped (with a message) when the identity changed. For a call file, "identity changed" means the `reads` identities or the `call_hash` differ. The loader refuses a taxon that two entries list (as `load_status_source` does).
- One call, one variant, one species-level taxon per command (as today).

### 3.5 Report

The core report adds a table "Call calibration" (call, variant, taxon, status, calibration set, positives, negatives, sensitivity, specificity) from valid call sources. Stale sources are listed with the reason. `calls.long.tsv.gz` keeps its columns. `status_basis` carries the new basis strings.

### 3.6 Failure modes

| case | behaviour |
|---|---|
| `reads` names a module the call no longer reads, or misses one it reads | stale, reason in basis |
| identity differs (version, params, artefact) | stale |
| `call_hash` differs (expression or threshold changed) | stale |
| same taxon in two entries | file refused at load |
| call name or variant unknown to the config | file refused at load |
| a read module is `unavailable` in the run | the call status is not written (the command refuses: a read module has no identity) |
| file for a call that has `per_variant` and no `variant` | refused |
| status cap by leakage | as today |

## 4. Tests (written before code)

1. Valid call entry gives the call's status and basis; a descendant taxon inherits; a sibling clade does not; the most specific wins.
2. OR call with two modules: call entry applies, module entries ignored for that call.
3. Each stale case in 3.6 gives `unvalidated` or the module logic, with the reason in the basis.
4. Backward compatibility: the stored R0 module files give the same `calls.long` as before (golden test).
5. `calibrate truth --call-status` on a fixture for `tandem_repeat_protein`: file written, then used by the engine in a second run.
6. Refusal tests: unavailable read module, wrong run identity, genus taxon, duplicate taxon.
7. Report test: the call table appears.

Mutation checks (each must make a named test fail): ignore `reads`; ignore `call_hash`; use the call entry for a sibling clade; let a stale file raise the status; drop the identity check in `--call-status`.

## 5. Files that change

`src/cellsurface_sorting_hat/status.py` (loader, validity, resolution), `engine.py` (call status hook), `cli.py` (load call sources, pass the hook, report), `calibration/cli.py` and `calibration/measure.py` (the `--call-status` mode), tests under `tests/cellsurface_sorting_hat/`, and `docs/paper/03-status-and-validation-rules.md` (rules section).

Out of scope: a status for a composite call measured as a whole; migration of module files; any change to the interval or `estimated` rules; new truth tables.

## 6. Questions for the reviewer

Q1. Section 3.3: should a leaf call status stand in for its modules inside a composite call? Without it, a composite call stays `unvalidated` even when its repeat leaf is `smoke`. With it, the composite takes a status from a measurement of a part. Decision 5 (weakest of the deciding modules) and decision 3 (status from a measurement of that call) pull in different directions.
Q2. `iuis_allergen_homolog` reads `pfam_allergen`, which is `unavailable` until a family is active. Is "refuse to write" right, or should the call be measurable on the modules that exist?
Q3. Is `call_hash` over the expression and the thresholds enough? Module parameters are covered by `params_hash`.
Q4. Should module files and call files be allowed together for one call? Section 3.2 says the call file wins. Is a conflict ever possible?
Q5. Anything in the failure table that is missing.
