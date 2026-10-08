"""Write calls.long.tsv.gz, calls.wide.tsv.gz, evidence.tsv.gz, proteins.tsv.gz, report.md, run.json."""

import csv
import gzip
import io
import json
from collections import Counter
from dataclasses import asdict, dataclass, field

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.logic import CALLED, NOT_ASSESSABLE, NOT_CALLED

LONG_COLUMNS = ["protein", "call", "variant", "value", "status", "status_basis", "other_basis"]
MAX_LISTED_INVALID = 20


def _known_limits(thresholds):
    percent = thresholds.get("antigen_percentile_max", "?")
    if isinstance(percent, float):
        percent = f"{percent:g}"
    return [
        "Research use only. This is not a regulatory allergenicity assessment and not a diagnostic "
        "result. No row is supported by an IgE, antibody or T-cell measurement.",
        "`signal_peptide_protein` means that SignalP calls a signal peptide (rule R0) and nothing more. "
        "Such proteins are secreted, wall-bound or GPI-anchored. They are not shown to be exposed at the "
        "cell surface, and glycosylation is not assessed. Plasma membrane mucins can be missed.",
        "GPI-anchored and secreted enzymes, and non-adhesive structural wall proteins, get no finer "
        "label in version 1. There is no `cell_wall_protein` call.",
        "`tandem_repeat_protein` and `wall_family_domain` are evidence. Repeat proteins include "
        "intracellular ones (ubiquitin, calmodulin, ankyrin proteins). A domain of a family that is "
        "linked to adhesion or wall function in at least one species (CFEM, Bys1, Als) is "
        "not shown to mediate adhesion here. The adhesion call needs a signal peptide. "
        "`hydrophobin_domain` (hydrophobin Pfam models, category 2c) and `hsba_domain` (HsbA, PF12296, "
        "which may not be a hydrophobin) are separate evidence calls. They do not set `wall_family_domain`.",
        f"`cocci_specificity_rank_top15` is the top {percent}% of a fixed Coccidioides immitis ranking "
        "(similarity to IEDB antigens, prevalence, absence of orthologs in confounder fungi). It is not "
        "epitope prediction. The cut was set after the four anchors were seen. The ranking prints NOT "
        "CALIBRATED. `serodiagnostic_marker_candidate` adds a signal peptide. Peptide level only; glycan "
        "epitopes are not assessed.",
        "`iuis_allergen_similarity` and `iuis_allergen_homolog` are sequence similarity to allergens in "
        "the WHO/IUIS fungal set (IgE binding in patients). They suggest possible IgE cross-reactivity "
        "at most. A protein with no hit is not thereby a non-allergen. WHO/IUIS lists no Coccidioides "
        "allergen.",
        "Cell wall integrity signaling, septation, polarized growth, polysaccharide chemistry, "
        "non-protein adhesins, moonlighting proteins and biofilm are not categories.",
        "The taxon you give is recorded as given. It is not checked against the sequences.",
        "Leakage: overlap between the Phase C positives and the SignalP 6 training data was not "
        "measured. The status of rule R0 does not account for it.",
    ]


@dataclass
class RunInfo:
    version: str
    taxa: dict  # taxon ID -> number of proteins
    taxonomy_sha256: str
    config_sha256: str
    default_gate: str
    thresholds: dict
    n_proteins: int
    n_invalid: int
    invalid: list  # [id, note] pairs
    n_trailing_stop: int
    module_states: dict  # module name -> run state
    module_notes: dict = field(default_factory=dict)  # module name -> text
    state_counts: dict = field(default_factory=dict)  # module -> {protein state: count}
    absent_modules: list = field(default_factory=list)  # needed by the rules, not found
    unavailable_variants: list = field(default_factory=list)
    inconsistent: dict = field(default_factory=dict)  # module -> groups of identical sequences
    map_ids_not_in_fasta: int = 0
    unmatched_module_ids: dict = field(default_factory=dict)  # module -> rows with no FASTA ID
    calibration: list = field(default_factory=list)  # rows from cli.calibration_rows
    call_status_sources: list = field(default_factory=list)  # rows from CallStatusResolver.rows()
    # one {name, version, params_hash, artefact_hash} per module that was loaded
    module_identities: list = field(default_factory=list)


def _gz(text):
    return gzip.compress(text.encode(), mtime=0)


def _tsv(header, rows):
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return _gz(buf.getvalue())


def column_name(call, variant):
    return f"{call}[{variant}]" if variant else call


def write_long(path, records):
    rows = [
        [r.protein, r.call, r.variant, r.value, r.status, r.status_basis, r.other_basis]
        for r in records
    ]
    write_atomic(path, _tsv(LONG_COLUMNS, rows))


def write_wide(path, records, protein_ids):
    columns, table = [], {pid: {} for pid in protein_ids}
    for r in records:
        name = column_name(r.call, r.variant)
        if name not in columns:
            columns.append(name)
        row = table[r.protein]
        row[name] = r.value
        row[name + "_status"] = r.status
        if r.call.startswith("other_"):
            row[name + "_basis"] = r.other_basis
    header = ["protein"]
    for name in columns:
        header += [name, name + "_status"]
        if name.startswith("other_"):
            header.append(name + "_basis")
    rows = [[pid] + [table[pid].get(h, "") for h in header[1:]] for pid in protein_ids]
    write_atomic(path, _tsv(header, rows))


def write_evidence(path, rows):
    write_atomic(path, _tsv(["protein", "module", "field", "value"], rows))


def write_proteins(path, proteins, taxa):
    header = ["id", "sha256", "taxon", "state", "note", "trailing_stop", "ambiguous_fraction"]
    rows = [
        [
            p.id,
            p.sha256,
            taxa[p.id],
            p.state,
            p.note,
            int(p.trailing_stop),
            f"{p.ambiguous_fraction:.4f}",
        ]
        for p in proteins
    ]
    write_atomic(path, _tsv(header, rows))


def write_run_json(path, info):
    write_atomic(path, (json.dumps(asdict(info), indent=2, sort_keys=True) + "\n").encode())


def _warnings(info):
    out = []
    if info.default_gate in info.unavailable_variants:
        out.append(
            f"The default gate {info.default_gate} is not available. Gated calls are not written for it."
        )
    for module, state in sorted(info.module_states.items()):
        if state in ("error", "partial"):
            note = info.module_notes.get(module, "")
            out.append(f"Module {module} is in state {state}. {note}".strip())
    for module, n in sorted(info.inconsistent.items()):
        out.append(
            f"Module {module} gives different rows to {n} group(s) of identical sequences. "
            "Check the module."
        )
    if info.absent_modules:
        out.append(
            "Modules the rules need and that were not found: " + ", ".join(info.absent_modules)
        )
    for module, n in sorted(info.unmatched_module_ids.items()):
        out.append(
            f"Module {module} has {n} row(s) whose ID is not in the FASTA. They are ignored."
        )
    if info.map_ids_not_in_fasta:
        out.append(f"{info.map_ids_not_in_fasta} ID(s) in --taxon-map are not in the FASTA.")
    return out


def render_report(info, records):
    counts = Counter((r.call, r.variant, r.value) for r in records)
    keys = sorted({(r.call, r.variant) for r in records})
    lines = [
        "# cellsurface_sorting_hat report",
        "",
        f"- version: {info.version}",
        f"- proteins: {info.n_proteins} ({info.n_invalid} invalid, excluded from all modules; "
        f"{info.n_trailing_stop} had a trailing `*`, removed)",
        f"- taxa (taxon ID: proteins): {', '.join(f'{t}: {n}' for t, n in sorted(info.taxa.items()))}",
        f"- taxonomy file sha256: {info.taxonomy_sha256}",
        f"- categories.yaml sha256: {info.config_sha256}",
        f"- default gate: {info.default_gate}",
        "- thresholds: " + ", ".join(f"{k} = {v}" for k, v in sorted(info.thresholds.items())),
        "",
    ]
    warnings = _warnings(info)
    if warnings:
        lines += ["## Warnings", ""] + [f"- **WARNING: {w}**" for w in warnings] + [""]
    if info.invalid:
        shown = info.invalid[:MAX_LISTED_INVALID]
        lines += ["## Invalid proteins", ""] + [f"- {i}: {note}" for i, note in shown]
        if len(info.invalid) > len(shown):
            lines.append(f"- ... and {len(info.invalid) - len(shown)} more (see proteins.tsv.gz)")
        lines.append("")
    lines += ["## Module run states", "", "| module | run state |", "|---|---|"]
    lines += [f"| {m} | {s} |" for m, s in sorted(info.module_states.items())]
    if info.unavailable_variants:
        lines += ["", "Unavailable step 1 variants: " + ", ".join(info.unavailable_variants)]
    if any(info.state_counts.values()):
        lines += [
            "",
            "## Proteins that are not `ok` in a module",
            "",
            "| module | state | proteins |",
            "|---|---|---|",
        ]
        for m, states in sorted(info.state_counts.items()):
            lines += [f"| {m} | {s} | {n} |" for s, n in sorted(states.items())]
    if info.calibration:
        lines += [
            "",
            "## Module calibration",
            "",
            "| module | taxon | status | measured on call | calibration set | positives "
            "| negatives | sensitivity [95% CI] | specificity [95% CI] |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for c in info.calibration:
            n_pos = c["n_pos"] if c["n_pos"] != "" else "-"
            n_neg = c["n_neg"] if c["n_neg"] != "" else "-"
            lines.append(
                f"| {c['module']} | {c['taxon']} | {c['status']} | {c.get('measured_call') or '-'} "
                f"| {c['calibration_set'] or '-'} | {n_pos} | {n_neg} "
                f"| {c['sensitivity']} | {c['specificity']} |"
            )
    if info.call_status_sources:
        lines += [
            "",
            "## Call calibration",
            "",
            "Calls that read several modules, measured as a whole (`status/calls/`). A valid file sets "
            "the status of its call. Statuses of composite calls (for example "
            "`cell_wall_adhesion_candidate`) are derived: the weakest status of the measured leaf "
            "calls that decided them. A stale file is listed with its reason and is not used.",
            "",
            "| call | variant | taxon | status | calibration set | positives | negatives "
            "| sensitivity [95% CI] | specificity [95% CI] | validity |",
            "|---|---|---|---|---|---|---|---|---|---|",
        ]
        for c in info.call_status_sources:
            if c["valid"]:
                lines.append(
                    f"| {c['call']} | {c['variant'] or '-'} | {c['taxon']} | {c['status']} "
                    f"| {c['calibration_set'] or '-'} | {c['n_pos'] if c['n_pos'] != '' else '-'} "
                    f"| {c['n_neg'] if c['n_neg'] != '' else '-'} | {c['sensitivity']} "
                    f"| {c['specificity']} | valid |"
                )
            else:
                lines.append(
                    f"| {c['call']} | {c['variant'] or '-'} | - | - | - | - | - | - | - "
                    f"| stale: {c['reason']} |"
                )
    lines += [
        "",
        "## Calls",
        "",
        "| call | variant | called | not_called | not_assessable |",
        "|---|---|---|---|---|",
    ]
    for call, variant in keys:
        c = [counts[(call, variant, v)] for v in (CALLED, NOT_CALLED, NOT_ASSESSABLE)]
        lines.append(f"| {call} | {variant or '-'} | {c[0]} | {c[1]} | {c[2]} |")
    basis = Counter(r.other_basis for r in records if r.call.startswith("other_") and r.other_basis)
    if basis:
        lines += ["", "## `other_basis` (categories left out because they were not assessable)", ""]
        lines += [f"- {b}: {n}" for b, n in sorted(basis.items())]
    lines += ["", "## Known limits", ""] + [
        f"{i}. {t}" for i, t in enumerate(_known_limits(info.thresholds), 1)
    ]
    return "\n".join(lines) + "\n"


def write_atomic_text(path, text):
    write_atomic(path, text.encode())


def write_report(path, info, records):
    write_atomic_text(path, render_report(info, records))
