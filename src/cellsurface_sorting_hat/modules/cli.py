"""Command ``cellsurface_sorting_hat_module``: turn one tool's output into a module table."""

import argparse
import json
import sys
from pathlib import Path

from cellsurface_sorting_hat.cli import InputError, RunError, assign_taxa, read_taxon_map
from cellsurface_sorting_hat.fasta import FastaError, read_fasta
from cellsurface_sorting_hat.modules import allergen, lookups, pfam, repeats, signalp
from cellsurface_sorting_hat.modules.base import ModuleSpec, invalid_row, write_module

APPLICABLE_RS = {
    246410
}  # C. immitis RS (taxon_id in analysis/step1_compare/phasec/tc_taxon_clades.tsv)


# blanks are real: the spherule table has rows without a padj or a day 8 value
EXPRESSION_NUMBERS = {
    "log2fc_48h": "finite_or_blank",
    "padj_48h": "nonneg_or_blank",
    "log2fc_8d": "finite_or_blank",
}


def _common(p):
    p.add_argument("--fasta", required=True)
    p.add_argument("--workdir", required=True)


def _taxon_id(text):
    value = int(text)
    if value < 2:
        raise argparse.ArgumentTypeError(f"taxon ID {value} is not a taxon of an organism")
    return value


def _blast_cap(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be an integer >= 1")
    return value


def _taxa_args(p):
    p.add_argument("--taxon", type=int)
    p.add_argument("--taxon-map")
    p.add_argument("--applicable-taxa", type=_taxon_id, nargs="+", default=sorted(APPLICABLE_RS))


def build_parser():
    ap = argparse.ArgumentParser(prog="cellsurface_sorting_hat_module")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("signalp", help="SignalP 6 prediction_results.txt -> step1_rule@R0")
    _common(p)
    p.add_argument("--results", required=True)
    p.add_argument("--signalp-version", required=True, help="for example 6.0h-gpu (recorded)")
    p.add_argument("--signalp-mode", default="fast")

    p = sub.add_parser("pfam", help="hmmsearch --domtblout -> pfam_adhesion and pfam_allergen")
    _common(p)
    p.add_argument("--domtbl", required=True)
    p.add_argument("--family-table", required=True)
    p.add_argument("--pfam-release", required=True)
    p.add_argument(
        "--pfam-sha256", required=True, help="sha256 of Pfam-A.hmm, from provenance.json"
    )
    p.add_argument("--hmmer-version", required=True, help="for example 3.4 (recorded)")
    p.add_argument(
        "--sp-module", help="step1_rule@R0 table in --workdir, for second_condition=signal_peptide"
    )
    p.add_argument("--tm-module", help="tm table in --workdir, for second_condition=no_tm")

    for name in ("repeat02", "repeat14"):
        p = sub.add_parser(name, help="repeat detector table -> " + name)
        _common(p)
        p.add_argument("--table", required=True)
        p.add_argument("--min-coverage", type=float, default=repeats.DEFAULT_MIN_COVERAGE)
        p.add_argument("--min-copies", type=float, default=repeats.DEFAULT_MIN_COPIES)
        p.add_argument(
            "--script",
            required=True,
            help="detector script that made --table; its sha256 is part of the module identity",
        )

    p = sub.add_parser(
        "allergen",
        help="BLASTP results against the IUIS fungal allergens; an empty BLAST file is refused "
        "because a finished run with no hit cannot be told apart from a failed run",
    )
    _common(p)
    p.add_argument("--blast", required=True, help="outfmt 6 with: " + allergen.BLAST_FIELDS)
    p.add_argument("--allergen-fasta", required=True)
    p.add_argument("--blast-version", required=True)
    p.add_argument("--evalue", required=True, help="BLAST -evalue used by the job (recorded)")
    p.add_argument("--seg", required=True, help="BLAST low-complexity masking used by the job")
    p.add_argument(
        "--max-target-seqs",
        type=_blast_cap,
        required=True,
        help="BLAST -max_target_seqs used by the job; must exceed the database size for recall tests",
    )

    p = sub.add_parser("antigen", help="antigen ranking lookup (Coccidioides)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--ranking", required=True)
    p.add_argument("--protein-map", required=True)

    p = sub.add_parser("cys", help="Cys-rich tiers (analysis/cys_candidates candidates.tsv.gz)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--candidates", required=True)

    p = sub.add_parser("expression", help="spherule table lookup (Coccidioides)")
    _common(p)
    _taxa_args(p)
    p.add_argument("--table", required=True)
    p.add_argument("--protein-map", required=True)

    p = sub.add_parser("tm", help="TMHMM table -> tm")
    _common(p)
    p.add_argument("--table", required=True)
    p.add_argument("--tmhmm-version", default="2.0c")
    return ap


def _taxa(args, proteins):
    tmap = read_taxon_map(args.taxon_map) if args.taxon_map else {}
    if args.taxon is None and not tmap:
        raise RunError("give --taxon or --taxon-map")
    foreign = sorted(set(tmap) - {p.id for p in proteins})
    if tmap and len(foreign) == len(tmap):
        raise RunError(
            f"{args.taxon_map}: none of its {len(tmap)} ID(s) is in the FASTA, "
            f"for example {foreign[0]!r}"
        )
    note = f"{len(foreign)} taxon map ID(s) are not in the FASTA" if foreign else ""
    taxa = assign_taxa(proteins, args.taxon, tmap)
    bad = sorted({t for t in taxa.values() if t < 2})
    if bad:  # 0 is unset and 1 is the root of the taxonomy: neither names an organism
        raise RunError(f"taxon ID {bad[0]} is not a taxon of an organism")
    return taxa, note


def _check_ids(label, found, ids):
    """Every ID in a result table must be a FASTA ID. ``label`` names the file."""
    unknown = sorted(set(found) - set(ids))
    if not unknown:
        return
    example = f"for example {unknown[0]!r}"
    if len(unknown) == len(set(found)):
        raise RunError(
            f"no protein has a result: all {len(unknown)} ID(s) in {label} are not in the FASTA, "
            f"{example}"
        )
    raise RunError(
        f"{label}: {len(unknown)} of {len(set(found))} ID(s) are not in the FASTA, {example}"
    )


def _write(workdir, spec, columns, rows, extra_note=""):
    """Write the module. Refuse a table without one usable result: every protein ``error`` (an ID
    mismatch) or no applicable protein ``ok`` (a lookup that matched nothing). Some ``error`` rows
    make the run state ``partial``."""
    counted = [r for r in rows if r["state"] not in ("na_invalid", "not_applicable")]
    if counted and not any(r["state"] == "ok" for r in counted):
        states = sorted({r["state"] for r in counted})
        raise RunError(
            f"{spec.name}: no protein has a result (states: {', '.join(states)}); "
            "check that the IDs agree"
        )
    errors = sum(1 for r in rows if r["state"] == "error")
    state = "partial" if errors else "ok"
    notes = [f"{errors} protein(s) have no result"] if errors else []
    if not counted and any(r["state"] == "not_applicable" for r in rows):
        notes.append("0 applicable proteins: every protein is not_applicable for this table")
    if extra_note:
        notes.append(extra_note)
    note = "; ".join(notes)
    return write_module(workdir, spec, columns, rows, run_state=state, note=note)


def _family_digest(families):
    """Hash of what decides a call (accession, module, condition, active), not of the free text."""
    import hashlib

    rows = sorted((f.pfam_acc, f.module, f.second_condition, f.active) for f in families)
    return hashlib.sha256(repr(rows).encode()).hexdigest()


def run(args):
    proteins = read_fasta(args.fasta)
    ids = [p.id for p in proteins]
    w = args.workdir
    if args.cmd == "signalp":
        spec = ModuleSpec(
            "step1_rule@R0",
            "1",
            {"rule": "R0", "mode": args.signalp_mode, "organism": "eukarya"},
            (),
            {"signalp": args.signalp_version},
            artefact_digest=_tool_digest("signalp", args.signalp_version),
        )
        parsed = signalp.parse_signalp(args.results)
        _check_ids(args.results, parsed, ids)
        rows = signalp.signalp_rows(proteins, parsed)
        return _write(w, spec, signalp.COLUMNS, rows)
    if args.cmd == "pfam":
        families = pfam.load_family_table(args.family_table)
        hits = pfam.parse_domtblout(args.domtbl)
        _check_ids(args.domtbl, {h["target"] for h in hits}, ids)
        needed = {f.second_condition for f in families if f.active}
        for cond, option, given in (
            ("signal_peptide", "--sp-module", args.sp_module),
            ("no_tm", "--tm-module", args.tm_module),
        ):
            if cond in needed and not given:
                raise RunError(f"an active family has second_condition={cond}: give {option}")
        conditions = {}
        sp_calls = tm_counts = None
        if args.sp_module:
            raw, conditions["sp_module"] = _condition_table(
                w, args.sp_module, "call", ids, "--sp-module"
            )
            sp_calls = {k: v for k, v in raw.items() if v in ("called", "not_called")}
        if args.tm_module:
            raw, conditions["tm_module"] = _condition_table(
                w, args.tm_module, "n_tm_mature", ids, "--tm-module"
            )
            tm_counts = {k: int(v) for k, v in raw.items() if v.isdigit()}
        out = []
        for module in pfam.MODULES:
            active = sorted(f.pfam_acc for f in families if f.module == module and f.active)
            params = {
                "pfam_release": args.pfam_release,
                "cut": "ga",
                "families": active,
                "conditions": conditions,
            }
            spec = ModuleSpec(
                module,
                "1",
                params,
                (),
                {"hmmer": args.hmmer_version},
                artefact_digest=args.pfam_sha256 + ":" + _family_digest(families),
            )
            if not active:  # no family has passed its specificity test: no domain test was made
                rows = [
                    invalid_row(p) if p.state != "ok" else {"id": p.id, "state": "unavailable"}
                    for p in proteins
                ]
                out.append(
                    write_module(
                        w,
                        spec,
                        pfam.COLUMNS,
                        rows,
                        run_state="unavailable",
                        note="no active family",
                    )
                )
                continue
            rows = pfam.pfam_rows(proteins, hits, families, module, sp_calls, tm_counts)
            out.append(_write(w, spec, pfam.COLUMNS, rows))
        return out
    if args.cmd in ("repeat02", "repeat14"):
        params = {
            "min_coverage": args.min_coverage,
            "min_copies": args.min_copies,
            "min_len": repeats.DETECTOR_MIN_LEN,
            "script_sha256": _file_digest(args.script),
        }
        spec = ModuleSpec(args.cmd, "1", params)
        parsed = repeats.parse_repeat_table(args.table)
        _check_ids(args.table, parsed, ids)
        rows = repeats.repeat_rows(proteins, parsed, args.min_coverage, args.min_copies)
        return _write(w, spec, repeats.COLUMNS, rows)
    if args.cmd == "allergen":
        if not any(line.strip() for line in Path(args.blast).read_text().splitlines()):
            raise RunError(f"{args.blast}: the BLAST table is empty (a failed or truncated run?)")
        best = allergen.parse_blast(args.blast)
        _check_ids(args.blast, best, ids)
        meta_path = Path(str(args.allergen_fasta) + ".meta.tsv")
        meta = allergen.read_meta(meta_path) if meta_path.exists() else {}
        spec = ModuleSpec(
            "allergen_homology",
            "1",
            {
                "evalue": args.evalue,
                "program": "blastp",
                "seg": args.seg,
                "max_target_seqs": args.max_target_seqs,
            },
            (args.allergen_fasta,),
            {"blast": args.blast_version},
        )
        rows = allergen.allergen_rows(proteins, best, meta)
        note = f"{len(best)} of {len(proteins)} FASTA protein(s) have a BLAST hit"
        return _write(w, spec, allergen.COLUMNS, rows, extra_note=note)
    if args.cmd == "tm":
        table, _ = lookups.read_table(args.table, "protein_id", unique=True)
        _check_ids(args.table, table, ids)
        spec = ModuleSpec(
            "tm",
            "1",
            {"signal_peptide_window": lookups.SIGNAL_PEPTIDE_WINDOW},
            (),
            {"tmhmm": args.tmhmm_version},
            artefact_digest=_tool_digest("tmhmm", args.tmhmm_version),
        )
        return _write(w, spec, lookups.TM_COLUMNS, lookups.tm_rows(proteins, table))
    taxa, taxa_note = _taxa(args, proteins)
    applicable = set(args.applicable_taxa)
    if args.cmd == "antigen":
        by_gene, several = lookups.ranking_by_gene(args.ranking)
        pmap = lookups.load_protein_map(args.protein_map)
        params = {"applicable_taxa": sorted(applicable), "genes_with_several_ranking_rows": several}
        spec = ModuleSpec("antigen_lookup", "1", params, (args.ranking, args.protein_map))
        rows = lookups.antigen_rows(proteins, taxa, pmap, by_gene, applicable)
        return _write(w, spec, lookups.ANTIGEN_COLUMNS, rows, extra_note=taxa_note)
    if args.cmd == "cys":
        table, _ = lookups.read_table(
            args.candidates, "protein_id", required=("tier", "cys_frac"), unique=True
        )
        spec = ModuleSpec(
            "cys_rich", "1", {"applicable_taxa": sorted(applicable)}, (args.candidates,)
        )
        rows = lookups.cys_rows(proteins, taxa, table, applicable)
        return _write(w, spec, lookups.CYS_COLUMNS, rows, extra_note=taxa_note)
    if args.cmd == "expression":
        table, _ = lookups.read_table(
            args.table,
            "gene_id",
            required=("log2fc_48h", "padj_48h", "log2fc_8d"),
            numeric=EXPRESSION_NUMBERS,
            unique=True,
        )
        pmap = lookups.load_protein_map(args.protein_map)
        spec = ModuleSpec(
            "expression",
            "1",
            {"applicable_taxa": sorted(applicable)},
            (args.table, args.protein_map),
        )
        rows = lookups.expression_rows(proteins, taxa, pmap, table, applicable)
        return _write(w, spec, lookups.EXPRESSION_COLUMNS, rows, extra_note=taxa_note)
    raise AssertionError(args.cmd)


def _tool_digest(tool, version):
    """Identity of a tool: the measurement of a rule belongs to the tool version, not to one input."""
    import hashlib

    return hashlib.sha256(f"{tool}:{version}".encode()).hexdigest()


def _condition_table(workdir, module, column, ids, option):
    """Read one column of a module table that a Pfam second condition needs.

    Refuse a missing table or run record, an unusable run state, a missing column and IDs that are
    not FASTA IDs. Returns ``({id: value}, ...)`` with the values of rows in state ``ok`` only, and
    ``{"params_hash", "artefact_hash"}``.
    """
    import csv
    import gzip

    base = Path(workdir) / "modules"
    table, record = base / f"{module}.tsv.gz", base / f"{module}.json"
    try:
        run = json.loads(record.read_text())
    except (OSError, ValueError) as err:
        raise RunError(f"{option} {module}: cannot read the run record {record}: {err}") from err
    state = run.get("run_state")
    if state in ("unavailable", "not_run", "error"):
        raise RunError(f"{option} {module}: run state of {record} is {state!r}, not usable")
    try:
        with gzip.open(table, "rt", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if column not in (reader.fieldnames or []):
                raise RunError(f"{option} {module}: {table} has no column {column!r}")
            rows = [(r["id"], r["state"], r[column]) for r in reader]
            values = {i: v for i, state, v in rows if state == "ok"}  # other states: no value
    except (EOFError, gzip.BadGzipFile, OSError, csv.Error, KeyError, UnicodeDecodeError) as err:
        raise RunError(f"{option} {module}: cannot read {table}: {err!r}") from err
    _check_ids(str(table), {i for i, _, _ in rows}, ids)
    identity = {k: run.get(k, "") for k in ("params_hash", "artefact_hash")}
    return values, {"module": module, **identity}


def _file_digest(path):
    from cellsurface_sorting_hat.modules.base import sha256_file

    return sha256_file(path)


def main(argv=None):
    try:
        records = run(build_parser().parse_args(argv))
    except (RunError, InputError, FastaError, ValueError, OSError) as err:
        print(f"cellsurface_sorting_hat_module: error: {err}", file=sys.stderr)
        return 2
    for record in records if isinstance(records, list) else [records]:
        print(json.dumps({k: record[k] for k in ("module", "run_state", "n_rows")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
