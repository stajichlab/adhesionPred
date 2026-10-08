"""Command ``cellsurface_sorting_hat_calibrate``: write status sources and check panels."""

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

from cellsurface_sorting_hat.calibration import phasec
from cellsurface_sorting_hat.calibration.call_files import write_call_status
from cellsurface_sorting_hat.calibration.measure import (
    build_measure,
    make_entry,
    write_status_source,
)
from cellsurface_sorting_hat.calibration.panel import panel_check
from cellsurface_sorting_hat.engine import call_eligible, load_config, reads_of_call
from cellsurface_sorting_hat.fasta import read_fasta
from cellsurface_sorting_hat.modules import allergen, pfam
from cellsurface_sorting_hat.taxonomy import TaxonError


def build_parser():
    ap = argparse.ArgumentParser(prog="cellsurface_sorting_hat_calibrate")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser(
        "phasec", help="status source for step1_rule@R0 from the Phase C metrics.json"
    )
    p.add_argument("--workdir", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument(
        "--clusters", required=True, help="Phase C clusters.tsv.gz (seq_sha256, cluster_id)"
    )
    p.add_argument(
        "--eval-table", required=True, help="Phase C eval_table.tsv.gz (origin, class, source_ids)"
    )
    p.add_argument(
        "--phasec-signalp-module",
        required=True,
        help="SignalP environment module that Phase C used, for example signalp/6-gpu",
    )
    p.add_argument(
        "--phasec-signalp-mode",
        required=True,
        help="SignalP mode that Phase C used, for example fast",
    )
    p.add_argument("--set-species", required=True, help="TSV: set_key, scientific_name")
    p.add_argument("--names-dmp", required=True)
    p.add_argument(
        "--nodes-dmp", required=True, help="NCBI nodes.dmp (ranks; species level is enforced)"
    )

    p = sub.add_parser(
        "truth", help="sensitivity and specificity of one call against a truth table"
    )
    p.add_argument("--workdir", required=True)
    target = p.add_mutually_exclusive_group(required=True)
    target.add_argument(
        "--module", help="module that receives the status entry (a call that reads one module)"
    )
    target.add_argument(
        "--call-status",
        action="store_true",
        help="write status/calls/<call>.json: the status of a call that reads two or more modules",
    )
    p.add_argument(
        "--config",
        help="categories.yaml of the run (default: the packaged file); with --call-status its hash "
        "must equal config_sha256 in run.json",
    )
    p.add_argument(
        "--calls-long", required=True, help="calls.long.tsv.gz of a run on the truth proteins"
    )
    p.add_argument("--call", required=True)
    p.add_argument("--variant", default="")
    p.add_argument("--truth", required=True, help="TSV: id, label (1 or 0), cluster")
    p.add_argument("--calibration-set", required=True)
    p.add_argument(
        "--taxa",
        type=int,
        nargs="+",
        required=True,
        help="the ONE tested taxon (species, or a strain or subspecies below a species); one call per species",
    )
    p.add_argument(
        "--nodes-dmp", required=True, help="NCBI nodes.dmp (ranks; species level is enforced)"
    )
    p.add_argument("--notes", default="")
    p.add_argument(
        "--leakage",
        required=True,
        choices=["none", "partial", "tuned_on_truth", "in_reference", "unknown"],
        help="did the truth proteins help to set the rule or its cutoffs? anything but 'none' caps the status at smoke",
    )
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--seed", type=int, default=1)

    p = sub.add_parser(
        "allergen-lso", help="leave-species-out recall of the allergen set (sensitivity only)"
    )
    p.add_argument("--blast", required=True, help="allergens against allergens, outfmt 6")
    p.add_argument(
        "--allergen-fasta", required=True, help="the searched FASTA; fixes the denominator"
    )

    p = sub.add_parser(
        "pfam-specificity", help="specificity test of one Pfam family on one proteome"
    )
    p.add_argument("--family", required=True, help="for example PF05730")
    p.add_argument("--domtbl", required=True, help="hmmsearch --domtblout of the proteome")
    p.add_argument(
        "--members",
        required=True,
        help="TSV: pfam_acc, protein_id (members known from curation, not from Pfam)",
    )
    p.add_argument("--universe-fasta", required=True, help="the proteome that was searched")
    p.add_argument(
        "--tm-table",
        help="TMHMM table (protein_id, pred_hel) to show helices beside each non-member hit",
    )

    p = sub.add_parser("panel", help="report-only check of calls against a panel")
    p.add_argument("--calls-long", required=True)
    p.add_argument("--panel", required=True)
    return ap


LEAKAGE = ("none", "partial", "tuned_on_truth", "in_reference", "unknown")
CALL_VALUES = ("called", "not_called", "not_assessable")


def _check_taxon(taxon, what):
    """A status entry never covers the root (1) or taxon 0."""
    if isinstance(taxon, bool) or not isinstance(taxon, int) or taxon < 2:
        raise TaxonError(
            f"{what}: taxon {taxon!r} is not allowed; give a species-level taxon ID (>= 2)"
        )
    return taxon


def read_nodes(path):
    """``nodes.dmp`` -> ``(parent, rank)`` dicts; parse errors name the path and line."""
    parent, rank = {}, {}
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            f = [x.strip() for x in line.rstrip("\n").split("|")]
            if len(f) < 3 or not f[0]:
                continue
            try:
                parent[int(f[0])], rank[int(f[0])] = int(f[1]), f[2]
            except ValueError:
                raise ValueError(f"{path}:{n}: taxon IDs must be integers") from None
    if not parent:
        raise ValueError(f"{path}: no nodes")
    return parent, rank


def _require_species(taxon, nodes, what):
    """The taxon must be in nodes.dmp and be a species or lie below a species node.

    A genus or any higher rank is refused, as is a ``no rank`` taxon that is not under a species."""
    _check_taxon(taxon, what)
    parent, rank = nodes
    if taxon not in parent:
        raise TaxonError(f"{what}: taxon {taxon} is not in nodes.dmp")
    t, seen = taxon, set()
    while t not in seen:
        seen.add(t)
        if rank[t] == "species":
            return taxon
        if parent[t] == t or parent[t] not in parent:
            break
        t = parent[t]
    raise TaxonError(
        f"{what}: taxon {taxon} has rank {rank[taxon]!r}; a status entry is one species "
        "(species, or below a species)"
    )


def _need_columns(reader, columns, path):
    missing = [c for c in columns if c not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"{path}: missing column(s) {missing}")


def _set_taxa(path, names_dmp, nodes):
    """Phase C set key -> ``[taxon]``. One species per set; the taxon comes from names.dmp only."""
    try:
        names = phasec.read_names(names_dmp)
    except ValueError as err:
        raise ValueError(f"{names_dmp}: cannot parse names.dmp: {err}") from err
    by_set = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        # a line that starts with "#" is a comment; keep the file line number of the other lines
        kept = [(n, line) for n, line in enumerate(fh, 1) if not line.startswith("#")]
        reader = csv.DictReader((line for _, line in kept), delimiter="\t")
        _need_columns(reader, ("set_key", "scientific_name"), path)
        for r in reader:
            where = f"{path}:{kept[reader.line_num - 1][0]}"
            key, name = (r["set_key"] or "").strip(), (r["scientific_name"] or "").strip()
            if not key or not name:
                raise ValueError(f"{where}: set_key and scientific_name must not be empty")
            if key not in phasec.SETS:
                raise ValueError(f"{where}: {key!r} is not a Phase C set ({sorted(phasec.SETS)})")
            if key in by_set:
                raise ValueError(f"{where}: set {key!r} appears twice; one species per set")
            taxon = phasec.species_taxid(names, name)
            by_set[key] = [_require_species(taxon, nodes, f"{where} {name!r}")]
    if not by_set:
        raise ValueError(f"{path}: no set rows")
    return by_set


def _require_run_record(workdir, module):
    path = Path(workdir) / "modules" / f"{module}.json"
    if not path.is_file():
        raise ValueError(f"{path}: module run record not found; run the module first")
    return path


def _merge_entries(path, new_entries, workdir=None, module=None):
    """Entries already in the status source stay, except those of the same calibration set.

    Old entries are dropped (with a message) when the module identity in the workdir differs from
    the identity in the old file: a measurement of an older version must not carry the new one.
    """
    if not Path(path).exists():
        return list(new_entries)
    try:
        old_file = json.loads(Path(path).read_text())
        old = old_file["entries"]
        for e in old:
            e["measure"]["calibration_set"]
    except (ValueError, KeyError, TypeError) as err:
        raise ValueError(f"{path}: cannot read the existing status source: {err!r}") from err
    if workdir is not None:
        now = json.loads(_require_run_record(workdir, module).read_text())
        keys = ("module", "version", "params_hash", "artefact_hash")
        if any(old_file.get(k) != now.get(k) for k in keys):
            print(
                f"{path}: module identity changed; {len(old)} old entr(ies) dropped",
                file=sys.stderr,
            )
            return list(new_entries)
    names = {e["measure"]["calibration_set"] for e in new_entries}
    return [e for e in old if e["measure"]["calibration_set"] not in names] + list(new_entries)


def _plain(text):
    """Text for the notes: no ``=`` or ``;``, so a ``call=`` in it cannot become a ``call=`` field."""
    return " ".join(str(text).replace("=", ":").replace(";", ",").split())


def _check_run_identity(calls_long, workdir, module):
    """``run.json`` next to ``calls_long`` must carry the identity of ``module`` that the work
    directory's module record carries now. A status belongs to the data it was measured on."""
    path = Path(calls_long).with_name("run.json")
    if not path.is_file():
        raise ValueError(f"{path}: not found; it must be next to --calls-long (the run output)")
    try:
        run_json = json.loads(path.read_text(encoding="utf-8-sig"))
        listed = run_json["module_identities"]
        mine = [m for m in listed if m.get("name") == module]
    except (ValueError, KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f"{path}: no usable module_identities ({exc!r})") from exc
    if len(mine) != 1:
        raise ValueError(f"{path}: module_identities has no entry for module {module!r}")
    record = json.loads(_require_run_record(workdir, module).read_text())
    for key in ("version", "params_hash", "artefact_hash"):
        if str(mine[0].get(key)) != str(record.get(key)):
            raise ValueError(
                f"module {module!r}: {key} in {path} ({mine[0].get(key)!r}) differs from "
                f"{Path(workdir) / 'modules' / (module + '.json')} ({record.get(key)!r})"
            )


def _read_call_values(path, call, variant):
    """protein -> value for one call and variant; a conflicting duplicate or an unknown value is refused."""
    out = {}
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        _need_columns(reader, ("protein", "call", "variant", "value"), path)
        for r in reader:
            if r["call"] != call or r["variant"] != variant:
                continue
            where = f"{path}:{reader.line_num}"
            if r["value"] not in CALL_VALUES:
                raise ValueError(
                    f"{where}: value must be one of {CALL_VALUES} (got {r['value']!r})"
                )
            if out.setdefault(r["protein"], r["value"]) != r["value"]:
                raise ValueError(f"{where}: {r['protein']} has two different values for this call")
    return out


def _read_protein_taxa(calls_long):
    """``proteins.tsv.gz`` next to ``calls_long`` -> ``{id: taxon}``."""
    path = Path(calls_long).with_name("proteins.tsv.gz")
    if not path.is_file():
        raise ValueError(
            f"{path}: not found; it must be next to --calls-long (the run output directory)"
        )
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        _need_columns(reader, ("id", "taxon"), path)
        out = {}
        for r in reader:
            try:
                out[r["id"]] = int(r["taxon"])
            except ValueError:
                raise ValueError(
                    f"{path}:{reader.line_num}: taxon {r['taxon']!r} is not an integer"
                ) from None
    return out


def _is_below(taxon, ancestor, parent):
    """True when ``taxon`` is ``ancestor`` or has it on its path to the root."""
    seen = set()
    while taxon not in seen:
        if taxon == ancestor:
            return True
        seen.add(taxon)
        up = parent.get(taxon)
        if up is None or up == taxon:
            return False
        taxon = up
    return False


def _check_truth_taxa(matched, protein_taxa, taxon, parent, proteins_path):
    """Every matched truth protein must be of the tested taxon or below it."""
    bad = []
    for pid in matched:
        t = protein_taxa.get(pid)
        if t is None or not _is_below(t, taxon, parent):
            bad.append((pid, t))
    if bad:
        pid, t = bad[0]
        raise ValueError(
            f"{len(bad)} matched truth protein(s) are not of --taxa {taxon} or below it "
            f"(see {proteins_path}); for example {pid!r} has taxon {'missing' if t is None else t}"
        )


def _read_truth(path, calls):
    """Yield ``(id, label, call value or None, cluster)`` for each truth row; refuses a bad row."""
    seen = set()
    rows = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        _need_columns(reader, ("id", "label", "cluster"), path)
        for r in reader:
            where = f"{path}:{reader.line_num}"
            pid, label, cluster = r["id"], r["label"], (r["cluster"] or "").strip()
            if not pid or not cluster:
                raise ValueError(f"{where}: id and cluster must not be empty")
            if label not in ("0", "1"):
                raise ValueError(f"{where}: label must be 1 or 0 (got {label!r})")
            if pid in seen:
                raise ValueError(f"{where}: id {pid!r} appears twice")
            seen.add(pid)
            rows.append((pid, int(label), calls.get(pid), cluster))
    return rows


def _check_call_run(calls_long, workdir, cfg, call, variant):
    """Checks for ``--call-status``; returns the modules the call reads.

    The call must be eligible. ``run.json`` next to ``calls_long`` must come from a run with this
    config and every read module in state ``ok`` (a call measured with a module missing is a
    different rule). Each module identity must be the one in the work directory now.
    """
    ok, reason = call_eligible(cfg, call)
    if not ok:
        raise ValueError(
            f"unknown call {call!r}"
            if reason == "unknown call"
            else f"call {call!r} is not eligible for a call status file ({reason}); "
            "it needs an expression with no ref, not kind other, reading two or more modules"
        )
    reads = reads_of_call(cfg, call, variant)
    path = Path(calls_long).with_name("run.json")
    if not path.is_file():
        raise ValueError(f"{path}: not found; it must be next to --calls-long (the run output)")
    try:
        run_json = json.loads(path.read_text(encoding="utf-8-sig"))
    except ValueError as exc:
        raise ValueError(f"{path}: not valid JSON ({exc})") from exc
    if "config_sha256" not in run_json:
        raise ValueError(f"{path}: no config_sha256; run the core command again")
    if run_json["config_sha256"] != cfg.sha256:
        raise ValueError(
            f"{path}: the run used another config (config_sha256 differs from --config)"
        )
    states = run_json.get("module_states")
    if not isinstance(states, dict):
        raise ValueError(f"{path}: no module_states; run the core command again")
    for module in reads:
        state = states.get(module, "missing")
        if state != "ok":
            raise ValueError(
                f"module {module!r} was not ok in the run (state: {state}); the call was not "
                "measured with all the modules it reads"
            )
    for module in reads:
        _require_run_record(workdir, module)
        _check_run_identity(calls_long, workdir, module)
    return reads


def run(args):
    if args.cmd == "phasec":
        record = _require_run_record(args.workdir, phasec.R0_MODULE)
        sp_version = phasec.check_signalp(
            record, args.phasec_signalp_module, args.phasec_signalp_mode
        )
        nodes = read_nodes(args.nodes_dmp)
        set_taxa = _set_taxa(args.set_species, args.names_dmp, nodes)
        counts = phasec.count_clusters(args.clusters, args.eval_table)
        entries = phasec.entries_from_phasec(
            args.metrics,
            set_taxa,
            counts,
            extra_notes=(
                f"signalp_module={args.phasec_signalp_module}; "
                f"signalp_mode={args.phasec_signalp_mode}; "
                f"signalp_record_version={_plain(sp_version)}"
            ),
            where_counts=f"{args.eval_table} joined to {args.clusters}",
        )
        path = Path(args.workdir) / "status" / f"{phasec.R0_MODULE}.json"
        return write_status_source(
            args.workdir,
            phasec.R0_MODULE,
            _merge_entries(path, entries, args.workdir, phasec.R0_MODULE),
        )
    if args.cmd == "truth":
        if args.n_boot < 1:
            raise ValueError("--n-boot must be at least 1")
        cfg = load_config(args.config) if args.config else load_config()
        if args.module:
            _require_run_record(args.workdir, args.module)
        if len(args.taxa) != 1:
            raise ValueError(
                "--taxa takes exactly one taxon: a status entry is one species; "
                "run one call per species"
            )
        nodes = read_nodes(args.nodes_dmp)
        _require_species(args.taxa[0], nodes, "--taxa")
        if args.call_status:
            reads = _check_call_run(args.calls_long, args.workdir, cfg, args.call, args.variant)
        else:
            reads = reads_of_call(cfg, args.call, args.variant)
            if args.module not in reads:
                raise ValueError(
                    f"--module {args.module!r} is not read by call {args.call!r}; "
                    f"the call reads: {', '.join(reads)}"
                )
            if len(reads) > 1:
                raise ValueError(
                    f"call {args.call!r} reads {len(reads)} modules ({', '.join(reads)}); a status "
                    "comes from a measurement of one module on a call that reads only that module. "
                    "Give a call that reads one module, or use --call-status"
                )
            _check_run_identity(args.calls_long, args.workdir, args.module)
        protein_taxa = _read_protein_taxa(args.calls_long)
        calls = _read_call_values(args.calls_long, args.call, args.variant)
        y, called, clusters, unmatched, unknown, unknown_pos = [], [], [], 0, 0, 0
        matched = []
        for pid, label, value, cluster in _read_truth(args.truth, calls):
            if value is not None:
                matched.append(pid)
            if value is None:
                unmatched += 1
            elif value == "not_assessable":
                unknown += 1
                unknown_pos += label == 1
            else:
                y.append(label)
                called.append(value == "called")
                clusters.append(cluster)
        if not y:
            raise ValueError("no truth protein has a call in --calls-long")
        _check_truth_taxa(
            matched,
            protein_taxa,
            args.taxa[0],
            nodes[0],
            Path(args.calls_long).with_name("proteins.tsv.gz"),
        )
        n_pos_called = sum(1 for lab, c in zip(y, called, strict=True) if lab == 1 and c)
        n_pos_all = sum(1 for lab in y if lab == 1) + unknown_pos
        bound = f"{n_pos_called / n_pos_all:.3f}" if n_pos_all else "NA"
        notes = (
            f"{args.notes} truth rows without a call: {unmatched}; not assessable: {unknown}; "
            f"sensitivity if not assessable positives count as missed: {bound}"
        ).strip()
        if not any(lab == 0 for lab in y):
            notes += " specificity not measured (no negatives)"
        if args.call_status:
            notes += f" call={args.call}; variant={args.variant}; reads={','.join(reads)}"
        else:
            notes += f" call={args.call}; module={args.module}; variant={args.variant}"
        measure = build_measure(
            args.calibration_set,
            str(args.truth),
            y,
            called,
            clusters,
            notes,
            args.n_boot,
            args.seed,
        )
        measure["notes"] = f"{measure['notes']} leakage: {args.leakage}".strip()
        cap = None if args.leakage == "none" else "smoke"
        entry = make_entry(args.taxa, measure, source=str(args.truth), cap=cap)
        if args.call_status:
            return write_call_status(
                args.workdir, cfg, args.call, args.variant, cfg.sha256, [entry]
            )
        path = Path(args.workdir) / "status" / f"{args.module}.json"
        return write_status_source(
            args.workdir, args.module, _merge_entries(path, [entry], args.workdir, args.module)
        )
    if args.cmd == "pfam-specificity":
        hits = {h["target"] for h in pfam.parse_domtblout(args.domtbl) if h["acc"] == args.family}
        with open(args.members, encoding="utf-8-sig", newline="") as fh:
            members = {
                r["protein_id"]
                for r in csv.DictReader(fh, delimiter="\t")
                if r["pfam_acc"] == args.family
            }
        universe = {p.id for p in read_fasta(args.universe_fasta)}
        outside = members - universe
        if outside:
            raise ValueError(
                f"{len(outside)} member(s) are not in the proteome, for example {sorted(outside)[0]}"
            )
        rep = pfam.specificity_report(hits & universe, members, universe)
        print(f"family\t{args.family}\thits\t{len(hits & universe)}\tmembers\t{len(members)}")
        for key in ("tp", "fp", "fn", "tn", "sensitivity", "specificity"):
            print(f"{key}\t{rep[key]}")
        tm = {}
        if args.tm_table:
            with open(args.tm_table, encoding="utf-8-sig", newline="") as fh:
                tm = {r["protein_id"]: r["pred_hel"] for r in csv.DictReader(fh, delimiter="\t")}
        for pid in rep["nonmember_hits"]:
            print(f"nonmember_hit\t{pid}\tn_tm={tm.get(pid, 'NA')}")
        for pid in rep["missed_members"]:
            print(f"missed_member\t{pid}")
        return None
    if args.cmd == "allergen-lso":
        ids = [p.id for p in read_fasta(args.allergen_fasta)]
        rep = allergen.lso_report(args.blast, ids)
        print(f"sequences\t{rep['n_sequences']}\tspecies\t{rep['n_species']}")
        for r in rep["recall"]:
            print(f"{r['rule']}\t{r['recovered']}/{r['n']}")
        return None
    rows, summary = panel_check(args.calls_long, args.panel)
    for r in rows:
        print(
            "\t".join(
                [r["protein"], r["call"], r["variant"], r["expected"], r["observed"], r["verdict"]]
            )
        )
    print(json.dumps(summary), file=sys.stderr)
    return None


def main(argv=None):
    try:
        run(build_parser().parse_args(argv))
    except (ValueError, KeyError, TaxonError, OSError) as err:
        print(f"cellsurface_sorting_hat_calibrate: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
