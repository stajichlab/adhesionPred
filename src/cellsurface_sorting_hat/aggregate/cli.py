"""Command ``cellsurface_sorting_hat_aggregate``: counts per call across the runs of many proteomes, grouped by taxa, with a report and figures."""

import argparse
import csv
import sys
from pathlib import Path

from cellsurface_sorting_hat.aggregate import figures, groups, inputs, report, tables
from cellsurface_sorting_hat.calibration.cli import read_nodes


def build_parser():
    ap = argparse.ArgumentParser(prog="cellsurface_sorting_hat_aggregate", description=__doc__)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--manifest",
        help="TSV with run_dir and an optional label (one engine output directory per proteome)",
    )
    src.add_argument(
        "--runs",
        nargs="+",
        help="engine output directories (the label is the name of the directory, or of its parent for `out`)",
    )
    ap.add_argument("--taxdump", required=True, help="nodes.dmp of the NCBI taxonomy")
    ap.add_argument("--names", help="names.dmp (readable group labels)")
    ap.add_argument(
        "--rank", default="class", help="lineage rank of the default groups (default: class)"
    )
    ap.add_argument(
        "--groups", help="TSV taxon_id, group[, note]; overrides the lineage group for those taxa"
    )
    ap.add_argument("--variant", help="step 1 variant to read (default: the gate of the runs)")
    ap.add_argument(
        "--measured-only",
        action="store_true",
        help="summaries over proteomes with a measured status only",
    )
    ap.add_argument("--allow-mixed-config", action="store_true")
    ap.add_argument("--calls", help="comma-separated calls for the figures (default: every call)")
    ap.add_argument("--out", required=True, help="new output directory")
    return ap


def _tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow(
                [
                    "NA"
                    if r.get(c) is None
                    else (f"{r[c]:.4f}" if isinstance(r.get(c), float) else r.get(c, ""))
                    for c in cols
                ]
            )


def run(a):
    out = Path(a.out)
    if out.exists():
        raise inputs.AggregateError(f"{out} exists; nothing is overwritten")
    entries = (
        inputs.read_manifest(a.manifest)
        if a.manifest
        else [(d, inputs.default_label(d)) for d in a.runs]
    )
    runs = inputs.read_runs(entries, a.variant, a.allow_mixed_config)
    parent, rank = read_nodes(a.taxdump)
    names = groups.read_names(a.names) if a.names else None
    table = groups.read_group_table(a.groups) if a.groups else None
    grp = groups.assign_groups(runs, parent, rank, a.rank, names, table)
    runs.sort(key=lambda r: (grp[r.label][0], r.label))
    prows = tables.proteome_rows(runs, grp)
    crows = tables.count_rows(runs, grp)
    grows = tables.group_rows(crows, a.measured_only)
    mixed = len({r.config_sha256 for r in runs}) > 1
    warn = tables.warnings(runs, crows, grows, mixed)
    calls = (
        [c.strip() for c in a.calls.split(",")]
        if a.calls
        else list(dict.fromkeys(c["call"] for c in crows))
    )
    unknown = [c for c in calls if c not in {r["call"] for r in crows}]
    if unknown:
        raise inputs.AggregateError(f"--calls names calls that are not in the runs: {unknown}")
    out.mkdir(parents=True)
    _tsv(out / "proteomes.tsv", tables.PROTEOME_COLS, prows)
    _tsv(out / "counts_by_proteome.tsv", tables.COUNT_COLS, crows)
    _tsv(out / "counts_by_group.tsv", tables.GROUP_COLS, grows)
    matrix, proteomes, all_calls = tables.status_matrix(crows)
    usable_calls = {
        r["call"]
        for r in crows
        if r["called_per_1000_assessable"] is not None and r["low_assessable"] == "no"
    }
    shown = [
        c for c in all_calls if c in calls and c in usable_calls
    ]  # a call with no usable rate in any proteome has nothing to draw
    figures.heatmap([r for r in crows if r["call"] in shown], proteomes, shown, out / "heatmap.png")
    figures.group_bars([r for r in crows if r["call"] in shown], shown, out / "groups.png")
    figures.status_panel(matrix, proteomes, shown, out / "status_coverage.png")
    params = {
        "rank": a.rank,
        "groups_table": bool(a.groups),
        "variant": runs[0].variant or "(none)",
        "measured_only": a.measured_only,
    }
    text = report.render_report(
        params,
        prows,
        crows,
        grows,
        warn,
        {"heatmap": "heatmap.png", "groups": "groups.png", "status": "status_coverage.png"},
    )
    (out / "report.md").write_text(text)
    print(
        f"wrote {out}: {len(runs)} proteomes, {len({p['group'] for p in prows})} groups, {len(warn)} warning(s)"
    )
    return 0


def main(argv=None):
    try:
        return run(build_parser().parse_args(argv))
    except (inputs.AggregateError, OSError, ValueError) as err:
        print(f"cellsurface_sorting_hat_aggregate: error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
