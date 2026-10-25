"""Counts per proteome and call, summaries per group, the status matrix and the warnings (standard library only).

Guard: a cell in which fewer than half of the proteins are assessable (``low_assessable``) keeps its raw rate in ``counts_by_proteome`` but is left out of the group
rate summaries and drawn as such in the figures, because the rate is then dominated by which proteins happen to be assessable (a call that is assessable only where it is
called has a rate of 1,000 by construction). Denominator rule (owner decision): rates are per 1,000 ASSESSABLE proteins for that call and proteome (called + not_called), never per all proteins; the number of proteins
is listed beside them. When no protein is assessable the rate is ``NA``, never 0. The status of a call in a proteome is the weakest status among the records that are
called (if any are called), else among the not-called records, else ``unvalidated``; a call with an ``unvalidated`` status is marked and is not "measured".
"""

import statistics

from cellsurface_sorting_hat.status import weakest

PROTEOME_COLS = [
    "proteome",
    "run_dir",
    "taxon",
    "group",
    "group_source",
    "proteins",
    "invalid",
    "config_sha256",
    "engine_version",
    "absent_modules",
]
COUNT_COLS = [
    "proteome", "group", "call", "proteins", "called", "not_called", "not_assessable", "not_applicable", "assessable", "called_per_1000_assessable",
    "assessable_fraction", "low_assessable", "status", "status_mixed", "measured", "unvalidated_marker",
]  # fmt: skip
GROUP_COLS = [
    "group", "call", "n_proteomes", "n_rate", "n_low_assessable", "n_measured", "total_called", "total_assessable", "pooled_per_1000", "mean_per_1000", "median_per_1000", "min_per_1000",
    "max_per_1000", "mean_per_1000_measured", "median_per_1000_measured",
]  # fmt: skip
MARK = "*"  # unvalidated marker in tables and figures
LOW_ASSESSABLE = 0.5  # a cell with fewer than this fraction of assessable proteins is flagged and left out of group rate summaries
NOT_ASSESSED_WARN = 0.05  # warn when more than 5% of the proteins are not assessable for a call


def fmt(x, digits=3):
    return "NA" if x is None else f"{x:.{digits}f}"


def proteome_rows(runs, groups):
    return [
        {"proteome": r.label, "run_dir": r.run_dir, "taxon": r.taxon, "group": groups[r.label][0], "group_source": groups[r.label][1], "proteins": r.n_proteins,
         "invalid": r.n_invalid, "config_sha256": r.config_sha256, "engine_version": r.version, "absent_modules": ",".join(r.absent_modules)}
        for r in runs
    ]  # fmt: skip


def count_rows(runs, groups):
    rows = []
    for r in runs:
        for call in r.call_order:
            recs = r.calls[call]
            n = {
                v: sum(1 for value, _ in recs if value == v)
                for v in ("called", "not_called", "not_assessable", "not_applicable")
            }
            assessable = n["called"] + n["not_called"]
            rate = None if assessable == 0 else 1000.0 * n["called"] / assessable
            frac = assessable / r.n_proteins if r.n_proteins else 0.0
            pool = [s for v, s in recs if v == "called"] or [
                s for v, s in recs if v == "not_called"
            ]
            status = weakest(pool) if pool else "unvalidated"
            rows.append({
                "proteome": r.label, "group": groups[r.label][0], "call": call, "proteins": r.n_proteins, **n, "assessable": assessable,
                "called_per_1000_assessable": rate, "assessable_fraction": frac,
                "low_assessable": "yes" if frac < LOW_ASSESSABLE else "no", "status": status, "status_mixed": "yes" if len(set(pool)) > 1 else "no",
                "measured": "yes" if status != "unvalidated" else "no", "unvalidated_marker": MARK if status == "unvalidated" else "",
            })  # fmt: skip
    return rows


def group_rows(rows, measured_only=False):
    """Summaries per group and call over the proteomes of the group. ``measured_only`` keeps only measured proteomes in every summary column."""
    out, order = [], []
    seen = set()
    for r in rows:
        key = (r["group"], r["call"])
        if key not in seen:
            seen.add(key)
            order.append(key)
    for group, call in order:
        part = [r for r in rows if r["group"] == group and r["call"] == call]
        base = [r for r in part if r["measured"] == "yes"] if measured_only else part
        usable = [
            r
            for r in base
            if r["called_per_1000_assessable"] is not None and r["low_assessable"] == "no"
        ]
        rates = [r["called_per_1000_assessable"] for r in usable]
        meas = [
            r["called_per_1000_assessable"]
            for r in part
            if r["measured"] == "yes"
            and r["called_per_1000_assessable"] is not None
            and r["low_assessable"] == "no"
        ]
        called, assessable = sum(r["called"] for r in usable), sum(r["assessable"] for r in usable)
        out.append({
            "group": group, "call": call, "n_proteomes": len(base), "n_rate": len(rates), "n_low_assessable": sum(1 for r in base if r["low_assessable"] == "yes"), "n_measured": sum(1 for r in part if r["measured"] == "yes"),
            "total_called": called, "total_assessable": assessable, "pooled_per_1000": None if assessable == 0 else 1000.0 * called / assessable,
            "mean_per_1000": statistics.fmean(rates) if rates else None, "median_per_1000": statistics.median(rates) if rates else None,
            "min_per_1000": min(rates) if rates else None, "max_per_1000": max(rates) if rates else None,
            "mean_per_1000_measured": statistics.fmean(meas) if meas else None, "median_per_1000_measured": statistics.median(meas) if meas else None,
        })  # fmt: skip
    return out


def status_matrix(rows):
    """({proteome: {call: status or None}}, proteome order, call order); None where the call has no record in the proteome."""
    proteomes, calls = [], []
    for r in rows:
        if r["proteome"] not in proteomes:
            proteomes.append(r["proteome"])
        if r["call"] not in calls:
            calls.append(r["call"])
    m = {p: dict.fromkeys(calls) for p in proteomes}
    for r in rows:
        m[r["proteome"]][r["call"]] = r["status"]
    return m, proteomes, calls


def warnings(runs, rows, grouped, mixed_config, invalid_warn=0.05):
    out = []
    if mixed_config:
        out.append(
            f"The runs were made with {len({r.config_sha256 for r in runs})} different categories.yaml files; their calls are not defined by the same rules."
        )
    for r in runs:
        if r.n_proteins and r.n_invalid / r.n_proteins > invalid_warn:
            out.append(
                f"{r.label}: {r.n_invalid} of {r.n_proteins} proteins are invalid and excluded from every call."
            )
    for r in runs:
        mine = [row for row in rows if row["proteome"] == r.label]
        none = [row["call"] for row in mine if row["assessable"] == 0]
        if none:
            absent = (
                f"; modules without a table in that run: {', '.join(r.absent_modules)}"
                if r.absent_modules
                else ""
            )
            out.append(
                f"{r.label}: no protein is assessable for {len(none)} call(s) (rate NA): {', '.join(f'`{c}`' for c in none)}{absent}."
            )
        part = [
            f"`{row['call']}` ({100 * row['not_assessable'] / row['proteins']:.0f}%)"
            for row in mine
            if row["assessable"] > 0
            and row["proteins"]
            and row["not_assessable"] / row["proteins"] > NOT_ASSESSED_WARN
        ]
        if part:
            out.append(
                f"{r.label}: more than {int(NOT_ASSESSED_WARN * 100)}% of the proteins are not assessable for {', '.join(part)}; the denominators are smaller."
            )
    for gname in sorted({g["group"] for g in grouped}):
        n = max(g["n_proteomes"] for g in grouped if g["group"] == gname)
        if n <= 2:
            out.append(
                f"Group `{gname}` has {n} proteome(s): a mean over so few is descriptive only."
            )
    return out
