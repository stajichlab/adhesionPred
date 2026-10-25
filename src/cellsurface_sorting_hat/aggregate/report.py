"""The Markdown report of an aggregation (text from the tables; no number is computed here that is not in a table)."""

from cellsurface_sorting_hat.aggregate.tables import MARK, fmt

LIMITS = [
    "A count is a count of calls. It is not a prevalence and not corrected for false positives; expected false calls need a measured specificity, which most calls do not have.",
    "A proteome-call cell with fewer than half of its proteins assessable is flagged (`low_assessable`) and left out of the group summaries: its rate depends on which proteins happen to be assessable "
    "(for example a call that is assessable only where it is called has a rate of 1,000 by construction).",
    f"Rates are per 1,000 assessable proteins for that call in that proteome (a protein that is not assessable for the call is left out of the denominator). Proteomes differ in size, "
    f"annotation quality and completeness; a mean over proteomes of unequal size is a descriptive number, not an estimate, and `n` is printed beside it. A call marked {MARK} has the status "
    "`unvalidated` in that proteome.",
    "The status of a call is for one species and taxon and only where its status file is installed and current. A call that reads a module without a table in the run is "
    "`not_assessable` for every protein and is shown as NA.",
    "A call that uses clade thresholds (for example `model_surface_candidate`) is calibrated only for the clades its thresholds table lists as calibrated; elsewhere its count is not comparable.",
]


def table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def render_report(params, proteomes, counts, grouped, warn, figures):
    L = ["# cellsurface_sorting_hat aggregate report", ""]
    L += [
        f"- proteomes: {len(proteomes)}; groups: {len({p['group'] for p in proteomes})} (rank `{params['rank']}`"
        + (", overridden by a group table" if params["groups_table"] else "")
        + ")"
    ]
    L += [
        f"- calls read at the step 1 variant: `{params['variant']}` (ungated calls have no variant)"
    ]
    L += [
        f"- categories.yaml sha256 of the runs: {', '.join(sorted({p['config_sha256'] for p in proteomes}))}"
    ]
    L += [
        f"- summaries over: {'measured proteomes only' if params['measured_only'] else 'all proteomes'}",
        "",
    ]
    L += (
        ["## Proteomes", ""]
        + table(
            ["proteome", "taxon", "group", "proteins", "invalid", "modules without a table"],
            [
                [
                    p["proteome"],
                    p["taxon"],
                    p["group"],
                    p["proteins"],
                    p["invalid"],
                    p["absent_modules"] or "none",
                ]
                for p in proteomes
            ],
        )
        + [""]
    )
    L += [
        "## Group summaries",
        "",
        f"Called per 1,000 assessable proteins. `n` = proteomes of the group, `meas` = proteomes in which the call has a measured status; {MARK} = no measured status in any of them.",
        "",
    ]
    calls = list(dict.fromkeys(c["call"] for c in counts))
    for call in calls:
        rows = [g for g in grouped if g["call"] == call]
        L += [f"### `{call}`", ""]
        L += table(["group", "n", "low", "meas", "mean [min, max]", "median", "pooled", "called / assessable"], [
            [g["group"] + ("" if g["n_measured"] else MARK), g["n_proteomes"], g["n_low_assessable"], g["n_measured"], f"{fmt(g['mean_per_1000'], 2)} [{fmt(g['min_per_1000'], 2)}, {fmt(g['max_per_1000'], 2)}]",
             fmt(g["median_per_1000"], 2), fmt(g["pooled_per_1000"], 2), f"{g['total_called']} / {g['total_assessable']}"] for g in rows
        ]) + [""]  # fmt: skip
    L += [
        "## Status coverage",
        "",
        f"{sum(1 for c in counts if c['measured'] == 'yes')} of {len(counts)} proteome-call cells have a measured status.",
        "",
    ]
    L += [f"![status of each call in each proteome]({figures['status']})", ""]
    L += [
        "## Figures",
        "",
        f"![called per 1,000 assessable proteins]({figures['heatmap']})",
        "",
        f"![mean per group]({figures['groups']})",
        "",
    ]
    L += (
        ["## Warnings", ""] + ([f"- **WARNING: {w}**" for w in warn] if warn else ["- none"]) + [""]
    )
    L += ["## Known limits", ""] + [f"- {x}" for x in LIMITS] + [""]
    return "\n".join(L)
