"""Static figures of the aggregation (matplotlib, Agg backend; imported only here). PNG bytes are reproducible: no Software metadata, fixed fonts and size."""

import statistics

STATUS_COLORS = {
    "unvalidated": "#bbbbbb",
    "smoke": "#f0ad4e",
    "estimated": "#5cb85c",
    None: "#ffffff",
}
DPI = 120


def _mpl():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "savefig.dpi": DPI,
            "figure.dpi": DPI,
            "axes.unicode_minus": False,
        }
    )
    return matplotlib, plt


def _save(fig, path, plt):
    fig.savefig(path, format="png", metadata={"Software": None}, bbox_inches="tight")
    plt.close(fig)


def heatmap(rows, proteomes, calls, path):
    """Called per 1,000 assessable proteins: rows = proteomes (in the given order), columns = calls; NA grey; unvalidated cells hatched and marked *."""
    import numpy as np
    from matplotlib.patches import Rectangle  # noqa: I001  (after the backend is chosen in _mpl)

    _, plt = _mpl()
    by = {(r["proteome"], r["call"]): r for r in rows}
    z = np.full((len(proteomes), len(calls)), np.nan)
    for i, p in enumerate(proteomes):
        for j, c in enumerate(calls):
            r = by.get((p, c))
            if (
                r is not None
                and r["called_per_1000_assessable"] is not None
                and r["low_assessable"] == "no"
            ):
                z[i, j] = np.log10(1.0 + r["called_per_1000_assessable"])
    fig, ax = plt.subplots(
        figsize=(max(6.0, 0.55 * len(calls) + 3), max(3.0, 0.38 * len(proteomes) + 2))
    )
    cmap = plt.get_cmap("viridis").with_extremes(bad="#dddddd")
    im = ax.imshow(np.ma.masked_invalid(z), aspect="auto", cmap=cmap)
    for i, p in enumerate(proteomes):
        for j, c in enumerate(calls):
            r = by.get((p, c))
            if r is None or r["called_per_1000_assessable"] is None:
                ax.text(j, i, "NA", ha="center", va="center", fontsize=6, color="#555555")
                continue
            if r["low_assessable"] == "yes":
                z[i, j] = np.nan
                ax.text(j, i, "<50%", ha="center", va="center", fontsize=5, color="#555555")
                continue
            ax.text(
                j,
                i,
                f"{r['called_per_1000_assessable']:.1f}",
                ha="center",
                va="center",
                fontsize=6,
                color="white" if z[i, j] < 0.6 * np.nanmax(z) else "black",
            )
            if r["measured"] == "no":
                ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        fill=False,
                        hatch="///",
                        edgecolor="#ffffff",
                        linewidth=0,
                        alpha=0.55,
                    )
                )
    ax.set_xticks(range(len(calls)), calls, rotation=90, fontsize=7)
    ax.set_yticks(range(len(proteomes)), proteomes, fontsize=7)
    ax.set_title(
        "Called per 1,000 assessable proteins (hatched = status unvalidated; <50% = fewer than half of the proteins assessable)",
        fontsize=8,
    )
    fig.colorbar(im, ax=ax, label="log10(1 + called per 1,000 assessable)", shrink=0.8)
    _save(fig, path, plt)


def group_bars(rows, calls, path):
    """One panel per call: mean over the proteomes of each group, the number of proteomes above each bar, +/- SD only where a group has 3 or more proteomes."""
    import math

    _, plt = _mpl()
    groups = []
    for r in rows:
        if r["group"] not in groups:
            groups.append(r["group"])
    ncol = min(5, len(calls))
    nrow = math.ceil(len(calls) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.3 * nrow + 0.8), squeeze=False)
    for k, call in enumerate(calls):
        ax = axes[k // ncol][k % ncol]
        means, errs, ns = [], [], []
        for g in groups:
            rates = [
                r["called_per_1000_assessable"]
                for r in rows
                if r["group"] == g
                and r["call"] == call
                and r["called_per_1000_assessable"] is not None
            ]
            means.append(statistics.fmean(rates) if rates else 0.0)
            errs.append(statistics.stdev(rates) if len(rates) >= 3 else 0.0)
            ns.append(len(rates))
        ax.bar(
            range(len(groups)),
            means,
            yerr=[e if n >= 3 else 0 for e, n in zip(errs, ns, strict=True)],
            color="#4a7ebb",
            capsize=2,
        )
        for x, (m, n) in enumerate(zip(means, ns, strict=True)):
            ax.text(x, m, f"n={n}", ha="center", va="bottom", fontsize=5)
        ax.set_xticks(range(len(groups)), groups, rotation=60, ha="right", fontsize=5)
        ax.set_title(call, fontsize=6)
        ax.tick_params(axis="y", labelsize=6)
    for k in range(len(calls), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle(
        "Mean called per 1,000 assessable proteins by group (n = proteomes with a rate; bars with SD only for n >= 3)",
        fontsize=8,
    )
    fig.tight_layout()
    _save(fig, path, plt)


def status_panel(matrix, proteomes, calls, path):
    """Colour = status of each call in each proteome (grey unvalidated, amber smoke, green estimated, white no record)."""
    import numpy as np
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    _, plt = _mpl()
    order = [None, "unvalidated", "smoke", "estimated"]
    z = np.array([[order.index(matrix[p][c]) for c in calls] for p in proteomes])
    cmap = ListedColormap([STATUS_COLORS[s] for s in order])
    fig, ax = plt.subplots(
        figsize=(max(6.0, 0.5 * len(calls) + 3), max(3.0, 0.35 * len(proteomes) + 2))
    )
    ax.imshow(z, aspect="auto", cmap=cmap, vmin=0, vmax=3)
    ax.set_xticks(range(len(calls)), calls, rotation=90, fontsize=7)
    ax.set_yticks(range(len(proteomes)), proteomes, fontsize=7)
    ax.set_title("Status of each call in each proteome", fontsize=9)
    ax.legend(
        handles=[Patch(color=STATUS_COLORS[s], label=str(s) if s else "no record") for s in order],
        fontsize=6,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
    )
    _save(fig, path, plt)
