"""Groups of taxa for the aggregation: NCBI lineage rank by default, an owner table that overrides it.

``--rank RANK`` walks from the taxon of a proteome up the lineage (``nodes.dmp``) to the first node of that rank; the label is its scientific name when ``names.dmp``
is given, else ``taxon:<id>``. A group table (TSV ``taxon_id``, ``group``, optional ``note``) overrides the lineage group for the taxon it lists and for taxa below it that
have no more specific row. A taxon that is not in ``nodes.dmp`` stops the aggregation; a proteome with no node of the rank above it goes to ``no <rank> above``.
"""

import csv

from cellsurface_sorting_hat.aggregate.inputs import AggregateError


def read_names(path):
    """``names.dmp`` -> {taxon: scientific name}."""
    out = {}
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            f = [x.strip() for x in line.rstrip("\n").split("|")]
            if len(f) >= 4 and f[3] == "scientific name":
                try:
                    out[int(f[0])] = f[1]
                except ValueError:
                    continue
    return out


def read_group_table(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not {"taxon_id", "group"} <= set(reader.fieldnames):
            raise AggregateError(f"{path}: needs the columns 'taxon_id' and 'group'")
        table = {}
        for n, r in enumerate(reader, 2):
            try:
                t = int(r["taxon_id"])
            except (TypeError, ValueError):
                raise AggregateError(f"{path}:{n}: taxon_id is not an integer") from None
            g = (r["group"] or "").strip()
            if not g:
                raise AggregateError(f"{path}:{n}: empty group")
            if t in table:
                raise AggregateError(f"{path}:{n}: taxon {t} is listed twice")
            table[t] = g
    if not table:
        raise AggregateError(f"{path}: no rows")
    return table


def _chain(taxon, parent):
    if taxon not in parent:
        raise AggregateError(f"taxon {taxon} is not in nodes.dmp")
    chain, seen = [taxon], {taxon}
    while (
        parent[chain[-1]] != chain[-1]
        and parent[chain[-1]] in parent
        and parent[chain[-1]] not in seen
    ):
        chain.append(parent[chain[-1]])
        seen.add(chain[-1])
    return chain


def assign_groups(runs, parent, rank, rank_name, names=None, table=None):
    """{label: (group, source)} where source is ``table`` or ``lineage``."""
    names = names or {}
    out = {}
    for r in runs:
        chain = _chain(r.taxon, parent)
        hit = next((t for t in chain if table and t in table), None)
        if hit is not None:
            out[r.label] = (table[hit], "table")
            continue
        node = next((t for t in chain if rank.get(t) == rank_name), None)
        out[r.label] = (
            names.get(node, f"taxon:{node}") if node is not None else f"no {rank_name} above",
            "lineage",
        )
    return out
