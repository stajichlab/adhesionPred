"""Summarize round 1 of the cross-clade domain search (search.sbatch output) into small tables without sequences.

Declared before the run: a hit is E <= 1e-3 (primary cut-off); E <= 1e-2 and E <= 0.1 are reported beside it. Query hits to any query
accession are dropped (self and paralog queries, so the H99 panel cannot score its own paralogs). Decoys (one reversal and --n-shuffles
local-window shuffles per query, from build_inputs.py) are searched like the queries; all their hits are kept and counted as false hits.

  python3 -I analysis/crossclade_domain_search/summarize.py --inputs /path/build_inputs_out --results /path/search_out \
      --flags docs/reports/data/aga1_nonyeast_test_2026-10-26/all_proteins_flags.tsv --out /path/NEW_DIR

Writes hits.tsv (every query hit at E <= 0.1, with target coordinates, coverage and surface flags), decoy_hits.tsv (every decoy hit at
E <= 0.1), summary.tsv (method x panel x cut-off: distinct targets hit, decoy hits and decoys with a hit, targets with a signal peptide or
GPI, PF21671 members recovered among non-query proteins) and per_query.tsv (targets hit per species, number of decoys, best decoy E-value).
"""

import argparse
import csv
import sys
from pathlib import Path

CUTOFFS = (1e-3, 1e-2, 0.1)
METHODS = ("foldseek", "mmseqs")
COLS = (
    "query",
    "target",
    "fident",
    "alnlen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "bits",
    "qlen",
    "tlen",
)
FLAG_COLS = ("sp_yes", "gpi_yes", "disorder_longest_run", "broad_flag", "narrow_flag")


class SummaryError(Exception):
    pass


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_hits(path):
    """Read a headerless search table with COLS; fail on a short row."""
    out = []
    with open(path, newline="") as fh:
        for n, row in enumerate(csv.reader(fh, delimiter="\t"), 1):
            if not row:
                continue
            if len(row) != len(COLS):
                raise SummaryError(f"{path}:{n}: {len(row)} columns, expected {len(COLS)}")
            d = dict(zip(COLS, row, strict=True))
            d["evalue"] = float(d["evalue"])
            out.append(d)
    return out


def query_accession(qid):
    """`ACC_PFxxxxx_start_end` or a decoy id `decNN_ACC_...` -> ACC."""
    return decoy_query(qid).rsplit("_", 3)[0] if is_decoy(qid) else qid.rsplit("_", 3)[0]


def is_decoy(qid):
    return qid.startswith("dec") and qid[3:5].isdigit() and qid[5:6] == "_"


def decoy_query(did):
    """`decNN_<query_id>` -> `<query_id>`."""
    if not is_decoy(did):
        raise SummaryError(f"not a decoy id: {did}")
    return did.split("_", 1)[1]


def _cov(start, end, length):
    return f"{(int(end) - int(start) + 1) / int(length):.2f}"


def summarize(domains, targets, truth, flags, results, decoy_ids):
    """Query hits to ANY query accession are dropped (self and paralog queries); decoy hits are all kept (a decoy hitting its own
    forward protein is a real false hit). PF21671 recovery counts only members that are not themselves queries."""
    dom = {r["query_id"]: r for r in domains}
    qaccs = {r["accession"] for r in domains}
    tsp = {r["accession"]: r["species"] for r in targets}
    species = sorted(set(tsp.values()))
    pf = {r["accession"] for r in truth if r["in_targets"] == "yes"} - qaccs
    n_dec = {}
    for did in decoy_ids:
        qid = decoy_query(did)
        if qid not in dom:
            raise SummaryError(f"decoy {did} has no query in query_domains.tsv")
        n_dec[qid] = n_dec.get(qid, 0) + 1
    hits_out, decoy_out, summary, per_query = [], [], [], []
    for method in METHODS:
        q = [
            h
            for h in read_hits(Path(results) / f"{method}_queries.tsv")
            if h["target"] not in qaccs
        ]
        d = read_hits(Path(results) / f"{method}_decoys.tsv")
        for h in q + d:
            if h["target"] not in tsp:
                raise SummaryError(f"{method}: target {h['target']} not in target_proteins.tsv")
        for h in q:
            if h["query"] not in dom:
                raise SummaryError(f"{method}: query {h['query']} not in query_domains.tsv")
        for h in d:
            if h["query"] not in decoy_ids:
                raise SummaryError(f"{method}: decoy {h['query']} not in decoys.faa")
        for h in q:
            if h["evalue"] <= max(CUTOFFS):
                f = flags.get(h["target"], {})
                hits_out.append(
                    [
                        method,
                        h["query"],
                        dom[h["query"]]["panel"],
                        dom[h["query"]]["gene"],
                        h["target"],
                        tsp[h["target"]],
                        f"{h['evalue']:.3g}",
                        h["bits"],
                        h["fident"],
                        h["alnlen"],
                        h["tstart"],
                        h["tend"],
                        _cov(h["qstart"], h["qend"], h["qlen"]),
                        _cov(h["tstart"], h["tend"], h["tlen"]),
                        *(f.get(c, "NA") for c in FLAG_COLS),
                        "yes" if h["target"] in pf else "no",
                    ]
                )
        for h in d:
            if h["evalue"] <= max(CUTOFFS):
                qid = decoy_query(h["query"])
                decoy_out.append(
                    [
                        method,
                        h["query"],
                        qid,
                        dom[qid]["panel"],
                        h["target"],
                        tsp[h["target"]],
                        f"{h['evalue']:.3g}",
                        h["bits"],
                        _cov(h["qstart"], h["qend"], h["qlen"]),
                    ]
                )
        for panel in sorted({r["panel"] for r in domains}):
            pq = [h for h in q if dom[h["query"]]["panel"] == panel]
            pd = [h for h in d if dom[decoy_query(h["query"])]["panel"] == panel]
            n_panel_dec = sum(n for qid, n in n_dec.items() if dom[qid]["panel"] == panel)
            for cut in CUTOFFS:
                tq = {h["target"] for h in pq if h["evalue"] <= cut}
                dh = [h for h in pd if h["evalue"] <= cut]
                row = {
                    "method": method,
                    "panel": panel,
                    "evalue_max": cut,
                    "targets_hit": len(tq),
                    "decoy_hits": len(dh),
                    "decoy_targets_hit": len({h["target"] for h in dh}),
                    "decoys_with_hit": f"{len({h['query'] for h in dh})}/{n_panel_dec}",
                    "targets_hit_with_sp": sum(
                        1 for t in tq if flags.get(t, {}).get("sp_yes") == "1"
                    ),
                    "targets_hit_with_gpi": sum(
                        1 for t in tq if flags.get(t, {}).get("gpi_yes") == "1"
                    ),
                }
                for sp in species:
                    members = {a for a in pf if tsp[a] == sp}
                    row[f"pf21671_recovered_{sp}"] = f"{len(tq & members)}/{len(members)}"
                    row[f"non_pf21671_hits_{sp}"] = len({t for t in tq if tsp[t] == sp} - members)
                summary.append(row)
        for qid in sorted(dom):
            de = [h["evalue"] for h in d if decoy_query(h["query"]) == qid]
            best = f"{min(de):.3g}" if de else "NA"
            for cut in CUTOFFS:
                tq = {h["target"] for h in q if h["query"] == qid and h["evalue"] <= cut}
                per_query.append(
                    [
                        method,
                        qid,
                        dom[qid]["panel"],
                        dom[qid]["gene"],
                        cut,
                        *(len({t for t in tq if tsp[t] == sp}) for sp in species),
                        n_dec.get(qid, 0),
                        best,
                    ]
                )
    return hits_out, decoy_out, summary, per_query


def write(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        out = Path(a.out)
        if out.exists():
            raise SummaryError(f"{out} exists; nothing is overwritten")
        inp = Path(a.inputs)
        domains, targets, truth = (
            read_tsv(inp / f)
            for f in ("query_domains.tsv", "target_proteins.tsv", "truth_pf21671.tsv")
        )
        flags = {r["accession"]: r for r in read_tsv(a.flags)}
        missing = [r["accession"] for r in targets if r["accession"] not in flags]
        if missing:
            raise SummaryError(
                f"{len(missing)} targets have no row in {a.flags} (first: {missing[0]})"
            )
        decoy_ids = [
            ln[1:].split()[0]
            for ln in (inp / "decoys.faa").read_text().splitlines()
            if ln.startswith(">")
        ]
        hits, decoys, summary, per_query = summarize(
            domains, targets, truth, flags, a.results, set(decoy_ids)
        )
        species = sorted({r["species"] for r in targets})
        out.mkdir(parents=True)
        write(
            out / "hits.tsv",
            [
                "method",
                "query_id",
                "panel",
                "gene",
                "target",
                "target_species",
                "evalue",
                "bits",
                "fident",
                "alnlen",
                "tstart",
                "tend",
                "qcov",
                "tcov",
                *FLAG_COLS,
                "in_pf21671",
            ],
            hits,
        )
        write(
            out / "decoy_hits.tsv",
            [
                "method",
                "decoy_id",
                "query_id",
                "panel",
                "target",
                "target_species",
                "evalue",
                "bits",
                "qcov",
            ],
            decoys,
        )
        write(out / "summary.tsv", list(summary[0]), [list(r.values()) for r in summary])
        write(
            out / "per_query.tsv",
            [
                "method",
                "query_id",
                "panel",
                "gene",
                "evalue_max",
                *(f"targets_{s}" for s in species),
                "n_decoys",
                "best_decoy_evalue",
            ],
            per_query,
        )
        print(f"wrote {out}: {len(hits)} hit rows")
    except SummaryError as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
