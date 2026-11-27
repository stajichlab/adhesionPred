"""Summarize round 2 of the cross-clade domain search: three structure methods and the 2-of-3 rule.

Methods, each put on ONE scale, the round 1 E-value (expected chance hits in a search of the N1 round 1 target proteins):
  prostt5   round 1 Foldseek with ProstT5-predicted 3Di; its E-value as is (round 1 hits.tsv, E <= 0.1 kept there)
  af_fs     Foldseek on the AlphaFold DB models; E-value x N1 / (number of target models searched)
  reseek    Reseek v3 on the AlphaFold DB models; P-value (per pair, independent of database size) x N1

Declared before the run: a method SUPPORTS a query-target pair when its scaled score is <= 1e-3 AND better than that query's null:
  prostt5: the best E-value of the query's 51 sequence decoys (round 1 per_query.tsv);
  af_fs, reseek: the best scaled score of the query against background half A, leaving out the target itself.
A target is CONFIRMED when at least one query pair has support from 2 of the 3 methods. Background half B is scored by the same
rule (null from half A) and gives the rate at which proteins drawn at random pass, per method and for the 2-of-3 rule.
Caveats: af_fs and reseek read the same AFDB model, and af_fs and prostt5 share Foldseek's scoring, so the three votes are not
independent; every scored target was selected by a round 1 hit at E <= 0.1. Half B measures the rule as a whole.

  python3 -I analysis/crossclade_domain_search/round2_summarize.py --round1-inputs /path/round1_in \
      --round1-summary /path/round1_summary --inputs2 /path/round2_in --results2 /path/round2_out --out /path/NEW_DIR

Writes pairs.tsv, targets.tsv, nulls.tsv, background_hits.tsv and summary.tsv (no sequences).
"""

import argparse
import csv
import sys
from pathlib import Path

CUT = 1e-3
METHODS = ("prostt5", "af_fs", "reseek")
FS_COLS = (
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
RS_COLS = ("query", "target", "qlo", "qhi", "ql", "tlo", "thi", "tl", "pctid", "pvalue")


class SummaryError(Exception):
    pass


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm_label(label, known):
    """Map a structure label from a search tool (path, file name, maybe a chain suffix) to a known id."""
    name = Path(label.split()[0]).name
    for ext in (".gz", ".pdb", ".cif"):
        if name.endswith(ext):
            name = name[: -len(ext)]
    for cand in (name, name.rsplit("_", 1)[0], name[:-1], name.rsplit("-", 1)[0]):
        if cand in known:
            return cand
    raise SummaryError(f"structure label {label!r} matches no input id")


def read_scores(path, cols, score_col, known, scale=1.0):
    """{(query, target): best score x scale} from a headerless (or `query`-headed) search table."""
    best = {}
    with open(path, newline="") as fh:
        for n, row in enumerate(csv.reader(fh, delimiter="\t"), 1):
            if not row or (n == 1 and row[0] == "query"):
                continue
            if len(row) != len(cols):
                raise SummaryError(f"{path}:{n}: {len(row)} columns, expected {len(cols)}")
            d = dict(zip(cols, row, strict=True))
            key = (norm_label(d["query"], known), norm_label(d["target"], known))
            s = float(d[score_col]) * scale
            if key not in best or s < best[key]:
                best[key] = s
    return best


def fmt(x):
    return "NA" if x is None else f"{x:.3g}"


def summarize(domains, models, r1_hits, r1_per_query, af, rs, truth):
    """af and rs are already on the round 1 E-value scale."""
    dom = {d["query_id"]: d for d in domains}
    qaccs = {d["accession"] for d in domains}
    mod = {m["id"]: m for m in models}
    tmod = {i: m for i, m in mod.items() if m["role"] != "query"}
    ok = {i for i, m in tmod.items() if m["status"] == "ok"}
    sp = {i: m["species"] for i, m in mod.items()}
    pf = {r["accession"] for r in truth if r["in_targets"] == "yes"} - qaccs
    half = {h: {i for i in ok if tmod[i]["bg_half"] == h} for h in ("A", "B")}
    scored = {i for i, m in tmod.items() if m["role"] in ("candidate", "pf21671")}
    for q, t in list(af) + list(rs):
        if q not in dom:
            raise SummaryError(f"query {q} not in query_domains.tsv")
        if t not in tmod:
            raise SummaryError(f"target {t} not in models.tsv")
    flags, prost, mm = {}, {}, {}
    for h in r1_hits:
        flags[h["target"]] = h
        tab = prost if h["method"] == "foldseek" else mm
        key = (h["query_id"], h["target"])
        e = float(h["evalue"])
        if key not in tab or e < tab[key]:
            tab[key] = e
    tabs = {"prostt5": prost, "af_fs": af, "reseek": rs}
    decoy = {}
    for r in r1_per_query:
        if r["method"] == "foldseek" and r["best_decoy_evalue"] != "NA":
            decoy[r["query_id"]] = float(r["best_decoy_evalue"])
    a_scores = {}  # (q, method) -> sorted [(score, acc)] over half A
    for m in ("af_fs", "reseek"):
        for (q, t), s in tabs[m].items():
            if t in half["A"]:
                a_scores.setdefault((q, m), []).append((s, t))
    for v in a_scores.values():
        v.sort()

    def null(q, m, t):
        if m == "prostt5":
            return decoy.get(q)
        for s, acc in a_scores.get((q, m), []):
            if acc != t:
                return s
        return None

    def pair(q, t):
        s = {m: tabs[m].get((q, t)) for m in METHODS}
        sup = {}
        for m in METHODS:
            n = null(q, m, t)
            sup[m] = s[m] is not None and s[m] <= CUT and (n is None or s[m] < n)
        return s, sup

    pairs = []
    keys = {k for tab in tabs.values() for k in tab if k[1] in scored}
    for q, t in sorted(keys):
        s, sup = pair(q, t)
        pairs.append(
            {
                "query_id": q,
                "panel": dom[q]["panel"],
                "gene": dom[q]["gene"],
                "target": t,
                "species": sp[t],
                "role": tmod[t]["role"],
                "in_pf21671": "yes" if t in pf else "no",
                "model": tmod[t]["status"],
                "mean_plddt": tmod[t]["mean_plddt"],
                "residues_kept": tmod[t]["residues_kept"],
                "prostt5_e": fmt(s["prostt5"]),
                "af_fs_e": fmt(s["af_fs"]),
                "reseek_e": fmt(s["reseek"]),
                "mmseqs_e": fmt(mm.get((q, t))),
                **{f"support_{m}": "yes" if sup[m] else "no" for m in METHODS},
                "n_support": sum(sup.values()),
            }
        )

    def best_scaled(p):
        return min(
            (float(p[c]) for c in ("prostt5_e", "af_fs_e", "reseek_e") if p[c] != "NA"), default=9.0
        )

    targets = []
    for t in sorted({p["target"] for p in pairs}):
        tp = [p for p in pairs if p["target"] == t]
        best = max(tp, key=lambda p: (p["n_support"], -best_scaled(p)))
        f = flags.get(t, {})
        targets.append(
            {
                "target": t,
                "species": sp[t],
                "role": tmod[t]["role"],
                "in_pf21671": "yes" if t in pf else "no",
                "panels": ",".join(sorted({p["panel"] for p in tp})),
                "model": tmod[t]["status"],
                "mean_plddt": tmod[t]["mean_plddt"],
                "sp_yes": f.get("sp_yes", "NA"),
                "gpi_yes": f.get("gpi_yes", "NA"),
                "disorder_longest_run": f.get("disorder_longest_run", "NA"),
                "queries_2of3": sum(1 for p in tp if p["n_support"] >= 2),
                "confirmed": "yes" if any(p["n_support"] >= 2 for p in tp) else "no",
                "best_query": best["query_id"],
                "best_prostt5_e": best["prostt5_e"],
                "best_af_fs_e": best["af_fs_e"],
                "best_reseek_e": best["reseek_e"],
            }
        )
    bg_hits = []
    for m in ("af_fs", "reseek"):
        for (q, t), s in sorted(tabs[m].items()):
            if t in half["A"] | half["B"] and s <= CUT:
                bg_hits.append(
                    {
                        "method": m,
                        "query_id": q,
                        "target": t,
                        "species": sp[t],
                        "bg_half": tmod[t]["bg_half"],
                        "role": tmod[t]["role"],
                        "scaled_e": fmt(s),
                    }
                )
    b_support = {}  # (q, b) -> sup for half B
    for q in dom:
        for b in half["B"]:
            if any((q, b) in tabs[m] for m in METHODS):
                b_support[(q, b)] = pair(q, b)[1]
    summary = []
    species = sorted({sp[i] for i in tmod})
    for panel in sorted({d["panel"] for d in domains}):
        qs = {q for q, d in dom.items() if d["panel"] == panel}
        pp = [p for p in pairs if p["query_id"] in qs]
        qmodel_missing = sum(1 for q in qs if mod.get(q, {}).get("status") != "ok")
        for m in (*METHODS, "2of3"):
            if m == "2of3":
                hit = {p["target"] for p in pp if p["n_support"] >= 2}
                bpass = {
                    b for (q, b), sup in b_support.items() if q in qs and sum(sup.values()) >= 2
                }
            else:
                hit = {p["target"] for p in pp if p[f"support_{m}"] == "yes"}
                bpass = {b for (q, b), sup in b_support.items() if q in qs and sup[m]}
            row = {"panel": panel, "method": m, "targets_supported": len(hit)}
            for s in species:
                members = {a for a in pf if sp.get(a) == s}
                row[f"pf21671_supported_{s}"] = f"{len(hit & members)}/{len(members)}"
                row[f"pf21671_no_model_{s}"] = len(members - ok)
                row[f"non_pf21671_supported_{s}"] = len({t for t in hit if sp[t] == s} - members)
            row["half_B_passing"] = f"{len(bpass)}/{len(half['B'])}"
            row["queries_without_model"] = "NA" if m == "prostt5" else f"{qmodel_missing}/{len(qs)}"
            summary.append(row)
    null_rows = []
    for q in sorted(dom):
        for m in METHODS:
            v = a_scores.get((q, m), [])
            null_rows.append(
                {
                    "query_id": q,
                    "method": m,
                    "query_model": mod.get(q, {}).get("status", "NA"),
                    "null_n": 51 if m == "prostt5" else f"{len(v)}/{len(half['A'])}",
                    "null_best": fmt(decoy.get(q) if m == "prostt5" else (v[0][0] if v else None)),
                }
            )
    return pairs, targets, null_rows, bg_hits, summary


def write(path, rows, header=None):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=list(rows[0]) if rows else (header or []),
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for opt in ("--round1-inputs", "--round1-summary", "--inputs2", "--results2", "--out"):
        ap.add_argument(opt, required=True)
    a = ap.parse_args(argv)
    try:
        out = Path(a.out)
        if out.exists():
            raise SummaryError(f"{out} exists; nothing is overwritten")
        r1, s1, i2, r2 = map(Path, (a.round1_inputs, a.round1_summary, a.inputs2, a.results2))
        domains = read_tsv(r1 / "query_domains.tsv")
        truth = read_tsv(r1 / "truth_pf21671.tsv")
        n1 = len(read_tsv(r1 / "target_proteins.tsv"))
        models = read_tsv(i2 / "models.tsv")
        known = {m["id"] for m in models}
        n_db = sum(1 for m in models if m["role"] != "query" and m["status"] == "ok")
        if n_db == 0:
            raise SummaryError("no usable target models")
        af = read_scores(r2 / "foldseek_af.tsv", FS_COLS, "evalue", known, n1 / n_db)
        rs = read_scores(r2 / "reseek.tsv", RS_COLS, "pvalue", known, n1)
        pairs, targets, nulls, bg_hits, summary = summarize(
            domains,
            models,
            read_tsv(s1 / "hits.tsv"),
            read_tsv(s1 / "per_query.tsv"),
            af,
            rs,
            truth,
        )
        out.mkdir(parents=True)
        write(out / "pairs.tsv", pairs)
        write(out / "targets.tsv", targets)
        write(out / "nulls.tsv", nulls)
        write(
            out / "background_hits.tsv",
            bg_hits,
            ["method", "query_id", "target", "species", "bg_half", "role", "scaled_e"],
        )
        write(out / "summary.tsv", summary)
        (out / "scales.txt").write_text(
            f"N1 (round 1 targets): {n1}\nN_db (target models searched): {n_db}\n"
            f"af_fs scale: {n1 / n_db:.4g}\nreseek scale: {n1}\n"
        )
        print(
            f"wrote {out}: {len(pairs)} pairs, {sum(t['confirmed'] == 'yes' for t in targets)} confirmed targets"
        )
    except SummaryError as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
