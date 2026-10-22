#!/usr/bin/env python3
"""Label MMseqs2 clusters `seed-like` or `seed-independent` (plan Task 5.3; spec D8, D20).

DRAFT, not run. Stdlib only. Run with `/usr/bin/python3.12 -I` on HPCC or `python3 -I` locally.

A cluster is `seed-like` when ANY member has a hit with at least 30% identity to a Pfam seed
sequence of one of the models (target coverage 0.5 or more, E <= 1e-3). This is the rule of
`analysis/lectin_validation/finalize_truth.py` and `measure.py` (the `SPLIT_SIM` line). Otherwise
the cluster is `seed-independent`. `seed-independent` does NOT exclude remote homology that
MMseqs2 at 30% misses.

Inputs (all made by cluster_split.sh, or by hand for a test):
  --members   TSV with columns accession, class (positive | negative | family_rep), slice and
              optionally more; one line per sequence that was clustered
  --clusters  MMseqs2 easy-cluster output `*_cluster.tsv` (representative<TAB>member)
  --seed-hits MMseqs2 easy-search output, columns
              query,target,pident,alnlen,qcov,tcov,evalue,bits
              (target = `<PFxxxxx>|<seq name>`). Identity may be a percentage or a fraction; a
              file whose largest value is 1.0 or less is read as a fraction (as finalize_truth.py does).
  --controls  optional: a control-set directory (data/controls/wall-fusion-kiss2019) whose
              controls.tsv and controls_local_sequence.tsv give the clade and stratum of each cluster

Outputs in --out:
  members_split.tsv   accession, cluster_id, class, slice, seed_class, best_identity, nearest_model,
                      member_seed_hit, leakage_suggested (in_reference | none)
  clusters_split.tsv  one line per cluster: members, classes, slices, seed_class, clade and stratum columns
  clusters.tsv        accession<TAB>cluster_id (the control-set file of COMMON-RULES)
  summary.tsv         clusters per class and slice, seed-like and seed-independent; never pooled

The two negative sources (data/curated/adhesins/hard_negative_seeds.tsv and the shared set of 161
rows) are NOT chosen here: the members file says which sequences were clustered.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

SEED_IDENTITY = 30.0
MIN_TCOV = 0.5
MAX_EVALUE = 1e-3
CLASSES = ("positive", "negative", "family_rep")
HIT_COLUMNS = ("query", "target", "pident", "alnlen", "qcov", "tcov", "evalue", "bits")


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def read_members(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for r in read_tsv(path):
        for col in ("accession", "class", "slice"):
            if col not in r:
                raise ValueError(f"{path.name}: missing column {col}")
        if r["class"] not in CLASSES:
            raise ValueError(f"{r['accession']}: class {r['class']!r} is not one of {CLASSES}")
        if r["accession"] in out:
            raise ValueError(f"{r['accession']}: listed twice in {path.name}")
        out[r["accession"]] = r
    return out


def read_clusters(path: Path) -> dict[str, list[str]]:
    """{representative: [members]} from the two-column MMseqs2 table."""
    out: dict[str, list[str]] = defaultdict(list)
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 2 and f[0] and f[1]:
            out[f[0]].append(f[1])
    return dict(out)


def read_hits(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < len(HIT_COLUMNS) or f[0] == "query":
            continue
        rec = dict(zip(HIT_COLUMNS, f, strict=False))
        rec["pident"], rec["tcov"], rec["evalue"] = float(f[2]), float(f[5]), float(f[6])
        rows.append(rec)
    if rows and max(r["pident"] for r in rows) <= 1.0:  # a fraction, not a percentage
        for r in rows:
            r["pident"] *= 100.0
    return rows


def best_hits(hits: list[dict]) -> dict[str, dict]:
    """{query: {identity, model, qualifying, any_identity}} for the best hit of each query.

    `qualifying` is the best hit that meets identity, target coverage and E-value together.
    `any_identity` is the best identity of any reported hit (for the report only).
    """
    best: dict[str, dict] = {}
    for h in hits:
        b = best.setdefault(
            h["query"], {"identity": None, "model": "", "any_identity": h["pident"]}
        )
        b["any_identity"] = max(b["any_identity"], h["pident"])
        ok = h["pident"] >= SEED_IDENTITY and h["tcov"] >= MIN_TCOV and h["evalue"] <= MAX_EVALUE
        if ok and (b["identity"] is None or h["pident"] > b["identity"]):
            b["identity"], b["model"] = h["pident"], h["target"].split("|")[0]
    return best


def split(members: dict[str, dict], clusters: dict[str, list[str]], hits: list[dict]):
    """Return (member rows, cluster rows) sorted deterministically."""
    in_clusters = [m for ms in clusters.values() for m in ms]
    if len(in_clusters) != len(set(in_clusters)):
        raise ValueError("a sequence is in two clusters")
    if set(in_clusters) != set(members):
        missing = sorted(set(members) - set(in_clusters))
        extra = sorted(set(in_clusters) - set(members))
        raise ValueError(
            f"clusters and members differ: {len(missing)} members not clustered "
            f"(e.g. {missing[:1]}), {len(extra)} clustered sequences not in members (e.g. {extra[:1]})"
        )
    best = best_hits(hits)
    unknown = sorted(set(best) - set(members))
    if unknown:
        raise ValueError(
            f"seed hits for {len(unknown)} sequences that are not members, e.g. {unknown[0]}"
        )

    ordered = sorted(clusters.values(), key=lambda ms: min(ms))
    member_rows, cluster_rows = [], []
    for i, ms in enumerate(ordered, 1):
        cid = f"KS{i:04d}"
        hit_members = [m for m in ms if best.get(m, {}).get("identity") is not None]
        seed_class = "seed-like" if hit_members else "seed-independent"
        classes = sorted({members[m]["class"] for m in ms})
        slices = sorted({members[m]["slice"] for m in ms if members[m]["slice"]})
        for m in sorted(ms):
            b = best.get(m, {})
            qual = b.get("identity") is not None
            member_rows.append(
                {
                    "accession": m,
                    "cluster_id": cid,
                    "class": members[m]["class"],
                    "slice": members[m]["slice"],
                    "seed_class": seed_class,
                    "best_identity": f"{b['identity']:.1f}" if qual else "",
                    "nearest_model": b.get("model", "") if qual else "",
                    "member_seed_hit": "yes" if qual else "no",
                    "leakage_suggested": "in_reference" if qual else "none",
                }
            )
        cluster_rows.append(
            {
                "cluster_id": cid,
                "n_members": len(ms),
                "classes": ",".join(classes),
                "slices": ",".join(slices),
                "mixed_pos_neg": "yes" if {"positive", "negative"} <= set(classes) else "no",
                "seed_class": seed_class,
                "n_members_with_seed_hit": len(hit_members),
                "best_identity": (
                    f"{max(best[m]['identity'] for m in hit_members):.1f}" if hit_members else ""
                ),
                "members": ",".join(sorted(ms)),
            }
        )
    return member_rows, cluster_rows


def join_controls(
    cluster_rows: list[dict], member_rows: list[dict], controls_dir: Path
) -> list[str]:
    """Add clade and stratum columns from a control-set directory; return the new column names."""
    meta: dict[str, dict] = {}
    for name in ("controls.tsv", "controls_local_sequence.tsv"):
        p = controls_dir / name
        if p.is_file():
            for r in read_tsv(p):
                meta[r["accession"]] = r
    cols = [
        "clades_fine",
        "clades_coarse",
        "strata_fine",
        "strata_coarse",
        "kiss_clusters",
        "jgi_accept",
    ]
    by_cluster: dict[str, list[str]] = defaultdict(list)
    for m in member_rows:
        by_cluster[m["cluster_id"]].append(m["accession"])
    src = {
        "clades_fine": "clade_fine",
        "clades_coarse": "clade_coarse",
        "strata_fine": "cluster_stratum_fine",
        "strata_coarse": "cluster_stratum_coarse",
        "kiss_clusters": "kiss_cluster_number",
        "jgi_accept": "jgi_accept",
    }
    for c in cluster_rows:
        rows = [meta[a] for a in by_cluster[c["cluster_id"]] if a in meta]
        for col in cols:
            c[col] = ",".join(sorted({r[src[col]] for r in rows if r.get(src[col])}))
    return cols


def summary(member_rows: list[dict], cluster_rows: list[dict]) -> list[dict]:
    """Clusters per (class, slice) and seed class. A cluster counts in each slice it touches.

    Nothing is pooled across slices (D20). Family representatives are listed apart.
    """
    cls_of = {m["cluster_id"]: set() for m in member_rows}
    for m in member_rows:
        cls_of[m["cluster_id"]].add((m["class"], m["slice"]))
    count: Counter = Counter()
    for c in cluster_rows:
        for cl, sl in cls_of[c["cluster_id"]]:
            count[(cl, sl, c["seed_class"])] += 1
    out = []
    for cl, sl in sorted({(k[0], k[1]) for k in count}):
        out.append(
            {
                "class": cl,
                "slice": sl,
                "clusters_seed_like": count[(cl, sl, "seed-like")],
                "clusters_seed_independent": count[(cl, sl, "seed-independent")],
                "clusters": count[(cl, sl, "seed-like")] + count[(cl, sl, "seed-independent")],
            }
        )
    return out


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(c, "") for c in header])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--members", type=Path, required=True)
    ap.add_argument("--clusters", type=Path, required=True)
    ap.add_argument("--seed-hits", type=Path, required=True)
    ap.add_argument("--controls", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    members = read_members(a.members)
    member_rows, cluster_rows = split(members, read_clusters(a.clusters), read_hits(a.seed_hits))
    extra = join_controls(cluster_rows, member_rows, a.controls) if a.controls else []
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(
        a.out / "members_split.tsv",
        list(member_rows[0]),
        member_rows,
    )
    write_tsv(
        a.out / "clusters_split.tsv",
        [
            "cluster_id",
            "n_members",
            "classes",
            "slices",
            "mixed_pos_neg",
            "seed_class",
            "n_members_with_seed_hit",
            "best_identity",
            *extra,
            "members",
        ],
        cluster_rows,
    )
    (a.out / "clusters.tsv").write_text(
        "accession\tcluster_id\n"
        + "".join(f"{m['accession']}\t{m['cluster_id']}\n" for m in member_rows)
    )
    write_tsv(
        a.out / "summary.tsv",
        ["class", "slice", "clusters", "clusters_seed_like", "clusters_seed_independent"],
        summary(member_rows, cluster_rows),
    )
    like = sum(1 for c in cluster_rows if c["seed_class"] == "seed-like")
    print(
        f"{len(member_rows)} sequences in {len(cluster_rows)} clusters: {like} seed-like, "
        f"{len(cluster_rows) - like} seed-independent",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
