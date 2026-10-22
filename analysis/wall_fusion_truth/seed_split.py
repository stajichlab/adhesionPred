#!/usr/bin/env python3
"""Label MMseqs2 clusters `seed-like` or `seed-independent` (plan Task 5.3; spec D8, D20).

DRAFT, not run. Stdlib only. Run with `/usr/bin/python3.12 -I` on HPCC or `python3 -I` locally.

A cluster is `seed-like` when ANY member has a hit with at least 30% identity to a Pfam seed
sequence (target coverage 0.5 or more, E <= 1e-3). This is the rule of
`analysis/lectin_validation/finalize_truth.py` and `measure.py` (the `SPLIT_SIM` line). Otherwise
the cluster is `seed-independent`. `seed-independent` does NOT exclude remote homology that
MMseqs2 at 30% misses.

Two classes are written (owner decision 2026-10-22):
  seed_class_slice  the models of the row's own slice: wall and fusion rows against the 21 wall and
                    fusion models (analysis/wall_fusion_report/wall_fusion_models.tsv); adhesion rows
                    against the 18 active family rows of data/sorting_hat/family_table.tsv. Negatives
                    and family representatives have no slice and do not count here: a cluster of only
                    negatives is `not_applicable`. A mixed cluster (control and negative) is judged by
                    its control members.
  seed_class_all39  every member against all 39 models.
A negative is judged by all 39 for `leakage_suggested`; a control by its slice models.

Inputs (all made by cluster_split.sh, or by hand for a test):
  --members   TSV with columns accession, class (positive | negative | family_rep), slice and
              optionally more; one line per sequence that was clustered
  --clusters  MMseqs2 easy-cluster output `*_cluster.tsv` (representative<TAB>member)
  --seed-hits MMseqs2 easy-search output, columns
              query,target,pident,alnlen,qcov,tcov,evalue,bits
              (target = `<PFxxxxx>|<seq name>`). Identity may be a percentage or a fraction; a
              file whose largest value is 1.0 or less is read as a fraction (as finalize_truth.py does).
  --controls  optional: a control-set directory (data/controls/wall-fusion-kiss2019) whose
              controls.tsv and controls_local_sequence.tsv (BOTH are the set) give the clade and
              stratum of each cluster; family_reps.tsv gives the flagged group of representatives
  --wall-models, --family-table  the two model lists (default: the repository files)

Outputs in --out:
  members_split.tsv   accession, cluster_id, class, slice, overlap_with, in_mixed_cluster, both seed classes,
                      the best identity and nearest model for each, leakage_suggested (in_reference | none)
  clusters_split.tsv  one line per cluster: members, classes, slices, mixed_pos_neg (control and negative
                      in one cluster: filter on it), has_overlap_copy, both seed classes, clade and stratum columns
  clusters.tsv        accession<TAB>cluster_id (the control-set file of COMMON-RULES)
  summary.tsv         clusters per class and slice for both seed classes; never pooled; the family
                      representatives (not clustered) are a flagged group when --controls is given

The negatives are whatever the members file holds (the owner chose the shared set of 161 rows,
data/controls/hard-negatives-shared). A negative that is also a Kiss control is written as
`<accession>#neg` and has `overlap_with` set; its cluster is a mixed cluster.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WALL_MODELS = ROOT / "analysis/wall_fusion_report/wall_fusion_models.tsv"
FAMILY_TABLE = ROOT / "data/sorting_hat/family_table.tsv"
SEED_IDENTITY = 30.0
MIN_TCOV = 0.5
MAX_EVALUE = 1e-3
CLASSES = ("positive", "negative", "family_rep")
HIT_COLUMNS = ("query", "target", "pident", "alnlen", "qcov", "tcov", "evalue", "bits")
SLICE_SET = {"adhesion": "family", "wall": "wall_fusion", "fusion": "wall_fusion"}


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def read_model_sets(wall: Path = WALL_MODELS, family: Path = FAMILY_TABLE) -> dict[str, set[str]]:
    """{'wall_fusion': the 21 wall and fusion models, 'family': the active adhesion family rows}."""
    return {
        "wall_fusion": {r["id"] for r in read_tsv(wall)},
        "family": {
            r["pfam_acc"] for r in read_tsv(family) if r["active"] == "yes" and r["pfam_acc"]
        },
    }


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
        if r["slice"] and r["slice"] not in SLICE_SET:
            raise ValueError(
                f"{r['accession']}: slice {r['slice']!r} is not one of {sorted(SLICE_SET)}"
            )
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


def best_hits(hits: list[dict], model_sets: dict[str, set[str]]) -> dict[str, dict]:
    """{query: {'all39': (identity, model) or None, 'wall_fusion': ..., 'family': ...}}.

    Only a hit that meets identity, target coverage and E-value together counts. A model that is in
    neither model set is an error (the seed FASTA and the model lists must agree).
    """
    known = model_sets["wall_fusion"] | model_sets["family"]
    best: dict[str, dict] = {}
    for h in hits:
        ok = h["pident"] >= SEED_IDENTITY and h["tcov"] >= MIN_TCOV and h["evalue"] <= MAX_EVALUE
        b = best.setdefault(h["query"], {"all39": None, "wall_fusion": None, "family": None})
        if not ok:
            continue
        model = h["target"].split("|")[0]
        if model not in known:
            raise ValueError(f"seed hit to {model}, which is in neither model list")
        for key, member in (
            ("all39", True),
            ("wall_fusion", model in model_sets["wall_fusion"]),
            ("family", model in model_sets["family"]),
        ):
            if member and (b[key] is None or h["pident"] > b[key][0]):
                b[key] = (h["pident"], model)
    return best


def _fmt(hit) -> tuple[str, str]:
    return (f"{hit[0]:.1f}", hit[1]) if hit else ("", "")


def split(
    members: dict[str, dict],
    clusters: dict[str, list[str]],
    hits: list[dict],
    model_sets: dict[str, set[str]],
):
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
    best = best_hits(hits, model_sets)
    unknown = sorted(set(best) - set(members))
    if unknown:
        raise ValueError(
            f"seed hits for {len(unknown)} sequences that are not members, e.g. {unknown[0]}"
        )

    def slice_hit(m: str):
        """The best hit against the models of the member's own slice (None: no hit or no slice)."""
        key = SLICE_SET.get(members[m]["slice"])
        return best.get(m, {}).get(key) if key else None

    ordered = sorted(clusters.values(), key=lambda ms: min(ms))
    member_rows, cluster_rows = [], []
    for i, ms in enumerate(ordered, 1):
        cid = f"KS{i:04d}"
        with_slice = [m for m in ms if SLICE_SET.get(members[m]["slice"])]
        hit_slice = [m for m in with_slice if slice_hit(m)]
        hit_all = [m for m in ms if best.get(m, {}).get("all39")]
        class_slice = (
            "not_applicable"
            if not with_slice
            else ("seed-like" if hit_slice else "seed-independent")
        )
        class_all = "seed-like" if hit_all else "seed-independent"
        classes = sorted({members[m]["class"] for m in ms})
        slices = sorted({members[m]["slice"] for m in ms if members[m]["slice"]})
        mixed = "yes" if {"positive", "negative"} <= set(classes) else "no"
        for m in sorted(ms):
            sh, ah = slice_hit(m), best.get(m, {}).get("all39")
            id_s, mod_s = _fmt(sh)
            id_a, mod_a = _fmt(ah)
            # a control is judged by its slice models; a negative or representative by all 39
            leak = sh if SLICE_SET.get(members[m]["slice"]) else ah
            member_rows.append(
                {
                    "accession": m,
                    "cluster_id": cid,
                    "class": members[m]["class"],
                    "slice": members[m]["slice"],
                    "overlap_with": members[m].get("overlap_with", ""),
                    "in_mixed_cluster": mixed,
                    "seed_class_slice": class_slice,
                    "seed_class_all39": class_all,
                    "best_identity_slice": id_s,
                    "nearest_model_slice": mod_s,
                    "best_identity_all39": id_a,
                    "nearest_model_all39": mod_a,
                    "member_hit_slice": "yes" if sh else "no",
                    "member_hit_all39": "yes" if ah else "no",
                    "leakage_suggested": "in_reference" if leak else "none",
                }
            )
        cluster_rows.append(
            {
                "cluster_id": cid,
                "n_members": len(ms),
                "classes": ",".join(classes),
                "slices": ",".join(slices),
                "mixed_pos_neg": mixed,
                "has_overlap_copy": "yes"
                if any(members[m].get("overlap_with") for m in ms)
                else "no",
                "seed_class_slice": class_slice,
                "seed_class_all39": class_all,
                "n_members_with_hit_slice": len(hit_slice),
                "n_members_with_hit_all39": len(hit_all),
                "best_identity_slice": (
                    f"{max(slice_hit(m)[0] for m in hit_slice):.1f}" if hit_slice else ""
                ),
                "best_identity_all39": (
                    f"{max(best[m]['all39'][0] for m in hit_all):.1f}" if hit_all else ""
                ),
                "members": ",".join(sorted(ms)),
            }
        )
    return member_rows, cluster_rows


def join_controls(
    cluster_rows: list[dict], member_rows: list[dict], controls_dir: Path
) -> list[str]:
    """Add clade and stratum columns from a control-set directory; return the new column names.

    The set is two files, controls.tsv and controls_local_sequence.tsv: both are read.
    """
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


def summary(
    member_rows: list[dict], cluster_rows: list[dict], controls_dir: Path | None
) -> list[dict]:
    """Clusters per (class, slice) and seed class. A cluster counts in each slice it touches.

    Nothing is pooled across slices (D20). The family representatives are not clustered; with a
    control-set directory they are listed as a flagged group (distinct Kiss clusters, no seed class).
    """
    by_id = {c["cluster_id"]: c for c in cluster_rows}
    touched: dict[str, set] = defaultdict(set)
    for m in member_rows:
        touched[m["cluster_id"]].add((m["class"], m["slice"]))
    count: Counter = Counter()
    for cid, pairs in touched.items():
        for cl, sl in pairs:
            count[(cl, sl, by_id[cid]["seed_class_slice"], by_id[cid]["seed_class_all39"])] += 1
    out = []
    for cl, sl in sorted({(k[0], k[1]) for k in count}):
        keys = [k for k in count if k[:2] == (cl, sl)]

        def n(pos, value, keys=keys):
            return sum(count[k] for k in keys if k[pos] == value)

        out.append(
            {
                "class": cl,
                "slice": sl,
                "clusters": sum(count[k] for k in keys),
                "slice_seed_like": n(2, "seed-like"),
                "slice_seed_independent": n(2, "seed-independent"),
                "all39_seed_like": n(3, "seed-like"),
                "all39_seed_independent": n(3, "seed-independent"),
            }
        )
    if controls_dir is not None and (controls_dir / "family_reps.tsv").is_file():
        reps = defaultdict(set)
        for r in read_tsv(controls_dir / "family_reps.tsv"):
            reps[r["slice"]].add(r["kiss_cluster_number"])
        for sl in sorted(reps):
            out.append(
                {
                    "class": "family_rep (not clustered)",
                    "slice": sl,
                    "clusters": len(reps[sl]),
                    "slice_seed_like": "",
                    "slice_seed_independent": "",
                    "all39_seed_like": "",
                    "all39_seed_independent": "",
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
    ap.add_argument("--wall-models", type=Path, default=WALL_MODELS)
    ap.add_argument("--family-table", type=Path, default=FAMILY_TABLE)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    members = read_members(a.members)
    model_sets = read_model_sets(a.wall_models, a.family_table)
    member_rows, cluster_rows = split(
        members, read_clusters(a.clusters), read_hits(a.seed_hits), model_sets
    )
    extra = join_controls(cluster_rows, member_rows, a.controls) if a.controls else []
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "members_split.tsv", list(member_rows[0]), member_rows)
    write_tsv(
        a.out / "clusters_split.tsv",
        [
            "cluster_id",
            "n_members",
            "classes",
            "slices",
            "mixed_pos_neg",
            "has_overlap_copy",
            "seed_class_slice",
            "seed_class_all39",
            "n_members_with_hit_slice",
            "n_members_with_hit_all39",
            "best_identity_slice",
            "best_identity_all39",
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
        [
            "class",
            "slice",
            "clusters",
            "slice_seed_like",
            "slice_seed_independent",
            "all39_seed_like",
            "all39_seed_independent",
        ],
        summary(member_rows, cluster_rows, a.controls),
    )
    like = sum(1 for c in cluster_rows if c["seed_class_all39"] == "seed-like")
    mixed = sum(1 for c in cluster_rows if c["mixed_pos_neg"] == "yes")
    print(
        f"{len(member_rows)} sequences in {len(cluster_rows)} clusters: {like} seed-like (all 39), "
        f"{len(cluster_rows) - like} seed-independent; {mixed} mixed (control and negative)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
