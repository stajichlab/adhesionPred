#!/usr/bin/env python3
"""PR 6 for hyd1_like: the per-unit training manifest, the cluster units, the task list and the cutoff (DRAFT, not run on HPCC).

Stdlib only; run with `python3 -I`. It writes no call, no status, no tier and no sequence. Spec: docs/superpowers/specs/2026-10-23-pr6-hyd1-like-design.md.

Reuse, not copy: the independence helpers, the alignment checks, the cutoff rule and the freeze guard are IMPORTED from pr6_ecm33_build.py (loaded
by path; `python3 -I` does not put this folder on sys.path). The ecm33 files are not edited. What is new here is only (a) units that are a clade OR a
30 percent MMseqs2 cluster (`clade:<name>`, `cluster:<representative>`, `all`) and (b) the two readings of the training rule.

Subcommands
  cluster-units  units of the leave-cluster-out kind from the MMseqs2 `_cluster.tsv` over the members (after pr6_hyd1_extract.sbatch)
  manifest       the training set of ONE unit after the independence rule; the held-out ids are never in it
  tasks          the (unit, mode) task list of the build array, from units.tsv and units_cluster.tsv; refuses an unfrozen freeze
  check-alignment, cutoff   as in pr6_ecm33_build.py (same code); the cutoff refuses a missing or unfrozen freeze

For a `cluster:` unit the held-out members are the whole 30 percent cluster, so `cluster30` removes no further training member and `none` is the
same set: ONE column only (the task list gives cluster units the headline mode only).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


e33 = _load("pr6_ecm33_build")
mem = _load("pr6_hyd1_members")

# the same helpers, under the same names (the tests compare them with the ecm33 ones)
sha256_file = e33.sha256_file
read_tsv = e33.read_tsv
fasta_id = e33.fasta_id
read_cluster_tsv = e33.read_cluster_tsv
check_alignment_ids = e33.check_alignment_ids
occupancy_check = e33.occupancy_check
lowest_cutoff = e33.lowest_cutoff
run_cutoff = e33.run_cutoff
extra_call_limit_total = e33.extra_call_limit_total

EXTRA_FREEZE_KEYS = ("training_rule", "comparators", "heldout_truth", "evalue_db_size")


def load_freeze(path: Path) -> dict:
    """ecm33's guard (missing, incomplete, unfrozen, modes) plus the hyd1 keys."""
    fz = e33.load_freeze(path)
    missing = [k for k in EXTRA_FREEZE_KEYS if k not in fz]
    if missing:
        raise SystemExit(f"FATAL: freeze file {path} lacks {missing}")
    if fz["training_rule"].get("reading") not in ("total", "per_clade"):
        raise SystemExit(
            f"FATAL: freeze file {path}: training_rule.reading must be total or per_clade"
        )
    return fz


def held_out_ids(members: list[dict], unit: str, clusters: dict[str, str] | None) -> set[str]:
    kind, _, name = unit.partition(":")
    if unit == "all":
        return set()
    if kind == "clade":
        held = {fasta_id(m) for m in members if m["clade_fine"] == name}
    elif kind == "cluster":
        if clusters is None:
            raise ValueError("a cluster unit needs the MMseqs2 cluster table")
        held = {fasta_id(m) for m in members if clusters.get(fasta_id(m)) == name}
    else:
        raise ValueError(f"unknown unit {unit}")
    if not held:
        raise ValueError(f"no members in unit {unit}")
    return held


def make_manifest(
    members: list[dict],
    unit: str,
    clusters: dict[str, str] | None,
    mode: str = "cluster30",
    reading: str = "total",
) -> dict:
    if mode not in ("cluster30", "none"):
        raise ValueError(f"unknown mode {mode}")
    held = held_out_ids(members, unit, clusters)
    train = [m for m in members if fasta_id(m) not in held]
    excluded = []
    if mode == "cluster30" and held:
        if clusters is None:
            raise ValueError("mode cluster30 needs the MMseqs2 cluster table")
        missing = sorted(i for i in held if i not in clusters)
        if missing:
            raise ValueError(f"held-out ids missing from the cluster table: {missing[:5]}")
        reps = {clusters[i] for i in held}
        excluded = [m for m in train if clusters.get(fasta_id(m)) in reps]
        gone = {fasta_id(m) for m in excluded}
        train = [m for m in train if fasta_id(m) not in gone]
    by: dict[str, int] = {}
    for m in train:
        by[m["clade_fine"]] = by.get(m["clade_fine"], 0) + 1
    total, per = mem.training_rules(by)
    rule = total if reading == "total" else per
    reason = ""
    if not rule:
        reason = (
            f"training set after the independence rule: {len(train)} sequences in clades {dict(sorted(by.items()))} "
            f"(rule {reading}: {mem.MIN_TRAIN_SEQ} sequences, {mem.MIN_TRAIN_CLADES} clades)"
        )
    return {
        "unit": unit,
        "mode": mode,
        "reading": reading,
        "held_out_ids": sorted(held),
        "training": [
            {
                "id": fasta_id(m),
                "clade_fine": m["clade_fine"],
                "cluster": (clusters or {}).get(fasta_id(m), ""),
            }
            for m in sorted(train, key=fasta_id)
        ],
        "excluded_30pct_neighbours": sorted(fasta_id(m) for m in excluded),
        "n_held_out": len(held),
        "n_training": len(train),
        "n_training_clades": len(by),
        "n_excluded_neighbours": len(excluded),
        "training_rule_met": rule,
        "status": "buildable" if rule else "not_buildable",
        "reason": reason,
    }


CLUSTER_UNIT_HEADER = ["unit", "cluster", "n_members", "clades", "has_query", "status"]


def cluster_units(
    members: list[dict],
    clusters: dict[str, str],
    reading: str = "total",
    min_unit: int = mem.MIN_UNIT_MEMBERS,
) -> list[dict]:
    """One row per MMseqs2 30 percent cluster of the members. status: unit (>= min_unit members, training rule met after removing the cluster),
    kiss_single (the cluster of the query), too_few_members, not_buildable. Most clusters will be small; that is a result, not a bug."""
    by: dict[str, list[dict]] = {}
    for m in members:
        by.setdefault(clusters[fasta_id(m)], []).append(m)
    rows = []
    for rep in sorted(by):
        grp = by[rep]
        ids = {fasta_id(m) for m in grp}
        rest: dict[str, int] = {}
        for m in members:
            if fasta_id(m) not in ids:
                rest[m["clade_fine"]] = rest.get(m["clade_fine"], 0) + 1
        total, per = mem.training_rules(rest)
        met = total if reading == "total" else per
        has_q = any(m["is_query"] == "yes" for m in grp)
        if not met:
            status = "not_buildable"
        elif has_q:
            status = (
                "kiss_single"
                if len(grp) == 1
                else "unit"
                if len(grp) >= min_unit
                else "too_few_members"
            )
        else:
            status = "unit" if len(grp) >= min_unit else "too_few_members"
        rows.append(
            {
                "unit": f"cluster:{rep}",
                "cluster": rep,
                "n_members": len(grp),
                "clades": ",".join(sorted({m["clade_fine"] for m in grp})),
                "has_query": "yes" if has_q else "no",
                "status": status,
            }
        )
    return rows


def task_list(
    units: list[dict], cluster_rows: list[dict], head_mode: str, second_mode: str
) -> list[tuple[str, str]]:
    """Clade units (and the Kiss unit): both modes; cluster units: the headline mode only; then `all` once."""
    out: list[tuple[str, str]] = []
    for u in units:
        if u["status"] in ("unit", "kiss_single"):
            out += [(u["unit"], head_mode), (u["unit"], second_mode)]
    kiss_clade = any(u["status"] == "kiss_single" for u in units)
    for u in cluster_rows:
        # a cluster that is the query alone is the same single-protein unit as the query's clade: not run twice
        if u["status"] == "unit" or (u["status"] == "kiss_single" and not kiss_clade):
            out.append((u["unit"], head_mode))
    out.append(("all", "none"))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    cu = sub.add_parser("cluster-units")
    cu.add_argument("--members", required=True, type=Path)
    cu.add_argument("--clusters", required=True, type=Path)
    cu.add_argument("--reading", default="total", choices=["total", "per_clade"])
    cu.add_argument("--out", required=True, type=Path)
    m = sub.add_parser("manifest")
    m.add_argument("--members", required=True, type=Path)
    m.add_argument("--unit", required=True, help="clade:<name>, cluster:<representative> or all")
    m.add_argument("--clusters", type=Path)
    m.add_argument("--mode", required=True, choices=["cluster30", "none"])
    m.add_argument("--freeze", required=True, type=Path)
    m.add_argument("--out", required=True, type=Path)
    t = sub.add_parser("tasks")
    t.add_argument("--units", required=True, type=Path)
    t.add_argument("--units-cluster", required=True, type=Path)
    t.add_argument("--freeze", required=True, type=Path)
    c = sub.add_parser("check-alignment")
    c.add_argument("--manifest", required=True, type=Path)
    c.add_argument("--alignment", required=True, type=Path)
    c.add_argument("--freeze", type=Path)
    k = sub.add_parser("cutoff")
    for flag in ("proteome-hits", "labelled", "neg-groups", "neg-hits", "sizes"):
        k.add_argument(f"--{flag}", required=True, type=Path)
    k.add_argument("--freeze", required=True, type=Path)
    k.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)

    if a.cmd == "cluster-units":
        members = read_tsv(a.members)
        rows = cluster_units(members, read_cluster_tsv(a.clusters), a.reading)
        mem.write_tsv(a.out, CLUSTER_UNIT_HEADER, rows)
        print(f"{len(rows)} clusters; statuses {sorted({r['status'] for r in rows})}")
        return 0
    if a.cmd == "tasks":
        fz = load_freeze(a.freeze)
        for unit, mode in task_list(
            read_tsv(a.units),
            read_tsv(a.units_cluster),
            fz["independence_mode"],
            fz["independence_second_column"],
        ):
            print(f"{unit}\t{mode}")
        return 0
    if a.cmd == "manifest":
        fz = load_freeze(a.freeze)  # no manifest, so no build, before the freeze
        members = read_tsv(a.members)
        clusters = read_cluster_tsv(a.clusters) if a.clusters else None
        man = make_manifest(members, a.unit, clusters, a.mode, fz["training_rule"]["reading"])
        man["members_sha256"] = sha256_file(a.members)
        a.out.mkdir(parents=True, exist_ok=True)
        (a.out / "manifest.json").write_text(json.dumps(man, indent=1, sort_keys=True) + "\n")
        (a.out / "train_ids.txt").write_text("".join(t["id"] + "\n" for t in man["training"]))
        (a.out / "heldout_ids.txt").write_text("".join(i + "\n" for i in man["held_out_ids"]))
        print(
            f"unit {a.unit} ({a.mode}): held out {man['n_held_out']}, training {man['n_training']} in "
            f"{man['n_training_clades']} clades, excluded 30% neighbours {man['n_excluded_neighbours']}: {man['status']} {man['reason']}"
        )
        return 0 if man["status"] == "buildable" else 3
    # check-alignment and cutoff: the ecm33 code, with the hyd1 freeze guard in front of the cutoff
    if a.cmd == "check-alignment":
        argv2 = ["check-alignment", "--manifest", str(a.manifest), "--alignment", str(a.alignment)]
        if a.freeze:
            load_freeze(a.freeze)
            argv2 += ["--freeze", str(a.freeze)]
        return e33.main(argv2)
    load_freeze(a.freeze)
    return e33.main(
        [
            "cutoff",
            "--freeze",
            str(a.freeze),
            "--proteome-hits",
            str(a.proteome_hits),
            "--labelled",
            str(a.labelled),
            "--neg-groups",
            str(a.neg_groups),
            "--neg-hits",
            str(a.neg_hits),
            "--sizes",
            str(a.sizes),
            "--out",
            str(a.out),
        ]
    )


if __name__ == "__main__":
    sys.exit(main())
