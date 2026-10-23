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
  check-alignment   as in pr6_ecm33_build.py (same code)
  cutoff            the ecm33 rule (same `run_cutoff`) with `panel_ignored` (owner decision Q8) left out of the extra calls and the hard-negative rate AND the
                    same numbers with those ids counted; refuses a missing or unfrozen freeze
  unit `external`   the externally curated proteins of the freeze (heldout_truth.external_set), one column (cluster30), descriptive

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

EXTRA_FREEZE_KEYS = (
    "training_rule",
    "comparators",
    "heldout_truth",
    "evalue_db_size",
    "proteome_other_hydrophobin",
    "panel_ignored",
)


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
    if fz["proteome_other_hydrophobin"] not in ("ignore", "count"):
        raise SystemExit(
            f"FATAL: freeze file {path}: proteome_other_hydrophobin must be ignore or count"
        )
    if not isinstance(fz["panel_ignored"], list):
        raise SystemExit(f"FATAL: freeze file {path}: panel_ignored must be a list of ids")
    ext = fz["heldout_truth"].get("external_set")
    if ext is not None and (not isinstance(ext, list) or len(set(ext)) != len(ext)):
        raise SystemExit(
            f"FATAL: freeze file {path}: heldout_truth.external_set must be null or a list of distinct ids"
        )
    return fz


def external_ids_of(fz: dict) -> list[str]:
    return list(fz["heldout_truth"].get("external_set") or [])


def held_out_ids(
    members: list[dict],
    unit: str,
    clusters: dict[str, str] | None,
    external_ids: list[str] | None = None,
) -> set[str]:
    kind, _, name = unit.partition(":")
    if unit == "all":
        return set()
    if unit == "external":
        # the externally curated proteins; a training member is never a test protein (the owner's rule for A0A6V8R0V1)
        held = set(external_ids or [])
        if not held:
            raise ValueError("unit external needs freeze heldout_truth.external_set")
        both = held & {fasta_id(m) for m in members}
        if both:
            raise ValueError(f"external ids that are training members: {sorted(both)}")
        return held
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
    external_ids: list[str] | None = None,
) -> dict:
    if mode not in ("cluster30", "none"):
        raise ValueError(f"unknown mode {mode}")
    if unit == "external" and mode != "cluster30":
        raise ValueError("unit external has one column: cluster30")
    held = held_out_ids(members, unit, clusters, external_ids)
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
    units: list[dict],
    cluster_rows: list[dict],
    head_mode: str,
    second_mode: str,
    external: bool = False,
) -> list[tuple[str, str]]:
    """Clade units (and the Kiss unit): both modes; cluster units: the headline mode only; the `external` line once (cluster30), when the
    freeze names an external set; then `all` once."""
    out: list[tuple[str, str]] = []
    for u in units:
        if u["status"] in ("unit", "kiss_single"):
            out += [(u["unit"], head_mode), (u["unit"], second_mode)]
    kiss_clade = any(u["status"] == "kiss_single" for u in units)
    for u in cluster_rows:
        # a cluster that is the query alone is the same single-protein unit as the query's clade: not run twice
        if u["status"] == "unit" or (u["status"] == "kiss_single" and not kiss_clade):
            out.append((u["unit"], head_mode))
    if external:
        out.append(("external", "cluster30"))
    out.append(("all", "none"))
    return out


def counts_at(
    cut: float,
    proteome_hits: list[dict],
    labelled: set[str],
    neg_groups: dict[str, str],
    neg_hits: dict[str, float],
    sizes: dict[str, int],
    neg_ignored: set[str],
) -> dict:
    """Extra calls per proteome (unlabelled hits at >= cut) and the call rate per hard-negative group at one cutoff."""
    cands = [h for h in proteome_hits if h["id"] not in labelled]
    extra = {
        p: sum(1 for x in cands if x["proteome"] == p and float(x["score"]) >= cut) for p in sizes
    }
    hn = e33.hn_score_table(neg_groups, neg_hits, neg_ignored)
    rates = {g: (sum(1 for s in v if s >= cut) / len(v) if v else 0.0) for g, v in hn.items()}
    return {
        "extra_calls_by_proteome": extra,
        "extra_calls_total": sum(extra.values()),
        "hard_negative_call_rates": rates,
        "hard_negative_max_rate": max(rates.values()) if rates else 0.0,
    }


def run_cutoff_hyd1(
    fz: dict,
    proteome_hits: list[dict],
    labelled: set[str],
    neg_groups: dict[str, str],
    neg_hits: dict[str, float],
    sizes: dict[str, int],
) -> dict:
    """The ecm33 cutoff rule with the owner decision of 2026-10-23 (Q8): the hydrophobin-class panel hits that phmmer does not call family are
    listed by id in the freeze (`panel_ignored`). While proteome_other_hydrophobin is `ignore`, those ids are left out of the extra-call count
    and the hard-negative rate (the headline cutoff). The result ALSO holds the same two numbers with those ids COUNTED: the cutoff that would
    follow, and the extra calls and hard-negative rates the headline cutoff would give. negatives_ignored (the freeze list, empty) is separate."""
    panel = set(fz["panel_ignored"]) if fz["proteome_other_hydrophobin"] == "ignore" else set()
    base_ign = set(fz["negatives_ignored"])
    head_fz = {**fz, "negatives_ignored": sorted(base_ign | (panel & set(neg_groups)))}
    res = e33.run_cutoff(head_fz, proteome_hits, labelled | panel, neg_groups, neg_hits, sizes)
    res["negatives_ignored"] = sorted(base_ign)
    res["panel_ignored"] = sorted(panel)
    res["panel_ignored_negatives_left_out_of_rate"] = sorted(panel & set(neg_groups))
    res["hard_negative_max_rate"] = max(res["hard_negative_call_rates"].values(), default=0.0)
    counted_fz = {**fz, "negatives_ignored": sorted(base_ign)}
    alt = e33.run_cutoff(counted_fz, proteome_hits, labelled, neg_groups, neg_hits, sizes)
    at_head = counts_at(
        res["cutoff"], proteome_hits, labelled, neg_groups, neg_hits, sizes, base_ign
    )
    res["panel_ignored_counted"] = {
        "cutoff_if_counted": alt["cutoff"],
        "extra_calls_total_at_cutoff_if_counted": alt["extra_calls_total"],
        "at_headline_cutoff": at_head,
    }
    return res


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
            external=bool(external_ids_of(fz)),
        ):
            print(f"{unit}\t{mode}")
        return 0
    if a.cmd == "manifest":
        fz = load_freeze(a.freeze)  # no manifest, so no build, before the freeze
        members = read_tsv(a.members)
        clusters = read_cluster_tsv(a.clusters) if a.clusters else None
        man = make_manifest(
            members, a.unit, clusters, a.mode, fz["training_rule"]["reading"], external_ids_of(fz)
        )
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
    fz = load_freeze(a.freeze)
    sizes = {r[0]: int(r[1]) for r in e33._read_two_col(a.sizes)}
    if set(sizes) != set(fz["kept_proteomes"]):
        raise SystemExit(
            f"FATAL: proteomes {sorted(sizes)} are not the frozen kept set {sorted(fz['kept_proteomes'])}"
        )
    res = run_cutoff_hyd1(
        fz,
        [{"proteome": r[0], "id": r[1], "score": r[2]} for r in e33._read_two_col(a.proteome_hits)],
        {line.strip() for line in open(a.labelled) if line.strip()},
        {r[0]: r[1] for r in e33._read_two_col(a.neg_groups)},
        {r[0]: float(r[1]) for r in e33._read_two_col(a.neg_hits)},
        sizes,
    )
    res["freeze_sha256"] = sha256_file(a.freeze)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    alt = res["panel_ignored_counted"]
    print(
        f"cutoff {res['cutoff']} bits; extra calls {res['extra_calls_total']} (limit {res['extra_call_limit_total']:.1f}); "
        f"with the {len(res['panel_ignored'])} panel_ignored ids counted: cutoff {alt['cutoff_if_counted']}, "
        f"extra calls at the headline cutoff {alt['at_headline_cutoff']['extra_calls_total']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
