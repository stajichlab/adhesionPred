#!/usr/bin/env python3
"""PR 6 step 6.1 (build) for ecm33: the per-unit training manifest, the independence check and the cutoff rule.

DRAFT, not run on HPCC. Stdlib only; run with `python3 -I` (`/usr/bin/python3.12 -I` on HPCC). It reads tables and writes tables and
JSON. It writes no call, no status, no tier, and no sequence. Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md
(sections 4.3 to 4.5 and 4.8; decisions 7 and 8).

Subcommands
  manifest         the training set of ONE held-out unit (a fine clade) or of the `all` model, after the independence rule
  check-alignment  the build alignment contains exactly the manifest's training ids: no held-out id, no excluded 30 percent neighbour
  cutoff           the bit-score cutoff, from the FALSE-HIT side only (the four kept proteomes and the shared hard negatives);
                   refuses unless the freeze file is frozen

Independence (spec 4.4). A training member is excluded when it falls in the same MMseqs2 cluster (30 percent identity, 50 percent
coverage: the Task 5.3 rule, the `_cluster.tsv` of `mmseqs easy-cluster` over the members) as any held-out member of the unit.
`--mode none` skips that exclusion (kept only to report how much it removes: for a conserved family the 30 percent clusters can
span several clades, so the exclusion can empty the training set, and the manifest then says `not_buildable`).

Cutoff (decision 8, the rule of analysis/hydrophobin_truth/relaxed.py `lowest_cutoff`): the lowest bit score at or above the floor at
which (i) every kept proteome has at most per10k unlabelled extra calls per 10,000 proteins and (ii) every hard-negative group (except
the ignored accessions) has a call rate of at most hn_rate. NO positive score is an input: the function does not take one. The limits
and the floor are read from the freeze file, which must say `"frozen": true` before any held-out score exists.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

MIN_TRAIN_SEQ = 5  # decision 7, provisional
MIN_TRAIN_CLADES = 2


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fasta_id(m: dict) -> str:
    """The header the extract job gives a member: the query keeps its accession, a UniRef50 representative its cluster id."""
    return (
        m["query_accession"]
        if m.get("is_query") == "yes" and m.get("query_accession")
        else m["hit_id"]
    )


def read_cluster_tsv(path: Path) -> dict[str, str]:
    """{member: representative} from the `_cluster.tsv` of mmseqs easy-cluster (two columns: representative, member)."""
    out = {}
    for line in open(path):
        if not line.strip():
            continue
        rep, mem = line.rstrip("\n").split("\t")[:2]
        out[mem] = rep
    return out


# ---------------------------------------------------------------------------------------------------------------------------
# manifest
# ---------------------------------------------------------------------------------------------------------------------------


def make_manifest(
    members: list[dict],
    unit: str,
    clusters: dict[str, str] | None,
    mode: str = "cluster30",
    min_seq: int = MIN_TRAIN_SEQ,
    min_clades: int = MIN_TRAIN_CLADES,
) -> dict:
    """Training set of `unit` (a fine clade; `all` = no held-out unit). Held-out members are never in the training set; with
    mode `cluster30` neither is any training member in the same 30 percent cluster as a held-out member."""
    if mode not in ("cluster30", "none"):
        raise ValueError(f"unknown mode {mode}")
    held = [m for m in members if unit != "all" and m["clade_fine"] == unit]
    if unit != "all" and not held:
        raise ValueError(f"no members in clade {unit}")
    held_ids = {fasta_id(m) for m in held}
    train = [m for m in members if fasta_id(m) not in held_ids]
    excluded = []
    if mode == "cluster30" and held:
        if clusters is None:
            raise ValueError("mode cluster30 needs the MMseqs2 cluster table")
        held_reps = {clusters[i] for i in held_ids if i in clusters}
        missing = sorted(i for i in held_ids if i not in clusters)
        if missing:
            raise ValueError(f"held-out ids missing from the cluster table: {missing[:5]}")
        kept = []
        for m in train:
            if clusters.get(fasta_id(m)) in held_reps:
                excluded.append(m)
            else:
                kept.append(m)
        train = kept
    clades: dict[str, int] = {}
    for m in train:
        clades[m["clade_fine"]] = clades.get(m["clade_fine"], 0) + 1
    rule = len(train) >= min_seq and len(clades) >= min_clades
    reason = ""
    if not rule:
        reason = (
            f"training set after the independence rule has {len(train)} sequences in {len(clades)} clades "
            f"(minimum {min_seq} in {min_clades})"
        )
    return {
        "unit": unit,
        "mode": mode,
        "held_out_ids": sorted(held_ids),
        "training": [
            {
                "id": fasta_id(m),
                "clade_fine": m["clade_fine"],
                "cluster": (clusters or {}).get(fasta_id(m), ""),
            }
            for m in sorted(train, key=fasta_id)
        ],
        "excluded_30pct_neighbours": sorted(fasta_id(m) for m in excluded),
        "n_held_out": len(held_ids),
        "n_training": len(train),
        "n_training_clades": len(clades),
        "n_excluded_neighbours": len(excluded),
        "training_rule_met": rule,
        "status": "buildable" if rule else "not_buildable",
        "reason": reason,
    }


def check_alignment_ids(manifest: dict, aligned_ids: list[str]) -> list[str]:
    """Problems (empty = fine): a held-out or excluded id in the alignment, an id not in the manifest, a training id missing."""
    train = {t["id"] for t in manifest["training"]}
    held = set(manifest["held_out_ids"])
    excl = set(manifest["excluded_30pct_neighbours"])
    got = set(aligned_ids)
    problems = []
    if len(got) != len(aligned_ids):
        problems.append("duplicate ids in the alignment")
    for i in sorted(got & held):
        problems.append(f"held-out id {i} is in the alignment")
    for i in sorted(got & excl):
        problems.append(f"excluded 30 percent neighbour {i} is in the alignment")
    for i in sorted(got - train - held - excl):
        problems.append(f"id {i} is not in the manifest")
    for i in sorted(train - got):
        problems.append(f"training id {i} is missing from the alignment")
    return problems


def fasta_ids(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in open(path) if line.startswith(">")]


def read_aligned(path: Path) -> list[str]:
    seqs, buf = [], []
    for line in open(path):
        if line.startswith(">"):
            if buf:
                seqs.append("".join(buf))
            buf = []
        else:
            buf.append(line.strip())
    if buf:
        seqs.append("".join(buf))
    return seqs


def occupancy_check(
    aligned: list[str], occupancy: float, min_columns: int, min_median_coverage: float
) -> dict:
    """The pre-registered stop condition of the build (a PROPOSAL in the freeze file; the hydrophobin precedent stopped on a cysteine
    column check, which does not apply here): at least `min_columns` columns hold a residue in at least `occupancy` of the sequences, and
    the median sequence has a residue in at least `min_median_coverage` of those columns. It is a stop condition, not a tuning step:
    a failing alignment stops the unit and is reported; no other aligner is tried silently."""
    n = len(aligned)
    if n == 0:
        return {"ok": False, "n_columns": 0, "median_coverage": 0.0, "reason": "empty alignment"}
    width = len(aligned[0])
    cols = [j for j in range(width) if sum(1 for s in aligned if s[j] not in "-.") / n >= occupancy]
    if not cols:
        return {"ok": False, "n_columns": 0, "median_coverage": 0.0, "reason": "no occupied column"}
    cov = sorted(sum(1 for j in cols if s[j] not in "-.") / len(cols) for s in aligned)
    median = cov[n // 2] if n % 2 else (cov[n // 2 - 1] + cov[n // 2]) / 2
    ok = len(cols) >= min_columns and median >= min_median_coverage
    return {
        "ok": ok,
        "n_columns": len(cols),
        "median_coverage": median,
        "reason": "" if ok else "below the pre-registered thresholds",
    }


# ---------------------------------------------------------------------------------------------------------------------------
# freeze and cutoff
# ---------------------------------------------------------------------------------------------------------------------------

FREEZE_KEYS = (
    "frozen",
    "target",
    "cutoff_floor",
    "extra_calls_per_10000",
    "hard_negative_rate",
    "min_u",
    "min_unit_members",
    "independence_mode",
    "negatives_ignored",
    "kept_proteomes",
    "alignment_check",
)


def load_freeze(path: Path) -> dict:
    """The freeze file, or a refusal: scoring and the cutoff need it BEFORE any held-out score is looked at."""
    if not Path(path).is_file():
        raise SystemExit(
            f"FATAL: freeze file {path} is missing; nothing is scored before the freeze"
        )
    fz = json.load(open(path))
    missing = [k for k in FREEZE_KEYS if k not in fz]
    if missing:
        raise SystemExit(f"FATAL: freeze file {path} lacks {missing}")
    if fz["frozen"] is not True:
        raise SystemExit(
            f'FATAL: freeze file {path} is not frozen ("frozen": true is set by the owner\'s freeze commit)'
        )
    return fz


def _ok(c, cands, sizes, hn_scores, per10k, hn_rate, ignore):
    by: dict[str, int] = {}
    for x in cands:
        if x["score"] >= c:
            by[x["proteome"]] = by.get(x["proteome"], 0) + 1
    for p, n in by.items():
        if n > per10k * sizes[p] / 10000.0:
            return False
    for g, scores in hn_scores.items():
        if g in ignore or not scores:
            continue
        if sum(1 for s in scores if s >= c) / len(scores) > hn_rate:
            return False
    return True


def lowest_cutoff(cands, sizes, hn_scores, per10k, hn_rate, floor=0.0, ignore_groups=()):
    """Same rule as analysis/hydrophobin_truth/relaxed.py lowest_cutoff (a test compares the two). cands: dicts with id, proteome,
    score (unlabelled hits only); sizes: {proteome: n proteins}; hn_scores: {group: [score per negative, 0 when not hit]}."""
    ignore = set(ignore_groups)
    scores = sorted(
        {x["score"] for x in cands if x["score"] >= floor}
        | {s for g, v in hn_scores.items() if g not in ignore for s in v if s >= floor}
    )
    for c in [floor] + scores:
        if _ok(c, cands, sizes, hn_scores, per10k, hn_rate, ignore):
            included = [s for s in scores if s >= c]
            return included[0] if included else floor
    return scores[-1] + 0.1 if scores else floor


def extra_call_limit_total(sizes: dict[str, int], per10k: float) -> float:
    """5 per 10,000 over 32,491 proteins = 16.2455 (about 16 extra calls in all)."""
    return per10k * sum(sizes.values()) / 10000.0


def hn_score_table(
    neg_groups: dict[str, str], neg_hits: dict[str, float], ignored: set[str]
) -> dict[str, list[float]]:
    """{group: score per negative}; a negative that is not hit scores 0; the ignored accessions are left out of every group."""
    out: dict[str, list[float]] = {}
    for acc, g in neg_groups.items():
        if acc in ignored:
            continue
        out.setdefault(g, []).append(neg_hits.get(acc, 0.0))
    return out


def run_cutoff(
    fz: dict,
    proteome_hits: list[dict],
    labelled: set[str],
    neg_groups: dict,
    neg_hits: dict,
    sizes: dict,
) -> dict:
    cands = [
        {"id": h["id"], "proteome": h["proteome"], "score": float(h["score"])}
        for h in proteome_hits
        if h["id"] not in labelled
    ]
    ignored = set(fz["negatives_ignored"])
    hn = hn_score_table(neg_groups, neg_hits, ignored)
    per10k = float(fz["extra_calls_per_10000"])
    rate = float(fz["hard_negative_rate"])
    floor = float(fz["cutoff_floor"])
    cut = lowest_cutoff(cands, sizes, hn, per10k, rate, floor)
    extra = {p: sum(1 for x in cands if x["proteome"] == p and x["score"] >= cut) for p in sizes}
    return {
        "cutoff": cut,
        "floor": floor,
        "extra_calls_per_10000": per10k,
        "extra_calls_by_proteome": extra,
        "extra_calls_total": sum(extra.values()),
        "extra_call_limit_by_proteome": {p: per10k * n / 10000.0 for p, n in sizes.items()},
        "extra_call_limit_total": extra_call_limit_total(sizes, per10k),
        "hard_negative_rate_limit": rate,
        "hard_negative_call_rates": {
            g: (sum(1 for s in v if s >= cut) / len(v) if v else 0.0) for g, v in hn.items()
        },
        "negatives_ignored": sorted(ignored),
        "positives_used": False,
    }


# ---------------------------------------------------------------------------------------------------------------------------


def _read_two_col(path: Path, header: bool = True) -> list[list[str]]:
    rows = [line.rstrip("\n").split("\t") for line in open(path) if line.strip()]
    return rows[1:] if header else rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("manifest")
    m.add_argument("--members", required=True, type=Path)
    m.add_argument("--unit", required=True, help="a fine clade, or `all`")
    m.add_argument(
        "--clusters", type=Path, help="mmseqs easy-cluster _cluster.tsv over the members"
    )
    m.add_argument("--mode", default=None, choices=["cluster30", "none"])
    m.add_argument("--freeze", type=Path)
    m.add_argument("--out", required=True, type=Path)
    c = sub.add_parser("check-alignment")
    c.add_argument("--manifest", required=True, type=Path)
    c.add_argument("--alignment", required=True, type=Path)
    c.add_argument(
        "--freeze", type=Path, help="also run the pre-registered occupancy check of the freeze file"
    )
    k = sub.add_parser("cutoff")
    k.add_argument("--freeze", required=True, type=Path)
    k.add_argument(
        "--proteome-hits",
        required=True,
        type=Path,
        help="TSV proteome, id, score (every hit at or above the floor)",
    )
    k.add_argument(
        "--labelled",
        required=True,
        type=Path,
        help="ids of the proteome proteins that are family members (one per line)",
    )
    k.add_argument(
        "--neg-groups",
        required=True,
        type=Path,
        help="TSV accession, group (every shared negative)",
    )
    k.add_argument(
        "--neg-hits",
        required=True,
        type=Path,
        help="TSV accession, score (the negatives that were hit)",
    )
    k.add_argument("--sizes", required=True, type=Path, help="TSV proteome, n_proteins")
    k.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)

    if a.cmd == "manifest":
        members = read_tsv(a.members)
        mode = a.mode
        if mode is None:
            mode = load_freeze(a.freeze)["independence_mode"] if a.freeze else "cluster30"
        clusters = read_cluster_tsv(a.clusters) if a.clusters else None
        man = make_manifest(members, a.unit, clusters, mode)
        man["members_sha256"] = sha256_file(a.members)
        a.out.mkdir(parents=True, exist_ok=True)
        (a.out / "manifest.json").write_text(json.dumps(man, indent=1, sort_keys=True) + "\n")
        (a.out / "train_ids.txt").write_text("".join(t["id"] + "\n" for t in man["training"]))
        (a.out / "heldout_ids.txt").write_text("".join(i + "\n" for i in man["held_out_ids"]))
        print(
            f"unit {a.unit}: held out {man['n_held_out']}, training {man['n_training']} in {man['n_training_clades']} clades, "
            f"excluded 30% neighbours {man['n_excluded_neighbours']}: {man['status']} {man['reason']}"
        )
        return 0 if man["status"] == "buildable" else 3
    if a.cmd == "check-alignment":
        man = json.load(open(a.manifest))
        problems = check_alignment_ids(man, fasta_ids(a.alignment))
        for p in problems:
            print("FATAL:", p, file=sys.stderr)
        if problems:
            return 1
        if a.freeze:
            ac = load_freeze(a.freeze)["alignment_check"]
            r = occupancy_check(
                read_aligned(a.alignment),
                float(ac["occupancy"]),
                int(ac["min_columns"]),
                float(ac["min_median_coverage"]),
            )
            print(f"alignment check: {r}")
            if not r["ok"]:
                print(
                    "FATAL: the pre-registered alignment check failed; stop and report (no other aligner is tried)",
                    file=sys.stderr,
                )
                return 4
        return 0
    fz = load_freeze(a.freeze)
    sizes = {r[0]: int(r[1]) for r in _read_two_col(a.sizes)}
    kept = set(fz["kept_proteomes"])
    if set(sizes) != kept:
        raise SystemExit(
            f"FATAL: proteomes {sorted(sizes)} are not the frozen kept set {sorted(kept)}"
        )
    res = run_cutoff(
        fz,
        [{"proteome": r[0], "id": r[1], "score": r[2]} for r in _read_two_col(a.proteome_hits)],
        {line.strip() for line in open(a.labelled) if line.strip()},
        {r[0]: r[1] for r in _read_two_col(a.neg_groups)},
        {r[0]: float(r[1]) for r in _read_two_col(a.neg_hits)},
        sizes,
    )
    res["freeze_sha256"] = sha256_file(a.freeze)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(
        f"cutoff {res['cutoff']} bits; extra calls {res['extra_calls_total']} (limit {res['extra_call_limit_total']:.1f})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
