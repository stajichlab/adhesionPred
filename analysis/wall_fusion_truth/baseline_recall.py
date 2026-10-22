#!/usr/bin/env python3
"""Baseline recall of the Pfam models on the Kiss 2019 control set (plan Task 5.4; spec D8, D9 step 1, D20).

DRAFT, not run. Stdlib only. Run with `/usr/bin/python3.12 -I` on HPCC or `python3 -I` locally.
It reads tables and writes tables. It writes no call, no status and no tier (nothing of
`cellsurface_sorting_hat` is imported), and the files it writes are listed in OUTPUTS.

What a "hit" is: the sequence has a row in the `hmmsearch --cut_ga --domtblout` table of one model,
so it passed that model's GA threshold. No other threshold is used. A model without a GA line is
REPORTED (`no_ga_models.tsv`), never skipped silently and never folded into a recall number.

Per slice (never pooled across slices; there is no line that adds adhesion, wall and fusion):
  adhesion  recall of the 18 active adhesion family rows (family coverage: a control is recalled
            when any of the 18 models hits it; this is NOT the sensitivity of the sorting-hat call)
  wall      recall of the 21 wall and fusion models of analysis/wall_fusion_report/wall_fusion_models.tsv;
            and, as a descriptive false-hit panel, the hits of the wall rows on the 18 family models
            (`wall_rows_false_hits.tsv`): wall rows are never negatives and never truth for the adhesion call
  fusion    descriptive only: hits on the 21 models and on the 2 fusion models; never `estimated`
The shared hard negatives (class `negative`, six of them `<accession>#neg` copies of controls) give the
false hits per model (`negatives_false_hits.tsv`), reported per cluster composition (`negative_only`,
`mixed`). The 35 family representatives are a flagged descriptive group (`family_rep_descriptive.tsv`):
never test controls, never counted toward a floor.

Units and labels (owner decisions of 2026-10-22, report section 10): every number carries its n.
The unit of independence is the MMseqs2 cluster (clusters_split.tsv); a cluster is recalled when at
least one of its control members is hit (`k_clusters_any`), and `k_clusters_all` counts clusters whose
members are all hit. `seed-like` and `seed-independent` are always separate lines, under both
`seed_class_slice` and `seed_class_all39` (`seed_class_kind`). A cell is `estimated` only with at least
5 clusters of its own (the cell's `n_clusters`) and never for fusion; otherwise `descriptive`. The
interval is a Wilson 95% interval on the cluster counts (a proxy for independence, not a bootstrap).

Inputs:
  --domtbl           hmmsearch --cut_ga --domtblout of the 39 models against kiss_sequences.faa
  --members          members_split.tsv (Task 5.3)
  --clusters         clusters_split.tsv (Task 5.3)
  --controls         data/controls/wall-fusion-kiss2019 (controls.tsv AND controls_local_sequence.tsv
                     are the set; family_reps.tsv for the representatives)
  --domtbl-famrep    optional: the same search against family_rep_sequences.faa
  --models-no-ga     optional: one model accession per line that has no GA line (searched with an
                     E-value cut-off instead, in --domtbl-noga / --domtbl-noga-famrep)
  --wall-models, --family-table  the two model lists (default: the repository files)
  --out              a new or empty directory

Usage: python3 -I analysis/wall_fusion_truth/baseline_recall.py --domtbl D --members M --clusters C \
           --controls DIR --out OUT
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WALL_MODELS = ROOT / "analysis/wall_fusion_report/wall_fusion_models.tsv"
FAMILY_TABLE = ROOT / "data/sorting_hat/family_table.tsv"
FLOOR = 5
SLICES = ("adhesion", "wall", "fusion")
SEED_CLASSES = ("seed-like", "seed-independent")
KINDS = ("slice", "all39")
OUTPUTS = (
    "recall_by_slice.tsv",
    "recall_by_clade.tsv",
    "clusters_missed.tsv",
    "wall_rows_false_hits.tsv",
    "negatives_false_hits.tsv",
    "family_rep_descriptive.tsv",
    "no_ga_models.tsv",
    "members_hits.tsv",
    "provenance.tsv",
)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.write_text(
        "\t".join(header) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows)
    )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_model_lists(wall: Path = WALL_MODELS, family: Path = FAMILY_TABLE) -> dict[str, set[str]]:
    wall_rows = read_tsv(wall)
    return {
        "wall_fusion": {r["id"] for r in wall_rows},
        "fusion": {r["id"] for r in wall_rows if r["set"] == "cell_fusion_candidate_models"},
        "family": {
            r["pfam_acc"] for r in read_tsv(family) if r["active"] == "yes" and r["pfam_acc"]
        },
    }


def parse_domtbl(path: Path, require_cut_ga: bool = True) -> dict[str, set[str]]:
    """{sequence name: set of model accessions without version} from an hmmsearch --domtblout file.

    The file must end with the `# [ok]` trailer; with require_cut_ga its option line must show
    --cut_ga (a table made with another threshold is refused)."""
    text = path.read_text()
    lines = text.splitlines()
    if not lines or lines[-1].strip() != "# [ok]":
        raise ValueError(f"{path.name}: no '# [ok]' trailer (hmmsearch did not finish)")
    if require_cut_ga and not any(
        ln.startswith("# Option settings:") and "--cut_ga" in ln for ln in lines
    ):
        raise ValueError(f"{path.name}: the option line does not show --cut_ga")
    hits: dict[str, set[str]] = defaultdict(set)
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split()
        if len(f) < 22:
            raise ValueError(f"{path.name}: short domain line: {ln[:80]}")
        hits[f[0]].add(f[4].split(".")[0])  # query accession without the version
    return dict(hits)


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def f3(x: float) -> str:
    return "" if x != x else f"{x:.3f}"


def status(slice_: str, n_clusters: int) -> str:
    if n_clusters == 0:
        return "none"
    return "estimated" if n_clusters >= FLOOR and slice_ != "fusion" else "descriptive"


def load_controls(controls: Path) -> dict[str, dict]:
    """Both files are the set; keyed by accession (the local key jgi:<portal>:<id> for the second)."""
    out: dict[str, dict] = {}
    for name in ("controls.tsv", "controls_local_sequence.tsv"):
        for r in read_tsv(controls / name):
            if r["accession"] in out:
                raise ValueError(f"{r['accession']}: in both control files")
            out[r["accession"]] = r
    return out


def build_rows(
    members: list[dict],
    controls: dict[str, dict],
    clusters: dict[str, dict],
    hits: dict[str, set[str]],
    models: dict[str, set[str]],
) -> list[dict]:
    """One record per positive control. Kiss rows only: family representatives never enter here."""
    rows = []
    for m in members:
        if m["class"] != "positive":
            continue
        c = controls.get(m["accession"])
        if c is None:
            raise ValueError(f"{m['accession']}: a positive member that is not in the control set")
        sl = c["stratum"]
        own = models["family"] if sl == "adhesion" else models["wall_fusion"]
        cross = models["wall_fusion"] if sl == "adhesion" else models["family"]
        h = hits.get(m["accession"], set())
        cl = clusters[m["cluster_id"]]
        rows.append(
            {
                "accession": m["accession"],
                "gene": c["gene"],
                "species": c["species"],
                "slice": sl,
                "cluster_id": m["cluster_id"],
                "clade_fine": c["clade_fine"],
                "clade_coarse": c["clade_coarse"],
                "cluster_stratum_fine": c["cluster_stratum_fine"],
                "cluster_stratum_coarse": c["cluster_stratum_coarse"],
                "jgi_accept": c["jgi_accept"],
                "seed_class_slice": cl["seed_class_slice"],
                "seed_class_all39": cl["seed_class_all39"],
                "hits": sorted(h),
                "hit_own": bool(h & own),
                "hit_cross": bool(h & cross),
                "hit_fusion2": bool(h & models["fusion"]),
            }
        )
    return rows


def cell(rows: list[dict]) -> tuple[int, int, int, int, int]:
    """n_rows, k_rows, n_clusters, k_clusters_any, k_clusters_all."""
    by: dict[str, list[bool]] = defaultdict(list)
    for r in rows:
        by[r["cluster_id"]].append(r["hit_own"])
    return (
        len(rows),
        sum(r["hit_own"] for r in rows),
        len(by),
        sum(any(v) for v in by.values()),
        sum(all(v) for v in by.values()),
    )


def recall_by_slice(rows: list[dict]) -> list[list]:
    out = []
    for sl in SLICES:
        for kind in KINDS:
            for sc in SEED_CLASSES:
                sub = [r for r in rows if r["slice"] == sl and r[f"seed_class_{kind}"] == sc]
                n, k, nc, kany, kall = cell(sub)
                lo, hi = wilson(kany, nc)
                out.append(
                    [
                        sl,
                        kind,
                        sc,
                        n,
                        k,
                        nc,
                        kany,
                        kall,
                        f3(k / n if n else math.nan),
                        f3(kany / nc if nc else math.nan),
                        f3(lo),
                        f3(hi),
                        status(sl, nc),
                    ]
                )
    return out


def recall_by_clade(rows: list[dict]) -> list[list]:
    out = []
    for level in ("fine", "coarse"):
        for sl in SLICES:
            clades = sorted({r[f"clade_{level}"] for r in rows if r["slice"] == sl})
            for clade in clades:
                base = [r for r in rows if r["slice"] == sl and r[f"clade_{level}"] == clade]
                for stratum in ("any", "clade_only", "shared"):
                    sel = base
                    if stratum != "any":
                        sel = [r for r in base if r[f"cluster_stratum_{level}"] == stratum]
                    for kind in KINDS:
                        for sc in SEED_CLASSES:
                            sub = [r for r in sel if r[f"seed_class_{kind}"] == sc]
                            if not sub:
                                continue
                            n, k, nc, kany, _ = cell(sub)
                            out.append(
                                [
                                    level,
                                    sl,
                                    clade,
                                    stratum,
                                    kind,
                                    sc,
                                    n,
                                    k,
                                    nc,
                                    kany,
                                    f3(k / n),
                                    status(sl, nc),
                                ]
                            )
    return out


def clusters_missed(rows: list[dict], clusters: dict[str, dict]) -> list[list]:
    by: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by[r["cluster_id"]].append(r)
    out = []
    for cid, rs in sorted(by.items(), key=lambda kv: (kv[1][0]["slice"], kv[0])):
        if any(r["hit_own"] for r in rs):
            continue
        cl = clusters[cid]
        out.append(
            [
                rs[0]["slice"],
                cid,
                len(rs),
                ",".join(sorted(r["accession"] for r in rs)),
                ",".join(sorted({r["gene"] for r in rs})),
                ",".join(sorted({r["species"] for r in rs})),
                ",".join(sorted({r["clade_fine"] for r in rs})),
                cl["strata_fine"],
                cl["seed_class_slice"],
                cl["seed_class_all39"],
                cl.get("best_identity_all39", ""),
                ",".join(sorted({r["jgi_accept"] for r in rs})),
            ]
        )
    return out


def wall_panel(rows: list[dict], models: dict[str, set[str]]) -> list[list]:
    wall = [r for r in rows if r["slice"] == "wall"]
    out = []
    for m in sorted(models["family"]):
        hit = [r for r in wall if m in r["hits"]]
        out.append(
            [
                m,
                len(wall),
                len(hit),
                len({r["cluster_id"] for r in hit}),
                len({r["cluster_id"] for r in wall}),
                ",".join(sorted(r["accession"] for r in hit)),
                "descriptive: wall rows are never negatives and never truth for the adhesion call",
            ]
        )
    return out


def negatives_panel(
    members: list[dict], hits: dict[str, set[str]], models: dict[str, set[str]]
) -> list[list]:
    neg = [m for m in members if m["class"] == "negative"]
    groups = {
        "negative_only": [m for m in neg if m["in_mixed_cluster"] != "yes"],
        "mixed": [m for m in neg if m["in_mixed_cluster"] == "yes"],
    }
    out = []
    for model in sorted(models["family"] | models["wall_fusion"]):
        mset = "family" if model in models["family"] else "wall_fusion"
        for comp, ms in groups.items():
            hit = [m for m in ms if model in hits.get(m["accession"], set())]
            out.append(
                [
                    model,
                    mset,
                    comp,
                    len(ms),
                    len(hit),
                    len({m["cluster_id"] for m in hit}),
                    ",".join(sorted(m["accession"] for m in hit)),
                ]
            )
    return out


def famrep_panel(
    reps: list[dict], hits: dict[str, set[str]], models: dict[str, set[str]]
) -> list[list]:
    out = []
    for sl in SLICES:
        rs = [r for r in reps if r["slice"] == sl]
        own = models["family"] if sl == "adhesion" else models["wall_fusion"]
        seqs = sorted({r["accession"] for r in rs})
        k = sum(bool(hits.get(a, set()) & own) for a in seqs)
        out.append(
            [
                sl,
                len(rs),
                len(seqs),
                k,
                "family_rep",
                "no",
                "no",
                "descriptive: a representative is the Afu ortholog, not the Kiss protein",
            ]
        )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--domtbl", type=Path, required=True)
    ap.add_argument("--members", type=Path, required=True)
    ap.add_argument("--clusters", type=Path, required=True)
    ap.add_argument("--controls", type=Path, required=True)
    ap.add_argument("--domtbl-famrep", type=Path)
    ap.add_argument("--models-no-ga", type=Path)
    ap.add_argument("--domtbl-noga", type=Path)
    ap.add_argument("--domtbl-noga-famrep", type=Path)
    ap.add_argument("--wall-models", type=Path, default=WALL_MODELS)
    ap.add_argument("--family-table", type=Path, default=FAMILY_TABLE)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.out.exists() and any(a.out.iterdir()):
        raise SystemExit(f"{a.out} is not empty; nothing is overwritten")

    models = read_model_lists(a.wall_models, a.family_table)
    no_ga = set()
    if a.models_no_ga:
        no_ga = {ln.strip() for ln in a.models_no_ga.read_text().splitlines() if ln.strip()}
    if no_ga - (models["family"] | models["wall_fusion"]):
        raise SystemExit("models-no-ga names a model that is not one of the 39")
    members = read_tsv(a.members)
    clusters = {r["cluster_id"]: r for r in read_tsv(a.clusters)}
    controls = load_controls(a.controls)
    hits = parse_domtbl(a.domtbl)
    for acc in hits:
        if acc not in {m["accession"] for m in members}:
            raise SystemExit(f"{acc}: a hit for a sequence that is not in members.tsv")
    rows = build_rows(members, controls, clusters, hits, models)

    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(
        a.out / "recall_by_slice.tsv",
        [
            "slice",
            "seed_class_kind",
            "seed_class",
            "n_rows",
            "k_rows",
            "n_clusters",
            "k_clusters_any",
            "k_clusters_all",
            "row_recall",
            "cluster_recall_any",
            "wilson95_low",
            "wilson95_high",
            "status",
        ],
        recall_by_slice(rows),
    )
    write_tsv(
        a.out / "recall_by_clade.tsv",
        [
            "level",
            "slice",
            "clade",
            "stratum",
            "seed_class_kind",
            "seed_class",
            "n_rows",
            "k_rows",
            "n_clusters",
            "k_clusters_any",
            "row_recall",
            "status",
        ],
        recall_by_clade(rows),
    )
    write_tsv(
        a.out / "clusters_missed.tsv",
        [
            "slice",
            "cluster_id",
            "n_members",
            "accessions",
            "genes",
            "species",
            "clades_fine",
            "strata_fine",
            "seed_class_slice",
            "seed_class_all39",
            "best_identity_all39",
            "jgi_accept",
        ],
        clusters_missed(rows, clusters),
    )
    write_tsv(
        a.out / "wall_rows_false_hits.tsv",
        [
            "family_model",
            "n_wall_rows",
            "n_hit_rows",
            "n_hit_clusters",
            "n_wall_clusters",
            "hit_accessions",
            "note",
        ],
        wall_panel(rows, models),
    )
    write_tsv(
        a.out / "negatives_false_hits.tsv",
        [
            "model",
            "model_set",
            "cluster_composition",
            "n_negatives",
            "n_hit",
            "n_hit_clusters",
            "hit_accessions",
        ],
        negatives_panel(members, hits, models),
    )
    reps = read_tsv(a.controls / "family_reps.tsv")
    rep_hits = parse_domtbl(a.domtbl_famrep) if a.domtbl_famrep else {}
    write_tsv(
        a.out / "family_rep_descriptive.tsv",
        [
            "slice",
            "n_rep_rows",
            "n_rep_sequences",
            "k_hit_own_models",
            "group",
            "counts_toward_floor",
            "may_be_test_control",
            "note",
        ],
        famrep_panel(reps, rep_hits, models) if a.domtbl_famrep else [],
    )
    noga_rows = []
    noga_hits = parse_domtbl(a.domtbl_noga, require_cut_ga=False) if a.domtbl_noga else {}
    for m in sorted(no_ga):
        n_seq = sum(m in v for v in noga_hits.values())
        noga_rows.append(
            [m, "no GA line", "E-value cut-off, not GA", n_seq, "not in any recall number"]
        )
    write_tsv(
        a.out / "no_ga_models.tsv",
        ["model", "reason", "threshold_used", "n_sequences_hit", "note"],
        noga_rows,
    )
    write_tsv(
        a.out / "members_hits.tsv",
        [
            "accession",
            "gene",
            "species",
            "slice",
            "cluster_id",
            "clade_fine",
            "seed_class_slice",
            "seed_class_all39",
            "hit_own_slice_models",
            "hit_other_slice_models",
            "hit_fusion_models",
            "models_hit",
        ],
        [
            [
                r["accession"],
                r["gene"],
                r["species"],
                r["slice"],
                r["cluster_id"],
                r["clade_fine"],
                r["seed_class_slice"],
                r["seed_class_all39"],
                "yes" if r["hit_own"] else "no",
                "yes" if r["hit_cross"] else "no",
                "yes" if r["hit_fusion2"] else "no",
                ",".join(r["hits"]),
            ]
            for r in rows
        ],
    )
    prov = [
        ["domtbl_sha256", sha256(a.domtbl)],
        ["members_sha256", sha256(a.members)],
        ["clusters_sha256", sha256(a.clusters)],
        ["n_positive_rows", len(rows)],
        ["n_models_family", len(models["family"])],
        ["n_models_wall_fusion", len(models["wall_fusion"])],
        ["n_models_no_ga", len(no_ga)],
        ["floor", FLOOR],
    ]
    write_tsv(a.out / "provenance.tsv", ["key", "value"], prov)
    return 0


if __name__ == "__main__":
    sys.exit(main())
