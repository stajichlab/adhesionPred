#!/usr/bin/env python3
"""PR 6: specificity tests of the Pfam families that the step 6.0 scan names as candidates.

DRAFT, not run. Stdlib only. Run `python3 -I analysis/wall_fusion_truth/pr6_family_specificity.py <command>` (on HPCC:
`/usr/bin/python3.12 -I`). Owner decisions of 2026-10-22: the candidate list is the union of (a) a Pfam hit that covers at
least 0.5 of the sequence and (b) a best hit with an E-value of 1e-30 or better (both limits inclusive); every candidate
family gets a specificity test BEFORE any panel addition; nothing enters the active panel without a `family_measure.tsv`
row AND the owner's sign-off; `narrow_ok` stays `no` until the owner signs it.

THIS SCRIPT NEVER WRITES `data/sorting_hat/family_measure.tsv` OR `family_table.tsv` AND CHANGES NO TIER, NO PANEL AND NO
STATUS. It writes report-only tables. The only bridge to the measurement file is the `convert` command, which writes a
PROPOSED measurement file somewhere else, prints the tier that `cellsurface_sorting_hat.family_tier` derives from it
(never hand-set) and refuses the two committed paths.

Commands:
  candidates  scan_per_cluster.tsv (the step 6.0 record) + clusters_missed.tsv -> candidate_clusters.tsv,
              candidate_families.tsv, family_members.tsv (the family's own evidenced members), families.txt (the
              model NAMES for `hmmfetch -f`, one per line). Offline and deterministic.
  summarize   the candidate tables + `hmmsearch --cut_ga --domtblout` files of the Kiss positives, the shared hard
              negatives and the kept proteomes -> family_specificity.tsv (one row per family), positives_hits.tsv,
              negatives_hits.tsv, proteome_hits.tsv (every hit with the FASTA header text), provenance.tsv.
  convert     family_specificity.tsv -> a proposed `family_measure.tsv`-format file (`measure_basis=descriptive_shared_sets`,
              so the derived tier is `smoke` at most), a dry run unless --out is given.

Method notes (also in README_pr6_family_specificity.md):
  * The family's own evidenced members (D19) are all members of the candidate clusters in which the family has a
    qualifying hit. The own-member recall is CIRCULAR: the family was selected because it hits these members. It is
    reported as a count, never as a sensitivity, and the `leakage` field says `tuned_on_truth`.
  * Negatives are the shared hard-negative set of 161 rows (137 N1 and 24 N2; `data/controls/hard-negatives-shared`),
    scored under their plain accession (the `#neg` copies exist only in the Kiss clustering FASTA and are not used here;
    the six accessions that are both a Kiss control and a shared negative are positives in one table and negatives in
    the other, as decided under D20). Counts follow the 2026-10-09 family test: N1 rows and N1 clusters are the
    denominators, N2 hits are a separate count.
  * The four kept proteomes (one per species, decision 9b) give hits per proteome and per 10,000 proteins, with the
    FASTA header text of every hit. This is descriptive: the proteomes are not labelled.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCAN_RECORD = ROOT / "docs/reports/data/pr6_pfam_scan_2026-10-22/scan_per_cluster.tsv"
MISSED = ROOT / "docs/reports/data/kiss_baseline_2026-10-22/clusters_missed.tsv"
MEMBERS_SPLIT = ROOT / "docs/reports/data/kiss_cluster_split_2026-10-22/members_split.tsv"
NEG_DIR = ROOT / "data/controls/hard-negatives-shared"
FAMILY_TABLE = ROOT / "data/sorting_hat/family_table.tsv"
FAMILY_MEASURE = ROOT / "data/sorting_hat/family_measure.tsv"
INPUTS_TSV = ROOT / "scripts/sorting_hat/refresh_3_10/inputs.tsv"

COVERAGE_MIN = 0.5  # inclusive; owner decision 1 of 2026-10-22
EVALUE_MAX = 1e-30  # inclusive; owner decision 1 of 2026-10-22
KEPT_PROTEOMES = ("Afum_Af293_UniProt", "Calb_SC5314", "Cimm_RS", "Scer_S288C")  # decision 9b
LEAKAGE = "tuned_on_truth"  # the family was selected because it hits its own members
SLICES = ("adhesion", "wall", "fusion")
ANNOTATION_MAX = 160
PER = 10_000

HIT_RE = re.compile(r"^(PF\d+):([^:]+):frac=([0-9.]+):score=([0-9.]+):ievalue=([^:]+):seq=(.+)$")

CANDIDATE_OUTPUTS = (
    "candidate_clusters.tsv",
    "candidate_families.tsv",
    "family_members.tsv",
    "families.txt",
)
SUMMARY_OUTPUTS = (
    "family_specificity.tsv",
    "positives_hits.tsv",
    "negatives_hits.tsv",
    "proteome_hits.tsv",
    "provenance.tsv",
)


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---- the candidate rule ----------------------------------------------------------------------------


def parse_other_pfam(cell: str) -> list[dict]:
    """The hits of one `other_pfam_families` cell of scan_per_cluster.tsv (entries are separated by `;`)."""
    hits = []
    for entry in (cell or "").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        m = HIT_RE.match(entry)
        if not m:
            raise ValueError(f"unreadable other_pfam_families entry: {entry[:100]}")
        acc, name, frac, score, ievalue, seq = m.groups()
        hits.append(
            {
                "acc": acc,
                "name": name,
                "frac": float(frac),
                "score": float(score),
                "ievalue": float(ievalue),
                "seq": seq,
            }
        )
    return hits


def rule_met(hit: dict) -> tuple[bool, bool]:
    """(coverage rule, E-value rule) for one hit; both limits are inclusive."""
    return hit["frac"] >= COVERAGE_MIN, hit["ievalue"] <= EVALUE_MAX


def rule_label(hit: dict) -> str:
    cov, ev = rule_met(hit)
    return "coverage+evalue" if cov and ev else "coverage" if cov else "evalue"


def select_candidates(scan_rows: list[dict], missed_rows: list[dict]) -> dict:
    """The candidate clusters and families.

    A cluster is a candidate when its class is `hit_other_pfam` and at least one of its hits meets a rule. A family is a
    candidate when one of ITS hits meets a rule (a weaker hit of the same cluster does not make its family a
    candidate). `hit_panel_model` and `no_pfam_hit` clusters are not candidates (the panel families have their own
    measurement; the no-hit clusters are the gap search)."""
    members = {r["cluster_id"]: [a for a in r["accessions"].split(",") if a] for r in missed_rows}
    clusters, families = [], {}
    for row in sorted(scan_rows, key=lambda r: r["cluster_id"]):
        if row["class"] != "hit_other_pfam":
            continue
        if row["cluster_id"] not in members:
            raise ValueError(f"{row['cluster_id']} is not in clusters_missed.tsv")
        qualifying = [h for h in parse_other_pfam(row["other_pfam_families"]) if any(rule_met(h))]
        if not qualifying:
            continue
        by_family: dict[str, list[dict]] = defaultdict(list)
        for h in qualifying:
            by_family[h["acc"]].append(h)
        clusters.append(
            {
                "slice": row["slice"],
                "cluster_id": row["cluster_id"],
                "group": row["group"],
                "genes": row["genes"],
                "clades_fine": row["clades_fine"],
                "strata_fine": row["strata_fine"],
                "jgi_accept": row["jgi_accept"],
                "members": members[row["cluster_id"]],
                "families": sorted(by_family),
                "rules": sorted({rule_label(h) for h in qualifying}),
            }
        )
        for acc, hs in by_family.items():
            fam = families.setdefault(
                acc,
                {"acc": acc, "name": hs[0]["name"], "clusters": set(), "slices": set(), "hits": []},
            )
            fam["clusters"].add(row["cluster_id"])
            fam["slices"].add(row["slice"])
            fam["hits"].extend(hs)
    return {"clusters": clusters, "families": families, "members": members}


def write_candidates(
    scan: Path, missed: Path, out: Path, family_table: Path = FAMILY_TABLE
) -> dict:
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"{out} exists and is not empty; nothing is overwritten")
    sel = select_candidates(read_tsv(scan), read_tsv(missed))
    in_table = {r["pfam_acc"] for r in read_tsv(family_table) if r["pfam_acc"]}
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(
        out / "candidate_clusters.tsv",
        [
            "slice",
            "cluster_id",
            "group",
            "genes",
            "clades_fine",
            "strata_fine",
            "jgi_accept",
            "n_members",
            "candidate_families",
            "rule",
        ],
        [
            [
                c["slice"],
                c["cluster_id"],
                c["group"],
                c["genes"],
                c["clades_fine"],
                c["strata_fine"],
                c["jgi_accept"],
                len(c["members"]),
                ";".join(c["families"]),
                ";".join(c["rules"]),
            ]
            for c in sel["clusters"]
        ],
    )
    fam_rows, member_rows = [], []
    for acc in sorted(sel["families"]):
        f = sel["families"][acc]
        cov = any(rule_met(h)[0] for h in f["hits"])
        ev = any(rule_met(h)[1] for h in f["hits"])
        fam_rows.append(
            [
                acc,
                f["name"],
                len(f["clusters"]),
                ";".join(sorted(f["clusters"])),
                ";".join(s for s in SLICES if s in f["slices"]),
                len(f["hits"]),
                f"{max(h['frac'] for h in f['hits']):.3f}",
                f"{min(h['ievalue'] for h in f['hits']):.3g}",
                "yes" if cov else "no",
                "yes" if ev else "no",
                "yes" if acc in in_table else "no",
            ]
        )
        for cid in sorted(f["clusters"]):
            for a in sel["members"][cid]:
                member_rows.append([acc, cid, a])
    write_tsv(
        out / "candidate_families.tsv",
        [
            "pfam_acc",
            "name",
            "n_clusters",
            "clusters",
            "slices",
            "n_hits",
            "max_frac",
            "min_ievalue",
            "rule_coverage",
            "rule_evalue",
            "in_family_table",
        ],
        fam_rows,
    )
    write_tsv(out / "family_members.tsv", ["pfam_acc", "cluster_id", "accession"], member_rows)
    (out / "families.txt").write_text(
        "".join(f"{sel['families'][a]['name']}\n" for a in sorted(sel["families"]))
    )
    return sel


# ---- hmmsearch tables and FASTA headers ------------------------------------------------------------


def parse_domtbl(path: Path) -> dict[str, set[str]]:
    """{target sequence name: set of family accessions without version} from `hmmsearch --cut_ga --domtblout`.

    The file must end with the `# [ok]` trailer and its option line must show --cut_ga (a table made with another
    threshold is refused). The orientation must be hmmsearch's (column 5 is the Pfam accession of the model)."""
    lines = Path(path).read_text().splitlines()
    if not lines or lines[-1].strip() != "# [ok]":
        raise ValueError(f"{Path(path).name}: no '# [ok]' trailer (hmmsearch did not finish)")
    if not any(ln.startswith("# Option settings:") and "--cut_ga" in ln for ln in lines):
        raise ValueError(f"{Path(path).name}: the option line does not show --cut_ga")
    hits: dict[str, set[str]] = defaultdict(set)
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split()
        if len(f) < 22:
            raise ValueError(f"{Path(path).name}: short domain line: {ln[:80]}")
        if not f[4].startswith("PF"):
            raise ValueError(
                f"{Path(path).name}: column 5 is not a Pfam accession (hmmscan orientation?)"
            )
        hits[f[0]].add(f[4].split(".")[0])
    return dict(hits)


def fasta_headers(path: Path) -> dict[str, str]:
    """{first token of the header: the rest of the header}. Sequence lines are never read into the result."""
    out: dict[str, str] = {}
    with Path(path).open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                head = ln[1:].strip()
                name, _, rest = head.partition(" ")
                if name in out:
                    raise ValueError(f"{Path(path).name}: duplicate sequence name {name}")
                out[name] = rest
    return out


def kept_proteome_counts(inputs: Path = INPUTS_TSV) -> dict[str, int]:
    """{kept proteome: number of proteins} from the committed refresh inputs table."""
    table = {r["prot"]: int(r["n_proteins"]) for r in read_tsv(inputs) if r["n_proteins"].isdigit()}
    missing = [p for p in KEPT_PROTEOMES if p not in table]
    if missing:
        raise ValueError(f"{inputs}: no protein count for {missing}")
    return {p: table[p] for p in KEPT_PROTEOMES}


# ---- the specificity table -------------------------------------------------------------------------


def summarize(
    cand_dir: Path,
    members_split: Path,
    neg_dir: Path,
    domtbl_pos: Path,
    domtbl_neg: Path,
    proteomes: dict[str, tuple[Path, Path]],
    out: Path,
    expected_proteins: dict[str, int] | None = None,
) -> list[dict]:
    """One row per candidate family. `proteomes` is {name: (domtbl, fasta)}."""
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"{out} exists and is not empty; nothing is overwritten")
    cand_dir = Path(cand_dir)
    families = read_tsv(cand_dir / "candidate_families.tsv")
    own: dict[str, set[str]] = defaultdict(set)
    own_clusters: dict[str, set[str]] = defaultdict(set)
    for r in read_tsv(cand_dir / "family_members.tsv"):
        own[r["pfam_acc"]].add(r["accession"])
        own_clusters[r["pfam_acc"]].add(r["cluster_id"])
    split = {r["accession"]: r for r in read_tsv(members_split) if r["class"] == "positive"}
    # negatives: the committed shared set. N1 and N2 as in the 2026-10-09 family test.
    neg_level = {
        r["accession"]: r["evidence_level"] for r in read_tsv(Path(neg_dir) / "controls.tsv")
    }
    neg_cluster = {
        r["accession"]: r["cluster_id"]
        for r in read_tsv(Path(neg_dir) / "clusters.tsv")
        if r["set"] == "hard_negative"
    }
    n1 = {a for a, lv in neg_level.items() if lv == "N1"}
    n2 = {a for a, lv in neg_level.items() if lv == "N2"}
    if not n1 or set(neg_level) != n2 | n1:
        raise ValueError("the shared negative set must hold N1 and N2 rows only")
    n1_clusters = {neg_cluster[a] for a in n1}
    pos_hits, neg_hits = parse_domtbl(domtbl_pos), parse_domtbl(domtbl_neg)
    stray = set(pos_hits) - set(split)
    if stray:
        raise ValueError(
            f"{Path(domtbl_pos).name}: hits on sequences that are not Kiss positives: {sorted(stray)[:3]}"
        )
    stray = set(neg_hits) - set(neg_level)
    if stray:
        raise ValueError(
            f"{Path(domtbl_neg).name}: hits on sequences that are not shared negatives: {sorted(stray)[:3]}"
        )
    prot_hits, prot_head, prot_n = {}, {}, {}
    for name, (dom, fasta) in proteomes.items():
        prot_head[name] = fasta_headers(fasta)
        prot_n[name] = len(prot_head[name])
        if expected_proteins is not None and prot_n[name] != expected_proteins.get(name):
            raise ValueError(
                f"{name}: {prot_n[name]} proteins in the FASTA, expected {expected_proteins.get(name)}"
            )
        prot_hits[name] = parse_domtbl(dom)
        stray = set(prot_hits[name]) - set(prot_head[name])
        if stray:
            raise ValueError(
                f"{name}: hits on sequences that are not in its FASTA: {sorted(stray)[:3]}"
            )
    kept_total = sum(prot_n.values())
    names = list(proteomes)

    def hit_set(table: dict[str, set[str]], acc: str) -> set[str]:
        return {s for s, accs in table.items() if acc in accs}

    rows, pos_rows, neg_rows, prot_rows = [], [], [], []
    for fam in families:
        acc = fam["pfam_acc"]
        members = own[acc]
        absent = members - set(split)
        if absent:
            raise ValueError(
                f"{acc}: own members that are not Kiss positives: {sorted(absent)[:3]}"
            )
        pos = hit_set(pos_hits, acc)
        own_hit = pos & members
        other = pos - members
        other_by_slice = {s: len([p for p in other if split[p]["slice"] == s]) for s in SLICES}
        h1, h2 = hit_set(neg_hits, acc) & n1, hit_set(neg_hits, acc) & n2
        row = {
            "pfam_acc": acc,
            "name": fam["name"],
            "in_family_table": fam["in_family_table"],
            "n_candidate_clusters": len(own_clusters[acc]),
            "candidate_clusters": ";".join(sorted(own_clusters[acc])),
            "n_own_members": len(members),
            "own_hits": len(own_hit),
            "own_recall_circular": f"{len(own_hit) / len(members):.3f}",
            "leakage": LEAKAGE,
        }
        for s in SLICES:
            row[f"other_positive_hits_{s}"] = other_by_slice[s]
        row.update(
            {
                "n_neg_n1": len(n1),
                "neg_hits_n1": len(h1),
                "n_clusters_neg_n1": len(n1_clusters),
                "neg_clusters_hit_n1": len({neg_cluster[a] for a in h1}),
                "n_neg_n2": len(n2),
                "neg_hits_n2": len(h2),
            }
        )
        total = 0
        for name in names:
            hit = hit_set(prot_hits[name], acc)
            total += len(hit)
            row[f"{name}_hits"] = len(hit)
            row[f"{name}_per10k"] = f"{PER * len(hit) / prot_n[name]:.2f}"
            for p in sorted(hit):
                prot_rows.append([acc, fam["name"], name, p, prot_head[name][p][:ANNOTATION_MAX]])
        row["kept_hits_total"] = total
        row["kept_n_proteins"] = kept_total
        row["kept_per10k"] = f"{PER * total / kept_total:.2f}"
        rows.append(row)
        for p in sorted(pos):
            pos_rows.append(
                [
                    acc,
                    p,
                    split[p]["cluster_id"],
                    split[p]["slice"],
                    "own" if p in members else "other",
                ]
            )
        for a in sorted(hit_set(neg_hits, acc)):
            neg_rows.append([acc, a, neg_level[a], neg_cluster[a]])
    out.mkdir(parents=True, exist_ok=True)
    header = list(rows[0]) if rows else []
    write_tsv(out / "family_specificity.tsv", header, [[r[k] for k in header] for r in rows])
    write_tsv(
        out / "positives_hits.tsv",
        ["pfam_acc", "accession", "cluster_id", "slice", "own_or_other"],
        pos_rows,
    )
    write_tsv(
        out / "negatives_hits.tsv",
        ["pfam_acc", "accession", "evidence_level", "neg_cluster"],
        neg_rows,
    )
    write_tsv(
        out / "proteome_hits.tsv",
        ["pfam_acc", "name", "proteome", "protein", "annotation"],
        prot_rows,
    )
    write_tsv(
        out / "provenance.tsv",
        ["item", "value"],
        [
            ["script", "analysis/wall_fusion_truth/pr6_family_specificity.py"],
            ["coverage_min", COVERAGE_MIN],
            ["evalue_max", EVALUE_MAX],
            ["kept_proteomes", ";".join(f"{n}={prot_n[n]}" for n in names)],
            ["kept_total", kept_total],
            ["n_families", len(rows)],
            ["tier_changes", "none"],
            ["panel_changes", "none"],
        ],
    )
    return rows


# ---- the dry-run converter -------------------------------------------------------------------------


def convert(
    specificity: Path, out: Path | None, source: str, date: str, pfam_release: str, hmmer: str
) -> list[dict]:
    """Proposed `family_measure.tsv` rows from family_specificity.tsv. Writes nothing without `out`.

    `measure_basis=descriptive_shared_sets` (positives are the family's own members chosen by the family's hits, the
    negatives a shared set), so the derived tier is `smoke` at most. `n_pos_unknown` is the number of own members: nobody
    has checked whether they are in the family's Pfam seed (the leakage is `tuned_on_truth`). The sensitivity and the
    specificity columns stay empty: no interval was measured."""
    sys.path.insert(0, str(ROOT / "src"))
    from cellsurface_sorting_hat import family_tier

    forbidden = {Path(family_tier.MEASURE_PATH).resolve(), Path(family_tier.TABLE_PATH).resolve()}
    if out is not None and Path(out).resolve() in forbidden:
        raise ValueError(
            f"{out} is a committed measurement or family table; this script never writes them "
            "(a measurement row needs the owner's sign-off and its own PR)"
        )
    existing = {r["pfam_acc"] for r in read_tsv(FAMILY_MEASURE)}
    rows = []
    for r in read_tsv(specificity):
        n_own = int(r["n_own_members"])
        kept = {
            k[: -len("_hits")]: int(r[k])
            for k in r
            if k.endswith("_hits") and not k.startswith(("own", "neg", "other", "kept"))
        }
        row = dict.fromkeys(family_tier.MEASURE_COLUMNS, "")
        row.update(
            {
                "pfam_acc": r["pfam_acc"],
                "name": r["name"],
                "measure_basis": "descriptive_shared_sets",
                "denominator": "family_own_kiss_members",
                "n_pos": n_own,
                "pos_hits": int(r["own_hits"]),
                "pos_denominator": n_own,
                "coverage": f"{int(r['own_hits']) / n_own:.3f}",
                "n_neg": int(r["n_neg_n1"]),
                "n_clusters_pos": int(r["n_candidate_clusters"]),
                "n_clusters_neg": int(r["n_clusters_neg_n1"]),
                "neg_hits": int(r["neg_hits_n1"]),
                "neg_hits_n2": int(r["neg_hits_n2"]),
                "neg_clusters_hit": int(r["neg_clusters_hit_n1"]),
                "negative_basis": "curated",
                "n_pos_unknown": n_own,
                "source": source,
                "date": date,
                "pfam_release": pfam_release,
                "hmmer_version": hmmer,
                "notes": (
                    f"own-member hits {r['own_hits']} of {n_own} (circular: the family was selected because it hits "
                    f"them; leakage {r['leakage']}; not a sensitivity); N1 hits {r['neg_hits_n1']} in "
                    f"{r['neg_clusters_hit_n1']} of {r['n_clusters_neg_n1']} clusters; N2 hits {r['neg_hits_n2']}; "
                    f"kept proteomes {r['kept_hits_total']} hits in {r['kept_n_proteins']} proteins "
                    f"({r['kept_per10k']} per 10,000: "
                    + ", ".join(f"{k} {v}" for k, v in sorted(kept.items()))
                    + "); lists are agent-made and not owner-reviewed; proposal only"
                ),
            }
        )
        row = {
            k: str(v) for k, v in row.items()
        }  # the derivation reads strings, like the committed file
        row["_tier"] = family_tier.derive_tier(row)
        row["_already_measured"] = "yes" if r["pfam_acc"] in existing else "no"
        rows.append(row)
    if out is not None:
        family_tier.write_measure(
            out, [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
        )
    return rows


# ---- command line ----------------------------------------------------------------------------------


def parse_proteome_arg(text: str) -> tuple[str, Path, Path]:
    try:
        name, dom, fasta = text.split(":")
    except ValueError:
        raise argparse.ArgumentTypeError("--proteome NAME:DOMTBL:FASTA") from None
    return name, Path(dom), Path(fasta)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("candidates")
    c.add_argument("--scan", type=Path, default=SCAN_RECORD)
    c.add_argument("--missed", type=Path, default=MISSED)
    c.add_argument("--family-table", type=Path, default=FAMILY_TABLE)
    c.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("summarize")
    s.add_argument("--candidates", type=Path, required=True)
    s.add_argument("--members-split", type=Path, default=MEMBERS_SPLIT)
    s.add_argument("--neg-dir", type=Path, default=NEG_DIR)
    s.add_argument("--domtbl-pos", type=Path, required=True)
    s.add_argument("--domtbl-neg", type=Path, required=True)
    s.add_argument("--proteome", type=parse_proteome_arg, action="append", required=True)
    s.add_argument("--inputs", type=Path, default=INPUTS_TSV)
    s.add_argument("--out", type=Path, required=True)
    v = sub.add_parser("convert")
    v.add_argument("--specificity", type=Path, required=True)
    v.add_argument("--out", type=Path, help="without it: a dry run that prints the derived tiers")
    v.add_argument("--source", required=True)
    v.add_argument("--date", required=True)
    v.add_argument("--pfam-release", default="38.2")
    v.add_argument("--hmmer-version", default="HMMER 3.4")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "candidates":
            sel = write_candidates(a.scan, a.missed, a.out, a.family_table)
            print(
                f"{len(sel['clusters'])} candidate clusters, {len(sel['families'])} candidate families"
            )
        elif a.cmd == "summarize":
            counts = kept_proteome_counts(a.inputs)
            given = {n: (d, f) for n, d, f in a.proteome}
            if sorted(given) != sorted(KEPT_PROTEOMES):
                raise ValueError(
                    f"the proteomes must be exactly {list(KEPT_PROTEOMES)} (one per species)"
                )
            rows = summarize(
                a.candidates,
                a.members_split,
                a.neg_dir,
                a.domtbl_pos,
                a.domtbl_neg,
                {n: given[n] for n in KEPT_PROTEOMES},
                a.out,
                counts,
            )
            print(f"{len(rows)} families summarized")
        else:
            rows = convert(a.specificity, a.out, a.source, a.date, a.pfam_release, a.hmmer_version)
            for r in rows:
                print(
                    f"{r['pfam_acc']}\t{r['name']}\ttier={r['_tier']}\talready_measured={r['_already_measured']}"
                )
            print(
                "dry run: nothing written"
                if a.out is None
                else f"wrote {a.out} (a proposal; no committed file changed)"
            )
    except (ValueError, OSError) as err:
        print(f"pr6_family_specificity: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
