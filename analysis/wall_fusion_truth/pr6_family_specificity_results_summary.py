"""Combine the two PR 6 family-specificity runs into one table (stdlib, offline, deterministic).

Reads `<record>/run1` and `<record>/run2` (tables copied back from HPCC), the committed map of Kiss positives to
proteome ids (`kiss_positive_proteome_map.tsv`, see `pr6_kiss_proteome_map.py`) and `members_split.tsv`, and writes
`<record>/summary.tsv`: one row per family with the counts of `family_specificity.tsv` plus derived columns.

Every proteome hit is put in exactly one of three groups, by the id its own proteome uses:
- `own_in_proteomes`: the hit is an own member of the family (a Kiss positive of the family's candidate cluster);
- `other_kiss_positive_in_proteomes`: the hit is a Kiss positive of another cluster;
- `non_kiss_hits`: the hit is not a Kiss positive at all.
Non-member = not an own member: `kept_hits_excl_own` = other Kiss positives + non-Kiss hits. The limit (at most 5
per 10,000 proteins, hydrophobin precedent) is read on `kept_hits_excl_own`; `non_kiss_per10k` leaves the other
Kiss positives out as well and is shown beside it. `over_limit_excl_own` is `yes` when `kept_per10k_excl_own` is
above 5 (16 hits of 32,491 are 4.92; 17 are 5.23). `own_mapped_not_hit` counts own members that the map places in a
kept proteome but that the family model did not hit there (a sequence difference or a secondary id; the proteome
FASTA is on HPCC only, so this cannot be checked here). Nothing here is a sensitivity or a specificity.

Usage: python3 -I analysis/wall_fusion_truth/pr6_family_specificity_results_summary.py [record_dir]
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pr6_kiss_proteome_map as pm  # noqa: E402

RECORD = ROOT / "docs/reports/data/pr6_family_specificity_results_2026-10-24"
MEMBERS = ROOT / "docs/reports/data/kiss_cluster_split_2026-10-22/members_split.tsv"
MAP_NAME = "kiss_positive_proteome_map.tsv"
PER = 10_000
LIMIT_PER10K = 5.0
COLUMNS = [
    "run",
    "pfam_acc",
    "name",
    "in_family_table",
    "n_own_members",
    "own_hits",
    "neg_hits_n1",
    "neg_hits_n2",
    "other_positive_hits_adhesion",
    "other_positive_hits_wall",
    "other_positive_hits_fusion",
    "kept_hits_total",
    "kept_per10k",
    "own_in_proteomes",
    "kept_hits_excl_own",
    "kept_per10k_excl_own",
    "other_kiss_positive_in_proteomes",
    "non_kiss_hits",
    "non_kiss_per10k",
    "own_mapped_not_hit",
    "over_limit_excl_own",
]


def own_accessions(family: dict, members: list[dict]) -> set[str]:
    """Kiss positives of the family's candidate clusters (`candidate_clusters` is `;`-separated)."""
    clusters = set(family["candidate_clusters"].split(";"))
    return {
        m["accession"] for m in members if m["class"] == "positive" and m["cluster_id"] in clusters
    }


def _key(proteome: str, protein: str) -> tuple[str, str]:
    return (proteome, pm.normalize_proteome_id(proteome, protein))


def summarize(record: Path = RECORD, members_path: Path = MEMBERS) -> list[dict]:
    members = pm.read_tsv(members_path)
    lookup = pm.lookup_from_map(pm.read_tsv(record / MAP_NAME))
    out = []
    for run in ("run1", "run2"):
        hits = pm.read_tsv(record / run / "proteome_hits.tsv")
        for r in pm.read_tsv(record / run / "family_specificity.tsv"):
            own = own_accessions(r, members)
            counts = {"own": 0, "other_kiss_positive": 0, "non_kiss": 0}
            seen_own = set()
            for h in hits:
                if h["pfam_acc"] != r["pfam_acc"]:
                    continue
                kind = pm.classify_hit(h["proteome"], h["protein"], own, lookup)
                counts[kind] += 1
                if kind == "own":
                    seen_own.add(lookup[_key(h["proteome"], h["protein"])])
            if sum(counts.values()) != int(r["kept_hits_total"]):
                raise ValueError(
                    f"{run} {r['pfam_acc']}: hit rows do not add up to kept_hits_total"
                )
            n_prot = int(r["kept_n_proteins"])
            excl = counts["other_kiss_positive"] + counts["non_kiss"]
            mapped_own = {a for a in lookup.values() if a in own}
            row = {k: r[k] for k in COLUMNS if k in r}
            row.update(
                run=run,
                own_in_proteomes=counts["own"],
                kept_hits_excl_own=excl,
                kept_per10k_excl_own=f"{PER * excl / n_prot:.2f}",
                other_kiss_positive_in_proteomes=counts["other_kiss_positive"],
                non_kiss_hits=counts["non_kiss"],
                non_kiss_per10k=f"{PER * counts['non_kiss'] / n_prot:.2f}",
                own_mapped_not_hit=len(mapped_own - seen_own),
                over_limit_excl_own="yes" if PER * excl / n_prot > LIMIT_PER10K else "no",
            )
            out.append(row)
    out.sort(key=lambda x: (x["run"], -float(x["kept_per10k"]), x["pfam_acc"]))
    return out


def main(argv: list[str]) -> None:
    record = Path(argv[1]) if len(argv) > 1 else RECORD
    rows = summarize(record)
    with (record / "summary.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in COLUMNS])
    print(f"wrote {record / 'summary.tsv'}: {len(rows)} families")


if __name__ == "__main__":
    main(sys.argv)
