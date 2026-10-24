"""Write `kiss_positive_proteome_map.tsv` from the committed UniProt export (stdlib, offline, deterministic).

Inputs: `members_split.tsv` (the 132 Kiss positives) and the UniProt TSV of the 112 non-JGI positives, downloaded once
from rest.uniprot.org (release and query in `kiss_positive_proteome_map.provenance.txt`). The id rules are in
`pr6_kiss_proteome_map.py`.

Usage: python3 -I analysis/wall_fusion_truth/pr6_build_kiss_positive_proteome_map.py [record_dir]
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pr6_kiss_proteome_map as pm  # noqa: E402

RECORD = ROOT / "docs/reports/data/pr6_family_specificity_results_2026-10-24"
MEMBERS = ROOT / "docs/reports/data/kiss_cluster_split_2026-10-22/members_split.tsv"
UNIPROT_NAME = "kiss_positives_uniprot_2026_03.tsv"
MAP_NAME = "kiss_positive_proteome_map.tsv"


def build(record: Path = RECORD, members: Path = MEMBERS) -> list[dict]:
    return pm.build_map(pm.read_tsv(members), pm.read_tsv(record / UNIPROT_NAME))


def main(argv: list[str]) -> None:
    record = Path(argv[1]) if len(argv) > 1 else RECORD
    rows = build(record)
    pm.write_map(record / MAP_NAME, rows)
    counts = {s: sum(r["status"] == s for r in rows) for s in pm.STATUSES}
    print(f"wrote {record / MAP_NAME}: {len(rows)} positives, {counts}")


if __name__ == "__main__":
    main(sys.argv)
