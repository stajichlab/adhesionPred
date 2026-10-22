#!/usr/bin/env python3
"""PR 6 ecm33 pre-freeze check: which of the 161 shared hard negatives are members of the Ecm33 family? (DRAFT, nothing run on HPCC)

The freeze `pr6_ecm33_freeze.json` leaves the accessions of `negatives_ignored` out of the hard-negative rate, because a correct Ecm33
model should hit them: they are true family members that the shared set labels as negatives for the ADHESION call. The freeze lists
P38248 (S. cerevisiae Ecm33). Q4WNS8 (A. fumigatus ecm33) is a second one, and a member under another name cannot be excluded by
reading names. This script turns the output of the pre-check job (pr6_ecm33_negatives_precheck.sbatch: `phmmer` of the Kiss Ecm33
query Q5AGC4 against the 161 shared negatives) into the complete list, with the SAME definition of 'family member' that the freeze uses
for proteome proteins (`labelled_rule`): phmmer of Q5AGC4 finds it at E <= 1e-5.

    python3 -I pr6_ecm33_negatives_precheck.py summarize --strict phmmer_E1e-5.tbl --weak phmmer_E1e-3.tbl \\
        --controls data/controls/hard-negatives-shared/controls.tsv \\
        --negatives-fasta data/controls/hard-negatives-shared/sequences.faa --already-ignored P38248 --out OUT_DIR

Files written (ids, genes, scores and lengths only; no sequence):
  negatives_family_table.tsv   one row per negative hit at E <= 1e-5: accession, gene, species, stratum, evalue, score, length, listed_in_freeze
  weak_E1e-3.tsv              DESCRIPTIVE ONLY: hits at 1e-5 < E <= 1e-3 (never part of the proposed list)
  proposed_negatives_ignored.json   the exact one-line JSON edit for the freeze (the union of the hits and the already listed accessions)
  summary.txt                 counts, the per-stratum numbers, and a LOUD line if an expected member is not found

The script changes nothing in the freeze. The owner decides the list: it is a proposal made from the data.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

STRICT_E = 1e-5
WEAK_E = 1e-3
EXPECTED = (
    "P38248",
    "Q4WNS8",
)  # the two true Ecm33 proteins known among the negatives; the check must find both


def parse_tblout(path: Path) -> list[dict]:
    """phmmer `--tblout`: 18 columns then the description. The `# [ok]` trailer must be present (the search finished)."""
    lines = path.read_text().splitlines()
    if not lines or lines[-1].strip() != "# [ok]":
        raise ValueError(f"{path.name}: no '# [ok]' trailer (the search did not finish)")
    out = []
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split(None, 18)
        if len(f) < 18:
            raise ValueError(f"{path.name}: short line: {ln[:80]}")
        out.append({"target": f[0], "query": f[2], "evalue": float(f[4]), "score": float(f[5])})
    return out


def read_controls(path: Path) -> dict[str, dict]:
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}


def fasta_lengths(path: Path) -> dict[str, int]:
    """{first token of the header: residue count}. Only the lengths are kept."""
    out: dict[str, int] = {}
    name = None
    for ln in path.read_text().splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            if name in out:
                raise ValueError(f"duplicate header {name} in {path.name}")
            out[name] = 0
        elif name is not None:
            out[name] += len(ln.strip())
    return out


def build(
    strict: list[dict],
    weak: list[dict],
    controls: dict[str, dict],
    lengths: dict[str, int],
    already: list[str],
    expected: tuple[str, ...] = EXPECTED,
) -> dict:
    """The table, the weak table, the proposed list and the loud lines. Pure function (tested on fixtures)."""
    unknown = sorted({h["target"] for h in strict + weak if h["target"] not in controls})
    strict_ok = [h for h in strict if h["evalue"] <= STRICT_E]
    members = {h["target"]: h for h in strict_ok if h["target"] in controls}
    rows = []
    for acc in sorted(members, key=lambda a: (members[a]["evalue"], a)):
        c = controls[acc]
        rows.append(
            {
                "accession": acc,
                "gene": c.get("gene", ""),
                "species": c.get("species", ""),
                "stratum": c.get("stratum", ""),
                "evalue": members[acc]["evalue"],
                "score": members[acc]["score"],
                "length": lengths.get(acc, ""),
                "listed_in_freeze": "yes" if acc in already else "no",
            }
        )
    weak_rows = []
    for h in sorted(weak, key=lambda x: (x["evalue"], x["target"])):
        if (
            STRICT_E < h["evalue"] <= WEAK_E
            and h["target"] in controls
            and h["target"] not in members
        ):
            c = controls[h["target"]]
            weak_rows.append(
                {
                    "accession": h["target"],
                    "gene": c.get("gene", ""),
                    "species": c.get("species", ""),
                    "stratum": c.get("stratum", ""),
                    "evalue": h["evalue"],
                    "score": h["score"],
                    "length": lengths.get(h["target"], ""),
                }
            )
    proposed = sorted(set(already) | set(members))
    loud = [
        f"LOUD: expected family member {a} is NOT found at E <= {STRICT_E:g}"
        for a in expected
        if a not in members
    ]
    loud += [f"LOUD: hit accession {a} is not in controls.tsv" for a in unknown]
    loud += [
        f"LOUD: already listed in the freeze but not found by the search: {a}"
        for a in sorted(set(already) - set(members))
    ]
    strata: dict[str, list[str]] = {}
    for a in members:
        strata.setdefault(controls[a].get("stratum", ""), []).append(a)
    return {
        "rows": rows,
        "weak_rows": weak_rows,
        "proposed": proposed,
        "new": sorted(set(members) - set(already)),
        "loud": loud,
        "strata": {k: sorted(v) for k, v in sorted(strata.items())},
    }


def json_edit(proposed: list[str]) -> str:
    """The one line to put in the freeze file."""
    return '"negatives_ignored": ' + json.dumps(proposed) + ","


TABLE_HEADER = [
    "accession",
    "gene",
    "species",
    "stratum",
    "evalue",
    "score",
    "length",
    "listed_in_freeze",
]
WEAK_HEADER = ["accession", "gene", "species", "stratum", "evalue", "score", "length"]


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r[h] for h in header])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summarize")
    s.add_argument("--strict", required=True, type=Path, help="phmmer --tblout at E 1e-5")
    s.add_argument(
        "--weak", required=True, type=Path, help="phmmer --tblout at E 1e-3 (descriptive)"
    )
    s.add_argument("--controls", required=True, type=Path)
    s.add_argument("--negatives-fasta", required=True, type=Path)
    s.add_argument(
        "--already-ignored",
        nargs="*",
        default=["P38248"],
        help="the accessions the freeze lists now",
    )
    s.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    if a.out.exists() and any(a.out.iterdir()):
        print(f"FATAL: {a.out} is not empty; nothing is overwritten", file=sys.stderr)
        return 1
    res = build(
        parse_tblout(a.strict),
        parse_tblout(a.weak),
        read_controls(a.controls),
        fasta_lengths(a.negatives_fasta),
        list(a.already_ignored),
    )
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "negatives_family_table.tsv", TABLE_HEADER, res["rows"])
    write_tsv(a.out / "weak_E1e-3.tsv", WEAK_HEADER, res["weak_rows"])
    edit = json_edit(res["proposed"])
    (a.out / "proposed_negatives_ignored.json").write_text(edit + "\n")
    lines = [
        "PR 6 ecm33 pre-freeze check (DRAFT): which shared hard negatives are Ecm33 family members?",
        f"definition: phmmer of Q5AGC4 against the shared negatives finds the protein at E <= {STRICT_E:g} (the freeze's labelled_rule)",
        f"family members among the negatives: {len(res['rows'])}",
        f"already listed in the freeze: {', '.join(a.already_ignored) or '(none)'}",
        f"NEW (not yet listed): {', '.join(res['new']) or '(none)'}",
        "per stratum: " + "; ".join(f"{k}: {', '.join(v)}" for k, v in res["strata"].items()),
        f"weak hits (1e-5 < E <= 1e-3), descriptive only, never in the proposed list: {len(res['weak_rows'])}",
        "proposed edit of the freeze (the owner decides; this script changes nothing):",
        "  " + edit,
    ] + res["loud"]
    (a.out / "summary.txt").write_text("\n".join(lines) + "\n")
    print((a.out / "summary.txt").read_text())
    for ln in res["loud"]:
        print(ln, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
