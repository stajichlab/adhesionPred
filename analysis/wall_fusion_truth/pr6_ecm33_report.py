#!/usr/bin/env python3
"""PR 6 step 6.2 for ecm33: the held-out report (recall per held-out clade, the comparators side by side, the false hits).

DRAFT, not run on HPCC. Stdlib only; `python3 -I`. It reads the tables of the build job and writes tables. It writes no call, no status,
no tier, no sequence, and it does not touch data/sorting_hat/family_measure.tsv or family_table.tsv (a custom model is report-only,
decision 6). Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md (sections 4.2, 4.6, 4.7, 5).

It REFUSES to run unless the freeze file is frozen and the cutoff files record that very file (sha256): scoring comes after the freeze.

What n means. The units are fine clades; a unit's members are UniRef50 representatives that a phmmer search of Q5AGC4 (E 1e-5) found.
Members of a clade are related and are not independent clusters, and "truth" is defined by the same search that supplies the training
members. So every cell is `descriptive`: counts and a recall fraction, no interval, no `estimated`.

Inputs (--results DIR, one sub-directory per unit written by pr6_ecm33_build.sbatch)
  <unit>/manifest.json      the training manifest (status buildable | not_buildable)
  <unit>/cutoff.json        the custom model's cutoff (pr6_ecm33_build.py cutoff) and its freeze sha256
  <unit>/relaxed_cutoff.json  the relaxed-GA Pfam cutoff of PF12454 (the same rule, the same data)
  <unit>/hits_custom.tsv    id, score of every held-out member hit by the custom model at or above the floor
  <unit>/hits_pfam_ga.tsv   id of the held-out members hit by PF12454 at its GA (hmmsearch --cut_ga)
  <unit>/hits_pfam_relaxed.tsv  id, score of every held-out member hit by PF12454 at or above the floor
  <unit>/hits_iter.tsv      id, evalue of the held-out members found by the jackhmmer iteration seeded by the training alignment
  <unit>/neg_hits.tsv       accession, score of the shared negatives hit by the custom model
  <unit>/proteome_hits.tsv  proteome, id, score of the proteome proteins hit by the custom model
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

ITER_E = 1e-5  # the inclusion E-value of the searches (fixed before any result)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def load_frozen(path: Path) -> dict:
    if not Path(path).is_file():
        raise SystemExit(
            f"FATAL: freeze file {path} is missing; nothing is reported before the freeze"
        )
    fz = json.load(open(path))
    if fz.get("frozen") is not True:
        raise SystemExit(f"FATAL: freeze file {path} is not frozen")
    for k in ("min_u", "extra_calls_per_10000", "hard_negative_rate", "kept_proteomes"):
        if k not in fz:
            raise SystemExit(f"FATAL: freeze file lacks {k}")
    return fz


def read_tsv(path: Path) -> list[dict]:
    if not Path(path).is_file():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def recovered_by_score(rows: list[dict], cutoff: float) -> set[str]:
    return {r["id"] for r in rows if float(r["score"]) >= cutoff}


def recovered_by_evalue(rows: list[dict], e: float = ITER_E) -> set[str]:
    return {r["id"] for r in rows if float(r["evalue"]) <= e}


def unit_row(
    unit: str, held: list[str], custom: set, pfam: set, relaxed: set, it: set, min_u: int
) -> dict:
    """One held-out unit. u = held-out members that EVERY comparator misses (Pfam GA, relaxed-GA Pfam, iteration); the custom model
    adds recall only there (spec 4.6). The unit is testable when u >= min_u; the freeze's rule is at least ceil(u/2) of them
    recovered. Every cell is descriptive."""
    h = set(held)
    comparators = (pfam | relaxed | it) & h
    missed_by_all = h - comparators
    custom_h = custom & h
    new = custom_h & missed_by_all
    u = len(missed_by_all)
    need = math.ceil(u / 2) if u else 0
    return {
        "unit": unit,
        "n_held_out": len(h),
        "custom_hit": len(custom_h),
        "pfam_ga_hit": len(pfam & h),
        "pfam_relaxed_hit": len(relaxed & h),
        "iteration_hit": len(it & h),
        "any_comparator_hit": len(comparators),
        "u_missed_by_every_comparator": u,
        "custom_recovers_of_u": len(new),
        "min_u_freeze": min_u,
        "testable_u": "yes" if u >= min_u else "no",
        "ceil_u_half": need,
        "rule_ok": ("yes" if len(new) >= need else "no") if u >= min_u else "not_testable",
        "recall_custom": f"{len(custom_h)}/{len(h)}",
        "status": "descriptive",
    }


UNIT_HEADER = [
    "unit",
    "n_held_out",
    "custom_hit",
    "pfam_ga_hit",
    "pfam_relaxed_hit",
    "iteration_hit",
    "any_comparator_hit",
    "u_missed_by_every_comparator",
    "custom_recovers_of_u",
    "min_u_freeze",
    "testable_u",
    "ceil_u_half",
    "rule_ok",
    "recall_custom",
    "status",
]


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r[h] for h in header])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--freeze", required=True, type=Path)
    ap.add_argument("--units", required=True, type=Path, help="units.tsv of pr6_ecm33_members.py")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    fz = load_frozen(a.freeze)
    fz_sha = sha256_file(a.freeze)
    units = [u["clade_fine"] for u in read_tsv(a.units) if u["status"] == "unit"]
    rows, bad, extra_rows, hn_rows = [], [], [], []
    for unit in units:
        d = a.results / unit
        man_path = d / "manifest.json"
        if not man_path.is_file():
            bad.append(f"{unit}: manifest.json is missing")
            continue
        man = json.load(open(man_path))
        if man["status"] != "buildable":
            rows.append(
                {
                    **{k: "" for k in UNIT_HEADER},
                    "unit": unit,
                    "n_held_out": man["n_held_out"],
                    "rule_ok": "not_buildable",
                    "status": "descriptive",
                }
            )
            continue
        cut = json.load(open(d / "cutoff.json"))
        rel = json.load(open(d / "relaxed_cutoff.json"))
        for name, c in (("cutoff.json", cut), ("relaxed_cutoff.json", rel)):
            if c.get("freeze_sha256") != fz_sha:
                bad.append(
                    f"{unit}/{name} was not made with this freeze file (sha256 differs): scoring before or without the freeze"
                )
        if cut.get("positives_used") is not False:
            bad.append(f"{unit}: cutoff.json does not state that no positive score was used")
        custom = recovered_by_score(read_tsv(d / "hits_custom.tsv"), float(cut["cutoff"]))
        relaxed = recovered_by_score(read_tsv(d / "hits_pfam_relaxed.tsv"), float(rel["cutoff"]))
        pfam = {r["id"] for r in read_tsv(d / "hits_pfam_ga.tsv")}
        it = recovered_by_evalue(read_tsv(d / "hits_iter.tsv"))
        rows.append(
            unit_row(unit, man["held_out_ids"], custom, pfam, relaxed, it, int(fz["min_u"]))
        )
        extra_rows.append(
            {
                "unit": unit,
                "cutoff": cut["cutoff"],
                "extra_calls_total": cut["extra_calls_total"],
                "extra_call_limit_total": f"{cut['extra_call_limit_total']:.4f}",
                "within_limit": "yes"
                if cut["extra_calls_total"] <= cut["extra_call_limit_total"]
                else "no",
                "by_proteome": ";".join(
                    f"{p}={n}" for p, n in sorted(cut["extra_calls_by_proteome"].items())
                ),
            }
        )
        for g, rate in sorted(cut["hard_negative_call_rates"].items()):
            hn_rows.append(
                {
                    "unit": unit,
                    "group": g,
                    "call_rate": f"{rate:.4f}",
                    "frozen_limit": cut["hard_negative_rate_limit"],
                    "within_limit": "yes" if rate <= cut["hard_negative_rate_limit"] else "no",
                }
            )
    if bad:
        for b in bad:
            print("FATAL:", b, file=sys.stderr)
        return 1
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "recall_by_unit.tsv", UNIT_HEADER, rows)
    write_tsv(
        a.out / "extra_calls.tsv",
        [
            "unit",
            "cutoff",
            "extra_calls_total",
            "extra_call_limit_total",
            "within_limit",
            "by_proteome",
        ],
        extra_rows,
    )
    write_tsv(
        a.out / "hard_negative_calls.tsv",
        ["unit", "group", "call_rate", "frozen_limit", "within_limit"],
        hn_rows,
    )
    (a.out / "summary.txt").write_text(
        "\n".join(
            [
                "PR 6 ecm33 held-out report (DRAFT): every cell is descriptive; members are related UniRef50 representatives, not independent clusters.",
                f"freeze sha256 {fz_sha}",
                f"units reported {len(rows)}; testable by u (u >= {fz['min_u']}): {sum(1 for r in rows if r['testable_u'] == 'yes')}",
                "no call, no status, no tier, no family_measure row is written; a custom model is report-only (decision 6)",
            ]
        )
        + "\n"
    )
    print((a.out / "summary.txt").read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
