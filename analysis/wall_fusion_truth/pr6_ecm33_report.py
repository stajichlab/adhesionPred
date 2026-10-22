#!/usr/bin/env python3
"""PR 6 step 6.2 for ecm33: the held-out report (recall per held-out clade, the comparators side by side, the false hits).

DRAFT, not run on HPCC. Stdlib only; `python3 -I`. It reads the tables of the build job and writes tables. It writes no call, no status,
no tier, no sequence, and it does not touch data/sorting_hat/family_measure.tsv or family_table.tsv (a custom model is report-only,
decision 6). Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md (sections 4.2, 4.6, 4.7, 5).

It REFUSES to run unless the freeze file is frozen and the cutoff files record that very file (sha256): scoring comes after the freeze.

What n means. The units are fine clades; a unit's members are UniRef50 representatives that a phmmer search of Q5AGC4 (E 1e-5) found.
Members of a clade are related and are not independent clusters, and "truth" is defined by the same search that supplies the training
members. So every cell is `descriptive`: counts and a recall fraction, no interval, no `estimated`.

Independence modes (owner decision 2026-10-22). Every held-out unit is built and scored under TWO rules, side by side: the HEADLINE
`cluster30` (spec 4.4: no training member in the same 30 percent identity / 50 percent coverage cluster as a held-out member) and a
second, DESCRIPTIVE column `none` (only the held-out members themselves are removed). The modes come from the freeze file. When the
headline rule removed no neighbour, the build job does not repeat the unit (SAME_AS_HEADLINE) and the second column repeats the headline.

P38248 (S. cerevisiae Ecm33, a shared hard negative that is a true family member) is left out of the hard-negative rate by the freeze
(`negatives_ignored`) and is reported on its own line (p38248_line.tsv): its custom score, the cutoff and whether it is hit.

Inputs (--results DIR, <unit>/<mode>/ written by pr6_ecm33_build.sbatch; `all` is not scored here)
  <unit>/<mode>/manifest.json      the training manifest (status buildable | not_buildable)
  <unit>/<mode>/cutoff.json        the custom model's cutoff (pr6_ecm33_build.py cutoff) and its freeze sha256
  <unit>/<mode>/relaxed_cutoff.json  the relaxed-GA Pfam cutoff of PF12454 (the same rule, the same data)
  <unit>/<mode>/hits_custom.tsv    id, score of every held-out member hit by the custom model at or above the floor
  <unit>/<mode>/hits_pfam_ga.tsv   id of the held-out members hit by PF12454 at its GA (hmmsearch --cut_ga)
  <unit>/<mode>/hits_pfam_relaxed.tsv  id, score of every held-out member hit by PF12454 at or above the floor
  <unit>/<mode>/hits_iter.tsv      id, evalue of the held-out members found by the jackhmmer iteration seeded by the training alignment
  <unit>/<mode>/neg_hits.tsv       accession, score of the shared negatives hit by the custom model
  <unit>/<mode>/proteome_hits.tsv  proteome, id, score of the proteome proteins hit by the custom model
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
    "independence_mode",
    "role",
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
WIDE_HEADER = [
    "unit",
    "n_held_out",
    "n_excluded_neighbours_headline",
    "n_training_headline",
    "n_training_second",
    "custom_hit_headline",
    "custom_hit_second",
    "u_headline",
    "u_second",
    "rule_ok_headline",
    "rule_ok_second",
    "second_is_same_training_set",
]


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r[h] for h in header])


P38248 = "P38248"


def score_dir(
    unit: str, mode: str, role: str, d: Path, fz: dict, fz_sha: str, bad: list[str]
) -> tuple[dict, dict | None, list, list, dict | None]:
    """(recall row, manifest, extra-call rows, hard-negative rows, P38248 line) of one (unit, mode) directory."""
    man_path = d / "manifest.json"
    if not man_path.is_file():
        bad.append(f"{unit}/{mode}: manifest.json is missing")
        return {}, None, [], [], None
    man = json.load(open(man_path))
    base = {"independence_mode": mode, "role": role}
    if man["status"] != "buildable":
        row = {
            **{k: "" for k in UNIT_HEADER},
            **base,
            "unit": unit,
            "n_held_out": man["n_held_out"],
            "rule_ok": "not_buildable",
            "status": "descriptive",
        }
        return row, man, [], [], None
    cut = json.load(open(d / "cutoff.json"))
    rel = json.load(open(d / "relaxed_cutoff.json"))
    for name, c in (("cutoff.json", cut), ("relaxed_cutoff.json", rel)):
        if c.get("freeze_sha256") != fz_sha:
            bad.append(
                f"{unit}/{mode}/{name} was not made with this freeze file (sha256 differs): scoring before or without the freeze"
            )
    if cut.get("positives_used") is not False:
        bad.append(f"{unit}/{mode}: cutoff.json does not state that no positive score was used")
    custom = recovered_by_score(read_tsv(d / "hits_custom.tsv"), float(cut["cutoff"]))
    relaxed = recovered_by_score(read_tsv(d / "hits_pfam_relaxed.tsv"), float(rel["cutoff"]))
    pfam = {r["id"] for r in read_tsv(d / "hits_pfam_ga.tsv")}
    it = recovered_by_evalue(read_tsv(d / "hits_iter.tsv"))
    row = {
        **unit_row(unit, man["held_out_ids"], custom, pfam, relaxed, it, int(fz["min_u"])),
        **base,
    }
    extra = [
        {
            "unit": unit,
            "independence_mode": mode,
            "role": role,
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
    ]
    hn = [
        {
            "unit": unit,
            "independence_mode": mode,
            "role": role,
            "group": g,
            "call_rate": f"{rate:.4f}",
            "frozen_limit": cut["hard_negative_rate_limit"],
            "within_limit": "yes" if rate <= cut["hard_negative_rate_limit"] else "no",
        }
        for g, rate in sorted(cut["hard_negative_call_rates"].items())
    ]
    # P38248 on its own line: it is outside the hard-negative rate (freeze negatives_ignored) but a true family member
    neg = {r["accession"]: float(r["score"]) for r in read_tsv(d / "neg_hits.tsv")}
    negr = {r["accession"]: float(r["score"]) for r in read_tsv(d / "neg_hits_relaxed.tsv")}
    line = None
    if P38248 in fz["negatives_ignored"]:
        sc = neg.get(P38248)
        rs = negr.get(P38248)
        line = {
            "unit": unit,
            "independence_mode": mode,
            "role": role,
            "accession": P38248,
            "custom_score": "" if sc is None else sc,
            "custom_cutoff": cut["cutoff"],
            "custom_hit": "yes" if sc is not None and sc >= float(cut["cutoff"]) else "no",
            "pfam_relaxed_score": "" if rs is None else rs,
            "pfam_relaxed_cutoff": rel["cutoff"],
            "pfam_relaxed_hit": "yes" if rs is not None and rs >= float(rel["cutoff"]) else "no",
            "note": "true Ecm33; left out of the hard-negative rate by the freeze (negatives_ignored); the shared negative set is not changed",
        }
    return row, man, extra, hn, line


EXTRA_HEADER = [
    "unit",
    "independence_mode",
    "role",
    "cutoff",
    "extra_calls_total",
    "extra_call_limit_total",
    "within_limit",
    "by_proteome",
]
HN_HEADER = [
    "unit",
    "independence_mode",
    "role",
    "group",
    "call_rate",
    "frozen_limit",
    "within_limit",
]
P38248_HEADER = [
    "unit",
    "independence_mode",
    "role",
    "accession",
    "custom_score",
    "custom_cutoff",
    "custom_hit",
    "pfam_relaxed_score",
    "pfam_relaxed_cutoff",
    "pfam_relaxed_hit",
    "note",
]


def wide_row(
    unit: str, head: dict, second: dict, head_man: dict | None, second_man: dict | None, same: bool
) -> dict:
    return {
        "unit": unit,
        "n_held_out": head.get("n_held_out", ""),
        "n_excluded_neighbours_headline": head_man["n_excluded_neighbours"] if head_man else "",
        "n_training_headline": head_man["n_training"] if head_man else "",
        "n_training_second": second_man["n_training"] if second_man else "",
        "custom_hit_headline": head.get("custom_hit", ""),
        "custom_hit_second": second.get("custom_hit", ""),
        "u_headline": head.get("u_missed_by_every_comparator", ""),
        "u_second": second.get("u_missed_by_every_comparator", ""),
        "rule_ok_headline": head.get("rule_ok", ""),
        "rule_ok_second": second.get("rule_ok", ""),
        "second_is_same_training_set": "yes" if same else "no",
    }


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
    head_mode, second_mode = fz["independence_mode"], fz["independence_second_column"]
    units = [u["clade_fine"] for u in read_tsv(a.units) if u["status"] == "unit"]
    rows, wide, extra_rows, hn_rows, p_rows, bad = [], [], [], [], [], []
    for unit in units:
        head, head_man, ex, hn, pl = score_dir(
            unit, head_mode, "headline", a.results / unit / head_mode, fz, fz_sha, bad
        )
        if not head:
            continue
        rows.append(head)
        extra_rows += ex
        hn_rows += hn
        p_rows += [pl] if pl else []
        sdir = a.results / unit / second_mode
        same = (sdir / "SAME_AS_HEADLINE").is_file()
        if same:
            # the headline rule removed no neighbour: the second build is the same training set; the numbers are the headline's
            second = {**head, "independence_mode": second_mode, "role": "descriptive"}
            second_man = head_man
            rows.append(second)
            extra_rows += [
                {**r, "independence_mode": second_mode, "role": "descriptive"} for r in ex
            ]
            hn_rows += [{**r, "independence_mode": second_mode, "role": "descriptive"} for r in hn]
            p_rows += (
                [{**pl, "independence_mode": second_mode, "role": "descriptive"}] if pl else []
            )
        else:
            second, second_man, ex2, hn2, pl2 = score_dir(
                unit, second_mode, "descriptive", sdir, fz, fz_sha, bad
            )
            if not second:
                continue
            rows.append(second)
            extra_rows += ex2
            hn_rows += hn2
            p_rows += [pl2] if pl2 else []
        wide.append(wide_row(unit, head, second, head_man, second_man, same))
    if bad:
        for b in bad:
            print("FATAL:", b, file=sys.stderr)
        return 1
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "recall_by_unit.tsv", UNIT_HEADER, rows)
    write_tsv(a.out / "side_by_side.tsv", WIDE_HEADER, wide)
    write_tsv(a.out / "extra_calls.tsv", EXTRA_HEADER, extra_rows)
    write_tsv(a.out / "hard_negative_calls.tsv", HN_HEADER, hn_rows)
    write_tsv(a.out / "p38248_line.tsv", P38248_HEADER, p_rows)
    n_testable = sum(1 for r in rows if r["role"] == "headline" and r["testable_u"] == "yes")
    (a.out / "summary.txt").write_text(
        "\n".join(
            [
                "PR 6 ecm33 held-out report (DRAFT): every cell is descriptive; members are related UniRef50 representatives, not independent clusters.",
                f"freeze sha256 {fz_sha}",
                f"independence: headline {head_mode}; second, descriptive column {second_mode} (side_by_side.tsv)",
                f"units reported {len(wide)}; testable by u (u >= {fz['min_u']}) under the headline rule: {n_testable}",
                "P38248 (a true Ecm33 among the shared negatives) is outside the hard-negative rate and has its own line (p38248_line.tsv)",
                "no call, no status, no tier, no family_measure row is written; a custom model is report-only (decision 6)",
            ]
        )
        + "\n"
    )
    print((a.out / "summary.txt").read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
