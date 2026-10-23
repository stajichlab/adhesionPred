#!/usr/bin/env python3
"""PR 6 for hyd1_like: the held-out report (DRAFT, not run on HPCC). Stdlib only; `python3 -I`.

Writes tables only: no call, no status, no tier, no sequence, nothing in data/sorting_hat. REFUSES to run unless the freeze is frozen and every
cutoff.json records that very freeze file (sha256). ecm33 code is imported (loaded by path), not copied.

Comparators (spec section 5; the freeze names them). u is counted twice, because the ecm33 run taught that a comparator seeded by the training
alignment ties the custom model by construction (u = 0 in all 11 units):
  INDEPENDENT of the training alignment (they enter u_indep):
    pfam_ga           hydrophobin-class Pfam models at GA (hits_pfam_ga.tsv, ids)
    relaxed_panel     data/sorting_hat/hydrophobin_relaxed.hmm, `--nobias --cut_ga` (hits_relaxed_panel.tsv, ids)
    phmmer_seed       phmmer of the single Kiss protein, one pass (hits_phmmer_seed.tsv, id, evalue; the file is absent for the Kiss unit, whose
                      held-out member IS the seed: the comparator is then `seed_in_unit` and does not count)
  CEILING (reported, never in u_indep): jackhmmer_seed, jackhmmer -N 3 from the single Kiss protein (hits_jackhmmer_seed.tsv). The members were
    found by that very iteration, so it recovers them by construction; it measures what a search from the Kiss protein alone reaches.
A unit is testable when u_indep >= min_u (freeze); the custom model must then recover ceil(u/2) of those u members. Whatever the numbers, no cell
is `estimated`: every cell is `descriptive`, except that the Kiss unit row (leave-cluster-out of KS0017, n = 1 cluster) carries the label `partial`.

Rows (column `line`): headline (cluster30) and descriptive (none) for clade and cluster units; the Kiss unit in both columns (kiss_unit_*); the
`external` line (curated hydrophobins related to the Kiss protein, minus the training members; tiny, descriptive). EVERY row prints the exact two-sided
95 percent Clopper-Pearson interval for each recall and for recovered-of-u. A wide interval is not a result.
The cutoff columns print the extra calls and the hard-negative rate twice: with the panel_ignored ids left out (the headline rule, Q8 `ignore`) and
with them counted (columns *_panel_counted, evaluated at the headline cutoff; cutoff_if_panel_counted is the cutoff that rule alone would give).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


e33 = _load("pr6_ecm33_report")
ev = _load("pr6_ecm33_eval_report")  # the exact stdlib Clopper-Pearson function
clopper_pearson = ev.clopper_pearson
fmt_ci = ev.fmt_ci
read_tsv = e33.read_tsv
write_tsv = e33.write_tsv
sha256_file = e33.sha256_file
recovered_by_score = e33.recovered_by_score
ITER_E = e33.ITER_E


def load_frozen(path: Path) -> dict:
    fz = e33.load_frozen(path)
    for k in ("comparators", "training_rule", "heldout_truth"):
        if k not in fz:
            raise SystemExit(f"FATAL: freeze file lacks {k}")
    return fz


def unit_row(
    unit: str,
    held: list[str],
    custom: set,
    pfam: set,
    relaxed: set,
    seed: set | None,
    jack: set,
    min_u: int,
    kiss: bool,
) -> dict:
    """seed is None when the seed is the held-out member (the Kiss unit): that comparator does not count."""
    h = set(held)
    indep = (pfam | relaxed | (seed or set())) & h
    u_set = h - indep
    u = len(u_set)
    ceiling = (indep | jack) & h
    u_all = len(h - ceiling)
    new = custom & u_set
    need = math.ceil(u / 2) if u else 0
    testable = u >= min_u
    row = {
        "unit": unit,
        "n_held_out": len(h),
        "custom_hit": len(custom & h),
        "pfam_ga_hit": len(pfam & h),
        "relaxed_panel_hit": len(relaxed & h),
        "phmmer_seed_hit": "seed_in_unit" if seed is None else len(seed & h),
        "jackhmmer_seed_hit_ceiling": len(jack & h),
        "u_indep": u,
        "u_with_ceiling": u_all,
        "custom_recovers_of_u": len(new),
        "min_u_freeze": min_u,
        "testable_u": "yes" if testable else "no",
        "ceil_u_half": need,
        "rule_ok": ("yes" if len(new) >= need else "no") if testable else "not_testable",
        "recall_custom": f"{len(custom & h)}/{len(h)}",
        "label": "partial" if kiss else "descriptive",
    }
    row.update(intervals(row, len(h), len(u_set), len(new)))
    return row


CI_NOTE = (
    "Intervals are the exact two-sided 95 percent Clopper-Pearson interval of k of n (stdlib, same function as the ecm33 evaluation). "
    "A wide interval is not a result: with these n it is the honest statement of how little each row shows."
)


def intervals(row: dict, n_held: int, n_u: int, n_new: int) -> dict:
    """CI columns for every recall (custom, each comparator, the ceiling) and for recovered-of-u. NA when the denominator is 0 or the
    comparator was not run (seed_in_unit)."""
    out = {"recall_custom_ci95": fmt_ci(row["custom_hit"], n_held)}
    for col, key in (
        ("pfam_ga_hit", "pfam_ga_ci95"),
        ("relaxed_panel_hit", "relaxed_panel_ci95"),
        ("phmmer_seed_hit", "phmmer_seed_ci95"),
        ("jackhmmer_seed_hit_ceiling", "jackhmmer_seed_ceiling_ci95"),
    ):
        out[key] = "NA" if row[col] == "seed_in_unit" else fmt_ci(int(row[col]), n_held)
    out["recovered_of_u_ci95"] = fmt_ci(n_new, n_u)
    return out


HEADER = list(unit_row("x", [], set(), set(), set(), set(), set(), 3, False)) + [
    "line",
    "independence_mode",
    "cutoff",
    "extra_calls_total",
    "extra_call_limit_total",
    "hard_negative_max_rate",
    "panel_ignored_n",
    "extra_calls_total_panel_counted",
    "hard_negative_max_rate_panel_counted",
    "cutoff_if_panel_counted",
]


def line_of(unit: str, mode: str, fz: dict) -> str:
    """Which line of the report a row is: the headline column, the descriptive column, the external line or the Kiss unit."""
    if unit == "external":
        return "external_descriptive"
    if unit == "clade:" + fz["heldout_truth"]["kiss_unit_clade"]:
        return "kiss_unit_" + ("headline" if mode == fz["independence_mode"] else "descriptive")
    return "headline" if mode == fz["independence_mode"] else "descriptive"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--results", required=True, type=Path, help="<unit-dir>/<mode>/ from pr6_hyd1_build.sbatch"
    )
    ap.add_argument("--freeze", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    fz = load_frozen(a.freeze)
    fz_sha = sha256_file(a.freeze)
    rows, bad = [], []
    for man_path in sorted(a.results.glob("*/*/manifest.json")):
        d = man_path.parent
        man = json.load(open(man_path))
        unit, mode = man["unit"], man["mode"]
        if unit == "all":  # the all-members model holds nothing out: no held-out row
            continue
        if man["status"] != "buildable":
            rows.append(
                {
                    **{k: "" for k in HEADER},
                    "unit": unit,
                    "independence_mode": mode,
                    "line": line_of(unit, mode, fz),
                    "n_held_out": man["n_held_out"],
                    "rule_ok": "not_buildable",
                    "label": "descriptive",
                }
            )
            continue
        cut = json.load(open(d / "cutoff.json"))
        if cut.get("freeze_sha256") != fz_sha:
            bad.append(
                f"{unit}/{mode}: cutoff.json was not made with this freeze file (scoring before or without the freeze)"
            )
        if cut.get("positives_used") is not False:
            bad.append(f"{unit}/{mode}: cutoff.json does not state that no positive score was used")
        seed_file = d / "hits_phmmer_seed.tsv"
        seed = (
            {r["id"] for r in read_tsv(seed_file) if float(r["evalue"]) <= ITER_E}
            if seed_file.is_file()
            else None
        )
        r = unit_row(
            unit,
            man["held_out_ids"],
            recovered_by_score(read_tsv(d / "hits_custom.tsv"), float(cut["cutoff"])),
            {x["id"] for x in read_tsv(d / "hits_pfam_ga.tsv")},
            {x["id"] for x in read_tsv(d / "hits_relaxed_panel.tsv")},
            seed,
            {
                x["id"]
                for x in read_tsv(d / "hits_jackhmmer_seed.tsv")
                if float(x["evalue"]) <= ITER_E
            },
            int(fz["min_u"]),
            unit == "clade:" + fz["heldout_truth"]["kiss_unit_clade"],
        )
        alt = cut.get("panel_ignored_counted", {})
        at = alt.get("at_headline_cutoff", {})
        r.update(
            line=line_of(unit, mode, fz),
            independence_mode=mode,
            cutoff=cut["cutoff"],
            extra_calls_total=cut["extra_calls_total"],
            extra_call_limit_total=f"{cut['extra_call_limit_total']:.4f}",
            hard_negative_max_rate=cut.get("hard_negative_max_rate", ""),
            panel_ignored_n=len(cut.get("panel_ignored", [])),
            extra_calls_total_panel_counted=at.get("extra_calls_total", ""),
            hard_negative_max_rate_panel_counted=at.get("hard_negative_max_rate", ""),
            cutoff_if_panel_counted=alt.get("cutoff_if_counted", ""),
        )
        rows.append(r)
    if bad:
        for b in bad:
            print("FATAL:", b, file=sys.stderr)
        return 1
    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "recall_by_unit.tsv", HEADER, rows)
    (a.out / "notes.txt").write_text(CI_NOTE + "\n")
    print(
        f"{len(rows)} unit/mode rows written; every cell descriptive (the Kiss unit row is labelled partial)"
    )
    print(CI_NOTE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
