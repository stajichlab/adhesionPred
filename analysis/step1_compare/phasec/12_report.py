#!/usr/bin/env python3
"""Phase C step 5: report.md from metrics.json and findings.json (Phase C spec 4, E3).

Reads $STEP1_WORKDIR/phasec/metrics.json and findings.json (SHA-256 checked against
evaluate_run.json). Writes report.md and report_run.json to $STEP1_WORKDIR/phasec/.

Every number in report.md is a value of metrics.json or findings.json, printed by `fmt`
(integers with thousands separators, other numbers with 3 decimals). Before it writes, the
script checks this with `unsupported_numbers` and stops if a number has no source. The report
quotes the findings fields; it adds no free-text claim about accuracy. The literature test set
has no negatives, so its tables give recall only (spec 3.2). The fixed caveat lines (CAVEATS)
name known limits of the method.

STOP (exit 2, no output): stale inputs; a number in the report without a source.
"""

import argparse
import re
import sys
from pathlib import Path

import evalio
import paths
import runinfo

OUTPUT_NAMES = ("report.md", "report_run.json")
NUMBER = re.compile(r"(?<![\w.\-])-?\d[\d,]*(?:\.\d+)?(?![\w])")
HEADLINE = (
    "recall",
    "precision",
    "fpr",
    "roc_auc",
    "pr_auc",
    "precision_at_recall_0.8",
    "precision_at_recall_0.9",
    "recall_at_fpr_0.01",
)
LITERATURE_HEADLINE = ("recall",)  # positives only: recall is the one defined metric
CAVEATS = (
    "The pooled S1:all ROC-AUC and PR-AUC rank decision values from five fold models. Each "
    "fold model has its own score scale, so the pooled values mix scales (review M-3; the "
    "per-fold diagnostic is not computed).",
    "The ML decision threshold and the Platt scaling are fitted on the inner out-of-fold "
    "predictions and then applied to the model refitted on all training rows. The refitted "
    "model can have another score scale, so binary calls and calibrated probabilities can "
    "shift (review M-4; not measured).",
    "The report-number check is set membership: each printed number equals some value in "
    "metrics.json or findings.json. The check does not show that a number is in the correct "
    "cell (review M-7).",
)
CANNOT_SHOW = (
    "The whole-proteome prevalence of surface proteins. It is not measured; the prevalence "
    "table uses assumed values (ruling C-11).",
    "Basidiomycota performance. Cneo_H99_GOA and Umay_MYCMD are smoke tests, not the Q8 "
    "validation (owner decision on the Basidiomycota question).",
    "Whether SignalP under-calls C. immitis RS proteins. No truth set shows it; the literature "
    "rows give recall only.",
    "Performance on P-gpi proteins. P-gpi is a list until curated_gpi.tsv has literature rows "
    "(R-B).",
    "Variance from refitting the models. The intervals describe sampling of the test set only.",
)


def fmt(x) -> str:
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, int):
        return f"{x:,}"
    return f"{x:.3f}"


def ci(m: dict | None) -> str:
    if not m or m.get("value") is None:
        return "n/a"
    return f"{fmt(m['value'])} [{fmt(m['lo'])}, {fmt(m['hi'])}]"


def _leaves(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield k
            yield from _leaves(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _leaves(v)
    else:
        yield obj


def allowed_numbers(*objs) -> set[str]:
    out = set()
    for obj in objs:
        for x in _leaves(obj):
            if isinstance(x, bool) or x is None:
                continue
            if isinstance(x, str):
                out |= set(NUMBER.findall(x))
                continue
            if isinstance(x, int):
                out |= {f"{x:,}", str(x)}
            else:
                out.add(f"{x:.3f}")
    return out


def unsupported_numbers(text: str, allowed: set[str]) -> list[str]:
    return [t for t in NUMBER.findall(text) if t not in allowed]


def table(header, rows) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out + [""]


def render(m: dict, f: dict) -> str:
    st = m["settings"]
    cands, variants = st["candidates"], st["variants"]
    lines = [
        "# Step 1 Phase C: rule versus ML",
        "",
        f"Intervals: {fmt(st['ci_level'])}% cluster bootstrap, {fmt(st['n_resamples'])} "
        f"resamples, seed {st['seed']}. The intervals describe sampling of the test set only. "
        "The models are fitted once per split and are not refitted inside the bootstrap.",
        "Headline truth: direct evidence (`homology_only == no`). All non-IEA truth is beside it.",
        f"A test set is an estimate when the recall half-width is at most "
        f"{fmt(st['estimate_half_width'])} for R2 and every ML candidate (V-go, direct truth) "
        f"and the set has at least {fmt(st['estimate_min_direct_positives'])} direct-evidence "
        "positives; otherwise it is a smoke test.",
        "",
        "## Caveats",
        "",
        *[f"- {x}" for x in CAVEATS],
        "",
        "## Estimate or smoke test",
        "",
    ]
    rows = []
    for name, ts in m["test_sets"].items():
        n = ts["truth"]["direct"]["n"]["all"]
        rows.append([name, ts["label"], fmt(n["pos"]), fmt(n["neg"]), fmt(ts["floor_met"]),
                     fmt(ts["max_recall_half_width"]), fmt(ts["max_fpr_half_width"]),
                     ", ".join(ts["zero_width_recall_interval"]) or "none"])  # fmt: skip
    lines += table(["Test set", "Label", "Positives", "Negatives", "Count floor met",
                    "Max recall half-width", "Max FPR half-width", "Zero-width recall interval"],
                   rows)  # fmt: skip
    lines += ["A zero-width interval (every resample gives the same recall, for example 1 of 1) "
              "meets the half-width rule but carries no information about precision of the "
              "estimate. The count floor keeps such a small set a smoke test.", ""]  # fmt: skip
    for name, ts in m["test_sets"].items():
        lit = ts["kind"] == "literature"
        heads = LITERATURE_HEADLINE if lit else HEADLINE
        lines += [f"## {name} ({ts['label']})", ""]
        if lit:
            lines += [
                "Literature rows have positives only: the tables give recall only (Phase C spec). "
                "Context: the V-kw training of S1 and S2 keeps the C. immitis and "
                "C. posadasii T-c rows (counts in the context section).",
                "",
            ]
        for truth in ("direct", "all"):
            block = ts["truth"][truth]
            n = block["n"]["all"]
            lines += [f"Truth `{truth}`: {fmt(n['pos'])} positives, {fmt(n['neg'])} negatives.", ""]
            rows = []
            for v in variants:
                for c in cands:
                    mm = block["metrics"]["all"][v][c]
                    rows.append([c, v, ts["label"], *[ci(mm.get(k)) for k in heads]])
            lines += table(["Candidate", "Variant", "Label", *heads], rows)
        block = ts["truth"]["direct"]
        cols = LITERATURE_HEADLINE if lit else ("recall", "fpr", "roc_auc")
        rows = []
        for stratum, n in block["n"].items():
            for c in cands:
                mm = block["metrics"][stratum]["V-go"][c]
                rows.append([stratum, c, ts["label"], fmt(n["pos"]), fmt(n["neg"]),
                             *[ci(mm.get(k)) for k in cols]])  # fmt: skip
        lines += ["Strata (direct truth, V-go):", ""]
        lines += table(["Stratum", "Candidate", "Label", "Positives", "Negatives", *cols], rows)
        cols = LITERATURE_HEADLINE if lit else ("recall", "fpr", "roc_auc", "pr_auc")
        rows = []
        for c in cands:
            eff = block["variant_effect"]["all"][c]
            rows.append([c, ts["label"], *[ci(eff.get(k)) for k in cols]])
        lines += ["Variant effect, V-kw minus V-go (paired intervals; test-set sampling only):", ""]
        lines += table(["Candidate", "Label", *cols], rows)
        amb = ts["ambiguous"]
        lines += [f"Ambiguous genes: {fmt(amb['n'])} (score distribution only); with "
                  f"high-throughput-only internal evidence: {fmt(ts['ambiguous_htp_only']['n'])}.", ""]  # fmt: skip
        if lit:
            continue  # no negatives: no comparison at the rule's FPR, no prevalence table
        rows = []
        for v in variants:
            for c, vr in block["vs_rule"].get(v, {}).items():
                rows.append(
                    [
                        c,
                        v,
                        ts["label"],
                        ci(vr["precision_at_rule_recall"]),
                        ci(vr["recall_at_rule_fpr"]),
                    ]
                )
        if rows:
            lines += [
                "ML against the rule R2 (precision at the rule's recall, recall at the rule's FPR):",
                "",
            ]
            lines += table(
                ["Candidate", "Variant", "Label", "Precision at R2 recall", "Recall at R2 FPR"],
                rows,
            )
        rows = []
        for c in cands:
            pv = block["prevalence"]["V-go"][c]
            rows.append([c, ts["label"], *[ci(pv[k]) for k in pv]])
        prev = [fmt(p) for p in st["prevalences"]]
        lines += ["Precision at assumed prevalence (assumed, not measured), V-go:", ""]
        lines += table(["Candidate", "Label", *prev], rows)
        if "calibration" in ts:
            rows = []
            for v, per in ts["calibration"].items():
                for c, cal in per.items():
                    rows.append([c, v, ts["label"], fmt(cal["brier"])])
            lines += [
                "Calibration after Platt scaling on inner out-of-fold predictions (Brier score):",
                "",
            ]
            lines += table(["Candidate", "Variant", "Label", "Brier"], rows)
            rows = []
            for v, per in ts["calibration"].items():
                for c, cal in per.items():
                    cells = [f"{fmt(b['frac_pos'])} ({fmt(b['n'])})" for b in cal["bins"]]
                    rows.append([c, v, ts["label"], *cells])
            bins = next(iter(next(iter(ts["calibration"].values())).values()))["bins"]
            head = [f"{fmt(b['lo'])} to {fmt(b['hi'])}" for b in bins]
            lines += ["Reliability: observed positive fraction (count) per bin of the calibrated "
                      "probability:", ""]  # fmt: skip
            lines += table(["Candidate", "Variant", "Label", *head], rows)
    lines += ["## Agreement of R2 with ML (proteome_calls.tsv.gz holds the protein IDs)", ""]
    rows = []
    for set_id, per in m["agreement"]["proteomes"].items():
        for v, cc in per.items():
            for c, k in cc.items():
                rows.append([set_id, c, v, fmt(k["rule1_ml1"]), fmt(k["rule1_ml0"]),
                             fmt(k["rule0_ml1"]), fmt(k["rule0_ml0"])])  # fmt: skip
    lines += table(["Set", "ML", "Variant", "both", "rule only", "ML only", "neither"], rows)
    lines += ["The C. immitis RS row is read as held-out, but V-kw training of S1 and S2 keeps "
              "the C. immitis and C. posadasii T-c rows (context section).", ""]  # fmt: skip
    rows = []
    for cls, per in m["agreement"]["truth"].items():
        for v, cc in per.items():
            for c, k in cc.items():
                rows.append([cls, c, v, fmt(k["rule1_ml1"]), fmt(k["rule1_ml0"]),
                             fmt(k["rule0_ml1"]), fmt(k["rule0_ml0"])])  # fmt: skip
    lines += ["Agreement on the GO truth rows (all sources), by class:", ""]
    lines += table(["Class", "ML", "Variant", "both", "rule only", "ML only", "neither"], rows)
    lines += ["## Named panel", ""]
    rows = []
    for p in m["named_panel"]:
        calls = [fmt(p["calls"].get(f"{c}|V-go")) if p["found"] else "n/a" for c in cands]
        rows.append([p["name"], p["id"], fmt(p["found"]), p.get("label", ""), p.get("class", ""),
                     p.get("score_source", ""), *calls])  # fmt: skip
    lines += ["Calls under V-go (1 = called); scores and V-kw calls are in metrics.json.", ""]
    lines += table(["Protein", "ID", "Found", "Label", "Class", "Score source", *cands], rows)
    a, b, c3 = f["a_b1_not_saturated"], f["b_ml_beats_b1_and_r2_on_nsec_s1"], f["c_same_under_s2"]
    lines += ["## Findings (quoted from findings.json)", ""]
    lines += [f"- (a) B1 ROC-AUC below {fmt(a['threshold'])} on S1 and S2 (direct truth, V-go): "
              f"holds = {fmt(a['holds'])}."]  # fmt: skip
    lines += [f"  - {k}: {fmt(v)}" for k, v in a["b1_roc_auc"].items()]
    lines += ["Findings (b) and (c) are read per ML candidate; no headline ML candidate is named "
              "(ruling C-13). \"Any\" means at least one listed candidate.", ""]  # fmt: skip
    lines += [f"- (b) Any ML candidate beats B1 and R2 on N-sec FPR at the recall of R2 "
              f"({b['test_set']}): holds for {', '.join(b['holds_for']) or 'no ML candidate'}."]  # fmt: skip
    for cand, entry in b["candidates"].items():
        lines += [f"  - {cand}: B1 minus ML {ci(entry['B1'])}; R2 minus ML {ci(entry['R2'])}"]
    lines += [
        f"- (c) Any ML candidate, the same on every S2 test set: holds for "
        f"{', '.join(c3['holds_for']) or 'no ML candidate'}.",
        "",
    ]
    lines += ["## Context", ""]
    for key, n in m["context"]["onygenales_tc_rows_in_vkw_training"].items():
        lines += [f"- {key}: {fmt(n)} C. immitis and C. posadasii T-c rows in V-kw training."]
    lines += ["", "## What the data cannot show", ""]
    lines += [f"- {x}" for x in CANNOT_SHOW]
    return "\n".join(lines) + "\n"


def run(work: Path, arguments=()):
    out = evalio.out_dir(work)
    ev = evalio.require_current(out, "evaluate_run.json", ("metrics.json", "findings.json"), "11")
    m = evalio.read_json(out / "metrics.json")
    f = evalio.read_json(out / "findings.json")
    text = render(m, f)
    bad = unsupported_numbers(text, allowed_numbers(m, f))
    if bad:
        raise evalio.StopError(f"report numbers without a source in metrics.json: {bad[:5]}")
    log = {
        "all_sources": ev["all_sources"],
        "truth_set_sha256": ev["truth_set_sha256"],
        "input_sha256": {n: ev["outputs_sha256"][n] for n in ("metrics.json", "findings.json")},
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    evalio.write_outputs(
        out, {"report.md": lambda p: Path(p).write_text(text)}, "report_run.json", log
    )
    return text


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        text = run(work, list(argv) if argv is not None else sys.argv[1:])
    except (evalio.StopError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(f"report.md: {len(text.splitlines())} lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
