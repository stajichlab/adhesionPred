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

STOP (exit 2): stale inputs; a number in the report without a source; an input that has a null
or another wrong type where the report needs a value. A STOP removes an older report.md and
report_run.json, so that a stale report never stays beside a failed run.
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
STRATA_COLUMNS = ("recall", "fpr", "roc_auc")
EFFECT_COLUMNS = ("recall", "fpr", "roc_auc", "pr_auc")
DECISION_COLUMNS = ("recall", "fpr", "roc_auc")
SOURCES = ("oof", "final", "in_sample")
# Tables whose cells tests/step1_compare/test_phasec_report.py compares with metrics.json and
# findings.json one by one. The third caveat names exactly this list.
CELL_CHECKED = (
    "decision table",
    "metric tables of each test set (direct and all truth)",
    "strata tables",
    "variant-effect tables",
    "ML-against-rule tables",
    "prevalence tables",
    "Brier tables",
    "reliability tables",
    "agreement tables",
    "findings (b) and (c) tables",
    "named-panel call and score tables",
    "fitted-settings table",
)
OPERATING_POINT = (
    "Operating point: ML candidates at their own fitted call threshold, rules (R0, R1, R2) at "
    "their fitted setting, B0 and B1 at their own fitted call."
)
# final review I-2. R2 calls a subset of the proteins that R0 calls (tests/step1_compare/
# test_phasec_rules.py checks the subset and shows that J of R2 can still exceed J of R0).
RULE_SUBSET = (
    "R2 is R0 restricted by g and t: R0 calls every protein that R2 calls. So the recall and the "
    "FPR of R2 cannot exceed those of R0; its recall minus FPR can. When the fitted g and t are "
    "the loosest cell of the grid (g weakly and the lowest t; see Fitted settings), the rule was "
    "set to the loosest value the grid allows."
)
EDGE_SENTENCE = (
    "A setting at a grid edge means the data may prefer a value outside the grid; the report "
    "cannot say."
)
FITTED_HEADER = ["Unit", "Variant", "Candidate", "g", "t", "J", "C", "H variant", "Threshold",
                 "Platt a", "Platt b", "Convergence warnings", "At grid edge"]  # fmt: skip
EDGE_PARAMETERS = ("C", "g", "t")
ZERO_WIDTH = (
    "Zero-width recall interval lists the candidates whose recall interval has the same lower "
    "and upper bound."
)
CAVEATS = (
    "The pooled S1:all ROC-AUC and PR-AUC rank decision values from five fold models. Each "
    "fold model has its own score scale, so the pooled values mix scales (plan review item "
    "M-3; the per-fold diagnostic is not computed).",
    "The ML decision threshold and the Platt scaling are fitted on the inner out-of-fold "
    "predictions and then applied to the model refitted on all training rows. The refitted "
    "model can have another score scale, so binary calls and calibrated probabilities can "
    "shift (plan review item M-4; not measured).",
    "The report-number check is a typo guard: each printed number equals some value in "
    "metrics.json or findings.json. It does not show that a number is in the correct cell. "
    "Cell-level tests compare every cell of these tables with the JSON files: "
    + "; ".join(CELL_CHECKED)
    + ". Every other number is checked only by the typo guard (plan review item M-7).",
)
CANNOT_SHOW = (
    "The whole-proteome prevalence of surface proteins. It is not measured. The prevalence "
    "table uses assumed values.",
    "Basidiomycota performance. Cneo_H99_GOA and Umay_MYCMD are smoke tests. A validation "
    "needs curated Basidiomycota truth, which does not exist yet.",
    "Whether SignalP under-calls C. immitis RS proteins. No truth set shows it; the literature "
    "rows give recall only.",
    "Performance on P-gpi proteins. P-gpi is a list until curated_gpi.tsv has literature rows "
    "(the curation of GPI-anchored proteins from the literature, not yet done).",
    "Variance from refitting the models. The intervals describe sampling of the test set only.",
)
SIGN = (
    "The difference is the comparator minus the ML candidate. A value below zero means the ML "
    "candidate has the higher N-sec false-positive rate (the baseline is better). The finding "
    "holds only when the paired {level}% interval is entirely above zero."
)
# final review M-3; {fraction} is min_defined_fraction of findings.json
DEFINED_RULE = (
    "A comparator counts as beaten only when its difference is defined in at least a fraction "
    "{fraction} of the resamples. A cell shows (n_defined k of B) when its interval rests on "
    "k of the B resamples only."
)
# (code, meaning); "{long}" is the long-protein cut-off from settings
GLOSSARY = (
    ("B0", "baseline: logistic regression on log protein length."),
    ("B1", "baseline: amino-acid fractions and log length, logistic regression."),
    ("R0", "rule: SignalP calls a signal peptide (SP)."),
    ("R1", "rule: R0 and a GPI call at or above the class g."),
    ("R2", "the rule: SP and (GPI call at or above g, or Ser+Thr fraction at or above t)."),
    ("M8", "ESM-2 8M embedding of the N-terminal window, logistic regression."),
    ("M35", "ESM-2 35M embedding of the N-terminal window, logistic regression."),
    ("M8-C", "as M8, but proteins longer than {long} aa use the C-terminal window."),
    ("M35-C", "as M35, but proteins longer than {long} aa use the C-terminal window."),
    ("H", "hybrid: embedding, SignalP probability, GPI score and Ser+Thr in one regression."),
    ("ML", "any of M8, M35, M8-C, M35-C and H."),
    ("V-go", "training variant: GO-labelled proteins only."),
    (
        "V-kw",
        "training variant: V-go plus the keyword-only (T-c) rows that pass the leakage "
        "controls.",
    ),  # fmt: skip
    ("S1", "split: homology-grouped cross-validation (out-of-fold test)."),
    ("S2", "split: leave one species out."),
    ("S3", "split: leave one clade out."),
    ("FULL", "the model fitted on the whole training pool; it scores the proteomes."),
    ("N-int", "negative stratum: internal proteins."),
    ("N-sec", "negative stratum: proteins of the secretory pathway that are not surface proteins."),
    ("PM-TM", "negative stratum: plasma-membrane proteins with a transmembrane segment."),
    ("pos", "positive protein (wall or extracellular)."),
    ("neg", "negative protein."),
    ("FPR", "false-positive rate."),
    ("J", "Youden's J: recall minus FPR."),
    ("C", "inverse regularisation strength of a logistic regression, fitted on the inner folds."),
    ("grid edge", "the first or the last value of the grid of a fitted setting."),
    ("ROC-AUC", "area under the receiver operating characteristic curve."),
    ("PR-AUC", "area under the precision-recall curve."),
    ("oof", "score from a model that did not train on the protein or its cluster."),
    ("in_sample", "score from a model that trained on the protein."),
    ("final", "score from the FULL model for a protein outside every training table."),
    ("g", "cut-off class of the GPI call, fitted on the training rows."),
    ("t", "cut-off of the Ser+Thr fraction, fitted on the training rows."),
    ("T-c", "keyword-only training tier: proteins labelled surface by UniProt keywords, not GO."),
    ("P-gpi", "plasma-membrane protein with curated GPI evidence; a list, not scored."),
    ("homology_only", "yes when every supporting GO code is a homology-transfer code."),
    ("direct", "truth of the genes with homology_only = no; the headline truth."),
    ("all", "truth of all non-IEA labelled genes; shown beside direct."),
    ("Platt scaling", "logistic fit that turns a decision value into a probability."),
    ("Brier score", "mean squared difference of probability and label; lower is better."),
    ("half-width", "half of the width of an interval: (upper bound minus lower bound) / two."),
    (
        "paired interval",
        "interval of a difference, with both terms resampled on the same " "clusters.",
    ),  # fmt: skip
    ("S1:all", "S1 test set pooled over all folds; S1:<source> is one source of it."),
    ("S2-<source>:<source>", "S2 test set: the left-out species."),
    (
        "S3-<clade>:clade",
        "S3 test set: all sources of the left-out clade; S3-<clade>:<source> " "is one source.",
    ),  # fmt: skip
    ("literature", "S3-Eurotiomycetes:literature, curated literature rows (positives only)."),
    ("estimate", "test set with enough direct-evidence positives and narrow recall intervals."),
    ("smoke test", "test set that does not meet the estimate rule; it checks the pipeline only."),
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
        # the percent form of the assumed prevalences ("1% prevalence") is derived from them
        if isinstance(obj, dict) and isinstance(obj.get("settings"), dict):
            out |= {f"{p * 100:g}" for p in obj["settings"].get("prevalences", [])}
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


def cell_metric(block: dict, stratum: str, variant: str, cand: str, key: str):
    """One metric of one key (stratum, variant, candidate, metric) of a truth block."""
    return block["metrics"].get(stratum, {}).get(variant, {}).get(cand, {}).get(key)


def cell(block: dict, stratum: str, variant: str, cand: str, key: str) -> str:
    return ci(cell_metric(block, stratum, variant, cand, key))


def table(header, rows) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out + [""]


def unit_words(key: str) -> str:
    name, _, unit = key.partition("|")
    if name == "S1":
        return f"S1 fold {unit}"
    if name == "FULL":
        return "FULL model (whole training pool)"
    return f"{name} training set"


def ci_n(m: dict, d: dict | None) -> str:
    """ci() plus `(n_defined k of B)` when the interval rests on fewer than all B resamples."""
    text = ci(d)
    n_res = m["settings"]["n_resamples"]
    if text != "n/a" and d.get("n_defined") is not None and d["n_defined"] < n_res:
        text += f" (n_defined {fmt(d['n_defined'])} of {fmt(n_res)})"
    return text


def nsec_fpr(m: dict, name: str, cand: str) -> str:
    block = m["test_sets"][name]["truth"]["direct"]["nsec_fpr_at_rule_recall"]
    return ci_n(m, block.get("V-go", {}).get("fpr", {}).get(cand))


def edge_words(entry: dict) -> str:
    """The fitted settings of one candidate that sit at a grid edge, or `no`."""
    return ", ".join(k for k, v in entry["at_grid_edge"].items() if v) or "no"


def fitted_rows(m: dict) -> list[list[str]]:
    rows = []
    for key, per in m["fitted_settings"].items():
        unit, _, variant = key.rpartition("|")
        for c, e in per.items():
            text = "n/a" if e.get("g") is None else e["g"]
            rows.append([unit, variant, c, text, fmt(e.get("t")), fmt(e.get("j")),
                         fmt(e.get("C")), e.get("h_variant") or "n/a", fmt(e.get("threshold")),
                         fmt(e.get("platt_a")), fmt(e.get("platt_b")),
                         fmt(e.get("convergence_warnings")), edge_words(e)])  # fmt: skip
    return rows


def finding_rows(m: dict, name: str, entries: dict) -> list[list[str]]:
    label = m["test_sets"][name]["label"]
    return [[name, label, c, nsec_fpr(m, name, "B1"), nsec_fpr(m, name, "R2"),
             nsec_fpr(m, name, c), ci_n(m, e["B1"]), ci_n(m, e["R2"]), fmt(e["beats_B1"]),
             fmt(e["beats_R2"]), fmt(e["holds"])] for c, e in entries.items()]  # fmt: skip


FINDING_HEADER = ["Test set", "Label", "ML candidate", "B1 N-sec FPR", "R2 N-sec FPR",
                  "ML N-sec FPR", "B1 minus ML", "R2 minus ML", "Beats B1", "Beats R2",
                  "Holds"]  # fmt: skip


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
        "positives; otherwise it is a smoke test. Every number of a test set carries the "
        "label of that set.",
        "",
        "## Caveats",
        "",
        *[f"- {x}" for x in CAVEATS],
        "",
        "## Decision table",
        "",
        "V-go, direct truth, stratum all except the N-sec column (N-sec stratum). Every candidate "
        "is listed. The report names no best candidate: a choice among ML candidates on the "
        "test data would use the test data twice.",
        OPERATING_POINT,
        RULE_SUBSET,
        "",
    ]
    rows = []
    for name, ts in m["test_sets"].items():
        block = ts["truth"]["direct"]
        for c in cands:
            rows.append([name, ts["label"], c,
                         cell(block, "all", "V-go", c, "recall"),
                         cell(block, "all", "V-go", c, "fpr"),
                         cell(block, "N-sec", "V-go", c, "fpr"),
                         cell(block, "all", "V-go", c, "roc_auc")])  # fmt: skip
    lines += table(["Test set", "Label", "Candidate", "recall", "FPR", "N-sec FPR", "ROC-AUC"],
                   rows)  # fmt: skip
    lines += ["## Estimate or smoke test", "", ZERO_WIDTH, ""]
    rows = []
    for name, ts in m["test_sets"].items():
        n = ts["truth"]["direct"]["n"]["all"]
        rows.append([name, ts["label"], fmt(n["pos"]), fmt(n["neg"]), fmt(ts["floor_met"]),
                     fmt(ts["max_recall_half_width"]), fmt(ts["max_fpr_half_width"]),
                     ", ".join(ts["zero_width_recall_interval"]) or "none"])  # fmt: skip
    lines += table(["Test set", "Label", "Positives", "Negatives", "Count floor met",
                    "Max recall half-width", "Max FPR half-width", "Zero-width recall interval"],
                   rows)  # fmt: skip
    lines += ["A zero-width interval (every resample gives the same recall) meets the half-width "
              "rule but says nothing about the exactness of the estimate. The count floor keeps "
              "such a small set a smoke test.", ""]  # fmt: skip
    lines += ["## Glossary", ""]
    lines += [f"- {code}: {text.format(long=fmt(st['long_cutoff']))}" for code, text in GLOSSARY]
    lines += [""]
    for name, ts in m["test_sets"].items():
        lit = ts["kind"] == "literature"
        heads = LITERATURE_HEADLINE if lit else HEADLINE
        lab = ts["label"]
        lines += [f"## {name} ({lab})", ""]
        if lit:
            lines += [
                "Literature rows have positives only, so the tables give recall only. "
                "Context: the V-kw training of S1 and S2 keeps the C. immitis and "
                "C. posadasii T-c rows (counts in the context section).",
                "",
            ]
        for truth in ("direct", "all"):
            block = ts["truth"][truth]
            n = block["n"]["all"]
            lines += [f"Truth `{truth}` ({lab}): {fmt(n['pos'])} positives, "
                      f"{fmt(n['neg'])} negatives.", ""]  # fmt: skip
            rows = []
            for v in variants:
                for c in cands:
                    rows.append([c, v, lab, *[cell(block, "all", v, c, k) for k in heads]])
            lines += table(["Candidate", "Variant", "Label", *heads], rows)
        block = ts["truth"]["direct"]
        cols = LITERATURE_HEADLINE if lit else STRATA_COLUMNS
        rows = []
        for stratum, n in block["n"].items():
            for c in cands:
                rows.append([stratum, c, lab, fmt(n["pos"]), fmt(n["neg"]),
                             *[cell(block, stratum, "V-go", c, k) for k in cols]])  # fmt: skip
        lines += [f"Strata (direct truth, V-go; {lab}):", ""]
        lines += table(["Stratum", "Candidate", "Label", "Positives", "Negatives", *cols], rows)
        cols = LITERATURE_HEADLINE if lit else EFFECT_COLUMNS
        rows = []
        for c in cands:
            eff = block["variant_effect"]["all"][c]
            rows.append([c, lab, *[ci(eff.get(k)) for k in cols]])
        lines += ["Variant effect, V-kw minus V-go (paired intervals; test-set sampling only):", ""]
        lines += table(["Candidate", "Label", *cols], rows)
        amb = ts["ambiguous"]
        lines += [f"Ambiguous genes ({lab}): {fmt(amb['n'])} (score distribution only); with "
                  f"high-throughput-only internal evidence: {fmt(ts['ambiguous_htp_only']['n'])}.", ""]  # fmt: skip
        if lit:
            continue  # no negatives: no comparison at the rule's FPR, no prevalence table
        rows = []
        for v in variants:
            for c, vr in block["vs_rule"].get(v, {}).items():
                rows.append([c, v, lab, ci(vr["precision_at_rule_recall"]),
                             ci(vr["recall_at_rule_fpr"])])  # fmt: skip
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
        keys = list(next(iter(block["prevalence"]["V-go"].values())))
        for c in cands:
            pv = block["prevalence"]["V-go"][c]
            rows.append([c, lab, *[ci(pv[k]) for k in keys]])
        prev = [f"{float(k) * 100:g}% prevalence" for k in keys]
        lines += [f"Precision at assumed prevalence (assumed, not measured), V-go; {lab}:", ""]
        lines += table(["Candidate", "Label", *prev], rows)
        if "calibration" in ts:
            rows = []
            for v, per in ts["calibration"].items():
                for c, cal in per.items():
                    rows.append([c, v, lab, fmt(cal["brier"])])
            lines += [
                f"Calibration after Platt scaling on inner out-of-fold predictions "
                f"(Brier score; {lab}):",
                "",
            ]
            lines += table(["Candidate", "Variant", "Label", "Brier"], rows)
            rows = []
            for v, per in ts["calibration"].items():
                for c, cal in per.items():
                    cells = [f"{fmt(b['frac_pos'])} ({fmt(b['n'])})" for b in cal["bins"]]
                    rows.append([c, v, lab, *cells])
            bins = next(iter(next(iter(ts["calibration"].values())).values()))["bins"]
            head = [f"{fmt(b['lo'])} to {fmt(b['hi'])}" for b in bins]
            lines += [f"Reliability ({lab}): observed positive fraction (count) per bin of the "
                      "calibrated probability:", ""]  # fmt: skip
            lines += table(["Candidate", "Variant", "Label", *head], rows)
    counts = m["agreement"].get("score_source_counts", {})

    def source_cells(group: str, key: str) -> list[str]:
        have = counts.get(group, {}).get(key, {})
        return [fmt(have[s]) if s in have else "none" for s in SOURCES]

    lines += ["## Agreement of R2 with ML (proteome_calls.tsv.gz holds the protein IDs)", ""]
    lines += ["The columns oof, final and in_sample count the proteins of the set by the source "
              "of their score. A count that is not listed is zero (none).", ""]  # fmt: skip
    rows = []
    for set_id, per in m["agreement"]["proteomes"].items():
        for v, cc in per.items():
            for c, k in cc.items():
                rows.append([set_id, c, v, fmt(k["rule1_ml1"]), fmt(k["rule1_ml0"]),
                             fmt(k["rule0_ml1"]), fmt(k["rule0_ml0"]),
                             *source_cells("proteomes", set_id)])  # fmt: skip
    lines += table(["Set", "ML", "Variant", "both", "rule only", "ML only", "neither", *SOURCES],
                   rows)  # fmt: skip
    lines += ["The C. immitis RS row mixes out-of-fold and final-model scores, as the counts in "
              "the oof and final columns show. It is not purely held-out, because V-kw training "
              "of S1 and S2 keeps the C. immitis and C. posadasii T-c rows (context "
              "section).", ""]  # fmt: skip
    rows = []
    for cls, per in m["agreement"]["truth"].items():
        for v, cc in per.items():
            for c, k in cc.items():
                rows.append([cls, c, v, fmt(k["rule1_ml1"]), fmt(k["rule1_ml0"]),
                             fmt(k["rule0_ml1"]), fmt(k["rule0_ml0"]),
                             *source_cells("truth", cls)])  # fmt: skip
    lines += ["Agreement on the GO truth rows (all sources), by class:", ""]
    lines += table(["Class", "ML", "Variant", "both", "rule only", "ML only", "neither", *SOURCES],
                   rows)  # fmt: skip
    lines += ["## Named panel", ""]
    lines += ["Calls and scores of the named proteins that the panel found. Each cell reads "
              "V-go / V-kw. A call of 1 means called. Rules have no score (n/a).", ""]  # fmt: skip
    found = [p for p in m["named_panel"] if p["found"]]
    rows = []
    for p in m["named_panel"]:
        rows.append([p["name"], p["id"], fmt(p["found"]), p.get("label", ""), p.get("class", ""),
                     p.get("score_source", "")])  # fmt: skip
    lines += table(["Protein", "ID", "Found", "Label", "Class", "Score source"], rows)
    for what, field in (("Calls", "calls"), ("Scores", "scores")):
        rows = [[p["name"], *[f"{fmt(p[field].get(f'{c}|V-go'))} / "
                              f"{fmt(p[field].get(f'{c}|V-kw'))}" for c in cands]]
                for p in found]  # fmt: skip
        lines += [f"{what} (V-go / V-kw):", ""]
        lines += table(["Protein", *cands], rows)
    a, b, c3 = f["a_b1_not_saturated"], f["b_ml_beats_b1_and_r2_on_nsec_s1"], f["c_same_under_s2"]
    sign = SIGN.format(level=fmt(st["ci_level"]))
    frac = f["b_ml_beats_b1_and_r2_on_nsec_s1"]["min_defined_fraction"]
    defined_rule = DEFINED_RULE.format(fraction=fmt(frac))
    lines += ["## Findings (quoted from findings.json)", ""]
    lines += ["- (a) B1 ROC-AUC is below the saturation threshold on S1 and S2 "
              f"(direct truth, V-go): holds = {fmt(a['holds'])}."]  # fmt: skip
    lines += [f"  - {k} ({m['test_sets'][k]['label']}): B1 ROC-AUC {fmt(v)}, threshold "
              f"{fmt(a['threshold'])}" for k, v in a["b1_roc_auc"].items()]  # fmt: skip
    lines += ["Findings (b) and (c) are read per ML candidate. No headline ML candidate is named, "
              "so the choice is not made on test data. \"Any\" means at least one listed "
              "candidate.", ""]  # fmt: skip
    lab_b = m["test_sets"][b["test_set"]]["label"]
    lines += [f"- (b) Any ML candidate beats B1 and R2 on N-sec FPR at the recall of R2 "
              f"({b['test_set']}, {lab_b}): holds for "
              f"{', '.join(b['holds_for']) or 'no ML candidate'}.", sign, defined_rule, ""]  # fmt: skip
    lines += table(FINDING_HEADER, finding_rows(m, b["test_set"], b["candidates"]))
    n_found = f"{fmt(c3['n_s2'])} of {fmt(c3['expected_s2'])} expected S2 test sets"
    labelled = ", ".join(f"{k} ({m['test_sets'][k]['label']})" for k in c3["test_sets"])
    lines += [f"- (c) Any ML candidate, the same on every S2 test set: holds for "
              f"{', '.join(c3['holds_for']) or 'no ML candidate'}.",
              f"  Found {n_found}: {labelled}.", sign, defined_rule, ""]  # fmt: skip
    rows = []
    for name, entry in c3["test_sets"].items():
        rows += finding_rows(m, name, entry["candidates"])
    lines += table(FINDING_HEADER, rows)
    lines += ["## Fitted settings", ""]
    lines += ["Fitted settings per unit (from scores_run.json, copied by 11). g, t and J are for "
              "the rules; C, the H variant, the threshold and Platt a and b are for the "
              "logistic regressions. metrics.json also holds the inner out-of-fold PR-AUC of "
              "each C.", ""]  # fmt: skip
    lines += table(FITTED_HEADER, fitted_rows(m))
    ec = m["grid_edge_counts"]
    lines += ["Settings at a grid edge: " + "; ".join(
        f"{k} {fmt(ec[k]['at_edge'])} of {fmt(ec[k]['n'])}" for k in EDGE_PARAMETERS) + ".",
        EDGE_SENTENCE, ""]  # fmt: skip
    lines += ["## Context", ""]
    lines += ["Number of C. immitis and C. posadasii keyword-only (T-c) rows in the V-kw training "
              "of each unit. A unit is one fitted model. S1 fold k is the model without fold k "
              "of S1. FULL is the model fitted on the whole training pool. The other names are "
              "the model fitted for that split.", ""]  # fmt: skip
    for key, n in m["context"]["onygenales_tc_rows_in_vkw_training"].items():
        unit = "row" if n == 1 else "rows"
        lines += [f"- {unit_words(key)} (key {key}): {fmt(n)} {unit}."]
    lines += ["", "Each count covers the C. immitis and C. posadasii rows together."]
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


def remove_stale(work: Path) -> None:
    """A STOP removes an older report, so a stale report never stays beside a failed run."""
    for name in OUTPUT_NAMES:
        try:
            (evalio.out_dir(work) / name).unlink(missing_ok=True)
        except OSError:
            pass


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        text = run(work, list(argv) if argv is not None else sys.argv[1:])
    except (evalio.StopError, ValueError, OSError, KeyError) as exc:
        remove_stale(work)
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    except (TypeError, AttributeError, IndexError) as exc:
        remove_stale(work)
        print(f"STOP: an input has a null or a wrong type ({type(exc).__name__}: {exc})",
              file=sys.stderr)  # fmt: skip
        return 2
    print(f"report.md: {len(text.splitlines())} lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
