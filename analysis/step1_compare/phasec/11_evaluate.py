#!/usr/bin/env python3
"""Phase C step 4: metrics with cluster-bootstrap intervals, strata, agreement, findings.

Reads from $STEP1_WORKDIR/phasec/ (SHA-256 checked against the run JSONs): eval_table.tsv.gz,
eval_literature.tsv (08), split_members.tsv.gz, clusters.tsv.gz, max_identity.tsv.gz (09),
scores.tsv.gz (10). Reads phaseb/sequence_members.tsv.gz (its SHA-256 must equal
features_run.json) and sequence_sets.tsv (--sets) for the proteome sets.

Test sets: S1:all (out-of-fold, all folds pooled) and S1:<source>; S2-<source>:<source>;
S3-<clade>:clade, S3-<clade>:<source>, S3-Eurotiomycetes:literature. Truth `direct` =
homology_only no (headline); `all` = all non-IEA labels (beside it). Literature rows have
positives only: a truth block without negatives stores precision, PR-AUC, precision at recall
and precision at the rule's recall as null (spec 3.2: recall only). The bootstrap weights (bootstrap.py) are drawn once per test set and truth and
serve every candidate, variant and stratum (paired).

Writes to $STEP1_WORKDIR/phasec/:

  metrics.json           values with 95% intervals per test set, truth, stratum, variant,
                         candidate; estimate or smoke-test label; variant effect; comparison
                         with the rule; prevalence table; calibration; agreement; named panel
  findings.json          findings (a), (b), (c) of spec 4 as booleans with their numbers
  proteome_calls.tsv.gz  one row per protein of every proteome set: call and score of every
                         candidate and variant, score_source
  evaluate_run.json      input hashes, seed, resamples, git commit, library versions

STOP (exit 2, no output): a stale input; scores that do not cover a test row; a T-c row in a
test set.
"""

import argparse
import csv
import gzip
import math
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import agreement
import bootstrap
import evalio
import findings
import manifest
import metrics
import numpy as np
import paths
import runinfo
import splits
import truth_table

OUTPUT_NAMES = ("metrics.json", "findings.json", "proteome_calls.tsv.gz", "evaluate_run.json")
METRICS_SCHEMA = "step1-phasec-metrics/1"
FINDINGS_SCHEMA = "step1-phasec-findings/1"
BINARY = ("recall", "precision", "fpr")
RECALL_LEVELS = (0.8, 0.9)
FPR_LEVEL = 0.01
SCORED = (
    "roc_auc",
    "pr_auc",
    "precision_at_recall_0.8",
    "precision_at_recall_0.9",
    "recall_at_fpr_0.01",
)
# undefined without negatives: precision is 1 (or n/a) by construction (spec 3.2, review I-1)
PRECISION_METRICS = (
    "precision",
    "pr_auc",
    "precision_at_recall_0.8",
    "precision_at_recall_0.9",
)
STRATA = (
    "all",
    "wall",
    "extracellular-only",
    "N-int",
    "N-sec",
    "PM-TM",
    "long",
    "identity_below_0.3",
)
TRUTHS = ("direct", "all")
PREVALENCES = (0.01, 0.03, 0.05, 0.10)
CALIBRATION_BINS = 10
QUANTILES = (0.1, 0.25, 0.5, 0.75, 0.9)
PROTEOME_KINDS = ("download", "site")
ONYGENALES_TC_TAXA = ("246410", "443226")
PROTEOME_BASE = (
    "set_id",
    "source_id",
    "gene_id",
    "seq_sha256",
    "class",
    "label",
    "origin",
    "score_source",
)
SCORE_COLUMNS = (
    "split_id", "fold", "variant", "candidate", "seq_sha256", "part", "score", "prob", "call",
)  # fmt: skip
LONG_CUTOFF = 1022


def num(x):
    """JSON number or None (NaN and inf become None)."""
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def summary(arr) -> dict:
    """Point value (row 0) and the bootstrap interval (rows 1..B)."""
    return {"value": num(arr[0]), **bootstrap.interval(arr[1:])}


@dataclass
class TestSet:
    name: str
    split: str
    kind: str
    rows: list


def read_scores(path: Path) -> dict:
    raw: dict = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if tuple(next(reader)) != SCORE_COLUMNS:
            raise evalio.StopError(f"{path}: unexpected header")
        for split, fold, variant, cand, h, _part, score, prob, call in reader:
            d = raw.setdefault((split, fold, variant, cand), ([], [], [], []))
            d[0].append(h)
            d[1].append(float(score) if score else math.nan)
            d[2].append(float(prob) if prob else math.nan)
            d[3].append(call == "1")
    out = {}
    for key, (hs, s, p, c) in raw.items():
        out[key] = {
            "index": {h: i for i, h in enumerate(hs)},
            "score": np.array(s),
            "prob": np.array(p),
            "call": np.array(c, dtype=bool),
        }
    return out


def test_sets(members, table, literature, identity, species_rows) -> list[TestSet]:
    lit = {r["seq_sha256"]: r for r in literature if r["seq_sha256"]}
    ident = {(r["split_id"], r["seq_sha256"]): r["below_0.3"] for r in identity}
    train_sources = sorted(r["source_id"] for r in species_rows if r["role"] == "train")
    rows_by_split: dict[str, list] = {}
    for m in members:
        if m["part"] not in ("test", "test_lit"):
            continue
        if m["origin"] == "tc":
            raise evalio.StopError(f"T-c row {m['seq_sha256']} is in test set {m['split_id']}")
        h = m["seq_sha256"]
        if m["part"] == "test_lit":
            info = {"label": "literature", "subset": "", "stratum": "literature", "d8_class": "",
                    "homology_only": "no", "internal_evidence_htp_only": "no",
                    "source_ids": "literature", "gene_ids": lit[h]["accession"],
                    "length": lit[h]["length"]}  # fmt: skip
        else:
            info = table[h]
        rows_by_split.setdefault(m["split_id"], []).append(
            {
                **{
                    k: info[k]
                    for k in (
                        "label",
                        "subset",
                        "stratum",
                        "d8_class",
                        "homology_only",
                        "internal_evidence_htp_only",
                        "source_ids",
                        "gene_ids",
                    )
                },
                "length": int(info["length"]),
                "seq_sha256": h,
                "split": m["split_id"],
                "fold": m["fold"],
                "class": m["class"],
                "cluster": m["cluster_id"],
                "part": m["part"],
                "below": ident.get((m["split_id"], h), ""),
            }  # fmt: skip
        )
    out = []
    for split in sorted(rows_by_split):
        rows = rows_by_split[split]
        go_rows = [r for r in rows if r["part"] == "test"]
        lit_rows = [r for r in rows if r["part"] == "test_lit"]
        sources = sorted({s for r in go_rows for s in r["source_ids"].split(",")})
        if split == "S1":
            out.append(TestSet("S1:all", split, "pooled", go_rows))
            for s in train_sources:
                out.append(
                    TestSet(
                        f"S1:{s}",
                        split,
                        "species",
                        [r for r in go_rows if s in r["source_ids"].split(",")],
                    )
                )
        elif split.startswith("S2-"):
            out.append(TestSet(f"{split}:{split[3:]}", split, "species", go_rows))
        else:
            out.append(TestSet(f"{split}:clade", split, "clade", go_rows))
            for s in sources:
                out.append(
                    TestSet(
                        f"{split}:{s}",
                        split,
                        "species",
                        [r for r in go_rows if s in r["source_ids"].split(",")],
                    )
                )
            if lit_rows:
                out.append(TestSet(f"{split}:literature", split, "literature", lit_rows))
    return out


def gather(scores, rows, variant: str, cand: str):
    """score, prob, call arrays for the rows (each row knows its split and fold)."""
    n = len(rows)
    s, p, c = np.full(n, math.nan), np.full(n, math.nan), np.zeros(n, dtype=bool)
    for i, r in enumerate(rows):
        d = scores.get((r["split"], r["fold"], variant, cand))
        j = None if d is None else d["index"].get(r["seq_sha256"])
        if j is None:
            raise evalio.StopError(
                f"no {cand} {variant} score for {r['seq_sha256']} in {r['split']} fold {r['fold']}"
            )
        s[i], p[i], c[i] = d["score"][j], d["prob"][j], d["call"][j]
    return s, p, c


def strata_masks(rows, y) -> dict:
    def has(col, value):
        return np.array([value in r[col].split(",") for r in rows], dtype=bool)

    out = {
        "all": np.ones(len(rows), dtype=bool),
        "wall": y & has("subset", "wall"),
        "extracellular-only": y & has("subset", "extracellular-only"),
        "N-int": ~y & has("stratum", "N-int"),
        "N-sec": ~y & has("stratum", "N-sec"),
        "PM-TM": ~y & has("stratum", "PM-TM"),
        "long": np.array([r["length"] > LONG_CUTOFF for r in rows], dtype=bool),
    }
    if any(r["below"] for r in rows):
        out["identity_below_0.3"] = np.array([r["below"] == "yes" for r in rows], dtype=bool)
    return out


def evaluate_truth(ts, truth, scores, cands, n_resamples, seed):
    rows = [
        r for r in ts.rows
        if r["class"] in ("pos", "neg") and (truth == "all" or r["homology_only"] == "no")
    ]  # fmt: skip
    y = np.array([r["class"] == "pos" for r in rows], dtype=bool)
    W = bootstrap.cluster_weights([r["cluster"] for r in rows], n_resamples,
                                  bootstrap.seed_for(seed, f"{ts.name}|{truth}"))  # fmt: skip
    Wx = np.vstack([np.ones((1, len(rows))), W])
    no_negatives = not bool((~y).any())
    masks = strata_masks(rows, y)
    data = {(v, c): gather(scores, rows, v, c) for v in evalio.VARIANTS for c in cands}
    arrays: dict = {}
    block = {"n": {}, "metrics": {}, "variant_effect": {}}
    for stratum, mask in masks.items():
        Wm, ym = Wx[:, mask], y[mask]
        block["n"][stratum] = {"pos": int(ym.sum()), "neg": int((~ym).sum())}
        per = block["metrics"][stratum] = {}
        for v in evalio.VARIANTS:
            per[v] = {}
            for c in cands:
                s, _, call = data[(v, c)]
                vals = {"recall": metrics.recall(Wm, ym, call[mask]),
                        "precision": metrics.precision(Wm, ym, call[mask]),
                        "fpr": metrics.fpr(Wm, ym, call[mask])}  # fmt: skip
                if not np.isnan(s).any():
                    vals.update(metrics.scored_summary(Wm, ym, s[mask], RECALL_LEVELS, FPR_LEVEL))
                if no_negatives:
                    for name in PRECISION_METRICS:
                        if name in vals:
                            vals[name] = np.full(len(Wm), np.nan)
                for name, arr in vals.items():
                    arrays[(stratum, v, c, name)] = arr
                per[v][c] = {name: summary(arr) for name, arr in vals.items()}
        block["variant_effect"][stratum] = {
            c: {
                name: summary(
                    arrays[(stratum, "V-kw", c, name)] - arrays[(stratum, "V-go", c, name)]
                )
                for name in per["V-go"][c]
            }
            for c in cands
        }
    block["vs_rule"], block["prevalence"], block["nsec_fpr_at_rule_recall"] = {}, {}, {}
    for v in evalio.VARIANTS:
        rule_rec = arrays.get(("all", v, "R2", "recall"))
        rule_fpr = arrays.get(("all", v, "R2", "fpr"))
        block["vs_rule"][v] = {}
        block["prevalence"][v] = {
            c: {
                f"{pi}": summary(metrics.precision_at_prevalence(
                    arrays[("all", v, c, "recall")], arrays[("all", v, c, "fpr")], pi))
                for pi in PREVALENCES
            }
            for c in cands
        }  # fmt: skip
        if rule_rec is None:
            continue
        nsec = masks["N-sec"]
        comp = {"R2": arrays[("N-sec", v, "R2", "fpr")]}
        for c in cands:
            s = data[(v, c)][0]
            if c == "R2" or np.isnan(s).any():
                continue
            if c in findings.ML:
                prec = metrics.precision_at_recall(Wx, y, s, rule_rec)
                if no_negatives:
                    prec = np.full(len(Wx), np.nan)
                block["vs_rule"][v][c] = {
                    "precision_at_rule_recall": summary(prec),
                    "recall_at_rule_fpr": summary(metrics.recall_at_fpr(Wx, y, s, rule_fpr)),
                }  # fmt: skip
            comp[c] = metrics.fpr_at_recall(Wx, y, s, rule_rec, nsec)
        block["nsec_fpr_at_rule_recall"][v] = {
            "fpr": {c: summary(a) for c, a in comp.items()},
            "diff": {
                ml: {k: summary(comp[k] - comp[ml]) for k in findings.COMPARATORS if k in comp}
                for ml in findings.ML
                if ml in comp
            },
        }
    return block, rows, y, data


def calibration(rows, y, data, cands) -> dict:
    out = {}
    for (v, c), (_, prob, _) in data.items():
        if c not in findings.ML or np.isnan(prob).any():
            continue
        bins = np.clip((prob * CALIBRATION_BINS).astype(int), 0, CALIBRATION_BINS - 1)
        out.setdefault(v, {})[c] = {
            "brier": num(metrics.brier(np.ones(len(y)), y, prob)[0]),
            "bins": [
                {"lo": k / CALIBRATION_BINS, "hi": (k + 1) / CALIBRATION_BINS,
                 "n": int((bins == k).sum()),
                 "mean_prob": num(prob[bins == k].mean()) if (bins == k).any() else None,
                 "frac_pos": num(y[bins == k].mean()) if (bins == k).any() else None}
                for k in range(CALIBRATION_BINS)
            ],
        }  # fmt: skip
    return out


def excluded_views(ts, scores, cands) -> dict:
    amb = [r for r in ts.rows if r["class"] == "excluded" and "ambiguous" in r["label"].split(",")]
    out = {}
    for name, rows in (("ambiguous", amb),
                       ("ambiguous_htp_only", [r for r in amb if r["internal_evidence_htp_only"] == "yes"])):  # fmt: skip
        view = {"n": len(rows)}
        for v in evalio.VARIANTS:
            for c in cands:
                if not rows:
                    continue
                s, _, call = gather(scores, rows, v, c)
                entry = {"call_fraction": num(call.mean())}
                if not np.isnan(s).any():
                    entry["score_quantiles"] = {f"{q}": num(np.quantile(s, q)) for q in QUANTILES}
                view.setdefault(v, {})[c] = entry
        out[name] = view
    lists = {}
    for name, test in (("pm-unresolved", lambda r: "pm-unresolved" in r["stratum"].split(",")),
                       ("P-gpi", lambda r: "P-gpi" in r["d8_class"].split(","))):  # fmt: skip
        picked = [r for r in ts.rows if test(r)]
        entries = []
        for r in picked:
            calls = {}
            for v in evalio.VARIANTS:
                for c in cands:
                    calls[f"{c}|{v}"] = int(gather(scores, [r], v, c)[2][0])
            entries.append(
                {"gene_ids": r["gene_ids"], "seq_sha256": r["seq_sha256"], "calls": calls}
            )
        lists[name] = entries
    out["lists"] = lists
    return out


def lookup_tables(members):
    s1 = {m["seq_sha256"]: m["fold"] for m in members
          if m["split_id"] == "S1" and m["part"] in ("test", "test_tc")}  # fmt: skip
    full_train = {m["seq_sha256"] for m in members
                  if m["split_id"] == "FULL" and m["part"] in ("train", "train_tc")}  # fmt: skip
    return s1, full_train


def hash_score(scores, s1, h, v, c):
    key = ("S1", s1[h], v, c) if h in s1 else ("FULL", "0", v, c)
    d = scores[key]
    j = d["index"][h]
    return d["score"][j], d["call"][j]


def proteome_rows(seq_members, set_ids, table, scores, s1, full_train, cands):
    rows = []
    for m in seq_members:
        if m["set_id"] not in set_ids:
            continue
        h = m["seq_sha256"]
        t = table.get(h, {})
        row = {"set_id": m["set_id"], "source_id": m["source_id"], "gene_id": m["gene_id"],
               "seq_sha256": h, "class": t.get("class", ""), "label": t.get("label", ""),
               "origin": t.get("origin", ""),
               "score_source": agreement.score_source(h, s1, full_train)}  # fmt: skip
        for c in cands:
            for v in evalio.VARIANTS:
                s, call = hash_score(scores, s1, h, v, c)
                row[f"call_{c}_{v}"] = "1" if call else "0"
                row[f"score_{c}_{v}"] = "" if math.isnan(s) else format(float(s), ".10g")
        rows.append(row)
    return rows


def proteome_columns(cands) -> tuple:
    return PROTEOME_BASE + tuple(
        f"{kind}_{c}_{v}" for c in cands for v in evalio.VARIANTS for kind in ("call", "score")
    )


def agreement_block(prows, table, scores, s1, cands) -> dict:
    out = {"proteomes": {}, "truth": {}}
    ml = [c for c in cands if c in findings.ML]
    for set_id in sorted({r["set_id"] for r in prows}):
        rows = [r for r in prows if r["set_id"] == set_id]
        out["proteomes"][set_id] = {
            v: {c: agreement.agreement_counts([r[f"call_R2_{v}"] == "1" for r in rows],
                                              [r[f"call_{c}_{v}"] == "1" for r in rows])
                for c in ml}
            for v in evalio.VARIANTS
        }  # fmt: skip
    for cls in ("pos", "neg"):
        hashes = [h for h, t in table.items() if t["origin"] == "go" and t["class"] == cls]
        out["truth"][cls] = {}
        for v in evalio.VARIANTS:
            rule = [hash_score(scores, s1, h, v, "R2")[1] for h in hashes]
            out["truth"][cls][v] = {
                c: agreement.agreement_counts(rule, [hash_score(scores, s1, h, v, c)[1] for h in hashes])
                for c in ml
            }  # fmt: skip
    return out


def named_panel(seq_members, table, literature, scores, s1, full_train, cands) -> list:
    panel = list(agreement.PANEL)
    have = {p[1] for p in panel}
    for r in literature:
        if r["lit_class"] == "hard_negative" and r["accession"] not in have:
            panel.append((r["gene"], r["accession"], "literature hard_negative row (ruling C-9)"))
    found = agreement.find_panel_hashes(panel, seq_members)
    lit = {r["accession"]: r for r in literature}
    out = []
    for name, pid, why in panel:
        hit = found.get(pid)
        entry = {"name": name, "id": pid, "why": why, "found": hit is not None}
        if hit:
            h = hit["seq_sha256"]
            t = table.get(h, {})
            entry.update({"seq_sha256": h, "set_ids": sorted(set(hit["set_ids"])),
                          "class": t.get("class", ""), "label": t.get("label", ""),
                          "stratum": t.get("stratum", ""),
                          "literature_class": lit.get(pid, {}).get("lit_class", ""),
                          "score_source": agreement.score_source(h, s1, full_train),
                          "calls": {}, "scores": {}})  # fmt: skip
            for c in cands:
                for v in evalio.VARIANTS:
                    s, call = hash_score(scores, s1, h, v, c)
                    entry["calls"][f"{c}|{v}"] = int(call)
                    entry["scores"][f"{c}|{v}"] = num(s)
        out.append(entry)
    return out


def build_findings(test_blocks: dict) -> dict:
    s2 = sorted(n for n, b in test_blocks.items() if b["split"].startswith("S2-"))
    b1 = {}
    for name in ["S1:all", *s2]:
        if name in test_blocks:
            m = test_blocks[name]["truth"]["direct"]["metrics"]["all"]["V-go"].get("B1", {})
            b1[name] = m.get("roc_auc", {}).get("value")

    def diffs(name):
        block = test_blocks[name]["truth"]["direct"]["nsec_fpr_at_rule_recall"].get("V-go", {})
        return block.get("diff", {})

    return {
        "schema": FINDINGS_SCHEMA,
        "a_b1_not_saturated": findings.finding_a(b1),
        "b_ml_beats_b1_and_r2_on_nsec_s1": {
            "test_set": "S1:all",
            **findings.finding_b(diffs("S1:all") if "S1:all" in test_blocks else {}),
        },
        "c_same_under_s2": findings.finding_c({n: findings.finding_b(diffs(n)) for n in s2}),
    }


def run(work: Path, sets_path: Path, species_path: Path, n_resamples: int, arguments=()):
    work = Path(work)
    out = evalio.out_dir(work)
    build = evalio.require_current(
        out, "build_run.json", ("eval_table.tsv.gz", "eval_literature.tsv"), "08"
    )
    split_log = evalio.require_current(
        out,
        "splits_run.json",
        ("split_members.tsv.gz", "clusters.tsv.gz", "max_identity.tsv.gz"),
        "09",
    )
    score_log = evalio.require_current(
        out, "scores_run.json", ("scores.tsv.gz",), "10_fit_and_score.py"
    )
    for name in ("eval_table.tsv.gz", "eval_literature.tsv", "eval_sequences.fasta.gz"):
        if split_log["input_sha256"].get(name) != build["outputs_sha256"].get(name):
            raise evalio.StopError(
                f"splits_run.json was made from another {name} than build_run.json records; "
                "re-run 09_make_splits.py"
            )
    for name in ("split_members.tsv.gz", "clusters.tsv.gz"):
        if score_log["input_sha256"].get(name) != split_log["outputs_sha256"].get(name):
            raise evalio.StopError(
                f"scores.tsv.gz was made from another {name} than splits_run.json records; "
                "re-run 10"
            )
    fr = evalio.read_json(work / "phaseb" / "features_run.json")
    seq_members_path = work / "phaseb" / "sequence_members.tsv.gz"
    if manifest.sha256_file(seq_members_path) != fr["input_sha256"]["sequence_members.tsv.gz"]:
        raise evalio.StopError("phaseb/sequence_members.tsv.gz differs from the file 07 read")
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    literature = truth_table.read_tsv(out / "eval_literature.tsv")
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    splits.check_test_truth(members)
    identity = truth_table.read_tsv(out / "max_identity.tsv.gz")
    if manifest.sha256_file(species_path) != build["input_sha256"]["species.tsv"]:
        raise evalio.StopError(f"{species_path} differs from the species.tsv that 08 used")
    species = truth_table.read_tsv(species_path)
    scores = read_scores(out / "scores.tsv.gz")
    cands = tuple(score_log["candidates"])
    blocks = {}
    for ts in test_sets(members, table, literature, identity, species):
        res = {"split": ts.split, "kind": ts.kind, "truth": {}}
        for truth in TRUTHS:
            block, rows, y, data = evaluate_truth(
                ts, truth, scores, cands, n_resamples, evalio.SEED
            )
            res["truth"][truth] = block
            if truth == "direct" and ts.split.startswith("S2-"):
                res["calibration"] = calibration(rows, y, data, cands)
        direct = res["truth"]["direct"]["metrics"]["all"]["V-go"]
        hw = {
            c: bootstrap.half_width(direct[c]["recall"])
            for c in findings.LABEL_CANDIDATES
            if c in direct
        }
        res["recall_half_width"] = hw
        res["fpr_half_width"] = {c: bootstrap.half_width(direct[c]["fpr"]) for c in hw}
        n_pos = res["truth"]["direct"]["n"]["all"]["pos"]
        res["n_direct_positives"] = n_pos
        res["floor_met"] = findings.floor_met(n_pos)
        res["label"] = findings.estimate_label(hw, n_pos)
        for key, src in (
            ("max_recall_half_width", hw),
            ("max_fpr_half_width", res["fpr_half_width"]),
        ):
            vals = list(src.values())
            res[key] = None if not vals or None in vals else max(vals)
        res["zero_width_recall_interval"] = sorted(
            c for c in hw if direct[c]["recall"]["lo"] is not None
            and direct[c]["recall"]["lo"] == direct[c]["recall"]["hi"]
        )  # fmt: skip
        res.update(excluded_views(ts, scores, cands))
        blocks[ts.name] = res
    s1, full_train = lookup_tables(members)
    set_rows = truth_table.read_tsv(sets_path)
    proteome_sets = {r["set_id"] for r in set_rows if r["kind"] in PROTEOME_KINDS}
    seq_members = truth_table.read_tsv(seq_members_path)
    prows = proteome_rows(seq_members, proteome_sets, table, scores, s1, full_train, cands)
    onygenales = Counter(
        f"{m['split_id']}|{m['fold']}" for m in members
        if m["part"] == "train_tc"
        and set(table[m["seq_sha256"]]["taxon_ids"].split(",")) & set(ONYGENALES_TC_TAXA)
    )  # fmt: skip
    metrics_json = {
        "schema": METRICS_SCHEMA,
        "settings": {
            "n_resamples": n_resamples,
            "seed": evalio.SEED,
            "ci_level": 95,
            "estimate_half_width": findings.ESTIMATE_HALF_WIDTH,
            "estimate_min_direct_positives": findings.MIN_DIRECT_POSITIVES,
            "saturation_auc": findings.SATURATION_AUC,
            "prevalences": list(PREVALENCES),
            "recall_levels": list(RECALL_LEVELS),
            "fpr_level": FPR_LEVEL,
            "long_cutoff": LONG_CUTOFF,
            "identity_cutoff": splits.IDENTITY_CUT,
            "calibration_bins": CALIBRATION_BINS,
            "candidates": list(cands),
            "variants": list(evalio.VARIANTS),
        },  # fmt: skip
        "test_sets": blocks,
        "agreement": agreement_block(prows, table, scores, s1, cands),
        "score_sources": {
            s: dict(sorted(Counter(r["score_source"] for r in prows if r["set_id"] == s).items()))
            for s in sorted(proteome_sets)
        },
        "named_panel": named_panel(seq_members, table, literature, scores, s1, full_train, cands),
        "context": {"onygenales_tc_rows_in_vkw_training": dict(sorted(onygenales.items()))},
        "tc_removed": split_log["tc_removed_by_split_fold_rule"],
    }
    findings_json = build_findings(blocks)
    log = {
        "all_sources": build["all_sources"],
        "truth_set_sha256": build["truth_set_sha256"],
        "input_sha256": {
            "scores.tsv.gz": score_log["outputs_sha256"]["scores.tsv.gz"],
            "split_members.tsv.gz": split_log["outputs_sha256"]["split_members.tsv.gz"],
            "eval_table.tsv.gz": build["outputs_sha256"]["eval_table.tsv.gz"],
            "sequence_members.tsv.gz": fr["input_sha256"]["sequence_members.tsv.gz"],
            "sequence_sets.tsv": manifest.sha256_file(sets_path),
        },
        "seed": evalio.SEED,
        "n_resamples": n_resamples,
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    columns = proteome_columns(cands)

    evalio.write_outputs(
        out,
        {
            "metrics.json": lambda p: runinfo.write_json(p, metrics_json),
            "findings.json": lambda p: runinfo.write_json(p, findings_json),
            "proteome_calls.tsv.gz": lambda p: truth_table.write_tsv(p, columns, prows),
        },
        "evaluate_run.json",
        log,
    )
    return metrics_json, findings_json


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sets", default=str(paths.STEP1_DIR / "sequence_sets.tsv"))
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--n-resamples", type=int, default=bootstrap.N_RESAMPLES)
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        m, f = run(work, Path(args.sets), Path(args.species), args.n_resamples,
                   list(argv) if argv is not None else sys.argv[1:])  # fmt: skip
    except (evalio.StopError, splits.SplitError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for name, block in m["test_sets"].items():
        print(f"{name}\t{block['label']}")
    print(f"finding_a_holds={f['a_b1_not_saturated']['holds']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
