#!/usr/bin/env python3
"""PR 6 hyd1_like pre-freeze check: which of the 161 shared hard negatives are hyd1_like family members or other hydrophobin-class proteins?

DRAFT, nothing run on HPCC. The ecm33 pre-check found Q4WNS8 this way, after the freeze draft had named only P38248. Same method, imported (not copied)
from pr6_ecm33_negatives_precheck.py, with two additions: (1) the hydrophobin-class Pfam panel at GA over the same negatives, because a hydrophobin that
is not a homolog of the Kiss protein by phmmer is still a hydrophobin, and (2) an optional relatedness table of the Kiss protein against the committed
hydrophobin truth proteins (descriptive; it says how many curated hydrophobins the Kiss protein is related to at all).

    python3 -I pr6_hyd1_negatives_precheck.py summarize --strict phmmer_E1e-5.tbl --weak phmmer_E1e-3.tbl --panel panel_ga.tbl \\
        --controls data/controls/hard-negatives-shared/controls.tsv --negatives-fasta data/controls/hard-negatives-shared/sequences.faa \\
        [--related related_E1e-3.tbl] --out OUT_DIR

Files written (ids, genes, scores and lengths only; no sequence):
  negatives_family_table.tsv   negatives found by phmmer of A0ACE4KFX7 at E <= 1e-5 (the freeze's labelled_rule)
  negatives_panel_table.tsv    negatives hit by a hydrophobin-class Pfam model at GA (model, score), whether or not phmmer found them
  weak_E1e-3.tsv               DESCRIPTIVE ONLY
  related.tsv                  (with --related) truth proteins within E 1e-3 of the Kiss protein, descriptive
  proposed_negatives_ignored.json  the union of the phmmer hits (the family rule); the panel hits are NOT added: that is an owner question
  summary.txt                  counts and LOUD lines. The freeze file is not read or changed.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


e33 = _load("pr6_ecm33_negatives_precheck")
QUERY = "A0ACE4KFX7"
PANEL_HEADER = ["accession", "gene", "species", "stratum", "model", "score", "also_phmmer_family"]
RELATED_HEADER = ["target", "evalue", "score"]


def panel_rows(panel: list[dict], controls: dict[str, dict], family: set[str]) -> list[dict]:
    """One row per (negative, model) GA hit, best score first. In the real file `query` is the model name."""
    rows = []
    for h in sorted(panel, key=lambda x: (-x["score"], x["target"])):
        if h["target"] not in controls:
            continue
        c = controls[h["target"]]
        rows.append(
            {
                "accession": h["target"],
                "gene": c.get("gene", ""),
                "species": c.get("species", ""),
                "stratum": c.get("stratum", ""),
                "model": h["query"],
                "score": h["score"],
                "also_phmmer_family": "yes" if h["target"] in family else "no",
            }
        )
    return rows


def related_rows(related: list[dict], max_e: float = e33.WEAK_E) -> list[dict]:
    return [
        {"target": h["target"], "evalue": h["evalue"], "score": h["score"]}
        for h in sorted(related, key=lambda x: (x["evalue"], x["target"]))
        if h["evalue"] <= max_e
    ]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summarize")
    for flag in ("strict", "weak", "panel", "controls", "negatives-fasta"):
        s.add_argument(f"--{flag}", required=True, type=Path)
    s.add_argument("--related", type=Path)
    s.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    if a.out.exists() and any(a.out.iterdir()):
        print(f"FATAL: {a.out} is not empty; nothing is overwritten", file=sys.stderr)
        return 1
    controls = e33.read_controls(a.controls)
    res = e33.build(
        e33.parse_tblout(a.strict),
        e33.parse_tblout(a.weak),
        controls,
        e33.fasta_lengths(a.negatives_fasta),
        [],
        expected=(),
    )
    family = {r["accession"] for r in res["rows"]}
    prow = panel_rows(e33.parse_tblout(a.panel), controls, family)
    a.out.mkdir(parents=True, exist_ok=True)
    e33.write_tsv(a.out / "negatives_family_table.tsv", e33.TABLE_HEADER, res["rows"])
    e33.write_tsv(a.out / "negatives_panel_table.tsv", PANEL_HEADER, prow)
    e33.write_tsv(a.out / "weak_E1e-3.tsv", e33.WEAK_HEADER, res["weak_rows"])
    n_rel = ""
    if a.related:
        rel = related_rows(e33.parse_tblout(a.related))
        e33.write_tsv(a.out / "related.tsv", RELATED_HEADER, rel)
        n_rel = f"truth proteins within E 1e-3 of {QUERY} (descriptive): {len(rel)}"
    edit = e33.json_edit(res["proposed"])
    (a.out / "proposed_negatives_ignored.json").write_text(edit + "\n")
    only_panel = sorted({r["accession"] for r in prow} - family)
    lines = [
        "PR 6 hyd1_like pre-freeze check (DRAFT): shared hard negatives that are family members or hydrophobin-class",
        f"family definition: phmmer of {QUERY} finds the negative at E <= {e33.STRICT_E:g} (the freeze's labelled_rule)",
        f"family members among the negatives: {len(res['rows'])} ({', '.join(sorted(family)) or 'none'})",
        f"hydrophobin-class Pfam GA hits among the negatives: {len({r['accession'] for r in prow})} ({', '.join(sorted({r['accession'] for r in prow})) or 'none'})",
        f"panel hits that phmmer does NOT call family (owner question 8; not in the proposed list): {', '.join(only_panel) or 'none'}",
        f"weak hits (1e-5 < E <= 1e-3), descriptive: {len(res['weak_rows'])}",
        n_rel,
        "proposed edit of the freeze (the owner decides; this script changes nothing):",
        "  " + edit,
    ] + res["loud"]
    (a.out / "summary.txt").write_text("\n".join(x for x in lines if x) + "\n")
    print((a.out / "summary.txt").read_text())
    return 0


if __name__ == "__main__":
    sys.exit(main())
