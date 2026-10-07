"""Compare two independent fills of the repeat-mechanism curation table.

Input: two copies of `data/controls/repeat-mechanism/curation_table.tsv`, one filled by each agent.
Output (in --out): `comparison.tsv`, `review_queue.tsv`, `summary.md`.

What it does:
1. Joins the two tables row by row (accession, plus the order among equal accessions).
2. Reports raw agreement and Cohen's kappa for the mechanism label, and for the binary question
   that matters for the repeat call: class 2a or not.
3. Checks, for each agent, whether the quote contains words that fit the label it chose
   (`quote_check`: pass, fail, na). This is a keyword check. It cannot prove support. A fail means
   "a person must read this row". A pass does not mean the label is right.
4. Builds a review queue: every row that needs a person.
5. Counts independent clusters of 2a under three rules, from loosest to strictest.

It changes no input file. Extra columns that an agent added (for example `evidence_basis`) are
carried into `comparison.tsv` with the suffix `_A` or `_B`.

Usage:
    python3.12 compare_mechanism_labels.py --agent-a A.tsv --agent-b B.tsv --out DIR
"""

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

REQUIRED = [
    "accession",
    "cls",
    "evidence_level",
    "cluster_id",
    "mechanism_label",
    "label_confidence",
    "mechanism_source_pmid",
    "mechanism_quote",
]
LABELS = ["2a", "2b-i", "2b-ii", "2b-iii", "2c", "2d", "other", "unknown"]
CONFIDENCE = ["high", "medium", "low"]

# Words that fit each label. A keyword check, not a proof.
QUOTE_WORDS = {
    "2a": r"repeat|tandem|repetitive",
    "2b-i": r"CFEM|cysteine-rich extracellular membrane",
    "2b-ii": r"cys(teine)?[- ]knot|cystine[- ]knot",
    "2b-iii": r"Bys1|thaumatin",
    "2c": r"hydrophobin|rodlet|8[- ]Cys|eight cysteine",
    "2d": r"moonlight|chaperon|enolase|heat shock|Hsp\d+|mitochondri|cytoplasm|cytosol",
}


def read_table(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        header = reader.fieldnames or []
    missing = [c for c in REQUIRED if c not in header]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    return header, rows


def keyed(rows, path):
    """Key = (accession, n-th row with that accession). Refuses a row without an accession."""
    seen = Counter()
    out = {}
    for n, row in enumerate(rows, start=2):
        acc = row["accession"].strip()
        if not acc:
            raise ValueError(f"{path}:{n}: row without an accession")
        seen[acc] += 1
        out[(acc, seen[acc])] = row
    return out


def quote_check(label, quote):
    """`pass`, `fail` or `na` (no check exists for the label, or no label)."""
    label = label.strip()
    if label not in QUOTE_WORDS:
        return "na"
    return "pass" if re.search(QUOTE_WORDS[label], quote or "", re.IGNORECASE) else "fail"


def cohen_kappa(pairs):
    """Cohen's kappa for a list of (label_a, label_b). None if undefined."""
    n = len(pairs)
    if n == 0:
        return None
    observed = sum(1 for a, b in pairs if a == b) / n
    ca = Counter(a for a, _ in pairs)
    cb = Counter(b for _, b in pairs)
    expected = sum((ca[k] / n) * (cb[k] / n) for k in set(ca) | set(cb))
    if expected == 1.0:
        return None
    return (observed - expected) / (1 - expected)


def fmt(x):
    return "n/a" if x is None else f"{x:.3f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--agent-a", required=True)
    ap.add_argument("--agent-b", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    head_a, rows_a = read_table(args.agent_a)
    head_b, rows_b = read_table(args.agent_b)
    ka, kb = keyed(rows_a, args.agent_a), keyed(rows_b, args.agent_b)
    if set(ka) != set(kb):
        only_a, only_b = sorted(set(ka) - set(kb)), sorted(set(kb) - set(ka))
        raise ValueError(
            f"the two tables hold different rows: only in A {only_a[:5]}, only in B {only_b[:5]}"
        )
    extra = [c for c in dict.fromkeys(head_a + head_b) if c not in REQUIRED]
    extra = [c for c in extra if c not in ("reviewer", "review_date")]
    base = [c for c in head_a if c in ("gene", "organism", "family", "pmids") or c in REQUIRED[:4]]

    # Adhesin rows with evidence E1 or E2 are the ones that need a label.
    keys = sorted(ka, key=lambda k: (k[0], k[1]))
    target = [
        k for k in keys if ka[k]["cls"] == "adhesin" and ka[k]["evidence_level"] in ("E1", "E2")
    ]
    for k in target:
        for tag, tab in (("A", ka), ("B", kb)):
            lab = tab[k]["mechanism_label"].strip()
            if lab not in LABELS:
                raise ValueError(f"agent {tag}, {k[0]}: label '{lab}' is not one of {LABELS}")
            conf = tab[k]["label_confidence"].strip()
            if conf not in CONFIDENCE:
                raise ValueError(
                    f"agent {tag}, {k[0]}: confidence '{conf}' is not one of {CONFIDENCE}"
                )

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    comp_cols = (
        base
        + [
            f"{c}_{t}"
            for t in "AB"
            for c in (
                "mechanism_label",
                "label_confidence",
                "mechanism_source_pmid",
                "mechanism_quote",
            )
        ]
        + [f"quote_check_{t}" for t in "AB"]
        + ["label_agree", "is_2a_agree", "queue_reasons"]
        + [f"{c}_{t}" for c in extra for t in "AB"]
    )
    records, queue = [], []
    for k in target:
        ra, rb = ka[k], kb[k]
        la, lb = ra["mechanism_label"].strip(), rb["mechanism_label"].strip()
        ca, cb = ra["label_confidence"].strip(), rb["label_confidence"].strip()
        qa, qb = quote_check(la, ra["mechanism_quote"]), quote_check(lb, rb["mechanism_quote"])
        reasons = []
        if la != lb:
            reasons.append("labels_differ")
        if (la == "2a") != (lb == "2a"):
            reasons.append("2a_vs_not_2a")
        if "low" in (ca, cb):
            reasons.append("low_confidence")
        if qa == "fail":
            reasons.append("quote_check_fail_A")
        if qb == "fail":
            reasons.append("quote_check_fail_B")
        if la == "unknown" or lb == "unknown":
            reasons.append("unknown_label")
        for tag, tab in (("A", ra), ("B", rb)):
            if not tab["mechanism_source_pmid"].strip() or not tab["mechanism_quote"].strip():
                reasons.append(f"no_source_or_quote_{tag}")
        rec = {c: ra.get(c, "") for c in base}
        for tag, tab in (("A", ra), ("B", rb)):
            for c in (
                "mechanism_label",
                "label_confidence",
                "mechanism_source_pmid",
                "mechanism_quote",
            ):
                rec[f"{c}_{tag}"] = tab[c]
        rec.update(
            quote_check_A=qa,
            quote_check_B=qb,
            label_agree="yes" if la == lb else "no",
            is_2a_agree="yes" if (la == "2a") == (lb == "2a") else "no",
            queue_reasons=";".join(reasons),
        )
        for c in extra:
            rec[f"{c}_A"], rec[f"{c}_B"] = ra.get(c, ""), rb.get(c, "")
        records.append((rec, la, lb, qa, qb, ra["cluster_id"]))
        if reasons:
            queue.append(rec)

    def write(path, cols, rows):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    write(out / "comparison.tsv", comp_cols, [r[0] for r in records])
    write(out / "review_queue.tsv", comp_cols, queue)

    # Summary numbers.
    n = len(records)
    pairs = [(la, lb) for _, la, lb, _, _, _ in records]
    pairs_2a = [(la == "2a", lb == "2a") for _, la, lb, _, _, _ in records]
    agree = sum(1 for a, b in pairs if a == b)
    agree_2a = sum(1 for a, b in pairs_2a if a == b)
    conf = Counter((la, lb) for _, la, lb, _, _, _ in records)
    reasons = Counter(r for rec, *_ in records for r in rec["queue_reasons"].split(";") if r)

    def clusters(rule):
        c = defaultdict(list)
        for rec, la, lb, qa, qb, cl in records:
            c[cl].append((rec, la, lb, qa, qb))
        return sum(1 for rows in c.values() if rows and any(rule(*r) for r in rows))

    def both2a(rec, la, lb, qa, qb):
        return la == "2a" and lb == "2a"

    def both2a_q(rec, la, lb, qa, qb):
        return both2a(rec, la, lb, qa, qb) and qa == "pass" and qb == "pass"

    def any2a(rec, la, lb, qa, qb):
        return la == "2a" or lb == "2a"

    lines = [
        "# Comparison of two fills of the repeat-mechanism table",
        "",
        f"Agent A: `{args.agent_a}`  \nAgent B: `{args.agent_b}`  \nRows compared (adhesin, E1 or E2): {n}",
        "",
        "## Agreement",
        "",
        f"- Label (all classes): {agree} of {n} agree ({agree / n:.1%}); Cohen's kappa {fmt(cohen_kappa(pairs))}.",
        f"- Class 2a or not: {agree_2a} of {n} agree ({agree_2a / n:.1%}); Cohen's kappa {fmt(cohen_kappa(pairs_2a))}.",
        "- Kappa is unstable when one class holds most rows. Read it with the counts below.",
        "",
        "## Distribution",
        "",
        "| Label | Agent A | Agent B |",
        "|---|---|---|",
    ]
    ca_, cb_ = Counter(p[0] for p in pairs), Counter(p[1] for p in pairs)
    for lab in LABELS:
        lines.append(f"| {lab} | {ca_[lab]} | {cb_[lab]} |")
    lines += [
        "",
        "Confidence (A / B): "
        + "; ".join(
            f"{c} {sum(1 for r in records if r[0]['label_confidence_A'] == c)} / "
            f"{sum(1 for r in records if r[0]['label_confidence_B'] == c)}"
            for c in CONFIDENCE
        ),
    ]
    lines += [
        "",
        "## Confusion matrix (rows A, columns B)",
        "",
        "| A \\ B | " + " | ".join(LABELS) + " |",
        "|---|" + "---|" * len(LABELS),
    ]
    for la in LABELS:
        lines.append(f"| {la} | " + " | ".join(str(conf[(la, lb)]) for lb in LABELS) + " |")
    lines += [
        "",
        "## Quote check (keyword check; a fail means a person reads the row)",
        "",
        "| | pass | fail | na |",
        "|---|---|---|---|",
    ]
    for tag, idx in (("A", 3), ("B", 4)):
        c = Counter(r[idx] for r in records)
        lines.append(f"| Agent {tag} | {c['pass']} | {c['fail']} | {c['na']} |")
    lines += [
        "",
        "Rows labelled 2a with a quote that has no repeat, tandem or repetitive wording:",
        f"A {sum(1 for r in records if r[1] == '2a' and r[3] == 'fail')}, "
        f"B {sum(1 for r in records if r[2] == '2a' and r[4] == 'fail')}.",
        "Rows labelled 2a with confidence `high` and a failed quote check (the confidence does not "
        "match the quote): "
        f"A {sum(1 for r in records if r[1] == '2a' and r[3] == 'fail' and r[0]['label_confidence_A'] == 'high')}, "
        f"B {sum(1 for r in records if r[2] == '2a' and r[4] == 'fail' and r[0]['label_confidence_B'] == 'high')}.",
        "",
        "## Independent clusters of class 2a (30% identity, 50% coverage)",
        "",
        "| Rule | Clusters | Needed for `estimated` |",
        "|---|---|---|",
        f"| Either agent labels at least one row of the cluster 2a (loosest) | {clusters(any2a)} | 20 |",
        f"| Both agents label at least one row of the cluster 2a | {clusters(both2a)} | 20 |",
        f"| Both agree on 2a and both quotes pass the keyword check (strictest) | {clusters(both2a_q)} | 20 |",
        "",
        "The strictest row is the number to quote in the paper until a person has read the review queue.",
        "A cluster count says nothing about sensitivity. It only shows whether the 20-cluster floor "
        "can be met.",
        "",
        "## Review queue",
        "",
        f"{len(queue)} of {n} rows need a person. File: `review_queue.tsv`. Reasons:",
        "",
    ]
    lines += [f"- {k}: {v}" for k, v in reasons.most_common()]
    if extra:
        lines += ["", "## Extra columns carried through", "", ", ".join(f"`{c}`" for c in extra)]
        for c in extra:
            if c == "evidence_basis":
                lines += ["", "`evidence_basis` (A / B):"]
                for v in sorted(
                    {r[0].get(f"{c}_A", "") for r in records}
                    | {r[0].get(f"{c}_B", "") for r in records}
                ):
                    lines.append(
                        f"- `{v}`: {sum(1 for r in records if r[0][f'{c}_A'] == v)} / "
                        f"{sum(1 for r in records if r[0][f'{c}_B'] == v)}"
                    )
    (out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{n} rows; label agreement {agree / n:.1%}; queue {len(queue)}; written to {out}")


if __name__ == "__main__":
    main()
