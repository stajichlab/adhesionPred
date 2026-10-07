"""Adjudicate two independent fills of the repeat-mechanism table, using evidence tiers.

Inputs: agent-A table, agent-B table (same rows, filled independently), and the UniProt repeat-feature
table (`uniprot_repeat_features.py`). It changes none of them.

Rule. Each source gives a vote for a mechanism label at an evidence tier:

| Tier (rank) | Source |
|---|---|
| paper_stated (3) | an agent marks `stated_in_paper` |
| database_annotation (2) | UniProt has 2 or more `Repeat` features: votes `2a` |
| family_inference (1) | an agent marks `family_inference` |
| background_knowledge (0) | an agent marks `background_knowledge` |

A vote of `unknown` is no vote. The label at the highest tier wins, whatever the lower tiers say.
Two different labels at the same highest tier cannot be settled by the rule: the row goes to the
expert review packet and gets the final label `unresolved`.

Two questions are kept apart, because a repeat detector is tested on the first one:

1. **Has the protein a tandem repeat region?** Only `2a` votes count. A label of `other`, `2c` or
   `2d` from an agent is silence about repeats (the paper it read gives another mechanism). It is not
   evidence against repeats. The tier of this answer is the highest tier of any `2a` vote
   (`repeat_evidence_tier`). Cluster counts for the repeat call use this question.
2. **What mechanism does the evidence give?** The rule above (`final_label`).

This breaks ties by strength of evidence. It does not decide which agent was right, and it does not
read any paper. Rows that need a person are written to `expert_review_packet.tsv`.

Outputs in --out: `curation_table.consensus.tsv`, `expert_review_packet.tsv`, `adjudication_summary.md`.

Usage: python3.12 adjudicate_mechanism_labels.py --agent-a A.tsv --agent-b B.tsv \
         --uniprot uniprot_repeat_features.tsv --out DIR
"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

RANK = {
    "stated_in_paper": 3,
    "database_annotation": 2,
    "family_inference": 1,
    "background_knowledge": 0,
}
LABELS = ["2a", "2b-i", "2b-ii", "2b-iii", "2c", "2d", "other", "unknown"]
MIN_REPEATS = 2


def read(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def keyed(rows):
    seen, out = Counter(), {}
    for row in rows:
        acc = row["accession"].strip()
        seen[acc] += 1
        out[(acc, seen[acc])] = row
    return out


def votes(row_a, row_b, uniprot):
    """List of (source, label, rank). `unknown` gives no vote."""
    out = []
    for tag, row in (("A", row_a), ("B", row_b)):
        label = row["mechanism_label"].strip()
        basis = row.get("evidence_basis", "").strip()
        if label not in LABELS:
            raise ValueError(f"agent {tag}, {row['accession']}: label '{label}' is not in {LABELS}")
        if basis not in RANK:
            raise ValueError(
                f"agent {tag}, {row['accession']}: evidence_basis '{basis}' not in {sorted(RANK)}"
            )
        if label != "unknown":
            out.append((tag, label, RANK[basis]))
    u = uniprot.get(row_a["accession"].strip())
    if u and int(u["n_repeat_features"]) >= MIN_REPEATS:
        out.append(("UniProt", "2a", RANK["database_annotation"]))
    return out


def decide(vote_list):
    """Return (final_label, tier, rule, needs_expert)."""
    if not vote_list:
        return "unknown", "none", "no_claim", "no"
    top = max(rank for _, _, rank in vote_list)
    tier = [k for k, v in RANK.items() if v == top][0]
    top_labels = {label for _, label, rank in vote_list if rank == top}
    all_labels = {label for _, label, _ in vote_list}
    if len(top_labels) > 1:
        return "unresolved", tier, "conflict_at_same_tier", "yes"
    final = next(iter(top_labels))
    if len(all_labels) == 1:
        return final, tier, "agreement", "no"
    return final, tier, "higher_tier_wins", "no"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--agent-a", required=True)
    ap.add_argument("--agent-b", required=True)
    ap.add_argument("--uniprot", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows_a, rows_b = read(args.agent_a), read(args.agent_b)
    ka, kb = keyed(rows_a), keyed(rows_b)
    if set(ka) != set(kb):
        raise ValueError("the two tables hold different rows")
    uniprot = {r["accession"]: r for r in read(args.uniprot)}
    target = [
        k for k in ka if ka[k]["cls"] == "adhesin" and ka[k]["evidence_level"] in ("E1", "E2")
    ]

    out_rows, packet = [], []
    for k in target:
        ra, rb = ka[k], kb[k]
        vl = votes(ra, rb, uniprot)
        final, tier, rule, expert = decide(vl)
        rep_votes = [(s, r) for s, lab, r in vl if lab == "2a"]
        rep_rank = max((r for _, r in rep_votes), default=None)
        rep_tier = "none" if rep_rank is None else [k for k, v in RANK.items() if v == rep_rank][0]
        rep_sources = (
            ",".join(sorted({s for s, r in rep_votes if r == rep_rank})) if rep_votes else ""
        )
        u = uniprot.get(k[0], {})
        rec = {
            "accession": k[0],
            "gene": ra.get("gene", ""),
            "organism": ra.get("organism", ""),
            "cluster_id": ra["cluster_id"],
            "evidence_level": ra["evidence_level"],
            "label_A": ra["mechanism_label"],
            "basis_A": ra.get("evidence_basis", ""),
            "label_B": rb["mechanism_label"],
            "basis_B": rb.get("evidence_basis", ""),
            "uniprot_repeat_features": u.get("n_repeat_features", ""),
            "final_label": final,
            "final_tier": tier,
            "rule": rule,
            "needs_expert": expert,
            "repeat_evidence_tier": rep_tier,
            "repeat_evidence_sources": rep_sources,
        }
        out_rows.append(rec)
        if expert == "yes":
            packet.append(
                {
                    **rec,
                    "protein_name": ra.get("protein_name", ""),
                    "pmids_in_table": ra.get("pmids", ""),
                    "source_A": ra.get("mechanism_source_pmid", ""),
                    "quote_A": ra.get("mechanism_quote", ""),
                    "source_B": rb.get("mechanism_source_pmid", ""),
                    "quote_B": rb.get("mechanism_quote", ""),
                    "question": "Which class of the repository scheme (2a, 2b-i, 2b-ii, 2b-iii, 2c, 2d, other) does "
                    "the cited evidence support for THIS protein? Quote the sentence.",
                    "expert_label": "",
                    "expert_source_pmid": "",
                    "expert_quote": "",
                    "expert_name": "",
                }
            )

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    def write(name, rows):
        if not rows:
            (out / name).write_text("")
            return
        with open(out / name, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    write("curation_table.consensus.tsv", out_rows)
    write("expert_review_packet.tsv", packet)

    n = len(out_rows)
    rules = Counter(r["rule"] for r in out_rows)
    finals = Counter(r["final_label"] for r in out_rows)
    by_cluster = defaultdict(list)
    for r in out_rows:
        by_cluster[r["cluster_id"]].append(r)

    def clusters(tiers):
        return sum(
            1
            for rows in by_cluster.values()
            if any(r["repeat_evidence_tier"] in tiers for r in rows)
        )

    t1 = clusters({"stated_in_paper"})
    t2 = clusters({"stated_in_paper", "database_annotation"})
    t3 = clusters({"stated_in_paper", "database_annotation", "family_inference"})
    lines = [
        "# Adjudication of the repeat-mechanism labels",
        "",
        f"Rows (adhesin, E1 or E2): {n}. Rule: highest evidence tier wins; a same-tier conflict goes to an expert.",
        "",
        "## How rows were settled",
        "",
        "| Rule | Rows |",
        "|---|---|",
    ]
    lines += [f"| {k} | {v} |" for k, v in rules.most_common()]
    lines += ["", "## Final mechanism label", "", "| Label | Rows |", "|---|---|"]
    lines += [f"| {k} | {finals[k]} |" for k in LABELS + ["unresolved"] if finals[k]]
    lines += [
        "",
        "## Independent clusters with a repeat region (30% identity, 50% coverage), cumulative by evidence tier",
        "",
        "| Evidence tier included | Clusters | Needed |",
        "|---|---|---|",
        f"| paper_stated | {t1} | 20 |",
        f"| + database_annotation (UniProt repeat features) | {t2} | 20 |",
        f"| + family_inference | {t3} | 20 |",
        "",
        "Votes of `other`, `2c` or `2d` are not counted against repeats. They mean that the cited paper gives "
        "another mechanism or none.",
        "",
        "A cluster count does not make a status. `estimated` also needs a 95% half-width of at most 0.10, which "
        "needs about 61 clusters at sensitivity 0.8 (binomial approximation).",
        "",
        f"## Expert review packet\n\n{len(packet)} rows. File: `expert_review_packet.tsv`. Each row has both agents' "
        "labels, sources and quotes, and blank `expert_*` columns.",
    ]
    (out / "adjudication_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(
        f"{n} rows; rules {dict(rules)}; 2a clusters by tier {t1}/{t2}/{t3}; packet {len(packet)}"
    )


if __name__ == "__main__":
    main()
