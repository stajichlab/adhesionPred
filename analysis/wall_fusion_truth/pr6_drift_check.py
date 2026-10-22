#!/usr/bin/env python3
"""PR 6 drift check for the queries whose outside homologs come only from jackhmmer iteration (DRAFT, not run).

aga1 (KS0080, P32323), ece1 (KS0107, Q07730) and, after the clade-bucket fix, hyd1_like (KS0017, A0ACE4KFX7) are
`buildable_iterated_only`: the one-pass phmmer set does not meet the training rule, the jackhmmer set does. Iteration
can pull in unrelated proteins that share only a composition (Ser/Thr-rich, repeat-rich or short-motif proteins). The
raw output holds the final iteration only, so this check does not read iterations. It asks of every jackhmmer outside
fungal hit: does a direct pairwise alignment to the query support it, or is it explained by composition?

Four commands (standard library only; run with `python3 -I`; no sequence is ever written to a committed table):

  extract-ids  --record DIR --out DIR     the jackhmmer outside fungal hits of the queries -> ids.txt, hits_manifest.tsv
  split        --manifest TSV --hits-faa FA --queries-faa FA --out DIR [--seed N]
                                          per query: <q>.faa (query), <q>.hits.faa, <q>.null.faa (hit residues shuffled)
  summarize    --manifest TSV --queries-faa FA --hits-faa FA --real-dir DIR --null-dir DIR --out DIR
                                          classify every hit; the training rule on SUPPORTED hits only
  (the searches between them are `phmmer -E 10 --domE 10 --domtblout`, see pr6_drift_check.sbatch)

FREEZE (written before any result is read; every number is a PROPOSAL for the owner, not a decision):
  E_SUPPORT        1e-5   the best domain i-Evalue of the pairwise phmmer query-vs-hit must be <= this (the inclusion E-value
                          of the searches)
  COVER_MIN        0.5    the aligned part of the query OR of the hit must be >= this fraction of its length (a hit that
                          matches only a short stretch of both is not a homolog of the protein)
  NULL_RATIO_MAX   0.5    the same pairwise search against the hit with its residues SHUFFLED (same composition, order
                          destroyed; seed fixed) must score LESS than this fraction of the real full-sequence score;
                          otherwise the score is explained by composition (`bias_driven`)
  MIN_SEQUENCES / MIN_CLADES  5 / 2   the training rule of decision 7, applied to SUPPORTED hits only; the 2-clade condition
                          counts named clades only (hits without a fine clade are not a clade)
  Composition columns (fraction of S/T/N/Q/G/P, low-complexity fraction) are DESCRIPTIVE: the classifier does not use them.

Classes of one hit:
  supported_homolog  E_SUPPORT met AND COVER_MIN met AND the null score ratio < NULL_RATIO_MAX
  bias_driven        E_SUPPORT met AND the null score ratio >= NULL_RATIO_MAX (a shuffled hit scores about as well)
  unsupported        no pairwise hit at E_SUPPORT, or the coverage is below COVER_MIN
Verdict of a query: `supported_rule_met` when the supported hits meet the rule, else `not_supported` (the query stays not
buildable). A query is never promoted by this script; the owner decides.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

E_SUPPORT = 1e-5
COVER_MIN = 0.5
NULL_RATIO_MAX = 0.5
MIN_SEQUENCES = 5
MIN_CLADES = 2
NO_CLADE = "(no clade)"
DEFAULT_QUERIES = ("P32323", "Q07730", "A0ACE4KFX7")  # aga1, ece1, hyd1_like
LOW_COMPLEXITY_WINDOW = 12
LOW_COMPLEXITY_BITS = (
    2.2  # windows with Shannon entropy below this many bits count as low complexity (descriptive)
)
BIAS_RESIDUES = set("STNQGP")
AA = set("ACDEFGHIKLMNPQRSTVWY")
MANIFEST_HEADER = ["query", "hit", "clade_fine", "taxid", "evalue", "score", "in_phmmer"]
CLASS_HEADER = [
    "query",
    "hit",
    "clade_fine",
    "in_phmmer",
    "class",
    "best_domain_ievalue",
    "query_cover",
    "hit_cover",
    "real_score",
    "null_score",
    "null_ratio",
    "hit_len",
    "hit_stnqgp",
    "hit_low_complexity",
    "query_stnqgp",
    "query_low_complexity",
]
SUMMARY_HEADER = [
    "query",
    "n_jackhmmer_outside",
    "n_supported",
    "n_bias_driven",
    "n_unsupported",
    "supported_named_clades",
    "supported_unassigned",
    "verdict",
    "reason",
]


# ---------------------------------------------------------------------------------------------------- helpers


def read_tsv(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    """{first header token: sequence}. A repeated id is an error (a hit must be unique)."""
    out: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(parts)
                name = line[1:].split()[0]
                if name in out:
                    raise ValueError(f"{path}: id {name} occurs twice")
                parts = []
            elif name is not None:
                parts.append(line.strip())
    if name is not None:
        out[name] = "".join(parts)
    return out


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(v), 60):
                fh.write(v[i : i + 60] + "\n")


def stnqgp_fraction(seq: str) -> float:
    return sum(1 for c in seq if c in BIAS_RESIDUES) / len(seq) if seq else 0.0


def low_complexity_fraction(seq: str) -> float:
    """Fraction of residues inside at least one window of LOW_COMPLEXITY_WINDOW residues whose Shannon entropy is below
    LOW_COMPLEXITY_BITS (a simple SEG-like measure; descriptive only)."""
    n, w = len(seq), LOW_COMPLEXITY_WINDOW
    if n < w:
        return 0.0
    covered = [False] * n
    for i in range(n - w + 1):
        win = seq[i : i + w]
        counts = Counter(win)
        h = -sum(c / w * math.log2(c / w) for c in counts.values())
        if h < LOW_COMPLEXITY_BITS:
            for j in range(i, i + w):
                covered[j] = True
    return sum(covered) / n


def union_length(intervals: list[tuple[int, int]]) -> int:
    """Total length of the union of closed integer intervals (1-based, inclusive)."""
    total = 0
    cur_s = cur_e = None
    for s, e in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


# ---------------------------------------------------------------------------------------------------- extract-ids


def extract_ids(record: Path, out: Path, queries: tuple[str, ...]) -> tuple[list[list], list[str]]:
    """The jackhmmer outside fungal hits of `queries` from the COMMITTED record, with a flag for those the one-pass
    phmmer search also found."""
    hits = read_tsv(record / "gap_homologs.tsv.gz")
    phmmer = {
        (h["query"], h["hit_accession"])
        for h in hits
        if h["method"] == "phmmer" and h["where"] == "outside_clade"
    }
    rows = []
    for h in sorted(hits, key=lambda h: (h["query"], h["hit_accession"])):
        if h["method"] != "jackhmmer" or h["where"] != "outside_clade" or h["query"] not in queries:
            continue
        rows.append(
            [
                h["query"],
                h["hit_accession"],
                h["clade_fine"] or NO_CLADE,
                h["taxid"],
                h["evalue"],
                h["score"],
                "yes" if (h["query"], h["hit_accession"]) in phmmer else "no",
            ]
        )
    ids = sorted({r[1] for r in rows})
    return rows, ids


# ---------------------------------------------------------------------------------------------------- split / shuffle


def shuffled(seq: str, rng: random.Random) -> str:
    letters = list(seq)
    rng.shuffle(letters)
    return "".join(letters)


def split_inputs(
    manifest: Path, hits_faa: Path, queries_faa: Path, out: Path, seed: int
) -> dict[str, int]:
    rows = read_tsv(manifest)
    hits = read_fasta(hits_faa)
    queries = read_fasta(queries_faa)
    by_query: dict[str, list[str]] = defaultdict(list)
    for r in rows:
        by_query[r["query"]].append(r["hit"])
    out.mkdir(parents=True, exist_ok=True)
    counts = {}
    for q in sorted(by_query):
        if q not in queries:
            raise ValueError(f"query {q} is not in {queries_faa}")
        missing = [h for h in by_query[q] if h not in hits]
        if missing:
            raise ValueError(
                f"{len(missing)} hits of {q} are not in {hits_faa} (first: {missing[0]})"
            )
        rng = random.Random(
            f"{seed}:{q}"
        )  # one stream per query: the shuffle does not depend on the other queries
        real = {h: hits[h] for h in by_query[q]}
        null = {h: shuffled(s, rng) for h, s in real.items()}
        write_fasta(out / f"{q}.faa", {q: queries[q]})
        write_fasta(out / f"{q}.hits.faa", real)
        write_fasta(out / f"{q}.null.faa", null)
        counts[q] = len(real)
    return counts


# ---------------------------------------------------------------------------------------------------- domtblout


def parse_domtbl(path: Path) -> dict[str, dict]:
    """Per target of a `phmmer --domtblout` of ONE query: best domain i-Evalue, the query and target intervals of all
    domains at or below E 10 (union), the full-sequence score, qlen and tlen. The file must end with `# [ok]`."""
    lines = path.read_text().splitlines()
    if not lines or lines[-3:].count("# [ok]") != 1:
        raise ValueError(f"{path}: no '# [ok]' trailer (the search did not finish)")
    out: dict[str, dict] = {}
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 22:
            raise ValueError(f"{path}: a domain line has {len(f)} fields, expected at least 22")
        target, tlen, qlen = f[0], int(f[2]), int(f[5])
        full_score, i_e = float(f[7]), float(f[12])
        q_iv = (
            int(f[15]),
            int(f[16]),
        )  # hmm from, hmm to = positions on the QUERY (phmmer builds the model from it)
        t_iv = (int(f[17]), int(f[18]))  # ali from, ali to = positions on the target
        d = out.setdefault(
            target,
            {"best_ie": i_e, "q": [], "t": [], "score": full_score, "qlen": qlen, "tlen": tlen},
        )
        d["best_ie"] = min(d["best_ie"], i_e)
        d["q"].append(q_iv)
        d["t"].append(t_iv)
    return out


def classify(row: dict, real: dict | None, null: dict | None) -> dict:
    """The class of one hit under the FREEZE of the module docstring."""
    out = {
        "best_domain_ievalue": "",
        "query_cover": "",
        "hit_cover": "",
        "real_score": "",
        "null_score": "",
        "null_ratio": "",
    }
    if real is None or real["best_ie"] > E_SUPPORT:
        if real is not None:
            out["best_domain_ievalue"] = f"{real['best_ie']:.3g}"
            out["real_score"] = f"{real['score']:.1f}"
        out["class"] = "unsupported"
        return out
    q_cover = union_length(real["q"]) / real["qlen"]
    h_cover = union_length(real["t"]) / real["tlen"]
    null_score = null["score"] if null is not None else 0.0
    ratio = null_score / real["score"] if real["score"] > 0 else 1.0
    out.update(
        best_domain_ievalue=f"{real['best_ie']:.3g}",
        query_cover=f"{q_cover:.3f}",
        hit_cover=f"{h_cover:.3f}",
        real_score=f"{real['score']:.1f}",
        null_score=f"{null_score:.1f}",
        null_ratio=f"{ratio:.3f}",
    )
    if ratio >= NULL_RATIO_MAX:
        out["class"] = "bias_driven"
    elif max(q_cover, h_cover) < COVER_MIN:
        out["class"] = "unsupported"
    else:
        out["class"] = "supported_homolog"
    return out


def verdict(rows: list[dict]) -> tuple[str, str, int, int]:
    """The training rule on the supported hits: (verdict, reason, named clades, unassigned sequences)."""
    sup = [r for r in rows if r["class"] == "supported_homolog"]
    named = {r["clade_fine"] for r in sup if r["clade_fine"] != NO_CLADE}
    unassigned = sum(1 for r in sup if r["clade_fine"] == NO_CLADE)
    if len(sup) >= MIN_SEQUENCES and len(named) >= MIN_CLADES:
        return "supported_rule_met", "", len(named), unassigned
    if len(sup) < MIN_SEQUENCES:
        reason = f"only {len(sup)} supported outside sequences (minimum {MIN_SEQUENCES})"
    else:
        reason = f"supported hits fall in {len(named)} named clade (minimum {MIN_CLADES})"
    return "not_supported", reason, len(named), unassigned


def summarize(
    manifest: Path, queries_faa: Path, hits_faa: Path, real_dir: Path, null_dir: Path, out: Path
) -> list[list]:
    rows = read_tsv(manifest)
    queries = read_fasta(queries_faa)
    hits = read_fasta(hits_faa)
    out.mkdir(parents=True, exist_ok=True)
    qs = sorted({r["query"] for r in rows})
    classified: list[dict] = []
    summary = []
    for q in qs:
        if q not in queries:
            raise ValueError(f"query {q} is not in {queries_faa}")
        real = parse_domtbl(real_dir / f"{q}.dom")
        null = parse_domtbl(null_dir / f"{q}.dom")
        qseq = queries[q]
        mine = []
        for r in (r for r in rows if r["query"] == q):
            seq = hits.get(r["hit"])
            if seq is None:
                raise ValueError(f"hit {r['hit']} of {q} is not in {hits_faa}")
            c = classify(r, real.get(r["hit"]), null.get(r["hit"]))
            c.update(
                query=q,
                hit=r["hit"],
                clade_fine=r["clade_fine"],
                in_phmmer=r["in_phmmer"],
                hit_len=len(seq),
                hit_stnqgp=f"{stnqgp_fraction(seq):.3f}",
                hit_low_complexity=f"{low_complexity_fraction(seq):.3f}",
                query_stnqgp=f"{stnqgp_fraction(qseq):.3f}",
                query_low_complexity=f"{low_complexity_fraction(qseq):.3f}",
            )
            mine.append(c)
        classified += mine
        v, reason, n_named, n_un = verdict(mine)
        cnt = Counter(c["class"] for c in mine)
        summary.append(
            [
                q,
                len(mine),
                cnt["supported_homolog"],
                cnt["bias_driven"],
                cnt["unsupported"],
                n_named,
                n_un,
                v,
                reason,
            ]
        )
    write_tsv(
        out / "drift_hits.tsv", CLASS_HEADER, [[c[k] for k in CLASS_HEADER] for c in classified]
    )
    write_tsv(out / "drift_summary.tsv", SUMMARY_HEADER, summary)
    per_clade = defaultdict(Counter)
    for c in classified:
        per_clade[(c["query"], c["clade_fine"])][c["class"]] += 1
    write_tsv(
        out / "drift_per_clade.tsv",
        ["query", "clade_fine", "supported_homolog", "bias_driven", "unsupported"],
        [
            [q, cl, v["supported_homolog"], v["bias_driven"], v["unsupported"]]
            for (q, cl), v in sorted(per_clade.items())
        ],
    )
    lines = [
        "PR 6 drift check (DRAFT summary; descriptive only: no model, no call, no status, no tier)"
    ]
    lines.append(
        f"freeze: E <= {E_SUPPORT:g}, cover >= {COVER_MIN}, null ratio < {NULL_RATIO_MAX}, rule {MIN_SEQUENCES} sequences in "
        f"{MIN_CLADES} named clades (all PROPOSALS for the owner)"
    )
    for s in summary:
        lines.append(
            f"{s[0]}: {s[1]} jackhmmer outside hits; supported {s[2]}, bias_driven {s[3]}, unsupported {s[4]}; "
            f"supported named clades {s[5]}; {s[7]}" + (f" ({s[8]})" if s[8] else "")
        )
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    return summary


# ---------------------------------------------------------------------------------------------------- main


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract-ids")
    e.add_argument("--record", type=Path, required=True)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--queries", default=",".join(DEFAULT_QUERIES))
    s = sub.add_parser("split")
    s.add_argument("--manifest", type=Path, required=True)
    s.add_argument("--hits-faa", type=Path, required=True)
    s.add_argument("--queries-faa", type=Path, required=True)
    s.add_argument("--out", type=Path, required=True)
    s.add_argument("--seed", type=int, default=1)
    m = sub.add_parser("summarize")
    m.add_argument("--manifest", type=Path, required=True)
    m.add_argument("--queries-faa", type=Path, required=True)
    m.add_argument("--hits-faa", type=Path, required=True)
    m.add_argument("--real-dir", type=Path, required=True)
    m.add_argument("--null-dir", type=Path, required=True)
    m.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    if a.out.exists() and any(a.out.iterdir()):
        print(f"FATAL: {a.out} exists and is not empty; nothing is overwritten", file=sys.stderr)
        return 1
    if a.cmd == "extract-ids":
        rows, ids = extract_ids(a.record, a.out, tuple(x for x in a.queries.split(",") if x))
        a.out.mkdir(parents=True, exist_ok=True)
        write_tsv(a.out / "hits_manifest.tsv", MANIFEST_HEADER, rows)
        (a.out / "ids.txt").write_text("".join(i + "\n" for i in ids))
        per = Counter(r[0] for r in rows)
        print(
            f"{len(rows)} hits of {len(per)} queries, {len(ids)} distinct ids: {dict(sorted(per.items()))}"
        )
    elif a.cmd == "split":
        counts = split_inputs(a.manifest, a.hits_faa, a.queries_faa, a.out, a.seed)
        print(f"split: {counts}")
    else:
        summary = summarize(a.manifest, a.queries_faa, a.hits_faa, a.real_dir, a.null_dir, a.out)
        print((a.out / "summary.txt").read_text(), end="")
        return 0 if summary else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
