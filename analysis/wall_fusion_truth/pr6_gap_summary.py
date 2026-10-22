#!/usr/bin/env python3
"""PR 6 gap search: summary of the outside-homolog search for the clusters with no Pfam hit, and the Ecm33 check.

DRAFT, not run. Stdlib only. Run with `/usr/bin/python3.12 -I` on HPCC or `python3 -I` locally.
It reads tables and writes tables. It writes no call, no status and no tier (nothing of
`cellsurface_sorting_hat` is imported), and no sequence. Spec: docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md,
decisions 3, 5 and 7 (steps 6.1 and 6.2): outside homologs by sequence search, how many gap clusters have any, and the
minimum training set.

Inputs
  --targets      custom_targets.tsv of the step 6.0 record (the clusters with no Pfam hit: slice, cluster_id, genes,
                 clades_fine, accessions). The held-out unit of a target is its own fine clade (all four are clade-only).
  --phmmer-dir / --jackhmmer-dir
                 one `<accession>.tbl` per target query: the `--tblout` of phmmer (one pass) and of jackhmmer
                 (iterations) against the UniProt set. The description column of a UniProt hit carries `OX=<taxid>`.
  --nodes, --names  the NCBI taxdump (nodes.dmp, names.dmp) that turns `OX=` into a clade.
  --ecm33-domtbl, --ecm33-model, --ecm33-lengths  (optional) the Ecm33 check: `hmmsearch --domtblout` of the model
                 PF12454 against Q5AGC4 and P38248, the fetched model file (for GA and LENG) and a two-column table
                 sequence<TAB>length.

Clades are those of `map_kiss.clades_from_lineage` (fine level: Saccharomycotina for the yeasts, the class for
Pezizomycotina and Basidiomycota, the phylum otherwise). Only fungi (taxon 4751 among the ancestors) count toward the
minimum training set; hits from other kingdoms are counted apart and are descriptive.

Minimum training set (decision 7, provisional): AT LEAST 5 distinct outside sequences IN TOTAL, spread over AT LEAST 2
clades outside the held-out clade. (The stricter reading, 5 sequences in each of 2 clades, is a descriptive column.)
Status of a target:
  buildable                 the rule is met by the phmmer (one pass) hits
  buildable_iterated_only   met only by the jackhmmer hits (iteration may drift: a flag for the owner, not a decision)
  not_buildable             not met; the reason is written. It is a known gap, not a model.

Outputs (OUTPUTS): gap_per_query.tsv, gap_clades.tsv, gap_homologs.tsv (one line per hit: ids only), ecm33_check.tsv,
summary.txt, provenance.tsv.

Usage: python3 -I analysis/wall_fusion_truth/pr6_gap_summary.py --targets T --phmmer-dir P --jackhmmer-dir J \
           --nodes nodes.dmp --names names.dmp --out OUT [--ecm33-domtbl D --ecm33-model M --ecm33-lengths L]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

FUNGI = 4751
MIN_SEQUENCES = 5  # decision 7, provisional
MIN_CLADES = 2
METHODS = ("phmmer", "jackhmmer")
OUTPUTS = (
    "gap_per_query.tsv",
    "gap_clades.tsv",
    "gap_homologs.tsv",
    "ecm33_check.tsv",
    "summary.txt",
    "provenance.tsv",
)
ECM33_MODEL_ACC = "PF12454"
SCER_ECM33 = "P38248"  # S. cerevisiae Ecm33 (406 aa, report 07 row 22)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------- taxonomy


def load_nodes(path: Path) -> tuple[dict[int, int], dict[int, str]]:
    """{taxid: parent}, {taxid: rank} from nodes.dmp (fields separated by '\\t|\\t')."""
    parent: dict[int, int] = {}
    rank: dict[int, str] = {}
    with path.open(encoding="utf-8") as fh:
        for ln in fh:
            f = ln.split("\t|\t")
            if len(f) < 3:
                continue
            tid = int(f[0])
            parent[tid] = int(f[1])
            rank[tid] = f[2]
    if not parent:
        raise ValueError(f"{path.name}: no nodes")
    return parent, rank


def load_names(path: Path, wanted: set[int]) -> dict[int, str]:
    """Scientific names of the wanted taxids only (names.dmp is large; it is streamed)."""
    out: dict[int, str] = {}
    with path.open(encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t|\t")
            if len(f) < 4 or not f[3].startswith("scientific name"):
                continue
            tid = int(f[0])
            if tid in wanted:
                out[tid] = f[1]
    return out


def ancestry(tid: int, parent: dict[int, int]) -> list[int]:
    """The taxon and its ancestors up to the root (a taxon that is its own parent ends the walk)."""
    chain, seen = [], set()
    while tid in parent and tid not in seen:
        chain.append(tid)
        seen.add(tid)
        if parent[tid] == tid:
            break
        tid = parent[tid]
    return chain


def clades_from_lineage(lineage: dict) -> tuple[str, str]:
    """(coarse, fine) clade; the same rule as analysis/wall_fusion_truth/map_kiss.py (a test compares the two)."""
    phylum = lineage.get("phylum", ("", ""))[0]
    sub = lineage.get("subphylum", ("", ""))[0]
    cls = lineage.get("class", ("", ""))[0]
    if phylum == "Ascomycota":
        coarse = sub
        fine = cls if sub == "Pezizomycotina" else sub
    elif phylum == "Basidiomycota":
        coarse, fine = phylum, cls
    else:
        coarse = fine = phylum
    return coarse, fine


class Taxonomy:
    def __init__(self, nodes: Path, names: Path, taxids: set[int]):
        self.parent, self.rank = load_nodes(nodes)
        chains = {t: ancestry(t, self.parent) for t in taxids}
        need = {a for c in chains.values() for a in c}
        self.names = load_names(names, need)
        self.chains = chains

    def info(self, tid: int) -> dict:
        """{'fungus': bool, 'clade_coarse', 'clade_fine', 'species_taxid'} of a taxon ('' clade when unknown)."""
        chain = self.chains.get(tid)
        if not chain:
            return {"fungus": False, "clade_coarse": "", "clade_fine": "", "species_taxid": tid}
        lineage = {}
        species = tid
        for a in chain:  # from the taxon upward
            r = self.rank.get(a, "")
            if r in ("phylum", "subphylum", "class") and r not in lineage:
                lineage[r] = (self.names.get(a, ""), a)
            if r == "species" and species == tid:
                species = a
        coarse, fine = clades_from_lineage(lineage)
        return {
            "fungus": FUNGI in chain,
            "clade_coarse": coarse,
            "clade_fine": fine,
            "species_taxid": species,
        }


# ---------------------------------------------------------------------------------------------------- hits

OX = re.compile(r"\bOX=(\d+)")


def uniprot_accession(target: str) -> str:
    """`sp|P12345|NAME_SPECIES` or `tr|A0A..|NAME` -> the accession; any other name is returned unchanged."""
    parts = target.split("|")
    return parts[1] if len(parts) >= 3 and parts[0] in ("sp", "tr") else target


def parse_tblout(path: Path) -> list[dict]:
    """The hits of a phmmer/jackhmmer `--tblout`: target, E-value, score, taxid from `OX=` in the description.

    The file must end with the `# [ok]` trailer (HMMER did not finish otherwise). A line has 18 columns followed by
    the description of the target. A hit without `OX=` has taxid None."""
    lines = path.read_text().splitlines()
    if not lines or lines[-1].strip() != "# [ok]":
        raise ValueError(f"{path.name}: no '# [ok]' trailer (the search did not finish)")
    out = []
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split(None, 18)
        if len(f) < 18:
            raise ValueError(f"{path.name}: short line: {ln[:80]}")
        desc = f[18] if len(f) > 18 else ""
        m = OX.search(desc)
        out.append(
            {
                "target": f[0],
                "accession": uniprot_accession(f[0]),
                "evalue": float(f[4]),
                "score": float(f[5]),
                "taxid": int(m.group(1)) if m else None,
            }
        )
    return out


def classify_hits(query_acc: str, held_out: str, hits: list[dict], tax: Taxonomy) -> list[dict]:
    """One row per distinct hit accession (the best E-value kept), with its clade and where it falls.

    where: self (the query), no_taxid, non_fungal, inside_clade (the held-out clade), outside_clade."""
    best: dict[str, dict] = {}
    for h in hits:
        cur = best.get(h["accession"])
        if cur is None or h["evalue"] < cur["evalue"]:
            best[h["accession"]] = h
    rows = []
    for acc in sorted(best):
        h = best[acc]
        info = tax.info(h["taxid"]) if h["taxid"] is not None else None
        if acc == query_acc:
            where = "self"
        elif info is None:
            where = "no_taxid"
        elif not info["fungus"]:
            where = "non_fungal"
        elif info["clade_fine"] == held_out:
            where = "inside_clade"
        else:
            where = "outside_clade"
        rows.append(
            {
                "accession": acc,
                "evalue": h["evalue"],
                "score": h["score"],
                "taxid": h["taxid"],
                "species_taxid": info["species_taxid"] if info else None,
                "clade_fine": info["clade_fine"] if info else "",
                "where": where,
            }
        )
    return rows


def training_rule(rows: list[dict]) -> dict:
    """The minimum training set of decision 7 on the outside fungal hits of one search."""
    outside = [r for r in rows if r["where"] == "outside_clade"]
    per_clade: dict[str, list[dict]] = defaultdict(list)
    for r in outside:
        per_clade[r["clade_fine"] or "(no clade)"].append(r)
    n_seq = len(outside)
    n_clades = len(per_clade)
    n_species = len({r["species_taxid"] for r in outside})
    met = n_seq >= MIN_SEQUENCES and n_clades >= MIN_CLADES
    strict = sum(1 for v in per_clade.values() if len(v) >= MIN_SEQUENCES) >= MIN_CLADES
    if met:
        reason = ""
    elif n_seq == 0:
        reason = "no outside fungal homolog"
    elif n_seq < MIN_SEQUENCES:
        reason = f"only {n_seq} outside sequences (minimum {MIN_SEQUENCES})"
    else:
        reason = f"outside hits fall in {n_clades} clade (minimum {MIN_CLADES})"
    return {
        "n_outside_seq": n_seq,
        "n_outside_clades": n_clades,
        "n_outside_species": n_species,
        "rule_met": met,
        "strict_per_clade_met": strict,
        "reason": reason,
        "per_clade": {c: len(v) for c, v in sorted(per_clade.items())},
    }


def target_status(rules: dict[str, dict]) -> tuple[str, str]:
    p, j = rules["phmmer"], rules["jackhmmer"]
    if p["rule_met"]:
        return "buildable", ""
    if j["rule_met"]:
        return (
            "buildable_iterated_only",
            "met only after iteration; jackhmmer may have drifted (check before use)",
        )
    return "not_buildable", j["reason"] if j["n_outside_seq"] >= p["n_outside_seq"] else p["reason"]


# ---------------------------------------------------------------------------------------------------- Ecm33 check


def parse_model_header(path: Path) -> dict:
    """NAME, ACC, LENG, GA (sequence, domain) from the head of an HMMER3 model file."""
    out: dict = {}
    for ln in path.read_text().splitlines():
        if ln.startswith("HMM "):
            break
        f = ln.split(None, 1)
        if len(f) == 2:
            out.setdefault(f[0], f[1].strip())
    if "GA" not in out:
        raise ValueError(f"{path.name}: the model has no GA line")
    ga = out["GA"].rstrip(";").split()
    return {
        "name": out.get("NAME", ""),
        "acc": out.get("ACC", "").split(".")[0],
        "leng": int(out["LENG"]),
        "ga_seq": float(ga[0]),
        "ga_dom": float(ga[1]),
    }


def parse_hmmsearch_domtbl(path: Path, model_name: str) -> dict[str, dict]:
    """{sequence: {score_seq, evalue_seq, best_dom_score, best_dom_ievalue, hmm_span, ali_span, tlen}}.

    hmmsearch orientation: column 1 is the sequence and column 4 the model (hmmscan's is the reverse and is refused)."""
    lines = path.read_text().splitlines()
    if not lines or lines[-1].strip() != "# [ok]":
        raise ValueError(f"{path.name}: no '# [ok]' trailer (hmmsearch did not finish)")
    out: dict[str, dict] = {}
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split()
        if len(f) < 22:
            raise ValueError(f"{path.name}: short domain line: {ln[:80]}")
        if f[3] != model_name:
            raise ValueError(
                f"{path.name}: column 4 is {f[3]!r}, not the model {model_name!r} (hmmscan orientation?)"
            )
        seq = f[0]
        d = out.setdefault(
            seq,
            {
                "score_seq": float(f[7]),
                "evalue_seq": float(f[6]),
                "best_dom_score": float(f[13]),
                "best_dom_ievalue": float(f[12]),
                "hmm_from": int(f[15]),
                "hmm_to": int(f[16]),
                "ali_from": int(f[17]),
                "ali_to": int(f[18]),
                "tlen": int(f[2]),
            },
        )
        if float(f[13]) > d["best_dom_score"]:
            d.update(
                best_dom_score=float(f[13]),
                best_dom_ievalue=float(f[12]),
                hmm_from=int(f[15]),
                hmm_to=int(f[16]),
                ali_from=int(f[17]),
                ali_to=int(f[18]),
            )
    return out


def ecm33_rows(model: dict, hits: dict[str, dict], lengths: dict[str, int]) -> list[list]:
    """One row per queried sequence: the score against the GA thresholds of PF12454 and the plausibility of its length."""
    scer = lengths.get(SCER_ECM33)
    rows = []
    for seq in sorted(lengths):
        n = lengths[seq]
        ratio = f"{n / scer:.3f}" if scer else ""
        h = hits.get(seq)
        if h is None:
            rows.append(
                [
                    seq,
                    n,
                    ratio,
                    model["acc"],
                    model["ga_seq"],
                    model["ga_dom"],
                    "",
                    "",
                    "",
                    "",
                    "",
                    "no_hit_at_E10",
                ]
            )
            continue
        cover = (h["hmm_to"] - h["hmm_from"] + 1) / model["leng"]
        if h["score_seq"] >= model["ga_seq"] and h["best_dom_score"] >= model["ga_dom"]:
            verdict = "pass_GA"
        else:
            verdict = "below_GA"
        rows.append(
            [
                seq,
                n,
                ratio,
                model["acc"],
                model["ga_seq"],
                model["ga_dom"],
                h["score_seq"],
                h["evalue_seq"],
                h["best_dom_score"],
                h["best_dom_ievalue"],
                f"{cover:.3f}",
                verdict,
            ]
        )
    return rows


ECM33_HEADER = [
    "sequence",
    "length",
    "length_ratio_vs_P38248",
    "model",
    "ga_seq",
    "ga_dom",
    "score_seq",
    "evalue_seq",
    "best_dom_score",
    "best_dom_ievalue",
    "hmm_cover_of_model",
    "verdict",
]


# ---------------------------------------------------------------------------------------------------- main


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--targets", type=Path, required=True)
    ap.add_argument("--phmmer-dir", type=Path, required=True)
    ap.add_argument("--jackhmmer-dir", type=Path, required=True)
    ap.add_argument("--nodes", type=Path, required=True)
    ap.add_argument("--names", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ecm33-domtbl", type=Path)
    ap.add_argument("--ecm33-model", type=Path)
    ap.add_argument("--ecm33-lengths", type=Path)
    a = ap.parse_args(argv)
    if a.out.exists() and any(a.out.iterdir()):
        sys.exit(f"refusing: {a.out} exists and is not empty")
    targets = read_tsv(a.targets)
    if not targets:
        sys.exit("no targets")
    dirs = {"phmmer": a.phmmer_dir, "jackhmmer": a.jackhmmer_dir}
    hits: dict[tuple[str, str], list[dict]] = {}
    taxids: set[int] = set()
    for t in targets:
        for acc in t["accessions"].split(","):
            for m in METHODS:
                path = dirs[m] / f"{acc}.tbl"
                if not path.exists():
                    sys.exit(f"missing {path}")
                hits[(acc, m)] = parse_tblout(path)
                taxids |= {h["taxid"] for h in hits[(acc, m)] if h["taxid"] is not None}
    tax = Taxonomy(a.nodes, a.names, taxids)

    per_query, clade_rows, hom_rows = [], [], []
    counts = defaultdict(int)
    for t in targets:
        for acc in t["accessions"].split(","):
            held = t["clades_fine"]
            rules, classified = {}, {}
            for m in METHODS:
                classified[m] = classify_hits(acc, held, hits[(acc, m)], tax)
                rules[m] = training_rule(classified[m])
                for r in classified[m]:
                    hom_rows.append(
                        [
                            acc,
                            m,
                            r["accession"],
                            r["evalue"],
                            r["score"],
                            r["taxid"] if r["taxid"] is not None else "",
                            r["clade_fine"],
                            r["where"],
                        ]
                    )
                for clade, n in rules[m]["per_clade"].items():
                    clade_rows.append([acc, m, clade, n])
            status, reason = target_status(rules)
            counts[status] += 1
            mode = "leave-clade-out" if t["slice"] == "wall" else "leave-cluster-out (partial)"
            w = {m: defaultdict(int) for m in METHODS}
            for m in METHODS:
                for r in classified[m]:
                    w[m][r["where"]] += 1
            per_query.append(
                [
                    t["slice"],
                    t["cluster_id"],
                    t["genes"],
                    acc,
                    held,
                    mode,
                    *(
                        v
                        for m in METHODS
                        for v in (
                            rules[m]["n_outside_seq"],
                            rules[m]["n_outside_clades"],
                            rules[m]["n_outside_species"],
                            "yes" if rules[m]["rule_met"] else "no",
                            "yes" if rules[m]["strict_per_clade_met"] else "no",
                            w[m]["inside_clade"],
                            w[m]["non_fungal"],
                        )
                    ),
                    status,
                    reason,
                ]
            )
    header = ["slice", "cluster_id", "genes", "query", "held_out_clade", "mode"]
    for m in METHODS:
        header += [
            f"{m}_outside_seq",
            f"{m}_outside_clades",
            f"{m}_outside_species",
            f"{m}_rule_met",
            f"{m}_strict_5_in_2_clades",
            f"{m}_inside_clade_seq",
            f"{m}_non_fungal_seq",
        ]
    header += ["status", "reason"]

    a.out.mkdir(parents=True, exist_ok=True)
    write_tsv(a.out / "gap_per_query.tsv", header, per_query)
    write_tsv(
        a.out / "gap_clades.tsv", ["query", "method", "clade_fine", "n_outside_seq"], clade_rows
    )
    write_tsv(
        a.out / "gap_homologs.tsv",
        ["query", "method", "hit_accession", "evalue", "score", "taxid", "clade_fine", "where"],
        hom_rows,
    )
    ecm_lines: list[str] = []
    if a.ecm33_domtbl and a.ecm33_model and a.ecm33_lengths:
        model = parse_model_header(a.ecm33_model)
        if model["acc"] != ECM33_MODEL_ACC:
            sys.exit(f"the fetched model is {model['acc']}, not {ECM33_MODEL_ACC}")
        hits_e = parse_hmmsearch_domtbl(a.ecm33_domtbl, model["name"])
        lengths = {r[0]: int(r[1]) for r in csv.reader(a.ecm33_lengths.open(), delimiter="\t") if r}
        rows = ecm33_rows(model, hits_e, lengths)
        write_tsv(a.out / "ecm33_check.tsv", ECM33_HEADER, rows)
        for r in rows:
            ecm_lines.append(
                f"ecm33 check {r[0]} ({r[1]} aa): {r[11]}"
                + (f"; score {r[6]} vs GA {r[4]} (domain {r[8]} vs {r[5]})" if r[6] != "" else "")
            )
    else:
        write_tsv(a.out / "ecm33_check.tsv", ECM33_HEADER, [])
    lines = [
        "PR 6 gap search: outside homologs of the clusters with no Pfam hit (DRAFT summary)",
        f"targets: {len(per_query)} queries; "
        + "; ".join(f"{k} {v}" for k, v in sorted(counts.items())),
        f"rule (decision 7, provisional): >= {MIN_SEQUENCES} outside fungal sequences in total over >= {MIN_CLADES} clades",
        "descriptive only: no model is built, no call, no status, no tier; not_buildable is a known gap",
        *ecm_lines,
    ]
    for r in per_query:
        lines.append(
            f"{r[3]} ({r[2]}, {r[0]}, held-out {r[4]}): {r[-2]}" + (f" ({r[-1]})" if r[-1] else "")
        )
    (a.out / "summary.txt").write_text("\n".join(lines) + "\n")
    write_tsv(
        a.out / "provenance.tsv",
        ["item", "value"],
        [
            ["targets_sha256", sha256(a.targets)],
            ["nodes_sha256", sha256(a.nodes)],
            ["names_sha256", sha256(a.names)],
            ["min_sequences", MIN_SEQUENCES],
            ["min_clades", MIN_CLADES],
            ["script", Path(__file__).name],
        ],
    )
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
