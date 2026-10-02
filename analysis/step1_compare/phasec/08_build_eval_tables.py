#!/usr/bin/env python3
"""Phase C step 1: the evaluation table (Phase C spec 3.1, 3.2, 3.4 items 1, 2 and 5c).

Reads from $STEP1_WORKDIR: truth_set_triaged.tsv.gz (03), keyword_tier.tsv.gz (04),
phaseb/features.tsv.gz, phaseb/features_unique.tsv.gz (07), phaseb/unique_sequences.tsv.gz
(05), and the run JSONs d8_run.json, keyword_tier_run.json, phaseb/features_run.json,
phaseb/emb/embedding_run.json. Reads species.tsv, data/curated/adhesins/eurotiomycetes_seeds.tsv
and phasec/tc_taxon_clades.tsv. Writes to $STEP1_WORKDIR/phasec/:

  eval_table.tsv.gz        one row per unique sequence: GO truth genes of non-alternate sources
                           with a label other than `unlabelled`, and the T-c rows that survive
                           the precedence rule (ruling C-6)
  eval_literature.tsv      one row per literature seed with an accession (spec 3.2)
  eval_dedupe_log.tsv      one row per dropped or reclassified member
  eval_sequences.fasta.gz  the sequences of the table and the literature rows (MMseqs2 input)
  build_run.json           counts, input hashes, all_sources, git commit, library versions

STOP (exit 2, no output): a run JSON without all_sources: true; features_run.json
input_sha256["unique_sequences.tsv.gz"] differs from embedding_run.json
unique_sequences_sha256 or from the current file (Phase B review item M3); the truth set hash
differs between d8_run.json, features_run.json and keyword_tier_run.json; the current
truth_set_triaged.tsv.gz differs from the file 07 read; a labelled truth row whose only evidence
is IEA, or IEA in an evidence column; a T-c taxon_id without a row in tc_taxon_clades.tsv; a
source without a row in species.tsv or with another role there; a table or literature hash
without a valid feature row (empty or non-finite sp_prob, gpi_prob, ser_thr_frac; unknown
gpi_call) or without a sequence.
"""

import argparse
import gzip
import io
import sys
from collections import Counter
from pathlib import Path

import dedupe
import evalio
import labelmap
import manifest
import paths
import runinfo
import truth_table

OUTPUT_NAMES = (
    "eval_table.tsv.gz",
    "eval_literature.tsv",
    "eval_dedupe_log.tsv",
    "eval_sequences.fasta.gz",
    "build_run.json",
)
LITERATURE_COLUMNS = (
    "accession",
    "gene",
    "lit_class",
    "moonlighting",
    "species",
    "order",
    "seq_sha256",
    "length",
    "emb_row",
    "emb_cterm_row",
    "literature_positive",
)
GPI_CALLS = ("highly_probable", "probable", "weakly", "none", "too_short")
EVIDENCE_COLUMNS = ("surface_evidence", "internal_evidence", "secretory_evidence")
DEFAULT_SEEDS = Path("data/curated/adhesins/eurotiomycetes_seeds.tsv")


def check_chain(work: Path) -> dict:
    """Run JSON checks. Return the hashes for build_run.json."""
    work = Path(work)
    fr = evalio.read_json(work / "phaseb" / "features_run.json")
    er = evalio.read_json(work / "phaseb" / "emb" / "embedding_run.json")
    d8 = evalio.read_json(work / "d8_run.json")
    kt = evalio.read_json(work / "keyword_tier_run.json")
    for name, log in (
        ("features_run.json", fr),
        ("d8_run.json", d8),
        ("keyword_tier_run.json", kt),
    ):
        if log.get("all_sources") is not True:
            raise evalio.StopError(f"{name} does not say all_sources: true")
    unique = work / "phaseb" / "unique_sequences.tsv.gz"
    recorded = fr.get("input_sha256", {}).get("unique_sequences.tsv.gz")
    if recorded != er.get("unique_sequences_sha256"):
        raise evalio.StopError(
            f"features_run.json input_sha256['unique_sequences.tsv.gz'] ({recorded}) differs "
            f"from embedding_run.json unique_sequences_sha256 "
            f"({er.get('unique_sequences_sha256')}); features and embeddings come from "
            "different unique sequence sets; re-run 07 and the assembly"
        )
    current = manifest.sha256_file(unique)
    if recorded != current:
        raise evalio.StopError(
            f"phaseb/unique_sequences.tsv.gz ({current}) differs from the file 07 read "
            f"({recorded}); re-run 07 and the assembly"
        )
    truth_hashes = {
        "d8_run.json": d8.get("truth_set_sha256"),
        "features_run.json": fr.get("truth_set_sha256"),
        "keyword_tier_run.json": kt.get("truth_set_sha256"),
    }
    if len(set(truth_hashes.values())) != 1 or None in truth_hashes.values():
        raise evalio.StopError(f"truth_set_sha256 differs between run JSONs: {truth_hashes}")
    triaged = work / "truth_set_triaged.tsv.gz"
    triaged_sha = manifest.sha256_file(triaged)
    if fr.get("input_sha256", {}).get("truth_set_triaged.tsv.gz") != triaged_sha:
        raise evalio.StopError(
            "truth_set_triaged.tsv.gz differs from the file that 07 read; re-run 07"
        )
    return {
        "truth_set_sha256": truth_hashes["d8_run.json"],
        "unique_sequences_sha256": current,
        "truth_set_triaged_sha256": triaged_sha,
    }


def iea_problems(truth_rows) -> list[str]:
    """Labelled rows whose evidence is IEA only, and IEA in a label evidence column."""
    bad = []
    for r in truth_rows:
        key = f"{r['source_id']}:{r['gene_id']}"
        if r["label"] == "unlabelled":
            continue
        codes = set(r["evidence_codes"].split(",")) - {""}
        if not codes - {"IEA"}:
            bad.append(f"{key} has label {r['label']} with IEA evidence only")
        for col in EVIDENCE_COLUMNS:
            if "IEA" in r[col].split(","):
                bad.append(f"{key} has IEA in {col}")
    return bad


def read_tc_clades(path: Path) -> dict[str, str]:
    rows = truth_table.read_tsv(path)
    out = {}
    for r in rows:
        if r["taxon_id"] in out:
            raise evalio.StopError(f"{path}: taxon_id {r['taxon_id']} occurs twice")
        if not r["clade"]:
            raise evalio.StopError(f"{path}: taxon_id {r['taxon_id']} has an empty clade")
        out[r["taxon_id"]] = r["clade"]
    return out


def check_tc_clades(kw_rows, clades: dict[str, str]) -> None:
    missing = sorted({r["taxon_id"] for r in kw_rows} - set(clades))
    if missing:
        raise evalio.StopError(
            f"tc_taxon_clades.tsv has no row for taxon_id {', '.join(missing)}; add the clade "
            "(or `other`) for every taxon_id of keyword_tier.tsv.gz"
        )


def go_members(truth_rows, feature_rows, species_rows) -> tuple[list[dict], Counter]:
    """GO members of non-alternate sources with a label other than unlabelled.

    Return (members, count of labelled genes without a sequence per source)."""
    species = {r["source_id"]: r for r in species_rows}
    truth = {(r["source_id"], r["gene_id"]): r for r in truth_rows}
    for r in truth_rows:
        s = species.get(r["source_id"])
        if s is None:
            raise evalio.StopError(f"source {r['source_id']} has no row in species.tsv")
        if s["role"] != r["role"]:
            raise evalio.StopError(
                f"source {r['source_id']}: role {r['role']!r} in the truth set, "
                f"{s['role']!r} in species.tsv"
            )
    members, seen = [], set()
    for f in feature_rows:
        if f["set_id"] != "truth":
            continue
        key = (f["source_id"], f["gene_id"])
        t = truth.get(key)
        if t is None:
            raise evalio.StopError(f"truth member {key[0]}:{key[1]} has no truth row")
        seen.add(key)
        if t["role"] == "alternate_file" or t["label"] == "unlabelled":
            continue
        members.append(
            {
                "seq_sha256": f["seq_sha256"],
                "class": labelmap.class_of(t["label"], t["d8_class"]),
                "label": t["label"],
                "subset": t["subset"],
                "stratum": labelmap.stratum_of(t["label"], t["subset"], t["d8_class"]),
                "d8_class": t["d8_class"],
                "homology_only": t["homology_only"],
                "internal_evidence_htp_only": t["internal_evidence_htp_only"],
                "source_id": t["source_id"],
                "gene_id": t["gene_id"],
                "species": t["species"],
                "role": t["role"],
                "clade": species[t["source_id"]]["in_clade"],
                "taxon_id": t["taxon_id"],
                "length": f["length"],
                "emb_row": f["emb_row"],
                "emb_cterm_row": f["emb_cterm_row"],
            }
        )
    unmatched = Counter(
        r["source_id"]
        for key, r in truth.items()
        if key not in seen and r["role"] != "alternate_file" and r["label"] != "unlabelled"
    )
    return members, unmatched


def kw_features(feature_rows) -> dict[str, dict]:
    """uniprot_kw members by accession (gene_id)."""
    return {f["gene_id"]: f for f in feature_rows if f["set_id"] == "uniprot_kw"}


def tc_members(kw_rows, kw_by_acc: dict, clades: dict) -> list[dict]:
    out = []
    for r in kw_rows:
        f = kw_by_acc.get(r["accession"])
        if f is None or f["seq_sha256"] != r["seq_sha256"]:
            raise evalio.StopError(
                f"T-c row {r['accession']} has no uniprot_kw feature row with its seq_sha256; "
                "re-run 05 and 07 after 04"
            )
        out.append(
            {
                **{k: r[k] for k in ("accession", "genome", "taxon_id", "seq_sha256")},
                "clade": clades[r["taxon_id"]],
                "length": f["length"],
                "emb_row": f["emb_row"],
                "emb_cterm_row": f["emb_cterm_row"],
            }
        )
    return out


def read_seeds(path: Path) -> list[dict]:
    import csv

    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def literature_rows(seeds, kw_by_acc: dict) -> tuple[list[dict], list[str]]:
    """Rows with an accession. Positive: a sequence in the keyword set and moonlighting not YES.

    Return (rows, genes without an accession)."""
    rows, no_accession = [], []
    for s in seeds:
        query = (s.get("uniprot_query") or "").strip()
        if not query.startswith("accession:"):
            no_accession.append(s["gene"])
            continue
        acc = query.split(":", 1)[1].strip()
        f = kw_by_acc.get(acc)
        moon = (s.get("moonlighting") or "").strip().upper() == "YES"
        rows.append(
            {
                "accession": acc,
                "gene": s["gene"],
                "lit_class": s["class"],
                "moonlighting": s["moonlighting"],
                "species": s["species"],
                "order": s["order"],
                "seq_sha256": f["seq_sha256"] if f else "",
                "length": f["length"] if f else "",
                "emb_row": f["emb_row"] if f else "",
                "emb_cterm_row": f["emb_cterm_row"] if f else "",
                "literature_positive": "yes" if f and not moon else "no",
            }
        )
    return rows, no_accession


def check_features(hashes, features_by_hash: dict) -> None:
    for h in sorted(hashes):
        f = features_by_hash.get(h)
        if f is None:
            raise evalio.StopError(f"hash {h} has no row in features_unique.tsv.gz")
        for col in ("sp_prob", "gpi_prob", "ser_thr_frac"):
            evalio.float_or_stop(f[col], f"features_unique.tsv.gz {h} {col}")
        if not f["sp_prediction"]:
            raise evalio.StopError(f"features_unique.tsv.gz {h}: empty sp_prediction")
        if f["gpi_call"] not in GPI_CALLS:
            raise evalio.StopError(f"features_unique.tsv.gz {h}: gpi_call {f['gpi_call']!r}")


def fasta_bytes(hashes, seq_by_hash: dict) -> bytes:
    buf = io.StringIO()
    for h in sorted(hashes):
        seq = seq_by_hash.get(h)
        if seq is None:
            raise evalio.StopError(f"hash {h} has no row in unique_sequences.tsv.gz")
        buf.write(f">{h}\n{seq}\n")
    raw = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(buf.getvalue().encode())
    return raw.getvalue()


def run(work: Path, species_path: Path, seeds_path: Path, clades_path: Path, arguments=()):
    work = Path(work)
    chain = check_chain(work)
    truth_rows = truth_table.read_tsv(work / "truth_set_triaged.tsv.gz")
    problems = iea_problems(truth_rows)
    if problems:
        raise evalio.StopError(f"{len(problems)} IEA problems, first: {problems[0]}")
    features = truth_table.read_tsv(work / "phaseb" / "features.tsv.gz")
    fu = {
        r["seq_sha256"]: r for r in truth_table.read_tsv(work / "phaseb" / "features_unique.tsv.gz")
    }
    kw_rows = truth_table.read_tsv(work / "keyword_tier.tsv.gz")
    clades = read_tc_clades(clades_path)
    check_tc_clades(kw_rows, clades)
    species_rows = truth_table.read_tsv(species_path)
    members, unmatched = go_members(truth_rows, features, species_rows)
    go_rows, go_log = dedupe.merge_go(members)
    go_hashes: dict[str, set] = {}
    for m in members:
        go_hashes.setdefault(m["seq_sha256"], set()).add(m["class"])
    go_detail = {h: ",".join(sorted(v)) for h, v in go_hashes.items()}
    kw_by_acc = kw_features(features)
    tc_rows, tc_log = dedupe.apply_precedence(go_detail, tc_members(kw_rows, kw_by_acc, clades))
    table = sorted([*go_rows.values(), *tc_rows.values()], key=lambda r: r["seq_sha256"])
    lit, no_accession = literature_rows(read_seeds(seeds_path), kw_by_acc)
    hashes = {r["seq_sha256"] for r in table} | {r["seq_sha256"] for r in lit if r["seq_sha256"]}
    check_features(hashes, fu)
    seqs = {
        r["seq_sha256"]: r["sequence"]
        for r in truth_table.read_tsv(work / "phaseb" / "unique_sequences.tsv.gz")
        if r["seq_sha256"] in hashes
    }
    fasta = fasta_bytes(hashes, seqs)
    log_rows = go_log + tc_log
    shared = Counter(
        (m["class"], m["label"], m["d8_class"])
        for m in members
        if m["seq_sha256"] in {r["seq_sha256"] for r in kw_rows}
    )
    log = {
        "all_sources": True,
        "truth_set_sha256": chain["truth_set_sha256"],
        "input_sha256": {
            "truth_set_triaged.tsv.gz": chain["truth_set_triaged_sha256"],
            "features.tsv.gz": manifest.sha256_file(work / "phaseb" / "features.tsv.gz"),
            "features_unique.tsv.gz": manifest.sha256_file(
                work / "phaseb" / "features_unique.tsv.gz"
            ),
            "unique_sequences.tsv.gz": chain["unique_sequences_sha256"],
            "keyword_tier.tsv.gz": manifest.sha256_file(work / "keyword_tier.tsv.gz"),
            "species.tsv": manifest.sha256_file(species_path),
            "eurotiomycetes_seeds.tsv": manifest.sha256_file(seeds_path),
            "tc_taxon_clades.tsv": manifest.sha256_file(clades_path),
        },
        "go_members_by_source_class": {
            s: dict(sorted(Counter(m["class"] for m in members if m["source_id"] == s).items()))
            for s in sorted({m["source_id"] for m in members})
        },
        "labelled_genes_without_sequence": dict(sorted(unmatched.items())),
        "table_rows_by_origin_class": {
            f"{o}:{c}": n
            for (o, c), n in sorted(Counter((r["origin"], r["class"]) for r in table).items())
        },
        "log_rows_by_reason": dict(sorted(Counter(r["reason"] for r in log_rows).items())),
        "tc_rows": len(kw_rows),
        "tc_rows_dropped_by_go_class": dict(sorted(Counter(r["detail"] for r in tc_log).items())),
        "go_members_sharing_a_tc_hash": {
            f"{c}:{lab}:{d8}": n for (c, lab, d8), n in sorted(shared.items())
        },
        "literature_rows": len(lit),
        "literature_positives": sum(r["literature_positive"] == "yes" for r in lit),
        "literature_no_accession": sorted(no_accession),
        "sequences": len(hashes),
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    out = evalio.out_dir(work)
    evalio.write_outputs(
        out,
        {
            "eval_table.tsv.gz": lambda p: truth_table.write_tsv(p, dedupe.TABLE_COLUMNS, table),
            "eval_literature.tsv": lambda p: truth_table.write_tsv(p, LITERATURE_COLUMNS, lit),
            "eval_dedupe_log.tsv": lambda p: truth_table.write_tsv(p, dedupe.LOG_COLUMNS, log_rows),
            "eval_sequences.fasta.gz": lambda p: Path(p).write_bytes(fasta),
        },
        "build_run.json",
        log,
    )
    return table, lit, log_rows, log


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--seeds", default=None, help=f"default: $PROJ_ROOT/{DEFAULT_SEEDS}")
    parser.add_argument("--tc-clades", default=str(evalio.PHASEC_DIR / "tc_taxon_clades.tsv"))
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    seeds = Path(args.seeds) if args.seeds else paths.repo_root() / DEFAULT_SEEDS
    try:
        table, lit, log_rows, log = run(
            work,
            Path(args.species),
            seeds,
            Path(args.tc_clades),
            list(argv) if argv is not None else sys.argv[1:],
        )
    except (evalio.StopError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(f"table_rows={len(table)} sequences={log['sequences']} literature={len(lit)}")
    for key, n in log["table_rows_by_origin_class"].items():
        print(f"  {key}={n}")
    for key, n in log["tc_rows_dropped_by_go_class"].items():
        print(f"  tc_dropped_go_class_{key}={n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
