"""Build the D1 truth table and its counts from filtered GAF rows. Standard library only."""

import csv
import gzip
import io
from dataclasses import dataclass, field
from pathlib import Path

import labels

TRUTH_COLUMNS = (
    "source_id",
    "species",
    "taxon_id",
    "in_clade",
    "role",
    "gene_id",
    "symbol",
    "synonym1",
    "label",
    "subset",
    "stratum",
    "tier",
    "label_no_homology",
    "label_experimental",
    "homology_only",
    "pm_candidate",
    "evidence_codes",
    "surface_evidence",
    "internal_evidence",
    "internal_evidence_htp_only",
    "secretory_evidence",
    "source_file",
    "source_sha256",
    "source_date",
    "obo_sha256",
)

COUNT_COLUMNS = (
    "source_id",
    "primary_db",
    "genes_cc",
    "genes_noniea_cc",
    "p_ext",
    "p_ext_wall",
    "p_ext_extonly",
    "n_int",
    "n_sec",
    "ambiguous",
    "pm_candidates",
    "cc_iea_triples",
    "cc_triples",
    "cc_iea_frac",
    "all_aspects_iea_frac",
    "obsolete_rows",
    "unknown_term_rows",
    "exp_p_ext",
    "exp_n_int",
    "exp_n_sec",
    "exp_ambiguous",
    "nohom_p_ext",
    "nohom_n_int",
    "nohom_n_sec",
    "nohom_ambiguous",
    "direct_p_ext",
    "direct_n_int",
    "direct_n_sec",
    "direct_ambiguous",
    "ambiguous_htp_only",
)


@dataclass
class SourceInfo:
    source_id: str
    species: str
    taxon_id: str
    in_clade: str
    role: str
    source_file: str
    source_sha256: str
    source_date: str
    obo_sha256: str


@dataclass
class GeneRecord:
    gene_id: str
    symbol: str
    synonym1: str
    rows: set[tuple[str, str]] = field(default_factory=set)  # (term, evidence)


def collect_genes(filtered) -> dict[str, GeneRecord]:
    genes: dict[str, GeneRecord] = {}
    for r in filtered.cc_rows:
        rec = genes.get(r.gene_id)
        if rec is None:
            rec = GeneRecord(r.gene_id, r.symbol, r.synonyms.split("|")[0])
            genes[r.gene_id] = rec
        rec.rows.add((r.term, r.evidence))
    return genes


def _terms(rec: GeneRecord, ontology, policy: str | None) -> frozenset[str]:
    out: set[str] = set()
    for term, evidence in rec.rows:
        if policy is None or labels.policy_accepts(policy, evidence):
            out |= ontology.ancestors(term)
    return frozenset(out)


def _codes(rec: GeneRecord, ontology, targets: frozenset[str]) -> str:
    codes = {ev for term, ev in rec.rows if ev != "IEA" and ontology.ancestors(term) & targets}
    return ",".join(sorted(codes))


def build_truth_rows(filtered, ontology, info: SourceInfo) -> list[dict[str, str]]:
    rows = []
    for gene_id, rec in sorted(collect_genes(filtered).items()):
        any_terms = _terms(rec, ontology, None)
        by_policy = {p: _terms(rec, ontology, p) for p in labels.POLICIES}
        label = labels.classify(by_policy["non_iea"], any_terms)
        label_nohom = labels.classify(by_policy["no_homology"], any_terms)
        label_exp = labels.classify(by_policy["experimental"], any_terms)
        subset = labels.subset_of(label, by_policy["non_iea"])
        internal = _codes(rec, ontology, labels.INTERNAL)
        rows.append(
            {
                "source_id": info.source_id,
                "species": info.species,
                "taxon_id": info.taxon_id,
                "in_clade": info.in_clade,
                "role": info.role,
                "gene_id": gene_id,
                "symbol": rec.symbol,
                "synonym1": rec.synonym1,
                "label": label,
                "subset": subset,
                "stratum": subset or label,
                "tier": "T-a",
                "label_no_homology": label_nohom,
                "label_experimental": label_exp,
                "homology_only": "yes" if label != label_nohom else "no",
                "pm_candidate": "yes"
                if labels.is_pm_candidate(label, by_policy["non_iea"])
                else "no",
                "evidence_codes": ",".join(sorted({ev for _, ev in rec.rows})),
                "surface_evidence": _codes(rec, ontology, labels.SURFACE),
                "internal_evidence": internal,
                "internal_evidence_htp_only": "yes"
                if labels.htp_only(set(internal.split(",")) - {""})
                else "no",
                "secretory_evidence": _codes(rec, ontology, labels.SECRETORY),
                "source_file": info.source_file,
                "source_sha256": info.source_sha256,
                "source_date": info.source_date,
                "obo_sha256": info.obo_sha256,
            }
        )
    return rows


def in_direct_stratum(row: dict[str, str], label: str) -> bool:
    """Direct-evidence truth (spec 3.3): the label is the same with and without homology codes.

    This is an intersection. Labels recomputed without homology codes are NOT used, because
    that moves genes whose only internal term is IBA or ISS into P-ext (contradicts Q9).
    """
    return row["label"] == label and row["homology_only"] == "no"


def _n(rows, column, value) -> int:
    return sum(1 for r in rows if r[column] == value)


def count_rows(rows: list[dict[str, str]], filtered, source_id: str) -> dict[str, str]:
    triples = {(r.gene_id, r.term, r.evidence) for r in filtered.cc_rows}
    iea = sum(1 for t in triples if t[2] == "IEA")
    iea_all = sum(1 for t in filtered.all_aspect_triples if t[2] == "IEA")
    counts = {
        "source_id": source_id,
        "primary_db": filtered.primary_db,
        "genes_cc": len(rows),
        "genes_noniea_cc": sum(1 for r in rows if r["evidence_codes"] != "IEA"),
        "p_ext": _n(rows, "label", labels.P_EXT),
        "p_ext_wall": _n(rows, "subset", "wall"),
        "p_ext_extonly": _n(rows, "subset", "extracellular-only"),
        "n_int": _n(rows, "label", labels.N_INT),
        "n_sec": _n(rows, "label", labels.N_SEC),
        "ambiguous": _n(rows, "label", labels.AMBIGUOUS),
        "pm_candidates": _n(rows, "pm_candidate", "yes"),
        "cc_iea_triples": iea,
        "cc_triples": len(triples),
        "cc_iea_frac": f"{iea / len(triples):.3f}",
        "all_aspects_iea_frac": f"{iea_all / len(filtered.all_aspect_triples):.3f}",
        "obsolete_rows": filtered.dropped["obsolete_term"],
        "unknown_term_rows": filtered.unknown_term_rows,
        "exp_p_ext": _n(rows, "label_experimental", labels.P_EXT),
        "exp_n_int": _n(rows, "label_experimental", labels.N_INT),
        "exp_n_sec": _n(rows, "label_experimental", labels.N_SEC),
        "exp_ambiguous": _n(rows, "label_experimental", labels.AMBIGUOUS),
        "nohom_p_ext": _n(rows, "label_no_homology", labels.P_EXT),
        "nohom_n_int": _n(rows, "label_no_homology", labels.N_INT),
        "nohom_n_sec": _n(rows, "label_no_homology", labels.N_SEC),
        "nohom_ambiguous": _n(rows, "label_no_homology", labels.AMBIGUOUS),
        "direct_p_ext": sum(in_direct_stratum(r, labels.P_EXT) for r in rows),
        "direct_n_int": sum(in_direct_stratum(r, labels.N_INT) for r in rows),
        "direct_n_sec": sum(in_direct_stratum(r, labels.N_SEC) for r in rows),
        "direct_ambiguous": sum(in_direct_stratum(r, labels.AMBIGUOUS) for r in rows),
        "ambiguous_htp_only": sum(
            r["label"] == labels.AMBIGUOUS and r["internal_evidence_htp_only"] == "yes"
            for r in rows
        ),
    }
    return {k: str(v) for k, v in counts.items()}


def write_tsv(path: str | Path, columns, rows) -> None:
    """Write a TSV. A path ending in .gz is gzip-compressed with a fixed mtime (byte-stable)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = buffer.getvalue().encode("utf-8")
    if path.suffix == ".gz":
        with (
            open(path, "wb") as raw,
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz,
        ):
            gz.write(data)
    else:
        path.write_bytes(data)


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
