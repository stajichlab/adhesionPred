"""Build the D1 truth table from filtered GAF rows. Standard library only."""

from dataclasses import dataclass, field

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
