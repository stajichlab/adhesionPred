"""D8: split P-ext genes with a non-IEA plasma-membrane term into P-gpi, PM-TM, pm-unresolved.

Standard library only. Evidence comes from UniProtKB REST JSON (fields accession, reviewed,
ft_lipid, ft_transmem, xref_sgd, xref_cgd, xref_pombase) and from curated_gpi.tsv (literature).

Rules (spec 2.2, Q2 and Q3):
- P-gpi: a literature row in curated_gpi.tsv, or a reviewed UniProt entry with a Lipidation
  feature whose description starts with "GPI-anchor" and whose evidence includes a code in
  CURATED_GPI_ECO. Predictor output never counts.
- PM-TM: not P-gpi, and at least one UniProt Transmembrane feature (any evidence code; the
  codes are recorded because they are often ECO:0000255, sequence analysis).
- pm-unresolved: neither. The spec does not define this case; it stays P-ext and is listed.
P-gpi is reported as a list, not a scored stratum, until curated_gpi.tsv has literature rows.
"""

import urllib.parse
from dataclasses import dataclass, field

CURATED_GPI_ECO = frozenset({"ECO:0000269"})  # experimental evidence used in manual assertion
UNIPROT_FIELDS = "accession,reviewed,ft_lipid,ft_transmem,xref_sgd,xref_cgd,xref_pombase"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
XREF_DB = {"sgd": "SGD", "cgd": "CGD", "pombase": "PomBase"}
P_GPI, PM_TM, PM_UNRESOLVED = "P-gpi", "PM-TM", "pm-unresolved"


class UniprotError(RuntimeError):
    """A UniProt response lacks expected fields, or no candidate gene matched any entry."""


@dataclass
class UniprotEvidence:
    accession: str
    reviewed: bool
    gpi_eco: list[str] = field(default_factory=list)
    gpi_feature_count: int = 0  # GPI-anchor Lipidation features (with or without evidence)
    gpi_features_without_eco: int = 0  # of those, features with no evidence code
    tm_count: int = 0
    tm_eco: list[str] = field(default_factory=list)
    xrefs: dict[str, set[str]] = field(default_factory=dict)

    @property
    def curated_gpi(self) -> bool:
        return self.reviewed and bool(set(self.gpi_eco) & CURATED_GPI_ECO)


def _eco(feature: dict) -> list[str]:
    return [e["evidenceCode"] for e in feature.get("evidences", []) if "evidenceCode" in e]


def parse_entry(entry: dict) -> UniprotEvidence:
    missing = [k for k in ("primaryAccession", "entryType") if k not in entry]
    if missing:
        raise UniprotError(f"UniProt entry lacks {missing}: {str(entry)[:80]}")
    ev = UniprotEvidence(
        accession=entry["primaryAccession"],
        reviewed=entry.get("entryType", "").startswith("UniProtKB reviewed"),
    )
    gpi, tm = set(), set()
    for feature in entry.get("features", []):
        if feature.get("type") == "Lipidation" and feature.get("description", "").startswith(
            "GPI-anchor"
        ):
            ev.gpi_feature_count += 1
            codes = _eco(feature)
            if not codes:
                ev.gpi_features_without_eco += 1
            gpi.update(codes)
        elif feature.get("type") == "Transmembrane":
            ev.tm_count += 1
            tm.update(_eco(feature))
    ev.gpi_eco, ev.tm_eco = sorted(gpi), sorted(tm)
    for xref in entry.get("uniProtKBCrossReferences", []):
        if "database" not in xref or "id" not in xref:
            raise UniprotError(f"{ev.accession}: cross-reference lacks database or id: {xref}")
        ev.xrefs.setdefault(xref["database"], set()).add(xref["id"])
    return ev


def query_terms(rows: list[dict[str, str]], mapping: str) -> dict[str, str]:
    """Return UniProt query term -> gene_id for each truth row."""
    if mapping == "uniprot":
        return {f"accession:{r['gene_id']}": r["gene_id"] for r in rows}
    return {f"xref:{mapping}-{r['gene_id']}": r["gene_id"] for r in rows}


def search_urls(terms: list[str], batch: int = 50) -> list[str]:
    urls = []
    for start in range(0, len(terms), batch):
        query = " OR ".join(f"({t})" for t in terms[start : start + batch])
        params = {"query": query, "fields": UNIPROT_FIELDS, "format": "json", "size": "500"}
        urls.append(f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}")
    return urls


def organism_gpi_url(taxon_id: str) -> str:
    query = f"(organism_id:{taxon_id}) AND (reviewed:true) AND (ft_lipid:GPI-anchor)"
    params = {"query": query, "fields": UNIPROT_FIELDS, "format": "json", "size": "500"}
    return f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"


def entries_for_gene(gene_id: str, mapping: str, entries: list[UniprotEvidence]):
    if mapping == "uniprot":
        return [e for e in entries if e.accession == gene_id]
    db = XREF_DB[mapping]
    return [e for e in entries if gene_id in e.xrefs.get(db, set())]


def classify_pm(entries: list[UniprotEvidence], literature: bool) -> tuple[str, str]:
    if literature:
        return P_GPI, "literature row in curated_gpi.tsv"
    for e in entries:
        if e.curated_gpi:
            return P_GPI, f"{e.accession} reviewed GPI-anchor {','.join(e.gpi_eco)}"
    with_tm = [e for e in entries if e.tm_count > 0]
    if with_tm:
        e = with_tm[0]
        return PM_TM, f"{e.accession} {e.tm_count} TM {','.join(e.tm_eco) or 'no ECO'}"
    if not entries:
        return PM_UNRESOLVED, "no UniProt entry found"
    return PM_UNRESOLVED, "no curated GPI evidence and no TM feature"
