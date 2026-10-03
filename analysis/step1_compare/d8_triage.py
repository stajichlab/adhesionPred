"""D8: split P-ext genes with a non-IEA plasma-membrane term into P-gpi, PM-TM, pm-unresolved.

Standard library only. Evidence comes from UniProtKB REST JSON (fields accession, reviewed,
ft_lipid, ft_transmem, xref_sgd, xref_cgd, xref_pombase) and from curated_gpi.tsv (literature).

Rules (spec 2.2, Q2 and Q3):
- P-gpi: a literature row in curated_gpi.tsv, or a reviewed UniProt entry with a Lipidation
  feature whose description starts with "GPI-anchor" and whose evidence includes a code in
  CURATED_GPI_ECO. Predictor output never counts. A UniProt TM feature blocks a literature row
  unless the row has override_tm=yes or a reviewed UniProt entry has the GPI evidence
  (owner decision 8, 2026-10-02).
- PM-TM: not P-gpi, and at least one UniProt Transmembrane feature (any evidence code; the
  codes are recorded because they are often ECO:0000255, sequence analysis). This includes a
  literature row that a TM feature blocks.
- pm-unresolved: neither. The spec does not define this case; it stays P-ext and is listed.
P-gpi is reported as a list, not a scored stratum, until a test set has 20 direct P-gpi
positives (owner decision 12, 2026-10-02).
"""

import csv
import datetime
import re
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


def classify_pm(
    entries: list[UniprotEvidence], literature: bool, override_tm: bool = False
) -> tuple[str, str]:
    """Class and reason for one P-ext gene with a plasma-membrane term.

    A literature row gives P-gpi, except when a UniProt TM feature blocks it: the gene has a TM
    feature, `override_tm` is false, and no reviewed UniProt entry has experimental GPI evidence.
    A blocked gene stays PM-TM (owner decision 8, 2026-10-02). The reason says so. When a literature
    row and a reviewed UniProt entry both give P-gpi, the reason names both."""
    with_tm = [e for e in entries if e.tm_count > 0]
    curated = [e for e in entries if e.curated_gpi]
    blocked = literature and bool(with_tm) and not override_tm and not curated
    if literature and not blocked:
        reason = "literature row in curated_gpi.tsv"
        if override_tm:
            reason += "; override_tm=yes"
            if not with_tm:
                reason += " not needed (no TM feature)"
        if curated:
            reason += f"; {curated[0].accession} reviewed GPI-anchor {','.join(curated[0].gpi_eco)}"
        return P_GPI, reason
    if curated:
        e = curated[0]
        return P_GPI, f"{e.accession} reviewed GPI-anchor {','.join(e.gpi_eco)}"
    if with_tm:
        e = with_tm[0]
        reason = f"{e.accession} {e.tm_count} TM {','.join(e.tm_eco) or 'no ECO'}"
        if blocked:
            reason += "; literature row blocked by the TM feature (override_tm=no)"
        return PM_TM, reason
    if not entries:
        return PM_UNRESOLVED, "no UniProt entry found"
    return PM_UNRESOLVED, "no curated GPI evidence and no TM feature"


CURATED_GPI_COLUMNS = (
    "source_id",
    "gene_id",
    "symbol",
    "pmid",
    "note",
    "species",
    "uniprot_accession",
    "evidence_level",
    "evidence_note",
    "reviewer",
    "review_date",
    "override_tm",
)
CURATED_GPI_REQUIRED = ("pmid", "evidence_note", "reviewer", "review_date")
EVIDENCE_LEVELS = frozenset({"direct", "transfer"})
UNMATCHED_COLUMNS = ("source_id", "gene_id", "symbol", "reason", "label")
CONFLICT_COLUMNS = (
    "source_id",
    "gene_id",
    "symbol",
    "uniprot_accessions",
    "tm_count",
    "tm_eco",
    "d8_reason",
    "override_tm",
)


class CuratedGpiError(ValueError):
    """curated_gpi.tsv has a missing column, a bad value or a repeated gene."""


def read_curated_gpi(path) -> list[dict[str, str]]:
    """Read and check curated_gpi.tsv. A header-only file gives an empty list.

    The file is plain tab-separated text: a double quote is an ordinary character, so a quoted
    sentence in `evidence_note` stays as written. A UTF-8 byte order mark is accepted. Messages
    name the physical line of the file. The PMID is checked for its form only. The reviewer opens
    each PMID and checks that it resolves (spec section 4)."""
    rows, seen = [], set()
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        missing = [c for c in CURATED_GPI_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise CuratedGpiError(f"{path}: missing columns {missing}")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise CuratedGpiError(f"{path}: duplicate column names")
        for row in reader:
            where = f"{path} line {reader.line_num}"
            if None in row or None in row.values():
                raise CuratedGpiError(
                    f"{where}: wrong number of fields (the header has {len(reader.fieldnames)})"
                )
            key = (row["source_id"], row["gene_id"])
            if not all(key):
                raise CuratedGpiError(f"{where}: source_id and gene_id are required")
            if key in seen:
                raise CuratedGpiError(f"{where}: {key[0]} {key[1]} appears twice")
            seen.add(key)
            if row["override_tm"] not in ("yes", "no"):
                raise CuratedGpiError(
                    f"{where}: override_tm must be yes or no, not {row['override_tm']!r}"
                )
            if row["evidence_level"] not in EVIDENCE_LEVELS:
                raise CuratedGpiError(f"{where}: evidence_level must be direct or transfer")
            for column in CURATED_GPI_REQUIRED:
                if not row[column]:
                    raise CuratedGpiError(f"{where}: {column} is required")
            if not re.fullmatch(r"[0-9]+(;[0-9]+)*", row["pmid"]):
                raise CuratedGpiError(f"{where}: pmid must be digits, joined with ';'")
            if not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}", row["review_date"]):
                raise CuratedGpiError(f"{where}: review_date must be YYYY-MM-DD")
            try:
                datetime.date.fromisoformat(row["review_date"])
            except ValueError as exc:
                raise CuratedGpiError(f"{where}: review_date is not a date: {exc}") from exc
            rows.append(row)
    return rows


def check_curated_gpi(curated_rows, truth_rows) -> list[dict[str, str]]:
    """Rows of curated_gpi.tsv that cannot change a D8 class, with the reason.

    D8 reads a row only for a P-ext gene that is a plasma-membrane candidate."""
    truth = {(r["source_id"], r["gene_id"]): r for r in truth_rows}
    problems = []
    for row in curated_rows:
        gene = truth.get((row["source_id"], row["gene_id"]))
        if gene is None:
            reason, label = "no_truth_gene", ""
        elif gene["label"] != "P-ext":
            reason, label = "outside_p_ext", gene["label"]
        elif gene["pm_candidate"] != "yes":
            reason, label = "not_pm_candidate", gene["label"]
        else:
            continue
        problems.append(
            {
                "source_id": row["source_id"],
                "gene_id": row["gene_id"],
                "symbol": row["symbol"],
                "reason": reason,
                "label": label,
            }
        )
    return problems


def tm_conflicts(triage_rows, literature) -> list[dict[str, str]]:
    """Literature genes that a TM feature blocked. They stay PM-TM until the owner reviews them.

    `literature` maps (source_id, gene_id) to the row's override_tm flag. A PM-TM gene with a
    literature row is always a blocked gene, because override_tm=yes would have given P-gpi."""
    return [
        {
            "source_id": t["source_id"],
            "gene_id": t["gene_id"],
            "symbol": t["symbol"],
            "uniprot_accessions": t["uniprot_accessions"],
            "tm_count": t["tm_count"],
            "tm_eco": t["tm_eco"],
            "d8_reason": t["d8_reason"],
            "override_tm": "no",
        }
        for t in triage_rows
        if t["d8_class"] == PM_TM and (t["source_id"], t["gene_id"]) in literature
    ]
