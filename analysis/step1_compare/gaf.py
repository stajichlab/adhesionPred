"""Read GAF 2.x files and apply the D1 row filters (spec section 3.1, steps 2 to 4).

Standard library only. The filters copy /tmp/glyco_spec/d1_count.py in this order:
1. column 1 must be the most frequent database in the file, and column 12 must be
   `protein` or `gene_product`;
2. if a taxon is given, the first taxon in column 13 must equal it;
3. (gene, term, evidence, aspect) is recorded for the all-aspect IEA fraction;
4. aspect must be `C` and the qualifier must not contain `NOT`;
5. rows with an obsolete term are counted and dropped.
"""

import collections
import gzip
from dataclasses import dataclass, field
from pathlib import Path

HOMOLOGY_CODES = frozenset({"IBA", "IBD", "IKR", "IRD", "ISS", "ISO", "ISA", "ISM", "RCA"})
EXPERIMENTAL_CODES = frozenset(
    {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"}
)
PROTEIN_TYPES = frozenset({"protein", "gene_product"})
GZIP_MAGIC = b"\x1f\x8b"


class GafFormatError(ValueError):
    """A GAF file is empty, truncated or not a GAF file."""


@dataclass(frozen=True, slots=True)
class GafRow:
    db: str
    gene_id: str
    symbol: str
    qualifier: str
    term: str
    evidence: str
    aspect: str
    synonyms: str
    object_type: str
    taxon: str


def open_text(path: str | Path):
    """Open plain or gzip text. Gzip is detected by its magic bytes, not by the suffix."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == GZIP_MAGIC:
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def iter_gaf_rows(path: str | Path):
    with open_text(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            if raw.startswith("!") or not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) < 13:
                raise GafFormatError(f"{path}:{lineno}: {len(f)} columns, GAF needs at least 13")
            yield GafRow(f[0], f[1], f[2], f[3], f[4], f[6], f[8], f[10], f[11], f[12])


def header_value(path: str | Path, key: str) -> str:
    """Return the value of the first `!key:` header line, or '' if absent."""
    prefix = f"!{key}:"
    with open_text(path) as handle:
        for raw in handle:
            if not raw.startswith("!"):
                break
            if raw.startswith(prefix):
                return raw[len(prefix) :].strip()
    return ""


def taxon_matches(taxon_field: str, taxon_id: str) -> bool:
    """True if the first taxon in GAF column 13 has this NCBI id (`taxon:` or `NCBITaxon:`)."""
    first = taxon_field.split("|")[0]
    return first.split(":")[-1] == str(taxon_id)


@dataclass
class FilteredGaf:
    path: str
    primary_db: str
    cc_rows: list[GafRow]
    all_aspect_triples: set[tuple[str, str, str, str]]
    dropped: collections.Counter = field(default_factory=collections.Counter)
    unknown_term_rows: int = 0


def filter_gaf(path: str | Path, ontology, taxon_id: str | None = None) -> FilteredGaf:
    rows = list(iter_gaf_rows(path))
    if not rows:
        raise GafFormatError(f"{path}: no data rows")
    primary = collections.Counter(r.db for r in rows).most_common(1)[0][0]
    out = FilteredGaf(str(path), primary, [], set())
    for r in rows:
        if r.db != primary:
            out.dropped["other_db"] += 1
            continue
        if r.object_type not in PROTEIN_TYPES:
            out.dropped["object_type"] += 1
            continue
        if taxon_id and not taxon_matches(r.taxon, taxon_id):
            out.dropped["taxon"] += 1
            continue
        out.all_aspect_triples.add((r.gene_id, r.term, r.evidence, r.aspect))
        if r.aspect != "C":
            out.dropped["other_aspect"] += 1
            continue
        if "NOT" in r.qualifier:
            out.dropped["not_qualifier"] += 1
            continue
        if r.term in ontology.obsolete:
            out.dropped["obsolete_term"] += 1
            continue
        if not ontology.known(r.term):
            out.unknown_term_rows += 1
        out.cc_rows.append(r)
    return out
