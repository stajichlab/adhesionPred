"""WHO/IUIS fungal allergen sequences and BLASTP results -> module ``allergen_homology``.

Module fields: ``identity`` (percent identity of the best local alignment), ``aligned_length``
(alignment length, gaps included), ``coverage`` (alignment length as a percent of the allergen
sequence, capped at 100), ``allergen_name``, ``allergen_species``, ``exposure`` and ``evidence``
(the IUIS evidence text of the matched allergen, or empty). A protein with no hit has ``0``, ``0``,
``0`` and empty text (state ``ok``). The engine makes two calls from these fields:
``iuis_allergen_similarity`` (identity >= 35% and aligned length >= 80, ONE local alignment; this is
not the sliding 80-residue window of the Codex rule) and ``iuis_allergen_homolog`` (identity >= 70%
and coverage >= 80%).
"""

import csv
import math
import re
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

BLAST_FIELDS = "qseqid sseqid pident length qlen slen bitscore evalue"
COLUMNS = [
    "identity",
    "aligned_length",
    "coverage",
    "allergen_name",
    "allergen_species",
    "exposure",
    "evidence",
]
META_COLUMNS = ["id", "name", "species", "tax_order", "exposure", "evidence", "length"]
RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWYXBZUJO")
MIN_USABLE_LENGTH = (
    30  # shorter entries are peptide fragments; they cannot reach an 80 aa alignment
)
_ID_BAD = re.compile(r"[^A-Za-z0-9_.\-]+")
_N_FIELDS = len(BLAST_FIELDS.split())


def _clean(value):
    """One table cell as plain text: whitespace collapsed, a pandas ``nan`` cell empty."""
    text = " ".join(str(value or "").split())
    return "" if text.lower() == "nan" else text


def _require_columns(reader, path, needed):
    missing = [c for c in needed if c not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"{path}: missing column(s) {missing}")


def _number(text, where, name, cast=float, low=None, high=None, positive=False):
    """Parse one BLAST number. Refuse text, nan, inf, a negative value or one outside the range."""
    try:
        value = cast(text)
    except ValueError:
        raise ValueError(f"{where}: {name} {text!r} is not a number") from None
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{where}: {name} {text!r} is not finite")
    if value < 0 or (positive and value == 0):
        raise ValueError(
            f"{where}: {name} {text!r} must be {'positive' if positive else 'not negative'}"
        )
    if high is not None and value > high:
        raise ValueError(f"{where}: {name} {text!r} is above {high}")
    return value


def _parse_line(path, n, line):
    f = line.split("\t")
    where = f"{path}:{n}"
    if len(f) != _N_FIELDS:
        raise ValueError(f"{where}: expected {_N_FIELDS} fields, found {len(f)}")
    return f[0], {
        "subject": f[1],
        "identity": _number(f[2], where, "identity", high=100.0),
        "length": _number(f[3], where, "alignment length", int, positive=True),
        "qlen": _number(f[4], where, "query length", int, positive=True),
        "slen": _number(f[5], where, "subject length", int, positive=True),
        "bitscore": _number(f[6], where, "bit score"),
        "evalue": _number(f[7], where, "e-value"),
    }


def _blast_lines(path):
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if line.strip():
            yield n, line


def build_allergen_fasta(isoallergen_tsv, out_fasta, allergen_tsv=None):
    """Write the IUIS fungal sequences as FASTA and ``<out_fasta>.meta.tsv``.

    ID = ``<IsoAllergenName>|<AllergenID>``. All whitespace is removed from a sequence. An entry
    whose sequence has characters that are not residues (free text such as "N-terminal peptide: ...")
    or is shorter than 30 residues is not written and is listed in the returned ``skipped``. The
    meta table carries species, order, exposure and the IUIS evidence text of each allergen (from
    ``allergen_tsv``, the ``iuis_fungal_allergens.tsv`` table, when given).

    Returns ``(n_written, skipped)`` where ``skipped`` is a list of ``(id, reason)``.
    """
    allergen = {}
    if allergen_tsv:
        with open(allergen_tsv, encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            _require_columns(reader, allergen_tsv, ["AllergenID"])
            for r in reader:
                allergen[str(r["AllergenID"])] = r
    n, skipped, seen = 0, [], set()
    meta_path = Path(str(out_fasta) + ".meta.tsv")
    with (
        open(isoallergen_tsv, encoding="utf-8-sig", newline="") as fh,
        open(out_fasta, "w") as out,
        open(meta_path, "w") as meta,
    ):
        meta.write("\t".join(META_COLUMNS) + "\n")
        reader = csv.DictReader(fh, delimiter="\t")
        _require_columns(reader, isoallergen_tsv, ["AllergenID", "Sequence"])
        for r in reader:
            raw = (r.get("Sequence") or "").strip()
            name = _ID_BAD.sub("_", (r.get("IsoName") or r.get("Name") or "").strip())
            pid = f"{name}|{r['AllergenID']}"
            seq = "".join(raw.split()).upper()
            if not seq or seq == "NAN":
                continue
            if not set(seq) <= RESIDUES:
                skipped.append((pid, "not a protein sequence"))
                continue
            if len(seq) < MIN_USABLE_LENGTH:
                skipped.append((pid, f"peptide fragment of {len(seq)} residues"))
                continue
            if pid in seen:
                raise ValueError(f"duplicate allergen ID {pid}")
            seen.add(pid)
            out.write(f">{pid}\n{seq}\n")
            a = allergen.get(str(r["AllergenID"]), {})
            evidence = _clean(a.get("Allergenicity"))[:300]
            meta.write(
                "\t".join(
                    [
                        pid,
                        _clean(a.get("Name")),
                        _clean(a.get("Species")),
                        _clean(a.get("TaxOrder")),
                        _clean(a.get("Exposure")),
                        evidence,
                        str(len(seq)),
                    ]
                )
                + "\n"
            )
            n += 1
    return n, skipped


def read_meta(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def parse_blast(path):
    """Return the best local alignment per query: ``{query: {...}}`` (highest bit score)."""
    best = {}
    for n, line in _blast_lines(path):
        q, hit = _parse_line(path, n, line)
        cur = best.get(q)
        if cur is None or (hit["bitscore"], hit["identity"]) > (cur["bitscore"], cur["identity"]):
            best[q] = hit
    return best


def allergen_name(subject):
    """``Asp_f_1.0101|11`` -> ``Asp_f_1.0101``."""
    return subject.split("|", 1)[0]


def species_code(subject):
    """``Asp_f_1.0101|11`` -> ``Asp_f`` (genus and species letters of the allergen name)."""
    parts = allergen_name(subject).split("_")
    return "_".join(parts[:2]) if len(parts) >= 3 else parts[0]


def allergen_rows(proteins, best, meta=None):
    check_blast_ids(best, [p.id for p in proteins])
    meta = meta or {}
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        h = best.get(p.id)
        if h is None:
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "identity": "0",
                    "aligned_length": "0",
                    "coverage": "0",
                    "allergen_name": "",
                    "allergen_species": "",
                    "exposure": "",
                    "evidence": "",
                }
            )
            continue
        m = meta.get(h["subject"], {})
        cov = min(100.0, 100.0 * h["length"] / h["slen"]) if h["slen"] else 0.0
        rows.append(
            {
                "id": p.id,
                "state": "ok",
                "identity": f"{h['identity']:.2f}",
                "aligned_length": str(h["length"]),
                "coverage": f"{cov:.1f}",
                "allergen_name": allergen_name(h["subject"]),
                "allergen_species": m.get("species", ""),
                "exposure": m.get("exposure", ""),
                "evidence": m.get("evidence", ""),
            }
        )
    return rows


def check_blast_ids(best, fasta_ids):
    """A BLAST table whose queries are not FASTA IDs belongs to another proteome."""
    unknown = sorted(set(best) - set(fasta_ids))
    if unknown:
        raise ValueError(
            f"{len(unknown)} BLAST query ID(s) are not in the FASTA, for example {unknown[0]!r}"
        )


def parse_blast_hits(path):
    """All rows of a BLAST table as dicts."""
    hits = []
    for n, line in _blast_lines(path):
        q, hit = _parse_line(path, n, line)
        hits.append({"query": q, **hit})
    return hits


def best_hit_other_species(hits, ids):
    """Per query, the best hit (bit score) to an allergen of another species (see ``species_code``)."""
    best = {}
    species = {i: species_code(i) for i in ids}
    for h in hits:
        q, s = h["query"], h["subject"]
        if q == s or q not in species or s not in species or species[q] == species[s]:
            continue
        cur = best.get(q)
        if cur is None or (h["bitscore"], h["identity"]) > (cur["bitscore"], cur["identity"]):
            best[q] = h
    return best


def recall_by_rule(ids, best_other, rules):
    """Share of allergens whose best hit in ANOTHER species meets each rule.

    ``rules`` is a list of ``(name, min_identity, min_aligned_length, min_coverage)``. These are the
    rules the engine applies, so each number belongs to one call. The denominator is ``ids`` (all
    sequences of the FASTA); a sequence with no hit counts as not recovered. This is sensitivity only.
    """
    out = []
    for name, min_id, min_len, min_cov in rules:
        k = 0
        for i in ids:
            h = best_other.get(i)
            if (
                h
                and h["identity"] >= min_id
                and h["length"] >= min_len
                and 100.0 * h["length"] / h["slen"] >= min_cov
            ):
                k += 1
        out.append({"rule": name, "recovered": k, "n": len(ids)})
    return out


def lso_report(blast_path, fasta_ids, rules=None):
    """Leave-species-out recall of the allergen set against itself (all-against-all BLAST table).

    Every sequence of ``fasta_ids`` must appear as a query (a sequence with no self-hit means that
    BLAST masked or dropped it); otherwise the function refuses.
    """
    rules = rules or [
        ("iuis_allergen_similarity", 35.0, 80, 0.0),
        ("iuis_allergen_homolog", 70.0, 0, 80.0),
    ]
    hits = parse_blast_hits(blast_path)
    check_blast_ids({h["query"]: h for h in hits}, fasta_ids)
    check_blast_ids({h["subject"]: h for h in hits}, fasta_ids)
    queries = {h["query"] for h in hits}
    missing = sorted(set(fasta_ids) - queries)
    if missing:
        raise ValueError(
            f"{len(missing)} sequence(s) have no BLAST line, for example {missing[0]!r}"
        )
    best = best_hit_other_species(hits, fasta_ids)
    species = {species_code(i) for i in fasta_ids}
    return {
        "n_sequences": len(fasta_ids),
        "n_species": len(species),
        "recall": recall_by_rule(sorted(fasta_ids), best, rules),
    }
