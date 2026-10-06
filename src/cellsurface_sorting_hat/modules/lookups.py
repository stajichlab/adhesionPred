"""Lookups into precomputed tables: antigen ranking, Cys-rich tiers, spherule expression, TM helices.

These modules need an ID map. For *C. immitis* RS the RefSeq protein ID (``XP_...``) maps to a
gene ID (``CIMG_...``) through ``protein_map.tsv``. The antigen ranking and the spherule table are
keyed on that gene ID. A protein with no map entry or no table row is ``not_in_reference``. A protein
from a taxon where the table does not apply is ``not_applicable``.
"""

import csv
import gzip
import math
import re
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

ANTIGEN_COLUMNS = [
    "percentile",
    "antigenicity",
    "specificity",
    "prevalence",
    "max_crossreact",
    "rank",
    "idmap_method",
]
CYS_COLUMNS = ["tier", "cys_frac"]
EXPRESSION_COLUMNS = ["log2fc", "padj", "log2fc_8d"]
TM_COLUMNS = ["n_tm", "n_tm_mature", "topology"]
SIGNAL_PEPTIDE_WINDOW = (
    35  # a helix that starts at or before this residue may be the signal peptide
)


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def _check_number(path, line, column, text, mode):
    """Refuse a value that is not a finite number. ``mode`` is ``nonneg`` or ``finite``,
    optionally with the suffix ``_or_blank``."""
    base, _, blank = mode.partition("_or_")
    where = f"{path}:{line}: column {column!r}"
    if text.strip() == "":
        if blank:
            return
        raise ValueError(f"{where} is empty, a number is expected")
    try:
        value = float(text)
    except ValueError:
        raise ValueError(f"{where} value {text!r} is not a number") from None
    if not math.isfinite(value):
        raise ValueError(f"{where} value {text!r} is not finite")
    if base == "nonneg" and value < 0:
        raise ValueError(f"{where} value {text!r} is negative")


def _data_rows(path, delimiter, required, numeric):
    """Yield ``(line, row)``. Refuse a missing column, a short or long row and a bad number."""
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        fields = reader.fieldnames or []
        for col in required:
            if col not in fields:
                raise ValueError(f"{path}: missing column {col!r}")
        for r in reader:
            line = reader.line_num
            if None in r or any(v is None for v in r.values()):
                raise ValueError(f"{path}:{line}: expected {len(fields)} fields, row differs")
            for col, mode in numeric.items():
                _check_number(path, line, col, r[col], mode)
            yield line, r


def read_table(path, key, delimiter="\t", required=(), numeric=None):
    """Read a TSV into ``{row[key]: row}``. A repeated key keeps the first row and is counted.

    ``required`` lists more columns that must exist. ``numeric`` maps a column to ``nonneg`` or
    ``finite`` (add ``_or_blank`` to allow an empty value). Errors name the path and line.
    """
    numeric = numeric or {}
    rows, repeated = {}, 0
    for _, r in _data_rows(path, delimiter, (key, *required, *numeric), numeric):
        if r[key] in rows:
            repeated += 1
        else:
            rows[r[key]] = r
    return rows, repeated


def load_protein_map(path):
    rows, _ = read_table(path, "protein_id", required=("gene_id",))
    return {pid: r["gene_id"] for pid, r in rows.items()}


def gene_of_ranking_id(protein):
    """``CIMG_04613-t26_1-p1`` -> ``CIMG_04613``."""
    return protein.split("-", 1)[0]


RANKING_NUMBERS = {
    "percentile": "nonneg",
    "antigenicity": "nonneg",
    "specificity": "nonneg",
    "prevalence": "nonneg",
    "max_fungal_crossreact_pid": "nonneg",
}


def ranking_by_gene(ranking_tsv):
    """Best (lowest rank) ranking row per gene. Returns ``(rows, n_genes_with_several_rows)``.

    A tie in rank keeps the first row. Each number must be finite and not negative; ``rank`` must
    be a positive integer.
    """
    rows, several = {}, set()
    numeric = {"rank": "nonneg", **RANKING_NUMBERS}
    for line, r in _data_rows(ranking_tsv, "\t", ("protein", *numeric), numeric):
        if not r["rank"].strip().isdigit() or int(r["rank"]) < 1:
            raise ValueError(
                f"{ranking_tsv}:{line}: column 'rank' value {r['rank']!r} is not a positive integer"
            )
        gene = gene_of_ranking_id(r["protein"])
        cur = rows.get(gene)
        if cur is not None:
            several.add(gene)
        if cur is None or int(r["rank"]) < int(cur["rank"]):
            rows[gene] = r
    return rows, len(several)


def _check_taxa(proteins, taxa, module):
    """Every protein needs a taxon ID. A missing one is an input error, not a silent default."""
    for p in proteins:
        if p.id not in taxa:
            raise ValueError(f"{module}: no taxon ID for protein {p.id!r}")


def _applicable(taxa, pid, applicable_taxa):
    return taxa[pid] in applicable_taxa


def antigen_rows(proteins, taxa, protein_map, by_gene, applicable_taxa):
    _check_taxa(proteins, taxa, "antigen")
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        else:
            gene = protein_map.get(p.id) or (
                gene_of_ranking_id(p.id) if p.id.startswith("CIMG_") else None
            )
            r = by_gene.get(gene)
            if r is None:
                rows.append({"id": p.id, "state": "not_in_reference"})
            else:
                rows.append(
                    {
                        "id": p.id,
                        "state": "ok",
                        "percentile": r["percentile"],
                        "antigenicity": r["antigenicity"],
                        "specificity": r["specificity"],
                        "prevalence": r["prevalence"],
                        "max_crossreact": r["max_fungal_crossreact_pid"],
                        "rank": r["rank"],
                        "idmap_method": "gene_best_transcript",
                    }
                )
    return rows


def cys_rows(proteins, taxa, cys_table, applicable_taxa):
    """``cys_table``: ``{protein_id: row}`` from ``candidates.tsv.gz`` of ``analysis/cys_candidates``."""
    _check_taxa(proteins, taxa, "cys")
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        elif p.id not in cys_table:
            rows.append({"id": p.id, "state": "not_in_reference"})
        else:
            r = cys_table[p.id]
            rows.append({"id": p.id, "state": "ok", "tier": r["tier"], "cys_frac": r["cys_frac"]})
    return rows


def expression_rows(proteins, taxa, protein_map, spherule, applicable_taxa):
    """``spherule``: ``{gene_id: row}`` from ``spherule_surface_table.tsv.gz``."""
    _check_taxa(proteins, taxa, "expression")
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif not _applicable(taxa, p.id, applicable_taxa):
            rows.append({"id": p.id, "state": "not_applicable"})
        else:
            r = spherule.get(protein_map.get(p.id, ""))
            if r is None:
                rows.append({"id": p.id, "state": "not_in_reference"})
            else:
                rows.append(
                    {
                        "id": p.id,
                        "state": "ok",
                        "log2fc": r["log2fc_48h"],
                        "padj": r["padj_48h"],
                        "log2fc_8d": r["log2fc_8d"],
                    }
                )
    return rows


def tm_rows(proteins, tmhmm):
    """``tmhmm``: ``{protein_id: row}`` from the TMHMM table (columns ``pred_hel``, ``topology``)."""
    foreign = sorted(set(tmhmm) - {p.id for p in proteins})
    if foreign:
        raise ValueError(
            f"tm: {len(foreign)} ID(s) in the TMHMM table are not in the FASTA, for example {foreign[0]!r}"
        )
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in tmhmm:
            rows.append({"id": p.id, "state": "error"})
        else:
            r = tmhmm[p.id]
            n_tm = r["pred_hel"]
            if not n_tm.isdigit():
                raise ValueError(
                    f"tm: pred_hel {n_tm!r} for {p.id!r} is not a non-negative integer"
                )
            starts = [int(a) for a, _ in re.findall(r"(\d+)-(\d+)", r["topology"])]
            mature = sum(1 for a in starts if a > SIGNAL_PEPTIDE_WINDOW)
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "n_tm": n_tm,
                    "n_tm_mature": mature,
                    "topology": r["topology"],
                }
            )
    return rows
