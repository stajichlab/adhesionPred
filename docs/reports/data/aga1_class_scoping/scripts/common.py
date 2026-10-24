"""Shared helpers for the aga1-class scoping study (stdlib only). Reads committed tables; never writes sequences into the repo."""

import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
OUT = ROOT / "docs/reports/data/aga1_class_scoping"
CACHE = Path(
    os.environ.get("AGA1_SCOPING_CACHE", Path.home() / ".cache/adhesionPred/aga1_class_scoping")
)
RCT = ROOT / "analysis/calibration_truth/repeat_call_truth"
CTRL = ROOT / "data/controls"
AA = "ACDEFGHIKLMNPQRSTVWY"


def read_tsv(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
    head = rows[0]
    return [dict(zip(head, r + [""] * (len(head) - len(r)))) for r in rows[1:]]  # noqa: B905


def write_tsv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


def read_fasta(path):
    seqs, name = {}, None
    for ln in open(path):
        ln = ln.strip()
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(ln)
    return {k: "".join(v) for k, v in seqs.items()}


def clades_from_lineage_text(text):
    """(coarse, fine) clade from a UniProt 'Taxonomic lineage' string with '(rank)' suffixes.

    Same rule as analysis/wall_fusion_truth/map_kiss.py clades_from_lineage.
    """
    ranks = {}
    for part in text.split(", "):
        if part.endswith(")") and " (" in part:
            name, rank = part[:-1].rsplit(" (", 1)
            ranks.setdefault(rank, name)
    phylum, sub, cls = ranks.get("phylum", ""), ranks.get("subphylum", ""), ranks.get("class", "")
    if phylum == "Ascomycota":
        return sub, (cls if sub == "Pezizomycotina" else sub)
    if phylum == "Basidiomycota":
        return phylum, cls
    return phylum, phylum


def family_from_lineage_text(text):
    """The family-rank taxon name of a UniProt lineage string ('' if none)."""
    for part in text.split(", "):
        if part.endswith("(family)"):
            return part[: -len(" (family)")]
    return ""


def truth_tables():
    """{'Scer_S288C': rows, 'Calb_SC5314': rows} of the committed truth_v3 annotated tables."""
    return {p: read_tsv(RCT / f"truth_v3.{p}.annotated.tsv") for p in ("Scer_S288C", "Calb_SC5314")}


def all_named_accessions():
    """Accessions named in committed tables that this study needs a sequence for."""
    acc = set()
    for rows in truth_tables().values():
        acc |= {r["accession"] for r in rows if r["accession"]}
    for r in read_tsv(CTRL / "repeat-mechanism/adjudication/curation_table.consensus.tsv"):
        acc.add(r["accession"])
    for r in read_tsv(CTRL / "repeat-mechanism/curation_table.tsv"):
        acc.add(r["accession"])
    for r in read_tsv(CTRL / "hard-negatives-shared/controls.tsv"):
        acc.add(r["accession"])
    for f in ("curation_table.other_species.tsv", "hard_negatives_other_species.tsv"):
        acc |= {r["accession"] for r in read_tsv(CTRL / "repeat-adhesins-other-species" / f)}
    return acc


def uniprot_cache():
    """{accession: row} from every cached UniProt TSV batch (cache lives outside the repo)."""
    out = {}
    for f in sorted(CACHE.glob("uniprot_batch_*.tsv")):
        for r in read_tsv(f):
            out[r["Entry"]] = r
    return out


def local_sequences():
    """{accession: sequence} from data/controls/*/sequences.faa (committed; used in memory only)."""
    seqs = {}
    for f in sorted(CTRL.glob("*/sequences.faa")):
        for k, v in read_fasta(f).items():
            seqs.setdefault(k, v)
    return seqs


def all_sequences():
    """Local sequences first, then the UniProt cache (outside the repo). Returns ({acc: seq}, {acc: source})."""
    seqs, src = {}, {}
    for k, v in local_sequences().items():
        seqs[k], src[k] = v, "local_controls_faa"
    for k, r in uniprot_cache().items():
        if k not in seqs and r.get("Sequence"):
            seqs[k], src[k] = r["Sequence"], "uniprot_rest_cache"
    return seqs, src
