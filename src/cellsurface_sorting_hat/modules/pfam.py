"""Pfam family table and ``hmmsearch --domtblout`` results -> modules ``pfam_adhesion``, ``pfam_allergen``."""

import csv
from dataclasses import dataclass
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

FAMILY_COLUMNS = [
    "pfam_acc",
    "name",
    "class",
    "module",
    "source_pmid",
    "pfam_release",
    "specificity_note",
    "second_condition",
    "active",
    "active_by",
    "active_date",
]
MODULES = ("pfam_adhesion", "pfam_allergen")
SECOND_CONDITIONS = ("", "signal_peptide", "no_tm")
COLUMNS = ["hit", "families"]


class FamilyTableError(ValueError):
    """The family table is not valid."""


class NoActiveFamilyError(ValueError):
    """No family of the module has passed its specificity test and sign-off."""


@dataclass(frozen=True)
class Family:
    pfam_acc: str
    name: str
    cls: str
    module: str
    second_condition: str
    active: bool


def load_family_table(path):
    rows, seen = [], set()
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(FAMILY_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise FamilyTableError(f"{path}: missing column(s) {sorted(missing)}")
        for n, r in enumerate(reader, 2):
            where = f"{path}:{n}"
            acc = r["pfam_acc"].strip()
            if not acc.startswith("PF") or "." in acc:
                raise FamilyTableError(f"{where}: pfam_acc must look like PF05730 (no version)")
            if acc in seen:
                raise FamilyTableError(f"{where}: duplicate {acc}")
            seen.add(acc)
            if r["module"] not in MODULES:
                raise FamilyTableError(f"{where}: module must be one of {MODULES}")
            if r["second_condition"] not in SECOND_CONDITIONS:
                raise FamilyTableError(
                    f"{where}: second_condition must be one of {SECOND_CONDITIONS}"
                )
            if r["active"] not in ("yes", "no"):
                raise FamilyTableError(f"{where}: active must be yes or no")
            if r["active"] == "yes" and not (r["active_by"].strip() and r["active_date"].strip()):
                raise FamilyTableError(f"{where}: an active family needs active_by and active_date")
            rows.append(
                Family(
                    acc,
                    r["name"],
                    r["class"],
                    r["module"],
                    r["second_condition"],
                    r["active"] == "yes",
                )
            )
    return rows


def parse_domtblout(path):
    """Return a list of ``{"target", "acc", "ievalue", "score"}`` (accession without version).

    hmmsearch writes ``# [ok]`` as the last line of a finished run and ``--cut_ga`` in its option
    settings. A file without the trailer is truncated; a file without ``--cut_ga`` used another
    cutoff. Both are refused.
    """
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    if not any(line.strip() == "# [ok]" for line in lines[-3:]):
        raise ValueError(f"{path}: no '# [ok]' trailer; the hmmsearch run is not finished")
    if not any("--cut_ga" in line for line in lines if line.startswith("#")):
        raise ValueError(f"{path}: hmmsearch was not run with --cut_ga")
    hits = []
    for n, line in enumerate(lines, 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 22:
            raise ValueError(f"{path}:{n}: expected at least 22 fields, found {len(f)}")
        hits.append(
            {
                "target": f[0],
                "acc": f[4].split(".")[0],
                "ievalue": float(f[12]),
                "score": float(f[13]),
            }
        )
    return hits


def check_hit_ids(hits, fasta_ids):
    """A domain table whose targets are not FASTA IDs belongs to another proteome."""
    unknown = sorted({h["target"] for h in hits} - set(fasta_ids))
    if unknown:
        raise ValueError(
            f"{len(unknown)} domain-table target(s) are not in the FASTA, for example {unknown[0]!r}"
        )


def pfam_rows(proteins, hits, families, module, sp_calls=None, tm_counts=None):
    """Rows for ``module``. A protein is a hit if it has a domain of an active family of this module.

    ``second_condition = signal_peptide`` counts a hit only when ``sp_calls[id]`` is ``called``.
    ``second_condition = no_tm`` counts a hit only when ``tm_counts[id]`` is 0 (for example a CFEM
    domain in a receptor with transmembrane helices is not a cell wall CFEM protein). If the needed
    table is None, such a family never counts.
    """
    wanted = {f.pfam_acc: f for f in families if f.module == module and f.active}
    if not wanted:
        raise NoActiveFamilyError(f"{module}: no family of this module is active")
    found = {}
    for h in hits:
        fam = wanted.get(h["acc"])
        if fam is None:
            continue
        if (
            fam.second_condition == "signal_peptide"
            and (sp_calls or {}).get(h["target"]) != "called"
        ):
            continue
        if fam.second_condition == "no_tm" and (tm_counts or {}).get(h["target"]) != 0:
            continue
        found.setdefault(h["target"], set()).add(fam.pfam_acc)
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        accs = sorted(found.get(p.id, ()))
        rows.append(
            {"id": p.id, "state": "ok", "hit": "1" if accs else "0", "families": ",".join(accs)}
        )
    return rows


def specificity_report(hit_ids, member_ids, universe_ids):
    """Compare the proteins hit by a family with the known members (a specificity test).

    ``universe_ids`` are all proteins searched. Returns counts and the non-member hits, which a
    person must read before the family is made active.
    """
    hit, members, universe = set(hit_ids), set(member_ids), set(universe_ids)
    if not members <= universe:
        raise ValueError("member_ids must be a subset of universe_ids")
    tp, fp, fn = hit & members, hit - members, members - hit
    tn = universe - hit - members
    return {
        "tp": len(tp),
        "fp": len(fp),
        "fn": len(fn),
        "tn": len(tn),
        "sensitivity": len(tp) / len(members) if members else None,
        "specificity": len(tn) / (len(tn) + len(fp)) if (len(tn) + len(fp)) else None,
        "nonmember_hits": sorted(fp),
        "missed_members": sorted(fn),
    }
