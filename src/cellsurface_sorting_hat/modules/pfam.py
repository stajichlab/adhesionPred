"""Pfam family table and ``hmmsearch --domtblout`` results -> modules ``pfam_adhesion``, ``pfam_allergen``."""

import csv
from dataclasses import dataclass
from pathlib import Path

from cellsurface_sorting_hat.family_tier import SMOKE, TIERS, UNTESTED
from cellsurface_sorting_hat.modules.base import invalid_row

# Columns of data/sorting_hat/family_table.tsv.
# `active` means the owner enabled the family. It does not mean tested.
# The family tier describes the family model's own evidence; the call status describes the call's
# measurement; neither implies the other.
# A `narrow_ok=yes` needs a tier of smoke or better and a sign-off (who, when) with an evidence
# summary: the family's own members, how many the model hits, the leakage of those members and the
# hits on near-neighbour negatives. A count on a shared panel (coverage) is not enough on its own.
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
    "specificity_tier",  # derived from family_measure.tsv by family_tier.py; never edited by hand
    "narrow_ok",  # yes: the narrow call reads this family (owner sign-off); no: only the broad call
    "narrow_ok_by",
    "narrow_ok_date",
    "narrow_ok_evidence",  # members, hits, leakage and near-neighbour negative hits behind a yes
]
# Text that marks a ``narrow_ok_evidence`` as a pre-filled draft. A row with ``narrow_ok=yes`` must
# not carry it: the owner replaces the draft with the real evidence when signing.
DRAFT_EVIDENCE_MARKERS = ("draft", "not signed")
MODULES = ("pfam_adhesion", "pfam_allergen", "pfam_hydrophobin", "pfam_hsba")
SECOND_CONDITIONS = ("", "signal_peptide", "no_tm")
COLUMNS = ["hit", "families"]
NARROW_MODULES = ("pfam_adhesion",)  # the narrow call reads only this module


def columns_for(module):
    """Columns of a Pfam module table after ``id`` and ``state``. Only ``pfam_adhesion`` has
    ``hit_narrow``: no call reads a narrow flag of the other modules."""
    return ["hit", "hit_narrow", "families"] if module in NARROW_MODULES else list(COLUMNS)


class FamilyTableError(ValueError):
    """The family table is not valid."""


class NoActiveFamilyError(ValueError):
    """No family of the module is active (enabled by the owner; a row's specificity_note says whether a test was run)."""


@dataclass(frozen=True)
class Family:
    pfam_acc: str
    name: str
    cls: str
    module: str
    second_condition: str
    active: bool
    tier: str = UNTESTED  # evidence tier of the family model (family_tier.py)
    narrow_ok: bool = False  # signed by the owner: the narrow call reads this family


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
            if r["specificity_tier"] not in TIERS:
                raise FamilyTableError(f"{where}: specificity_tier must be one of {TIERS}")
            if r["narrow_ok"] not in ("yes", "no"):
                raise FamilyTableError(f"{where}: narrow_ok must be yes or no")
            if r["narrow_ok"] == "yes":
                if r["specificity_tier"] == UNTESTED:
                    raise FamilyTableError(
                        f"{where}: narrow_ok=yes needs a specificity_tier of {SMOKE} or better, "
                        f"not {UNTESTED}"
                    )
                for column in ("narrow_ok_by", "narrow_ok_date", "narrow_ok_evidence"):
                    if not r[column].strip():
                        raise FamilyTableError(f"{where}: narrow_ok=yes needs {column}")
                evidence = r["narrow_ok_evidence"].lower()
                if any(marker in evidence for marker in DRAFT_EVIDENCE_MARKERS):
                    raise FamilyTableError(
                        f"{where}: narrow_ok=yes but narrow_ok_evidence is still a draft "
                        f"(it contains one of {DRAFT_EVIDENCE_MARKERS}); replace it with the evidence"
                    )
            if r.get("narrow_ok_draft", "") not in ("", "yes", "no"):
                raise FamilyTableError(f"{where}: narrow_ok_draft must be yes, no or empty")
            rows.append(
                Family(
                    acc,
                    r["name"],
                    r["class"],
                    r["module"],
                    r["second_condition"],
                    r["active"] == "yes",
                    r["specificity_tier"],
                    r["narrow_ok"] == "yes",
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
                "seq_score": float(f[7]),
                "query": f[3],
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
    domain in a receptor with transmembrane helices is not a cell wall CFEM protein).
    ``sp_calls`` and ``tm_counts`` hold only proteins whose row in the condition module is ``ok``.
    A protein with a domain hit of such a family and no usable condition value is ``error``: its
    result is unknown, not 0. A protein with no domain hit is ``ok`` with hit 0.

    ``pfam_adhesion`` rows also have ``hit_narrow``: a hit in a family with tier ``smoke`` or better
    and ``narrow_ok`` (owner sign-off). It implies ``hit``: every narrow family is an active family.
    """
    wanted = {f.pfam_acc: f for f in families if f.module == module and f.active}
    if not wanted:
        raise NoActiveFamilyError(f"{module}: no family of this module is active")
    found, narrow, unknown = {}, set(), set()
    for h in hits:
        fam = wanted.get(h["acc"])
        if fam is None:
            continue
        cond = fam.second_condition
        if cond == "signal_peptide":
            value = (sp_calls or {}).get(h["target"])
            if value is None:
                unknown.add(h["target"])
                continue
            if value != "called":
                continue
        elif cond == "no_tm":
            value = (tm_counts or {}).get(h["target"])
            if value is None:
                unknown.add(h["target"])
                continue
            if value != 0:
                continue
        found.setdefault(h["target"], set()).add(fam.pfam_acc)
        if fam.narrow_ok and fam.tier != UNTESTED:
            narrow.add(h["target"])
    with_narrow = module in NARROW_MODULES
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        accs = sorted(found.get(p.id, ()))
        if p.id in unknown and not accs:
            rows.append({"id": p.id, "state": "error"})
            continue
        row = {"id": p.id, "state": "ok", "hit": "1" if accs else "0", "families": ",".join(accs)}
        if with_narrow:
            row["hit_narrow"] = "1" if p.id in narrow else "0"
        rows.append(row)
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
