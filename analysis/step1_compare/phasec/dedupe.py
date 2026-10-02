"""Dedupe and precedence for the Phase C table (Phase C spec 3.4 item 2, ruling C-6).

This module does not import the surface_glyco dedupe function; it has its own rules and tests.

GO members (truth genes) are grouped by seq_sha256:
- one class in the group: one table row; the members' values are merged (see merge_group).
- pos and neg in the group: the hash is dropped (parent spec 4 step 4); every member is logged
  with reason `both_classes`.
- pos or neg together with excluded: the row becomes excluded and every member is logged with
  reason `class_and_excluded` (plan decision: an excluded label on the same sequence makes the
  label uncertain; 2 hashes in the Phase A data).

T-c rows: a T-c row whose hash has any GO member (pos, neg or excluded, and also a hash dropped
as `both_classes`) is dropped with reason `go_label_wins`: the GO label wins (ruling C-6). The
other T-c rows are grouped by hash into one positive row each.
"""

from labelmap import EXCLUDED, NEG, POS

TABLE_COLUMNS = (
    "seq_sha256",
    "origin",
    "class",
    "label",
    "subset",
    "stratum",
    "d8_class",
    "homology_only",
    "internal_evidence_htp_only",
    "source_ids",
    "gene_ids",
    "species",
    "roles",
    "clades",
    "taxon_ids",
    "length",
    "emb_row",
    "emb_cterm_row",
)
LOG_COLUMNS = ("origin", "source_id", "gene_id", "seq_sha256", "class", "reason", "detail")


def _join(values) -> str:
    return ",".join(sorted({v for v in values if v != ""}))


def _check_same(group: list[dict], column: str, h: str) -> str:
    values = {m[column] for m in group}
    if len(values) != 1:
        raise ValueError(f"hash {h}: members disagree on {column} ({sorted(values)})")
    return values.pop()


def merge_group(h: str, group: list[dict], cls: str) -> dict:
    """One table row from GO members that share seq_sha256 `h`.

    homology_only is `no` when any member has direct evidence (the same sequence carries the
    label without homology codes); internal_evidence_htp_only is `yes` when any member says so;
    list columns hold the sorted unique values, comma separated."""
    return {
        "seq_sha256": h,
        "origin": "go",
        "class": cls,
        "label": _join(m["label"] for m in group),
        "subset": _join(m["subset"] for m in group),
        "stratum": _join(m["stratum"] for m in group),
        "d8_class": _join(m["d8_class"] for m in group),
        "homology_only": "no" if any(m["homology_only"] == "no" for m in group) else "yes",
        "internal_evidence_htp_only": "yes"
        if any(m["internal_evidence_htp_only"] == "yes" for m in group)
        else "no",
        "source_ids": _join(m["source_id"] for m in group),
        "gene_ids": _join(m["gene_id"] for m in group),
        "species": _join(m["species"] for m in group),
        "roles": _join(m["role"] for m in group),
        "clades": _join(m["clade"] for m in group),
        "taxon_ids": _join(m["taxon_id"] for m in group),
        "length": _check_same(group, "length", h),
        "emb_row": _check_same(group, "emb_row", h),
        "emb_cterm_row": _check_same(group, "emb_cterm_row", h),
    }


def merge_go(members: list[dict]) -> tuple[dict[str, dict], list[dict]]:
    """Group GO members by hash. Return (table rows by hash, log rows)."""
    groups: dict[str, list[dict]] = {}
    for m in members:
        groups.setdefault(m["seq_sha256"], []).append(m)
    rows, log = {}, []
    for h in sorted(groups):
        group = groups[h]
        classes = {m["class"] for m in group}
        if {POS, NEG} <= classes:
            for m in group:
                log.append(_log("go", m, m["class"], "both_classes", _join(classes)))
            continue
        if EXCLUDED in classes and len(classes) > 1:
            for m in group:
                log.append(_log("go", m, m["class"], "class_and_excluded", _join(classes)))
            rows[h] = merge_group(h, group, EXCLUDED)
            continue
        rows[h] = merge_group(h, group, classes.pop())
    return rows, log


def _log(origin, m, cls, reason, detail) -> dict:
    return {
        "origin": origin,
        "source_id": m.get("source_id", "T-c"),
        "gene_id": m.get("gene_id", m.get("accession", "")),
        "seq_sha256": m["seq_sha256"],
        "class": cls,
        "reason": reason,
        "detail": detail,
    }


def apply_precedence(
    go_hashes: dict[str, str], tc_rows: list[dict]
) -> tuple[dict[str, dict], list[dict]]:
    """Drop T-c rows whose hash has any GO member; merge the rest by hash.

    `go_hashes` maps every hash of a GO member (kept or dropped) to its class detail (for the
    log). Each T-c row needs seq_sha256, accession, taxon_id, clade, length, emb_row,
    emb_cterm_row. Return (T-c table rows by hash, log rows)."""
    kept: dict[str, list[dict]] = {}
    log = []
    for r in tc_rows:
        h = r["seq_sha256"]
        if h in go_hashes:
            log.append(_log("tc", r, POS, "go_label_wins", go_hashes[h]))
            continue
        kept.setdefault(h, []).append(r)
    rows = {}
    for h in sorted(kept):
        group = kept[h]
        rows[h] = {
            "seq_sha256": h,
            "origin": "tc",
            "class": POS,
            "label": "T-c",
            "subset": "",
            "stratum": "T-c",
            "d8_class": "",
            "homology_only": "",
            "internal_evidence_htp_only": "",
            "source_ids": "T-c",
            "gene_ids": _join(r["accession"] for r in group),
            "species": _join(r["genome"] for r in group),
            "roles": "tc",
            "clades": _join(r["clade"] for r in group),
            "taxon_ids": _join(r["taxon_id"] for r in group),
            "length": _check_same(group, "length", h),
            "emb_row": _check_same(group, "emb_row", h),
            "emb_cterm_row": _check_same(group, "emb_cterm_row", h),
        }
    return rows, log
