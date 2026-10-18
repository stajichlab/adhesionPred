"""Task 4.3 (L3): the generated ``docs/CLASSES.md`` names the families that the calls read.

``CLASSES.md`` is generated from ``docs/class_descriptions.yaml`` and the family table. These tests
read the committed text and the table. They do not import from ``analysis/`` (D22) and need no
network. Task 2.1 of the plan (a test of stale phrases) is not in this tree; only Task 4.3 is here.
"""

import csv
import re
from pathlib import Path

from cellsurface_sorting_hat.engine import load_config
from cellsurface_sorting_hat.family_tier import check, derive_tiers

ROOT = Path(__file__).resolve().parents[2]
CLASSES = (ROOT / "docs" / "CLASSES.md").read_text(encoding="utf-8")
TABLE = ROOT / "data" / "sorting_hat" / "family_table.tsv"
NARROW, BROAD = "wall_family_narrow_domain", "wall_family_domain"


def _families():
    with open(TABLE, encoding="utf-8", newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["module"] == "pfam_adhesion"]


def _row_of(call):
    """The table row of one call in CLASSES.md (the text between the first two pipes names it)."""
    rows = [ln for ln in CLASSES.splitlines() if ln.startswith(f"| `{call}` |")]
    assert len(rows) == 1, f"CLASSES.md has {len(rows)} rows for {call}"
    return rows[0]


def _named(text, family):
    return re.search(rf"(?<![A-Za-z0-9_]){re.escape(family['name'])}(?![A-Za-z0-9_])", text)


def test_there_are_narrow_and_active_families_to_check():
    """A guard: the checks below must not pass because the table has no such family."""
    fams = _families()
    assert any(f["narrow_ok"] == "yes" for f in fams)
    assert any(f["active"] == "yes" and f["narrow_ok"] == "no" for f in fams)


def test_every_narrow_family_is_named_for_the_narrow_call():
    row = _row_of(NARROW)
    missing = [f["name"] for f in _families() if f["narrow_ok"] == "yes" and not _named(row, f)]
    assert missing == [], f"narrow_ok=yes families not named in the {NARROW} row: {missing}"


def test_every_active_family_is_named_for_the_broad_call():
    row = _row_of(BROAD)
    missing = [f["name"] for f in _families() if f["active"] == "yes" and not _named(row, f)]
    assert missing == [], f"active families not named in the {BROAD} row: {missing}"


def test_a_family_that_is_not_narrow_is_not_named_for_the_narrow_call():
    row = _row_of(NARROW)
    wrong = [f["name"] for f in _families() if f["narrow_ok"] != "yes" and _named(row, f)]
    assert wrong == [], f"families without narrow_ok=yes named in the {NARROW} row: {wrong}"


def test_a_family_that_is_not_active_is_not_named_for_the_broad_call():
    row = _row_of(BROAD)
    wrong = [f["name"] for f in _families() if f["active"] != "yes" and _named(row, f)]
    assert wrong == [], f"inactive families named in the {BROAD} row: {wrong}"


def test_classes_md_mentions_each_call_of_categories_yaml():
    missing = [c["name"] for c in load_config().calls if f"`{c['name']}`" not in CLASSES]
    assert missing == []


def test_the_stored_tier_equals_the_derived_tier():
    assert check() == []


def test_a_narrow_family_derives_smoke_or_better():
    tiers = derive_tiers()
    weak = [
        f["name"]
        for f in _families()
        if f["narrow_ok"] == "yes" and tiers[f["pfam_acc"]] not in ("smoke", "estimated")
    ]
    assert weak == []


def test_the_name_check_notices_a_family_removed_from_a_text():
    """The helper is not vacuous: it finds a family in a text and does not find it when removed."""
    family = {"name": "Flo11"}
    assert _named("reads Flo11, GLEYA and PA14", family)
    assert not _named("reads GLEYA and PA14", family)
    assert not _named("reads Flo11b", family)  # a longer name is another family
