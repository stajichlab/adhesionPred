"""One derived label per protein: ``primary_label`` (D1, D14 of the 2026-10-17 spec).

The calls stay independent flags. The label is a summary for a reader, derived for the default gate
(``categories.yaml`` ``default_gate``). Precedence, first match wins:

1. ``cell_wall_adhesion_candidate`` is called.
2. ``surface_attachment_candidate`` is called (the basis names the evidence that held).
3. ``other_surface_no_mechanism`` is called.
4. ``other_not_surface`` is called.
5. ``mechanism_without_class``: a mechanism is called that neither composite explains, whatever the
   signal peptide says (no signal peptide, a mechanism no composite reads, or only the antigen rank).
6. ``not_assessable``.

The mechanism calls are declared in ``MECHANISM_CALLS``. Allergen calls and ``serodiagnostic_marker_candidate`` are not label sources. The antigen rank is a
mechanism of the ``other_*`` calls, so it counts for label 5. The function is total: every protein
gets exactly one label. Labels 3 and 4 cannot hold together with labels 1, 2 or each other.
"""

from cellsurface_sorting_hat.engine import variant_label
from cellsurface_sorting_hat.logic import CALLED

CWAC = "cell_wall_adhesion_candidate"
SAC = "surface_attachment_candidate"
OSN = "other_surface_no_mechanism"
ONS = "other_not_surface"
MWC = "mechanism_without_class"
NAS = "not_assessable"
LABELS = (CWAC, SAC, OSN, ONS, MWC, NAS)


# The mechanism calls of label 5. This is the one place to add a new mechanism call; it must also
# be added to the `mechanism` list of both `other_*` calls in categories.yaml. `_mechanisms`
# refuses a config whose `other_*` lists name a call that is not declared here, and
# test_labels.py fails when a call declared here is missing from those lists.
MECHANISM_CALLS = (
    "tandem_repeat_protein",
    "wall_family_domain",
    "hydrophobin_domain",
    "hsba_domain",
    "cocci_specificity_rank_top15",
)


def _mechanisms(cfg):
    """``MECHANISM_CALLS``. A call in an ``other_*`` mechanism list that is not declared there
    would never reach label 5; that is an error, not a silent gap."""
    listed = {m for c in cfg.calls if c.get("kind") == "other" for m in c["mechanism"]}
    missing = sorted(listed - set(MECHANISM_CALLS))
    if missing:
        raise ValueError(
            f"mechanism call(s) {missing} are in an other_* list of the config but not in "
            "labels.MECHANISM_CALLS; add them there"
        )
    return list(MECHANISM_CALLS)


def primary_labels(records, cfg):
    """``{protein: (label, basis)}`` for the default gate, in the order the proteins first appear.

    The basis of ``surface_attachment_candidate`` is its ``other_basis`` (the evidence that held);
    of ``mechanism_without_class`` the called mechanisms; of ``not_assessable`` the mechanism
    inputs that were not assessable (``other_basis`` of an ``other_*`` call). It is empty otherwise."""
    gate = variant_label(cfg.default_gate)
    mechanisms = _mechanisms(cfg)
    values, basis, proteins = {}, {}, []
    for r in records:
        if r.protein not in basis:
            proteins.append(r.protein)
            basis[r.protein] = {}
            values[r.protein] = {}
        if r.variant in ("", gate):
            values[r.protein][r.call] = r.value
            basis[r.protein][r.call] = r.other_basis
    out = {}
    for pid in proteins:
        v, b = values[pid], basis[pid]
        called = [m for m in mechanisms if v.get(m) == CALLED]
        if v.get(CWAC) == CALLED:
            out[pid] = (CWAC, "")
        elif v.get(SAC) == CALLED:
            out[pid] = (SAC, b.get(SAC, ""))
        elif v.get(OSN) == CALLED:
            out[pid] = (OSN, "")
        elif v.get(ONS) == CALLED:
            out[pid] = (ONS, "")
        elif called:
            out[pid] = (MWC, ",".join(called))
        else:
            unknown = b.get(OSN) or b.get(ONS) or ""
            out[pid] = (NAS, unknown)
    return out


def label_counts(labels):
    """``{label: number of proteins}`` for every label, also the ones with no protein."""
    counts = dict.fromkeys(LABELS, 0)
    for label, _ in labels.values():
        counts[label] += 1
    return counts
