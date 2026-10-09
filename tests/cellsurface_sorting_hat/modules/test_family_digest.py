"""The artefact digest of a Pfam module covers the rows of that module only."""

from cellsurface_sorting_hat.modules import pfam
from cellsurface_sorting_hat.modules.cli import _family_digest


def _fam(acc, module, active=True, cond=""):
    return pfam.Family(acc, "n", "c", module, cond, active)


def test_digest_ignores_rows_of_other_modules():
    a = [_fam("PF00001", "pfam_adhesion"), _fam("PF00002", "pfam_allergen")]
    b = [_fam("PF00001", "pfam_adhesion"), _fam("PF00002", "pfam_allergen", active=False)]
    assert _family_digest(a, "pfam_adhesion") == _family_digest(b, "pfam_adhesion")
    assert _family_digest(a, "pfam_allergen") != _family_digest(b, "pfam_allergen")


def test_digest_changes_with_own_rows():
    a = [_fam("PF00001", "pfam_adhesion")]
    b = [_fam("PF00001", "pfam_adhesion"), _fam("PF00003", "pfam_adhesion")]
    c = [_fam("PF00001", "pfam_adhesion", cond="no_tm")]
    assert _family_digest(a, "pfam_adhesion") != _family_digest(b, "pfam_adhesion")
    assert _family_digest(a, "pfam_adhesion") != _family_digest(c, "pfam_adhesion")
