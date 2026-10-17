import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/fungi_scan/summarize_scan.py"
spec = importlib.util.spec_from_file_location("summarize_scan", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def rows():
    return [
        {
            "protein": "a",
            "call": "signal_peptide_protein",
            "variant": "R0",
            "value": "called",
            "other_basis": "",
        },
        {
            "protein": "b",
            "call": "signal_peptide_protein",
            "variant": "R0",
            "value": "not_called",
            "other_basis": "",
        },
        {
            "protein": "a",
            "call": "tandem_repeat_protein",
            "variant": "",
            "value": "called",
            "other_basis": "",
        },
        {
            "protein": "a",
            "call": "surface_attachment_candidate",
            "variant": "R0",
            "value": "called",
            "other_basis": "tandem_repeat_protein,hydrophobin_domain",
        },
        {
            "protein": "c",
            "call": "surface_attachment_candidate",
            "variant": "R0",
            "value": "called",
            "other_basis": "hsba_domain",
        },
        {
            "protein": "d",
            "call": "surface_attachment_candidate",
            "variant": "R0",
            "value": "not_assessable",
            "other_basis": "",
        },
    ]


def test_count_calls_counts_called_per_call_and_variant():
    c = m.count_calls(rows())
    assert c[("signal_peptide_protein", "R0")]["called"] == 1
    assert c[("signal_peptide_protein", "R0")]["not_called"] == 1
    assert c[("surface_attachment_candidate", "R0")]["called"] == 2
    assert c[("surface_attachment_candidate", "R0")]["not_assessable"] == 1


def test_basis_combinations():
    out = m.basis_counts(rows(), "surface_attachment_candidate")
    assert out["tandem_repeat_protein,hydrophobin_domain"] == 1 and out["hsba_domain"] == 1


def test_family_counts_split_the_comma_list():
    mod = [
        {"hit": "1", "families": "PF05730,PF07691"},
        {"hit": "1", "families": "PF05730"},
        {"hit": "0", "families": ""},
    ]
    assert m.family_counts(mod) == {"PF05730": 2, "PF07691": 1}


def test_rate_per_1000():
    assert m.per_1000(5, 10000) == 0.5
    assert m.per_1000(0, 0) == 0.0
