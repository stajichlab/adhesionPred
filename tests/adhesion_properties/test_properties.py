"""Tests for analysis/adhesion_properties/properties.py."""

import math

from properties import sequence_properties


def test_sequence_properties_pure_ser_thr_pro():
    result = sequence_properties("SSSTTTPPP")
    assert math.isclose(result["pct_ser"], 100 / 3, rel_tol=1e-6)
    assert math.isclose(result["pct_thr"], 100 / 3, rel_tol=1e-6)
    assert math.isclose(result["pct_pro"], 100 / 3, rel_tol=1e-6)
    assert math.isclose(result["pct_ser_thr_pro"], 100.0, rel_tol=1e-6)
    assert result["aromaticity"] == 0.0
    assert result["net_charge_ph7"] == 0.0


def test_sequence_properties_aromaticity():
    result = sequence_properties("FWYA")
    assert math.isclose(result["aromaticity"], 75.0, rel_tol=1e-6)


def test_sequence_properties_net_charge():
    result = sequence_properties("KRDE")
    assert result["net_charge_ph7"] == 0.0
    result2 = sequence_properties("KKRR")
    assert result2["net_charge_ph7"] == 4.0


def test_sequence_properties_empty_string_returns_nan():
    result = sequence_properties("")
    assert math.isnan(result["pct_ser"])
    assert math.isnan(result["mean_hydrophobicity"])


def test_sequence_properties_handles_unknown_residue():
    """An unusual/ambiguous residue code (X) should not raise, and should
    still be counted in length for percentage denominators."""
    result = sequence_properties("SSXX")
    assert math.isclose(result["pct_ser"], 50.0, rel_tol=1e-6)
