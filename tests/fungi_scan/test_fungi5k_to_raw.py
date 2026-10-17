import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/fungi_scan/run/fungi5k_to_raw.py"
spec = importlib.util.spec_from_file_location("fungi5k_to_raw", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_signalp_rows_follow_the_scan_layout():
    gff = ["## gff-version 3", "P1 P1\tSignalP-6.0\tsignal_peptide\t1\t18\t0.75\t.\t.\t."]
    rows = m.signalp_rows(["P1", "P2"], gff)
    assert rows[0].startswith("P1\tSP\t0.250000\t0.750000\tCS pos: 18-19")
    assert rows[1] == "P2\tOTHER\t1.000000\t0.000000\t"


def test_fasta_ids_take_the_first_token(tmp_path):
    f = tmp_path / "x.fa"
    f.write_text(">A desc\nMK\n>B\nMM\n")
    assert m.fasta_ids(f) == ["A", "B"]
