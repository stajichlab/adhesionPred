import gzip
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/secretome_discovery/pools.py"
spec = importlib.util.spec_from_file_location("pools", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_composition_counts_percent_and_sequons():
    c = m.composition("CPSTNASNPT")  # N-A-S sequon at 5; N-P-T is not one
    assert c["length"] == 10
    assert c["pct_cys"] == 10.0 and c["pct_pro"] == 20.0 and c["pct_st"] == 40.0
    assert c["sequons_per_100aa"] == 10.0


def test_sequon_needs_ns_or_nt_and_no_proline_in_the_middle():
    assert m.composition("NASNGTNPS")["sequons_per_100aa"] == round(100 * 2 / 9, 2)


def test_er_retention_signal_is_a_c_terminal_hdel_or_kdel():
    assert m.er_retained("MKAAAHDEL")
    assert m.er_retained("MKAAAKDEL")
    assert not m.er_retained("MKAHDELAA")


def test_pfam_reader_gives_protein_to_unversioned_accessions(tmp_path):
    f = tmp_path / "d.domtblout.gz"
    line = "P1 - 100 CFEM PF05730.17 70 1e-9 30 0 1 1 1e-9 1e-9 30 0 1 60 1 60 1 60 0.9 desc\n"
    with gzip.open(f, "wt") as fh:
        fh.write(
            "# comment\n" + line + line.replace("P1", "P2").replace("PF05730.17", "PF00624.24")
        )
    got = m.read_pfam(f)
    assert got == {"P1": {"PF05730"}, "P2": {"PF00624"}}


def test_predgpi_reader_keeps_only_gpi_anchor_features(tmp_path):
    f = tmp_path / "g.gff3.gz"
    with gzip.open(f, "wt") as fh:
        fh.write(
            "P1\tPredGPI\tChain\t1\t9\t1.0\t.\t.\te\nP2\tPredGPI\tGPI-anchor\t5\t5\t1.0\t.\t.\te\n"
        )
    assert m.read_gpi(f) == {"P2"}


def test_pools_exclude_membrane_and_er_proteins_and_split_by_mechanism():
    seqs = {"a": "MKKK", "b": "MKKK", "c": "MKKKHDEL", "d": "MKKK", "e": "MKKK"}
    sp = {"a", "b", "c", "d"}
    no_mech = {"a", "c", "e"}
    n_tm_mature = {"a": 0, "b": 0, "c": 0, "d": 1, "e": 0}
    s, u = m.pools(seqs, sp, no_mech, n_tm_mature)
    assert s == {"a", "b"}  # c is ER-retained, d has a mature TM helix, e has no signal peptide
    assert u == {"a"}  # b has mechanism evidence
