import gzip
import os
import shutil
import subprocess

import feature_parsers as fp
import pytest

A, B, C, D = "a" * 64, "b" * 64, "c" * 64, "d" * 64


@pytest.fixture
def phaseb_fixtures(fixtures_dir):
    return fixtures_dir / "phaseb"


def test_signalp_rows_parse_by_header_name(phaseb_fixtures):
    calls = fp.parse_signalp(phaseb_fixtures / "signalp_prediction_results.txt")
    assert set(calls) == {A, B, C, D}
    assert calls[A] == fp.SignalPCall("OTHER", 1.0, 0.0, None, None)
    assert calls[B].sp_prob == pytest.approx(0.000014)
    assert calls[C] == fp.SignalPCall("SP", 0.000293, 0.999663, 24, 0.5487)
    assert calls[D].cs_end == 22 and calls[D].cs_prob == pytest.approx(0.7016)


def test_signalp_gff_agrees_with_prediction_results(phaseb_fixtures):
    calls = fp.parse_signalp(phaseb_fixtures / "signalp_prediction_results.txt")
    ends = fp.parse_signalp_gff(phaseb_fixtures / "signalp_output.gff3")
    assert ends == {C: 24, D: 22}
    fp.check_signalp_consistency(calls, ends)
    with pytest.raises(fp.OutputFormatError, match="disagree"):
        fp.check_signalp_consistency(calls, {C: 25, D: 22})


def test_signalp_reads_gzip(phaseb_fixtures, tmp_path):
    src = phaseb_fixtures / "signalp_prediction_results.txt"
    gz = tmp_path / "prediction_results.txt.gz"
    gz.write_bytes(gzip.compress(src.read_bytes()))
    assert fp.parse_signalp(gz) == fp.parse_signalp(src)


@pytest.mark.parametrize(
    "line, message",
    [
        (f"{A}\tSP\t0.1\t0.9\t\n", "SP row with CS field"),
        (f"{A}\tOTHER\t0.9\t0.1\tCS pos: 3-4. Pr: 0.2\n", "OTHER row with a CS field"),
        (f"{A}\tLIPO\t0.1\t0.9\t\n", "unexpected prediction"),
        (f"{A} OTHER 1.0 0.0\n", "fields"),
    ],
)
def test_signalp_malformed_rows_stop(tmp_path, line, message):
    path = tmp_path / "p.txt"
    path.write_text("# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n" + line)
    with pytest.raises(fp.OutputFormatError, match=message):
        fp.parse_signalp(path)


def test_signalp_missing_header_or_duplicate_id_stops(tmp_path):
    path = tmp_path / "p.txt"
    path.write_text(f"{A}\tOTHER\t1.0\t0.0\t\n")
    with pytest.raises(fp.OutputFormatError, match="before the '# ID' header"):
        fp.parse_signalp(path)
    path.write_text(
        "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        f"{A}\tOTHER\t1.0\t0.0\t\n{A}\tOTHER\t1.0\t0.0\t\n"
    )
    with pytest.raises(fp.OutputFormatError, match="duplicate id"):
        fp.parse_signalp(path)


def test_predgpi_scores_parse(phaseb_fixtures):
    calls = fp.parse_predgpi_scores(phaseb_fixtures / "predgpi_scores.tsv")
    assert calls[A] == fp.GpiCall("weakly", 0.55, 80, 0.0095988, -0.708166)
    assert calls[B] == fp.GpiCall("none", 0.0, None, 0.878876, -1.68145)
    assert calls[C].omega == 1512 and calls[C].prob == 1.0
    assert calls[D] == fp.GpiCall("too_short", 0.0, None, None, None)


def test_predgpi_wrong_columns_stop(tmp_path):
    path = tmp_path / "g.tsv"
    path.write_text("id\tgpi\n" + f"{A}\tyes\n")
    with pytest.raises(fp.OutputFormatError, match="columns"):
        fp.parse_predgpi_scores(path)


def test_ser_thr_fraction():
    assert fp.ser_thr_fraction("STSTAAAA") == 0.5
    assert fp.ser_thr_fraction("") == 0.0


PREDGPI_HOME = os.environ.get("PREDGPI_HOME")


@pytest.mark.skipif(not PREDGPI_HOME, reason="needs `module load predgpi/202001`")
def test_predgpi_wrapper_matches_cli(tmp_path):
    from conftest import STEP1_DIR

    python = shutil.which("python")  # the module's Python 3.9 with numpy
    fasta = tmp_path / "in.fasta"
    text = open(os.path.join(PREDGPI_HOME, "testdata", "test.fasta")).read()
    fasta.write_text(text + ">short1\nMKVLAAGIVALLLAAG\n")
    cli_out, wrap_out = tmp_path / "cli.gff3", tmp_path / "scores.tsv"
    cli = os.path.join(PREDGPI_HOME, "predgpi.py")
    subprocess.run([python, cli, "-f", str(fasta), "-o", str(cli_out), "-m", "gff3"], check=True)
    wrapper = STEP1_DIR / "jobs" / "predgpi_scores.py"
    subprocess.run(
        [python, str(wrapper), "--fasta", str(fasta), "--out", str(wrap_out)], check=True
    )
    scores = fp.parse_predgpi_scores(wrap_out)
    n = 0
    for line in cli_out.read_text().splitlines():
        f = line.split("\t")
        call = scores[f[0]]
        if f[2] == "GPI-anchor":
            assert call.call in ("highly_probable", "probable", "weakly")
            assert call.omega == int(f[3]) and call.prob == float(f[5])
        else:
            assert call.call in ("none", "too_short")
        n += 1
    assert n == len(scores) >= 3
    assert scores["short1"].call == "too_short"
