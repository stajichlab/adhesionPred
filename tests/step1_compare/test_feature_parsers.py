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


GPI_HEAD = "id\tlength\tgpi_call\tgpi_prob\tomega\tfpr\tsvm\n"


def _gpi(tmp_path, line, name="g.tsv"):
    path = tmp_path / name
    path.write_text(GPI_HEAD + line)
    return path


@pytest.mark.parametrize(
    "line, message",
    [
        (f"{A}\t90\tweakly\t0.55\t8x\t0.007\t-1\n", r"g.tsv:2: omega '8x'"),
        (f"{A}\t90\tweakly\tabc\t80\t0.007\t-1\n", r"g.tsv:2: gpi_prob 'abc'"),
        (f"{A}\t90\tweakly\t0.55\t80\t0.007\tnan?\n", r"g.tsv:2: svm 'nan\?'"),
        (f"{A}\t90\tweakly\t0.55\t80\t1e\t-1\n", r"g.tsv:2: fpr '1e'"),
        (f"{A}\t90\tweakly\t0.70\t80\t0.007\t-1\n", "expected 0.55"),
        (f"{A}\t90\tprobable\t0.70\t80\t0.007\t-1\n", "expected 0.0015 to 0.005"),
        (f"{A}\t90\tweakly\t0.55\t\t0.007\t-1\n", "omega"),
        (f"{A}\t90\tnone\t0\t80\t0.5\t-1\n", "none row has omega"),
        (f"{A}\t90\tnone\t0\t\t0.005\t-1\n", "expected 0.01 to 1.0"),
        (f"{A}\t30\ttoo_short\t0\t\t\t\n{A}\t30\ttoo_short\t0\t\t\t\n", "g.tsv:3: duplicate"),
        (f"{A}\t30\ttoo_short\t0.55\t\t\t\n", "expected 0.0"),
        (f"{A}\t30\ttoo_short\t0\t10\t\t\n", "too_short row has omega"),
    ],
)
def test_predgpi_bad_rows_stop_with_file_and_line(tmp_path, line, message):
    with pytest.raises(fp.OutputFormatError, match=message):
        fp.parse_predgpi_scores(_gpi(tmp_path, line))


def test_predgpi_fpr_at_the_class_bounds_is_accepted(tmp_path):
    rows = (
        f"{A}\t90\thighly_probable\t1.0\t5\t0.0015\t1\n"
        f"{B}\t90\tprobable\t0.70\t5\t0.0015\t1\n"
        f"{C}\t90\tweakly\t0.55\t5\t0.01\t1\n"
        f"{D}\t90\tnone\t0\t\t0.01\t1\n"
    )
    assert len(fp.parse_predgpi_scores(_gpi(tmp_path, rows))) == 4


@pytest.mark.parametrize(
    "line, message",
    [
        (
            "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n" + f"{A}\tOTHER\tx\t0.1\t\n",
            "p.txt:2",
        ),
        (f"{A}\tSP\t0.1\t0.9\tCS pos: 3-4. Pr: 1.2.3\n", "CS probability"),
    ],
)
def test_signalp_bad_numbers_stop_with_file_and_line(tmp_path, line, message):
    path = tmp_path / "p.txt"
    head = "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
    path.write_text(line if line.startswith("# ID") else head + line)
    with pytest.raises(fp.OutputFormatError, match=message):
        fp.parse_signalp(path)


def test_signalp_gff_bad_end_and_duplicate_id_stop(tmp_path):
    path = tmp_path / "o.gff3"
    row = f"{A}\tSignalP-6.0\tsignal_peptide\t1\t{{end}}\t0.9\t.\t.\t.\n"
    path.write_text("##gff-version 3\n" + row.format(end="2x"))
    with pytest.raises(fp.OutputFormatError, match=r"o.gff3:2: GFF end '2x'"):
        fp.parse_signalp_gff(path)
    path.write_text("##gff-version 3\n" + row.format(end=24) + row.format(end=25))
    with pytest.raises(fp.OutputFormatError, match=r"o.gff3:3: duplicate id"):
        fp.parse_signalp_gff(path)


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


def _wrapper():
    import predgpi_scores  # jobs/predgpi_scores.py; PredGPI itself is imported lazily

    return predgpi_scores


@pytest.mark.parametrize(
    "fpr, call, prob",
    [
        (0.0, "highly_probable", "1.0"),
        (0.0015, "highly_probable", "1.0"),
        (0.0015000001, "probable", "0.70"),
        (0.005, "probable", "0.70"),
        (0.0050000001, "weakly", "0.55"),
        (0.01, "weakly", "0.55"),
        (0.0100000001, "none", "0"),
    ],
)
def test_wrapper_classify_boundaries(fpr, call, prob):
    assert _wrapper().classify(fpr) == (call, prob)


def test_wrapper_length_boundary_is_40_and_41():
    wrapper = _wrapper()
    calls = []

    def predict(seq_t):
        calls.append(seq_t)
        return 5, 0.5, 0.2

    short = wrapper.score_row("s", "A" * 40, predict).split("\t")
    assert short == ["s", "40", "too_short", "0", "", "", ""] and calls == []
    long = wrapper.score_row("l", "A" * 41, predict).split("\t")
    assert long == ["l", "41", "none", "0", "", "0.2", "0.5"] and len(calls) == 1


def test_wrapper_applies_the_residue_substitutions_and_computes_omega():
    wrapper = _wrapper()
    seen = []

    def predict(seq_t):
        seen.append(seq_t)
        return 30, 1.5, 0.001

    row = wrapper.score_row("p", "UZBX" + "K" * 46, predict).split("\t")
    assert seen == ["CAAA" + "K" * 46]
    assert row[2:5] == ["highly_probable", "1.0", "20"] and row[5] == "0.001"


def test_wrapper_read_fasta_removes_inner_spaces(tmp_path):
    path = tmp_path / "in.fasta"
    path.write_text(">p1 desc\nMK LV\nAA G\n>p2\nGG\n")
    assert list(_wrapper().read_fasta(str(path))) == [("p1", "MKLVAAG"), ("p2", "GG")]


def test_wrapper_error_removes_tmp_and_stops(tmp_path, capsys):
    wrapper = _wrapper()
    fasta, out = tmp_path / "in.fasta", tmp_path / "scores.tsv"
    fasta.write_text(f">a\n{'A' * 50}\n>b\n{'A' * 50}\n")

    def boom(seq_t):
        raise RuntimeError("model failed")

    assert wrapper.write_scores(str(fasta), str(out), boom) == 2
    assert "STOP: RuntimeError: model failed" in capsys.readouterr().err
    assert list(tmp_path.glob("scores.tsv*")) == []


def test_wrapper_duplicate_id_removes_tmp_and_stops(tmp_path, capsys):
    wrapper = _wrapper()
    fasta, out = tmp_path / "in.fasta", tmp_path / "scores.tsv"
    fasta.write_text(">a\nMKV\n>a\nMKV\n")
    assert wrapper.write_scores(str(fasta), str(out), lambda s: (1, 1.0, 1.0)) == 2
    assert "duplicate id a" in capsys.readouterr().err
    assert list(tmp_path.glob("scores.tsv*")) == []


def test_wrapper_output_parses_with_the_parser(tmp_path):
    wrapper = _wrapper()
    fasta, out = tmp_path / "in.fasta", tmp_path / "scores.tsv"
    fasta.write_text(f">a\nMKV\n>b\n{'A' * 60}\n>c\n{'S' * 60}\n")
    fprs = iter([0.0001, 0.3])
    assert wrapper.write_scores(str(fasta), str(out), lambda s: (10, 0.7, next(fprs))) == 0
    calls = fp.parse_predgpi_scores(out)
    assert [calls[k].call for k in "abc"] == ["too_short", "highly_probable", "none"]
