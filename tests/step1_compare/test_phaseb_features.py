import gzip
import json

import pytest
import seqsets
import truth_table
from conftest import load_script

A, B, C, D = "a" * 64, "b" * 64, "c" * 64, "d" * 64
SEQS = {A: "MKSTST" * 20, B: "M" * 84, C: "MS" * 600, D: "MKT" * 10}


def _work(tmp_path, fixtures_dir, drop_gpi=None):
    work = tmp_path / "w"
    out = work / "phaseb"
    unique = []
    crow = 0
    for i, h in enumerate(sorted(SEQS)):
        seq = SEQS[h]
        long = len(seq) > 1022
        unique.append(
            {
                "row": str(i),
                "seq_sha256": h,
                "length": str(len(seq)),
                "cterm_row": str(crow) if long else "",
                "sequence": seq,
            }
        )
        crow += long
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    members = [
        {"set_id": "truth", "source_id": "Scer_SGD", "gene_id": "S1", "seq_sha256": A},
        {"set_id": "truth", "source_id": "Calb_CGD", "gene_id": "C1", "seq_sha256": A},
        {"set_id": "truth", "source_id": "Scer_SGD", "gene_id": "S2", "seq_sha256": C},
        {"set_id": "Cimm", "source_id": "Cimm", "gene_id": "CIMG_1", "seq_sha256": B},
        {"set_id": "Cimm", "source_id": "Cimm", "gene_id": "CIMG_2", "seq_sha256": D},
    ]
    for m in members:
        m["length"] = str(len(SEQS[m["seq_sha256"]]))
    truth_table.write_tsv(out / "sequence_members.tsv.gz", seqsets.MEMBER_COLUMNS, members)
    truth = [
        {"source_id": "Scer_SGD", "gene_id": "S1", "label": "P-ext", "subset": "wall",
         "stratum": "P-gpi", "d8_class": "P-gpi", "homology_only": "no", "role": "train"},
        {"source_id": "Calb_CGD", "gene_id": "C1", "label": "N-int", "subset": "",
         "stratum": "N-int", "d8_class": "", "homology_only": "no", "role": "train"},
        {"source_id": "Scer_SGD", "gene_id": "S2", "label": "N-sec", "subset": "",
         "stratum": "N-sec", "d8_class": "", "homology_only": "yes", "role": "train"},
    ]  # fmt: skip
    truth_table.write_tsv(work / "truth_set_triaged.tsv.gz", list(truth[0]), truth)
    (work / "d8_run.json").write_text(json.dumps({"all_sources": True}))
    (out / "prepare_run.json").write_text(json.dumps({"all_sources": True}))
    fx = fixtures_dir / "phaseb"
    sp = out / "signalp" / "part_000"
    sp.mkdir(parents=True)
    for name, dest in (("signalp_prediction_results.txt", "prediction_results.txt.gz"),
                       ("signalp_output.gff3", "output.gff3.gz")):  # fmt: skip
        (sp / dest).write_bytes(gzip.compress((fx / name).read_bytes()))
    gpi_lines = (fx / "predgpi_scores.tsv").read_text().splitlines(keepends=True)
    if drop_gpi:
        gpi_lines = [x for x in gpi_lines if not x.startswith(drop_gpi)]
    (out / "predgpi").mkdir()
    (out / "predgpi" / "part_000.tsv.gz").write_bytes(gzip.compress("".join(gpi_lines).encode()))
    return work, out


def test_features_join_members_truth_and_embedding_rows(tmp_path, fixtures_dir):
    work, out = _work(tmp_path, fixtures_dir)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 0
    feats = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "features_unique.tsv.gz")}
    assert len(feats) == 4
    assert feats[C]["sp_prediction"] == "SP" and feats[C]["sp_cs_end"] == "24"
    assert feats[A]["gpi_call"] == "weakly" and feats[A]["gpi_omega"] == "80"
    assert feats[A]["ser_thr_frac"] == "0.666667"
    assert feats[C]["cterm_row"] == "0" and feats[A]["cterm_row"] == ""
    rows = truth_table.read_tsv(out / "features.tsv.gz")
    assert len(rows) == 5
    s1 = next(r for r in rows if r["gene_id"] == "S1")
    c1 = next(r for r in rows if r["gene_id"] == "C1")
    assert s1["label"] == "P-ext" and s1["d8_class"] == "P-gpi" and c1["label"] == "N-int"
    assert s1["emb_row"] == c1["emb_row"] == feats[A]["row"]  # one embedding per sequence
    cimm = next(r for r in rows if r["gene_id"] == "CIMG_1")
    assert cimm["label"] == "" and cimm["sp_prediction"] == "OTHER"
    cov = {r["set_id"]: r for r in truth_table.read_tsv(out / "feature_coverage.tsv")}
    assert cov["truth"] == {
        "set_id": "truth",
        "members": "3",
        "unique_sequences": "2",
        "no_signalp": "0",
        "no_predgpi": "0",
        "over_1022": "1",
        "gpi_too_short": "0",
    }
    assert cov["Cimm"]["gpi_too_short"] == "1"
    log = json.loads((out / "features_run.json").read_text())
    assert log["all_sources"] is True and log["missing_signalp"] == 0
    assert log["sp_predictions"] == {"OTHER": 2, "SP": 2}


def test_missing_predgpi_call_stops_unless_allowed(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir, drop_gpi=B)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "1 no PredGPI call" in capsys.readouterr().err
    assert not (out / "features.tsv.gz").exists()
    assert build.main(["--work-dir", str(work), "--allow-missing-calls"]) == 0
    cov = {r["set_id"]: r for r in truth_table.read_tsv(out / "feature_coverage.tsv")}
    assert cov["Cimm"]["no_predgpi"] == "1" and cov["truth"]["no_predgpi"] == "0"


def test_stale_tool_output_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique[:3])
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work), "--allow-missing-calls"]) == 2
    assert "not unique sequences" in capsys.readouterr().err


def test_id_in_two_parts_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    src = out / "predgpi" / "part_000.tsv.gz"
    (out / "predgpi" / "part_001.tsv.gz").write_bytes(src.read_bytes())
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "two parts" in capsys.readouterr().err


def test_truth_member_without_truth_row_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    rows = truth_table.read_tsv(work / "truth_set_triaged.tsv.gz")[:2]
    truth_table.write_tsv(work / "truth_set_triaged.tsv.gz", list(rows[0]), rows)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "Scer_SGD:S2" in capsys.readouterr().err


def test_partial_triage_is_refused(tmp_path, fixtures_dir, capsys):
    work, _ = _work(tmp_path, fixtures_dir)
    (work / "d8_run.json").write_text(json.dumps({"all_sources": False}))
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "03_triage_pm.py" in capsys.readouterr().err


@pytest.mark.parametrize("name", ["signalp", "predgpi"])
def test_missing_tool_output_folder_stops(tmp_path, fixtures_dir, capsys, name):
    import shutil

    work, out = _work(tmp_path, fixtures_dir)
    shutil.rmtree(out / name)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "no J1 output" in capsys.readouterr().err


def test_unreferenced_stale_id_stops_with_its_message(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    part = out / "predgpi" / "part_000.tsv.gz"
    text = gzip.decompress(part.read_bytes()).decode()
    text += "e" * 64 + "\t50\tnone\t0\t\t0.9\t-1.7\n"  # an id that is in no member
    part.write_bytes(gzip.compress(text.encode()))
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    err = capsys.readouterr().err
    assert "PredGPI output has 1 ids that are not unique sequences" in err
    assert not (out / "features.tsv.gz").exists()
