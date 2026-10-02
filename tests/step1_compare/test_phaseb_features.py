import gzip
import json

import manifest
import pytest
import seqhash
import seqsets
import truth_table
from conftest import load_script

_SEQ_TEXT = {"a": "MKSTST" * 20, "b": "M" * 84, "c": "MS" * 600, "d": "MKT" * 10}
# the fixture files use the ids a*64 .. d*64; the tests use the real hashes of the sequences
PLACEHOLDER = {k * 64: seqhash.seq_sha256(v) for k, v in _SEQ_TEXT.items()}
A, B, C, D = (PLACEHOLDER[k * 64] for k in "abcd")
SEQS = {A: _SEQ_TEXT["a"], B: _SEQ_TEXT["b"], C: _SEQ_TEXT["c"], D: _SEQ_TEXT["d"]}
TRUTH_SHA = "7" * 64


def _real_ids(text: str) -> str:
    for old, new in PLACEHOLDER.items():
        text = text.replace(old, new)
    return text


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
    (work / "d8_run.json").write_text(
        json.dumps({"all_sources": True, "truth_set_sha256": TRUTH_SHA})
    )
    (work / "sequence_run.json").write_text(
        json.dumps({"all_sources": True, "truth_set_sha256": TRUTH_SHA})
    )
    truth_file = work / "truth_sequences.tsv.gz"
    truth_table.write_tsv(
        truth_file, ["source_id", "gene_id"], [{"source_id": "s", "gene_id": "g"}]
    )
    (out / "prepare_run.json").write_text(
        json.dumps(
            {
                "all_sources": True,
                "truth_set_sha256": TRUTH_SHA,
                "inputs": {"truth": {"sha256": manifest.sha256_file(truth_file)}},
            }
        )
    )
    fx = fixtures_dir / "phaseb"
    sp = out / "signalp" / "part_000"
    sp.mkdir(parents=True)
    for name, dest in (("signalp_prediction_results.txt", "prediction_results.txt.gz"),
                       ("signalp_output.gff3", "output.gff3.gz")):  # fmt: skip
        (sp / dest).write_bytes(gzip.compress(_real_ids((fx / name).read_text()).encode()))
    gpi_lines = _real_ids((fx / "predgpi_scores.tsv").read_text()).splitlines(keepends=True)
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
    assert log["truth_set_sha256"] == TRUTH_SHA
    for name, path in (
        ("truth_set_triaged.tsv.gz", work / "truth_set_triaged.tsv.gz"),
        ("sequence_members.tsv.gz", out / "sequence_members.tsv.gz"),
        ("unique_sequences.tsv.gz", out / "unique_sequences.tsv.gz"),
    ):
        assert log["input_sha256"][name] == manifest.sha256_file(path)
    for rel in (
        "signalp/part_000/prediction_results.txt.gz",
        "signalp/part_000/output.gff3.gz",
        "predgpi/part_000.tsv.gz",
    ):
        assert log["tool_outputs_sha256"][rel] == manifest.sha256_file(out / rel)
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
    # unique_sequences and members no longer hold D, but the J1 output still does
    work, out = _work(tmp_path, fixtures_dir)
    unique = [
        r for r in truth_table.read_tsv(out / "unique_sequences.tsv.gz") if r["seq_sha256"] != D
    ]
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    members = [
        r for r in truth_table.read_tsv(out / "sequence_members.tsv.gz") if r["seq_sha256"] != D
    ]
    truth_table.write_tsv(out / "sequence_members.tsv.gz", seqsets.MEMBER_COLUMNS, members)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work), "--allow-missing-calls"]) == 2
    assert "not unique sequences" in capsys.readouterr().err


def _rewrite(path, edit):
    rows = truth_table.read_tsv(path)
    cols = list(rows[0])
    truth_table.write_tsv(path, cols, edit(rows))


def test_member_hash_absent_from_unique_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    _rewrite(
        out / "unique_sequences.tsv.gz", lambda rows: [r for r in rows if r["seq_sha256"] != D]
    )
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work), "--allow-missing-calls"]) == 2
    assert "CIMG_2 has seq_sha256" in capsys.readouterr().err


def test_sequence_not_matching_its_id_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)

    def mutate(rows):
        for r in rows:
            if r["seq_sha256"] == B:
                r["sequence"] = "M" * 83 + "K"  # same length, same id, other sequence
        return rows

    _rewrite(out / "unique_sequences.tsv.gz", mutate)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert f"sequence {B} does not match its seq_sha256" in capsys.readouterr().err
    assert not (out / "features.tsv.gz").exists()


def test_lower_case_sequence_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)

    def mutate(rows):
        for r in rows:
            if r["seq_sha256"] == B:
                r["sequence"] = r["sequence"].lower()
        return rows

    _rewrite(out / "unique_sequences.tsv.gz", mutate)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "not cleaned upper case" in capsys.readouterr().err


def test_duplicate_unique_id_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    _rewrite(out / "unique_sequences.tsv.gz", lambda rows: rows + [dict(rows[0])])
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "twice" in capsys.readouterr().err


def test_duplicate_member_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    _rewrite(out / "sequence_members.tsv.gz", lambda rows: rows + [dict(rows[0])])
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "truth:Scer_SGD:S1 occurs twice" in capsys.readouterr().err


def test_duplicate_truth_row_stops(tmp_path, fixtures_dir, capsys):
    work, _ = _work(tmp_path, fixtures_dir)
    _rewrite(work / "truth_set_triaged.tsv.gz", lambda rows: rows + [dict(rows[0])])
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "Scer_SGD:S1 twice" in capsys.readouterr().err


def test_stale_truth_set_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    (work / "sequence_run.json").write_text(
        json.dumps({"all_sources": True, "truth_set_sha256": "8" * 64})
    )
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "truth_set_sha256" in capsys.readouterr().err
    assert not (out / "features.tsv.gz").exists()


def _prepare_edit(out, **changes):
    path = out / "prepare_run.json"
    log = json.loads(path.read_text())
    log.update(changes)
    path.write_text(json.dumps(log))
    return log


def test_prepare_from_another_truth_set_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    _prepare_edit(out, truth_set_sha256="8" * 64)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    err = capsys.readouterr().err
    assert "prepare_run.json truth_set_sha256" in err and "re-run 05" in err
    assert not (out / "features.tsv.gz").exists()
    _prepare_edit(out, truth_set_sha256="")  # 05 could not read sequence_run.json
    assert build.main(["--work-dir", str(work)]) == 2


def test_truth_sequences_file_changed_after_prepare_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    truth_table.write_tsv(
        work / "truth_sequences.tsv.gz",
        ["source_id", "gene_id"],
        [{"source_id": "s", "gene_id": "other"}],
    )
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    err = capsys.readouterr().err
    assert "truth sequences" in err and "re-run 05" in err
    assert not (out / "features.tsv.gz").exists()


def test_prepare_without_a_record_of_the_truth_file_stops(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    _prepare_edit(out, inputs={})
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "truth sequences" in capsys.readouterr().err


def test_truth_set_id_comes_from_sequence_sets(tmp_path, fixtures_dir):
    work, out = _work(tmp_path, fixtures_dir)
    _rewrite(
        out / "sequence_members.tsv.gz",
        lambda rows: [
            {**r, "set_id": "gold" if r["set_id"] == "truth" else r["set_id"]} for r in rows
        ],
    )
    log = _prepare_edit(out)
    _prepare_edit(out, inputs={"gold": log["inputs"]["truth"]})
    sets = tmp_path / "sets.tsv"
    sets.write_text(
        "set_id\tkind\tlocation\tnote\n"
        "gold\ttruth\ttruth_sequences.tsv.gz\t\nCimm\tsite\tx:y\t\n"
    )
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work), "--sets", str(sets)]) == 0
    rows = truth_table.read_tsv(out / "features.tsv.gz")
    s1 = next(r for r in rows if r["gene_id"] == "S1")
    assert s1["set_id"] == "gold" and s1["label"] == "P-ext" and s1["role"] == "train"


def test_partial_prepare_is_refused(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    (out / "prepare_run.json").write_text(json.dumps({"all_sources": False}))
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 2
    assert "prepare_run.json" in capsys.readouterr().err
    (out / "prepare_run.json").write_text("{}")  # all_sources missing
    assert build.main(["--work-dir", str(work)]) == 2


def test_tool_ids_are_the_join_key(tmp_path, fixtures_dir, capsys):
    work, out = _work(tmp_path, fixtures_dir)
    build = load_script("07_build_features")
    assert build.main(["--work-dir", str(work)]) == 0
    base = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "features_unique.tsv.gz")}

    def swap(text):
        return text.replace(C, "@").replace(D, C).replace("@", D)

    for rel in ("signalp/part_000/prediction_results.txt.gz", "signalp/part_000/output.gff3.gz"):
        path = out / rel
        path.write_bytes(gzip.compress(swap(gzip.decompress(path.read_bytes()).decode()).encode()))
    assert build.main(["--work-dir", str(work)]) == 0
    now = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "features_unique.tsv.gz")}
    for col in ("sp_prediction", "sp_prob", "sp_cs_end"):
        assert now[C][col] == base[D][col] and now[D][col] == base[C][col]
    assert base[C]["sp_cs_end"] != base[D]["sp_cs_end"]  # the swap changes the values
    # an id that no longer matches any unique sequence stops
    for rel in ("signalp/part_000/prediction_results.txt.gz", "signalp/part_000/output.gff3.gz"):
        path = out / rel
        path.write_bytes(
            gzip.compress(gzip.decompress(path.read_bytes()).replace(D.encode(), b"f" * 64))
        )
    assert build.main(["--work-dir", str(work)]) == 2
    assert "not unique sequences" in capsys.readouterr().err


def test_outputs_are_byte_stable(tmp_path, fixtures_dir):
    from hashlib import sha256

    work, out = _work(tmp_path, fixtures_dir)
    build = load_script("07_build_features")
    digests = []
    for _ in range(2):
        assert build.main(["--work-dir", str(work)]) == 0
        digests.append({n: sha256((out / n).read_bytes()).hexdigest() for n in build.OUTPUT_NAMES})
    assert digests[0] == digests[1]


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
