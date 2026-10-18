"""refresh_checks.py models-searched, and the provenance written by pfam_hmmsearch.sbatch (stub tools, no HMMER)."""

import csv
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SH = ROOT / "scripts" / "sorting_hat"
TABLE = ROOT / "data" / "sorting_hat" / "family_table.tsv"
spec = importlib.util.spec_from_file_location(
    "refresh_checks", SH / "refresh_3_10" / "refresh_checks.py"
)
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)

# the real old format (the 14 work directories of 2026-10-08): no model list
OLD_PROVENANCE = {
    "pfam_hmm": "/db/pfam/current/Pfam-A.hmm",
    "resolved": "/db/pfam/Pfam38.2/Pfam-A.hmm",
    "release": "38.2",
    "sha256": "4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e",
    "hmmer": "# HMMER 3.4 (Aug 2023); http://hmmer.org/",
}


def table_models():
    with open(TABLE, encoding="utf-8") as fh:
        return [
            {"acc": r["pfam_acc"], "name": r["name"]} for r in csv.DictReader(fh, delimiter="\t")
        ]


def write(tmp_path, record):
    p = tmp_path / "provenance.json"
    p.write_text(json.dumps(record))
    return str(p)


def test_old_provenance_without_a_list_gives_exit_4(tmp_path, capsys):
    assert rc.models_searched(write(tmp_path, OLD_PROVENANCE), str(TABLE)) == 4
    assert "NO RECORDED MODEL LIST" in capsys.readouterr().out


def test_new_provenance_with_every_table_row_passes(tmp_path):
    new = dict(OLD_PROVENANCE, family_table_sha256="x", models=table_models())
    assert rc.models_searched(write(tmp_path, new), str(TABLE)) == 0


def test_a_model_in_the_record_but_not_in_the_table_is_allowed(tmp_path):
    models = table_models() + [{"acc": "PF99999", "name": "Extra"}]
    assert rc.models_searched(write(tmp_path, dict(OLD_PROVENANCE, models=models)), str(TABLE)) == 0


def test_new_provenance_missing_a_table_row_fails_and_names_it(tmp_path, capsys):
    models = table_models()
    dropped = models.pop(0)
    assert rc.models_searched(write(tmp_path, dict(OLD_PROVENANCE, models=models)), str(TABLE)) == 1
    assert f"NOT SEARCHED {dropped['acc']} {dropped['name']}" in capsys.readouterr().out


def test_a_renamed_model_counts_as_missing(tmp_path):
    models = table_models()
    models[0] = dict(models[0], name="Other")
    assert rc.models_searched(write(tmp_path, dict(OLD_PROVENANCE, models=models)), str(TABLE)) == 1


def test_subcommand_exit_codes(tmp_path):
    cmd = [sys.executable, str(SH / "refresh_3_10" / "refresh_checks.py"), "models-searched"]
    old = subprocess.run(cmd + [write(tmp_path, OLD_PROVENANCE), str(TABLE)], capture_output=True)
    assert old.returncode == 4


def _run_pfam_job(tmp_path):
    names = {m["name"]: m["acc"] for m in table_models()}
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (tmp_path / "Pfam38.2").mkdir()
    (tmp_path / "Pfam38.2" / "Pfam-A.hmm").write_text("db")
    (tmp_path / "q.faa").write_text(">q1\nMKV\n")
    stubs = {
        "hmmfetch": 'while read k; do a=$(awk -F"\\t" -v n="$k" \'$2==n{print $1}\' '
        + f"{TABLE}"
        + '); printf "NAME  %s\\nACC   %s.1\\n//\\n" "$k" "$a"; done < "$3"\n',
        "hmmsearch": 'if [ "$1" = -h ]; then echo "# hmmsearch :: x"; echo "# HMMER 3.4 (Aug 2023)"; exit 0; fi\n'
        'while [ $# -gt 0 ]; do [ "$1" = --domtblout ] && o="$2"; shift; done\n'
        'printf "# a\\n# --cut_ga\\n# [ok]\\n" > "$o"\n',
    }
    for name, body in stubs.items():
        (bindir / name).write_text("#!/bin/bash\n" + body)
        (bindir / name).chmod(0o755)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    # /usr/bin/python3.12 is the cluster interpreter; use the test interpreter in the copy that is run
    script = tmp_path / "pfam.sbatch"
    script.write_text(
        (SH / "pfam_hmmsearch.sbatch").read_text().replace("/usr/bin/python3.12", sys.executable)
    )
    env = dict(
        os.environ,
        PATH=f"{bindir}:{os.environ['PATH']}",
        SCRATCH=str(scratch),
        FASTA=str(tmp_path / "q.faa"),
        WORKDIR=str(tmp_path / "w"),
        FAMILY_TABLE=str(TABLE),
        PFAM_RELEASE="38.2",
        PFAM_HMM=str(tmp_path / "Pfam38.2" / "Pfam-A.hmm"),
    )
    cmd = f"module() {{ :; }}; export -f module; bash {script}"
    return subprocess.run(
        ["bash", "-c", cmd], env=env, cwd=tmp_path, capture_output=True, text=True
    ), names


def test_pfam_job_writes_the_model_list_and_the_table_hash_and_old_readers_still_work(tmp_path):
    r, names = _run_pfam_job(tmp_path)
    assert r.returncode == 0, r.stderr
    prov = tmp_path / "w" / "raw" / "pfam" / "provenance.json"
    rec = json.loads(prov.read_text())
    import hashlib

    assert rec["family_table_sha256"] == hashlib.sha256(TABLE.read_bytes()).hexdigest()
    assert {(m["acc"], m["name"]) for m in rec["models"]} == {(a, n) for n, a in names.items()}
    # the keys read by analysis/sorting_hat_run/run_converters.sh and by job 0
    assert (
        rec["hmmer"].split()[2] == "3.4" and rec["release"] == "38.2" and len(rec["sha256"]) == 64
    )
    assert (
        f'"sha256": "{rec["sha256"]}"' in prov.read_text()
        and '"release": "38.2"' in prov.read_text()
    )
    # and the comparison accepts the file that the job wrote
    assert rc.models_searched(str(prov), str(TABLE)) == 0
    assert not list(prov.parent.glob(".tmp.*"))


R310 = SH / "refresh_3_10"


def test_every_refresh_script_is_valid_bash_and_has_no_draft_header():
    for script in [*R310.glob("*.sbatch"), R310 / "lib.sh"]:
        assert subprocess.run(["bash", "-n", str(script)], capture_output=True).returncode == 0
        head = "\n".join(script.read_text().splitlines()[:3])
        assert "DRAFT" not in head and "not submitted" not in head, script.name
        assert "run on HPCC on 2026-10-17" in head, script.name
    assert "pyyaml==6.0.3 numpy==2.5.4" in (R310 / "lib.sh").read_text().split("else", 1)[1]


def _models_case_fragment():
    text = (R310 / "02_pfam_modules.sbatch").read_text()
    start = text.index("  rc=0\n")
    return text[start : text.index("  check_fasta", start)]


def _run_fragment(tmp_path, exit_code, allow):
    stub = tmp_path / "py"
    stub.write_text(f"#!/bin/bash\nexit {exit_code}\n")
    stub.chmod(0o755)
    env = dict(os.environ, PY=str(stub), R="r", W="w", FAMILY_TABLE="t", p="P1")
    env.pop("ALLOW_UNRECORDED_MODELS", None)
    if allow:
        env["ALLOW_UNRECORDED_MODELS"] = "1"
    script = (
        'fail() { echo "FATAL: $*" >&2; exit 1; }\n' + _models_case_fragment() + "echo continued\n"
    )
    return subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)


def test_job_2_with_provenance_without_a_list_fails_unless_the_override_is_set(tmp_path):
    r = _run_fragment(tmp_path, 4, allow=False)
    assert r.returncode == 1 and "FATAL" in r.stderr and "continued" not in r.stdout
    r = _run_fragment(tmp_path, 4, allow=True)
    assert r.returncode == 0 and "continued" in r.stdout
    assert r.stdout.count("WARNING: P1:") == 1


def test_job_2_fails_on_a_missing_table_row_even_with_the_override(tmp_path):
    r = _run_fragment(tmp_path, 1, allow=True)
    assert r.returncode == 1 and "FATAL" in r.stderr and "continued" not in r.stdout
    assert _run_fragment(tmp_path, 0, allow=False).stdout.strip() == "continued"


def test_job_4_does_not_write_the_final_name_when_the_values_differ():
    text = (R310 / "04_truth_repeat_call.sbatch").read_text()
    body = text.split('mkdir -p "$W/status/calls"', 1)[1]
    rejected, final = body.split('if [ "$rc" -ne 0 ]; then', 1)[1].split("\nfi\n", 1)
    assert ".rejected" in rejected and 'exit "$rc"' in rejected
    assert 'mv "$t" "$W/status/calls/tandem_repeat_protein.$ext.rejected"' in rejected
    assert 'mv "$t" "$W/status/calls/tandem_repeat_protein.$ext"' not in rejected
    assert 'mv "$t" "$W/status/calls/tandem_repeat_protein.$ext"' in final
