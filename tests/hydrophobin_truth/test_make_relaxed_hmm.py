import importlib.util
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/make_relaxed_hmm.py"
spec = importlib.util.spec_from_file_location("make_relaxed_hmm", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

HMM = """HMMER3/f [3.4 | Aug 2023]
NAME  Toy
ACC   PF99999.1
LENG  2
GA    23.00 22.00;
TC    23.10 22.10;
NC    22.90 21.90;
HMM          A
//
HMMER3/f [3.4 | Aug 2023]
NAME  Toy2
ACC   PF99998.1
LENG  2
GA    30.00 29.00;
HMM          A
//
"""


def test_rewrite_ga_sets_sequence_cutoff_and_a_very_low_domain_cutoff_for_every_model():
    out = m.rewrite_ga(HMM, 8.5)
    assert out.count("GA    8.50 -1000.00;") == 2
    assert "23.00 22.00" not in out and "30.00 29.00" not in out
    assert "TC    23.10 22.10;" in out  # other cutoff lines are untouched


def test_rewrite_ga_refuses_a_file_without_ga_lines():
    with pytest.raises(ValueError):
        m.rewrite_ga("HMMER3/f\nNAME x\n//\n", 5.0)


@pytest.mark.skipif(shutil.which("hmmsearch") is None, reason="hmmer is not on PATH")
def test_cut_ga_run_equals_offline_threshold_including_a_protein_with_all_domains_below_zero(
    tmp_path,
):
    # tool test: a real Pfam model, hmmsearch on sequences with a range of scores
    pfam = Path("/bigdata/operations/pkgadmin/srv/projects/db/pfam/current/Pfam-A.hmm")
    if not pfam.exists():
        pytest.skip("Pfam-A.hmm not available")
    hmm = tmp_path / "one.hmm"
    with open(hmm, "w") as f:
        subprocess.run(["hmmfetch", str(pfam), "Hydrophobin"], stdout=f, check=True)
    # a family-like sequence (from the model consensus) and a decoy
    cons = subprocess.run(
        ["hmmemit", "-c", str(hmm)], capture_output=True, text=True, check=True
    ).stdout
    faa = tmp_path / "q.faa"
    seq = "".join(cons.splitlines()[1:])
    faa.write_text(
        f">cons\n{seq}\n>half\n{seq[: len(seq) // 2]}\nMKKLLAAGGSSTTNNPPQQ\n>decoy\nMKKLLAAGGSSTTNNPPQQ\n"
    )
    tbl = tmp_path / "t.tbl"
    subprocess.run(
        [
            "hmmsearch",
            "-T",
            "-1000",
            "--domT",
            "-1000",
            "--tblout",
            str(tbl),
            "-o",
            "/dev/null",
            str(hmm),
            str(faa),
        ],
        check=True,
    )
    scores = m.score_table(tbl.read_text().splitlines())
    for cutoff in (
        sorted(scores.values())[len(scores) // 2],
        min(scores.values()),
    ):  # the lowest score has domains below 0
        rew = tmp_path / "rew.hmm"
        rew.write_text(m.rewrite_ga(hmm.read_text(), cutoff))
        dom = tmp_path / "d.tbl"
        subprocess.run(
            [
                "hmmsearch",
                "--cut_ga",
                "--domtblout",
                str(dom),
                "-o",
                "/dev/null",
                str(rew),
                str(faa),
            ],
            check=True,
        )
        hits = {
            line.split()[0] for line in dom.read_text().splitlines() if not line.startswith("#")
        }
        assert hits == {k for k, v in scores.items() if v >= cutoff}, cutoff
