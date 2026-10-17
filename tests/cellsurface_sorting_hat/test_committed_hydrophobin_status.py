"""The 8 committed ``pfam_hydrophobin`` status files must stay valid after the Pfam rebuild.

A rebuilt run record must keep the identity of the modules that have no narrow call. The status
files hold ``params_hash`` and ``artefact_hash``; the real ``params_hash`` also holds the identity
of the R0 table of the real run, which a test cannot rebuild. So the test takes the hash from the
file and checks that nothing else in the identity moved: name, version and artefact hash (Pfam
sha256 plus the family digest of the real table), and that the params carry no new key.
"""

import json
from pathlib import Path

import pytest

from cellsurface_sorting_hat.modules.cli import main
from cellsurface_sorting_hat.status import ModuleIdentity, load_status_source, resolve_entry

ROOT = Path(__file__).resolve().parents[2]
FILES = sorted((ROOT / "docs/reports/data/sorting_hat/call_status/hydrophobin").glob("pfam_*.json"))
FAMILY_TABLE = ROOT / "data/sorting_hat/family_table.tsv"
OPTIONS = "# Option settings:     hmmsearch --cut_ga --cpu 2 --noali fam.hmm in.fasta\n"


def test_there_are_hydrophobin_status_files():
    assert FILES  # the parametrised test below must not be empty


class SameTaxon:
    """A lineage in which a taxon is a descendant of itself only."""

    def is_descendant_or_self(self, taxon, tested):
        return taxon == tested

    def depth(self, taxon):
        return 1


@pytest.fixture
def rebuilt(tmp_path, write_module):
    """Run the ``pfam`` command on the real family table; return the records by module."""
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKTAYIAKQRQISFVKSHFSRQ\n>B\nMNLLPQWERTYIPASDFG\n")
    domtbl = tmp_path / "d.domtbl"
    domtbl.write_text(OPTIONS + "# [ok]\n")
    wd = tmp_path / "wd"
    for name, col in (("sp", "call"), ("tm", "n_tm_mature")):
        value = "called" if col == "call" else "0"
        write_module(wd, name, [{"id": i, "state": "ok", col: value} for i in "AB"])

    def build(sha):
        argv = ["pfam", "--fasta", str(fasta), "--workdir", str(wd), "--domtbl", str(domtbl)]
        argv += ["--family-table", str(FAMILY_TABLE), "--pfam-release", "38.2"]
        argv += ["--pfam-sha256", sha, "--hmmer-version", "3.4", "--sp-module", "sp"]
        argv += ["--tm-module", "tm"]
        assert main(argv) == 0
        return {
            m: json.loads((wd / "modules" / f"{m}.json").read_text())
            for m in ("pfam_hydrophobin", "pfam_hsba", "pfam_allergen")
        }

    return build


@pytest.mark.parametrize("path", FILES, ids=lambda p: p.name)
def test_a_committed_hydrophobin_status_file_is_not_stale_after_a_rebuild(path, rebuilt):
    source = load_status_source(path)
    sha = source.identity.artefact_hash.split(":")[0]
    record = rebuilt(sha)["pfam_hydrophobin"]
    # the params hash of the real run holds the identity of its real R0 table; take it from the file
    running = ModuleIdentity(
        record["module"],
        record["version"],
        source.identity.params_hash,
        record["artefact_hash"],
    )
    entry, _, reason = resolve_entry(source, running, source.entries[0].taxa[0], SameTaxon())
    assert "stale" not in reason, reason
    assert entry is not None


def test_the_params_of_the_modules_without_a_narrow_call_hold_no_snapshot_keys(rebuilt):
    for module, record in rebuilt("ab").items():
        for cond in record["params"]["conditions"].values():
            assert set(cond) == {"module", "params_hash", "artefact_hash"}, module
        assert "narrow_families" not in record["params"], module
        assert set(record["condition_snapshots"]) == {"sp_module", "tm_module"}, module
