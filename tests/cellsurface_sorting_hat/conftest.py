"""Shared fixtures: a toy taxonomy and a helper that writes module tables."""

import csv
import gzip
import json

import pytest

# taxid, parent. 1 root; 10 Fungi; 20 Eurotiomycetes; 30 Eurotiales; 31 Onygenales;
# 40 Aspergillus fumigatus (in 30); 41 Coccidioides immitis (in 31); 42 Aspergillus nidulans (in 30)
TOY_NODES = {1: 1, 10: 1, 20: 10, 30: 20, 31: 20, 40: 30, 41: 31, 42: 30}


@pytest.fixture
def nodes_dmp(tmp_path):
    path = tmp_path / "nodes.dmp"
    path.write_text("".join(f"{t}\t|\t{p}\t|\tno rank\t|\n" for t, p in TOY_NODES.items()))
    return path


@pytest.fixture
def write_module():
    def _write(workdir, name, rows, meta=None, status=None):
        folder = workdir / "modules"
        folder.mkdir(parents=True, exist_ok=True)
        columns = sorted({k for r in rows for k in r} - {"id"})
        with gzip.open(folder / f"{name}.tsv.gz", "wt") as fh:
            writer = csv.DictWriter(fh, ["id"] + columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
        record = {"module": name, "version": "1", "params_hash": "p", "artefact_hash": "a"}
        record.update(meta or {})
        (folder / f"{name}.json").write_text(json.dumps(record))
        if status is not None:
            sfolder = workdir / "status"
            sfolder.mkdir(exist_ok=True)
            (sfolder / f"{name}.json").write_text(
                json.dumps(
                    {
                        "module": name,
                        "version": "1",
                        "params_hash": "p",
                        "artefact_hash": "a",
                        "entries": status,
                    }
                )
            )

    return _write
