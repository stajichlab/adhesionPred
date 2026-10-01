import ast
import sys
from pathlib import Path

import paths


def test_workdir_from_env(monkeypatch, tmp_path):
    monkeypatch.setenv("STEP1_WORKDIR", str(tmp_path))
    assert paths.workdir() == tmp_path
    assert paths.downloads_dir() == tmp_path / "downloads"


def test_workdir_from_site_yaml(monkeypatch, tmp_path):
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "site.yaml").write_text("fungi5k_input: /a\nworkdir: /x/y\n")
    monkeypatch.delenv("STEP1_WORKDIR", raising=False)
    monkeypatch.setenv("PROJ_ROOT", str(tmp_path))
    assert paths.repo_root() == tmp_path.resolve()
    assert paths.workdir() == Path("/x/y/step1_compare")


def test_relative_site_yaml_workdir_is_resolved_against_the_repo_root(monkeypatch, tmp_path):
    root = tmp_path / "root"
    (root / "config").mkdir(parents=True)
    (root / "config" / "site.yaml").write_text("workdir: _wd\n")
    monkeypatch.delenv("STEP1_WORKDIR", raising=False)
    monkeypatch.setenv("PROJ_ROOT", str(root))
    monkeypatch.chdir(paths.STEP1_DIR)
    assert paths.workdir() == root.resolve() / "_wd" / "step1_compare"


def test_running_from_the_script_folder_writes_nothing_under_the_repo(
    monkeypatch, tmp_path, fixtures_dir
):
    import shutil

    import manifest
    import truth_table
    from conftest import load_script
    from gaf_fixture import write_golden_gaf

    extract = load_script("01_extract_go_truth")
    root = tmp_path / "root"
    (root / "config").mkdir(parents=True)
    (root / "config" / "site.yaml").write_text("workdir: _wd\n")
    inputs = tmp_path / "in"
    inputs.mkdir()
    shutil.copy(fixtures_dir / "mini.obo", inputs / "go-basic.obo")
    write_golden_gaf(inputs / "golden.gaf")
    rows = [
        {c: "" for c in manifest.MANIFEST_COLUMNS}
        | {"file": f, "mode": "strict", "sha256": manifest.sha256_file(inputs / f)}
        for f in ("go-basic.obo", "golden.gaf")
    ]
    manifest.write_manifest(tmp_path / "manifest.tsv", rows)
    species = [
        {
            "source_id": "Fix_SGD",
            "species": "Fixture yeast",
            "taxon_id": "559292",
            "taxon_filter": "",
            "in_clade": "Saccharomycotina",
            "role": "train",
            "gaf_file": "golden.gaf",
        }
    ]
    truth_table.write_tsv(tmp_path / "species.tsv", list(species[0]), species)
    repo = paths.STEP1_DIR.parents[1]
    before = {p for p in paths.STEP1_DIR.rglob("*") if "__pycache__" not in p.parts}
    monkeypatch.delenv("STEP1_WORKDIR", raising=False)
    monkeypatch.setenv("PROJ_ROOT", str(root))
    monkeypatch.chdir(paths.STEP1_DIR)
    stray = paths.STEP1_DIR / "_wd"
    try:
        argv = ["--species", str(tmp_path / "species.tsv")]
        argv += ["--manifest", str(tmp_path / "manifest.tsv"), "--input-dir", str(inputs)]
        assert extract.main(argv) == 0
        assert (root / "_wd" / "step1_compare" / "truth_set.tsv.gz").exists()
        after = {p for p in paths.STEP1_DIR.rglob("*") if "__pycache__" not in p.parts}
        assert after == before, sorted(str(p.relative_to(repo)) for p in after - before)
    finally:
        if stray.exists() and not any(p == stray for p in before):
            shutil.rmtree(stray)


def test_every_module_imports_only_stdlib_or_local():
    modules = sorted(paths.STEP1_DIR.glob("*.py"))
    assert len(modules) >= 15, modules  # the scan must see the real modules
    local = {p.stem for p in modules}
    scanned = 0
    for path in modules:
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [(node.module or "").split(".")[0]]
            else:
                continue
            for name in names:
                scanned += 1
                assert name in sys.stdlib_module_names or name in local, f"{path.name}: {name}"
    assert scanned >= 50, scanned


def test_import_scan_catches_a_third_party_import(tmp_path):
    tree = ast.parse("import numpy\nfrom pandas import DataFrame\nimport json\n")
    names = [
        (a.name if isinstance(node, ast.Import) else node.module).split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.Import | ast.ImportFrom)
        for a in node.names
    ]
    assert [n for n in names if n not in sys.stdlib_module_names] == ["numpy", "pandas"]


def test_no_shell_script_uses_bash_source():
    # The folder has no shell script; the Python scripts are scanned too, so the test is not empty.
    scripts = sorted(paths.STEP1_DIR.rglob("*.sh")) + sorted(paths.STEP1_DIR.glob("*.py"))
    assert len(scripts) >= 15, scripts
    for script in scripts:
        assert "BASH_SOURCE" not in script.read_text(), script


def test_relative_step1_workdir_resolves_to_an_absolute_path(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("STEP1_WORKDIR", "rel/work")
    got = paths.workdir()
    assert got.is_absolute()
    assert got == (tmp_path / "rel" / "work").resolve()
    assert paths.downloads_dir() == got / "downloads"


def test_every_output_file_has_a_columns_md_heading():
    from conftest import load_script

    names = []
    for script in (
        "01_extract_go_truth",
        "02_attach_sequences",
        "03_triage_pm",
        "04_build_keyword_tier",
    ):
        names += load_script(script).OUTPUT_NAMES
    names += ["keyword_sequences.json", "keyword_sequences.fasta.gz"]  # D10 sequence cache
    assert len(names) >= 16
    headings = [
        line[3:].split()[0]
        for line in (paths.STEP1_DIR / "COLUMNS.md").read_text().splitlines()
        if line.startswith("## ")
    ]
    missing = [n for n in names if n not in headings]
    assert not missing, f"COLUMNS.md has no heading for {missing}"
