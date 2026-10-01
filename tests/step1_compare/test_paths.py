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


def test_every_module_imports_only_stdlib_or_local():
    modules = sorted(paths.STEP1_DIR.glob("*.py"))
    local = {p.stem for p in modules}
    for path in modules:
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [(node.module or "").split(".")[0]]
            else:
                continue
            for name in names:
                assert name in sys.stdlib_module_names or name in local, f"{path.name}: {name}"


def test_no_shell_script_uses_bash_source():
    for script in paths.STEP1_DIR.rglob("*.sh"):
        assert "BASH_SOURCE" not in script.read_text(), script
