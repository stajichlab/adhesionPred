"""Model directory lookup, file names, and the no-model case."""

import sys

import pytest

from surface_glyco import config
from surface_glyco.scripts import predict as predict_mod


def test_model_filename():
    assert config.model_filename("esm2_t6_8M_UR50D") == "surface_glyco_model_esm2_t6_8M_UR50D.pkl"


def test_packaged_dir_has_readme_and_no_pickles():
    pkg = config.PACKAGED_MODELS_DIR
    assert (pkg / "README.md").is_file()
    assert list(pkg.glob("*.pkl")) == []


def test_stray_models_dir_in_cwd_does_not_override(tmp_path, monkeypatch):
    (tmp_path / "models").mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("SURFACE_GLYCO_MODELS_DIR", raising=False)
    assert config.get_models_dir() == config.PACKAGED_MODELS_DIR


def test_env_var_overrides(tmp_path, monkeypatch):
    monkeypatch.setenv("SURFACE_GLYCO_MODELS_DIR", str(tmp_path))
    assert config.get_models_dir() == tmp_path


def test_predict_without_a_model_explains_what_to_do(tmp_path, monkeypatch, capsys):
    fa = tmp_path / "in.fa"
    fa.write_text(">a\nMKT\n")
    monkeypatch.setenv("SURFACE_GLYCO_MODELS_DIR", str(tmp_path / "empty"))
    monkeypatch.setattr(sys, "argv", ["surface_glyco_predict", "--input", str(fa)])
    with pytest.raises(SystemExit) as e:
        predict_mod.cli()
    assert e.value.code == 1
    captured = capsys.readouterr()  # read once: a second call returns empty
    assert "surface_glyco_train" in captured.err and "--model" in captured.err


def test_train_default_output_is_not_inside_the_package(tmp_path, monkeypatch):
    from surface_glyco.scripts.train import default_output_path

    monkeypatch.chdir(tmp_path)
    out = default_output_path("esm2_t6_8M_UR50D")
    assert out == tmp_path / "models" / "surface_glyco_model_esm2_t6_8M_UR50D.pkl"
    assert config.PACKAGED_MODELS_DIR not in out.parents
