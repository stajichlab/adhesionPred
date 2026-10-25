"""Tests of analysis/aga1_class/feature_classifier.py (synthetic data only) and of the committed README hashes."""

import hashlib
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis/aga1_class/feature_classifier.py"
DATA = ROOT / "docs/reports/data/aga1_feature_classifier_2026-10-25"
spec = importlib.util.spec_from_file_location("feature_classifier", SCRIPT)
fc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fc)


# ---------------------------------------------------------------- feature functions
def test_sequon_counts_overlaps_and_excludes_proline():
    assert fc.sequon_per100("NAS") == pytest.approx(100 / 3)
    assert fc.sequon_per100("NPS") == 0
    assert fc.sequon_per100("NNST") == pytest.approx(200 / 4)  # N-N-S and N-S-T
    assert fc.sequon_per100("AAAA") == 0
    assert np.isnan(fc.sequon_per100(""))


def test_composition_fractions():
    c = fc.composition("SSTPNQGACDEKRAVILMFW")
    assert c["frac_S"] == pytest.approx(2 / 20)
    assert c["frac_ST"] == pytest.approx(3 / 20)
    assert c["frac_charged"] == pytest.approx(4 / 20)  # D E K R
    assert c["frac_hydrophobic"] == pytest.approx(8 / 20)  # A (twice) V I L M F W
    with pytest.raises(fc.FeatureError):
        fc.composition("")


def test_entropy_low_for_low_complexity():
    assert fc.aa_entropy("S" * 50) == 0
    assert fc.aa_entropy(fc.AA20) == pytest.approx(np.log2(20))


def test_st_in_disordered_uses_only_disordered_residues():
    assert fc.st_in_disordered("SSAA", [0.9, 0.9, 0.1, 0.1]) == 1.0
    assert np.isnan(fc.st_in_disordered("SSAA", [0.1] * 4))
    with pytest.raises(fc.FeatureError):
        fc.st_in_disordered("SS", [0.9])


def test_gpi_na_is_missing_with_indicator_never_zero():
    yes, miss, conf = fc.gpi_encoding("NA", "NA")
    assert np.isnan(yes) and np.isnan(conf) and miss == 1.0
    assert fc.gpi_encoding("no", "NA") == (0.0, 0.0, 0.0)
    assert fc.gpi_encoding("yes", "high") == (1.0, 0.0, 3.0)
    with pytest.raises(fc.FeatureError):
        fc.gpi_encoding("maybe", "high")
    with pytest.raises(fc.FeatureError):
        fc.gpi_encoding("yes", "NA")


def _synthetic_inputs():
    seqs = {"P1": "NASTSTSTSTS" * 10, "P2": "MKLVNPSAAAGGG" * 8, "P3": "MKLV" * 5}  # P3 is 20 aa
    proteins = pd.DataFrame(
        {"accession": list(seqs), "gene": ["a", "b", ""], "role": ["class_positive", "hard_negative", "hard_negative"],
         "pos_tier": ["A", "", ""], "clade_coarse": ["X"] * 3, "cluster_id": ["C1", "C2", ""]}
    )  # fmt: skip
    preds = pd.DataFrame(
        {
            "accession": list(seqs),
            "signal_peptide_pred": ["yes", "no", "no"],
            "sp_prob": ["0.9", "0.1", "0.2"],
            "gpi_pred": ["yes", "no", "NA"],
            "gpi_confidence_class": ["high", "NA", "NA"],
            "disorder_n_res": [len(s) for s in seqs.values()],
            "disorder_frac": ["0.8", "0.1", "0.0"],
            "disorder_longest_run": ["90", "5", "0"],
            "disorder_mean": ["0.7", "0.2", "0.1"],
        }
    )
    stalk = pd.DataFrame(
        {"accession": list(seqs), "n_res": [len(s) for s in seqs.values()], "stalk_run_len": ["80", "0", "0"],
         "max_st_window_frac": ["0.9", "0", "0"], "ctail_hydro": ["-1.0", "0.5", "NA"]}
    )  # fmt: skip
    rep = pd.DataFrame(
        {
            "protein": ["P1", "P2"],
            "rep_coverage": [0.5, 0.0],
            "rep_n_copies": [10, 0],
            "rep_z_seq": [5.0, 1.0],
        }
    )  # P3 < 80 aa is absent from the detector table
    raw = {a: (list(s), [0.9 if i < 50 else 0.1 for i in range(len(s))]) for a, s in seqs.items()}
    return proteins, preds, stalk, rep, raw


def test_build_features_missing_values_and_no_text_columns():
    df = fc.build_features(*_synthetic_inputs())
    p3 = df[df.accession == "P3"].iloc[0]
    assert p3["rep_missing"] == 1.0 and np.isnan(p3["rep_coverage"])  # < 80 aa: missing, not 0
    assert p3["gpi_missing"] == 1.0 and np.isnan(p3["gpi_yes"])  # GPI NA: missing, not 0
    assert np.isnan(p3["ctail_hydro"])
    assert p3["cluster_id"] == "SINGLETON_P3"  # cluster-less protein gets its own cluster
    p1 = df[df.accession == "P1"].iloc[0]
    assert p1["rep_missing"] == 0.0 and p1["gpi_yes"] == 1.0 and p1["length"] == 110
    assert list(df.columns) == fc.META + fc.FEATURES
    assert all(df[c].dtype != object for c in fc.FEATURES)


def test_build_features_refuses_mismatched_tables():
    proteins, preds, stalk, rep, raw = _synthetic_inputs()
    with pytest.raises(fc.FeatureError):
        fc.build_features(proteins, preds.iloc[:2], stalk, rep, raw)
    preds2 = preds.copy()
    preds2.loc[0, "disorder_n_res"] = 5
    with pytest.raises(fc.FeatureError):
        fc.build_features(proteins, preds2, stalk, rep, raw)


# ---------------------------------------------------------------- labels and rules
def _table(n_pos=12, n_neg=60, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_pos + n_neg):
        pos = i < n_pos
        r = {c: rng.normal() for c in fc.FEATURES}
        r.update(
            accession=f"A{i}", gene="", role="class_positive" if pos else "hard_negative", pos_tier="A" if pos else "",
            clade_coarse="X", cluster_id=f"C{i // 3}",
        )  # fmt: skip
        r["frac_ST"] += 3 * pos  # a planted signal
        r["length"] = abs(r["length"]) * 100 + 50
        r["disorder_longest_run"] = 250.0 if pos else 10.0
        r["sp_yes"], r["gpi_yes"] = 1.0, 1.0
        r["rep_coverage"], r["rep_n_copies"] = (0.5, 4.0) if i % 2 == 0 else (np.nan, np.nan)
        rows.append(r)
    return pd.DataFrame(rows)[fc.META + fc.FEATURES]


def test_label_sets_exclude_tier_b_and_set_aside():
    df = _table(4, 6)
    df.loc[0, "pos_tier"] = "B"
    df.loc[5, "role"] = "excluded_indirect_regulator"
    df.loc[6, "accession"] = "Q01589"
    elig, pos, neg, extra = fc.label_sets(df)
    assert not elig[0] and extra[0] and not elig[5] and extra[5]
    assert elig[6] and neg[6]
    elig2, _, _, extra2 = fc.label_sets(df, drop=fc.NEW_COL_SED)
    assert not elig2[6] and not extra2[6]  # left out entirely


def test_rules_missing_repeats_never_call_narrow():
    df = _table(4, 6)
    broad, narrow = fc.rule_flags(df)
    assert broad[:4].all() and not broad[4:].any()
    assert narrow[0] and not narrow[1]  # row 1 has a missing repeat value
    assert not (narrow & ~broad).any()


# ---------------------------------------------------------------- grouping and in-fold fitting
def test_grouped_split_never_mixes_a_cluster():
    df = _table(20, 100)
    y = (df.role == "class_positive").astype(int).to_numpy()
    g = df.cluster_id.to_numpy()
    for seed in range(5):
        for tr, te in fc.grouped_splits(y, g, 5, seed):
            assert not set(g[tr]) & set(g[te])
            assert y[te].sum() > 0


def test_imputer_and_scaler_fit_inside_training_rows_only():
    X = np.array([[1.0, 10.0], [2.0, np.nan], [3.0, 30.0], [4.0, 20.0], [5.0, 10.0], [6.0, 40.0]])
    y = np.array([0, 0, 1, 0, 1, 1])
    cols = ["frac_S", "frac_T"]
    train = np.array([0, 1, 2, 3])
    pipe = fc.build_pipeline("logreg", cols).fit(X[train], y[train])
    assert pipe.named_steps["imp"].statistics_[1] == pytest.approx(np.median([10.0, 30.0, 20.0]))
    assert pipe.named_steps["sc"].mean_[0] == pytest.approx(
        X[train, 0].mean()
    )  # test rows 5, 6 do not enter
    # a test row with a NaN is imputed with the training median
    out = pipe.predict_proba(np.array([[3.0, np.nan]]))
    assert out.shape == (1, 2) and np.isfinite(out).all()


def test_cv_one_seed_runs_with_missing_values_and_extras():
    df = _table(15, 75)
    df.loc[::7, "frac_ST"] = np.nan
    elig = np.arange(len(df)) < 80
    y = (df.role.to_numpy() == "class_positive").astype(int)
    X = df[fc.FEATURES].to_numpy(float)
    res = fc.cv_one_seed(
        "logreg",
        fc.FEATURES,
        0,
        X[elig],
        y[elig],
        df.cluster_id.to_numpy()[elig],
        X[~elig],
        df.cluster_id.to_numpy()[~elig],
        (0.05,),
    )
    assert np.isfinite(res["oof"]).all()
    assert fc.auc(res["oof"], y[elig]) > 0.9  # planted signal
    assert res["called"].shape == (1, 80)
    assert (res["ex_n"] > 0).all()  # extras are scored by models that did not see their cluster


# ---------------------------------------------------------------- statistics
def test_auc_and_bootstrap():
    y = np.array([0, 0, 0, 1, 1])
    assert fc.auc([0.1, 0.2, 0.3, 0.9, 0.8], y) == 1.0
    assert fc.auc([0.9, 0.9, 0.9, 0.9, 0.9], y) == 0.5
    g = np.array(list("aabbcdeeff")[:5])
    point, lo, hi = fc.cluster_bootstrap(lambda i: float(i.size), g, 20, 0)
    assert point == 5 and lo <= hi


# ---------------------------------------------------------------- apply CLI
def _write(df, path):
    df.to_csv(path, sep="\t", index=False)


def test_cli_apply_scores_refuses_overwrite_and_fails_on_missing_columns(tmp_path, capsys):
    train = _table(15, 75)
    new = _table(3, 3, seed=5)
    _write(train, tmp_path / "train.tsv")
    _write(new, tmp_path / "new.tsv")
    out = tmp_path / "scores.tsv"
    args = [
        "apply",
        "--train-features",
        str(tmp_path / "train.tsv"),
        "--new-features",
        str(tmp_path / "new.tsv"),
        "--model",
        "logreg",
        "--out",
        str(out),
    ]
    assert fc.main(args) == 0
    s = pd.read_csv(out, sep="\t")
    assert (
        list(s.columns) == ["accession", "score_logreg"]
        and len(s) == 6
        and s.score_logreg.between(0, 1).all()
    )
    assert fc.main(args) == 1  # refuses the existing output
    assert "FATAL:" in capsys.readouterr().err
    _write(new.drop(columns=["frac_ST"]), tmp_path / "bad.tsv")
    bad = [
        "apply",
        "--train-features",
        str(tmp_path / "train.tsv"),
        "--new-features",
        str(tmp_path / "bad.tsv"),
        "--out",
        str(tmp_path / "o2.tsv"),
    ]
    assert fc.main(bad) == 1
    err = capsys.readouterr().err
    assert "FATAL:" in err and "frac_ST" in err
    assert not (tmp_path / "o2.tsv").exists()


def test_score_new_never_fits_on_tier_b_or_excluded_roles():
    train = _table(15, 75)
    base = fc.score_new(train, train.iloc[:5], "logreg")
    extra = train.iloc[:5].copy()
    extra["role"] = "excluded_indirect_regulator"
    extra["pos_tier"] = ""
    both = pd.concat([train, extra.assign(accession=[f"X{i}" for i in range(5)], frac_ST=100.0)])
    again = fc.score_new(both, train.iloc[:5], "logreg")
    assert np.allclose(base.iloc[:, 1], again.iloc[:, 1])


# ---------------------------------------------------------------- committed record
def test_committed_readme_hashes_match():
    text = (DATA / "README.md").read_text()
    pairs = re.findall(r"^\| `([^`]+)` \| ([0-9a-f]{64}) \|", text, flags=re.M)
    assert len(pairs) >= 5
    for name, h in pairs:
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == h, name


def test_committed_features_are_numbers_only_and_complete():
    df = pd.read_csv(DATA / "features.tsv", sep="\t", keep_default_na=False, na_values=[""])
    assert list(df.columns) == fc.META + fc.FEATURES
    assert len(df) == 850 and df.accession.is_unique
    assert all(pd.api.types.is_numeric_dtype(df[c]) for c in fc.FEATURES)
    manifest = (ROOT / "analysis/aga1_class/manifest_accessions.txt").read_text().split()
    assert sorted(df.accession) == sorted(manifest)
    assert (df.role.eq("class_positive") & df.pos_tier.eq("A")).sum() == 35
    assert df.loc[df.gpi_missing == 1, "gpi_yes"].isna().all()
    assert df.loc[df.rep_missing == 1, "rep_coverage"].isna().all()
    assert not any("rep_unit" in c for c in df.columns)
