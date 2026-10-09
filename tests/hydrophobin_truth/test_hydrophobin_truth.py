import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/build_truth.py"
spec = importlib.util.spec_from_file_location("build_truth", SCRIPT)
bt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bt)


def entry(name, comments=None, reviewed=True, acc="P00001", seq="MKC"):
    return {
        "entryType": "UniProtKB reviewed (Swiss-Prot)"
        if reviewed
        else "UniProtKB unreviewed (TrEMBL)",
        "primaryAccession": acc,
        "uniProtkbId": "X_TEST",
        "organism": {"scientificName": "Test sp", "taxonId": 1},
        "proteinDescription": {"recommendedName": {"fullName": {"value": name}}},
        "comments": comments or [],
        "sequence": {"value": seq},
    }


def comment(kind, code, source="PubMed", pmid="123"):
    ev = {"evidenceCode": code}
    if source:
        ev.update({"source": source, "id": pmid})
    return {"commentType": kind, "texts": [{"value": "text", "evidences": [ev]}]}


def test_name_match_requires_hydrophobin_or_rodlet():
    assert bt.is_named(entry("Class I hydrophobin rodA"))
    assert bt.is_named(entry("Rodlet protein A"))
    assert not bt.is_named(entry("C2H2 type master regulator brlA"))
    assert not bt.is_named(entry("ATP synthase subunit 9"))


def test_t2_needs_experimental_pubmed_evidence_on_function():
    e = entry("Class I hydrophobin rodA", [comment("FUNCTION", "ECO:0000269")])
    assert bt.assign_tier(e) == ("T2", ["123"])


def test_t3_when_only_rule_or_similarity_evidence():
    e = entry("Class I hydrophobin 2", [comment("FUNCTION", "ECO:0000255", source=None)])
    assert bt.assign_tier(e)[0] == "T3"
    e = entry("Class I hydrophobin 2", [comment("SIMILARITY", "ECO:0000305", source=None)])
    assert bt.assign_tier(e)[0] == "T3"
    e = entry(
        "Class I hydrophobin 2", [comment("FUNCTION", "ECO:0000250", source="UniProtKB", pmid="Q1")]
    )
    assert bt.assign_tier(e)[0] == "T3"


def test_experimental_evidence_on_unrelated_comment_type_does_not_make_t2():
    e = entry("Class I hydrophobin x", [comment("INDUCTION", "ECO:0000269")])
    assert bt.assign_tier(e)[0] == "T3"


def test_unnamed_entry_is_not_a_positive():
    e = entry("C2H2 type master regulator brlA", [comment("FUNCTION", "ECO:0000269")])
    assert bt.assign_tier(e) == (None, [])


def test_unreviewed_entry_is_never_t2():
    e = entry("Class I hydrophobin", [comment("FUNCTION", "ECO:0000269")], reviewed=False)
    assert bt.assign_tier(e)[0] == "T3"


def test_split_keeps_clusters_whole_and_is_reproducible():
    clusters = {f"p{i}": f"c{i // 3}" for i in range(30)}
    a = bt.split_clusters(clusters, seed=7, dev_fraction=0.3)
    b = bt.split_clusters(clusters, seed=7, dev_fraction=0.3)
    assert a == b
    for c in set(clusters.values()):
        members = [p for p, cl in clusters.items() if cl == c]
        assert len({a[p] for p in members}) == 1
    assert set(a.values()) == {"dev", "test"}
