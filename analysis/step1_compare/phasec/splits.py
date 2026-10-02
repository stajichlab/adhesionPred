"""Clusters, split schemes S1, S2, S3, the T-c leakage controls (Phase C spec 3.4 items 3-6).

Split ids:
- S1: homology-grouped 5-fold CV over the GO rows of the training sources (role `train`).
  StratifiedGroupKFold(5, shuffle, seed) on the pos and neg rows, groups = clusters. Every
  other row of the training sources and every T-c row takes the fold of its cluster; a cluster
  without a fold gets int(sha256(cluster_id)[:16], 16) mod 5. So T-c rows get out-of-fold
  scores too (spec 4, score_source).
- S2-<source>: leave-species-out. For each training source B: train on the other training
  sources, test B. For each `test_species` source: train on all training sources.
- S3-<clade>: leave-clade-out. Train on all training sources; test every `test_clade` or
  `undecided` source of the clade. S3-Eurotiomycetes also tests the literature rows.
- FULL: train on all training sources and all T-c rows; no test. Its models score the
  proteins that are in no training table (score_source `final`).

Parts: train (GO pos/neg training rows), train_tc (T-c rows, V-kw only), test (GO rows of the
test set, all classes), test_tc (S1: T-c rows of the test fold; scored, never truth), test_lit
(literature rows). T-c rows are removed from train_tc, first match wins:
a_test_protein (accession or hash of a test protein), b_cluster_mate (cluster holds a test
protein, ruling C-5), c_test_taxon (S2: taxon_id of the test species; S3: taxon in the test
clade, ruling C-7).

In S1 the guarantee of ruling C-5 comes from assigning T-c rows to the fold of their cluster, so
rules a and b are expected to remove few or no rows there (b_cluster_mate counts can be 0 in S1).

In S2 and S3 a GO training protein may share a cluster with a test protein: that is homology
across species, which the maximum-identity stratum measures (ruling C-4). Only S1 folds and
the T-c rows must not share clusters with the test set.
"""

import hashlib
from dataclasses import dataclass, field

import numpy as np
from sklearn.model_selection import StratifiedGroupKFold

S1_FOLDS = 5
LITERATURE_CLADE = "Eurotiomycetes"
PARTS = ("train", "train_tc", "test", "test_tc", "test_lit")
TRAIN_PARTS = ("train", "train_tc")
TEST_PARTS = ("test", "test_tc", "test_lit")
MEMBER_COLUMNS = ("split_id", "fold", "seq_sha256", "part", "origin", "class", "cluster_id")
REMOVED_COLUMNS = ("split_id", "fold", "seq_sha256", "gene_ids", "taxon_ids", "rule")
IDENTITY_COLUMNS = ("split_id", "seq_sha256", "max_identity", "below_0.3")
CLUSTER_COLUMNS = ("seq_sha256", "cluster_id")
IDENTITY_CUT = 0.3


class SplitError(ValueError):
    """A split breaks a leakage rule or the cluster table does not fit the sequences."""


@dataclass(frozen=True)
class SplitDef:
    split_id: str
    train_sources: tuple
    test_sources: tuple = ()
    test_taxa: frozenset = field(default_factory=frozenset)
    test_clades: frozenset = field(default_factory=frozenset)
    literature: bool = False
    folds: int = 1


def read_cluster_tsv(lines) -> dict[str, str]:
    """MMseqs2 `<prefix>_cluster.tsv` lines (representative TAB member) -> member: cluster."""
    out = {}
    for n, line in enumerate(lines, 1):
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            raise SplitError(f"cluster table line {n}: expected 2 fields, got {len(parts)}")
        rep, member = parts
        if member in out:
            raise SplitError(f"cluster table: {member} is in two clusters")
        out[member] = rep
    return out


def check_clusters(cluster_of: dict, hashes) -> None:
    hashes = set(hashes)
    missing = sorted(hashes - set(cluster_of))
    extra = sorted(set(cluster_of) - hashes)
    if missing or extra:
        raise SplitError(
            f"cluster table: {len(missing)} sequences without a cluster, {len(extra)} unknown "
            f"members (for example {(missing or extra)[0]})"
        )


def split_definitions(species_rows) -> list[SplitDef]:
    train = tuple(sorted(r["source_id"] for r in species_rows if r["role"] == "train"))
    taxon = {r["source_id"]: r["taxon_id"] for r in species_rows}
    defs = [SplitDef("S1", train, folds=S1_FOLDS)]
    if len(train) > 1:
        for b in train:
            rest = tuple(s for s in train if s != b)
            defs.append(SplitDef(f"S2-{b}", rest, (b,), frozenset({taxon[b]})))
    for r in sorted(species_rows, key=lambda r: r["source_id"]):
        if r["role"] == "test_species":
            s = r["source_id"]
            defs.append(SplitDef(f"S2-{s}", train, (s,), frozenset({taxon[s]})))
    clades: dict[str, list[str]] = {}
    for r in species_rows:
        if r["role"] in ("test_clade", "undecided"):
            clades.setdefault(r["in_clade"], []).append(r["source_id"])
    for clade in sorted(clades):
        defs.append(
            SplitDef(
                f"S3-{clade}",
                train,
                tuple(sorted(clades[clade])),
                test_clades=frozenset({clade}),
                literature=clade == LITERATURE_CLADE,
            )
        )
    defs.append(SplitDef("FULL", train))
    return defs


def _sources(row) -> set[str]:
    return set(row["source_ids"].split(","))


def _hash_fold(cluster_id: str, n_folds: int) -> int:
    return int(hashlib.sha256(cluster_id.encode()).hexdigest()[:16], 16) % n_folds


def s1_folds(table, cluster_of: dict, train_sources, seed: int, n_folds=S1_FOLDS):
    """Fold of every S1 row (GO rows of the training sources, T-c rows), by cluster."""
    train_sources = set(train_sources)
    fit = sorted(
        (r for r in table if r["origin"] == "go" and r["class"] in ("pos", "neg")
         and _sources(r) <= train_sources),
        key=lambda r: r["seq_sha256"],
    )  # fmt: skip
    y = np.array([r["class"] == "pos" for r in fit])
    groups = np.array([cluster_of[r["seq_sha256"]] for r in fit])
    cv = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    fold_of_cluster = {}
    for k, (_, test) in enumerate(cv.split(np.zeros(len(fit)), y, groups)):
        for c in groups[test]:
            fold_of_cluster[c] = k
    out = {}
    for r in table:
        in_s1 = r["origin"] == "tc" or (r["origin"] == "go" and _sources(r) <= train_sources)
        if not in_s1:
            continue
        c = cluster_of[r["seq_sha256"]]
        out[r["seq_sha256"]] = fold_of_cluster.get(c, _hash_fold(c, n_folds))
    return out


def tc_removals(tc_rows, test_rows, cluster_of: dict, taxa=frozenset(), clades=frozenset()):
    """Split T-c rows into (kept, [(row, rule)]) for one test set (rules a, b, c)."""
    test_acc = {a for r in test_rows for a in r["gene_ids"].split(",") if a}
    test_hash = {r["seq_sha256"] for r in test_rows}
    test_clusters = {cluster_of[r["seq_sha256"]] for r in test_rows}
    kept, removed = [], []
    for r in tc_rows:
        accs = set(r["gene_ids"].split(","))
        if accs & test_acc or r["seq_sha256"] in test_hash:
            removed.append((r, "a_test_protein"))
        elif cluster_of[r["seq_sha256"]] in test_clusters:
            removed.append((r, "b_cluster_mate"))
        elif set(r["taxon_ids"].split(",")) & taxa or set(r["clades"].split(",")) & clades:
            removed.append((r, "c_test_taxon"))
        else:
            kept.append(r)
    return kept, removed


def _member(split_id, fold, row, part, cls, cluster_of):
    return {
        "split_id": split_id,
        "fold": str(fold),
        "seq_sha256": row["seq_sha256"],
        "part": part,
        "origin": row.get("origin", "lit"),
        "class": cls,
        "cluster_id": cluster_of[row["seq_sha256"]],
    }


def _lit_as_rows(literature):
    rows = []
    for r in literature:
        if r["seq_sha256"]:
            rows.append(
                {
                    "seq_sha256": r["seq_sha256"],
                    "origin": "lit",
                    "class": "pos" if r["literature_positive"] == "yes" else "excluded",
                    "gene_ids": r["accession"],
                }
            )
    return rows


def build(table, literature, cluster_of: dict, defs, seed: int):
    """Return (members, removed) for every split definition."""
    go = [r for r in table if r["origin"] == "go"]
    tc = [r for r in table if r["origin"] == "tc"]
    lit = _lit_as_rows(literature)
    members, removed = [], []
    for d in defs:
        train_sources = set(d.train_sources)
        if d.split_id == "S1":
            fold = s1_folds(table, cluster_of, d.train_sources, seed, d.folds)
            s1_go = [r for r in go if _sources(r) <= train_sources]
            for f in range(d.folds):
                test = [r for r in s1_go if fold[r["seq_sha256"]] == f]
                train = [
                    r for r in s1_go if fold[r["seq_sha256"]] != f and r["class"] != "excluded"
                ]
                tc_test = [r for r in tc if fold[r["seq_sha256"]] == f]
                tc_kept, tc_out = tc_removals(
                    [r for r in tc if fold[r["seq_sha256"]] != f], test, cluster_of
                )
                members += _parts(d.split_id, f, train, tc_kept, test, tc_test, [], cluster_of)
                removed += _removed(d.split_id, f, tc_out)
            continue
        if d.split_id == "FULL":
            train = [r for r in go if _sources(r) <= train_sources and r["class"] != "excluded"]
            members += _parts("FULL", 0, train, tc, [], [], [], cluster_of)
            continue
        test_sources = set(d.test_sources)
        for r in go:
            if _sources(r) & test_sources and _sources(r) & train_sources:
                raise SplitError(
                    f"{d.split_id}: sequence {r['seq_sha256']} belongs to a training source and "
                    f"a test source ({r['source_ids']})"
                )
        train = [r for r in go if _sources(r) <= train_sources and r["class"] != "excluded"]
        test = [r for r in go if _sources(r) & test_sources]
        test_lit = lit if d.literature else []
        tc_kept, tc_out = tc_removals(tc, test + test_lit, cluster_of, d.test_taxa, d.test_clades)
        members += _parts(d.split_id, 0, train, tc_kept, test, [], test_lit, cluster_of)
        removed += _removed(d.split_id, 0, tc_out)
    return members, removed


def _parts(split_id, fold, train, train_tc, test, test_tc, test_lit, cluster_of):
    out = [_member(split_id, fold, r, "train", r["class"], cluster_of) for r in train]
    out += [_member(split_id, fold, r, "train_tc", "pos", cluster_of) for r in train_tc]
    out += [_member(split_id, fold, r, "test", r["class"], cluster_of) for r in test]
    out += [_member(split_id, fold, r, "test_tc", "pos", cluster_of) for r in test_tc]
    out += [_member(split_id, fold, r, "test_lit", r["class"], cluster_of) for r in test_lit]
    return out


def _removed(split_id, fold, tc_out):
    return [
        {
            "split_id": split_id,
            "fold": str(fold),
            "seq_sha256": r["seq_sha256"],
            "gene_ids": r["gene_ids"],
            "taxon_ids": r["taxon_ids"],
            "rule": rule,
        }
        for r, rule in tc_out
    ]


def _groups(members):
    out: dict[tuple, list] = {}
    for m in members:
        out.setdefault((m["split_id"], m["fold"]), []).append(m)
    return out


def check_test_truth(members) -> None:
    """T-c rows are never test truth (spec 2.3, Q4)."""
    bad = [m for m in members if m["origin"] == "tc" and m["part"] in ("test", "test_lit")]
    if bad:
        raise SplitError(f"{bad[0]['split_id']}: T-c row {bad[0]['seq_sha256']} is in a test set")


def check_no_shared_hash(members) -> None:
    for (split_id, fold), rows in _groups(members).items():
        train = {m["seq_sha256"] for m in rows if m["part"] in TRAIN_PARTS}
        test = {m["seq_sha256"] for m in rows if m["part"] in TEST_PARTS}
        both = train & test
        if both:
            raise SplitError(f"{split_id} fold {fold}: {sorted(both)[0]} is in train and test")


def check_no_cluster_spans(members) -> None:
    """S1: no cluster in both train and test of a fold. S2, S3: no T-c training row shares a
    cluster with the test set (ruling C-5)."""
    for (split_id, fold), rows in _groups(members).items():
        if split_id == "S1":
            train = {m["cluster_id"] for m in rows if m["part"] in TRAIN_PARTS}
        else:
            train = {m["cluster_id"] for m in rows if m["part"] == "train_tc"}
        test = {m["cluster_id"] for m in rows if m["part"] in TEST_PARTS}
        both = train & test
        if both:
            raise SplitError(
                f"{split_id} fold {fold}: cluster {sorted(both)[0]} spans train and test"
            )


def parse_hits(lines) -> dict[str, float]:
    """Highest fident per query from `query TAB target TAB fident` lines; self hits skipped."""
    best: dict[str, float] = {}
    for n, line in enumerate(lines, 1):
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 3:
            raise SplitError(f"search result line {n}: expected 3 fields, got {len(parts)}")
        query, target, fident = parts
        if query == target:
            continue
        value = float(fident)
        if value > best.get(query, -1.0):
            best[query] = value
    return best


def identity_rows(split_id: str, test_hashes, best: dict) -> list[dict]:
    rows = []
    for h in sorted(set(test_hashes)):
        value = best.get(h)
        rows.append(
            {
                "split_id": split_id,
                "seq_sha256": h,
                "max_identity": "" if value is None else f"{value:.4f}",
                "below_0.3": "yes" if value is None or value < IDENTITY_CUT else "no",
            }
        )
    return rows
