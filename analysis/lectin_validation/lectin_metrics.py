#!/usr/bin/env python3.12
"""Metric functions for the lectin validation (pure functions, no file access).

The independent unit is the cluster (MMseqs2, 30% identity, 50% coverage). A cluster counts as a
hit when at least one of its members in the group has a hit. Wilson 95% intervals are computed
on clusters. The sequence-level interval is given too; it ignores the cluster structure and
is too narrow.
"""

import math

Z95 = 1.959963984540054

VALID_LEAKAGE = {"none", "partial", "tuned_on_truth", "in_reference", "unknown"}
POSITIVE_CLASSES = {
    "structure_backed",
    "host_binding_h1a",
    "chitin_masking_h1b",
    "literature_binding",
    "known_gap",
    "pa14_owned_elsewhere",
}
NEGATIVE_CLASSES = {"er_housekeeping", "enzyme_lectin_like_domain", "n1_shared", "assumed"}


def wilson(k, n, z=Z95):
    """Wilson score interval for k successes in n trials. Returns (lo, hi); (nan, nan) if n == 0."""
    if n <= 0:
        return (float("nan"), float("nan"))
    if k < 0 or k > n:
        raise ValueError("k must be between 0 and n")
    p = k / n
    z2 = z * z
    denom = 1 + z2 / n
    centre = (p + z2 / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def hit_accessions(hit_rows, models, floors=None):
    """Accessions with a hit to any model in `models`.

    hit_rows: iterable of (accession, model, score). floors: {model: minimum score}. A model
    without a floor uses the GA (hits are already at GA when the table comes from hmmsearch --cut_ga).
    """
    floors = floors or {}
    models = set(models)
    out = set()
    for acc, model, score in hit_rows:
        if model in models and float(score) >= floors.get(model, float("-inf")):
            out.add(acc)
    return out


def summarize(members, hit_accs):
    """Counts and intervals for a group.

    members: list of (accession, cluster_id). hit_accs: set of accessions with a hit.
    Returns a dict: n_seq, k_seq, n_clu, k_clu, seq_lo, seq_hi, clu_lo, clu_hi (nan when n is 0).
    """
    n_seq = len(members)
    k_seq = sum(1 for a, _ in members if a in hit_accs)
    clusters = {}
    for a, c in members:
        clusters.setdefault(c, False)
        if a in hit_accs:
            clusters[c] = True
    n_clu = len(clusters)
    k_clu = sum(1 for v in clusters.values() if v)
    sl, sh = wilson(k_seq, n_seq)
    cl, ch = wilson(k_clu, n_clu)
    return {
        "n_seq": n_seq,
        "k_seq": k_seq,
        "n_clu": n_clu,
        "k_clu": k_clu,
        "seq_lo": sl,
        "seq_hi": sh,
        "clu_lo": cl,
        "clu_hi": ch,
    }


def fmt(k, n, lo, hi):
    """'k/n (lo-hi)' with two decimals, or 'n=0'."""
    if n == 0:
        return "n=0"
    return f"{k}/{n} ({lo:.2f}-{hi:.2f})"


def meets_floor(n_pos_clusters, n_neg_clusters, floor=20):
    """True when both classes have at least `floor` independent clusters (COMMON-RULES rule 7, first part)."""
    return n_pos_clusters >= floor and n_neg_clusters >= floor


def cluster_seed_split(members, flagged):
    """Split a group at the cluster level by a flag (for example "in_reference").

    members: list of (accession, cluster_id). flagged: set of accessions with the flag.
    A cluster is `similar` when at least one of its members in the group is flagged. Every
    member of a `similar` cluster goes to the `similar` list, flagged or not. The `independent`
    list holds the members of clusters with no flagged member. A cluster is in exactly one
    list, so the cluster counts of the two lists add up to the number of clusters of the group.
    Returns (similar, independent), both lists of (accession, cluster_id).
    """
    similar_clusters = {c for a, c in members if a in flagged}
    similar = [(a, c) for a, c in members if c in similar_clusters]
    independent = [(a, c) for a, c in members if c not in similar_clusters]
    return similar, independent


def identity_bin(pident):
    """Bin for the best identity to a seed sequence (None = no seed hit at 30% or more)."""
    if pident is None:
        return "no seed hit (<30%)"
    if pident >= 98.0:
        return "98% or more"
    if pident >= 90.0:
        return "90 to 97.9%"
    if pident >= 50.0:
        return "50 to 89.9%"
    return "30 to 49.9%"
