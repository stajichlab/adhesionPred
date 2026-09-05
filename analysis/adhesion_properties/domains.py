"""Functional domain flag computation for the protein-properties analysis."""

import re
import sys
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from db import fetch_domain_hits  # noqa: E402

_CAZY_FAM_RE = re.compile(r"^([^(]+)")
MIN_ADHESION_COUNT_FOR_ENRICHMENT = 5


def clean_cazy_family(cazyme_fam: str) -> str:
    """Strip the alignment-range suffix from a CAZy family string,
    e.g. 'GH114(59-261)' -> 'GH114'."""
    match = _CAZY_FAM_RE.match(cazyme_fam)
    return match.group(1) if match else cazyme_fam


def compute_domain_flags(con, protein_ids: Iterable[str]) -> pd.DataFrame:
    """Return one row per protein_id with boolean/categorical domain flags."""
    protein_ids = list(protein_ids)
    base = pd.DataFrame({"protein_id": protein_ids})

    pfam_hits = set(fetch_domain_hits(con, "pfam", protein_ids)["protein_id"])
    cazy_hits = set(fetch_domain_hits(con, "cazy_overview", protein_ids)["protein_id"])
    merops_hits = set(fetch_domain_hits(con, "merops", protein_ids)["protein_id"])
    signalp_hits = set(fetch_domain_hits(con, "signalp", protein_ids)["protein_id"])

    tmhmm = fetch_domain_hits(con, "tmhmm", protein_ids)
    tm_hits = set(tmhmm.loc[tmhmm["PredHel"] >= 1, "protein_id"]) if not tmhmm.empty else set()

    targetp = fetch_domain_hits(con, "targetp", protein_ids)
    targetp_map = (
        dict(zip(targetp["protein_id"], targetp["prediction"])) if not targetp.empty else {}
    )

    base["has_pfam"] = base["protein_id"].isin(pfam_hits)
    base["has_cazy"] = base["protein_id"].isin(cazy_hits)
    base["has_merops"] = base["protein_id"].isin(merops_hits)
    base["has_signal_peptide"] = base["protein_id"].isin(signalp_hits)
    base["has_tm_helix"] = base["protein_id"].isin(tm_hits)
    base["targetp_category"] = base["protein_id"].map(targetp_map).fillna("none")

    return base


def top_domain_table(
    con,
    table: str,
    id_field: str,
    adhesion_ids: Iterable[str],
    background_ids: Iterable[str],
    top_n: int = 20,
) -> pd.DataFrame:
    """Distinct-protein frequency of each domain identity among adhesion vs.
    background proteins, restricted to domains with at least
    MIN_ADHESION_COUNT_FOR_ENRICHMENT adhesion hits (avoids a single rare
    hit with zero background hits producing a meaningless infinite ratio),
    sorted by enrichment ratio (adhesion_rate / background_rate).
    """
    adhesion_ids = list(adhesion_ids)
    background_ids = list(background_ids)
    n_adhesion = len(adhesion_ids)
    n_background = len(background_ids)

    adhesion_hits = fetch_domain_hits(con, table, adhesion_ids)
    background_hits = fetch_domain_hits(con, table, background_ids)

    if table == "cazy_overview":
        if not adhesion_hits.empty:
            adhesion_hits = adhesion_hits.assign(
                **{id_field: adhesion_hits[id_field].map(clean_cazy_family)}
            )
        if not background_hits.empty:
            background_hits = background_hits.assign(
                **{id_field: background_hits[id_field].map(clean_cazy_family)}
            )

    def _distinct_counts(hits: pd.DataFrame) -> pd.Series:
        if hits.empty:
            return pd.Series(dtype=int)
        return (
            hits.drop_duplicates(["protein_id", id_field]).groupby(id_field)["protein_id"].nunique()
        )

    adhesion_counts = _distinct_counts(adhesion_hits)
    background_counts = _distinct_counts(background_hits)

    all_ids = sorted(set(adhesion_counts.index) | set(background_counts.index))
    rows = []
    for domain_id in all_ids:
        a_count = int(adhesion_counts.get(domain_id, 0))
        b_count = int(background_counts.get(domain_id, 0))
        a_rate = a_count / n_adhesion if n_adhesion else 0.0
        b_rate = b_count / n_background if n_background else 0.0
        enrichment = (a_rate / b_rate) if b_rate > 0 else float("inf")
        rows.append(
            {
                "domain_id": domain_id,
                "adhesion_count": a_count,
                "adhesion_rate": a_rate,
                "background_count": b_count,
                "background_rate": b_rate,
                "enrichment_ratio": enrichment,
            }
        )

    result = pd.DataFrame(rows)
    if result.empty:
        return result
    frequent = result[result["adhesion_count"] >= MIN_ADHESION_COUNT_FOR_ENRICHMENT]
    return frequent.sort_values("enrichment_ratio", ascending=False).head(top_n)
