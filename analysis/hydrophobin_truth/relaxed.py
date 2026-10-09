"""Pure functions for the relaxed Pfam level and the cutoff rule of the hydrophobin spec (section 6.1).

No positive protein's score is an input to lowest_cutoff. Positive scores are used only by choose_option, to choose
between the two search options (default filters, --nobias).
"""

FLOOR = 0.0


def parse_tblout(lines):
    """Best full-sequence score per target from hmmsearch --tblout lines (the score is column 6)."""
    best = {}
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        c = line.split()
        score = float(c[5])
        if c[0] not in best or score > best[c[0]]:
            best[c[0]] = score
    return best


def conditioned_extra(items):
    """Ids that count as extra calls: R0 called, at least 8 cysteines, not called by the strict level, not labelled."""
    return sorted(
        i
        for i, x in items.items()
        if x["r0"] == "called" and x["n_cys"] >= 8 and not x["strict"] and not x["labelled"]
    )


def _ok(c, cands, sizes, hn_scores, per10k, hn_rate, ignore):
    by = {}
    for x in cands:
        if x["score"] >= c:
            by[x["proteome"]] = by.get(x["proteome"], 0) + 1
    for p, n in by.items():
        if n > per10k * sizes[p] / 10000.0:
            return False
    for g, scores in hn_scores.items():
        if g in ignore or not scores:
            continue
        if sum(1 for s in scores if s >= c) / len(scores) > hn_rate:
            return False
    return True


def lowest_cutoff(cands, sizes, hn_scores, per10k, hn_rate, floor=FLOOR, ignore_groups=()):
    """Lowest bit score at or above the floor at which (i) each tuning proteome has at most per10k extra calls per
    10,000 proteins and (ii) each hard-negative group (except ignored ones) has a call rate at most hn_rate.
    Ties go to the higher cutoff: when several cutoffs give the same call set, the lowest included score is returned.
    Returns the floor when no tuning score is at or above it. cands: dicts with id, proteome, score."""
    ignore = set(ignore_groups)
    scores = sorted(
        {x["score"] for x in cands if x["score"] >= floor}
        | {s for g, v in hn_scores.items() if g not in ignore for s in v if s >= floor}
    )
    for c in [floor] + scores:
        if _ok(c, cands, sizes, hn_scores, per10k, hn_rate, ignore):
            included = [s for s in scores if s >= c]
            return included[0] if included else floor
    return scores[-1] + 0.1 if scores else floor


def choose_option(recovered):
    """recovered: {option: set of recovered Pfam-missed clusters}. More clusters wins; a tie goes to default filters."""
    best = max(recovered.values(), key=len)
    for opt in ("default", "nobias"):
        if opt in recovered and len(recovered[opt]) == len(best):
            return opt
