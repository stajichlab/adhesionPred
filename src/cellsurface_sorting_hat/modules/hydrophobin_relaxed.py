"""Relaxed hydrophobin level: the seven hydrophobin-class Pfam models with a frozen full-sequence cutoff -> module ``hydrophobin_relaxed``.

The shipped ``hydrophobin_relaxed.hmm`` holds the models with the GA line rewritten to the frozen full-sequence cutoff and a domain
cutoff of -1000, so ``hmmsearch --cut_ga`` decides a hit by the full-sequence score. ``hit`` is 1 when the protein has a row in the
domain table, at least 8 cysteines in the full sequence, and the R0 call is ``called``. It is 0 when any of these fails. It is empty
(not assessable) when the protein would be a hit except that its R0 row is missing or not usable.
"""

import re
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row, sha256_file

COLUMNS = ["hit", "score", "n_cys", "query"]
MIN_CYS = 8
DOMAIN_GA = -1000.0


def model_info(path):
    """Model names with lengths, and the single GA pair of the file. Refuses models with different GA lines."""
    text = Path(path).read_text()
    names, ga = {}, set()
    for block in text.split("//"):
        n = re.search(r"^NAME\s+(\S+)", block, re.M)
        if not n:
            continue
        ln = re.search(r"^LENG\s+(\d+)", block, re.M)
        g = re.search(r"^GA\s+(\S+)\s+(\S+);", block, re.M)
        if not ln or not g:
            raise ValueError(f"{path}: model {n.group(1)} has no LENG or GA line")
        names[n.group(1)] = int(ln.group(1))
        ga.add((float(g.group(1)), float(g.group(2))))
    if not names:
        raise ValueError(f"{path}: no models")
    if len(ga) != 1:
        raise ValueError(f"{path}: the models have different GA lines: {sorted(ga)}")
    return {"names": names, "ga": next(iter(ga)), "path": str(path), "sha256": sha256_file(path)}


def check_table(path, info, search_option):
    """The domain table must come from this file, with --cut_ga and the stated search option. Raises ValueError otherwise."""
    from cellsurface_sorting_hat.modules import pfam

    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    head = [ln for ln in lines if ln.startswith("#")]
    opt_line = next((ln for ln in head if "Option settings" in ln), "")
    if "--cut_ga" not in opt_line:
        raise ValueError(f"{path}: hmmsearch was not run with --cut_ga")
    has_nobias = "--nobias" in opt_line
    if has_nobias != (search_option == "nobias"):
        raise ValueError(
            f"{path}: search option is {search_option!r} but the run {'used' if has_nobias else 'did not use'} --nobias"
        )
    q = next((ln for ln in head if ln.startswith("# Query file:")), "")
    if Path(q.split(":", 1)[1].strip()).name != Path(info["path"]).name:
        raise ValueError(
            f"{path}: query file {q.split(':', 1)[-1].strip()!r} is not {Path(info['path']).name!r}"
        )
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        if f[3] not in info["names"]:
            raise ValueError(f"{path}: model {f[3]!r} is not in {Path(info['path']).name}")
        if int(f[5]) != info["names"][f[3]]:
            raise ValueError(
                f"{path}: model {f[3]} length {f[5]} differs from {info['names'][f[3]]} in the file"
            )
    return pfam.parse_domtblout(path)


def relaxed_rows(proteins, hits, sp_calls, info):
    best = {}
    for h in hits:
        t = h["target"]
        if t not in best or h["seq_score"] > best[t]["seq_score"]:
            best[t] = h
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        n_cys = p.sequence.upper().count("C")
        h = best.get(p.id)
        if h is None or n_cys < MIN_CYS:
            hit = "0"
        else:
            sp = sp_calls.get(p.id)
            hit = "" if sp is None else ("1" if sp == "called" else "0")
        rows.append(
            {
                "id": p.id,
                "state": "ok",
                "hit": hit,
                "score": "" if h is None else f"{h['seq_score']}",
                "n_cys": str(n_cys),
                "query": "" if h is None else h["query"],
            }
        )
    return rows


def module_params(info, search_option, condition):
    return {
        "cutoff_seq_bits": info["ga"][0],
        "domain_ga": info["ga"][1],
        "score_type": "full-sequence bit score (--cut_ga with the rewritten GA line)",
        "search_option": search_option,
        "min_cysteines": MIN_CYS,
        "models": sorted(info["names"]),
        "conditions": {"sp_module": condition},
    }
