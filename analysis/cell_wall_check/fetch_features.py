#!/usr/bin/env python3.12
"""Fetch UniProt feature annotations (signal peptide, GPI anchor, TM, propeptide) for the accessions in test_set.tsv.
Output: uniprot_features.tsv. Data is untrusted text; only feature type/begin/end are parsed."""

import csv
import json
import urllib.request

rows = list(csv.DictReader(open("test_set.tsv"), delimiter="\t"))
with open("uniprot_features.tsv", "w") as out:
    out.write("gene\tuniprot\tlength\tsignal\tgpi_site\tn_tm\ttm_ranges\tpropeptide_or_chain\n")
    for r in rows:
        d = json.load(
            urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{r['uniprot']}.json")
        )
        sig, gpi, tm, ch = [], [], [], []
        for f in d.get("features", []):
            t = f["type"]
            b = f["location"]["start"]["value"]
            e = f["location"]["end"]["value"]
            if t == "Signal":
                sig.append(f"{b}-{e}")
            elif t == "Lipidation" and "GPI" in f.get("description", ""):
                gpi.append(f"{b}")
            elif t == "Transmembrane":
                tm.append(f"{b}-{e}")
            elif t in ("Propeptide", "Chain"):
                ch.append(f"{t}:{b}-{e}")
        out.write(
            "\t".join(
                [
                    r["gene"],
                    r["uniprot"],
                    r["length"],
                    ",".join(sig) or "-",
                    ",".join(gpi) or "-",
                    str(len(tm)),
                    ",".join(tm) or "-",
                    ";".join(ch) or "-",
                ]
            )
            + "\n"
        )
