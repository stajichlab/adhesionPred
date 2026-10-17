# ruff: noqa
"""Fetch PubMed abstracts (title + abstract) via E-utilities; cache in dl/abs/PMID.txt; print them.
Usage: abs.py PMID [PMID ...]"""

import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dl", "abs")
os.makedirs(D, exist_ok=True)
need = [p for p in sys.argv[1:] if not os.path.exists(f"{D}/{p}.txt")]
for i in range(0, len(need), 50):
    ids = ",".join(need[i : i + 50])
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id={ids}"
    )
    x = subprocess.run(["curl", "-s", "-m", "120", url], capture_output=True, text=True).stdout
    root = ET.fromstring(x)
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//PMID")
        title = (
            "".join(art.find(".//ArticleTitle").itertext())
            if art.find(".//ArticleTitle") is not None
            else ""
        )
        ab = " ".join("".join(a.itertext()) for a in art.findall(".//AbstractText"))
        jr = art.findtext(".//Journal/ISOAbbreviation") or ""
        yr = art.findtext(".//PubDate/Year") or ""
        open(f"{D}/{pmid}.txt", "w").write(f"{title}\n[{jr} {yr}]\n{ab}\n")
    time.sleep(0.4)
for p in sys.argv[1:]:
    f = f"{D}/{p}.txt"
    print(f"=== {p}")
    print(open(f).read() if os.path.exists(f) else "NOT FOUND")
