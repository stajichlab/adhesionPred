import json
import re
import sys
import urllib.parse
import urllib.request

from lxml import etree


def get(url):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "curation/1.0"}), timeout=60
    ).read()


def pmid2pmc(pmid):
    j = json.loads(
        get(
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json&resultType=core"
        )
    )
    r = j["resultList"]["result"]
    return r[0] if r else None


def fulltext(pmc):
    try:
        x = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML")
    except Exception:
        return None
    return x if len(x) > 1000 else None


def xml2txt(x):
    root = etree.fromstring(x)
    out = []
    for el in root.iter("p", "title", "article-title", "abstract", "caption", "td", "th", "label"):
        t = "".join(el.itertext())
        out.append(re.sub(r"\s+", " ", t).strip())
    return "\n".join(out)


if __name__ == "__main__":
    for pmid in sys.argv[1:]:
        r = pmid2pmc(pmid)
        if not r:
            print(pmid, "NOT FOUND")
            continue
        pmc = r.get("pmcid")
        print(pmid, pmc, r.get("title"), r.get("journalTitle"), r.get("pubYear"))
        ab = r.get("abstractText", "")
        ab = re.sub(r"<[^>]+>", "", ab) if False else ab
        # keep abstract with tags stripped (no space insertion)
        abx = re.sub(r"<[^>]+>", "", ab)
        open(f"texts/PMID{pmid}_abstract.txt", "w").write(f"{r.get('title')}\n{abx}\n")
        if pmc:
            x = fulltext(pmc)
            if x:
                open(f"texts/{pmc}.xml", "wb").write(x)
                open(f"texts/{pmc}.txt", "w").write(xml2txt(x))
                print("  fulltext saved", pmc, len(x))
            else:
                print("  no EuropePMC full text", pmc)


def efetch(pmc):
    x = get(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmc}&rettype=xml"
    )
    if b"<body" not in x:
        return None
    open(f"texts/{pmc}.xml", "wb").write(x)
    open(f"texts/{pmc}.txt", "w").write(xml2txt(x))
    return len(x)


def html2txt(h):
    import html as _h

    h = re.sub(r"<(script|style)[^>]*>.*?</\1>", "", h, flags=re.S)
    h = re.sub(r"</?(p|div|h\d|li|tr|br|section|table|figcaption|caption)[^>]*>", "\n", h)
    h = re.sub(r"</t[dh]>", " | ", h)
    h = re.sub(r"<[^>]+>", "", h)
    h = _h.unescape(h)
    return "\n".join(re.sub(r"[ \t\xa0]+", " ", ln).strip() for ln in h.split("\n") if ln.strip())
